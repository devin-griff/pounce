# Covariance and information: design notes

This document records the design of pyomo-pounce's parameter
uncertainty subsystem as it stands: `covariance()` and
`information()`, the declaration surface that feeds them, the `wrt=`
block selection both accessors take, and `retain_kkt()`. Everything
here is computed from ONE ordinary solve. The solver factorizes the
KKT matrix to solve the NLP; the subsystem keeps that factorization
and answers every uncertainty question against it, so no second
solve, no finite differencing, and no perturbed re-solves appear
anywhere below.

User-facing documentation lives in `docs/src/sensitivity.md`; demo
notebooks 31 (information and identifiability) and 32 (one solve,
many questions) walk the surface end to end. This file records the
constructions and the reasons they are what they are.

## The declaration surface

`declare_fitted(vars)` flags free variables as the fitted parameters
of a least-squares problem. `declare_residual(container, group=)`
flags indexed variables holding the residuals, one member per data
point; the residual count and the SSR are derived from them, so no
data counts are passed. Groups partition residuals into noise
populations: containers sharing a group pool into one estimated
variance, distinct groups get their own and switch the covariance to
the heteroscedastic sandwich form. `retain_kkt(model)` keeps the
factorization with nothing declared at all, for problems where no
subset is THE fitted set and every question arrives as an explicit
`wrt=` block.

The retention policy in one place: the factorization is kept if
anything is declared or `retain_kkt()` was called; a `Covariance` or
`Information` result whose `conditioned_on` has not been read keeps
the session, and with it the factorization, alive until first
access; and an undeclared solve without the call takes the ordinary
subprocess path and pays nothing.

Explicit call-time forms (`solve(m, fitted=..., residuals=...)`)
mirror the declarations and are deliberately solve-local: they are
not written back into the registry, so repeated solves of one model
(the receding-horizon pattern) do not accumulate duplicate
declarations. The registry is deepcopy-aware, so `model.clone()`
carries declarations and the retain flag while the session, which
holds handles into one converged factorization, stays behind.

## Two spaces of rows

Every variable index the session hands around is user-space, the
`.col` file order over ALL of the model's variables (full-x). The
factor's x block is the algorithm's space (var-x): a variable whose
bounds are equal is removed from the solve under the default
`fixed_variable_treatment = make_parameter`, and every later column
moves up. The two spaces coincide exactly when the model has no
fixed variable, which makes confusing them invisible until it is
not: indexing the factor with a full-x row silently reads a
NEIGHBORING variable's numbers. `Solver.primal_rows` translates
full-x to factor rows (with `None` for removed variables, which the
accessors refuse with an actionable message), and
`session.scatter_x` translates factor-space vectors back. Every
factor index in the subsystem routes through this pair; the activity
report and `row_normal` read in full-x, and the two spaces are kept
in separately named variables (`rows` vs `krows`) wherever both
appear.

## The block and its covariance

`wrt=` selects the block: a Var (scalar or indexed, every member),
an indexed slice, a `(Var, iterable)` pair (a tuple of two Vars is a
two-member block, not a pair), data objects, or a list mixing these,
normalized to an ordered duplicate-free list. Omitted, the block is
the declared fitted set, and the accessors behave exactly as they
did before `wrt=` existed.

For block rows `krows`, unit backsolves against the held factor give
the K-inverse columns, and

    M[i, j] = (K^-1)[krows[i], krows[j]]

is the block of the inverse KKT matrix, symmetrized. `kkt_solve`
returns natural units (the backsolver unscales when NLP scaling was
active), so M and everything built from it composes without scale
factors. The homoscedastic Lagrangian covariance is `2 sigma^2 M`
with the membership corrections below; the factor 2 comes from the
objective being a plain sum of squares rather than the statistician's
half.

Sigma is a property of the FIT, never of the block: estimated from
declared residuals as `SSR / (n - n_fit)` per group, taken from
`sigma_sq=` when known (a float, or a per-group dict), or derived
from `n_data=` against the solve-time objective value (recorded at
the solve, so writing measurements into the model afterwards cannot
silently rescale the covariance). A sub-block's marginal therefore
equals the corresponding entries of the default answer exactly.

Marginal is the operative word: everything outside the block adjusts
rather than being held fixed, because the K-inverse block IS the
marginal covariance. One solve serves as many blocks as are asked
about, each call re-reducing onto its own argument.

## Gauss-Newton

`hessian="gauss-newton"` replaces the exact reduced Hessian with the
expected information `2 J^T J`, where J is the residual Jacobian
with respect to the block. J is recovered exactly from the same
backsolves by an identity: the residual rows of the K-inverse
columns are `J M`, so `Z_r inv(M) = J` and the factor's barrier
weight cancels regardless of its size. The product is formed over
the whole block and sliced afterwards, so pinned members still have
rows when their disposition is assembled. Gauss-Newton is
structurally positive semidefinite, which is what MHE arrival costs
and scipy-parity consumers need; the Lagrangian default is the
observed information, exact at the solution. The two agree for
linear models and in the small-residual limit and differ by
O(residual x curvature) otherwise. The heteroscedastic sandwich uses
the same recovered per-group Jacobians in both Hessian modes.

## Tangent recovery: the exact reduced Hessian

`information()` returns the reduced Hessian over the block, natural
units, no sigma anywhere, so for a homoscedastic fit covariance
equals `2 sigma^2 inv(information)` on the free block. The
Lagrangian form is built by tangent recovery. The x blocks of the
K-inverse columns are `T M`: each column satisfies the linearized
equalities and has the block as its own coordinates, so

    T = Zx inv(M),    R = T^T H T,

with H applied through `Solver.hessian_vec` (the exact Lagrangian
Hessian times a user-space vector, natural units, one product per
block column; factor-space tangents are scattered to full-x first).
The factor's barrier weight W = H + Sigma cancels multiplicatively
inside the recovery instead of being subtracted off.

Accuracy boundary, measured: machine-exact for equality and
variable-bound activity, including pinned variables at
`Sigma/q ~ 3e10` where a subtraction loses ten digits. A binding
INEQUALITY row is the exception: it couples through its slack
barrier with a large but finite weight, leaving roughly 1e-6
relative residue at practical barrier parameters and degrading as mu
tightens, because the pinned combination drives M toward
singularity. The construction requires the square estimation
structure: the equality constraints determine the non-fitted
variables given the block, checked as
`n_var - n_eq == n_params` before any tangent is formed.

Blocks that are not the fitted set route as follows. A block that
parameterizes the constraint manifold (size equal to the degrees of
freedom) gets the direct tangent. A proper sub-block of the fitted
set gets its marginal information by Schur complement of the exact
tangent R over the fitted block: free fitted variables outside the
block are profiled out, pinned ones are conditioned on (rows
dropped), and no covariance is ever inverted, so a pinned member
costs no digits. Any other within-count block reduces off the held
factor with the classification corrections, which is benign for free
coordinates since their slice carries no barrier term. Fitted-level
binding rows decline the Schur route with a warning, since their
projection does not compose simply with marginalization.

## Activity classification

The Rust core classifies every bounded variable and inequality row
of the converged solve. With `Sigma = z/s` summed over whichever
bound sides exist and q the curvature reference, the ratio
`r = Sigma/q` is compared against mu-dependent edges: inactive below
`sqrt(mu)`, strongly active above `1/sqrt(mu)`, weakly active inside
the fixed band `[1e-1, 1e1]`, ambiguous in the gaps. Above
`mu = 1e-4` the edges have thinned into the band and only the two
clear calls are made. For a variable, q is the raw Lagrangian
diagonal `|H_ii|`; for a row it is `|grad_d^T H grad_d| / ||grad_d||^4`,
the fourth power making r invariant to row rescaling and absorbing
the solver's per-row scale. A q below the floor
`sqrt(eps) * max(1, max_j |H_jj|)` is unidentified. An inactive
variable with `Sigma > 100 mu` is flagged contaminated, the
mu-relative threshold being the only one an inactive entry can
meaningfully exceed. The report is user-space indexed, and its
exports (`var_sigma`, `row_sigma`, `row_normal(j)`,
`hessian_vec(v)`) follow the natural-units contract; classification
itself runs on the solver's scaled quantities, where the ratio is
invariant.

Block members are classified at the REDUCED level, because in the
residual idiom a fitted parameter has zero raw curvature: the
effective curvature is `q_red = |diag(inv(M)) - Sigma|`, clamped to
a cancellation floor rather than refused (a huge Sigma cancelling
inside q_red would otherwise misfile a strongly active entry), and
the same ratio edges make the call. Variables outside the block are
classified per candidate as a singleton block by the identical rule,
one backsolve giving `(K^-1)_ii`, behind a cheap `Sigma > sqrt(mu)`
prefilter so only near-bound variables pay. This is scale-invariant
where any absolute Sigma threshold is not.

## Dispositions

The two accessors share one classification pass, so membership, row
handling, and their warnings cannot drift. For a block member:

- Free (inactive or unbounded): an ordinary row of the answer. A
  weakly active member is KEPT free at its true variance, its own
  barrier weight subtracted from the reduced block so it reports the
  curvature q rather than the factor's 2q; the warning states that
  boundary asymptotics are nonstandard.
- Strongly active (pinned): `covariance()` embeds a zero row, the
  variance conditional on the bound. `information()` returns S, the
  Schur reduction onto the pinned set, because zero variance is not
  zero information; cross blocks between free and pinned are zero,
  and S is conditional on the rest of the pinned set.

A strongly active general row whose normal is supported on the block
pins a DIRECTION rather than a coordinate: the free block is reduced
on the null space of the binding normals and pushed back, singular
by the number of binding rows, identically in both accessors, and
the projection annihilates the row's barrier weight exactly. The
conditional information along the pinned combination is reported in
the warning, computed by the tangent route inside the square
structure and by factor subtraction outside it. A row whose support
leaves the block cannot be represented by a restricted normal, so it
is warned about and left unprojected, with the raw classification as
its honest status. An indefinite Lagrangian information block is
returned as computed with a warning naming Gauss-Newton as the
positive semidefinite alternative: refusing would withhold the
finding that the point is not a minimum or the model is
over-parameterized.

Diagnostics name what they touch: on the default block they speak of
fitted parameters, under an explicit `wrt=` of block members and
variables outside the block, and every message is prefixed by the
accessor that produced it.

## Rank and singularity guards

Whether LAPACK raises on a structurally singular system is
build-dependent: the same matrix raises on one BLAS and returns
garbage on another, so no structural condition in the subsystem is
guarded by catching `LinAlgError` alone. Explicit blocks are gated
by count (more coordinates than the fit has degrees of freedom) and
then by a rank test on M (fp-detectable dependence, a duplicated
design point being the canonical case), each path with its own
message. A rank-deficient block is the prediction-band case:
`covariance()` returns the homoscedastic Lagrangian marginal
`2 sigma^2 M` with membership handling bypassed, exactly the ribbon
around a fitted trajectory, while `information()` raises toward
`covariance()`, per-group noise and Gauss-Newton raise because their
profiled Jacobians need `inv(M)`, and the singular free block inside
the S computation is rank-gated the same way. `_SingularBlock`, a
dedicated exception from the shared inversion helper, stands as the
last resort behind all of these so no rescue path does control flow
on message text. A held factor carrying inertia-correction
perturbations warns that the answer is regularized rather than
exact, since the isotropic delta_w lands on the free block and
survives projection.

## The reporting surface

`Covariance` carries the matrix keyed by the block's data objects
(either index order), `std_err`, `correlation` (entries with exactly
zero variance report correlation 0), `sigma_sq` as used, and
`eigen()`. `Information` carries the matrix, `params`, and
`eigen()`, eigenvalues ascending, for identifiability diagnosis: a
near-zero eigenvalue is a direction the data does not determine and
its eigenvector names the combination. Both carry `conditioned_on`,
the strongly active variables OUTSIDE the block: their barrier
weight stays in the held factor and drives the coupling through them
to zero as mu falls, so the block's numbers are conditional on those
bounds rather than marginal over them, and the conditioning travels
with the answer. The list is computed on first access and cached, so
calls that never read it pay nothing; until then the pending
computation keeps the session alive.

## Validation

The suites anchor to closed-form values rather than to the
implementation: the linear model's information is exactly `2 X^T X`,
restricted least squares gives the pinned dispositions, the hat
matrix gives the prediction band, and the inverse identity ties the
two accessors together. Load-bearing constructions are
mutation-verified (a broken tangent, a slice-first Gauss-Newton, a
block-sized degrees of freedom, a re-widened tuple guard each fail a
named test). Two fixture axes are always exercised because their
absence has repeatedly hidden bugs: objective scaling engaged with
`df != 1` asserted, and an inert fixed variable ahead of the block
in `.col` order. Structural refusals are tested through fixtures
that reach them deterministically, bit-identical coordinates rather
than near-singularity, so the tests do not inherit LAPACK's
build-dependence.
