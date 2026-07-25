# Covariance and information: v0.10 roadmap

**Status: roadmap proposal for pyomo-pounce, targeting v0.10.** This note
scopes the post-solve second-order surface in pyomo-pounce: `covariance()`,
shipping in v0.9, and the additions v0.10 should make around it. Nothing here
changes an existing signature. Companion to the active-set sensitivity roadmap
(`sensitivity-roadmap.md`), which extends `estimate()`; this note extends
`covariance()`.

## State of the art

Parameter covariance from the reduced Hessian is standard: at an estimation
optimum the covariance is the scaled inverse of the reduced Hessian of the
Lagrangian. The reduced Hessian is the object the sensitivity engines form
first. sIPOPT and k_aug both compute it, the information matrix, and hand it
back directly: k_aug writes the reduced Hessian to a file, and sIPOPT exposes
it through `rh_eigendecomp`. The parameter covariance is a downstream step on
top of it, invert and scale. The covariance-first tools sit a level up, at the
estimation frontend: Pyomo's `parmest` inverts the reduced Hessian to report a
covariance, and scipy's `curve_fit` returns the covariance (`pcov`) in the
Gauss-Newton form.

## Benefit hypothesis

The reduced Hessian and the covariance recipe are established, and pounce ships
them in its core, QP, and `curve_fit` interfaces (see Related reduced-Hessian
work below). What no pounce interface offers together, on a Pyomo model, is the
un-inverted object and a reduce-onto block chosen after the solve rather than
at declaration time. The next section states each against the shipped code.

Moving horizon estimation wants both at once, and is the worked case: its
information-form arrival cost is `Γ(x0) = ½ (x0 − x̂)ᵀ Π⁻¹ (x0 − x̂)`,
where `Π⁻¹` is the reduced Hessian of the one-step subproblem (the previous
arrival cost plus the stage leaving the window) marginalized onto the arrival
state. That is Lagrangian information, not covariance, on a block the solve did
not know about in advance.

## Where we are

v0.9 ships `covariance()` (in `pyomo-pounce/pyomo_pounce/sens.py`). You
`declare_fitted` a set of free variables, solve, and `covariance(model)`
returns their asymptotic covariance: the scaled parameter block of the
inverse KKT matrix, $2\sigma^2 (K^{-1})_{pp}$, where $\sigma^2$ is the residual
variance and the 2 is the sum-of-squares convention that makes the reduced
Hessian $2J^\top J$. Under the hood that block
is `inv(d2f*/dp2)`, the inverse reduced Hessian of the eliminated problem,
obtained by one backsolve per declared variable against the held
factorization. `hessian=` selects the Lagrangian (observed) or Gauss-Newton
(expected) form.

Two limits matter for what comes next:

1. **There is no un-inverted accessor.** `covariance()` returns the inverse
   reduced Hessian. A caller who wants the reduced Hessian itself, the
   information matrix, has to invert the covariance back. This is a
   pyomo-interface gap, not a pounce one: the core and QP interfaces already
   expose the reduced Hessian directly (see Related reduced-Hessian work
   below).
2. **The reduce-onto block is fixed at declaration time.** `covariance()`
   reports over the whole `declare_fitted` set. Asking about a different
   block means re-declaring and re-solving.

## Preconditions

Four conditions underwrite the whole surface, and the roadmap handles two of
them explicitly.

Strict complementarity failing is the weakly active case, which item 0 detects
and reports rather than assuming away. An active set that changes under
perturbation is item 2's conditional-versus-marginal distinction, which is
stated with the matrix.

The other two are assumed. LICQ is what makes the active-set KKT nonsingular,
and it is the precondition for the roadmap's own validation, since matching a
free block against the same model solved with a variable fixed is a
bounds-to-equalities substitution. Second-order sufficiency is what makes the
reduced Hessian positive definite on the free space; when it fails the result
is the indefinite free block item 1 returns with a warning, which is a
diagnostic rather than an error but is not a covariance.

## Related reduced-Hessian work in other pounce interfaces

pounce already surfaces the reduced Hessian outside pyomo-pounce. The core
`Problem.solve_with_sens` returns it in natural (unscaled) units, so
`-inv(reduced_hessian)` is directly the covariance, with an eigendecomposition
(pounce#128, mirroring sIPOPT's `rh_eigendecomp`).
`QpSensitivity.reduced_hessian` mirrors it on the convex-QP side, and
`pounce.curve_fit` is the scipy-style covariance frontend for callable models,
of which `covariance()` is the Pyomo-model sibling.

## Roadmap

Item 0 is Rust core work and everything below consumes it. Items 1 to 3 are
additive pyomo surface and can land before it, on the shipped classifier,
which is correct for interior solutions. Item 4 cannot.

**0. The bound classifier → Rust core.** Which regime each bounded variable is in,
decided by the ratio of the barrier curvature on it to the objective's own
curvature there, `Σ_i / |H_ii|`. That ratio is `O(μ)` when the bound is
inactive, `O(1)` when weakly active, and `O(1/μ)` when strongly active, so the
three separate by orders of magnitude at any `μ` and the weakly active case
sits at a fixed place rather than moving with the solve.

Thresholding `s` and `z` separately does not work at any constant. Both are
`O(sqrt(μ))` at weak activity, so a fixed threshold reports the regime only
when `μ` lands in its band. `sqrt(μ)` as the threshold fails the other way, by
putting the case the test exists to find exactly on it.

Cutoffs, for `μ ≤ 1e-4`: `inactive` below `sqrt(μ)`, `weakly active` in
`[1e-1, 1e1]`, `strongly active` above `1/sqrt(μ)`, `ambiguous` in the two
gaps. The outer edges are the geometric midpoints between the regimes'
scalings, so the weakly active case sits at the centre with
`0.5·|log μ|` decades of margin either side. Above `μ = 1e-4` the outer
edges close inside the inner band and everything not clearly outside is
`ambiguous`, which is the honest answer at that tolerance.

The inner band is fixed and the outer edges are not, because they absorb
different things. `Σ = H_ii` exactly at weak activity in the decoupled case,
so the ratio is `1` there whatever `μ` and whatever the problem scaling, and
what the band has to tolerate is the coupling drift
`Σ_i = H_ii + \sum_{j \ne i} H_{ij} x_j / x_i`, which is `O(1)` and does not
move with `μ`. Scaling the band with `μ` would widen it as the solve tightens,
swallowing more of the other two regimes exactly where the classification
should be sharpest.

The denominator has to be checked first. `|H_ii|` is a curvature scale only
while the objective actually curves in that direction; on a poorly identified
one it is noise, and the ratio inherits that. An inactive bound a full unit
away from its variable classifies as `ambiguous` once `H_ii` reaches `1e-6`
and as `strongly active` by `1e-13`, which projects the variable out and
reports zero variance for a parameter the data barely constrains. Tightening
`μ` relocates the misfile rather than removing it. So `|H_ii|` is compared to
the free block's own curvature scale before anything else, and below it the
variable returns `unidentified` with no bound regime, since the bound is not
what is wrong. Its disposition is $F$: whatever else is uncertain, a direction
the objective barely curves in must not come back with zero variance. The sign
of `H_ii` is reported alongside, since the absolute value would otherwise hide
the indefinite case item 1 returns.

`ambiguous` means the regime is undetermined at this `μ`, and the fix is to
re-solve at a tighter tolerance: the drift that puts an inactive or strongly
active variable inside the band is `O(1)`, while the edges move as
`sqrt(μ)`, so tightening separates them. The classification is returned to
the caller in every case, since which region a variable falls in is not
stable near a transition.

The classifier requires `bound_relax_factor = 0`, and checks it rather than
documenting it. The default lets a converged primal sit outside its bound by
`factor · max(1, |bound|)`, so `s` is not the distance to the user's bound and
`Σ = z/s` is scaled wrong, badly when the bound is large in magnitude. A user
hits that by doing nothing, so a stated precondition the shipped default
violates is not enough.

Two conditions are checked at every call and warned on, not only in the tests.
`s·z` away from `μ` means the point is off the central path or the bound was
relaxed, so the slack being classified is not the user's slack. `Σ_i/|H_ii|`
non-negligible on a variable the reduction kept means barrier curvature
survived the projection into a block that is supposed to be free of it.

This is new code in the Rust core, alongside `classify_working_set`
(`crates/pounce-sensitivity/src/convenience.rs`, exposed through
`crates/pounce-py`), which answers the membership question on fixed
`mult_tol` and `primal_tol` and is the only thing in this area that ships
today. It needs the solver to expose the bound multipliers, `μ`, and the
Hessian diagonal: the pyomo session carries only the primal and the bounds,
and the held factor is barrier-augmented (`sigma_x` enters the augmented
system as `d_x` in `kkt/pd_full_space_solver.rs`), so `kkt_solve` inverts $W$
and neither `Σ` nor $H$ is recoverable above it.

`diagnose_bounds` (`python/notebooks/barrier_curvature_sensitivity.ipynb` §5)
is the reference implementation of the shape and the oracle to validate
against, not a primitive to call. Its `curv_ratio` is the quantity above; its
status field, set from fixed `s` and `z` thresholds, is the part done
differently here.

**1. `information()`, the un-inverted sibling of `covariance()`.** Returns the
information matrix over the block: the reduced Hessian, formed as the Schur
complement onto the block's rows off the held factor, not by inverting the
covariance.

Natural units, the core's convention; pyomo `covariance()` carries the `2σ²`
scaling on top. Same `hessian=` selector, Lagrangian (default) or Gauss-Newton
(PSD). The Gauss-Newton path has to form `JᵀJ` over all fitted variables and
slice to the free block afterwards: it slices first today, so the pinned rows
are gone before the matrix exists and item 4's `S` has nothing to build from.

An indefinite Lagrangian free block is returned as computed, with a warning
naming Gauss-Newton as the PSD alternative. It is the honest curvature and
itself a finding, that the point is not a minimum or the model is
over-parameterized, so refusing withholds a diagnostic. Substituting
Gauss-Newton silently would return something other than the `hessian=` the
caller asked for, and a consumer that needs PSD, such as an arrival cost, can
ask for it.

Pinned variables are projected out: the information matrix is restricted to the
free block, the square block over the variables that remain, which is
`covariance()`'s existing construction. The embedding differs.
`covariance()` embeds a pinned parameter as a zero row, reading as zero
variance; the same zeros in an information matrix read as zero information.
So `information()` returns the free block plus the reduction onto the pinned
set, the finite weight on each pinned variable as it leaves its bound,
computable only from the held factor. Item 4 decides membership and gives the
expression.

It carries `covariance()`'s inertia-correction guardrail (`sens.py:814-820`).
`δ_w I` is isotropic, so unlike `Σ` it lands on the free block and survives the
projection.

**2. `wrt=` block selection.** `covariance(model, wrt=block)` and
`information(model, wrt=block)` reduce onto the given block, any free
variables, off the held factor, post-solve. The factor captured at the
solution covers every free variable, so the block is a call argument, not a
fixed declaration. `declare_fitted` becomes the default block when `wrt=` is
omitted, which keeps `covariance(model)` behaving exactly as in v0.9. The
block, whether declared or passed to `wrt=`, accepts a slice (`m.x[t, :]`) or
a `(Var, time)` pair, not just a hand-listed VarData set, so an MHE arrival
state at one time point is one call rather than an enumeration.

Each call re-reduces onto its own argument, giving that block's marginal, so
one solve serves as many blocks as you ask about.

A strongly active variable outside the block is not deleted: its `Σ` stays in
the held factor, and as `μ` falls that growing diagonal drives the coupling
through it to zero, so the block converges to the value conditional on that
bound rather than the marginal over it. The active set is returned with the
matrix, and the block's numbers move with `μ` on the way there.

**3. `retain_kkt()`, a factor-retention switch decoupled from
declarations.** The solve factors the KKT to solve the NLP; the only question
is whether that factor is kept for post-solve queries. Today any declaration
keeps it (`declare_sens_param`, `declare_fitted`, or `declare_residual`).
`retain_kkt()` keeps it with no declaration at all, which is what item 2's
declaration-free flow needs: the MHE case, where the arrival state and the
parameters are each queried by `wrt=` and neither is a default. It defaults
off, so a solve with no sensitivity pays nothing.

| setup | factor kept | `covariance(model)` | `covariance(model, wrt=T)` |
|---|---|---|---|
| nothing | no | error | error |
| `declare_fitted(S)` | yes | over S | over T |
| `retain_kkt()` only | yes | error, no default | over T |
| `retain_kkt()` + `declare_fitted(S)` | yes | over S | over T |

The columns show `covariance()`; `information()` follows the same rows, since
factor retention and the default block are accessor-agnostic. `gradient()` and
`estimate()` are unaffected, since the `declare_sens_param` they require
already keeps the factor. A block queried under that declaration alone comes
out conditional on the pinned parameter, since fixing an input conditions
rather than marginalizes.

**4. Membership and dispositions.** What `covariance()` and `information()`
each return for a variable, given item 0's classification. The barrier diagonal
it is classified on sums over both bounds,
`Σ_i = z^L_i/s^L_i + z^U_i/s^U_i`.

Write $W$ for the primal Hessian block the held factor carries, $H = W - \Sigma$
for the Lagrangian one, $F$ for the variables the reduction keeps and $A$ for
the ones it projects out, and $S = H_{AA} - H_{AF} H_{FF}^{-1} H_{FA}$ for the
reduction onto the pinned set. Each accessor returns a matrix over the whole
block; the columns are the row a fitted variable $i$ gets in each:

| bound regime | `s` | `z` | `Σ` as `μ → 0` | $i$ in | `covariance()` row | `information()` row |
|---|---|---|---|---|---|---|
| inactive | `O(1)` | `→ 0` | `μ/s² → 0` | $F$ | $2\sigma^2 (H_{FF}^{-1})_{iF}$ | $H_{iF}$ |
| strongly active | `→ 0` | `O(1)` | `z²/μ → ∞` | $A$ | $0$ | $S_{iA}$ |
| weakly active | `→ 0` | `→ 0` | finite, `O(1)` | $F$ | $2\sigma^2 (H_{FF}^{-1})_{iF}$ | $H_{iF}$ |

The `Σ` column is what skipping the subtraction in $H$ would cost.

$S$ is conditional on the rest of $A$: with more than one pinned variable,
$S_{ii}$ holds the others at their bounds rather than marginalizing over them.

`ambiguous` and `unidentified` both go to $F$ with the weakly active row, and
both warn: `ambiguous` that the regime is undetermined at this `μ` and that
re-solving tighter will settle it, `unidentified` that the objective barely
curves in that direction so the bound question does not arise and the variance
is large rather than small. $F$ is the conservative side for `covariance()`,
which reports a variance rather than asserting zero, and the
anti-conservative side for `information()`, which reports full information on
a variable that may not have it, so the warning is doing real work rather than
annotating a number that is fine.

`covariance()` ships a slack-only test today (`sens.py:826-827`,
`tol = 1e-6 * (1.0 + abs(xv))`), which pins a weakly active variable and
deletes its information.

## Scope and compatibility

Items 1 to 3 are additive: `information()` is a new function, `wrt=` (with its
slice and `(Var, time)` block forms) is a new optional keyword, and
`retain_kkt()` is new surface. No signature changes, and v0.9
`covariance(model)` with no `wrt=` reduces onto the declared set, which is
exactly the v0.10 no-argument default, so the v0.9 surface is a
forward-compatible subset.

Items 0 and 4 are the exception. Item 0 is Rust core work, since the multipliers,
`μ` and the barrier diagonal all have to reach Python before anything can
classify. Item 4 then changes which variables `covariance()` projects out, so
a model with a weakly active bound gets different numbers than v0.9 returns.

Until they land, `information()` inherits the shipped slack-only
classification, so items 1 to 3 are complete for interior solutions and misfile
a weakly active bound exactly as `covariance()` does now.

## Validation

- A `μ` sweep over all three regimes, checking together: item 0's ratio goes
  as `O(μ)`, `O(1)` and `O(1/μ)`; the weakly active value holds at `1` and its
  free-block diagonal matches the objective's curvature rather than twice it;
  every other free-block number moves by `O(μ)` and no more; and `Σ_i/|H_ii|`
  stays `O(μ)` on the variables the reduction keeps, since a non-negligible
  ratio there is barrier curvature surviving the projection. Against
  `diagnose_bounds` on the same points, including one it calls `ambiguous`.
- The classification itself across the `μ` sweep, not just the numbers: record
  every label change. An `ambiguous` variable must settle into a regime as `μ`
  tightens, which is what item 0's re-solve instruction promises, and a
  variable crossing between $F$ and $A$ changes the returned matrix's rank,
  which a caller holding a block across two solves has to see.
- Item 0 under variable scaling: hold the regime and move the objective's
  curvature away from the bound's scale. The weakly active ratio stays at `1`;
  the other two drift toward the band, and at loose `μ` they enter it, which
  must come back as `ambiguous` rather than as a confident wrong regime.
- One-sided finite differences at a transition, against the classification and
  against `covariance()`'s free block. A symmetric difference there returns the
  barrier's smoothing value rather than either true one-sided derivative, so it
  passes while hiding the active-set change; only the one-sided pair shows the
  two disagree. Run with `bound_relax_factor = 0`, or the slacks being
  classified are not distances to the user's bounds.
- `information()` on a problem where the solve reports non-zero
  `kkt_perturbations`: the inertia-correction guardrail fires, and the returned
  block differs from the unregularized one by the isotropic `δ_w` and nothing
  else.
- `information(...)` against `inv(covariance(...))` to tolerance on a
  well-conditioned block with no active bound and pooled residuals; the
  conditioning advantage on a deliberately ill-identified one. The identity
  holds only there: a bound active makes `covariance()` singular (the pinned
  row is zero) and grouped residuals of unequal variance make the two objects
  different.
- A strongly active fitted variable: the free block matches the same model solved
  with that variable fixed (a bounds-to-equalities substitution, so LICQ is
  assumed), and the pinned variable reports $S_{ii}$ with the activity
  classification. $S_{ii}$ is checkable directly: it is the second derivative
  of the objective minimized over the free variables with that variable held at
  a value and its own bound dropped, and the first derivative of that same
  function is the bound multiplier.
- Two strongly active fitted variables: the off-diagonal $S_{ij}$ matches the
  same construction, and $S_{ii}$ differs from the same variable's value when
  the other is free.
- The `μ` sweep is necessary, not sufficient, on its own: the weakly-active
  case is `μ`-invariant and barrier-inflated at once, so it only certifies
  anything alongside item 0's classification. A block conditioned on a strongly
  active variable outside it does move with `μ`, which is correct.
- `s·z` against `μ` on every bounded variable. A mismatch means the point is
  off the central path or the solver relaxed the bound, and the slack being
  classified is not the true slack.
- An indefinite Lagrangian free block is returned as computed, with the
  warning, and is not silently replaced by the Gauss-Newton form.
- The marginal identity: `inv(state block of covariance(wrt={state,
  params}))` against `information(wrt=state)`, both the parameter-marginal
  state information.
- The conditional identity: the state block of `information(wrt={state,
  params})` against `information(wrt=state)` computed with the parameters
  fixed.
- Lagrangian versus Gauss-Newton agree on a linear model and in the
  small-residual limit; the Lagrangian can go indefinite where Gauss-Newton
  stays PSD.
