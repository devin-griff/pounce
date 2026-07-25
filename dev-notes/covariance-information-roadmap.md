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

Two things the pyomo-pounce interface lacks, which no pounce interface offers
together:

- an `information()` accessor consistent with `covariance()` and the core's
  natural-units reduced Hessian, so a Pyomo model gets the un-inverted
  object without the invert-then-reinvert round trip; and
- post-solve `wrt=` block selection off one retained factor, reducing onto
  arbitrary free-variable blocks from a single solve, the
  one-solve-two-blocks flow the MHE arrival cost needs, with `retain_kkt()`
  as the declaration-free enabler.

The reduced Hessian and the covariance recipe are established, and pounce ships
them in its core, QP, and `curve_fit` interfaces (see Related reduced-Hessian
work below).

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

## Related reduced-Hessian work in other pounce interfaces

pounce already surfaces the reduced Hessian outside pyomo-pounce. The core
`Problem.solve_with_sens` returns it in natural (unscaled) units, so
`-inv(reduced_hessian)` is directly the covariance, with an eigendecomposition
(pounce#128, mirroring sIPOPT's `rh_eigendecomp`).
`QpSensitivity.reduced_hessian` mirrors it on the convex-QP side, and
`pounce.curve_fit` is the scipy-style covariance frontend for callable models,
of which `covariance()` is the Pyomo-model sibling.

## The MHE arrival cost

The motivating consumer is moving horizon estimation. Its information-form
arrival cost is `Γ(x0) = ½ (x0 − x̂)ᵀ Π⁻¹ (x0 − x̂)`, where the
weighting `Π⁻¹` is the reduced Hessian of the one-step subproblem (the
previous arrival cost plus the stage leaving the window) marginalized onto the
arrival state, the Lagrangian information, not the covariance. That
un-inverted, per-block object is what the roadmap below adds.

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

The inner band is wider than the ideal case needs, since `Σ = H_ii` exactly
at weak activity and the ratio is `1` there whatever the problem scaling.
It is set for the other two, which carry scale: `μ/(s² H_ii)` and
`z²/(μ H_ii)` drift toward the centre when the objective's curvature is far
from the bound's own scale, and the drift is what the gaps absorb.

The classification is returned to the caller in every case, since which
region a variable falls in is not stable near a transition.

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

Rust rather than the pyomo frontend because the consumers are plural: item 4
below, the active-set sensitivity roadmap's item 0 report and item 3 trigger,
and the jax and torch paths, whose active-set-margin monitor asks a related
question.

**1. `information()`, the un-inverted sibling of `covariance()`.** Returns the
information matrix over the block: the reduced Hessian, formed as the Schur
complement onto the block's rows off the held factor, not by inverting the
covariance.

Natural units, the core's convention; pyomo `covariance()` carries the `2σ²`
scaling on top. Same `hessian=` selector, Lagrangian (default) or Gauss-Newton
(PSD).

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

`ambiguous` goes to $F$ with the weakly active row. That is the conservative
side for `covariance()`, which reports a variance rather than asserting zero,
and the anti-conservative side for `information()`, which reports full
information on a variable that may not have it. Item 0's classification is
what carries the doubt either way.

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

- Item 0's ratio across a `μ` sweep on all three regimes: `O(μ)`, `O(1)` and
  `O(1/μ)`, with the weakly active value fixed as `μ` falls. Against
  `diagnose_bounds` on the same points, including one it calls `ambiguous`.
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
- Refining the solver's `μ` moves the free-block numbers by `O(μ)` and no
  more. Necessary, not sufficient: the weakly-active case is `μ`-invariant and
  barrier-inflated at once, so it pairs with item 0's
  classification. A block conditioned on a strongly active variable outside it
  does move with `μ`, which is correct.
- A weakly-active fitted variable (slack and multiplier both near zero) stays
  in the free block with its diagonal matching the objective's curvature in
  that variable, not twice it, and is reported as weakly active. Both hold
  across a sweep in `μ`.
- `Σ_i / |H_ii|` on every variable the reduction keeps. It is `O(μ)` on a
  genuinely inactive one, so a non-negligible ratio there is barrier curvature
  surviving into the free block, which the `μ` sweep alone does not catch.
- `s·z` against `μ` on every bounded variable. A mismatch means the point is
  off the central path or the solver relaxed the bound, and the slack being
  classified is not the true slack.
- An indefinite Lagrangian free block returns the settled outcome, refusal or a
  Gauss-Newton fallback, not a matrix that makes a downstream quadratic
  unbounded below.
- The marginal identity: `inv(state block of covariance(wrt={state,
  params}))` against `information(wrt=state)`, both the parameter-marginal
  state information.
- The conditional identity: the state block of `information(wrt={state,
  params})` against `information(wrt=state)` computed with the parameters
  fixed.
- Lagrangian versus Gauss-Newton agree on a linear model and in the
  small-residual limit; the Lagrangian can go indefinite where Gauss-Newton
  stays PSD.
- MHE sanity: reducing a linear-Gaussian one-step subproblem onto its arrival
  state matches the closed-form update.
