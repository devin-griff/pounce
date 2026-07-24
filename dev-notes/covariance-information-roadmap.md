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
Gauss-Newton form. The information form is the natural one for an
information-form arrival cost.

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
work below). This is an interface contribution on the pyomo side, layered on
them.

## Where we are

v0.9 ships `covariance()` (in `pyomo-pounce/pyomo_pounce/sens.py`). You
`declare_fitted` a set of free variables, solve, and `covariance(model)`
returns their asymptotic covariance: the scaled parameter block of the
inverse KKT matrix, `2 * sigma_sq * (K^-1)_pp`. Under the hood that block
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
arrival cost is `Gamma(x0) = 0.5 (x0 - xhat)^T Pi^{-1} (x0 - xhat)`, where the
weighting `Pi^{-1}` is the reduced Hessian marginalized onto the arrival
state, the Lagrangian information, not the covariance. That un-inverted,
per-block object is what the roadmap below adds.

## Roadmap

**1. `information()`, the un-inverted sibling of `covariance()`.** Returns
the information matrix over the block: the reduced Hessian, formed directly
from the held factor (the Schur complement onto the block's rows) rather than
by inverting the covariance, which skips the round trip and stays well-scaled
in the poorly-identified directions.

It returns natural units, the core's convention, since a consumer whose
objective already carries its own inverse-covariance weights needs the
unscaled object; pyomo `covariance()` carries the `2 * sigma_sq` scaling on
top. `inv(covariance(...))` recovers it only for pooled residuals: with
labeled residual groups of unequal variance `covariance()` is a
heteroscedastic sandwich (`sens.py:995` Lagrangian, `sens.py:969`
Gauss-Newton), whose inverse is no scalar multiple of a reduced Hessian.
Same `hessian=` selector: Lagrangian (default, the exact reduced Hessian,
what the information-form arrival cost wants) and Gauss-Newton (PSD).

Bound-active directions are projected out: the information matrix is
restricted to the free block, which is `covariance()`'s existing construction,
never inverted first and then restricted. The embedding differs.
`covariance()` embeds a pinned parameter as a zero row, reading as zero
variance; the same zeros in an information matrix read as zero information,
the opposite claim. So `information()` returns the free block plus, for each
pinned direction, the retained row with `Sigma` removed: the finite weight,
describing how the objective curves as that variable leaves its bound, and
computable only from the held factor. Item 4 gives the per-regime dispositions.

It carries `covariance()`'s inertia-correction guardrail (`sens.py:814-820`).
`Sigma` is rank-structured onto the deleted directions, so the projection
removes it; `delta_w * I` is isotropic, lands on the free block, and survives.
pounce bakes it into the held factor (`kkt_perturbations()` in `solver.rs`),
and it is injected precisely where the Hessian is indefinite or near-singular,
which is the poorly-identified regime.

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

A bound-active variable outside the block is not deleted: its `Sigma` stays in
the held factor, and as `mu` falls that growing diagonal drives the coupling
through it to zero, so the block converges to the value conditional on that
bound rather than the marginal over it. The active set is returned with the
matrix, and the block's numbers move with `mu` on the way there.

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

**4. Joint activity classification.** `covariance()` and `information()`
classify a bound as active on slack **and** multiplier, with tolerances tied to
`mu` (compare `s` to `sqrt(mu)`, and `s*z` to `mu`), giving three outcomes.
The barrier's diagonal `Sigma_i = z_i / s_i = mu / s_i^2` separates them:

| regime | slack `s` | multiplier `z` | `Sigma` as `mu -> 0` | reported as |
|---|---|---|---|---|
| inactive | `O(1)` | `-> 0` | `mu / s^2 -> 0` | in the free block, unadjusted |
| strongly active | `-> 0` | `O(1)` | `z^2 / mu -> infinity` | out of the free block, its retained row less `Sigma` |
| weakly active | `-> 0` | `-> 0` | finite, `O(1)` | in the free block, its diagonal less `Sigma` |

`Sigma` comes off in the last two rows for the same reason and by the same
subtraction: it is the barrier's, not the objective's. It is left in the first
because `O(mu)` is below the tolerance anything else is computed to. A weakly
active variable sits within `O(sqrt(mu))` of its bound while carrying finite
information of the same order as the objective's own curvature, and its
multiplier is near zero, so the bound is not holding it; leaving `Sigma` in
there would roughly double the curvature reported for that direction, an error
that does not shrink with `mu`. The classification is returned with the matrix
in every regime, since which side of a tolerance a variable falls on is not
stable.

`covariance()` ships a slack-only test today (`sens.py:826-827`,
`tol = 1e-6 * (1.0 + abs(xv))`), which pins such a variable and deletes its
information, so this is the one item that changes `covariance()`'s numbers
rather than only adding surface. A solver that relaxed the bound reports a
slack that is not the true slack.

## Scope and compatibility

pyomo-pounce only. Items 1 to 3 are additive: `information()` is a new
function, `wrt=` (with its slice and `(Var, time)` block forms) is a new
optional keyword, and `retain_kkt()` is new surface. No signature changes, and
v0.9 `covariance(model)` with no `wrt=` reduces onto the declared set, which is
exactly the v0.10 no-argument default, so the v0.9 surface is a
forward-compatible subset.

Item 4 is the exception: the joint activity test changes which directions
`covariance()` projects out, so a model with a weakly active bound gets
different numbers than v0.9 returns.

## Validation

- `information(...)` against `inv(covariance(...))` to tolerance on a
  well-conditioned block with no active bound and pooled residuals; the
  conditioning advantage on a deliberately ill-identified one. The identity
  holds only there: a bound active makes `covariance()` singular (the pinned
  row is zero) and grouped residuals of unequal variance make the two objects
  different.
- A bound-active fitted variable: the free block matches the same model solved
  with that variable fixed (a bounds-to-equalities substitution, so LICQ is
  assumed), and the pinned direction reports its retained-row value with the
  activity classification.
- Refining the solver's `mu` moves the free-block numbers by `O(mu)` and no
  more. Necessary, not sufficient: the weakly-active case is `mu`-invariant and
  barrier-inflated at once, so it pairs with the slack-and-multiplier
  classification. A block conditioned on a bound-active variable outside it
  does move with `mu`, which is correct.
- A weakly-active fitted variable (slack and multiplier both near zero) stays
  in the free block with its diagonal matching the objective's curvature in
  that direction, not twice it, and is reported as weakly active. Both hold
  across a sweep in `mu`.
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
