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
Gauss-Newton form. So across the field the un-inverted information matrix is
not the rare quantity but the native engine output, with the covariance
derived from it. The information form is the natural one for an
information-form arrival cost.

## Benefit hypothesis

The contribution is not the reduced Hessian or the covariance recipe. Both
are established, and pounce already ships them in its core, QP, and
`curve_fit` interfaces (see Related reduced-Hessian work below). It is two
things the
pyomo-pounce interface lacks and that no pounce interface offers together:

- an `information()` accessor consistent with `covariance()` and the core's
  natural-units reduced Hessian, so a Pyomo model gets the un-inverted
  object without the invert-then-reinvert round trip; and
- post-solve `wrt=` block selection off one retained factor, reducing onto
  arbitrary free-variable blocks from a single solve, the
  one-solve-two-blocks flow the MHE arrival cost needs, with `retain_kkt()`
  as the declaration-free enabler.

So this is an interface and ergonomics contribution on the pyomo side,
layered on the core's existing reduced Hessian, not new numerics.

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

So the object and its natural-units convention already exist in pounce. This
roadmap brings them to the pyomo-pounce interface, which today exposes only
`covariance()`, the inverse, over a fixed declared set, with no reduced-Hessian
accessor and no per-call block.

## The MHE arrival cost

The motivating consumer is moving horizon estimation. Its information-form
arrival cost is `Gamma(x0) = 0.5 (x0 - xhat)^T Pi^{-1} (x0 - xhat)`, where the
weighting `Pi^{-1}` is the reduced Hessian marginalized onto the arrival
state, the Lagrangian information, not the covariance. That un-inverted,
per-block object is what the roadmap below adds.

The estimator solves with `retain_kkt()` set and then reduces: onto the arrival
state, one time slice, for `Pi^{-1}`, and onto the parameter block for their
covariance. Those two calls off one factor are the whole demand this note has
to meet. Which problem it reduces is the estimator's, and it is a one-step
subproblem rather than the whole window, since overlapping windows would
otherwise count the overlap's measurements twice.

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
computable only from the held factor. The activity classification is returned
with them.

It carries `covariance()`'s inertia-correction guardrail (`sens.py:814-820`).
`Sigma` is rank-structured onto the deleted directions, so the projection
removes it; `delta_w * I` is isotropic, lands on the free block, and survives.
pounce bakes it into the held factor (`kkt_perturbations()` in `solver.rs`),
and it is injected precisely where the Hessian is indefinite or near-singular,
which is the poorly-identified regime.

**2. `wrt=` block selection on both.** `covariance(model, wrt=block)` and
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
declarations.** The solve factors the KKT to solve the NLP; the only
question is whether that factor is kept for post-solve queries. Today it is
kept whenever a declaration is present (`declare_sens_param`,
`declare_fitted`, or `declare_residual`). Item 2 lets the block move to a
call argument, so a caller may want no declaration at all. `retain_kkt()`
keeps the factor without committing to a block, param, or residual. It
defaults off, so a solve with no sensitivity pays nothing.

`retain_kkt()` is not specific to this surface. Keeping the factor is the
substrate every sensitivity feature rests on, and `gradient()` and
`estimate()` already get it as a side effect of the `declare_sens_param`
they require, so they never call `retain_kkt()` and it has no user-facing
effect on them. The only flow that needs a declaration-free retain is
`covariance` / `information` driven purely by `wrt=`.

`wrt=` itself is not gated by `retain_kkt()`. It needs the factor kept, and
`declare_fitted` already keeps it, so `wrt=` works with `declare_fitted`
alone. `retain_kkt()` earns its place only when you want the factor kept
with no default block at all: the declaration-free MHE case, where the
arrival state and the parameters are each queried by `wrt=` and neither is a
default.

| setup | factor kept | `covariance(model)` | `covariance(model, wrt=T)` |
|---|---|---|---|
| nothing | no | error | error |
| `declare_fitted(S)` | yes | over S | over T |
| `retain_kkt()` only | yes | error, no default | over T |
| `retain_kkt()` + `declare_fitted(S)` | yes | over S | over T |

The columns show `covariance()` for concreteness; `information()` follows the
same rows, since factor retention and the default block are accessor-agnostic.

Any declaration keeps the factor, not only `declare_fitted`, so `retain_kkt()`
is needed only when nothing at all is declared. In particular
`declare_sens_param` alone (no `declare_fitted`, no `retain_kkt()`) does
support `covariance(model, wrt=T)` and `information(model, wrt=T)` off the same
solve. It just carries no default block, so a bare `covariance(model)` errors,
exactly the `retain_kkt()`-only row. The block `T` then comes out conditional
on the pinned parameter, since fixing an input conditions rather than
marginalizes.

**4. Joint activity classification.** Both accessors classify a bound as active
on slack **and** multiplier, with tolerances tied to `mu` (compare `s` to
`sqrt(mu)`, and `s*z` to `mu`), giving three outcomes: free, pinned, and weakly
active, which is flagged rather than forced into either. The barrier's
diagonal `Sigma_i = z_i / s_i = mu / s_i^2` is what separates them:

| regime | slack `s` | multiplier `z` | `Sigma` as `mu -> 0` |
|---|---|---|---|
| inactive | `O(1)` | `-> 0` | `mu / s^2 -> 0` |
| strongly active | `-> 0` | `O(1)` | `z^2 / mu -> infinity` |
| weakly active | `-> 0` | `-> 0` | finite, `O(1)` |

A weakly active variable sits within `O(sqrt(mu))` of its bound while carrying
finite information of the same order as the objective's own curvature. Treated
as free its block carries roughly twice that curvature, an error that does not
shrink with `mu`; treated as pinned, finite information is deleted.
`covariance()` ships a slack-only test today (`sens.py:826-827`,
`tol = 1e-6 * (1.0 + abs(xv))`), which pins it, so this is the one item that
changes `covariance()`'s numbers rather than only adding surface. A solver that
relaxed the bound reports a slack that is not the true slack.

## Scope and compatibility

pyomo-pounce only. Items 1 to 3 are additive: `information()` is a new
function, `wrt=` (with its slice and `(Var, time)` block forms) is a new
optional keyword, and `retain_kkt()` is new surface. No signature changes, and
v0.9 `covariance(model)` with no `wrt=` reduces onto the declared set, which is
exactly the v0.10 no-argument default, so the v0.9 surface is a
forward-compatible subset.

Item 4 is the exception. The joint activity test changes which directions
`covariance()` projects out, so a model with a weakly active bound gets
different numbers than v0.9 returns. That is a fix rather than a break, but it
is a behavior change and wants its own release note. Nothing here needs to be
rushed into v0.9.

## Validation

- `information(...)` against `inv(covariance(...))` to tolerance on a
  well-conditioned block with no active bound and pooled residuals; the
  conditioning advantage on a deliberately ill-identified one. The identity
  holds only there: a bound active makes `covariance()` singular (the pinned
  row is zero) and grouped residuals of unequal variance make the two objects
  genuinely different.
- A bound-active fitted variable: the free block matches the same model solved
  with that variable fixed (a bounds-to-equalities substitution, so LICQ is
  assumed), and the pinned direction reports its retained-row value with the
  activity classification.
- Refining the solver's `mu` moves the free-block numbers by `O(mu)` and no
  more. Necessary, not sufficient: the weakly-active case is `mu`-invariant and
  barrier-inflated at once, so it pairs with the slack-and-multiplier
  classification. A block conditioned on a bound-active variable outside it
  does move with `mu`, which is correct.
- A weakly-active fitted variable (slack and multiplier both near zero) is
  flagged, not classified free (carrying `2q`) or pinned (dropping finite
  information), and the flag survives a sweep in `mu`.
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
