# Covariance and information: v0.10 roadmap

**Status: roadmap proposal for pyomo-pounce, targeting v0.10.** This note
scopes the post-solve second-order surface in pyomo-pounce: `covariance()`,
shipping in v0.9, and the additions v0.10 should make around it. Everything
here is additive to the v0.9 `covariance()` surface; nothing changes an
existing signature. Companion to the active-set sensitivity roadmap
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
per-block object is what the roadmap below adds; the concrete per-window loop
is in the MHE section.

## Roadmap

**1. `information()`, the un-inverted sibling of `covariance()`.** Returns
the information matrix over the block: the reduced Hessian, formed directly
from the held factor (the Schur complement onto the block's rows) rather than
by inverting the covariance, which skips the round trip and stays well-scaled
in the poorly-identified directions.

The definition is the reduced Hessian, not `inv(covariance(...))`. The two
agree only in the homoscedastic case: `covariance()` returns a heteroscedastic
sandwich whenever labeled residual groups carry unequal variances
(`sens.py:995` Lagrangian, `sens.py:969` Gauss-Newton), and inverting a
sandwich is not a scaled reduced Hessian by any scalar, however well
conditioned the block. Units follow the same care: pyomo `covariance()` carries
the `2 * sigma_sq` scaling while the core's reduced Hessian is natural-units,
and a consumer whose objective already carries its own inverse-covariance
weights wants the unscaled one, so the note must say which `information()`
returns rather than promising consistency with both. Same
`hessian=` selector: Lagrangian (default, the exact reduced Hessian, what
the information-form arrival cost wants) and Gauss-Newton (PSD). Same
scaling conventions as `covariance()`; the bound handling is shared only in
part, see Active bounds below.

**2. `wrt=` block selection on both.** `covariance(model, wrt=block)` and
`information(model, wrt=block)` reduce onto the given block, any free
variables, off the held factor, post-solve. The factor captured at the
solution covers every free variable, so the block is a call argument, not a
fixed declaration. `declare_fitted` becomes the default block when `wrt=` is
omitted, which keeps `covariance(model)` behaving exactly as in v0.9. Each
call reduces onto its own argument, so one solve serves as many blocks as
you ask about. The block, whether declared or passed to `wrt=`, accepts a
slice (`m.x[t, :]`) or a `(Var, time)` pair, not just a hand-listed VarData
set, so an MHE arrival state at one time point is one call rather than an
enumeration.

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

## Marginal versus conditional: the one semantic to get right

Each call reduces onto its argument, and that yields the block's
**marginal**: everything not in the block, other states, parameters, is
integrated out through the KKT reduction. This is what the arrival cost
wants.

The trap is the asymmetry between the two objects. Slicing a covariance
gives a marginal; slicing an information matrix gives a **conditional**. So
`information(wrt=T)` must re-reduce onto `T`, not slice a joint reduction
over some larger set. If it sliced, the arrival-state block of a joint
`{state, params}` information would be the information conditional on the
parameters held fixed, not the marginal that carries their uncertainty.
Making `wrt=` mean "reduce onto this" gives the marginal directly, and the
answer does not depend on what else was declared.

A useful consequence: `information(wrt=arrival_block)` marginalizes the
parameters for free, because they are simply not in the block. You reach for
the conditional only by deliberately putting the parameters in the block and
slicing.

## Active bounds

Three things have to be right: which directions count as bound-active, what the
projection does with them, and what gets reported. The barrier curvature study
(`python/notebooks/barrier_curvature_sensitivity.ipynb`) derives the regimes;
its operative sentence is that an active bound deletes a direction and adds no
intrinsic curvature.

### Three regimes, not two

The barrier adds a diagonal `Sigma_i = z_i / s_i = mu / s_i^2` to the primal
Hessian. Its size depends on the activity regime, not merely on whether a bound
exists:

| regime | slack `s` | multiplier `z` | `Sigma` as `mu -> 0` |
|---|---|---|---|
| inactive | `O(1)` | `-> 0` | `mu / s^2 -> 0` |
| strongly active | `-> 0` | `O(1)` | `z^2 / mu -> infinity` |
| weakly active | `-> 0` | `-> 0` | finite, `-> q` |

The first two the projection handles cleanly. The large term sits exactly on
the direction a projection deletes, and the vanishing one is `O(mu)`, so under
strict complementarity the free block is the objective's own curvature and no
separate scrubbing step is needed.

The third breaks that. At a weakly active bound `s = sqrt(mu/q)` and
`z = sqrt(mu*q)`, so `Sigma = q` exactly, at every `mu`: finite, independent of
`mu`, and the same order as the objective's own curvature, on a direction whose
multiplier is near zero. Classified free, the block carries `2q` instead of
`q`, so the information is doubled and the error does not shrink as `mu` falls.
Classified pinned, a direction carrying genuine finite information is deleted.
Neither is right, so activity needs a third outcome rather than a free/pinned
split, and the strict-complementarity qualifier above is load bearing.

### Activity detection

Classification is on slack **and** multiplier, never distance to the bound
alone, with tolerances tied to `mu` (compare `s` to `sqrt(mu)`, and `s*z` to
`mu`) rather than fixed constants. Both small is weakly active, which is
flagged: neither silently included nor silently deleted.

The rule that ships is slack-only (`sens.py:826-827`, `tol = 1e-6 * (1.0 +
abs(xv))`), which misclassifies exactly that case. The study's Ipopt run at a
weakly active point returns `x = 5.64e-7` with `z_L = 5.64e-7`, inside that
tolerance, so a direction with finite information is projected out. Upgrading
the test changes `covariance()` as well as `information()`, so it is its own
change rather than something `information()` inherits quietly. A solver that
relaxed the bound reports a slack that is not the true slack, which the test
has to account for.

### The projection is shared, the embedding is not

A fitted variable sitting at one of its bounds at the optimum is projected
out: both accessors work in the free (off-bound) directions. `covariance()`
restricts the information matrix to the free block and inverts that, rather
than inverting first and restricting, which is the construction it already
ships. `information()` restricts the same way.

What cannot carry over is the embedding. `covariance()` embeds a pinned
parameter as a zero row and column, which reads as zero variance and is
correct for something the active set holds fixed. The same zeros in an
information matrix read as zero information, that is, infinite variance: the
same array making the opposite claim. So `information()` returns the reduced
Hessian over the free block and reports the pinned directions as active-set
metadata, not as matrix entries. A pinned direction has no finite information
to report; conditional on the bound staying active it is unbounded, and
unbounded is not something a caller can multiply.

The failure mode under strict complementarity is skipping the projection, not
picking the wrong Hessian: the bound-active entry of the full
`(W + Sigma)^{-1}` is a variance collapsing to zero, a constraint artifact read
as precision.

### What the projection does not remove

`Sigma` is rank-structured onto the deleted directions, so the projection
handles it. Inertia correction is not. `delta_w * I` is isotropic: it lands on
the free block and survives the projection. pounce bakes it into the held
factor (`kkt_perturbations()` in `solver.rs`) and `covariance()` already guards
on it (`sens.py:814-820`); `information()` inherits the same guardrail. This
matters most where item 1 claims the headline advantage, since `delta` is
injected precisely when the Hessian is indefinite or near-singular, which is
the poorly-identified regime.

### Bound-active variables outside the block

The above concerns bound-active members of the block. One outside it is not
deleted: its `Sigma` lives in the held factor and enters the Schur complement
onto the block, so the result is conditional on that bound rather than marginal
over it. On the study's coupled example `information(wrt=x2)` converges to 2.0,
conditional on `x1` pinned, where the bound-free marginal is 1.5. Conditional
information always dominates the marginal, so the caller gets a silently
over-confident matrix, and unlike a clean free block this number does move with
`mu` (1.9938 at `mu=1e-2`, 2.0000 at `mu=1e-8`). Both facts qualify the
Marginal versus conditional section above: the reduction marginalizes over free
directions and conditions on bound-active ones, so the answer does depend on
the active set even though it does not depend on the declarations.

## MHE, the motivating consumer

Per window the estimator solves its NLP once with `retain_kkt()` set, then
reduces onto the arrival state for the next window's `Pi^{-1}` and onto the
parameter block for their covariance. The arrival block is the components of
the state that become the next window's start, one time slice, not the whole
window.

The reduction is the mechanism; which problem it is applied to is the
estimator's business, and it is not the whole window. The arrival-cost update
is one step,

    Gamma_{k+1}(z) = min over the dropped state of
                     [ Gamma_k(x_dropped) + the one stage cost leaving the window ]

subject to the arrival state equalling `z`. Reducing the entire solved window
onto the arrival state instead would fold the overlap's stage costs into
`Pi^{-1}` while those same residuals stay live in the next window's objective,
counting every measurement in the overlap twice. On a scalar linear-Gaussian
three-point window that is 6.28x over-confident, and the factor grows with the
horizon. So the consumer builds the one-step subproblem and reduces that.

Whether the arrival cost carries parameter uncertainty (marginal) or treats
parameters as known (conditional) is a modeling choice, set by whether the
parameters are in the block.

A bound-active arrival state is where the reporting convention matters, because
the consumer cannot reconstruct what the accessor discards. Three candidate
weights exist for a pinned direction: zero, the barrier's `z^2/mu`, and the
retained row with `Sigma` removed. On the study's coupled example those are 0,
1.6e8, and 1.5, and only the third is both finite and meaningful. Only pounce
can compute it, so reporting zero or dropping the row is lossy in a way the
caller cannot undo. `information()` should expose the retained-row value
alongside the activity classification and let the estimator decide what its
arrival cost does with it.

## Scope boundary: mechanism in pounce, policy in the caller

`information()`, `covariance()`, `wrt=`, and `retain_kkt()` are mechanisms:
stateless queries against the held factorization. What a consumer does with the
numbers is policy and belongs downstream, the same split the active-set
sensitivity roadmap draws. Which block to ask about, how to build the
subproblem being reduced, what an arrival cost weights, and what to do when the
active set churns between windows are all the estimator's, not pounce's. The
MHE material above is motivation for the mechanism and a check that the
mechanism suffices, not a specification of the estimator.

The one place the boundary is not clean is reporting: an accessor that discards
a number the caller cannot recompute has made a policy decision by omission.
That is why the pinned-direction convention is settled here rather than
downstream.

## Scope and compatibility

pyomo-pounce only. All three items are additive to v0.9: `information()` is
a new function, `wrt=` (with its slice and `(Var, time)` block forms) is a
new optional keyword, and `retain_kkt()` is new surface. Nothing
changes an existing signature. v0.9 `covariance(model)` with no `wrt=` reduces onto the declared
set, which is exactly the v0.10 no-argument default, so the v0.9 surface is
a forward-compatible subset. Nothing here needs to be rushed into v0.9.

## Validation

- `information(...)` against `inv(covariance(...))` to tolerance on a
  well-conditioned block with no active bound and pooled (homoscedastic)
  residuals; the conditioning advantage on a deliberately ill-identified one.
  The identity is scoped twice over: with a bound active `covariance()` is
  singular by construction (the pinned row is zero), so the round trip is
  undefined rather than ill-conditioned, and with grouped residuals of unequal
  variance the two objects genuinely differ, so a separate statement of what
  they satisfy there is owed.
- A bound-active fitted variable: the free block matches the same model solved
  with that variable fixed (which presumes LICQ, since that is a
  bounds-to-equalities substitution), and the pinned direction is reported as
  active with its retained-row value rather than as a zero or an inflated
  entry.
- Refining the solver's `mu` moves the free-block numbers by `O(mu)` and no
  more. This is necessary, not a certificate: the weakly-active case is
  `mu`-invariant and barrier-inflated at once, so pair it with the
  slack-and-multiplier classification. A block whose reduction conditions on a
  bound-active variable outside it does move with `mu`, and that is correct
  rather than a failure.
- A weakly-active fitted variable (slack and multiplier both near zero): it is
  flagged as degenerate rather than silently classified free (which would carry
  `2q`) or pinned (which would drop finite information), and the flag survives
  a sweep in `mu`.
- An indefinite Lagrangian free block: `information()` does whatever the note
  settles on, refuse or fall back to Gauss-Newton, rather than returning a
  matrix that would make a downstream quadratic unbounded below.
- The marginal identity: `inv(state block of covariance(wrt={state,
  params}))` against `information(wrt=state)`, both the parameter-marginal
  state information.
- The conditional identity: the state block of `information(wrt={state,
  params})` against `information(wrt=state)` computed with the parameters
  fixed.
- Lagrangian versus Gauss-Newton agree on a linear model and in the
  small-residual limit; the Lagrangian can go indefinite where Gauss-Newton
  stays PSD.
- MHE recursion sanity: the posterior information equals prior plus data on a
  linear-Gaussian window.
