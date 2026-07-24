# Review: the covariance and information roadmap

**Status: working review, not for PR #262.** A read of
`covariance-information-roadmap.md` against John Kitchin's barrier-curvature
study (`python/notebooks/barrier_curvature_sensitivity.ipynb`) and against the
shipped `covariance()` implementation. Sixteen findings survived an adversarial
check; fourteen more were raised and rejected as churn. Every numerical claim
below was recomputed rather than taken from the roadmap.

## Bottom line

The roadmap's central mechanism is wrong, and it is wrong independently of the
bounds discussion that prompted this review. The MHE recursion double-counts by
a factor that grows with the horizon. The bounds material added after reading
the study fixed one real defect (the embedding asymmetry) and introduced
several new ones, because it reasons from a two-regime picture where the study
establishes three.

Nothing here is on PR #262. The MHE defect predates today's edits and is in the
version the PR currently carries.

The four clusters, worst first.

## 1. The MHE section

**1.1 The recursion double-counts.** The roadmap says
`information(model, wrt=arrival_block)` on the solved window is `Pi^{-1}` for
the next window, "giving prior plus window data, the recursion." It is not.
MHE windows overlap in N-1 time points, so reducing the *whole* window folds
the overlap's stage costs into `Pi^{-1}`, and those same measurement residuals
and process-noise terms appear again as live terms in the next window's
objective. Every measurement in the overlap is counted twice, once as a weight
and once as a residual.

Checked on a scalar linear-Gaussian window (a=0.9, process information 1.0,
measurement information 4.0, prior information 2.0, three points, reduced onto
x1): the full-window reduction gives arrival information 5.529 where the
one-step recursion gives 0.881. That is 6.28x over-confident on a three-point
window, and the error grows with horizon length. It is benign only at N=1 or
for non-overlapping windows, neither of which the roadmap mentions.

The exact recursion is one step:
`Gamma_{k+1}(z) = min over the dropped state of [Gamma_k + the one dropped
stage cost]`, subject to the arrival state equalling z. So `Pi^{-1}` is the
reduction of the dropped stage's subproblem, not of the full window.

The roadmap's validation item is phrased in the same wrong terms ("the
posterior information equals prior plus data on a linear-Gaussian window"), so
it would ratify the double count rather than catch it.

**1.2 "The bound does its own work" is false.** The roadmap argues that a
bound-active arrival state needs no information entry because the bound is
still a bound in the next window. A bound is one-sided. It blocks motion below
`l` and does nothing above it. An absent information entry licenses motion in
both directions at zero arrival-cost penalty, so on the feasible half-line the
next window has neither a penalty nor a constraint and can revise the arrival
state arbitrarily far into the interior on weak new data. That is precisely the
failure the roadmap names one sentence earlier and then declares unnecessary.

"Double-counts" is also the wrong description of the alternative. A hard
constraint and a stiff quadratic coincide only exactly at the bound; off it the
quadratic penalizes feasible motion the constraint permits. The conclusion (do
not carry barrier curvature) is right, but for the study's reason, that
`Sigma ~ z^2/mu` is a homotopy artifact, not for the roadmap's reason.

**1.3 The pinned entry is finite, not unbounded.** The roadmap says a pinned
direction is unbounded and therefore unreportable. That holds only under the
conditioning it silently chose, conditioning on the bound staying exactly
active. The MHE consumer wants a different and well-posed question: how much
does the past cost curve as the arrival state moves off the bound into the
interior.

On the study's own Section 4 example (W=[[2,1],[1,2]], b=[-1,0.5], x1 at its
lower bound with z1=1.25), define `Gamma(z)` as the past problem with the
arrival state pinned to z and its own bound dropped. Then `Gamma'(0) = 1.2500`,
exactly the bound multiplier, and `Gamma''(0) = 1.5000`, exactly
`W11 - W12 W22^{-1} W21`. Both finite. The three candidate weights for that
arrival state are therefore 0 (the roadmap's design), 1.6e8 (barrier), and 1.5
(correct). The roadmap rejects the second and ships the first.

The linear term has no home at all: `0.5 (x0-xhat)^T Pi^{-1} (x0-xhat)` has
zero gradient at `xhat` by construction, so the quadratic arrival-cost form
cannot represent a value function whose slope at the solution is a nonzero
multiplier. The roadmap also never says what `xhat` is or where the loop gets
it, which matters most in exactly this case.

**1.4 It is a conditional, not the promised marginal.** The "Marginal versus
conditional" section is emphatic that reducing onto a block integrates
everything else out and that the answer does not depend on what else was
declared. The bounds section then projects bound-active variables out. Those
are incompatible: a projected-out variable is held fixed, not integrated out.

Quantified on the study's coupled 2x2: information on x2 conditional on x1
pinned is `W22 = 2.0`, where the bound-free marginal is
`1/(W^{-1})_22 = 1.5`. So the projection reports 33% more information (25%
smaller variance) on the free block than a marginal would. In MHE this is not
an edge case, since interior window states sit on bounds routinely
(nonnegative concentrations, saturated valves, envelope limits), and the
over-confidence is wrong the moment the bound releases in a later window.

**1.5 Active-set churn is unaddressed.** `Pi^{-1}` changes rank window to
window as bounds activate and release, and near a degenerate bound the
classification itself is unstable. The roadmap does not say what the loop does
when the free block changes dimension between windows.

**1.6 No validation covers any of this.** There is no MHE validation item with
an active bound, at the arrival state or at an interior state.

## 2. The barrier regimes

**2.1 The third regime is missing.** The bounds section has two regimes,
strongly active (`Sigma ~ z^2/mu`, unbounded) and inactive (`Sigma ~ mu/s^2`,
vanishing). The study's Section 1.3 has three. The missing one is weakly
active, where slack and multiplier vanish together: `s = sqrt(mu/q)`,
`z = sqrt(mu*q)`, and

    Sigma = mu/x^2 = mu/(mu/q) = q, exactly, at every mu.

Recomputed independently at q=2.7 across mu from 1e-2 to 1e-12: `Sigma` printed
2.700000 every time, and a free block containing that direction carries
`q + Sigma = 5.400000` every time. The study's own diagnostic output agrees:
"WEAKLY ACTIVE (p=0), mu=1e-08 ... barrier_c 1.00e+00  H_ii 1.00e+00  ratio
1.00e+00".

So the barrier curvature is finite, mu-independent, and of the same order as
the objective's own curvature. It sits on a direction whose multiplier is near
zero, which the projection has no clear claim on.

**2.2 The central claim fails there.** "On the free space they agree to
`O(mu)`" is false in that regime: classified free, the block carries `2q`
instead of `q`, so information is doubled and variance halved, silently, and
the error does not shrink as `mu` falls. Classified pinned, the direction is
dropped entirely and the next window gets no information on a direction the
bound is not actually holding. Neither of the roadmap's two dispositions
returns the correct `q`. The study's position is that the object is set-valued
there and must be flagged, not silently included or silently deleted.

**2.3 The mu-independence check has a false pass.** The roadmap's validation
bullet says refining `mu` must not move the free-block numbers, "which is what
separates information from barrier curvature." In the weakly-active regime the
contaminated number is exactly mu-invariant, so the check passes while the
block is 2x wrong. It is a necessary check, not a certificate.

## 3. Activity detection

**3.1 No criterion is given.** The whole bounds section rests on knowing which
directions are pinned, and the roadmap never states the rule.

**3.2 The rule that ships is the one the study forbids.** `covariance()`
classifies on slack alone, with no multiplier test (`sens.py:826-827`):

    tol = 1e-6 * (1.0 + abs(xv))
    if xv - lo[r] < tol or hi[r] - xv < tol:

The study's Section 5 is explicit that classification must use both slack and
multiplier, "never distance to the bound alone." Its own Ipopt run at the
weakly-active point returns `x = 5.64e-7` with `z_L = 5.64e-7`. That slack
falls under the shipped 1e-6 tolerance, so the shipped rule projects out a
direction whose information is finite. The study also shows its recommended
classifier returning "ambiguous" for a cleanly inactive case and a cleanly
strongly-active one at `mu=1e-8`, so tolerances need tying to `mu` (compare `s`
to `sqrt(mu)` and `s*z` to `mu`) rather than fixed constants.

This is a defect in shipped code, not only in the roadmap. The roadmap inherits
it by saying the projection is shared.

## 4. Scope and specification gaps

**4.1 Inertia correction contaminates the free block.** The bounds section
argues no scrubbing is needed because `Sigma` lives on the deleted direction.
That argument is `Sigma`-specific. The study flags a second solver artifact:
inertia correction adds `delta_w * I`, which is isotropic, lands on the free
block, and survives the projection. This is not hypothetical here.
`solver.rs` documents `kkt_perturbations()` as inertia-correction perturbations
baked into the held KKT factor, and `covariance()` already guards on it
(`sens.py:814-820`). The sharp point: `delta` is injected precisely when the
Hessian is indefinite or near-singular, which is both the "Lagrangian can go
indefinite" case the roadmap names and the "poorly-identified directions" where
item 1 claims `information()`'s headline advantage. So the roadmap advertises
better scaling in exactly the regime where the retained factor is least likely
to be the exact KKT matrix, and never says `information()` carries the same
guardrail.

**4.2 Bound-active variables outside the block are unspecified.** The bounds
section is written entirely about bound-active members of the block. For a
bound-active variable outside it, nothing deletes it: its `Sigma ~ z^2/mu`
lives in the held factor and enters the Schur complement onto the block. The
study's Section 4 is exactly this configuration. Working its numbers, the
reduced information on x2 is `W22 - W21 (W11+Sigma1)^{-1} W12`, which is 1.9938
at `mu=1e-2` and 2.0000 at `mu=1e-8`, against the true marginal 1.5. Two
consequences the roadmap never states: `information(wrt=x2)` is not the
marginal its own semantics section promises, and this free-block number does
move with `mu`, so the mu-independence check would flag a correct
implementation.

**4.3 Item 1 defines `information()` twice, incompatibly.** It is "the scaled
reduced Hessian" and also numerically `inv(covariance(...))`. Those coincide
only in the homoscedastic case. `covariance()` returns a heteroscedastic
sandwich whenever labeled residual groups have unequal variances
(`sens.py:995` Lagrangian, `sens.py:969` Gauss-Newton), and inverting a
sandwich is not a scaled reduced Hessian by any scalar. So on a grouped model
the two definitions give different matrices and validation bullet 1 fails for a
correct implementation of either. "Well-conditioned block" does not cover this.

A quieter version sits between sections: the benefit hypothesis promises an
accessor consistent with both `covariance()` and the core's natural-units
reduced Hessian, but those are different conventions. The core's is unscaled;
pyomo's `covariance()` carries `2 * sigma_sq`. The MHE consumer needs the
unscaled one, since its objective already carries its own inverse-covariance
weights.

**4.4 Preconditions are unstated.** The study closes with an explicit
assumptions list; the roadmap has none. Two bite. LICQ is what makes the
active-set KKT nonsingular and is the precondition for the roadmap's own
validation test (the free block matching the model solved with that variable
fixed, which is a bounds-to-equalities substitution). Second-order sufficiency
matters because the roadmap makes Lagrangian the default for `information()`
while noting elsewhere that the Lagrangian can go indefinite. An indefinite
`Pi^{-1}` makes the next window's NLP unbounded below along the
negative-curvature direction. The roadmap never says what `information()` does
when the free block is not positive definite.

## What this implies

The bounds section can be repaired in place: add the third regime, state the
activity criterion, scope the claims to in-block variables, demote the
mu check, and note the inertia-correction guardrail.

The MHE section cannot. Its recursion is wrong, its bound argument is wrong,
and the object it wants to hand forward has a linear term the arrival-cost form
it writes cannot carry. That section wants rewriting around the one-step
recursion, and the question of what to do at a bound-active arrival state wants
deciding rather than patching.

Two findings are about shipped code rather than the roadmap: the slack-only
activity test (3.2), and whether `information()` inherits the
`kkt_perturbations` guardrail (4.1). Both want their own issues.
