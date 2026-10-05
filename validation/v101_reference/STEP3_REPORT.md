# Step 3 — Design lifting, scale anchoring and finite-noise audit on the genericized Synechocystis benchmark

## Inputs

This step uses **only** the discrete-consistent Fréchet Jacobians frozen in Step 2:

- `C_PG = 11`, `C_Gc = 9`;
- F6P/GAP homotopy `epsilon = 1, 1/2, ..., 1/512`;
- 64 internal subdivisions per experimental interval;
- full observed Jacobian `215 x 60` with the global flux-pool scale as the only exact null.

No result from the earlier approximate-J Synechocystis finite-noise campaign is reused.

## 1. Asymptotic design lifting

The direct log-pool rows for F6P and GAP are algebraic order-zero observations.  A log `co2in` flux row is used as an independent scale anchor.  Constant statistical weights do not alter the formal order count.

### Scale conditioned / known

| Design | m=0 | m=1 | m=2 |
|---|---:|---:|---:|
| labeling baseline | 57 | 2 | 0 |
| + F6P pool | 58 | 1 | 0 |
| + GAP pool | 58 | 1 | 0 |
| + F6P + GAP | **59** | **0** | **0** |

Thus each absolute fast-pool measurement removes one of the two first-order losses, and both measurements lift the complete 59-dimensional scale-quotiented problem to order zero.

### Full 60-parameter problem

| Design | m=0 | m=1 | m=2 |
|---|---:|---:|---:|
| + F6P pool | 58 | 2 | 0 |
| + GAP pool | 58 | 2 | 0 |
| + F6P + GAP | 59 | 1 | 0 |
| + F6P + GAP + independent scale anchor | **60** | **0** | **0** |

This is the rank accounting that was missing from the older benchmark interpretation.  In the full problem, a single pool row is first needed to break the exact global scale degeneracy, so one pool alone does not reduce the number of m=1 modes.  Two pool rows leave one first-order scale-dominated mode.  A third independent order-zero scale anchor removes it.

At epsilon=1/512 the remaining m=1 mode after `F6P+GAP` has slope `1.00009` over the deepest four points and overlap `0.98895` with the normalized global scale direction.  After adding the scale anchor, all 60 fitted four-point slopes are consistent with m=0; the largest absolute slope is about `3.4e-3`.

## 2. Fixed-pseudocount ILR weighting

The same classification is obtained after transforming every MID by the regularized ILR map with

- pseudocount `delta = 1e-4`,
- `sigma_ILR = 0.08`,
- log-pool uncertainty `sigma_pool = 0.10`,
- log-scale-anchor uncertainty `sigma_scale = 0.10`.

The weighted counts are identical to the raw counts.  This is a direct numerical check, on the revised benchmark, of the theorem that a fixed-pseudocount ILR map preserves asymptotic orders.

## 3. Finite-noise local uncertainty at epsilon=1/64

All covariance matrices are formed directly from the sensitivity SVD,

`Cov = V diag(sigma_i^{-2}) V^T`,

using the same rank decision as the Moore-Penrose inverse.  No Gram matrix `A^T A` is formed.

### Scale-conditioned comparison

| Design | SE(log F6P) | SE(log GAP) | SE(fast contrast) | SE(fast sum) |
|---|---:|---:|---:|---:|
| baseline | 3.922e4 | 3.769e4 | 5.374e4 | 8.436e3 |
| + F6P | 0.100 | 1.146e4 | 8.103e3 | 8.103e3 |
| + GAP | 1.193e4 | 0.100 | 8.433e3 | 8.433e3 |
| + F6P + GAP | **0.100** | **0.100** | **0.100** | **0.100** |

Thus one pool measurement constrains the measured coordinate but leaves the independent partner direction extremely uncertain.  Both measurements constrain the two targeted fast-pool coordinates at the assumed 10% log-scale precision.

A 5000-replicate matched linearized Gaussian Monte Carlo reproduces these local predictions (for example baseline RMSEs are 3.866e4 and 3.710e4 for F6P and GAP, and the two-pool RMSEs are 0.0998 and 0.0985).  As in the revised FreeFlux analysis, this Monte Carlo is an implementation check, **not** a nonlinear coverage study.

## 4. Practical effect of the scale anchor

In the full 60-parameter problem, `F6P+GAP` leaves one first-order mode that is 98.75% aligned with the global scale direction at epsilon=1/64.  Adding the `co2in` log-flux anchor converts that mode to m=0.

This is an **asymptotic** lifting statement, not a guarantee of good finite precision.  With a working 10% log-error on `co2in`, the local standard deviation of the global log-scale coefficient falls from about `5.91e3` without the anchor to about `23.7` with it.  The reduction is enormous but the remaining finite-noise uncertainty is still poor because order-zero directions with very small prefactors and strong correlations remain.

After both pool measurements and the scale anchor, the weakest finite-noise singular vector is no longer a fast-pool or scale mode: it is 99.96% aligned with the lifted PG-minus-Gc direction.  This is an explicit example of the paper's central caveat that `m=0` does not imply practical identifiability.

## 5. Correct regularized-ILR raw-space calibration

For `a = p + delta 1`, `s = 1 + n delta`, the inverse tangent derivative used for raw-MID calibration is

`G_delta = [diag(a) - a a^T / s] H^T`.

At the genericized nominal point the induced raw-component standard deviations have

- minimum: about `9e-6`,
- median: `8.342e-3`,
- maximum: `2.5611e-2`.

The old delta=0 approximation produced an artificial near-zero minimum (`~5e-12`).  The median and maximum change little, but the low-probability tail is now calibrated consistently with the pseudocount-regularized likelihood.

## Step-3 verdict

- `GENERICIZED EXACT-J DESIGN LIFTING = PASS`
- `CONDITIONED: 57+2(m=1) -> 58+1 -> 59+0 = PASS`
- `FULL: F6P+GAP LEAVES ONE m=1 SCALE MODE = PASS`
- `FULL: + INDEPENDENT SCALE ANCHOR -> 60 m=0 = PASS`
- `FIXED-PSEUDOCOUNT ILR PRESERVES ORDER COUNTS = PASS`
- `SVD-BASED FINITE-NOISE FAST-POOL CONTRACTION = PASS`
- `REGULARIZED ILR RAW-SPACE CALIBRATION = PASS`

## Consequence for the manuscript

All numerical Synechocystis values in the current v12 finite-noise section that came from the old approximate Jacobian must be replaced.  The robust qualitative statement is now sharper:

1. the genericized benchmark has two resolved m=1 F6P/GAP modes;
2. with scale conditioned, independent F6P and GAP pool measurements lift those two modes to m=0;
3. in the full 60-parameter problem the same two measurements leave one m=1 scale-dominated mode;
4. an additional independent order-zero scale anchor lifts that final mode, although finite-noise precision can remain poor because of small order-zero prefactors.

The next pre-release task is the focused FreeFlux nonlinear-profile audit (common physical domain, per-point optimizer diagnostics, and refined threshold crossings), followed by the Synechocystis local-null audit and integration into the clean v1.0.1 reproduction path.
