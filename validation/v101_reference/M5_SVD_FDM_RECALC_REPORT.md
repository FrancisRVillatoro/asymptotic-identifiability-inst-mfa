# M5 SVD/FDM recalculation audit

Date: 2026-09-16
Reference reproducibility release: GitHub/Zenodo v1.0.0
Purpose: resolve the numerical-method inconsistency identified in peer-review item M5.

## 1. SVD covariance convention

All local covariance calculations should use a single SVD of the weighted sensitivity matrix

    A = U diag(s) V^T

with cutoff

    tol = max(A.shape) * eps_machine * s_max.

For retained singular values, the Moore-Penrose inverse and covariance are

    A^+ = V diag(1/s) U^T,
    Cov = V diag(1/s^2) V^T.

This avoids forming A^T A and ensures that covariance and linearized Monte Carlo use the same numerical rank.

## 2. FreeFlux toy finite-noise calculation

The finite-noise toy calculation was reconstructed from the archived v1.0.0 equations and observation model. Centered finite differences were evaluated at

    h_theta = 8e-5, 4e-5, 2e-5, 1e-5, 5e-6.

The manuscript value h_theta = 2e-5 is in the stable plateau.

At h_theta = 2e-5, the SVD-based values are:

| design | sigma_min | SE(q1) | SE(q2) |
|---|---:|---:|---:|
| baseline | 0.09378294738 | 1.087726091 | 10.662919304 |
| +Fum | 1.022068887 | 0.9625280884 | 0.7534234592 |
| +Fum + Asp-fast | 3.314999894 | 0.2041054077 | 0.1682224529 |
| +Fum + Asp-fast + C_AKG | 1.202345448 | 0.2041054076 | 0.1682224529 |

The corresponding pinv(A^T A) values agree with the direct-SVD covariance to about 1e-11 relative or better in this small toy problem. Thus the toy local-SE numbers are unaffected by replacing the Gram pseudoinverse.

Finite-difference stability is much stronger than the precision displayed in the paper. Across the full h sweep, the relative variation of SE(q2) is about 1.1e-8 at baseline, 2.5e-9 after Fum, 7.6e-10 after Asp-fast, and 1.2e-9 after the pool measurement. The baseline q2 direction changes by less than 1.5e-6 degrees.

The step-doubling Jacobian differences exhibit the expected second-order regime before the roundoff plateau. From 8e-5 -> 4e-5 and 4e-5 -> 2e-5, the observed orders are approximately 1.94--1.95 for all four designs. Refinement below 2e-5 reaches the numerical floor, so h_theta = 2e-5 is an appropriate working value.

## 3. Toy linearized Monte Carlo

The linearized Monte Carlo was regenerated using the same SVD pseudoinverse and the archived seed 20260905, after reproducing the original RNG consumption by the shared rich-design noisy observation.

For q2 the reproducible values are:

| design | SVD local SE | linearized MC RMSE | nominal 95% coverage |
|---|---:|---:|---:|
| baseline | 10.6629193 | 10.4857676 | 0.95625 |
| +Fum | 0.75342346 | 0.75276173 | 0.95175 |
| +Fum + Asp-fast | 0.16822245 | 0.16455809 | 0.95200 |
| +pool | 0.16822245 | 0.17007411 | 0.94800 |

These exactly reproduce the archived v1.0.0 `toy_compositional_monte_carlo.csv` values. The manuscript v11 numbers 10.672, 0.760 and 0.171 are therefore stale and should not be retained. Since M1 already demotes this Monte Carlo to an implementation check, the cleanest main-text revision is to report the local SE contraction and leave the linearized Monte Carlo values to the supplement.

## 4. Synechocystis SVD covariance correction

The archived finite-noise script used SVD for singular values and `pinv(A)` for Monte Carlo, but used `pinv(A.T @ A)` for the reported modal covariance. This causes a rank inconsistency only in the most ill-conditioned baseline case.

Using SD_i = 1/sigma_i from the already archived direct SVD gives:

| design | corrected largest modal SD | corrected second modal SD | old Gram largest SD |
|---|---:|---:|---:|
| baseline | 204190.1285 | 33851.9169 | 33851.9187 |
| +F6P | 47493.5673 | 3254.81694 | 47493.6845 |
| +GAP | 46295.5935 | 3255.28671 | 46295.7286 |
| +F6P,+GAP | 3255.47280 | 2683.33207 | 3255.47274 |

Thus `pinv(A.T @ A)` dropped the weakest baseline singular mode. The direct SVD and `pinv(A)` calculations do not. The F6P, GAP and both-pool designs are unchanged to the displayed precision.

The archived F6P/GAP coordinate Monte Carlo values were already computed using `pinv(A)`, not the Gram matrix, so they do not inherit this dropped-mode error. In particular the baseline RMSE values 1.449e5 and 1.490e5 are consistent with the corrected very weak modal scale and are much larger than the erroneous old Gram `largest_sd`.

## 5. FreeFlux exact positional-isotopomer FDM audit

The exact positional-isotopomer response used in the source-faithful runtime comparison was independently differentiated with centered differences at

    h_theta = 4e-5, 2e-5, 1e-5, 5e-6, 2.5e-6.

Successive relative Frobenius differences are

    5.16e-10, 1.00e-9, 1.63e-9, 3.29e-9.

The non-monotone increase at the smallest steps is the expected floating-point floor; the entire variation is orders of magnitude below the archived runtime-vs-exact sensitivity discrepancy (3.5736e-6 at internal propagation step 0.00625). Thus the parameter finite-difference step h_theta = 1e-5 is not the limiting error in that cross-check. The source-faithful runtime itself was not rerun here at multiple parameter-FD steps; the archived time-grid refinement remains the direct check of its dominant propagation error.

## 6. Consequences for the revised paper

1. Replace every numerical covariance based on `pinv(A.T @ A)` by the SVD formula above.
2. State explicitly which sensitivity calculations use centered finite differences and which propagate differentiated EMU equations.
3. Retain the toy local SE values: they are unchanged.
4. Correct the Synechocystis baseline modal SD from 3.385e4 to 2.042e5 wherever that diagnostic is reported. The order-lifting claims and direct singular-value spectra are unchanged.
5. Replace/remove the stale manuscript Monte Carlo RMSE values. The archived/recomputed q2 RMSEs are 10.486, 0.753 and 0.165 for baseline, Fum and Asp-fast.
6. Keep the source-faithful FreeFlux sensitivity comparison 3.5736e-6; the exact-response FDM audit shows parameter-FD error is negligible relative to that propagation discrepancy.

## 7. Verdict

The recalculation confirms that the central scientific spectra and the toy finite-noise local uncertainty contraction are stable. One auxiliary Synechocystis covariance diagnostic was wrong because the Gram pseudoinverse discarded the weakest baseline mode; it is corrected by using the sensitivity SVD consistently. No asymptotic order changes.

M5 can be closed after these corrected outputs and method descriptions are integrated into the revised manuscript/supplement and the reproducibility repository is updated in a new post-v1.0.0 release.
