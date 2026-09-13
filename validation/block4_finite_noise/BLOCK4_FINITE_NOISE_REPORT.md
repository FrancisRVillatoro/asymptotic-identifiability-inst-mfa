# Block 4 — finite-noise compositional validation

## Purpose

This block tests whether the asymptotic order hierarchy survives a finite-noise observation model and whether order lifting produces measurable contraction of uncertainty. It deliberately separates the algebraic order m, defined for an epsilon-independent weighting with a nonsingular limit, from finite-noise effective slopes when the likelihood-induced weighting depends on the simulated composition.

## Statistical observation model

For each complete MID p in the simplex, the observation is transformed to pseudocount-regularized ILR coordinates

    z(p) = H log(p + delta),

where H is a Helmert contrast matrix and delta = 1e-4. We use

    z_obs = z_true + eta,   eta ~ N(0, sigma_ILR^2 I),
    sigma_ILR = 0.08.

Absolute pool measurements are modeled as

    log C_obs = log C_true + xi,   xi ~ N(0, 0.1^2).

The random seed is 20260905.

The delta-method raw-space standard deviations induced by sigma_ILR=0.08 in the Synechocystis benchmark have median 8.06e-3, 75th percentile 1.85e-2, 95th percentile 2.42e-2 and maximum 2.56e-2. These values lie in the same range as example component uncertainties 0.005, 0.013 and 0.028 distributed with the published Synechocystis FluxML data. Thus the statistical model is compositional and heteroscedastic in raw MID space, rather than Euclidean with one arbitrarily deleted component.

## FreeFlux toy benchmark

At the published point, after quotienting the exact structural fibers, the finite-noise local standard errors along the two asymptotic weak directions are:

| design | SE(q1, m=1) | SE(q2, baseline m=2) |
|---|---:|---:|
| baseline | 1.0877 | 10.6629 |
| + Fum at slow times | 0.9625 | 0.7534 |
| + Fum + Asp-fast | 0.2041 | 0.1682 |
| + above + absolute AKG pool | 0.2041 | 0.1682 |

Hence the q2 uncertainty contracts by a factor 14.15 from baseline to the slow-Fum design and by a factor 63.39 from baseline to the fast-lifted design.

A 4000-replicate local Monte Carlo in transformed observation space agrees with the linear covariance calculation. For q2 the RMSE changes

    10.6719 -> 0.7603 -> 0.1706,

with empirical 95% coverage 0.947, 0.951 and 0.948, respectively.

For one perturbed data realization, a nonlinear profile along the baseline q2 direction (with Gauss-Newton nuisance profiling) does not cross Delta chi^2=3.84 on [-1.5,1.5] under the baseline design. Adding slow Fum gives an approximate 95% interval

    [-1.4334, 0.9032],  width 2.3366,

and adding the fast Asp experiment contracts this to

    [-0.3346, 0.2954],  width 0.6300.

This is direct finite-noise evidence that lowering the order can turn a broad nonlinear valley into a substantially narrower profile.

## Synechocystis benchmark

The same ILR likelihood was applied to the automatic F6P-GAP fast partition. At epsilon=1/64 the labeling-only problem remains globally very sloppy even after the exact scale is quotiented. The smallest singular values are

    baseline: 4.90e-6,
    + C_F6P: 2.11e-5,
    + C_GAP: 2.16e-5,
    + both pools: 3.07e-4.

The large remaining covariance modes are not a contradiction of order lifting: m describes the power of epsilon, not the O(1) prefactor or unrelated sloppy combinations.

A 5000-replicate local Monte Carlo makes the pool-specific effect explicit. At epsilon=1/64:

- baseline: RMSE(log C_F6P)=1.45e5 and RMSE(log C_GAP)=1.49e5;
- +C_F6P: RMSE(log C_F6P)=0.1009 while GAP remains at 4.77e4;
- +C_GAP: RMSE(log C_GAP)=0.0996 while F6P remains at 4.51e4;
- +both: both RMSEs are approximately 0.10, i.e. the imposed absolute-pool measurement scale.

The baseline weak singular vector is almost exactly the fast-pool contrast: its two largest coefficients are log C_GAP = -0.7162 and log C_F6P = +0.6977.

Profiling the physically transparent contrast

    a = (log C_F6P - log C_GAP)/sqrt(2)

shows that baseline, +F6P and +GAP do not reach Delta chi^2=3.84 over [-1.5,1.5], whereas measuring both pools gives

    [-0.19384, 0.18583], width 0.37967.

Thus both independent pool measurements are required to close this two-dimensional first-order fast sector, exactly as predicted by the order-lifting calculation.

## Epsilon-dependent weighting

Because the ILR Jacobian contains diag(1/(p+delta)), the statistical weighting varies with epsilon as the MIDs move along the homotopy. The resulting finite-noise slopes are therefore not the formal Smith orders and are reported only as effective preasymptotic slopes.

For Synechocystis, over epsilon <= 1/4, the two baseline weak singular values have effective slopes approximately 1.026 and 0.755. A single pool measurement leaves one weak effective slope of 0.813 (+F6P) or 0.860 (+GAP), while measuring both gives slopes 0.0046 and 0.0069 for the two weakest modes, i.e. effectively O(1). The design hierarchy survives even though the numerical slopes are no longer quantized exactly.

The toy benchmark displays the same phenomenon: under ILR weighting the weakest baseline preasymptotic slopes over the tested range are approximately 1.68 and 0.78 rather than the formal 2 and 1. This is precisely why the manuscript distinguishes asymptotic identifiability order from a finite-noise, composition-dependent Fisher slope.

## Scope and caveat

The Synechocystis finite-noise experiment here intentionally isolates the labeling likelihood plus candidate pool measurements after structural quotienting. The official FreeFlux synthetic fitting tutorial also supplies measured fluxes, so the extremely large nuisance uncertainties reported here must not be interpreted as uncertainty estimates for the complete FreeFlux fit. They are a controlled test of whether the fast-pool information loss and its repair survive a realistic compositional observation model.

## Conclusion

Block 4 supports three claims.

1. The asymptotic design ranking survives finite, heteroscedastic, compositional MID noise.
2. In the toy benchmark, m=2 -> 1 -> 0 lifting produces large local, Monte-Carlo and nonlinear-profile contraction.
3. In the larger benchmark, lifting the fast-pool sector pins the targeted pools but does not automatically cure unrelated O(1) sloppiness; information order is therefore a hierarchy of asymptotic degeneracy, not a complete uncertainty metric.

A full Bayesian posterior analysis remains optional strengthening work rather than a blocking requirement for the present manuscript.
