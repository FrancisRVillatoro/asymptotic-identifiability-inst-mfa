# Final audit of the nonlinear FreeFlux finite-noise profiles

## Verdict

**PASS after numerical corrections to the reported crossing locations.**

The nonlinear reoptimization itself is supported. The earlier grid-interpolated
crossings were too coarse and should be replaced.

## Common physical domain

For audit purposes the common domain is the componentwise log-parameter box

`|theta_i-theta0_i| <= 2.5`

after fixing the exact structural gauge(s). None of the points supporting the
threshold claims reaches this physical boundary:

- Fum crossing neighborhood: max physical displacement = 2.016.
- Asp lower crossing: 0.493.
- Asp upper crossing: 0.385.
- Fum at a=2: 2.476.
- Baseline endpoints a=-2,+2: 1.353, 1.632.

Thus the original nuisance-coordinate box does not determine the threshold
crossings or the tested-range conclusions. The Fum point at a=2 is close to,
but still inside, the common physical boundary.

## Refined threshold crossings

Using direct nonlinear reoptimization near the threshold, enlarged nuisance
bounds, and two independent starting points at each refined point:

- Fum lower crossing: `a ~= -1.280247` -> report **-1.280**.
- Fum has no upper crossing through `a=2`; there `Delta chi2 = 0.356199`.
- Fum + Asp-fast lower crossing: `a ~= -0.403851` -> report **-0.404**.
- Fum + Asp-fast upper crossing: `a ~= 0.212244` -> report **0.212**.
- Corresponding two-sided width: **0.616**.

The two starting points at the refined threshold evaluations agree in chi-square
to approximately 1e-9 or better.

## Baseline range

At `a=-2` and `a=2` the fully nonlinear profile gives

- `Delta chi2 = 0.050934`,
- `Delta chi2 = 0.250630`,

well below 3.84, with both optima interior to the common physical domain.

## Optimizer diagnostics

All targeted optimizations returned `success=True` and status 2 (objective/ftol
termination). The first-order optimality measures are of order 1e-4--1e-3.
Although this is not machine-precision stationarity, the objective values are
stable under enlarged nuisance bounds and alternate starts to far better than
the precision needed for three-decimal crossing locations.

## Figure and table consequence

The old values `-1.32` and `[-0.39,0.21]` must be replaced by

- **Fum: lower crossing about -1.280; no upper crossing through 2.**
- **Fum + Asp-fast: about [-0.404, 0.212], width about 0.616.**

Do not use the Fum `a=-2` point, or Asp large-|a| points that were marked
`nuisance_hit_bound=True`, as evidence of the unconstrained profile. They should
be omitted or visibly marked as bound-limited in the plotting data.

## Scientific conclusion

The qualitative conclusion is unchanged and strengthened:

1. `m=2 -> m=1` produces a large local uncertainty contraction but does not
   yield a bounded two-sided profile over the tested range.
2. The subsequent `m=1 -> m=0` design produces a bounded two-sided nonlinear
   profile in this example.
3. Order lifting and practical identifiability remain distinct notions.
