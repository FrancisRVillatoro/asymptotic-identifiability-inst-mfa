# Block 2 — Synechocystis INST-13C-MFA benchmark

## Benchmark dimensions

- 38 net reactions
- 60 directional fluxes
- 31 balanced intracellular pools
- rank(S_total) = 31
- 29 steady-state flux degrees of freedom
- 31 log-pool parameters
- joint local parameter dimension = 60
- 215 independent MID observations after removing one compositional component per MID and t=0 rows

The distributed synthetic flux vector has a small GAP steady-state residual 0.0326324. Projection to the nearest steady-state directional-flux vector changes the flux vector by relative norm 8.44e-5 and leaves all directional fluxes positive.

## Nominal identifiability

The 215 x 60 joint sensitivity matrix has numerical rank 59. The single structural null is the global simultaneous flux/pool scale symmetry.

## Fast-pool control

The raw top turnover group contains GAP, TA and F6P. Because TA is an atom-mapping auxiliary pool, a biochemical-pool control was performed with only F6P and GAP.

Homotopy:
C_F6P(epsilon)=epsilon C_F6P(1),
C_GAP(epsilon)=epsilon C_GAP(1),
with all fluxes and other pools fixed.

## Asymptotic result

For the biochemical fast2 partition:
- sigma_58 exponent = 0.999381
- sigma_59 exponent = 1.002753
- sigma_57 exponent = 0.000072

Thus there are exactly two m=1 directions, giving two Fisher eigenvalues Theta(epsilon^2), and no m=2 / Theta(epsilon^4) direction.

The fast3 control {F6P, TA, GAP} gives three m=1 singular values with exponents 0.9883, 0.9992 and 1.0034, again with no m=2 direction.

## Weak-direction coordinates

At epsilon=1/512 the two m=1 directions project onto (log C_F6P, log C_GAP) as

[[ 0.65839593 -0.75260257]
 [-0.7453093  -0.65003417]]

with det = -0.988901545 and cond_2 = 1.011220037.

## Experimental lifting

After quotienting the exact global scale symmetry:
- baseline: two m=1 directions;
- + log C_F6P: one m=1 direction;
- + log C_GAP: one m=1 direction;
- + both: all 59 identifiable directions are m=0.

Weakest-singular-value asymptotic exponents:
- baseline: 1.002753 (the other weak mode: 0.999381)
- +F6P: 0.999475
- +GAP: 1.000498
- +F6P+GAP: 0.000071

## Scientific conclusion

This second benchmark does not reproduce the m=2 phenomenon of the FreeFlux toy TCA network. Instead it shows a robust first-order hierarchy: each selected fast biochemical pool contributes one m=1 direction, while the remaining identifiable directions stay O(1). This demonstrates that m=2 is not automatic; it requires an additional moment degeneracy such as the OAA-Fum first-moment cancellation in the toy benchmark.
