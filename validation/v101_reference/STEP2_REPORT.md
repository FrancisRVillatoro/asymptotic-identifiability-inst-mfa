# Step 2 — Genericized Synechocystis benchmark with the exact discrete Jacobian

## Objective

Remove the accidental PG/Gc branch-coalescence singularity found in Step 1 without materially changing the synthetic benchmark, then repeat the fast F6P/GAP asymptotic analysis with the discrete-consistent Fréchet Jacobian.

## Genericized nominal point

The distributed nominal point had `C_PG=C_Gc=10`, producing an exact swap symmetry and a second local null. I replaced only those two pools by

- `C_PG = 11`
- `C_Gc = 9`

leaving their arithmetic mean equal to 10. All fluxes, all other concentrations, tracer settings and observation times are unchanged.

This perturbation changes the nominal observed MID vector only by

- relative 2-norm: `1.85449e-6`
- maximum absolute component: `4.15705e-6`
- RMS component change: `6.00064e-7`.

Thus the benchmark is numerically almost unchanged at the observation level while moving off the accidental equality surface.

## Jacobian and time resolution

All reported final spectra use the exact derivative of the piecewise-linear discrete propagator developed in Step 1, with 64 internal subdivisions per experimental interval.

At the genericized point the full 215x60 Jacobian has numerical rank 59 at epsilon=1, 1/64 and 1/512. The only remaining numerical null is the known global flux-pool scale direction.

For the two fast F6P/GAP modes, a 16->32->64 check gives:

| epsilon | mode | sub16 | sub32 | sub64 | p=2 Richardson | rel. sub64->Richardson |
|---:|---|---:|---:|---:|---:|---:|
| 1/64 | weaker | 1.82500e-8 | 1.74785e-8 | 1.72364e-8 | 1.71557e-8 | 4.70e-3 |
| 1/64 | stronger | 1.89999e-7 | 1.92106e-7 | 1.91740e-7 | 1.91619e-7 | 6.36e-4 |
| 1/512 | weaker | 2.28408e-9 | 2.18880e-9 | 2.15888e-9 | 2.14891e-9 | 4.64e-3 |
| 1/512 | stronger | 2.39875e-8 | 2.42319e-8 | 2.41816e-8 | 2.41649e-8 | 6.94e-4 |

The sub32-to-sub64 weak-subspace angles are small; the individual weakest-vector angle is larger near epsilon=1/64 because of a nearby O(1) PG/Gc-derived mode, but the singular values themselves show sub-percent Richardson residuals at sub64.

## Asymptotic spectrum

The full sweep uses

`epsilon = 1, 1/2, ..., 1/512`.

Right singular vectors were tracked across epsilon by maximum-overlap matching, avoiding mode-order errors at crossings.

At the deepest point the two F6P/GAP modes are:

1. **contrast mode**
   - F6P/GAP projection: `0.999999`
   - normalized fast contrast `(F6P-GAP)/sqrt(2)`: `0.999993`
   - sigma at epsilon=1/512: `2.15888e-9`
   - fitted slope over deepest 4 points: `0.99905`

2. **sum mode**
   - F6P/GAP projection: `0.988814`
   - normalized fast sum `(F6P+GAP)/sqrt(2)`: `0.988808`
   - sigma at epsilon=1/512: `2.41816e-8`
   - fitted slope over deepest 4 points: `0.99593`

The corresponding 3/4/5-point slopes are:

- contrast: `0.99939 / 0.99905 / 0.99847`
- sum: `0.99861 / 0.99593 / 0.99986`.

Hence the resolved classification is

`m = 1, 1`

for the two fast-pool modes.

Among the other 57 scale-quotiented modes, the largest absolute 4-point slope is only `0.00372`; all classify as `m=0`. No resolved `m=2` mode is present in this validated genericized sweep.

## PG/Gc lifted mode

The former exact PG/Gc null becomes a genuine O(1) mode:

- sigma at epsilon=1/512: `3.58417e-7`
- projection on `(log C_PG-log C_Gc)/sqrt(2)`: `0.998679`
- 4-point slope: `1.23e-4`.

Thus it is cleanly separated from the two shrinking F6P/GAP modes in the asymptotic regime.

## Robustness to the amount of genericization

Using symmetric off-equality perturbations

`C_PG=10(1+d), C_Gc=10(1-d)`

with `d=0.08, 0.10, 0.12`, the two fast singular values and their F6P/GAP projections are essentially unchanged at epsilon=1/64 and 1/512. For example, the weakest fast sigma at epsilon=1/512 is `2.28407e-9` to the shown precision for all three choices.

Therefore the m=1,m=1 result is not an artifact of selecting d=0.10.

## Step-2 verdict

- `ACCIDENTAL PG/Gc NULL REMOVED = PASS`
- `ONLY GLOBAL SCALE NULL REMAINS = PASS`
- `EXACT DISCRETE JACOBIAN USED = PASS`
- `TWO F6P/GAP MODES m=1,1 = PASS`
- `NO RESOLVED m=2 MODE IN VALIDATED GENERICIZED SWEEP = PASS`
- `GENERICIZATION ROBUSTNESS d=0.08..0.12 = PASS`

The next step is to rerun design lifting / scale anchoring and finite-noise calculations on this genericized, exact-J benchmark. Those downstream calculations should not reuse the old rank-59 approximate Jacobian outputs.
