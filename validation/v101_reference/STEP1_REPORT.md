# Step 1 — Discrete-consistent Synechocystis Jacobian and 16→32→64 audit

## Scope

This package closes Step 1 of the post-referee numerical audit. It replaces the previous approximate sensitivity propagation by a Jacobian that is the exact derivative of the same piecewise-linear discrete propagator used for the forward map, up to floating-point error.

The original forward step is

\[
X_1=\Phi X_0-\Gamma G_0-\Omega(G_1-G_0),
\]
with
\[
\Phi=e^{F\Delta t},\qquad
\Gamma=(\Phi-I)F^{-1},\qquad
\Omega=(\Gamma/\Delta t-I)F^{-1}.
\]
For each parameter direction, the new code uses

\[
D\Phi=L_{\exp}(F\Delta t,DF\,\Delta t),
\]
computed by `scipy.linalg.expm_frechet`, together with

\[
D(M^{-1})=-M^{-1}(DM)M^{-1},\qquad
D(F^{-1})=-F^{-1}(DF)F^{-1},
\]
\[
D\Gamma=D\Phi F^{-1}+(\Phi-I)D(F^{-1}),
\]
\[
D\Omega=(D\Gamma/\Delta t)F^{-1}+(\Gamma/\Delta t-I)D(F^{-1}).
\]
The source-EMU derivatives are propagated recursively through the same convolution graph as the nominal forward solution.

All dense linear-algebra runs used one BLAS thread (`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) to avoid oversubscription artifacts.

## 1. Same discrete forward map

At ε=1/64 with 16 internal subdivisions per experimental interval, the refactored forward map agrees with the original source-faithful implementation to

- relative 2-norm difference: `6.31e-13`
- maximum absolute difference: `9.30e-13`.

Thus the Fréchet Jacobian differentiates the intended discrete propagator, not a replacement model.

## 2. Independent derivative check

A complex-compatible algebraic copy of the same propagator was used only as an independent directional derivative check. At ε=1/64:

- random normalized directions agree with the Fréchet Jacobian at roughly `1e-12` relative error;
- the two principal weak F6P–GAP modes agree at `1e-5` to `1e-6` relative error, with cosine equal to 1 to displayed precision;
- the scale and PG–Gc null directions are both at numerical-zero level, so relative errors there are not meaningful.

The centered-difference audit is also included. As expected, finite subtraction loses accuracy in the weakest cancellation mode even when it works well on ordinary columns; this is why the Fréchet/complex-step pair is the relevant certification.

## 3. 16→32→64 convergence

The exact discrete Jacobian was recomputed at subdivisions 16, 32 and 64 for ε=1, 1/64 and 1/512.

Global convergence is essentially second order at all three ε values:

| ε | forward order | Jacobian order | rel. `J32-J64` / `J64` |
|---:|---:|---:|---:|
| 1 | 1.979 | 1.945 | 4.62e-4 |
| 1/64 | 1.979 | 1.945 | 4.52e-4 |
| 1/512 | 1.979 | 1.945 | 4.52e-4 |

For ε=1/64, the two F6P–GAP weak singular values (after the two exact first-order null directions discussed below) converge as

| mode | sub16 | sub32 | sub64 | p=2 Richardson | rel. sub64→Richardson |
|---|---:|---:|---:|---:|---:|
| weaker | 1.8309e-8 | 1.7565e-8 | 1.7339e-8 | 1.7263e-8 | 4.38e-3 |
| stronger | 1.9192e-7 | 1.9350e-7 | 1.9305e-7 | 1.9290e-7 | 7.81e-4 |

The corresponding two-dimensional right-singular subspace changes by only about `0.039°` and `0.013°` from sub32 to sub64 at ε=1/64.

Hence 64 subdivisions are adequate for the Step-1 certification at the sub-percent level in the weakest resolved mode. Any final exponent fit should nevertheless be based on the validated Jacobian, not on the old approximate sensitivity propagation.

## 4. New finding: a second exact first-order null at the distributed nominal point

The discrete-consistent Jacobian has numerical rank 58, not 59, at the distributed nominal concentration point. The two null directions are:

1. the known global flux–pool scale direction;
2. the antisymmetric pool direction

\[
q_{PG-Gc}=\frac{1}{\sqrt2}(e_{\log C_{PG}}-e_{\log C_{Gc}}).
\]

The alignment of the two smallest right singular vectors with these analytical directions is 1.0 to displayed precision at ε=1, 1/64 and 1/512, for subdivisions 16, 32 and 64.

This second null is not a continuous scale fiber. It is caused by the exact discrete symmetry of the serial PG→Gc block at the distributed nominal point, where

\[
C_{PG}=C_{Gc}=10
\]

and the steady-state fluxes through `rbc2`, `pgp` and `gld` coincide. The input–output map is invariant under interchange of the two pool sizes. Direct swap tests give relative output differences of `4e-16`–`1e-15` for examples `(PG,Gc)=(9,11),(5,20),(1,30)` and their swaps.

Moving away from the symmetric point lifts this first-order singularity. At subdivisions 16 and ε=1, multiplying only `C_PG` by 1.001, 1.01 and 1.1 raises the formerly null singular value to approximately

`1.52e-9`, `1.50e-8`, and `1.29e-7`, respectively.

This must be addressed before the final Synechocystis benchmark is frozen. The clean options are either (a) move the synthetic nominal point slightly off the artificial equality `C_PG=C_Gc`, then rerun the benchmark, or (b) retain the symmetric point but explicitly classify the PG/Gc branch-coalescence singularity separately from structural fibers and from the fast–slow Smith orders. Option (a) is likely cleaner for a generic network benchmark.

## 5. Consequence for the previous referee concern

Step 1 now distinguishes two separate errors that were conflated before:

- the old sensitivity propagation was not the exact derivative of the discrete forward step;
- the forward/sensitivity discretization also required a time-resolution convergence check.

The new Fréchet Jacobian removes the first issue exactly. The 16→32→64 study then addresses the second issue directly and shows near-second-order convergence.

## Files

- `scripts/syn_discrete_frechet_jacobian.py`: exact derivative of the discrete propagator.
- `scripts/syn_discrete_complex_step.py`: independent complex-step directional validator.
- `scripts/run_exact_jacobian.py`: single-point runner.
- `scripts/analyze_convergence.py`: regenerates the convergence summary from the NPZ results.
- `results/syn_exactJ_*.npz`: full 215×60 exact discrete Jacobians at ε=1, 1/64, 1/512 and subdivisions 16, 32, 64.
- `results/syn_exactJ_time_convergence_summary.csv`: time-refinement metrics and Richardson estimates.
- `results/frechet_complex_step_directional_validation.csv`: independent derivative validation.
- `results/pg_gc_symmetry_audit.csv`: scale/null alignments, PG↔Gc swap tests and genericity lift.
- `results/clean_smoke_match.json`: clean-package reproducibility check.

## Step-1 verdict

`DISCRETE-CONSISTENT JACOBIAN = PASS`

`16→32→64 TIME CONVERGENCE = PASS`

`OLD RANK-59 CLAIM AT THE DISTRIBUTED NOMINAL POINT = FAIL (rank 58 after exact differentiation)`

The next scientific step should therefore use this validated Jacobian and decide how to remove or explicitly treat the accidental PG/Gc symmetry before fitting the final Synechocystis asymptotic orders.
