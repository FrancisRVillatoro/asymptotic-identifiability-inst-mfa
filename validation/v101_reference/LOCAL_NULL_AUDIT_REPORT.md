# Synechocystis local-null audit after genericization and exact-J correction

## Question

Does the genericized 60-parameter benchmark show any additional **local continuous null direction** beyond the known global flux--pool scale symmetry?

This audit does **not** claim global structural identifiability and does not exclude isolated/discrete global symmetries. It addresses the local statement needed by the manuscript.

## Exact-J nominal checks

The discrete-consistent Frechet Jacobian is 215 x 60. At the genericized point `C_PG=11`, `C_Gc=9`, the numerical rank is 59 at three widely separated homotopy locations:

| epsilon | rank | second-smallest sigma | null sigma | |<v_null,q_scale>| | min sigma after quotienting scale |
|---:|---:|---:|---:|---:|---:|
| 1 | 59 | 3.110264e-07 | 2.071362e-14 | 1.000000000000000 | 3.110264e-07 |
| 0.015625 | 59 | 1.723641e-08 | 2.719839e-14 | 0.999999999999999 | 1.723641e-08 |
| 0.001953125 | 59 | 2.158885e-09 | 2.001978e-14 | 0.999999999999967 | 2.158885e-09 |


The residual of the analytic scale direction, normalized by the operator norm of J, is at most `3.585e-14`. Thus the only numerical null at the frozen generic point is the known scale direction.

## Independent nearby-point check

A deterministic random perturbation was applied to **all 31 pool concentrations**, with maximum absolute log shift `0.050` (about 5%), while retaining `C_PG != C_Gc`. The observed vector changes by `4.084e-03` in relative 2-norm. An independent exact discrete Jacobian was then recomputed at that perturbed point.

Result:

- rank = `59`;
- second-smallest singular value = `3.491001e-07`;
- null singular value = `1.331310e-14`;
- alignment with the analytic scale direction = `1.000000000000000`;
- relative scale residual = `1.969e-14`;
- smallest singular value after quotienting scale = `3.491001e-07`.

Hence the rank-59 result is not an isolated consequence of the exact nominal concentrations.

## Supporting PG/Gc genericity scan

Step 2 independently varied the genericization to

`C_PG=10(1+d), C_Gc=10(1-d)` with `d=0.08, 0.10, 0.12`

at epsilon `1/64` and `1/512`. Across this scan the two F6P/GAP weak singular values and their fast-subspace projections remain essentially unchanged. This supports that the selected `d=0.10` point is not a special fine-tuned replacement point.

## Interpretation

A second smooth structural fiber would generate an additional local null direction throughout a neighborhood. The exact-J rank-59 result at the generic point, its persistence under an independent nearby concentration perturbation, and the PG/Gc off-equality robustness scan provide strong numerical evidence that **no additional local continuous null direction was found** beyond global scale.

The correct manuscript wording is therefore:

> "At the genericized nominal point, the 215 x 60 sensitivity matrix has numerical rank 59. The single numerical null is aligned with the known global flux--pool scale symmetry; no additional local continuous null direction was detected in nearby-point audits. This is a local numerical statement, not a proof of global structural identifiability."

## Verdict

- `GENERIC POINT RANK 59 AT eps=1,1/64,1/512 = PASS`
- `NULL ALIGNED WITH GLOBAL SCALE = PASS`
- `INDEPENDENT NEARBY 5% CONCENTRATION PERTURBATION RANK 59 = PASS`
- `PG/Gc GENERICIZATION ROBUSTNESS d=0.08..0.12 = PASS`
- `GLOBAL STRUCTURAL IDENTIFIABILITY CLAIM = NOT MADE`

**The Synechocystis local-null audit is closed.**

The remaining pre-release task is integration of Steps 1--3 and the final FreeFlux profile corrections into the manuscript/reproducibility tree, followed by a clean end-to-end reproduction check before freezing `v1.0.1`.
