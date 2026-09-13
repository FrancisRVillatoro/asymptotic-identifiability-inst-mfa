# Block 1 - FreeFlux 0.3.8 runtime/source validation

## Scope and provenance

The audit is fixed to FreeFlux 0.3.8 and Git commit
`ec05c47bbc2e4ac58bb85d39408ff4ef4016a15c`.  The current execution
environment is Python 3.13, whereas the package metadata advertises Python
3.7-3.10 and pins legacy NumPy/SciPy ranges.  Therefore the simulation path
was executed through a source-faithful extraction of the official EMU
implementation.  The only compatibility changes were:

- `np.float` -> `float`;
- `scipy.linalg.pinv2` -> `scipy.linalg.pinv`;
- omission of fitting, optimization and result-display layers not used by
  INST simulation.

The EMU decomposition, equivalent-EMU handling, construction of the A, B and
M matrices, natural-abundance model, tracer MDVs, initialization and the
official piecewise-linear INST propagator are retained.

This is stronger than a reimplementation from documentation but is not
claimed to be an untouched-wheel execution.  A Python-3.10 runner and pinned
requirements are included so that the unmodified wheel can be executed in a
supported environment.

## Benchmark reproduction

The source-faithful runtime used the official toy reaction file, fluxes,
concentrations, AcCoA tracer strategy, target EMUs `Glu_123`, `Glu_12345` and
`Cit_2345`, and times 0, 0.1, 0.2, 0.5, 1 and 2.

The simulated `Glu_12345` values on the six-time tutorial grid round exactly
to every value in the official `measured_inst_MDVs.tsv` fixture.

## Comparison with exact full-isotopomer propagation

FreeFlux advances each EMU block over each user-supplied interval using a
piecewise-linear approximation to the lower-size forcing.  Consequently the
six-time tutorial grid is not identical to exact propagation of the coupled
full-isotopomer linear system.  Subdivision of the same intervals removes
this discretization difference quadratically:

| internal maximum step | maximum error | RMS error | relative Frobenius |
|---:|---:|---:|---:|
| 1 | 4.583550e-03 | 7.834499e-04 | 2.118305e-03 |
| 0.1 | 1.385403e-04 | 3.068939e-05 | 8.297850e-05 |
| 0.05 | 3.415810e-05 | 7.579971e-06 | 2.049486e-05 |
| 0.025 | 8.534241e-06 | 1.895068e-06 | 5.123919e-06 |
| 0.0125 | 2.133477e-06 | 4.738333e-07 | 1.281159e-06 |
| 0.00625 | 5.333545e-07 | 1.184604e-07 | 3.202952e-07 |

The fitted orders over the four finest grids are 2.000302 for the
maximum error and 1.999895 for the RMS error.

## Sensitivity comparison

Central finite differences in the nine log-parameters were evaluated with
identical parameterization in both engines.  At internal maximum step
0.00625, the source-faithful runtime and exact full-isotopomer sensitivities
satisfy:

- global maximum absolute difference: 1.156630e-06;
- global RMS difference: 1.137860e-07;
- relative Frobenius difference: 3.573638e-06;
- worst relative column error: 1.766312e-05
  for `log v6b`.

Halving the internal step from 0.0125 to 0.00625 gives observed orders
2.000007, 1.999932 and 1.999932 for the
maximum, RMS and relative-Frobenius sensitivity errors, respectively.

## Conclusion

The official rounded tutorial fixture is reproduced exactly.  The remaining
runtime-versus-exact differences are a second-order time-grid effect of the
FreeFlux piecewise-linear EMU propagator, not a disagreement in atom mapping,
natural abundance, flux conventions, pool dilution, or the terminal Glu
balance.  The same conclusion holds for the nine-parameter sensitivity
matrix under grid refinement.
