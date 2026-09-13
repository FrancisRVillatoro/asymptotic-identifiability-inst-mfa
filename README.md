# Asymptotic identifiability — reproducibility code and data

This repository contains **code, data, machine-readable numerical outputs, audit logs, and reproducible figure files** supporting the study

**Asymptotic Identifiability Orders and Order Lifting by Experimental Design in Isotopically Nonstationary 13C Metabolic Flux Analysis**

by Francisco R. Villatoro.

The manuscript itself is **not** part of this repository. Paper sources and article PDFs are distributed only through the journal/preprint channel.

## Scientific reproduction

The compute-only scientific pipeline is

```bash
bash reproduce_science.sh
```

It reruns:

1. FreeFlux source-faithful runtime validation and sensitivity refinement;
2. the network-based synthetic Synechocystis sensitivity/order analyses;
3. the fast/slow partition sweep and diagnostics;
4. the finite-noise/profile calculations;
5. the Python finite-noise figure refresh.

For provenance capture, use

```bash
bash run_full_science_audit.sh
```

which records the software environment, stdout/stderr, timing, before/after SHA-256 manifests and a machine-readable completion record under `audit/full_runs/<UTC timestamp>/`.

## Tested Picasso environment

The final post-patch end-to-end scientific rerun (job 2279016) completed with exit code 0 on the UMA SCBI Picasso cluster using:

- Python 3.9.13
- NumPy 1.24.4
- SciPy 1.9.1
- pandas 1.4.4
- Matplotlib 3.5.2
- openpyxl 3.0.10
- SymPy 1.10.1

The runtime compatibility layer also supports newer pandas by using `DataFrame.map` when available and `DataFrame.applymap` on pandas 1.4.x.

## Repository contents

- `data/`: frozen numerical inputs and summary data;
- `validation/freeflux_runtime/`: Block 1 FreeFlux validation;
- `validation/synechocystis_block2/`: network-based synthetic Synechocystis benchmark;
- `validation/partition_block3/`: partition sweep and weak-subspace diagnostics;
- `validation/block4_finite_noise/`: finite-noise and profile calculations;
- `validation/block5_feasibility/`: experimental-feasibility support;
- `figure_sources/`: Python sources for reproducible scientific figures;
- `figures/`: rendered reproducible figures;
- `scripts/`: repository/audit helper scripts;
- `audit/`: provenance manifests and final-run evidence.

## Reproducibility policy

The canonical v1.0.0 release baseline is the **post-patch successful full run**, not the older pre-patch Block-6 package. The release-assembly script `scripts/assemble_postpatch_release.py` combines this code tree with the archived scientific snapshot from job 2279016, freezes the new canonical SHA-256 manifest, verifies that no manuscript files are present, and creates the final code-and-data archive.

No LaTeX compilation is part of the scientific reproduction workflow.

## Strict content policy

This repository deliberately excludes manuscript files and manuscript-build
artifacts. In particular, it contains no `.tex`, `.bib`, `.bbl`, article PDF,
Biber/LaTeX build logs, or manuscript-verification scripts. Rendered scientific
figures are included only as reproducibility outputs; their scientific source
code lives under `figure_sources/` or the corresponding validation directory.
Run `python scripts/guard_no_paper_artifacts.py` to verify this policy.
