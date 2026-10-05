# Asymptotic identifiability — reproducibility code and data

This repository contains **code, data, machine-readable numerical outputs, audit records, and reproducible scientific figures** supporting the study

**Asymptotic Identifiability Orders and Order Lifting by Experimental Design in Isotopically Nonstationary 13C Metabolic Flux Analysis**

by Francisco R. Villatoro.

The manuscript itself is **not** part of this repository. Paper sources and article PDFs are distributed only through the journal/preprint channel.

## Version v1.0.1

Version `v1.0.1` incorporates the post-review scientific corrections verified by the clean Picasso end-to-end rerun (job `2570248`, run stamp `20261002T045228Z`):

- discrete-consistent Fréchet Jacobian for the Synechocystis benchmark;
- genericization from `C_PG=C_Gc=10` to `C_PG=11`, `C_Gc=9`;
- 16/32/64 internal-subdivision convergence checks and overlap-based mode tracking;
- scale-conditioned and full 60-parameter design lifting;
- regularized-ILR finite-noise calculations with direct-SVD covariance;
- nearby-point local-null audit;
- physical-point (`epsilon=1`) finite-noise control;
- final FreeFlux SVD/FDM and nonlinear-profile refinement audits;
- regenerated synthetic validation and automated checks against the frozen manuscript claims.

The successful run records `EXIT_CODE=0`, `COMPLETION_MARKER_FOUND=true`, `V15_CLAIM_AUDIT_PASS=true`, and `V101_CANDIDATE=VERIFIED`.

## Full scientific reproduction

Run

```bash
bash reproduce_science.sh
```

For provenance capture use

```bash
bash run_full_science_audit.sh
```

On Picasso use `run_full_science_picasso.slurm`.

## Strict content policy

The scientific repository contains no article `.tex`, `.bib`, article PDF, or manuscript-build artifacts. Scientific figure outputs are allowed; manuscript files remain in the journal/preprint channel.

## Archived release

- GitHub release: https://github.com/FrancisRVillatoro/asymptotic-identifiability-inst-mfa/releases/tag/v1.0.1
- Zenodo concept DOI (all versions): https://doi.org/10.5281/zenodo.22737808

The Zenodo version DOI for `v1.0.1` is added to `main` immediately after the GitHub-triggered Zenodo archive is published.
