#!/usr/bin/env bash
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT" || false
FAIL=0
stage() {
  label="$1"; shift
  t0=$(date +%s)
  echo "=== START ${label} | $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  "$@"
  rc=$?
  t1=$(date +%s)
  echo "=== END   ${label} | rc=${rc} elapsed=$((t1-t0)) s ==="
  if [ "$rc" -ne 0 ]; then FAIL=1; fi
  return 0
}

stage "v1.0.1 candidate preflight" python scripts/v101_preflight.py
if [ "$FAIL" -eq 0 ]; then stage "Inherited v1.0.0 canonical-data hash check" python scripts/verify_inherited_v100.py; fi
if [ "$FAIL" -eq 0 ]; then stage "Synthetic validation regeneration" bash -c 'cd validation/synthetic_validation_v101 && python validate_asymptotic_identifiability.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux v1.0.0 source-faithful runtime cross-check" bash -c 'cd validation/freeflux_runtime && python compare_runtime_to_exact.py && python refine_sensitivity.py && python make_block1_deliverables.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux finite-noise and nonlinear profiles" bash -c 'cd validation/freeflux_v101 && python run_block4_toy.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux refined profile verification" bash -c 'cd validation/freeflux_v101 && python verify_profile_refinement.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux SVD/FDM toy audit" bash -c 'cd validation/freeflux_v101 && python m5_toy_svd_fdm_audit.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux exact positional-isotopomer FDM audit" bash -c 'cd validation/freeflux_v101 && python m5_freeflux_exact_fdm_audit.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "FreeFlux OAA-Fum symbolic check" bash -c 'cd validation/freeflux_v101 && python freeflux_OAA_Fum_symbolic.py > freeflux_OAA_Fum_symbolic_result.txt'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis exact discrete Jacobian + genericization" bash -c 'cd validation/synechocystis_v101 && python run_exact_generic.py --mode full'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis exact-J analysis" bash -c 'cd validation/synechocystis_v101 && python analyze_exact_generic.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis Step-3 metadata" bash -c 'cd validation/synechocystis_v101 && python prepare_step3_metadata.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis design lifting + finite noise" bash -c 'cd validation/synechocystis_v101 && python run_step3.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis local-null audit" bash -c 'cd validation/synechocystis_v101 && python run_local_null_audit.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "Synechocystis physical epsilon=1 finite-noise control" bash -c 'cd validation/synechocystis_v101 && python run_physical_eps1_control.py'; fi
if [ "$FAIL" -eq 0 ]; then stage "Key scientific figure refresh" python figure_sources/plot_v101_key_figures.py; fi
if [ "$FAIL" -eq 0 ]; then stage "Frozen-v15 numerical claim audit" python scripts/check_v15_claims.py; fi
if [ "$FAIL" -eq 0 ] && [ -f scripts/guard_no_paper_artifacts.py ]; then stage "No-paper-artifacts guard" python scripts/guard_no_paper_artifacts.py; fi

if [ "$FAIL" -eq 0 ]; then
  echo "FULL_SCIENTIFIC_REPRODUCTION_V101_COMPLETED"
else
  echo "FULL_SCIENTIFIC_REPRODUCTION_V101_FAILED"
fi
test "$FAIL" -eq 0
