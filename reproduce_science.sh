#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

stage() {
  local label="$1"; shift
  local t0 t1
  t0=$(date +%s)
  echo "=== START ${label} | $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  "$@"
  t1=$(date +%s)
  echo "=== END   ${label} | elapsed=$((t1-t0)) s ==="
}

stage "Block 1A FreeFlux runtime vs exact" bash -c 'cd validation/freeflux_runtime && python compare_runtime_to_exact.py'
stage "Block 1B FreeFlux sensitivity refinement" bash -c 'cd validation/freeflux_runtime && python refine_sensitivity.py && python make_block1_deliverables.py'

stage "Block 2A Synechocystis nominal" bash -c 'cd validation/synechocystis_block2 && python run_syn_nominal.py'
stage "Block 2B Synechocystis joint sensitivity" bash -c 'cd validation/synechocystis_block2 && python joint_sensitivity_cached.py'
stage "Block 2C Synechocystis fast2 sweep" bash -c 'cd validation/synechocystis_block2 && python asymptotic_sweep_fast2_extended.py'
stage "Block 2D Synechocystis fast3 sweep" bash -c 'cd validation/synechocystis_block2 && python asymptotic_sweep_fast3_extended.py'
stage "Block 2E Synechocystis design lifting" bash -c 'cd validation/synechocystis_block2 && python fast2_design_lifting.py'

stage "Block 3A partition sweeps" bash -c 'cd validation/partition_block3 && python run_partition_sweeps.py'
stage "Block 3B partition diagnostics" bash -c 'cd validation/partition_block3 && python finalize_block3.py'

stage "Block 4A toy finite-noise" bash -c 'cd validation/block4_finite_noise && python run_block4_toy.py'
stage "Block 4B toy profiles" bash -c 'cd validation/block4_finite_noise && python run_toy_profiles_only.py'
stage "Block 4C Synechocystis finite-noise" bash -c 'cd validation/block4_finite_noise && python run_block4_syn.py'
stage "Block 4D Synechocystis pool-contrast profiles" bash -c 'cd validation/block4_finite_noise && python run_block4_syn_poolcontrast_profiles.py'

stage "Figure 7 refresh (Python only)" python figure_sources/fig7_finite_noise_validation.py

echo "FULL_SCIENTIFIC_REPRODUCTION_COMPLETED"
