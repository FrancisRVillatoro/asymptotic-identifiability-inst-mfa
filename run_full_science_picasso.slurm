#!/bin/bash -l
#SBATCH --job-name=ident_science_repro
#SBATCH --partition=cpu_parti
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=audit/science-slurm-%j.out
#SBATCH --error=audit/science-slurm-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

module purge
module load anaconda/2022.10

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p audit

echo "=== PICASSO COMPUTE PREFLIGHT ==="
command -v python
python --version
python - <<'PY'
import numpy, scipy, pandas, matplotlib, openpyxl, sympy
print("numpy      ", numpy.__version__)
print("scipy      ", scipy.__version__)
print("pandas      ", pandas.__version__)
print("matplotlib ", matplotlib.__version__)
print("openpyxl   ", openpyxl.__version__)
print("sympy      ", sympy.__version__)
PY
for c in sha256sum tee awk grep find sort xargs; do
  command -v "$c" >/dev/null
  echo "$c=$(command -v "$c")"
done

echo "=== BLOCK 3 STRESS SMOKE ==="
python scripts/block3_stress_smoke.py
echo "=== START SCIENTIFIC AUDIT ==="
bash run_full_science_audit.sh
