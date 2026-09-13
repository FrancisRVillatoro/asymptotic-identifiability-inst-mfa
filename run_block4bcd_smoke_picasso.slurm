#!/bin/bash -l
#SBATCH --job-name=ident_b4bcd_smoke
#SBATCH --partition=cpu_parti
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=audit/block4bcd-smoke-%j.out
#SBATCH --error=audit/block4bcd-smoke-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module purge
module load anaconda/2022.10
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python --version
python - <<'PY'
import numpy, scipy, pandas
print('numpy ', numpy.__version__)
print('scipy ', scipy.__version__)
print('pandas', pandas.__version__)
PY

echo '=== START BLOCK4B ==='
(cd validation/block4_finite_noise && python run_toy_profiles_only.py)
echo 'BLOCK4B=PASS'

echo '=== START BLOCK4C ==='
(cd validation/block4_finite_noise && python run_block4_syn.py)
echo 'BLOCK4C=PASS'

echo '=== START BLOCK4D ==='
(cd validation/block4_finite_noise && python run_block4_syn_poolcontrast_profiles.py)
echo 'BLOCK4D=PASS'

echo 'BLOCK4BCD_SMOKE=PASS'
