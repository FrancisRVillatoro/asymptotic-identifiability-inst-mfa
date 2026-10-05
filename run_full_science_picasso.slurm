#!/bin/bash -l
#SBATCH --job-name=ident_v101_rc1
#SBATCH --partition=cpu_parti
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=audit/v101-science-slurm-%j.out
#SBATCH --error=audit/v101-science-slurm-%j.err

cd "${SLURM_SUBMIT_DIR}" || false
module purge
module load anaconda/2022.10
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
mkdir -p audit

echo "=== PICASSO v1.0.1 RC1 PREFLIGHT ==="
command -v python
python --version
python - <<'PY'
import numpy,scipy,pandas,matplotlib,openpyxl,sympy
for m in (numpy,scipy,pandas,matplotlib,openpyxl,sympy): print(m.__name__,m.__version__)
PY
FAIL=0
python validation/synechocystis_v101/run_exact_generic.py --mode smoke || FAIL=1
if [ "$FAIL" -eq 0 ]; then bash run_full_science_audit.sh || FAIL=1; fi
if [ "$FAIL" -eq 0 ]; then echo "PICASSO_V101_RC1_JOB_PASS"; else echo "PICASSO_V101_RC1_JOB_FAIL"; fi
test "$FAIL" -eq 0
