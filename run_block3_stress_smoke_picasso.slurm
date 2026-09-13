#!/bin/bash -l
#SBATCH --job-name=ident_b3smoke
#SBATCH --partition=cpu_parti
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=audit/block3-smoke-%j.out
#SBATCH --error=audit/block3-smoke-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module purge
module load anaconda/2022.10
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
python scripts/block3_stress_smoke.py
