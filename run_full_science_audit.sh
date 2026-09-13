#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

STAMP="${FULL_RUN_STAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="audit/full_runs/${STAMP}"
mkdir -p "$RUN_DIR"
rm -f audit/full_science_success.json audit/full_run_success.json

stable_manifest() {
  local out="$1"
  find . -type f \
    ! -path './audit/full_runs/*' \
    ! -path '*/__pycache__/*' \
    ! -name '*.pyc' \
    ! -name '*.aux' ! -name '*.bcf' ! -name '*.blg' ! -name '*.log' \
    ! -name '*.out' ! -name '*.run.xml' ! -name '*.toc' ! -name '*.synctex.gz' \
    ! -name 'full_run_success.json' ! -name 'full_science_success.json' \
    -print0 | sort -z | xargs -0 sha256sum > "$out"
}

outputs_manifest() {
  local out="$1"
  {
    find data validation figures -type f \
      \( -name '*.csv' -o -name '*.json' -o -name '*.npz' -o -name '*.npy' -o -name '*.pdf' -o -name '*.png' -o -name '*.txt' \) \
      ! -path '*/__pycache__/*' -print0
  } | sort -z | xargs -0 sha256sum > "$out"
}

{
  echo "phase=scientific_compute"
  echo "full_run_stamp_utc=$STAMP"
  echo "start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host=$(hostname)"
  echo "pwd=$ROOT"
  echo "uname=$(uname -a)"
  echo "shell=${BASH_VERSION:-unknown}"
  echo "python=$(command -v python)"
  python --version 2>&1
  python - <<'PY'
mods=['numpy','scipy','pandas','matplotlib','openpyxl','sympy']
for m in mods:
    mod=__import__(m)
    print(f'{m}={getattr(mod,"__version__","unknown")}')
PY
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
  echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-unset}"
  echo "MKL_NUM_THREADS=${MKL_NUM_THREADS:-unset}"
} > "$RUN_DIR/environment_compute.txt"

stable_manifest "$RUN_DIR/manifest_before_sha256.txt"
outputs_manifest "$RUN_DIR/outputs_before_sha256.txt"

set +e
if [[ -x /usr/bin/time ]]; then
  /usr/bin/time -v -o "$RUN_DIR/time_science_verbose.txt" \
    bash reproduce_science.sh \
    > >(tee "$RUN_DIR/science_stdout.log") \
    2> >(tee "$RUN_DIR/science_stderr.log" >&2)
  RC=$?
else
  echo "/usr/bin/time not available" > "$RUN_DIR/time_science_verbose.txt"
  bash reproduce_science.sh \
    > >(tee "$RUN_DIR/science_stdout.log") \
    2> >(tee "$RUN_DIR/science_stderr.log" >&2)
  RC=$?
fi
set -e

stable_manifest "$RUN_DIR/manifest_after_science_sha256.txt"
outputs_manifest "$RUN_DIR/outputs_after_science_sha256.txt"

MARKER=false
if grep -Fq 'FULL_SCIENTIFIC_REPRODUCTION_COMPLETED' "$RUN_DIR/science_stdout.log"; then
  MARKER=true
fi

STDOUT_SHA=$(sha256sum "$RUN_DIR/science_stdout.log" | awk '{print $1}')
STDERR_SHA=$(sha256sum "$RUN_DIR/science_stderr.log" | awk '{print $1}')
OUTPUTS_SHA=$(sha256sum "$RUN_DIR/outputs_after_science_sha256.txt" | awk '{print $1}')
END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)

python - "$STAMP" "$RUN_DIR" "$RC" "$MARKER" "$END_UTC" "$STDOUT_SHA" "$STDERR_SHA" "$OUTPUTS_SHA" <<'PY'
import json,sys
stamp,run_dir,rc,marker,end,so,se,oa=sys.argv[1:]
ev={
  'phase':'scientific_compute',
  'stamp_utc':stamp,
  'run_dir':run_dir,
  'exit_code':int(rc),
  'completion_marker_found':marker.lower()=='true',
  'end_utc':end,
  'stdout_sha256':so,
  'stderr_sha256':se,
  'outputs_after_science_sha256_file_sha256':oa,
}
open(f'{run_dir}/science_run_result.json','w').write(json.dumps(ev,indent=2)+'\n')
if ev['exit_code']==0 and ev['completion_marker_found']:
    open('audit/full_science_success.json','w').write(json.dumps(ev,indent=2)+'\n')
PY

{
  echo "END_UTC=$END_UTC"
  echo "EXIT_CODE=$RC"
  echo "SCIENCE_COMPLETION_MARKER_FOUND=$MARKER"
  echo "SCIENCE_STDOUT_SHA256=$STDOUT_SHA"
  echo "SCIENCE_STDERR_SHA256=$STDERR_SHA"
  echo "SCIENCE_OUTPUTS_AFTER_FILE_SHA256=$OUTPUTS_SHA"
  echo "RUN_DIR=$RUN_DIR"
  if [[ "$RC" -eq 0 && "$MARKER" == true ]]; then
    echo "BLOCK6_SCIENCE=VERIFIED_AWAITING_FINAL_BUILD"
  else
    echo "BLOCK6_SCIENCE=FAILED"
  fi
} | tee "$RUN_DIR/SCIENCE_RUN_SUMMARY.txt"

[[ "$RC" -eq 0 ]]
[[ "$MARKER" == true ]]
