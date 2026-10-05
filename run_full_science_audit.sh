#!/usr/bin/env bash
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT" || false
STAMP="${FULL_RUN_STAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="audit/full_runs_v101/${STAMP}"
mkdir -p "$RUN_DIR"
rm -f audit/full_science_v101_success.json

stable_manifest() {
  out="$1"
  find . -type f ! -path './audit/full_runs_v101/*' ! -path '*/__pycache__/*' ! -name '*.pyc' -print0 | sort -z | xargs -0 sha256sum > "$out"
}
outputs_manifest() {
  out="$1"
  find data validation figures figures_v101 -type f \( -name '*.csv' -o -name '*.json' -o -name '*.npz' -o -name '*.npy' -o -name '*.pdf' -o -name '*.png' -o -name '*.txt' \) ! -path '*/__pycache__/*' -print0 2>/dev/null | sort -z | xargs -0 sha256sum > "$out"
}
{
 echo "phase=v1.0.1_candidate_scientific_compute"
 echo "stamp_utc=$STAMP"
 echo "start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
 echo "host=$(hostname)"
 echo "pwd=$ROOT"
 echo "uname=$(uname -a)"
 echo "python=$(command -v python)"
 python --version 2>&1
 python - <<'PY'
mods=['numpy','scipy','pandas','matplotlib','openpyxl','sympy']
for m in mods:
    mod=__import__(m); print(f'{m}={getattr(mod,"__version__","unknown")}')
PY
 echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
 echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-unset}"
 echo "MKL_NUM_THREADS=${MKL_NUM_THREADS:-unset}"
} > "$RUN_DIR/environment_compute.txt"
stable_manifest "$RUN_DIR/manifest_before_sha256.txt"
outputs_manifest "$RUN_DIR/outputs_before_sha256.txt"
if [ -x /usr/bin/time ]; then
  /usr/bin/time -v -o "$RUN_DIR/time_science_verbose.txt" bash reproduce_science.sh > >(tee "$RUN_DIR/science_stdout.log") 2> >(tee "$RUN_DIR/science_stderr.log" >&2)
  RC=$?
else
  echo "/usr/bin/time not available" > "$RUN_DIR/time_science_verbose.txt"
  bash reproduce_science.sh > >(tee "$RUN_DIR/science_stdout.log") 2> >(tee "$RUN_DIR/science_stderr.log" >&2)
  RC=$?
fi
stable_manifest "$RUN_DIR/manifest_after_sha256.txt"
outputs_manifest "$RUN_DIR/outputs_after_sha256.txt"
MARKER=false
if grep -Fq 'FULL_SCIENTIFIC_REPRODUCTION_V101_COMPLETED' "$RUN_DIR/science_stdout.log"; then MARKER=true; fi
CLAIM_PASS=false
if grep -Fq 'V15_CLAIM_AUDIT=PASS' "$RUN_DIR/science_stdout.log"; then CLAIM_PASS=true; fi
END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
python - "$STAMP" "$RUN_DIR" "$RC" "$MARKER" "$CLAIM_PASS" "$END_UTC" <<'PY'
import json,sys
stamp,run_dir,rc,marker,claim,end=sys.argv[1:]
ev={'phase':'v1.0.1_candidate_scientific_compute','stamp_utc':stamp,'run_dir':run_dir,'exit_code':int(rc),'completion_marker_found':marker.lower()=='true','v15_claim_audit_pass':claim.lower()=='true','end_utc':end}
open(f'{run_dir}/science_run_result.json','w').write(json.dumps(ev,indent=2)+'\n')
if ev['exit_code']==0 and ev['completion_marker_found'] and ev['v15_claim_audit_pass']:
    open('audit/full_science_v101_success.json','w').write(json.dumps(ev,indent=2)+'\n')
PY
{
 echo "END_UTC=$END_UTC"
 echo "EXIT_CODE=$RC"
 echo "COMPLETION_MARKER_FOUND=$MARKER"
 echo "V15_CLAIM_AUDIT_PASS=$CLAIM_PASS"
 echo "RUN_DIR=$RUN_DIR"
 if [ "$RC" -eq 0 ] && [ "$MARKER" = true ] && [ "$CLAIM_PASS" = true ]; then echo "V101_CANDIDATE=VERIFIED"; else echo "V101_CANDIDATE=FAILED"; fi
} | tee "$RUN_DIR/SCIENCE_RUN_SUMMARY.txt"
if [ "$RC" -eq 0 ] && [ "$MARKER" = true ] && [ "$CLAIM_PASS" = true ]; then true; else false; fi
