#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
required=[
 'validation/synechocystis_block2/asymptotic_sweep_fast2_extended.py',
 'validation/synechocystis_block2/source_faithful_freeflux_runtime.py',
 'validation/freeflux_runtime/compare_runtime_to_exact.py',
 'data/freeflux_spectrum_extended.csv',
 'data/freeflux_design_lifting.csv',
]
missing=[x for x in required if not (ROOT/x).exists()]
print('ROOT',ROOT)
if missing:
    print('V101_PREFLIGHT=FAIL missing base-v1.0.0 files')
    for x in missing: print('MISSING',x)
    raise SystemExit(1)
try:
    sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    print('GIT_HEAD',sha)
    if sha!='9c886060635688f97c7129402d787fac65602ae9':
        print('NOTE: HEAD is not the original v1.0.0 commit; this is acceptable only after the RC overlay itself has been committed locally.')
except Exception as e:
    print('GIT_HEAD unavailable',e)
print('V101_PREFLIGHT=PASS')
