#!/usr/bin/env python3
"""Non-writing stress smoke for the stiffest Block-3 partitions.

Exercises the epsilon=1/512 endpoint for the two broad partitions before the
full sweep.  It is intentionally read-only with respect to scientific outputs.
"""
from pathlib import Path
import numpy as np
import scipy.linalg as la

ROOT=Path(__file__).resolve().parents[1]
runner=ROOT/'validation'/'partition_block3'/'run_partition_sweeps.py'
code=runner.read_text()
prefix=code.split('for name,fast in PARTITIONS.items():')[0]
ns={'__file__':str(runner),'__name__':'block3_stress_smoke'}
exec(compile(prefix,str(runner),'exec'),ns)
compute_J=ns['compute_J']; concs0=ns['concs0']

cases={
    'auto_secondary_phys':['G2P','S7P','E4P','G3P','F6P','GAP'],
    'auto_secondary_raw':['G2P','S7P','E4P','G3P','F6P','TA','GAP'],
}
eps=1.0/512.0
for name,fast in cases.items():
    cc=dict(concs0)
    for m in fast:
        cc[m]=concs0[m]*eps
    y,J,_=compute_J(cc)
    ny=int(np.size(y)-np.count_nonzero(np.isfinite(y)))
    nJ=int(np.size(J)-np.count_nonzero(np.isfinite(J)))
    if ny or nJ:
        bad=np.argwhere(~np.isfinite(J))[:8].tolist()
        raise SystemExit(f'FAIL {name}: nonfinite_y={ny} nonfinite_J={nJ} first_bad={bad}')
    s=la.svdvals(J)
    if not np.isfinite(s).all():
        raise SystemExit(f'FAIL {name}: non-finite singular values')
    stats=compute_J.__globals__.get('_EXPM_FALLBACK_STATS',{})
    print(f'PASS {name} eps={eps:.12g} J_shape={J.shape} sigma_min_nonstruct={s[-2]:.12e} sigma_struct={s[-1]:.12e} expm_stats={stats}',flush=True)
print('BLOCK3_STRESS_SMOKE=PASS')
