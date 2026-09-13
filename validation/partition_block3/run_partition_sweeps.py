#!/usr/bin/env python3
"""Regenerate the seven Block-3 epsilon sweeps from the Block-2 Synechocystis engine.

This script is intentionally separate from the manuscript build because the full sweep is
computationally heavier. It uses the same variational sensitivity engine as Block 2 and writes
partition_sweep_*.csv, weak_subspaces_epsmin_*.npz and block3_full_summary_*.json.
"""
from pathlib import Path
import importlib.util, json, time
import numpy as np
import pandas as pd
import scipy.linalg as la

OUT=Path(__file__).resolve().parent
B2=OUT.parent/'synechocystis_block2'
modpath=B2/'asymptotic_sweep_fast2_extended.py'
spec=importlib.util.spec_from_file_location('b2sweep',modpath)
mod=importlib.util.module_from_spec(spec)
# Execute only the common model/compute_J definitions, not the original fast2 sweep.
code=modpath.read_text()
prefix=code.split("partitions={'fast2'")[0]
ns={'__file__':str(modpath),'__name__':'b2_common'}
exec(compile(prefix,str(modpath),'exec'),ns)
compute_J=ns['compute_J']; concs0=ns['concs0']; pool_ids=ns['pool_ids']; model=ns['model']

PARTITIONS={
 'auto_fast1_raw':['GAP'],
 'auto_primary_phys':['F6P','GAP'],
 'auto_primary_raw':['F6P','TA','GAP'],
 'primary_plus1_phys':['G3P','F6P','GAP'],
 'primary_plus1_raw':['G3P','F6P','TA','GAP'],
 'auto_secondary_phys':['G2P','S7P','E4P','G3P','F6P','GAP'],
 'auto_secondary_raw':['G2P','S7P','E4P','G3P','F6P','TA','GAP'],
}
EPS=np.array([1.,.5,.25,.125,.0625,.03125,.015625,.0078125,.00390625,.001953125])

def classify(slopes):
    m0=[];m1=[];m2=[];other=[]
    # Last singular value is the exact global-scale structural null.
    for j,p in enumerate(slopes[:-1],1):
        if abs(p)<.15: m0.append(j)
        elif .75<p<1.25: m1.append(j)
        elif 1.75<p<2.25: m2.append(j)
        else: other.append(j)
    return m0,m1,m2,other

for name,fast in PARTITIONS.items():
    tstart=time.time(); rows=[]; specs=[]; ys=[]; Vlast=None
    for eps in EPS:
        cc=dict(concs0)
        for m in fast: cc[m]=concs0[m]*eps
        y,J,labels=compute_J(cc)
        ny=int(np.size(y)-np.count_nonzero(np.isfinite(y)))
        nJ=int(np.size(J)-np.count_nonzero(np.isfinite(J)))
        print(f"{name} eps={eps:.12g} finite_y={ny==0} finite_J={nJ==0}", flush=True)
        if ny or nJ:
            bad=np.argwhere(~np.isfinite(J))
            sample=bad[:8].tolist()
            raise FloatingPointError(
                f"non-finite Block-3 Jacobian: partition={name}, epsilon={eps:.12g}, "
                f"nonfinite_y={ny}, nonfinite_J={nJ}, first_bad_indices={sample}"
            )
        U,s,Vh=la.svd(J,full_matrices=False)
        specs.append(s); ys.append(y); Vlast=Vh
        for j,z in enumerate(s,1): rows.append((name,eps,j,float(z),float(z*z)))
    specs=np.asarray(specs)
    slopes=np.array([np.polyfit(np.log(EPS[-4:]),np.log(specs[-4:,j]),1)[0] if np.all(specs[-4:,j]>1e-13) else np.nan for j in range(specs.shape[1])])
    m0,m1,m2,other=classify(slopes)
    # First-order convergence of observables toward the epsilon->0 reduced limit is estimated
    # from successive differences, matching the Block-3 audit's purpose.
    dy=np.array([la.norm(ys[i]-ys[i+1]) for i in range(len(ys)-1)])
    rel=np.array([dy[i]/max(la.norm(ys[i+1]),1e-300) for i in range(len(dy))])
    xo=np.log(EPS[1:][-4:]); yo=np.log(np.maximum(dy[-4:],1e-300))
    out_order=float(np.polyfit(xo,yo,1)[0])
    summary={name:{
      'fast_pools':fast,'n_fast':len(fast),'sigma_exponents':slopes.tolist(),
      'm0_indices_desc':m0,'m1_indices_desc':m1,'m1_exponents':[float(slopes[i-1]) for i in m1],
      'm2_indices_desc':m2,'m2_exponents':[float(slopes[i-1]) for i in m2],
      'other_indices_desc':other,'other_exponents':[float(slopes[i-1]) for i in other],
      'structural_count':1,'output_convergence_order':out_order,
      'output_last_relative_delta':float(rel[-1]),'elapsed_seconds':float(time.time()-tstart)
    }}
    pd.DataFrame(rows,columns=['partition','epsilon','index_desc','sigma','fim_eigenvalue']).to_csv(OUT/f'partition_sweep_{name}.csv',index=False)
    np.savez_compressed(OUT/f'weak_subspaces_epsmin_{name}.npz',**{name:Vlast})
    (OUT/f'block3_full_summary_{name}.json').write_text(json.dumps({'partitions':summary},indent=2,allow_nan=True))
    print(name, 'm1=',len(m1),'m2=',len(m2),'other=',len(other),'seconds=',time.time()-tstart,flush=True)

print('Run finalize_block3.py next to regenerate combined diagnostics and summaries.')
