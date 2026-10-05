#!/usr/bin/env python3
from pathlib import Path
import argparse, json, time
import numpy as np
import scipy.linalg as la
import syn_discrete_frechet_jacobian as sj

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
OUT.mkdir(parents=True,exist_ok=True)
EPS=np.array([1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125])

def qscale_basis():
    alpha=sj.ns['Brel'].T@np.ones(len(sj.total_ids))
    qscale=np.r_[alpha,np.ones(sj.nc)]
    qscale/=la.norm(qscale)
    Q=la.null_space(qscale.reshape(1,-1))
    return qscale,Q

def concentrations(eps, d=0.10):
    cc=dict(sj.concs0)
    cc['PG']=10.0*(1.0+d)
    cc['Gc']=10.0*(1.0-d)
    cc['F6P']*=eps
    cc['GAP']*=eps
    return cc

def fname(eps,sub,prefix='generic_exactJ'):
    return OUT/f'{prefix}_eps{eps:.12g}_sub{sub}.npz'

def run_one(eps,sub,d=0.10,prefix='generic_exactJ'):
    qscale,Q=qscale_basis()
    cc=concentrations(eps,d)
    t=time.time(); y,J,setup=sj.compute(cc,sub=sub)
    A=J@Q
    U,s,Vh=la.svd(A,full_matrices=False)
    Vphys=Q@Vh.T
    p=fname(eps,sub,prefix)
    np.savez_compressed(p,y=y,J=J,s=s,Vphys=Vphys,qscale=qscale,epsilon=eps,C_PG=cc['PG'],C_Gc=cc['Gc'],subdivisions=sub,setup_seconds=setup,elapsed_seconds=time.time()-t)
    print(f'{prefix} eps={eps:.12g} sub={sub} d={d:.3f} rank_full={np.linalg.matrix_rank(J)} sigma_tail={s[-3:]} elapsed={time.time()-t:.2f}s',flush=True)
    return y,J,s

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=['smoke','full'],default='full')
    args=ap.parse_args()
    if args.mode=='smoke':
        run_one(1/64,16)
        print('V101_EXACT_J_SMOKE_COMPLETED')
        return
    # Canonical sub64 sweep.
    for eps in EPS:
        run_one(float(eps),64)
    # Time-resolution audit at widely separated epsilon values.
    for eps in (1.0,1/64,1/512):
        for sub in (16,32):
            run_one(float(eps),sub)
    # Original distributed equality point for the size of the genericization perturbation.
    cc=dict(sj.concs0)
    t=time.time(); y0,J0,_=sj.compute(cc,sub=64)
    y1=np.load(fname(1.0,64))['y']
    change={
        'C_PG_original':float(sj.concs0['PG']),'C_Gc_original':float(sj.concs0['Gc']),
        'C_PG_generic':11.0,'C_Gc_generic':9.0,
        'relative_2norm':float(la.norm(y1-y0)/la.norm(y0)),
        'max_abs_component':float(np.max(np.abs(y1-y0))),
        'rms_component':float(np.sqrt(np.mean((y1-y0)**2))),
        'elapsed_original_seconds':float(time.time()-t),
    }
    (OUT/'genericization_output_change.json').write_text(json.dumps(change,indent=2)+'\n')
    print('genericization change',json.dumps(change),flush=True)
    # Robustness scan around d=0.10. Sub16 is sufficient for the local genericity audit;
    # canonical spectra above remain sub64.
    for d in (0.08,0.10,0.12):
        for eps in (1/64,1/512):
            run_one(float(eps),16,d=d,prefix=f'genericity_d{d:.2f}')
    print('V101_EXACT_J_FULL_COMPLETED')

if __name__=='__main__':
    main()
