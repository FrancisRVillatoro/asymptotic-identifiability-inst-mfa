#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.linalg import helmert
HERE=Path(__file__).resolve().parent; OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
DELTA=1e-4; SIGMA_ILR=.08; SIGMA_POOL=.10; SIGMA_SCALE=.10
m=np.load(OUT/'step3_metadata.npz',allow_pickle=True)
nf=int(m['nf']); P=int(m['P']); pool_ids=[str(x) for x in m['pool_ids']]; pidx={x:i for i,x in enumerate(pool_ids)}; qscale=m['qscale']; Q=m['Qscale']; rF=m['rF']; rG=m['rG']; rco=m['rco']
idxF=nf+pidx['F6P']; idxG=nf+pidx['GAP']; idxPG=nf+pidx['PG']; idxGc=nf+pidx['Gc']
unitF=np.eye(P)[idxF]; unitG=np.eye(P)[idxG]; qpg=np.zeros(P); qpg[idxPG]=1/np.sqrt(2);qpg[idxGc]=-1/np.sqrt(2)
TARGETS=[('G3P_23',2),('G3P_123',3),('DHAP_123',3),('PEP_123',3),('Fum_1234',4),('R5P_12345',5),('RuBP_12345',5),('Cit_12345',5),('Cit_123456',6),('S7P_1234567',7)]; TIMES=[10,30,60,120,240]
GROUPS=[]; row=0
for name,nred in TARGETS:
    for t in TIMES:
        idx=np.arange(row,row+nred); GROUPS.append((f'{name}|t={t}',idx)); row+=nred
assert row==215
d=np.load(OUT/'generic_exactJ_eps1_sub64.npz'); y,J=d['y'],d['J']
Js=[]
for pref,idx in GROUPS:
    pr=y[idx]; p=np.maximum(np.r_[pr,1-pr.sum()],0.0); H=helmert(len(p),full=False); Jr=J[idx,:]; Jfull=np.vstack([Jr,-Jr.sum(axis=0)]); Js.append(H@np.diag(1/(p+DELTA))@Jfull)
Jz=np.vstack(Js); A0=Jz/SIGMA_ILR

def cov_full(A,basis):
    U,s,Vh=la.svd(A,full_matrices=False); tol=max(A.shape)*np.finfo(float).eps*s[0]; keep=s>tol; Vr=Vh[keep].T; cov=(Vr*(1/s[keep]**2))@Vr.T; return s,basis@cov@basis.T,basis@Vh[-1]
DES=[
 ('cond_baseline',A0@Q,Q),
 ('cond_F6P',np.vstack([A0,rF/SIGMA_POOL])@Q,Q),
 ('cond_GAP',np.vstack([A0,rG/SIGMA_POOL])@Q,Q),
 ('cond_both',np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL])@Q,Q),
 ('full_both',np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL]),np.eye(P)),
 ('full_both_scale',np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL,rco/SIGMA_SCALE]),np.eye(P)),
]
rows=[]
for name,A,basis in DES:
    s,cov,vweak=cov_full(A,basis)
    def se(q): return float(np.sqrt(max(0.0,q@cov@q)))
    rows.append((name,float(s[-1]),se(unitF),se(unitG),se(qscale),float(abs(vweak@qscale)),float(abs(vweak@qpg)),float(la.norm(vweak[[idxF,idxG]]))))
df=pd.DataFrame(rows,columns=['design','sigma_min','se_logC_F6P','se_logC_GAP','se_normalized_scale_direction','weakest_proj_scale','weakest_proj_PGminusGc','weakest_proj_F6P_GAP'])
df.to_csv(OUT/'finite_noise_local_eps1.csv',index=False)
summary={
 'epsilon':1.0,
 'conditioned_baseline_F6P_GAP_SE':[float(df[df.design=='cond_baseline'].se_logC_F6P.iloc[0]),float(df[df.design=='cond_baseline'].se_logC_GAP.iloc[0])],
 'conditioned_both_F6P_GAP_SE':[float(df[df.design=='cond_both'].se_logC_F6P.iloc[0]),float(df[df.design=='cond_both'].se_logC_GAP.iloc[0])],
 'scale_anchor_observation_log_sd':SIGMA_SCALE,
 'note':'The 0.10 value is the assumed log-scale SD of the directly measured anchor; Euclidean coefficients along normalized parameter-space scale vectors are coordinate projections, not the anchor measurement uncertainty.',
}
(OUT/'physical_eps1_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(df.to_string(index=False)); print(json.dumps(summary,indent=2)); print('V101_PHYSICAL_EPS1_CONTROL_COMPLETED')
