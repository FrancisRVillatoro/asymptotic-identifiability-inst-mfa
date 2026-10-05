#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment
import syn_discrete_frechet_jacobian as sj

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
EPS=np.array([1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125])
alpha=sj.ns['Brel'].T@np.ones(len(sj.total_ids)); qscale=np.r_[alpha,np.ones(sj.nc)]; qscale/=la.norm(qscale); Q=la.null_space(qscale.reshape(1,-1))
fast_idx=[sj.nf+sj.pidx['F6P'],sj.nf+sj.pidx['GAP']]
pg_idx=sj.nf+sj.pidx['PG']; gc_idx=sj.nf+sj.pidx['Gc']; qpg=np.zeros(sj.P); qpg[pg_idx]=1; qpg[gc_idx]=-1; qpg/=la.norm(qpg)

def fpath(eps,sub=64,prefix='generic_exactJ'):
    return OUT/f'{prefix}_eps{eps:.12g}_sub{sub}.npz'

vecs=[];ss=[];spec=[];runs=[]
for eps in EPS:
    d=np.load(fpath(float(eps)))
    s=d['s'];V=d['Vphys'];ss.append(s);vecs.append(V);runs.append((eps,s[0],s[-1]))
    for k,z in enumerate(s,1): spec.append((eps,k,z))
pd.DataFrame(runs,columns=['epsilon','sigma_max','sigma_min_quotient']).to_csv(OUT/'generic_run_summary.csv',index=False)
pd.DataFrame(spec,columns=['epsilon','index_desc','sigma']).to_csv(OUT/'generic_spectrum.csv',index=False)
# Track modes from deep epsilon to shallow by overlap.
tracked=[vecs[-1].copy()]; prev=vecs[-1]
for Vcur in vecs[-2::-1]:
    r,c=linear_sum_assignment(-np.abs(prev.T@Vcur)); order=np.empty(len(c),int); order[r]=c
    Vord=Vcur[:,order]; sign=np.sign(np.sum(prev*Vord,axis=0)); sign[sign==0]=1; Vord*=sign
    tracked.append(Vord); prev=Vord
tracked=tracked[::-1]
rows=[]
for ie,eps in enumerate(EPS):
    Vorig=vecs[ie]; Vt=tracked[ie]; s=ss[ie]; ov=np.abs(Vorig.T@Vt); match=np.argmax(ov,axis=0)
    for mid,v in enumerate(Vt.T,1):
        rows.append((eps,mid,match[mid-1]+1,s[match[mid-1]],la.norm(v[fast_idx]),abs(v@qpg),abs(v@qscale),ov[match[mid-1],mid-1]))
t=pd.DataFrame(rows,columns=['epsilon','mode_id','index_desc_current','sigma','proj_F6P_GAP','proj_PGminusGc','proj_scale','tracking_overlap'])
t.to_csv(OUT/'generic_mode_tracking.csv',index=False)
fits=[]
for mid,g in t.groupby('mode_id'):
    # smallest epsilon rows are at tail after sorting epsilon descending, use explicit sort ascending.
    g=g.sort_values('epsilon')
    vals=[]
    for n in (3,4,5,6):
        gg=g.head(n); vals.append(float(np.polyfit(np.log(gg.epsilon),np.log(gg.sigma),1)[0]))
    deep=g.iloc[0]
    fits.append((mid,*vals,deep.proj_F6P_GAP,deep.proj_PGminusGc,deep.sigma))
f=pd.DataFrame(fits,columns=['mode_id','slope3','slope4','slope5','slope6','deep_proj_F6P_GAP','deep_proj_PGminusGc','sigma_deep']).sort_values('sigma_deep')
f.to_csv(OUT/'generic_mode_exponents.csv',index=False)
# Numerical rank/null audit of the full J at canonical points.
ranks=[]
for eps in (1.0,1/64,1/512):
    d=np.load(fpath(eps)); J=d['J']; u,s,Vh=la.svd(J,full_matrices=False); tol=max(J.shape)*np.finfo(float).eps*s[0]
    ranks.append((eps,int(np.sum(s>tol)),s[-2],s[-1],abs(float(Vh[-1]@qscale)),float(la.norm(J@qscale)/la.norm(J,2)),float(la.svdvals(J@Q)[-1])))
pd.DataFrame(ranks,columns=['epsilon','rank','sigma_second_smallest','sigma_null','null_scale_alignment','scale_residual_rel','quotient_sigma_min']).to_csv(OUT/'generic_rank_audit.csv',index=False)
# Richardson-style 16/32/64 convergence of the two deepest fast-pool modes.
conv=[]
for eps in (1.0,1/64,1/512):
    vals=[]
    for sub in (16,32,64):
        d=np.load(fpath(eps,sub)); s=d['s']; V=d['Vphys']; proj=np.linalg.norm(V[fast_idx,:],axis=0); ids=np.argsort(proj)[-2:]
        z=sorted([float(s[i]) for i in ids])
        vals.append((sub,z[0],z[1]))
    # p=2 Richardson from 32/64
    w16,w32,w64=vals[0][1],vals[1][1],vals[2][1]; s16,s32,s64=vals[0][2],vals[1][2],vals[2][2]
    conv.append((eps,w16,w32,w64,w64+(w64-w32)/3,abs((w64+(w64-w32)/3-w64)/(w64+(w64-w32)/3)),s16,s32,s64,s64+(s64-s32)/3,abs((s64+(s64-s32)/3-s64)/(s64+(s64-s32)/3))))
pd.DataFrame(conv,columns=['epsilon','weak_sub16','weak_sub32','weak_sub64','weak_richardson_p2','weak_rel64_richardson','strong_sub16','strong_sub32','strong_sub64','strong_richardson_p2','strong_rel64_richardson']).to_csv(OUT/'generic_time_convergence_fast_modes.csv',index=False)
# Genericity d scan: record two modes with largest F6P/GAP projection.
gen=[]
for dval in (0.08,0.10,0.12):
    pref=f'genericity_d{dval:.2f}'
    for eps in (1/64,1/512):
        d=np.load(fpath(eps,16,pref)); s=d['s']; V=d['Vphys']; proj=np.linalg.norm(V[fast_idx,:],axis=0); ids=np.argsort(proj)[-2:]
        ids=ids[np.argsort(s[ids])]
        for j,i in enumerate(ids,1): gen.append((dval,eps,j,float(s[i]),float(proj[i]),float(abs(V[:,i]@qpg))))
pd.DataFrame(gen,columns=['d','epsilon','fast_mode','sigma','proj_F6P_GAP','proj_PGminusGc']).to_csv(OUT/'genericity_robustness.csv',index=False)
summary={
  'rank_audit':pd.read_csv(OUT/'generic_rank_audit.csv').to_dict(orient='records'),
  'fast_mode_exponents':f.sort_values('deep_proj_F6P_GAP',ascending=False).head(2).to_dict(orient='records'),
  'genericization_output_change':json.loads((OUT/'genericization_output_change.json').read_text()),
}
(OUT/'generic_exactJ_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(f.sort_values('deep_proj_F6P_GAP',ascending=False).head(4).to_string(index=False))
print(pd.read_csv(OUT/'generic_rank_audit.csv').to_string(index=False))
print('V101_GENERIC_ANALYSIS_COMPLETED')
