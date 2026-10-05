#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.linalg import helmert
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
INP = ROOT / 'results'
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)

EPS = np.array([1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125])
DELTA = 1e-4
SIGMA_ILR = 0.08
SIGMA_POOL = 0.10
SIGMA_SCALE = 0.10
EPS0 = 1.0/64.0
SEED = 20260905
NMC = 5000

m = np.load(INP/'step3_metadata.npz',allow_pickle=True)
nf = int(m['nf']); nc = int(m['nc']); P = int(m['P'])
pool_ids = [str(x) for x in m['pool_ids']]
pidx = {x:i for i,x in enumerate(pool_ids)}
qscale = m['qscale']; Q = m['Qscale']; rF = m['rF']; rG = m['rG']; rco = m['rco']
idxF = nf+pidx['F6P']; idxG = nf+pidx['GAP']; idxPG = nf+pidx['PG']; idxGc = nf+pidx['Gc']
qcontrast=np.zeros(P);qcontrast[idxF]=1/np.sqrt(2);qcontrast[idxG]=-1/np.sqrt(2)
qsum=np.zeros(P);qsum[idxF]=1/np.sqrt(2);qsum[idxG]=1/np.sqrt(2)
qpg=np.zeros(P);qpg[idxPG]=1/np.sqrt(2);qpg[idxGc]=-1/np.sqrt(2)
SCALE_NORM_RAW=1.0/float(rco@qscale)  # raw global log-scale vector has rco dot qscale_raw = 1

TARGETS=[('G3P_23',2),('G3P_123',3),('DHAP_123',3),('PEP_123',3),('Fum_1234',4),('R5P_12345',5),('RuBP_12345',5),('Cit_12345',5),('Cit_123456',6),('S7P_1234567',7)]
TIMES=[10,30,60,120,240]
GROUPS=[]; row=0
for name,nred in TARGETS:
    for t in TIMES:
        idx=np.arange(row,row+nred)
        GROUPS.append((f'{name}|t={t}',idx))
        row += nred
assert row == 215

def eps_file(eps):
    cands=[INP/f'generic_exactJ_eps{eps:.12g}_sub64.npz',INP/f'generic_exactJ_eps{eps}_sub64.npz']
    for p in cands:
        if p.exists(): return p
    raise FileNotFoundError(eps)

def load_eps(eps):
    d=np.load(eps_file(eps))
    assert abs(float(d['C_PG'])-11.0)<1e-12
    assert abs(float(d['C_Gc'])-9.0)<1e-12
    assert int(d['subdivisions'])==64
    return d['y'], d['J']

def transform_ilr(y,J):
    zs=[];Js=[];ps=[];meta=[]
    for pref,idx in GROUPS:
        pr=y[idx]; p=np.r_[pr,1.0-pr.sum()]
        if p.min() < -1e-9:
            raise RuntimeError(f'Negative MID component {pref}: {p.min()}')
        p=np.maximum(p,0.0)
        H=helmert(len(p),full=False)
        a=p+DELTA
        z=H@np.log(a)
        Jr=J[idx,:]; Jfull=np.vstack([Jr,-Jr.sum(axis=0)])
        Jz=H@np.diag(1.0/a)@Jfull
        zs.append(z);Js.append(Jz);ps.append(p)
        meta.append((pref,float(p.min()),float(p.max())))
    return np.concatenate(zs),np.vstack(Js),ps,meta

def track_modes(mats,basis):
    ss=[]; Vs=[]
    for A in mats:
        _,s,Vh=la.svd(A,full_matrices=False)
        ss.append(s)
        Vs.append(Vh.T if basis is None else basis@Vh.T)
    tracked=[Vs[-1].copy()]; saligned=[ss[-1].copy()]
    prev=Vs[-1]
    for k in range(len(Vs)-2,-1,-1):
        Vcur=Vs[k]
        r,c=linear_sum_assignment(-np.abs(prev.T@Vcur))
        order=np.empty(len(c),dtype=int); order[r]=c
        Vord=Vcur[:,order]
        sg=np.sign(np.sum(prev*Vord,axis=0)); sg[sg==0]=1
        Vord*=sg
        tracked.append(Vord); saligned.append(ss[k][order]); prev=Vord
    return np.asarray(saligned[::-1]), tracked[::-1]

def write_modes(family,designs,basis,prefix):
    sweep=[]; fits=[]; classification=[]
    for dname,mats in designs.items():
        sal,Vs=track_modes(mats,basis)
        for ie,eps in enumerate(EPS):
            for j,(sig,v) in enumerate(zip(sal[ie],Vs[ie].T),1):
                sweep.append((family,dname,eps,j,float(sig),float(la.norm(v[[idxF,idxG]])),float(abs(v@qpg)),float(abs(v@qscale))))
        for j in range(sal.shape[1]):
            slopes=[]
            for n in (3,4,5):
                slopes.append(float(np.polyfit(np.log(EPS[-n:]),np.log(sal[-n:,j]),1)[0]))
            v=Vs[-1][:,j]
            s4=slopes[1]
            if abs(s4)<0.2: cls='m=0'
            elif 0.8<s4<1.2: cls='m=1'
            elif 1.7<s4<2.3: cls='m=2'
            else: cls='other'
            fits.append((family,dname,j+1,*slopes,float(sal[-1,j]),float(la.norm(v[[idxF,idxG]])),float(abs(v@qpg)),float(abs(v@qscale)),cls))
        fd=[r for r in fits if r[0]==family and r[1]==dname]
        counts={c:sum(x[-1]==c for x in fd) for c in ('m=0','m=1','m=2','other')}
        classification.append((family,dname,counts['m=0'],counts['m=1'],counts['m=2'],counts['other']))
    cols=['family','design','epsilon','mode_id','sigma','proj_F6P_GAP','proj_PGminusGc','proj_scale']
    pd.DataFrame(sweep,columns=cols).to_csv(OUT/f'{prefix}_mode_sweep.csv',index=False)
    cols2=['family','design','mode_id','slope3','slope4','slope5','sigma_deep','proj_F6P_GAP_deep','proj_PGminusGc_deep','proj_scale_deep','class']
    pd.DataFrame(fits,columns=cols2).to_csv(OUT/f'{prefix}_mode_exponents.csv',index=False)
    pd.DataFrame(classification,columns=['family','design','n_m0','n_m1','n_m2','n_other']).to_csv(OUT/f'{prefix}_classification.csv',index=False)

# Cache exact raw and ILR Jacobians.
cache={}
for eps in EPS:
    y,J=load_eps(eps); z,Jz,ps,meta=transform_ilr(y,J); cache[eps]=(y,J,z,Jz,ps,meta)

# Unweighted asymptotic design lifting. Constant measurement weights do not affect the formal orders.
raw_cond={k:[] for k in ('baseline','F6P','GAP','both')}
raw_full={k:[] for k in ('F6P','GAP','both','both+scale')}
for eps in EPS:
    J=cache[eps][1]
    raw_cond['baseline'].append(J@Q)
    raw_cond['F6P'].append(np.vstack([J,rF])@Q)
    raw_cond['GAP'].append(np.vstack([J,rG])@Q)
    raw_cond['both'].append(np.vstack([J,rF,rG])@Q)
    raw_full['F6P'].append(np.vstack([J,rF]))
    raw_full['GAP'].append(np.vstack([J,rG]))
    raw_full['both'].append(np.vstack([J,rF,rG]))
    raw_full['both+scale'].append(np.vstack([J,rF,rG,rco]))
write_modes('conditioned_scale_known',raw_cond,Q,'raw_design')
write_modes('full_60_parameter',raw_full,None,'raw_full_design')

# ILR-weighted finite-range sweep; this should preserve the same formal order counts.
w_cond={k:[] for k in ('baseline','F6P','GAP','both')}
w_full={k:[] for k in ('both','both+scale')}
for eps in EPS:
    Jz=cache[eps][3]
    A0=Jz/SIGMA_ILR
    w_cond['baseline'].append(A0@Q)
    w_cond['F6P'].append(np.vstack([A0,rF/SIGMA_POOL])@Q)
    w_cond['GAP'].append(np.vstack([A0,rG/SIGMA_POOL])@Q)
    w_cond['both'].append(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL])@Q)
    w_full['both'].append(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL]))
    w_full['both+scale'].append(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL,rco/SIGMA_SCALE]))
write_modes('conditioned_scale_known',w_cond,Q,'weighted_design')
write_modes('full_60_parameter',w_full,None,'weighted_full_design')

def svd_inverse(A):
    U,s,Vh=la.svd(A,full_matrices=False)
    tol=max(A.shape)*np.finfo(float).eps*s[0]
    keep=s>tol
    Ur=U[:,keep]; sr=s[keep]; Vr=Vh[keep,:].T
    Aplus=(Vr*(1.0/sr))@Ur.T
    cov=(Vr*(1.0/sr**2))@Vr.T
    return s,Vh,keep,Aplus,cov,tol

# Finite-noise local uncertainty at epsilon=1/64.
y,J,z,Jz,ps,meta=cache[EPS0]
A0=Jz/SIGMA_ILR
DES={
 'cond_baseline':(A0@Q,Q,['mid']),
 'cond_F6P':(np.vstack([A0,rF/SIGMA_POOL])@Q,Q,['mid','F6P']),
 'cond_GAP':(np.vstack([A0,rG/SIGMA_POOL])@Q,Q,['mid','GAP']),
 'cond_both':(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL])@Q,Q,['mid','F6P','GAP']),
 'full_both':(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL]),np.eye(P),['mid','F6P','GAP']),
 'full_both_scale':(np.vstack([A0,rF/SIGMA_POOL,rG/SIGMA_POOL,rco/SIGMA_SCALE]),np.eye(P),['mid','F6P','GAP','scale']),
}
unitF=np.eye(P)[idxF]; unitG=np.eye(P)[idxG]
local=[]; inv={}; weakdiag=[]
for name,(A,basis,parts) in DES.items():
    s,Vh,keep,Aplus,cov,tol=svd_inverse(A); covp=basis@cov@basis.T; inv[name]=(Aplus,basis,parts)
    vweak=basis@Vh[-1]
    weakdiag.append((name,float(s[-1]),float(la.norm(vweak[[idxF,idxG]])),float(abs(vweak@qpg)),float(abs(vweak@qscale))))
    def se(q): return float(np.sqrt(max(0.0,q@covp@q)))
    local.append((name,A.shape[1],int(keep.sum()),float(s[-1]),float(s[-2]),float(1/s[-1]),float(1/s[-2]),se(unitF),se(unitG),se(qcontrast),se(qsum),se(qscale),se(qscale)/SCALE_NORM_RAW))
localdf=pd.DataFrame(local,columns=['design','nparam','rank','sigma_min','sigma_second','sd_mode_max','sd_mode_second','se_logC_F6P','se_logC_GAP','se_fast_contrast','se_fast_sum','se_scale_normalized','se_global_log_scale_beta'])
localdf.to_csv(OUT/'finite_noise_local_eps1_64.csv',index=False)
pd.DataFrame(weakdiag,columns=['design','sigma_min','weakest_proj_F6P_GAP','weakest_proj_PGminusGc','weakest_proj_scale']).to_csv(OUT/'finite_noise_weakest_mode_eps1_64.csv',index=False)

# Matched linearized Gaussian Monte Carlo for implementation checking.
rng=np.random.default_rng(SEED)
e_mid=rng.standard_normal((Jz.shape[0],NMC)); eF=rng.standard_normal((1,NMC)); eG=rng.standard_normal((1,NMC)); eS=rng.standard_normal((1,NMC))
mc=[]
for name,(Aplus,basis,parts) in inv.items():
    blocks=[e_mid]
    if 'F6P' in parts: blocks.append(eF)
    if 'GAP' in parts: blocks.append(eG)
    if 'scale' in parts: blocks.append(eS)
    dx=Aplus@np.vstack(blocks); dth=basis@dx
    for label,q in [('F6P',unitF),('GAP',unitG),('fast_contrast',qcontrast),('fast_sum',qsum),('scale_normalized',qscale),('global_log_scale_beta',qscale/SCALE_NORM_RAW)]:
        a=q@dth
        mc.append((name,label,float(np.sqrt(np.mean(a*a))),float(np.quantile(a,.025)),float(np.quantile(a,.975))))
pd.DataFrame(mc,columns=['design','quantity','rmse','q025','q975']).to_csv(OUT/'finite_noise_linear_mc_eps1_64.csv',index=False)

# Correct raw-MID delta-method calibration for the fixed-pseudocount ILR model.
cal=[]; old=[]
for (pref,_),p in zip(GROUPS,cache[1.0][4]):
    p=np.maximum(p,0.0); n=len(p); H=helmert(n,full=False); a=p+DELTA; ss=float(a.sum())
    G=(np.diag(a)-np.outer(a,a)/ss)@H.T
    sd=SIGMA_ILR*np.sqrt(np.diag(G@G.T))
    G0=(np.diag(p)-np.outer(p,p))@H.T
    sd0=SIGMA_ILR*np.sqrt(np.diag(G0@G0.T))
    for k in range(n): cal.append((pref,k,float(p[k]),float(sd[k]),float(sd0[k])))
caldf=pd.DataFrame(cal,columns=['group','component','p','regularized_delta_method_sd','delta0_approx_sd'])
caldf.to_csv(OUT/'ilr_raw_sd_calibration_regularized.csv',index=False)

# Compact summary.
def class_counts(path):
    return pd.read_csv(path).to_dict(orient='records')
summary={
 'benchmark':{'C_PG':11.0,'C_Gc':9.0,'subdivisions':64,'epsilon_grid':EPS.tolist()},
 'noise_model':{'delta':DELTA,'sigma_ilr':SIGMA_ILR,'sigma_log_pool':SIGMA_POOL,'sigma_log_scale_anchor':SIGMA_SCALE,'mc_seed':SEED,'mc_replicates':NMC},
 'raw_conditioned_classification':class_counts(OUT/'raw_design_classification.csv'),
 'raw_full_classification':class_counts(OUT/'raw_full_design_classification.csv'),
 'weighted_conditioned_classification':class_counts(OUT/'weighted_design_classification.csv'),
 'weighted_full_classification':class_counts(OUT/'weighted_full_design_classification.csv'),
 'finite_noise_local':localdf.to_dict(orient='records'),
 'calibration':{
    'regularized_min':float(caldf.regularized_delta_method_sd.min()),
    'regularized_median':float(caldf.regularized_delta_method_sd.median()),
    'regularized_max':float(caldf.regularized_delta_method_sd.max()),
    'delta0_min':float(caldf.delta0_approx_sd.min()),
    'delta0_median':float(caldf.delta0_approx_sd.median()),
    'delta0_max':float(caldf.delta0_approx_sd.max()),
 }
}
(OUT/'step3_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary['raw_conditioned_classification'],indent=2))
print(json.dumps(summary['raw_full_classification'],indent=2))
print(localdf.to_string(index=False))
