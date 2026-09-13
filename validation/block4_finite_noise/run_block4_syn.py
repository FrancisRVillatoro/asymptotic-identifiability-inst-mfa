#!/usr/bin/env python3
from pathlib import Path
import json, time, warnings
import numpy as np, pandas as pd
import scipy.linalg as la
from scipy.linalg import helmert, null_space, pinv, svd
warnings.filterwarnings('ignore')

OUT=Path(__file__).resolve().parent; OUT.mkdir(exist_ok=True)
SRC=Path(__file__).resolve().parent.parent/'synechocystis_block2'/'asymptotic_sweep_fast2_extended.py'
code=SRC.read_text(); prefix=code.split("partitions={'fast2'")[0]
ns={'__file__':str(SRC),'__name__':'block4_syn_common'}; exec(compile(prefix,str(SRC),'exec'),ns)
compute_J=ns['compute_J']; concs0=ns['concs0']; pool_ids=ns['pool_ids']; pidx=ns['pidx']; nf=ns['nf']; nc=ns['nc']; P=ns['P']
Brel=ns['Brel']; L=ns['L']; fss=ns['fss']; total_ids=ns['total_ids']; model=ns['model']; sim=ns['sim']

DELTA=1e-4; SIGMA_ILR=0.08; SIGMA_POOL=0.10; FAST=['F6P','GAP']
alpha=Brel.T@np.ones(len(total_ids)); qscale=np.r_[alpha,np.ones(nc)]; qscale/=la.norm(qscale)
Q=null_space(qscale.reshape(1,-1)) # 60 x 59

# transform reduced raw MID rows to regularized ILR rows and sensitivities
def transform_ilr(y,J,labels):
    # labels are contiguous by MID/time group
    groups=[]; i=0
    while i<len(labels):
        pref=labels[i].rsplit('|M',1)[0]; idx=[]
        while i<len(labels) and labels[i].rsplit('|M',1)[0]==pref:
            idx.append(i); i+=1
        groups.append((pref,idx))
    zs=[]; Js=[]; outmeta=[]
    for pref,idx in groups:
        pr=y[idx]; plast=1.0-pr.sum(); p=np.r_[pr,plast]
        Jr=J[idx,:]; Jfull=np.vstack([Jr,-Jr.sum(axis=0)])
        H=helmert(len(p),full=False)
        z=H@np.log(np.maximum(p,0)+DELTA)
        Dz=H@np.diag(1.0/(np.maximum(p,0)+DELTA))
        Jz=Dz@Jfull
        zs.append(z); Js.append(Jz)
        for k in range(len(z)): outmeta.append((pref,k,float(p.min()),float(p.max())))
    return np.concatenate(zs),np.vstack(Js),outmeta

def pool_row(pool):
    r=np.zeros(P); r[nf+pidx[pool]]=1.; return r

def design_matrix(Jz,design):
    A=[Jz/SIGMA_ILR]
    if design in ('F6P','both'): A.append(pool_row('F6P')[None,:]/SIGMA_POOL)
    if design in ('GAP','both'): A.append(pool_row('GAP')[None,:]/SIGMA_POOL)
    return np.vstack(A)@Q

# Compute four homotopy points.
epsvals=np.array([1.0,0.25,0.0625,0.015625])
DES=['baseline','F6P','GAP','both']
records=[]; caches={}
for eps in epsvals:
    cc=dict(concs0)
    for m in FAST: cc[m]=concs0[m]*eps
    st=time.time(); y,J,labels=compute_J(cc); z,Jz,meta=transform_ilr(y,J,labels)
    caches[eps]=(cc,y,J,labels,z,Jz,meta)
    print('eps',eps,'compute sec',time.time()-st,flush=True)
    for d in DES:
        A=design_matrix(Jz,d); s=svd(A,compute_uv=False)
        records += [(eps,d,j+1,float(x),float(1/x)) for j,x in enumerate(s)]

df=pd.DataFrame(records,columns=['epsilon','design','index_desc','sigma','local_sd_mode'])
df.to_csv(OUT/'syn_compositional_local_sweep.csv',index=False)
# exponent fits on last three points
fits=[]
for d in DES:
    piv=df[df.design==d].pivot(index='epsilon',columns='index_desc',values='sigma').sort_index()
    ee=piv.index.to_numpy(); use=ee<=0.25; x=np.log(ee[use])
    for j in piv.columns:
        yy=piv.loc[ee[use],j].to_numpy(); p=np.polyfit(x,np.log(yy),1)[0]
        fits.append((d,int(j),float(p),float(-p)))
fitdf=pd.DataFrame(fits,columns=['design','index_desc','sigma_exponent','uncertainty_sd_exponent'])
fitdf.to_csv(OUT/'syn_compositional_exponents.csv',index=False)

# Monte Carlo at epsilon=1/64
eps0=0.015625; cc,y,J,labels,z,Jz,meta=caches[eps0]
rng=np.random.default_rng(20260905); mc=[]; local=[]
for d in DES:
    A=design_matrix(Jz,d); cov=pinv(A.T@A); vals,vecs=la.eigh(cov); order=np.argsort(vals)[::-1]; vals=vals[order]; vecs=vecs[:,order]
    local.append((d,float(np.sqrt(vals[0])),float(np.sqrt(vals[1])),float(svd(A,compute_uv=False)[-1]),float(svd(A,compute_uv=False)[-2])))
    nmc=5000; e=rng.standard_normal((A.shape[0],nmc)); dx=pinv(A)@e; dth=Q@dx
    for pool in FAST:
        arr=dth[nf+pidx[pool],:]; lo,hi=np.quantile(arr,[.025,.975])
        mc.append((d,pool,float(np.sqrt(np.mean(arr**2))),float(lo),float(hi)))
pd.DataFrame(local,columns=['design','largest_sd','second_largest_sd','sigma_min','sigma_second_min']).to_csv(OUT/'syn_compositional_local_uncertainty_eps1_64.csv',index=False)
pd.DataFrame(mc,columns=['design','pool','rmse_logC','q025','q975']).to_csv(OUT/'syn_compositional_monte_carlo_eps1_64.csv',index=False)

# Nonlinear weak-coordinate + linearized nuisance profile at epsilon=1/64.
Abase=design_matrix(Jz,'baseline'); U,s,Vh=svd(Abase,full_matrices=False); qweak=Q@Vh[-1]; qweak/=la.norm(qweak)
# save direction composition
qrows=[]
for i,x in enumerate(qweak[:nf]): qrows.append(('flux_tangent_'+str(i+1),float(x)))
for m,x in zip(pool_ids,qweak[nf:]): qrows.append(('logC_'+m,float(x)))
pd.DataFrame(qrows,columns=['parameter','qweak']).to_csv(OUT/'syn_compositional_weak_direction_eps1_64.csv',index=False)

# true transformed prediction from cache z; matched noise
zobs=z+SIGMA_ILR*rng.standard_normal(len(z))
poolobs={m:np.log(cc[m])+SIGMA_POOL*rng.standard_normal() for m in FAST}

# helper nonlinear prediction along qweak; keep same prepared EMU decomposition
def nonlinear_z(a):
    u=a*qweak[:nf]; f=fss+L@u
    if np.min(f)<=0: return None
    for k,v in zip(total_ids,f): model.total_fluxes[k]=float(v)
    for m in pool_ids: model.concentrations[m]=float(cc[m]*np.exp(a*qweak[nf+pidx[m]]))
    raw=sim.simulate_raw()
    zz=[]
    # exact order matching compute_J labels
    for emuid in model.target_EMUs:
        for t in [10.,30.,60.,120.,240.]:
            p=np.asarray(raw[emuid][t],float); H=helmert(len(p),full=False); zz.append(H@np.log(np.maximum(p,0)+DELTA))
    return np.concatenate(zz)

Uprof=null_space(np.vstack([qscale,qweak]))
agrid=np.array([-1.5,-1.0,-.75,-.5,-.3,-.2,-.1,0,.1,.2,.3,.5,.75,1.0,1.5])
profiles=[]
for a in agrid:
    zp=nonlinear_z(float(a))
    if zp is None: continue
    rmid=(zp-zobs)/SIGMA_ILR
    for d in DES:
        rows=[Jz/SIGMA_ILR]; rr=[rmid]
        if d in ('F6P','both'):
            rows.append(pool_row('F6P')[None,:]/SIGMA_POOL)
            pred=np.log(cc['F6P'])+a*qweak[nf+pidx['F6P']]; rr.append(np.array([(pred-poolobs['F6P'])/SIGMA_POOL]))
        if d in ('GAP','both'):
            rows.append(pool_row('GAP')[None,:]/SIGMA_POOL)
            pred=np.log(cc['GAP'])+a*qweak[nf+pidx['GAP']]; rr.append(np.array([(pred-poolobs['GAP'])/SIGMA_POOL]))
        Afull=np.vstack(rows); r=np.concatenate(rr); AU=Afull@Uprof
        res=r-AU@(pinv(AU)@r); chi=float(res@res)
        profiles.append((d,float(a),chi))
prof=pd.DataFrame(profiles,columns=['design','a_qweak','chi2_gn_profile'])
prof['delta_chi2']=prof.groupby('design')['chi2_gn_profile'].transform(lambda x:x-x.min())
prof.to_csv(OUT/'syn_compositional_noisy_gn_profiles_eps1_64.csv',index=False)

def interval(d,thr=3.84):
    d=d.sort_values('a_qweak'); x=d.a_qweak.to_numpy(); y=d.delta_chi2.to_numpy(); im=np.argmin(y); lo=hi=None
    for i in range(im-1,-1,-1):
        if y[i]>=thr and y[i+1]<thr:
            lo=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    for i in range(im,len(x)-1):
        if y[i]<thr and y[i+1]>=thr:
            hi=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    return lo,hi
ints=[]
for d in DES:
    lo,hi=interval(prof[prof.design==d]); ints.append((d,lo,hi,None if lo is None or hi is None else hi-lo))
pd.DataFrame(ints,columns=['design','lo95','hi95','width95']).to_csv(OUT/'syn_compositional_noisy_gn_profile_intervals_eps1_64.csv',index=False)

# raw-space SD calibration under logistic-normal delta method at published eps=1
cc1,y1,J1,l1,z1,Jz1,meta1=caches[1.0]
# reconstruct full groups and compute raw component sd distribution
cal=[]; i=0
while i<len(l1):
    pref=l1[i].rsplit('|M',1)[0]; idx=[]
    while i<len(l1) and l1[i].rsplit('|M',1)[0]==pref: idx.append(i); i+=1
    pr=y1[idx]; p=np.r_[pr,1-pr.sum()]; H=helmert(len(p),full=False); G=(np.diag(p)-np.outer(p,p))@H.T
    sd=SIGMA_ILR*np.sqrt(np.diag(G@G.T))
    for k,v in enumerate(sd): cal.append((pref,k,float(p[k]),float(v)))
pd.DataFrame(cal,columns=['group','component','p','delta_method_sd']).to_csv(OUT/'syn_ilr_raw_sd_calibration.csv',index=False)

summary={
 'noise_model':{'family':'pseudocount-regularized ILR Gaussian','delta':DELTA,'sigma_ilr':SIGMA_ILR,'sigma_log_pool':SIGMA_POOL,'seed':20260905},
 'epsilon_profile':eps0,
 'local':pd.read_csv(OUT/'syn_compositional_local_uncertainty_eps1_64.csv').to_dict(orient='records'),
 'profile_intervals':pd.read_csv(OUT/'syn_compositional_noisy_gn_profile_intervals_eps1_64.csv').replace({np.nan:None}).to_dict(orient='records'),
 'weak_direction_top':sorted(qrows,key=lambda x:abs(x[1]),reverse=True)[:10]
}
(OUT/'syn_block4_summary.json').write_text(json.dumps(summary,indent=2))
print('\nLOCAL EPS=1/64'); print(pd.read_csv(OUT/'syn_compositional_local_uncertainty_eps1_64.csv').to_string(index=False))
print('\nMC'); print(pd.read_csv(OUT/'syn_compositional_monte_carlo_eps1_64.csv').to_string(index=False))
print('\nPROFILE'); print(pd.read_csv(OUT/'syn_compositional_noisy_gn_profile_intervals_eps1_64.csv').to_string(index=False))
print('\nBOTTOM EXPONENTS')
for d in DES:
    print(d); print(fitdf[fitdf.design==d].tail(5).to_string(index=False))
