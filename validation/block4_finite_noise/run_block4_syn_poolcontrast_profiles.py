#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd, scipy.linalg as la, warnings
from scipy.linalg import helmert, null_space, pinv, svd
warnings.filterwarnings('ignore')
OUT=Path(__file__).resolve().parent
src=OUT/'run_block4_syn.py'
code=src.read_text(); prefix=code.split('# Compute four homotopy points.')[0]
ns={'__file__':str(src),'__name__':'block4_syn_profile_common'}
exec(compile(prefix,str(src),'exec'),ns)
compute_J=ns['compute_J']; concs0=ns['concs0']; pool_ids=ns['pool_ids']; pidx=ns['pidx']; nf=ns['nf']; nc=ns['nc']; P=ns['P']; qscale=ns['qscale']
L=ns['L']; fss=ns['fss']; total_ids=ns['total_ids']; model=ns['model']; sim=ns['sim']; DELTA=ns['DELTA']; SIGMA_ILR=ns['SIGMA_ILR']; SIGMA_POOL=ns['SIGMA_POOL']; transform_ilr=ns['transform_ilr']; pool_row=ns['pool_row']
FAST=['F6P','GAP']; eps0=1/64; cc=dict(concs0)
for m in FAST: cc[m]=concs0[m]*eps0
y,J,labels=compute_J(cc); z,Jz,meta=transform_ilr(y,J,labels)
# pure fast-pool contrast qC = (log F6P - log GAP)/sqrt(2)
qC=np.zeros(P); qC[nf+pidx['F6P']]=1/np.sqrt(2); qC[nf+pidx['GAP']]=-1/np.sqrt(2)
# nuisance: remove exact global scale and hold contrast fixed
Uprof=null_space(np.vstack([qscale,qC]))
rng=np.random.default_rng(20260905)
zobs=z+SIGMA_ILR*rng.standard_normal(len(z)); poolobs={m:np.log(cc[m])+SIGMA_POOL*rng.standard_normal() for m in FAST}
DES=['baseline','F6P','GAP','both']

def nonlinear_z(a):
    # fluxes unchanged; only fast pool contrast perturbed
    for k,v in zip(total_ids,fss): model.total_fluxes[k]=float(v)
    for m in pool_ids: model.concentrations[m]=float(cc[m]*np.exp(a*qC[nf+pidx[m]]))
    raw=sim.simulate_raw(); zz=[]
    for emuid in model.target_EMUs:
        for t in [10.,30.,60.,120.,240.]:
            p=np.asarray(raw[emuid][t],float); H=helmert(len(p),full=False); zz.append(H@np.log(np.maximum(p,0)+DELTA))
    return np.concatenate(zz)

def rows_for(d):
    rows=[Jz/SIGMA_ILR]
    if d in ('F6P','both'): rows.append(pool_row('F6P')[None,:]/SIGMA_POOL)
    if d in ('GAP','both'): rows.append(pool_row('GAP')[None,:]/SIGMA_POOL)
    return np.vstack(rows)

agrid=np.array([-1.5,-1.0,-.75,-.5,-.35,-.25,-.15,-.08,0,.08,.15,.25,.35,.5,.75,1.0,1.5])
profiles=[]
for a in agrid:
    zp=nonlinear_z(float(a)); rmid=(zp-zobs)/SIGMA_ILR
    for d in DES:
        rr=[rmid]
        if d in ('F6P','both'):
            pred=np.log(cc['F6P'])+a*qC[nf+pidx['F6P']]; rr.append(np.array([(pred-poolobs['F6P'])/SIGMA_POOL]))
        if d in ('GAP','both'):
            pred=np.log(cc['GAP'])+a*qC[nf+pidx['GAP']]; rr.append(np.array([(pred-poolobs['GAP'])/SIGMA_POOL]))
        A=rows_for(d); r=np.concatenate(rr); AU=A@Uprof; res=r-AU@(pinv(AU)@r); chi=float(res@res); profiles.append((d,float(a),chi))
prof=pd.DataFrame(profiles,columns=['design','a_contrast','chi2_gn_profile']); prof['delta_chi2']=prof.groupby('design')['chi2_gn_profile'].transform(lambda x:x-x.min()); prof.to_csv(OUT/'syn_poolcontrast_noisy_gn_profiles_eps1_64.csv',index=False)

def interval(d,thr=3.84):
    d=d.sort_values('a_contrast'); x=d.a_contrast.to_numpy(); y=d.delta_chi2.to_numpy(); im=np.argmin(y); lo=hi=None
    for i in range(im-1,-1,-1):
        if y[i]>=thr and y[i+1]<thr:
            lo=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    for i in range(im,len(x)-1):
        if y[i]<thr and y[i+1]>=thr:
            hi=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    return lo,hi
ints=[]
for d in DES:
    lo,hi=interval(prof[prof.design==d]); ints.append((d,lo,hi,None if lo is None or hi is None else hi-lo,float(prof[prof.design==d].delta_chi2.max())))
out=pd.DataFrame(ints,columns=['design','lo95','hi95','width95','max_delta_chi2']); out.to_csv(OUT/'syn_poolcontrast_noisy_gn_profile_intervals_eps1_64.csv',index=False)
print(out.to_string(index=False))
for d in DES:
 print('\n',d); print(prof[prof.design==d][['a_contrast','delta_chi2']].to_string(index=False))
