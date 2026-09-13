#!/usr/bin/env python3
from pathlib import Path
import csv, json, warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import null_space, helmert, svd, pinv
from scipy.optimize import least_squares
warnings.filterwarnings('ignore')

OUT=Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)

# ---------------- Model ----------------
POOLS=['OAA','Cit','AKG','Suc','Fum','Glu']
NC={'OAA':4,'Cit':6,'AKG':5,'Suc':4,'Fum':4,'Glu':5}
NISO={m:2**NC[m] for m in POOLS}
offset={}; nstate=0
for m in POOLS:
    offset[m]=nstate; nstate+=NISO[m]
CONST=nstate; N=nstate+1
SLOW_TIMES=np.array([0.1,0.2,0.5,1.0,2.0])

def bits(i,n): return tuple((i>>(n-1-j))&1 for j in range(n))
def ibits(b):
    x=0
    for q in b: x=(x<<1)|int(q)
    return x

def trans_matrix(dst_n,src_n,mapper):
    T=np.zeros((2**dst_n,2**src_n))
    for j in range(2**src_n):
        for out,wgt in mapper(bits(j,src_n)): T[ibits(out),j]+=wgt
    return T

def make_v1(accoa_mix):
    def mapper(o):
        a,b,c,d=o
        return [((d,c,b,f,e,a),w) for (e,f),w in accoa_mix]
    return trans_matrix(6,4,mapper)
T_v1_base=make_v1([((0,0),.5),((0,1),.25),((1,1),.25)])
T_v1_unlab=make_v1([((0,0),1.)])
T_v2=trans_matrix(5,6,lambda x:[(x[:5],1.)])
T_v3=np.eye(32)
def map_v4(x):
    s=x[1:5]; return [(s,.5),(s[::-1],.5)]
T_v4=trans_matrix(4,5,map_v4)
def map_sym4(x): return [(x,.5),(x[::-1],.5)]
T_sym=trans_matrix(4,4,map_sym4)

def rs(m): return slice(offset[m],offset[m]+NISO[m])
def fluxes(u,w,r): return {'v1':u+w,'v2':u+w,'v3':u,'v4':w,'v5':w,'v6f':w+r,'v6b':r,'v7':u}

def K_values(u,w,r,C,acc='baseline',asp='unlabeled'):
    F=fluxes(u,w,r); K=np.zeros((N,N))
    T_v1=T_v1_base if acc=='baseline' else T_v1_unlab
    ext_oaa=np.zeros(16); ext_oaa[0 if asp=='unlabeled' else 15]=1.
    def tr(dst,src,v,T): K[rs(dst),rs(src)] += (v/C[dst])*T
    def out(pool,v): K[rs(pool),rs(pool)] -= (v/C[pool])*np.eye(NISO[pool])
    def ext(dst,v,d): K[rs(dst),CONST] += (v/C[dst])*d
    tr('OAA','Fum',F['v6f'],T_sym); ext('OAA',F['v7'],ext_oaa); out('OAA',F['v1']+F['v6b'])
    tr('Cit','OAA',F['v1'],T_v1); out('Cit',F['v2'])
    tr('AKG','Cit',F['v2'],T_v2); out('AKG',F['v3']+F['v4'])
    tr('Suc','AKG',F['v4'],T_v4); out('Suc',F['v5'])
    tr('Fum','Suc',F['v5'],T_sym); tr('Fum','OAA',F['v6b'],T_sym); out('Fum',F['v6f'])
    tr('Glu','AKG',F['v3'],T_v3); out('Glu',F['v3'])
    return K

def initial_state():
    q=np.zeros(N)
    for m in POOLS: q[offset[m]]=1.
    q[CONST]=1.; return q
Q0=initial_state()

def mdv_full(pool,pos):
    n=NC[pool]; pos=[p-1 for p in pos]
    M=np.zeros((len(pos)+1,N))
    for i in range(2**n):
        b=bits(i,n); M[sum(b[j] for j in pos),offset[pool]+i]+=1.
    return M
M_GLU123=mdv_full('Glu',[1,2,3]); M_GLU5=mdv_full('Glu',[1,2,3,4,5]); M_CIT=mdv_full('Cit',[2,3,4,5])
M_FUM=mdv_full('Fum',[1,2,3,4]); M_OAA=mdv_full('OAA',[1,2,3,4])

names=['log v3','log v4','log v6b']+[f'log C_{m}' for m in POOLS]
u0,w0,r0=5.,5.,7.5
C0={'OAA':.1,'Cit':5.,'AKG':.3,'Suc':1.,'Fum':.2,'Glu':.5}
theta0=np.array([np.log(u0),np.log(w0),np.log(r0)]+[np.log(C0[m]) for m in POOLS])

def unpack(th):
    u,w,r=np.exp(th[:3]); C={m:np.exp(th[3+i]) for i,m in enumerate(POOLS)}
    return u,w,r,C

def trajectories(th,times,acc='baseline',asp='unlabeled'):
    u,w,r,C=unpack(th); K=sparse.csr_matrix(K_values(u,w,r,C,acc=acc,asp=asp))
    return [expm_multiply(K*float(t),Q0) for t in times]

def groups_for_design(th,design):
    groups=[]
    slow=trajectories(th,SLOW_TIMES,'baseline','unlabeled')
    for t,q in zip(SLOW_TIMES,slow):
        groups.append((f'Glu123@{t:g}',M_GLU123@q))
        groups.append((f'Glu12345@{t:g}',M_GLU5@q))
        groups.append((f'Cit2345@{t:g}',M_CIT@q))
        if design in ('fum','asp','pool'):
            groups.append((f'Fum1234@{t:g}',M_FUM@q))
    if design in ('asp','pool'):
        # Nominal physical sampling times; protocol clock is fixed, not parameter dependent.
        tau0=C0['OAA']/(u0+w0+r0)
        ftimes=tau0*np.array([0.25,0.5,1.,2.,4.])
        fast=trajectories(th,ftimes,'unlabeled','labeled')
        for t,q in zip(ftimes,fast):
            groups.append((f'AspFast_OAA@{t:.8f}',M_OAA@q))
            groups.append((f'AspFast_Fum@{t:.8f}',M_FUM@q))
    return groups

# ---------------- Statistical model ----------------
DELTA=1e-4
SIGMA_ILR=0.08
SIGMA_POOL=0.10

def ilr(p):
    p=np.asarray(p,float); H=helmert(len(p),full=False)
    return H@np.log(np.maximum(p,0)+DELTA)

def prediction(th,design):
    parts=[]; labels=[]
    for label,p in groups_for_design(th,design):
        z=ilr(p); parts.append(z); labels += [f'{label}|ilr{i+1}' for i in range(len(z))]
    if design=='pool':
        parts.append(np.array([th[3+POOLS.index('AKG')]])); labels.append('logC_AKG_pool')
    return np.concatenate(parts),labels

def jac_fd(th,design,h=2e-5):
    y,_=prediction(th,design); J=np.empty((len(y),len(th)))
    for j in range(len(th)):
        tp=th.copy(); tm=th.copy(); tp[j]+=h; tm[j]-=h
        J[:,j]=(prediction(tp,design)[0]-prediction(tm,design)[0])/(2*h)
    # whiten rows: all ILR have sigma_ilr; final pool row sigma_pool
    scale=np.full(len(y),SIGMA_ILR)
    if design=='pool': scale[-1]=SIGMA_POOL
    return J/scale[:,None], y, scale

def whitened_residual(th,design,yobs):
    yp,_=prediction(th,design)
    scale=np.full(len(yp),SIGMA_ILR)
    if design=='pool': scale[-1]=SIGMA_POOL
    return (yp-yobs)/scale

# Structural directions
# global simultaneous log flux/pool scale
g=np.ones(9); g/=np.linalg.norm(g)
# exact OAA-Fum-exchange symmetry tangent (unnormalized) from the analytic invariant family
s2=np.array([0.,0.,15.,10.,0.,0.,0.,-1.,0.]); s2/=np.linalg.norm(s2)
Qs,_=np.linalg.qr(np.column_stack([g,s2])); Qs=Qs[:,:2]
Nbase=null_space(Qs.T)

# Baseline weak modes under the compositional likelihood
Jb, ybase, scaleb = jac_fd(theta0,'baseline')
Ub,sb,Vhb=svd(Jb@Nbase,full_matrices=False)
q2=Nbase@Vhb[-1]; q2/=np.linalg.norm(q2)
q1=Nbase@Vhb[-2]; q1/=np.linalg.norm(q1)

# Master synthetic noisy observations, shared by nested designs.
rng=np.random.default_rng(20260905)
# Generate the richest design once and subset by labels to guarantee matched perturbations.
yrich,lab_rich=prediction(theta0,'pool')
noise=np.zeros_like(yrich)
noise[:-1]=SIGMA_ILR*rng.standard_normal(len(yrich)-1)
noise[-1]=SIGMA_POOL*rng.standard_normal()
yrich_obs=yrich+noise
master=dict(zip(lab_rich,yrich_obs))

def noisy_obs(design):
    y,labs=prediction(theta0,design)
    return np.array([master[x] for x in labs]),labs

# Local information and Monte Carlo uncertainty
DESIGNS=['baseline','fum','asp','pool']
pretty={'baseline':'baseline','fum':'+Fum slow','asp':'+Fum + Asp-fast','pool':'+Fum + Asp-fast + pool AKG'}
local_rows=[]; mc_rows=[]
for design in DESIGNS:
    Jw,y,scale=jac_fd(theta0,design)
    if design=='baseline': Q=Nbase
    elif design in ('fum','asp'): Q=null_space(g.reshape(1,-1))
    else: Q=np.eye(9)
    A=Jw@Q
    U,s,Vh=svd(A,full_matrices=False)
    cov=pinv(A.T@A)
    # variances of physical q1/q2 coordinates
    b1=Q.T@q1; b2=Q.T@q2
    se1=float(np.sqrt(max(0,b1@cov@b1))); se2=float(np.sqrt(max(0,b2@cov@b2)))
    local_rows.append((design,A.shape[1],s[-1],s[-2] if len(s)>1 else np.nan,se1,se2))
    # local Monte Carlo under transformed Gaussian likelihood
    nmc=4000
    e=rng.standard_normal((A.shape[0],nmc))
    dxi=pinv(A)@e
    e1=b1@dxi; e2=b2@dxi
    for label,arr,se in [('q1',e1,se1),('q2',e2,se2)]:
        qlo,qhi=np.quantile(arr,[0.025,0.975])
        mc_rows.append((design,label,float(np.sqrt(np.mean(arr**2))),float(qlo),float(qhi),float(np.mean(np.abs(arr)<=1.96*se))))

pd.DataFrame(local_rows,columns=['design','identifiable_dim','sigma_min','sigma_second_min','se_q1','se_q2']).to_csv(OUT/'toy_compositional_local_uncertainty.csv',index=False)
pd.DataFrame(mc_rows,columns=['design','direction','rmse','q025','q975','nominal95_coverage']).to_csv(OUT/'toy_compositional_monte_carlo.csv',index=False)

# Full nonlinear MLE and profile along the original baseline m=2 direction.
agrid=np.array([-2.0,-1.5,-1.0,-0.75,-0.5,-0.30,-0.20,-0.10,-0.05,0.,0.05,0.10,0.20,0.30,0.5,0.75,1.0,1.5,2.0])
prof_rows=[]; fit_summary=[]
for design in DESIGNS:
    yobs,_=noisy_obs(design)
    if design=='baseline': Q=Nbase
    elif design in ('fum','asp'): Q=null_space(g.reshape(1,-1))
    else: Q=np.eye(9)
    # global MLE in identifiable coordinates
    def rfull(x): return whitened_residual(theta0+Q@x,design,yobs)
    fit=least_squares(rfull,np.zeros(Q.shape[1]),bounds=(-2.5,2.5),max_nfev=250,xtol=2e-9,ftol=2e-9,gtol=2e-9)
    chi_min=float(fit.fun@fit.fun); thhat=theta0+Q@fit.x; ahat=float(q2@(thhat-theta0))
    # profile nuisance basis
    if design=='baseline': Uprof=null_space(np.vstack([g,s2,q2]))
    elif design in ('fum','asp'): Uprof=null_space(np.vstack([g,q2]))
    else: Uprof=null_space(q2.reshape(1,-1))
    zprev=np.zeros(Uprof.shape[1])
    vals=[]
    # start near 0 and warm-start outward separately
    for sign in [1,-1]:
        z=np.zeros(Uprof.shape[1])
        grid=sorted([a for a in agrid if (a>=0 if sign==1 else a<=0)],key=abs)
        for a in grid:
            def rp(zv): return whitened_residual(theta0+a*q2+Uprof@zv,design,yobs)
            fp=least_squares(rp,z,bounds=(-2.5,2.5),max_nfev=120,xtol=3e-8,ftol=3e-8,gtol=3e-8)
            z=fp.x; chi=float(fp.fun@fp.fun)
            vals.append((float(a),chi,max(0.,chi-chi_min),bool(np.any(np.isclose(np.abs(z),2.5,atol=2e-4)))))
    # dedupe a=0
    vv={a:(chi,d,hit) for a,chi,d,hit in vals}
    for a in sorted(vv):
        chi,d,hit=vv[a]; prof_rows.append((design,a,chi,d,hit))
    fit_summary.append((design,chi_min,ahat,float(np.linalg.norm(fit.fun)),bool(fit.success),int(fit.nfev)))

profdf=pd.DataFrame(prof_rows,columns=['design','a_q2','chi2','delta_chi2','nuisance_hit_bound'])
profdf.to_csv(OUT/'toy_compositional_noisy_profiles.csv',index=False)
pd.DataFrame(fit_summary,columns=['design','chi2_min','a_hat_q2','resid_norm','success','nfev']).to_csv(OUT/'toy_compositional_fit_summary.csv',index=False)

# Profile interval interpolation around minimum on grid.
def interval_from_grid(df,thr=3.84):
    d=df.sort_values('a_q2'); x=d.a_q2.to_numpy(); y=d.delta_chi2.to_numpy()
    imin=np.argmin(y); lo=None; hi=None
    for i in range(imin-1,-1,-1):
        if y[i]>=thr and y[i+1]<thr:
            lo=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    for i in range(imin,len(x)-1):
        if y[i]<thr and y[i+1]>=thr:
            hi=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
    return lo,hi
ints=[]
for d in DESIGNS:
    lo,hi=interval_from_grid(profdf[profdf.design==d]); ints.append((d,lo,hi,None if lo is None or hi is None else hi-lo))
pd.DataFrame(ints,columns=['design','profile95_lo','profile95_hi','profile95_width']).to_csv(OUT/'toy_compositional_profile_intervals.csv',index=False)

# Noise calibration in original MID coordinates via delta method at truth.
raws=[]
for label,p in groups_for_design(theta0,'pool'):
    D=len(p); H=helmert(D,full=False)
    # standard logistic-normal (delta ignored in leading-order calibration)
    G=(np.diag(p)-np.outer(p,p))@H.T
    sd=SIGMA_ILR*np.sqrt(np.diag(G@G.T))
    for i,z in enumerate(sd): raws.append((label,i,float(p[i]),float(z)))
pd.DataFrame(raws,columns=['group','component','p','delta_method_sd']).to_csv(OUT/'toy_ilr_raw_sd_calibration.csv',index=False)

summary={
    'noise_model':{'delta':DELTA,'sigma_ilr':SIGMA_ILR,'sigma_log_pool':SIGMA_POOL,'seed':20260905},
    'q1':dict(zip(names,q1.tolist())),'q2':dict(zip(names,q2.tolist())),
    'local':pd.DataFrame(local_rows,columns=['design','identifiable_dim','sigma_min','sigma_second_min','se_q1','se_q2']).to_dict(orient='records'),
    'profile_intervals':pd.DataFrame(ints,columns=['design','profile95_lo','profile95_hi','profile95_width']).replace({np.nan:None}).to_dict(orient='records')
}
(OUT/'toy_compositional_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
