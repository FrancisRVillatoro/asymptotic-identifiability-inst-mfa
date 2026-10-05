#!/usr/bin/env python3
from pathlib import Path
import json, math, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import null_space, helmert, svd

OUT=Path(__file__).resolve().parent
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
        for out,wgt in mapper(bits(j,src_n)):
            T[ibits(out),j]+=wgt
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
    s=x[1:5]
    return [(s,.5),(s[::-1],.5)]
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
    q[CONST]=1.
    return q
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
    u,w,r,C=unpack(th)
    K=sparse.csr_matrix(K_values(u,w,r,C,acc=acc,asp=asp))
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
        tau0=C0['OAA']/(u0+w0+r0)
        ftimes=tau0*np.array([0.25,0.5,1.,2.,4.])
        fast=trajectories(th,ftimes,'unlabeled','labeled')
        for t,q in zip(ftimes,fast):
            groups.append((f'AspFast_OAA@{t:.8f}',M_OAA@q))
            groups.append((f'AspFast_Fum@{t:.8f}',M_FUM@q))
    return groups

DELTA=1e-4; SIGMA_ILR=0.08; SIGMA_POOL=0.10

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

def jac_fd(th,design,h):
    y,_=prediction(th,design); J=np.empty((len(y),len(th)))
    for j in range(len(th)):
        tp=th.copy(); tm=th.copy(); tp[j]+=h; tm[j]-=h
        J[:,j]=(prediction(tp,design)[0]-prediction(tm,design)[0])/(2*h)
    scale=np.full(len(y),SIGMA_ILR)
    if design=='pool': scale[-1]=SIGMA_POOL
    return J/scale[:,None]

g=np.ones(9); g/=np.linalg.norm(g)
s2=np.array([0.,0.,15.,10.,0.,0.,0.,-1.,0.]); s2/=np.linalg.norm(s2)
Qs,_=np.linalg.qr(np.column_stack([g,s2])); Qs=Qs[:,:2]
Nbase=null_space(Qs.T)

def quotient(design):
    if design=='baseline': return Nbase
    if design in ('fum','asp'): return null_space(g.reshape(1,-1))
    return np.eye(9)

def svd_objects(A):
    U,s,Vh=svd(A,full_matrices=False)
    tol=max(A.shape)*np.finfo(float).eps*s[0]
    keep=s>tol
    Ur=U[:,keep]; sr=s[keep]; Vhr=Vh[keep,:]
    Aplus=(Vhr.T*(1.0/sr))@Ur.T
    cov=(Vhr.T*(1.0/sr**2))@Vhr
    return U,s,Vh,keep,Aplus,cov,tol

# Reference q1/q2 at manuscript step h=2e-5.
Jref=jac_fd(theta0,'baseline',2e-5)
Ur,sr,Vhr=svd(Jref@Nbase,full_matrices=False)
q2_ref=Nbase@Vhr[-1]; q2_ref/=np.linalg.norm(q2_ref)
q1_ref=Nbase@Vhr[-2]; q1_ref/=np.linalg.norm(q1_ref)

hs=np.array([8e-5,4e-5,2e-5,1e-5,5e-6])
designs=['baseline','fum','asp','pool']
cache={}
rows=[]
start=time.time()
for h in hs:
    for d in designs:
        t0=time.time(); J=jac_fd(theta0,d,float(h)); Q=quotient(d); A=J@Q
        U,s,Vh,keep,Aplus,cov,tol=svd_objects(A)
        b1=Q.T@q1_ref; b2=Q.T@q2_ref
        se1=float(np.sqrt(max(0.,b1@cov@b1))); se2=float(np.sqrt(max(0.,b2@cov@b2)))
        # bottom singular direction for baseline only, lifted to full coordinates
        if d=='baseline':
            q=Nbase@Vh[-1]; q/=np.linalg.norm(q)
            dot=abs(float(np.dot(q,q2_ref))); dot=min(1.0,max(-1.0,dot))
            angle_deg=float(np.degrees(np.arccos(dot)))
        else:
            angle_deg=np.nan
        cache[(float(h),d)]=(J,s,se1,se2,angle_deg,tol,int(keep.sum()))
        rows.append((h,d,J.shape[0],A.shape[1],float(s[-1]),float(s[-2]),se1,se2,angle_deg,tol,int(keep.sum()),time.time()-t0))
        print(f'h={h:.1e} {d:8s} smin={s[-1]:.12g} se2={se2:.12g} angle={angle_deg:.3e} sec={time.time()-t0:.2f}',flush=True)

main=pd.DataFrame(rows,columns=['h_theta','design','nobs','q_dim','sigma_min','sigma_second_min','se_q1','se_q2','baseline_q2_angle_deg_vs_h2e-5','svd_tol','rank_kept','seconds'])
main.to_csv(OUT/'toy_svd_fdm_audit.csv',index=False)

# Step-doubling differences J(h) versus J(h/2).
rd=[]
for h1,h2 in zip(hs[:-1],hs[1:]):
    for d in designs:
        J1=cache[(float(h1),d)][0]; J2=cache[(float(h2),d)][0]
        rel=float(np.linalg.norm(J1-J2)/np.linalg.norm(J2))
        maxabs=float(np.max(np.abs(J1-J2)))
        rd.append((h1,h2,d,rel,maxabs))
pd.DataFrame(rd,columns=['h','h_over_2','design','rel_Frobenius_Jh_minus_Jh2','max_abs_Jh_minus_Jh2']).to_csv(OUT/'toy_fd_step_differences.csv',index=False)

# SVD-vs-Gram comparison at manuscript h.
comp=[]
for d in designs:
    h=2e-5; J,s,se1_svd,se2_svd,ang,tol,rk=cache[(h,d)]; Q=quotient(d); A=J@Q
    G=A.T@A
    covg=np.linalg.pinv(G)
    b1=Q.T@q1_ref; b2=Q.T@q2_ref
    se1g=float(np.sqrt(max(0.,b1@covg@b1))); se2g=float(np.sqrt(max(0.,b2@covg@b2)))
    comp.append((d,se1_svd,se1g,se2_svd,se2g,se1g/se1_svd,se2g/se2_svd))
pd.DataFrame(comp,columns=['design','se_q1_svd','se_q1_gram_pinv','se_q2_svd','se_q2_gram_pinv','gram_over_svd_q1','gram_over_svd_q2']).to_csv(OUT/'toy_svd_vs_gram.csv',index=False)

summary={
    'runtime_seconds':time.time()-start,
    'reference_h_theta':2e-5,
    'svd_cutoff':'max(A.shape)*eps_machine*sigma_max',
    'q_reference':'baseline two weakest singular directions at h_theta=2e-5',
    'files':['toy_svd_fdm_audit.csv','toy_fd_step_differences.csv','toy_svd_vs_gram.csv']
}
(OUT/'toy_svd_fdm_summary.json').write_text(json.dumps(summary,indent=2))
print('DONE',json.dumps(summary,indent=2))
