#!/usr/bin/env python3
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.linalg import expm, svd

OUT=Path(__file__).resolve().parent
TARGET_TIMES=np.array([0.0,0.1,0.2,0.5,1.0,2.0])
P_NAT=0.0107
POOLS=['OAA','Cit','AKG','Suc','Fum','Glu']
NC={'OAA':4,'Cit':6,'AKG':5,'Suc':4,'Fum':4,'Glu':5}
NISO={m:2**NC[m] for m in POOLS}
offset={}; nstate=0
for m in POOLS:
    offset[m]=nstate; nstate+=NISO[m]
CONST=nstate; N=nstate+1

def bits(i,n): return tuple((i>>(n-1-j))&1 for j in range(n))
def ibits(b):
    x=0
    for q in b: x=(x<<1)|int(q)
    return x

def trans_matrix(dst_n,src_n,mapper):
    T=np.zeros((2**dst_n,2**src_n))
    for j in range(2**src_n):
        for out,w in mapper(bits(j,src_n)):
            T[ibits(out),j]+=w
    return T

def bernoulli_iso(n,p=P_NAT):
    arr=np.zeros(2**n)
    for i in range(2**n):
        b=bits(i,n)
        arr[i]=np.prod([p if z else 1-p for z in b])
    return arr

_accoa=[]
for e in (0,1):
    pe=P_NAT if e else 1-P_NAT
    for f in (0,1):
        pf=P_NAT if f else 1-P_NAT
        _accoa.append(((e,f),0.5*pe*pf))
    _accoa.append(((e,1),0.25*pe))
_accoa.append(((1,1),0.25))
_acc={}
for b,w in _accoa: _acc[b]=_acc.get(b,0.0)+w
ACCOA_MIX=list(_acc.items())

def map_v1(o):
    a,b,c,d=o
    return [((d,c,b,f,e,a),w) for (e,f),w in ACCOA_MIX]
T_v1=trans_matrix(6,4,map_v1)
T_v2=trans_matrix(5,6,lambda x:[(x[:5],1.)])
T_v3=np.eye(32)
def map_v4(x):
    s=x[1:5]; return [(s,.5),(s[::-1],.5)]
T_v4=trans_matrix(4,5,map_v4)
def map_sym4(x): return [(x,.5),(x[::-1],.5)]
T_sym=trans_matrix(4,4,map_sym4)
EXT_OAA=bernoulli_iso(4)

def rs(m): return slice(offset[m],offset[m]+NISO[m])
def fluxes_from_theta(theta):
    u,w,r=np.exp(theta[:3]); return {'v1':u+w,'v2':u+w,'v3':u,'v4':w,'v5':w,'v6_f':w+r,'v6_b':r,'v7':u}
def concentrations_from_theta(theta): return {m:np.exp(theta[3+i]) for i,m in enumerate(POOLS)}
def build_K(theta):
    F=fluxes_from_theta(theta); C=concentrations_from_theta(theta); K=np.zeros((N,N))
    def tr(dst,src,v,T): K[rs(dst),rs(src)] += (v/C[dst])*T
    def out(pool,v): K[rs(pool),rs(pool)] -= (v/C[pool])*np.eye(NISO[pool])
    def ext(dst,v,d): K[rs(dst),CONST] += (v/C[dst])*d
    tr('OAA','Fum',F['v6_f'],T_sym); ext('OAA',F['v7'],EXT_OAA); out('OAA',F['v1']+F['v6_b'])
    tr('Cit','OAA',F['v1'],T_v1); out('Cit',F['v2'])
    tr('AKG','Cit',F['v2'],T_v2); out('AKG',F['v3']+F['v4'])
    tr('Suc','AKG',F['v4'],T_v4); out('Suc',F['v5'])
    tr('Fum','Suc',F['v5'],T_sym); tr('Fum','OAA',F['v6_b'],T_sym); out('Fum',F['v6_f'])
    tr('Glu','AKG',F['v3'],T_v3); out('Glu',F['v3'])
    return K

def initial_state():
    q=np.zeros(N)
    for m in POOLS: q[rs(m)]=bernoulli_iso(NC[m])
    q[CONST]=1.; return q
Q0=initial_state()
def mdv_matrix(pool,pos):
    pos=[p-1 for p in pos]; M=np.zeros((len(pos)+1,N))
    for i in range(2**NC[pool]):
        b=bits(i,NC[pool]); M[sum(b[j] for j in pos),offset[pool]+i]+=1
    return M
OBS={
'Glu_123':mdv_matrix('Glu',[1,2,3]),
'Glu_12345':mdv_matrix('Glu',[1,2,3,4,5]),
'Cit_2345':mdv_matrix('Cit',[2,3,4,5])}
theta0=np.array([np.log(5.),np.log(5.),np.log(7.5)]+[np.log(x) for x in [.1,5.,.3,1.,.2,.5]])

def exact_stack(theta):
    K=build_K(theta); states={float(t):expm(K*t)@Q0 for t in TARGET_TIMES}; vals=[]
    for emu,M in OBS.items():
        for t in TARGET_TIMES:
            x=M@states[float(t)]; vals.extend(x[:-1])
    return np.asarray(vals)

def jac_fd(theta,h):
    y0=exact_stack(theta); J=np.empty((len(y0),len(theta)))
    for j in range(len(theta)):
        tp=theta.copy(); tm=theta.copy(); tp[j]+=h; tm[j]-=h
        J[:,j]=(exact_stack(tp)-exact_stack(tm))/(2*h)
    return J

hs=[4e-5,2e-5,1e-5,5e-6,2.5e-6]
cache={}; rows=[]; tstart=time.time()
for h in hs:
    t0=time.time(); J=jac_fd(theta0,h); s=svd(J,compute_uv=False); cache[h]=J
    rows.append((h,np.linalg.norm(J),s[0],s[-1],time.time()-t0))
    print('h',h,'norm',np.linalg.norm(J),'smin',s[-1],'sec',time.time()-t0,flush=True)
pd.DataFrame(rows,columns=['h_theta','J_frobenius','sigma_max','sigma_min_raw','seconds']).to_csv(OUT/'freeflux_exact_fdm_audit.csv',index=False)
rd=[]
for h1,h2 in zip(hs[:-1],hs[1:]):
    J1,J2=cache[h1],cache[h2]
    rd.append((h1,h2,float(np.linalg.norm(J1-J2)/np.linalg.norm(J2)),float(np.max(np.abs(J1-J2)))))
pd.DataFrame(rd,columns=['h','h_over_2','rel_Frobenius_Jh_minus_Jh2','max_abs_Jh_minus_Jh2']).to_csv(OUT/'freeflux_exact_fd_step_differences.csv',index=False)
(OUT/'freeflux_exact_fdm_summary.json').write_text(json.dumps({'runtime_seconds':time.time()-tstart,'note':'Exact positional-isotopomer simulator only; source-faithful runtime derivative not rerun.'},indent=2))
print(pd.read_csv(OUT/'freeflux_exact_fd_step_differences.csv').to_string(index=False))
