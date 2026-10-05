#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import least_squares
from scipy.linalg import null_space

OUT=Path(__file__).resolve().parent
POOLS=['OAA','Cit','AKG','Suc','Fum','Glu']
NC={'OAA':4,'Cit':6,'AKG':5,'Suc':4,'Fum':4,'Glu':5}
NISO={m:2**NC[m] for m in POOLS}
offset={}; nstate=0
for m in POOLS:
    offset[m]=nstate; nstate += NISO[m]
CONST=nstate; N=nstate+1
TIMES=np.array([0.,0.1,0.2,0.5,1.,2.])

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
accoa_mix=[((0,0),.5),((0,1),.25),((1,1),.25)]
def map_v1(o):
    a,b,c,d=o
    return [((d,c,b,f,e,a),w) for (e,f),w in accoa_mix]
T_v1=trans_matrix(6,4,map_v1)
T_v2=trans_matrix(5,6,lambda x:[(x[:5],1.)])
T_v3=np.eye(32)
def map_v4(x):
    s=x[1:5]; return [(s,.5),(s[::-1],.5)]
T_v4=trans_matrix(4,5,map_v4)
def map_sym4(x): return [(x,.5),(x[::-1],.5)]
T_sym=trans_matrix(4,4,map_sym4)
ext_oaa=np.zeros(16); ext_oaa[0]=1.

def rs(m): return slice(offset[m],offset[m]+NISO[m])
def fluxes(u,w,r):
    return {'v1':u+w,'v2':u+w,'v3':u,'v4':w,'v5':w,'v6f':w+r,'v6b':r,'v7':u}
def K_values(u,w,r,C):
    F=fluxes(u,w,r)
    K=np.zeros((N,N))
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

def mdv(pool,pos):
    n=NC[pool]; pos=[p-1 for p in pos]
    M=np.zeros((len(pos)+1,N))
    for i in range(2**n):
        b=bits(i,n); M[sum(b[j] for j in pos),offset[pool]+i]+=1.
    return M
OBS=np.vstack([mdv('Glu',[1,2,3])[:-1],mdv('Glu',[1,2,3,4,5])[:-1],mdv('Cit',[2,3,4,5])[:-1]])

names=['log v3','log v4','log v6b']+[f'log C_{m}' for m in POOLS]
u0,w0,r0=5.,5.,7.5
C0={'OAA':.1,'Cit':5.,'AKG':.3,'Suc':1.,'Fum':.2,'Glu':.5}
theta0=np.array([np.log(u0),np.log(w0),np.log(r0)]+[np.log(C0[m]) for m in POOLS])
eps_pub=2/175
FAST={'OAA','AKG','Fum'}

def theta_for_eps(eps):
    th=theta0.copy(); scale=eps/eps_pub
    for i,m in enumerate(POOLS):
        if m in FAST: th[3+i]+=np.log(scale)
    return th

def unpack(th):
    u,w,r=np.exp(th[:3]); C={m:np.exp(th[3+i]) for i,m in enumerate(POOLS)}
    return u,w,r,C

def output(th):
    u,w,r,C=unpack(th); K=sparse.csr_matrix(K_values(u,w,r,C))
    traj=expm_multiply(K,Q0,start=0.0,stop=2.0,num=21,endpoint=True)
    idx=[0,1,2,5,10,20]
    return np.concatenate([OBS@traj[i] for i in idx])

def jac_fd(th,h=2e-5):
    y=output(th); J=np.empty((len(y),len(th)))
    for j in range(len(th)):
        tp=th.copy(); tm=th.copy(); tp[j]+=h; tm[j]-=h
        J[:,j]=(output(tp)-output(tm))/(2*h)
    return J

# Exact finite structural family in the quotient by global scale.
Ainv=C0['OAA']*C0['Fum']/(w0+r0)
Binv=C0['OAA']+C0['Fum']*(u0+w0+r0)/(w0+r0)
def structural_family(r):
    p=w0+r
    disc=Binv**2-4*Ainv*(p+u0)
    CO=(Binv-np.sqrt(disc))/2
    CF=Ainv*p/CO
    C=dict(C0); C['OAA']=CO; C['Fum']=CF
    th=np.array([np.log(u0),np.log(w0),np.log(r)]+[np.log(C[m]) for m in POOLS])
    return th,CO,CF
ynom=output(theta0)
inv_rows=[]
for r in [5.,6.,7.5,9.,10.,12.]:
    th,CO,CF=structural_family(r)
    inv_rows.append((r,CO,CF,float(np.max(np.abs(output(th)-ynom)))))

# Structural quotient basis.
g=np.ones(9)
s2=np.array([0.,0.,15.,10.,0.,0.,0.,-1.,0.])
Qs,_=np.linalg.qr(np.column_stack([g,s2])); Qs=Qs[:,:2]
Nq=null_space(Qs.T)  # 9x7

def quotient_modes(eps):
    J=jac_fd(theta_for_eps(eps))
    Jr=J@Nq
    F=Jr.T@Jr
    vals,vecs=np.linalg.eigh(F)
    # ascending: m2, m1, m1, O(1)...
    q2=Nq@vecs[:,0]; q1=Nq@vecs[:,1]
    q2/=np.linalg.norm(q2); q1/=np.linalg.norm(q1)
    return vals,q1,q2

vals_pub,q1_pub,q2_pub=quotient_modes(eps_pub)

def nuisance_basis(q):
    return null_space(np.vstack([Qs.T,q.reshape(1,-1)]))

sigma=.01
out_cache={}
def profile_at(eps,q,a,z0=None,max_nfev=60):
    th0=theta_for_eps(eps); y0=output(th0); U=nuisance_basis(q)
    if z0 is None: z0=np.zeros(U.shape[1])
    def resid(z): return (output(th0+a*q+U@z)-y0)/sigma
    fit=least_squares(resid,z0,bounds=(-2.,2.),max_nfev=max_nfev,xtol=1e-9,ftol=1e-9,gtol=1e-9)
    return float(fit.fun@fit.fun),fit.x,bool(np.any(np.isclose(np.abs(fit.x),2.,atol=2e-4)))

# Full nonlinear profiles at the published epsilon.
agrid=np.array([-0.5,-0.25,0.0,0.25,0.5])
prof=[]
for label,q in [('m1',q1_pub),('m2',q2_pub)]:
    vals={}; hits={}
    for sign in [1,-1]:
        z=np.zeros(6)
        grid=sorted([a for a in agrid if (a>=0 if sign==1 else a<=0)],key=abs)
        for a in grid:
            L,z,hit=profile_at(eps_pub,q,float(a),z,max_nfev=25)
            vals[float(a)]=L; hits[float(a)]=hit
    for a in agrid: prof.append((eps_pub,label,float(a),vals[float(a)],hits[float(a)]))

with open(OUT/'freeflux_structural_invariance.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['v6b','C_OAA','C_Fum','max_output_difference']); w.writerows(inv_rows)
with open(OUT/'freeflux_profile_likelihood.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['epsilon','direction','profile_coordinate','delta_chi2','nuisance_hit_bound']); w.writerows(prof)

plt.figure(figsize=(7.2,5.2))
for label,marker in [('m1','o'),('m2','s')]:
    d=np.array([(a,L) for e,lab,a,L,hb in prof if lab==label]); d=d[np.argsort(d[:,0])]
    plt.plot(d[:,0],d[:,1],marker=marker,label=label)
plt.axhline(3.84,linestyle='--',label=r'$\Delta\chi^2=3.84$')
plt.xlabel('Coordenada perfilada (norma log-paramétrica)')
plt.ylabel(r'$\Delta\chi^2$ esperado, $\sigma=0.01$')
plt.legend(); plt.tight_layout(); plt.savefig(OUT/'freeflux_profile_likelihood.png',dpi=200); plt.close()

summary=[]
summary += ['FreeFlux structural symmetry + nonlinear profile audit','='*60]
summary += [f'A={Ainv:.12g}',f'B={Binv:.12g}','Exact tangent modulo global scale: (d log v6b,d log C_OAA,d log C_Fum)=(15,10,-1)','']
summary += ['Finite exact-family checks:']
for r,CO,CF,err in inv_rows: summary.append(f' r={r:4.1f} CO={CO:.9f} CF={CF:.9f} max|dy|={err:.3e}')
summary += ['','Published-point quotient FIM eigenvalues ascending: '+np.array2string(vals_pub,precision=8)]
summary += ['q1 dominant: '+str(sorted(zip(names,q1_pub),key=lambda z:abs(z[1]),reverse=True)[:5])]
summary += ['q2 dominant: '+str(sorted(zip(names,q2_pub),key=lambda z:abs(z[1]),reverse=True)[:5])]
summary += ['', 'Noise assumption for profile: independent sigma=0.01 MID units, synthetic noiseless truth.']
for label in ['m1','m2']:
    summary.append(label+' profile: '+str([(a,L) for e,lab,a,L,hb in prof if lab==label]))
(OUT/'freeflux_structural_profile_summary.txt').write_text('\n'.join(summary))
print('\n'.join(summary))
