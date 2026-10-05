#!/usr/bin/env python3
from pathlib import Path
import time
from collections.abc import Iterable
import numpy as np
import scipy.linalg as la

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'synechocystis_block2'
src=ROOT/'asymptotic_sweep_fast2_extended.py'
prefix=src.read_text().split("partitions={'fast2'")[0]
ns={'__file__':str(src),'__name__':'discrete_vec'}
exec(compile(prefix,str(src),'exec'),ns)
model=ns['model']; concs0=ns['concs0']; pool_ids=ns['pool_ids']; pidx=ns['pidx']; nf=ns['nf']; nc=ns['nc']; P=ns['P']
L=ns['L']; fss=ns['fss']; total_ids=ns['total_ids']; tidx=ns['tidx']; Aall_flux=ns['Aall_flux']; Ball_flux=ns['Ball_flux']
OBS=[10.,30.,60.,120.,240.]; BASE=np.array([0.,10.,30.,60.,120.,240.])

def make_grid(sub):
 g=[0.]
 for a,b in zip(BASE[:-1],BASE[1:]): g.extend(np.linspace(a,b,sub+1)[1:].tolist())
 return [float(round(x,12)) for x in g]

def conv_vals(vs):
 v=np.asarray(vs[0],float)
 for z in vs[1:]: v=np.convolve(v,np.asarray(z,float))
 return v

def conv_der_all(vs,ds):
 v=np.asarray(vs[0],float); d=np.asarray(ds[0],float)
 for v2,d2 in zip(vs[1:],ds[1:]):
  v2=np.asarray(v2,float); d2=np.asarray(d2,float)
  out=np.empty((P,len(v)+len(v2)-1))
  for q in range(P): out[q]=np.convolve(d[q],v2)+np.convolve(v,d2[q])
  d=out; v=np.convolve(v,v2)
 return d

def setup(concs,sub):
 for k,v in zip(total_ids,fss): model.total_fluxes[k]=float(v)
 for k,v in concs.items(): model.concentrations[k]=float(v)
 dts=sorted(set(round((b-a)/sub,12) for a,b in zip(BASE[:-1],BASE[1:])))
 mats={}; props={}; dprops={}
 for size in sorted(model.matrix_As):
  lambA,fa,pros=model.matrix_As[size]; lambB,fb,srcs=model.matrix_Bs[size]; lambM,mids=model.matrix_Ms[size]
  A=np.asarray(lambA(*model.total_fluxes[fa]),float); B=np.asarray(lambB(*model.total_fluxes[fb]),float); M=np.asarray(lambM(*model.concentrations[mids]),float)
  Minv=la.inv(M); F=Minv@A; Finv=la.inv(F); I=np.eye(F.shape[0])
  Af=Aall_flux[size]; Bf=Ball_flux[size]
  Ad=np.concatenate([Af,np.zeros((nc,*Af.shape[1:]))])
  Bd=np.concatenate([Bf,np.zeros((nc,*Bf.shape[1:]))])
  Md=np.zeros((P,len(pros),len(pros)))
  for i,e in enumerate(pros):
   if e.metabolite_id in pidx: Md[nf+pidx[e.metabolite_id],i,i]=concs[e.metabolite_id]
  Mid=-np.einsum('ij,qjk,kl->qil',Minv,Md,Minv,optimize=True)
  Fd=np.einsum('qij,jk->qik',Mid,A,optimize=True)+np.einsum('ij,qjk->qik',Minv,Ad,optimize=True)
  Finvd=-np.einsum('ij,qjk,kl->qil',Finv,Fd,Finv,optimize=True)
  mats[size]=(A,B,Minv,F,Finv,srcs,pros,Ad,Bd,Mid)
  for dt in dts:
   Phi=la.expm(F*dt); Gamma=(Phi-I)@Finv; Omega=(Gamma/dt-I)@Finv
   props[(size,dt)]=(Phi,Gamma,Omega)
   dp=np.empty((P,3,F.shape[0],F.shape[1]))
   for q in range(P):
    dPhi=la.expm_frechet(F*dt,Fd[q]*dt,compute_expm=False)
    dGamma=dPhi@Finv+(Phi-I)@Finvd[q]
    dOmega=(dGamma/dt)@Finv+(Gamma/dt-I)@Finvd[q]
    dp[q,0]=dPhi;dp[q,1]=dGamma;dp[q,2]=dOmega
   dprops[(size,dt)]=dp
 return mats,props,dprops

def compute(concs,sub=16):
 tsetup=time.time();mats,props,dprops=setup(concs,sub); setup_sec=time.time()-tsetup
 nom={}; Ys={0.:{}}; Xs={0.:{}}; ders={}; Yd={0.:{}}; Xd={0.:{}}
 for size in mats:
  Ys[0.][size]=np.asarray(model.initial_matrix_Ys[size]); Xs[0.][size]=np.asarray(model.initial_matrix_Xs[size])
  Yd[0.][size]=np.zeros((P,*Ys[0.][size].shape)); Xd[0.][size]=np.zeros((P,*Xs[0.][size].shape))
 def get_nom(e,t): return model.substrate_MDVs[e].value.copy() if e in model.substrate_MDVs else nom[e][t]
 def get_der(e,t): return np.zeros((P,len(get_nom(e,t)))) if e in model.substrate_MDVs else ders[e][t]
 tp=0.; grid=make_grid(sub)
 for t in grid:
  if t==0: continue
  dt=float(round(t-tp,12));Ys[t]={};Xs[t]={};Yd[t]={};Xd[t]={}
  for size in mats:
   A,B,Minv,F,Finv,srcs,pros,Ad,Bd,Mid=mats[size];Phi,Gamma,Omega=props[(size,dt)];dp=dprops[(size,dt)]
   dPhi=dp[:,0];dGamma=dp[:,1];dOmega=dp[:,2]
   Xp=Xs[tp][size];Yp=Ys[tp][size];Xdp=Xd[tp][size];Ydp=Yd[tp][size]
   Gp=Minv@B@Yp
   Gpd=(np.einsum('qij,jk,km->qim',Mid,B,Yp,optimize=True)+np.einsum('ij,qjk,km->qim',Minv,Bd,Yp,optimize=True)+np.einsum('ij,jk,qkm->qim',Minv,B,Ydp,optimize=True))
   yl=[]; ydl=[]
   for src0 in srcs:
    if not isinstance(src0,Iterable):
     yl.append(get_nom(src0,t));ydl.append(get_der(src0,t))
    else:
     vs=[get_nom(e,t) for e in src0];ds=[get_der(e,t) for e in src0]
     yl.append(conv_vals(vs));ydl.append(conv_der_all(vs,ds))
   Y=np.asarray(yl); Yq=np.stack(ydl,axis=1) # P, nsrc, mdv
   G=Minv@B@Y
   Gd=(np.einsum('qij,jk,km->qim',Mid,B,Y,optimize=True)+np.einsum('ij,qjk,km->qim',Minv,Bd,Y,optimize=True)+np.einsum('ij,jk,qkm->qim',Minv,B,Yq,optimize=True))
   DG=G-Gp
   X=Phi@Xp-Gamma@Gp-Omega@DG
   Xq=(np.einsum('qij,jm->qim',dPhi,Xp,optimize=True)+np.einsum('ij,qjm->qim',Phi,Xdp,optimize=True)-np.einsum('qij,jm->qim',dGamma,Gp,optimize=True)-np.einsum('ij,qjm->qim',Gamma,Gpd,optimize=True)-np.einsum('qij,jm->qim',dOmega,DG,optimize=True)-np.einsum('ij,qjm->qim',Omega,Gd-Gpd,optimize=True))
   Ys[t][size]=Y;Xs[t][size]=X;Yd[t][size]=Yq;Xd[t][size]=Xq
   for i,e in enumerate(pros): nom.setdefault(e,{})[t]=X[i]; ders.setdefault(e,{})[t]=Xq[:,i,:]
  tp=t
 id2={e.id:e for e in nom};targets=[id2[x] for x in model.target_EMUs]
 y=[];J=np.zeros((215,P)); row=0
 for e in targets:
  for t in OBS:
   vals=nom[e][t][:-1]; dd=ders[e][t][:,:-1]; n=len(vals)
   y.extend(vals); J[row:row+n,:]=dd.T; row+=n
 return np.asarray(y),J,setup_sec

if __name__=='__main__':
 cc=dict(concs0);cc['F6P']*=1/64;cc['GAP']*=1/64
 t=time.time();y,J,ss=compute(cc,16);print('elapsed',time.time()-t,'setup',ss,'sv',la.svdvals(J)[-6:])
