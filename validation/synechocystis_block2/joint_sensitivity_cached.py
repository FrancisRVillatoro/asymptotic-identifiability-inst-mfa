#!/usr/bin/env python3
import ast, importlib.util, json, time, warnings
from collections.abc import Iterable
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.linalg as la
import openpyxl
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parent

spec=importlib.util.spec_from_file_location('ff',ROOT/'source_faithful_freeflux_runtime.py')
ff=importlib.util.module_from_spec(spec); spec.loader.exec_module(ff)
tr=ast.parse((ROOT/'run_syn_nominal.py').read_text()); fluxes=None
for node in tr.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='fluxes' for t in node.targets):
        fluxes=ast.literal_eval(node.value); break
wb=openpyxl.load_workbook(ROOT/'concentrations.xlsx',data_only=True)
concs={r[0]:float(r[1]) for r in list(wb.active.iter_rows(values_only=True))[1:] if r[0] and r[1] is not None}
pool_ids=sorted(concs); pidx={m:i for i,m in enumerate(pool_ids)}

model=ff.Model('syn'); model.read_from_file(ROOT/'reactions.tsv')
Stot=model._full_total_stoichiometric_matrix(tuple(model.metabolites),tuple(model.reactions)).drop(index=model.end_substrates+model.end_products)
total_ids=list(model.totalfluxids); tidx={k:i for i,k in enumerate(total_ids)}
f0=np.array([fluxes[k] for k in total_ids],float)
Nabs=la.null_space(Stot.values); fss=Nabs@(Nabs.T@f0)
Brel=la.null_space(Stot.values@np.diag(fss)); L=np.diag(fss)@Brel
nf=L.shape[1]; nc=len(pool_ids); P=nf+nc

sim=ff.InstSimulator(model)
sim.set_target_EMUs({'G3P':['23','123'],'DHAP':'123','PEP':'123','Fum':'1234','R5P':'12345','RuBP':'12345','Cit':['12345','123456'],'S7P':'1234567'})
sim.set_labeling_strategy('CO2.ex',labeling_pattern=['1'],percentage=[0.5],purity=[0.997],label_atom='C')
for k,v in zip(total_ids,fss): sim.set_flux(k,float(v))
for k,v in concs.items(): sim.set_concentration(k,float(v))
sim.set_timepoints([0,10,30,60,120,240])
t0=time.time(); sim.prepare(); prep=time.time()-t0

# parameter derivatives of A/B/M (all parameter columns, cheap and exact)
Aall={}; Ball={}; Mall={}
for size,EAM in model.EAMs.items():
    lambA,fluxidsA,productEMUs=model.matrix_As[size]
    lambB,fluxidsB,sourceEMUs=model.matrix_Bs[size]
    cA=[]
    for i in range(len(fluxidsA)):
        args=np.zeros(len(fluxidsA)); args[i]=1.; cA.append(np.asarray(lambA(*args),float))
    cB=[]
    for i in range(len(fluxidsB)):
        args=np.zeros(len(fluxidsB)); args[i]=1.; cB.append(np.asarray(lambB(*args),float))
    cA=np.asarray(cA); cB=np.asarray(cB)
    LA=L[[tidx[x] for x in fluxidsA],:] if fluxidsA else np.zeros((0,nf))
    LB=L[[tidx[x] for x in fluxidsB],:] if fluxidsB else np.zeros((0,nf))
    Af=np.einsum('iab,ik->kab',cA,LA) if len(fluxidsA) else np.zeros((nf,len(productEMUs),len(productEMUs)))
    Bf=np.einsum('iab,ik->kab',cB,LB) if len(fluxidsB) else np.zeros((nf,len(productEMUs),len(sourceEMUs)))
    Aall[size]=np.concatenate([Af,np.zeros((nc,*Af.shape[1:]))],axis=0)
    Ball[size]=np.concatenate([Bf,np.zeros((nc,*Bf.shape[1:]))],axis=0)
    Md=np.zeros((P,len(productEMUs),len(productEMUs)))
    for i,emu in enumerate(productEMUs):
        if emu.metabolite_id in pidx: Md[nf+pidx[emu.metabolite_id],i,i]=concs[emu.metabolite_id]
    Mall[size]=Md

# Nominal forward pass and cached propagators.
nomvals={}; Ys={0.0:{}}; Xs={0.0:{}}; cache={}
for size in sorted(model.matrix_As):
    Ys[0.0][size]=model.initial_matrix_Ys[size]
    Xs[0.0][size]=model.initial_matrix_Xs[size]

def get_nom(emu,t):
    if emu in model.substrate_MDVs: return model.substrate_MDVs[emu].value.copy()
    return nomvals[emu][t]

def conv_vals(vs):
    v=np.asarray(vs[0],float)
    for z in vs[1:]: v=np.convolve(v,np.asarray(z,float))
    return v

tprev=0.0
for t in model.timepoints:
    if t==0: continue
    dt=t-tprev; Ys[t]={}; Xs[t]={}
    for size in sorted(model.matrix_As):
        lambA,fluxidsA,productEMUs=model.matrix_As[size]
        lambB,fluxidsB,sourceEMUs=model.matrix_Bs[size]
        lambM,metabids=model.matrix_Ms[size]
        A=np.asarray(lambA(*model.total_fluxes[fluxidsA]),float)
        B=np.asarray(lambB(*model.total_fluxes[fluxidsB]),float)
        M=np.asarray(lambM(*model.concentrations[metabids]),float)
        Minv=la.pinv(M); F=Minv@A; Finv=la.pinv(F); I=np.eye(F.shape[0])
        Phi=la.expm(F*dt); Gamma=(Phi-I)@Finv; Omega=(Gamma/dt-I)@Finv
        Xprev=Xs[tprev][size]; Yprev=Ys[tprev][size]; Gprev=Minv@B@Yprev
        ylist=[]
        for src in sourceEMUs:
            if not isinstance(src,Iterable): v=get_nom(src,t)
            else: v=conv_vals([get_nom(e,t) for e in src])
            ylist.append(v)
        Y=np.asarray(ylist,float); G=Minv@B@Y
        X=Phi@Xprev-Gamma@Gprev-Omega@(G-Gprev)
        Ys[t][size]=Y; Xs[t][size]=X
        cache[(t,size)]={'A':A,'B':B,'Minv':Minv,'Phi':Phi,'Gamma':Gamma,'Omega':Omega,
                         'Xprev':Xprev,'Yprev':Yprev,'X':X,'Y':Y,'sourceEMUs':sourceEMUs,'productEMUs':productEMUs}
        for i,e in enumerate(productEMUs): nomvals.setdefault(e,{})[t]=X[i]
    tprev=t

# map target id -> EMU object
id2emu={e.id:e for e in nomvals}
targets=[id2emu[x] for x in model.target_EMUs]
obs_labels=[]; y=[]
for e in targets:
    for t in [10.,30.,60.,120.,240.]:
        v=nomvals[e][t]
        for m in range(len(v)-1): obs_labels.append(f'{e.id}|t={t:g}|M{m}'); y.append(v[m])
y=np.asarray(y,float); J=np.zeros((len(y),P))

# scalar derivative convolution
def conv_der(vs,ds):
    v=np.asarray(vs[0],float); d=np.asarray(ds[0],float)
    for v2,d2 in zip(vs[1:],ds[1:]):
        v2=np.asarray(v2,float); d2=np.asarray(d2,float)
        d=np.convolve(d,v2)+np.convolve(v,d2); v=np.convolve(v,v2)
    return d

# One parameter at a time, reusing nominal propagators: no repeated expm.
start=time.time()
for q in range(P):
    ders={}; Yd={0.0:{}}; Xd={0.0:{}}
    for size in sorted(model.matrix_As):
        Yd[0.0][size]=np.zeros_like(model.initial_matrix_Ys[size])
        Xd[0.0][size]=np.zeros_like(model.initial_matrix_Xs[size])
    tprev=0.0
    for t in model.timepoints:
        if t==0: continue
        Yd[t]={}; Xd[t]={}
        for size in sorted(model.matrix_As):
            C=cache[(t,size)]; A=C['A']; B=C['B']; Minv=C['Minv']; Phi=C['Phi']; Gamma=C['Gamma']; Omega=C['Omega']
            Xprev=C['Xprev']; Yprev=C['Yprev']; X=C['X']; Y=C['Y']; srcs=C['sourceEMUs']; pros=C['productEMUs']
            Ad=Aall[size][q]; Bd=Ball[size][q]; Md=Mall[size][q]; Minvd=-Minv@Md@Minv
            Xdp=Xd[tprev][size]; Ydp=Yd[tprev][size]
            Hprev=(Minvd@A@Xprev + Minv@Ad@Xprev - Minv@B@Ydp - Minvd@B@Yprev - Minv@Bd@Yprev)
            ydlist=[]
            for src in srcs:
                if not isinstance(src,Iterable): d=np.zeros_like(get_nom(src,t))
                else:
                    vs=[]; ds=[]
                    for e in src:
                        vs.append(get_nom(e,t)); ds.append(ders[e][t] if e not in model.substrate_MDVs else np.zeros_like(get_nom(e,t)))
                    d=conv_der(vs,ds)
                ydlist.append(d)
            Ydq=np.asarray(ydlist,float)
            H=(Minvd@A@X + Minv@Ad@X - Minv@B@Ydq - Minvd@B@Y - Minv@Bd@Y)
            Xdq=Phi@Xdp+Gamma@Hprev+Omega@(H-Hprev)
            Yd[t][size]=Ydq; Xd[t][size]=Xdq
            for i,e in enumerate(pros): ders.setdefault(e,{})[t]=Xdq[i]
        tprev=t
    rr=[]
    for e in targets:
        for t in [10.,30.,60.,120.,240.]: rr.extend(ders[e][t][:-1])
    J[:,q]=rr
    if (q+1)%10==0 or q==P-1: print('param',q+1,'/',P,flush=True)
deriv_s=time.time()-start

U,s,Vh=la.svd(J,full_matrices=False)
alpha=Brel.T@np.ones(len(total_ids)); qscale=np.r_[alpha,np.ones(nc)]; qscale/=la.norm(qscale)
scale_res=la.norm(J@qscale)
weak=[]
for kk in range(min(15,len(s))):
    idx=len(s)-1-kk; qv=Vh[idx]; qf=Brel@qv[:nf]; qc=qv[nf:]
    weak.append({'asc':kk+1,'sigma':float(s[idx]),
                 'top_flux_rel':[(a,float(b)) for a,b in sorted(zip(total_ids,qf),key=lambda z:abs(z[1]),reverse=True)[:8]],
                 'top_logC':[(a,float(b)) for a,b in sorted(zip(pool_ids,qc),key=lambda z:abs(z[1]),reverse=True)[:8]]})
summary={'stoich_shape':list(Stot.shape),'stoich_rank':int(np.linalg.matrix_rank(Stot.values)),'directional_fluxes':len(total_ids),'flux_dof':nf,'pools':nc,'joint_dim':P,'obs_dim':J.shape[0],
'nominal_balance_residual_norm':float(la.norm(Stot.values@f0)),'nominal_balance_residual_max':float(np.max(np.abs(Stot.values@f0))),'nominal_balance_residual_metabolite':str(Stot.index[np.argmax(np.abs(Stot.values@f0))]),
'steady_state_projection_rel':float(la.norm(fss-f0)/la.norm(f0)),'steady_state_projection_max_abs':float(np.max(np.abs(fss-f0))),'projected_balance_residual_norm':float(la.norm(Stot.values@fss)),
'prepare_seconds':prep,'derivative_seconds':deriv_s,'sigma_max':float(s[0]),'sigma_min':float(s[-1]),'rank_rel_1e8':int(np.sum(s>s[0]*1e-8)),'rank_rel_1e10':int(np.sum(s>s[0]*1e-10)),'rank_rel_1e12':int(np.sum(s>s[0]*1e-12)),'rank_rel_1e14':int(np.sum(s>s[0]*1e-14)),'global_scale_residual':float(scale_res),'weak_modes':weak}
pd.DataFrame(Stot).to_csv(ROOT/'stoichiometric_total_balanced.csv')
pd.DataFrame(Brel,index=total_ids,columns=[f'u{i+1:02d}' for i in range(nf)]).to_csv(ROOT/'flux_relative_null_basis.csv')
pd.DataFrame({'flux_id':total_ids,'distributed_nominal':f0,'steady_state_projected':fss,'correction':fss-f0}).to_csv(ROOT/'steady_state_projection.csv',index=False)
pd.DataFrame(J,index=obs_labels,columns=[f'flux_tangent_{i+1:02d}' for i in range(nf)]+[f'logC_{m}' for m in pool_ids]).to_csv(ROOT/'joint_sensitivity_nominal.csv')
pd.DataFrame({'index_desc':np.arange(1,len(s)+1),'sigma':s,'fim_eigenvalue':s*s}).to_csv(ROOT/'joint_spectrum_nominal.csv',index=False)
np.savez_compressed(ROOT/'joint_sensitivity_nominal.npz',J=J,y=y,singular_values=s,Vh=Vh,Brel=Brel,L=L,f0=f0,fss=fss,total_ids=np.array(total_ids,dtype=object),pool_ids=np.array(pool_ids,dtype=object),labels=np.array(obs_labels,dtype=object))
(ROOT/'joint_sensitivity_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='weak_modes'},indent=2))
print('weakest',s[-15:])
for w in weak[:8]: print(w)
