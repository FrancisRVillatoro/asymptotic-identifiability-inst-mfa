#!/usr/bin/env python3
import ast, importlib.util, json, time, warnings
from collections.abc import Iterable
from pathlib import Path
import numpy as np, pandas as pd, scipy.linalg as la, openpyxl
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ff',ROOT/'source_faithful_freeflux_runtime.py')
ff=importlib.util.module_from_spec(spec); spec.loader.exec_module(ff)
# inputs
tr=ast.parse((ROOT/'run_syn_nominal.py').read_text()); fluxes=None
for node in tr.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='fluxes' for t in node.targets): fluxes=ast.literal_eval(node.value); break
wb=openpyxl.load_workbook(ROOT/'concentrations.xlsx',data_only=True)
concs0={r[0]:float(r[1]) for r in list(wb.active.iter_rows(values_only=True))[1:] if r[0] and r[1] is not None}
pool_ids=sorted(concs0); pidx={m:i for i,m in enumerate(pool_ids)}
model=ff.Model('syn'); model.read_from_file(ROOT/'reactions.tsv')
Stot=model._full_total_stoichiometric_matrix(tuple(model.metabolites),tuple(model.reactions)).drop(index=model.end_substrates+model.end_products)
total_ids=list(model.totalfluxids); tidx={k:i for i,k in enumerate(total_ids)}
f0=np.array([fluxes[k] for k in total_ids]); Nabs=la.null_space(Stot.values); fss=Nabs@(Nabs.T@f0)
Brel=la.null_space(Stot.values@np.diag(fss)); L=np.diag(fss)@Brel
nf=L.shape[1]; nc=len(pool_ids); P=nf+nc
# turnover table
outflow=(-np.minimum(Stot.values,0))@fss
turn=outflow/np.array([concs0[m] for m in Stot.index])
tdf=pd.DataFrame({'metabolite':Stot.index,'outflow':outflow,'turnover':turn}).sort_values('turnover').reset_index(drop=True)
tdf['log10_turnover']=np.log10(tdf.turnover); tdf['gap_to_next']=tdf.log10_turnover.shift(-1)-tdf.log10_turnover
tdf.to_csv(ROOT/'turnover_spectrum.csv',index=False)
# simulator once
sim=ff.InstSimulator(model); sim.set_target_EMUs({'G3P':['23','123'],'DHAP':'123','PEP':'123','Fum':'1234','R5P':'12345','RuBP':'12345','Cit':['12345','123456'],'S7P':'1234567'})
sim.set_labeling_strategy('CO2.ex',labeling_pattern=['1'],percentage=[0.5],purity=[0.997],label_atom='C')
for k,v in zip(total_ids,fss): sim.set_flux(k,float(v))
for k,v in concs0.items(): sim.set_concentration(k,float(v))
sim.set_timepoints([0,10,30,60,120,240]); sim.prepare()
# A/B derivatives independent of concentrations
Aall_flux={}; Ball_flux={}
for size,EAM in model.EAMs.items():
    lambA,fluxidsA,pros=model.matrix_As[size]; lambB,fluxidsB,srcs=model.matrix_Bs[size]
    cA=[]
    for i in range(len(fluxidsA)):
        a=np.zeros(len(fluxidsA)); a[i]=1.; cA.append(np.asarray(lambA(*a),float))
    cB=[]
    for i in range(len(fluxidsB)):
        a=np.zeros(len(fluxidsB)); a[i]=1.; cB.append(np.asarray(lambB(*a),float))
    cA=np.asarray(cA); cB=np.asarray(cB)
    LA=L[[tidx[x] for x in fluxidsA],:] if fluxidsA else np.zeros((0,nf)); LB=L[[tidx[x] for x in fluxidsB],:] if fluxidsB else np.zeros((0,nf))
    Aall_flux[size]=np.einsum('iab,ik->kab',cA,LA) if len(fluxidsA) else np.zeros((nf,len(pros),len(pros)))
    Ball_flux[size]=np.einsum('iab,ik->kab',cB,LB) if len(fluxidsB) else np.zeros((nf,len(pros),len(srcs)))

def conv_vals(vs):
    v=np.asarray(vs[0],float)
    for z in vs[1:]: v=np.convolve(v,np.asarray(z,float))
    return v

def conv_der(vs,ds):
    v=np.asarray(vs[0],float); d=np.asarray(ds[0],float)
    for v2,d2 in zip(vs[1:],ds[1:]):
        v2=np.asarray(v2,float); d2=np.asarray(d2,float); d=np.convolve(d,v2)+np.convolve(v,d2); v=np.convolve(v,v2)
    return d

def compute_J(concs):
    for k,v in concs.items(): model.concentrations[k]=float(v)
    # parameter derivatives for current log concentrations
    Aall={}; Ball={}; Mall={}
    for size,EAM in model.EAMs.items():
        pros=model.matrix_As[size][2]; srcs=model.matrix_Bs[size][2]
        Af=Aall_flux[size]; Bf=Ball_flux[size]
        Aall[size]=np.concatenate([Af,np.zeros((nc,*Af.shape[1:]))])
        Ball[size]=np.concatenate([Bf,np.zeros((nc,*Bf.shape[1:]))])
        Md=np.zeros((P,len(pros),len(pros)))
        for i,e in enumerate(pros):
            if e.metabolite_id in pidx: Md[nf+pidx[e.metabolite_id],i,i]=concs[e.metabolite_id]
        Mall[size]=Md
    # nominal cache
    nom={}; Ys={0.:{}}; Xs={0.:{}}; cache={}
    for size in sorted(model.matrix_As): Ys[0.][size]=model.initial_matrix_Ys[size]; Xs[0.][size]=model.initial_matrix_Xs[size]
    def get_nom(e,t): return model.substrate_MDVs[e].value.copy() if e in model.substrate_MDVs else nom[e][t]
    tp=0.
    for t in model.timepoints:
        if t==0: continue
        dt=t-tp; Ys[t]={}; Xs[t]={}
        for size in sorted(model.matrix_As):
            lambA,fa,pros=model.matrix_As[size]; lambB,fb,srcs=model.matrix_Bs[size]; lambM,mids=model.matrix_Ms[size]
            A=np.asarray(lambA(*model.total_fluxes[fa]),float); B=np.asarray(lambB(*model.total_fluxes[fb]),float); M=np.asarray(lambM(*model.concentrations[mids]),float)
            Minv=la.pinv(M); F=Minv@A; Finv=la.pinv(F); I=np.eye(F.shape[0]); Phi=la.expm(F*dt); Gamma=(Phi-I)@Finv; Omega=(Gamma/dt-I)@Finv
            Xp=Xs[tp][size]; Yp=Ys[tp][size]; Gp=Minv@B@Yp
            yl=[]
            for src in srcs: yl.append(get_nom(src,t) if not isinstance(src,Iterable) else conv_vals([get_nom(e,t) for e in src]))
            Y=np.asarray(yl); G=Minv@B@Y; X=Phi@Xp-Gamma@Gp-Omega@(G-Gp)
            Ys[t][size]=Y; Xs[t][size]=X; cache[(t,size)]=(A,B,Minv,Phi,Gamma,Omega,Xp,Yp,X,Y,srcs,pros)
            for i,e in enumerate(pros): nom.setdefault(e,{})[t]=X[i]
        tp=t
    id2={e.id:e for e in nom}; targets=[id2[x] for x in model.target_EMUs]
    labels=[]; y=[]
    for e in targets:
        for t in [10.,30.,60.,120.,240.]:
            for m,z in enumerate(nom[e][t][:-1]): labels.append(f'{e.id}|t={t:g}|M{m}'); y.append(z)
    J=np.zeros((len(y),P))
    for q in range(P):
        ders={}; Yd={0.:{}}; Xd={0.:{}}
        for size in sorted(model.matrix_As): Yd[0.][size]=np.zeros_like(model.initial_matrix_Ys[size]); Xd[0.][size]=np.zeros_like(model.initial_matrix_Xs[size])
        tp=0.
        for t in model.timepoints:
            if t==0: continue
            Yd[t]={}; Xd[t]={}
            for size in sorted(model.matrix_As):
                A,B,Minv,Phi,Gamma,Omega,Xp,Yp,X,Y,srcs,pros=cache[(t,size)]
                Ad=Aall[size][q]; Bd=Ball[size][q]; Md=Mall[size][q]; Mid=-Minv@Md@Minv
                Xdp=Xd[tp][size]; Ydp=Yd[tp][size]
                Hp=Mid@A@Xp+Minv@Ad@Xp-Minv@B@Ydp-Mid@B@Yp-Minv@Bd@Yp
                ydl=[]
                for src in srcs:
                    if not isinstance(src,Iterable): d=np.zeros_like(get_nom(src,t))
                    else:
                        vs=[]; ds=[]
                        for e in src: vs.append(get_nom(e,t)); ds.append(np.zeros_like(get_nom(e,t)) if e in model.substrate_MDVs else ders[e][t])
                        d=conv_der(vs,ds)
                    ydl.append(d)
                Yq=np.asarray(ydl); H=Mid@A@X+Minv@Ad@X-Minv@B@Yq-Mid@B@Y-Minv@Bd@Y
                Xq=Phi@Xdp+Gamma@Hp+Omega@(H-Hp); Yd[t][size]=Yq; Xd[t][size]=Xq
                for i,e in enumerate(pros): ders.setdefault(e,{})[t]=Xq[i]
            tp=t
        rr=[]
        for e in targets:
            for t in [10.,30.,60.,120.,240.]: rr.extend(ders[e][t][:-1])
        J[:,q]=rr
    return np.asarray(y),J,labels


fast=['F6P','GAP']
epsvals=np.array([1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125])
# quotient out exact global scaling direction
alpha=Brel.T@np.ones(len(total_ids))
qscale=np.r_[alpha,np.ones(nc)]; qscale=qscale/la.norm(qscale)
Q=la.null_space(qscale.reshape(1,-1))  # 60 x 59
# algebraic log-pool measurement rows
rF=np.zeros(P); rF[nf+pidx['F6P']]=1.0
rG=np.zeros(P); rG[nf+pidx['GAP']]=1.0
designs={'baseline':[], '+F6P':[], '+GAP':[], '+F6P+GAP':[]}
rows=[]
for eps in epsvals:
    cc=dict(concs0)
    for m in fast: cc[m]=concs0[m]*eps
    y,J,labels=compute_J(cc)
    mats={
      'baseline':J@Q,
      '+F6P':np.vstack([J,rF])@Q,
      '+GAP':np.vstack([J,rG])@Q,
      '+F6P+GAP':np.vstack([J,rF,rG])@Q,
    }
    for name,A in mats.items():
        sv=la.svdvals(A); designs[name].append(sv)
        for i,z in enumerate(sv): rows.append((eps,name,i+1,z,z*z))
    print('eps',eps,'weak', {k:np.asarray(v[-1])[-4:].tolist() for k,v in designs.items()}, flush=True)
# fits last 4 eps
fits=[]
x=np.log(epsvals[-4:])
for name,arr in designs.items():
    A=np.asarray(arr)
    for j in range(A.shape[1]):
        slope=np.polyfit(x,np.log(A[-4:,j]),1)[0]
        fits.append((name,j+1,float(slope),float(2*slope),float(A[0,j]),float(A[-1,j])))
pd.DataFrame(rows,columns=['epsilon','design','index_desc','sigma','fim_eigenvalue']).to_csv(ROOT/'fast2_design_lifting_sweep.csv',index=False)
pd.DataFrame(fits,columns=['design','index_desc','sigma_exponent','fim_exponent','sigma_eps1','sigma_epsmin']).to_csv(ROOT/'fast2_design_lifting_exponents.csv',index=False)
summary={}
for name,arr in designs.items():
    A=np.asarray(arr); df=pd.DataFrame([r for r in fits if r[0]==name],columns=['design','index_desc','sigma_exponent','fim_exponent','sigma_eps1','sigma_epsmin'])
    summary[name]={'weakest_eps1':A[0,-4:].tolist(),'weakest_epsmin':A[-1,-4:].tolist(),'bottom_fits':df.tail(6).to_dict(orient='records')}
(ROOT/'fast2_design_lifting_summary.json').write_text(json.dumps(summary,indent=2))
for name in designs:
    print('\n',name)
    df=pd.DataFrame([r for r in fits if r[0]==name],columns=['design','index_desc','sigma_exponent','fim_exponent','sigma_eps1','sigma_epsmin'])
    print(df.tail(8).to_string(index=False))
