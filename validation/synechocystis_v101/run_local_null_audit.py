#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy.linalg as la
import syn_discrete_frechet_jacobian as sj
HERE=Path(__file__).resolve().parent; OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
alpha=sj.ns['Brel'].T@np.ones(len(sj.total_ids)); qscale=np.r_[alpha,np.ones(sj.nc)]; qscale/=la.norm(qscale); Q=la.null_space(qscale.reshape(1,-1))

def pfile(eps,sub=64): return OUT/f'generic_exactJ_eps{eps:.12g}_sub{sub}.npz'
rows=[]
for label,eps in [('generic_nominal_eps1',1.0),('generic_nominal_eps1_64',1/64),('generic_nominal_eps1_512',1/512)]:
    d=np.load(pfile(eps)); J=d['J']; y=d['y']; U,s,Vh=la.svd(J,full_matrices=False); tol=max(J.shape)*np.finfo(float).eps*s[0]; sq=la.svdvals(J@Q)
    rows.append(dict(label=label,epsilon=eps,subdivisions=64,rank=int(np.sum(s>tol)),sigma_max=s[0],sigma_second_smallest=s[-2],sigma_null=s[-1],null_to_second=s[-1]/s[-2],null_scale_alignment=abs(float(Vh[-1]@qscale)),scale_residual_rel=float(la.norm(J@qscale)/la.norm(J,2)),quotient_sigma_min=sq[-1],standard_rank=int(np.linalg.matrix_rank(J)),max_abs_log_concentration_shift=np.nan,output_relative_change_vs_nominal=np.nan,output_max_abs_change_vs_nominal=np.nan))
# Independent deterministic 5% multipool concentration perturbation at epsilon=1.
nom=dict(sj.concs0); nom['PG']=11.; nom['Gc']=9.
rng=np.random.default_rng(20262026); z=rng.uniform(-1,1,len(sj.pool_ids)); z*=0.05/np.max(np.abs(z)); cc=dict(nom)
for i,m in enumerate(sj.pool_ids): cc[m]*=float(np.exp(z[i]))
y,J,_=sj.compute(cc,sub=16); U,s,Vh=la.svd(J,full_matrices=False); tol=max(J.shape)*np.finfo(float).eps*s[0]; sq=la.svdvals(J@Q)
ref=np.load(pfile(1.0))['y']
row=dict(label='nearby_random_concentration_5pct_eps1',epsilon=1.0,subdivisions=16,rank=int(np.sum(s>tol)),sigma_max=s[0],sigma_second_smallest=s[-2],sigma_null=s[-1],null_to_second=s[-1]/s[-2],null_scale_alignment=abs(float(Vh[-1]@qscale)),scale_residual_rel=float(la.norm(J@qscale)/la.norm(J,2)),quotient_sigma_min=sq[-1],standard_rank=int(np.linalg.matrix_rank(J)),max_abs_log_concentration_shift=float(np.max(np.abs(z))),output_relative_change_vs_nominal=float(la.norm(y-ref)/la.norm(ref)),output_max_abs_change_vs_nominal=float(np.max(np.abs(y-ref))))
rows.append(row)
np.savez_compressed(OUT/'nearby_random_concentration_5pct_eps1_sub16.npz',y=y,J=J,s=s,Vh=Vh,qscale=qscale,pool_ids=np.asarray(sj.pool_ids,dtype=object),concentrations=np.asarray([cc[m] for m in sj.pool_ids]),log_shifts=z)
df=pd.DataFrame(rows); df.to_csv(OUT/'local_null_fullrank_checks.csv',index=False)
summary={
 'all_rank_59':bool((df['rank']==59).all()),
 'min_null_scale_alignment':float(df.null_scale_alignment.min()),
 'max_scale_residual_rel':float(df.scale_residual_rel.max()),
 'nearby_output_relative_change':float(row['output_relative_change_vs_nominal']),
 'nearby_max_abs_log_concentration_shift':float(row['max_abs_log_concentration_shift']),
 'global_structural_identifiability_claim':False,
}
(OUT/'local_null_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(df.to_string(index=False)); print(json.dumps(summary,indent=2)); print('V101_LOCAL_NULL_AUDIT_COMPLETED')
