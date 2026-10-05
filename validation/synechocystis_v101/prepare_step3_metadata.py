#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import scipy.linalg as la
import syn_discrete_frechet_jacobian as sj
HERE=Path(__file__).resolve().parent; OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
alpha=sj.ns['Brel'].T@np.ones(len(sj.total_ids)); qscale=np.r_[alpha,np.ones(sj.nc)]; qscale/=la.norm(qscale); Q=la.null_space(qscale.reshape(1,-1))
rF=np.zeros(sj.P); rG=np.zeros(sj.P); rF[sj.nf+sj.pidx['F6P']]=1.0; rG[sj.nf+sj.pidx['GAP']]=1.0
rco=np.zeros(sj.P); rco[:sj.nf]=sj.L[sj.tidx['co2in'],:]/sj.fss[sj.tidx['co2in']]
np.savez_compressed(OUT/'step3_metadata.npz',nf=sj.nf,nc=sj.nc,P=sj.P,pool_ids=np.asarray(sj.pool_ids,dtype=object),qscale=qscale,Qscale=Q,rF=rF,rG=rG,rco=rco,total_ids=np.asarray(sj.total_ids,dtype=object))
print('STEP3_METADATA_WRITTEN',OUT/'step3_metadata.npz')
