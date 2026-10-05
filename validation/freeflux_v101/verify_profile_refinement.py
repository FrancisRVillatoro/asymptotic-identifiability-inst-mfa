#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
SRC=HERE/'run_block4_toy.py'
# Load model, structural directions, q2 and the matched fixed noise realization,
# but do not execute the coarse Monte Carlo/profile loops below this marker.
code=SRC.read_text().split('# Local information and Monte Carlo uncertainty')[0]
ns={'__file__':str(SRC),'__name__':'freeflux_profile_refinement'}
exec(compile(code,str(SRC),'exec'),ns)
theta0=ns['theta0']; q2=ns['q2']; g=ns['g']; s2=ns['s2']; Nbase=ns['Nbase']; noisy_obs=ns['noisy_obs']; whitened_residual=ns['whitened_residual']

TARGETS={
 'baseline':[-2.0,2.0],
 'fum':[-1.280247,2.0],
 'asp':[-0.403851,0.212244],
}

def bases(design):
    if design=='baseline':
        Q=Nbase; U=null_space(np.vstack([g,s2,q2]))
    else:
        Q=null_space(g.reshape(1,-1)); U=null_space(np.vstack([g,q2]))
    return Q,U

def global_min(design,yobs,Q):
    def r(x): return whitened_residual(theta0+Q@x,design,yobs)
    fit=least_squares(r,np.zeros(Q.shape[1]),bounds=(-2.5,2.5),max_nfev=500,xtol=2e-10,ftol=2e-10,gtol=2e-10)
    return float(fit.fun@fit.fun),fit

def eval_profile(design,a,yobs,U,start):
    def r(z): return whitened_residual(theta0+a*q2+U@z,design,yobs)
    fit=least_squares(r,start,bounds=(-5.0,5.0),max_nfev=700,xtol=3e-10,ftol=3e-10,gtol=3e-10)
    chi=float(fit.fun@fit.fun); physical=float(np.max(np.abs(a*q2+U@fit.x)))
    return chi,fit,physical

rows=[]; summaries=[]
for design,alist in TARGETS.items():
    yobs,_=noisy_obs(design); Q,U=bases(design); chi_min,gfit=global_min(design,yobs,Q)
    for a in alist:
        starts=[np.zeros(U.shape[1]),0.1*np.sin(np.arange(U.shape[1])+1.0)]
        vals=[]
        for si,start in enumerate(starts):
            chi,fit,physical=eval_profile(design,a,yobs,U,start)
            vals.append(chi)
            rows.append((design,a,si,chi,max(0.0,chi-chi_min),bool(fit.success),int(fit.status),float(fit.optimality),int(fit.nfev),float(np.max(np.abs(fit.x))),physical))
        summaries.append((design,a,chi_min,float(min(vals)-chi_min),float(abs(vals[0]-vals[1]))))

df=pd.DataFrame(rows,columns=['design','a','start','chi2','delta_chi2','success','status','optimality','nfev','z_inf','physical_inf'])
df.to_csv(HERE/'freeflux_profile_refine_points.csv',index=False)
sdf=pd.DataFrame(summaries,columns=['design','a','chi2_min','delta_chi2_best','two_start_abs_chi2_difference'])
sdf.to_csv(HERE/'freeflux_profile_refinement_summary.csv',index=False)
# Machine-readable interpretation used in v15/v9.
out={
 'baseline':{'a_minus2':float(sdf[(sdf.design=='baseline')&(sdf.a==-2.0)].delta_chi2_best.iloc[0]),'a_plus2':float(sdf[(sdf.design=='baseline')&(sdf.a==2.0)].delta_chi2_best.iloc[0])},
 'fum':{'lower_root_test_a':-1.280247,'lower_delta_chi2':float(sdf[(sdf.design=='fum')&(np.isclose(sdf.a,-1.280247))].delta_chi2_best.iloc[0]),'a_plus2_delta_chi2':float(sdf[(sdf.design=='fum')&(sdf.a==2.0)].delta_chi2_best.iloc[0])},
 'asp':{'lower_root_test_a':-0.403851,'lower_delta_chi2':float(sdf[(sdf.design=='asp')&(np.isclose(sdf.a,-0.403851))].delta_chi2_best.iloc[0]),'upper_root_test_a':0.212244,'upper_delta_chi2':float(sdf[(sdf.design=='asp')&(np.isclose(sdf.a,0.212244))].delta_chi2_best.iloc[0]),'width':0.616095},
}
(HERE/'freeflux_profile_refinement_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(sdf.to_string(index=False)); print(json.dumps(out,indent=2)); print('V101_FREEFLUX_PROFILE_REFINEMENT_COMPLETED')
