from pathlib import Path
import numpy as np, pandas as pd
from scipy.linalg import null_space, pinv
OUT=Path(__file__).resolve().parent
src=OUT/'run_block4_toy.py'
code=src.read_text()
ns={'__file__':str(src),'__name__':'block4_toy_common'}
exec(compile(code.split('# Master synthetic noisy observations')[0],str(src),'exec'),ns)
prediction=ns['prediction']; jac_fd=ns['jac_fd']; theta0=ns['theta0']; q2=ns['q2']; g=ns['g']; s2=ns['s2']; SIGMA_ILR=ns['SIGMA_ILR']; SIGMA_POOL=ns['SIGMA_POOL']
rng=np.random.default_rng(20260905)
DESIGNS=['baseline','fum','asp']
obs={}
for d in DESIGNS:
 y,_=prediction(theta0,d); obs[d]=y+SIGMA_ILR*rng.standard_normal(len(y))
agrid=np.array([-1.5,-1.0,-.75,-.5,-.3,-.2,-.1,0,.1,.2,.3,.5,.75,1.0,1.5])
rows=[]
for d in DESIGNS:
 if d=='baseline': U=null_space(np.vstack([g,s2,q2]))
 else: U=null_space(np.vstack([g,q2]))
 Jw,_,_=jac_fd(theta0,d); AU=Jw@U; P=np.eye(AU.shape[0])-AU@pinv(AU)
 vals=[]
 for a in agrid:
  yp,_=prediction(theta0+a*q2,d); r=(yp-obs[d])/SIGMA_ILR; chi=float(np.linalg.norm(P@r)**2); vals.append((a,chi))
 cmin=min(c for a,c in vals)
 for a,c in vals: rows.append((d,a,c,c-cmin))
 print(d,'min',min(vals,key=lambda z:z[1]),'end',vals[0],vals[-1])
df=pd.DataFrame(rows,columns=['design','a_q2','chi2_gn_profile','delta_chi2']); df.to_csv(OUT/'toy_compositional_noisy_gn_profiles.csv',index=False)
# intervals
def interval(d,thr=3.84):
 x=d.a_q2.to_numpy(); y=d.delta_chi2.to_numpy(); im=np.argmin(y); lo=hi=None
 for i in range(im-1,-1,-1):
  if y[i]>=thr and y[i+1]<thr:
   lo=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
 for i in range(im,len(x)-1):
  if y[i]<thr and y[i+1]>=thr:
   hi=x[i]+(thr-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]); break
 return lo,hi
ints=[]
for design in DESIGNS:
 d=df[df.design==design].sort_values('a_q2'); lo,hi=interval(d); ints.append((design,lo,hi,None if lo is None or hi is None else hi-lo))
pd.DataFrame(ints,columns=['design','lo95','hi95','width95']).to_csv(OUT/'toy_compositional_noisy_gn_profile_intervals.csv',index=False)
print(pd.DataFrame(ints,columns=['design','lo95','hi95','width95']))
