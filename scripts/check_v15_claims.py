#!/usr/bin/env python3
from pathlib import Path
import json,sys,math
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SYN=ROOT/'validation/synechocystis_v101/results'
FF=ROOT/'validation/freeflux_v101'
issues=[]
def check(cond,msg):
    if cond: print('PASS',msg)
    else: print('FAIL',msg); issues.append(msg)
def close(x,y,rtol=0.02,atol=0): return abs(x-y)<=atol+rtol*abs(y)
# Synechocystis exact-J genericization and nullspace.
chg=json.loads((SYN/'genericization_output_change.json').read_text())
check(close(chg['relative_2norm'],1.85449e-6,rtol=0.05),'genericization observable change ~1.85e-6')
r=pd.read_csv(SYN/'generic_rank_audit.csv')
check((r['rank']==59).all(),'rank 59 at eps=1,1/64,1/512')
check((r['null_scale_alignment']>0.999999999).all(),'null aligns with global scale')
check((r['scale_residual_rel']<1e-11).all(),'analytic scale residual near machine precision')
# Unweighted deep slopes used in main text.
e=pd.read_csv(SYN/'generic_mode_exponents.csv').sort_values('deep_proj_F6P_GAP',ascending=False).head(2)
slopes=sorted(e['slope4'].to_list())
check(all(0.985<s<1.015 for s in slopes),'two F6P/GAP modes have m=1 slopes')
# Weighted finite-noise order check.
we=pd.read_csv(SYN/'weighted_design_mode_exponents.csv')
wb=we[(we.design=='baseline')&(we['class']=='m=1')]
check(len(wb)==2 and np.all((wb.slope4>0.98)&(wb.slope4<1.02)),'weighted baseline retains two m=1 modes')
# Local-null independent nearby audit.
ln=pd.read_csv(SYN/'local_null_fullrank_checks.csv')
near=ln[ln.label=='nearby_random_concentration_5pct_eps1'].iloc[0]
check(int(near['rank'])==59,'nearby 5% concentration perturbation rank 59')
check(close(float(near['output_relative_change_vs_nominal']),4.083601e-3,rtol=0.03),'nearby output change ~4.08e-3')
check(float(near['null_scale_alignment'])>0.999999999,'nearby null remains scale-aligned')
# Physical epsilon=1 finite-noise control.
pd1=pd.read_csv(SYN/'finite_noise_local_eps1.csv').set_index('design')
check(close(pd1.loc['cond_baseline','se_logC_F6P'],683.37,rtol=0.03),'physical eps=1 baseline F6P SE ~683')
check(close(pd1.loc['cond_baseline','se_logC_GAP'],634.42,rtol=0.03),'physical eps=1 baseline GAP SE ~634')
check(close(pd1.loc['cond_both','se_logC_F6P'],0.10,rtol=0.02) and close(pd1.loc['cond_both','se_logC_GAP'],0.10,rtol=0.02),'physical eps=1 both pools pin F6P/GAP to ~0.10')
# ILR calibration should no longer contain artificial near-zero tail.
cal=pd.read_csv(SYN/'ilr_raw_sd_calibration_regularized.csv')
check(cal['regularized_delta_method_sd'].min()>5e-6,'regularized ILR calibration removes artificial near-zero tail')
# FreeFlux local covariance / MC values.
loc=pd.read_csv(FF/'toy_compositional_local_uncertainty.csv').set_index('design')
check(close(loc.loc['baseline','se_q2'],10.663,rtol=0.01),'FreeFlux baseline local SE(q2)')
check(close(loc.loc['fum','se_q2'],0.753,rtol=0.02),'FreeFlux +Fum local SE(q2)')
check(close(loc.loc['asp','se_q2'],0.168,rtol=0.02),'FreeFlux +Fum+Asp local SE(q2)')
mc=pd.read_csv(FF/'toy_compositional_monte_carlo.csv')
for d,target in [('baseline',10.486),('fum',0.753),('asp',0.165)]:
    x=float(mc[(mc.design==d)&(mc.direction=='q2')].rmse.iloc[0]); check(close(x,target,rtol=0.03),f'FreeFlux linearized MC {d} q2')
# Refined fixed-realization profile verification.
pr=pd.read_csv(FF/'freeflux_profile_refinement_summary.csv')
for d,a in [('fum',-1.280247),('asp',-0.403851),('asp',0.212244)]:
    row=pr[(pr.design==d)&(np.isclose(pr.a,a))].iloc[0]
    check(abs(float(row.delta_chi2_best)-3.84)<0.08,f'{d} profile threshold at a={a}')
    check(float(row.two_start_abs_chi2_difference)<1e-5,f'{d} two-start agreement at a={a}')
check(float(pr[(pr.design=='fum')&(pr.a==2.0)].delta_chi2_best.iloc[0])<0.6,'Fum has no upper crossing through a=2')
check(float(pr[(pr.design=='baseline')&(pr.a==-2.0)].delta_chi2_best.iloc[0])<0.5 and float(pr[(pr.design=='baseline')&(pr.a==2.0)].delta_chi2_best.iloc[0])<0.5,'baseline has no crossing on [-2,2]')
# M5 numerical stability.
fd=pd.read_csv(FF/'toy_fd_step_differences.csv')
check(float(fd['rel_Frobenius_Jh_minus_Jh2'].max())<5e-5,'toy centered-difference refinement stable')
ex=pd.read_csv(FF/'freeflux_exact_fd_step_differences.csv')
check(float(ex['rel_Frobenius_Jh_minus_Jh2'].max())<1e-7,'exact FreeFlux FDM refinement below runtime discretization error')
# Synthetic validation generator.
sv=pd.read_csv(ROOT/'validation/synthetic_validation_v101/asymptotic_identifiability_validation.csv')
check((sv.status=='PASS').all(),'synthetic validation cases all PASS')
status={'n_checks':None,'n_fail':len(issues),'issues':issues}
(ROOT/'audit/v101_claim_check.json').write_text(json.dumps(status,indent=2)+'\n')
if issues:
    print('V15_CLAIM_AUDIT=FAIL'); raise SystemExit(1)
print('V15_CLAIM_AUDIT=PASS')
