#!/usr/bin/env python3
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures_v101'; OUT.mkdir(exist_ok=True)
SV=ROOT/'validation/synthetic_validation_v101'; FF=ROOT/'validation/freeflux_v101'; SYN=ROOT/'validation/synechocystis_v101/results'
# Synthetic validation panel A + Guppy panel B.
val=pd.read_csv(SV/'asymptotic_identifiability_validation.csv')
# Parse semicolon-separated exponents; use model rows only.
pred=[];obs=[]
for _,r in val.iloc[:-1].iterrows():
    ps=[x.strip() for x in str(r.predicted_fim_exponents).replace(';',',').split(',') if x.strip() not in ('','inf','∞','None','structural')]
    os=[x.strip() for x in str(r.observed_fim_exponents).replace(';',',').split(',') if x.strip() not in ('','inf','∞','None','structural')]
    for a,b in zip(ps,os):
        try: pred.append(float(a));obs.append(float(b))
        except: pass
fig,ax=plt.subplots(figsize=(7.2,5.4)); ax.scatter(pred,obs); lo=min(pred+obs)-.2;hi=max(pred+obs)+.2;ax.plot([lo,hi],[lo,hi]);ax.set_xlabel('Predicted FIM exponent');ax.set_ylabel('Measured log-log exponent');fig.tight_layout();fig.savefig(OUT/'fig2a_validation.pdf');plt.close(fig)
# Guppy information eigenvalues from direct singular values.
import sys
sys.path.insert(0, str(SV))
spec=importlib.util.spec_from_file_location('validate',SV/'validate_asymptotic_identifiability.py'); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
eps=np.logspace(-3,-5,13); curves=[]
for e in eps:
    s=np.linalg.svd(mod.guppy_matrix(float(e)),compute_uv=False); curves.append(s*s)
curves=np.asarray(curves)
fig,ax=plt.subplots(figsize=(7.2,5.4))
for j in range(curves.shape[1]): ax.loglog(eps,curves[:,j],marker='o',markersize=2)
ax.set_xlabel(r'$\varepsilon$');ax.set_ylabel('Information eigenvalue');fig.tight_layout();fig.savefig(OUT/'fig2b_guppy.pdf');plt.close(fig)
# FreeFlux fixed-realization nonlinear profile.
# The combined paper figure uses the same audited profile data in two panels:
# a full nonlinear view and a threshold-scale view. Bound-limited exploratory
# points are shown only when they lie inside the displayed vertical range and
# are never connected to the profile curves.
prof=pd.read_csv(FF/'toy_compositional_noisy_profiles.csv')
refine=pd.read_csv(FF/'freeflux_profile_refinement_summary.csv')
prof=prof[prof.design.isin(['baseline','fum','asp'])].copy()
refine=refine[refine.design.isin(['baseline','fum','asp'])].copy()
labels={'baseline':'baseline','fum':'+ Fum slow','asp':'+ Fum + Asp-fast'}
markers={'baseline':'o','fum':'s','asp':'^'}
colors={'baseline':'#1f77b4','fum':'#d95f02','asp':'#2ca02c'}
plot_order=['fum','asp','baseline']
zorders={'fum':3,'asp':4,'baseline':6}
root_pts=refine[refine.delta_chi2_best.sub(3.84).abs()<0.02]
fig,axes=plt.subplots(1,2,figsize=(10.5,4.6),constrained_layout=True)
for ax,ylim,panel,title in zip(
    axes,[(0,60),(0,8)],['(a)','(b)'],
    ['Fully nonlinear profiles','Threshold-scale view']):
    ymax=ylim[1]
    for d in plot_order:
        q=prof[prof.design==d].sort_values('a_q2')
        free=q[~q.nuisance_hit_bound.astype(bool)]
        bound=q[q.nuisance_hit_bound.astype(bool)]
        ax.plot(free.a_q2,free.delta_chi2,
                marker=markers[d],markersize=5.4 if d=='baseline' else 5.0,
                linewidth=2.2 if d=='baseline' else 1.8,
                color=colors[d],markerfacecolor=colors[d],
                markeredgecolor='white' if d=='baseline' else colors[d],
                markeredgewidth=0.6 if d=='baseline' else 0.4,
                label=labels[d],zorder=zorders[d])
        bound_visible=bound[bound.delta_chi2<ymax]
        if len(bound_visible):
            ax.scatter(bound_visible.a_q2,bound_visible.delta_chi2,
                       marker=markers[d],s=52 if d=='baseline' else 46,
                       facecolors='none',edgecolors=colors[d],linewidths=1.5,
                       zorder=zorders[d]+0.2)
    ax.axhline(3.84,color='0.35',linestyle='--',linewidth=1.2,zorder=1)
    if panel=='(b)':
        ax.text(1.94,3.84+0.11,r'$\Delta\chi^2 = 3.84$',ha='right',va='bottom',color='0.25',fontsize=9)
    for _,r in root_pts[root_pts.design!='baseline'].iterrows():
        ax.scatter([r['a']],[3.84],marker='x',s=70,linewidths=1.8,color='black',zorder=8)
    ax.set_xlim(-2.05,2.05); ax.set_ylim(*ylim)
    ax.set_xlabel(r'profile coordinate $a$'); ax.set_ylabel(r'$\Delta\chi^2$')
    ax.set_title(title); ax.text(0.02,0.98,panel,transform=ax.transAxes,ha='left',va='top',fontweight='bold')
    ax.tick_params(direction='out',length=4,width=0.8)
    ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]); ax.set_axisbelow(True)
handles,labs=axes[1].get_legend_handles_labels()
order=[labs.index('baseline'),labs.index('+ Fum slow'),labs.index('+ Fum + Asp-fast')]
axes[1].legend([handles[i] for i in order],[labs[i] for i in order],loc='upper right',frameon=False)
fig.savefig(OUT/'fig8_freeflux_profiles.pdf',bbox_inches='tight')
fig.savefig(OUT/'fig8_freeflux_profiles.png',dpi=240,bbox_inches='tight')
plt.close(fig)
# Preserve the historical single-panel outputs for backward compatibility.
fig,ax=plt.subplots(figsize=(7.2,5.4))
for d,label in [('baseline','baseline'),('fum','+ Fum slow'),('asp','+ Fum + Asp-fast')]:
    q=prof[(prof.design==d)&(~prof.nuisance_hit_bound.astype(bool))].sort_values('a_q2');ax.plot(q.a_q2,q.delta_chi2,marker='o',label=label)
ax.axhline(3.84,linestyle='--');ax.set_xlabel('profile coordinate a');ax.set_ylabel(r'$\Delta\chi^2$');ax.legend();fig.tight_layout();fig.savefig(OUT/'fig8a_freeflux_profile.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.2,5.4))
for d,label in [('baseline','baseline'),('fum','+ Fum slow'),('asp','+ Fum + Asp-fast')]:
    q=prof[(prof.design==d)&(~prof.nuisance_hit_bound.astype(bool))].sort_values('a_q2');ax.plot(q.a_q2,q.delta_chi2,marker='o',label=label)
ax.axhline(3.84,linestyle='--');ax.set_ylim(0,8);ax.set_xlabel('profile coordinate a');ax.set_ylabel(r'$\Delta\chi^2$');ax.legend();fig.tight_layout();fig.savefig(OUT/'fig8b_freeflux_profile_zoom.pdf');plt.close(fig)
# Regularized ILR raw-space calibration.
cal=pd.read_csv(SYN/'ilr_raw_sd_calibration_regularized.csv'); v=np.sort(cal.regularized_delta_method_sd.to_numpy())
fig,ax=plt.subplots(figsize=(7.2,5.4));ax.plot(np.arange(1,len(v)+1),v);ax.set_yscale('log');ax.axhline(.005,linestyle='--');ax.axhline(.028,linestyle=':');ax.set_xlabel('MID component rank');ax.set_ylabel('raw-space SD');fig.tight_layout();fig.savefig(OUT/'figS11a_ilr_calibration.pdf');plt.close(fig)
# Exact-J ILR-weighted scaling.
sw=pd.read_csv(SYN/'weighted_design_mode_sweep.csv'); ex=pd.read_csv(SYN/'weighted_design_mode_exponents.csv')
fig,ax=plt.subplots(figsize=(7.2,5.4))
for d,label in [('baseline','baseline weak 1'),('baseline','baseline weak 2'),('F6P','+F6P remaining weak'),('both','+F6P,+GAP weakest m=0')]:
    ee=ex[ex.design==d].sort_values('sigma_deep')
    if d=='baseline': mode=int(ee.iloc[0 if '1' in label else 1].mode_id)
    else: mode=int(ee.iloc[0].mode_id)
    q=sw[(sw.design==d)&(sw.mode_id==mode)].sort_values('epsilon');ax.loglog(q.epsilon,q.sigma,marker='o',label=label)
ax.set_xlabel(r'$\varepsilon$');ax.set_ylabel('whitened sensitivity singular value');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'figS11c_synech_scaling.pdf');plt.close(fig)
print('V101_KEY_FIGURES_COMPLETED')
