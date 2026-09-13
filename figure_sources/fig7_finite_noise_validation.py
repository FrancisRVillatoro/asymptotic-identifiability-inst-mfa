from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
data = root / 'validation' / 'block4_finite_noise'
out = root / 'figures'

cal = pd.read_csv(data / 'syn_ilr_raw_sd_calibration.csv')
toyprof = pd.read_csv(data / 'toy_compositional_noisy_gn_profiles.csv')
synsweep = pd.read_csv(data / 'syn_compositional_local_sweep.csv')
synprof = pd.read_csv(data / 'syn_poolcontrast_noisy_gn_profiles_eps1_64.csv')

fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2))
ax = axes[0,0]
v=np.sort(cal['delta_method_sd'].to_numpy()); ax.plot(np.arange(1,len(v)+1),v)
ax.axhline(0.005,linestyle='--',label='published example SD = 0.005')
ax.axhline(0.028,linestyle=':',label='published example SD = 0.028')
ax.set_yscale('log'); ax.set_xlabel('MID component rank'); ax.set_ylabel('Raw-space SD (delta method)'); ax.set_title('(a) Logistic-normal noise calibration'); ax.legend(fontsize=8)

ax=axes[0,1]
for d,m,l in [('baseline','o','baseline'),('fum','s','+ Fum slow'),('asp','^','+ Fum + Asp-fast')]:
    rr=toyprof[toyprof.design==d].sort_values('a_q2'); ax.plot(rr.a_q2,rr.delta_chi2,marker=m,label=l)
ax.axhline(3.84,linestyle='--',label=r'$\Delta\chi^2=3.84$'); ax.set_xlabel(r'profile coordinate along baseline $m=2$ direction'); ax.set_ylabel(r'$\Delta\chi^2$'); ax.set_title('(b) Toy benchmark: noisy profile contraction'); ax.legend(fontsize=8)

ax=axes[1,0]
for d,idx,m,l in [('baseline',59,'o',r'baseline $\sigma_{59}$'),('baseline',58,'s',r'baseline $\sigma_{58}$'),('F6P',59,'^',r'$+C_{\rm F6P}$, weakest'),('both',59,'D',r'$+C_{\rm F6P},C_{\rm GAP}$, weakest')]:
    rr=synsweep[(synsweep.design==d)&(synsweep.index_desc==idx)].sort_values('epsilon'); ax.loglog(rr.epsilon,rr.sigma,marker=m,label=l)
ax.set_xlabel(r'$\varepsilon$'); ax.set_ylabel('Whitened sensitivity singular value'); ax.set_title('(c) Synechocystis: finite-noise scaling'); ax.legend(fontsize=8)

ax=axes[1,1]
for d,m,l in [('baseline','o','baseline'),('F6P','s',r'$+C_{\rm F6P}$'),('GAP','^',r'$+C_{\rm GAP}$'),('both','D',r'$+C_{\rm F6P},C_{\rm GAP}$')]:
    rr=synprof[synprof.design==d].sort_values('a_contrast'); ax.plot(rr.a_contrast,rr.delta_chi2,marker=m,label=l)
ax.axhline(3.84,linestyle='--',label=r'$\Delta\chi^2=3.84$'); ax.set_xlabel(r'fast-pool contrast $(\log C_{\rm F6P}-\log C_{\rm GAP})/\sqrt{2}$'); ax.set_ylabel(r'$\Delta\chi^2$'); ax.set_title(r'(d) Synechocystis at $\varepsilon=1/64$'); ax.legend(fontsize=7)

fig.tight_layout(); fig.savefig(out/'fig7_finite_noise_validation.pdf',bbox_inches='tight'); fig.savefig(out/'fig7_finite_noise_validation.png',dpi=220,bbox_inches='tight'); plt.close(fig)
