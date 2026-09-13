from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
turn=pd.read_csv(ROOT/'data/synechocystis_turnover_spectrum.csv')
out=ROOT/'figures'
fig, ax = plt.subplots(figsize=(8.6,4.6))
t=turn.sort_values('turnover').reset_index(drop=True)
x=np.arange(len(t))
ax.semilogy(x,t['turnover'],marker='o',markersize=3)
ax.set_xticks(x)
ax.set_xticklabels(t['metabolite'],rotation=90,fontsize=7)
ax.set_ylabel('Turnover rate')
ax.set_xlabel('Metabolite (in increasing turnover order)')
ax.set_title('Automatic fast-block candidates from turnover gaps')
for label in ['G3P','Mal','TA']:
    i=int(t.index[t.metabolite==label][0])
    ax.axvline(i+0.5,linestyle='--')
    ax.text(i+0.55,t['turnover'].max()/1.6,f'cut after {label}',rotation=90,va='top',fontsize=8)
fig.tight_layout()
fig.savefig(out/'fig6_synechocystis_turnover_candidates.pdf',bbox_inches='tight')
fig.savefig(out/'fig6_synechocystis_turnover_candidates.png',dpi=220,bbox_inches='tight')
plt.close(fig)
