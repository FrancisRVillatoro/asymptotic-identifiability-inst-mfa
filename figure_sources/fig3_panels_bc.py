from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / 'figures'
DATADIR = ROOT / 'data'

with (DATADIR / 'freeflux_spectrum_extended.csv').open(newline='', encoding='utf-8') as f:
    rows = list(csv.DictReader(f))

eps = np.array([float(r['epsilon']) for r in rows])
lam = np.array([[float(r[f'lambda_desc_{j}']) for j in range(1, 10)] for r in rows])
eps_pub = 2 / 175

# Panel (b)
fig, ax = plt.subplots(figsize=(4.4, 4.2))
names = ['OAA', 'Cit', 'AKG', 'Suc', 'Fum', 'Glu']
values = [175.0, 2.0, 100.0 / 3.0, 5.0, 62.5, 10.0]
x = np.arange(len(names))
ax.bar(x, values)
ax.set_yscale('log')
ax.set_xticks(x, names, rotation=35, ha='right', fontsize=13)
ax.tick_params(axis='y', labelsize=13)
ax.set_ylabel('Turnover rate', fontsize=16)
ax.set_title('(b) Published turnover hierarchy', fontsize=20, fontweight='bold', pad=10)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig3_panel_b.pdf', bbox_inches='tight')
fig.savefig(FIGDIR / 'fig3_panel_b.png', dpi=240, bbox_inches='tight')
plt.close(fig)

# Panel (c)
fig, ax = plt.subplots(figsize=(5.2, 4.2))
for j in range(7):
    ax.loglog(eps, lam[:, j], marker='o', markersize=4,
              label=rf'$\lambda_{{{j+1}}}$')
ax.axvline(eps_pub, linestyle='--',
           label=r'$\varepsilon_{\mathrm{pub}}=2/175$')
ax.set_xlim(eps.min() * 0.9, eps.max() * 1.08)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel(r'$\varepsilon$', fontsize=16)
ax.set_ylabel('Information eigenvalue', fontsize=16)
ax.set_title('(c) Joint flux--pool spectrum', fontsize=20, fontweight='bold', pad=10)
ax.legend(fontsize=10, ncol=2)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig3_panel_c.pdf', bbox_inches='tight')
fig.savefig(FIGDIR / 'fig3_panel_c.png', dpi=240, bbox_inches='tight')
plt.close(fig)
