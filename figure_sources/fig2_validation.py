from pathlib import Path
import ast
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


validation = read_csv(DATA / "asymptotic_identifiability_validation.csv")
predicted, measured = [], []
for row in validation:
    try:
        pvals = ast.literal_eval(row["predicted_fim_exponents"])
        ovals = ast.literal_eval(row["observed_fim_exponents"])
    except (SyntaxError, ValueError):
        continue
    for pval, oval in zip(pvals, ovals):
        if pval is None or oval is None:
            continue
        predicted.append(float(pval))
        oval = float(oval)
        measured.append(0.0 if abs(oval) < 0.15 else oval)

guppy = read_csv(DATA / "phaseI_guppy_fim_spectrum.csv")
epsilon = np.array([float(row["epsilon"]) for row in guppy])
eigenvalues = np.array(
    [[float(row[f"lambda{i}"]) for i in range(1, 6)] for row in guppy]
)

fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)

ax = axes[0]
ax.scatter(predicted, measured, s=28)
lo = min(predicted + measured) - 0.15
hi = max(predicted + measured) + 0.15
ax.plot([lo, hi], [lo, hi], linestyle="--")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("Predicted FIM exponent")
ax.set_ylabel("Measured log-log exponent")
ax.set_title("(a) All validation cases", pad=8)

ax = axes[1]
for index in range(5):
    ax.loglog(
        epsilon,
        eigenvalues[:, index],
        marker="o",
        markersize=3,
        label=rf"$\lambda_{{{index + 1}}}$",
    )
ax.set_xlabel(r"$\varepsilon$")
ax.set_ylabel("FIM eigenvalue")
ax.set_title("(b) Reversible KFP benchmark", pad=8)
ax.legend(fontsize=8, ncol=2)

pdf = OUT / "fig2_validation.pdf"
png = OUT / "fig2_validation.png"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=240, bbox_inches="tight")
plt.close(fig)
print(pdf)
print(png)
