from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"
EPS_PUB = 2.0 / 175.0


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


design_rows = read_csv(DATA / "freeflux_design_lifting.csv")
close_rows = read_csv(DATA / "freeflux_close_AKG_mode.csv")
full9_rows = read_csv(DATA / "full9_final_scaling.csv")

published_points = {
    "baseline": 2.46959157e-10,
    "baseline+OAA": 1.31180578e-8,
    "baseline+Fum": 1.20117421e-6,
    "baseline+Asp-fast": 5.522284e-5,
}

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.45), constrained_layout=True)

# Panel (a)
ax = axes[0]
for design, marker, label in [
    ("baseline", "o", "baseline"),
    ("baseline+OAA", "s", "+ OAA"),
    ("baseline+Fum", "^", "+ Fum"),
    ("baseline+Asp-fast", "D", "+ Asp-fast"),
]:
    rows = [row for row in design_rows if row["design"] == design]
    epsilon = np.array([float(row["epsilon"]) for row in rows])
    column = "lambda_desc_7" if design in {"baseline", "baseline+OAA"} else "lambda_desc_8"
    values = np.array([float(row[column]) for row in rows])
    epsilon = np.append(epsilon, EPS_PUB)
    values = np.append(values, published_points[design])
    order = np.argsort(epsilon)
    ax.loglog(epsilon[order], values[order], marker=marker, markersize=4, label=label)

ax.axvline(EPS_PUB, linestyle="--", label=r"$\varepsilon_{\mathrm{pub}}$")
ax.set_xlabel(r"$\varepsilon$")
ax.set_ylabel("Weakest nonstructural eigenvalue")
ax.set_title("(a) Observable/tracer lifting", pad=8)
# Low-right, inside the axes and clear of the curves.
ax.legend(fontsize=7, loc="lower right", bbox_to_anchor=(0.98, 0.03))

# Panel (b)
ax = axes[1]
for design, marker, label in [
    ("Asp-fast", "o", "scale fixed: Asp-fast"),
    ("Asp-fast+pool-AKG", "s", r"scale fixed: + $C_{\mathrm{AKG}}$"),
]:
    rows = [row for row in close_rows if row["design"] == design]
    epsilon = np.array([float(row["epsilon"]) for row in rows])
    values = np.array([float(row["lambda_desc_8"]) for row in rows])
    order = np.argsort(epsilon)
    ax.loglog(epsilon[order], values[order], marker=marker, markersize=4, label=label)

for design, marker, label in [
    ("Asp-fast+CAKG", "^", r"full 9D: + $C_{\mathrm{AKG}}$"),
    ("Asp-fast+CAKG+CCit", "D", r"full 9D: + $C_{\mathrm{AKG}},C_{\mathrm{Cit}}$"),
]:
    rows = [row for row in full9_rows if row["design"] == design]
    epsilon = np.array([float(row["epsilon"]) for row in rows])
    values = np.array([float(row["lambda_min"]) for row in rows])
    order = np.argsort(epsilon)
    ax.loglog(epsilon[order], values[order], marker=marker, markersize=4, label=label)

ax.axvline(EPS_PUB, linestyle="--")
ax.set_xlabel(r"$\varepsilon$")
ax.set_ylabel("Smallest eigenvalue")
ax.set_title("(b) Global-scale accounting", pad=8)
# Mid-left near 1e-4, as requested, without covering the blue/green trends.
ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(0.035, 0.58))

pdf = OUT / "fig5_design_lifting.pdf"
png = OUT / "fig5_design_lifting.png"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=240, bbox_inches="tight")
plt.close(fig)
print(pdf)
print(png)
