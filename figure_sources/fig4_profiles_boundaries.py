from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


profiles = read_csv(DATA / "freeflux_profile_likelihood.csv")
profile_m1 = [row for row in profiles if row.get("direction") == "m1"]
profile_m2 = [row for row in profiles if row.get("direction") == "m2"]

m2_rows = read_csv(DATA / "freeflux_m2_branch_continuation.csv")
m2_rows = [row for row in m2_rows if row.get("branch") in {"plus", "minus"}]

m1_coordinate, m1_objective = [], []
with (DATA / "freeflux_m1_continuation_trace.csv").open(encoding="utf-8") as stream:
    for line in stream:
        line = line.strip()
        if not line or line.startswith("far_branch"):
            break
        if line.startswith("a,"):
            continue
        coordinate, objective = line.split(",")
        m1_coordinate.append(float(coordinate))
        m1_objective.append(float(objective))

fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.1), constrained_layout=True)

ax = axes[0]
for label, rows, marker in [
    ("m=1", profile_m1, "o"),
    ("m=2", profile_m2, "s"),
]:
    x = np.array([float(row["profile_coordinate"]) for row in rows])
    y = np.array([float(row["delta_chi2"]) for row in rows])
    order = np.argsort(x)
    ax.semilogy(x[order], np.maximum(y[order], 1e-12), marker=marker, label=label)
ax.axhline(3.84, linestyle="--", label=r"$\Delta\chi^2=3.84$")
ax.set_xlabel("Profile coordinate")
ax.set_ylabel(r"$\Delta\chi^2$")
ax.set_title("(a) Nonlinear profiles", pad=8)
ax.legend(fontsize=8)

ax = axes[1]
for branch, marker in [("plus", "o"), ("minus", "s")]:
    rows = [row for row in m2_rows if row["branch"] == branch]
    coordinate = np.array([float(row["a"]) for row in rows])
    oaa = np.array([float(row["C_OAA_multiplier"]) for row in rows])
    fum = np.array([float(row["C_Fum_multiplier"]) for row in rows])
    ax.plot(coordinate, oaa, marker=marker, label=f"{branch}: OAA")
    ax.plot(coordinate, fum, marker=marker, linestyle="--", label=f"{branch}: Fum")
ax.set_xlabel(r"Continuation coordinate $a$")
ax.set_ylabel("Pool-size multiplier")
ax.set_title(r"(b) Two ends of the $m=2$ valley", pad=8)
ax.legend(fontsize=7)

ax = axes[2]
ax.semilogy(
    m1_coordinate,
    np.maximum(m1_objective, 1e-12),
    marker="o",
    label=r"$m=1$ branch",
)
ax.axhline(
    5.369302360156972e-6,
    linestyle="--",
    label=r"$C_{\mathrm{Fum}}=0$ boundary",
)
ax.set_xlabel(r"Continuation coordinate $a$")
ax.set_ylabel(r"$\Delta\chi^2$")
ax.set_title("(c) Approach to reduced boundary", pad=8)
ax.legend(fontsize=8)

pdf = OUT / "fig4_profiles_boundaries.pdf"
png = OUT / "fig4_profiles_boundaries.png"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=240, bbox_inches="tight")
plt.close(fig)
print(pdf)
print(png)
