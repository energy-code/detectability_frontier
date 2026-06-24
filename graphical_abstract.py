"""
graphical_abstract.py  --  Entropy graphical abstract (standalone; not a body figure)

An original conceptual schematic for the GA slot: the minimum detectable RMS dose
(the "possibility floor") as study size N grows from a small cohort to the entire
world population. The floor flattens onto a hard limit that routine nuclear-plant
public doses fall at or below — i.e. more people barely lowers it, and below the
world floor no study can resolve the effect at all.

Meets MDPI GA spec: >= 1100 x 560 px, PNG, Arial, decimal points, no "Graphical
Abstract" heading in the image, not a copy/superposition of the body figures.

Run:  py graphical_abstract.py   ->   graphical_abstract.png
"""
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Arial if available (MDPI-recommended), else a clean sans fallback
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

# constants (identical to detectability_frontier.py)
K = stats.norm.ppf(0.95) + stats.norm.ppf(0.80)   # ~2.486
P_BASE = 0.20
BETA_LNT = 0.5
N_US, N_WORLD = 3.35e8, 8e9

def zfloor_mGy(N, beta=BETA_LNT, p=P_BASE, k=K):
    return k / (beta * np.sqrt(N * p)) * 1e3

fig, ax = plt.subplots(figsize=(11.0, 5.8))   # -> 2200 x 1160 px at dpi=200

Nn = np.logspace(4, np.log10(N_WORLD), 500)   # curve stops at the world population
floor = zfloor_mGy(Nn)
world_floor = zfloor_mGy(N_WORLD)   # 0.124
us_floor = zfloor_mGy(N_US)         # 0.608

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e4, 9.5e9); ax.set_ylim(0.03, 200)

# (1) the absolute-impossibility zone: below the world floor, nothing resolves
ax.axhspan(0.03, world_floor, color="#d62728", alpha=0.13, zorder=0)
ax.text(1.6e4, 0.052, "Undetectable by ANY study\n— even the entire world population",
        fontsize=11, color="#7c1414", va="center", ha="left", fontweight="bold")

# (2) routine nuclear-plant public dose band (mostly trapped at/below the floor)
ax.axhspan(0.07, 3.5, color="#1f77b4", alpha=0.10, zorder=0)
ax.text(1e5, 0.5, "routine NPP public dose (0.07–3.5 mGy)", fontsize=9.5,
        color="#11507e", va="center", ha="center")

# (3) the floor itself
ax.plot(Nn, floor, color="black", lw=3.2, zorder=4,
        label="minimum detectable dose (LNT-scale effect)")

# (4) study-size markers along the floor
pts = [(1e5, "100k\ncohort"), (1e6, "1M"), (N_US, "US\n335M"), (N_WORLD, "world\n8B")]
for N, lab in pts:
    f = zfloor_mGy(N)
    ax.plot(N, f, "o", ms=9, color="black", zorder=5)
    ax.annotate(f"{f:.3g} mGy", (N, f), textcoords="offset points", xytext=(0, 11),
                ha="center", fontsize=10.5, fontweight="bold")
    ax.annotate(lab, (N, f), textcoords="offset points", xytext=(0, -26),
                ha="center", fontsize=9.5, color="#333333")

# (5) the floor: even at world-population scale the minimum detectable dose
#     bottoms out here; below it nothing is resolvable
ax.plot([1e4, N_WORLD], [world_floor, world_floor], ls="--", lw=1.4,
        color="#7c1414", zorder=3)
ax.text(1.6e4, world_floor * 1.28, "absolute floor ≈ 0.12 mGy", fontsize=10.5,
        color="#7c1414", fontweight="bold", ha="left", va="bottom")

ax.set_xlabel("study size  N  (people)", fontsize=12.5)
ax.set_ylabel("RMS cumulative dose (mGy)", fontsize=12.5)
ax.set_title("Low-dose radiation: a detectability floor no study can cross",
             fontsize=15.5, fontweight="bold", pad=12)
ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.4)
ax.tick_params(labelsize=10.5)

fig.tight_layout()
fig.savefig("graphical_abstract.png", dpi=200)
print("wrote graphical_abstract.png")
w, h = fig.get_size_inches() * fig.dpi
print(f"size: {int(w)} x {int(h)} px (need >= 1100 x 560)")
print(f"US floor {us_floor:.3f} mGy; world floor {world_floor:.3f} mGy")
