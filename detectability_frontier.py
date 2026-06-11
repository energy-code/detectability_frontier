#!/usr/bin/env python
"""
The Detectability Frontier of Low-Dose Radiation Epidemiology
=============================================================
Single-file analysis reproducing every numerical result, table, and figure in
"The Detectability Frontier of Low-Dose Radiation Epidemiology: A Source-Agnostic
Floor in the Linear-ERR Framework" (A. Stein).

Run:
    python detectability_frontier.py

Outputs (to the working directory):
    - console tables: INWORKS calibration (Table 1), validity frontier (Table 2),
      minimum detectable dose / possibility floor (Table 3), bounded consequences
      (Table 4), parameter sensitivity (Table 5), and the Section 2.5 / 3.1 checks;
    - figures: detectability_frontier.png, possibility_floor.png,
      inworks_calibration.png.

Requires: numpy, scipy, matplotlib. Everything is closed-form; no data files.

Method
------
Under the linear excess-relative-rate (ERR) stratified Poisson model
    rate = lambda0 * (1 + beta * Z),
the low-dose Fisher information for beta is I(beta) ~= D * <Z^2> (expected deaths
times the death-weighted mean-square dose), so the standard error is the
Cramer-Rao bound SE(beta) = 1 / sqrt(D * <Z^2>). A study is valid (significant at
one-sided level alpha and powered at 1-eta) iff MDE = (z_alpha + z_pow)*SE <= beta,
i.e. on the boundary

    (z_alpha + z_pow)^2 = N * p * Zrms^2 * beta^2.

INWORKS (Richardson et al., BMJ 2023) is used only to confirm that this SE
reproduces its published confidence intervals (calibration); beta is free
throughout and no INWORKS conclusion is adopted.
"""
import numpy as np
from scipy import stats

# --- conventions ------------------------------------------------------------
ALPHA_1SIDED, POWER = 0.05, 0.80
K = stats.norm.ppf(1 - ALPHA_1SIDED) + stats.norm.ppf(POWER)   # ~= 2.487
P_BASE = 0.20              # baseline lifetime solid-cancer mortality
LIFETIME_YEARS = 70
N_US, N_WORLD = 3.35e8, 8.0e9
BETA_LNT = 0.5            # nominal LNT-scale ERR per Gy (reference magnitude only)
Z95 = stats.norm.ppf(0.95)

# INWORKS restricted-range results (Richardson et al., BMJ 2023, suppl. Table C):
# (label, deaths D, ERR/Gy, lower 90% CL, upper 90% CL)
INWORKS_ROWS = [
    ("Full",     28089, 0.52,  0.27, 0.77),
    ("<400 mGy", 27960, 0.63,  0.34, 0.92),
    ("<200 mGy", 27429, 0.97,  0.55, 1.39),
    ("<100 mGy", 26283, 1.12,  0.45, 1.80),
    ("<50 mGy",  24518, 1.38,  0.20, 2.60),
    ("<20 mGy",  21293, 1.30, -1.33, 4.06),
]

# --- closed forms (Zrms in Gy unless noted) ---------------------------------
def n_required(zrms, beta):   return (K / (zrms * beta)) ** 2 / P_BASE
def min_dose_mGy(N, beta):    return 1e3 * K / (beta * np.sqrt(N * P_BASE))
def max_hidden_rate(N):       return K / np.sqrt(N * P_BASE)
def K_of(alpha, power):       return stats.norm.ppf(1 - alpha) + stats.norm.ppf(power)
def floor_mGy(beta, N, p, Kv):return 1e3 * Kv / (beta * np.sqrt(N * p))

def se_from_ci(lo, hi):       return (hi - lo) / (2 * Z95)
def zrms_from_ci(D, lo, hi):  return 1.0 / (se_from_ci(lo, hi) * np.sqrt(D))   # Gy


# === Table 1: INWORKS calibration ===========================================
def table1_calibration():
    print("TABLE 1  INWORKS calibration (recover Zrms from published 90% CIs)")
    print(f"{'Range':<10}{'Deaths':>8}{'ERR':>6}{'SE':>8}{'Zrms(mGy)':>11}{'%deaths':>9}")
    full = INWORKS_ROWS[0][1]
    for label, D, err, lo, hi in INWORKS_ROWS:
        se = se_from_ci(lo, hi)
        zrms = zrms_from_ci(D, lo, hi) * 1e3
        print(f"{label:<10}{D:>8d}{err:>6.2f}{se:>8.3f}{zrms:>11.1f}{100*D/full:>8.1f}%")
    print()


# === Table 2: validity frontier (people required) ==========================
def table2_frontier():
    doses = [0.07, 0.2, 0.7, 3.5, 20.0, 100.0]
    betas = [0.1, 0.5, 1.0, 5.0]
    print("TABLE 2  Validity frontier: people required (rows = mGy, cols = beta/Gy)")
    print(f"{'dose\\beta':<10}" + "".join(f"{b:>12g}" for b in betas))
    for d in doses:
        print(f"{d:<10g}" + "".join(f"{n_required(d*1e-3, b):>12.2g}" for b in betas))
    print()


# === Table 3: minimum detectable RMS dose by effect x study size ===========
def table3_floor():
    betas = [0.1, 0.5, 1.0, 5.0]
    sizes = [("1e5", 1e5), ("1e6", 1e6), ("US", N_US), ("World", N_WORLD)]
    print("TABLE 3  Min detectable RMS dose (mGy): rows = beta/Gy, cols = study size")
    print(f"{'beta':<8}" + "".join(f"{lab:>12}" for lab, _ in sizes))
    for b in betas:
        print(f"{b:<8g}" + "".join(f"{min_dose_mGy(N, b):>12.3g}" for _, N in sizes))
    print()


# === Table 4: bounded consequences =========================================
def table4_consequences():
    sizes = [1e5, 3e5, 1e6, 3e7, N_US, N_WORLD]
    print("TABLE 4  Bounded consequences (largest undetectable effect by study size)")
    print(f"{'N':>12}{'maxRate':>10}{'deaths/M/yr':>14}{'ind.risk':>12}")
    for N in sizes:
        r = max_hidden_rate(N)
        deaths_per_M_yr = r * (P_BASE / LIFETIME_YEARS) * 1e6
        ind_risk = r * P_BASE
        print(f"{N:>12.2g}{100*r:>9.3g}%{deaths_per_M_yr:>14.3g}{ind_risk:>12.2g}")
    print()


# === Table 5: sensitivity of the US-ceiling LNT floor to p and K ===========
def table5_sensitivity():
    base = floor_mGy(BETA_LNT, N_US, 0.20, K)
    print(f"TABLE 5  US-ceiling LNT-scale floor sensitivity (base {base:.3f} mGy)")
    rows = [
        ("p = 0.10",                0.10, K),
        ("p = 0.20 (base)",         0.20, K),
        ("p = 0.30",                0.30, K),
        ("power 0.70",              0.20, K_of(0.05, 0.70)),
        ("two-sided 0.05",          0.20, K_of(0.025, 0.80)),
        ("power 0.95",              0.20, K_of(0.05, 0.95)),
        ("alpha 0.01, power 0.95",  0.20, K_of(0.01, 0.95)),
        ("joint worst case",        0.10, K_of(0.01, 0.95)),
    ]
    for lab, p, Kv in rows:
        f = floor_mGy(BETA_LNT, N_US, p, Kv)
        print(f"  {lab:<24}{f:>8.2f} mGy   x{f/base:>5.2f}")
    print()


# === Section 3.1: significance != power (INWORKS <50 mGy stratum) ==========
def section31_stratum():
    D, ERR, lo, hi = 24518, 1.38, 0.20, 2.60
    se = se_from_ci(lo, hi)
    mde = K * se
    power_at_lnt = stats.norm.cdf(BETA_LNT / se - Z95)
    print("SECTION 3.1  INWORKS <50 mGy: statistically significant yet underpowered")
    print(f"  ERR={ERR}/Gy (90% CI {lo}-{hi}); SE={se:.3f}/Gy; "
          f"MDE(80% power)={mde:.2f}/Gy ({mde/BETA_LNT:.1f}x LNT slope); "
          f"power at beta=0.5: {power_at_lnt:.2f}")
    print()


# === Section 2.5: MDE = published upper CL at <20 mGy ======================
def section25_check():
    D, lo, hi = 21293, -1.33, 4.06
    se = se_from_ci(lo, hi)
    mde = K * se
    print(f"SECTION 2.5  <20 mGy: SE={se:.3f}/Gy, MDE=K*SE={mde:.2f}/Gy "
          f"(= published upper 90% CL {hi}/Gy; ~{mde/BETA_LNT:.0f}x LNT slope)")
    print()


# === Figures ===============================================================
def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter
    betas = [0.1, 0.5, 1.0, 5.0]

    # Fig 1: validity frontier, dose vs N
    doses = np.logspace(np.log10(0.02), np.log10(200), 400) * 1e-3   # Gy
    fig, ax = plt.subplots(figsize=(8, 6))
    for b in betas:
        ax.plot(doses * 1e3, n_required(doses, b), lw=2, label=f"beta={b}/Gy")
    ax.axhline(N_US, ls="--", color="grey", lw=1.2)
    ax.axhline(N_WORLD, ls="--", color="black", lw=1.2)
    ax.text(180, N_US * 1.4, "US population (3.35e8)", fontsize=8, color="grey", ha="right")
    ax.text(180, N_WORLD * 1.4, "world population (8e9)", fontsize=8, color="black", ha="right")
    ax.axhspan(1e5, 1e6, color="tab:orange", alpha=0.12)
    ax.text(0.025, 3e5, "typical 10-mile NPP cohort (1e5-1e6)", fontsize=8, color="tab:orange")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("RMS cumulative dose (mGy)")
    ax.set_ylabel("People required for a valid result")
    ax.set_title("Validity frontier: dose vs N, by effect size\n"
                 f"(below/left of each line = not recoverable; "
                 f"alpha={ALPHA_1SIDED} 1-sided, power={POWER})")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=20))
        axis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10)), numticks=80))
        axis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", ls="-", lw=0.7, alpha=0.55)
    ax.grid(which="minor", ls=":", lw=0.4, alpha=0.35)
    ax.set_ylim(1e2, 1e13)
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout(); fig.savefig("detectability_frontier.png", dpi=140); plt.close(fig)

    # Fig 2: possibility floor vs effect size
    bb = np.logspace(np.log10(0.05), np.log10(10.0), 400)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(bb, floor_mGy(bb, N_US, P_BASE, K), lw=2, color="tab:blue",
            label="US population ceiling (3.35e8)")
    ax.plot(bb, floor_mGy(bb, N_WORLD, P_BASE, K), lw=2, color="tab:red",
            label="world population ceiling (8e9)")
    for Nc, col in [(N_US, "tab:blue"), (N_WORLD, "tab:red")]:
        f = floor_mGy(BETA_LNT, Nc, P_BASE, K)
        ax.plot(0.5, f, "o", color=col, ms=7, zorder=5)
        ax.annotate(f"{f:.3g} mGy", (0.5, f), textcoords="offset points",
                    xytext=(8, 6), fontsize=9, color=col)
    ax.axvline(0.5, ls=":", color="grey", lw=1.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("True effect size  beta  (ERR per Gy)")
    ax.set_ylabel("Possibility floor: RMS cumulative dose (mGy)")
    ax.set_title("Possibility floor: dose below which detection is impossible")
    ax.grid(which="both", ls=":", lw=0.4, alpha=0.4)
    ax.legend(fontsize=9); fig.tight_layout()
    fig.savefig("possibility_floor.png", dpi=140); plt.close(fig)

    # Fig 3: INWORKS calibration collapse
    labels = [r[0] for r in INWORKS_ROWS]
    zrms = [zrms_from_ci(D, lo, hi) * 1e3 for _, D, _, lo, hi in INWORKS_ROWS]
    half = [(hi - lo) / 2.0 for _, _, _, lo, hi in INWORKS_ROWS]
    x = np.arange(len(labels))
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax1.plot(x, zrms, "o-", color="tab:blue", lw=2, ms=7, label="recovered Zrms (mGy)")
    ax1.set_xticks(x); ax1.set_xticklabels(labels)
    ax1.set_xlabel("Restricted dose range (progressively narrower)")
    ax1.set_ylabel("recovered Zrms (mGy)", color="tab:blue"); ax1.set_ylim(0, 45)
    ax2 = ax1.twinx()
    ax2.plot(x, half, "s--", color="tab:red", lw=2, ms=7, label="90% CI half-width (/Gy)")
    ax2.set_ylabel("90% CI half-width (ERR per Gy)", color="tab:red"); ax2.set_ylim(0, 3.0)
    ax1.set_title("INWORKS calibration: information collapses with dose range")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], fontsize=9, loc="upper center")
    fig.tight_layout(); fig.savefig("inworks_calibration.png", dpi=140); plt.close(fig)


def main():
    print(f"K = z(0.95) + z({POWER}) = {K:.3f};  p = {P_BASE};  "
          f"ceilings: US = {N_US:g}, world = {N_WORLD:g}\n")
    table1_calibration()
    table2_frontier()
    table3_floor()
    table4_consequences()
    table5_sensitivity()
    section31_stratum()
    section25_check()
    figures()
    print("figures written: detectability_frontier.png, possibility_floor.png, "
          "inworks_calibration.png")


if __name__ == "__main__":
    main()
