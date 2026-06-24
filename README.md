# Detectability Frontier — Low-Dose Radiation Epidemiology

Code for **"The Detectability Frontier: An Information Floor for Low-Dose
Radiation Epidemiology"** (A. Stein).

The analysis reproduces every numerical result, table, and figure in the paper
from closed-form expressions; there are no input data files. Three scripts:

- `detectability_frontier.py` (Python) — prints all tables and the in-text
  checks and writes **Figure 2** plus the auxiliary figures.
- `floor_figure.R` (R) — draws **Figure 1**, the possibility-floor heatmap.
- `graphical_abstract.py` (Python) — draws the journal graphical abstract.

## Run

```bash
pip install numpy scipy matplotlib
python detectability_frontier.py
python graphical_abstract.py
```

```r
install.packages(c("ggplot2", "geomtextpath", "scales"))
Rscript floor_figure.R
```

`detectability_frontier.py` prints:

- **Table 1** — INWORKS calibration (recovered RMS dose from published 90% CIs)
- **Table 2** — validity frontier (people required by dose and effect)
- **Table 3** — minimum detectable RMS dose by effect size and study size
  (the possibility floor at the US and world population ceilings)
- **Table 4** — bounded consequences (largest undetectable effect by study size)
- **Table 5** — sensitivity of the floor to baseline mortality `p` and test
  stringency `K = z_alpha + z_pow`
- two in-text checks: the MDE-equals-published-upper-CL check (<20 mGy) and the
  significance-without-power check (<50 mGy)

and writes:

- `detectability_frontier.png` — **Figure 2**, the validity frontier
  (people required vs dose, one curve per effect size)
- `possibility_floor.png` — auxiliary (floor vs effect size at the US/world
  ceilings; superseded in the manuscript by Figure 1)
- `inworks_calibration.png` — auxiliary visualization of the Table 1 INWORKS
  collapse (not a numbered figure in the manuscript)

`floor_figure.R` writes:

- `floor_heatmap.png` — **Figure 1**, the possibility floor as a field
  (people required across dose and effect size, with the population floors marked)

`graphical_abstract.py` writes:

- `graphical_abstract.png` — the journal graphical abstract (not a numbered
  manuscript figure)

## Method

Under the linear excess-relative-rate (ERR) stratified Poisson model
`rate = λ₀·(1 + βZ)`, the low-dose Fisher information for the slope is
`I(β) ≈ D·⟨Z²⟩` (expected deaths × death-weighted mean-square dose), so the
standard error is the Cramér–Rao bound `SE(β) = 1/√(D·⟨Z²⟩)`. A study is valid
(significant at one-sided `α`, powered at `1−η`) iff `MDE = (z_α+z_pow)·SE ≤ β`,
giving the closed-form frontier `(z_α+z_pow)² = N·p·Z_rms²·β²`.

## Data

No new data are used. The calibration relies on published 90% confidence
intervals from the International Nuclear Workers Study (Richardson et al., *BMJ*
2023;382:e074520, Supplementary Table C), hard-coded in the script; `β` is a
free parameter throughout and no INWORKS conclusion is adopted.

## License

MIT — see [LICENSE.md](LICENSE.md).
