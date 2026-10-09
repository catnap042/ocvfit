# ocvfit

[![CI](https://github.com/catnap042/ocvfit/actions/workflows/ci.yml/badge.svg)](https://github.com/catnap042/ocvfit/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Fit full-cell open-circuit voltage (OCV) curves from half-cell OCP curves and
quantify lithium-ion degradation modes (LLI, LAM_PE, LAM_NE).**

`ocvfit` is a small, dependency-light Python package (NumPy + SciPy; Matplotlib
optional) that implements the classic "electrode alignment" approach to
non-destructive battery diagnostics: a low-rate full-cell voltage curve is
reconstructed from the open-circuit potentials (OCPs) of the positive and
negative electrodes, and the fitted electrode capacities and lithium inventory
reveal how the cell has aged.

![Example: BOL vs. aged fit](docs/images/quickstart.png)

*Synthetic LFP/graphite cell built from published OCP functions (LFP: Afshar
et al. 2017; graphite: Chen et al. 2020, the pairing of the PyBaMM `Prada2013`
parameter set). Panels: fitted voltage, residuals, aligned electrode potentials
and differential voltage, for a begin-of-life and an aged curve. The graphite
staging peaks in dV/dQ (panel d) are what anchor the fit on the flat LFP
plateau. Produced by [`examples/quickstart.py`](examples/quickstart.py).*

## Features

- **Three model parameterisations**
  - `CathodeFixedModel` - 3 parameters (`Kn`, `Sn`, `V_offset`), the positive
    window is fixed and only the negative curve is scaled and shifted;
  - `WindowModel` - free electrode windows (`y0`, `y100`, `x0`, `x100`,
    `V_offset`), optionally with a blended positive electrode;
  - `CapacityModel` - absolute electrode capacities `Q_PE`, `Q_NE`, cyclable
    lithium `n_Li` and a polarisation offset `R`, optionally anchored to the
    rested full-charge voltage. Maps directly onto LLI / LAM_PE / LAM_NE.
- **Blended positive electrodes** (e.g. NMC + LFP) with one capacity per
  component, using iso-potential capacity addition.
- **Robust fitting**: Huber loss, edge and slope-aware weights, end-point
  anchoring and a valid-region penalty for SOC-axis models; grid pre-screening
  and a two-stage (`soft_l1` then least-squares) fit for the capacity model.
- **Degradation modes** and the feasible range of side-reaction lithium loss.
- **Identifiability checks**: correlation-corrected profile-likelihood
  intervals plus a reference-curve bias budget.
- **Half-cell curve utilities**: shape-preserving PCHIP interpolation without
  extrapolation, PAVA monotonisation, adaptive PCHIP compression of dense curves.
- **ICA/DVA**: dQ/dV and dV/dQ with Savitzky-Golay filtering, peak detection.
- **Literature OCP functions** for examples and testing (no measured data are bundled).

## Installation

```bash
pip install git+https://github.com/catnap042/ocvfit
# with plotting support
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit"
```

Requires Python >= 3.10, NumPy and SciPy.

## Quickstart

```python
import ocvfit as of

# 1. Half-cell OCP curves on a lithiation-fraction axis (here: literature functions;
#    use of.OCPCurve(x, U) or of.OCPCurve.from_soc(...) for your own half-cell data)
cell = of.demo_cells("lfp")               # synthetic LFP/graphite BOL + aged discharge curves
model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)  # 3.45 V rest at full charge

# 2. Fit each low-rate discharge curve (q = Ah discharged from full charge)
bol = of.fit_capacity(model, cell.q["BOL"], cell.voltage["BOL"])
aged = of.fit_capacity(model, cell.q["aged"], cell.voltage["aged"])
print(aged.summary())

# 3. Degradation modes
print(of.degradation_modes(bol, aged))    # {'LLI': 0.080, 'LAM_PE': 0.048, 'LAM_NE': 0.060, ...}

# 4. Is each quantity actually determined by the data? (profile interval, add the bias budget)
print(of.mode_interval(bol, aged, "LAM_PE"))

# 5. Plot
from ocvfit.plotting import plot_fit_report
plot_fit_report({"BOL": bol, "aged": aged}).savefig("fit.png")
```

SOC-axis models use `fit_ocv`:

```python
pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
neg = of.literature_ocp("graphite_chen2020")
fit = of.fit_ocv(of.CathodeFixedModel(pos, neg), soc, voltage, loss="v2")
```

More examples:

- [`examples/quickstart.py`](examples/quickstart.py) - LFP/graphite, full workflow with
  identifiability check (headline figure). `--chemistry nmc532` runs the same workflow on an
  NMC532/graphite cell ([figure](docs/images/nmc532_graphite.png)).
- [`examples/cathode_fixed.py`](examples/cathode_fixed.py) - the three-parameter cathode-fixed
  method on NMC811/graphite ([figure](docs/images/cathode_fixed.png)).
- [`examples/blended_cathode.py`](examples/blended_cathode.py) - NMC811 + LFP blended positive electrode.

### Example results

Ground truth of the synthetic cells: LLI 8 %, LAM_PE 5 %, LAM_NE 6 %, 1 mV
voltage noise, 300 points per curve, end of discharge limited by the negative
electrode. Final interval = correlation-corrected profile interval +/- bias
budget (half-cell curves perturbed by +/-3 mV offset and 2 mV shape error, 8
draws), from `python examples/quickstart.py [--chemistry nmc532]`:

| Cell | Mode | Fitted | True | Final interval | Verdict |
|---|---|---|---|---|---|
| LFP/graphite | LLI | 8.01 % | 8.00 % | [7.87, 8.16] % | identifiable |
| LFP/graphite | LAM_PE | 4.76 % | 5.00 % | [-11.7, 21.1] % | **not identifiable** |
| LFP/graphite | LAM_NE | 6.05 % | 6.00 % | [5.53, 6.58] % | identifiable |
| NMC532/graphite | LLI | 8.00 % | 8.00 % | [7.95, 8.05] % | identifiable |
| NMC532/graphite | LAM_PE | 4.97 % | 5.00 % | [4.67, 5.26] % | identifiable |
| NMC532/graphite | LAM_NE | 6.00 % | 6.00 % | [5.49, 6.51] % | identifiable |

**LAM_PE is weakly identifiable for LFP/graphite in this design.** With
noise only, it is recovered well: over 10 noise seeds LAM_PE = 5.06 +/- 0.26 %
(LLI 8.00 +/- 0.02 %, LAM_NE 6.02 +/- 0.09 %), and the profile interval alone
is about +/-0.7 pp. But the LFP plateau is nearly flat (about 20 mV over the
whole lithiation range), and the cell is negative-limited, so the steep
lithiated end of the LFP curve is never reached. The positive capacity is then
constrained only by the plateau slope, which a few mV of reference-curve error
can mimic. The bias budget for LAM_PE is about 16 pp. LLI and LAM_NE stay well
determined because the graphite staging transitions (the dV/dQ peaks) fix the
negative electrode scale and position. If the discharge reaches the LFP
lithiation end (a positive-limited cell, see
`tests/test_lfp_identifiability.py`), LAM_PE becomes identifiable again (bias
budget < 0.1 pp). For LFP cells, report LAM_PE as an interval unless such an
end point is visible.

## Method

### Conventions

Every electrode curve is a function `U(x)` of the **lithiation fraction** `x`
in `[0, 1]` (potential in V vs Li/Li+), so both potentials decrease with `x`.
`x` denotes the negative and `y` the positive electrode. Half-cell data
recorded against a delithiation SOC are converted with
`OCPCurve.from_soc(soc, U, coordinate="delithiation")`. Curves are interpolated
with PCHIP (shape preserving, no overshoot between plateaus and steep regions)
and **never extrapolated**: points that would require an electrode outside its
measured range are flagged invalid.

### Full-cell voltage

All models evaluate both electrodes **on the same full-cell coordinate** and
subtract:

$$V = U_{PE}(y) - U_{NE}(x) + V_{offset}$$

**Cathode-fixed model** (full-cell SOC $s$, 0 = discharged). The positive
window is fixed to the whole reference curve, and the negative axis is
transformed before it is compared:

$$y(s) = y_{max} - s\,(y_{max}-y_{min}),\qquad u_n(s) = S_n + s/K_n$$

where $u_n$ is the normalised coordinate of the negative reference curve.
$K_n$ is the negative/positive capacity ratio in the window ($K_n>1$ compresses
the negative curve) and $S_n$ the slippage of its start point. A rising $S_n$
indicates lithium inventory loss, a falling $K_n$ loss of negative active
material. Because the positive window is fixed, LAM_PE cannot be detected and
is absorbed by the other parameters, so use this model as a fast baseline.

**Window model.** Both windows are free:

$$y(s) = y_0 + (y_{100}-y_0)\,s,\qquad x(s) = x_0 + (x_{100}-x_0)\,s$$

**Capacity model** (absolute capacity $q$ discharged from full charge, Ah):

$$L_{NE}(q) = x_{100} Q_{NE} - q,\qquad L_{PE}(q) = n_{Li} - x_{100} Q_{NE} + q$$

$$V(q) = U_{PE}\!\left(L_{PE}/Q_{PE}\right) - U_{NE}\!\left(L_{NE}/Q_{NE}\right) - R$$

with $n_{Li} = y_{100} Q_{PE} + x_{100} Q_{NE}$. If the rested full-charge
voltage $V_{full}$ is known, $x_{100}$ is solved from
$U_{PE}(y_{100}) - U_{NE}(x_{100}) = V_{full}$ (full-charge anchor), leaving
$Q_{PE}, Q_{NE}, n_{Li}, R$ free. Otherwise $x_{100}$ is fitted as well. $R$
lumps polarisation and reference mismatch. It is a fitted offset and should
not be read as a DC resistance.

### Blended electrodes

All active materials of an electrode share one potential at equilibrium. The
stored lithium is the sum over components,

$$L_{PE}(U) = \sum_k Q_k\, x_k(U),$$

and the blended potential is its numerical inverse. Each component has its own
capacity $Q_k$ and therefore its own LAM. Outside its measured potential range
a component is held at its end lithiation (phase exhausted or full), and
non-monotone measured segments are repaired with PAVA before inversion.
Potentials or SOCs are never averaged.

### Objectives and optimisation

For SOC-axis models (`fit_ocv`) the recommended `"v2"` objective is

$$J = \frac{\sum_i w_i\,\rho_\delta(r_i)}{\sum_i w_i} + \lambda_{edge}\left[\bar r_{low}^2 + \bar r_{high}^2\right] + \lambda_{valid}\max(0, f_{min} - f_{valid})$$

with Huber loss $\rho_\delta$ ($\delta$ = 2 mV), weights $w_i$ = edge weight
(x10 within 5 % of either end) times $1 + 0.5\,|dV/ds|/\mathrm{median}|dV/ds|$,
the mean residuals $\bar r$ of the two end regions, and a penalty on the
fraction of points outside the valid region (full coverage required by
default, so the optimiser cannot drop poorly fitted end points). Optimisation:
grid or random pre-screening, then L-BFGS-B and a Nelder-Mead polish, or
differential evolution.

For the capacity model (`fit_capacity`), starting points are pre-screened on a
grid of $(Q_{PE}, Q_{NE}, n_{Li})$ (and blend fractions), and the 10 best are
refined with a robust `soft_l1` least-squares fit (5 mV scale) followed by
ordinary least squares. Random starts often land in a solution where the wrong
electrode limits the end of discharge, which is why the grid pre-screen is
used. The few large residuals of the steep end of discharge otherwise make the
valley around the true solution too narrow for plain least squares.

### Degradation modes

$$LLI = 1 - \frac{n_{Li}^{aged}}{n_{Li}^{BOL}},\quad LAM_{PE} = 1 - \frac{Q_{PE}^{aged}}{Q_{PE}^{BOL}},\quad LAM_{NE} = 1 - \frac{Q_{NE}^{aged}}{Q_{NE}^{BOL}}$$

These three quantities are all that OCV data can determine (given unchanged
half-cell shapes, fully disconnected inactive material and uniform ageing).
The common five-mode split (side-reaction LLI plus lithiated/delithiated LAM of
each electrode) has two extra degrees of freedom that leave no trace in the
OCV. `side_reaction_interval` therefore reports the feasible range

$$\left[\max\left(0,\ \Delta n_{Li} - C_{max}\right),\ \Delta n_{Li} - C_{min}\right],\quad
C_{min} = x_0 \Delta Q_{NE} + y_{100} \Delta Q_{PE},\quad C_{max} = x_{100} \Delta Q_{NE} + y_0 \Delta Q_{PE}$$

and a point value only under an explicit assumption (material lost at full charge).

### Identifiability

Residuals along an OCV curve are strongly correlated, so naive confidence
intervals are far too narrow. `profile_interval` fixes one parameter, refits
the others and accepts values with

$$\Delta\chi^2 \le 3.84\,\frac{1+\rho}{1-\rho}$$

where $\rho$ is the lag-1 autocorrelation of the residuals (inflation capped at
100). `bias_budget` perturbs the half-cell curves (default +/-3 mV offset and
2 mV smooth shape error), refits, and takes the largest deviation. A quantity
whose final interval (profile interval +/- bias budget) has a half-width above
1 percentage point should be reported as an interval only (`is_identifiable`).

Whether an electrode or component capacity can be resolved depends on whether
the data show where its plateau starts and ends. For example, an LFP plateau at
the end of discharge in a negative-limited cell is only partially visible (see
[Example results](#example-results)).

## Data sources and citations

No measured data are included in this repository. Examples and tests use
synthetic curves generated from published analytic OCP functions, re-typed
from the [PyBaMM](https://github.com/pybamm-team/PyBaMM) parameter library
(BSD-3-Clause, (c) the PyBaMM team):

| Function | Material | Source |
|---|---|---|
| `graphite_chen2020`, `nmc811_chen2020` | Graphite, NMC811 (LG M50) | C.-H. Chen et al., *J. Electrochem. Soc.* 167 (2020) 080534, [doi:10.1149/1945-7111/ab9050](https://doi.org/10.1149/1945-7111/ab9050) |
| `graphite_mohtat2020`, `nmc532_mohtat2020` | Graphite, NMC532 | P. Mohtat et al., *J. Electrochem. Soc.* 167 (2020) 110561, [doi:10.1149/1945-7111/aba5d1](https://doi.org/10.1149/1945-7111/aba5d1) |
| `lfp_afshar2017` | LFP | S. Afshar, K. Morris, A. Khajepour, "Efficient electrochemical model for lithium-ion cells", arXiv:1709.03970 (2017) |

The LFP/graphite demo uses the same pairing as the PyBaMM `Prada2013`
parameter set: Afshar 2017 LFP OCP + Chen 2020 graphite OCP (the set's cell
parameters come from E. Prada et al., *J. Electrochem. Soc.* 160 (2013) A616,
[doi:10.1149/2.053304jes](https://doi.org/10.1149/2.053304jes), which does not
itself provide OCP functions).

Method background:

- C. R. Birkl, M. R. Roberts, E. McTurk, P. G. Bruce, D. A. Howey, "Degradation diagnostics for lithium ion cells", *J. Power Sources* 341 (2017) 373-386, [doi:10.1016/j.jpowsour.2016.12.011](https://doi.org/10.1016/j.jpowsour.2016.12.011).
- M. Dubarry, C. Truchot, B. Y. Liaw, "Synthesize battery degradation modes via a diagnostic and prognostic model", *J. Power Sources* 219 (2012) 204-216, [doi:10.1016/j.jpowsour.2012.07.016](https://doi.org/10.1016/j.jpowsour.2012.07.016).
- D. Beck, M. Dubarry, "Electrode Blending Simulations Using the Mechanistic Degradation Modes Modeling Approach", *Batteries* 10 (2024) 159, [doi:10.3390/batteries10050159](https://doi.org/10.3390/batteries10050159).
- J. Schmitt et al., "Determination of degradation modes of lithium-ion batteries considering aging-induced changes in the half-cell open-circuit potential curve of silicon-graphite", *J. Power Sources* 532 (2022) 231296, [doi:10.1016/j.jpowsour.2022.231296](https://doi.org/10.1016/j.jpowsour.2022.231296).

Openly licensed measured datasets that are suitable for trying the package on
real cells (not redistributed here; check each licence): Kirkaldy et al. 2024
LG M50T ageing data ([Zenodo 10637534](https://zenodo.org/records/10637534),
CC BY 4.0) and the Dubarry Gr//LFP synthetic library
([Mendeley Data](https://data.mendeley.com/datasets/bs2j56pn7y/3), CC BY 4.0).

## Limitations

- For LFP/graphite cells whose discharge is negative-limited, LAM_PE is
  typically **not identifiable** (flat LFP plateau, LFP end point not
  reached). LLI and LAM_NE remain well determined. See
  [Example results](#example-results).
- Half-cell curve shapes are assumed unchanged by ageing. This fails for
  silicon-containing negatives (see Schmitt et al. 2022), where LAM_NE tends to
  be overestimated.
- Low-rate (e.g. C/20) curves contain polarisation and hysteresis. They are
  *pseudo*-OCV. Pair discharge data with discharge-direction half-cell curves
  (important for LFP), and treat the fitted offset `R` as a nuisance parameter.
- Half-cell axes normalised by their own capacity are operational coordinates,
  not verified absolute stoichiometries.
- Partial voltage windows (e.g. only 100 -> 50 % SOC) can leave LLI and LAM_NE
  unidentifiable; always check intervals.
- Blending assumes equilibrium between components (valid at low rate) and a
  constant component ratio within each curve.
- The five-mode split is not determined by OCV alone; only its feasible range
  is reported.
- Window end points used in `side_reaction_interval` are BOL values; their
  drift during ageing is ignored.

## Development

```bash
pip install -e ".[test]"
pytest
python examples/quickstart.py
```

## License

[MIT](LICENSE) (c) 2026 catnap42. The literature OCP expressions are taken from
the cited publications via the PyBaMM parameter library (BSD-3-Clause).
