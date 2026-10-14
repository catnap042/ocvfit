[README](../README.md) · [Skill](skill.md) · [MCP](mcp.md) · [CLI](cli.md) · [Python API](python-api.md) · [Examples](examples.md) · [Method](method.md) · [Glossary](../GLOSSARY.md)

# Python API

All public functions are available from the top-level `ocvfit` namespace. Each
snippet below is self-contained.

## Installation

```bash
pip install git+https://github.com/catnap042/ocvfit
# with plotting support (Matplotlib; recommended for the CLI figures)
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit"
# with the optional MCP server
pip install "ocvfit[mcp] @ git+https://github.com/catnap042/ocvfit"
```

Development (editable) install with test dependencies:

```bash
git clone https://github.com/catnap042/ocvfit
cd ocvfit
pip install -e ".[test,mcp]"
pytest
```

Requires Python 3.10 or newer, NumPy 1.22 or newer and SciPy 1.9 or newer. CI runs the
tests on Linux (Python 3.10 to 3.13) and on Windows and macOS (Python 3.12).


## Quickstart

```python
import ocvfit as of

# Synthetic LFP/graphite cell: BOL and aged low-rate discharge curves with known truth
cell = of.demo_cells("lfp")
model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)  # 3.45 V rest at full charge

# Fit each curve (q = Ah discharged from full charge)
bol = of.fit_capacity(model, cell.q["BOL"], cell.voltage["BOL"])
aged = of.fit_capacity(model, cell.q["aged"], cell.voltage["aged"])
print(aged.summary())

modes = of.degradation_modes(bol, aged)
print({k: round(100 * modes[k], 2) for k in ("LLI", "LAM_PE", "LAM_NE")})  # percent
print("true:", {k: round(100 * v, 2) for k, v in cell.true_modes.items()})

# Profile interval of a mode (add the bias budget for the final interval, see the usage guide)
lo, hi = of.mode_interval(bol, aged, "LAM_NE")
print(f"LAM_NE profile interval: [{100 * lo:.2f}, {100 * hi:.2f}] %")
```


### 1. Half-cell OCP curves

Electrode curves are `OCPCurve` objects. Each is a function `U(x)` of the
**lithiation fraction** `x` (0 = delithiated, 1 = lithiated), in V vs Li/Li+,
so both electrodes' potentials decrease with `x`. Interpolation is PCHIP and
returns `nan` outside the data range.

```python
import numpy as np
import ocvfit as of

# From arrays on a lithiation axis (order and duplicates do not matter)
x = np.linspace(0.0, 1.0, 101)
graphite = of.OCPCurve(x, of.LITERATURE_OCP["graphite_chen2020"](x), name="graphite")

# From a half-cell recorded against SOC that grows with delithiation (x = 1 - soc)
soc = np.linspace(0.0, 1.0, 101)
u_pos = of.LITERATURE_OCP["lfp_afshar2017"](1.0 - soc)
lfp = of.OCPCurve.from_soc(soc, u_pos, coordinate="delithiation", name="LFP", smooth_window=1)

# Published literature curves
nmc = of.literature_ocp("nmc811_chen2020", x_min=0.27, x_max=0.90)
print(sorted(of.LITERATURE_OCP))

print(graphite(0.5), graphite(1.2))   # second value is nan: no extrapolation
print(lfp.x_range, lfp.voltage_range)
print(lfp.lithiation(3.40))           # inverse curve x(U), saturated at the ends

# Repair a slightly non-monotone measured segment (pool-adjacent-violators)
print(of.monotonize([3.50, 3.45, 3.46, 3.40], decreasing=True))

# Compress a dense curve to a few PCHIP knots within 1 mV and rebuild it
xd = np.linspace(0, 1, 3000)
knots, report = of.compress_pchip(xd, of.LITERATURE_OCP["graphite_chen2020"](xd), tol=1e-3)
print(report["n_knots"], report["max_error"])
graphite_small = of.OCPCurve.from_knots(knots, name="graphite")
```

### 2. Full-cell models

Every model maps a parameter vector and an abscissa (SOC, or discharged capacity
for `CapacityModel`) to the full-cell voltage. `components()` also returns
the electrode potentials, the lithiation fractions and the valid mask.

```python
import numpy as np
import ocvfit as of

pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
neg = of.literature_ocp("graphite_chen2020")
soc = np.linspace(0, 1, 50)

# (a) Cathode-fixed: params [Kn, Sn, v_offset]
cf = of.CathodeFixedModel(pos, neg)
v = cf([1.1, 0.05, 0.0], soc)
print(cf.param_names, cf.valid_soc_range([1.1, 0.05, 0.0]))

# (b) Free windows: params [y0, y100, x0, x100, v_offset]
wm = of.WindowModel(of.literature_ocp("nmc811_chen2020"), neg)
v = wm([0.90, 0.27, 0.03, 0.91, 0.0], soc)

# (b') Blended positive: an extra capacity-fraction parameter per additional component
lfp = of.literature_ocp("lfp_afshar2017")
lfp.name = "LFP"
wb = of.WindowModel([of.literature_ocp("nmc811_chen2020"), lfp], neg)
print(wb.param_names)  # [..., 'frac_LFP']

# (c) Capacity model on a capacity axis (Ah): params [Q_pe, Q_ne, n_li, R]
cm = of.CapacityModel(of.literature_ocp("lfp_afshar2017"), neg, v_full=3.45)
q = np.linspace(0, 2.0, 50)
comp = cm.components([2.5, 2.75, 2.4125, 0.005], q)
print(comp.voltage[:3], comp.positive[:3], comp.negative[:3], comp.valid.all())
print(cm.state([2.5, 2.75, 2.4125, 0.005], q_end=2.0))  # x100, y100, x0, y0, N/P ...

# Without a rest-voltage anchor, x100 becomes a free parameter
print(of.CapacityModel(of.literature_ocp("lfp_afshar2017"), neg).param_names)
```

### 3. Fitting

`fit_ocv` fits SOC-axis models. `fit_capacity` fits `CapacityModel`. Both
return a `FitResult` with `params`, `param_dict`, `fitted`, `residuals`,
`rmse`, `mae`, `max_abs_error`, `valid_fraction`, `summary()` and `to_dict()`.

```python
import numpy as np
import ocvfit as of

pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
neg = of.literature_ocp("graphite_chen2020")
model = of.CathodeFixedModel(pos, neg)
soc = np.linspace(0, 1, 120)
v = model([1.1, 0.05, 0.01], soc) + np.random.default_rng(0).normal(0, 1e-3, soc.size)

# loss: "v2" (Huber + edge/slope weights + end anchor), "rmse" (weighted), "mse"
# method: "fast" (pre-screen + L-BFGS-B + Nelder-Mead), "differential_evolution", "local"
fit = of.fit_ocv(model, soc, v, loss="v2", method="fast",
                 huber_delta=2e-3, lambda_edge=50.0, edge_frac=0.05, slope_weight=0.5,
                 min_valid_ratio=1.0)
print(fit.summary())

fit_de = of.fit_ocv(model, soc, v, loss="mse", method="differential_evolution",
                    bounds=[(0.8, 1.3), (-0.2, 0.3), (-0.1, 0.1)], seed=1)
print(fit_de.param_dict, f"{1e3 * fit_de.rmse:.2f} mV")

# Capacity model: grid pre-screen -> best n_starts -> soft_l1 (f_scale) -> least squares
cell = of.demo_cells("nmc532")
cm = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
cfit = of.fit_capacity(cm, cell.q["BOL"], cell.voltage["BOL"], n_starts=10, f_scale=5e-3)
print(cfit.to_dict()["state"])
```

Tuning tips for the `"v2"` objective: increase `lambda_edge` (50, then 100,
then 200) if the curve ends fit poorly. Increase `huber_delta` (2 mV to 3 mV)
if the ends contain outliers. Adjust `edge_frac` (0.03 to 0.08) to define the
end regions.

### 4. Degradation modes and identifiability

```python
import ocvfit as of

cell = of.demo_cells("lfp")
model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
bol = of.fit_capacity(model, cell.q["BOL"], cell.voltage["BOL"])
aged = of.fit_capacity(model, cell.q["aged"], cell.voltage["aged"])
modes = of.degradation_modes(bol, aged)  # fractions; also dn_li_Ah, dQ_pe_Ah, dQ_ne_Ah

# Correlation-corrected profile interval of each mode (from both fits)
intervals = {k: of.mode_interval(bol, aged, k) for k in ("LLI", "LAM_PE", "LAM_NE")}

# Reference-curve bias budget: perturb half-cells (+/-3 mV offset, 2 mV shape), refit
def refit(positives, negative):
    m = of.CapacityModel(positives[0], negative, v_full=cell.v_full)
    b = of.fit_capacity(m, cell.q["BOL"], cell.voltage["BOL"], n_starts=3)
    a = of.fit_capacity(m, cell.q["aged"], cell.voltage["aged"], n_starts=3)
    return of.degradation_modes(b, a)

budget = of.bias_budget(refit, [cell.positive], cell.negative,
                        {k: modes[k] for k in intervals}, n_draws=8, seed=1)

for k, (lo, hi) in intervals.items():
    lo, hi = lo - budget[k], hi + budget[k]
    verdict = "identifiable" if of.is_identifiable(0.5 * (hi - lo)) else "interval only"
    print(f"{k}: {100 * modes[k]:.2f} % in [{100 * lo:.2f}, {100 * hi:.2f}] % -> {verdict}")

# Lower-level: profile one parameter of one fit
pr = of.profile_interval(aged, "n_li")
print(pr.interval, pr.rho, pr.threshold)

# Side-reaction share of the total lithium loss (feasible range, not a point value)
print(of.side_reaction_interval(bol, aged)["side_reaction_interval"])
```

For the 3-parameter cathode-fixed model, `cathode_fixed_indicators(bol_fit,
eol_fit, q_nominal)` returns the first-order indicators `LLI ~ dSn * Q` and
`LAM_NE ~ -dKn * Q`. Use them as indicators only, because the model assumes
LAM_PE = 0.

### 5. dQ/dV and dV/dQ

```python
import ocvfit as of

cell = of.demo_cells("nmc532", noise=0.0)
q, v = cell.q["BOL"], cell.voltage["BOL"]

V, dqdv = of.incremental_capacity(v, q, dv=2e-3, window=11, polyorder=3)  # uniform voltage grid
Q, dvdq = of.differential_voltage(q, v, window=21, polyorder=3)           # uniform capacity grid
peaks = of.find_ic_peaks(V, dqdv, height_ratio=0.1)
print(peaks["voltage"], peaks["height"])
```

Use low-rate data (C/20 to C/25) and a 2 to 5 mV voltage step. Too large a
window broadens peaks or shifts them, so check settings on reference curves.

### 6. Plotting

```python
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import ocvfit as of
from ocvfit.plotting import plot_fit_report

cell = of.demo_cells("lfp")
model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
fits = {k: of.fit_capacity(model, cell.q[k], cell.voltage[k]) for k in ("BOL", "aged")}
fig = plot_fit_report(fits, title="My cell")  # voltage, residuals, electrode potentials, |dV/dQ|
fig.savefig(Path(tempfile.mkdtemp()) / "report.png", dpi=120)
```

### 7. Using your own CSV data

Expected format: comma-separated, one header row, one curve per file.

| File | Columns | Notes |
|---|---|---|
| Half-cell curve | `x,voltage_V` | `x` = lithiation fraction 0-1 (or an SOC column with `coordinate="delithiation"`), voltage in V vs Li/Li+. Use the same current direction as the full-cell curve (important for LFP hysteresis). |
| Full-cell curve (capacity model) | `capacity_Ah,voltage_V` | Capacity discharged **from full charge** (starts at 0, increases), one low-rate discharge branch. |
| Full-cell curve (SOC models) | `soc,voltage_V` | Full-cell SOC 0-1 (0 = discharged). |

Column names are arbitrary, because you pass them explicitly. Example files
(synthetic, made by `examples/make_example_csv.py`) are in
[`examples/data/`](../examples/data). Run this from the repository root:

```python
import ocvfit as of

pos = of.io.load_ocp("examples/data/lfp_halfcell.csv", "x", "voltage_V", name="LFP")
neg = of.io.load_ocp("examples/data/graphite_halfcell.csv", "x", "voltage_V", name="graphite")
q, v = of.io.load_xy("examples/data/lfp_fullcell_aged.csv", "capacity_Ah", "voltage_V")

# v_full: measured rest voltage after a full charge. Pass None to fit x100 instead.
model = of.CapacityModel(pos, neg, v_full=3.45)
fit = of.fit_capacity(model, q, v)
print(fit.summary())

# Optional: store a compressed half-cell curve
knots, report = of.compress_pchip(neg.x, neg.voltage, tol=1e-3)
print(report["n_knots"], "knots")
```

Data that come as capacity in mAh against a delithiation axis can be turned
into a curve by normalising first, for example `soc = Q / Q_max`, then calling
`OCPCurve.from_soc(soc, U, coordinate="delithiation")`.


## Using open datasets

No external data are redistributed here. To try real cells, download an openly
licensed dataset yourself, check its licence, and convert one low-rate
discharge branch plus matching half-cell curves to the CSV format above.

| Dataset | Chemistry | What to use | Licence |
|---|---|---|---|
| Kirkaldy et al. 2024, LG M50T ageing ([Zenodo 10637534](https://zenodo.org/records/10637534)) | NMC811 / graphite-SiOx | C/10 discharge curves from the RPTs as full-cell data. The dataset also ships its own degradation-mode results for comparison. | CC BY 4.0 |
| Dubarry et al., Gr//LFP synthetic library ([Mendeley Data](https://data.mendeley.com/datasets/bs2j56pn7y/3)) | LFP / graphite | Simulated V-Q curves with known LLI/LAM labels. Good for benchmarking. | CC BY 4.0 |
| PyBaMM parameter sets ([docs](https://docs.pybamm.org/en/latest/source/api/parameters/parameter_sets.html)) | several | Analytic half-cell OCPs, already built in via `of.literature_ocp`. | BSD-3-Clause |

Practical notes:

- Pair each full-cell branch with half-cell curves recorded in the same
  direction and at a similarly low rate.
- The Kirkaldy cells have a silicon-containing negative electrode. Graphite-only
  reference curves will then bias LAM_NE (see FAQ).
- Convert capacity to "Ah discharged from full charge". If a rested voltage
  after full charge is available, pass it as `v_full`. Otherwise use
  `v_full=None`, which makes `x100` a free parameter.
