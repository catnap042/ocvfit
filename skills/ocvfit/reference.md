# ocvfit skill reference

Long-form reference for [SKILL.md](SKILL.md). Everything here describes the fixed
behaviour of the `ocvfit` command-line tool (module `ocvfit.workflow`); the MCP server
(`ocvfit-mcp`) calls the same code and returns the same JSON.

## Commands

| command | needs | writes |
|---|---|---|
| `check` | any of `--pe`, `--ne`, `--bol`, `--aged` | `check.json`, `check.md`, `check.png` (input curves) |
| `fit` | `--pe`, `--ne`, `--bol` (optional `--aged`) | `fit.json`, `fit.md`, `fit.png` |
| `modes` | `--pe`, `--ne`, `--bol`, `--aged`; capacity model | `modes.json`, `modes.md` |
| `ica` | `--bol` and/or `--aged` | `ica.json`, `ica.md`, `ica.png`, `ica_<label>.csv` |
| `report` | as `fit` (modes only if `--aged`) | `report.json`, `report.md`, `report.png`, `report_ica.png`, `ica_<label>.csv` |
| `references` | - | prints the built-in references (JSON) |
| `settings` | - | prints the frozen settings (JSON) |

Common options: `--out DIR` (default `ocvfit_out`), `--title`, `--no-figures`, `--json`
(print the full JSON to stdout).

## Exit codes

| code | meaning | what to do |
|---|---|---|
| 0 | success | read `warnings` |
| 1 | internal error | report the message; do not work around it with custom code |
| 2 | usage error (missing or incompatible options) | fix the command line |
| 3 | invalid input data | relay `errors`; fix the file or pass the right option |
| 4 | fit failed (no finite model over the data) | check half-cell coverage / chemistry pairing |

The JSON and Markdown files are written for every exit code, except when the command line
itself cannot be parsed (argparse prints the usage and exits with 2).

## Input checks performed by `check` (and before every fit)

Half-cell curves: file and columns exist; at least 10 numeric points; voltages in V (median
above 10 is treated as mV); x within [0, 1] (percent detected); potential must fall with
lithiation after the coordinate conversion (otherwise: use `--pe-coordinate delithiation`);
non-monotone steps above 2 mV, short x coverage, unusual potential levels and duplicates
give warnings. Built-in references are checked against the electrode role.

Full-cell curves: file and columns exist; at least 20 numeric points; voltages in V;
capacity monotonically non-decreasing and starting near 0 (charge discharged from full
charge); voltage falling with capacity (charge curves are rejected); voltage span, capacity
units and consistency with `--v-full` give warnings.

## JSON document

Top level (all commands): `tool`, `schema_version`, `command`, `status`
(`ok`/`ok_with_warnings`/`error`), `exit_code`, `errors[]`, `warnings[]`, `versions`
(ocvfit, numpy, scipy, python, matplotlib), `options` (the inputs you gave), `settings`
(frozen settings), `inputs` (per-file diagnostics incl. SHA-256), `files` (relative names
in `--out`). The MCP server adds `output_dir` and `absolute_files`.

`fits.<label>` (`label` = `BOL`, `aged`):

| key | meaning |
|---|---|
| `model` | `CapacityModel`, `WindowModel` or `CathodeFixedModel` |
| `params` | fitted parameters (capacity model: `Q_pe` or `Q_<component>`, `Q_ne`, `n_li` in Ah, `R` in V, `x100` if no `--v-full`) |
| `state` | capacity model only: `x100`, `y100`, `x0`, `y0` (electrode lithiation at 100 % / 0 % SOC), `np_ratio` |
| `rmse_mV`, `max_abs_error_mV`, `valid_fraction` | fit quality (unweighted) |
| `parameters_at_bounds` | parameters that ended on a bound (unreliable) |
| `success`, `message` | optimiser status |

`degradation` (`modes` and `report` with `--aged`):

| key | meaning |
|---|---|
| `modes.<mode>.value` | fraction (0.05 = 5 %); `LLI = 1 - n_li(aged)/n_li(BOL)`, `LAM_PE = 1 - Q_pe(aged)/Q_pe(BOL)`, `LAM_NE = 1 - Q_ne(aged)/Q_ne(BOL)` |
| `modes.<mode>.profile_interval` | correlation-corrected profile-likelihood interval |
| `modes.<mode>.bias_budget` | largest change under perturbed half-cell curves (+/-3 mV offset, 2 mV shape error, 8 draws) |
| `modes.<mode>.interval` | final interval = profile interval +/- bias budget |
| `modes.<mode>.half_width` | half of the final interval width (fraction) |
| `modes.<mode>.verdict` | `identifiable` (half-width <= 0.01), `not_identifiable`, `not_assessed` |
| `modes.<mode>.report_as` | how to report it |
| `modes.<mode>.interval_contains_zero` | no detectable change if true |
| `absolute_Ah` | `dn_li_Ah`, `dQ_pe_Ah`, `dQ_ne_Ah` |
| `side_reaction_interval` | feasible range of the side-reaction share of lithium loss (fraction of BOL lithium) |

`differential.<label>`: `dqdv_peaks` (`voltage_V`, `height_Ah_per_V`, `prominence`),
`dvdq_peaks` (`capacity`, `height_V_per_Ah`, `prominence`), `curve_file`;
`differential.dqdv_peak_shifts`: matched BOL/aged peaks within 50 mV and the shift in mV.

## Frozen settings

| area | setting |
|---|---|
| capacity model fit | `fit_capacity`: grid prescreen, best 10 starts, robust soft-L1 least squares (f_scale 5 mV), then linear least squares |
| window / cathode-fixed fit | `fit_ocv`: weighted RMSE loss, fast method (300 random/grid starts, L-BFGS-B, Nelder-Mead), seed 42; abscissa SOC = 1 - q/q_max |
| profile interval | relative step 1e-3, threshold 3.84 (1+rho)/(1-rho), inflation capped at 100 |
| bias budget | 8 draws, +/-3 mV offset, 2 mV smooth shape error, seed 1, 3 starts per refit |
| identifiable | final half-width <= 1 percentage point |
| ICA | 500-point capacity grid, Savitzky-Golay pre-smoothing (31, order 3), dQ/dV on a 5 mV grid (window 9, order 3), dV/dQ window 21; peaks with prominence >= 10 % of max, >= 20 mV apart; dV/dQ peaks searched in 5-95 % of the capacity span |
| literature references | sampled on 2001 points over x in [0, 1] |

`ocvfit settings` prints the exact values. The Python API allows other settings; the CLI,
the skill and the MCP server do not.

## Method in brief

The capacity model writes the full-cell voltage at discharged capacity q as
`V(q) = U_PE((n_li - x100 Q_ne + q)/Q_pe) - U_NE(x100 - q/Q_ne) - R`, with electrode
capacities `Q_pe`, `Q_ne`, cyclable lithium `n_li` and the negative lithiation at full charge
`x100` (solved from `--v-full` if given). Half-cell curves are PCHIP-interpolated and never
extrapolated. Degradation modes compare the BOL and aged fits. Assumptions: half-cell curve
shapes do not change with ageing, ageing is uniform, the curves are near-equilibrium.

## FAQ

- *Why is LAM_PE not identifiable for my LFP cell?* The LFP potential is flat over most of the
  window and the end of discharge is limited by graphite, so the positive capacity barely
  changes the curve. The tool reports an interval; that is the correct result.
- *Can I use charge curves?* Not directly; supply a low-rate discharge curve, or convert a
  charge curve to discharged capacity (q = q_max - q, reversed) and state that you did.
- *No half-cell data?* Use built-in references and state that results carry the
  reference-curve error.
- *Different results on another machine?* The JSON is deterministic for the same inputs and
  library versions; `versions` records them.
- *Blended positive electrodes?* Repeat `--pe`. LAM per component and total LAM_PE are
  reported as indicative (`not_assessed`); LLI and LAM_NE keep their verdicts.

## Example (repository example data)

```bash
ocvfit report --pe examples/data/lfp_halfcell.csv --ne examples/data/graphite_halfcell.csv \
  --bol examples/data/lfp_fullcell_bol.csv --aged examples/data/lfp_fullcell_aged.csv \
  --v-full 3.45 --out results
```

Synthetic ground truth: LLI 8 %, LAM_PE 5 %, LAM_NE 6 %. Expected verdicts: LLI identifiable,
LAM_PE not identifiable (LFP), LAM_NE just above the 1 pp threshold with these short
150-point curves (not identifiable).
