---
name: ocvfit
description: Fit lithium-ion full-cell OCV curves from half-cell OCP curves and quantify degradation modes (LLI, LAM_PE, LAM_NE) with identifiability verdicts, plus dQ/dV and dV/dQ peaks, using the fixed, deterministic ocvfit command-line tool. Use when the user has low-rate (pseudo-OCV) discharge curves of a cell (begin-of-life and/or aged) and wants electrode alignment, degradation modes, or ICA/DVA. Never write custom fitting code; run the fixed commands.
---

# ocvfit: fixed OCV fitting and degradation-mode analysis

Run the fixed `ocvfit` commands and read their JSON output. **Do not write your own
fitting, optimisation or ICA code**, and do not change settings: every setting is frozen
in the tool and echoed in the output.

## Install (once)

```bash
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit.git"
ocvfit --version
```

If `ocvfit` is not on PATH, use `python -m ocvfit` or `python <skill dir>/scripts/run_ocvfit.py`
(same arguments).

To make this skill available to an agent, copy the folder from the repository:
`cp -r skills/ocvfit ~/.claude/skills/` (all projects) or `cp -r skills/ocvfit .claude/skills/`
(one project). Agents that read `.agents/skills/` work the same way.

## Inputs

| input | option | CSV columns (default names) |
|---|---|---|
| positive half-cell | `--pe` (repeat for a blend) | `x` = lithiation fraction 0-1, `voltage_V` (V vs Li/Li+) |
| negative half-cell | `--ne` | same as above |
| BOL full-cell curve | `--bol` | `capacity_Ah` = charge discharged from full charge (starts at 0, increasing), `voltage_V` |
| aged full-cell curve | `--aged` | same as BOL |
| rested full-charge voltage | `--v-full` (and `--v-full-aged`) | number in V; recommended |

- No half-cell data? Use a built-in reference instead of a path: `--pe lfp|nmc811|nmc532`,
  `--ne graphite` (or `graphite_mohtat2020` with `nmc532`). Say so in the report.
- Half-cell x axis is SOC/delithiation (potential rises with x)? Add `--pe-coordinate delithiation`
  (or `--ne-coordinate`); never edit the data yourself.
- Different column names? Pass `--halfcell-x-col`, `--halfcell-v-col`, `--q-col`, `--v-col`.
- Curves must be low-rate (C/20 or slower) discharge branches, voltages in V.

## Commands, in order

```bash
ocvfit check  --pe PE.csv --ne NE.csv --bol BOL.csv --aged AGED.csv --v-full 3.45 --out results
ocvfit report --pe PE.csv --ne NE.csv --bol BOL.csv --aged AGED.csv --v-full 3.45 --out results
```

1. `check` first. Exit code 3 means invalid input: read `errors` in `results/check.json`,
   tell the user exactly what to fix (or rerun with the right option, e.g. a coordinate or
   column flag). Do not "repair" data silently.
2. `report` does check + fit + modes (needs `--aged`) + ICA in one go. The single steps
   are `fit`, `modes`, `ica` with the same options.
3. Read `results/report.json` (machine-readable) and `results/report.md` (summary); show
   `results/report.png` (fit) and `results/report_ica.png` (dQ/dV, dV/dQ) to the user.

Exit codes: 0 ok (check `warnings`), 1 internal error, 2 usage error, 3 invalid input,
4 fit failed. Use `--model window` or `--model cathode-fixed` only when the user asks.

## Interpreting the JSON

- `status`: `ok`, `ok_with_warnings` or `error`. Always relay `warnings`.
- `fits.<BOL|aged>`: `params`, `state` (x100, y100, x0, y0, np_ratio), `rmse_mV`,
  `parameters_at_bounds`. RMSE above about 5 mV, points outside the half-cell range, or
  parameters at bounds mean the result is unreliable: say so.
- `R` is a lumped voltage offset (polarisation + reference mismatch) at the test current,
  **not** a DC resistance.
- `degradation.modes.<LLI|LAM_PE|LAM_NE>`: `value`, `interval` (profile interval plus
  reference-curve bias budget), `half_width`, `verdict`, `report_as`.
  - `verdict == "identifiable"` (half-width at most 1 percentage point): report value and interval.
  - `"not_identifiable"`: report **the interval only**, no point value.
  - `"not_assessed"` (blends): indicative only.
- LFP: LAM_PE is usually not identifiable (flat LFP plateau, negative-limited end of
  discharge). This is expected, not a failure.
- `side_reaction_interval`: OCV data give only a range for the side-reaction share of LLI.
- Cathode-fixed results (`cathode_fixed_indicators`) are a first-order baseline that assumes
  no LAM_PE; never present them as degradation modes.
- `differential`: dQ/dV peak voltages, dV/dQ peak capacities, BOL-to-aged peak shifts.

## Do not

- write or run ad-hoc fitting code, change frozen settings, or post-process numbers
  yourself (no re-fitting, no averaging of runs);
- claim a point value for a mode whose verdict is not `identifiable`;
- call `R` a resistance, or treat literature-reference results as cell-specific truth;
- invent missing inputs (for example a `--v-full` value) - ask the user instead.

## Reporting template

1. Inputs used (files or built-in references, v_full) and any warnings.
2. Fit quality: RMSE (mV) for BOL and aged.
3. Degradation modes table: mode, value (identifiable only), final interval, verdict.
4. Figures: `report.png`, `report_ica.png`.
5. Caveats: identifiability, reference-curve assumptions, R meaning.

See [reference.md](reference.md) for the JSON schema, frozen settings, method and FAQ.
