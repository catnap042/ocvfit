[README](../README.md) · [Skill](skill.md) · [MCP](mcp.md) · [CLI](cli.md) · [Python API](python-api.md) · [Examples](examples.md) · [Method](method.md)

# Command-line tool

The CLI, the [Agent Skill](skill.md) and the [MCP server](mcp.md) all run the **same
fixed, deterministic code path** (`ocvfit.workflow`): settings are frozen, seeds are
fixed, and every result file echoes the settings, the input hashes and the library
versions, so two runs on the same inputs give identical JSON.

Installed with the package as the `ocvfit` console script (or `python -m ocvfit`):

```bash
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit"

ocvfit check  --pe examples/data/lfp_halfcell.csv --ne examples/data/graphite_halfcell.csv \
              --bol examples/data/lfp_fullcell_bol.csv --aged examples/data/lfp_fullcell_aged.csv \
              --v-full 3.45 --out results
ocvfit report --pe examples/data/lfp_halfcell.csv --ne examples/data/graphite_halfcell.csv \
              --bol examples/data/lfp_fullcell_bol.csv --aged examples/data/lfp_fullcell_aged.csv \
              --v-full 3.45 --out results
```

| command | does | writes to `--out` |
|---|---|---|
| `check` | validates columns, units (V vs mV), monotonicity, capacity origin, axis direction, coverage | `check.json`, `check.md`, `check.png` |
| `fit` | fits BOL (and aged) curves; capacity model by default (`--model window` / `cathode-fixed` optional) | `fit.json`, `fit.md`, `fit.png` |
| `modes` | LLI / LAM_PE / LAM_NE with profile interval + bias budget and an identifiability verdict | `modes.json`, `modes.md` |
| `ica` | dQ/dV and dV/dQ, peak positions, BOL-to-aged peak shifts | `ica.json`, `ica.md`, `ica.png`, `ica_<label>.csv` |
| `report` | check + fit + modes (with `--aged`) + ica | `report.json`, `report.md`, `report.png`, `report_ica.png` |

Inputs use the [CSV format above](python-api.md#7-using-your-own-csv-data) (other column names via
`--halfcell-x-col`, `--halfcell-v-col`, `--q-col`, `--v-col`; SOC-axis half-cells via
`--pe-coordinate delithiation`). Without half-cell data, use a built-in reference in
place of a path: `--pe lfp|nmc811|nmc532`, `--ne graphite|graphite_mohtat2020`
(`ocvfit references` lists them, `ocvfit settings` prints the frozen settings).
Exit codes: 0 success (check `warnings`), 1 internal error, 2 usage error, 3 invalid
input, 4 fit failed. Errors are reported as messages, without a Python traceback (add
`--debug` to see one). Input files must be `.csv`, `.tsv`, `.txt` or `.dat` and are limited to 25 MiB and
100 000 data rows each.
Results overwrite only ocvfit's own file names in `--out`; the run is refused if one of
them is a symbolic link.

Excerpt of `results/report.md` for the example data (synthetic truth: LLI 8 %,
LAM_PE 5 %, LAM_NE 6 %; 150 points per curve):

```text
| mode | value | final interval | half-width | verdict |
|---|---|---|---|---|
| LLI | 7.98 % | [7.65, 8.29] % | 0.32 pp | identifiable |
| LAM_PE | (5.75 %) | [-23.54, 34.95] % | 29.24 pp | not_identifiable |
| LAM_NE | (5.84 %) | [4.26, 7.20] % | 1.47 pp | not_identifiable |
```

With these shorter, noisier example curves, LAM_NE ends up just above the 1 pp
threshold, while it is identifiable in the 300-point worked example below. The verdict
follows the data, not the chemistry label.


## Exit codes

| code | meaning |
|---|---|
| 0 | success (still read `warnings`) |
| 1 | internal error |
| 2 | usage error (bad options) |
| 3 | invalid input (see `errors` in the JSON) |
| 4 | fit failed |

Other commands: `ocvfit references` (built-in literature curves), `ocvfit settings`
(frozen settings), `ocvfit --version`. Output options: `--out DIR` (default
`ocvfit_out`), `--title`, `--no-figures`, `--json` (print the full JSON to stdout),
`--debug`.
