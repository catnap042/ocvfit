[README](../README.md) · [Skill](skill.md) · [MCP](mcp.md) · [CLI](cli.md) · [Python API](python-api.md) · [Examples](examples.md) · [Method](method.md) · [Glossary](../GLOSSARY.md)

# Agent Skill guide

[`skills/ocvfit/`](../skills/ocvfit) is an [Agent Skill](https://docs.claude.com/en/docs/claude-code/skills):
a folder with a `SKILL.md` (YAML front matter + instructions), a
[`reference.md`](../skills/ocvfit/reference.md) (JSON schema, frozen settings, method, FAQ)
and one fixed wrapper script. It lets an AI agent run a complete OCV fit and
degradation-mode analysis on your CSV files **without writing any fitting code**.

## Install

```bash
# 1. the tool (Matplotlib included for the figures)
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit.git"
ocvfit --version

# 2. the skill folder
git clone https://github.com/catnap042/ocvfit
mkdir -p ~/.claude/skills && cp -r ocvfit/skills/ocvfit ~/.claude/skills/   # all projects
# or, for one project only:
mkdir -p .claude/skills && cp -r ocvfit/skills/ocvfit .claude/skills/
```

Claude Code picks the skill up from these folders. Other agents that read the Agent
Skills format (for example from `.agents/skills/`) use the same folder unchanged.

## What the skill tells the agent

- **When to use it**: low-rate (pseudo-OCV, C/20 or slower) discharge curves of a cell
  (BOL and/or aged), and a request for electrode alignment, degradation modes, or ICA/DVA.
- **Which inputs it needs**: positive and negative half-cell CSVs (or built-in literature
  references `--pe lfp|nmc811|nmc532`, `--ne graphite|graphite_mohtat2020`), BOL and aged
  full-cell CSVs, and the rested full-charge voltage `--v-full`. If an input is missing,
  the agent asks you instead of inventing it.
- **The exact commands, in order**: `ocvfit check` first, then `ocvfit report`
  (or the single steps `fit`, `modes`, `ica`).
- **How to read the JSON**: fit RMSE and parameters at bounds, identifiability verdicts,
  the LFP LAM_PE caveat, `R` is a lumped offset and not a DC resistance, cathode-fixed is
  only a baseline.
- **What not to do**: no hand-written fitting code, no changed defaults, no averaging of
  runs, no point value for a mode that is not identifiable.
- **How to report**: inputs and warnings, fit quality, a degradation-mode table, the
  figures, and caveats.

## The scripts

The skill contains a single script, `scripts/run_ocvfit.py`. It holds **no analysis
code**: it imports `ocvfit.cli.main` and forwards all arguments, so

```bash
python skills/ocvfit/scripts/run_ocvfit.py report --pe lfp --ne graphite --bol bol.csv --aged aged.csv --out results
```

is the same as `ocvfit report ...`. It exists for agents whose shell does not have the
`ocvfit` console script on `PATH` (`python -m ocvfit` works too). If the package is not
installed, it prints the install command and exits with code 1.

## Determinism

All settings (loss, optimiser, number of starts, seeds, profile and bias-budget settings,
ICA filtering) are frozen inside `ocvfit.workflow`; `ocvfit settings` prints them.
Every JSON file echoes the settings, SHA-256 hashes of the inputs and the library
versions, so two runs on the same inputs and versions give identical JSON. Only
`--model window|cathode-fixed` and the input/column options change what is run.

## Outputs

Written to `--out` (default `ocvfit_out`):

| file | content |
|---|---|
| `check.json` / `check.md` / `check.png` | input validation: columns, V vs mV, monotonicity, capacity origin, axis direction, coverage |
| `fit.json` / `fit.md` / `fit.png` | fitted parameters, electrode windows, RMSE |
| `modes.json` / `modes.md` | LLI / LAM_PE / LAM_NE with final intervals and verdicts |
| `ica.json` / `ica.md` / `ica.png` / `ica_<label>.csv` | dQ/dV and dV/dQ curves and peaks |
| `report.json` / `report.md` / `report.png` / `report_ica.png` | everything above in one run |

See [CLI](cli.md) for every option and [`reference.md`](../skills/ocvfit/reference.md)
for the JSON schema.

## Troubleshooting

| symptom | cause and fix |
|---|---|
| `ocvfit: command not found` | use `python -m ocvfit` or the wrapper script; or point to the venv's `bin/` |
| `ocvfit is not installed` from the wrapper | run the `pip install` line above in the agent's Python |
| exit code 3 | invalid input; read `errors` in `check.json` (wrong column names → `--q-col` etc.; SOC-axis half-cell → `--pe-coordinate delithiation`; file type must be `.csv/.tsv/.txt/.dat`, ≤ 25 MiB, ≤ 100 000 rows) |
| exit code 4 | the fit failed; check coverage, units and `--v-full` |
| warning: figures skipped | install the `plot` extra (Matplotlib) |
| RMSE > ~5 mV or parameters at bounds | reference curves do not match the cell; treat results as unreliable |
| LAM_PE `not_identifiable` on LFP | expected (flat LFP plateau); report the interval only |
| refused: symbolic link in `--out` | ocvfit will not overwrite symlinked result files; pick another directory |
| skill not picked up | check the folder is `~/.claude/skills/ocvfit/SKILL.md` (not nested twice) and restart the agent |
