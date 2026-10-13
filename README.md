# ocvfit

[![CI](https://github.com/catnap042/ocvfit/actions/workflows/ci.yml/badge.svg)](https://github.com/catnap042/ocvfit/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://github.com/catnap042/ocvfit/blob/main/pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Low-rate full-cell discharge curves + half-cell OCP curves in → electrode alignment,
degradation modes (LLI, LAM_PE, LAM_NE) with identifiability verdicts, and dQ/dV / dV/dQ
out, as JSON, Markdown and PNG.**

Ask your AI agent to *"check my CSV and decompose the ageing"* and it runs a fixed,
deterministic fitting pipeline instead of writing its own code. The same pipeline is
available as an [Agent Skill](#quickstart-with-an-ai-agent-skill), an
[MCP server](#mcp-server), a [command-line tool](#command-line-tool) and a
[Python API](docs/python-api.md).

| LFP/graphite (headline example) | NMC532/graphite |
|---|---|
| ![LFP/graphite](docs/images/quickstart.png) | ![NMC532/graphite](docs/images/nmc532_graphite.png) |

## What it is and isn't

| it is | it isn't |
|---|---|
| non-destructive diagnostics by *electrode alignment* of pseudo-OCV curves | a physics-based (P2D/SPM) cell simulator |
| a fixed pipeline with frozen settings, so results are reproducible | a tool that tunes itself to give "nice" numbers |
| honest about uncertainty: a mode that the data do not determine is reported as an interval only | a source of cell-specific truth when literature reference curves are used |
| single-material and blended positive electrodes | blended *negative* electrodes (not yet supported) |
| local only: no telemetry, no network connections | a cloud service |

> [!IMPORTANT]
> All bundled examples are **synthetic data built from published OCP functions** (no
> measured data are included). Results on your own cells depend on how well the half-cell
> reference curves match your electrodes; validate them (for example against
> harvested-electrode half-cells or post-mortem analysis) before drawing conclusions.

## Quickstart with an AI agent (Skill)

**1. Install the tool and the skill** (once):

```bash
pip install "ocvfit[plot] @ git+https://github.com/catnap042/ocvfit.git"
git clone https://github.com/catnap042/ocvfit
mkdir -p ~/.claude/skills && cp -r ocvfit/skills/ocvfit ~/.claude/skills/   # all projects
# or one project: mkdir -p .claude/skills && cp -r ocvfit/skills/ocvfit .claude/skills/
```

Claude Code loads the skill from `~/.claude/skills/` or `.claude/skills/`; other agents
that read the Agent Skills format use the same folder (for example in `.agents/skills/`).

**2. Ask your agent.** Copy-paste one of these (the example files are in
[`examples/data/`](examples/data)):

```text
Use ocvfit to check my data: positive half-cell examples/data/lfp_halfcell.csv,
negative half-cell examples/data/graphite_halfcell.csv, BOL curve
examples/data/lfp_fullcell_bol.csv, rested full-charge voltage 3.45 V.
```

```text
Fit the BOL curve examples/data/lfp_fullcell_bol.csv and the aged curve
examples/data/lfp_fullcell_aged.csv with ocvfit (LFP and graphite half-cells in
examples/data, v_full 3.45 V) and decompose the ageing into LLI, LAM_PE and LAM_NE.
Tell me which modes are identifiable.
```

```text
Run an ocvfit dQ/dV and dV/dQ analysis on my BOL and aged curves and list the peak
shifts. I have no half-cell data, so use the built-in lfp and graphite references.
```

```text
Make a full ocvfit report for my cell and summarise it in a table with the figures.
```

**3. What you get** (in `results/` or the folder you name):

| output | content |
|---|---|
| `report.json` | everything machine-readable: settings, input hashes, fits, modes, peaks, warnings |
| `report.md` | human summary: fit RMSE, a degradation-mode table with intervals and verdicts |
| `report.png` | four-panel fit report (BOL and aged) |
| `report_ica.png` | dQ/dV and dV/dQ, BOL vs aged |

The agent first runs `ocvfit check`, then `ocvfit report`, reads the JSON, and reports
point values only for identifiable modes. Example from `report.md` (synthetic truth: LLI 8 %,
LAM_PE 5 %, LAM_NE 6 %):

```text
| mode | value | final interval | half-width | verdict |
|---|---|---|---|---|
| LLI | 7.98 % | [7.65, 8.29] % | 0.32 pp | identifiable |
| LAM_PE | (5.75 %) | [-23.54, 34.95] % | 29.24 pp | not_identifiable |
| LAM_NE | (5.84 %) | [4.26, 7.20] % | 1.47 pp | not_identifiable |
```

Full guide (what the skill does, determinism, outputs, troubleshooting): **[docs/skill.md](docs/skill.md)**.

## MCP server

For MCP clients (Claude Desktop, Cursor, Claude Code, ...):

```bash
pip install "ocvfit[mcp] @ git+https://github.com/catnap042/ocvfit"
claude mcp add --transport stdio --scope user ocvfit -- ocvfit-mcp   # Claude Code
```

Claude Desktop (`claude_desktop_config.json`) or Cursor (`~/.cursor/mcp.json`):

```json
{
  "mcpServers": {
    "ocvfit": {
      "command": "/path/to/venv/bin/ocvfit-mcp",
      "args": [],
      "env": {"OCVFIT_MCP_ROOT": "/path/to/data"}
    }
  }
}
```

`OCVFIT_MCP_ROOT` is optional and confines reads and writes to that directory tree.
Use the absolute path from `which ocvfit-mcp`.

| tool | CLI equivalent |
|---|---|
| `check_data` | `ocvfit check` |
| `fit_cell` | `ocvfit fit` |
| `degradation_modes` | `ocvfit modes` |
| `differential_analysis` | `ocvfit ica` |
| `full_report` | `ocvfit report` |
| `literature_references` | `ocvfit references` |

Details: [docs/mcp.md](docs/mcp.md).

## Command-line tool

```bash
ocvfit report --pe examples/data/lfp_halfcell.csv --ne examples/data/graphite_halfcell.csv \
              --bol examples/data/lfp_fullcell_bol.csv --aged examples/data/lfp_fullcell_aged.csv \
              --v-full 3.45 --out results
```

| command | does |
|---|---|
| `ocvfit check` | validate inputs (columns, units, monotonicity, coverage) |
| `ocvfit fit` | fit BOL (and aged) curves |
| `ocvfit modes` | LLI / LAM_PE / LAM_NE with intervals and verdicts |
| `ocvfit ica` | dQ/dV and dV/dQ peaks and shifts |
| `ocvfit report` | all of the above in one run |
| `ocvfit references` / `ocvfit settings` | built-in reference curves / frozen settings |

Exit codes 0 ok, 1 internal, 2 usage, 3 invalid input, 4 fit failed. Details:
[docs/cli.md](docs/cli.md).

## Python API

```python
import ocvfit as of

cell = of.demo_cells("lfp")  # synthetic LFP/graphite cell with known truth
model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
bol = of.fit_capacity(model, cell.q["BOL"], cell.voltage["BOL"])
aged = of.fit_capacity(model, cell.q["aged"], cell.voltage["aged"])
print(of.degradation_modes(bol, aged))
```

Half-cell curves, the three models, fitting options, identifiability, ICA/DVA, plotting
and your own CSV data: **[docs/python-api.md](docs/python-api.md)**.

## Features

- **Half-cell OCP curves**: shape-preserving PCHIP interpolation that never
  extrapolates, conversion from delithiation-SOC axes, PAVA monotonisation,
  inverse curves, and adaptive PCHIP compression of dense data.
- **Three full-cell models**:
  - `CathodeFixedModel`: 3 parameters (`Kn`, `Sn`, `V_offset`). A fast baseline.
  - `WindowModel`: free electrode windows (`y0`, `y100`, `x0`, `x100`, `V_offset`) and optional blend fractions.
  - `CapacityModel`: absolute capacities `Q_PE`, `Q_NE`, cyclable lithium `n_Li` and offset `R`, with an optional full-charge voltage anchor.
- **Blended positive electrodes** (for example NMC + LFP) with one capacity per
  component. Lithium is shared between components at a common potential.
- **Robust fitting**: Huber loss, edge and slope-aware weights, end-point
  anchoring and a valid-region penalty. Optimisers are grid pre-screening,
  L-BFGS-B, Nelder-Mead and differential evolution. The capacity model uses a
  two-stage `soft_l1` then least-squares fit.
- **Degradation modes**: LLI, LAM_PE (per component) and LAM_NE, plus the
  feasible range of side-reaction lithium loss.
- **Identifiability**: profile-likelihood intervals corrected for residual
  correlation, a bias budget from perturbed reference curves, and a 1
  percentage-point verdict.
- **ICA/DVA**: dQ/dV and dV/dQ with Savitzky-Golay filtering, and peak detection.
- **Plotting**: a four-panel fit report.
- **Literature OCP functions**: LFP, NMC811, NMC532 and two graphites, for
  examples, tests and synthetic studies.
- **Agent-ready**: a deterministic `ocvfit` CLI with frozen settings and JSON/PNG/Markdown
  output, an Agent Skill, and an optional MCP server
  (see the [Skill quickstart](#quickstart-with-an-ai-agent-skill)).

## Privacy and network

| | |
|---|---|
| telemetry, analytics, update checks | none |
| network connections | none (library, CLI and MCP server) |
| shell commands from the MCP server | none |
| files read | only the input files you pass |
| files written | only ocvfit's own result names in `--out` / `out_dir`; refused if one is a symbolic link |
| MCP sandbox | optional `OCVFIT_MCP_ROOT`: paths outside it are refused |
| input limits | `.csv`, `.tsv`, `.txt`, `.dat` only; 25 MiB and 100 000 data rows per file |

## Documentation

| page | content |
|---|---|
| [docs/skill.md](docs/skill.md) | Agent Skill: install, what it does, determinism, outputs, troubleshooting |
| [docs/mcp.md](docs/mcp.md) | MCP server: tools, client configuration, sandbox |
| [docs/cli.md](docs/cli.md) | command-line tool: commands, options, outputs, exit codes |
| [docs/python-api.md](docs/python-api.md) | Python API: installation, usage guide, own data, open datasets |
| [docs/examples.md](docs/examples.md) | worked LFP/graphite and NMC532/graphite examples |
| [docs/method.md](docs/method.md) | method and theory, FAQ and limitations |

## Citing and references

If you use `ocvfit`, please cite the repository
(`https://github.com/catnap042/ocvfit`) and the publications below that your
analysis depends on.

OCP functions (re-typed from the [PyBaMM](https://github.com/pybamm-team/PyBaMM)
parameter library, BSD-3-Clause, (c) the PyBaMM team):

| Function | Material | Source |
|---|---|---|
| `graphite_chen2020`, `nmc811_chen2020` | Graphite, NMC811 (LG M50) | C.-H. Chen et al., *J. Electrochem. Soc.* 167 (2020) 080534, [doi:10.1149/1945-7111/ab9050](https://doi.org/10.1149/1945-7111/ab9050) |
| `graphite_mohtat2020`, `nmc532_mohtat2020` | Graphite, NMC532 | P. Mohtat et al., *J. Electrochem. Soc.* 167 (2020) 110561, [doi:10.1149/1945-7111/aba5d1](https://doi.org/10.1149/1945-7111/aba5d1) |
| `lfp_afshar2017` | LFP | S. Afshar, K. Morris, A. Khajepour, "Efficient electrochemical model for lithium-ion cells", arXiv:1709.03970 (2017) |

The LFP/graphite example uses the pairing of the PyBaMM `Prada2013` parameter
set: Afshar 2017 LFP plus Chen 2020 graphite. The set's cell parameters come
from E. Prada et al., *J. Electrochem. Soc.* 160 (2013) A616,
[doi:10.1149/2.053304jes](https://doi.org/10.1149/2.053304jes), which does not
itself provide OCP functions.

Method background:

- C. R. Birkl, M. R. Roberts, E. McTurk, P. G. Bruce, D. A. Howey, "Degradation diagnostics for lithium ion cells", *J. Power Sources* 341 (2017) 373-386, [doi:10.1016/j.jpowsour.2016.12.011](https://doi.org/10.1016/j.jpowsour.2016.12.011).
- M. Dubarry, C. Truchot, B. Y. Liaw, "Synthesize battery degradation modes via a diagnostic and prognostic model", *J. Power Sources* 219 (2012) 204-216, [doi:10.1016/j.jpowsour.2012.07.016](https://doi.org/10.1016/j.jpowsour.2012.07.016).
- D. Beck, M. Dubarry, "Electrode Blending Simulations Using the Mechanistic Degradation Modes Modeling Approach", *Batteries* 10 (2024) 159, [doi:10.3390/batteries10050159](https://doi.org/10.3390/batteries10050159).
- J. Schmitt et al., "Determination of degradation modes of lithium-ion batteries considering aging-induced changes in the half-cell open-circuit potential curve of silicon-graphite", *J. Power Sources* 532 (2022) 231296, [doi:10.1016/j.jpowsour.2022.231296](https://doi.org/10.1016/j.jpowsour.2022.231296).

## Contributing

Questions, bug reports and pull requests are welcome via
[GitHub Issues](https://github.com/catnap042/ocvfit/issues); see
[CONTRIBUTING.md](CONTRIBUTING.md) for the development setup and
[SECURITY.md](SECURITY.md) for reporting security problems privately. Changes are listed
in [CHANGELOG.md](CHANGELOG.md).

## License

[MIT](LICENSE) (c) 2026 catnap042. The literature OCP expressions come from the
cited publications via the PyBaMM parameter library (BSD-3-Clause); its copyright and
license notice is in [`LICENSES/PyBaMM-BSD-3-Clause.txt`](LICENSES/PyBaMM-BSD-3-Clause.txt)
and ships with the package. The CSV files in `examples/data/` are synthetic and covered by
the MIT licence. No measured or third-party data are included.
