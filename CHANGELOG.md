# Changelog

All notable changes are listed here. The project follows
[Semantic Versioning](https://semver.org/); before 1.0, minor versions may change the API.

## Unreleased

### Added
- `ocvfit` command-line tool with frozen, deterministic workflows: `check`, `fit`,
  `modes`, `ica`, `report` (JSON + Markdown + PNG outputs, documented exit codes).
- Optional MCP server `ocvfit-mcp` (`ocvfit[mcp]`) with the same workflows as tools.
- Agent Skill in `skills/ocvfit/` (SKILL.md, reference, wrapper script).
- `python -m ocvfit`.
- Input limits (25 MiB, 100 000 rows per file), output-directory guards, `--debug`, and the
  optional `OCVFIT_MCP_ROOT` sandbox.
- PyBaMM BSD-3-Clause notice in `LICENSES/`; SECURITY, CONTRIBUTING and issue templates.

### Changed
- Example full-cell data split into `lfp_fullcell_bol.csv` and `lfp_fullcell_aged.csv`.

## 0.1.0 - 2026-10-09

First public version: OCP curves (PCHIP, no extrapolation, PAVA monotonisation),
cathode-fixed / window / capacity full-cell models, blended electrodes, fitting
(`v2`, `rmse`, `mse` objectives), degradation modes with profile intervals, bias budget and
identifiability verdicts, dQ/dV and dV/dQ, plotting helpers, examples and tests.
