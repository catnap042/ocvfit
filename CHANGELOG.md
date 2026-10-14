# Changelog

All notable changes are listed here. The project follows
[Semantic Versioning](https://semver.org/); before 1.0, minor versions may change the API.

## Unreleased

### Docs
- README opens with the problem it solves and a 30-second install + quickstart on the
  bundled synthetic data.
- New `GLOSSARY.md` (terms matched to code names) and `SCOPE.md` (data policy,
  supported vs not yet, declined ideas). Declined ideas are recorded in `.out-of-scope/`.

### Tests
- New `tests/test_known_answers.py`: closed-form toy OCPs and full-cell curves built
  without ocvfit code; LLI/LAM_PE/LAM_NE recovered within 0.5 pp over 3 noise seeds,
  plus a gate test that must flag LAM_PE as not identifiable for a flat positive OCP.

## 0.2.1 - 2026-10-13

### Docs
- README restructured around the Agent Skill and MCP server quickstarts; detailed
  usage moved to `docs/` (skill, mcp, cli, python-api, examples, method).
- Author name unified to `catnap042` (LICENSE, pyproject).

## 0.2.0 - 2026-10-09

### Added
- `ocvfit` command-line tool with frozen, deterministic workflows: `check`, `fit`,
  `modes`, `ica`, `report` (JSON + Markdown + PNG outputs, documented exit codes).
- Optional MCP server `ocvfit-mcp` (`ocvfit[mcp]`) with the same workflows as tools.
- Agent Skill in `skills/ocvfit/` (SKILL.md, reference, wrapper script).
- `python -m ocvfit`.
- Input limits (.csv/.tsv/.txt/.dat, 25 MiB, 100 000 rows per file), output-directory
  guards, a warning when figures are skipped, `--debug`, and the
  optional `OCVFIT_MCP_ROOT` sandbox.
- PyBaMM BSD-3-Clause notice in `LICENSES/`; SECURITY, CONTRIBUTING and issue templates.

### Changed
- Example full-cell data split into `lfp_fullcell_bol.csv` and `lfp_fullcell_aged.csv`.

## 0.1.0 - 2026-10-09

First public version: OCP curves (PCHIP, no extrapolation, PAVA monotonisation),
cathode-fixed / window / capacity full-cell models, blended electrodes, fitting
(`v2`, `rmse`, `mse` objectives), degradation modes with profile intervals, bias budget and
identifiability verdicts, dQ/dV and dV/dQ, plotting helpers, examples and tests.
