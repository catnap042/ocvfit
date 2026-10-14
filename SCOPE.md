# Scope

What ocvfit is for, what it accepts, and what it declines. Ideas that were considered
and declined each have a short record in [`.out-of-scope/`](.out-of-scope), so the same
discussion does not have to happen twice.

## The bar

A change belongs in ocvfit if it makes **non-destructive electrode-alignment
diagnostics from low-rate OCV curves** more correct, more reproducible, or easier for a
person or an AI agent to run, **without** hiding uncertainty and **without** adding
network access or proprietary data.

## Data policy (fixed)

- The repository contains only **synthetic data** generated from published OCP
  functions or closed-form toy functions, and references to **public** datasets.
- It does **not** accept, store or publish company, customer or any other
  real/proprietary measurement data. This applies to issues, pull requests, test
  fixtures and examples too. If a bug only shows up on your data, please reduce it
  to a synthetic or public example.
- The tool runs locally. No telemetry, no network calls, no upload of user files.

## In scope

| area | status |
|---|---|
| single-material positive and negative electrodes | supported |
| blended **positive** electrodes (shared potential) | supported |
| capacity model anchored to `v_full`; window and cathode-fixed models | supported (cathode-fixed is flagged as unreliable on LFP) |
| LLI / LAM_PE / LAM_NE with profile intervals, bias budget and identifiability verdicts | supported |
| dQ/dV and dV/dQ peaks and shifts | supported |
| CLI, MCP server, Agent Skill, Python API | supported |
| blended **negative** electrodes (for example graphite + Si) | not yet; contributions welcome |
| more literature reference OCPs with a clear licence | welcome |

## Out of scope

| idea | why | record |
|---|---|---|
| bundling real or proprietary cell data | data policy above | [real-cell-data.md](.out-of-scope/real-cell-data.md) |
| physics-based P2D/SPM simulation or rate-dependent fitting | different method; mature tools exist | [p2d-simulation.md](.out-of-scope/p2d-simulation.md) |
| reporting a single number when the mode is not identifiable | would hide the main result | [force-point-estimate.md](.out-of-scope/force-point-estimate.md) |
| hosted/cloud service or usage telemetry | local-only by design | [cloud-and-telemetry.md](.out-of-scope/cloud-and-telemetry.md) |
| auto-tuning settings per dataset to minimise RMSE | breaks reproducibility | [auto-tuned-settings.md](.out-of-scope/auto-tuned-settings.md) |

## Already decided

- The identifiability threshold is a 1 pp half-width and is part of the frozen settings
  (`ocvfit settings`).
- Settings are frozen per release. Changing them is a versioned, CHANGELOG-listed change.

## Unclear requests

Open an issue with a synthetic or public example and the result you expected. If it
is declined, it gets a record in `.out-of-scope/`.
