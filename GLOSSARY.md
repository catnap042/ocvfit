# Glossary

The terms used in the ocvfit code, docs, reports and Skill, and the names they go by
in the code. Use these words in issues and pull requests. Words under *Avoid* mean
the same thing but cause confusion.

## Curves and electrodes

**OCP (open-circuit potential)**:
Equilibrium potential of a single electrode material against Li/Li+, as a function
of its lithiation. Measured on a **half-cell**. Code: `OCPCurve`, `--pe` / `--ne`.
_Avoid_: half-cell OCV, electrode OCV

**OCV (open-circuit voltage)**:
Equilibrium voltage of the **full-cell**: positive OCP minus negative OCP. ocvfit
works with *pseudo-OCV*, a low-rate (for example C/20 or slower) discharge curve.
The small remaining polarisation is absorbed by `R`.
_Avoid_: OCP (for the full cell)

**Half-cell**:
One electrode material against a lithium-metal counter electrode. It gives an
**OCP** curve.

**Full-cell**:
Positive electrode against negative electrode, the cell you actually age. It gives
the **OCV** / discharge curves (`--bol`, `--aged`).

**Positive electrode (PE)** / **Negative electrode (NE)**:
Used instead of cathode/anode so that the meaning stays the same on charge and
discharge. Code suffixes `_pe` / `_ne`, mode names `LAM_PE` / `LAM_NE`.
_Avoid_: cathode, anode (except in the historical model name `CathodeFixedModel`)

**Blended electrode**:
An electrode made of several active materials that share one potential (for example
LFP + NMC). Supported for the **positive** electrode (`BlendedElectrode`, one
`Q_<name>` per component). Blended negative electrodes are not supported yet; see
[SCOPE.md](SCOPE.md).

**LFP**:
Lithium iron phosphate positive electrode. Its OCP is almost flat across most of
the lithiation range, which is why LAM_PE is often **not identifiable** for LFP cells.

## Coordinates and capacities

**Lithiation x (negative), y (positive)**:
Fraction of the lithium sites occupied in the material, from 0 to 1. `x` is always
the negative electrode, `y` always the positive. `x100`, `y100` are the values at
full charge and `x0`, `y0` the values at the end of the measured discharge.
_Avoid_: stoichiometry (fine in prose, but the code says lithiation), SOC (for an electrode)

**SOC (state of charge)**:
Full-cell state of charge, 1 at full charge. Only the window and cathode-fixed
models use an SOC axis. The capacity model uses discharged capacity `q` (Ah).

**q**:
Capacity discharged from the fully charged state, in Ah, increasing during discharge.
It is the x-axis of `--bol` / `--aged` for the capacity model.

**Q_pe, Q_ne**:
Active capacity of each electrode (Ah): the capacity it would deliver from lithiation
0 to 1.

**n_li**:
Cyclable lithium inventory (Ah) shared by the two electrodes.

**np_ratio**:
`Q_ne / Q_pe`.

**v_full**:
Rested full-cell voltage at full charge (V), given with `--v-full`. It anchors `x100`
so that `x100` is solved rather than fitted.

**R**:
Constant voltage offset (V) that lumps low-rate polarisation and reference mismatch.
It is not a cell resistance.

## Ageing

**BOL (beginning of life)**:
The reference state that ageing is measured against (`--bol`). Every mode is relative
to the BOL fit.

**Aged**:
Any later state of the same cell (`--aged`).

**LLI (loss of lithium inventory)**:
`1 - n_li(aged) / n_li(BOL)`. Reported as a fraction in code and as % in reports.

**LAM_PE / LAM_NE (loss of active material)**:
`1 - Q_pe(aged) / Q_pe(BOL)` and the same for `Q_ne`. Blends add `LAM_<component>`.

**Degradation modes**:
LLI, LAM_PE and LAM_NE together. Code: `degradation_modes()`, `ocvfit modes`.

**Side-reaction interval**:
The range of LLI that may come from side reactions rather than lithium trapped in
lost active material. OCV data alone cannot pin it down.

## Analysis and trust

**dQ/dV (ICA)** / **dV/dQ (DVA)**:
Incremental-capacity and differential-voltage curves. Their peaks are fingerprints of
electrode phase transitions. Code: `ocvfit ica`, `differential_analysis`.

**Identifiability**:
Whether the data actually determine a mode. The half-width of the final interval
(profile interval plus reference-curve **bias budget**) must be 1 percentage point
(pp) or less (`is_identifiable(half_width, threshold=0.01)`). Modes that fail are
reported as intervals only, with the value in parentheses.
_Avoid_: confidence (it is not a probability statement)

**Profile interval**:
Interval of one parameter found by re-optimising all the others
(`profile_interval`, `mode_interval`).

**Bias budget**:
How much a mode moves when the half-cell reference curves are perturbed within a
realistic range (`bias_budget`). It covers reference-curve mismatch, which noise alone
does not.

**pp (percentage point)**:
Unit for differences and half-widths of modes given in %.

**Synthetic data**:
Curves generated from published OCP functions or closed-form toy functions with a
known answer. They are the only data that ship with this repository.
