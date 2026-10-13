[README](../README.md) · [Skill](skill.md) · [MCP](mcp.md) · [CLI](cli.md) · [Python API](python-api.md) · [Examples](examples.md) · [Method](method.md)

# Method and theory

### Conventions

Every electrode curve is a function `U(x)` of the **lithiation fraction** `x`
in `[0, 1]` (potential in V vs Li/Li+), so both potentials decrease with `x`.
`x` denotes the negative and `y` the positive electrode. Half-cell data
recorded against a delithiation SOC are converted with
`OCPCurve.from_soc(soc, U, coordinate="delithiation")`. Curves are interpolated
with PCHIP (shape preserving, no overshoot between plateaus and steep regions)
and **never extrapolated**: points that would require an electrode outside its
measured range are flagged invalid.

### Full-cell voltage

All models evaluate both electrodes **on the same full-cell coordinate** and
subtract:

```math
V = U_{PE}(y) - U_{NE}(x) + V_{offset}
```

**Cathode-fixed model** (full-cell SOC $`s`$, 0 = discharged). The positive
window is fixed to the whole reference curve, and the negative axis is
transformed before it is compared:

```math
y(s) = y_{max} - s\,(y_{max}-y_{min}),\qquad u_n(s) = S_n + s/K_n
```

where $`u_n`$ is the normalised coordinate of the negative reference curve.
$`K_n`$ is the negative/positive capacity ratio in the window ($`K_n>1`$ compresses
the negative curve) and $`S_n`$ the slippage of its start point. A rising $`S_n`$
indicates lithium inventory loss, a falling $`K_n`$ loss of negative active
material. Because the positive window is fixed, LAM_PE cannot be detected and
is absorbed by the other parameters, so use this model as a fast baseline.

**Window model.** Both windows are free:

```math
y(s) = y_0 + (y_{100}-y_0)\,s,\qquad x(s) = x_0 + (x_{100}-x_0)\,s
```

**Capacity model** (absolute capacity $`q`$ discharged from full charge, Ah):

```math
L_{NE}(q) = x_{100} Q_{NE} - q,\qquad L_{PE}(q) = n_{Li} - x_{100} Q_{NE} + q
```

```math
V(q) = U_{PE}\!\left(L_{PE}/Q_{PE}\right) - U_{NE}\!\left(L_{NE}/Q_{NE}\right) - R
```

with $`n_{Li} = y_{100} Q_{PE} + x_{100} Q_{NE}`$. If the rested full-charge
voltage $`V_{full}`$ is known, $`x_{100}`$ is solved from
$`U_{PE}(y_{100}) - U_{NE}(x_{100}) = V_{full}`$ (full-charge anchor), leaving
$`Q_{PE}, Q_{NE}, n_{Li}, R`$ free. Otherwise $`x_{100}`$ is fitted as well. $`R`$
lumps polarisation and reference mismatch. It is a fitted offset and should
not be read as a DC resistance.

### Blended electrodes

All active materials of an electrode share one potential at equilibrium. The
stored lithium is the sum over components,

```math
L_{PE}(U) = \sum_k Q_k\, x_k(U),
```

and the blended potential is its numerical inverse. Each component has its own
capacity $`Q_k`$ and therefore its own LAM. Outside its measured potential range
a component is held at its end lithiation (phase exhausted or full), and
non-monotone measured segments are repaired with PAVA before inversion.
Potentials or SOCs are never averaged.

### Objectives and optimisation

For SOC-axis models (`fit_ocv`) the recommended `"v2"` objective is

```math
J = \frac{\sum_i w_i\,\rho_\delta(r_i)}{\sum_i w_i} + \lambda_{edge}\left[\bar r_{low}^2 + \bar r_{high}^2\right] + \lambda_{valid}\max(0, f_{min} - f_{valid})
```

with Huber loss $`\rho_\delta`$ ($`\delta`$ = 2 mV), weights $`w_i`$ = edge weight
(x10 within 5 % of either end) times $`1 + 0.5\,|dV/ds|/\mathrm{median}|dV/ds|`$,
the mean residuals $`\bar r`$ of the two end regions, and a penalty on the
fraction of points outside the valid region (full coverage required by
default, so the optimiser cannot drop poorly fitted end points). Optimisation:
grid or random pre-screening, then L-BFGS-B and a Nelder-Mead polish, or
differential evolution.

For the capacity model (`fit_capacity`), starting points are pre-screened on a
grid of $`(Q_{PE}, Q_{NE}, n_{Li})`$ (and blend fractions), and the 10 best are
refined with a robust `soft_l1` least-squares fit (5 mV scale) followed by
ordinary least squares. Random starts often land in a solution where the wrong
electrode limits the end of discharge, which is why the grid pre-screen is
used. The few large residuals of the steep end of discharge otherwise make the
valley around the true solution too narrow for plain least squares.

### Degradation modes

```math
LLI = 1 - \frac{n_{Li}^{aged}}{n_{Li}^{BOL}},\quad LAM_{PE} = 1 - \frac{Q_{PE}^{aged}}{Q_{PE}^{BOL}},\quad LAM_{NE} = 1 - \frac{Q_{NE}^{aged}}{Q_{NE}^{BOL}}
```

These three quantities are all that OCV data can determine (given unchanged
half-cell shapes, fully disconnected inactive material and uniform ageing).
The common five-mode split (side-reaction LLI plus lithiated/delithiated LAM of
each electrode) has two extra degrees of freedom that leave no trace in the
OCV. `side_reaction_interval` therefore reports the feasible range

```math
\left[\max\left(0,\ \Delta n_{Li} - C_{max}\right),\ \Delta n_{Li} - C_{min}\right],\quad C_{min} = x_0 \Delta Q_{NE} + y_{100} \Delta Q_{PE},\quad C_{max} = x_{100} \Delta Q_{NE} + y_0 \Delta Q_{PE}
```

and a point value only under an explicit assumption (material lost at full charge).

### Identifiability

Residuals along an OCV curve are strongly correlated, so naive confidence
intervals are far too narrow. `profile_interval` fixes one parameter, refits
the others and accepts values with

```math
\Delta\chi^2 \le 3.84\,\frac{1+\rho}{1-\rho}
```

where $`\rho`$ is the lag-1 autocorrelation of the residuals (inflation capped at
100). `bias_budget` perturbs the half-cell curves (default +/-3 mV offset and
2 mV smooth shape error), refits, and takes the largest deviation. A quantity
whose final interval (profile interval +/- bias budget) has a half-width above
1 percentage point should be reported as an interval only (`is_identifiable`).

Whether an electrode or component capacity can be resolved depends on whether
the data show where its plateau starts and ends. For example, an LFP plateau at
the end of discharge in a negative-limited cell is only partially visible (see
[Worked examples](examples.md)).


## FAQ and limitations

**Which model should I use?** Use `CapacityModel` for degradation modes. Use
`WindowModel` when you only have SOC-normalised data. Use `CathodeFixedModel` as
a fast baseline: it assumes LAM_PE = 0.

**Why is LAM_PE not identifiable for my LFP cell?** The LFP plateau is almost
flat. Unless the measured window reaches the steep lithiated end of the LFP
curve (which happens in a positive-limited cell), the positive capacity is only
weakly constrained. Reference-curve errors of a few mV can then shift it by
many percent. Report it as an interval, extend the voltage window, or measure
the positive electrode directly (for example a half-cell from a harvested
electrode).

**Can I get side-reaction LLI versus lithium trapped in dead material?** Not
from OCV data alone. Those two parts are extra degrees of freedom that leave no
trace in the curve. `side_reaction_interval` gives the feasible range and a
point value only under an explicit assumption.

**What does `R` / `v_offset` mean?** It is a lumped offset covering
polarisation at low rate and mismatch between the reference curves. It is not a
DC resistance.

Limitations:

- Half-cell curve shapes are assumed unchanged by ageing. This fails for
  silicon-containing negatives, where LAM_NE tends to be overestimated
  (Schmitt et al. 2022).
- Low-rate curves contain polarisation and hysteresis, so they are pseudo-OCV.
- Half-cell axes normalised to their own capacity are operational coordinates,
  not verified absolute stoichiometries.
- Partial windows (for example only 100 to 50 % SOC) can make LLI and LAM_NE
  unidentifiable. Always check the intervals.
- Blending assumes equilibrium between components (valid at low rate).
- `side_reaction_interval` uses BOL window end points and ignores their drift
  during ageing.
