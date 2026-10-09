"""LFP/graphite: identifiability of the degradation modes.

The LFP plateau is almost flat, so the positive capacity is only determined if
the steep lithiation end of the LFP curve lies inside the measured window. In
the headline (negative-limited) design it does not: LLI and LAM_NE remain
identifiable thanks to the graphite staging features, but LAM_PE is dominated
by the reference-curve bias budget.
"""

import numpy as np
import pytest

import ocvfit as of

MODES = ("LLI", "LAM_PE", "LAM_NE")


def _budget(cell, q, v, ref, n_draws=8):
    def refit(pos, neg):
        m = of.CapacityModel(pos[0], neg, v_full=cell.v_full)
        b = of.fit_capacity(m, q["BOL"], v["BOL"], n_starts=3)
        a = of.fit_capacity(m, q["aged"], v["aged"], n_starts=3)
        return of.degradation_modes(b, a)

    return of.bias_budget(refit, [cell.positive], cell.negative, {k: ref[k] for k in MODES}, n_draws, seed=1)


def _half_width(bol, aged, budget, mode):
    lo, hi = of.mode_interval(bol, aged, mode)
    return 0.5 * (hi - lo) + budget[mode]


def test_lfp_negative_limited_lam_pe_weakly_identifiable(demo, demo_fits):
    bol, aged = demo_fits["BOL"], demo_fits["aged"]
    modes = of.degradation_modes(bol, aged)
    budget = _budget(demo, demo.q, demo.voltage, modes)
    half = {k: _half_width(bol, aged, budget, k) for k in MODES}
    assert of.is_identifiable(half["LLI"])
    assert of.is_identifiable(half["LAM_NE"])
    assert not of.is_identifiable(half["LAM_PE"])  # > 1 pp: report an interval only
    lo, hi = of.mode_interval(bol, aged, "LAM_PE")
    assert lo - budget["LAM_PE"] <= demo.true_modes["LAM_PE"] <= hi + budget["LAM_PE"]


def test_lfp_positive_limited_lam_pe_identifiable(demo):
    # Larger negative electrode and lithium inventory: the LFP lithiation end is reached
    model = of.CapacityModel(demo.positive, demo.negative, v_full=demo.v_full)
    truth = {"BOL": [2.5, 3.0, 2.62, 0.005], "aged": [2.375, 2.82, 2.62 * 0.97, 0.012]}
    q, v, fits = {}, {}, {}
    for i, (k, p) in enumerate(truth.items()):
        q[k], v[k] = of.synthesize_discharge(model, p, v_min=2.5, n=300, noise=1e-3, seed=i)
        fits[k] = of.fit_capacity(model, q[k], v[k])
    modes = of.degradation_modes(fits["BOL"], fits["aged"])
    assert modes["LAM_PE"] == pytest.approx(0.05, abs=0.005)
    assert modes["LLI"] == pytest.approx(0.03, abs=0.005)
    budget = _budget(demo, q, v, modes)
    assert budget["LAM_PE"] < 0.005
    assert np.isfinite(budget["LLI"])
