import numpy as np
import pytest

import ocvfit as of


def test_side_reaction_interval(demo_fits):
    out = of.side_reaction_interval(demo_fits["BOL"], demo_fits["aged"])
    lo, hi = out["side_reaction_interval"]
    assert 0.0 <= lo <= hi <= out["total_LLI"]
    assert lo <= out["if_lost_at_full_charge"] <= hi
    t = out["trapped_lithium_Ah"]
    assert t["min"] <= t["full_charge"] <= t["max"]


def test_cathode_fixed_indicators():
    pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
    neg = of.literature_ocp("graphite_chen2020")
    m = of.CathodeFixedModel(pos, neg)
    s = np.linspace(0, 1, 80)
    bol = of.fit_ocv(m, s, m([1.15, 0.03, 0.0], s))
    eol = of.fit_ocv(m, s, m([1.05, 0.08, 0.0], s))
    ind = of.cathode_fixed_indicators(bol, eol, q_nominal=5.0)
    assert ind["delta_Sn"] == pytest.approx(0.05, abs=2e-3)
    assert ind["LLI_Ah"] > 0 and ind["LAM_NE_Ah"] > 0


def test_lag1_autocorrelation():
    rng = np.random.default_rng(0)
    assert abs(of.lag1_autocorrelation(rng.normal(size=5000))) < 0.05
    e = rng.normal(size=5000)
    ar = np.zeros_like(e)
    for i in range(1, e.size):
        ar[i] = 0.8 * ar[i - 1] + e[i]
    assert of.lag1_autocorrelation(ar) == pytest.approx(0.8, abs=0.05)


def test_profile_interval_covers_truth(demo, demo_fits):
    fit = demo_fits["aged"]
    for name in ("n_li", "Q_pe", "Q_ne"):
        pr = of.profile_interval(fit, name)
        lo, hi = pr.interval
        assert lo < fit.param_dict[name] < hi
        assert lo <= demo.truth["aged"][name] <= hi, name
        assert pr.threshold >= 3.84


def test_bias_budget_runs(demo, demo_fits):
    ref = of.degradation_modes(demo_fits["BOL"], demo_fits["aged"])

    def refit(pos, neg):
        m = of.CapacityModel(pos[0], neg, v_full=demo.v_full)
        b = of.fit_capacity(m, demo.q["BOL"], demo.voltage["BOL"], n_starts=3)
        a = of.fit_capacity(m, demo.q["aged"], demo.voltage["aged"], n_starts=3)
        return of.degradation_modes(b, a)

    dev = of.bias_budget(refit, [demo.positive], demo.negative, {"LLI": ref["LLI"]}, n_draws=2)
    assert 0.0 <= dev["LLI"] < 0.05


def test_is_identifiable():
    assert of.is_identifiable(0.005)
    assert not of.is_identifiable(0.02)
    assert not of.is_identifiable(float("nan"))
