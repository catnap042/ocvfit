import numpy as np
import pytest

import ocvfit as of


@pytest.fixture
def curves():
    return of.literature_ocp("nmc811_chen2020", 0.27, 0.90), of.literature_ocp("graphite_chen2020")


def test_cathode_fixed_valid_range(curves):
    pos, neg = curves
    m = of.CathodeFixedModel(pos, neg)
    s = np.linspace(0, 1, 101)
    comp = m.components([1.2, -0.1, 0.0], s)
    lo, hi = m.valid_soc_range([1.2, -0.1, 0.0])
    assert lo == pytest.approx(0.12)
    assert np.all(~comp.valid[s < lo - 1e-9])
    assert np.all(comp.valid[(s > lo + 1e-9) & (s <= hi)])


def test_cathode_fixed_endpoints(curves):
    pos, neg = curves
    m = of.CathodeFixedModel(pos, neg)
    v = m([1.0, 0.0, 0.01], np.array([0.0, 1.0]))
    assert v[0] == pytest.approx(pos(0.90) - neg(0.0) + 0.01)
    assert v[1] == pytest.approx(pos(0.27) - neg(1.0) + 0.01)


def test_window_model_formula(curves):
    pos, neg = of.literature_ocp("nmc811_chen2020"), curves[1]
    m = of.WindowModel(pos, neg)
    p = [0.9, 0.27, 0.03, 0.91, 0.005]
    s = np.array([0.0, 0.5, 1.0])
    y = 0.9 + (0.27 - 0.9) * s
    x = 0.03 + (0.91 - 0.03) * s
    np.testing.assert_allclose(m(p, s), pos(y) - neg(x) + 0.005)


def test_single_component_blend_matches_curve():
    c = of.literature_ocp("nmc532_mohtat2020")
    b = of.BlendedElectrode([c])
    y = np.linspace(0.1, 0.9, 9)
    np.testing.assert_allclose(b.potential(y * 2.0, [2.0]), c(y), atol=5e-4)


def test_blend_iso_potential_lithium_balance():
    a = of.literature_ocp("nmc811_chen2020")
    b = of.literature_ocp("lfp_afshar2017")
    blend = of.BlendedElectrode([a, b])
    q = np.array([3.0, 2.0])
    lith = np.array([1.5, 2.5, 3.5])
    u = blend.potential(lith, q)
    np.testing.assert_allclose(blend.lithium(u, q), lith, atol=2e-3)


def test_capacity_model_anchor(demo):
    m = of.CapacityModel(demo.positive, demo.negative, v_full=demo.v_full)
    t = demo.truth["BOL"]
    p = [t["Q_pe"], t["Q_ne"], t["n_li"], t["R"]]
    st = m.state(p)
    rest = demo.positive(st["y100"]) - demo.negative(st["x100"])
    assert rest == pytest.approx(demo.v_full, abs=1e-6)
    assert m(p, np.array([0.0]))[0] == pytest.approx(demo.v_full - t["R"], abs=1e-6)
    assert st["n_li"] == pytest.approx(st["y100"] * st["Q_pe"] + st["x100"] * st["Q_ne"])


def test_capacity_model_free_x100(demo):
    m = of.CapacityModel(demo.positive, demo.negative)
    assert m.param_names[-1] == "x100"
    v = m([5.0, 5.5, 4.875, 0.0, 0.85], np.array([0.0, 1.0]))
    assert np.all(np.isfinite(v))
