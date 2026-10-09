"""Parameter-recovery tests on synthetic data generated from known parameters."""

import numpy as np
import pytest

import ocvfit as of


def assert_close(got, want, atol):
    err = np.abs(np.asarray(got) - np.asarray(want))
    assert np.all(err <= np.asarray(atol)), f"got {got}, want {want}, err {err}"


@pytest.mark.parametrize("which", ["lfp", "nmc532"])
def test_capacity_model_recovers_parameters(which, request):
    cell = request.getfixturevalue("demo" if which == "lfp" else "demo_nmc")
    fits = request.getfixturevalue("demo_fits" if which == "lfp" else "demo_nmc_fits")
    for label, fit in fits.items():
        truth = cell.truth[label]
        got = fit.param_dict
        for k in ("Q_pe", "Q_ne", "n_li"):
            assert got[k] == pytest.approx(truth[k], rel=0.01), (which, label, k)
        assert got["R"] == pytest.approx(truth["R"], abs=2e-3)
        assert fit.rmse < 1.5e-3  # noise level is 1 mV
        assert fit.valid_fraction == 1.0


@pytest.mark.parametrize("which", ["lfp", "nmc532"])
def test_degradation_modes_recovered(which, request):
    cell = request.getfixturevalue("demo" if which == "lfp" else "demo_nmc")
    fits = request.getfixturevalue("demo_fits" if which == "lfp" else "demo_nmc_fits")
    modes = of.degradation_modes(fits["BOL"], fits["aged"])
    true = cell.true_modes
    # LLI and LAM_NE: within 0.5 percentage point for both chemistries
    for k in ("LLI", "LAM_NE"):
        assert modes[k] == pytest.approx(true[k], abs=0.005), (which, k)
    # LAM_PE: noise-only recovery within 1 pp (see the LFP identifiability tests)
    assert modes["LAM_PE"] == pytest.approx(true["LAM_PE"], abs=0.01), which


@pytest.mark.parametrize("loss", ["v2", "rmse", "mse"])
def test_cathode_fixed_recovery(loss):
    pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
    neg = of.literature_ocp("graphite_chen2020")
    m = of.CathodeFixedModel(pos, neg)
    s = np.linspace(0, 1, 120)
    true = np.array([1.1, 0.05, 0.01])
    v = m(true, s) + np.random.default_rng(1).normal(0, 1e-3, s.size)
    fit = of.fit_ocv(m, s, v, loss=loss)
    assert_close(fit.params, true, [0.02, 0.003, 0.002])
    assert fit.rmse < 1.5e-3


def test_cathode_fixed_differential_evolution():
    pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
    neg = of.literature_ocp("graphite_chen2020")
    m = of.CathodeFixedModel(pos, neg)
    s = np.linspace(0, 1, 80)
    true = np.array([1.15, 0.08, -0.005])
    fit = of.fit_ocv(m, s, m(true, s), method="differential_evolution")
    assert_close(fit.params, true, [5e-3, 1e-3, 1e-3])


def test_window_model_recovery():
    m = of.WindowModel(of.literature_ocp("nmc811_chen2020"), of.literature_ocp("graphite_chen2020"))
    s = np.linspace(0, 1, 150)
    true = np.array([0.90, 0.27, 0.03, 0.91, 0.0])
    v = m(true, s) + np.random.default_rng(2).normal(0, 1e-3, s.size)
    fit = of.fit_ocv(m, s, v, loss="mse")
    assert_close(fit.params, true, [0.01, 0.01, 0.01, 0.01, 3e-3])


def test_blended_capacity_model_recovery():
    nmc = of.literature_ocp("nmc811_chen2020")
    nmc.name = "NMC811"
    lfp = of.literature_ocp("lfp_afshar2017")
    lfp.name = "LFP"
    m = of.CapacityModel([nmc, lfp], of.literature_ocp("graphite_chen2020"), v_full=4.2)
    assert m.param_names == ["Q_NMC811", "Q_LFP", "Q_ne", "n_li", "R"]
    true = np.array([3.0, 2.0, 5.5, 4.9, 0.005])
    q, v = of.synthesize_discharge(m, true, v_min=2.7, n=250, noise=1e-3, seed=3)
    fit = of.fit_capacity(m, q, v)
    np.testing.assert_allclose(fit.params[:4], true[:4], rtol=0.03)


def test_fit_result_reporting(demo_fits):
    fit = demo_fits["BOL"]
    d = fit.to_dict()
    assert d["model"] == "CapacityModel"
    assert {"x100", "y100", "x0", "y0", "np_ratio"} <= set(d["state"])
    assert "RMSE" in fit.summary()


def test_objective_penalises_shrinking_valid_region():
    pos = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
    neg = of.literature_ocp("graphite_chen2020")
    m = of.CathodeFixedModel(pos, neg)
    s = np.linspace(0, 1, 50)
    v = m([1.1, 0.05, 0.0], s)
    full = of.objective([1.1, 0.05, 0.0], m, s, v)
    shrunk = of.objective([0.9, 0.05, 0.0], m, s, v)  # loses the high-SOC end
    assert shrunk > 1.0 > full


def test_huber_and_weights():
    r = np.array([0.0, 1e-3, 1e-2])
    h = of.huber_loss(r, delta=2e-3)
    assert h[1] == pytest.approx(0.5e-6)
    assert h[2] == pytest.approx(2e-3 * (1e-2 - 1e-3))
    w = of.build_weights(np.linspace(0, 1, 101), edge_frac=0.05, edge_weight=10.0)
    assert w[0] == 10 and w[50] == 1 and w[-1] == 10
