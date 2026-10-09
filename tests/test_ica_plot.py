import numpy as np
import pytest

import ocvfit as of


def test_incremental_capacity_linear():
    v = np.linspace(3.0, 4.0, 400)
    q = 2.5 * (4.0 - v)  # discharge: capacity grows as voltage falls
    vg, d = of.incremental_capacity(v, q, dv=2e-3)
    np.testing.assert_allclose(d, -2.5, atol=1e-8)
    assert vg[0] == pytest.approx(3.0)


def test_differential_voltage_linear():
    q = np.linspace(0, 5, 300)
    v = 4.2 - 0.2 * q
    qg, d = of.differential_voltage(q, v)
    np.testing.assert_allclose(d, -0.2, atol=1e-8)


def test_find_ic_peaks():
    v = np.linspace(3.0, 4.0, 1001)
    y = np.exp(-(((v - 3.4) / 0.02) ** 2)) + 0.5 * np.exp(-(((v - 3.8) / 0.03) ** 2))
    peaks = of.find_ic_peaks(v, y, height_ratio=0.2)
    np.testing.assert_allclose(peaks["voltage"], [3.4, 3.8], atol=2e-3)


def test_plot_report(demo_fits, tmp_path):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    from ocvfit.plotting import plot_fit_report

    fig = plot_fit_report(demo_fits, title="test")
    out = tmp_path / "fig.png"
    fig.savefig(out, dpi=60)
    assert out.stat().st_size > 1000
