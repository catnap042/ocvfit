import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

import ocvfit as of


def test_clean_curve_sorts_and_averages_duplicates():
    x, y = of.clean_curve([0.5, 0.1, 0.5, np.nan, 0.9], [2.0, 1.0, 4.0, 5.0, 3.0])
    np.testing.assert_allclose(x, [0.1, 0.5, 0.9])
    np.testing.assert_allclose(y, [1.0, 3.0, 3.0])


def test_monotonize_pava_decreasing():
    y = np.array([4.0, 3.9, 3.95, 3.8, 3.7])
    out = of.monotonize(y, decreasing=True)
    assert np.all(np.diff(out) <= 0)
    np.testing.assert_allclose(out[1:3], [3.925, 3.925])
    strict = of.monotonize(y, decreasing=True, min_step=1e-9)
    assert np.all(np.diff(strict) < 0)


def test_monotonize_keeps_monotone_input():
    y = np.linspace(4.2, 3.0, 50)
    np.testing.assert_allclose(of.monotonize(y), y)


def test_ocp_curve_never_extrapolates():
    c = of.OCPCurve(np.linspace(0.2, 0.8, 20), np.linspace(4.0, 3.5, 20))
    v = c([0.1, 0.5, 0.9])
    assert np.isnan(v[0]) and np.isnan(v[2]) and np.isfinite(v[1])


def test_from_soc_delithiation_axis():
    soc = np.linspace(0, 1, 11)
    volt = 3.5 + 0.7 * soc  # potential rises with delithiation
    c = of.OCPCurve.from_soc(soc, volt, coordinate="delithiation")
    # lithiation x = 1 - soc, so U decreases with x
    assert c(0.0) == pytest.approx(4.2)
    assert c(1.0) == pytest.approx(3.5)


def test_lithiation_inverse_round_trip():
    c = of.literature_ocp("nmc532_mohtat2020", 0.05, 0.95)
    x = np.linspace(0.1, 0.9, 25)
    np.testing.assert_allclose(c.lithiation(c(x)), x, atol=2e-4)
    lo, hi = c.x_range
    assert c.lithiation(10.0) == pytest.approx(lo)
    assert c.lithiation(-10.0) == pytest.approx(hi)


def test_compress_pchip_meets_tolerance():
    x = np.linspace(0, 1, 3000)
    v = of.LITERATURE_OCP["graphite_chen2020"](x)
    knots, rep = of.compress_pchip(x, v, tol=1e-3)
    assert rep["max_error"] <= 1e-3
    assert rep["n_knots"] < 100
    rebuilt = PchipInterpolator(knots[:, 0], knots[:, 1])(x)
    assert np.max(np.abs(rebuilt - v)) <= 1e-3 + 1e-12
    curve = of.OCPCurve.from_knots(knots)
    np.testing.assert_allclose(curve(knots[:, 0]), knots[:, 1])


def test_compress_pchip_needs_two_points():
    with pytest.raises(ValueError):
        of.compress_pchip([0.5, 0.5], [1.0, 2.0])
