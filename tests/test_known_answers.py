"""Known-answer tests whose expected values do not come from ocvfit itself.

Why the answers are independent of the code under test:

* The half-cell OCP functions below are small closed-form toys written only for
  these tests. They are not the literature functions in ``ocvfit.literature``.
* The full-cell curves are generated here with NumPy and SciPy, straight from
  the electrode-balance equations: x = x100 - q/Q_ne, y = (n_Li - x100*Q_ne + q)/Q_pe,
  V = U_pe(y) - U_ne(x) - R, with x100 solved by our own root find from v_full.
  They do not go through ``CapacityModel``, ``synthesize_discharge``, the PCHIP
  interpolation in ``OCPCurve`` or any other ocvfit code.
* The expected degradation modes are the percentages we set by hand
  (LLI 5 %, LAM_PE 3 %, LAM_NE 4 %), not values taken from an earlier fit.

ocvfit only sees the sampled OCP tables and the noisy voltage curves, so a bug in
its model would show up as a recovery error rather than cancel out.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

import ocvfit as of

TRUE_MODES = {"LLI": 0.05, "LAM_PE": 0.03, "LAM_NE": 0.04}
NOISE_V = 1e-3  # 1 mV Gaussian voltage noise
SEEDS = (0, 1, 2)


def u_sloped_pe(y):
    """Toy layered-oxide positive OCP: a sloped curve with one weak step."""
    y = np.asarray(y, dtype=float)
    return 4.25 - 0.9 * y - 0.08 * np.tanh((y - 0.55) / 0.06) - 0.35 * np.exp(-(1 - y) / 0.03)


def u_flat_pe(y):
    """Toy LFP-like positive OCP: almost flat, with steep ends only near y=0 and y=1."""
    y = np.asarray(y, dtype=float)
    return 3.42 - 0.01 * (y - 0.5) + 0.5 * np.exp(-y / 0.02) - 0.6 * np.exp(-(1 - y) / 0.02)


def u_ne(x):
    """Toy graphite-like negative OCP with two staging steps."""
    x = np.asarray(x, dtype=float)
    return (
        0.08
        + 0.6 * np.exp(-x / 0.03)
        + 0.05 * (1 - np.tanh((x - 0.2) / 0.03))
        + 0.04 * (1 - np.tanh((x - 0.55) / 0.03))
    )


def full_cell(u_pe, q_pe, q_ne, n_li, r, v_full, depth, seed, n=300):
    """Discharge curve built from the balance equations (no ocvfit code)."""
    x100 = brentq(lambda x: float(u_pe((n_li - x * q_ne) / q_pe) - u_ne(x)) - v_full, 0.3, 0.999)
    y100 = (n_li - x100 * q_ne) / q_pe
    q_max = depth * min((x100 - 0.01) * q_ne, (0.99 - y100) * q_pe)
    q = np.linspace(0.0, q_max, n)
    v = u_pe(y100 + q / q_pe) - u_ne(x100 - q / q_ne) - r
    return q, v + np.random.default_rng(seed).normal(0.0, NOISE_V, n)


def aged_params(bol):
    q_pe, q_ne, n_li, r = bol
    return (
        q_pe * (1 - TRUE_MODES["LAM_PE"]),
        q_ne * (1 - TRUE_MODES["LAM_NE"]),
        n_li * (1 - TRUE_MODES["LLI"]),
        r + 0.003,
    )


def fit_pair(u_pe, bol, v_full, seed, depth):
    pos = of.OCPCurve.from_function(u_pe, 0.0, 1.0, 2001, name="pe")
    neg = of.OCPCurve.from_function(u_ne, 0.0, 1.0, 2001, name="ne")
    model = of.CapacityModel(pos, neg, v_full=v_full)
    fb = of.fit_capacity(model, *full_cell(u_pe, *bol, v_full, depth, seed))
    fa = of.fit_capacity(model, *full_cell(u_pe, *aged_params(bol), v_full, depth, seed + 100))
    return fb, fa


@pytest.mark.parametrize("seed", SEEDS)
def test_sloped_cell_recovers_known_modes(seed):
    """All three modes are visible in the curve, so each must be recovered.

    Tolerance: 0.5 pp per mode (observed error <= 0.15 pp at 1 mV noise), and
    the profile interval of each mode, widened by 0.2 pp, must contain the truth.
    """
    bol = (3.0, 3.3, 2.7, 0.005)  # Q_pe, Q_ne, n_Li (Ah), R (V)
    fb, fa = fit_pair(u_sloped_pe, bol, v_full=4.10, seed=seed, depth=0.97)
    assert fb.rmse < 2 * NOISE_V and fa.rmse < 2 * NOISE_V
    modes = of.degradation_modes(fb, fa)
    for name, truth in TRUE_MODES.items():
        assert modes[name] == pytest.approx(truth, abs=0.005), name
        lo, hi = of.mode_interval(fb, fa, name)
        assert lo - 0.002 <= truth <= hi + 0.002, name
        assert of.is_identifiable(0.5 * (hi - lo)), name
    # Absolute capacities too, not just ratios (1 % relative)
    got = fb.param_dict
    assert got["Q_pe"] == pytest.approx(bol[0], rel=0.01)
    assert got["Q_ne"] == pytest.approx(bol[1], rel=0.01)
    assert got["n_li"] == pytest.approx(bol[2], rel=0.01)
    # Aged absolute capacities against the hand-set aged truth
    aged = aged_params(bol)
    got_a = fa.param_dict
    for k, i in (("Q_pe", 0), ("Q_ne", 1), ("n_li", 2)):
        assert got_a[k] == pytest.approx(aged[i], rel=0.01), k


@pytest.mark.parametrize("seed", SEEDS)
def test_flat_positive_gate_fires(seed):
    """Gate test: this case is built so that LAM_PE cannot be identified.

    The positive OCP is flat across the whole measured window (the cell is
    negative-limited, so its steep lithiation end is never reached). The data
    carry no LAM_PE information, so the identifiability check must return
    False. LLI and LAM_NE stay identifiable through the negative staging steps.
    If this test starts passing with a narrow LAM_PE interval, the uncertainty
    code has become over-confident.
    """
    bol = (3.6, 3.0, 2.8, 0.005)
    fb, fa = fit_pair(u_flat_pe, bol, v_full=3.45, seed=seed, depth=0.97)
    lo, hi = of.mode_interval(fb, fa, "LAM_PE")
    assert not of.is_identifiable(0.5 * (hi - lo))
    assert lo <= TRUE_MODES["LAM_PE"] <= hi  # wide, but still honest
    for name in ("LLI", "LAM_NE"):
        lo, hi = of.mode_interval(fb, fa, name)
        assert of.is_identifiable(0.5 * (hi - lo)), name
