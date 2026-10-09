"""Incremental capacity (dQ/dV) and differential voltage (dV/dQ) analysis.

Derivatives are computed with a Savitzky-Golay filter on a *uniform* grid
(the data are first resampled by linear interpolation), which smooths and
differentiates in one step. Typical settings: 2-5 mV voltage step, odd window
of 9-15 samples, polynomial order 2-3. Over-smoothing broadens peaks, shifts
their position or removes them; validate settings on reference curves.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import find_peaks, savgol_filter

from .ocp import clean_curve

FloatArray = NDArray[np.float64]

__all__ = ["resample_uniform", "incremental_capacity", "differential_voltage", "find_ic_peaks"]


def resample_uniform(x: ArrayLike, y: ArrayLike, step: float) -> tuple[FloatArray, FloatArray]:
    """Linearly resample ``y(x)`` on a uniform grid of spacing ``step``."""
    xs, ys = clean_curve(x, y)
    n = int(np.floor((xs[-1] - xs[0]) / step)) + 1
    xg = xs[0] + step * np.arange(n)
    return xg, np.interp(xg, xs, ys)


def incremental_capacity(
    voltage: ArrayLike,
    capacity: ArrayLike,
    dv: float = 2e-3,
    window: int = 11,
    polyorder: int = 3,
    clip_percentile: float | None = None,
) -> tuple[FloatArray, FloatArray]:
    """dQ/dV on a uniform voltage grid.

    Parameters
    ----------
    voltage, capacity
        One monotone branch (charge *or* discharge) of a low-rate curve.
    dv
        Voltage step of the resampling grid (V).
    window, polyorder
        Savitzky-Golay window (odd, in samples) and polynomial order.
    clip_percentile
        If set (e.g. 99), values whose magnitude exceeds this percentile are
        replaced by ``nan`` to suppress spikes at flat plateaus.

    Returns
    -------
    V, dQdV
        Voltage grid and dQ/dV (Ah/V). The sign follows the data: a
        discharge curve with capacity counted positive gives negative values.
    """
    vg, qg = resample_uniform(voltage, capacity, dv)
    if len(vg) < window:
        raise ValueError("not enough points for the requested window")
    d = savgol_filter(qg, window, polyorder, deriv=1, delta=dv)
    if clip_percentile is not None:
        lim = np.nanpercentile(np.abs(d), clip_percentile)
        d = np.where(np.abs(d) > lim, np.nan, d)
    return vg, d


def differential_voltage(
    capacity: ArrayLike,
    voltage: ArrayLike,
    dq: float | None = None,
    window: int = 11,
    polyorder: int = 3,
) -> tuple[FloatArray, FloatArray]:
    """dV/dQ on a uniform capacity grid (default step: 1/500 of the span)."""
    q = np.asarray(capacity, dtype=float)
    if dq is None:
        dq = float(np.nanmax(q) - np.nanmin(q)) / 500.0
    qg, vg = resample_uniform(q, voltage, dq)
    if len(qg) < window:
        raise ValueError("not enough points for the requested window")
    return qg, savgol_filter(vg, window, polyorder, deriv=1, delta=dq)


def find_ic_peaks(
    voltage: ArrayLike,
    dqdv: ArrayLike,
    height_ratio: float = 0.1,
    prominence: float | None = None,
) -> dict[str, FloatArray]:
    """Locate peaks of ``|dQ/dV|``.

    Peaks lower than ``height_ratio`` times the maximum are ignored.
    Returns peak ``voltage``, ``height`` and ``prominence`` arrays.
    """
    v = np.asarray(voltage, dtype=float)
    y = np.abs(np.asarray(dqdv, dtype=float))
    ok = np.isfinite(y)
    v, y = v[ok], y[ok]
    idx, props = find_peaks(y, height=height_ratio * np.max(y), prominence=prominence or 0.0)
    return {"voltage": v[idx], "height": y[idx], "prominence": props["prominences"]}
