"""Half-cell open-circuit potential (OCP) curves.

All curves in this package are expressed against the *lithiation fraction*
(stoichiometry) of the electrode material, ``x`` in ``[0, 1]``, with the
potential in volts versus Li/Li+. This is the convention used by most
physics-based battery models: the potential of both the positive and the
negative electrode decreases as the lithiation fraction increases.

Measured half-cell data are often stored against a state of charge (SOC) that
increases with *delithiation* instead. Use :meth:`OCPCurve.from_soc` with
``coordinate="delithiation"`` to convert such data.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d

FloatArray = NDArray[np.float64]

__all__ = [
    "OCPCurve",
    "clean_curve",
    "smooth",
    "monotonize",
    "compress_pchip",
]


def clean_curve(x: ArrayLike, y: ArrayLike) -> tuple[FloatArray, FloatArray]:
    """Remove non-finite points, sort by ``x`` and average duplicate ``x`` values.

    Parameters
    ----------
    x, y
        Abscissa and ordinate of a sampled curve.

    Returns
    -------
    tuple of ndarray
        Strictly increasing ``x`` and the matching ``y``. When several samples
        share the same ``x``, their ``y`` values are averaged.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.shape != y.shape:
        raise ValueError("x and y must have the same length")
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    ux, inv = np.unique(x, return_inverse=True)
    uy = np.bincount(inv, weights=y) / np.bincount(inv)
    return ux, uy


def smooth(y: ArrayLike, window: int = 5) -> FloatArray:
    """Moving-average smoothing (edge values are repeated at the boundaries).

    Returns the input unchanged when it is shorter than ``window`` or when
    ``window <= 1``.
    """
    y = np.asarray(y, dtype=float)
    if window <= 1 or len(y) < window:
        return y.copy()
    return uniform_filter1d(y, size=window, mode="nearest")


def monotonize(
    y: ArrayLike,
    weights: ArrayLike | None = None,
    decreasing: bool = True,
    min_step: float = 0.0,
) -> FloatArray:
    """Make a sequence monotone with the pool-adjacent-violators algorithm (PAVA).

    This is a weighted least-squares isotonic regression: neighbouring samples
    that violate the requested order are merged into their weighted mean. It is
    used to repair small non-monotone segments of measured OCP curves before
    the curve is inverted (potential -> lithiation fraction).

    Parameters
    ----------
    y
        Values ordered by increasing abscissa.
    weights
        Optional positive weights (default: equal weights).
    decreasing
        If True (default) return a non-increasing sequence, else non-decreasing.
    min_step
        Optional small step (in the units of ``y``) added to break exact ties so
        that the result is strictly monotone. Use e.g. ``1e-9`` V.
    """
    y = np.asarray(y, dtype=float)
    w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    sign = -1.0 if decreasing else 1.0
    z = sign * y  # solve the non-decreasing problem
    means: list[float] = []
    wts: list[float] = []
    counts: list[int] = []
    for zi, wi in zip(z, w):
        means.append(float(zi))
        wts.append(float(wi))
        counts.append(1)
        while len(means) > 1 and means[-2] > means[-1]:
            m2, w2, c2 = means.pop(), wts.pop(), counts.pop()
            m1, w1, c1 = means.pop(), wts.pop(), counts.pop()
            wt = w1 + w2
            means.append((m1 * w1 + m2 * w2) / wt)
            wts.append(wt)
            counts.append(c1 + c2)
    out = np.repeat(np.asarray(means), counts)
    if min_step > 0:
        # enforce strict monotonicity by tiny increments inside tied blocks
        for i in range(1, len(out)):
            if out[i] <= out[i - 1]:
                out[i] = out[i - 1] + min_step
    return sign * out


def compress_pchip(
    x: ArrayLike, y: ArrayLike, tol: float = 1.2e-3, max_knots: int | None = None
) -> tuple[FloatArray, dict[str, float]]:
    """Compress a densely sampled curve into a small set of PCHIP knots.

    Greedy adaptive knot insertion: start from the two end points, rebuild a
    shape-preserving PCHIP interpolant, add the sample with the largest
    absolute residual, and repeat until the maximum residual on all samples is
    below ``tol``.

    Parameters
    ----------
    x, y
        Sampled curve. Duplicate ``x`` values are averaged first.
    tol
        Maximum absolute reconstruction error, in the units of ``y``
        (default 1.2 mV for a curve in volts).
    max_knots
        Optional hard limit on the number of knots.

    Returns
    -------
    knots : ndarray, shape (n, 2)
        Retained ``(x, y)`` control points, sorted by ``x``. The curve is
        reconstructed with ``PchipInterpolator(knots[:, 0], knots[:, 1])``.
    report : dict
        ``n_input``, ``n_unique``, ``n_knots``, ``max_error`` and ``rmse``
        evaluated on the (deduplicated) input samples.
    """
    x_in = np.asarray(x, dtype=float)
    ux, uy = clean_curve(x, y)
    if len(ux) < 2:
        raise ValueError("at least two distinct x values are required")
    keep = {0, len(ux) - 1}
    limit = len(ux) if max_knots is None else max(2, int(max_knots))
    while True:
        idx = np.array(sorted(keep))
        f = PchipInterpolator(ux[idx], uy[idx], extrapolate=False)
        residual = np.abs(f(ux) - uy)
        worst = int(np.nanargmax(residual))
        if residual[worst] <= tol or len(keep) >= limit:
            break
        keep.add(worst)
    knots = np.column_stack([ux[idx], uy[idx]])
    report = {
        "n_input": float(x_in.size),
        "n_unique": float(len(ux)),
        "n_knots": float(len(idx)),
        "max_error": float(np.nanmax(residual)),
        "rmse": float(np.sqrt(np.nanmean(residual**2))),
    }
    return knots, report


@dataclass
class OCPCurve:
    """A half-cell OCP curve ``U(x)`` interpolated with PCHIP.

    The interpolant is shape preserving (no spurious oscillations between
    plateaus and steep regions) and is **never extrapolated**: evaluating the
    curve outside its data range returns ``nan``.

    Parameters
    ----------
    x
        Lithiation fraction of the electrode material (any order; duplicates
        are averaged).
    voltage
        Potential in V vs Li/Li+.
    name
        Optional label used in plots and reports.
    smooth_window
        Optional moving-average window applied to the voltage before
        interpolation (``1`` = no smoothing). Useful for noisy measured curves.
    """

    x: FloatArray
    voltage: FloatArray
    name: str = "electrode"
    smooth_window: int = 1
    _interp: PchipInterpolator = field(init=False, repr=False)

    def __post_init__(self) -> None:
        x, v = clean_curve(self.x, self.voltage)
        if len(x) < 2:
            raise ValueError("an OCP curve needs at least two distinct points")
        v = smooth(v, self.smooth_window)
        self.x, self.voltage = x, v
        self._interp = PchipInterpolator(x, v, extrapolate=False)
        self._inverse_cache: tuple[FloatArray, FloatArray] | None = None

    # ------------------------------------------------------------------ build
    @classmethod
    def from_function(
        cls,
        func: Callable[[FloatArray], FloatArray],
        x_min: float = 0.0,
        x_max: float = 1.0,
        n: int = 1001,
        name: str = "electrode",
    ) -> OCPCurve:
        """Sample an analytic OCP function ``func(x)`` on ``[x_min, x_max]``."""
        x = np.linspace(x_min, x_max, n)
        return cls(x, np.asarray(func(x), dtype=float), name=name)

    @classmethod
    def from_soc(
        cls,
        soc: ArrayLike,
        voltage: ArrayLike,
        coordinate: Literal["lithiation", "delithiation"] = "lithiation",
        name: str = "electrode",
        smooth_window: int = 1,
    ) -> OCPCurve:
        """Build a curve from half-cell data stored against an SOC axis.

        ``coordinate="delithiation"`` means that the SOC grows as lithium is
        removed (typical for a positive-electrode charge curve); the axis is
        converted with ``x = 1 - soc``.
        """
        s = np.asarray(soc, dtype=float)
        x = 1.0 - s if coordinate == "delithiation" else s
        return cls(x, np.asarray(voltage, dtype=float), name=name, smooth_window=smooth_window)

    @classmethod
    def from_knots(cls, knots: ArrayLike, name: str = "electrode") -> OCPCurve:
        """Rebuild a curve from ``(x, U)`` knots returned by :func:`compress_pchip`."""
        k = np.asarray(knots, dtype=float)
        return cls(k[:, 0], k[:, 1], name=name)

    # --------------------------------------------------------------- evaluate
    @property
    def x_range(self) -> tuple[float, float]:
        """Lithiation-fraction range covered by the data."""
        return float(self.x[0]), float(self.x[-1])

    @property
    def voltage_range(self) -> tuple[float, float]:
        """(min, max) potential covered by the data."""
        return float(np.min(self.voltage)), float(np.max(self.voltage))

    def __call__(self, x: ArrayLike) -> FloatArray:
        """Potential at lithiation fraction ``x`` (``nan`` outside the data range)."""
        return np.asarray(self._interp(np.asarray(x, dtype=float)), dtype=float)

    def derivative(self, x: ArrayLike) -> FloatArray:
        """dU/dx of the interpolant (``nan`` outside the data range)."""
        return np.asarray(self._interp.derivative()(np.asarray(x, dtype=float)), dtype=float)

    def lithiation(self, potential: ArrayLike, n_grid: int = 4001) -> FloatArray:
        """Inverse curve ``x(U)``, saturated at the ends of the data range.

        A potential above the highest measured value returns the lowest
        lithiation fraction and vice versa ("the phase is exhausted"); the OCP
        itself is never extrapolated. Non-monotone segments are repaired with
        :func:`monotonize` before inversion.
        """
        if self._inverse_cache is None:
            xg = np.linspace(self.x[0], self.x[-1], n_grid)
            vg = monotonize(self(xg), decreasing=True)
            # strictly decreasing for np.interp
            vg = vg - np.arange(n_grid) * 1e-12
            self._inverse_cache = (vg[::-1].copy(), xg[::-1].copy())
        vp, xp = self._inverse_cache
        return np.interp(np.asarray(potential, dtype=float), vp, xp)

    def shifted(self, dv: float = 0.0, shape: Callable[[FloatArray], FloatArray] | None = None) -> OCPCurve:
        """Return a copy with a constant offset ``dv`` and an optional additive shape error."""
        v = self.voltage + dv
        if shape is not None:
            v = v + np.asarray(shape(self.x), dtype=float)
        return OCPCurve(self.x.copy(), v, name=self.name)
