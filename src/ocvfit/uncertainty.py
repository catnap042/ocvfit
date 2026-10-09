"""Uncertainty and identifiability checks for capacity-model fits.

Residuals of OCV fits are strongly correlated along the curve, so treating
every sample as independent gives intervals that are far too narrow. The
procedure implemented here is:

1. **Correlation-corrected profile interval** - fix one parameter on a grid,
   refit all others, and accept values whose chi-square increase is below
   ``3.84 * (1 + rho) / (1 - rho)`` (inflation capped at 100), where ``rho`` is
   the lag-1 autocorrelation of the best-fit residuals.
2. **Reference-curve bias budget** - perturb each half-cell curve with a
   random offset (default +/-3 mV) and a smooth shape error (default 2 mV),
   refit, and take the largest deviation of each quantity.
3. **Decision** - final interval = profile interval +/- bias budget; a
   quantity whose half-width exceeds a threshold (default 1 percentage point)
   is reported as *not identifiable* (interval only, no point value).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .fitting import FitResult
from .models import CapacityModel
from .ocp import OCPCurve

FloatArray = NDArray[np.float64]

__all__ = [
    "lag1_autocorrelation",
    "ProfileResult",
    "profile_interval",
    "perturb_curve",
    "bias_budget",
    "is_identifiable",
    "mode_interval",
]

#: Model parameter behind each degradation mode (``mode = 1 - aged / BOL``).
MODE_PARAMETERS = {"LLI": "n_li", "LAM_PE": "Q_pe", "LAM_NE": "Q_ne"}


def lag1_autocorrelation(residuals: ArrayLike) -> float:
    """Lag-1 autocorrelation coefficient of a residual sequence (nan ignored)."""
    r = np.asarray(residuals, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 3:
        return 0.0
    r = r - r.mean()
    den = float(np.sum(r * r))
    return float(np.sum(r[1:] * r[:-1]) / den) if den > 0 else 0.0


@dataclass
class ProfileResult:
    """Profile of one parameter."""

    name: str
    values: FloatArray
    delta_chi2: FloatArray
    threshold: float
    rho: float
    interval: tuple[float, float]


def profile_interval(
    fit: FitResult,
    param: str,
    values: ArrayLike | None = None,
    rel_step: float = 1e-3,
    max_steps: int = 200,
    max_inflation: float = 100.0,
    invalid_residual: float = 1.0,
) -> ProfileResult:
    """Correlation-corrected profile-likelihood interval for one parameter.

    Parameters
    ----------
    fit
        Best fit of a :class:`~ocvfit.models.CapacityModel`.
    param
        Parameter name (e.g. ``"n_li"``, ``"Q_pe"``, ``"Q_ne"``).
    values
        Optional explicit grid of fixed values. By default the profile walks
        outwards from the best value in relative steps of ``rel_step`` until
        the threshold is crossed on both sides (at most ``max_steps`` per
        side); the crossing points are linearly interpolated.
    """
    model = fit.model
    if not isinstance(model, CapacityModel):
        raise TypeError("profile_interval requires a CapacityModel fit")
    names = model.param_names
    j = names.index(param)
    best = fit.params.copy()
    q, v = fit.t, fit.voltage
    bounds = model.default_bounds(q)
    lo = np.array([b[0] for i, b in enumerate(bounds) if i != j])
    hi = np.array([b[1] for i, b in enumerate(bounds) if i != j])

    def resid_full(p: FloatArray) -> FloatArray:
        r = model(p, q) - v
        return np.where(np.isfinite(r), r, invalid_residual)

    ssr_min = float(np.sum(resid_full(best) ** 2))
    dof = max(1, len(q) - len(best))
    sigma2 = ssr_min / dof
    rho = lag1_autocorrelation(fit.residuals)
    inflation = min(max_inflation, (1 + rho) / (1 - rho)) if rho < 1 else max_inflation
    threshold = 3.84 * max(1.0, inflation)

    start = np.delete(best, j)

    def profile_at(val: float) -> float:
        def resid(p_free: FloatArray) -> FloatArray:
            return resid_full(np.insert(p_free, j, val))

        sol = least_squares(resid, np.clip(start, lo + 1e-12, hi - 1e-12), bounds=(lo, hi), x_scale="jac")
        return float((2 * sol.cost - ssr_min) / sigma2)

    if values is not None:
        vals = np.asarray(values, dtype=float)
        dchi = np.array([profile_at(float(v_)) for v_ in vals])
        inside = vals[dchi <= threshold]
        interval = (float(inside.min()), float(inside.max())) if len(inside) else (float(best[j]),) * 2
        return ProfileResult(param, vals, dchi, threshold, rho, interval)

    base_step = rel_step * abs(best[j]) if best[j] != 0 else rel_step
    pts: dict[float, float] = {float(best[j]): 0.0}
    ends = []
    for direction in (-1.0, 1.0):
        # shrink the step until the first step stays below the threshold
        step = base_step
        for _ in range(12):
            d1 = profile_at(float(best[j] + direction * step))
            if d1 <= threshold:
                break
            pts[float(best[j] + direction * step)] = d1
            step /= 4.0
        prev_v, prev_d = float(best[j]), 0.0
        end = prev_v
        for k in range(1, max_steps + 1):
            val = float(best[j] + direction * k * step)
            if not (bounds[j][0] <= val <= bounds[j][1]):
                end = prev_v
                break
            d = profile_at(val)
            pts[val] = d
            if d > threshold:
                if k == 1:  # quadratic profile around the optimum
                    end = best[j] + (val - best[j]) * np.sqrt(threshold / d)
                else:
                    end = prev_v + (val - prev_v) * (threshold - prev_d) / (d - prev_d)
                break
            prev_v, prev_d, end = val, d, val
        ends.append(float(end))
    vals = np.array(sorted(pts))
    dchi = np.array([pts[v_] for v_ in vals])
    interval = (float(min(ends)), float(max(ends)))
    return ProfileResult(param, vals, dchi, threshold, rho, interval)


def perturb_curve(
    curve: OCPCurve, rng: np.random.Generator, offset: float = 3e-3, shape: float = 2e-3
) -> OCPCurve:
    """Random reference-curve error: uniform offset plus a smooth shape distortion."""
    dv = rng.uniform(-offset, offset)
    a = rng.uniform(-1, 1, 3)
    a = shape * a / max(1.0, np.max(np.abs(a)))
    lo, hi = curve.x_range

    def distortion(x: FloatArray) -> FloatArray:
        u = (x - lo) / (hi - lo)
        return a[0] * np.sin(np.pi * u) + a[1] * np.sin(2 * np.pi * u) + a[2] * np.sin(3 * np.pi * u)

    return curve.shifted(dv, distortion)


def bias_budget(
    refit: Callable[[Sequence[OCPCurve], OCPCurve], dict[str, float]],
    positive: Sequence[OCPCurve],
    negative: OCPCurve,
    reference: dict[str, float],
    n_draws: int = 8,
    offset: float = 3e-3,
    shape: float = 2e-3,
    seed: int = 0,
) -> dict[str, float]:
    """Largest deviation of each quantity under perturbed half-cell curves.

    Parameters
    ----------
    refit
        Callable ``refit(positive_curves, negative_curve) -> {name: value}``
        that rebuilds the model with the given curves, refits the data and
        returns the quantities of interest (e.g. degradation modes).
    positive, negative
        Nominal half-cell curves.
    reference
        Quantities obtained with the nominal curves.
    """
    rng = np.random.default_rng(seed)
    dev = {k: 0.0 for k in reference}
    for _ in range(n_draws):
        pos = [perturb_curve(c, rng, offset, shape) for c in positive]
        neg = perturb_curve(negative, rng, offset, shape)
        out = refit(pos, neg)
        for k in dev:
            if k in out and np.isfinite(out[k]):
                dev[k] = max(dev[k], abs(out[k] - reference[k]))
    return dev


def is_identifiable(half_width: float, threshold: float = 0.01) -> bool:
    """True if an interval half-width (fraction) is within ``threshold``."""
    return bool(np.isfinite(half_width) and half_width <= threshold)


def mode_interval(bol: FitResult, aged: FitResult, mode: str) -> tuple[float, float]:
    """Conservative profile interval of a degradation mode (fraction).

    ``mode`` is ``"LLI"``, ``"LAM_PE"`` or ``"LAM_NE"`` (single-material
    positive electrode). The profile intervals of the underlying parameter in
    the BOL and the aged fit are combined as
    ``[1 - aged_hi / bol_lo, 1 - aged_lo / bol_hi]``. Add the bias budget to
    obtain the final interval.
    """
    param = MODE_PARAMETERS[mode]
    b_lo, b_hi = profile_interval(bol, param).interval
    a_lo, a_hi = profile_interval(aged, param).interval
    return 1 - a_hi / b_lo, 1 - a_lo / b_hi
