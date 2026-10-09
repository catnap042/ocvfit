"""Parameter estimation for the full-cell OCV models.

Two fitting front-ends are provided:

* :func:`fit_ocv` - for SOC-axis models (:class:`~ocvfit.models.CathodeFixedModel`
  and :class:`~ocvfit.models.WindowModel`). It minimises a scalar objective
  that is designed to keep the curve ends well fitted:

  - ``"v2"`` (default): Huber loss, edge and slope-aware weights, an end-point
    anchor penalty on the mean residual of each end region, and a penalty that
    prevents the optimiser from shrinking the valid region to escape errors;
  - ``"rmse"``: weighted RMSE with fixed end-region weights plus the same
    valid-region penalty;
  - ``"mse"``: plain mean squared error plus the valid-region penalty.

* :func:`fit_capacity` - for :class:`~ocvfit.models.CapacityModel`. Starting
  points are pre-screened on a grid of capacities, then a two-stage
  least-squares fit is run from the best candidates: first with a robust
  ``soft_l1`` loss (scale 5 mV) that is insensitive to the few large residuals
  of the steep end of discharge, then with ordinary least squares.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import differential_evolution, least_squares, minimize

from .models import CapacityModel, CathodeFixedModel, OCVModel

FloatArray = NDArray[np.float64]

__all__ = [
    "FitResult",
    "huber_loss",
    "build_weights",
    "objective",
    "fit_ocv",
    "fit_capacity",
]


@dataclass
class FitResult:
    """Outcome of a fit.

    ``rmse``, ``mae`` and ``max_abs_error`` are always the *unweighted*
    statistics of ``model - data`` over the valid points (V), independent of
    the loss function that was minimised (``loss``).
    """

    model: OCVModel
    params: FloatArray
    t: FloatArray
    voltage: FloatArray
    loss: float
    loss_type: str
    success: bool
    message: str
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def param_dict(self) -> dict[str, float]:
        return self.model.as_dict(self.params)

    @property
    def fitted(self) -> FloatArray:
        return self.model(self.params, self.t)

    @property
    def residuals(self) -> FloatArray:
        """Model minus data (V); ``nan`` where the model is invalid."""
        return self.fitted - self.voltage

    @property
    def valid_fraction(self) -> float:
        return float(np.mean(np.isfinite(self.residuals)))

    @property
    def rmse(self) -> float:
        r = self.residuals
        return float(np.sqrt(np.nanmean(r**2)))

    @property
    def mae(self) -> float:
        return float(np.nanmean(np.abs(self.residuals)))

    @property
    def max_abs_error(self) -> float:
        return float(np.nanmax(np.abs(self.residuals)))

    def summary(self) -> str:
        """Human readable one-block summary."""
        lines = [f"{type(self.model).__name__} fit ({self.loss_type} loss, success={self.success})"]
        for k, v in self.param_dict.items():
            lines.append(f"  {k:>12s} = {v: .6g}")
        if isinstance(self.model, CapacityModel):
            st = self.model.state(self.params, q_end=float(np.nanmax(self.t)))
            for k in ("x100", "y100", "x0", "y0", "np_ratio"):
                lines.append(f"  {k:>12s} = {st[k]: .6g}")
        lines.append(
            f"  RMSE = {1e3 * self.rmse:.3f} mV, MAE = {1e3 * self.mae:.3f} mV, "
            f"max |err| = {1e3 * self.max_abs_error:.3f} mV, valid = {100 * self.valid_fraction:.1f} %"
        )
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "model": type(self.model).__name__,
            "params": self.param_dict,
            "loss_type": self.loss_type,
            "loss": self.loss,
            "rmse_V": self.rmse,
            "mae_V": self.mae,
            "max_abs_error_V": self.max_abs_error,
            "valid_fraction": self.valid_fraction,
            "success": self.success,
            "message": self.message,
        }
        if isinstance(self.model, CapacityModel):
            out["state"] = self.model.state(self.params, q_end=float(np.nanmax(self.t)))
        return out


# ---------------------------------------------------------------- objectives
def huber_loss(r: ArrayLike, delta: float = 2e-3) -> FloatArray:
    """Huber loss: quadratic for ``|r| <= delta``, linear beyond."""
    r = np.asarray(r, dtype=float)
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * r**2, delta * (a - 0.5 * delta))


def _normalised(t: FloatArray) -> FloatArray:
    lo, hi = np.nanmin(t), np.nanmax(t)
    return (t - lo) / (hi - lo) if hi > lo else np.zeros_like(t)


def build_weights(
    t: ArrayLike,
    voltage: ArrayLike | None = None,
    slope_weight: float = 0.5,
    edge_frac: float = 0.05,
    edge_weight: float = 10.0,
) -> FloatArray:
    """Point weights emphasising the curve ends and information-rich regions.

    Points within ``edge_frac`` of either end of the (normalised) abscissa get
    ``edge_weight``. If ``voltage`` is given, weights are further multiplied by
    ``1 + slope_weight * |dV/dt| / median(|dV/dt|)``.
    """
    s = _normalised(np.asarray(t, dtype=float))
    w = np.ones_like(s)
    w[s <= edge_frac] *= edge_weight
    w[s >= 1.0 - edge_frac] *= edge_weight
    if voltage is not None and slope_weight > 0 and len(s) > 2:
        dv = np.abs(np.gradient(np.asarray(voltage, dtype=float), s))
        dv = dv / (np.nanmedian(dv) + 1e-12)
        w = w * (1.0 + slope_weight * np.nan_to_num(dv))
    return w


def objective(
    params: ArrayLike,
    model: OCVModel,
    t: FloatArray,
    voltage: FloatArray,
    loss: Literal["v2", "rmse", "mse"] = "v2",
    use_weights: bool = True,
    huber_delta: float = 2e-3,
    slope_weight: float = 0.5,
    edge_frac: float = 0.05,
    edge_weight: float = 10.0,
    lambda_edge: float = 50.0,
    min_valid_ratio: float = 1.0,
    lambda_valid: float = 1e3,
    anchor: Literal["both", "low", "high", "none"] = "both",
) -> float:
    """Scalar objective used by :func:`fit_ocv` (see module docstring).

    ``min_valid_ratio`` is the fraction of data points that must lie inside
    the model's valid region; each missing fraction below it costs
    ``lambda_valid``. The default (1.0) requires full coverage so that the
    optimiser cannot discard poorly fitted end points by shrinking the valid
    region. Lower it (e.g. 0.85) only if a reference curve genuinely does not
    cover the full-cell window.
    """
    v_model = model(params, t)
    valid = np.isfinite(v_model)
    ratio = valid.mean() if len(valid) else 0.0
    penalty = lambda_valid * (min_valid_ratio - ratio) if ratio < min_valid_ratio else 0.0
    min_points = 5 if loss == "rmse" else 8
    if valid.sum() < min_points:
        return 1e6 + penalty
    tv, r = t[valid], v_model[valid] - voltage[valid]
    s = _normalised(t)[valid]

    if loss == "mse":
        return float(np.mean(r**2) + penalty)

    if loss == "rmse":
        if use_weights:
            w = np.ones_like(s)
            w[s < 0.05] = 10.0
            w[(s >= 0.05) & (s < 0.15)] = 5.0
            w[s > 0.90] = 3.0
        else:
            w = np.ones_like(s)
        return float(np.sqrt(np.average(r**2, weights=w)) + penalty)

    if loss != "v2":
        raise ValueError("loss must be 'v2', 'rmse' or 'mse'")
    w = (
        build_weights(tv, voltage[valid], slope_weight, edge_frac, edge_weight)
        if use_weights
        else np.ones_like(s)
    )
    main = np.average(huber_loss(r, huber_delta), weights=w)
    edge = 0.0
    lo, hi = s <= edge_frac, s >= 1.0 - edge_frac
    if anchor in ("both", "low") and lo.any():
        edge += np.mean(r[lo]) ** 2
    if anchor in ("both", "high") and hi.any():
        edge += np.mean(r[hi]) ** 2
    return float(main + lambda_edge * edge + penalty)


def _start_points(
    model: OCVModel, bounds: Sequence[tuple[float, float]], n_random: int, seed: int
) -> list[FloatArray]:
    if isinstance(model, CathodeFixedModel):
        grid = itertools.product([0.9, 1.0, 1.1, 1.2], [0.0, 0.05, 0.1, 0.15], [-0.02, 0.0, 0.02])
        return [np.array(g, dtype=float) for g in grid]
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    pts = [0.5 * (lo + hi)]
    pts += list(lo + (hi - lo) * rng.random((n_random, len(bounds))))
    return pts


def fit_ocv(
    model: OCVModel,
    t: ArrayLike,
    voltage: ArrayLike,
    loss: Literal["v2", "rmse", "mse"] = "v2",
    method: Literal["fast", "differential_evolution", "local"] = "fast",
    bounds: Sequence[tuple[float, float]] | None = None,
    x0: ArrayLike | None = None,
    seed: int = 42,
    n_random: int = 300,
    **loss_kwargs: Any,
) -> FitResult:
    """Fit an SOC-axis OCV model to a measured full-cell curve.

    Parameters
    ----------
    model
        A :class:`~ocvfit.models.CathodeFixedModel` or
        :class:`~ocvfit.models.WindowModel` (any :class:`OCVModel` works).
    t, voltage
        Full-cell SOC (0-1) and voltage (V). Non-finite samples are dropped.
    loss
        ``"v2"`` (recommended), ``"rmse"`` or ``"mse"``.
    method
        ``"fast"``: pre-screen start points (a fixed grid for the
        cathode-fixed model, random samples within the bounds otherwise), then
        L-BFGS-B followed by a Nelder-Mead polish. ``"differential_evolution"``:
        global search with polishing. ``"local"``: L-BFGS-B from ``x0``.
    bounds
        Parameter bounds (default: ``model.default_bounds(t)``).
    x0
        Initial guess for ``method="local"``.
    seed
        Random seed (pre-screening and differential evolution).
    **loss_kwargs
        Forwarded to :func:`objective` (e.g. ``lambda_edge``, ``huber_delta``,
        ``edge_frac``, ``slope_weight``, ``min_valid_ratio``).
    """
    t = np.asarray(t, dtype=float)
    voltage = np.asarray(voltage, dtype=float)
    if t.shape != voltage.shape:
        raise ValueError("t and voltage must have the same shape")
    keep = np.isfinite(t) & np.isfinite(voltage)
    t, voltage = t[keep], voltage[keep]
    bounds = list(bounds) if bounds is not None else model.default_bounds(t)

    def f(p: FloatArray) -> float:
        return objective(p, model, t, voltage, loss=loss, **loss_kwargs)

    if method == "differential_evolution":
        opt = differential_evolution(f, bounds, seed=seed, polish=True, maxiter=250, tol=1e-7)
    else:
        if method == "fast":
            starts = _start_points(model, bounds, n_random, seed)
            start = min(starts, key=f)
        elif method == "local":
            if x0 is None:
                lo = np.array([b[0] for b in bounds])
                hi = np.array([b[1] for b in bounds])
                start = 0.5 * (lo + hi)
            else:
                start = np.asarray(x0, dtype=float)
        else:
            raise ValueError("method must be 'fast', 'differential_evolution' or 'local'")
        opt = minimize(f, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 800, "ftol": 1e-12})
        if method == "fast":
            nm = minimize(
                f,
                opt.x,
                method="Nelder-Mead",
                bounds=bounds,
                options={"maxiter": 4000, "xatol": 1e-9, "fatol": 1e-14},
            )
            if nm.fun < opt.fun:
                opt = nm
    return FitResult(
        model=model,
        params=np.asarray(opt.x, dtype=float),
        t=t,
        voltage=voltage,
        loss=float(opt.fun),
        loss_type=loss,
        success=bool(opt.success),
        message=str(opt.message),
    )


# ---------------------------------------------------------- capacity model
def _capacity_candidates(model: CapacityModel, span: float) -> list[FloatArray]:
    q_tot = span * np.array([1.0, 1.1, 1.25, 1.4, 1.6])
    q_ne = span * np.array([1.0, 1.1, 1.25, 1.4, 1.6, 1.8])
    n_li = span * np.array([0.9, 1.0, 1.1, 1.2, 1.35, 1.5])
    k = model.n_positive
    if k == 1:
        splits: list[FloatArray] = [np.array([1.0])]
    else:
        fr = [0.2, 0.4, 0.6, 0.8]
        splits = []
        for combo in itertools.product(fr, repeat=k - 1):
            if sum(combo) < 1.0:
                splits.append(np.array([1.0 - sum(combo), *combo]))
    x100s = [None] if model.v_full is not None else [0.6, 0.75, 0.9]
    out = []
    for qt, qn, nl, sp, xh in itertools.product(q_tot, q_ne, n_li, splits, x100s):
        p = [*(qt * sp), qn, nl, 0.0]
        if xh is not None:
            p.append(xh)
        out.append(np.array(p, dtype=float))
    return out


def fit_capacity(
    model: CapacityModel,
    q: ArrayLike,
    voltage: ArrayLike,
    n_starts: int = 10,
    f_scale: float = 5e-3,
    bounds: Sequence[tuple[float, float]] | None = None,
    candidates: Sequence[ArrayLike] | None = None,
    invalid_residual: float = 1.0,
) -> FitResult:
    """Fit a :class:`~ocvfit.models.CapacityModel` to a low-rate discharge curve.

    Parameters
    ----------
    model
        The capacity model (single or blended positive electrode).
    q
        Capacity discharged from the fully charged state (Ah).
    voltage
        Cell voltage (V).
    n_starts
        Number of best grid candidates used as starting points.
    f_scale
        Scale of the robust ``soft_l1`` stage (V).
    bounds
        Parameter bounds (default: ``model.default_bounds(q)``).
    candidates
        Optional explicit list of start vectors (overrides the grid).
    invalid_residual
        Residual (V) assigned to points where the model is undefined.
    """
    q = np.asarray(q, dtype=float)
    voltage = np.asarray(voltage, dtype=float)
    keep = np.isfinite(q) & np.isfinite(voltage)
    q, voltage = q[keep], voltage[keep]
    bounds = list(bounds) if bounds is not None else model.default_bounds(q)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    span = float(np.nanmax(q) - np.nanmin(q))

    def resid(p: FloatArray) -> FloatArray:
        r = model(p, q) - voltage
        return np.where(np.isfinite(r), r, invalid_residual)

    cands = [np.asarray(c, dtype=float) for c in (candidates or _capacity_candidates(model, span))]
    cands = [np.clip(c, lo + 1e-9, hi - 1e-9) for c in cands]
    scored = sorted(cands, key=lambda c: float(np.sum(resid(c) ** 2)))[:n_starts]

    best = None
    for c in scored:
        s1 = least_squares(resid, c, bounds=(lo, hi), loss="soft_l1", f_scale=f_scale, x_scale="jac")
        s2 = least_squares(resid, s1.x, bounds=(lo, hi), loss="linear", x_scale="jac")
        if best is None or s2.cost < best.cost:
            best = s2
    assert best is not None
    return FitResult(
        model=model,
        params=np.asarray(best.x, dtype=float),
        t=q,
        voltage=voltage,
        loss=float(2 * best.cost / len(q)),
        loss_type="two-stage soft_l1 -> least squares (reported loss = MSE)",
        success=bool(best.success),
        message=str(best.message),
        extra={"n_starts": len(scored)},
    )
