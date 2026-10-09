"""Full-cell open-circuit voltage (OCV) models built from half-cell OCP curves.

Three parameterisations are provided, from the simplest to the most physical:

* :class:`CathodeFixedModel` - the positive electrode window is fixed to the
  full reference curve and only the negative electrode is scaled and shifted
  (``Kn``, ``Sn``) plus a voltage offset. Three parameters, fast and robust,
  but it cannot detect loss of positive active material.
* :class:`WindowModel` - both electrode windows are free (stoichiometry at
  0 % and 100 % cell SOC for each electrode) plus a voltage offset. Optional
  blended positive electrode with free capacity fractions.
* :class:`CapacityModel` - electrode-capacity parameterisation on an absolute
  capacity axis (Ah): positive and negative electrode capacities, cyclable
  lithium and a polarisation offset, optionally anchored to a measured
  full-charge rest voltage. Its parameters map directly onto the degradation
  modes LLI, LAM_PE and LAM_NE.

All models return ``nan`` where an electrode would have to be evaluated outside
its measured OCP range (no extrapolation).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

from .blend import BlendedElectrode
from .ocp import OCPCurve

FloatArray = NDArray[np.float64]

__all__ = ["ModelComponents", "OCVModel", "CathodeFixedModel", "WindowModel", "CapacityModel"]


@dataclass
class ModelComponents:
    """Detailed model output.

    Attributes
    ----------
    voltage
        Full-cell voltage (``nan`` where invalid).
    valid
        Boolean mask of points where both electrodes are inside their data range.
    positive, negative
        Electrode potentials (V vs Li/Li+).
    y, x
        Lithiation fractions of the positive and negative electrodes.
    """

    voltage: FloatArray
    valid: NDArray[np.bool_]
    positive: FloatArray
    negative: FloatArray
    y: FloatArray
    x: FloatArray


class OCVModel:
    """Common interface of all full-cell models."""

    param_names: list[str]

    @property
    def n_params(self) -> int:
        return len(self.param_names)

    def default_bounds(self, t: ArrayLike | None = None) -> list[tuple[float, float]]:
        raise NotImplementedError

    def components(self, params: ArrayLike, t: ArrayLike) -> ModelComponents:
        raise NotImplementedError

    def __call__(self, params: ArrayLike, t: ArrayLike) -> FloatArray:
        """Full-cell voltage at abscissa ``t`` (SOC or capacity, model dependent)."""
        return self.components(params, t).voltage

    def as_dict(self, params: ArrayLike) -> dict[str, float]:
        return {k: float(v) for k, v in zip(self.param_names, np.asarray(params, dtype=float))}


def _assemble(y: FloatArray, x: FloatArray, up: FloatArray, un: FloatArray, offset: float) -> ModelComponents:
    v = up - un + offset
    valid = np.isfinite(v)
    return ModelComponents(
        voltage=np.where(valid, v, np.nan), valid=valid, positive=up, negative=un, y=y, x=x
    )


class CathodeFixedModel(OCVModel):
    """Cathode-fixed electrode alignment (three parameters).

    The positive electrode is assumed to be capacity limiting and its window
    is fixed to the whole reference curve. On the full-cell SOC axis ``s``
    (0 = discharged, 1 = fully charged):

    .. math::

        y(s) = y_{max} - s\\,(y_{max} - y_{min})

        u_n(s) = S_n + s / K_n

        V(s) = U_p(y(s)) - U_n(x(u_n(s))) + V_{offset}

    where ``u_n`` is the *normalised* lithiation coordinate of the negative
    reference curve (``x = x_min + u_n (x_max - x_min)``). The negative axis is
    transformed **before** it is compared with the positive electrode so that
    both potentials are subtracted on the same full-cell SOC coordinate.

    Parameters: ``Kn`` (negative/positive capacity ratio within the window;
    ``Kn > 1`` compresses the negative curve), ``Sn`` (slippage of the negative
    electrode start point; an increase indicates loss of lithium inventory)
    and ``v_offset`` (V).
    """

    param_names = ["Kn", "Sn", "v_offset"]

    def __init__(self, positive: OCPCurve, negative: OCPCurve) -> None:
        self.positive = positive
        self.negative = negative

    def default_bounds(self, t: ArrayLike | None = None) -> list[tuple[float, float]]:
        return [(0.8, 1.3), (-0.2, 0.3), (-0.1, 0.1)]

    def components(self, params: ArrayLike, t: ArrayLike) -> ModelComponents:
        kn, sn, off = (float(p) for p in np.asarray(params, dtype=float))
        s = np.asarray(t, dtype=float)
        ylo, yhi = self.positive.x_range
        xlo, xhi = self.negative.x_range
        y = yhi - s * (yhi - ylo)
        un = sn + s / kn
        x = xlo + un * (xhi - xlo)
        inside = (s >= 0) & (s <= 1) & (un >= 0) & (un <= 1)
        up = np.where(inside, self.positive(np.clip(y, ylo, yhi)), np.nan)
        uneg = np.where(inside, self.negative(np.clip(x, xlo, xhi)), np.nan)
        return _assemble(y, x, up, uneg, off)

    def valid_soc_range(self, params: ArrayLike) -> tuple[float, float]:
        """Full-cell SOC interval in which the transformed negative curve is defined."""
        kn, sn, _ = (float(p) for p in np.asarray(params, dtype=float))
        return max(0.0, -sn * kn), min(1.0, (1.0 - sn) * kn)


class WindowModel(OCVModel):
    """Free electrode windows (five parameters, plus blend fractions).

    .. math::

        y(s) = y_0 + (y_{100} - y_0)\\,s, \\qquad x(s) = x_0 + (x_{100} - x_0)\\,s

        V(s) = U_p(y(s)) - U_n(x(s)) + V_{offset}

    ``y0``/``y100`` are the positive lithiation fractions at 0 % / 100 % cell
    SOC and ``x0``/``x100`` the negative ones. If ``positive`` is a list of
    curves, the positive electrode is a blend (see :mod:`ocvfit.blend`) and
    ``K - 1`` capacity-fraction parameters ``frac_<name>`` are appended; the
    first component takes the remaining fraction. ``y`` is then the
    capacity-weighted lithiation of the blend.
    """

    def __init__(self, positive: OCPCurve | Sequence[OCPCurve], negative: OCPCurve) -> None:
        self.negative = negative
        if isinstance(positive, OCPCurve):
            self.positive_components = [positive]
        else:
            self.positive_components = list(positive)
        self.blend = BlendedElectrode(self.positive_components) if len(self.positive_components) > 1 else None
        self.param_names = ["y0", "y100", "x0", "x100", "v_offset"] + [
            f"frac_{c.name}" for c in self.positive_components[1:]
        ]

    def default_bounds(self, t: ArrayLike | None = None) -> list[tuple[float, float]]:
        b = [(0.5, 1.0), (0.0, 0.6), (0.0, 0.4), (0.5, 1.0), (-0.1, 0.1)]
        return b + [(0.0, 1.0)] * (len(self.positive_components) - 1)

    def positive_curve(self, params: ArrayLike) -> OCPCurve:
        """The (possibly blended) positive OCP for the given parameter vector."""
        p = np.asarray(params, dtype=float)
        if self.blend is None:
            return self.positive_components[0]
        rest = p[5:]
        fr = np.r_[max(0.0, 1.0 - rest.sum()), rest]
        return self.blend.as_curve(fr)

    def components(self, params: ArrayLike, t: ArrayLike) -> ModelComponents:
        p = np.asarray(params, dtype=float)
        y0, y100, x0, x100, off = p[:5]
        s = np.asarray(t, dtype=float)
        if self.blend is not None and p[5:].sum() > 1.0:
            nan = np.full_like(s, np.nan)
            return ModelComponents(nan, np.zeros_like(s, bool), nan, nan, nan, nan)
        y = y0 + (y100 - y0) * s
        x = x0 + (x100 - x0) * s
        if self.blend is None:
            up = self.positive_components[0](y)
        else:
            rest = p[5:]
            fr = np.r_[1.0 - rest.sum(), rest]
            up = self.blend.potential(y, fr)
        return _assemble(y, x, up, self.negative(x), off)


class CapacityModel(OCVModel):
    """Electrode-capacity model on an absolute capacity axis.

    The abscissa ``q`` is the capacity discharged from the fully charged state
    (Ah, increasing during discharge). With ``Q_pe`` (or one ``Q_<name>`` per
    blend component), ``Q_ne``, the cyclable lithium ``n_li`` (all in Ah) and
    the negative lithiation ``x100`` at full charge:

    .. math::

        L_{NE}(q) = x_{100} Q_{NE} - q, \\qquad L_{PE}(q) = n_{Li} - x_{100} Q_{NE} + q

        V(q) = U_{PE}(L_{PE}) - U_{NE}(L_{NE}/Q_{NE}) - R

    where ``U_PE(L)`` is ``U_p(L / Q_pe)`` for a single material and the
    iso-potential blend inverse for several materials. ``R`` (V) lumps
    polarisation and reference mismatch; it is positive for a discharge curve
    measured at a small current.

    Full-charge anchor: if ``v_full`` (the rested voltage at full charge) is
    given, ``x100`` is not a free parameter but is solved from
    ``U_PE(L_PE(0)) - U_NE(x100) = v_full``. Otherwise ``x100`` is fitted.
    """

    def __init__(
        self,
        positive: OCPCurve | Sequence[OCPCurve],
        negative: OCPCurve,
        v_full: float | None = None,
    ) -> None:
        self.negative = negative
        if isinstance(positive, OCPCurve):
            self.positive_components = [positive]
        else:
            self.positive_components = list(positive)
        self.blend = BlendedElectrode(self.positive_components) if len(self.positive_components) > 1 else None
        self.v_full = v_full
        if len(self.positive_components) == 1:
            cap_names = ["Q_pe"]
        else:
            cap_names = [f"Q_{c.name}" for c in self.positive_components]
        self.param_names = cap_names + ["Q_ne", "n_li", "R"] + ([] if v_full is not None else ["x100"])

    # -------------------------------------------------------------- helpers
    @property
    def n_positive(self) -> int:
        return len(self.positive_components)

    def unpack(self, params: ArrayLike) -> dict[str, float | FloatArray]:
        """Split a parameter vector into named quantities (``Q_pe`` = total)."""
        p = np.asarray(params, dtype=float)
        k = self.n_positive
        q_k = p[:k]
        out: dict[str, float | FloatArray] = {
            "Q_k": q_k,
            "Q_pe": float(np.sum(q_k)),
            "Q_ne": float(p[k]),
            "n_li": float(p[k + 1]),
            "R": float(p[k + 2]),
        }
        if self.v_full is None:
            out["x100"] = float(p[k + 3])
        return out

    def positive_potential(self, lithium: ArrayLike, q_k: ArrayLike) -> FloatArray:
        """Positive potential as a function of stored lithium (Ah)."""
        q_k = np.asarray(q_k, dtype=float)
        lith = np.asarray(lithium, dtype=float)
        if self.blend is None:
            return self.positive_components[0](lith / q_k[0])
        return self.blend.potential(lith, q_k)

    def solve_x100(self, q_k: ArrayLike, q_ne: float, n_li: float, v_full: float) -> float:
        """Negative lithiation at full charge that reproduces ``v_full`` (``nan`` if none)."""
        xlo, xhi = self.negative.x_range

        def g(x100: float) -> float:
            up = self.positive_potential(n_li - x100 * q_ne, q_k)
            return float(up - self.negative(x100) - v_full)

        grid = np.linspace(xlo, xhi, 201)
        vals = self.positive_potential(n_li - grid * q_ne, q_k) - self.negative(grid) - v_full
        ok = np.isfinite(vals)
        if ok.sum() < 2:
            return float("nan")
        gx, gv = grid[ok], vals[ok]
        sign = np.where(np.diff(np.sign(gv)) != 0)[0]
        if len(sign) == 0:
            return float("nan")
        i = sign[0]
        try:
            return float(brentq(g, gx[i], gx[i + 1], xtol=1e-12))
        except ValueError:
            return float("nan")

    def x100(self, params: ArrayLike) -> float:
        d = self.unpack(params)
        if self.v_full is None:
            return float(d["x100"])
        return self.solve_x100(d["Q_k"], float(d["Q_ne"]), float(d["n_li"]), float(self.v_full))

    def default_bounds(self, t: ArrayLike | None = None) -> list[tuple[float, float]]:
        span = 1.0 if t is None else float(np.nanmax(t) - np.nanmin(t))
        k = self.n_positive
        b = [(0.0 if k > 1 else 0.5 * span, 3.0 * span)] * k
        b += [(0.5 * span, 3.0 * span), (0.5 * span, 3.0 * span), (-0.1, 0.1)]
        if self.v_full is None:
            b += [(0.0, 1.0)]
        return b

    def components(self, params: ArrayLike, t: ArrayLike) -> ModelComponents:
        d = self.unpack(params)
        q = np.asarray(t, dtype=float)
        q_k = np.asarray(d["Q_k"], dtype=float)
        q_ne, n_li, r = float(d["Q_ne"]), float(d["n_li"]), float(d["R"])
        x100 = self.x100(params)
        if not np.isfinite(x100) or q_ne <= 0 or np.any(q_k < 0) or q_k.sum() <= 0:
            nan = np.full_like(q, np.nan)
            return ModelComponents(nan, np.zeros_like(q, bool), nan, nan, nan, nan)
        l_ne = x100 * q_ne - q
        l_pe = n_li - x100 * q_ne + q
        x = l_ne / q_ne
        y = l_pe / q_k.sum()
        up = self.positive_potential(l_pe, q_k)
        un = self.negative(x)
        return _assemble(y, x, up, un, -r)

    def state(self, params: ArrayLike, q_end: float | None = None) -> dict[str, float]:
        """Derived electrode state: window end points and capacity ratios.

        ``q_end`` is the discharged capacity at the end of the measured curve;
        if given, the 0 %-SOC end points ``x0``/``y0`` are reported for it.
        """
        d = self.unpack(params)
        q_k = np.asarray(d["Q_k"], dtype=float)
        q_pe, q_ne, n_li = float(d["Q_pe"]), float(d["Q_ne"]), float(d["n_li"])
        x100 = self.x100(params)
        out = {
            "Q_pe": q_pe,
            "Q_ne": q_ne,
            "n_li": n_li,
            "R": float(d["R"]),
            "x100": x100,
            "y100": (n_li - x100 * q_ne) / q_pe,
            "np_ratio": q_ne / q_pe,
        }
        for c, qk in zip(self.positive_components, q_k):
            if self.n_positive > 1:
                out[f"Q_{c.name}"] = float(qk)
        if q_end is not None:
            out["x0"] = x100 - q_end / q_ne
            out["y0"] = (n_li - x100 * q_ne + q_end) / q_pe
        return out
