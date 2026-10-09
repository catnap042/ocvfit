"""Degradation-mode quantification from fitted models.

Low-rate OCV data determine three quantities (assuming the half-cell curve
shapes do not change with ageing, inactive material is fully disconnected and
ageing is uniform): the positive electrode capacity ``Q_PE``, the negative
electrode capacity ``Q_NE`` and the cyclable lithium ``n_Li``. Hence

.. math::

    LLI = 1 - n_{Li}^{aged} / n_{Li}^{BOL},\\quad
    LAM_{PE} = 1 - Q_{PE}^{aged} / Q_{PE}^{BOL},\\quad
    LAM_{NE} = 1 - Q_{NE}^{aged} / Q_{NE}^{BOL}.

Here LLI is the *total* loss of lithium inventory. How much of it was consumed
by side reactions and how much was trapped in material that became inactive
cannot be separated from OCV data alone; :func:`side_reaction_interval` returns
the feasible range instead of a point value.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .fitting import FitResult
from .models import CapacityModel, CathodeFixedModel

__all__ = ["degradation_modes", "side_reaction_interval", "cathode_fixed_indicators"]


def _state(fit: FitResult) -> dict[str, float]:
    if not isinstance(fit.model, CapacityModel):
        raise TypeError("degradation_modes requires CapacityModel fits")
    return fit.model.state(fit.params, q_end=float(np.nanmax(fit.t)))


def degradation_modes(bol: FitResult, aged: FitResult) -> dict[str, float]:
    """LLI, LAM_PE (total and per blend component) and LAM_NE as fractions.

    Both fits must use :class:`~ocvfit.models.CapacityModel` with the same
    electrode definitions.
    """
    b, a = _state(bol), _state(aged)
    out = {
        "LLI": 1.0 - a["n_li"] / b["n_li"],
        "LAM_PE": 1.0 - a["Q_pe"] / b["Q_pe"],
        "LAM_NE": 1.0 - a["Q_ne"] / b["Q_ne"],
        "dn_li_Ah": b["n_li"] - a["n_li"],
        "dQ_pe_Ah": b["Q_pe"] - a["Q_pe"],
        "dQ_ne_Ah": b["Q_ne"] - a["Q_ne"],
    }
    assert isinstance(bol.model, CapacityModel)
    if bol.model.n_positive > 1:
        qb = bol.model.unpack(bol.params)["Q_k"]
        qa = aged.model.unpack(aged.params)["Q_k"]  # type: ignore[union-attr]
        for c, x, y in zip(bol.model.positive_components, np.atleast_1d(qb), np.atleast_1d(qa)):
            out[f"LAM_{c.name}"] = 1.0 - float(y) / float(x) if x > 0 else float("nan")
    return out


def side_reaction_interval(bol: FitResult, aged: FitResult) -> dict[str, Any]:
    """Split the total lithium loss into side-reaction loss and trapped lithium.

    With the BOL window end points (``y100``/``y0`` positive, ``x0``/``x100``
    negative lithiation at 100 % / 0 % SOC), lost electrode capacity carries
    away between

    * ``C_min = x0 * dQ_NE + y100 * dQ_PE`` (both electrodes lost in the
      delithiated state) and
    * ``C_max = x100 * dQ_NE + y0 * dQ_PE`` (both lost in the lithiated state)

    of lithium. The side-reaction loss therefore lies in
    ``[max(0, dn_Li - C_max), dn_Li - C_min]`` (lower bound clipped at zero
    because side reactions only consume lithium). ``if_lost_at_full_charge``
    gives the point value under the explicit assumption that material was lost
    at full charge (negative lithiated, positive delithiated).

    Approximation: window end points are taken at BOL; their drift during
    ageing is ignored. All values are fractions of the BOL cyclable lithium.
    """
    b = _state(bol)
    m = degradation_modes(bol, aged)
    dn, dpe, dne = m["dn_li_Ah"], m["dQ_pe_Ah"], m["dQ_ne_Ah"]
    c_min = b["x0"] * dne + b["y100"] * dpe
    c_max = b["x100"] * dne + b["y0"] * dpe
    c_full = b["x100"] * dne + b["y100"] * dpe
    n0 = b["n_li"]
    return {
        "total_LLI": dn / n0,
        "side_reaction_interval": (max(0.0, (dn - c_max) / n0), (dn - c_min) / n0),
        "side_reaction_raw_lower": (dn - c_max) / n0,
        "if_lost_at_full_charge": (dn - c_full) / n0,
        "trapped_lithium_Ah": {"min": c_min, "max": c_max, "full_charge": c_full},
    }


def cathode_fixed_indicators(bol: FitResult, aged: FitResult, q_nominal: float) -> dict[str, float]:
    """Indicative ageing metrics from two :class:`CathodeFixedModel` fits.

    ``LAM_NE ~ -dKn * Q_nominal`` and ``LLI ~ dSn * Q_nominal``. These are
    first-order indicators only: the cathode-fixed parameterisation assumes no
    loss of positive active material, so any LAM_PE is attributed to the other
    modes. Prefer :func:`degradation_modes` with a capacity model for
    quantitative work.
    """
    if not (isinstance(bol.model, CathodeFixedModel) and isinstance(aged.model, CathodeFixedModel)):
        raise TypeError("cathode_fixed_indicators requires CathodeFixedModel fits")
    pb, pa = bol.param_dict, aged.param_dict
    d_kn = pa["Kn"] - pb["Kn"]
    d_sn = pa["Sn"] - pb["Sn"]
    lam = -d_kn * q_nominal
    lli = d_sn * q_nominal
    return {
        "delta_Kn": d_kn,
        "delta_Sn": d_sn,
        "LAM_NE_Ah": lam,
        "LAM_NE_pct": 100 * lam / q_nominal,
        "LLI_Ah": lli,
        "LLI_pct": 100 * lli / q_nominal,
        "delta_v_offset_mV": 1e3 * (pa["v_offset"] - pb["v_offset"]),
    }
