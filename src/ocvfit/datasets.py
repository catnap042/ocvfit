"""Synthetic, literature-based data for examples and tests.

No measured data are bundled with ``ocvfit``. Half-cell curves are generated
from published analytic OCP functions (see :mod:`ocvfit.literature`) and
full-cell curves are synthesised with the models of this package, optionally
with Gaussian noise. Because the true parameters are known, these data are
suited to checking that a fitting procedure recovers them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from .literature import LITERATURE_OCP
from .models import CapacityModel
from .ocp import OCPCurve

FloatArray = NDArray[np.float64]

__all__ = ["sample_halfcell", "synthesize_discharge", "DemoCell", "DEMO_DESIGNS", "demo_cells"]


def sample_halfcell(
    name: str,
    x_min: float = 0.0,
    x_max: float = 1.0,
    n: int = 200,
    noise: float = 0.0,
    seed: int | None = None,
) -> OCPCurve:
    """Sample a literature OCP at ``n`` points, optionally adding noise (V)."""
    func = LITERATURE_OCP[name]
    rng = np.random.default_rng(seed)
    x = np.linspace(x_min, x_max, n)
    v = func(x) + (rng.normal(0.0, noise, n) if noise > 0 else 0.0)
    return OCPCurve(x, v, name=name)


def synthesize_discharge(
    model: CapacityModel,
    params: NDArray[np.float64] | list[float],
    v_min: float = 3.0,
    n: int = 200,
    noise: float = 1e-3,
    seed: int | None = 0,
) -> tuple[FloatArray, FloatArray]:
    """Discharge curve ``(q, V)`` from full charge down to ``v_min``.

    The end capacity is found by root finding on the noise-free model; ``n``
    points are spaced uniformly in capacity.
    """
    p = np.asarray(params, dtype=float)
    q_hi = 3.0 * float(model.unpack(p)["Q_pe"])
    grid = np.linspace(0.0, q_hi, 4001)
    v = model(p, grid)
    ok = np.isfinite(v)
    below = np.where(ok & (v <= v_min))[0]
    if len(below):
        i = int(below[0])
        q_end = float(brentq(lambda z: float(model(p, np.array([z]))[0]) - v_min, grid[i - 1], grid[i]))
    else:
        q_end = float(grid[np.where(ok)[0][-1]])
    q = np.linspace(0.0, q_end, n)
    rng = np.random.default_rng(seed)
    volt = model(p, q) + (rng.normal(0.0, noise, n) if noise > 0 else 0.0)
    return q, volt


@dataclass
class DemoCell:
    """A synthetic BOL/aged cell pair with known ground truth."""

    positive: OCPCurve
    negative: OCPCurve
    v_full: float
    truth: dict[str, dict[str, float]]
    q: dict[str, FloatArray] = field(default_factory=dict)
    voltage: dict[str, FloatArray] = field(default_factory=dict)

    @property
    def true_modes(self) -> dict[str, float]:
        b, a = self.truth["BOL"], self.truth["aged"]
        return {
            "LLI": 1 - a["n_li"] / b["n_li"],
            "LAM_PE": 1 - a["Q_pe"] / b["Q_pe"],
            "LAM_NE": 1 - a["Q_ne"] / b["Q_ne"],
        }


#: Demo cell designs. Capacities in Ah, R in V, voltages in V.
DEMO_DESIGNS: dict[str, dict] = {
    "lfp": {
        "description": "LFP/graphite: Afshar2017 LFP + Chen2020 graphite (PyBaMM Prada2013 set)",
        "positive": ("lfp_afshar2017", "LFP"),
        "negative": ("graphite_chen2020", "graphite"),
        "v_full": 3.45,
        "v_min": 2.5,
        "truth": {
            "BOL": {"Q_pe": 2.5, "Q_ne": 2.75, "n_li": 2.4125, "R": 0.005},
            "aged": {"Q_pe": 2.5 * 0.95, "Q_ne": 2.75 * 0.94, "n_li": 2.4125 * 0.92, "R": 0.012},
        },
    },
    "nmc532": {
        "description": "NMC532/graphite: Mohtat2020 OCP functions",
        "positive": ("nmc532_mohtat2020", "NMC532"),
        "negative": ("graphite_mohtat2020", "graphite"),
        "v_full": 4.2,
        "v_min": 3.0,
        "truth": {
            "BOL": {"Q_pe": 5.0, "Q_ne": 5.5, "n_li": 4.875, "R": 0.005},
            "aged": {"Q_pe": 5.0 * 0.95, "Q_ne": 5.5 * 0.94, "n_li": 4.875 * 0.92, "R": 0.012},
        },
    },
}


def demo_cells(chemistry: str = "lfp", noise: float = 1e-3, n: int = 300, seed: int = 0) -> DemoCell:
    """Synthetic BOL and aged low-rate discharge curves with known ground truth.

    Parameters
    ----------
    chemistry
        ``"lfp"`` (default): LFP/graphite with the LFP OCP of Afshar et al. 2017
        and the graphite OCP of Chen et al. 2020, the pairing used by the
        PyBaMM ``Prada2013`` parameter set. BOL ``Q_pe=2.5, Q_ne=2.75,
        n_li=2.4125 Ah, R=5 mV``; rested full-charge voltage 3.45 V; discharge
        to 2.5 V. ``"nmc532"``: NMC532/graphite (Mohtat et al. 2020 OCPs), BOL
        ``Q_pe=5.0, Q_ne=5.5, n_li=4.875 Ah, R=5 mV``; 4.2 V to 3.0 V.
        In both designs the aged cell has LLI 8 %, LAM_PE 5 %, LAM_NE 6 % and
        R = 12 mV, and the end of discharge is limited by the negative electrode.
    noise
        Standard deviation of the Gaussian voltage noise (V).
    n
        Points per curve (uniform in capacity).
    seed
        Random seed (BOL uses ``seed``, aged ``seed + 1``).
    """
    try:
        d = DEMO_DESIGNS[chemistry]
    except KeyError as exc:
        raise KeyError(f"unknown chemistry {chemistry!r}; choose from {sorted(DEMO_DESIGNS)}") from exc
    pos = OCPCurve.from_function(LITERATURE_OCP[d["positive"][0]], 0.0, 1.0, 2001, name=d["positive"][1])
    neg = OCPCurve.from_function(LITERATURE_OCP[d["negative"][0]], 0.0, 1.0, 2001, name=d["negative"][1])
    model = CapacityModel(pos, neg, v_full=d["v_full"])
    truth = {k: dict(v) for k, v in d["truth"].items()}
    cell = DemoCell(pos, neg, d["v_full"], truth)
    for i, (label, t) in enumerate(truth.items()):
        p = [t["Q_pe"], t["Q_ne"], t["n_li"], t["R"]]
        q, v = synthesize_discharge(model, p, v_min=d["v_min"], n=n, noise=noise, seed=seed + i)
        cell.q[label], cell.voltage[label] = q, v
    return cell
