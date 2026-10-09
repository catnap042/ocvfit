"""Published half-cell OCP functions used for examples, tests and synthetic data.

The analytic expressions below are taken from the peer-reviewed literature as
implemented in the PyBaMM parameter library (BSD-3-Clause,
https://github.com/pybamm-team/PyBaMM). They are re-typed here so that
``ocvfit`` does not depend on PyBaMM. All functions take the lithiation
fraction (stoichiometry) ``x`` and return the potential in V vs Li/Li+.

References
----------
[Chen2020]
    C.-H. Chen, F. Brosa Planella, K. O'Regan, D. Gastol, W. D. Widanage,
    E. Kendrick, "Development of Experimental Techniques for Parameterization
    of Multi-scale Lithium-ion Battery Models", J. Electrochem. Soc. 167
    (2020) 080534. https://doi.org/10.1149/1945-7111/ab9050
[Mohtat2020]
    P. Mohtat, S. Lee, V. Sulzer, J. B. Siegel, A. G. Stefanopoulou,
    "Differential Expansion and Voltage Model for Li-ion Batteries at
    Practical Charging Rates", J. Electrochem. Soc. 167 (2020) 110561.
    https://doi.org/10.1149/1945-7111/aba5d1
[Afshar2017]
    S. Afshar, K. Morris, A. Khajepour, "Efficient electrochemical model for
    lithium-ion cells", arXiv:1709.03970 (2017). Used by the PyBaMM
    ``Prada2013`` LFP parameter set.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .ocp import OCPCurve

__all__ = [
    "graphite_chen2020",
    "nmc811_chen2020",
    "graphite_mohtat2020",
    "nmc532_mohtat2020",
    "lfp_afshar2017",
    "LITERATURE_OCP",
    "literature_ocp",
]


def graphite_chen2020(x: ArrayLike) -> NDArray[np.float64]:
    """Graphite (LG M50) OCP from [Chen2020]."""
    s = np.asarray(x, dtype=float)
    return (
        1.9793 * np.exp(-39.3631 * s)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (s - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (s - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (s - 0.6103))
    )


def nmc811_chen2020(x: ArrayLike) -> NDArray[np.float64]:
    """NMC811 (LG M50) OCP from [Chen2020]."""
    s = np.asarray(x, dtype=float)
    return (
        -0.8090 * s
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (s - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (s - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (s - 0.3120))
    )


def graphite_mohtat2020(x: ArrayLike) -> NDArray[np.float64]:
    """Graphite OCP from the [Mohtat2020] parameter set."""
    s = np.asarray(x, dtype=float)
    return (
        0.063
        + 0.8 * np.exp(-75 * (s + 0.001))
        - 0.0120 * np.tanh((s - 0.127) / 0.016)
        - 0.0118 * np.tanh((s - 0.155) / 0.016)
        - 0.0035 * np.tanh((s - 0.220) / 0.020)
        - 0.0095 * np.tanh((s - 0.190) / 0.013)
        - 0.0145 * np.tanh((s - 0.490) / 0.020)
        - 0.0800 * np.tanh((s - 1.030) / 0.055)
    )


def nmc532_mohtat2020(x: ArrayLike) -> NDArray[np.float64]:
    """NMC532 OCP from the [Mohtat2020] parameter set."""
    s = np.asarray(x, dtype=float)
    return (
        4.3452
        - 1.6518 * s
        + 1.6225 * s**2
        - 2.0843 * s**3
        + 3.5146 * s**4
        - 2.2166 * s**5
        - 0.5623e-4 * np.exp(109.451 * s - 100.006)
    )


def lfp_afshar2017(x: ArrayLike) -> NDArray[np.float64]:
    """LFP OCP from [Afshar2017] (PyBaMM ``Prada2013`` parameter set)."""
    s = np.asarray(x, dtype=float)
    return 3.4077 - 0.020269 * s + 0.5 * np.exp(-150 * s) - 0.9 * np.exp(-30 * (1 - s))


#: Registry of the available literature OCP functions.
LITERATURE_OCP: dict[str, Callable[[ArrayLike], NDArray[np.float64]]] = {
    "graphite_chen2020": graphite_chen2020,
    "nmc811_chen2020": nmc811_chen2020,
    "graphite_mohtat2020": graphite_mohtat2020,
    "nmc532_mohtat2020": nmc532_mohtat2020,
    "lfp_afshar2017": lfp_afshar2017,
}


def literature_ocp(name: str, x_min: float = 0.0, x_max: float = 1.0, n: int = 1001) -> OCPCurve:
    """Return a literature OCP sampled on ``[x_min, x_max]`` as an :class:`OCPCurve`."""
    try:
        func = LITERATURE_OCP[name]
    except KeyError as exc:
        raise KeyError(f"unknown OCP {name!r}; choose from {sorted(LITERATURE_OCP)}") from exc
    return OCPCurve.from_function(func, x_min, x_max, n, name=name)
