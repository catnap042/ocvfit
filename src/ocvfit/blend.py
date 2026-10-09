"""Blended (multi-component) electrodes.

Physical assumption (iso-potential capacity addition): at equilibrium all
active materials in one electrode share the same potential ``U`` and lithium
distributes among them according to their own half-cell curves. The total
lithium stored in the electrode is therefore

.. math::

    L(U) = \\sum_k Q_k \\, x_k(U)

where ``Q_k`` is the capacity of component ``k`` and ``x_k(U)`` its lithiation
fraction obtained by inverting its OCP curve. The blended potential is the
numerical inverse ``U(L)``. Potentials or capacities are never averaged
directly, which would have no physical basis.

Outside the measured potential range of a component its lithiation fraction is
held at the end value (the phase is exhausted or full); the OCP curve itself is
never extrapolated.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .ocp import OCPCurve

__all__ = ["BlendedElectrode"]


class BlendedElectrode:
    """Equilibrium potential of an electrode made of several active materials.

    Parameters
    ----------
    components
        OCP curves of the individual materials (lithiation-fraction axis).
    n_grid
        Number of potential grid points used for the numerical inversion.
        The default (8001) gives a potential step well below 0.1 mV for a
        typical 2 V wide window.
    """

    def __init__(self, components: Sequence[OCPCurve], n_grid: int = 8001) -> None:
        if len(components) == 0:
            raise ValueError("at least one component is required")
        self.components = list(components)
        vmin = min(c.voltage_range[0] for c in self.components)
        vmax = max(c.voltage_range[1] for c in self.components)
        self.potential_grid = np.linspace(vmin, vmax, n_grid)
        # theta[k, j]: lithiation fraction of component k at potential_grid[j]
        self.theta = np.vstack([c.lithiation(self.potential_grid) for c in self.components])

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.components]

    def lithium(self, potential: ArrayLike, capacities: ArrayLike) -> NDArray[np.float64]:
        """Total stored lithium ``sum_k Q_k x_k(U)`` (same unit as ``capacities``)."""
        q = np.asarray(capacities, dtype=float)
        u = np.asarray(potential, dtype=float)
        return np.sum([qk * c.lithiation(u) for qk, c in zip(q, self.components)], axis=0)

    def potential(self, lithium: ArrayLike, capacities: ArrayLike) -> NDArray[np.float64]:
        """Equilibrium potential for a given total lithium content.

        Returns ``nan`` where the requested lithium is outside the range that
        the components can hold.
        """
        q = np.asarray(capacities, dtype=float)
        total = q @ self.theta  # decreasing with potential
        lith = np.asarray(lithium, dtype=float)
        lo, hi = total[-1], total[0]
        out = np.interp(lith, total[::-1], self.potential_grid[::-1])
        return np.where((lith >= lo) & (lith <= hi), out, np.nan)

    def as_curve(self, fractions: ArrayLike, name: str = "blend") -> OCPCurve:
        """Blended OCP on a normalised lithiation axis ``y = sum_k f_k x_k(U)``.

        ``fractions`` are relative capacities (normalised to sum to one).
        """
        f = np.asarray(fractions, dtype=float)
        f = f / np.sum(f)
        y = f @ self.theta
        # keep only strictly varying samples to avoid flat saturated ends
        keep = np.r_[True, np.abs(np.diff(y)) > 1e-12]
        return OCPCurve(y[keep], self.potential_grid[keep], name=name)
