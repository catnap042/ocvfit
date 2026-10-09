"""Minimal CSV helpers (NumPy only)."""

from __future__ import annotations

from os import PathLike

import numpy as np
from numpy.typing import NDArray

from .ocp import OCPCurve

__all__ = ["load_xy", "load_ocp", "save_knots"]


def load_xy(
    path: str | PathLike[str], x_col: str, y_col: str, delimiter: str = ","
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Read two named columns from a delimited text file with a header row."""
    data = np.genfromtxt(path, delimiter=delimiter, names=True, dtype=float, encoding="utf-8")
    return np.asarray(data[x_col], dtype=float), np.asarray(data[y_col], dtype=float)


def load_ocp(
    path: str | PathLike[str],
    x_col: str,
    v_col: str,
    coordinate: str = "lithiation",
    name: str = "electrode",
    delimiter: str = ",",
    smooth_window: int = 1,
) -> OCPCurve:
    """Load a half-cell curve from CSV (see :meth:`OCPCurve.from_soc`)."""
    x, v = load_xy(path, x_col, v_col, delimiter)
    return OCPCurve.from_soc(x, v, coordinate=coordinate, name=name, smooth_window=smooth_window)  # type: ignore[arg-type]


def save_knots(path: str | PathLike[str], knots: NDArray[np.float64], header: str = "x,voltage_V") -> None:
    """Write PCHIP knots (e.g. from :func:`~ocvfit.ocp.compress_pchip`) to CSV."""
    np.savetxt(path, knots, delimiter=",", header=header, comments="", fmt="%.10g")
