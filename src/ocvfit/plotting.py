"""Matplotlib figures (requires the ``plot`` extra: ``pip install ocvfit[plot]``)."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .fitting import FitResult
from .ica import differential_voltage
from .models import CapacityModel

__all__ = ["plot_fit_report"]


def plot_fit_report(fits: Mapping[str, FitResult], title: str | None = None):  # type: ignore[no-untyped-def]
    """Four-panel report for one or more fits on a common abscissa.

    Panels: (a) data and fitted model, (b) residuals in mV, (c) aligned
    electrode potentials of the fitted model (negative on the right axis),
    (d) |dV/dQ| (or |dV/dSOC|) of data and model. Returns the matplotlib ``Figure``.
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    (ax1, ax2), (ax3, ax4) = axes
    ax3n = ax3.twinx()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    capacity_axis = all(isinstance(f.model, CapacityModel) for f in fits.values())
    xlabel = "Discharged capacity (Ah)" if capacity_axis else "Full-cell SOC"
    dlabel = "|dV/dQ| (V/Ah)" if capacity_axis else "|dV/dSOC| (V)"

    model_slopes = []
    for i, (label, f) in enumerate(fits.items()):
        c = colors[i % len(colors)]
        t, v = f.t, f.voltage
        tt = np.linspace(np.nanmin(t), np.nanmax(t), 600)
        comp = f.model.components(f.params, tt)
        ax1.plot(t, v, ".", color=c, ms=3, alpha=0.6, label=f"{label} data")
        ax1.plot(tt, comp.voltage, "-", color=c, lw=1.6, label=f"{label} fit")
        ax2.plot(t, 1e3 * f.residuals, "-", color=c, lw=1, label=f"{label}: RMSE {1e3 * f.rmse:.2f} mV")
        ax3.plot(tt, comp.positive, "-", color=c, lw=1.6, label=f"{label} positive (left)")
        ax3n.plot(tt, comp.negative, "--", color=c, lw=1.6, label=f"{label} negative (right)")
        ok = np.isfinite(comp.voltage)
        qd, dvd = differential_voltage(t, v, window=21, polyorder=3)
        qm, dvm = differential_voltage(tt[ok], comp.voltage[ok], window=21, polyorder=3)
        model_slopes.append(np.abs(dvm))
        ax4.plot(qd, np.abs(dvd), ".", color=c, ms=2, alpha=0.5)
        ax4.plot(qm, np.abs(dvm), "-", color=c, lw=1.4, label=label)

    ax1.set(xlabel=xlabel, ylabel="Voltage (V)", title="(a) Full-cell voltage")
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set(xlabel=xlabel, ylabel="Model - data (mV)", title="(b) Residuals")
    ax3.set(
        xlabel=xlabel, ylabel="Positive potential (V vs Li/Li+)", title="(c) Aligned electrode potentials"
    )
    ax3n.set_ylabel("Negative potential (V vs Li/Li+)")
    ax4.set(xlabel=xlabel, ylabel=dlabel, title="(d) Differential voltage (dots: data, lines: model)")
    ax4.set_ylim(0, 3.0 * np.nanpercentile(np.concatenate(model_slopes), 80))
    for ax in (ax1, ax2, ax4):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax3.grid(alpha=0.3)
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3n.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, fontsize=8, loc="center left")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
