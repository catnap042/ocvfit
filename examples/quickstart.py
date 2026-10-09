"""End-to-end example: degradation modes of a synthetic LFP/graphite cell.

Half-cell curves are published analytic OCP functions: LFP from Afshar et al.
(2017) and graphite from Chen et al. (2020), the pairing used by the PyBaMM
``Prada2013`` parameter set. A begin-of-life (BOL) and an aged low-rate
discharge curve are synthesised with known electrode capacities and cyclable
lithium plus 1 mV Gaussian noise. Both curves are fitted with the
electrode-capacity model, the degradation modes are compared with the ground
truth, and each mode is checked for identifiability.

Run:   python examples/quickstart.py                       (LFP/graphite, headline figure)
       python examples/quickstart.py --chemistry nmc532    (NMC532/graphite, secondary)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import ocvfit as of  # noqa: E402
from ocvfit.plotting import plot_fit_report  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIGURES = {"lfp": "quickstart.png", "nmc532": "nmc532_graphite.png"}
TITLES = {
    "lfp": "LFP/graphite (synthetic; Afshar 2017 LFP + Chen 2020 graphite OCPs): BOL vs. aged",
    "nmc532": "NMC532/graphite (synthetic; Mohtat 2020 OCPs): BOL vs. aged",
}
MODE_PARAM = {"LLI": "n_li", "LAM_PE": "Q_pe", "LAM_NE": "Q_ne"}


def main(chemistry: str = "lfp") -> dict[str, dict[str, float]]:
    cell = of.demo_cells(chemistry, noise=1e-3, seed=0)
    model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
    fits = {label: of.fit_capacity(model, cell.q[label], cell.voltage[label]) for label in ("BOL", "aged")}

    print(f"== {of.datasets.DEMO_DESIGNS[chemistry]['description']}")
    print(f"{'':6s}{'param':>6s} {'fitted':>9s} {'true':>9s}")
    for label, fit in fits.items():
        for k, v in fit.param_dict.items():
            print(f"{label:6s}{k:>6s} {v:9.4f} {cell.truth[label][k]:9.4f}")
        print(f"{label:6s} RMSE {1e3 * fit.rmse:.2f} mV, max |err| {1e3 * fit.max_abs_error:.2f} mV")

    modes = of.degradation_modes(fits["BOL"], fits["aged"])
    truth = cell.true_modes

    def refit(pos, neg):
        m = of.CapacityModel(pos[0], neg, v_full=cell.v_full)
        b = of.fit_capacity(m, cell.q["BOL"], cell.voltage["BOL"], n_starts=3)
        a = of.fit_capacity(m, cell.q["aged"], cell.voltage["aged"], n_starts=3)
        return of.degradation_modes(b, a)

    budget = of.bias_budget(
        refit, [cell.positive], cell.negative, {k: modes[k] for k in MODE_PARAM}, n_draws=8, seed=1
    )
    print("\nDegradation modes: fitted vs. true, final interval = profile interval +/- bias budget")
    print("(bias budget: half-cell curves perturbed by +/-3 mV offset and 2 mV shape error)")
    report: dict[str, dict[str, float]] = {}
    for k in MODE_PARAM:
        lo, hi = of.mode_interval(fits["BOL"], fits["aged"], k)
        lo, hi = lo - budget[k], hi + budget[k]
        half = 0.5 * (hi - lo)
        verdict = "identifiable" if of.is_identifiable(half) else "NOT identifiable (interval only)"
        report[k] = {"fitted": modes[k], "true": truth[k], "low": lo, "high": hi, "bias": budget[k]}
        print(
            f"  {k:7s} {100 * modes[k]:5.2f} % (true {100 * truth[k]:4.2f} %)  "
            f"[{100 * lo:5.2f}, {100 * hi:5.2f}] %  half-width {100 * half:.2f} pp "
            f"(bias {100 * budget[k]:.2f} pp) -> {verdict}"
        )

    side = of.side_reaction_interval(fits["BOL"], fits["aged"])
    lo, hi = side["side_reaction_interval"]
    print(
        f"\nSide-reaction share of the lithium loss: [{100 * lo:.2f}, {100 * hi:.2f}] % "
        "(not determined by OCV data alone)."
    )

    fig = plot_fit_report(fits, title=TITLES[chemistry])
    path = ROOT / "docs" / "images" / FIGURES[chemistry]
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    print(f"\nFigure written to {path.relative_to(ROOT)}")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--chemistry", choices=sorted(FIGURES), default="lfp")
    main(parser.parse_args().chemistry)
