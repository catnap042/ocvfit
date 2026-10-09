"""The three-parameter cathode-fixed method on a synthetic NMC811/graphite cell.

The positive window is fixed to the reference curve and only the negative
curve is scaled (Kn) and shifted (Sn), plus a voltage offset. The enhanced
objective ("v2": Huber loss, edge/slope weights, end-point anchor, valid-region
penalty) is compared with a plain weighted RMSE.

Run:  python examples/cathode_fixed.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402

import ocvfit as of  # noqa: E402
from ocvfit.plotting import plot_fit_report  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "images" / "cathode_fixed.png"


def main() -> None:
    # Chen et al. 2020 OCPs; positive reference measured over the cell window only
    positive = of.literature_ocp("nmc811_chen2020", 0.27, 0.90)
    negative = of.literature_ocp("graphite_chen2020")
    model = of.CathodeFixedModel(positive, negative)

    soc = np.linspace(0, 1, 120)
    true = {"BOL": [1.15, 0.03, 0.005], "aged": [1.08, 0.09, 0.012]}
    rng = np.random.default_rng(7)
    data = {k: model(p, soc) + rng.normal(0, 2e-3, soc.size) for k, p in true.items()}

    fits = {}
    for label in ("BOL", "aged"):
        for loss in ("v2", "rmse"):
            fit = of.fit_ocv(model, soc, data[label], loss=loss, method="fast")
            p = fit.param_dict
            print(
                f"[{label:4s} | {loss:4s}] Kn={p['Kn']:.4f} Sn={p['Sn']:.4f} "
                f"V_offset={1e3 * p['v_offset']:.2f} mV  RMSE={1e3 * fit.rmse:.2f} mV   (true {true[label]})"
            )
            if loss == "v2":
                fits[label] = fit
    print("\nIndicative ageing metrics (5 Ah nominal):")
    for k, v in of.cathode_fixed_indicators(fits["BOL"], fits["aged"], q_nominal=5.0).items():
        print(f"  {k:18s} {v: .4f}")

    fig = plot_fit_report(fits, title="Cathode-fixed method (synthetic NMC811/graphite, Chen 2020 OCPs)")
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=110)
    print(f"\nFigure written to {FIG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
