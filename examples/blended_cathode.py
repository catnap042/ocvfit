"""Blended positive electrode (NMC811 + LFP) with one capacity per component.

Lithium is distributed between the two materials at a common potential
(iso-potential capacity addition). The fit returns a separate capacity for each
component, so their losses can be tracked independently - whether they are
*identifiable* depends on which electrode limits the end of discharge.

Run:  python examples/blended_cathode.py
"""

from __future__ import annotations

import ocvfit as of


def main() -> None:
    nmc = of.literature_ocp("nmc811_chen2020")
    nmc.name = "NMC811"
    lfp = of.literature_ocp("lfp_afshar2017")
    lfp.name = "LFP"
    graphite = of.literature_ocp("graphite_chen2020")

    model = of.CapacityModel([nmc, lfp], graphite, v_full=4.2)
    truth = [3.0, 2.0, 5.5, 4.9, 0.005]  # Q_NMC811, Q_LFP, Q_ne, n_li (Ah), R (V)
    q, v = of.synthesize_discharge(model, truth, v_min=2.7, n=250, noise=1e-3, seed=3)

    fit = of.fit_capacity(model, q, v)
    print(fit.summary())
    print("truth:", dict(zip(model.param_names, truth)))

    for name in ("Q_NMC811", "Q_LFP"):
        pr = of.profile_interval(fit, name)
        lo, hi = pr.interval
        print(f"{name}: profile interval [{lo:.3f}, {hi:.3f}] Ah (rho = {pr.rho:.2f})")


if __name__ == "__main__":
    main()
