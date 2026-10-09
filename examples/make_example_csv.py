"""Write small example CSV files in the format expected by ``ocvfit.io``.

All data are synthetic, generated from published OCP functions (LFP: Afshar et
al. 2017; graphite: Chen et al. 2020) and the package's own capacity model, so
they can be redistributed under the repository's MIT licence.

Run:  python examples/make_example_csv.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import ocvfit as of

OUT = Path(__file__).resolve().parent / "data"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(42)

    # Half-cell curves on a lithiation-fraction axis (x = 0 delithiated, 1 lithiated)
    x = np.linspace(0.0, 1.0, 201)
    for name, key in (("lfp_halfcell.csv", "lfp_afshar2017"), ("graphite_halfcell.csv", "graphite_chen2020")):
        u = of.LITERATURE_OCP[key](x) + rng.normal(0, 0.5e-3, x.size)
        np.savetxt(
            OUT / name, np.column_stack([x, u]), delimiter=",", header="x,voltage_V", comments="", fmt="%.6f"
        )

    # Full-cell low-rate discharge: capacity discharged from full charge (Ah) and voltage (V)
    cell = of.demo_cells("lfp", noise=1e-3, n=150, seed=5)
    np.savetxt(
        OUT / "lfp_fullcell_discharge.csv",
        np.column_stack([cell.q["aged"], cell.voltage["aged"]]),
        delimiter=",",
        header="capacity_Ah,voltage_V",
        comments="",
        fmt="%.6f",
    )
    print(f"wrote CSV files to {OUT.name}/")


if __name__ == "__main__":
    main()
