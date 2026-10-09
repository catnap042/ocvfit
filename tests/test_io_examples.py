from pathlib import Path

import numpy as np
import pytest

import ocvfit as of

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"


@pytest.mark.skipif(not DATA.exists(), reason="example data not present")
def test_example_csvs_load_and_fit():
    pe = of.io.load_ocp(DATA / "lfp_halfcell.csv", "x", "voltage_V", name="LFP")
    ne = of.io.load_ocp(DATA / "graphite_halfcell.csv", "x", "voltage_V", name="graphite")
    q, v = of.io.load_xy(DATA / "lfp_fullcell_discharge.csv", "capacity_Ah", "voltage_V")
    assert len(q) == 150 and np.all(np.diff(q) > 0)
    fit = of.fit_capacity(of.CapacityModel(pe, ne, v_full=3.45), q, v)
    assert fit.rmse < 3e-3
