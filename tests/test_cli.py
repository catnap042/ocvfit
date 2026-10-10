"""Tests of the fixed command-line workflows (``ocvfit`` console script)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import ocvfit as of
from ocvfit.cli import main

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"
PE = str(DATA / "lfp_halfcell.csv")
NE = str(DATA / "graphite_halfcell.csv")
BOL = str(DATA / "lfp_fullcell_bol.csv")
AGED = str(DATA / "lfp_fullcell_aged.csv")
CELL = ["--pe", PE, "--ne", NE, "--bol", BOL]


def run(args: list[str], out: Path) -> tuple[int, dict]:
    code = main([*args, "--out", str(out)])
    command = args[0]
    doc = json.loads((out / f"{command}.json").read_text())
    assert (out / f"{command}.md").exists()
    assert doc["exit_code"] == code
    return code, doc


def write_curve(path: Path, header: str, a: np.ndarray, b: np.ndarray) -> str:
    np.savetxt(path, np.column_stack([a, b]), delimiter=",", header=header, comments="", fmt="%.6f")
    return str(path)


# ----------------------------------------------------------------------- check
def test_check_example_data(tmp_path):
    code, doc = run(["check", *CELL, "--aged", AGED, "--v-full", "3.45"], tmp_path)
    assert code == 0 and doc["status"] == "ok"
    assert doc["inputs"]["bol"]["n_valid"] == 150
    assert len(doc["inputs"]["bol"]["sha256"]) == 64
    assert (tmp_path / "check.png").exists()


@pytest.mark.parametrize(
    "case, message",
    [
        ("mV", "look like mV"),
        ("charge", "charge curve"),
        ("columns", "not found"),
        ("percent", "percent"),
        ("soc_axis", "--pe-coordinate delithiation"),
    ],
)
def test_check_reports_problems(tmp_path, case, message):
    q, v = of.io.load_xy(BOL, "capacity_Ah", "voltage_V")
    x, u = of.io.load_xy(PE, "x", "voltage_V")
    args = ["check"]
    if case == "mV":
        args += ["--bol", write_curve(tmp_path / "b.csv", "capacity_Ah,voltage_V", q, 1e3 * v)]
    elif case == "charge":
        args += ["--bol", write_curve(tmp_path / "b.csv", "capacity_Ah,voltage_V", q, v[::-1])]
    elif case == "columns":
        args += ["--bol", write_curve(tmp_path / "b.csv", "Q,V", q, v)]
    elif case == "percent":
        args += ["--pe", write_curve(tmp_path / "p.csv", "x,voltage_V", 100 * x, u)]
    else:
        args += ["--pe", write_curve(tmp_path / "p.csv", "x,voltage_V", 1 - x, u)]
    code, doc = run(args, tmp_path / "out")
    assert code == 3 and doc["status"] == "error"
    assert any(message in e for e in doc["errors"]), doc["errors"]


def test_check_soc_axis_ok_with_coordinate_flag(tmp_path):
    x, u = of.io.load_xy(PE, "x", "voltage_V")
    p = write_curve(tmp_path / "p.csv", "soc,voltage_V", 1 - x, u)
    code, doc = run(
        ["check", "--pe", p, "--pe-coordinate", "delithiation", "--halfcell-x-col", "soc"], tmp_path / "o"
    )
    assert code == 0, doc["errors"]


def test_usage_errors(tmp_path):
    code, doc = run(["modes", *CELL], tmp_path / "a")  # no aged curve
    assert code == 2 and "aged" in doc["errors"][0]
    code, doc = run(["modes", *CELL, "--aged", AGED, "--model", "window"], tmp_path / "b")
    assert code == 2
    code, doc = run(["fit", "--pe", "no_such_file.csv", "--ne", "graphite", "--bol", BOL], tmp_path / "c")
    assert code == 3 and "built-in reference" in doc["errors"][0]


# ------------------------------------------------------------------------- fit
def test_fit_example_data(tmp_path):
    code, doc = run(["fit", *CELL, "--aged", AGED, "--v-full", "3.45"], tmp_path)
    assert code == 0
    bol, aged = doc["fits"]["BOL"], doc["fits"]["aged"]
    assert bol["rmse_mV"] < 2.0 and aged["rmse_mV"] < 2.0
    assert bol["params"]["n_li"] == pytest.approx(2.4125, rel=0.01)
    assert doc["settings"]["fit"]["capacity"]["n_starts"] == 10
    assert doc["versions"]["ocvfit"] == of.__version__
    assert (tmp_path / "fit.png").exists()


def test_fit_with_literature_references(tmp_path):
    code, doc = run(["fit", "--pe", "lfp", "--ne", "graphite", "--bol", BOL, "--no-figures"], tmp_path)
    assert code == 0
    assert doc["inputs"]["positive"][0]["kind"] == "literature"
    assert "x100" in doc["fits"]["BOL"]["params"]  # no --v-full: x100 fitted
    assert any("v_full not given" in w for w in doc["warnings"])
    assert not (tmp_path / "fit.png").exists()


@pytest.mark.parametrize("model", ["window", "cathode-fixed"])
def test_fit_soc_models(tmp_path, model):
    cell = of.demo_cells("nmc532", noise=1e-3, n=200, seed=0)
    b = write_curve(tmp_path / "b.csv", "capacity_Ah,voltage_V", cell.q["BOL"], cell.voltage["BOL"])
    args = ["fit", "--model", model, "--pe", "nmc532", "--ne", "graphite_mohtat2020", "--bol", b]
    code, doc = run([*args, "--no-figures"], tmp_path / "o")
    assert code == 0
    assert doc["fits"]["BOL"]["model"] in ("WindowModel", "CathodeFixedModel")
    if model == "window":
        assert doc["fits"]["BOL"]["rmse_mV"] < 2.0


# ----------------------------------------------------------------------- modes
@pytest.fixture(scope="module")
def nmc_csv(tmp_path_factory):
    d = tmp_path_factory.mktemp("nmc")
    cell = of.demo_cells("nmc532", noise=1e-3, n=300, seed=0)
    paths = {}
    for label in ("BOL", "aged"):
        paths[label] = write_curve(
            d / f"{label}.csv", "capacity_Ah,voltage_V", cell.q[label], cell.voltage[label]
        )
    return paths, cell.true_modes


def test_modes_nmc532_identifiable(tmp_path, nmc_csv):
    paths, truth = nmc_csv
    args = ["modes", "--pe", "nmc532", "--ne", "graphite_mohtat2020", "--v-full", "4.2"]
    code, doc = run([*args, "--bol", paths["BOL"], "--aged", paths["aged"]], tmp_path)
    assert code == 0
    modes = doc["degradation"]["modes"]
    for k in ("LLI", "LAM_PE", "LAM_NE"):
        assert modes[k]["verdict"] == "identifiable"
        lo, hi = modes[k]["interval"]
        assert lo <= truth[k] <= hi
        assert modes[k]["value"] == pytest.approx(truth[k], abs=0.01)


def test_modes_lfp_lam_pe_not_identifiable(tmp_path):
    code, doc = run(["modes", *CELL, "--aged", AGED, "--v-full", "3.45"], tmp_path)
    assert code == 0 and doc["status"] == "ok_with_warnings"
    modes = doc["degradation"]["modes"]
    assert modes["LLI"]["verdict"] == "identifiable"
    assert modes["LLI"]["interval"][0] <= 0.08 <= modes["LLI"]["interval"][1]
    assert modes["LAM_PE"]["verdict"] == "not_identifiable"
    assert modes["LAM_PE"]["report_as"].startswith("interval only")
    assert any("LAM_PE is not identifiable" in w for w in doc["warnings"])


# ------------------------------------------------------------------------- ica
def test_ica(tmp_path):
    code, doc = run(["ica", "--bol", BOL, "--aged", AGED], tmp_path)
    assert code == 0
    peaks = doc["differential"]["BOL"]["dqdv_peaks"]["voltage_V"]
    assert len(peaks) >= 1 and all(3.1 < p < 3.4 for p in peaks)
    assert (tmp_path / "ica.png").exists() and (tmp_path / "ica_BOL.csv").exists()
    assert doc["differential"]["dqdv_peak_shifts"]


# ---------------------------------------------------------------------- report
def test_report_and_determinism(tmp_path):
    args = ["report", *CELL, "--aged", AGED, "--v-full", "3.45"]
    code1, _ = run(args, tmp_path / "r1")
    code2, _ = run(args, tmp_path / "r2")
    assert code1 == code2 == 0
    j1 = (tmp_path / "r1" / "report.json").read_text()
    j2 = (tmp_path / "r2" / "report.json").read_text()
    assert j1 == j2
    doc = json.loads(j1)
    assert {"fits", "degradation", "differential", "settings", "versions", "inputs"} <= set(doc)
    for f in ("report.png", "report_ica.png", "report.md"):
        assert (tmp_path / "r1" / f).exists()


def test_references_and_version(capsys):
    assert main(["references"]) == 0
    refs = json.loads(capsys.readouterr().out)
    assert "lfp" in refs["lfp_afshar2017"]["aliases"]
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
