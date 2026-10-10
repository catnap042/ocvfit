"""Input limits and output-path guards of the CLI and MCP workflows."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from ocvfit import mcp_server, workflow
from ocvfit.cli import main

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"
BOL = str(DATA / "lfp_fullcell_bol.csv")
CELL = ["--pe", "lfp", "--ne", "graphite", "--bol", BOL]


def _check(args: list[str], out: Path) -> tuple[int, dict]:
    code = main(["check", *args, "--out", str(out), "--no-figures"])
    return code, json.loads((out / "check.json").read_text())


def test_too_many_rows_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(workflow, "MAX_INPUT_ROWS", 100)
    code, doc = _check(CELL, tmp_path)
    assert code == 3
    assert any("limit is 100" in e for e in doc["errors"])


def test_too_large_file_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(workflow, "MAX_INPUT_BYTES", 1000)
    code, doc = _check(CELL, tmp_path)
    assert code == 3
    assert any("MiB" in e for e in doc["errors"])


def test_binary_file_rejected(tmp_path):
    f = tmp_path / "blob.csv"
    f.write_bytes(b"capacity_Ah,voltage_V\n\x00\x01\x02\n")
    code, doc = _check(["--pe", "lfp", "--ne", "graphite", "--bol", str(f)], tmp_path / "o")
    assert code == 3
    assert any("not a text" in e for e in doc["errors"])


def test_wrong_file_header_not_echoed_in_full(tmp_path):
    f = tmp_path / "wide.csv"
    header = ",".join(f"{'secret' * 20}{i}" for i in range(50))
    f.write_text(header + "\n" + ",".join("1" for _ in range(50)) + "\n")
    code, doc = _check(["--pe", "lfp", "--ne", "graphite", "--bol", str(f)], tmp_path / "o")
    assert code == 3
    msg = " ".join(doc["errors"])
    assert "38 more" in msg and "secret" * 20 not in msg


def test_nan_inf_and_text_cells_are_dropped(tmp_path):
    q, v = np.loadtxt(BOL, delimiter=",", skiprows=1, unpack=True)
    lines = ["capacity_Ah,voltage_V"] + [f"{a},{b}" for a, b in zip(q, v)]
    lines[5] = "nan,3.3"
    lines[6] = "0.5,inf"
    lines[7] = "abc,-inf"
    f = tmp_path / "dirty.csv"
    f.write_text("\n".join(lines) + "\n")
    code, doc = _check(["--pe", "lfp", "--ne", "graphite", "--bol", str(f)], tmp_path / "o")
    assert code == 0
    assert doc["inputs"]["bol"]["n_valid"] == len(q) - 3


def test_output_dir_that_is_a_file(tmp_path, capsys):
    f = tmp_path / "taken"
    f.write_text("keep me\n")
    code = main(["check", *CELL, "--out", str(f)])
    assert code == 2
    assert f.read_text() == "keep me\n"
    assert "not a directory" in capsys.readouterr().err


@pytest.mark.skipif(os.name == "nt", reason="symlinks need extra privileges on Windows")
def test_symlinked_output_file_refused(tmp_path, capsys):
    victim = tmp_path / "victim.txt"
    victim.write_text("keep me\n")
    out = tmp_path / "out"
    out.mkdir()
    (out / "check.json").symlink_to(victim)
    code = main(["check", *CELL, "--out", str(out)])
    assert code == 2
    assert victim.read_text() == "keep me\n"
    assert "symlinks" in capsys.readouterr().err


def test_internal_error_has_no_traceback_unless_debug(tmp_path, monkeypatch, capsys):
    def boom(opts):
        raise RuntimeError("kaboom")

    monkeypatch.setitem(workflow._RUNNERS, "check", boom)
    assert main(["check", *CELL, "--out", str(tmp_path)]) == 1
    err = capsys.readouterr().err
    assert "kaboom" in err and "Traceback" not in err
    with pytest.raises(RuntimeError, match="kaboom"):
        main(["check", *CELL, "--out", str(tmp_path), "--debug"])


def test_mcp_root_sandbox(tmp_path, monkeypatch):
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setenv(mcp_server.ROOT_ENV, str(root))
    doc = mcp_server.check_data(pe=["lfp"], ne="graphite", bol=BOL, out_dir=str(root / "o"))
    assert doc["exit_code"] == 2 and "outside the allowed root" in doc["errors"][0]
    outside = str(tmp_path / "o")
    doc = mcp_server.check_data(pe=["lfp"], ne="graphite", bol=str(root / "x.csv"), out_dir=outside)
    assert doc["exit_code"] == 2
    data = root / "bol.csv"
    data.write_text(Path(BOL).read_text())
    doc = mcp_server.check_data(pe=["lfp"], ne="graphite", bol=str(data), out_dir=str(root / "o"))
    assert doc["exit_code"] == 0


def test_non_text_extension_rejected(tmp_path):
    f = tmp_path / "curve.json"
    f.write_text(Path(BOL).read_text())
    code, doc = _check(["--pe", "lfp", "--ne", "graphite", "--bol", str(f)], tmp_path / "o")
    assert code == 3
    assert any("not a delimited text file" in e for e in doc["errors"])


def test_missing_matplotlib_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr(workflow, "_has_matplotlib", lambda: False)
    code = main(["check", *CELL, "--out", str(tmp_path)])
    doc = json.loads((tmp_path / "check.json").read_text())
    assert code == 0
    assert any("matplotlib is not installed" in w for w in doc["warnings"])
