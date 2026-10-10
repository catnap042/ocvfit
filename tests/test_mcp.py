"""MCP server: tool functions match the CLI, and a stdio round trip works."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from ocvfit import mcp_server
from ocvfit.cli import main

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"
FILES = {
    "pe": [str(DATA / "lfp_halfcell.csv")],
    "ne": str(DATA / "graphite_halfcell.csv"),
    "bol": str(DATA / "lfp_fullcell_bol.csv"),
    "aged": str(DATA / "lfp_fullcell_aged.csv"),
}


def _strip(doc: dict) -> dict:
    doc = dict(doc)
    for k in ("output_dir", "absolute_files"):
        doc.pop(k, None)
    return doc


def test_tool_functions_match_cli(tmp_path):
    doc = mcp_server.fit_cell(**FILES, v_full=3.45, out_dir=str(tmp_path / "mcp"))
    assert doc["exit_code"] == 0
    assert Path(doc["absolute_files"]["figure"]).is_file()
    cli_args = ["fit", "--pe", FILES["pe"][0], "--ne", FILES["ne"], "--bol", FILES["bol"]]
    main([*cli_args, "--aged", FILES["aged"], "--v-full", "3.45", "--out", str(tmp_path / "cli")])
    cli = json.loads((tmp_path / "cli" / "fit.json").read_text())
    assert json.loads((tmp_path / "mcp" / "fit.json").read_text()) == cli
    assert _strip(doc) == cli


def test_check_and_ica_tools(tmp_path):
    doc = mcp_server.check_data(**FILES, out_dir=str(tmp_path))
    assert doc["exit_code"] == 0
    doc = mcp_server.differential_analysis(bol=FILES["bol"], out_dir=str(tmp_path))
    assert doc["exit_code"] == 0 and doc["differential"]["BOL"]["dqdv_peaks"]["voltage_V"]
    assert "lfp_afshar2017" in mcp_server.literature_references()


def test_stdio_round_trip(tmp_path):
    pytest.importorskip("mcp")
    import anyio
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    params = StdioServerParameters(command=sys.executable, args=["-m", "ocvfit.mcp_server"])

    async def go() -> tuple[list[str], dict]:
        async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
            await session.initialize()
            tools = sorted(t.name for t in (await session.list_tools()).tools)
            res = await session.call_tool(
                "full_report", {**FILES, "v_full": 3.45, "out_dir": str(tmp_path / "rep")}
            )
            return tools, json.loads(res.content[0].text)

    tools, doc = anyio.run(go)
    assert tools == sorted(
        ["check_data", "fit_cell", "degradation_modes", "differential_analysis", "full_report",
         "literature_references"]
    )  # fmt: skip
    assert doc["exit_code"] == 0
    assert doc["degradation"]["modes"]["LLI"]["verdict"] == "identifiable"
    assert (tmp_path / "rep" / "report.png").is_file()
