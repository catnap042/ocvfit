"""Optional MCP server exposing the fixed ``ocvfit`` workflows as tools.

Install with ``pip install "ocvfit[mcp]"`` and run ``ocvfit-mcp`` (stdio
transport). The tools call exactly the same functions as the command-line
interface (:mod:`ocvfit.workflow`), so results are identical. Inputs are file
paths (absolute paths are recommended, because relative paths are resolved
against the server's working directory); each tool returns the JSON document
written to ``out_dir`` plus absolute paths of the files it produced.

Works with the official ``mcp`` Python SDK 1.x (``FastMCP``) and 2.x
(``MCPServer``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import __version__
from .workflow import Options, list_references, run_command

INSTRUCTIONS = (
    "Deterministic full-cell OCV fitting and degradation-mode analysis for lithium-ion cells. "
    "Inputs are CSV paths: half-cell curves with columns x (lithiation fraction 0-1) and "
    "voltage_V, full-cell low-rate discharge curves with columns capacity_Ah (discharged from "
    "full charge) and voltage_V. Built-in references: lfp, nmc811, nmc532, graphite, "
    "graphite_mohtat2020. Run check_data first, then full_report (or fit_cell / "
    "degradation_modes / differential_analysis). Settings are frozen; do not write your own "
    "fitting code. Report a degradation mode as a value only when its verdict is "
    "'identifiable'; otherwise report the interval only."
)


def _run(command: str, out_dir: str, **kwargs: Any) -> dict[str, Any]:
    out = Path(out_dir).expanduser().resolve()
    pe = kwargs.pop("pe", None) or []
    given = {k: v for k, v in kwargs.items() if v is not None}
    opts = Options(pe=[str(p) for p in pe], out=str(out), **given)
    doc = run_command(command, opts)
    doc["absolute_files"] = {k: str(out / v) for k, v in doc.get("files", {}).items()}
    return doc


def check_data(
    pe: list[str] | None = None,
    ne: str | None = None,
    bol: str | None = None,
    aged: str | None = None,
    v_full: float | None = None,
    v_full_aged: float | None = None,
    pe_coordinate: str = "lithiation",
    ne_coordinate: str = "lithiation",
    halfcell_x_col: str = "x",
    halfcell_v_col: str = "voltage_V",
    q_col: str = "capacity_Ah",
    v_col: str = "voltage_V",
    out_dir: str = "ocvfit_out",
) -> dict[str, Any]:
    """Validate input files (columns, units, monotonicity, coverage) without fitting.

    pe: positive half-cell CSV paths or built-in references (several entries = blend).
    ne: negative half-cell CSV path or built-in reference. bol/aged: full-cell CSV paths.
    Returns status, errors, warnings and per-input diagnostics; exit_code 0 means usable.
    """
    return _run(
        "check",
        out_dir,
        pe=pe,
        ne=ne,
        bol=bol,
        aged=aged,
        v_full=v_full,
        v_full_aged=v_full_aged,
        pe_coordinate=pe_coordinate,
        ne_coordinate=ne_coordinate,
        halfcell_x_col=halfcell_x_col,
        halfcell_v_col=halfcell_v_col,
        q_col=q_col,
        v_col=v_col,
    )


def fit_cell(
    pe: list[str],
    ne: str,
    bol: str,
    aged: str | None = None,
    model: str = "capacity",
    v_full: float | None = None,
    v_full_aged: float | None = None,
    pe_coordinate: str = "lithiation",
    ne_coordinate: str = "lithiation",
    halfcell_x_col: str = "x",
    halfcell_v_col: str = "voltage_V",
    q_col: str = "capacity_Ah",
    v_col: str = "voltage_V",
    out_dir: str = "ocvfit_out",
) -> dict[str, Any]:
    """Fit the BOL (and optional aged) full-cell curve with frozen settings.

    model: 'capacity' (default, recommended), 'window' or 'cathode-fixed' (baseline only).
    v_full: rested full-charge voltage of the BOL cell (recommended).
    Returns fitted parameters, electrode state, RMSE, warnings and the figure path.
    """
    return _run(
        "fit",
        out_dir,
        pe=pe,
        ne=ne,
        bol=bol,
        aged=aged,
        model=model,
        v_full=v_full,
        v_full_aged=v_full_aged,
        pe_coordinate=pe_coordinate,
        ne_coordinate=ne_coordinate,
        halfcell_x_col=halfcell_x_col,
        halfcell_v_col=halfcell_v_col,
        q_col=q_col,
        v_col=v_col,
    )


def degradation_modes(
    pe: list[str],
    ne: str,
    bol: str,
    aged: str,
    v_full: float | None = None,
    v_full_aged: float | None = None,
    pe_coordinate: str = "lithiation",
    ne_coordinate: str = "lithiation",
    halfcell_x_col: str = "x",
    halfcell_v_col: str = "voltage_V",
    q_col: str = "capacity_Ah",
    v_col: str = "voltage_V",
    out_dir: str = "ocvfit_out",
) -> dict[str, Any]:
    """LLI, LAM_PE and LAM_NE with final intervals and identifiability verdicts.

    Uses the capacity model. Each mode has value, interval (profile +/- bias budget),
    half_width and verdict; report a point value only when verdict == 'identifiable'.
    """
    return _run(
        "modes",
        out_dir,
        pe=pe,
        ne=ne,
        bol=bol,
        aged=aged,
        v_full=v_full,
        v_full_aged=v_full_aged,
        pe_coordinate=pe_coordinate,
        ne_coordinate=ne_coordinate,
        halfcell_x_col=halfcell_x_col,
        halfcell_v_col=halfcell_v_col,
        q_col=q_col,
        v_col=v_col,
    )


def differential_analysis(
    bol: str | None = None,
    aged: str | None = None,
    q_col: str = "capacity_Ah",
    v_col: str = "voltage_V",
    out_dir: str = "ocvfit_out",
) -> dict[str, Any]:
    """dQ/dV and dV/dQ of the full-cell curves with peak positions (and BOL->aged peak shifts)."""
    return _run("ica", out_dir, bol=bol, aged=aged, q_col=q_col, v_col=v_col)


def full_report(
    pe: list[str],
    ne: str,
    bol: str,
    aged: str | None = None,
    model: str = "capacity",
    v_full: float | None = None,
    v_full_aged: float | None = None,
    pe_coordinate: str = "lithiation",
    ne_coordinate: str = "lithiation",
    halfcell_x_col: str = "x",
    halfcell_v_col: str = "voltage_V",
    q_col: str = "capacity_Ah",
    v_col: str = "voltage_V",
    title: str | None = None,
    out_dir: str = "ocvfit_out",
) -> dict[str, Any]:
    """One shot: check + fit + degradation modes (if aged) + dQ/dV, with figures and a summary."""
    return _run(
        "report",
        out_dir,
        pe=pe,
        ne=ne,
        bol=bol,
        aged=aged,
        model=model,
        v_full=v_full,
        v_full_aged=v_full_aged,
        pe_coordinate=pe_coordinate,
        ne_coordinate=ne_coordinate,
        halfcell_x_col=halfcell_x_col,
        halfcell_v_col=halfcell_v_col,
        q_col=q_col,
        v_col=v_col,
        title=title,
    )


def literature_references() -> dict[str, Any]:
    """Built-in literature half-cell references usable in place of a half-cell CSV."""
    return list_references()


TOOLS = (check_data, fit_cell, degradation_modes, differential_analysis, full_report, literature_references)


def build_server() -> Any:
    """Create the MCP server (``MCPServer`` on mcp 2.x, ``FastMCP`` on mcp 1.x)."""
    try:
        from mcp.server.mcpserver import MCPServer

        server = MCPServer("ocvfit", instructions=INSTRUCTIONS, version=__version__)
    except ImportError:
        try:
            from mcp.server.fastmcp import FastMCP
        except ImportError as exc:  # pragma: no cover - depends on the environment
            raise SystemExit('the MCP server needs the optional extra: pip install "ocvfit[mcp]"') from exc
        server = FastMCP("ocvfit", instructions=INSTRUCTIONS)
    for func in TOOLS:
        server.tool()(func)
    return server


def main() -> None:
    """Console entry point ``ocvfit-mcp`` (stdio transport)."""
    build_server().run(transport="stdio")


if __name__ == "__main__":  # pragma: no cover
    main()
