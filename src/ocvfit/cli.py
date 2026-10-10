"""Command-line interface: ``ocvfit <command> [options]``.

All fitting and analysis settings are frozen (see :data:`ocvfit.workflow.SETTINGS`);
the options only describe the inputs. Every command writes ``<command>.json``,
``<command>.md`` and figures to ``--out`` and exits with a documented code::

    0  success (warnings may be present)
    1  unexpected internal error
    2  usage error: missing or incompatible options
    3  input data invalid
    4  fit failed

Examples::

    ocvfit check  --pe lfp_halfcell.csv --ne graphite_halfcell.csv --bol bol.csv --aged aged.csv
    ocvfit report --pe lfp --ne graphite --bol bol.csv --aged aged.csv --v-full 3.45 --out results
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .workflow import EXIT_CODES, MODELS, Options, list_references, run_command

COMMANDS = {
    "check": "validate input CSV files (columns, units, monotonicity, coverage)",
    "fit": "fit the BOL (and optional aged) full-cell curve",
    "modes": "degradation modes LLI / LAM_PE / LAM_NE with intervals and identifiability verdicts",
    "ica": "dQ/dV and dV/dQ curves with peak positions",
    "report": "one shot: check + fit + modes (if --aged) + ica, with figures and summary",
}

EPILOG = "exit codes: " + "; ".join(f"{k} = {v}" for k, v in EXIT_CODES.items())


def _add_inputs(p: argparse.ArgumentParser, command: str) -> None:
    g = p.add_argument_group("inputs")
    if command != "ica":
        g.add_argument(
            "--pe",
            action="append",
            default=[],
            metavar="CSV|REF",
            help="positive half-cell CSV or built-in reference (lfp, nmc811, nmc532); "
            "repeat for blend components",
        )
        g.add_argument(
            "--ne", metavar="CSV|REF", help="negative half-cell CSV or built-in reference (graphite, ...)"
        )
    g.add_argument("--bol", metavar="CSV", help="begin-of-life full-cell low-rate discharge curve")
    g.add_argument("--aged", metavar="CSV", help="aged full-cell low-rate discharge curve")
    if command != "ica":
        g.add_argument(
            "--v-full",
            type=float,
            metavar="V",
            help="rested full-charge voltage of the BOL cell (recommended; else x100 is fitted)",
        )
        g.add_argument(
            "--v-full-aged",
            type=float,
            metavar="V",
            help="rested full-charge voltage of the aged cell (default: --v-full)",
        )
        g.add_argument(
            "--model", choices=MODELS, default="capacity", help="full-cell model (default: capacity)"
        )
        g.add_argument(
            "--pe-coordinate",
            choices=("lithiation", "delithiation"),
            default="lithiation",
            help="x axis of positive half-cell CSVs (delithiation: x = 1 - soc)",
        )
        g.add_argument("--ne-coordinate", choices=("lithiation", "delithiation"), default="lithiation")
        g.add_argument("--halfcell-x-col", default="x", help="half-cell x column (default: x)")
        g.add_argument(
            "--halfcell-v-col", default="voltage_V", help="half-cell voltage column (default: voltage_V)"
        )
    g.add_argument("--q-col", default="capacity_Ah", help="full-cell capacity column (default: capacity_Ah)")
    g.add_argument("--v-col", default="voltage_V", help="full-cell voltage column (default: voltage_V)")
    o = p.add_argument_group("output")
    o.add_argument("--out", default="ocvfit_out", help="output directory (default: ocvfit_out)")
    o.add_argument("--title", help="figure title")
    o.add_argument("--no-figures", action="store_true", help="skip PNG figures")
    o.add_argument("--json", action="store_true", help="print the full JSON document to stdout")
    o.add_argument(
        "--debug", action="store_true", help="show the Python traceback of unexpected internal errors"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ocvfit",
        description="Deterministic full-cell OCV fitting and degradation-mode analysis (frozen settings).",
        epilog=EPILOG,
    )
    parser.add_argument("--version", action="version", version=f"ocvfit {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="command")
    for name, help_ in COMMANDS.items():
        sp = sub.add_parser(name, help=help_, description=help_, epilog=EPILOG)
        _add_inputs(sp, name)
    sub.add_parser("references", help="list built-in literature half-cell references (JSON)")
    sub.add_parser("settings", help="print the frozen settings (JSON)")
    return parser


def options_from_args(args: argparse.Namespace) -> Options:
    return Options(
        pe=list(getattr(args, "pe", []) or []),
        ne=getattr(args, "ne", None),
        bol=args.bol,
        aged=args.aged,
        model=getattr(args, "model", "capacity"),
        v_full=getattr(args, "v_full", None),
        v_full_aged=getattr(args, "v_full_aged", None),
        pe_coordinate=getattr(args, "pe_coordinate", "lithiation"),
        ne_coordinate=getattr(args, "ne_coordinate", "lithiation"),
        halfcell_x_col=getattr(args, "halfcell_x_col", "x"),
        halfcell_v_col=getattr(args, "halfcell_v_col", "voltage_V"),
        q_col=args.q_col,
        v_col=args.v_col,
        out=args.out,
        title=args.title,
        figures=not args.no_figures,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 2
    if args.command == "references":
        print(json.dumps(list_references(), indent=2))
        return 0
    if args.command == "settings":
        from .workflow import SETTINGS

        print(json.dumps(SETTINGS, indent=2))
        return 0
    try:
        doc = run_command(args.command, options_from_args(args), debug=args.debug)
    except OSError as exc:  # e.g. output directory not writable
        if args.debug:
            raise
        print(f"ocvfit {args.command}: cannot write results: {exc.strerror or exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(doc, indent=2))
    elif doc.get("output_dir") is None:
        pass  # nothing was written; the errors are printed below
    else:
        out = doc["output_dir"]
        md = (Path(out) / doc["files"]["summary"]).read_text(encoding="utf-8")
        print(md.rstrip())
        print(f"\nResults: {out}/{doc['files']['json']}")
        for k in ("figure", "ica_figure"):
            if k in doc.get("files", {}):
                print(f"Figure:  {out}/{doc['files'][k]}")
    if doc["exit_code"]:
        print(f"ocvfit {args.command}: failed with exit code {doc['exit_code']}", file=sys.stderr)
        for e in doc["errors"]:
            print(f"  error: {e}", file=sys.stderr)
    return int(doc["exit_code"])


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
