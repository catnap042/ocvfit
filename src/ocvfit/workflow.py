"""Fixed, deterministic analysis workflows.

This module is the single code path behind the ``ocvfit`` command-line tool
and the optional MCP server (``ocvfit-mcp``). Every fitting and analysis
setting is frozen in :data:`SETTINGS`, random seeds are fixed, and every result
document echoes the settings, the inputs (with SHA-256 hashes) and the package
versions, so two runs on the same inputs give identical JSON.

Commands
--------
``check``  validate input files (columns, units, monotonicity, coverage).
``fit``    fit the BOL and, optionally, the aged full-cell curve.
``modes``  degradation modes (LLI, LAM_PE, LAM_NE) with intervals and an
           identifiability verdict.
``ica``    dQ/dV and dV/dQ curves with peak positions.
``report`` check + fit + modes + ica in one go, with figures and a summary.

Each command writes ``<command>.json``, ``<command>.md`` and (where useful)
``<command>.png`` into the output directory and returns the JSON document as a
dict. Exit codes are listed in :data:`EXIT_CODES`.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import platform
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from numpy.typing import NDArray
from scipy.signal import find_peaks, savgol_filter

from . import __version__
from .degradation import cathode_fixed_indicators, degradation_modes, side_reaction_interval
from .fitting import FitResult, fit_capacity, fit_ocv
from .ica import differential_voltage, incremental_capacity, resample_uniform
from .literature import LITERATURE_OCP
from .models import CapacityModel, CathodeFixedModel, WindowModel
from .ocp import OCPCurve
from .uncertainty import bias_budget, is_identifiable, mode_interval

FloatArray = NDArray[np.float64]

__all__ = [
    "SETTINGS",
    "EXIT_CODES",
    "LITERATURE_ALIASES",
    "Options",
    "WorkflowError",
    "run_check",
    "run_fit",
    "run_modes",
    "run_ica",
    "run_report",
    "run_command",
    "list_references",
]

SCHEMA_VERSION = "1"

#: Frozen settings used by every command. They are echoed in each JSON output.
SETTINGS: dict[str, Any] = {
    "fit": {
        "capacity": {
            "function": "ocvfit.fit_capacity",
            "n_starts": 10,
            "f_scale_V": 5e-3,
            "invalid_residual_V": 1.0,
            "abscissa": "discharged capacity from full charge (Ah)",
        },
        "soc_models": {
            "function": "ocvfit.fit_ocv",
            "loss": "rmse",
            "method": "fast",
            "seed": 42,
            "n_random": 300,
            "abscissa": "full-cell SOC = 1 - q / q_max of the curve",
        },
        "quality_warning": {"rmse_V": 5e-3, "max_abs_error_V": 20e-3, "valid_fraction": 1.0},
    },
    "modes": {
        "profile": {"function": "ocvfit.mode_interval", "rel_step": 1e-3, "chi2_base": 3.84},
        "bias_budget": {
            "function": "ocvfit.bias_budget",
            "n_draws": 8,
            "offset_V": 3e-3,
            "shape_V": 2e-3,
            "seed": 1,
            "refit_n_starts": 3,
        },
        "identifiable_if_half_width_le": 0.01,
        "final_interval": "profile interval +/- bias budget",
    },
    "ica": {
        "resample_points": 500,
        "presmooth_window": 31,
        "presmooth_polyorder": 3,
        "dqdv_step_V": 5e-3,
        "dqdv_window": 9,
        "dqdv_polyorder": 3,
        "dvdq_window": 21,
        "dvdq_polyorder": 3,
        "peak_prominence_ratio": 0.1,
        "peak_min_separation_V": 0.02,
        "dvdq_peak_min_separation_frac": 0.02,
        "dvdq_interior": [0.05, 0.95],
    },
    "literature_grid_points": 2001,
    "check": {
        "min_points_halfcell": 10,
        "min_points_fullcell": 20,
        "mV_suspect_if_median_above_V": 10.0,
        "nonmonotone_step_tolerance_V": 2e-3,
    },
}

#: Process exit codes of the CLI (also reported as ``exit_code`` in JSON).
EXIT_CODES: dict[int, str] = {
    0: "success (warnings may be present)",
    1: "unexpected internal error",
    2: "usage error: missing or incompatible options",
    3: "input data invalid (see errors)",
    4: "fit failed (no valid solution)",
}

#: Built-in literature references usable in place of a half-cell CSV.
LITERATURE_ALIASES: dict[str, str] = {
    "lfp": "lfp_afshar2017",
    "nmc811": "nmc811_chen2020",
    "nmc532": "nmc532_mohtat2020",
    "graphite": "graphite_chen2020",
    "graphite_mohtat2020": "graphite_mohtat2020",
    **{k: k for k in LITERATURE_OCP},
}

REFERENCE_NOTES = {
    "lfp_afshar2017": "LFP, Afshar et al. 2017 (PyBaMM Prada2013 set); pair with 'graphite'",
    "nmc811_chen2020": "NMC811 (LG M50), Chen et al. 2020; pair with 'graphite'",
    "nmc532_mohtat2020": "NMC532, Mohtat et al. 2020; pair with 'graphite_mohtat2020'",
    "graphite_chen2020": "graphite (LG M50), Chen et al. 2020",
    "graphite_mohtat2020": "graphite, Mohtat et al. 2020",
}

MODES = ("LLI", "LAM_PE", "LAM_NE")
MODE_PARAMETER = {"LLI": "n_li", "LAM_PE": "Q_pe", "LAM_NE": "Q_ne"}
MODELS = ("capacity", "window", "cathode-fixed")


class WorkflowError(Exception):
    """An error with a CLI exit code and a list of messages."""

    def __init__(self, code: int, messages: list[str] | str) -> None:
        self.code = code
        self.messages = [messages] if isinstance(messages, str) else list(messages)
        super().__init__("; ".join(self.messages))


@dataclass
class Options:
    """User options shared by all commands (everything else is frozen)."""

    pe: list[str] = field(default_factory=list)
    ne: str | None = None
    bol: str | None = None
    aged: str | None = None
    model: str = "capacity"
    v_full: float | None = None
    v_full_aged: float | None = None
    pe_coordinate: str = "lithiation"
    ne_coordinate: str = "lithiation"
    halfcell_x_col: str = "x"
    halfcell_v_col: str = "voltage_V"
    q_col: str = "capacity_Ah"
    v_col: str = "voltage_V"
    out: str = "ocvfit_out"
    title: str | None = None
    figures: bool = True


# --------------------------------------------------------------------- utilities
def _versions() -> dict[str, str]:
    out = {
        "ocvfit": __version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "python": platform.python_version(),
    }
    try:
        import matplotlib

        out["matplotlib"] = matplotlib.__version__
    except ImportError:  # pragma: no cover - optional
        pass
    return out


def _clean(obj: Any) -> Any:
    """JSON-safe copy: numpy -> python, floats to 10 significant digits, nan/inf -> None."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_clean(v) for v in obj.tolist()]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        x = float(obj)
        if not math.isfinite(x):
            return None
        return float(f"{x:.10g}")
    return obj


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe_name(text: str) -> str:
    name = re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_")
    return name or "electrode"


def read_table(path: str, label: str) -> tuple[list[str], dict[str, FloatArray], dict[str, int]]:
    """Read a delimited text file with one header row.

    Comma, semicolon and tab delimiters are detected; lines starting with
    ``#`` are skipped. Returns the column names, float arrays (``nan`` where a
    cell is not a number) and counts of skipped/bad cells.
    """
    p = Path(path)
    if not p.is_file():
        raise WorkflowError(3, f"{label}: file not found: {path}")
    text = p.read_text(encoding="utf-8-sig", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if len(lines) < 2:
        raise WorkflowError(3, f"{label}: {path} needs a header row and at least one data row")
    head = lines[0]
    delim = max((",", ";", "\t"), key=head.count)
    rows = list(csv.reader(lines, delimiter=delim))
    header = [h.strip() for h in rows[0]]
    cols: dict[str, list[float]] = {h: [] for h in header}
    bad = 0
    short = 0
    for r in rows[1:]:
        if len(r) < len(header):
            short += 1
        for i, h in enumerate(header):
            cell = r[i].strip() if i < len(r) else ""
            try:
                cols[h].append(float(cell))
            except ValueError:
                cols[h].append(float("nan"))
                bad += 1
    data = {h: np.asarray(v, dtype=float) for h, v in cols.items()}
    return header, data, {"rows": len(rows) - 1, "non_numeric_cells": bad, "short_rows": short}


def _get_columns(
    path: str, label: str, x_col: str, y_col: str
) -> tuple[FloatArray, FloatArray, dict[str, Any]]:
    header, data, counts = read_table(path, label)
    missing = [c for c in (x_col, y_col) if c not in data]
    if missing:
        raise WorkflowError(
            3,
            f"{label}: column(s) {missing} not found in {path}; available columns: {header}. "
            "Rename the columns or pass the column-name options.",
        )
    return data[x_col], data[y_col], counts


# ------------------------------------------------------------------------- check
def _is_reference(spec: str) -> bool:
    return spec in LITERATURE_ALIASES and not Path(spec).is_file()


def _check_halfcell(spec: str, role: str, opts: Options) -> dict[str, Any]:
    """Validate one half-cell input (literature key or CSV path)."""
    errors: list[str] = []
    warnings: list[str] = []
    label = f"{role} half-cell '{spec}'"
    if _is_reference(spec):
        key = LITERATURE_ALIASES[spec]
        info: dict[str, Any] = {"source": spec, "kind": "literature", "reference": key}
        info["note"] = REFERENCE_NOTES.get(key, "")
        if role == "positive" and key.startswith("graphite"):
            errors.append(f"{label}: '{spec}' is a negative-electrode (graphite) reference")
        if role == "negative" and not key.startswith("graphite"):
            errors.append(f"{label}: '{spec}' is a positive-electrode reference")
        warnings.append(
            f"{label}: literature curve used instead of a measured half-cell; absolute capacities "
            "and modes then carry the reference-curve error (covered only partly by the bias budget)"
        )
        info.update(errors=errors, warnings=warnings)
        return info
    if not Path(spec).is_file():
        errors.append(
            f"{label}: not an existing file and not a built-in reference "
            f"(choose from {sorted(set(LITERATURE_ALIASES))})"
        )
        return {"source": spec, "kind": "missing", "errors": errors, "warnings": warnings}
    coordinate = opts.pe_coordinate if role == "positive" else opts.ne_coordinate
    try:
        x, v, counts = _get_columns(spec, label, opts.halfcell_x_col, opts.halfcell_v_col)
    except WorkflowError as exc:
        return {"source": spec, "kind": "csv", "errors": exc.messages, "warnings": warnings}
    info = {
        "source": spec,
        "kind": "csv",
        "sha256": _sha256(Path(spec)),
        "columns": {"x": opts.halfcell_x_col, "voltage": opts.halfcell_v_col},
        "coordinate": coordinate,
        **counts,
    }
    ok = np.isfinite(x) & np.isfinite(v)
    n_bad = int(np.sum(~ok))
    if n_bad:
        warnings.append(f"{label}: {n_bad} row(s) with missing/non-numeric values are ignored")
    x, v = x[ok], v[ok]
    s = SETTINGS["check"]
    if len(x) < s["min_points_halfcell"]:
        errors.append(f"{label}: only {len(x)} valid points (need >= {s['min_points_halfcell']})")
        info.update(errors=errors, warnings=warnings, n_valid=len(x))
        return info
    info["n_valid"] = int(len(x))
    info["x_range"] = [float(x.min()), float(x.max())]
    info["voltage_range_V"] = [float(v.min()), float(v.max())]
    if np.median(v) > s["mV_suspect_if_median_above_V"]:
        errors.append(f"{label}: voltages look like mV (median {np.median(v):.0f}); convert to V")
    if x.min() < -0.02 or x.max() > 1.02:
        hint = " (looks like percent: divide by 100)" if x.max() <= 100.5 and x.max() > 1.5 else ""
        errors.append(
            f"{label}: x must be a lithiation (or SOC) fraction in [0, 1]; got "
            f"[{x.min():.4g}, {x.max():.4g}]{hint}"
        )
    n_dup = int(len(x) - len(np.unique(x)))
    if n_dup:
        warnings.append(f"{label}: {n_dup} duplicate x value(s) are averaged")
    xl = 1.0 - x if coordinate == "delithiation" else x
    order = np.argsort(xl, kind="stable")
    vs = v[order]
    slope = np.polyfit(xl[order], vs, 1)[0] if np.ptp(xl) > 0 else 0.0
    if slope > 0:
        other = "delithiation" if coordinate == "lithiation" else "lithiation"
        flag = "--pe-coordinate" if role == "positive" else "--ne-coordinate"
        errors.append(
            f"{label}: potential increases with lithiation; the axis is probably an SOC/"
            f"delithiation axis - rerun with {flag} {other}"
        )
    rises = int(np.sum(np.diff(vs) > s["nonmonotone_step_tolerance_V"]))
    if rises > 0.05 * len(vs):
        warnings.append(
            f"{label}: {rises} steps rise by more than {1e3 * s['nonmonotone_step_tolerance_V']:.0f} mV "
            "(noisy or non-monotone curve; consider a smoother half-cell curve)"
        )
    if np.ptp(xl) < 0.5:
        warnings.append(f"{label}: x covers only {np.ptp(xl):.2f} of the lithiation range")
    med = float(np.median(v))
    if role == "positive" and not (2.0 <= med <= 5.0):
        warnings.append(f"{label}: median potential {med:.3f} V is unusual for a positive electrode")
    if role == "negative" and med > 1.5:
        warnings.append(f"{label}: median potential {med:.3f} V is unusual for a negative electrode")
    info.update(errors=errors, warnings=warnings)
    return info


def _check_fullcell(path: str, label: str, opts: Options, v_full: float | None) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    name = f"{label} full-cell curve '{path}'"
    try:
        q, v, counts = _get_columns(path, name, opts.q_col, opts.v_col)
    except WorkflowError as exc:
        return {"source": path, "errors": exc.messages, "warnings": warnings}
    info: dict[str, Any] = {
        "source": path,
        "sha256": _sha256(Path(path)),
        "columns": {"capacity": opts.q_col, "voltage": opts.v_col},
        **counts,
    }
    ok = np.isfinite(q) & np.isfinite(v)
    if np.sum(~ok):
        warnings.append(f"{name}: {int(np.sum(~ok))} row(s) with missing/non-numeric values are ignored")
    q, v = q[ok], v[ok]
    s = SETTINGS["check"]
    info["n_valid"] = int(len(q))
    if len(q) < s["min_points_fullcell"]:
        errors.append(f"{name}: only {len(q)} valid points (need >= {s['min_points_fullcell']})")
        info.update(errors=errors, warnings=warnings)
        return info
    info["capacity_range"] = [float(q.min()), float(q.max())]
    info["voltage_range_V"] = [float(v.min()), float(v.max())]
    if np.median(v) > s["mV_suspect_if_median_above_V"]:
        errors.append(f"{name}: voltages look like mV (median {np.median(v):.0f}); convert to V")
    dq = np.diff(q)
    if q[-1] < q[0]:
        errors.append(
            f"{name}: capacity decreases along the file; the capacity column must be the charge "
            "discharged from full charge (0 at full charge, increasing). If you have remaining "
            "capacity, convert with q = Q_total - q_remaining"
        )
    elif np.any(dq < 0):
        errors.append(
            f"{name}: capacity is not monotonically non-decreasing ({int(np.sum(dq < 0))} backward "
            "step(s)); supply one discharge branch only, sorted by capacity"
        )
    if np.sum(dq == 0):
        warnings.append(f"{name}: {int(np.sum(dq == 0))} repeated capacity value(s)")
    span = float(np.ptp(q))
    if span <= 0:
        errors.append(f"{name}: capacity column is constant")
    elif q.min() > 0.02 * span or q.min() < -0.02 * span:
        warnings.append(
            f"{name}: capacity starts at {q.min():.4g}, not 0; the capacity model assumes q = 0 is "
            "the fully charged state (subtract the first value if the curve starts at full charge)"
        )
    if len(v) > 2 and v[-1] >= v[0]:
        errors.append(
            f"{name}: voltage does not fall with capacity - this looks like a charge curve. Use a "
            "low-rate discharge curve (or convert: q = q_max - q, reversed)"
        )
    if np.ptp(v) < 0.2:
        warnings.append(f"{name}: voltage span is only {1e3 * np.ptp(v):.0f} mV")
    if q.max() > 1000:
        warnings.append(f"{name}: capacities above 1000 - expected Ah (values are used as given)")
    if v_full is not None and len(v):
        if v[0] > v_full + 0.02:
            warnings.append(
                f"{name}: first voltage {v[0]:.4f} V is above v_full {v_full:.4f} V; check v_full"
            )
        elif v_full - v[0] > 0.3:
            warnings.append(
                f"{name}: v_full {v_full:.4f} V is far above the first voltage {v[0]:.4f} V; check v_full"
            )
    info.update(errors=errors, warnings=warnings)
    return info


def _collect(items: dict[str, Any]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    for item in items.values():
        if isinstance(item, list):
            for sub in item:
                errors += sub.get("errors", [])
                warnings += sub.get("warnings", [])
        elif isinstance(item, dict):
            errors += item.get("errors", [])
            warnings += item.get("warnings", [])
    return errors, warnings


def _check_inputs(opts: Options, need_halfcells: bool, need_bol: bool) -> dict[str, Any]:
    usage: list[str] = []
    if opts.model not in MODELS:
        usage.append(f"unknown model {opts.model!r}; choose from {list(MODELS)}")
    if opts.pe_coordinate not in ("lithiation", "delithiation"):
        usage.append("pe_coordinate must be 'lithiation' or 'delithiation'")
    if opts.ne_coordinate not in ("lithiation", "delithiation"):
        usage.append("ne_coordinate must be 'lithiation' or 'delithiation'")
    if need_halfcells and (not opts.pe or not opts.ne):
        usage.append("positive (--pe) and negative (--ne) half-cell inputs are required")
    if need_bol and not opts.bol:
        usage.append("a BOL full-cell curve (--bol) is required")
    if opts.model == "cathode-fixed" and len(opts.pe) > 1:
        usage.append("the cathode-fixed model supports a single positive material only")
    if usage:
        raise WorkflowError(2, usage)
    items: dict[str, Any] = {}
    if opts.pe:
        items["positive"] = [_check_halfcell(s, "positive", opts) for s in opts.pe]
    if opts.ne:
        items["negative"] = _check_halfcell(opts.ne, "negative", opts)
    if opts.bol:
        items["bol"] = _check_fullcell(opts.bol, "BOL", opts, opts.v_full)
    if opts.aged:
        vfa = opts.v_full_aged if opts.v_full_aged is not None else opts.v_full
        items["aged"] = _check_fullcell(opts.aged, "aged", opts, vfa)
    if opts.v_full is None and opts.model == "capacity" and opts.bol:
        items["note"] = {
            "warnings": [
                "v_full not given: x100 is fitted as a free parameter. Supplying the rested "
                "full-charge voltage (--v-full) usually improves conditioning."
            ]
        }
    return items


# -------------------------------------------------------------------- building
def _load_curve(spec: str, role: str, opts: Options) -> OCPCurve:
    if _is_reference(spec):
        key = LITERATURE_ALIASES[spec]
        return OCPCurve.from_function(
            LITERATURE_OCP[key], 0.0, 1.0, SETTINGS["literature_grid_points"], name=_safe_name(spec)
        )
    x, v, _ = _get_columns(spec, role, opts.halfcell_x_col, opts.halfcell_v_col)
    coordinate = opts.pe_coordinate if role == "positive" else opts.ne_coordinate
    return OCPCurve.from_soc(x, v, coordinate=coordinate, name=_safe_name(Path(spec).stem))  # type: ignore[arg-type]


def _electrodes(opts: Options) -> tuple[list[OCPCurve], OCPCurve]:
    pos = [_load_curve(s, "positive", opts) for s in opts.pe]
    seen: dict[str, int] = {}
    for c in pos:
        n = seen.get(c.name, 0)
        seen[c.name] = n + 1
        if n:
            c.name = f"{c.name}_{n + 1}"
    neg = _load_curve(opts.ne, "negative", opts)  # type: ignore[arg-type]
    return pos, neg


def _fullcell(path: str, opts: Options) -> tuple[FloatArray, FloatArray]:
    q, v, _ = _get_columns(path, "full-cell", opts.q_col, opts.v_col)
    ok = np.isfinite(q) & np.isfinite(v)
    return q[ok], v[ok]


def _capacity_model(pos: list[OCPCurve], neg: OCPCurve, v_full: float | None) -> CapacityModel:
    return CapacityModel(pos[0] if len(pos) == 1 else pos, neg, v_full=v_full)


def _fit_one(
    opts: Options, pos: list[OCPCurve], neg: OCPCurve, q: FloatArray, v: FloatArray, v_full: float | None
) -> FitResult:
    if opts.model == "capacity":
        s = SETTINGS["fit"]["capacity"]
        model = _capacity_model(pos, neg, v_full)
        fit = fit_capacity(
            model,
            q,
            v,
            n_starts=s["n_starts"],
            f_scale=s["f_scale_V"],
            invalid_residual=s["invalid_residual_V"],
        )
    else:
        s = SETTINGS["fit"]["soc_models"]
        soc = 1.0 - q / float(np.max(q))
        m: CathodeFixedModel | WindowModel
        m = (
            CathodeFixedModel(pos[0], neg)
            if opts.model == "cathode-fixed"
            else WindowModel(pos[0] if len(pos) == 1 else pos, neg)
        )
        fit = fit_ocv(m, soc, v, loss=s["loss"], method=s["method"], seed=s["seed"], n_random=s["n_random"])
    if not np.all(np.isfinite(fit.params)) or not np.isfinite(fit.rmse):
        raise WorkflowError(4, "fit failed: no parameter set gives a finite model over the data")
    return fit


def _fit_summary(fit: FitResult, label: str) -> tuple[dict[str, Any], list[str]]:
    d = fit.to_dict()
    warnings: list[str] = []
    qw = SETTINGS["fit"]["quality_warning"]
    if not fit.success:
        warnings.append(f"{label}: optimiser reported no convergence ({fit.message})")
    if fit.rmse > qw["rmse_V"]:
        warnings.append(
            f"{label}: RMSE {1e3 * fit.rmse:.2f} mV exceeds {1e3 * qw['rmse_V']:.0f} mV - half-cell "
            "curves may not match this cell, or the curve is not near-equilibrium"
        )
    if fit.max_abs_error > qw["max_abs_error_V"]:
        warnings.append(f"{label}: max |error| {1e3 * fit.max_abs_error:.1f} mV")
    if fit.valid_fraction < qw["valid_fraction"]:
        warnings.append(
            f"{label}: only {100 * fit.valid_fraction:.1f} % of points are inside the half-cell data "
            "range (the model is undefined outside it); extend the half-cell curves"
        )
    at_bound = []
    for name, val, (lo, hi) in zip(fit.model.param_names, fit.params, fit.model.default_bounds(fit.t)):
        tol = 1e-4 * (hi - lo)
        if val <= lo + tol or val >= hi - tol:
            at_bound.append(name)
    if at_bound:
        hint = (
            " (the cathode-fixed model needs the positive window to span the whole reference curve)"
            if isinstance(fit.model, CathodeFixedModel)
            else ""
        )
        warnings.append(f"{label}: parameter(s) {at_bound} at their bounds - result unreliable{hint}")
    d["parameters_at_bounds"] = at_bound
    d["rmse_mV"] = 1e3 * fit.rmse
    d["max_abs_error_mV"] = 1e3 * fit.max_abs_error
    d["n_points"] = int(len(fit.t))
    if isinstance(fit.model, CapacityModel):
        d["notes"] = [
            "R is a lumped voltage offset (polarisation + reference mismatch) at the test current, "
            "not a DC resistance"
        ]
    elif isinstance(fit.model, CathodeFixedModel):
        d["notes"] = ["cathode-fixed model assumes no LAM_PE; use it only as a baseline"]
    return d, warnings


# --------------------------------------------------------------------- outputs
def _base(command: str, opts: Options) -> dict[str, Any]:
    return {
        "tool": "ocvfit",
        "schema_version": SCHEMA_VERSION,
        "command": command,
        "status": "ok",
        "exit_code": 0,
        "errors": [],
        "warnings": [],
        "versions": _versions(),
        "options": {k: v for k, v in asdict(opts).items() if k not in ("out", "figures")},
        "settings": copy.deepcopy(SETTINGS),
    }


def _finish(doc: dict[str, Any], out: Path, md: str) -> dict[str, Any]:
    _set_status(doc)
    doc.setdefault("files", {})
    doc["files"]["json"] = f"{doc['command']}.json"
    doc["files"]["summary"] = f"{doc['command']}.md"
    clean = _clean(doc)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{doc['command']}.json").write_text(json.dumps(clean, indent=2, sort_keys=False) + "\n")
    (out / f"{doc['command']}.md").write_text(md)
    clean["output_dir"] = str(out)
    return clean


def _pct(x: float | None) -> str:
    return "n/a" if x is None or not np.isfinite(x) else f"{100 * x:.2f} %"


def _set_status(doc: dict[str, Any]) -> None:
    doc["warnings"] = list(dict.fromkeys(doc["warnings"]))
    if doc["errors"]:
        doc["status"] = "error"
    elif doc["warnings"]:
        doc["status"] = "ok_with_warnings"
    else:
        doc["status"] = "ok"


def _md_header(doc: dict[str, Any], title: str) -> list[str]:
    _set_status(doc)
    lines = [f"# ocvfit {title}", "", f"- status: **{doc['status']}** (exit code {doc['exit_code']})"]
    lines.append(f"- ocvfit {doc['versions']['ocvfit']}, model: {doc['options']['model']}")
    if doc["errors"]:
        lines += ["", "## Errors", ""] + [f"- {e}" for e in doc["errors"]]
    if doc["warnings"]:
        lines += ["", "## Warnings", ""] + [f"- {w}" for w in doc["warnings"]]
    return lines


def _md_fits(fits: dict[str, Any]) -> list[str]:
    if not fits:
        return []
    labels = list(fits)
    keys = list(fits[labels[0]]["params"])
    if "state" in fits[labels[0]]:
        keys += [k for k in ("x100", "y100", "x0", "y0", "np_ratio") if k in fits[labels[0]]["state"]]
    lines = ["", "## Fit", "", "| quantity | " + " | ".join(labels) + " |", "|---|" + "---|" * len(labels)]
    for k in keys:
        row = []
        for lab in labels:
            val = fits[lab]["params"].get(k, fits[lab].get("state", {}).get(k))
            row.append("" if val is None else f"{val:.5g}")
        lines.append(f"| {k} | " + " | ".join(row) + " |")
    lines.append("| RMSE (mV) | " + " | ".join(f"{fits[lab]['rmse_mV']:.3f}" for lab in labels) + " |")
    lines.append(
        "| valid points | "
        + " | ".join(f"{100 * fits[lab]['valid_fraction']:.1f} %" for lab in labels)
        + " |"
    )
    return lines


def _md_modes(modes: dict[str, Any]) -> list[str]:
    if not modes:
        return []
    lines = [
        "",
        "## Degradation modes",
        "",
        "| mode | value | final interval | half-width | verdict |",
        "|---|---|---|---|---|",
    ]
    for k, m in modes["modes"].items():
        iv = m.get("interval")
        ivs = "n/a" if iv is None or None in iv else f"[{100 * iv[0]:.2f}, {100 * iv[1]:.2f}] %"
        hw = m.get("half_width")
        hws = "n/a" if hw is None else f"{100 * hw:.2f} pp"
        val = _pct(m["value"]) if m["verdict"] == "identifiable" else f"({_pct(m['value'])})"
        lines.append(f"| {k} | {val} | {ivs} | {hws} | {m['verdict']} |")
    lines += [
        "",
        "Values in parentheses are not identifiable: report the interval only.",
    ]
    sr = modes.get("side_reaction_interval")
    if sr:
        lines.append(
            f"Side-reaction share of the lithium loss: [{_pct(sr[0])}, {_pct(sr[1])}] of BOL lithium "
            "(not determined by OCV data alone)."
        )
    return lines


def _md_ica(ica: dict[str, Any]) -> list[str]:
    if not ica:
        return []
    lines = ["", "## dQ/dV and dV/dQ peaks", ""]
    for lab, r in ica.items():
        if "error" in r:
            lines.append(f"- {lab}: {r['error']}")
            continue
        pk = ", ".join(f"{p:.3f}" for p in r["dqdv_peaks"]["voltage_V"]) or "none"
        dv = ", ".join(f"{p:.3f}" for p in r["dvdq_peaks"]["capacity"]) or "none"
        lines.append(f"- {lab}: dQ/dV peaks at {pk} V; dV/dQ peaks at capacity {dv}")
    return lines


# --------------------------------------------------------------------- figures
def _plt():  # type: ignore[no-untyped-def]
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def _fit_figure(fits: dict[str, FitResult], path: Path, title: str | None) -> bool:
    plt = _plt()
    if plt is None:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    from .plotting import plot_fit_report

    fig = plot_fit_report(fits, title=title)
    fig.savefig(path, dpi=110, metadata={"Software": None})
    plt.close(fig)
    return True


def _check_figure(opts: Options, path: Path) -> bool:
    plt = _plt()
    if plt is None:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for role, specs in (("positive", opts.pe), ("negative", [opts.ne] if opts.ne else [])):
        for s in specs:
            try:
                c = _load_curve(s, role, opts)
            except WorkflowError:
                continue
            a1.plot(c.x, c.voltage, "-", lw=1.3, label=f"{role}: {c.name}")
    a1.set(xlabel="Lithiation fraction x", ylabel="Potential (V vs Li/Li+)", title="Half-cell curves")
    for lab, p in (("BOL", opts.bol), ("aged", opts.aged)):
        if p:
            try:
                q, v = _fullcell(p, opts)
            except WorkflowError:
                continue
            a2.plot(q, v, ".", ms=2.5, label=lab)
    a2.set(xlabel=f"{opts.q_col}", ylabel="Voltage (V)", title="Full-cell curves")
    for a in (a1, a2):
        a.grid(alpha=0.3)
        if a.get_legend_handles_labels()[0]:
            a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110, metadata={"Software": None})
    plt.close(fig)
    return True


def _ica_figure(curves: dict[str, dict[str, FloatArray]], peaks: dict[str, Any], path: Path) -> bool:
    plt = _plt()
    if plt is None:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for lab, c in curves.items():
        (ln,) = a1.plot(c["V"], c["dqdv"], "-", lw=1.3, label=lab)
        pv = np.asarray(peaks[lab]["dqdv_peaks"]["voltage_V"], dtype=float)
        if len(pv):
            a1.plot(pv, np.interp(pv, c["V"], c["dqdv"]), "v", color=ln.get_color())
        a2.plot(c["Q"], c["dvdq"], "-", lw=1.3, label=lab, color=ln.get_color())
        pq = np.asarray(peaks[lab]["dvdq_peaks"]["capacity"], dtype=float)
        if len(pq):
            a2.plot(pq, np.interp(pq, c["Q"], c["dvdq"]), "v", color=ln.get_color())
    a1.set(xlabel="Voltage (V)", ylabel="|dQ/dV| (Ah/V)", title="Incremental capacity")
    a2.set(xlabel="Discharged capacity", ylabel="|dV/dQ| (V/Ah)", title="Differential voltage")
    allv = np.concatenate([c["dvdq"][c["interior"]] for c in curves.values()])
    if len(allv) and np.isfinite(allv).any():
        a2.set_ylim(0, 1.5 * np.nanmax(allv))
    for a in (a1, a2):
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110, metadata={"Software": None})
    plt.close(fig)
    return True


# -------------------------------------------------------------------- commands
def _usage_guard(doc: dict[str, Any], exc: WorkflowError, out: Path, title: str) -> dict[str, Any]:
    doc["errors"] += exc.messages
    doc["exit_code"] = exc.code
    doc["status"] = "error"
    return _finish(doc, out, "\n".join(_md_header(doc, title)) + "\n")


def run_check(opts: Options) -> dict[str, Any]:
    """Validate the input files without fitting."""
    out = Path(opts.out)
    doc = _base("check", opts)
    try:
        if not (opts.pe or opts.ne or opts.bol or opts.aged):
            raise WorkflowError(2, "nothing to check: give at least one of --pe, --ne, --bol, --aged")
        items = _check_inputs(opts, need_halfcells=False, need_bol=False)
    except WorkflowError as exc:
        return _usage_guard(doc, exc, out, "check")
    doc["inputs"] = items
    doc["errors"], doc["warnings"] = _collect(items)
    if doc["errors"]:
        doc["exit_code"] = 3
    if opts.figures and not doc["errors"] and _check_figure(opts, out / "check.png"):
        doc["files"] = {"figure": "check.png"}
    lines = _md_header(doc, "check")
    lines += ["", "## Inputs", ""]
    for k, item in items.items():
        if k == "note":
            continue
        for it in item if isinstance(item, list) else [item]:
            rng = it.get("voltage_range_V")
            extra = f", {it.get('n_valid', 0)} points, V in [{rng[0]:.3f}, {rng[1]:.3f}]" if rng else ""
            lines.append(f"- {k}: `{it['source']}` ({it.get('kind', 'csv')}{extra})")
    return _finish(doc, out, "\n".join(lines) + "\n")


def _prepare(opts: Options, doc: dict[str, Any], need_aged: bool) -> dict[str, Any]:
    items = _check_inputs(opts, need_halfcells=True, need_bol=True)
    if need_aged and not opts.aged:
        raise WorkflowError(2, "an aged full-cell curve (--aged) is required for degradation modes")
    if need_aged and opts.model != "capacity":
        raise WorkflowError(2, "degradation modes require --model capacity")
    doc["inputs"] = items
    errors, warnings = _collect(items)
    doc["warnings"] += warnings
    if errors:
        raise WorkflowError(3, errors)
    return items


def _fits(opts: Options) -> tuple[dict[str, FitResult], list[OCPCurve], OCPCurve]:
    pos, neg = _electrodes(opts)
    fits: dict[str, FitResult] = {}
    q, v = _fullcell(opts.bol, opts)  # type: ignore[arg-type]
    fits["BOL"] = _fit_one(opts, pos, neg, q, v, opts.v_full)
    if opts.aged:
        vfa = opts.v_full_aged if opts.v_full_aged is not None else opts.v_full
        q, v = _fullcell(opts.aged, opts)
        fits["aged"] = _fit_one(opts, pos, neg, q, v, vfa)
    return fits, pos, neg


def _fit_section(opts: Options, fits: dict[str, FitResult], doc: dict[str, Any]) -> dict[str, Any]:
    section: dict[str, Any] = {}
    for lab, f in fits.items():
        section[lab], w = _fit_summary(f, lab)
        doc["warnings"] += w
    if opts.model == "cathode-fixed" and "aged" in fits:
        q_nom = float(np.max(_fullcell(opts.bol, opts)[0]))  # type: ignore[arg-type]
        section["cathode_fixed_indicators"] = cathode_fixed_indicators(fits["BOL"], fits["aged"], q_nom)
        section["cathode_fixed_indicators"]["note"] = (
            "first-order baseline indicators; assume no LAM_PE. Use --model capacity for modes."
        )
    return section


def run_fit(opts: Options) -> dict[str, Any]:
    """Fit the BOL (and optional aged) full-cell curve with frozen settings."""
    out = Path(opts.out)
    doc = _base("fit", opts)
    try:
        _prepare(opts, doc, need_aged=False)
        fits, _, _ = _fits(opts)
    except WorkflowError as exc:
        return _usage_guard(doc, exc, out, "fit")
    doc["fits"] = _fit_section(opts, fits, doc)
    if opts.figures and _fit_figure(fits, out / "fit.png", opts.title):
        doc["files"] = {"figure": "fit.png"}
    lines = _md_header(doc, "fit") + _md_fits({k: v for k, v in doc["fits"].items() if k in fits})
    return _finish(doc, out, "\n".join(lines) + "\n")


def _modes_section(
    opts: Options, fits: dict[str, FitResult], pos: list[OCPCurve], neg: OCPCurve
) -> tuple[dict[str, Any], list[str]]:
    bol, aged = fits["BOL"], fits["aged"]
    values = degradation_modes(bol, aged)
    s = SETTINGS["modes"]
    bs = s["bias_budget"]
    vfa = opts.v_full_aged if opts.v_full_aged is not None else opts.v_full
    blend = len(pos) > 1

    def refit(p: list[OCPCurve], n: OCPCurve) -> dict[str, float]:
        b = fit_capacity(
            _capacity_model(list(p), n, opts.v_full), bol.t, bol.voltage, n_starts=bs["refit_n_starts"]
        )
        a = fit_capacity(
            _capacity_model(list(p), n, vfa), aged.t, aged.voltage, n_starts=bs["refit_n_starts"]
        )
        return degradation_modes(b, a)

    budget = bias_budget(
        refit,
        pos,
        neg,
        {k: values[k] for k in MODES},
        n_draws=bs["n_draws"],
        offset=bs["offset_V"],
        shape=bs["shape_V"],
        seed=bs["seed"],
    )
    warnings: list[str] = []
    modes: dict[str, Any] = {}
    for k in MODES:
        entry: dict[str, Any] = {"value": values[k], "parameter": MODE_PARAMETER[k], "bias_budget": budget[k]}
        if blend and k == "LAM_PE":
            entry.update(
                profile_interval=None,
                interval=None,
                half_width=None,
                identifiable=False,
                verdict="not_assessed",
                report_as="not assessed for blended positive electrodes; report as indicative only",
            )
        else:
            lo, hi = mode_interval(bol, aged, k)
            flo, fhi = lo - budget[k], hi + budget[k]
            half = 0.5 * (fhi - flo)
            ident = is_identifiable(half, s["identifiable_if_half_width_le"])
            entry.update(
                profile_interval=[lo, hi],
                interval=[flo, fhi],
                half_width=half,
                identifiable=ident,
                verdict="identifiable" if ident else "not_identifiable",
                report_as="value with interval" if ident else "interval only (no point value)",
                interval_contains_zero=bool(flo <= 0.0 <= fhi),
            )
            if not ident:
                warnings.append(
                    f"{k} is not identifiable (half-width {100 * half:.2f} pp > "
                    f"{100 * s['identifiable_if_half_width_le']:.0f} pp): report the interval only"
                )
        modes[k] = entry
    if blend:
        for key, val in values.items():
            if key.startswith("LAM_") and key not in MODES:
                modes[key] = {
                    "value": val,
                    "verdict": "not_assessed",
                    "identifiable": False,
                    "report_as": "per-component LAM of a blend; indicative only",
                }
    side = side_reaction_interval(bol, aged)
    section = {
        "modes": modes,
        "absolute_Ah": {k: values[k] for k in ("dn_li_Ah", "dQ_pe_Ah", "dQ_ne_Ah")},
        "side_reaction_interval": list(side["side_reaction_interval"]),
        "side_reaction_note": "share of BOL lithium lost to side reactions; OCV data give a range only",
        "notes": [
            "LAM_PE is weakly constrained when the positive potential is flat over the window (e.g. "
            "LFP) and the end of discharge is limited by the negative electrode; trust the verdict",
            "the final interval adds a reference-curve bias budget to the profile interval",
        ],
        "definitions": {
            "LLI": "1 - n_li(aged) / n_li(BOL)",
            "LAM_PE": "1 - Q_pe(aged) / Q_pe(BOL)",
            "LAM_NE": "1 - Q_ne(aged) / Q_ne(BOL)",
        },
    }
    return section, warnings


def run_modes(opts: Options) -> dict[str, Any]:
    """Degradation modes with profile intervals, bias budget and verdicts."""
    out = Path(opts.out)
    doc = _base("modes", opts)
    try:
        _prepare(opts, doc, need_aged=True)
        fits, pos, neg = _fits(opts)
    except WorkflowError as exc:
        return _usage_guard(doc, exc, out, "modes")
    doc["fits"] = _fit_section(opts, fits, doc)
    doc["degradation"], w = _modes_section(opts, fits, pos, neg)
    doc["warnings"] += w
    lines = _md_header(doc, "degradation modes") + _md_fits(doc["fits"]) + _md_modes(doc["degradation"])
    return _finish(doc, out, "\n".join(lines) + "\n")


def _differential(q: FloatArray, v: FloatArray) -> tuple[dict[str, Any], dict[str, FloatArray]]:
    s = SETTINGS["ica"]
    span = float(np.max(q) - np.min(q))
    qg, vg = resample_uniform(q, v, span / (s["resample_points"] - 1))
    vs = savgol_filter(vg, s["presmooth_window"], s["presmooth_polyorder"])
    V, d = incremental_capacity(
        vs, qg, dv=s["dqdv_step_V"], window=s["dqdv_window"], polyorder=s["dqdv_polyorder"]
    )
    y = np.abs(d)
    dist = max(1, int(round(s["peak_min_separation_V"] / s["dqdv_step_V"])))
    idx, props = find_peaks(
        np.nan_to_num(y), prominence=s["peak_prominence_ratio"] * np.nanmax(y), distance=dist
    )
    Q, dvq = differential_voltage(qg, vs, window=s["dvdq_window"], polyorder=s["dvdq_polyorder"])
    a = np.abs(dvq)
    lo, hi = s["dvdq_interior"]
    u = (Q - Q[0]) / (Q[-1] - Q[0])
    interior = (u >= lo) & (u <= hi)
    qi, ai = Q[interior], a[interior]
    dq = max(1, int(round(s["dvdq_peak_min_separation_frac"] * len(Q))))
    j, pj = find_peaks(ai, prominence=s["peak_prominence_ratio"] * np.max(ai), distance=dq)
    res = {
        "dqdv_peaks": {"voltage_V": V[idx], "height_Ah_per_V": y[idx], "prominence": props["prominences"]},
        "dvdq_peaks": {"capacity": qi[j], "height_V_per_Ah": ai[j], "prominence": pj["prominences"]},
        "n_points": int(len(q)),
    }
    return res, {"V": V, "dqdv": y, "Q": Q, "dvdq": a, "interior": interior}


def _ica_section(opts: Options, doc: dict[str, Any], out: Path) -> dict[str, Any]:
    section: dict[str, Any] = {}
    curves: dict[str, dict[str, FloatArray]] = {}
    for lab, p in (("BOL", opts.bol), ("aged", opts.aged)):
        if not p:
            continue
        q, v = _fullcell(p, opts)
        try:
            section[lab], curves[lab] = _differential(q, v)
        except ValueError as exc:
            section[lab] = {"error": f"differential analysis failed: {exc}"}
            doc["warnings"].append(f"{lab}: differential analysis failed ({exc})")
            continue
        c = curves[lab]
        name = f"ica_{lab}.csv"
        n = max(len(c["V"]), len(c["Q"]))
        cols = [
            np.pad(c["V"], (0, n - len(c["V"])), constant_values=np.nan),
            np.pad(c["dqdv"], (0, n - len(c["dqdv"])), constant_values=np.nan),
            np.pad(c["Q"], (0, n - len(c["Q"])), constant_values=np.nan),
            np.pad(c["dvdq"], (0, n - len(c["dvdq"])), constant_values=np.nan),
        ]
        out.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            out / name,
            np.column_stack(cols),
            delimiter=",",
            header="voltage_V,abs_dQdV_Ah_per_V,capacity,abs_dVdQ_V_per_Ah",
            comments="",
            fmt="%.8g",
        )
        section[lab]["curve_file"] = name
    if (
        "BOL" in section
        and "aged" in section
        and "error" not in section["BOL"]
        and "error" not in section["aged"]
    ):
        shifts = []
        for pv in section["BOL"]["dqdv_peaks"]["voltage_V"]:
            av = np.asarray(section["aged"]["dqdv_peaks"]["voltage_V"])
            if len(av):
                k = int(np.argmin(np.abs(av - pv)))
                if abs(av[k] - pv) <= 0.05:
                    shifts.append({"BOL_V": pv, "aged_V": av[k], "shift_mV": 1e3 * (av[k] - pv)})
        section["dqdv_peak_shifts"] = shifts
    if opts.figures and curves:
        ok = _ica_figure(
            curves, section, out / f"{doc['command']}_ica.png" if doc["command"] != "ica" else out / "ica.png"
        )
        if ok:
            doc.setdefault("files", {})["ica_figure"] = (
                "ica.png" if doc["command"] == "ica" else f"{doc['command']}_ica.png"
            )
    return section


def run_ica(opts: Options) -> dict[str, Any]:
    """dQ/dV and dV/dQ of the full-cell curves with peak positions."""
    out = Path(opts.out)
    doc = _base("ica", opts)
    try:
        if not opts.bol and not opts.aged:
            raise WorkflowError(2, "give at least one full-cell curve (--bol and/or --aged)")
        items: dict[str, Any] = {}
        if opts.bol:
            items["bol"] = _check_fullcell(opts.bol, "BOL", opts, None)
        if opts.aged:
            items["aged"] = _check_fullcell(opts.aged, "aged", opts, None)
        doc["inputs"] = items
        errors, warnings = _collect(items)
        doc["warnings"] += warnings
        if errors:
            raise WorkflowError(3, errors)
    except WorkflowError as exc:
        return _usage_guard(doc, exc, out, "ica")
    doc["differential"] = _ica_section(opts, doc, out)
    lines = _md_header(doc, "differential analysis") + _md_ica(
        {k: v for k, v in doc["differential"].items() if k in ("BOL", "aged")}
    )
    return _finish(doc, out, "\n".join(lines) + "\n")


def run_report(opts: Options) -> dict[str, Any]:
    """One-shot: check, fit, degradation modes (if an aged curve is given) and ICA."""
    out = Path(opts.out)
    doc = _base("report", opts)
    try:
        _prepare(opts, doc, need_aged=False)
        fits, pos, neg = _fits(opts)
    except WorkflowError as exc:
        return _usage_guard(doc, exc, out, "report")
    doc["fits"] = _fit_section(opts, fits, doc)
    if "aged" in fits and opts.model == "capacity":
        doc["degradation"], w = _modes_section(opts, fits, pos, neg)
        doc["warnings"] += w
    elif "aged" in fits:
        doc["warnings"].append("degradation modes need --model capacity; only fit indicators reported")
    doc["differential"] = _ica_section(opts, doc, out)
    if opts.figures and _fit_figure(fits, out / "report.png", opts.title):
        doc.setdefault("files", {})["figure"] = "report.png"
    lines = (
        _md_header(doc, "report")
        + _md_fits({k: v for k, v in doc["fits"].items() if k in fits})
        + _md_modes(doc.get("degradation", {}))
        + _md_ica({k: v for k, v in doc["differential"].items() if k in ("BOL", "aged")})
    )
    return _finish(doc, out, "\n".join(lines) + "\n")


_RUNNERS = {"check": run_check, "fit": run_fit, "modes": run_modes, "ica": run_ica, "report": run_report}


def run_command(command: str, opts: Options) -> dict[str, Any]:
    """Run one command by name; unexpected exceptions become exit code 1."""
    if command not in _RUNNERS:
        raise ValueError(f"unknown command {command!r}; choose from {sorted(_RUNNERS)}")
    try:
        return _RUNNERS[command](opts)
    except Exception as exc:  # noqa: BLE001 - reported to the caller as JSON
        doc = _base(command, opts)
        doc["errors"].append(f"internal error: {type(exc).__name__}: {exc}")
        doc["exit_code"] = 1
        return _finish(doc, Path(opts.out), "\n".join(_md_header(doc, command)) + "\n")


def list_references() -> dict[str, Any]:
    """Built-in literature references and their aliases."""
    aliases: dict[str, list[str]] = {}
    for alias, key in LITERATURE_ALIASES.items():
        aliases.setdefault(key, []).append(alias)
    return {
        key: {"aliases": sorted(a), "description": REFERENCE_NOTES.get(key, ""), "x_range": [0.0, 1.0]}
        for key, a in sorted(aliases.items())
    }
