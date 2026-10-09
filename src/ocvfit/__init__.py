"""ocvfit - half-cell based full-cell OCV fitting and degradation-mode analysis."""

from . import io
from .blend import BlendedElectrode
from .datasets import demo_cells, sample_halfcell, synthesize_discharge
from .degradation import cathode_fixed_indicators, degradation_modes, side_reaction_interval
from .fitting import FitResult, build_weights, fit_capacity, fit_ocv, huber_loss, objective
from .ica import differential_voltage, find_ic_peaks, incremental_capacity
from .literature import LITERATURE_OCP, literature_ocp
from .models import CapacityModel, CathodeFixedModel, ModelComponents, WindowModel
from .ocp import OCPCurve, clean_curve, compress_pchip, monotonize, smooth
from .uncertainty import bias_budget, is_identifiable, lag1_autocorrelation, mode_interval, profile_interval

__version__ = "0.1.0"

__all__ = [
    "BlendedElectrode",
    "CapacityModel",
    "CathodeFixedModel",
    "FitResult",
    "LITERATURE_OCP",
    "ModelComponents",
    "OCPCurve",
    "WindowModel",
    "bias_budget",
    "build_weights",
    "cathode_fixed_indicators",
    "clean_curve",
    "compress_pchip",
    "degradation_modes",
    "demo_cells",
    "differential_voltage",
    "find_ic_peaks",
    "fit_capacity",
    "fit_ocv",
    "huber_loss",
    "incremental_capacity",
    "is_identifiable",
    "lag1_autocorrelation",
    "literature_ocp",
    "mode_interval",
    "io",
    "monotonize",
    "objective",
    "profile_interval",
    "sample_halfcell",
    "side_reaction_interval",
    "smooth",
    "synthesize_discharge",
    "__version__",
]
