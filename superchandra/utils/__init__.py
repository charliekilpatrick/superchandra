"""
SuperChandra utilities: CLI options, coordinates, extinction, logging, math, spectrum, and I/O.

This package provides helper modules used by the main pipeline:

- options : CLI argument parser (get_parser).
- validation : Input validation (is_number).
- math : Numerical integration (integrate_trapezoid).
- coord : RA/Dec parsing (parse_coord).
- extinction : A_V and nH (Av_to_nH, get_MW_Av).
- logging : Package logging (get_logger, configure_logging).
- spectrum : Weight files and bolometric correction (analyze_weightfile, get_bolometric_correction).
- io : File I/O (gunzip).
- global_defaults : Pipeline defaults and spectral presets (see that module).
"""
from .options import get_parser
from .validation import is_number
from .math import integrate_trapezoid
from .coord import parse_coord
from .extinction import Av_to_nH, get_MW_Av
from .spectrum import analyze_weightfile, get_bolometric_correction
from .io import gunzip
from .logging import get_logger, configure_logging, PACKAGE_LOG_NAME

__all__ = [
    "get_parser",
    "is_number",
    "integrate_trapezoid",
    "parse_coord",
    "Av_to_nH",
    "get_MW_Av",
    "analyze_weightfile",
    "get_bolometric_correction",
    "gunzip",
    "get_logger",
    "configure_logging",
    "PACKAGE_LOG_NAME",
]
