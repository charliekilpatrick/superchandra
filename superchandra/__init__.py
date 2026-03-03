"""
SuperChandra: download and analyze Chandra/ACIS data for a sky position.

This package provides the main pipeline class and entry point for Chandra/ACIS
observation lookup (Chandra Data Archive), data download (evt2, asol, bpix, mask),
CIAO merge_obs, and aperture photometry with configurable spectral models for
flux and upper limits.

Subpackages / modules
--------------------
superchandra : Main module; defines `Chandra` and `main`.
utils : Helpers for options, coordinates, extinction, logging, math, spectrum, I/O.

See Also
--------
superchandra.superchandra.Chandra : Pipeline class.
superchandra.superchandra.main : CLI entry point.
"""
from ._version import __version__
from .superchandra import Chandra, main

__all__ = ["__version__", "Chandra", "main"]
