"""
Numerical utilities for the superchandra pipeline.

Provides trapezoid-rule integration used for spectra and CIAO weight-file
quantities (e.g. mean energy, bolometric correction).
"""
import numpy as np


def integrate_trapezoid(x, y):
    """
    Integrate y over x with the trapezoid rule.

    Uses `numpy.trapezoid` (or `numpy.trapz` for older NumPy).

    Parameters
    ----------
    x : array-like
        Abscissa values (same length as y).
    y : array-like
        Ordinate values (same length as x).

    Returns
    -------
    float or None
        Integral of y over x, or None if lengths differ.

    Raises
    ------
    RuntimeError
        If neither numpy.trapezoid nor numpy.trapz is available.
    """
    if len(x) != len(y):
        return None
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    if _trapz is None:
        raise RuntimeError("numpy.trapezoid or numpy.trapz required")
    return float(_trapz(y, x))
