"""
Spectral utilities for the superchandra pipeline.

Provides mean photon energy from CIAO instrument map weight files and
band-to-bolometric flux correction for blackbody (and power-law) models.
"""
import numpy as np

from .math import integrate_trapezoid


def analyze_weightfile(weightfile):
    """
    Compute mean photon energy from a CIAO weight file (energy vs weight).

    The weight file is a two-column ASCII file: energy in keV and weight.
    Mean energy = integral(E * w) / integral(w) over the file range.

    Parameters
    ----------
    weightfile : str
        Path to two-column file: energy (keV), weight.

    Returns
    -------
    float
        Average energy in keV (integral of E*w / integral of w).
    """
    energy, weights = np.loadtxt(weightfile, unpack=True)
    energy_int = integrate_trapezoid(energy, energy * weights)
    weight_int = integrate_trapezoid(energy, weights)
    return energy_int / weight_int


def get_bolometric_correction(model, spectype='blackbody'):
    """
    Bolometric correction: band flux / bolometric flux for the model.

    For blackbody, integrates the Planck function over the band (model
    "weights" emin, emax) and over full wavelength; returns the ratio so
    that band flux / ratio = bolometric flux. For power law (or any
    non-blackbody), returns 1.0.

    Parameters
    ----------
    model : dict
        Must have 'temperature' and 'weights' with 'emin' and 'emax' (keV)
        for blackbody.
    spectype : str, optional
        'blackbody' or other (returns 1.0). Default is 'blackbody'.

    Returns
    -------
    float
        Band-to-bolometric flux ratio (dimensionless).
    """
    if spectype != 'blackbody':
        return 1.0
    # h*c in keV*Angstrom for Planck function in wavelength.
    conversion = 12.3986006
    wavelengths = 1.0e-2 + 1.0e-2 * np.arange(60000)
    temperature = model['temperature']
    idx1 = (np.abs(wavelengths - (conversion / model['weights']['emax']))).argmin()
    idx2 = (np.abs(wavelengths - (conversion / model['weights']['emin']))).argmin()
    in_band_wave = wavelengths[idx1:idx2]
    in_band = integrate_trapezoid(
        in_band_wave,
        in_band_wave ** (-5) / (np.exp(conversion / (in_band_wave * temperature)) - 1)
    )
    bolometric = integrate_trapezoid(
        wavelengths,
        wavelengths ** (-5) / (np.exp(conversion / (wavelengths * temperature)) - 1)
    )
    return in_band / bolometric
