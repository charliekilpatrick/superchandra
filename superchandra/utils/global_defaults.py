"""
Central definitions for superchandra default parameters and spectral presets.

All hard-coded pipeline defaults (CDA URI, search radius, photometry aperture,
spectral model grids) live here so they can be adjusted and documented in one place.

Attributes
----------
CIAO_DIR : str
    Path to CIAO installation (from environment variable CIAO_DIR).
GLOBAL_DEFAULTS : dict
    Default URI, radius, phot_radius, bin, ciao for the pipeline.
MIN_EXPOSURE_KS : float
    Minimum exposure in ks for an observation to be kept (TE mode).
sss, bh, frb : dict
    Spectral presets: "type" and "models" for make_instmap_weights and flux.
"""
import os

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
# Path to CIAO installation; used when building commands (make_instmap_weights,
# merge_obs, etc.). Set via environment variable CIAO_DIR.
CIAO_DIR = os.environ.get("CIAO_DIR", "")

# ---------------------------------------------------------------------------
# Global pipeline defaults (observation lookup, photometry, merge)
# ---------------------------------------------------------------------------
# uri: Chandra Data Archive footprint service URL; returns VO table of ACIS observations.
# radius: Default search radius in degrees for finding observations covering the target.
# phot_radius: Default aperture radius in arcsec (background annulus scales from this).
# bin: Image binning factor (reserved for future use).
# ciao: CIAO installation path (same as CIAO_DIR), stored for downstream commands.
# ---------------------------------------------------------------------------

GLOBAL_DEFAULTS = {
    "ciao": CIAO_DIR,
    "uri": "https://cxcfps.cfa.harvard.edu/cgi-bin/cda/footprint/get_vo_table.pl",
    "radius": 0.2,
    "phot_radius": 2.216,
    "bin": 1,
}

# ---------------------------------------------------------------------------
# Observation filtering
# ---------------------------------------------------------------------------
# Minimum exposure in ks for an observation to be kept (TE mode only).
MIN_EXPOSURE_KS = 5.0

# ---------------------------------------------------------------------------
# Spectral model presets (for make_instmap_weights and flux calculation)
# ---------------------------------------------------------------------------
# Each preset has "type" (blackbody, blackhole, frb) and "models": list of dicts.
# Each model dict has "name", "weights" (emin, emax, ewidth in keV), and either
# "temperature" (blackbody) or "gamma" (power-law index).
# ---------------------------------------------------------------------------


def _sss_models():
    """
    Build blackbody grid for soft X-ray sources (e.g. supernova shock cooling).

    Temperature range 0.02--0.20 keV in 0.005 keV steps; band 0.3--1.0 keV.

    Returns
    -------
    list of dict
        Each dict has "name", "temperature", and "weights" (emin, emax, ewidth).
    """
    return [
        {
            "name": "t{:7.4f}".format(0.02 + 0.005 * t),
            "temperature": 0.02 + 0.005 * t,
            "weights": {"emin": 0.3, "emax": 1.0, "ewidth": 0.02},
        }
        for t in range(37)
    ]


# Soft X-ray blackbody grid (e.g. shock cooling); band 0.3--1.0 keV.
sss = {"type": "blackbody", "models": _sss_models()}

# Soft/hard black hole accretion disk presets; band 0.3--9.98 keV.
bh = {
    "type": "blackhole",
    "models": [
        {
            "name": "soft",
            "temperature": 1.0,
            "weights": {"emin": 0.3, "emax": 9.98, "ewidth": 0.02},
        },
        {
            "name": "hard",
            "gamma": 0.5,
            "weights": {"emin": 0.3, "emax": 9.98, "ewidth": 0.02},
        },
    ],
}

# Power-law preset for fast transients (e.g. FRBs); Gamma=2, band 2--10 keV.
frb = {
    "type": "frb",
    "models": [
        {
            "name": "frb",
            "gamma": 2.0,
            "weights": {"emin": 2.0, "emax": 10.0, "ewidth": 0.02},
        }
    ],
}
