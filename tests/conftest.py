"""
Pytest configuration and shared fixtures for superchandra tests.

Registers markers (ciao, network) and provides fixtures: chandra instance,
coordinates, temp dir, minimal weight file, minimal FITS image, and minimal expmap.
"""
import os
import tempfile
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

from superchandra.superchandra import Chandra


def pytest_configure(config):
    """Register custom pytest markers (ciao, network)."""
    config.addinivalue_line("markers", "ciao: tests that require CIAO tools (make_instmap_weights, merge_obs, etc.)")
    config.addinivalue_line("markers", "network: tests that require network access (CDA, IRSA)")


@pytest.fixture
def chandra():
    """Chandra pipeline instance with default options."""
    return Chandra()


@pytest.fixture
def coord_deg():
    """SkyCoord in degrees (ICRS)."""
    return SkyCoord(12.5, -5.2, unit=(u.deg, u.deg), frame='icrs')


@pytest.fixture
def coord_hms():
    """SkyCoord from colon-separated RA/Dec."""
    return SkyCoord("00:50:00", "-05:12:00", unit=(u.hour, u.deg), frame='icrs')


@pytest.fixture
def temp_dir(tmp_path):
    """Temporary directory for test outputs."""
    return tmp_path


@pytest.fixture
def minimal_weight_file(tmp_path):
    """Minimal two-column weight file (energy keV, weight) for analyze_weightfile."""
    path = tmp_path / "weights.txt"
    # Linear ramp: energy 0.5--2 keV, weight 1.
    energy = np.linspace(0.5, 2.0, 16)
    weight = np.ones_like(energy)
    np.savetxt(path, np.column_stack([energy, weight]))
    return str(path)


@pytest.fixture
def minimal_fits_image(tmp_path, coord_deg):
    """Minimal FITS image with WCS for do_photometry / get_expmap_value."""
    path = str(tmp_path / "test_image.fits")
    naxis = 50
    # Pixel scale ~1 arcsec (CDELT in deg).
    cdelt = 1.0 / 3600.0
    # Reference pixel at image center; CRVAL set to coord_deg.
    crval1 = coord_deg.ra.deg
    crval2 = coord_deg.dec.deg
    crpix1 = crpix2 = naxis / 2.0 + 0.5
    w = WCS(naxis=2)
    w.wcs.crval = [crval1, crval2]
    w.wcs.crpix = [crpix1, crpix2]
    w.wcs.cdelt = [cdelt, cdelt]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    data = np.zeros((naxis, naxis), dtype=np.float32)
    # Small source at center (few counts) for photometry tests.
    data[naxis // 2 - 2:naxis // 2 + 3, naxis // 2 - 2:naxis // 2 + 3] = 1.0
    header = w.to_header()
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(path, overwrite=True)
    return path


@pytest.fixture
def minimal_expmap(tmp_path, coord_deg):
    """Minimal exposure map FITS (same WCS as minimal_fits_image)."""
    path = str(tmp_path / "test_expmap.fits")
    naxis = 50
    # Same pixel scale as minimal_fits_image (~1 arcsec).
    cdelt = 1.0 / 3600.0
    crval1 = coord_deg.ra.deg
    crval2 = coord_deg.dec.deg
    crpix1 = crpix2 = naxis / 2.0 + 0.5
    w = WCS(naxis=2)
    w.wcs.crval = [crval1, crval2]
    w.wcs.crpix = [crpix1, crpix2]
    w.wcs.cdelt = [cdelt, cdelt]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    # Uniform exposure (constant value per pixel).
    data = np.ones((naxis, naxis), dtype=np.float32) * 1000.0
    header = w.to_header()
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(path, overwrite=True)
    return path
