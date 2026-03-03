"""
Coordinate parsing for the superchandra pipeline.

Parses RA/Dec (string or float) into an astropy SkyCoord in ICRS for
CLI arguments and pipeline use.
"""
from astropy.coordinates import SkyCoord
from astropy import units as u


def parse_coord(ra, dec):
    """
    Parse RA and Dec into an astropy SkyCoord (ICRS).

    If both ra and dec contain a colon (e.g. "00:50:00", "-05:12:00"),
    they are interpreted as hours and degrees (HM(S)S / DM(S)S).
    Otherwise they are interpreted as decimal degrees.

    Parameters
    ----------
    ra : str or float
        Right ascension: hours if string contains ':', else degrees.
    dec : str or float
        Declination in degrees (or DMS if string contains ':').

    Returns
    -------
    SkyCoord
        ICRS sky coordinates.
    """
    if ':' in str(ra) and ':' in str(dec):
        return SkyCoord(ra, dec, unit=(u.hour, u.deg), frame='icrs')
    return SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
