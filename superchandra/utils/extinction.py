"""
Extinction and column density utilities for the superchandra pipeline.

Provides A_V to nH conversion (MW-like dust) and Milky Way A_V/nH from
the IRSA DUST service for a given sky position.
"""
import requests
import xmltodict


def Av_to_nH(Av):
    """
    Convert extinction A_V to hydrogen column density nH (MW-like dust).

    Uses the conversion from Güver & Özel (2009, MNRAS, 400, 2050):
    nH = 2.21e21 * A_V (cm^-2).

    Parameters
    ----------
    Av : float
        V-band extinction in magnitudes.

    Returns
    -------
    float
        Column density in cm^-2.
    """
    return 2.21e21 * float(Av)


def get_MW_Av(coord, nH=False):
    """
    Get Milky Way A_V (or nH) toward a sky position from IRSA DUST.

    Queries the IRSA DUST web service; E(B-V) is converted to A_V using
    R_V = 3.1. Optionally returns nH via Av_to_nH.

    Parameters
    ----------
    coord : SkyCoord
        ICRS sky position (used for query).
    nH : bool, optional
        If True, return column density nH (cm^-2) instead of A_V.
        Default is False.

    Returns
    -------
    float
        A_V in magnitudes, or nH in cm^-2 if nH=True.

    Notes
    -----
    Requires network access to IRSA.
    """
    url = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?'
    url += 'locstr={ra}+{dec}+equ+j2000'
    url = url.format(ra=coord.ra.degree, dec=coord.dec.degree)
    r = requests.get(url)
    dictionary = xmltodict.parse(r.text)
    value = dictionary['results']['result'][0]['statistics']['meanValueSandF']
    # IRSA returns E(B-V); convert to A_V using R_V = 3.1.
    Av = float(value.split()[0]) * 3.1
    if nH:
        return Av_to_nH(Av)
    return Av
