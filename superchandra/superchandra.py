#!/usr/bin/env python3
"""
SuperChandra: Chandra/ACIS observation lookup, download, merge, and aperture
photometry with configurable spectral models for flux and upper limits.

Orchestrates Chandra Data Archive (CDA) queries, CIAO-based data download
and merge_obs, instrument map weight files, and aperture photometry; outputs
flux, errors, and upper limits per spectral model (e.g. blackbody, power law).

Notes
-----
Requires CIAO (Chandra Interactive Analysis of Observations) for
make_instmap_weights, merge_obs, download_chandra_obsid, and dmcopy.
Pipeline output uses the ``superchandra.utils.logging`` module (info, warning, error).

References
----------
.. [1] Kilpatrick et al. 2018, MNRAS, 481, 4123 (X-ray limits on SN Ia progenitor).

Author
------
C. D. Kilpatrick
"""
import os
import sys
import shutil
import glob
import tempfile
import warnings
import requests
import xmltodict
import numpy as np
from astropy.io.votable import parse as parse_votable
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import unique, Table
from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
from dateutil.parser import parse as dateparse
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus, aperture_photometry

warnings.filterwarnings("ignore")

from .utils.global_defaults import GLOBAL_DEFAULTS, MIN_EXPOSURE_KS, bh, frb, sss
from .utils.logging import get_logger, configure_logging
from .utils.validation import is_number
from .utils.math import integrate_trapezoid
from .utils.coord import parse_coord
from .utils.extinction import Av_to_nH, get_MW_Av
from .utils.spectrum import analyze_weightfile as _analyze_weightfile
from .utils.spectrum import get_bolometric_correction as _get_bolometric_correction

_log = get_logger("superchandra")

# Re-export utility functions for backward compatibility (e.g. tests, external imports).
__all__ = [
    "is_number",
    "integrate_trapezoid",
    "parse_coord",
    "Av_to_nH",
    "get_MW_Av",
    "Chandra",
    "main",
]


class Chandra:
    """
    Chandra/ACIS download and analysis for source photometry and flux limits.

    Manages observation lookup (CDA), data download (evt2, asol, bpix, mask),
    CIAO merge_obs, instrument map weight files, and aperture photometry for
    a set of spectral models (e.g. blackbody, power law). Holds pipeline state
    (coordinates, date filters, file lists, options) and output tables.
    """

    def __init__(self):
        """
        Initialize pipeline state: search criteria, paths, and output tables.

        Attributes
        ----------
        coord : SkyCoord or None
            Target sky position (set before get_obstable).
        before, after : datetime or None
            Optional date window for observation filtering.
        rawdir, weightdir : str
            Directories for raw data and weight files.
        outtable : str
            Output photometry table path.
        obstable : astropy.table.Table or None
            Observation table from CDA (filled by get_obstable).
        mwg_column, hst_column : float
            Galactic and host column densities (10^22 cm^-2).
        names : list of str
            Column names for final_phot.
        final_phot : astropy.table.Table
            Output photometry (empty until pipeline run).
        evt2files, asolfiles, bpixfiles, maskfiles : list of str
            Paths to downloaded/processed files.
        options : dict
            Pipeline options (global defaults and spectrum preset).
        """
        # Target position and optional date window for observation search.
        self.coord = None
        self.before = None
        self.after = None
        self.rawdir = 'raw/'
        self.weightdir = 'weights/'
        self.outtable = 'output.phot'
        self.obstable = None

        self.keepshort = False

        # Column densities in 10^22 cm^-2 (Galactic and host).
        self.mwg_column = 0.
        self.hst_column = 0.

        # Column names and empty table for final photometry output.
        self.names = ['temperature', 'flux','fluxerr', 'flux_limit', 's_to_n', 'total',
            'background', 'energy', 'bccorr', 'exposure', 
            'count_limit_poisson', 'flux_limit_poisson']
        self.final_phot = Table([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],
            [0.],[0.]],
            names=self.names)[:0].copy()

        # Lists of downloaded file paths (filled by collect_images).
        self.evt2files = []
        self.bpixfiles = []
        self.asolfiles = []
        self.maskfiles = []

        # Pipeline options: global defaults and spectral preset (e.g. sss, bh, frb).
        self.options = {"global": dict(GLOBAL_DEFAULTS), "spectrum": sss}

    def add_options(self, parser=None, usage=None):
        """
        Add or extend a CLI parser with all superchandra options.

        Delegates to `superchandra.utils.options.get_parser`.

        Parameters
        ----------
        parser : argparse.ArgumentParser, optional
            Parser to extend; if None, a new one is created.
        usage : str, optional
            Usage string for the new parser.

        Returns
        -------
        argparse.ArgumentParser
            Parser with superchandra options (ra, dec, --model, --download, etc.).
        """
        from .utils.options import get_parser
        return get_parser(parser=parser, usage=usage)

    def set_spectrum(self, spectrum_type):
        """
        Set the spectral model set (sss, bh, or frb).

        Parameters
        ----------
        spectrum_type : str
            One of 'sss', 'bh', 'frb'. Updates options['spectrum'] and
            output column names (e.g. temperature vs gamma) as needed.
        """
        spectra = {"sss": sss, "bh": bh, "frb": frb}
        self.options["spectrum"] = spectra[spectrum_type]
        if spectrum_type == "frb":
            self.options['spectrum'] = frb
            self.names[0] = 'gamma'
            self.final_phot = Table(
                [[0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.]],
                names=self.names
            )[:0].copy()

    def get_obstable(self, coord, radius, obsid=True):
        """
        Fetch Chandra/ACIS observation table from CDA for a given position and radius.

        Applies cuts: unique ObsId, off-axis < radius, TE mode only, optional
        date (before/after) and minimum exposure.

        Parameters
        ----------
        coord : SkyCoord
            Search center (ICRS).
        radius : float
            Search radius in degrees.
        obsid : bool, optional
            If True, keep unique ObsId only. Default True.

        Returns
        -------
        astropy.table.Table
            Table of matching observations (possibly empty).
        """
        r = requests.get(
            self.options["global"]["uri"],
            params={
                "pos": "{},{}".format(coord.ra.degree, coord.dec.degree),
                "size": radius,
                "grating": "NONE",
                "inst": "ACIS-I,ACIS-S",
            },
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xml", delete=False) as f:
            f.write(r.text)
            fname = f.name
        try:
            votable = parse_votable(fname)
            tbdata = votable.get_first_table().to_table()
        finally:
            os.unlink(fname)

        if len(tbdata) == 0:
            return tbdata

        if obsid:
            tbdata = unique(tbdata, keys='ObsId')

        # Apply off-axis cut: keep only observations whose aimpoint is within radius.
        # (radius is in degrees; separation is compared in arcmin.)
        remove = []
        for i, row in enumerate(tbdata):
            row_coord = SkyCoord(row['RA'], row['Dec'],
                                unit=(u.deg, u.deg), frame='icrs')
            if self.coord.separation(row_coord).arcmin > radius * 60.0:
                remove.append(i)
        if remove:
            tbdata.remove_rows(sorted(set(remove), reverse=True))

        # Keep only TE (timed exposure) mode; drop CC (continuous clocking) and others.
        remove = []
        for i, row in enumerate(tbdata):
            mode = row['SIMode'] if isinstance(row['SIMode'], str) else row['SIMode'].decode('utf-8')
            if not mode.startswith('TE'):
                remove.append(i)
        if remove:
            tbdata.remove_rows(sorted(set(remove), reverse=True))

        # Apply optional date filters: drop observations after self.before or before self.after.
        remove = []
        if self.before is not None:
            mjd_before = Time(self.before).mjd
            for i, row in enumerate(tbdata):
                mjd_obs = Time(row['obs_date']).mjd
                if mjd_obs > mjd_before:
                    _log.warning(
                        '%s after the input before: %s',
                        row['ObsId'], self.before.strftime('%Y-%m-%d'),
                    )
                    remove.append(i)
        if self.after is not None:
            mjd_after = Time(self.after).mjd
            for i, row in enumerate(tbdata):
                mjd_obs = Time(row['obs_date']).mjd
                if mjd_obs < mjd_after:
                    _log.warning(
                        '%s before the input after: %s',
                        row['ObsId'], self.after.strftime('%Y-%m-%d'),
                    )
                    remove.append(i)
        if remove:
            tbdata.remove_rows(sorted(set(remove), reverse=True))

        if not self.keepshort:
            mask = tbdata["Exposure"] > MIN_EXPOSURE_KS
            tbdata = tbdata[mask]

        return tbdata

    def make_banner(self, message):
        """
        Log a section banner: message followed by a line of hash characters.

        Parameters
        ----------
        message : str
            Banner text to log (e.g. pipeline phase description).
        """
        _log.info("\n\n%s\n%s\n", message, "#" * 80)

    def collect_images(self, obstable):
        """
        Download evt2, asol, bpix, and mask files for each ObsId and populate
        self.evt2files, self.asolfiles, self.bpixfiles, self.maskfiles.

        Uses CIAO download_chandra_obsid; copies/patches files into rawdir.

        Parameters
        ----------
        obstable : astropy.table.Table
            Table of observations (must have 'ObsId' column).

        Returns
        -------
        bool
            True if asol and evt2 lists match and are correctly paired.
        """
        if not os.path.exists(self.rawdir):
            os.makedirs(self.rawdir)

        for row in obstable:
            obsid = str(row['ObsId'])

            cmd = f'download_chandra_obsid {obsid} evt2,bpix,asol,msk'

            os.system(cmd)

            files = glob.glob(os.path.join(obsid, 'primary','*.fits.gz'))
            files += glob.glob(os.path.join(obsid, 'secondary','*.fits.gz'))

            # Match downloaded files to evt2, asol1, bpix1, msk1; copy or patch into rawdir.
            for file in files:
                if 'evt2' in file:
                    filename = f'raw/acis_{obsid}_evt2.fits.gz'
                    if filename not in self.evt2files:
                        # Ensure EVENTS header has required CIAO keywords; set defaults if missing.
                        hdu = fits.open(file, mode='update')
                        if 'SUM_2X2' not in hdu['EVENTS'].header.keys():
                            hdu['EVENTS'].header['SUM_2X2'] = 0
                        if 'ORC_MODE' not in hdu['EVENTS'].header.keys():
                            hdu['EVENTS'].header['ORC_MODE'] = 0
                        if 'OCLKPAIR' not in hdu['EVENTS'].header.keys():
                            hdu['EVENTS'].header['OCLKPAIR'] = 8
                        if 'FEP_CCD' not in hdu['EVENTS'].header.keys():
                            hdu['EVENTS'].header['FEP_CCD'] = 275638
                        asolfile = f'acis_{obsid}_asol1.fits'
                        hdu['EVENTS'].header['ASOLFILE'] = asolfile
                        hdu.writeto(filename, overwrite=True,
                                    output_verify='silentfix')
                        hdu.close()
                        self.evt2files.append(filename)
                if 'asol1' in file:
                    filename = f'raw/acis_{obsid}_asol1.fits.gz'
                    if filename not in self.asolfiles:
                        hdu = fits.open(file, mode='update')
                        hdu['ASPSOL'].header['OBI_NUM'] = 0
                        hdu.writeto(filename, overwrite=True,
                                    output_verify='silentfix')
                        hdu.close()
                        self.asolfiles.append(filename)
                if 'bpix1' in file:
                    filename = f'raw/acis_{obsid}_bpix1.fits.gz'
                    if filename not in self.bpixfiles:
                        shutil.copyfile(file, filename)
                        self.bpixfiles.append(filename)
                if 'msk1' in file:
                    filename = f'raw/acis_{obsid}_msk1.fits.gz'
                    if filename not in self.maskfiles:
                        shutil.copyfile(file, filename)
                        self.maskfiles.append(filename)

            shutil.rmtree(obsid)

        # Verify asol and evt2 lists have same length; sort and check 1:1 correspondence by ObsId.
        if len(self.asolfiles) == len(self.evt2files):
            self.asolfiles = sorted(self.asolfiles)
            self.evt2files = sorted(self.evt2files)
            self.bpixfiles = sorted(self.bpixfiles)
            self.maskfiles = sorted(self.maskfiles)
            for i,file in enumerate(self.asolfiles):
                cmpfile = file.replace('asol1', 'evt2')
                if cmpfile != self.evt2files[i]:
                    _log.warning("%s is not %s", cmpfile, self.evt2files[i])
                    return False
            return True
        return False


    def make_instmap_weights(self, weightfile, model, galnh, hostnh, z,
                            stdout=False):
        """
        Generate a CIAO instrument map weight file for a given spectral model.

        Parameters
        ----------
        weightfile : str
            Output path for the weight file.
        model : dict
            Model dict with 'weights' (emin, emax, ewidth) and either
            'temperature' (blackbody) or 'gamma' (power law).
        galnh : float
            Galactic nH (e.g. 10^22 cm^-2).
        hostnh : float
            Host nH (same units).
        z : float
            Redshift.
        stdout : bool, optional
            If False, suppress make_instmap_weights stdout/stderr. Default False.

        Returns
        -------
        bool
            True if weightfile was created or already existed.
        """
        if os.path.isfile(weightfile):
            return True
        w = model["weights"]
        emin, emax, ewid = w["emin"], min(10.0, w["emax"]), w["ewidth"]
        if "temperature" in model:
            spec, spec_str = "xszphabs.host*xsphabs.gal*xsbbody.p1", "p1.kt"
            val = "%7.4f" % model["temperature"]
        else:
            spec, spec_str = "xszphabs.host*xsphabs.gal*powlaw1d.p1", "p1.gamma"
            val = "%7.4f" % model["gamma"]
        cmd = (
            'make_instmap_weights {} "{}" '
            'paramvals="host.nh={};host.redshift={};gal.nh={};{}={}" '
            "emin={} emax={} ewidth={}"
        ).format(weightfile, spec, hostnh, z, galnh, spec_str, val, emin, emax, ewid)
        if not stdout:
            cmd += " > /dev/null 2> /dev/null"
        _log.info("%s", cmd)
        os.system(cmd)
        return os.path.isfile(weightfile)

    def merge_obs(self, files, asolfiles, bpixfiles, maskfiles,
                  outdir, weightfile, coord=None, stdout=False):
        """
        Run CIAO merge_obs to combine event lists and build merged image/expmap.

        Parameters
        ----------
        files : list of str
            evt2 file paths.
        asolfiles : list of str
            asol file paths.
        bpixfiles : list of str
            Bad pixel file paths.
        maskfiles : list of str
            Mask file paths.
        outdir : str
            Output directory (will contain merged_evt.fits, merged.fits, etc.).
        weightfile : str
            Path to instrument map weight file (bands=).
        coord : SkyCoord, optional
            If set, refcoord for regridding.
        stdout : bool, optional
            Whether to show merge_obs output. Default False.

        Returns
        -------
        bool
            True if merged.fits was produced.
        """
        # Write file lists for merge_obs (evt2, asol, bpix, mask).
        os.makedirs("tmpdir/", exist_ok=True)
        list_files = ("evt2.lis", "asol.lis", "bpix.lis", "mask.lis")
        list_data = (files, asolfiles, bpixfiles, maskfiles)
        for lpath, paths in zip(list_files, list_data):
            with open(lpath, "w") as f:
                f.writelines(p + "\n" for p in paths)

        cmd = "merge_obs @evt2.lis {} bands={} --clobber=yes ".format(outdir, weightfile)
        cmd += "asolfiles=@asol.lis badpixfiles=@bpix.lis maskfiles=@mask.lis --cleanup=yes tmpdir=\"tmpdir/\" "

        if coord:
            ra = coord.to_string(style='hmsdms', sep=':').split()[0]
            dec = coord.to_string(style='hmsdms', sep=':').split()[1]
            cmd += 'refcoord="{ra} {dec}" '.format(ra=ra, dec=dec)

        # Run merge_obs (up to 3 tries) until merged_evt and expmap exist.
        tries = 0
        merged_evt = os.path.join(outdir.rstrip('/'), 'merged_evt.fits')
        expmap = os.path.join(outdir.rstrip('/'), 'band1_thresh.expmap')
        merged_fits = os.path.join(outdir.rstrip('/'), 'merged.fits')
        while ((not os.path.exists(merged_evt) or not os.path.exists(expmap)) and tries < 3):
            _log.info("%s", cmd)
            os.system(cmd)
            tries += 1

        # If merge_obs produced merged_evt but not merged.fits, create image with dmcopy.
        if os.path.isfile(merged_evt):
            if not os.path.isfile(merged_fits):
                cmd_merge = 'dmcopy {0}[EVENTS] {1} option=image'.format(
                    merged_evt, merged_fits)
                os.system(cmd_merge)
            return os.path.isfile(merged_fits)
        return False

    def do_photometry(self, image, coord, radius=2.216):
        """
        Aperture photometry (circle) and background (annulus) at coord.

        Parameters
        ----------
        image : str
            Path to a FITS image (e.g. merged.fits).
        coord : SkyCoord
            Sky position for aperture center.
        radius : float, optional
            Aperture radius in arcsec. Background annulus is 2*radius--4*radius.
            Default 2.216.

        Returns
        -------
        tuple of (float, float)
            (source counts, background counts rescaled to aperture area).
        """
        with fits.open(image) as hdu:
            wcs = WCS(hdu[0].header)
            aperture = SkyCircularAperture(coord, radius * u.arcsec)
            background = SkyCircularAnnulus(coord, (2 * radius) * u.arcsec,
                                           (4 * radius) * u.arcsec)
            phot_table = aperture_photometry(hdu[0].data, aperture, wcs=wcs)
            back_table = aperture_photometry(hdu[0].data, background, wcs=wcs)
        phot = phot_table['aperture_sum'][0]
        # Rescale background to aperture area (annulus area ratio ~12 for 2r--4r vs radius r).
        back = back_table['aperture_sum'][0] / 12.0
        return phot, back

    def get_expmap_value(self, expmap, coord, radius=2.216):
        """
        Mean exposure (effective area × time) in an aperture at coord.

        Parameters
        ----------
        expmap : str
            Path to exposure map FITS.
        coord : SkyCoord
            Sky position.
        radius : float, optional
            Aperture radius in arcsec. Default 2.216.

        Returns
        -------
        float
            Exposure value (rescaled to aperture area).
        """
        with fits.open(expmap) as hdu:
            pixscale = np.abs(hdu[0].header['CDELT1'])
            wcs = WCS(hdu[0].header)
            aperture = SkyCircularAperture(coord, radius * u.arcsec)
            phot_table = aperture_photometry(hdu[0].data, aperture, wcs=wcs)
        phot = phot_table['aperture_sum'][0]
        phot_per_pixel = phot * (pixscale * 3600.0) ** 2 / (np.pi * radius ** 2)
        return phot_per_pixel

    def analyze_weightfile(self, weightfile):
        """
        Compute mean photon energy from a CIAO weight file (energy vs weight).

        Parameters
        ----------
        weightfile : str
            Path to two-column file: energy (keV), weight.

        Returns
        -------
        float
            Average energy (integral of E*w / integral of w).
        """
        return _analyze_weightfile(weightfile)

    def get_bolometric_correction(self, model, spectype='blackbody'):
        """
        Bolometric correction: band flux / bolometric flux for the model.

        Only implemented for blackbody; power law returns 1.0.

        Parameters
        ----------
        model : dict
            Must have 'temperature' and 'weights' (emin, emax) for blackbody.
        spectype : str, optional
            'blackbody' or other (returns 1.0). Default 'blackbody'.

        Returns
        -------
        float
            Band-to-bolometric flux ratio.
        """
        return _get_bolometric_correction(model, spectype=spectype)


def main():
    """
    Entry point for the superchandra CLI (ra dec and options).

    Parses command-line arguments, configures logging, runs the full pipeline:
    CDA lookup, download, merge_obs, aperture photometry, and flux/upper-limit
    calculation for each spectral model. Writes final photometry to the
    output table.

    Returns
    -------
    None

    Notes
    -----
    Exits with non-zero status on invalid ra/dec, no observations, or
    download failure. Call `configure_logging()` at start so all output
    uses the package logger.
    """
    configure_logging()
    usagestring = "superchandra ra dec"
    chandra = Chandra()
    parser = chandra.add_options(usage=usagestring)
    options = parser.parse_args()  # --help exits with 0.

    ra, dec = options.ra, options.dec
    if (not (is_number(ra) and is_number(dec))) and (":" not in str(ra) and ":" not in str(dec)):
        _log.error("Cannot interpret ra=%s, dec=%s.", ra, dec)
        sys.exit(1)
    chandra.coord = parse_coord(ra, dec)

    # Apply user-selected spectral preset (sss, bh, or frb).
    chandra.set_spectrum(options.model)

    if options.before is not None:
        chandra.before = dateparse(options.before)
    if options.after is not None:
        chandra.after = dateparse(options.after)

    if options.quick:
        models = chandra.options['spectrum']['models']
        # In quick mode, run only a single representative model per preset.
        if options.model=='sss':
            chandra.options['spectrum']['models']=[models[12]]
        elif options.model=='bh':
            chandra.options['spectrum']['models']=[models[0]]

    if options.search_radius:
        chandra.options['global']['radius']=options.search_radius

    if options.radius:
        chandra.options['global']['phot_radius'] = options.radius

    # Print pipeline start banner.
    message = 'Starting superchandra.py'
    chandra.make_banner(message)

    # Query CDA for observations covering the target position.
    message = 'Getting obstable for ra={ra}, dec={dec}'
    chandra.make_banner(message.format(ra=chandra.coord.ra.degree,
        dec=chandra.coord.dec.degree))
    chandra.obstable = chandra.get_obstable(chandra.coord,
        chandra.options['global']['radius'])

    if not chandra.obstable or len(chandra.obstable)==0:
        _log.warning("No observations to download for input coordinates. Exiting.")
        sys.exit()

    # Summarize obstable (number of ObsIds and total exposure).
    _log.info(
        "There are %d unique observation IDs in the obstable\nTOTAL EXPTIME: %s ks",
        len(chandra.obstable), "%7.2f" % np.sum(chandra.obstable['Exposure']),
    )

    # Download evt2, asol, bpix, and mask files for each ObsId.
    message = 'Downloading images for unique observation IDs...'
    chandra.make_banner(message)
    if not chandra.collect_images(chandra.obstable):
        _log.error("Could not properly download images.")
        sys.exit(1)

    # Print banner for analysis phase.
    message = 'Starting analysis of Chandra/ACIS images...'
    chandra.make_banner(message)

    # Set Galactic and host column densities (nH in 10^22 cm^-2) from IRSA and --extinction.
    chandra.mwg_column = get_MW_Av(chandra.coord, nH=True) / 1.0e22
    chandra.hst_column = Av_to_nH(options.extinction) / 1.0e22
    _log.info(
        "Galactic nH toward ra=%s, dec=%s is nH=%s",
        chandra.coord.ra.degree, chandra.coord.dec.degree, chandra.mwg_column,
    )

    # Ensure weight directory exists before building weight files.
    if not os.path.exists(chandra.weightdir):
        os.makedirs(chandra.weightdir)

    for model in chandra.options['spectrum']['models']:

        weightfile = chandra.weightdir + 'weights.'
        outdir = model['name'] + '/'
        spectrum = 0.
        spectype = ''

        if 'temperature' in model.keys():

            temp = np.round(model['temperature'], decimals=4)
            spectrum = temp
            spectype = 'blackbody'

            message = 'Starting analysis with '
            message += 'temperature={temp}, nH={nH}, host_nH={host_nH}, '
            message += 'host_z={z}'
            chandra.make_banner(message.format(temp=temp,
                nH=chandra.mwg_column, host_nH=chandra.hst_column,
                z=options.redshift))

            weightfile += 't{0}.nH{1}'.format('%7.4f'%temp,
                '%7.4f'%chandra.mwg_column)

        elif 'gamma' in model.keys():

            gamma = np.round(model['gamma'], decimals=4)
            spectrum = gamma
            spectype = 'powerlaw'

            message = 'Starting analysis with '
            message += 'gamma={gamma}, nH={nH}, host_nH={host_nH}, '
            message += 'host_z={z}'
            chandra.make_banner(message.format(gamma=gamma,
                nH=chandra.mwg_column, host_nH=chandra.hst_column,
                z=options.redshift))

            weightfile += 'g{0}.nH{1}'.format('%7.4f'%gamma,
                '%7.4f'%chandra.mwg_column)


        weightfile = weightfile.replace(' ','')
        outdir = outdir.replace(' ','')

        # Build CIAO instrument map weight file and run merge_obs for this model.
        banner = 'Making weight file: {weights}'
        chandra.make_banner(banner.format(weights=weightfile))
        chandra.make_instmap_weights(weightfile, model, chandra.mwg_column,
            chandra.hst_column, options.redshift, stdout=True)

        banner = 'Running merge_obs with weight file: {weights}'
        chandra.make_banner(banner.format(weights=weightfile))

        check = chandra.merge_obs(chandra.evt2files, chandra.asolfiles,
            chandra.bpixfiles, chandra.maskfiles,
            outdir, weightfile, coord=chandra.coord)

        if not check:
            _log.warning("merge_obs did not successfully complete for params=%s.", spectrum)
            continue

        # Post-processing: extract counts, exposure, and compute flux and limits.
        # Total source and background counts in merged image (aperture photometry).
        merged_path = os.path.join(outdir, 'merged.fits')
        (total, background) = chandra.do_photometry(merged_path, chandra.coord,
                                                    radius=chandra.options['global']['phot_radius'])

        # Mean exposure (effective area × time) in the aperture from exposure map.
        expmap_path = os.path.join(outdir, 'band1_thresh.expmap')
        exposure = chandra.get_expmap_value(expmap_path, chandra.coord,
                                            radius=chandra.options['global']['phot_radius'])

        # Mean photon energy and band-to-bolometric correction from weight file and model.
        energy = chandra.analyze_weightfile(weightfile)
        bccorr = chandra.get_bolometric_correction(model, spectype=spectype)

        # Net counts to flux (erg/s/cm^2) using energy, exposure, and bolometric correction.
        flux = (total - background) * energy / exposure / bccorr * 1.60218e-9

        # Signal-to-noise (approximate; background-dominated regime).
        sn = (total - background) / np.sqrt(total + background)

        # Upper limit: 3-sigma count upper bound (Gehrels 1986, eq. 9).
        lam = (total + 1) * (1 - 1/(9*(total+1)) + 3.0/(3*np.sqrt(total+1)))**3

        # Background-dominated 3-sigma count and corresponding flux for limit.
        lam_bkg = 3.0 * np.sqrt(total+background)

        lam_flux = (lam - background) * energy / exposure / bccorr * 1.60218e-9
        lam_flux_bkg = lam_bkg * energy / exposure / bccorr * 1.60218e-9

        chandra.final_phot.add_row((spectrum, flux, flux/sn, lam_flux, sn, total,
            background, energy, bccorr, exposure, lam_bkg, lam_flux_bkg))

    # Write final photometry table to disk.
    chandra.make_banner("Printing final output photometry to table={}".format(chandra.outtable))
    _log.info("\n%s", chandra.final_phot)
    chandra.final_phot.write(chandra.outtable, format='ascii.fixed_width',
                             overwrite=True)


if __name__ == '__main__':
    main()
