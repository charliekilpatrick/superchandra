#!/usr/bin/env python3

"""
By C. D. Kilpatrick 2019-05-25

v1.00: 2019-05-25. Base Chandra download, running ciao scripts
v1.01: 2019-06-02. Fixed issues with CC images, missing ASOL, off axis images
v1.02: 2019-06-09. Added options for host extinction and redshift.  Modified
                   script to stop if no images are found.

chandra.py: A set of scripts for downloading and analyzing Chandra/ACIS image
files for detecting or placing upper limits on the presence of emission at a
set of input coordinates.
"""

import warnings
warnings.filterwarnings('ignore')
import os,requests,sys,io,shutil,time,gzip,ftplib,xmltodict,copy,glob
from astropy.io.votable import parse
from astropy import utils
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import unique,Table
from astropy.time import Time
from astropy.io import ascii,fits
from astropy.wcs import WCS
from dateutil.parser import parse as dateparse
import numpy as np
from photutils import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry
from datetime import datetime

if 'CIAO_DIR' in os.environ.keys():
    ciao = os.environ['CIAO_DIR']
else:
    ciao = ''

global_defaults = {
    'ciao': ciao,
    'uri': 'https://cxcfps.cfa.harvard.edu/cgi-bin/'+\
        'cda/footprint/get_vo_table.pl',
    'ftp': 'cda.cfa.harvard.edu',
    'ftp_dir': 'pub/byobsid/{0}/{1}/',
    'radius': 0.2,
    'phot_radius': 2.216,
    'bin': 1
}
sss = {
    'type': 'blackbody',
    'models': []}

# blackbody with T=0.02-2 keV in 0.005 steps
for temp in np.arange(37):
    sss['models'].append({
        'name': 't{0}'.format('%7.4f'%(0.02+0.005*temp)),
        'temperature': 0.02+0.005*temp,
        'weights': {
            'emin': 0.3,
            'emax': 1.0,
            'ewidth': 0.02
        }})

# Soft and hard state black hole
bh = {
    'type': 'blackhole',
    'models':[
        {'name': 'soft',
         'temperature': 1.0,
         'weights': {
            'emin': 0.3,
            'emax': 9.98,
            'ewidth': 0.02
        }},
        {'name': 'hard',
         'gamma': 0.5,
         'weights': {
            'emin': 0.3,
            'emax': 9.98,
            'ewidth': 0.02
        }}]
}

# Color strings for download messages
green = '\033[1;32;40m'
red = '\033[1;31;40m'
end = '\033[0;0m'

# Check if num is a number
def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

# Integrate two arrays using the trapezoid rule
def integrate_trapezoid(x, y):
    # Check arrays are same length
    if len(x)!=len(y):
        return(None)
    else:
        total = 0
        for i in (np.arange(len(x)-1)+1):
            total += 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
        return(total)
    return(None)

def parse_coord(ra, dec):
    # Check if ':' in ra/dec
    if (':' in str(ra) and ':' in str(dec)):
        coord = SkyCoord(ra, dec, unit=(u.hour, u.deg), frame='icrs')
        return(coord)
    else:
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        return(coord)
    return(None)

def gunzip(file):
    newfile = file.replace('.gz','')
    os.system('gunzip -f '+file)
    return(newfile)

def Av_to_nH(Av):

    # Convert Av to nH using the conversion factor in
    # Guver & Ozel, 2009, MNRAS, 400, 2050
    return(2.21e21 * float(Av))

# Get MW Av based on input ra and dec
def get_MW_Av(coord, nH=False):

    url = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?'
    url += 'locstr={ra}+{dec}+equ+j2000'
    url = url.format(ra=coord.ra.degree, dec=coord.dec.degree)
    r = requests.get(url)

    dictionary = xmltodict.parse(r.text)
    value = dictionary['results']['result'][0]['statistics']['meanValueSandF']
    Av = float(value.split()[0]) * 3.1
    # Convert to nH if necessary
    if nH:
        return(Av_to_nH(Av))
    return(Av)

class chandra():
    def __init__(self):

        # Transient information/search criteria
        self.coord = None
        self.before = None
        self.after = None
        self.rawdir = 'raw/'
        self.weightdir = 'weights/'
        self.outtable = 'output.phot'
        self.obstable = None

        self.keepshort = False

        self.mwg_column = 0.
        self.hst_column = 0.

        # Dictionary with data output from post-processing analysis
        # Names for output photometry table
        names = ['temperature', 'flux', 'flux_limit', 's_to_n', 'total',
            'background', 'energy', 'bccorr', 'exposure']
        self.final_phot = Table([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]],
            names=names)[:0].copy()

        self.evt2files = []
        self.bpixfiles = []
        self.asolfiles = []
        self.maskfiles = []

        # Get options
        self.modeltype = 'sss'
        self.options = {'global': global_defaults}

        if self.modeltype.strip() == 'sss':
            self.options['spectrum']=sss

    def add_options(self, parser=None, usage=None):
        import argparse
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,
                conflict_handler='resolve')
        parser.add_argument('ra', default=None, type=str, help='RA')
        parser.add_argument('dec', default=None, type=str, help='DEC')
        parser.add_argument('--makeclean', default=False, action='store_true',
            help='Clean up all output files from previous runs then exit.')
        parser.add_argument('--download', default=False, action='store_true',
            help='Download the raw data files given input ra and dec.')
        parser.add_argument('--before', default=None, type=str,
            metavar='YYYY-MM-DD', help='Date after which we should reject '+\
            'all Chandra/ACIS observations for reduction.')
        parser.add_argument('--after', default=None, type=str,
            metavar='YYYY-MM-DD', help='Date before which we should reject '+\
            'all Chandra/ACIS observations for reduction.')
        parser.add_argument('--clobber', default=False, action='store_true',
            help='Overwrite files when using download mode.')
        parser.add_argument('--extinction', default=0., type=float,
            help='Host exinction (Av) for the observed SN.  MW-like dust is'+\
            ' assumed and relationship in Guver & Ozel 2009.')
        parser.add_argument('--redshift','-z', default=0., type=float,
            help='Host redshift for evaluating extinction.')
        parser.add_argument('--quick','-q', default=False, action='store_true',
            help='Only calculate the first model in a given set.')
        parser.add_argument('--model','-m', default='sss', type=str,
            help='Model spectrum to use (options are sss|bh).')
        parser.add_argument('--search_radius','-s', default=0.1667, type=float,
            help='Override search radius with input value (units are degrees).')
        return(parser)

    def set_spectrum(self, spectrum_type):

        if spectrum_type=='sss':
            self.options['spectrum']=sss
        elif spectrum_type=='bh':
            self.options['spectrum']=bh

    def get_obstable(self, coord, radius, obsid=True):
        url = self.options['global']['uri']
        ra  = coord.ra.degree
        dec = coord.dec.degree
        params = {
            'pos': '{0},{1}'.format(ra, dec),
            'size': radius,
            'grating': 'NONE',
            'inst': 'ACIS-I,ACIS-S'
        }
        if True:
            # Try to get a response from the URL
            r = requests.get(url, params=params)

            # Create a temporary file and then read it into a votable object.
            # This is terrible, but best way I've found to complete this part
            fname = 'test.xml'
            f = open(fname, 'w')
            f.write(r.text)
            f.close()
            votable = parse(fname)
            os.remove(fname)
        #except:
        #    return(None)

        tbdata = votable.get_first_table().to_table()

        if len(tbdata)==0:
            return(tbdata)

        # Typically only care about unique obsid, so default is to clip on that
        if obsid:
            tbdata = unique(tbdata, keys='ObsId')

        coords = SkyCoord(tbdata['RA'], tbdata['Dec'], unit=(u.deg, u.deg),
            frame='icrs')
        separation = coord.separation(coords)

        # Although using a large search radius, check to make sure that the
        # input coordinates are < 5 arcmin off axis
        remove = []
        for i,row in enumerate(copy.copy(tbdata)):
            row_coord = SkyCoord(row['RA'],row['Dec'],
                unit=(u.deg, u.deg), frame='icrs')
            if self.coord.separation(row_coord).arcmin > radius * 60.:
                remove.append(i)
        tbdata.remove_rows(remove)

        # Require that all observations have SI Mode = TE (timed exposure)
        # rather than CC = continuous clocking (and thus not imaging)
        remove = []
        for i,row in enumerate(copy.copy(tbdata)):
            if not isinstance(row['SIMode'], str):
                mode = row['SIMode'].decode('utf-8')
            else:
                mode = row['SIMode']
            if not mode.startswith('TE'):
                remove.append(i)
        tbdata.remove_rows(remove)

        # Check if we need before/after cuts
        remove = []
        if self.before is not None:
            mjd_before = Time(self.before).mjd
            for i,row in enumerate(copy.copy(tbdata)):
                mjd_obs = Time(row['obs_date']).mjd
                if mjd_obs > mjd_before:
                    warning = 'WARNING: {obsid} after the input before: {date}'
                    print(warning.format(obsid=row['ObsId'],
                        date=self.before.strftime('%Y-%m-%d')))
                    remove.append(i)
        if self.after is not None:
            mjd_after = Time(self.after).mjd
            for i,row in enumerate(copy.copy(tbdata)):
                mjd_obs = Time(row['obs_date']).mjd
                if mjd_obs < mjd_after:
                    warning = 'WARNING: {obsid} before the input after: {date}'
                    print(warning.format(obsid=row['ObsId'],
                        date=self.after.strftime('%Y-%m-%d')))
                    remove.append(i)


        tbdata.remove_rows(remove)

        if not self.keepshort:
            mask = tbdata['Exposure']>5.0
            tbdata = tbdata[mask]

        return(tbdata)

    # Make sure all standard output is formatted in the same way with banner
    # messages for each module
    def make_banner(self, message):
        n = 80
        print('')
        print('')
        print(message)
        # Print part below the message in the banner
        print('#' * n)
        print('#' * n)
        print('')
        print('')

    def download_image(self, url, filename, outdir='./', clobber=True):
        fullfilename = outdir + filename
        unpackfullfilename = fullfilename.replace('.gz','')
        message = 'Trying to download {image: <30}'
        sys.stdout.write(message.format(image=unpackfullfilename))
        sys.stdout.flush()
        if os.path.isfile(unpackfullfilename) and not clobber:
            message = '\r' + message
            message += green+' [SUCCESS - ALREADY EXISTS]'+end+'\n'
            sys.stdout.write(message.format(image=unpackfullfilename))
            return(True, unpackfullfilename)
        try:
            # utils.data.download_file can get buggy if the cache is
            # full.  Clear the cache even though we aren't using caching
            # to prevent download method from choking
            utils.data.clear_download_cache()
        except RuntimeError:
            pass
        try:
            dat = utils.data.download_file(url, cache=False,
                show_progress=False, timeout=120)
            shutil.move(dat, fullfilename)
            gunzip(fullfilename)
            #os.chmod(fullfilename, 0775)
            message = '\r' + message
            message += green+' [SUCCESS]'+end+'\n'
            sys.stdout.write(message.format(image=unpackfullfilename))
            return(True, unpackfullfilename)
        except:
            message = '\r' + message
            message += red+' [FAILURE]'+end+'\n'
            sys.stdout.write(message.format(image=unpackfullfilename))
            return(False, None)

    def collect_images(self, obstable):
        # Define the search coord/radius and grab all that
        # correspond to this region

        if not os.path.exists(self.rawdir):
            os.makedirs(self.rawdir)

        for row in obstable:
            # Get obsid from obstable
            obsid = str(row['ObsId'])

            cmd = f'download_chandra_obsid {obsid} evt2,bpix,asol,msk'

            os.system(cmd)

            files = glob.glob(os.path.join(obsid, 'primary','*.fits.gz'))
            files += glob.glob(os.path.join(obsid, 'secondary','*.fits.gz'))

            # Now iterate over files and grab the ones we want
            for file in files:
                if 'evt2' in file:
                    filename = f'raw/acis_{obsid}_evt2.fits.gz'
                    if filename not in self.evt2files:
                        # Validate that the evt2 file has SUM_2X2, ORC_MODE
                        # OCLKPAIR, FEP_CCD variables.
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
                        self.evt2files.append(filename)
                if 'asol1' in file:
                    filename = f'raw/acis_{obsid}_asol1.fits.gz'
                    if filename not in self.asolfiles:
                        # We want to treat every observation as
                        # separate/independent, so set OBI_NUM=0
                        hdu = fits.open(file, mode='update')
                        hdu['ASPSOL'].header['OBI_NUM'] = 0
                        hdu.writeto(filename, overwrite=True,
                            output_verify='silentfix')
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

        # Make sure asol and evt2 files are sorted and matched to each other
        if len(self.asolfiles) == len(self.evt2files):
            self.asolfiles = sorted(self.asolfiles)
            self.evt2files = sorted(self.evt2files)
            self.bpixfiles = sorted(self.bpixfiles)
            self.maskfiles = sorted(self.maskfiles)
            for i,file in enumerate(self.asolfiles):
                cmpfile = file.replace('asol1','evt2')
                if not (cmpfile == self.evt2files[i]):
                    print(cmpfile,'is not',self.evt2files[i])
                    return(False)
            return(True)
        else:
            return(False)


    # Generate a weight file using the ciao make_instmap_weights method
    def make_instmap_weights(self, weightfile, model, galnh, hostnh, z,
        stdout=False):

        if os.path.isfile(weightfile):
            return(True)

        # Get the parameters for instmap from options
        emin = model['weights']['emin']
        emax = model['weights']['emax']
        ewid = model['weights']['ewidth']

        if emax>10.0:
            emax = str(10.)

        # Construct the make_instmap_weights command
        # This is for a blackbody
        spectrum = ''
        if 'temperature' in model.keys():
            cmd = 'make_instmap_weights {weightfile} '
            cmd += '\"xszphabs.host*xsphabs.gal*xsbbody.p1\" '
            cmd += 'paramvals=\"host.nh={hostnh};host.redshift={z};'
            cmd += 'gal.nh={galnh};p1.kt={spectrum}\" '
            cmd += 'emin={emin} emax={emax} ewidth={ewid}'

            spectrum = '%7.4f'%model['temperature']

        # Use power law if gamma is in keys
        elif 'gamma' in model.keys():
            cmd = 'make_instmap_weights {weightfile} '
            cmd += '\"xszphabs.host*xsphabs.gal*powlaw1d.p1\" '
            cmd += 'paramvals=\"host.nh={hostnh};host.redshift={z};'
            cmd += 'gal.nh={galnh};p1.gamma={spectrum}\" '
            cmd += 'emin={emin} emax={emax} ewidth={ewid}'

            spectrum = '%7.4f'%model['gamma']

        # Format with remaining variable parameters
        cmd = cmd.format(weightfile=weightfile, hostnh=hostnh, z=z,
                spectrum=spectrum, galnh=galnh, emin=emin, emax=emax, ewid=ewid)

        if not stdout:
            cmd += ' > /dev/null 2> /dev/null'

        # Now execute the command
        print(cmd)
        os.system(cmd)

        # Check that the file was generated successfully
        if os.path.isfile(weightfile):
            return(True)
        else:
            return(False)

    def merge_obs(self, files, asolfiles, bpixfiles, maskfiles,
        outdir, weightfile, coord=None, stdout=False):

        # Make input file lists
        evt2_files = 'evt2.lis'
        asol_files = 'asol.lis'
        bpix_files = 'bpix.lis'
        mask_files = 'mask.lis'

        if not os.path.exists('tmpdir/'):
            os.makedirs('tmpdir/')

        with open(evt2_files, 'w') as f:
            for file in files:
                f.write("%s\n" % file)

        with open(asol_files, 'w') as f:
            for file in asolfiles:
                f.write("%s\n" % file)

        with open(bpix_files, 'w') as f:
            for file in bpixfiles:
                f.write("%s\n" % file)

        with open(mask_files, 'w') as f:
            for file in maskfiles:
                f.write("%s\n" % file)

        # Basic merge_obs command just requires input files and output dir
        cmd = 'merge_obs @{files} {outdir} bands={bands} --clobber=yes '
        cmd += 'asolfiles=@{asol} badpixfiles=@{bpix} maskfiles=@{mask} '
        cmd += '--cleanup=yes tmpdir=\"tmpdir/\" '
        cmd = cmd.format(files=evt2_files, outdir=outdir, bands=weightfile,
            asol=asol_files, bpix=bpix_files, mask=mask_files)

        # Check if we're regridding the files to a specific reference coord
        if coord:
            cmd += 'refcoord=\"{ra} {dec}\" '
            ra  = coord.to_string(style='hmsdms',sep=':').split()[0]
            dec = coord.to_string(style='hmsdms',sep=':').split()[1]
            cmd = cmd.format(ra=ra, dec=dec)

        if not stdout:
            cmd += ' '

        # Now run command
        tries = 0
        while ((not os.path.exists(outdir + '/merged_evt.fits') or
            not os.path.exists(outdir + 'band1_thresh.expmap')) and tries < 3):
            print(cmd)
            os.system(cmd)
            tries += 1

        # Check that we output the correct files
        if os.path.isfile(outdir + '/merged_evt.fits'):
            # Filter merged_evt.fits into merged file
            if not os.path.isfile(outdir + '/merged.fits'):
                cmd_merge = 'dmcopy {outdir}merged_evt.fits[EVENTS] '
                cmd_merge += '{outdir}merged.fits option=image'
                cmd_merge = cmd_merge.format(outdir=outdir)
                os.system(cmd_merge)
            if not os.path.isfile(outdir +'/merged.fits'):
                return(False)
            else:
                return(True)
        else:
            return(False)

    def do_photometry(self, image, coord, radius=2.216):

        hdu = fits.open(image)
        pixscale = np.abs(hdu[0].header['CDELT1'])

        wcs = WCS(hdu[0].header)

        # Construct apertures for photometry and background
        aperture = SkyCircularAperture(coord, radius * u.arcsec)
        background = SkyCircularAnnulus(coord, (2*radius) * u.arcsec,
            (4*radius) * u.arcsec)

        # Do photometry for counts and backgroundn
        phot_table = aperture_photometry(hdu[0].data, aperture, wcs=wcs)
        back_table = aperture_photometry(hdu[0].data, background, wcs=wcs)

        # Get photometry data and rescale background to area of photometry
        phot = phot_table['aperture_sum'][0]
        back = back_table['aperture_sum'][0] / 12.

        return(phot,back)

    def get_expmap_value(self, expmap, coord, radius=2.216):

        # Open the expmap HDU and get the pixel scale
        hdu = fits.open(expmap)
        pixscale = np.abs(hdu[0].header['CDELT1'])

        wcs = WCS(hdu[0].header)

        # Construct an aperture from the coordinate
        aperture = SkyCircularAperture(coord, radius * u.arcsec)

        # Do photometry on the image
        phot_table = aperture_photometry(hdu[0].data, aperture, wcs=wcs)

        # Get photometry value and rescale to per pixel
        phot = phot_table['aperture_sum'][0]
        print(phot)
        phot_per_pixel = phot * (pixscale * 3600.0)**2 / (np.pi * radius**2)

        return(phot_per_pixel)

    def analyze_weightfile(self, weightfile):

        energy, weights = np.loadtxt(weightfile, unpack=True)

        # Get average energy integrated over weights
        energy_int = integrate_trapezoid(energy, energy*weights)
        weight_int = integrate_trapezoid(energy, weights)
        average_energy = energy_int / weight_int

        return(average_energy)

    def get_bolometric_correction(self, model, spectype='blackbody'):

        if spectype=='blackbody':
            # Assuming that the input temperature is in keV, which is default
            # input for ciao/make_instmap_weights
            conversion = 12.3986006
            wavelengths = 1.0e-2 + 1.0e-2 * np.arange(60000)
            temperature = model['temperature']

            # Get index range corresponding to emin->emax
            # Get index range corresponding to emin->emax
            idx1 = (np.abs(wavelengths -
                (conversion / model['weights']['emax']))).argmin()
            idx2 = (np.abs(wavelengths -
                (conversion / model['weights']['emin']))).argmin()

            # Calculate Planck's law for temp in keV, wavelength in angstrom
            in_band_wave = wavelengths[idx1:idx2]
            in_band = integrate_trapezoid(in_band_wave, in_band_wave**(-5) /
                (np.exp(conversion/(in_band_wave * temperature)) - 1))
            bolometric = integrate_trapezoid(wavelengths, wavelengths**(-5) /
                (np.exp(conversion/(wavelengths * temperature)) - 1))

            bccorr = in_band / bolometric
            return(bccorr)

if __name__ == '__main__':
    # Start timer, create hst123 class obj, parse args
    start = time.time()
    usagestring='superchandra.py ra dec'
    chandra = chandra()

    # Handle ra/dec in command line so there's no ambiguity about declination
    if len(sys.argv) < 3:
        print(usagestring)
        sys.exit(1)
    else:
        ra = sys.argv[1]
        dec = sys.argv[2]
        if (not (is_number(ra) and is_number(dec)) and
           (':' not in ra and ':' not in dec)):
            error = 'ERROR: cannot interpret ra={ra}, dec={dec}.'
            print(error.format(ra=ra, dec=dec))
            sys.exit(1)
        else:
            chandra.coord = parse_coord(ra, dec)
            sys.argv[1] = str(chandra.coord.ra.degree)
            sys.argv[2] = str(chandra.coord.dec.degree)

    parser = chandra.add_options(usage=usagestring)
    options = parser.parse_args()

    # Set the model spectrum type
    chandra.set_spectrum(options.model)

    if options.before is not None:
        chandra.before = dateparse(options.before)
    if options.after is not None:
        chandra.after = dateparse(options.after)

    if options.quick:
        models = chandra.options['spectrum']['models']

        # For quick mode, only execute one model depending on type
        if options.model=='sss':
            chandra.options['spectrum']['models']=[models[12]]
        elif options.model=='bh':
            chandra.options['spectrum']['models']=[models[0]]

    if options.search_radius:
        chandra.options['global']['radius']=options.search_radius

    # Starting banner
    message = 'Starting superchandra.py'
    chandra.make_banner(message)

    # Grab obs table for input coordinate
    message = 'Getting obstable for ra={ra}, dec={dec}'
    chandra.make_banner(message.format(ra=chandra.coord.ra.degree,
        dec=chandra.coord.dec.degree))
    chandra.obstable = chandra.get_obstable(chandra.coord,
        chandra.options['global']['radius'])

    if not chandra.obstable or len(chandra.obstable)==0:
        print('WARNING: no observations to download for input coordinates.')
        print('Exiting...')
        sys.exit()

    # Analysis on observation IDs
    message = 'There are {n} unique observation IDs in the obstable\n'+\
              'TOTAL EXPTIME: {time} ks\n'
    print(message.format(n=len(chandra.obstable),
        time='%7.2f' % np.sum(chandra.obstable['Exposure'])))

    # Download the image files for each unique observation ID
    message = 'Downloading images for unique observation IDs...'
    chandra.make_banner(message)
    if not chandra.collect_images(chandra.obstable):
        message = 'ERROR: could not properly download images'
        print(message)
        sys.exit(1)

    if len(chandra.obstable) == 0:
        message = 'WARNING: there are no images for ra={ra}, dec={dec}.  '
        message += 'Exiting...'
        print(message.format(ra=chandra.coord.ra.degree,
            dec=chandra.coord.dec.degree))
        sys.exit(0)

    # Get metadata
    message = 'Starting analysis of Chandra/ACIS images...'
    chandra.make_banner(message)

    # Galactic nH
    chandra.mwg_column = get_MW_Av(chandra.coord, nH=True) / 1.0e22
    chandra.hst_column = Av_to_nH(options.extinction) / 1.0e22
    message = 'Galactic nH toward ra={ra}, dec={dec} is nH={nH}'
    print(message.format(ra=chandra.coord.ra.degree,
        dec=chandra.coord.dec.degree, nH=chandra.mwg_column,
        host_nH=chandra.mwg_column))

    # Now iterate over spectral models that we want to use
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

        # Make instmap weight file and merge observations
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
            warning = 'WARNING: merge_obs did not successfully complete '
            warning += 'for params={0}'.format(temp)
            print(warning)
            continue

        # Do post-processing analysis on the data
        # 1) Calculate total number of counts and background in merged event map
        (total, background) = chandra.do_photometry(outdir + 'merged.fits',
            chandra.coord, radius=chandra.options['global']['phot_radius'])

        # 2) Get value of exposure map at this location
        exposure = chandra.get_expmap_value(outdir + 'band1_thresh.expmap',
            chandra.coord, radius=chandra.options['global']['phot_radius'])

        # 3) Get the average energy per photon and bolometric correction from
        # the weight file
        energy = chandra.analyze_weightfile(weightfile)

        # 4) Get bolometric correction based on the input temperature
        bccorr = chandra.get_bolometric_correction(model, spectype=spectype)

        # 5) Calculate bolometric flux in erg/s/cm2 from all data
        flux = (total - background) * energy / exposure / bccorr * 1.60218e-9

        # 6) Determine S/N of flux (this is rough, assumes background limited)
        sn = (total - background)/np.sqrt(total + background)

        # 7) Determine limiting 3-sigma flux using Gehrels 1986 eq. 9
        lam = (total + 1) * (1 - 1/(9*(total+1)) + 3.0/(3*np.sqrt(total+1)))**3
        lam_flux = (lam - background) * energy / exposure / bccorr * 1.60218e-9

        chandra.final_phot.add_row((spectrum, flux, lam_flux, sn, total,
            background, energy, bccorr, exposure))

    message = 'Printing final output photometry to table={table}'
    chandra.make_banner(message.format(table=chandra.outtable))
    print(chandra.final_phot)
    chandra.final_phot.write(chandra.outtable,format='ascii.fixed_width',
        overwrite=True)
