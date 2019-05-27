#!/usr/bin/env python3

# chandra.py
#
# CDK v1.00: 2019-05-25. Base Chandra download, running ciao scripts.
#
# DESCRIPTION:
# A set of scripts for downloading and analyzing Chandra/ACIS image files for
# detecting or placing upper limits on the presence of emission at a set of
# input coordinates.

# REQUIREMENTS:
# In addition to the requirements in requirements.txt, this script requires the
# latest version (4.11) of ciao (to use ciao/make_instmap_weights and merge_obs,
# specifically).

import warnings,os,requests,sys,io,shutil,time,gzip,ftplib,xmltodict
from astropy.io.votable import parse
from astropy import utils
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import unique
from astropy.io import ascii
import numpy as np
from datetime import datetime
warnings.filterwarnings('ignore')

global_defaults = {
    'ciao': '/Users/ckilpatrick/scripts/ciao',
    'uri': 'https://cxcfps.cfa.harvard.edu/cgi-bin/'+\
        'cda/footprint/get_vo_table.pl',
    'ftp': 'cda.cfa.harvard.edu',
    'ftp_dir': 'pub/byobsid/{0}/{1}/primary/',
    'radius': 0.2
}
weights = {
    'emin': 0.3,
    'emax': 1.0,
    'ewidth': 0.02
}

# Color strings for download messages
green = '\033[1;32;40m'
red = '\033[1;31;40m'
end = '\033[0;0m'

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

# Get MW Av based on input ra and dec.  Option to calculate it in terms of
# nH using the conversion factor in Guver & Ozel, 2009, MNRAS, 400, 2050
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
        return(2.21e21 * float(Av))
    return(Av)

class chandra():
    def __init__(self):

        # Transient information/search criteria
        self.coord = None
        self.rawdir = 'raw/'
        self.weightdir = 'weights/'
        self.obstable = None

        self.evtfiles = []
        self.asolfiles = []

        # Get options
        self.options = {'global': global_defaults, 'weights': weights}


    def add_options(self, parser=None, usage=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage,
                conflict_handler='resolve')
        parser.add_option('--makeclean', default=False, action='store_true',
            help='Clean up all output files from previous runs then exit.')
        parser.add_option('--download', default=False, action='store_true',
            help='Download the raw data files given input ra and dec.')
        parser.add_option('--before', default=None, type='string',
            metavar='YYYY-MM-DD', help='Date after which we should reject '+\
            'all Chandra/ACIS observations for reduction.')
        parser.add_option('--after', default=None, type='string',
            metavar='YYYY-MM-DD', help='Date before which we should reject '+\
            'all Chandra/ACIS observations for reduction.')
        parser.add_option('--clobber', default=False, action='store_true',
            help='Overwrite files when using download mode.')
        parser.add_option('--ra', default=None, metavar='deg/HH:MM:SS',
            type='string', help='RA of interest.')
        parser.add_option('--dec', default=None, metavar='deg/DD:MM:SS',
            type='string', help='DEC of interest.')
        return(parser)

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
        try:
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
        except:
            return(None)

        tbdata = votable.get_first_table().to_table()

        # Typically only care about unique obsid, so default is to clip on that
        if obsid:
            tbdata = unique(tbdata, keys='ObsId')

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

        # First login to Chandra FTP server
        ftp_url = self.options['global']['ftp']
        ftp = ftplib.FTP(ftp_url)
        ftp.login('anonymous','')

        if not os.path.exists(self.rawdir):
            os.makedir(self.rawdir)

        for row in obstable:
            # Get obsid from obstable
            obsid = row['ObsId']

            # Get the first digit of the obsid
            digit = int(obsid) % 10

            # Construct ftp url
            listdir = self.options['global']['ftp_dir']
            listdir = listdir.format(digit, obsid)

            # Get complete list of files
            files = ftp.nlst(listdir)

            # Now iterate over files and grab the ones we want
            for file in files:
                if 'evt2' in file:
                    url = 'ftp://' + ftp_url + '/' + file
                    filename = 'acis_' + str(obsid) + '_evt2.fits.gz'
                    (success, ufile) = self.download_image(url, filename,
                        outdir=self.rawdir, clobber=False)
                    if success:
                        self.evtfiles.append(ufile)
                if 'asol1' in file:
                    url = 'ftp://' + ftp_url + '/' + file
                    filename = 'acis_' + str(obsid) + '_asol1.fits.gz'
                    (success, ufile) = self.download_image(url, filename,
                        outdir=self.rawdir, clobber=False)
                    if success:
                        self.asolfiles.append(ufile)


    # Generate a weight file using the ciao make_instmap_weights method
    def make_instmap_weights(self, weightfile, temperature, galnh,
        stdout=False):

        # Get the parameters for instmap from options
        emin = self.options['weights']['emin']
        emax = self.options['weights']['emax']
        ewid = self.options['weights']['ewidth']

        # Construct the make_instmap_weights command
        cmd = 'make_instmap_weights {weightfile} \"xsphabs.gal*xsbbody.p1\" '
        cmd += 'paramvals=\"gal.nh={galnh};p1.kt={temp}\" '
        cmd += 'emin={emin} emax={emax} ewidth={ewid} clobber=True'

        # Format with variable parameters
        cmd = cmd.format(weightfile=weightfile, galnh=galnh, temp=temperature,
            emin=emin, emax=emax, ewid=ewid)

        if not stdout:
            cmd += ' > /dev/null 2> /dev/null'

        # Now execute the command
        os.system(cmd)
        print(cmd)

        # Check that the file was generated successfully
        if os.path.isfile(weightfile):
            return(True)
        else:
            return(False)

    def merge_obs(self, files, asolfiles, outdir, weightfile,
        coord=None, stdout=False):

        # Get input file string
        input_files = ','.join(files)
        asol_files = ','.join(asolfiles)

        # Basic merge_obs command just requires input files and output dir
        cmd = 'merge_obs {files} {outdir} bands={bands} asolfiles={asol} '
        cmd = cmd.format(files=input_files, outdir=outdir, bands=weightfile,
            asol=asol_files)

        # Check if we're regridding the files to a specific reference coord
        if coord:
            cmd += 'refcoord=\"{ra} {dec}\" '
            ra  = coord.to_string(style='hmsdms',sep=':').split()[0]
            dec = coord.to_string(style='hmsdms',sep=':').split()[1]
            cmd = cmd.format(ra=ra, dec=dec)

        if not stdout:
            cmd += ' > /dev/null 2> /dev/null'

        # Now run command
        os.system(cmd)
        print(cmd)

        # Check that we output the correct files
        if (os.path.isfile(outdir + '/merged_evt.fits') and
            os.path.isfile(outdir + '/band1_thresh.expmap')):
            return(True)
        else:
            return(False)


if __name__ == '__main__':
    # Start timer, create hst123 class obj, parse args
    start = time.time()
    usagestring='USAGE: chandra.py'
    chandra = chandra()
    parser = chandra.add_options(usage=usagestring)
    options, args = parser.parse_args()

    # Starting banner
    message = 'Starting chandra.py'
    chandra.make_banner(message)

    # Parse the input coordinate
    chandra.coord = parse_coord(options.ra, options.dec)

    # Grab obs table for input coordinate
    message = 'Getting obstable for ra={ra}, dec={dec}'
    chandra.make_banner(message.format(ra=chandra.coord.ra.degree,
        dec=chandra.coord.dec.degree))
    chandra.obstable = chandra.get_obstable(chandra.coord,
        chandra.options['global']['radius'])

    # Analysis on observation IDs
    message = 'There are {n} unique observation IDs in the obstable'
    print(message.format(n=len(chandra.obstable)))

    # Download the image files for each unique observation ID
    message = 'Downloading images for unique observation IDs...'
    chandra.make_banner(message)
    chandra.collect_images(chandra.obstable)

    # Get metadata
    message = 'Starting analysis of Chandra/ACIS images...'
    chandra.make_banner(message)

    # Galactic nH
    galnh = get_MW_Av(chandra.coord, nH=True) / 1.0e22
    message = 'Galactic nH toward ra={ra}, dec={dec} is nH={nH}'
    print(message.format(ra=chandra.coord.ra.degree,
        dec=chandra.coord.dec.degree, nH=galnh))

    # Now iterate over spectral models that we want to use
    temperature = 0.020
    weightfile = chandra.weightdir + 'weights.t{temp}.nH{nH}'
    outdir = 't{temp}/'
    outdir = outdir.format(temp=temperature, nH=galnh)
    weightfile = weightfile.format(temp=temperature, nH=galnh)
    chandra.make_instmap_weights(weightfile, temperature, galnh, stdout=True)
    chandra.merge_obs(chandra.evtfiles, chandra.asolfiles,
        outdir, weightfile, coord=chandra.coord)
