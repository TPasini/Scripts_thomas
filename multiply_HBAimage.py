import sys
import os
import argparse
from astropy.io import fits

parser = argparse.ArgumentParser(description='Multiply fits image by a correction factor. This script was created to align the flux scale of LOFAR images to LoTSS-DR2.')
parser.add_argument('infile', help='Input image.', type=str)
parser.add_argument('factor', help='Multiplicative correction factor.', type=float)
parser.add_argument('-o','--outfile', help='Optional output image name. Default: infile.DR2corrected.fits', type=str, default=None)
parser.add_argument('--compress', help='Compress with gzip the output file. Default: False.', default=False, action='store_true')

args = vars(parser.parse_args())

infile  = args['infile']
factor  = args['factor']
outfile = args['outfile']

hdul = fits.open(infile)
data = hdul[0].data
hdul[0].data = data * factor

if infile.endswith('.fits') or infile.endswith('.gz'):
    prefix = infile.replace('.fits','').replace('.gz','')
else:
    print('ERROR: input file should be .fits or .fits.gz')
    sys.exit()

if outfile == None:
    outfile = prefix + '.DR2corrected.fits'

if args['compress'] == True:
    outfile += '.gz'

hdr = hdul[0].header
hdr['history'] = ''
hdr['history'] = 'Image multiplied by ' + str(factor) + ':'
hdr['history'] = infile + ' * ' + str(factor) + ' --> ' + outfile

hdul.writeto(outfile, overwrite=True)

hdul.close()

print('The multiplied file is: ' + outfile)