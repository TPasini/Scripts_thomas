import os, sys, argparse, logging
from astropy.io import fits
import numpy as np

parser = argparse.ArgumentParser(description='Make spectral curvature map by summing 2 spectral index maps between different frequencies')
parser.add_argument('images', nargs='+', help='List of spectral index maps to use')
parser.add_argument('--output', dest='output', default='curvature.fits', help='Name of output')

args = parser.parse_args()

# check input
if len(args.images) < 2:
    logging.error('Requires 2 images.')
    sys.exit(1)

allimages = args.images

image1 = fits.open(allimages[0])
image2 = fits.open(allimages[1])

data1 = image1[0].data
#data1 = np.nan_to_num(data1, nan=0)
data2 = image2[0].data
#data2 = np.nan_to_num(data2, nan=0)
head1 = image1[0].header
head2 = image2[0].header

freq1 = head1['FREQLO']
freq2 = head2['FREQLO']

beam1 = [head1['BMAJ'], head1['BMIN']]
beam2 = [head2['BMAJ'], head2['BMIN']]

beamdiff = [x - y for x, y in zip(beam1, beam2)]

if any(x*3600 > 1 for x in beamdiff):
    logging.error('Images must have the same resolution.')
    sys.exit(1)

if freq1 < freq2:
    maplow = data1
    maphigh = data2
else:
    maplow = data2
    maphigh = data1

curvature = -maplow + maphigh

outname = args.output
print('Writing final image...')
fits.writeto(outname, curvature, overwrite=True, header = head1)