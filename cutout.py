import os, sys
import numpy as np
from lib_fits import flatten
from astropy import wcs as pywcs
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.io.fits as fits
from astropy import units as u
import argparse

parser = argparse.ArgumentParser(description='Make cutout image of radio map. \n cutout.py <fits image> <RA> <DEC> <size>')
parser.add_argument('image', default=[], help='Fits image to start from.')
parser.add_argument('RA', type=float, help='RA where to center cutout')
parser.add_argument('DEC', type=float, help='DEC where to center cutout')
parser.add_argument('size', default=1, type=float, help='Size of the cutout in deg. Default is 1 deg.')

args = parser.parse_args()

image = args.image
RA = args.RA
DEC = args.DEC
size = args.size

def cutout(image, RA, DEC, size):
    """
    Creates cutout of image centered in RA, DEC with size in deg.
    The header gets rewritten correctly.
    """

    header_surv, data_surv = flatten(image)
    
    pixsize=size*u.deg/float(abs(header_surv['CDELT1']))

    if os.path.exists(image+'-cut.fits'):
        os.remove(image+'-cut.fits')

    w = pywcs.WCS(header_surv)
    coord_source = SkyCoord((RA*u.deg).astype(np.float64), (DEC*u.deg).astype(np.float64), frame='fk5')

    data_cut = Cutout2D(data_surv, coord_source, size=[size*u.deg,size*u.deg], wcs=w, mode='partial', fill_value=np.nan)

    header_surv['CRPIX1'] = int(pixsize.value/2)
    header_surv['CRPIX2'] = int(pixsize.value/2)
    header_surv['CRVAL1'] = RA
    header_surv['CRVAL2'] = DEC

    fits.writeto(image[0:-5]+'-cut.fits', data_cut.data, header = header_surv)
    

cutout(image, RA, DEC, size)
