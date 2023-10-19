import os, sys
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import lib_coordinates_mode as cm
from lib_fits import flatten
from astropy import wcs as pywcs
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as const
import astropy.io.fits as fits
import lib_radio as radiolib
from astropy.wcs import WCS
import math

file = os.path.basename('efeds_mosaic_flx.fits') #IMAGE NAME

file_cat = 'Coordinates.fits' #COORDINATE LIST NAME

cat = Table.read(file_cat)

#CHANGE DEPENDING ON COORDINATES COLUMNS NAMES:
racat = cat['BCG_RA']
deccat = cat['BCG_DEC']

def calc_noise(data):
    """
    Return the rms of the mask_rms region
    """
    
    rms = np.mean(data)  #IMPORTANT: USE NP.MEAN IF YOU HAVE THE RMS IMAGE, NP.NANSTD IF YOU HAVE THE NORMAL IMAGE
    check = math.isnan(rms)
    if check == False:
        print(rms)
        return rms


head, datafits = flatten(file)
w = WCS(head)
    
for i, coord in enumerate(racat):

    ra, dec = cm.getCoord(racat[i],deccat[i])
    
    # CONVERT COORDINATES INTO PHYSICAL
    rapix, decpix = w.wcs_world2pix(ra, dec, 0)
    
    sizemask=100 #half size of the mask
    
    min_x = int(rapix) - sizemask
    max_x = int(rapix) + sizemask
    min_y = int(decpix) - sizemask
    max_y = int(decpix) + sizemask
    
    restrict = datafits[min_y:max_y, min_x:max_x]

    rms = calc_noise(restrict)

    
    


