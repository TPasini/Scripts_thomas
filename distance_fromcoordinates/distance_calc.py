import numpy as np
import astropy.io.fits as fits
import astropy.coordinates as coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import csv
import multiprocessing
from astropy.coordinates import match_coordinates_sky
import os, sys

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

cat = Table.read('fordistance.txt', format='ascii')
    
ra1 = cat['ra']
dec1 = cat['dec']
ra2 = cat['RAdeg']
dec2 = cat['DEdeg']

coord_1 = SkyCoord(ra1 * u.deg, dec1 * u.deg, frame='fk5')

coord_2 = SkyCoord(ra2*u.deg, dec2*u.deg)

sep = coord_1.separation(coord_2)

with open ('sep.txt', 'w') as f:
    writer=csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
    writer.writerows(zip(sep*u.deg))
