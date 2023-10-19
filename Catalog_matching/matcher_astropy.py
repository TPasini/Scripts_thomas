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

cat1 = Table.read('COSMOSBGGs.fits')
cat2 = Table.read('vlacosmos.txt', format='ascii')
    
z1 = cat1['redshift_cl']
ra1 = cat1['ra_cl']
dec1 = cat1['dec_cl']
r200 = cat1['R200_deg_cl']

coord_1 = SkyCoord(ra1 * u.deg, dec1 * u.deg, frame='fk5')

conv = cosmo.kpc_proper_per_arcmin(z1)

septhreshold = 1.2*r200 * u.deg

ra2 = cat2['ra']
dec2 = cat2['dec']

coord_2 = SkyCoord(ra2*u.deg, dec2*u.deg)

idx, d2d, d3d = coord_1.match_to_catalog_sky(coord_2)

sep_constraint = d3d < septhreshold.value

c_matches = coord_1[sep_constraint]

catalog_matches = coord_2[idx[sep_constraint]]

separation = coord_1.separation(catalog_matches)

print(separation, septhreshold*3600*u.arcsec/u.deg)
