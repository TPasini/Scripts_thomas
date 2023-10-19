import os, sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import lib_coordinates_mode as cm
from lib_fits import flatten
import lib_fits as libfits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.io.fits as aif
import lib_radio as radiolib
import aplpy

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

if len(sys.argv) == 1:
    print('Usage: makeimage.py fitsname RA DEC redshift')
    sys.exit(0)

image = sys.argv[1]
RA = float(sys.argv[2])
DEC = float(sys.argv[3])
z = float(sys.argv[4])

conv = cosmo.kpc_proper_per_arcmin(z)/60 #kpc to arcsec conversion

head, data = flatten(image)
w = WCS(head)

fig = aplpy.FITSFigure(data)

fig.show_colorscale(7e-3, stretch='log',  cmap='hot')

rapix, decpix = w.wcs_world2pix(RA, DEC, 0)

fig.recenter(rapix, decpix, radius=200)
fig.add_colorbar()
fig.colorbar.show()
fig.colorbar.set_location('right')
fig.set_xaxis_coord_type('latitude')
#length = (20/conv.value) #Change value depending on N kpc length of the scalebar
#fig.add_scalebar(length)
#fig.scalebar.set_corner('bottom right')
#fig.scalebar.set_color('white')
#fig.scalebar.show()
#fig.add_beam()
#fig.beam.show()
#fig.beam.set_corner('bottom left')
fig.axis_labels.show()
fig.axis_labels.set_xtext('RA')
fig.axis_labels.set_ytext('DEC')
fig.tick_labels.show()
fig.ticks.show()

fig.save('myplot.pdf', dpi=300)

