import aplpy
import numpy
import matplotlib.pyplot as mpl
from astropy.io import fits
from astropy import units as u
import glob
import os
import sys
import numpy as np
from pylab import meshgrid, imshow
from astropy.wcs import WCS
from astropy.wcs import WCS as pywcs
from astropy.cosmology import FlatLambdaCDM
from astroquery.simbad import Simbad
from lib_fits import flatten


def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs

"""
input fits file, coordinates, redshift, zoom radius
"""

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

if len(sys.argv) == 1:
    print('Usage: makeimage.py fitsname RA DEC redshift [zoom radius(arcsec)]')
    sys.exit(0)

image = sys.argv[1]
RA = float(sys.argv[2])
DEC = float(sys.argv[3])
z = float(sys.argv[4])
try: radiusasec = float(sys.argv[5])
except: radiusasec = 50 #zoom on central 50" if radius not given

"""
estimate rms
"""
    
head, datafits = flatten(image)
w = WCS(head)

def calc_noise(data, niter=100, eps=1e-6):
    """
    Return the rms
    """

    rms = 1.;
    oldrms = 1.
    for i in range(niter):
        rms = np.mean(data)  #IMPORTANT: USE NP.MEAN IF YOU HAVE THE RMS IMAGE, NP.NANSTD IF YOU HAVE THE NORMAL IMAGE
        if np.abs(oldrms - rms)/rms < eps:
            return rms

        oldrms = rms

    raise Exception('Noise estimation failed to converge.')

# CONVERT COORDINATES INTO PHYSICAL
rapix, decpix = w.wcs_world2pix(RA, DEC, 0)

sizemask=50 #half size of the mask

min_x = int(rapix) - sizemask
max_x = int(rapix) - int(sizemask/2)
min_y = int(decpix) - sizemask
max_y = int(decpix) - int(sizemask/2)

restrict = datafits[min_y:max_y, min_x:max_x]

rms = calc_noise(restrict)


"""
apply flatten
"""

f = aplpy.FITSFigure(image)
#f.show_rgb('frame-irg-003698-4-0180.jpg')
fix_aplpy_fits(f)

"""
get conversions and calculate radius
"""

arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(z).value/60 # 1" = N1 kpc
#kpc_to_arcsec = 20/arcsec_to_kpc # 1 Mpc = N2 arcsec
#kpc_to_degrees = kpc_to_arcsec/3600 # 1 Mpc = N3 degrees

radius = radiusasec*arcsec_to_kpc

"""
set color
"""

f.show_colorscale(cmap='afmhot', stretch='power', vmax=0.02) #radio
#f.show_colorscale(cmap='viridis', stretch='power', smooth=7, vmax=1.5e-8) #X
f.add_colorbar()
f.colorbar.set_width(0.15)
f.colorbar.set_axis_label_text('[Jy beam$^{-1}$]')
#f.colorbar.set_axis_label_text(r'$\alpha_{54 \rm MHz}^{144 \rm MHz}$')
#f.colorbar.set_axis_label_text('[erg s$^{-1}$ cm$^{-2}$]')
f.colorbar.set_font(size=12)
f.colorbar.set_axis_label_font(size=12)

"""
recenter image and crop with a certain radius (all parameters in degrees!)
"""

f.recenter(RA, DEC, radius=(radiusasec/3600.))

"""
plot beam
"""

f.add_beam()
f.beam.set_color('black')
f.beam.set_frame(True)

"""
ticks
"""

f.ticks.set_color('white')
f.ticks.set_linewidth(2)
f.tick_labels.set_font(size='xx-large')
f.axis_labels.set_font(size='xx-large')

"""
scale bar
"""

#bar_length = round((radius/2.5)/10)*10 #default scalebar is 1/5 of the entire image
bar_length = 500
f.add_scalebar((bar_length/arcsec_to_kpc)/3600, '500 kpc', color='white')
f.scalebar.set_font(size='xx-large')

"""
Scale
"""

#text=str(arcsec_to_kpc)[:4]
#f.add_label(0.5,0.97, " 1'' / "+text+" kpc",relative=True, bbox=dict(facecolor='white', edgecolor='black',  pad=5.0))


"""
contours
"""

rms = 1.2e-3 #Jy/beam #input by hand if you desire
levels = np.array([-3,3,6,12,24,48,96,192,384,768,1536,3072])
levels = levels[:]*rms

f.show_contour(levels=levels, colors='white', linewidths=0.6)

"""
plot DS9 region and save plot
"""

f.show_regions('grid_8pix_fil.reg')
#f.show_regions('LBAbr0.reg')

f.save('plot.pdf', dpi=300)

