#!/usr/bin/env python
# coding: utf-8

import aplpy
from astropy.io import fits
from astropy import units as u
import os, sys
import numpy as np
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from lib_fits import flatten
import argparse
import logging
import warnings
import lib_log
from termcolor import colored

logger = lib_log.logger
#logger_obj = lib_log.Logger('plotting.logger')

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Make pdf image of fits file. Add contours, regions etc. \n    makeimage.py <fits image> *args')
parser.add_argument('image', default=[], help='Fits image to plot.')
parser.add_argument('--type', default='radio', help='radio/X')
parser.add_argument('--cmap', default='gnuplot2', help='Colormap to use.')
parser.add_argument('-z', '--z', type=float, help='Source redshift. Defaults to A1550.')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center image in deg (if not given, get full extent.)')
parser.add_argument('-s', '--size', default=600, type=float, help='Size of the final image in arcsec. Default to 600 arcsec')
parser.add_argument('-n', '--noise', type=float, help='Manually input noise (Jy/beam).')
parser.add_argument('--show_beam', action='store_true', help='Show image beam.')
parser.add_argument('--sbar', help='Length of size bar in kpc.')
parser.add_argument('--show_scale', action='store_true', help='Show arcsec to kpc scale.')
parser.add_argument('--show_contours', action='store_true', help='Show image contours.')
parser.add_argument('--interval', nargs='+', help='Manually input vmin and vmax.')
parser.add_argument('--stretch', default='linear', help='Stretch to apply to image. Default: linear.')
parser.add_argument('--smooth', default=7, help='Smoothing to apply to the X-ray image. Default: 7')
parser.add_argument('--region', nargs='+', help='DS9 region file to plot.')
parser.add_argument('-o', '--outfile', default='plot', help='Prefix of output image.')


def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs


args = parser.parse_args()

if args.radec == None:
    print(colored('ERROR:', 'red'), 'No input coordinates.')
    sys.exit()

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

filename = args.image
regions = args.region
z = args.z
plottype = args.type
plotcmap = args.cmap
center = args.radec
noise = args.noise
size = args.size
beam = args.show_beam
sbar = args.sbar
scale = args.show_scale
contours = args.show_contours
output = args.outfile
interval = args.interval
imstretch = args.stretch
imsmooth = args.smooth


f = aplpy.FITSFigure(filename)

if plottype=='radio':
    fix_aplpy_fits(f)
    f.recenter(center[0], center[1], radius=(size/3600.))
    if interval:
        if imstretch != 'log':
            f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(interval[0]), vmax=float(interval[1]))
        else:
            f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(interval[0]), vmax=float(interval[1]), vmid=float(interval[2]))
    else:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch)
    # f.add_colorbar()
    # f.colorbar.set_width(0.15)
    # f.colorbar.set_axis_label_text('[Jy beam$^{-1}$]')
    # #f.colorbar.set_axis_label_text('Polarisation degree')
    # f.colorbar.set_font(size=18)
    # f.colorbar.set_axis_label_font(size=18)
    # #f.colorbar.set_location('top')

    if contours:
        if not noise:

            head, datafits = flatten(filename)
            w = WCS(head)

            def calc_noise(data, niter=100, eps=1e-6):
                """
                Return the rms
                """

                rms = 1.;
                oldrms = 1.
                for i in range(niter):
                    rms = np.nanstd(data)  #IMPORTANT: USE NP.MEAN IF YOU HAVE THE RMS IMAGE, NP.NANSTD IF YOU HAVE THE NORMAL IMAGE
                    if np.abs(oldrms - rms)/rms < eps:
                        return rms

                    oldrms = rms

                raise Exception('Noise estimation failed to converge.')

            # CONVERT COORDINATES INTO PHYSICAL
            rapix, decpix = w.wcs_world2pix(center[0], center[1], 0)

            sizemask=400 #half size of the mask

            min_x = int(rapix) - sizemask
            max_x = int(rapix) - int(sizemask/2)
            min_y = int(decpix) - sizemask
            max_y = int(decpix) - int(sizemask/2)

            restrict = datafits[min_y:max_y, min_x:max_x]

            rms = calc_noise(restrict)

        else:
            rms = noise

        levels = np.array([3,6,12,24,48,96,192,384,768,1536,3072])
        levels = levels[:]*rms
        logger.info('Plotting contours..')
        f.show_contour(levels=levels, colors='white', linewidths=0.6)

    if beam:
        f.add_beam()
        f.beam.set_color('black')
        f.beam.set_frame(True)

elif plottype=='X':
    f.recenter(center[0], center[1], radius=(size/3600.))
    if interval:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(interval[0]), vmax=float(interval[1]), smooth=float(imsmooth))
    else:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch, smooth=float(imsmooth))
    f.add_colorbar()
    f.colorbar.set_width(0.15)
    f.colorbar.set_axis_label_text('[erg s$^{-1}$ cm$^{-2}$]')
    f.colorbar.set_font(size=12)
    f.colorbar.set_axis_label_font(size=12)

else:
    print(colored('ERROR:', 'red'), 'Plottype unrecognized.')
    sys.exit()


if scale:
    if not z:
        print(colored('ERROR:', 'red'), 'Redshift needed to plot scale.')
        sys.exit()
    arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(z).value/60
    text=str(arcsec_to_kpc)[:4]
    f.add_label(0.5,0.97, " 1'' / "+text+" kpc",relative=True, bbox=dict(facecolor='white', edgecolor='black',  pad=5.0))

if sbar:
    if not z:
        print(colored('ERROR:', 'red'), 'Redshift needed to plot scale.')
        sys.exit()
    bar_length = sbar
    arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(z).value/60
    if float(bar_length) >= 1000:
        f.add_scalebar((float(bar_length) / arcsec_to_kpc) / 3600, str(float(bar_length)/1000) + ' Mpc', color='white')
    else:
        f.add_scalebar((float(bar_length)/arcsec_to_kpc)/3600, bar_length+' kpc', color='white')
    f.scalebar.set_font(size='xx-large')


if regions is not None:
    if isinstance(regions, str):
        regions = [regions]
    for region in regions:
        logger.info('Plotting regions..')
        f.show_regions(region)

f.ticks.set_color('white')
f.ticks.set_linewidth(2)
f.tick_labels.set_font(size='xx-large')
f.axis_labels.set_font(size='xx-large')
f.axis_labels.set_ytext('Declination',)
f.axis_labels.set_xtext('Right Ascension')

if output:
    f.save(output+'.pdf', dpi=300)
else:
    f.save('plot.pdf', dpi=300)

