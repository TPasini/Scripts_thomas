#Script to create a mask.fits file from a ds9 .reg file.
#Makes use of lib_img in LiLF.

import lib_img, lib_log
import argparse
import os
import numpy as np

logger = lib_log.logger
#logger_obj = lib_log.Logger('plotting.logger')

parser = argparse.ArgumentParser(description='Produce .fits mask file from DS9 region. \n  mask_fromregion.py <fits image> <DS9 .reg file>')
parser.add_argument('image', default=[], help='Fits image to start from.')
parser.add_argument('mask', default=[], help='DS9 .reg file to use as mask. All pixels inside the mask will be set to 1, everything outside will be 0.')
parser.add_argument('--invert', default=False, help='If False, pixels inside masks will be 1 and everything outside will be 0. If True, invert.')
parser.add_argument('--partial', default=False, help='Set to false if you want to create a proper .mask.fits file. Set True if you only need to mask a region, and leave the rest unchanged.')

args = parser.parse_args()

filename = args.image
region = args.mask
invert = args.invert
partial = args.partial

os.system(f'cp {filename} mask.fits')
tomask = 'mask.fits'

if not args.partial:
    if not args.invert:
        logger.info('Setting pixels in mask = 1...')
        firstmask = lib_img.blank_image_reg(tomask, region, inverse=True)
        logger.info('Setting pixels outside mask = 0...')
        secondmask = lib_img.blank_image_reg(tomask, region, blankval=1)
    else:
        logger.info('Setting pixels in mask = 0...')
        firstmask = lib_img.blank_image_reg(tomask, region, inverse=True, blankval=1)
        logger.info('Setting pixels outside mask = 1...')
        secondmask = lib_img.blank_image_reg(tomask, region)
else:
    logger.info('Setting pixels inside mask = 0...')
    secondmask = lib_img.blank_image_reg(tomask, region, blankval=np.nan)

logger.info('Done.')

