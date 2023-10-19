from astropy.io import fits
from astropy.table import Table
import numpy as np
import scipy


t = Table.read('fitsnames.txt', format='ascii.commented_header')

def run(imagename):
  hdulist = fits.open(imagename)
  RA = hdulist[0].header['RA']
  DEC = hdulist[0].header['DEC']
  print(image['NAME'], RA, DEC)
  
  hdulist.close()
  
for image in t:
  imagename = image['NAME']+'.fits'
  run(imagename)
  
