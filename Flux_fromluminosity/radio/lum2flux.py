import numpy as np
from astropy.cosmology import FlatLambdaCDM
import csv
import sys, os

if len(sys.argv) == 1:
    print('flux2lum.py redshift lum(in W/Hz) [alpha=1]')
    sys.exit(0)

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
print('')
print('Currently using Flat LmbdaCDM H0=0.70 Om0=0.3.')


z = float(sys.argv[1]) # redshift
lum = float(sys.argv[2]) # in W/Hz
try: alpha = float(sys.argv[3])
except: alpha = -1

print('z=',z)
print('lum=',lum,' W/Hz')
print('alpha=',alpha)

dist = cosmo.luminosity_distance(z).value # in Mpc
print('Distance: ', dist, 'Mpc')

flux= (((float(lum))*(1e+7))/(4*3.1415*((1+z)**(alpha-1))*(dist*3.086*1e+24)**2))/(1e-23)

print('Flux: ',flux*1e+3, 'mJy')
