import numpy as np
from astropy.cosmology import FlatLambdaCDM
import csv
import sys, os

cosmo = FlatLambdaCDM(H0=71, Om0=0.27)
print('')
print('Currently using Flat LmbdaCDM H0=0.71 Om0=0.27.')

#def kcorr(flux, z, alpha):
# return flux*(1+z)**-(1+alpha)

def calc(flux, z):
  return flux*4*np.pi*(dist*3e+24)**2

print('')
print('Data are read from data.txt, which must contain, in this order, flux in erg/s/cm^2, redshift and, if needed, flux error.')
print('')
print('Results are stored in Results.txt.')
print('')
data=np.loadtxt('data.txt')

flux=data[:,0]
z=data[:,1]
try: flux_err=data[:,2]
except: flux_err=0

dist=cosmo.luminosity_distance(z).value #in Mpc

os.system('rm Results.txt')

with open ('Results.txt', 'a') as f:
 writer=csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
 writer.writerow(['Lum(erg/s)', '', 'Lum_err(W/Hz)'])
 writer.writerows(zip(calc(flux, z), calc(flux_err, z)))
 
os.system('open Results.txt')
