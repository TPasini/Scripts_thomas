import numpy as np
from astropy.cosmology import FlatLambdaCDM
import csv
import sys, os

#if len(sys.argv) == 0:
#    print('flux2lum.py alpha (defined positive)')
#    sys.exit(0)

cosmo = FlatLambdaCDM(H0=71, Om0=0.27)
print('')
print('Currently using Flat LmbdaCDM H0=0.71 Om0=0.27.')

def kcorr(flux, z, alpha):
 return flux*(1+z)**(alpha-1)

print('')
print('Data are read from data.txt, which must contain, in this order, flux in Jy, redshift and, if needed, flux error.')
print('')
print('Results are stored in Results.txt.')
print('')
data=np.loadtxt('data.txt')

flux=data[:,0]*1e-3
z=data[:,1]
# try: alpha=float(sys.argv[1])
# except: alpha=1
try: alpha=-float(sys.argv[3])
except: alpha=1
try: flux_err=data[:,2]*1e-3
except: flux_err=0

dist=cosmo.luminosity_distance(z).value #in Mpc

os.system('rm Results.txt')

with open ('Results.txt', 'a') as f:
 writer=csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
 writer.writerow(['Lum(W/Hz)', '', 'Lum_err(W/Hz)', '','Lum(erg/s/Hz)', '','Lum_err(erg/s/Hz)'])
 writer.writerows(zip(kcorr(flux,z,alpha) * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux_err,z,alpha) * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux,z,alpha) * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux_err,z,alpha) * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2)))
 
os.system('open Results.txt')
