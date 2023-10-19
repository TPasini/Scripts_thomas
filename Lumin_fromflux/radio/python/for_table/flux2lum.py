import numpy as np
from astropy.cosmology import FlatLambdaCDM
import csv
import sys, os
import astropy.io.fits as fits

cosmo = FlatLambdaCDM(H0=70, Om0=0.27)
print('')
print('Currently using Flat LmbdaCDM H0=0.70 Om0=0.27.')

if os.path.exists('Results.txt'):
    os.remove('Results.txt')

def kcorr(flux, z, alpha):
 return flux*(1+z)**(alpha-1)

print('')
print('Data are read from match.txt, which must contain flux in Jy, redshift and, if needed, flux error.')
print('')
print('Results are stored in Results.fits.')
print('')
data=np.loadtxt('match.txt')

flux=data[:,2]
z=data[:,1]
ra=data[:,0]
try: alpha=float(sys.argv[1]) #MUST BE DEFINED POSITIVE! SO E.G. ALPHA=0.7
except: alpha=1
print(alpha)
try: flux_err=data[:,3]
except: flux_err=0

#Comment if there are no flux errors available
for n, er in enumerate(flux_err):
   if er == -0.099:
       flux_err[n] = flux[n]/10 #10% of the total flux if error is not available

dist=cosmo.luminosity_distance(z).value #in Mpc

#os.system('rm Results.txt')

with open ('Results.txt', 'a') as f:
 writer=csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
 #writer.writerow(['Lum(W/Hz)', '', 'Lum_err(W/Hz)', '','Lum(erg/s/Hz)', '','Lum_err(erg/s/Hz)'])
 writer.writerows(zip(kcorr(flux,z,alpha) * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux_err,z,alpha) * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux,z,alpha) * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), kcorr(flux_err,z,alpha) * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2), ra, z))
 
#os.system('open Results.txt')
res=np.loadtxt('Results.txt')

if os.path.exists('Results.fits'):
    os.remove('Results.fits')

a1 = np.array(res[:,0])
a2 = np.array(res[:,1])
a3 = np.array(res[:,2])
a4 = np.array(res[:,3])
a5 = np.array(z)
a6 = np.array(ra)
col1 = fits.Column(name='Lum', array=a1, format='E', unit = 'W/Hz')
col2 = fits.Column(name='Lum_err', array=a2, format='E', unit = 'W/Hz')
col3 = fits.Column(name='Lum_erg', array=a3, format='E', unit = 'W/Hz')
col4 = fits.Column(name='Lum_err_erg', array=a4, format='E', unit = 'W/Hz')
col5 = fits.Column(name='z', array=a5, format='F')
col6 = fits.Column(name='ra', array=a6, format='F')

cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('Results.fits')
