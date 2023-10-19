import numpy as np
from astropy.cosmology import FlatLambdaCDM
import csv
import sys, os

cosmo = FlatLambdaCDM(H0=71, Om0=0.27)
print('')
print('Currently using Flat LmbdaCDM H0=0.71 Om0=0.27.')

data=np.loadtxt('lls_deg.txt')

angs = data[:,0]*3600
z=data[:,1]
try: angs_err=data[:,2]*3600
except: angs_err=0


conv=cosmo.kpc_proper_per_arcmin(z).value #in Mpc

os.system('rm lls.txt')

with open ('lls.txt', 'a') as f:
 writer=csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
 writer.writerow(['LLS(kpc)', '', '','','','LLS_err(kpc)'])
 writer.writerows(zip((angs/60)*conv, (angs_err/60)*conv))
os.system('open lls.txt')
