import numpy as np
import astropy.io.fits as fits
import astropy.coordinates as coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import csv
import multiprocessing
import os, sys

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

cat1 = Table.read('ared.txt', format='ascii')
cat2 = Table.read('bred.txt', format='ascii')

def multiprocessing_func():
    with open('matches.txt', 'a') as g:
        writer = csv.writer(g, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
        for i, c in enumerate(cat1):
            z1 = c['z']
            ra1 = c['RAJ2000_X']
            dec1 = c['DECJ2000_X']
            r200 = c['R200']
            
            coord_1 = SkyCoord(ra1 * u.deg, dec1 * u.deg, frame='fk5')
            
            conv = cosmo.kpc_proper_per_arcmin(z1)
            
            septhreshold = 0.7*r200
            
            with open('temporary.txt', 'a') as f:
                writer2 = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
                for j, g in enumerate(cat2):
                    z2 = g['z']
                    ra2 = g['Ra']
                    dec2 = g['Dec']
                    print(i,j)
                    
                    coord_2 = SkyCoord(ra2 * u.deg, dec2 * u.deg, frame='fk5')
                    
                    sep = coord_1.separation(coord_2)
                    
                    sepkpc =  sep.arcminute*conv.value
                 
                    if sepkpc < septhreshold and abs(z1-z2) < 0.02:
                        writer2.writerow([ra1, dec1, ra2, dec2, z1, z2, r200, sepkpc])
                
                findmin = np.loadtxt('temporary.txt')
                    
                if len(findmin)!=0:
                
                    minsep = np.min(findmin[:,7])
                    
                    for n, obj in enumerate(findmin[:,7]):
                        if obj > minsep:
                            index = 0
                            writer.writerow([ra1, dec1, findmin[n,2], findmin[n,3], z1, findmin[n,5], r200, findmin[n,7], index])
                        else:
                            index = 1
                            writer.writerow([ra1, dec1, findmin[n,2], findmin[n,3], z1, findmin[n,5], r200, findmin[n,7], index])
                
                
        
        if os.path.exists('temporary.txt'):
            os.remove('temporary.txt')
                        
 
if __name__ == '__main__':
    processes = []
    for i in range(1,4):
        p = multiprocessing.Process(target=multiprocessing_func)
        processes.append(p)
        p.start()
        
    for process in processes:
        process.join()

temp=np.loadtxt('matches.txt')
            
if os.path.exists('Matched_catalog.fits'):
    os.remove('Matched_catalog.fits')

a1 = np.array(temp[:,0])
a2 = np.array(temp[:,1])
a3 = np.array(temp[:,2])
a4 = np.array(temp[:,3])
a5 = np.array(temp[:,4])
a6 = np.array(temp[:,5])
a7 = np.array(temp[:,6])
a8 = np.array(temp[:,7])
a9 = np.array(temp[:,8])
col1 = fits.Column(name='RA1', array=a1, format='D', unit = 'deg')
col2 = fits.Column(name='DEC1', array=a2, format='D', unit = 'deg')
col3 = fits.Column(name='RA2', array=a3, format='D', unit = 'deg')
col4 = fits.Column(name='DEC2', array=a4, format='D', unit = 'deg')
col5 = fits.Column(name='z1', array=a5, format='D', unit = '')
col6 = fits.Column(name='z2', array=a6, format='D', unit = '')
col7 = fits.Column(name='R200', array=a7, format='D', unit='kpc')
col8 = fits.Column(name='Sep', array=a8, format='D', unit='kpc')
col9 = fits.Column(name='Index', array=a9, format='D')
cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('Matched_catalog.fits')


