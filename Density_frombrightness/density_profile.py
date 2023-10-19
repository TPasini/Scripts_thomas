#DEPROJECTED DENSITY PROFILE FROM SURFACE BRIGTHNESS. YOU NEED THE IMAGE, THE EXPOSURE MAP, THE RMF MATRIX AND, IF YOU WANT, THE BACKGROUND. YOU HAVE TO INPUT TEMPERATURE AND NH.

import numpy as np
import pyproffit
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import os, sys
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

cosmo = FlatLambdaCDM(H0=73, Om0=0.3)

if __name__ == '__main__':
   
    dat=pyproffit.Data(imglink='A1668_05_2.fits',
                       explink='perprofilo_broad_thresh.expmap')
                       
    bin=2. #Bin dimension in arcsec
                       
    #dat.dmfilth()
                       
    prof=pyproffit.Profile(dat,center_choice='peak',maxrad=2.,binsize=bin)

    prof.SBprofile()
    
#    prof.Plot()
#    plt.show()

    redshift=0.06355

    cf = prof.Emissivity(z=redshift,
                         kt=2.22,
                         nh=2.2e-2,
                         rmf='allcluster.rmf',
                         elow=0.5,
                         ehigh=2.)

    depr = pyproffit.Deproject(z=redshift, cf=cf, profile=prof)
    
    conv = (cosmo.kpc_proper_per_arcmin(redshift)/60)*u.arcmin/u.arcsec

    depr.Multiscale(nmcmc=1000, bkglim=30.)

    depr.Density()
    
    depr.SaveAll('dens.fits')
    
    data = Table.read('dens.fits', hdu=2)
    
    #depr.PlotDensity()
    
    plt.scatter(data['RADIUS']*60*conv.value, data['DENSITY'], alpha=0.8)
    plt.errorbar(data['RADIUS']*60*conv.value, data['DENSITY'], xerr=(bin/60)*30, yerr=[data['DENSITY_L'],data['DENSITY_H']], ls=' ', alpha=0.8)

    plt.xlabel('Radius [kpc]', fontsize=20)
    plt.ylabel(r'n$_e$ [cm$^{-3}]$', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')

    plt.show()
