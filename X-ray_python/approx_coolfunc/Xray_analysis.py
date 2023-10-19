#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021 - Thomas Pasini
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy as np
import os, sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy import units as u
import astropy.io.fits as fits
import csv
from matplotlib.ticker import FormatStrFormatter

H0=73
Om0=0.3
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

mu = 0.61
mh = 1.67e-24

print('')
print('########################')
print('#### X-RAY ANALYSIS ####')
print('########################')
print('')

print('Currently using cosmology with H0=', H0, 'OmegaM=', Om0)
print('')

print('The first input file, named norm.txt, must contain:')
print(' - min radius in arcsec;')
print(' - max radius in arcsec;')
print(' - radius error in arcsec;')
print(' - value of norm;')
print(' - lower error on norm;')
print(' - upper error on norm;')
print('')

print('The second input file, named kt.txt, must contain:')
print(' - kT in keV;')
print(' - lower error on kT;')
print(' - upper error on kT;')
print('')

datafile = np.loadtxt('norm.txt')
kt = np.loadtxt('kt.txt')

redshift_str = input("When the files are set, please enter the object redshift: ")
redshift = float(redshift_str)

dl = cosmo.luminosity_distance(redshift)
conv = (cosmo.kpc_proper_per_arcmin(redshift)/60)*u.arcmin/u.arcsec
distang = (dl*(3e+24)*u.cm/u.Mpc)/((1+redshift)**2)
print('The luminosity distance is', "%3.3f" % dl.value, 'Mpc')
print('The angular distance is', "%3.3e" % distang.value, 'cm')
print('The kpc/arcsec conversion is', "%3.3f" % conv.value, 'kpc/arcsec')

with open('elec_density.txt', 'w') as f:
    writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
    writer.writerow(['#rmin(kpc)', 'rmax(kpc)', 'Num_density(cm^-3)', 'lower_numdensity','upper_numdensity', 'Density(g/cm^3)'])
    for i,reg in enumerate(datafile):
        rmin = datafile[i,0]*u.arcsec*conv
        rmax = datafile[i,1]*u.arcsec*conv
        rmin_cm = rmin*3e+21*u.cm/u.kpc
        rmax_cm = rmax*3e+21*u.cm/u.kpc
        
        V=(((4*np.pi)/3)*(rmax_cm**3))-(((4*np.pi)/3)*(rmin_cm**3))
        errvol = 2*(datafile[i,2]**3)
        
        ne=np.sqrt(1e+14*(4*np.pi*(datafile[i,3]/u.cm**5)*((distang*(1+redshift))**2))/(0.82*V))
        rho=ne*mu*mh

        relvol=errvol/V
        relnormsu=datafile[i,5]/datafile[i,3]
        relnormgiu=datafile[i,4]/datafile[i,3]

        reltotsu=relnormsu+relvol.value
        reltotgiu=relnormgiu+relvol.value

        absnsu=reltotsu*datafile[i,3]
        absngiu=reltotgiu*datafile[i,3]
        
        writer.writerow(["%3.3f" % rmin.value, "%3.3f" % rmax.value, "%.5f" % ne.value,"%.5f" % absngiu, "%.5f" % absnsu,"%.3e" % rho.value])
        

dens = np.loadtxt('elec_density.txt')

with open('Results.txt', 'w') as g:
    writer = csv.writer(g, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
    writer.writerow(['#rmin(kpc)', 'rmax(kpc)', 'Num_density(cm^-3)', 'lower_numdensity','upper_numdensity', 'Density(g/cm^3)', 'norm', 'Pressure(dy/cm^2)', 'down_press', 'up_press', 't_cool(yr)', 'down_tcool', 'up_tcool'])
    for i, density in enumerate(dens):
        Temp = ((kt[i,0])*1000)*(1.6e-12)/(1.38e-16)
        press = 1.83*dens[i,2]*kt[i,0]*1000*(1.6e-12)
        t_cool=8.5e+10*(dens[i,2]/(1e-3))**(-1)*(Temp/1e+8)**(0.5)
        entropy=kt[i,0]*((dens[i,2])**(-2/3))
        
        errTsu=(kt[i,2]*1e+3*1.6e-12)/(1.38e-16)
        errTgiu=(kt[i,1]*1e+3*1.6e-12)/(1.38e-16)
        
        #Pressure errors
        errdown_press=np.sqrt((((1.83*3.38e-16*Temp)**2)*(dens[i,3]**2))+(((1.83*3.38e-16*dens[i,2])**2)*(errTgiu)**2))
        errup_press=np.sqrt((((1.83*3.38e-16*Temp)**2)*(dens[i,4]**2))+(((1.83*3.38e-16*dens[i,2])**2)*(errTsu)**2))

        #Cooling time errors
        errcoolsu=np.sqrt((((8.5e+10*((dens[i,2]/(1e-3))**-1)*((1/1e+4)*0.5*((Temp**(-0.5)))))**2)*(errTsu**2)) + (((8.5e+10*1e-3*(-1/(dens[i,2]**2))*((Temp/1e+8)**(-0.5)))**2)*(dens[i,4]**2)))
        errcoolgiu=np.sqrt((((8.5e+10*((dens[i,2]/(1e-3))**-1)*((1/1e+4)*0.5*((Temp**(-0.5)))))**2)*(errTgiu**2)) + (((8.5e+10*1e-3*(-1/(dens[i,2]**2))*((Temp/1e+8)**(-0.5)))**2)*(dens[i,3]**2)))
        
        #Entropy errors
        errentropy_up = np.sqrt(((dens[i,2]**(-2/3))*((kt[i,2])**2))  +(((((-2/3)*kt[i,0])*(dens[i,2]**(-5/3)))**2)*((dens[i,4])**2)))
        errentropy_down = np.sqrt(((dens[i,2]**(-2/3))*((kt[i,1])**2))+(((((-2/3)*kt[i,0])*(dens[i,2]**(-5/3)))**2)*((dens[i,3])**2)))
        
        writer.writerow([dens[i,0], dens[i,1], dens[i,2], dens[i,3], dens[i,4],dens[i,5], datafile[i,3], "%3.3e" % press, "%3.3e" % errdown_press, "%3.3e" % errup_press, "%3.3e" % t_cool, "%3.3e" % errcoolgiu, "%3.3e" % errcoolsu, "%3.3e" % entropy, "%3.3e" % errentropy_down, "%3.3e" % errentropy_up])
        
res = np.loadtxt('Results.txt')

#WRITE RESULTS TO FITS FILE FOR BETTER HANDLING
if os.path.exists('Results.fits'):
    os.remove('Results.fits')

a1 = np.array(res[:,0])
a2 = np.array(res[:,1])
a3 = np.array(res[:,2])
a4 = np.array(res[:,3])
a5 = np.array(res[:,4])
a6 = np.array(res[:,5])
a7 = np.array(res[:,6])
a8 = np.array(res[:,7])
a9 = np.array(res[:,8])
a10 = np.array(res[:,9])
a11 = np.array(res[:,10])
a12 = np.array(res[:,11])
a13 = np.array(res[:,12])
a14 = np.array(res[:,13])
a15 = np.array(res[:,14])
a16 = np.array(res[:,15])
col1 = fits.Column(name='min_r', array=a1, format='E', unit = 'kpc')
col2 = fits.Column(name='max_r', array=a2, format='E', unit = 'kpc')
col3 = fits.Column(name='elec_density', array=a3, format='E', unit = 'cm-3')
col4 = fits.Column(name='lower_ed', array=a4, format='E', unit = 'cm-3')
col5 = fits.Column(name='upper_ed', array=a5, format='E', unit = 'cm-3')
col6 = fits.Column(name='density', array=a6, format='E', unit = 'g/cm3')
col7 = fits.Column(name='norm', array=a7, format='E')
col8 = fits.Column(name='pressure', array=a8, format='E', unit = 'dy/cm2')
col9 = fits.Column(name='lower_p', array=a9, format='E', unit = 'dy/cm2')
col10 = fits.Column(name='upper_p', array=a10, format='E', unit = 'dy/cm2')
col11 = fits.Column(name='t_cool', array=a11, format='E', unit = 'yr')
col12 = fits.Column(name='lower_t', array=a12, format='E', unit = 'yr')
col13 = fits.Column(name='upper_t', array=a13, format='E', unit = 'yr')
col14 = fits.Column(name='entropy', array=a14, format='E', unit = 'keV cm2')
col15 = fits.Column(name='lower_ent', array=a15, format='E', unit = 'keV cm2')
col16 = fits.Column(name='upper_ent', array=a16, format='E', unit = 'keV cm2')
cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('Results.fits')


#PLOTS
print('')
askplot = input("Do you want to plot the results? (y/n)")
print('')
if askplot == 'y':
    
    fig = plt.figure(figsize=(10,10))
    fig2 = plt.figure(figsize=(11,11))
    fig3 = plt.figure(figsize=(10,10))
    fig4 = plt.figure(figsize=(10,10))
    fig5 = plt.figure(figsize=(10,10))
    
    
    ###
    ax1 = fig.add_subplot(111)
    #ax1.set_title(r'Temperature profile', fontsize=25)
    ax1.set_xlabel(r'r (kpc)', fontsize=22)
    ax1.set_ylabel(r'kT (keV)', fontsize=22)
    med = (res[:,0]+res[:,1])/2
    ax1.errorbar(med, kt[:,0], xerr = ((med - res[:,0]), (res[:,1] - med)), yerr = (kt[:,1],kt[:,2]), fmt = 'o')
    ax1.set_xscale('log')
    ax1.tick_params(labelsize=22)
#    ax11 = ax1.twiny()
#    ax11.errorbar(med*conv.value, kt[:,0], xerr = ((med*conv.value - datafile[:,0]*conv.value), (datafile[:,1]*conv.value - med*conv.value)), yerr = (kt[:,1],kt[:,2]), fmt = 'o') # Create a dummy plot
#    ax11.set_xscale('log')
#    ax11.set_xlabel(r'r (kpc)', fontsize=20)
#    ax11.tick_params(labelsize=20)
    #ax11.cla()
    fig.savefig('kT.pdf')
    
    ###
    ax2 = fig2.add_subplot(111)
    #ax2.set_title(r'Density profile', fontsize=22)
    ax2.set_xlabel(r'r (kpc)', fontsize=22)
    ax2.set_ylabel(r'Electron density (cm-3)', fontsize=22)
    med = (res[:,0]+res[:,1])/2
    ax2.errorbar(med, res[:,2], xerr = ((med - res[:,0]), (res[:,1] - med)), yerr = (res[:,3],res[:,4]), fmt = 'o')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
#    ax11 = ax2.twiny()
#    ax11.errorbar(med*conv.value, res[:,2], xerr = ((med*conv.value - datafile[:,0]*conv.value), (datafile[:,1]*conv.value - med*conv.value)), yerr = (res[:,3],res[:,4]), fmt = 'o') # Create a dummy plot
#    ax11.set_xscale('log')
#    ax11.set_xlabel(r'r (kpc)', fontsize=20)
#    ax11.tick_params(labelsize=20, which='both')
    ax2.tick_params(labelsize=22)
    fig2.savefig('density.pdf')
    
    ax3 = fig3.add_subplot(111)
    #ax3.set_title(r'Cooling time profile', fontsize=20)
    ax3.set_xlabel(r'r (kpc)', fontsize=22)
    ax3.set_ylabel(r'Cooling time (Gyr)', fontsize=22)
    med = (res[:,0]+res[:,1])/2
    ax3.errorbar(med, res[:,10]/1e+9, xerr = ((med - res[:,0]), (res[:,1] - med)), yerr = (res[:,11]/1e+9,res[:,12]/1e+9), fmt = 'o')
    ax3.hlines(7.7, 0, 300, colors='red', linestyles='dashed')
    ax3.text(2, 8.1, '7.7 Gyr', fontsize = 21)
    #Lin. reg.
    xclem= np.linspace(1.4,250,100)
    yclem= 0.572693*(xclem**0.734208)
    ax3.plot(xclem,yclem, label='Best fit', color = 'blue')
    ax3.set_xlim([1.3, 210])
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.tick_params(labelsize=22)
    fig3.savefig('tcool.pdf')
    
    ax4 = fig4.add_subplot(111)
    #ax4.set_title(r'Pressure profile', fontsize=20)
    ax4.set_xlabel(r'r (kpc)', fontsize=22)
    ax4.set_ylabel(r'Pressure (10$^{-11}$ dy/cm$^2$)', fontsize=22)
    med = (res[:,0]+res[:,1])/2
    ax4.errorbar(med, res[:,7]/1e-11, xerr = ((med - res[:,0]), (res[:,1] - med)), yerr = (res[:,8]/1e-11,res[:,9]/1e-11), fmt = 'o')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
#    ax11 = ax4.twiny()
#    ax11.errorbar(med*conv.value, res[:,7]/1e-11, xerr = ((med*conv.value - datafile[:,0]*conv.value), (datafile[:,1]*conv.value - med*conv.value)), yerr = (res[:,8]/1e-11,res[:,9]/1e-11), fmt = 'o') # Create a dummy plot
#    ax11.set_xscale('log')
#    ax11.set_xlabel(r'r (kpc)', fontsize=22)
#    ax11.tick_params(labelsize=20)
    ax4.tick_params(labelsize=20)
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax4.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
    fig4.savefig('press.pdf')
    
    ax5 = fig5.add_subplot(111)
    ax5.set_xlabel(r'r (kpc)', fontsize=22)
    ax5.set_ylabel(r'Entropy (keV cm$^2$)', fontsize=22)
    med = (res[:,0]+res[:,1])/2
    ax5.errorbar(med, res[:,13], xerr = ((med - res[:,0]), (res[:,1] - med)), yerr = (res[:,14],res[:,15]), fmt = 'o')
    ax5.set_xscale('log')
    ax5.tick_params(labelsize=22)
#    ax11 = ax1.twiny()
#    ax11.errorbar(med*conv.value, res[:,13], xerr = ((med*conv.value - datafile[:,0]*conv.value), (datafile[:,1]*conv.value - med*conv.value)), yerr = (res[:,14],res[:,15]), fmt = 'o') # Create a dummy plot
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.tick_params(labelsize=22, which='major')
    ax5.tick_params(labelsize=15, which='minor')
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax5.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
    #ax5.set_ylim(9,200)
    #ax5.set_xlim(0.9,200)
    #ax11.cla()
    fig5.savefig('Entropy.pdf')

    

os.remove('Results.txt')
os.remove('elec_density.txt')
