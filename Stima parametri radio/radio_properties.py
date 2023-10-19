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
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

lightvel=3e+10 #cm s-1
k=1.38e-16 #boltzmann cost

print('')
print('######################################################################')
print('############### RADIO SOURCE PROPERTIES ESTIMATION ###################')
print('######################################################################')
print('')

H0=73
Om0=0.3
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

print('Currently using cosmology with H0=', H0, 'OmegaM=', Om0)
print('')

redshift_str = input("Input redshift: ")
redshift = float(redshift_str)

lumdist=cosmo.luminosity_distance(redshift)

convfactor=cosmo.kpc_proper_per_arcmin(redshift)

flux_str = input("Input flux in mJy: ")
flux = float(flux_str)

fluxerr_str = input("Error: ")
fluxerr = float(fluxerr_str)

majarc_str = input("Input major axis in arcsec: ")
majarc = float(majarc_str)
maj=(majarc/60)*convfactor.value

majarcerr_str = input("Error: ")
majarcerr = float(majarcerr_str)
majerr=(majarcerr/60)*convfactor.value

minarc_str = input("Input minor axis in arcsec: ")
minarc = float(minarc_str)
min=(minarc/60)*convfactor.value

minarcerr_str = input("Error: ")
minarcerr = float(minarcerr_str)
minerr=(minarcerr/60)*convfactor.value

spindex_str = input("Input spectral index (defined positive): ")
spindex = float(spindex_str)

spindexerr_str = input("Error: ")
spindexerr = float(spindexerr_str)

freq_str = input("Frequency in Hz: ")
freq = float(freq_str)

shape_str = input("Oblate (1) or prolate (2) shape? ")
shape = float(shape_str)

if shape == 1:
    third=min
    thirderr=minerr
elif shape == 2:
    third=maj
    thirderr=majerr
    
L=4*np.pi*((lumdist.value*3.086e+24)**2)*(flux*1e-26)*((1+redshift)**(spindex-1))

bright=4*(flux*1e-26)*((lumdist.value*3.086e+24)**2)/(np.pi*(maj*3e+21)*(min*3e+21))

T=bright*(1/(2*k))*((lightvel**2)/(freq**2))

V=(1./6.)*np.pi*(maj*3e+21)*(min*3e+21)*(third*3e+21)/(27e+63)

E=2e+41*((L/1e7)**(4./7.))*(V**(3./7.))

epsilon=E/(V*27e+63)

H=((24./7.)*np.pi*epsilon)**(1./2.)

#ERRORS

deltaL=L*((((fluxerr/flux)**2)+((np.log10(1+redshift)*(np.log10(1+redshift)))*spindexerr**2))**0.5)

deltabright=((((4*((lumdist.value*3.086e+24)**2))/(np.pi*(maj*3e+21)*(min*3e+21)))**2)*((fluxerr*1e-26)**2)+
((((-4*(flux*1e-26)*((lumdist.value*3.086e+24)**2)/(np.pi*((maj*3e+21)**2)*(min*3e+21)))**2)*((majerr*3e+21)**2)) +
((((-4*(flux*1e-26)*((lumdist.value*3.086e+24)**2)/(np.pi*(maj*3e+21)*((min*3e+21)**2)))**2)*((minerr*3e+21)**2)))))**0.5

deltaT=((((1/(2*k))*((lightvel**2)/(freq**2)))**2)*deltabright**2)**0.5

deltaV=((4/3)*np.pi)*(((((min/2)*(third/2)**2)*((majerr/2)**2)) +
((((maj/2)*(min/2)**2)*((thirderr/2)**2))) + ((((maj/2)*(third/2)**2)*((minerr/2)**2))))**0.5)

deltaE=(2e+41)*((((((4./7.)*((L/1e+7)**(-3./7.))*(V**(3./7.)))**2)*((deltaL/1e+7)**2)) +
((((3./7.)*((L/1e+7)**(4./7.))*(V**(-4./7.)))**2)*(deltaV**2)))**0.5)

deltaepsilon=(((deltaE/(V*27e+63))**2)+((-E/(((V*27e+63)**2)))**2)*((deltaV*27e+63)**2))**0.5

deltaH=(((12./7.)*np.pi*epsilon)**(-0.5))*deltaepsilon

print('')
print('############################################')
print('#################  RESULTS  ################')
print('############################################')
print('')

rescaled_bright=bright/1e-23
rescaled_brighterr=deltabright/1e-23
rescaled_field=H/1e-6
rescaled_fielderr=deltaH/1e-6

print('Luminosity:', "%10.3e" % L, '+-', "%10.3e" % deltaL ,'erg/s Hz')
print('Brightness:', "%10.3e" % rescaled_bright, '+-', "%10.3e" % rescaled_brighterr,'Jy')
print('Brightness Temperature:', "%5.2f" % T, '+-', "%5.2f" % deltaT , 'K')
print('Volume:', "%5.2e" % V, '+-', "%5.2e" % deltaV , 'kpc^3')
print('Total energy:', "%5.3e" % E, '+-', "%5.3e" % deltaE ,'erg')
print('Energy density:', "%10.3e" % epsilon, '+-', "%10.3e" % deltaepsilon ,'erg/cm^3')
print('Equipartition Field:', "%5.2f" % rescaled_field, '+-', "%5.2f" % rescaled_fielderr , 'muG')

print('')
print('############################################')
