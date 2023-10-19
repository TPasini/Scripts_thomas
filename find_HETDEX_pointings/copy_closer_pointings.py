#!/usr/bin/env python3

# given coordinates, get all the pointings within a certain distance

import os, sys, argparse
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, FK5
from astropy import units as u
import shutil
from pathlib import Path
import pathlib
import astropy.io.fits as pyfits

gridfile = 'allsky-grid_minor.fits'

def get_data(fname,colname):
    data=pyfits.open(fname)
    data=data[1].data
    return data[colname]

def copy_pointings(cl_list,  newdir='Pointings', obsonly=False, maxdist=2):
    """
    Find pointings which cover a given RA-DEC and copy them in a new directory
    """

    # load all pointings
    t = Table.read(gridfile, format='fits')
    if obsonly:
        t = t[t['hrs'] > 0]
        
    clusters = cl_list
    cl_name = get_data(clusters, 'Name')
    ra = get_data(clusters, 'RA')
    dec = get_data(clusters, 'DEC')

    ras = t['ra']
    decs = t['dec']
    coord = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    # find distances
    t['dist'] = coord.separation(SkyCoord(ras*u.deg,decs*u.deg,frame='fk5'))
    print('')
    print('The following pointings cover the listed clusters:')
    print('')
    for pointing in t[t['dist'] < maxdist*u.deg]:
        print('%s (dist: %.2f deg) ---> Copying into target dir..' % (pointing['name'], pointing['dist']))
        
        paths=[Path('/home/bax8338/node36'), Path('/home/bax8338/node37'), Path('/home/bax8338/node38'), Path('/home/bax8338/node39')]
        for path in paths:
            list_dirs = [_d for _d in path.iterdir() if _d.is_dir()]
    	
            for dir in list_dirs:
            	str_dir = dir.as_posix()
            	if pointing['name'] == str_dir[-7:]:
                	if not os.path.exists(str(newdir)+'/'+pointing['name']):
                		shutil.copytree(str(path)+'/'+pointing['name']+'/ddcal', str(newdir)+'/'+pointing['name']+'/ddcal')
                		shutil.copytree(str(path)+'/'+pointing['name']+'/mss-avg', str(newdir)+'/'+pointing['name']+'/mss-avg')
                    
                    

parser = argparse.ArgumentParser(description='Get pointings close to certain coords')
#parser.add_argument('-r','--ra', dest='ra', action='store', default='', type=float, help='Ra in degrees.')
#parser.add_argument('-d','--dec', dest='dec', action='store', default='', type=float, help='Dec in degrees.')
parser.add_argument('-l','--list', dest='list', action='store', default='', type=str, help='Name of .fits file which lists Cluster Name, RA and DEC.')
parser.add_argument('-n','--newdir', dest='newdir', action='store', default='Pointings', type=pathlib.Path, help='Path and name of the new directory to which pointings are copied (default: "Pointings" in working directory).')
parser.add_argument('-o','--obsonly', dest='obsonly', action='store_true', help='Restrict to only observed pointings (default: False).')
parser.add_argument('-m','--maxdist', dest='maxdist', action='store', default='2', type=int, help='Max distance in degrees (defailt: 2 deg).')
args = parser.parse_args()

copy_pointings(args.list, args.newdir, args.obsonly, args.maxdist)

