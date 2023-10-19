import numpy as np
import csv
import sys, os, shutil
from astropy.table import Table

os.remove('eFEDS_groups.fits')
os.remove('provisional.txt')

def count(column):
    """Return the number of data in one column that respect the condition"""
    num=0
    for i,value in enumerate(column):
        if column[i] < 1:
            num = num + 1
    return num

#def restrict(column):
#    """Return a file containing only the wanted columns of the starting file that respect the condition"""
#    num=0
#    with open('eFEDS_groups.csv', 'w') as f:
#        writer = csv.writer(f, delimiter = " ")
#        writer.writerow(['RAJ2000', 'DECJ2000', 'z', 'T (keV)', 'M500 (Msun)', 'L500 (10^44 erg/s)'])
#        for i,value in enumerate(column):
#            if column[i] < 1 and column[i] > 0:
#                writer.writerow(["%3.5f" % data[i,0], "%3.5f" % data[i,1], data[i,2], data[i,3], "%4.3e" % data[i,6], "%4.3f" % data[i,9]])
#    return os.system('open eFEDS_groups.csv')
          
def restrict(column):
    """Return a file containing only the wanted columns of the starting file that respect the condition"""
    num=0
    with open('provisional.txt', 'w') as f:
        writer = csv.writer(f, delimiter = " ")
        for i,value in enumerate(column):
            if column[i] < 1 and column[i] > 0:
                writer.writerow(["%3.5f" % data[i,0], "%3.5f" % data[i,1], data[i,2], data[i,3], "%4.3e" % data[i,6], "%4.3f" % data[i,9]])
    provisional = np.loadtxt('provisional.txt')
    t = Table([[provisional[:,0]], [provisional[:,1]], [provisional[:,2]], [provisional[:,3]], [provisional[:,4]], [provisional[:,5]]], names=('RAJ2000', 'DECJ2000', 'z', 'T (keV)', 'M500 (Msun)', 'L500 (10^44 erg/s)'))
    t.write('eFEDS_groups.fits', format='fits')
            
 
data = np.loadtxt('eFEDS.txt', skiprows = 1)
col = data[:,10]
restrict(col)
