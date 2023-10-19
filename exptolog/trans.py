import numpy as np
import csv
import os

if os.path.exists('data'):
    os.remove('data')

data=np.loadtxt('match.txt')

xlum=data[:,5]*1e+44
radiolum=data[:,23]
z=data[:,2]

for i, num in enumerate(xlum):
    with open('data', 'a') as f:
        writer = csv.writer(f, delimiter=" ")
        writer.writerow(["%2.7f" % np.log10(xlum[i]), '1', '', "%2.7f" % np.log10(radiolum[i]), '1', '', "%1.8f" % z[i], '1'])


