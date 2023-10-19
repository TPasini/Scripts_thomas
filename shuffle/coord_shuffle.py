import numpy as np
from numpy import random
import csv

cat = np.loadtxt('vajtxt.txt')
z = cat[:,0]
lum= cat[:,1]
rai= cat[:,2]
deci= cat[:,3]
rabcg= cat[:,4]
decbcg= cat[:,5]


ras=np.random.shuffle(rai)

des=np.random.shuffle(deci)

with open ('shuffle.txt', 'w') as f:
    writer = csv.writer(f, delimiter=" ")
    writer.writerow(['#z', 'L_X(10^42)', 'RA', 'DEC'])
    for i, n in enumerate(z):
        writer.writerow([z[i],lum[i],rai[i],deci[i]])


