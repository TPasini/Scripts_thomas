import os, sys
import numpy as np

if len(sys.argv) == 1:
    print('Spindex.py flux1 flux2 freq1 freq2 [err1] [err2]')
    sys.exit(0)

flux1 = float(sys.argv[1])
flux2 = float(sys.argv[2])
freq1 = float(sys.argv[3])
freq2 = float(sys.argv[4])
try: err1 = float(sys.argv[5])
except: err1 = 0
try: err2 = float(sys.argv[6])
except: err2 = 0

alpha = (np.log10(flux1/flux2))/(np.log10(freq1/freq2))
erralpha=(1/(np.log(freq2)-np.log(freq1)))*(np.sqrt(((err1/flux1)**2)+((err2/flux2)**2)))

print('alpha=',alpha)
print('error=',erralpha)
