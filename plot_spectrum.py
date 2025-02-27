import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from matplotlib.ticker import ScalarFormatter


parser = argparse.ArgumentParser(description='Plot synchrotron spectrum given flux density and frequency')
parser.add_argument('-flux', nargs='+', type=float, help='Flux density values in Jy.')
parser.add_argument('-freq', nargs='+', type=float, help='Frequency in same order of flux density.')
parser.add_argument('--fluxerr', nargs='+', type=float, help='Flux density errors in Jy.')

args = parser.parse_args()
flux = args.flux
freq = args.freq
fluxerr = args.fluxerr

alpha=[]
alphaerr=[]
for i, n in enumerate(freq):
    if i < len(freq)-1:
        spindex = (np.log10(flux[i]/flux[i+1]))/(np.log10(freq[i]/freq[i+1]))
        alpha.append(spindex)
        if args.fluxerr:
            errspindex = (1 / (np.log10(freq[i + 1]) - np.log10(freq[i]))) * (np.sqrt(((fluxerr[i] / flux[i]) ** 2) + ((fluxerr[i + 1] / flux[i]) ** 2)))
            alphaerr.append(errspindex)
            plt.text((freq[i+1]+freq[i])/2, (flux[i+1]+flux[i])/2.5, f'$\\alpha = {spindex:.2f} \pm {errspindex:.2f}$', fontsize=12)
        else:
            plt.text((freq[i+1]+freq[i])/2, (flux[i+1]+flux[i])/2.5, f'$\\alpha = {spindex:.2f}$', fontsize=12)

plt.errorbar(freq, flux, yerr=fluxerr, elinewidth=0, capsize=0,capthick=0, marker='o', markersize=13)
if args.fluxerr:
    plt.fill_between(freq, np.subtract(flux, fluxerr), np.add(flux, fluxerr), alpha=0.3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (Hz)', fontsize=15)
plt.ylabel('Flux density (Jy)', fontsize=15)
plt.tick_params(labelsize=14)
plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))

plt.show()

