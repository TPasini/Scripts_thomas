#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import astropy.units as u
import pylab as pl
import matplotlib.pyplot as plt
import aplpy
from matplotlib import rc
from matplotlib.colors import ListedColormap
from matplotlib import colors
from matplotlib import cm
rc('text', usetex=False)
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from astropy.cosmology import LambdaCDM
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import os
import glob
import linmix
from PIL import Image

#my_cosmo = LambdaCDM(H0=70, Om0=0.3)
my_cosmo = LambdaCDM(73.0, 0.27, 0.73)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### This first step needs to occur outside Python

# First, use flux_extract.py to measure the X-ray flux density in the boxes
#python flux_extract.py --region Regions/zwcl1451_minihalo_boxes_20asec_core_excised.reg  \
#    --output FluxExtract/zwcl1451_minihalo_surfacebrightness_xray_boxes_20asec_coreexcised.txt \
#    FITS_XRay/ZwCl1455.0+2232_0.5-2.0_flux_3arcsecsmooth.fits
# Second, use flux_extract.py to measure the X-ray flux density in the boxes
#python flux_extract.py --region Regions/zwcl1451_minihalo_boxes_20asec_core_excised.reg \
#    --output FluxExtract/zwcl1451_minihalo_surfacebrightness_radio_boxes_20asec_coreexcised.txt \
#    FITS_Radio/zwcl1454_lotss-subtracted_wsclean_robust-05_multiscale_i_wgridder_clean_5sigmask_uvtaper_15asec-MFS-image.2ax.img.regrid.fits

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


# Define a quick-use linear correlation plot function, the classic "y=mx+b" format
def radio_xray(xray, m, b):
    return m * xray + b

# Define the beam solid angle
def beam_solid_angle(maj, min):
    return (np.pi / (4 * np.log(2))) * maj * min

# Get info from the radio FITS header to derive the beam solid angle (so you can compare units in per-arcsec^2)
radio_hdr   = fits.getheader('briggs_-0.25_taper20_masked_SUB2_restbeam35_PSZ2G226.18+76.79_image_9_-MFS-image.fits')
majax       = (radio_hdr['BMAJ']*u.degree).to(u.arcsec)
minax       = (radio_hdr['BMIN']*u.degree).to(u.arcsec)

beam_solang = beam_solid_angle( majax, minax )



# Load the text files with the measurements in
radio_meas = np.loadtxt('FluxExtract/A1413_SB_radio_noAGN_XrayArea_noMH.txt')
xray_meas  = np.loadtxt('FluxExtract/A1413_SB_Xray_noAGN_XrayArea_noMH.txt')


# Define your I_r, dI_r, I_x, dI_x...
region_no       = radio_meas[:,0]
radio_flux      = radio_meas[:,3] * 1.149 / beam_solang.value  # Radio flux density in Jy per arcsec^2
radio_flux_err  = radio_meas[:,4] / beam_solang.value          # Radio flux density in Jy per arcsec^2
log_radio_err   = 0.434*radio_flux_err/radio_flux

xray_flux       = xray_meas[:,3]
xray_flux_err   = xray_meas[:,4]
log_xray_err    = 0.434*xray_flux_err/xray_flux

# Run linmix to do the MCMC
lm = linmix.LinMix( np.log10(xray_flux), np.log10(radio_flux), log_xray_err, log_radio_err, K=2, seed=42)
lm.run_mcmc(silent=True)

print("** FIT RESULTS: SLOPE {0} +{1}/-{2} ".format( '%.03f'%(np.nanpercentile(lm.chain['beta'], [50])[0]),\
    '%.03f'%(np.nanpercentile(lm.chain['beta'], [84])[0] - np.nanpercentile(lm.chain['beta'], [50])[0]),\
    '%.03f'%(np.nanpercentile(lm.chain['beta'], [50])[0] - np.nanpercentile(lm.chain['beta'], [16])[0]) ))


# Use the linmix chain to populate a full set of models for plotting
log_xray_sb_model  = np.arange(-10, -5, 0.25)
log_radio_sb_model = []
# This is not efficiently coded, but I couldn't remember how to code it up more efficiently. It's fast enough for now...
for i in range(len(lm.chain)):
    log_radio_sb_model.append(radio_xray(log_xray_sb_model, lm.chain[i]['beta'], lm.chain[i]['alpha'] ))

log_radio_sb_model = np.array(log_radio_sb_model)




# Generate the I_x / I_r plot:
fig, (s1) = pl.subplots(nrows=1, ncols=1, figsize=(8,12))
# Plot your data
s1.errorbar( np.log10(xray_flux), pl.log10(radio_flux), xerr=log_xray_err, yerr=log_radio_err,\
    mec='k', mfc='dodgerblue', ecolor='k', marker='o', markersize=5, ls='', capsize=2, elinewidth=0.75,\
    zorder=1 )

# Define the "fill range" according to user-defined percentiles of the radio SB model.
# Percentages defined in the square brackets -- using 16 and 84 so this plots the 1-sigma uncertainties:
fill_range = np.nanpercentile(log_radio_sb_model, [16, 84], axis=0)
s1.fill_between(log_xray_sb_model, fill_range[0], fill_range[1], color='green', alpha=0.25, zorder=2)

# Now plot the 50th percentile to show the median fit:
log_radio_sb_median = np.nanpercentile(log_radio_sb_model, [50], axis=0)
s1.plot(log_xray_sb_model, log_radio_sb_median[0], ls='--', color='green', linewidth=1, zorder=2)

# The rest of this is largely formatting... Edit away to your heart's content
s1.set_xlabel(r'log$_{10} ~ I_{\rm{x}}$ [counts s$^{-1}$ arcsec$^{-2}$]',fontsize=13)
s1.set_ylabel(r'log$_{10} ~ I_{\rm{radio}}$ [Jy arcsec$^{-2}$]',fontsize=13)
s1.minorticks_on()
#s1.set_xticks(size=11)
#s1.set_yticks(size=11)
s1.tick_params(direction='in', color='k', which='both', top=True, bottom=True, left=True, right=True)
#s1.set_xlim(-8.5, -6.0)
#s1.set_ylim(-6.5, -3.75)

### Change these values to suit your frequency and fit results (hence the print statement earlier)
s1.text(-8.4, -4.1, r'145 MHz', color='k', weight='bold', fontsize=16,\
    horizontalalignment='left', verticalalignment='bottom' )
s1.text(-8.4, -4.3, r'$I_{\rm{R}} \propto I_{\rm{X}}^{0.782^{+0.17}_{-0.15}}$', color='k', weight='bold', fontsize=14,\
    horizontalalignment='left', verticalalignment='bottom' )

plt.savefig('MiniHalo_Xray_Radio_correlation.png')


