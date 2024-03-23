#!/usr/bin/env python
# coding: utf-8

import os, argparse, sys
import numpy as np
import glob
from LiLF import lib_util, lib_log, lib_ms
import astropy.io.fits as pyfits
from astropy.cosmology import FlatLambdaCDM
from astropy.wcs import WCS
from termcolor import colored
from lib_fits import flatten

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
logger_obj = lib_log.Logger('radio_workflow')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('radio_workflow.walker')

def get_data(fname, colname):
    data = pyfits.open(fname)
    data = data[1].data
    return data[colname]

def get_imgparameters(myimage):
    med = np.median(myimage)
    mean = np.mean(myimage)
    if mean / med > 1:
        min = -np.nanstd(myimage)
        max = 5 * np.nanstd(myimage)
    elif mean / med < -5:
        min = -np.nanstd(myimage)
        max = 2 * np.nanstd(myimage)
    else:
        min = -5 * np.nanstd(myimage)
        max = 20 * np.nanstd(myimage)
    return min, max

def get_center_coordinates(fits_filename):

    head, data = flatten(fits_filename)
    header = head

    wcs = WCS(header)

    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']

    center_pixel_x = naxis1 / 2.0
    center_pixel_y = naxis2 / 2.0

    center_world = wcs.all_pix2world(center_pixel_x, center_pixel_y, 0)
    ra, dec = center_world[0], center_world[1]

    return ra, dec

def calc_noise(data, niter=100, eps=1e-6):
    """
    Return the rms
    """

    rms = 1.;
    oldrms = 1.
    for i in range(niter):
        rms = np.nanstd(
            data)  # IMPORTANT: USE NP.MEAN IF YOU HAVE THE RMS IMAGE, NP.NANSTD IF YOU HAVE THE NORMAL IMAGE
        if np.abs(oldrms - rms) / rms < eps:
            return rms

        oldrms = rms

    raise Exception('Noise estimation failed to converge.')


parser = argparse.ArgumentParser(description='Automatically produce .fits file, spectral index maps and pdf of radio datasets.')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Working directory. Must contain subdirectories named "LBA, "HBA", "GMRT" (not necessarily all of them)')
parser.add_argument('--name', dest='name', type=str, default='target', help='Name of the target.')
parser.add_argument('--z', dest='redshift', type=float, help='Redshift of the target.')
parser.add_argument('--dosourcesub', dest='dosourcesub', help='Perform subtraction of compact sources and create DIFFUSE_SUB column. Default to False.', action='store_true')

args = parser.parse_args()
pathdir = args.path
targetname = args.name
redshift = args.redshift
redshift = float(redshift)

do_LBA = False
do_HBA = False

LBA_dir = f'{pathdir}/LBA'
if os.path.exists(LBA_dir):
    do_LBA = True
    print(colored(f'Found LBA data.', 'green'))
    if not os.path.exists(f'{pathdir}/LBA/images'):
        os.makedirs(f'{pathdir}/LBA/images')
    if not os.path.exists(f'{pathdir}/LBA/datasets'):
        os.makedirs(f'{pathdir}/LBA/datasets')
        os.system(f'mv {pathdir}/LBA/*MS* {pathdir}/LBA/datasets')
    MSs_LBA = lib_ms.AllMSs(glob.glob(f'{pathdir}/LBA/datasets/*MS*'), s)
    if not os.path.exists(f'{pathdir}/LBA/logs'):
        os.makedirs(f'{pathdir}/LBA/logs')
    if args.dosourcesub == True:
        with w.if_todo(f'subtraction_LBA'):
            if not os.path.exists(f'{pathdir}/LBA/subtraction'):
                os.makedirs(f'{pathdir}/LBA/subtraction')
            logger.info('Performing LBA compact-sources subtraction...')
            subtract_LBA = f'subtract_sources.py --onlydosub -i {pathdir}/LBA/subtraction/{targetname}_forsubtract --imsize 1024 --z {redshift} --instr LBA {pathdir}/LBA/datasets/*MS*'
            s.add(subtract_LBA, log='LBAsub.log', commandType="python")
            s.run()
            os.system(f'rm *.log')
            os.system('rm last_MyCasapy2BBS.obj')
    if not os.path.exists(f'{pathdir}/LBA/pdf'):
        os.makedirs(f'{pathdir}/LBA/pdf')


HBA_dir = f'{pathdir}/HBA'
if os.path.exists(HBA_dir):
    do_HBA = True
    print(colored(f'Found HBA data.', 'green'))
    if not os.path.exists(f'{pathdir}/HBA/images'):
        os.makedirs(f'{pathdir}/HBA/images')
    if not os.path.exists(f'{pathdir}/HBA/datasets'):
        os.makedirs(f'{pathdir}/HBA/datasets')
        os.system(f'mv {pathdir}/HBA/*calibrated {pathdir}/HBA/datasets')
        os.system(f'mv {pathdir}/HBA/*archive* {pathdir}/HBA/datasets')
    MSs_HBA = lib_ms.AllMSs(glob.glob(f'{pathdir}/HBA/datasets/*calibrated'),s) #+ glob.glob(f"{pathdir}/HBA/datasets/*archive*"), s)
    if not os.path.exists(f'{pathdir}/HBA/logs'):
        os.makedirs(f'{pathdir}/HBA/logs')
    if args.dosourcesub == True:
        with w.if_todo(f'subtraction_HBA'):
            if not os.path.exists(f'{pathdir}/HBA/subtraction'):
                os.makedirs(f'{pathdir}/HBA/subtraction')
            logger.info('Performing HBA compact-sources subtraction...')
            subtract_LBA = f'subtract_sources.py --onlydosub -i {pathdir}/HBA/subtraction/{targetname}_forsubtract --imsize 1024 --z {redshift} --instr HBA {pathdir}/HBA/datasets/*calibrated'
            s.add(subtract_LBA, log='HBAsub.log', commandType="python")
            s.run()
            os.system(f'rm *.log')
            os.system('rm last_MyCasapy2BBS.obj')
    if not os.path.exists(f'{pathdir}/HBA/pdf'):
        os.makedirs(f'{pathdir}/HBA/pdf')

if os.path.exists(LBA_dir) and os.path.exists(HBA_dir):
    if not os.path.exists(f'{pathdir}/spidx'):
        os.makedirs(f'{pathdir}/spidx')

if do_HBA == False and do_LBA == False and do_GMRT == False:
    print(colored(f'No HBA, LBA or GMRT directories found in {pathdir}. Stopping...', 'red'))
    sys.exit()

if redshift:
    taper_kpc = 50  # kpc
    conversion = float((cosmo.kpc_proper_per_arcmin(redshift) / 60).value)
    taper_asec = taper_kpc / conversion

    taper_kpclow = 100  # kpc
    taper_aseclow = taper_kpclow / conversion

if do_HBA == True and do_LBA == True:
    do_spidx = True
else:
    do_spidx = False

if do_HBA == True:

    with w.if_todo(f'wsclean_HBA_{targetname}'):
        logger.info(f'Imaging HBA {targetname}...')

        lib_util.run_wsclean(s, 'nominalHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_briggs_HBA', size=2000, scale='1.5arcsec',
                             no_update_model_required='', padding=1.4,
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, 'highHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_highres_HBA', size=2000, scale='1.2arcsec',
                             no_update_model_required='', padding=1.4,
                             weight='briggs -1', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper30_HBA', size=2000, scale='6arcsec',
                             no_update_model_required='', padding=1.4, taper_gaussian='30asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper90_HBA', size=2000, scale='18arcsec',
                             no_update_model_required='', padding=1.4, taper_gaussian='90asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowkpcHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper50kpc_HBA', size=2000, scale=f'{float(taper_asec / 5)}asec',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_asec)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowkpcSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper50kpcSUB_HBA', size=2000,
                             scale=f'{float(taper_asec / 5)}asec', data_column ='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_asec)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowkpcHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper100kpc_HBA', size=2000,
                             scale=f'{float(taper_aseclow / 5)}asec',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_aseclow)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper30SUB_HBA', size=2000, data_column ='DIFFUSE_SUB',
                             scale='6asec', no_update_model_required='', padding=1.4, taper_gaussian='30asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper90SUB_HBA', size=2000, data_column ='DIFFUSE_SUB',
                             scale='18asec', no_update_model_required='', padding=1.4, taper_gaussian='90asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowkpcSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper100kpcSUB_HBA', size=2000,
                             scale=f'{float(taper_aseclow / 5)}asec', data_column ='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_aseclow)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_spidx_HBA', size=2000,
                             scale='3asec', minuv_l = 80, maxuv_l = 14000, no_update_model_required='', padding=1.4,
                             weight='briggs -1', niter=20000, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_spidxSUB_HBA', size=2000, data_column ='DIFFUSE_SUB',
                             scale='6asec', no_update_model_required='', padding=1.4, taper_gaussian='30asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxlowSUBHBA.log', MSs_HBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/HBA/images/{targetname}_taper50kpc_spidxSUB_HBA', size=2000,
                             scale=f'{float(taper_aseclow/5)}asec', data_column ='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_aseclow)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        os.system(f'rm {pathdir}/HBA/images/*-00*-*.fits')
        os.system(f'rm {pathdir}/HBA/images/*dirty*')
        os.system(f'rm {pathdir}/HBA/images/*psf*')
        if not os.path.exists(f'{pathdir}/HBA/images/cleanfiles'):
            os.makedirs(f'{pathdir}/HBA/images/cleanfiles')
        os.system(f'mv {pathdir}/HBA/images/*model* {pathdir}/HBA/images/cleanfiles')
        os.system(f'mv {pathdir}/HBA/images/*residual* {pathdir}/HBA/images/cleanfiles')

if do_LBA  == True:
    with w.if_todo(f'wsclean_LBA_{targetname}'):
        logger.info(f'Imaging LBA {targetname}...')

        lib_util.run_wsclean(s, f'nominalLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_briggs_LBA', size=2000, scale='3arcsec',
                             no_update_model_required='', padding=1.4,
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'highLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_highres_LBA', size=2000, scale='2.5arcsec',
                             no_update_model_required='', padding=1.4,
                             weight='briggs -1', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper30_LBA', size=2000, scale='6arcsec',
                             no_update_model_required='', padding=1.4, taper_gaussian='30asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper90_LBA', size=2000, scale='18arcsec',
                             no_update_model_required='', padding=1.4, taper_gaussian='90asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowkpcLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper50kpc_LBA', size=2000,
                             scale=f'{float(taper_asec / 5)}asec',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_asec)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'/splowkpcLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper100kpc_LBA', size=2000,
                             scale=f'{float(taper_aseclow / 5)}asec',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_aseclow)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowkpcSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper50kpcSUB_LBA', size=2000,
                             scale=f'{float(taper_asec / 5)}asec', data_column='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_asec)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowkpcSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper100kpcSUB_LBA', size=2000,
                             scale=f'{float(taper_aseclow / 5)}asec', data_column='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_aseclow)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'lowSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper30SUB_LBA', size=2000,
                             data_column='DIFFUSE_SUB',
                             scale='6asec', no_update_model_required='', padding=1.4, taper_gaussian='30asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'splowSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper90SUB_LBA', size=2000,
                             data_column='DIFFUSE_SUB',
                             scale='18asec', no_update_model_required='', padding=1.4, taper_gaussian='90asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_spidx_LBA', size=2000,
                             scale='3asec', minuv_l=80, maxuv_l=14000, no_update_model_required='', padding=1.4,
                             weight='briggs -1', niter=20000, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_spidxSUB_LBA', size=2000,
                             data_column='DIFFUSE_SUB',
                             scale='5asec', no_update_model_required='', padding=1.4, taper_gaussian='25asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        lib_util.run_wsclean(s, f'spidxlowSUBLBA.log', MSs_LBA.getStrWsclean(), concat_mss=False,
                             name=f'{pathdir}/LBA/images/{targetname}_taper50kpc_spidxSUB_LBA', size=2000,
                             scale=f'{float(taper_asec/5)}asec', data_column='DIFFUSE_SUB',
                             no_update_model_required='', padding=1.4, taper_gaussian=f'{float(taper_asec)}asec',
                             weight='briggs -0.3', niter=20000, minuv_l=80, mgain=0.8, join_channels='',
                             fit_spectral_pol=3, channels_out=6, multiscale='', pol='i', reorder='',
                             baseline_averaging=10.0)

        os.system(f'rm {pathdir}/LBA/images/*-00*-*.fits')
        os.system(f'rm {pathdir}/LBA/images/*dirty*')
        os.system(f'rm {pathdir}/LBA/images/*psf*')
        if not os.path.exists(f'{pathdir}/LBA/images/cleanfiles'):
            os.makedirs(f'{pathdir}/LBA/images/cleanfiles')
        os.system(f'mv {pathdir}/LBA/images/*model* {pathdir}/LBA/images/cleanfiles')
        os.system(f'mv {pathdir}/LBA/images/*residual* {pathdir}/LBA/images/cleanfiles')

do_plots = True
if do_plots == True:
    with w.if_todo(f'Plots_{targetname}'):
        logger.info(f'Producing pdf files for {targetname}...')

        if do_HBA == True:
            ##HBA
            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_briggs_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            #Assuming that the centre is the same at all resolution and frequencies
            ra, dec = get_center_coordinates(f'{pathdir}/HBA/images/{targetname}_briggs_HBA-MFS-image.fits')
            cmd_plot_nom_HBA = f'makeimage.py {pathdir}/HBA/images/{targetname}_briggs_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} ' \
                               f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_briggs_HBA --interval {min} {max}'
            s.add(cmd_plot_nom_HBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_highres_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_high_HBA = f'makeimage.py {pathdir}/HBA/images/{targetname}_highres_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} ' \
                                f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_highres_HBA --interval {min} {max}'
            s.add(cmd_plot_high_HBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_taper30_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_low_HBA = f'makeimage.py {pathdir}/HBA/images/{targetname}_taper30_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                               f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_taper30_HBA --interval {min} {max}'
            s.add(cmd_plot_low_HBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_taper50kpc_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_lowkpc_HBA = f'makeimage.py {pathdir}/HBA/images/{targetname}_taper50kpc_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                  f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_taper50kpc_HBA --interval {min} {max}'
            s.add(cmd_plot_lowkpc_HBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_taper30SUB_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_low_HBA_sub = f'makeimage.py {pathdir}/HBA/images/{targetname}_taper30SUB_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                   f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_taper30SUB_HBA --interval {min} {max}'
            s.add(cmd_plot_low_HBA_sub, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/HBA/images/{targetname}_taper100kpcSUB_HBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_splow_HBA_sub = f'makeimage.py {pathdir}/HBA/images/{targetname}_taper100kpcSUB_HBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                        f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/HBA/pdf/{targetname}_taper100kpcSUB_HBA --interval {min} {max}'
            s.add(cmd_plot_splow_HBA_sub, log='plots.log', commandType="python")
            s.run()


        if do_LBA == True:
            ##LBA
            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_briggs_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            ra, dec = get_center_coordinates(f'{pathdir}/LBA/images/{targetname}_briggs_LBA-MFS-image.fits')
            cmd_plot_nom_LBA = f'makeimage.py {pathdir}/LBA/images/{targetname}_briggs_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} ' \
                               f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_briggs_LBA --interval {min} {max * 2}'
            s.add(cmd_plot_nom_LBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_highres_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_high_LBA = f'makeimage.py {pathdir}/LBA/images/{targetname}_highres_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} ' \
                                f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_highres_LBA --interval {min} {max}'
            s.add(cmd_plot_high_LBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_taper30_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_low_LBA = f'makeimage.py {pathdir}/LBA/images/{targetname}_taper30_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                               f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_taper30_LBA --interval {min} {max}'
            s.add(cmd_plot_low_LBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_taper50kpc_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_lowkpc_LBA = f'makeimage.py {pathdir}/LBA/images/{targetname}_taper50kpc_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                  f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_taper50kpc_LBA --interval {min} {max}'
            s.add(cmd_plot_lowkpc_LBA, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_taper30SUB_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_low_LBA_sub = f'makeimage.py {pathdir}/LBA/images/{targetname}_taper30SUB_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                   f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_taper30SUB_LBA --interval {min} {max}'
            s.add(cmd_plot_low_LBA_sub, log='plots.log', commandType="python")
            s.run()

            data, header = pyfits.getdata(f'{pathdir}/LBA/images/{targetname}_taper100kpcSUB_LBA-MFS-image.fits', header=True)
            min, max = get_imgparameters(data)
            cmd_plot_splow_LBA_sub = f'makeimage.py {pathdir}/LBA/images/{targetname}_taper100kpcSUB_LBA-MFS-image.fits --type radio -z {redshift} --radec {ra} {dec} -s 1000 ' \
                                        f'--show_beam --sbar 500 --show_contours --stretch power -o {pathdir}/LBA/pdf/{targetname}_taper100kpcSUB_LBA --interval {min} {max}'
            s.add(cmd_plot_splow_LBA_sub, log='plots.log', commandType="python")
            s.run()

if do_spidx == True:
    with w.if_todo(f'Spidx_{targetname}'):
        logger.info(f'Producing spidx maps for {targetname}...')
        splog = 'spidxlog.log'
        ra, dec = get_center_coordinates(f'{pathdir}/LBA/images/{targetname}_briggs_LBA-MFS-image.fits')

        cmd_spidx = f'spidxmap_wuplims.py --size 0.35 --radec {ra} {dec} --noise --sigma 2.5 --fluxscaleerr 0.1 --output {pathdir}/spidx/{targetname}_spidx.fits ' \
                    f'{pathdir}/HBA/images/{targetname}_spidx_HBA-MFS-image.fits {pathdir}/LBA/images/{targetname}_spidx_LBA-MFS-image.fits'
        s.add(cmd_spidx, log=splog, commandType="python")
        s.run()

        cmd_plot_spidx = f'plot_image.py --type si -c {ra} {dec} -s 20 20 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidx --interval 0 -3 ' \
                         f'--sbar_kpc 100 {pathdir}/spidx/{targetname}_spidx.fits'
        s.add(cmd_plot_spidx, log=splog, commandType="python")
        s.run()

        cmd_plot_spidxerr = f'plot_image.py --type sierr -c {ra} {dec} -s 20 20 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidxerr ' \
                            f'--sbar_kpc 100 {pathdir}/spidx/{targetname}_spidx-err.fits'
        s.add(cmd_plot_spidxerr, log=splog, commandType="python")
        s.run()

        cmd_spidx_sub = f'spidxmap_wuplims.py --size 0.65 --radec {ra} {dec} --noise --sigma 2.5 --fluxscaleerr 0.1 --beam 45 45 0 --output {pathdir}/spidx/{targetname}_spidx_sourcesub.fits ' \
                        f'{pathdir}/HBA/images/{targetname}_spidxSUB_HBA-MFS-image.fits {pathdir}/LBA/images/{targetname}_spidxSUB_LBA-MFS-image.fits'
        s.add(cmd_spidx_sub, log=splog, commandType="python")
        s.run()

        cmd_plot_spidx_sub = f'plot_image.py --type si -c {ra} {dec} -s 30 30 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidx_sourcesub --interval 0 -3 ' \
                             f'--sbar_kpc 100 spidx/{targetname}_spidx_sourcesub.fits'
        s.add(cmd_plot_spidx_sub, log=splog, commandType="python")
        s.run()

        cmd_plot_spidxerr_sub = f'plot_image.py --type sierr -c {ra} {dec} -s 30 30 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidxerr_sourcesub ' \
                                f'--sbar_kpc 100 spidx/{targetname}_spidx_sourcesub-err.fits'
        s.add(cmd_plot_spidxerr_sub, log=splog, commandType="python")
        s.run()

        cmd_spidx_superlow_sub = f'spidxmap_wuplims.py --size 0.65 --radec {ra} {dec} --noise --sigma 2.5 --fluxscaleerr 0.1 --beam 130 130 0 --output {pathdir}/spidx/{targetname}_spidx_taper50kpc_sourcesub.fits ' \
                                 f'{pathdir}/HBA/images/{targetname}_taper50kpc_spidxSUB_HBA-MFS-image.fits {pathdir}/LBA/images/{targetname}_taper50kpc_spidxSUB_LBA-MFS-image.fits'
        s.add(cmd_spidx_superlow_sub, log=splog, commandType="python")
        s.run()

        cmd_plot_spidx_splow_sub = f'plot_image.py --type si -c {ra} {dec} -s 40 40 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidx_taper50kpc_sourcesub --interval 0 -3 ' \
                                      f'--sbar_kpc 100 spidx/{targetname}_spidx_taper50kpc_sourcesub.fits'
        s.add(cmd_plot_spidx_splow_sub, log=splog, commandType="python")
        s.run()

        cmd_plot_spidxerr_splow_sub = f'plot_image.py --type sierr -c {ra} {dec} -s 40 40 --redshift {redshift} -o {pathdir}/spidx/{targetname}_spidxerr_taper50kpc_sourcesub ' \
                                         f'--sbar_kpc 100 spidx/{targetname}_spidx_taper50kpc_sourcesub-err.fits'
        s.add(cmd_plot_spidxerr_splow_sub, log=splog, commandType="python")
        s.run()

logger.info('Done.')

logfile = sorted(glob.glob(f'radio_workflow*.logger'))[-1]
with open(logfile, 'r') as f:
    last_line = f.readlines()[-1]
    if not "Done" in last_line:
        logger.error(f'Something went wrong - check the logfile {logfile}.')