import sys, os, glob, re, argparse
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.wcs
import warnings
import pyrap.tables as pt
from astropy.cosmology import FlatLambdaCDM
from LiLF import lib_ms, lib_img, lib_util, lib_log
import astropy.io.fits as pyfits

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)

dataset = 'LBA' #OR LBA

def get_data(fname,colname):
    data=pyfits.open(fname)
    data=data[1].data
    return data[colname]

################################
##These two functions are to avoid excess printing from pyrap.tables.
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__
################################

def clean(p, MSs, res='normal', size=[1, 1], empty=False, userReg=None, apply_beam=False, do_predict=False, datacol='DATA', minuv=30, numiter=100000, fitsmask=None):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.resolution

    if res == 'normal':
        pixscale = float('%.1f' % (pixscale / 2.5))
    elif res == 'high':
        pixscale = float('%.1f' % (pixscale / 3.5))
    elif res == 'ultrahigh':
        pixscale = float('%.1f' % (pixscale / 3.5))
    elif res == 'low':
        pass  # no change

    imsize = [int(size[0] * 1.5 / (pixscale / 3600.)), int(size[1] * 1.5 / (pixscale / 3600.))]  # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

    logger.debug('Image size: ' + str(imsize) + ' - Pixel scale: ' + str(pixscale))

    if res == 'normal':
        weight = 'briggs -0.3'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.6'
        maxuv_l = None
    elif res == 'ultrahigh':
        weight = 'briggs -1'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    if empty:
        logger.info('Cleaning empty (' + str(p) + ')...')
        if dataset == 'HBA':
            imagename = f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img/empty-' + str(p)
        else:
            imagename = f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img/empty-' + str(p)
        lib_util.run_wsclean(s, 'wscleanE-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
                                 data_column='SUBTRACTED_DATA',
                                 size=imsize, scale=str(pixscale) + 'arcsec',
                                 weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0,
                                 baseline_averaging='')
    else:
        arg_dict = dict()
        # in case userReg is provided -> shallow clean, mask, merge mask with userReg, deep clean with mask
        if userReg:
            # clean 1
            logger.info('Cleaning (' + str(p) + ')...')
            if dataset == 'HBA':
                imagename = f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img/extract-' + str(p)
            else:
                imagename = f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img/empty-' + str(p)
            lib_util.run_wsclean(s, 'wscleanA-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
                                 size=imsize, scale=str(pixscale) + 'arcsec',
                                 weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                                 mgain=0.4, parallel_deconvolution=512, auto_threshold=5, join_channels='', data_column=datacol,
                                 fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3, baseline_averaging='', fits_mask=fitsmask,
                                 **arg_dict)

            # New mask method using Makemask.py
            mask = imagename + '-MFS-image.fits'
            try:
                os.system(f'MakeMask.py --RestoredIm {mask} --Th 3 --Box 150,5')

            except:
                logger.warning('Fail to create mask for %s.' % imagename + '-MFS-image.fits')
                return
            lib_img.blank_image_reg(mask + '.mask.fits', userReg, inverse=False, blankval=1.)


        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        if userReg:
            logger.info('Cleaning w/ mask (' + str(p) + ')...')
        else:
            logger.info('Cleaning (' + str(p) + ')...')
        if dataset == 'HBA':
            imagenameM = f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img/extractM-' + str(p)
        else:
            imagenameM = f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img/extractM-' + str(p)

        if apply_beam:
            arg_dict['use_idg'] = ''
            arg_dict['idg_mode'] = 'cpu'
            arg_dict['grid_with_beam'] = ''
            arg_dict['beam_aterm_update'] = 800
        else:
            arg_dict['baseline_averaging'] = ''
            if userReg:
                arg_dict['reuse_psf'] = imagename
                arg_dict['reuse_dirty'] = imagename
                arg_dict['fits_mask'] = mask + '.mask.fits'

        lib_util.run_wsclean(s, 'wscleanB-' + str(p) + '.log', MSs.getStrWsclean(), name=imagenameM, do_predict=True,
                             size=imsize, scale=str(pixscale) + 'arcsec', weight=weight, niter=numiter, local_rms='',
                             no_update_model_required='', minuv_l=minuv, maxuv_l=maxuv_l, mgain=0.4, multiscale='',
                             parallel_deconvolution=512, auto_threshold=0.5, auto_mask=3.0, save_source_list='',
                             join_channels='', fit_spectral_pol=3, channels_out=6, data_column=datacol, fits_mask=fitsmask,
                             **arg_dict)  # , deconvolution_channels=3)
        os.system('cat /beegfs/bax8338/analysis/logs/wscleanB-' + str(p) + '.log | grep "background noise"')

    if do_predict:
        if dataset == 'HBA':
            imagename= f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img/extract-forpredict'
        else:
            imagename = f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img/extract-forpredict'

        lib_util.run_wsclean(s, 'wscleanS-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename, do_predict=True,
                             size=imsize, scale=str(pixscale) + 'arcsec', weight=weight, niter=100000, local_rms='',
                             no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.4, multiscale='',
                             parallel_deconvolution=512, auto_threshold=0.5, auto_mask=3.0, save_source_list='',
                             join_channels='', fit_spectral_pol=3, channels_out=6, data_column=datacol, **arg_dict)


if dataset == 'HBA':
    logger_obj = lib_log.Logger('/beegfs/bax8338/analysis/sourcesubHBA.logger', '/beegfs/bax8338/analysis/logs')
    logger = lib_log.logger
    s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry=False)
    w = lib_util.Walker('/beegfs/bax8338/analysis/sourcesubHBA.walker')
else:
    logger_obj = lib_log.Logger('/beegfs/bax8338/analysis/sourcesubLBA.logger', '/beegfs/bax8338/analysis/logs')
    logger = lib_log.logger
    s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry=False)
    w = lib_util.Walker('/beegfs/bax8338/analysis/sourcesubLBA.walker')

parser = argparse.ArgumentParser(description='Subtraction of compact sources from ms files.')
parser.add_argument('-l', '--list', dest='list', action='store', default='', type=str)

args = parser.parse_args()
cluster_list = args.list

cl_name = get_data(cluster_list,'Name')
# cl_ra = get_data(cluster_list,'RA')
# cl_dec = get_data(cluster_list,'DEC')
cl_z = get_data(cluster_list, 'z')

for n, cluster in enumerate(cl_name):

    if dataset == 'HBA':
        MSs_extract = lib_ms.AllMSs(glob.glob(f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/*.calibrated'), s)

        os.system(f'cd /beegfs/bax8338/analysis/{str(cluster)}/HBA')

        redo = False
        if os.path.exists(f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img') and redo == True:

            os.system(f'rm -rf /beegfs/bax8338/analysis/{str(cluster)}/HBA/img')
            os.makedirs(f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img')

        os.chdir(f'/beegfs/bax8338/analysis/{str(cluster)}/HBA/img')

    else:
        MSs_extract = lib_ms.AllMSs(glob.glob(f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/mss-extract/shiftavg/*.MS-extract'), s)

        os.system(f'cd /beegfs/bax8338/analysis/{str(cluster)}/LBA')

        if os.path.exists(f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img'):
            os.system(f'rm -rf /beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img')

        os.makedirs(f'/beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img')

    sourceLLS=0.25 #Mpc. TODO Should this be tuned for each cluster? Not sure how..
    oneradinmpc = cosmo.angular_diameter_distance(cl_z[n]) / (360. / (2. * np.pi))
    scalebarlengthdeg = sourceLLS / oneradinmpc.value
    minuv_forsub = 1./(scalebarlengthdeg*np.pi/180.)

    with w.if_todo(f'find_compact_sources {str(cluster)}'):
        clean('sub-highres', MSs_extract, res='ultrahigh', minuv = minuv_forsub)

    with w.if_todo(f'produce_mask {str(cluster)}'):

        makemask='MakeMask.py'
        logger.info('Subtracting compact sources...')
        highimagename  = 'extractM-sub-highres-MFS-image.fits'
        if dataset == 'HBA':
            os.system(f'MakeMask.py --RestoredIm /beegfs/bax8338/analysis/{str(cluster)}/HBA/img/{highimagename} --Th 3')
        else:
            os.system(f'MakeMask.py --RestoredIm /beegfs/bax8338/analysis/{str(cluster)}/LBA/sub/img/{highimagename} --Th 3')

        fits_mask = highimagename + '.mask.fits'
        clean('compactmask', MSs_extract, do_predict=True, minuv = minuv_forsub, res='ultrahigh')

    with w.if_todo(f'source_subtraction {str(cluster)}'):
        logger.info('Adding DIFFUSE_SUB column to datasets...')

        mslist = MSs_extract.mssListStr

        blockPrint()

        outcolumn='DIFFUSE_SUB'
        for ms in mslist:
            ts  = pt.table(ms, readonly=False)
            colnames = ts.colnames()
            if outcolumn not in colnames:
                desc = ts.getcoldesc('DATA')
                desc['name']=outcolumn
                ts.addcols(desc)
                ts.close()
            else:
                ts.close()

        for ms in mslist:
            ts  = pt.table(ms, readonly=False)
            colnames = ts.colnames()
            if 'CORRECTED_DATA' in colnames:
                data = ts.getcol('CORRECTED_DATA')
            else:
                data = ts.getcol('DATA')
            model = ts.getcol('MODEL_DATA')
            ts.putcol(outcolumn,data-model)
            ts.close()

        enablePrint()

        logger.info('Final imaging with compact sources subtracted...')
        clean('sourcesubtracted', MSs_extract, apply_beam=False, datacol='DIFFUSE_SUB', res='low', fitsmask=highimagename + '.mask.fits')