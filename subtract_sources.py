#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import os, sys
from astroquery.simbad import Simbad
import numpy as np
import pyrap.tables as pt
import os.path
import pyregion
import argparse
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def cleanup():
    os.system('rm -rf *000*.fits') # remove channel maps
    os.system('rm -rf *_subROBUST*.fits') # remove non-masked images
    os.system('rm -rf *_ROBUST*fits') # remove non-masked images
    os.system('rm -rf *-dirty.fits') # remove all dirty images
    return


def getimsize(boxfile, cellsize=1.5):
   """
   find imsize need to image a DS9 boxfile region
   """
   r = pyregion.open(boxfile)
   
   xs = np.ceil((r[0].coord_list[2])*1.6*3600./cellsize)
   ys = np.ceil((r[0].coord_list[3])*1.6*3600./cellsize)

   imsize = np.ceil(xs) # // Round up decimals to an integer
   if(imsize % 2 == 1): 
       imsize = imsize + 1
   return int(imsize)

def compute_uvmin(redshift,sourceLLS=1.0):
    '''
    sourceLLS in units of Mpc# 
    taper output for WSClean in arcsec
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
    scalebarlengthdeg    = sourceLLS/oneradinmpc.value
    
    return 1./(scalebarlengthdeg*np.pi/180.)

def compute_taper(redshift,taperscale):
    '''
    taperscale in units of kpc# 
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
    taper    = 1e-3*taperscale/(oneradinmpc.value)
    
    return taper*3600

def adjustniter_for_taper(taper, niter):
  if taper < 5:
    return int(niter)
  if taper >= 5 and taper < 15:
    return int(niter/2)
  if taper >= 15:
    return int(niter/4)


def makeimage(mslist, imageout, pixsize, imsize, channelsout=6, niter=15000, robust=-0.5, minuv=80, uvtaper=None, multiscale=False, predict=True,fitsmask=None, deepmultiscale=False, cluster_redshift=None, column=None, log='wscleanlog', concat_mss=False):

    wsclean = 'wsclean'

    if predict == False and freq == 'LBA' and concat_mss == True:
        from LiLF import lib_ms
        from itertools import groupby

        keyfunct = lambda x: ' '.join(sorted(lib_ms.MS(x).getAntennas()))
        MSs_list = sorted(mslist, key=keyfunct)  # needs to be sorted
        groups = []
        for k, g in groupby(MSs_list, keyfunct):
            groups.append(list(g))
        print(f"Found {len(groups)} groups of datasets with same antennas.")
        for i, group in enumerate(groups, start=1):
            antennas = ', '.join(lib_ms.MS(group[0]).getAntennas())
            print(f"WSClean MS group {i}: {group}")
            print(f"List of antennas: {antennas}")

        MSs_files_clean = []
        for g, group in enumerate(groups):
            os.system(f'taql select from {group} giving wsclean_concat_{g}.MS as plain')
            MSs_files_clean.append(f'wsclean_concat_{g}.MS')
        mslist = MSs_files_clean#' '.join(MSs_files_clean)

    if uvtaper != None:
       if uvtaper < 0:
           print('Not supported uvtaper', uvtaper)
       else:
           #print(imsize, pixsize, uvtaper)
           imsizein  = int(float(imsize)*(pixsize/(uvtaper/5.)))
           pixsizein = int(uvtaper/5.)
           if float(pixsizein) < pixsize: # to deal with rounding issues which cause a 0arcsec pixelsize
              pixsizein = pixsize
              imsizein  = imsize
              
    else:
       imsizein  = imsize 
       pixsizein = pixsize

    if int(imsizein) < 511: # otherwise images get too small for multiscales
        imsizein = 512
    
    baselineav = 2.5e3*60000.*2.*np.pi *float(pixsizein)/(24.*60.*60*float(imsizein))

    # limit baseline averaging to 10, fixes prob
    if baselineav > 10.0:
        baselineav = 10.
    
    baselineav = str (baselineav)
 
    # few simple checks to make sure we have useful data
    msliststring = ' '.join(map(str, mslist))
    os.system('rm -f ' + imageout + '-*.fits')
    imcol = 'CORRECTED_DATA'
    t = pt.table(mslist[0],readonly=True) # just test for first ms in mslist
    colnames =t.colnames()
    if 'CORRECTED_DATA' not in colnames: # check which column to image
      imcol = 'DATA' 
    t.close()
    
    if column != None:
      imcol = column


    # build wsclean command      
    cmd = wsclean + ' '
    cmd += '-no-update-model-required -minuv-l ' + str(minuv) + ' '
    cmd += '-size ' + str(imsizein) + ' ' + str(imsizein) + ' -reorder '
    cmd += '-weight briggs ' + str(robust) + ' -weighting-rank-filter 3 -clean-border 1 -use-wgridder '
    cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol +' -join-channels -channels-out '
    cmd += str(channelsout) + ' -padding 1.4 '
    #cmd += '-parallel-deconvolution ' + str(int(imsizein)/2) + ' '
    if multiscale:
       if predict:
         cmd += '-multiscale '+' -multiscale-scales 0,2,4,8,16 '  
       else:
         cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '

    if fitsmask != None:
      if os.path.isfile(fitsmask): 
        cmd += '-fits-mask '+ fitsmask + ' '
      else:
        print('fitsmask: ', fitsmask, 'does not exist')
        sys.exit()
    else:
        cmd += '-auto-mask 2.5 -auto-threshold 1.0 '
    
    if uvtaper != None:
       cmd += '-taper-gaussian ' +  str(uvtaper) + 'arcsec '
    
    cmd += '-fit-spectral-pol 3 '# -beam-shape 6arcsec 6arcsec 0deg '   
    cmd += '-pol i '
    cmd += '-baseline-averaging ' + baselineav + ' '
      
    
    cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' 

    print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
    #logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
    #logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
    os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring + f'>> {log}.log')

    if deepmultiscale:
        
      # predict first to fill MODEL_DATA so we can continue with clean
      cmdp = wsclean + ' -size ' 
      cmdp += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict ' 

      
      cmdp += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
      print('PREDICT STEP for continue: ', cmdp)
      os.system(cmdp)
       
      # NOW continue cleaning  
      cmd += '-niter ' + str(niter/15) + ' -multiscale -continue ' + msliststring
      print('WSCLEAN continue: ', cmd)
      os.system(cmd)


    if predict:
      cmd = wsclean + ' -size ' 
      cmd += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict ' 

      cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
      print('PREDICT STEP: ', cmd)
      os.system(cmd)


def subtractcompact(mslist, imageout, pixsize, imsize, minuv, channelsout=6, niter=15000, robust=-0.5, outcolumn='DIFFUSE_SUB'):

    makemask = 'MakeMask.py'
     
    makeimage(mslist, imageout +'_compact', pixsize, imsize, channelsout=channelsout, niter=niter, robust=robust, minuv=minuv, predict=False)

   # make a mask
    imagename  = imageout +'_compact' + '-MFS-image.fits'
    cmdm  = makemask + ' --Th=2.0 --RestoredIm=' + imagename
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

   # re-image with mask
    makeimage(mslist, imageout +'_compactmask', pixsize, imsize, channelsout=channelsout, niter=niter, robust=-0.5, minuv=minuv, multiscale=True, predict=True, fitsmask=fitsmask, deepmultiscale=False)
    
   # now subtract the columns 
    
    for ms in mslist:
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if outcolumn not in colnames:
            desc = ts.getcoldesc('DATA')
            desc['name']=outcolumn
            ts.addcols(desc)
            ts.close() # to write results
        else:
            print(outcolumn, ' already exists')
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
    return


makemask = 'MakeMask.py'


parser = argparse.ArgumentParser(description='Make images from extraction run.')
parser.add_argument('-b','--boxfile', help='optional boxfile to set imsize automatically', type=str)
#parser.add_argument('--fitsmask', help='fitsmask for deconvolution, if not provided use automasking', type=str)
parser.add_argument('--imsize', help='image size', type=int)
parser.add_argument('-instr','--instr', help='Telescope that performed the observation. Just useful to assign the right name to images. Default to LBA.', default='LBA', type=str)
parser.add_argument('-n', '--niter', help='niter, default=25000', default=25000, type=int)
parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
parser.add_argument('--minuv', help='inner uv-cut for image in lambda, default=80', default=80., type=float)
parser.add_argument('--nodosub', help='do not subtract compact sources', action='store_true')
parser.add_argument('--onlydosub', help='only produce images subtracting compact sources', action='store_true')
parser.add_argument('--pixelscale', help='pixels size in arcsec, deafult=1.5', default=1.5, type=float)
parser.add_argument('--sourceLLS', help='size in Mpc of diffuse emission for uvcut, default=0.4', default=0.4, type=float)
parser.add_argument('--z', help='redshift of cluster, not required if --nodosub is used', default=-1.0, type=float)
parser.add_argument('-i','--imagename', help='imagename', default='image', required=True, type=str)
parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=3.0', default=3.0, type=int)
parser.add_argument('ms', nargs='*', help='msfile(s)')

args = vars(parser.parse_args())

minuv    = args['minuv']  
pixsize  = args['pixelscale']
niter    = args['niter']  
mslist   = sorted(args['ms'])
imageout = args['imagename']
freq = str(args['instr'])
#imsize   = args['imsize']

if args['boxfile'] == None and args['imsize'] == None:
  print('Incomplete input detected, either boxfile or imsize is required')
  sys.exit()
if args['boxfile'] != None and args['imsize'] != None:
  print('Wrong input detected, both boxfile and imsize are set')
  sys.exit()

if args['boxfile'] != None:
  imsize   = str(getimsize(args['boxfile'], args['pixelscale']))
if args['imsize'] != None:
  imsize = str(args['imsize']) 


#if args['imsize'] == None:
#    print 'Error: imsize not provided'
#    sys.exit()


if args['z'] < 0: # if no redshift provided try to find it automatically
  customSimbad = Simbad()
  customSimbad.add_votable_fields('z_value')
  #obj=customSimbad.query_object( args['imagename'] )
  #print obj
  #sys.exit()
  try:
    obj=customSimbad.query_object( args['imagename'] )
  #except: # in this case the cluster name is recognized by Simbad as no error was raised above

    if obj['Z_VALUE'][0] > 0.0:
        print('Found redshift',  obj['Z_VALUE'][0])
        args['z'] = obj['Z_VALUE'][0]
        os.system('echo ' + str(args['z']) + ' > redshift-used.log')
    else:
        print('Warning: Cluster known but redshift not found, assuming z=0.2')
        args['z'] = 0.2
        os.system('echo ' + str(args['z']) + ' > redshift-used.log')
  except:
    print('Cluster name is not known by Simbad' )

if not args['nodosub']:
  
    if args['z'] < 0:
        print('You need to provide a redshift, none was given')
        sys.exit()

    minuv_forsub = compute_uvmin(args['z'], sourceLLS=args['sourceLLS'])

    subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'],\
                  niter=int(niter/1.25), robust=-0.5, outcolumn='DIFFUSE_SUB')

    #  -----------------------------------------------------------------
    #  --- make the taper 50 kpc image, compact source subtracted ----
    #  -----------------------------------------------------------------

    makeimage(mslist, imageout +f'_subROBUST-0.5TAPER50kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, predict=False,
              column='DIFFUSE_SUB', uvtaper=compute_taper(args['z'],50.))

    # make a mask
    imagename  = imageout +f'_subROBUST-0.5TAPER50kpc_{freq}' + '-MFS-image.fits'
    cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print('Producing source-subtracted image tapered at 50 kpc...')
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout +f'_masksubROBUST-0.5TAPER50kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv,
              multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column='DIFFUSE_SUB',
              uvtaper=compute_taper(args['z'],50.), log='50kpc')

    #  -----------------------------------------------------------------
    #  --- make the taper 100 kpc image, compact source subtracted ----
    #  -----------------------------------------------------------------

    makeimage(mslist, imageout + f'_subROBUST-0.5TAPER100kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'], 100.), niter), robust=-0.5, minuv=minuv,
              predict=False, column='DIFFUSE_SUB', uvtaper=compute_taper(args['z'], 100.))

    # make a mask
    imagename = imageout + f'_subROBUST-0.5TAPER100kpc_{freq}' + '-MFS-image.fits'
    cmdm = makemask + ' --Th=' + str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print('Producing source-subtracted image tapered at 100 kpc...')
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout + f'_masksubROBUST-0.5TAPER100kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'], 100.), niter), robust=-0.5, minuv=minuv,
              multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column='DIFFUSE_SUB',
              uvtaper=compute_taper(args['z'], 100.), log='100kpc')

    #  -----------------------------------------------------------------
    #  --- make the taper 150 kpc image, compact source subtracted ----
    #  -----------------------------------------------------------------

    makeimage(mslist, imageout + f'_subROBUST-0.5TAPER150kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'], 150.), niter), robust=-0.5, minuv=minuv,
              predict=False, column='DIFFUSE_SUB', uvtaper=compute_taper(args['z'], 150.))

    # make a mask
    imagename = imageout + f'_subROBUST-0.5TAPER150kpc_{freq}' + '-MFS-image.fits'
    cmdm = makemask + ' --Th=' + str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print('Producing source-subtracted image tapered at 150 kpc...')
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout + f'_masksubROBUST-0.5TAPER150kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'],
              niter=adjustniter_for_taper(compute_taper(args['z'], 150.), niter), robust=-0.5, minuv=minuv,
              multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column='DIFFUSE_SUB',
              uvtaper=compute_taper(args['z'], 150.), log='150kpc')


    #  -----------------------------------------------------------------
    #  --- make the taper 90 image, compact source subtracted ----
    #  -----------------------------------------------------------------

    makeimage(mslist, imageout + f'_subROBUST-0.5TAPER90_{freq}', 18, 2000, channelsout=args['channelsout'],
              niter=1000000, robust=-0.5, minuv=minuv,
              predict=False, column='DIFFUSE_SUB', uvtaper=90)

    # make a mask
    imagename = imageout + f'_subROBUST-0.5TAPER90_{freq}' + '-MFS-image.fits'
    cmdm = makemask + ' --Th=' + str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print('Producing source-subtracted image tapered at 90 arcsec...')
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout + f'_masksubROBUST-0.5TAPER90_{freq}', 18, 2000, channelsout=args['channelsout'],
              niter=1000000, robust=-0.5, minuv=minuv, predict=False, fitsmask=fitsmask, deepmultiscale=False,
              column='DIFFUSE_SUB', uvtaper=90, log='90arcsec')



if (not args['onlydosub'] == True) and (args['z'] >0.): # otherwise cannot computer compute_taper(z)


    #  --------------------------------------------------
    #  --- make the taper 50 kpc image ----
    #  --------------------------------------------------
    makeimage(mslist, imageout + f'_ROBUST-0.5TAPER50kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, predict=False, uvtaper=compute_taper(args['z'],50.))

    # make a mask
    imagename  = imageout +  f'_ROBUST-0.5TAPER50kpc_{freq}' +'-MFS-image.fits'
    cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print(cmdm)
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout +f'_maskROBUST-0.5TAPER50kpc_{freq}', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False,uvtaper=compute_taper(args['z'],50.))


if os.path_exists('wsclean_concat_*.MS'):
    os.remove('wsclean_concat_*')
cleanup()
