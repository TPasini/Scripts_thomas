import numpy
import os, sys
import glob
import subprocess

# fix for rficonsole LD_LIBRARY_PATH (cannot handle libgfortran provided by casapy)
# for the AOFlagger in the Reinout environment
# my_env = os.environ.copy()
# my_env["LD_LIBRARY_PATH"]="/lib:/usr/lib:/home/rvweeren/software/lofar/opt/LofIm/lib64:/home/rvweeren/software/root/lib:/home/rvweeren/software/heasoft-6.16/x86_64-unknown-linux-gnu-libc2.19-0/lib"
# GDG: no my_env in the subprocess for AOFlagger

# parang
# setjy modimage
# interp nearest
# timeragne
# check that modimage is right band (twice in scripts)
# no tfcrop flaggin, problem for calibrator scans with more than 50% flagged due to slews for example

# Nov 2015
# fixed tfcrop problem by unflag and flag action with tb.open
# linearflag in applycal last step!!
# changed to Taylor Butler 2013
# fixed aoflagger call with LD_LIBRARY_PATH and new rfis

basevis = '/lofar2/bba1916/D-array/17B-298.sb35113940.eb35590563.58363.5770291088.ms'
spwlist = [i for i in range(0,16)]

##test before on spw=6,12 -- good spw to look for bad antenna
#spwlist = ['6']

targetname = "A3411L"
phasecalname = "J0902-1415"

allspw_dataset   = 'Darray_allspw.ms'
allspw_target    = 'Darray.ms'

ref_ant = 'ea05'
num_ant = 26
num_plots = (num_ant/4)

pi = numpy.pi


# FLAGGING PERCENTAGE IN THE END OF THE CALIBRATION
def flag_stat(ms):
  default('flagdata')
  t = flagdata(vis=ms, mode='summary', field='', scan='', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False)
  log = 'Flag statistics:'
  log += '\nAntenna, '
  for k in sorted(t['antenna']):
      log += k +': %d.2%% - ' % (100.*t['antenna'][k]['flagged']/t['antenna'][k]['total'])
  log += '\nCorrelation, '
  for k, v in t['correlation'].items():
      log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
  log += '\nSpw, '
  for k, v in t['spw'].items():
      log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
  log += '\nTotal: %d.2%%' % (100.*t['flagged']/t['total'])
  
  print log.replace(' - \n','\n')



for spw_id in spwlist:
  path_plot='./plot'  
  if not os.path.exists(path_plot+str(spw_id)+'spw'):
    os.makedirs(path_plot+str(spw_id)+'spw')
  
  path_plot_spw=path_plot+str(spw_id)+'spw/'

  
  #uvlimit = "<300klambda"
  #uvlimitp= "<50klambda"
  uvlimit = "<1000klambda"
  uvlimitp= "<1000klambda"
  
  
  outname1         = 'spw' + str(spw_id) + '.ms'
  outname_3C147    = 'spw' + str(spw_id) + '.3C147.ms'
  #outname_3C286    = 'spw' + str(spw_id) + '.3C286.ms'
  outname_3C138    = 'spw' + str(spw_id) + '.3C138.ms'
  #outname_3C48     = 'spw' + str(spw_id) + '.3C48.ms'
  #outname_leak     = 'spw' + str(spw_id) + '.leakage.ms'
  hsmooth_file     = 'hsmooth_' + outname1

  outname_phasecal = 'spw' + str(spw_id) + '.phasecal.ms'
  outname_bandpass = 'spw' + str(spw_id) + '.bandpass.ms'

  
  if numpy.int(spw_id) >= 10:
    outname_target     = 'spw' + str(spw_id) + '.target.ms'  
    outname_targetcorr = 'spw' + str(spw_id) + '.targetcorr.ms'
  else:
    outname_target     = 'spw0' + str(spw_id) + '.target.ms'  
    outname_targetcorr = 'spw0' + str(spw_id) + '.targetcorr.ms'
    
 
  calgceff       = 'cal.gceff.spw.' + str(spw_id)
  calantpos      = 'cal.antpos.spw.'+ str(spw_id)
  calbp1         = 'cal.bcal1.spw.' + str(spw_id)
  calgcal1       = 'cal.gcal1.spw.' + str(spw_id)  
  calgcal2       = 'cal.gcal2.spw.' + str(spw_id) 
  calgcal3       = 'cal.gcal3.spw.' + str(spw_id)
  calgcalall     = 'cal.gcalall.spw.' + str(spw_id)
  calgcal4       = 'cal.gcal4.spw.' + str(spw_id)
  calbp2         = 'cal.bcal2.spw.' + str(spw_id)
  caldf1         = 'cal.dcal1.spw.' + str(spw_id)
  calxf1         = 'cal.xcal1.spw.' + str(spw_id)
  calk1          = 'cal.kcal1.spw.' + str(spw_id)
  calkross1      = 'cal.kross1.spw.'+ str(spw_id)
  #calgcal1_3C286 = 'cal.gcal3C286.spw.' + str(spw_id)
  calgcal1_3C147 = 'cal.gcal3C147.spw.' + str(spw_id)
  calgcal1_3C138 = 'cal.gcal3C138.spw.' + str(spw_id)  
  calflux        = 'cal.fluxcal.spw.' + str(spw_id)  
  calgcal_phasecal  = 'cal.gcalphasecal.spw.' + str(spw_id) 

  os.system('rm -rf ' + outname1)
  os.system('rm -rf ' + outname_3C147)
  os.system('rm -rf ' + outname_3C138)  
  #os.system('rm -rf ' + outname_3C286)
  #os.system('rm -rf ' + outname_leak)
  os.system('rm -rf ' + hsmooth_file) 
  os.system('rm -rf ' + outname_phasecal)
  os.system('rm -rf ' + outname_bandpass)
  os.system('rm -rf ' + outname_target)
  os.system('rm -rf ' + outname_targetcorr)
  os.system('rm -rf ' + calgceff)
  os.system('rm -rf ' + calflux)
  os.system('rm -rf ' + calgcalall)
  os.system('rm -rf '+calbp1+' '+calbp2+' '+calgcal1+' '+calgcal2+' '+caldf1+' '+ calgcal4+' '+\
             calxf1+' '+calantpos+' '+calk1+' '+calkross1+' '+calgcal3)
  os.system('rm -rf all_calibrators_spw.' + str(spw_id) + '.ms')

  os.system('rm -rf ' + allspw_dataset)
  os.system('rm -rf ' + allspw_target)
  
  
  #split of the spw  
  default("split")
  split(vis=basevis,outputvis=outname1,datacolumn="data",field="",spw=spw_id)
  
  #apply the hanning smooth because of the RFI
  hanningsmooth(vis=outname1,datacolumn="data",outputvis=hsmooth_file)
  
  os.system('rm -rf ' + outname1)
  
  # FLAG SHADOW
  flagdata(vis=hsmooth_file,mode="shadow",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
          uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",           \
           clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0, \
           quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
           ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
           maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
           freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
           growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
           maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
           flagbackup=False,savepars=False,cmdreason="",outfile="")        

  # FLAG OF THE FIRST SCAN
  flagdata(vis=hsmooth_file,mode="manual",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
        uvrange="",timerange="",correlation="",scan="1",intent="",array="",observation="",feed="",           \
         clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0, \
         quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
         ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
         maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
         freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
         growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
         maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
         flagbackup=False,savepars=False,cmdreason="",outfile="")        


  print 'make copy of FLAG table'
  tb.open(hsmooth_file, nomodify=False)
  flags = numpy.copy(tb.getcol('FLAG'))
 
  flags_new = numpy.copy(flags)
  idx = numpy.where(flags_new == True)
  flags_new  [idx] = False # remove all flags
  tb.putcol('FLAG', flags_new)
  tb.flush()
  tb.close()
  
  # TFCROP FLAG
  flagdata(vis=hsmooth_file,mode="tfcrop",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   
           uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",                       
           clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,
           quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  
           ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        
           maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         
           freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, 
           growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      
           maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",
           flagbackup=False,savepars=False,cmdreason="",outfile="")
  
  #######################################################################################################################################
  #manual flagging of the ea24 -- after looking at polarisation solutions
  #badantenna='ea24'    
  #flagdata(vis=hsmooth_file,mode="manual",autocorr=False,inpfile="",reason="any",spw="",field="",antenna=badantenna,   
           #uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",                       
           #clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,
           #quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  
           #ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        
           #maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         
           #freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, 
           #growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      
           #maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",
           #flagbackup=False,savepars=False,cmdreason="",outfile="") 
   ####################################################################################################################################### 
  

  print 'restore FLAG table'
  tb.open(hsmooth_file, nomodify=False)
  flags = numpy.copy(tb.getcol('FLAG'))
  flags[idx] = True
  tb.putcol('FLAG', flags)
  tb.flush()
  tb.close()
  
  # GENCAL GAIN CURVE 
  gencal(vis=hsmooth_file,caltable=calgceff,caltype="gceff",spw="",antenna="",pol="",parameter=[])
  
  
	# ANTPOS
  gencal(vis=hsmooth_file,caltable=calantpos,caltype='antpos',spw="",antenna="",pol="",parameter=[])
  
  
  if os.path.isdir(calantpos):  
     print 'Antenna position corrections found'
     applycal(vis=hsmooth_file,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",  \
             observation="",msselect="",gaintable=[calgceff,calantpos],gainfield=[''],interp=['',''],spwmap=[],\
             parang=False,calwt=False,applymode="",flagbackup=False)
  else:  
     print 'Antenna position corrections NOT found'
     applycal(vis=hsmooth_file,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="",  \
             observation="",msselect="",gaintable=calgceff,gainfield=[''],interp=['',''],spwmap=[],\
             parang=False,calwt=False,applymode="",flagbackup=False)    

  ## bandpass [3C147 + 3C286]
  #default("split")
  #split(vis=hsmooth_file,outputvis=outname_bandpass,datacolumn="corrected",field="1331+305=3C286,0542+498=3C147")
     
  # bandpass [3C286]
  #default("split")
  #split(vis=hsmooth_file,outputvis=outname_bandpass,datacolumn="corrected",field="1331+305=3C286")
  
  # bandpass [3C147 & 3C138]
  default("split")
  split(vis=hsmooth_file,outputvis=outname_bandpass,datacolumn="corrected",field="3C138,0542+498=3C147")

  # [3C138]
  default("split")
  split(vis=hsmooth_file,outputvis=outname_3C138,datacolumn="corrected",field='3C138')
        
        
  # 3C147
  default("split")
  split(vis=hsmooth_file,outputvis=outname_3C147,datacolumn="corrected",field="0542+498=3C147")
 
  # 3C286 and bandpass
  #default("split")
  #split(vis=hsmooth_file,outputvis=outname_3C286,datacolumn="corrected",field="1331+305=3C286", scan="")

  
  # leakage calibrator J1407+2827
  #default("split")
  #split(vis=hsmooth_file,outputvis=outname_leak,datacolumn="corrected",field="J1407+2827", scan="")

  
  #Phasecal, FIELD_ID TAKES CARE OF DUMMY
  default("split")
  split(vis=hsmooth_file,outputvis=outname_phasecal,datacolumn="corrected",field=phasecalname, scan="")
    
  # target
  default("split")
  split(vis=hsmooth_file,outputvis=outname_target,datacolumn="corrected",field=targetname) 
 
  
  os.system('rm -rf ' + hsmooth_file)   
 
  #SETJY BANDPASS
  setjy(vis=outname_bandpass, field="0542+498=3C147",spw="",selectdata=False,timerange="",scan="", 
        observation="",model="3C147_L.im",listmodels=False,scalebychan=True,
        fluxdensity=-1,spix=0,reffreq="1GHz",standard="Perley-Butler 2013",
        useephemdir=False,usescratch=True)
 
  setjy(vis=outname_bandpass, field='3C138',spw="",selectdata=False,timerange="",scan="", 
       observation="",model="3C138_L.im",listmodels=False,scalebychan=True,
       fluxdensity=-1,spix=0,reffreq="1GHz",standard="Perley-Butler 2013",
       useephemdir=False,usescratch=True)

  #SETJY 3C147
  setjy(vis=outname_3C147, field="0542+498=3C147",spw="",selectdata=False,timerange="",scan="", \
        observation="",model="3C147_L.im",listmodels=False,scalebychan=True,    \
        fluxdensity=-1,spix=0,reffreq="1GHz",standard="Perley-Butler 2013",\
        useephemdir=False,usescratch=True)

  #SETJY 3C138
  setjy(vis=outname_3C138, field='3C138',spw="",selectdata=False,timerange="",scan="", \
        observation="",model="3C138_L.im",listmodels=False,scalebychan=True,    \
        fluxdensity=-1,spix=0,reffreq="1GHz",standard="Perley-Butler 2013",\
        useephemdir=False,usescratch=True)


  ##SETJY leakage TBC
  #setjy(vis=outname_leak, field="J1407+2827",spw="",selectdata=False,timerange="",scan="", \
  #      observation="",model="",listmodels=False,scalebychan=True,    \
  #      fluxdensity=[0.218,0,0,0],spix=[-0.755],reffreq="1GHz",standard="manual",\
  #      useephemdir=False,usescratch=True)

  #SETJY 3C286
  #setjy(vis=outname_3C286, field="1331+305=3C286",spw="",selectdata=False,timerange="",scan="", \
  #      observation="",model="3C286_L.im",listmodels=False,scalebychan=True,    \
  #      fluxdensity=-1,spix=0,reffreq="1GHz",standard="Perley-Butler 2013",\
  #      useephemdir=False,usescratch=True) # usescratch=T for now !!!

  #Going back to original data
  tb.open(outname_3C138, nomodify=False)
  data = tb.getcol('MODEL_DATA')
  
  RLphase = -11.0*2.0 #Polarisation angle times 2
  polfrac = 0.075 # at 1.45 GHz
  
 
  StokesI = (0.5*(numpy.copy(data[0,:,:]) + numpy.copy(data[3,:,:])))
  Q = StokesI*polfrac*numpy.cos(RLphase*pi/180.0)
  U = StokesI*polfrac*numpy.sin(RLphase*pi/180.0)
  RL = (Q + (complex(0,1)*U))
  LR = (Q - (complex(0,1)*U))
  
  data[1,:,:] = RL
  data[2,:,:] = LR
  tb.putcol('MODEL_DATA', data)
  tb.flush()
  tb.close()

  #PRE-CALIBRATION:
  #initial phase calibration
  gaincal(vis=outname_bandpass,caltable=calgcal1,field="",spw="0:25~35",intent="",      \
          selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",\
          solint="4min",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable="",gainfield=[''],interp=['linear'],spwmap=[],\
          parang=True)

  
  figure_name='initph1_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal1, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
    plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  #initial bandpass calibration                         
  bandpass(vis=outname_bandpass, caltable=calbp1,field="",spw="",intent="",selectdata=False,timerange="",\
           uvrange="",antenna="",scan="",observation="",msselect="",solint="inf",combine="scan,field",     \
           refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,bandtype="B",smodel=[],append=False,\
           fillgaps=0,degamp=3,degphase=3,visnorm=False,maskcenter=0,maskedge=5,gaintable=[calgcal1],  \
           gainfield=[''],interp=['linear,linear'],spwmap=[],parang=True)


  figure_name='bpassA1_'
  for i in range(num_plots+1):
      os.system('rm -rf '+figure_name+str(i+1))
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp1, xaxis='freq', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', 
                plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
    
  figure_name='bpassP1_'
  for i in range(num_plots+1):
      os.system('rm -rf '+figure_name+str(i+1))
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp1, xaxis='freq', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
              plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
              
 	##apply above solution to the individual sources for flagging purposes (rficonsole)
  applycal(vis=outname_3C147,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
           observation="",msselect="",gaintable=[calgcal1,calbp1],gainfield=[''],interp=['linear','linear,linear'],spwmap=[],\
           parang=False,calwt=False,applymode="",flagbackup=False)


  applycal(vis=outname_3C138,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
           observation="",msselect="",gaintable=[calgcal1,calbp1],gainfield=[''],interp=['linear','linear,linear'],spwmap=[],\
           parang=False,calwt=False,applymode="",flagbackup=False)

  #concat 3C147+3C138
  os.system('rm -rf ' + outname_bandpass)
  concat(vis=[outname_3C147, outname_3C138],concatvis=outname_bandpass,freqtol="",dirtol="",\
         respectname=False,timesort=True,copypointing=True,visweightscale=[])


  #applycal(vis=outname_3C286,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
  #         observation="",msselect="",gaintable=[calgcal1,calbp1],gainfield=[''],interp=['linear','linear,linear'],spwmap=[],\
  #         parang=False,calwt=False,applymode="",flagbackup=False)

  ##concat 3C147+3C286
  #os.system('rm -rf ' + outname_bandpass)
  #concat(vis=[outname_3C147, outname_3C286],concatvis=outname_bandpass,freqtol="",dirtol="",\
  #       respectname=False,timesort=True,copypointing=True,visweightscale=[])
  


  # NOW START THE ACTUAL CALIBRATION
  #phase calibration on the central channels
  gaincal(vis=outname_bandpass,caltable=calgcal2,field="",spw="0:25~35",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
          solint="5min",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable="",gainfield=[''],interp=[''],spwmap=[],\
          parang=True)
  
  figure_name='initph2_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal2, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # dynamic delay
  # use 3C147
  gaincal(vis=outname_bandpass,caltable=calk1,field="3C138, 0542+498=3C147", spw="0:2~61",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
          solint="inf",combine="scan,field",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="K",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable=[calgcal2],gainfield=[''],interp=['linear'],spwmap=[], \
          parang=True)
         
  figure_name='dynamic_delay'
  os.system('rm -rf '+figure_name)
  plotcal(caltable=calk1, xaxis='antenna', yaxis='delay', figfile=path_plot_spw+figure_name+'.png')

  # BANDPASS  
  gt1 = [calgcal2,calk1]
  bandpass(vis=outname_bandpass, caltable=calbp2,field="",spw="",intent="",selectdata=True,timerange="",\
           uvrange=uvlimit,antenna="",scan="",observation="",msselect="",solint="inf",combine="scan,field",     \
           refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,bandtype="B",smodel=[],append=False,\
           fillgaps=0,degamp=3,degphase=3,visnorm=False,maskcenter=0,maskedge=5,gaintable=gt1,  \
           gainfield="",interp=['linear','nearest'],spwmap=[],parang=True)

  figure_name='bpassA2_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calbp2, xaxis='freq', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', 
            plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
    
  figure_name='bpassP2_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calbp2, xaxis='freq', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  applycal(vis=outname_phasecal,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
            observation="",msselect="",gaintable=[calbp2],gainfield=[''],interp=['linear,linearflag'],spwmap=[],       \
            parang=False,calwt=False,applymode="",flagbackup=False)
                
           
  # GAINCAL 3C147 + 3C138
  #phase calibration on the whole bandwith (except for the first channels)
  gt2 = [calbp2,calk1]
  gaincal(vis=outname_bandpass,caltable=calgcal3,field="",spw="0:2~61",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
          solint="5min",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable=gt2,gainfield=[''],interp=['linear,linear','nearest'],spwmap=[], \
          parang=True)

  
  figure_name='finalphase'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal3, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='finalamp'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal3, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
    
  os.system('rm -rf ' + outname_bandpass)
  
  #GAINCAL 3C138
  gt3 = [calbp2,calk1]
  gaincal(vis=outname_3C138,caltable=calgcal1_3C138,field="",spw="0:2~61",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
          solint="30s",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable=gt3,gainfield=[''],interp=['linear,linear','nearest'],spwmap=[], \
          parang=True)
          
  figure_name='finalphase_3C138_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal1_3C138, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='finalamp_3C138_'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal1_3C138, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  #KCROSS 3C138 (need polarized source)
  #cross-hand delay: residual differences on the reference antenna

  gt31 = [calgcal1_3C138,calbp2,calk1]
  gaincal(vis=outname_3C138,caltable=calkross1,field="",spw="0:2~61",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimitp,antenna="",scan="",observation="",msselect="",\
          solint="inf",combine="scan",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="KCROSS",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable=gt31,gainfield=[''],interp=['linear','linear,linear','nearest'],spwmap=[], \
          parang=True)
          
  figure_name='KCROSS_3C138.png'
  os.system('rm -rf '+path_plot_spw+figure_name)
  plotcal(caltable=calkross1, xaxis='antenna', yaxis='delay',figfile=path_plot_spw+figure_name)
  
  # use J1407+2827 [3C147] for Df terms (if using phasecal solve with Df+QU)
  gt4 = [calgcal3,calbp2,calk1,calkross1]
  polcal(vis=outname_3C147,caltable=caldf1,field="",spw="",intent="",selectdata=True,timerange="",\
         uvrange=uvlimit,antenna="",scan="",observation="",msselect="",solint="inf",combine="scan",     \
         preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,poltype="Df",smodel=[],append=False,\
         gaintable=gt4,gainfield=[''],interp=['linear','linear,linear','nearest','nearest'],spwmap=[])


  figure_name='Df_real'
  for i in range(num_plots+1):
    os.system('rm -rf '+path_plot_spw+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=caldf1, xaxis='freq', yaxis='real',
            subplot=221, antenna=ant_plot, iteration='antenna',
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='Df_imag'
  for i in range(num_plots+1):
    os.system('rm -rf '+path_plot_spw+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=caldf1, xaxis='freq', yaxis='imag', subplot=221,
            antenna=ant_plot, iteration='antenna',
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # Local cal (amp for elevation dependent problems) 
  gt4_a = [calbp2,calk1,calkross1,caldf1]
  gaincal(vis=outname_phasecal,caltable=calgcal_phasecal,field="",spw="0:2~61",intent="",      \
          selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
          solint="5min",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
          gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
          phasewrap=180.0,gaintable=gt4_a,gainfield=[''],interp=['linear,linear','nearest','nearest','nearest,linear'],spwmap=[], \
          parang=True)
          
  figure_name='loc_cal_amp'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal_phasecal, xaxis='time', yaxis='amp',
            subplot=221, antenna=ant_plot, iteration='antenna',
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='loc_cal_phase'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcal_phasecal, xaxis='time', yaxis='phase', subplot=221,
            antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180],
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')  
  
  # Xf
  gt5 = [calgcal1_3C138,calbp2,calk1,calkross1,caldf1]
  polcal(vis=outname_3C138,caltable=calxf1,field="",spw="",intent="",selectdata=True,timerange="",\
         uvrange=uvlimitp,antenna="",scan="",observation="",msselect="",solint="inf",combine="scan",     \
         preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,poltype="Xf",smodel=[],append=False,\
         gaintable=gt5,gainfield=[''],interp=['linear','linear,linear','nearest','nearest','nearest,linear'],\
         spwmap=[])

  figure_name='Xf.png'
  os.system('rm -rf '+path_plot_spw+figure_name)
  plotcal(caltable=calxf1, xaxis='chan', yaxis='phase',figfile=path_plot_spw+figure_name)

  # NOW CONCAT ALL CALIBRATION SOLUTIONS
  concat(vis=[outname_phasecal,outname_3C147,outname_3C138],concatvis="all_calibrators_spw." + str(spw_id) + ".ms",\
               freqtol="",dirtol="",respectname=False,timesort=True,copypointing=True,visweightscale=[])

  #os.system('rm -rf ' + outname_3C147)
  #os.system('rm -rf ' + outname_3C138)
  #os.system('rm -rf ' + outname_3C286)
  #os.system('rm -rf ' + outname_phasecal)
 
  # GAINCAL phasecal, 3C147, 3C286    
  gaincal(vis="all_calibrators_spw." +str(spw_id) + ".ms",caltable=calgcalall,field="",spw="0:2~61",intent="",      \
         selectdata=True,timerange="",uvrange=uvlimit,antenna="",scan="",observation="",msselect="",\
         solint="5min",combine="",preavg=-1.0,refant=ref_ant,minblperant=4,minsnr=3.0,solnorm=False,\
         gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,npointaver=3,\
         phasewrap=180.0,gaintable=[calbp2,calk1,calkross1,caldf1,calxf1],gainfield=[''],       \
         interp=['linear,linear','nearest','nearest','nearest,linear','nearest,linear'],spwmap=[], parang=True)
  
  
  figure_name='finalphase_tot'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcalall, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
            plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='finalamp_tot'
  for i in range(num_plots+1):
    os.system('rm -rf '+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calgcalall, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna',
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')
          
  #"0137+331=3C48"
  #"1331+305=3C138"
  #"0542+498=3C147"

  # FLUXSCALE FOR phasecal (fitorder not "used" as we have only 1 spw)
  fluxscale(vis="all_calibrators_spw." +str(spw_id) + ".ms",caltable=calgcalall,fluxtable=calflux,reference=['3C138, 0542+498=3C147'],\
            transfer=phasecalname,listfile="",append=False,refspwmap=[-1],incremental=False,fitorder=1)
            

  figure_name='flux_amp'
  for i in range(num_plots+1):
    os.system('rm -rf '+path_plot_spw+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calflux, xaxis='time', yaxis='amp',
            subplot=221, antenna=ant_plot, iteration='antenna',
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='flux_phase'
  for i in range(num_plots+1):
    os.system('rm -rf '+path_plot_spw+figure_name+str(i+1))
    ant_plot=str(i*4)+'~'+str(i*4+3)
    plotcal(caltable=calflux, xaxis='time', yaxis='phase', subplot=221,       
            antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180],
            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  os.system('rm -rf all_calibrators_spw.' + str(spw_id) + '.ms')
  
  
  # STOP HERE WHEN TESTING THE GOOD SPW -- look for bad antennas in the polarisation plots
  
  gt7 = [calbp2,calk1,calkross1,caldf1,calxf1,calgcalall]
  
  if all(os.path.exists(gt) for gt in gt7):
    print "spw", str(spw_id), "ok"

    applycal(vis=outname_target,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
           observation="",msselect="",gaintable=gt7,gainfield=[''],\
           interp=['nearest','nearest','nearest','nearest','nearest','nearest'],spwmap=[],\
           parang=True,calwt=False,applymode="",flagbackup=False)

    #TFCROP FLAG
    flagdata(vis=outname_target,mode="tfcrop",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
           uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",           \
           clipminmax=[],datacolumn="corrected",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,\
           quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
           ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
           maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
           freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
           growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
           maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
           flagbackup=False,savepars=False,cmdreason="",outfile="")
           
    #RFLAG FLAG   
    flagdata(vis=outname_target,mode="rflag",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
        uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",  \
        clipminmax=[],datacolumn="corrected",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,\
        quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
        ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
        maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
        freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
        growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
        maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
        flagbackup=False,savepars=False,cmdreason="",outfile="")
      
  
     
    print 'Apply again, but now with the calibration flagging on'
    print 'Do not skip this final applycal!'
    applycal(vis=outname_target,field="",spw="",intent="",selectdata=False,timerange="",uvrange="",antenna="",scan="", \
           observation="",msselect="",gaintable=gt7,gainfield=[''],\
           interp=['linear','linear,linearflag','nearest','nearest','nearest,linearflag','nearest,linearflag'],spwmap=[],\
           parang=True,calwt=False,applymode="",flagbackup=False)

  else:
    print "spw", str(spw_id), "needs to be flagged"

    flagdata(vis=outname_target,mode="manual",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
        uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",  \
        clipminmax=[],datacolumn="corrected",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,\
        quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
        ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
        maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
        freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
        growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
        maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
        flagbackup=False,savepars=False,cmdreason="",outfile="")
  


#os.system('rm -rf ' + basevis)       

#put together all the spw in one dataset (that is still divided into 16 spw)
all_spwlist = sorted(glob.glob('spw??.target.ms'))
concat(vis=all_spwlist, concatvis=allspw_dataset)

#flag_stat(allspw_dataset)

#for spw_id in spwlist:
  #if numpy.int(spw_id) >= 10:
    #os.system('rm -rf spw' + str(spw_id) + '.target.ms')
    #os.system('rm -rf spw' + str(spw_id) + '.targetcorr.ms')
  #else:
    #os.system('rm -rf spw0' + str(spw_id) + '.target.ms')
    #os.system('rm -rf spw0' + str(spw_id) + '.targetcorr.ms')
        
#split and avareging; we ignore the initial and final channel because they are bad (total of 16 channels)
split(vis=allspw_dataset, outputvis=allspw_target, datacolumn='corrected', field=targetname, spw='*:7~54', width=4, \
                timebin='10s')

#os.system('rm -rf ' + allspw_dataset)
#run AOFlagger to better remove the RFI
AOFstrategy = '/lofar2/bba3185/flagonavg_rlandlr.rfis'
#subprocess.call('aoflagger -strategy ' + AOFstrategy + ' -column DATA ' + allspw_target, shell=True)
subprocess.call('/home/lofar/opt3/aoflagger/bin/aoflagger -strategy ' + AOFstrategy + ' -column DATA ' + allspw_target, shell=True)


#further flagging after rl and lr correlation inspection -- on spw
#plotms(vis=allspw_target, xaxis='uvdist',yaxis='amp',correlation='rl,lr',coloraxis='antenna1')


  
'''
#######################################################################################################################################
#manual flagging of the ea24 -- after looking at polarisation solutions
badspw='8'    
flagdata(vis=allspw_target,mode="manual",autocorr=False,inpfile="",reason="any",spw=badspw,field="",antenna="",   
         uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",                       
         clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,
         quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  
         ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        
         maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         
         freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, 
         growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      
         maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",
         flagbackup=False,savepars=False,cmdreason="",outfile="") 
####################################################################################################################################### 

'''
flag_stat(allspw_target)    
         
