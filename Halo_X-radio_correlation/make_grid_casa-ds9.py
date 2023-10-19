#Put boxes in casa or ds9 region format on an image. Boxes are squares, the strating points is p1 (blc, trc) dx and dy are the increment for the next box. WARNING: it assumes images have NAXIS=4 as LOFAR images do (in the fct pix2world)
#### p are intended as vertex of the box as follows:
#
#   p2            p3
#
#
#   p1           p4
#
# It starts from bottom left (p1) and adds boxes down to bottom right if the image taking p4[0] > pf_3[0] as condition to stop
# The number of boxes put in the first line is recorded, it restarts from a box above p2 and goes on along the same line
# (from bottom left to bottom right) putting the same number of boxes as in the first line.
# It stops when p1[1] > p3_f[1]

from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy
import os 
from astropy.wcs import WCS


#*************************start inputs*******************************

#Image you want to put boxes on (boxes will be written in J2000 coordinate system)
fits_in='image_allMS_SHI_FTh_Scale1000_no0scale_4_mask_35asec.int.restored.fits'

#Directory where you want the region files to be written. It writes one region file per box boxX_Y.reg where X is the line, Y is the column. It also writes a Allregion_casa.reg or Allregion_ds9.reg where all the boxes are written
dir_out='PROVA_ds9'


#starting points (xblc,yblc):
if dir_out=='PROVA_ds9':
  p1=(1264, 1438)
  p_f3=(2075,1689)
  dx=100.
  dy=71.


#format of region file either ds9 or CASA
format='ds9'


#*************************end inputs*******************************



p2=[p1[0],p1[1]+dy]


w=WCS(fits_in)




p3=[p1[0]+dx,p1[1]+dy]
p4=[p1[0]+dx,p1[1]]


p1_up=p2
p4_up=p3

p2_up=[p1_up[0],p1_up[1]]
p3_up=[p1_up[0]+dx,p1_up[1]+dy]
   

ii=0
jj=0
count=0


if os.path.exists(dir_out+"/Allregion_CASA.reg" or dir_out+"/Allregion_ds9.reg"):
           print("Files with regions already exists. I'll stop here",dir_out+"/Allregion.reg")
           sys.exit()

if format=='CASA':
  file_in=fits_in.replace('fits','im')

  if not os.path.exists(file_in):
    importfits(fitsimage=fits_in,imagename=file_in)

  fa=open(dir_out+"/Allregion_CASA.reg", "w")
  fa.write('#CRTFv0 CASA Region Text Format version 0 \n')



if format=='ds9':
  intab="hmds" 
  outtab="::: "
  trantab = str.maketrans(intab, outtab)
  fa=open(dir_out+"/Allregion_ds9.reg", "w")
  fa.write('#  Region file format: DS9 version 4.1 \n')
  fa.write('#  Filename '+fits_in+' \n')
  fa.write('global color=green dash=0 dashlist=8 3 delete=1 edit=1 fixed=0 font="helvetica 10 normal roman" highlite=1 include=1 move=1 select=1 source=1 width=1 \n')
  fa.write('fk5 \n')


while (p1_up[1] <= p_f3[1])  :

    ii=0
       
    while (p4[0] <= p_f3[0]) or (jj  >0 and ii < count):

       p1_d=w.wcs_pix2world(p1[0],p1[1],0,0,0)
       p2_d=w.wcs_pix2world(p2[0],p2[1],0,0,0)
       p3_d=w.wcs_pix2world(p3[0],p3[1],0,0,0)
       p4_d=w.wcs_pix2world(p4[0],p4[1],0,0,0)

       p1w=SkyCoord(p1_d[0],p1_d[1],unit='deg')
       p2w=SkyCoord(p2_d[0],p2_d[1],unit='deg')
       p3w=SkyCoord(p3_d[0],p3_d[1],unit='deg')
       p4w=SkyCoord(p4_d[0],p4_d[1],unit='deg')


       p1hms=p1w.to_string('hmsdms')
       p1_l=''.join(p1hms)
       p1_reg=p1_l.replace('+',', +')

       p2hms=p2w.to_string('hmsdms')
       p2_l=''.join(p2hms)
       p2_reg=p2_l.replace('+',', +')

       p3hms=p3w.to_string('hmsdms')
       p3_l=''.join(p3hms)
       p3_reg=p3_l.replace('+',', +')

       p4hms=p4w.to_string('hmsdms')
       p4_l=''.join(p4hms)
       p4_reg=p4_l.replace('+',', +')
   
       
       if format=='CASA':
         string_reg="poly [["+p1_reg+"]," + "["+p2_reg+"]," + "["+p3_reg+"]," + "["+p4_reg+"]]"
         fo=open(dir_out+"/box"+str(jj)+'_'+str(ii)+"_casa.reg", "w")
         fo.write('#CRTFv0 CASA Region Text Format version 0 \n')
         fo.write(string_reg)
         fo.write (" coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=green, font=Helvetica, fontsize=11, fontstyle=normal, usetex=false \n")
         fo.close()
         
         fa.write(string_reg)
         fa.write (" coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=green, font=Helvetica, fontsize=11, fontstyle=normal, usetex=false \n")

       
       if format=='ds9':         
         string_reg="polygon("+p1_reg.translate(trantab)+","+p2_reg.translate(trantab)+","+p3_reg.translate(trantab)+"," +p4_reg.translate(trantab)+")"
         fo=open(dir_out+"/box"+str(jj)+'_'+str(ii)+"_ds9.reg", "w")
         fo.write('#  Region file format: DS9 version 4.1 \n')
         fo.write('#  Filename '+fits_in+' \n')
         fo.write('global color=green dash=0 dashlist=8 3 delete=1 edit=1 fixed=0 font="helvetica 10 normal roman" highlite=1 include=1 move=1 select=1 source=1 width=1 \n')
         fo.write('fk5 \n')
         fo.write(string_reg)
         fo.close()


         fa.write(string_reg+' \n')

       p1=p4
       p2=p3

       p4=[p1[0]+dx,p1[1]]
       p3=[p1[0]+dx,p1[1]+dy]

       print('writing box ', "box"+str(jj)+'_'+str(ii)+"_"+format+".reg in "+dir_out)

       if (jj==0):
           count=count+1
           
       ii=ii+1
    
   
    p1=p1_up
    p2=[p1_up[0],p1_up[1]+dy]
    p3=[p1_up[0]+dx,p1_up[1]+dy]
    p4=[p1_up[0]+dx,p1_up[1]]


    p1_up=p2_up
    p4_up=p3_up
    p2_up=[p1_up[0],p1_up[1]+dy]
    p3_up=[p1_up[0]+dx,p1_up[1]+dy]
    
    jj=jj+1


fa.close()
