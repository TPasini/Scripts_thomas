import numpy as np
from scipy import stats
from numpy import array
import itsample
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import csv
import random
from astropy.cosmology import FlatLambdaCDM
import os, sys, shutil
import itertools

#####################
###ITSAMPLE SCRIPT###
#####################
    
N=8000  #Number of data to sample
cosmo = FlatLambdaCDM(H0=71, Om0=0.3)
print('')
print('Using Flat LmbdaCDM H0=0.71 Om0=0.3.')
print('Sampling', N, 'values from the provided function. This could take a while..')

#Files containing the results are stored in the 'Results' directory
if os.path.exists('Results'):
    shutil.rmtree('Results')
os.makedirs('Results')


with open('Results/lx_z.txt', 'a') as f:
 writer = csv.writer(f, delimiter=" ")
 with open('Results/simulation.txt', 'a') as s:
  writer2 = csv.writer(s, delimiter=" ")
  writer2.writerow(["L_X(10^44 erg s-1)", "log L_R(erg s-1 Hz-1)", 'z'])
  writer2.writerow([])
  with open('Results/upperlimits.txt', 'a') as g:
   writer3 = csv.writer(g, delimiter=" ")
   writer3.writerow(["L_X(10^44 erg s-1)", "log L_R(erg s-1 Hz-1)", 'z'])
   writer3.writerow([])
   with open('Results/flux.txt', 'a') as h:
    writer4 = csv.writer(h, delimiter=" ")
    writer4.writerow(["F_X(erg s-1 cm-2)", "F_R(erg s-1 cm-2 Hz-1)", 'z'])
    writer4.writerow([])
    with open('Results/fluxupper.txt', 'a') as j:
      writer5 = csv.writer(j, delimiter=" ")
      writer5.writerow(["F_X(erg s-1 cm-2)", "F_R(erg s-1 cm-2 Hz-1)", 'z'])
      writer5.writerow([])
 
 
 
      for _ in itertools.repeat(None, N):
      
      ##SAME REDSHIFT SAMPLING THAN REAL DATA
      #Get the redshift distribution of real data
      #data=np.loadtxt('zlist.txt')
      #fit=stats.gaussian_kde(data)

      ##UNIFORM SAMPLING
      #redshift = np.random.uniform(0.01, 0.4, size=N)

       def fit(z):
        return (4*3.1415*3e+10*((cosmo.luminosity_distance(z).value)**2)) / ((1+z)*(cosmo.H0.value)*((cosmo.H(z)/cosmo.H(0))))
       redshift = itsample.sample(fit, 1, 0.01, 1.7)


##X-ray luminosity function clusters
       #def my_pdf(x):
        # return 5.06e-7*np.exp(-x/9.1)*x**(-1.85) #From BCS Ebeling+98, 0.1-2.4 keV
#
##X-ray luminosity function groups

       def my_pdf(x):
        A = 2.94e-7*((1+redshift)**(-1.2))
        B = 2.64*((1+redshift)**(-2))
        return A*np.exp(-x/B)*x**(-1.69)#From WARPS Koenns+18

#Plot of the function
#Uncomment if you want the plot of the given function
#
#       xvals = np.linspace(0.001,1000)
#       plt.plot(xvals,[my_pdf(x) for x in xvals], color='blue')
#       #plt.plot(xvals,[my_pdf2(x) for x in xvals], color='red')
#       plt.xscale('log')
#       plt.yscale('log')
#       plt.show()
        

##Sampling N luminosities from the given probability function (normalized and inverted) and association with z
       lx = itsample.sample(my_pdf, 1, 0.0008, 5.)  #0.0008 5 for groups

       writer.writerows(zip(lx,redshift))


#######################################
#######ONLY X-RAY FLUX LIMITED#########
#######################################

##From luminosity to flux
#data=np.loadtxt('Results/lx_z.txt')
#cosmo = FlatLambdaCDM(H0=71, Om0=0.3)
#xlum=data[:,0]
#z=data[:,1]
#dist = cosmo.luminosity_distance(z).value # in Mpc
#flux=(xlum*1e+44)/(4*math.pi*(dist*3.086e+24)**2)
##print(flux)
#
##Application of the flux limitations and association with random radio luminosity from uniform distribution
#limit = 2e-15 #erg s-1 cm-2
#with open('Results/simulation.txt', 'a') as f:
# for i,x in enumerate(flux):
#  if flux[i] > limit:
#   radio = np.random.uniform(27., 35.)
#   writer = csv.writer(f, delimiter=" ")
#   writer.writerow([lx[i], radio, redshift[i]])
#   #print(lx[i], radio, redshift[i])

################################################
######RADIO AND X-RAY FLUX LIMITED, NO CORR#####
################################################

##From luminosity to flux
#       xlum=lx
#       z=redshift
#       dist = cosmo.luminosity_distance(z).value # in Mpc
#       flux=(xlum*1e+44)/(4*math.pi*(dist*3.086e+24)**2)
#
#       limit = 2e-15 #erg s-1 cm-2   X-ray limit
#       limitradio=18e-6*1e-23 #erg s-1 cm-2 Hz-1   Radio limit
#
#
#       if flux > limit:
#          radio = np.random.uniform(27., 35.)
#          fluxr=(10**radio)/(4*math.pi*(dist*3.086e+24)**2)
#          if fluxr > limitradio:
#           writer2.writerow([xlum[0], radio, z[0]])
#           writer4.writerow([flux[0], fluxr, z[0]])
#          if fluxr < limitradio:
#           lumupper=limitradio*(4*math.pi*((dist*3.086e+24)**2))
#           writer3.writerow([xlum[0], lumupper[0], z[0]])
#           writer5.writerow([flux[0], fluxr[0], z[0]])
#
#print('')
#print("All the files are stored in the 'Results' directory")


################################################
######RADIO AND X-RAY FLUX LIMITED, CORR#####
################################################

#From luminosity to flux
       xlum=lx
       z=redshift
       dist = cosmo.luminosity_distance(z).value # in Mpc
       flux=(xlum*1e+44)/(4*math.pi*(dist*3.086e+24)**2)

#Groups
       limit = 2e-15 #erg s-1 cm-2   X-ray limit
       limitradio=18e-6*1e-23 #erg s-1 cm-2 Hz-1   Radio limit

#Clusters
       #limit = 4.4e-12 #erg s-1 cm-2   X-ray limit
       #limitradio=1.35e-3*1e-23 #erg s-1 cm-2 Hz-1   Radio limit


       if flux > limit:
        radioexact=1.0703*(np.log10(xlum*1e+44))-15.9  #Groups
        #radioexact=1.258*(np.log10(xlum*1e+44))-25.8014  #Clusters
        radio=np.random.normal(radioexact, scale=0.8)
        fluxr=(10**radio)/(4*math.pi*(dist*3.086e+24)**2)
        if fluxr > limitradio:
         writer2.writerow([lx[0], radio[0], z[0]])
         writer4.writerow([flux[0], fluxr[0], z[0]])
        if fluxr < limitradio:
         lumupper=limitradio*(4*math.pi*((dist*3.086e+24)**2))
         writer3.writerow([xlum[0], lumupper[0], z[0]])
         writer5.writerow([flux[0], fluxr[0], z[0]])

print('')
print("All the files are stored in the 'Results' directory")



##############
#####PLOT#####
##############

simulation=np.loadtxt('Results/simulation.txt', skiprows=2)
uppers=np.loadtxt('Results/upperlimits.txt', skiprows=2)
cm = plt.cm.get_cmap('jet') #jet
arrowleft=  u'$\u2190$'

plt.scatter(10**(simulation[:,1]), simulation[:,0]*1e+44, c=simulation[:,2], cmap=cm, vmin = 0.01, vmax = 2, label='Simulation')
plt.scatter(uppers[:,1], uppers[:,0]*1e+44, c=uppers[:,2], cmap=cm, vmin = 0.01, vmax = 2, label='Upper limits', marker=arrowleft)
plt.xlim(5e+26,1.2e+35)
plt.ylim(1e+41,1e+46)
cbar=plt.colorbar()
cbar.set_label('Redshift', labelpad=20, fontsize=17)
plt.clim([0.01,2])

plt.tick_params(labelsize=19)
cbar.ax.tick_params(labelsize=17)
plt.ylabel(r'L$_{X}$ (ergs s$^{-1}$)', fontsize=20)
plt.xlabel(r'P$_{1.4GHz}$ (ergs s$^{-1}$ Hz$^{-1}$)', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.legend(prop={'size': 15})
plt.show()
