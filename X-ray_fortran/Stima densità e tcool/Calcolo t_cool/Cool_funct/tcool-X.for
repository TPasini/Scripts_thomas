      program t_coolX
c     by Myriam Gitti
c
      real*8 ab(4),c(2,4),tem(2),da,kev,g,mu,mp,const,krc,ne,kt
      real*8 tcool
      integer cosmo
c
      data ab/0.,0.1,0.31623,1./
c
c
      print*,'==================================================' 
      print*,'This program calculates the cooling time as:'
      print*,'t_cool = (5/2 * kT) / (mu * X * n-e * Lambda),' 
      print*,'where mu=0.6, X=0.7, and Lambda is the cooling function'
      print*,'=================================================='

 10   continue

      print*,'Insert temperature (keV)'
      read*,kt
      print*,'Insert metallicity (solar)'
      read*,abund
      print*,'Insert electron density (cm^-3)'
      read*,ne 
      print*,'Insert redshift'
      read*,z
      open(10,file='cool-tab.dat',status='old')

      temp = kt*1.16e7
c
c-----  Interpolate emissivity from tables of Sutherland & Dopita 1993, ApJS,
c-----  Lambda_net/ne^2 (erg cm^3/s) in file "cool-tab.dat"
c-----  for metallicities Z between 0.0 - 1.0 and log(T) [K] between 6.5 - 8.5 
c
c  
c-----  Find abundance indices
c
      if((abund.gt.1.).or.(abund.lt.0.)) 
     #                        stop "metallicity out of range"
      if(abund.gt..31623) then
         imlow=3
         imhig=4
         goto 100
      endif
      if (abund.gt..1) then
         imlow=2
         imhig=3
      else 
         imlow=1
         imhig=2
      endif
 100  continue
c
c-----  Find temperature indices
c
      if((alog10(temp).gt.8.5).or.(alog10(temp).lt.6.5))
     #                             stop "temperature out of range"
      do 200 i=1,41
         read(10,*) tem(2),c(2,1),c(2,2),c(2,3),c(2,4)
         if (10**tem(2).gt.temp) goto 210
         tem(1) = tem(2)
         do 205 j=1,4
            c(1,j) = c(2,j)
 205     continue
 200  continue
 210  continue
      close (10)
c
c-----  Interpolation
c
      do 250 i=1,2 
         tem(i) = 10**tem(i)
         c(i,imlow) = 10**c(i,imlow) 
         c(i,imhig) = 10**c(i,imhig)
 250  continue
c
      par1 = (c(2,imlow) - c(1,imlow)) / (tem(2) - tem(1))
      par2 = c(1,imlow) - par1*tem(1)
      clow = par1*temp + par2
c
      par1 = (c(2,imhig) - c(1,imhig)) / (tem(2) - tem(1))
      par2 = c(1,imhig) - par1*tem(1)
      chig = par1*temp + par2
c
      par1 = (chig - clow) / (ab(imhig) - ab(imlow))
      par2 = clow - par1*ab(imlow)
c
c-----  Emissivity bolometric
c
      xcbol = par1*abund + par2 
      print*,'===> Cooling function from Sutherland & Dopita 1993:' 
      print*,'Lambda (erg cm^3 s^-1) = ',xcbol    
      print*,' '               
c
c        tempo di cooling (formula Morris & Fabian 2005)
       tcool = ((2.5*1.60217653d-9*kt)/(0.61*0.71*ne*xcbol))/3.1536d7
c         tempo di cooling (formula Pratt et al. 2002)
c       tcool2 = (2.9d10 * sqrt(kt))/(0.82*ne/1d-3)

999   continue                                                          
      print*,'===> The cooling time (Gyr) is:',tcool/1e9
      print*,'=================================================='

c
      stop
      end                                                               
