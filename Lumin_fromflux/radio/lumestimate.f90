!This program estimates the luminosity starting from the lum distance in Mpc and the flux in mJy (no k-correction)s

PROGRAM lumestimate
  IMPLICIT NONE
  REAL*8:: pi, freq
  REAL*8, ALLOCATABLE:: dl(:), flux(:), lum(:), z(:)
  INTEGER:: ndati, i, ans
  
  pi=3.1415

  Print*, 'Print 1 for K-correction, 2 for no K-correction'
  read(*,*) ans
  if(ans==1) then

   print*, 'File must contain lum. distance, flux and redshift'
   OPEN(30, file='data.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(dl(ndati), flux(ndati), lum(ndati), z(ndati))

    REWIND(30)

    OPEN(40, file='Lum_results.dat')
                                                     
 DO i=1,ndati
    READ(30,*) dl(i), flux(i), z(i)
    lum(i)=(((flux(i)*1e-3)*1d-23)*4*pi*((dl(i)*3d+24)**2)*((1+z(i))**(-0.2)))
    WRITE(40,'(F10.5,A)') lum(i)/(1d+32), ' *10^32 erg/s/Hz'
 END DO

else if(ans==2) then

      OPEN(30, file='data.dat')
     i=0

      DO
        READ(30,*,end=110)
        i=i+1
      END DO

      110 CONTINUE
      ndati=i

    ALLOCATE(dl(ndati), flux(ndati), lum(ndati))

    REWIND(30)

    OPEN(40, file='Lum_results.dat')
                                                     
 DO i=1,ndati
    READ(30,*) dl(i), flux(i)
    lum(i)=((flux(i)*1e-3)*1d-23)*4*pi*((dl(i)*3d+24)**2)
    WRITE(40,'(F10.5,A)') lum(i)/(1d+32), ' *10^32 erg/s/Hz'
 END DO


 end if
 

 END PROGRAM lumestimate
