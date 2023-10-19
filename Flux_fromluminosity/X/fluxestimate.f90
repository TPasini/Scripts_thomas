!This program estimates the X-ray flux starting from the luminosity in erg/s and lum distance

PROGRAM fluxestimate
  IMPLICIT NONE
  REAL*8:: pi, freq
  REAL*8, ALLOCATABLE:: dl(:), flux(:), lum(:)
  INTEGER:: ndati, i
  
  pi=3.1415

   OPEN(30, file='data.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(dl(ndati), flux(ndati), lum(ndati))

    REWIND(30)

    OPEN(40, file='Fluxes_results.dat')
                                                     
 DO i=1,ndati
    READ(30,*) dl(i), lum(i)
    flux(i)=lum(i)/(4*pi*((dl(i)*3d+24)**2))
    WRITE(40,*) flux(i)
 END DO

 END PROGRAM fluxestimate
