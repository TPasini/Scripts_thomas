
PROGRAM kcorr
  IMPLICIT NONE
  REAL*8:: pi, freq
  REAL*8, ALLOCATABLE:: dl(:), flux(:), lum(:), z(:), kcor(:)
  INTEGER:: ndati, i, ans
  
   print*, 'File must contain luminosity and redshift'
   OPEN(30, file='data.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(dl(ndati), flux(ndati), lum(ndati), z(ndati), kcor(ndati))

    REWIND(30)

    OPEN(40, file='kcorrected.dat')
                                                     
 DO i=1,ndati
    READ(30,*) lum(i), z(i)
    kcor(i)=lum(i)*((1+z(i))**(-0.2))
    WRITE(40,'(F10.5,A)') kcor(i)
 END DO

 END PROGRAM kcorr
