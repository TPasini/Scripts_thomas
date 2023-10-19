PROGRAM spindex_est
IMPLICIT NONE 
REAL*8, ALLOCATABLE:: f1(:), f2(:), err1(:), err2(:), alpha(:), erralpha(:)
REAL*8, ALLOCATABLE:: lum1(:), lum2(:)
INTEGER:: answer, i, ndati, answer2
REAL*8:: nu1, nu2, spindex, freq1, freq2

10 print*, 'What do you want to do?'
print*, '1. Estimate the synchrotron index'
print*, '2. Estimate a flux/luminosity assuming a sp. index'
read(*,*) answer2

IF(answer2==1) then

print*,''
print*, '##################################'
print*, '###ESTIMATE OF SYNCHROTRON INDEX##'
print*, '##################################'
print*, ''

print*, 'The spectral index is here defined as 𝛼 in a spectrum S∝𝜈^(-𝛼), where 𝜈 is the frequency.' 
print*, 'The software reads the flux measurements from a file that must be called forspindex.dat.'
print*, 'The columns must be flux1, err, flux2, err.'
print*,''

PRINT*, 'What is the frequency in Hz to which the first column of the file is referred to?'
read(*,*) nu1
PRINT*, 'What is the other frequency?'
read(*,*) nu2


OPEN(30, file='forspindex.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(f1(ndati), f2(ndati),err1(ndati), err2(ndati), alpha(ndati), erralpha(ndati))

    REWIND(30)

    OPEN(40, file='spindex_results.dat')
    WRITE(40,*) '# Sp.index    error'
                                                     
 DO i=1,ndati
   READ(30,*) f1(i), err1(i), f2(i), err2(i)
   alpha(i)=(log10(f1(i)/f2(i)))/(log10(nu2/nu1))
   erralpha(i)=(1/(log(nu2)-log(nu1)))*(sqrt(((err1(i)/f1(i))**2)+((err2(i)/f2(i))**2)))
   write(40,'(F10.2,F10.2)') alpha(i), abs(erralpha(i))
   
END DO


else if(answer2==2) then

print*,''
print*, '###################################'
print*, '##ESTIMATE OF FLUXES/LUMINOSITIES##'
print*, '###################################'
print*, ''

   print*, 'The spectral index is here defined as 𝛼 in a spectrum S∝𝜈^(-𝛼), where 𝜈 is the frequency.'
   print*, 'What is the assumed spectral index of the radio sources?'
   read(*,*) spindex

   print*, 'The fluxes/luminosities are read from the file datalum.dat'
   print*, 'The file must contain the fluxes/luminosities at the starting frequency'
   print*, 'The result will be in the same units as the given data'

   print*, 'What is the frequency in Hz of the given data?'
   read(*,*) freq1
   print*, 'What is the frequency in Hz you want to estimate the luminosity of?'
   read(*,*) freq2

   OPEN(30, file='datalum.dat')
     i=0

      DO
        READ(30,*,end=200)
        i=i+1
      END DO

      200 CONTINUE
      ndati=i

    ALLOCATE(lum1(ndati), lum2(ndati))

    REWIND(30)

    OPEN(40, file='lum_results.dat')
                                                     
 DO i=1,ndati
   READ(30,*) lum1(i)
   lum2(i)=((freq2/freq1)**(-spindex))*lum1(i)
   write(40,'(E10.5)') lum2(i) 
   
END DO

ELSE
   GO TO 10

END IF

   
   

END PROGRAM spindex_est
