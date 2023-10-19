PROGRAM columnmulti
IMPLICIT NONE 
REAL*8:: constant
INTEGER:: i, ndati, answer
REAL*8, ALLOCATABLE:: number(:), numberfinal(:), number2(:)

15 print*, 'What do you want to do?'
print*, '1. Multiply an array for a constant'
print*, '2. Multiply an array for another array'
read(*,*) answer

if(answer==1) then
   
print*, 'What constant do you want your array multiplied for?'
read(*,*) constant

OPEN(30, file='datacol.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(number(ndati), numberfinal(ndati))

    REWIND(30)

    OPEN(40, file='results.dat')
                                                     
 DO i=1,ndati
    READ(30,*) number(i)
    numberfinal(i)=number(i)*constant
    WRITE(40,*) numberfinal(i)
 END DO

else if(answer==2) then
   print*, 'The arrays are read from cols.dat'

   OPEN(30, file='cols.dat')
     i=0

      DO
        READ(30,*,end=110)
        i=i+1
      END DO

      110 CONTINUE
      ndati=i

    ALLOCATE(number(ndati), numberfinal(ndati), number2(ndati))

    REWIND(30)
    
   OPEN(40, file='results.dat')
                                                     
 DO i=1,ndati
    READ(30,*) number(i), number2(i)
    numberfinal(i)=number(i)*number2(i)
    WRITE(40,*) numberfinal(i)
 END DO

else
   print*, 'Print 1 or 2'
   go to 15

end if


END PROGRAM columnmulti
