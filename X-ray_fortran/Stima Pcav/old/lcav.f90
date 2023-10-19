PROGRAM lcav
IMPLICIT NONE
REAL*8:: pi
INTEGER:: i, ndati
REAL*8, ALLOCATABLE:: pv(:), t(:), l(:), ts(:), pvexp(:), errpvsu(:), errpvgiu(:)
REAL*8, ALLOCATABLE:: errelsu(:), errelgiu(:), errabssu(:), errabsgiu(:)

PRINT*, "Il programma legge i dati dal file pv.dat, in cui l'entalpia deve essere in unità di 10^57 e t in yr."
PRINT*, 'I risultati si trovano nel file res.dat in unità di 10^42.'
OPEN(10, file='pv.dat')

pi=3.1415
i=0

      DO
        READ(10,*,end=100)
        i=i+1
      END DO
 
100 CONTINUE
      ndati=i


ALLOCATE(pv(ndati), t(ndati), l(ndati), ts(ndati), pvexp(ndati))
ALLOCATE(errpvsu(ndati), errpvgiu(ndati), errelsu(ndati), errelgiu(ndati))
ALLOCATE(errabssu(ndati), errabsgiu(ndati))

REWIND(10)
OPEN(20, file='res.dat')

DO i=1, ndati
 READ(10,*) pv(i), t(i), errpvsu(i), errpvgiu(i)
 errelsu(i)=errpvsu(i)/pv(i)
 errelgiu(i)=errpvgiu(i)/pv(i)
 pvexp(i)=pv(i)*(1d+57)
 ts(i)=t(i)*(1e+7)*pi*(1e+7)
 l(i)=(4*pvexp(i))/ts(i)
 errabssu(i)=(errelsu(i))*l(i)
 errabsgiu(i)=(errelgiu(i))*l(i)
 WRITE(20,*) l(i)/(1d+42), 'erg/s', errabssu(i)/(1d+42), errabsgiu(i)/(1d+42)
END DO

END PROGRAM Lcav
