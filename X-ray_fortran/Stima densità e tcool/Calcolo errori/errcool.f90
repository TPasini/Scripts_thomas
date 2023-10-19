PROGRAM errors
  IMPLICIT NONE
  REAL*8, ALLOCATABLE:: errTsu(:), T(:), errnsu(:), n(:), tcool(:),relTsu(:), relnsu(:), relcoolsu(:), errcoolsu(:), kT(:)
  REAL*8, ALLOCATABLE:: errTgiu(:), errngiu(:),relTgiu(:), relngiu(:), relcoolgiu(:), errcoolgiu(:), errkTsu(:), errkTgiu(:)
  INTEGER:: i, j, ndati

  PRINT*, 'I dati vengono letti dai file errsu.dat e errgiu.dat.'
  PRINT*, 'Devono contenere:'
  PRINT*, '- Errore di kT in keV;'
  PRINT*, '- Valore di kT in keV;'
  PRINT*, '- Errore di n in cm-3;'
  PRINT*, '- Valore di n in cm-3;'
  PRINT*, '- Stima di t_cool in yr.'
  PRINT*, 'I risultati ordinati si trovano nel file Results_err.dat.'
  PRINT*, 'Il file Printerr.dat contiene di seguito tutti gli errori superiori e poi tutti quelli inferiori in Gyr.'

  OPEN(150, file='errsu.dat')

  i=0

      DO
        READ(150,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

      ALLOCATE(errkTsu(ndati), kT(ndati), errnsu(ndati), n(ndati), tcool(ndati), T(ndati), errTsu(ndati))
      ALLOCATE(relTsu(ndati), relnsu(ndati), relcoolsu(ndati), errcoolsu(ndati))

      REWIND(150)

      OPEN(5, file='Results_err.dat')
      OPEN(6, file='Printerr.dat')

      DO i=1,ndati
         READ(150,*) errkTsu(i), kT(i), errnsu(i), n(i), tcool(i)

         T(i)=(kT(i)*1e+3*1.6e-12)/(1.38e-16)
         errTsu(i)=(errkTsu(i)*1e+3*1.6e-12)/(1.38e-16)
         relTsu(i)=errTsu(i)/T(i)
         relnsu(i)=errnsu(i)/n(i)
         relcoolsu(i)=relTsu(i)+relnsu(i)
         errcoolsu(i)=relcoolsu(i)*tcool(i)

      WRITE(5,*) "L'errore superiore sul tempo di cooling della region",i, 'è:'
      WRITE(5,*) errcoolsu(i)/1e+9, 'Gyr'
      WRITE(6,*) errcoolsu(i)/1e+9

      END DO
      
         
OPEN(160, file='errgiu.dat')

  !j=0

     ! DO
      !  READ(160,*,end=100)
      !  j=j+1
      !END DO

      !ndati=j

      ALLOCATE(errkTgiu(ndati),errngiu(ndati), errTgiu(ndati))
      ALLOCATE(relTgiu(ndati), relngiu(ndati), relcoolgiu(ndati), errcoolgiu(ndati))

!REWIND(160)

      i=0

      DO i=1,ndati
         READ(160,*) errkTgiu(i), T(i), errngiu(i), n(i), tcool(i)

         T(i)=(kT(i)*1e+3*1.6e-12)/(1.38e-16)
         errTgiu(i)=(errkTgiu(i)*1e+3*1.6e-12)/(1.38e-16)
         relTgiu(i)=errTgiu(i)/T(i)
         relngiu(i)=errngiu(i)/n(i)
         relcoolgiu(i)=relTgiu(i)+relngiu(i)
         errcoolgiu(i)=relcoolgiu(i)*tcool(i)

    WRITE(5,*) "L'errore inferiore sul tempo di cooling della region",i, "è:"
    WRITE(5,*) errcoolgiu(i)/1e+9, 'Gyr'
    WRITE(6,*) errcoolgiu(i)/1e+9

   END DO
   
      

END PROGRAM errors
