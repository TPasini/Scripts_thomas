PROGRAM density
  IMPLICIT NONE
  REAL*8:: z, pi, DL, DA, mu, mh, prop
  REAL*8, ALLOCATABLE:: ne(:), A(:), r_min(:), V(:), rho(:), r_max(:), rarcmin(:), rarcmax(:)
  INTEGER:: ans, ndati, i

70 PRINT*, 'Vuoi calcolare una serie di valori da un file o uno solo da terminale?'
  PRINT*,       '1. File (modificare DL e redshift a seconda del cluster)'
  PRINT*,       '2. Terminale'

  READ(*,*) ans

  IF(ans==1) THEN
     PRINT*, 'Ricordarsi di modificare distanza di Luminosità, redshift e scala arcsec/kpc.'
     PRINT*, 'I dati vengono letti dal file norm.dat. I raggi devono essere dati in arcsec.'
     PRINT*, 'I risultati sono visualizzati nel file Dens_results.dat.'
     PRINT*, 'Per stima di tempo di cooling, si utilizzi rho.dat, dove viene riportata la densità in g/cm^3, aggiungendo T.'

!Cambiare Distanza di Luminosità, redshift e proporzione a seconda del cluster
DL=273.7
z=0.06355
prop=1.173

  pi=3.1415
  DA=(DL*3e+24)/((1+z)**2)
  mu=0.61
  mh=1.67e-24
     
     OPEN(30, file='norm.dat')
     i=0

      DO
        READ(30,*,end=100)
        i=i+1
      END DO

      100 CONTINUE
      ndati=i

    ALLOCATE(ne(ndati), A(ndati), r_min(ndati), V(ndati), rho(ndati), r_max(ndati), rarcmin(ndati), rarcmax(ndati))

    REWIND(30)

    OPEN(40, file='Dens_results.dat')
    OPEN(120, file='rho.dat')

    WRITE(40,*) '     Min Rad:   ', ' Max Rad: ',' Density:      ', '   Numeric Density:'
    WRITE(40,*) '                                                      '
 DO i=1,ndati
    READ(30,*) rarcmin(i), rarcmax(i), A(i)
    r_max(i)=rarcmax(i)*prop*3e+21
    r_min(i)=rarcmin(i)*prop*3e+21
    V(i)=(((4*pi)/3)*(r_max(i)**3))-(((4*pi)/3)*(r_min(i)**3))

    ne(i)=sqrt(1e+14*(4*pi*A(i)*((DA*(1+z))**2))/(0.82*V(i)))
    rho(i)=ne(i)*mu*mh

 WRITE(40,'(F10.2,A,F10.2,A,E10.3,A,E10.3,A)') r_min(i)/(3e+21), 'kpc', r_max(i)/(3e+21), 'kpc', rho(i), 'g/cm^3', ne(i), 'cm^-3'
    WRITE(120,*) rho(i)
    
 END DO


ELSE IF(ans==2) THEN
   ALLOCATE(ne(1), A(1), r_min(1), V(1), r_max(1), rarcmin(1), rarcmax(1))
   PRINT*, 'Fattore di normalizzazione?'
   READ(*,*) A
   PRINT*, 'Raggio min in arcsec?'
   READ(*,*) rarcmin
   PRINT*, 'Raggio max in arcsec?'
   READ(*,*) rarcmax
   PRINT*, 'Redshift?'
   READ(*,*) z
   PRINT*, 'Distanza di luminosità in Mpc?'
   READ(*,*) DL
   PRINT*, 'A quanto equivale 1 arcsec in kpc?'
   READ(*,*) prop

  pi=3.1415
  DA=(DL*3e+24)/((1+z)**2)
  mu=0.61
  mh=1.67e-24

   r_max=rarcmax*prop*3e+21
   r_min=rarcmin*prop*3e+21
   V=(((4*pi)/3)*(r_max**3))-(((4*pi)/3)*(r_min**3))
   ne=sqrt(1e+14*(4*pi*A*((DA*(1+z))**2))/(0.82*V))
   PRINT*, 'Numeric Density:', ne, 'cm^-3'


ELSE
   PRINT*, 'Digita 1 o 2'
   GO TO 70
END IF

END PROGRAM density
