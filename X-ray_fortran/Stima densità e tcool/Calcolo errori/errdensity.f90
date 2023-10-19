PROGRAM errden
  IMPLICIT NONE
  REAL*8:: r_min, r_max, V, errradius, errvol, norm, errnorm_su, errnorm_giu, pi
  REAL*8:: relvol, relnormsu, relnormgiu, reltotsu, reltotgiu, absnsu, absngiu

pi=3.1415
  
  PRINT*, 'Raggio min e max in arcsec?'
  READ(*,*) r_min, r_max

  V=(((4*pi)/3)*(r_max**3))-(((4*pi)/3)*(r_min**3))

  PRINT*, 'Errore sul raggio (arcsec)?'
  READ(*,*) errradius

  errvol=2*(errradius**3)

  PRINT*, 'Valore di norm?'
  READ(*,*) norm

  PRINT*, 'Errore superiore e inferiore su norm?'
  READ(*,*) errnorm_su, errnorm_giu

  relvol=errvol/V
  relnormsu=errnorm_su/norm
  relnormgiu=errnorm_giu/norm

  reltotsu=relnormsu+relvol
  reltotgiu=relnormgiu+relvol

  absnsu=reltotsu*norm
  absngiu=reltotgiu*norm

  PRINT*, 'Errore superiore di n:', absnsu
  PRINT*, 'Errore inferiore di n:', absngiu
  



END PROGRAM errden


