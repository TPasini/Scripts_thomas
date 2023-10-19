PROGRAM parametri
IMPLICIT NONE 
REAL*8:: z, B, B_CMB, t1, t2, t3, t4, t5, trad, alfa, gamma, freq1, freq2
REAL*8:: dalfa, dB, j, k, delta1, delta2, delta3, dt,a , m, n

print*, 'inserire indice spettrale'
read(*,*) alfa
print*, 'inserire errore indice spettrale'

read(*,*) dalfa

B=1.8 !in microgauss
z=0.137
B_CMB=3.25*((1.+z)**2) !in microgauss
gamma=0.55 !indice di iniezione


freq1=6. !in gigahertz
freq2=1.7 !in gigahertz
dB=0.02 !errore campo magnetico


!scomponiamo la formula per semplicità

!t1=1610.*(B**0.5)
t1=1590.*(B**0.5)
t2=B**2+B_CMB**2
t3=((alfa-gamma)/(freq1-freq2))**0.5
t4=(log(freq1/freq2))**0.5
t5=(1.+z)**0.5


trad=t1*t3*t4/(t2*t5)


!delta1=1620.*((log(freq1/freq2))**0.5)/t5
!delta2=2.*t2
!k=((alfa-(gamma-1.)/2.)*(dB**2)*(B_CMB**2-3*B**2)/(B*((B_CMB**2+B**2)**2)))
!j=B*(dalfa**2)/((alfa-(gamma-1.)/2.)**4)
!delta3=sqrt(k**2+j**2)
!dt=delta1*delta3/delta2


!delta1=trad/t3
!delta2=0.5/t3
!dt=delta1*delta2*dalfa

dt=0.5*trad*dalfa/(alfa-gamma)

print*, 'eta radiativa', trad, '+-', dt, 'Myr'








END PROGRAM
