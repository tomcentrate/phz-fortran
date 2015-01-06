IMPLICIT NONE
INTEGER, PARAMETER::Prec14 = SELECTED_REAL_KIND(14)
INTEGER :: nsteps,n
REAL(Kind=Prec14) :: vxj,vyj,vxa,vya ! velocities of jupiter, asteroid
REAL(Kind=Prec14) :: xj,yj,xa,ya ! positions of jupiter, asteroid
REAL(Kind=Prec14) :: xja,yja ! x,y separation between asteroid and jupiter
REAL(Kind=Prec14) :: time
REAL(Kind=Prec14) :: rj,ra,rja ! separation between jupiter-sun, asteroid-sun, and jupiter-asteroid
REAL(Kind=Prec14) :: ekin,epot,etot,amom

REAL(Kind=Prec14), PARAMETER :: PI=4.D0*DATAN(1.D0)
REAL(Kind=Prec14), PARAMETER :: dt=0.005d0 ! time step in units of years
REAL(Kind=Prec14), PARAMETER :: ms=1.0d0,mj=9.55d-4,ma=4.39d-10 ! asteroid mass is that of 1 Ceres
REAL(Kind=Prec14), PARAMETER :: xji=0.0d0,yji=5.200d0 ! initial position of jupiter in AU
REAL(Kind=Prec14), PARAMETER :: vxji=2.75535903d0,vyji=0.0d0 ! initial velocity of jupiter in AU/yr

REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.276d0 ! initial position of asteroid in AU
REAL(Kind=Prec14), PARAMETER :: vxai=3.471d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr

!REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.000d0 ! initial position of asteroid in AU
!REAL(Kind=Prec14), PARAMETER :: vxai=3.628d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr

!REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.700d0 ! initial position of asteroid in AU
!REAL(Kind=Prec14), PARAMETER :: vxai=3.267d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr

nsteps = 10000
! set initial positions and velocities
xj=xji
yj=yji
xa=xai
ya=yai
vxj=vxji
vyj=vyji
vxa=vxai
vya=vyai

time = 0.0d0

open(unit=10,file='positions.out')
open(unit=11,file='conserved.out')

do n=1,nsteps 
	time=time+dt
	!Calculate the distance among Jupiter, Asteroid, and the Sun
	ra = SQRT(xa*xa + ya*ya)
	rj = SQRT(xj*xj + yj*yj)
	rja = SQRT((xa-xj)*(xa-xj) + (ya-yj)*(ya-yj))
	
	!New jupiter velocities
	vxj = vxj - 4*PI*PI*xj*dt / (rj*rj*rj) - 4*PI*PI*(mj)*(xj-xa)*dt / (rja*rja*rja)
	vyj = vyj - 4*PI*PI*yj*dt / (rj*rj*rj) - 4*PI*PI*(mj)*(yj-ya)*dt / (rja*rja*rja)
	
	!New Asteroid velocities
	vxa = vxa - 4*PI*PI*xa*dt / (ra*ra*ra) - 4*PI*PI*(ma)*(xa-xj)*dt / (rja*rja*rja)
	vya = vya - 4*PI*PI*ya*dt / (ra*ra*ra) - 4*PI*PI*(ma)*(ya-yj)*dt / (rja*rja*rja)
	
	!Euler-Cromer for new Positions
	
	xj = xj + vxj*dt
	yj = yj + vyj*dt
	
	xa = xa + vxa*dt
	ya = ya + vya*dt	
	
	!if(mod(n,10)==0) then
	write(10,100) time,xj,yj,xa,ya,rj,ra
	write(11,200) time,amom,ekin,epot,etot	
    !end if
enddo

100 format(7f12.6)
200 format(5f12.6)

stop
end