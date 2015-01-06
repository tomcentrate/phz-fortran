		IMPLICIT NONE
		INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
		INTEGER :: n,nsteps
		INTEGER, PARAMETER :: NSTEPS=10000
		REAL(Kind=Prec14) :: vxj,vyj,vxa,vya ! velocities of jupiter, asteroid
		REAL(Kind=Prec14) :: vxji,vyji,vxai,vyai ! initial velocities of jupiter, asteroid
		REAL(Kind=Prec14) :: xj,yj,xa,ya ! positions of jupiter, asteroid
		REAL(Kind=Prec14) :: xja,yja ! x,y separation between asteroid and jupiter
		REAL(Kind=Prec14) :: xji,yji,xai,yai ! initial positions of jupiter, asteroid
		REAL(Kind=Prec14) :: alpha,alpha_j,alpha_a,dt,time
		REAL(Kind=Prec14) :: rj,ra,rja ! separation between jupiter-sun, asteroid-sun, and jupiter-asteroid
		REAL(Kind=Prec14) :: pi
		REAL(Kind=Prec14) :: ms,mj,ma,ekin,epot,etot,amom
		REAL(Kind=Prec14), PARAMETER :: dt=0.05d0 ! time step in units of years
		REAL(Kind=Prec14), PARAMETER :: ms=1.0d0,mj=9.55d-4,ma=4.39d-10 ! asteroid mass is that of 1 Cere
		REAL(Kind=Prec14), PARAMETER :: xji=0.0d0,yji=5.200d0 ! initial position of jupiter in AU
		REAL(Kind=Prec14), PARAMETER :: vxji=2.75535903d0,vyji=0.0d0 ! initial velocity of jupiter in AU/yr
		! REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.000d0 ! initial position of asteroid in AU
		! REAL(Kind=Prec14), PARAMETER :: vxai=3.628d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr
		REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.276d0 ! initial position of asteroid in AU
		REAL(Kind=Prec14), PARAMETER :: vxai=3.471d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr
		! REAL(Kind=Prec14), PARAMETER :: xai=0.0d0,yai=3.700d0 ! initial position of asteroid in AU
		! REAL(Kind=Prec14), PARAMETER :: vxai=3.267d0,vyai=0.0d0 ! initial velocity of asteroid in AU/yr
		open(unit=10,file='positions.out')
		open(unit=11,file='conserved.out')

		do n=1,nsteps

			ra = sqr(xa**2 + ya**2)
			rj = sqr(xj**2 + yj**2)
			rja = sqrt((xa-xj)**2 + (ya-yj)**2 )
		
			vxa = vxa - 4 * pi**2 * xa * dt /ra**3 - 4 * pi**2 * (mj/ms)  * (xa - xj) * dt / rja**3
			vya = vya - 4 * pi**2 * ya * dt /ra**3 - 4 * pi**2 * (mj/ms)  * (ya - yj) * dt / rja**3
			vxj = vxj - 4 * pi**2 * xj * dt /rj**3 - 4 * pi**2 * (ma/ms)  * (xj - xa) * dt / rja**3
			vyj = vyj - 4 * pi**2 * yj * dt /rj**3 - 4 * pi**2 * (ma/ms)  * (yj - ya) * dt / rja**3

			amom=mj*(xj*vyj-yj*vxj)+ma*(xa*vya-ya*vxa)
			ekin=0.5d0*mj*(vxj**2+vyj**2)+0.5d0*ma*(vxa**2+vya**2)
			epot=-4.0d0*pi**2*(mj/rj+ma/ra+(ma*mj/ms)/rja)
			etot=ekin+epot
			write(10,100) time,xj,yj,xa,ya,rj,ra
			write(11,200) time,amom,ekin,epot,etot

		enddo
100 	format(7f12.6)
200 	format(5f12.6)