	IMPLICIT NONE
	INTEGER,PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
	REAL(Kind=Prec14), PARAMETER :: dt=0.5d0 ! timestep
	REAL(Kind=Prec14), PARAMETER :: time=50.0d0 ! total simulation time
	REAL(Kind=Prec14), PARAMETER :: y0=1.0d0 ! initial displacement
	REAL(Kind=Prec14), PARAMETER :: w0=1.0d0 ! initial angular velocity
	REAL(Kind=Prec14) :: y,yn ! current position y, next position y
	REAL(Kind=Prec14) :: w, wn ! current angular vel w, next ang vel w
	REAL(Kind=Prec14) :: energy ! computed energy
	REAL(Kind=Prec14) :: energy0 ! computed energy at t = 0
	REAL(Kind=Prec14) :: t ! current time
	REAL(Kind=Prec14) :: A,B ! coefficients for analytical solutions.
	REAL(Kind=Prec14) :: yexact, wexact
	INTEGER :: nsteps,n
	nsteps= anint(time/dt) !
	print *, "Printing Stuff Here", nsteps
	y = y0
	w = w0
!Determine A,B assuming y(t) = Acos(omega t) + B sin(omega t )
	A = y
	B = w 
!Determine initial energy
	energy0 = 0.5d0*(w**2 + y**2)
	t = 0.0d0
	do n=1, nsteps
		t = t+dt
		wn= w-y*dt
		yn= y+w*dt
		w = wn
		y = yn
		energy = 0.5d0*(w**2 + y**2)
		yexact = A*dcos(t) + B*dsin(t)
		wexact = A*dsin(t) + B*dcos(t)

		print *, t,y,yexact,w,wexact,energy,energy0
	enddo

100	format(f8.4,6f12.6)
	stop
	end

