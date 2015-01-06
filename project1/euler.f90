	IMPLICIT NONE
	INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
	REAL(Kind=Prec14), PARAMETER :: dt=0.05d0 ! timestep
	REAL(Kind=Prec14), PARAMETER :: time=50.0d0 ! total simulation time
	REAL(Kind=Prec14), PARAMETER :: y0=1.0d0 ! initial displacement
	REAL(Kind=Prec14), PARAMETER :: w0=1.0d0 ! initial ang velocity
	REAL(Kind=Prec14) :: y,yn ! current position y, next position yn
	REAL(Kind=Prec14) :: w,wn ! current ang vel w, next ang vel wn
	REAL(Kind=Prec14) :: energy ! computed energy
	REAL(Kind=Prec14) :: energy0 ! computed energy at t=0
	REAL(Kind=Prec14) :: t ! current time
	REAL(Kind=Prec14) :: A,B ! coefficients for analytical solution, from initial conditions
	REAL(Kind=Prec14) :: yexact,wexact ! position, angular velocity from analytical solution
	INTEGER :: nsteps,n
	nsteps=anint(time/dt) ! number of integration steps, integer value
	y=y0 ; w=w0 ! set initial position, angular velocity
! Determine A,B assuming y(t)=Acos(omega t) + Bsin(omega t)
	A=y
	B=w
! Determine the initial energy
	energy0=0.5d0*(w**2+y**2)
	t=0.0d0

	do n=1,nsteps
		t=t+dt ! step forward time t
		! Integrate step using Euler
		wn=w-y*dt ! step forward angular velocity
		yn=y+w*dt ! step forward position
		w=wn
		y=yn
		! Compute current energy
		energy=0.5d0*(w**2+y**2)
		! Find exact result
		yexact=A*dcos(t)+B*dsin(t) ! Double precision sine, cosine
		wexact=-A*dsin(t)+B*dcos(t)
		! Output position, energy including analytical
		write(*,100) t,y,yexact,w,wexact,energy,energy0
	enddo
100 	 format(f8.4,6f12.6) ! Format statement
	stop
	end
