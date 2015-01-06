	IMPLICIT NONE
	INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
	INTEGER, PARAMETER :: nsteps=21
	INTEGER :: n
	REAL(Kind=Prec14),DIMENSION(nsteps) :: X,Y
	REAL(Kind=Prec14) :: s,h,t
	OPEN(1, FILE='data/results.dat')
	t=0
	h=0
	s=0
	do n=1,nsteps ! read in all the points
		READ(1,*) X(n),Y(n)
	enddo
! Trapezoidal Rule
	h=((X(nsteps) - X(1))/((nsteps-1))) ! end multiplier
	t = Y(1)+Y(nsteps) ! adds the first and last points
	do n=2,nsteps-1
		t=t+(2*Y(n)) ! adds the middle points
	enddo
	t = h*t/2
! End of trapezoidal rule
! Simpsons Rule
	do n=1,nsteps-2,2
		s = s+(Y(n)+(4*Y(n+1))+(Y(n+2)))
	enddo
	s =( s*h / 3)
! End Simpsons Rule
	write(6,*) "Simpsons Rule:    ", s
	write(6,*) "Trapezoidal Rule: ",t
	CLOSE(1)
	END
	

