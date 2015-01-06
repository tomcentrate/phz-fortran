! Program to perform a random walk in two dimensions
! repeats for many realizations and accmulates r**1 statistics
! Also compute entropy vs. time
	IMPLICIT NONE
	INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
	INTEGER :: i,j,is,ix,iy,it
	INTEGER :: inc !! How much to increment the walk by.
	INTEGER :: ne,nt
	REAL(Kind=Prec14) :: rnd,rnds,r2,entropy
	INTEGER, PARAMETER :: nx=200, ny=200, ixmin=-99, ixmax=100, iymin=-99, iymax=100 ! Lattice s
	INTEGER, PARAMETER :: nxg=8, nyg=8 ! grid sites
	INTEGER, PARAMETER :: ntmax=1000 ! maximum number of time steps
	INTEGER, PARAMETER :: nemax=100000 ! number of walks in ensemble
	INTEGER, PARAMETER :: ient=10 ! Number of steps to update p array
	INTEGER, PARAMETER :: ipmax=2000 ! Maximum number of time elements in p
	REAL(Kind=Prec14), DIMENSION(ntmax) ::r2a(ntmax) ! accumulated average of r**2
	REAL(Kind=Prec14), DIMENSION(nxg,nyg,ipmax) :: p ! Probability of visiting a site
	REAL(Kind=Prec14) :: t
	! set r2a=0, p=0
	r2a=0.0d0
	p=0.0d0
	! initialize random number generator
	call RANDOM_SEED
	do ne=1,nemax ! Nemax realizations of random walk
		i = 0
		j = 0
		do nt=1,ntmax
		call RANDOM_NUMBER(rnd)
	        call RANDOM_NUMBER(rnds)

			inc = 0
		! Calculates whether they go forwards or backwards
		if (rnd > 0.5d0) then 
			inc = 1
			else
			 inc = -1
		endif ! end calc.
		! Calculates if its either x or y that goes.
		if (rnds > 0.5d0) then
			 i = i + inc
			else 
			j = j + inc
		endif
		! Prevent oob condition for x
		!if (i < ixmin) i = 1
		!if (i > ixmax) i = ixmax
		! Prevent oob condition for y
		!if (j < iymin) j = 1
		!if (j > iymax) j = iymax

		!if(mod(nt,ient).eq.0) then ! update array
		!	it=nt/ient+1
		!	ix=int(dfloat(i+99)*dfloat(nxg)/dfloat(nx))+1
		!	iy=int(dfloat(j+99)*dfloat(nyg)/dfloat(ny))+1
		!	p(ix,iy,it)=p(ix,iy,it)+1.0d0
		! endif
		r2=real(i)**2+real(j)**2
		r2a(nt)=r2a(nt)+r2
		enddo ! Loop over step
	enddo ! Loop over random walks
	
	do nt=1, ntmax
		print *,nt,  r2a(nt)/dfloat(nemax)
	enddo
	end
