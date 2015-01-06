	subroutine four1(data,nn,isign)
	integer isign,nn
	double precision data(2*nn)
!	 replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or
!	 replaces data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign
!	 is input as -1. data must bbe a complex array of length nn or, equivalently, a real
!	 array of legnth 2*nn. nn must be an integer power of 2 (this is not checked for!)
	integer i,istep,j,m,mmax,n
	double precision tempi,tempr
	double precision theta,wi,wpi,wpr,wr,wtemp

	pi=4.0d0*datan(1.0d0)
	n=2*nn
	j=1
	do i=1,n,2
	 if(j.gt.i) then
	  tempr=data(j)
	  tempi=data(j+1)
	  data(j)=data(i)
	  data(j+1)=data(i+1)
	  data(i)=tempr
	  data(i+1)=tempi
	 endif
         m=nn
1	 if((m.ge.2).and.(j.gt.m)) then
	  j=j-m
	  m=m/2
	  goto 1
	 endif
	 j=j+m
	enddo

	mmax=2
2	if(n.gt.mmax) then
	 istep=2*mmax
	 theta=2.0d0*pi/(isign*mmax)
	 wpr=-2.0d0*dsin(0.5d0*theta)**2
	 wpi=dsin(theta)
	 wr=1.0d0
	 wi=0.0d0
	 do m=1,mmax,2
	 do i=m,n,istep
	  j=i+mmax
	  tempr=wr*data(j)-wi*data(j+1)
	  tempi=wr*data(j+1)+wi*data(j)
	  data(j)=data(i)-tempr
	  data(j+1)=data(i+1)-tempi
	  data(i)=data(i)+tempr
	  data(i+1)=data(i+1)+tempi
	 enddo
	wtemp=wr
	wr=wr*wpr-wi*wpi+wr
	wi=wi*wpr+wtemp*wpi+wi
	enddo
	mmax=istep
	goto 2
	endif
	return
	end
