    SUBROUTINE MAPDER

	use setup3d
	implicit none
	include 'mpif.h'

	integer :: is,js,ks,ifp,jfp,kfp,ii,jj
	double precision :: xi,eta,beta
	double precision,allocatable :: mdfunc(:,:)

	!COMPUTE THE DERIVATIVES AT THE SOLUTION POINTS
	!ALL QUADRATIC ELEMENTS

	if (CURVE_WALL==1) then
	    allocate(mdfunc(3,20))
	else
		allocate(mdfunc(3,8))
	end if

	do ks=1,N
	do js=1,N
	do is=1,N

	    xi   = Xs(is)
		eta  = Xs(js)
		beta = Xs(ks)

		if (CURVE_WALL==1) then
			call DERBASEFUNC20(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,20
				dmdxs(jj,ii,is,js,ks) = mdfunc(ii,jj)   
			end do 
			end do
		else
			call DERBASEFUNC(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,8
				dmdxs(jj,ii,is,js,ks) = mdfunc(ii,jj)
			end do
			end do
		end if

	end do
	end do
	end do


	do kfp=1,N
	do jfp=1,N
	do ifp=1,N+1
		
		xi   = Xf(ifp)
		eta  = Xs(jfp)
		beta = Xs(kfp)

		if (CURVE_WALL==1) then
			call DERBASEFUNC20(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,20
				dmdxf1(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)   
			end do 
			end do
		else
			call DERBASEFUNC(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,8
				dmdxf1(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)
			end do
			end do
		end if

	end do
	end do
	end do


	do kfp=1,N
	do jfp=1,N+1
	do ifp=1,N

		xi   = Xs(ifp)
		eta  = Xf(jfp)
		beta = Xs(kfp)
			
		if (CURVE_WALL==1) then
			call DERBASEFUNC20(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,20
				dmdxf2(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)   
			end do 
			end do
		else
			call DERBASEFUNC(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,8
				dmdxf2(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)
			end do
			end do
		end if

	end do
	end do
	end do


	do kfp=1,N+1
	do jfp=1,N
	do ifp=1,N
		
		xi   = Xs(ifp)
		eta  = Xs(jfp)
		beta = Xf(kfp)
			
		if (CURVE_WALL==1) then
			call DERBASEFUNC20(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,20
				dmdxf3(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)   
			end do 
			end do
		else
			call DERBASEFUNC(mdfunc,xi,eta,beta)
			do ii=1,3
			do jj=1,8
				dmdxf3(jj,ii,ifp,jfp,kfp) = mdfunc(ii,jj)
			end do
			end do
		end if

	end do
	end do
	end do
	

    END SUBROUTINE MAPDER

!###############################################################################
 
SUBROUTINE MAPBASEFUNC(func,xi,eta,beta)
!MAPPING FUNCTION FOR 3D ELEMENTS.

	implicit none

	double precision :: func(8),xi,eta,beta

	func(1) = (1-xi) * (1-eta) * (1-beta)
	func(2) = xi * (1-eta) * (1-beta)
	func(3) = xi * eta * (1-beta)
	func(4) = (1-xi) * eta * (1-beta)
	func(5) = (1-xi) * (1-eta) * beta
	func(6) = xi * (1-eta) * beta
	func(7) = xi * eta * beta
	func(8) = (1-xi) * eta * beta

	return
	
END SUBROUTINE MAPBASEFUNC
	
!###############################################################################
	
    SUBROUTINE DERBASEFUNC(dfunc,xi,eta,beta)
    !DERIVATIVE OF MAPPING FUNCTION FOR 3D ELEMENTS.

	!dfunc(1,:)   	derivative wrt xi
	!dfunc(2,:)   	derivative wrt eta
	!dfunc(3,:)   	derivative wrt beta
    
	implicit none

	double precision:: dfunc(3,8),xi,eta,beta

	dfunc(1,1) = (eta-1)*(1-beta)
	dfunc(2,1) = (xi-1)*(1-beta) 
	dfunc(3,1) = (xi-1)*(1-eta)
	
	dfunc(1,2) = (1-eta)*(1-beta)
	dfunc(2,2) = xi*(beta-1) 
	dfunc(3,2) = xi*(eta-1) 
			
	dfunc(1,3) = eta*(1-beta)
	dfunc(2,3) = xi*(1-beta) 
	dfunc(3,3) = -xi*eta 
			
	dfunc(1,4) = eta*(beta-1)
	dfunc(2,4) = (1-xi)*(1-beta)
	dfunc(3,4) = (xi-1)*eta

	dfunc(1,5) = (eta-1)*beta
	dfunc(2,5) = (xi-1)*beta
	dfunc(3,5) = (1-xi)*(1-eta) 

	dfunc(1,6) = (1-eta)*beta
	dfunc(2,6) = -xi*beta
	dfunc(3,6) = xi*(1-eta)

	dfunc(1,7) = eta*beta
	dfunc(2,7) = xi*beta
	dfunc(3,7) = xi*eta

	dfunc(1,8) = -eta*beta
	dfunc(2,8) = (1-xi)*beta
	dfunc(3,8) = (1-xi)*eta

	return
	
    END SUBROUTINE DERBASEFUNC
	
!###############################################################################

SUBROUTINE MAPBASEFUNC20(func,xi,eta,beta)
!MAPPING FUNCTION FOR 3D QUADRATIC ELEMENTS.

	implicit none

	double precision:: func(20),xi,eta,beta

	func(1) = (1-xi)*(1-eta)*(1-beta)*(-2*xi-2*eta-2*beta+1)
	func(4) = (1-xi)*eta*(1-beta)*(-2*xi+2*eta-2*beta-1)
	func(8) = (1-xi)*eta*beta*(-2*xi+2*eta+2*beta-3)
	func(5) = (1-xi)*(1-eta)*beta*(-2*xi-2*eta+2*beta-1)
	func(2) = xi*(1-eta)*(1-beta)*(2*xi-2*eta-2*beta-1)
	func(3) = xi*eta*(1-beta)*(2*xi+2*eta-2*beta-3)
	func(7) = xi*eta*beta*(2*xi+2*eta+2*beta-5)
	func(6) = xi*(1-eta)*beta*(2*xi-2*eta+2*beta-3)

	func(9)  = 4*xi*(1-xi)*(1-eta)*(1-beta)
	func(11) = 4*xi*(1-xi)*eta*(1-beta)
	func(15) = 4*xi*(1-xi)*eta*beta
	func(13) = 4*xi*(1-xi)*(1-eta)*beta
	func(12) = 4*(1-xi)*eta*(1-eta)*(1-beta)
	func(16) = 4*(1-xi)*eta*(1-eta)*beta
	func(14) = 4*xi*eta*(1-eta)*beta
	func(10) = 4*xi*eta*(1-eta)*(1-beta)

	func(18) = 4*(1-xi)*(1-eta)*beta*(1-beta)
	func(17) = 4*xi*(1-eta)*beta*(1-beta)
	func(19) = 4*xi*eta*beta*(1-beta)
	func(20) = 4*(1-xi)*eta*beta*(1-beta)

	return
	
END SUBROUTINE MAPBASEFUNC20

!###############################################################################

SUBROUTINE DERBASEFUNC20(dfunc,xi,eta,beta)
!DERIVATIVE OF MAPPING FUNCTION FOR 3D QUADRATIC ELEMENTS.

	!dfunc(1,:)   	derivative wrt xi
	!dfunc(2,:)   	derivative wrt eta
	!dfunc(3,:)   	derivative wrt beta
	!dfunc(4,:)		derivative wrt tau

	implicit none

	double precision :: dfunc(3,20),xi,eta,beta

	dfunc(1,1) = -2*(1-eta)*(1-xi)*(1-beta)-(1-eta)*(1-2*eta-2*xi-2*beta)*(1-beta)
	dfunc(2,1) = -2*(1-eta)*(1-xi)*(1-beta)-(1-xi)*(1-2*eta-2*xi-2*beta)*(1-beta)
	dfunc(3,1) = -(1-eta)*(1-xi)*(1-2*eta-2*xi-2*beta)-2*(1-eta)*(1-xi)*(1-beta)
	dfunc(1,4) = -2*eta*(1-xi)*(1-beta)- eta*(-1+2*eta-2*xi-2*beta)*(1-beta)
	dfunc(2,4) = 2*eta*(1-xi)*(1-beta)+(1-xi)*(-1+2*eta-2*xi-2*beta)*(1-beta)
	dfunc(3,4) = -eta*(1-xi)*(-1+2*eta-2*xi-2*beta)-2*eta*(1-xi)*(1-beta)
	dfunc(1,8) = -2*eta*(1-xi)*beta-eta*beta*(-3+2*eta-2*xi+2*beta)
	dfunc(2,8) = 2*eta*(1-xi)*beta+(1-xi)*beta*(-3+2*eta-2*xi+2*beta)
	dfunc(3,8) = 2*eta*(1-xi)*beta+eta*(1-xi)*(-3+2*eta-2*xi+2*beta)
	dfunc(1,5) = -2*(1-eta)*(1-xi)*beta-(1-eta)*beta*(-1-2*eta-2*xi+2*beta)
	dfunc(2,5) = -2*(1-eta)*(1-xi)*beta-(1-xi)*beta*(-1-2*eta-2*xi+2*beta)
	dfunc(3,5) = 2*(1-eta)*(1-xi)*beta+(1-eta)*(1-xi)*(-1-2*eta-2*xi+2*beta)
	dfunc(1,2) = 2*(1-eta)*xi*(1-beta)+(1-eta)*(-1-2*eta+2*xi-2*beta)*(1-beta)
	dfunc(2,2) = -2*(1-eta)*xi*(1-beta)-xi*(-1-2*eta+2*xi-2*beta)*(1-beta)
	dfunc(3,2) = -(1-eta)*xi*(-1-2*eta+2*xi-2*beta)-2*(1-eta)*xi*(1-beta)
	dfunc(1,3) = 2*eta*xi*(1-beta)+eta*(-3+2*eta+2*xi-2*beta)*(1-beta)
	dfunc(2,3) = 2*eta*xi*(1-beta)+xi*(-3+2*eta+2*xi-2*beta)*(1-beta)
	dfunc(3,3) = -eta*xi*(-3+2*eta+2*xi-2*beta)-2*eta*xi*(1-beta)
	dfunc(1,7) = 2*eta*xi*beta+eta*beta*(-5+2*eta+2*xi+2*beta)
	dfunc(2,7) = 2*eta*xi*beta+xi*beta*(-5+2*eta+2*xi+2*beta)
	dfunc(3,7) = 2*eta*xi*beta+eta*xi*(-5+2*eta+2*xi+2*beta)
	dfunc(1,6) = 2*(1-eta)*xi*beta+(1-eta)*beta*(-3-2*eta+2*xi+2*beta)
	dfunc(2,6) = -2*(1-eta)*xi*beta-xi*beta*(-3-2*eta+2*xi+2*beta)
	dfunc(3,6) = 2*(1-eta)*xi*beta+(1-eta)*xi*(-3-2*eta+2*xi+2*beta)

	dfunc(1,9) = 4*(1-eta)*(1-xi)*(1-beta)-4*(1-eta)*xi*(1-beta)
	dfunc(2,9) = -4*(1-xi)*xi*(1-beta)
	dfunc(3,9) = -4*(1-eta)*(1-xi)*xi

	dfunc(1,11) = 4*eta*(1-xi)*(1-beta)-4*eta*xi*(1-beta)
	dfunc(2,11) = 4*(1-xi)*xi*(1-beta)
	dfunc(3,11) = -4*eta*(1-xi)*xi
	dfunc(1,15) = 4*eta*(1-xi)*beta-4*eta*xi*beta
	dfunc(2,15) = 4*(1-xi)*xi*beta
	dfunc(3,15) = 4*eta*(1-xi)*xi
	dfunc(1,13) = 4*(1-eta)*(1-xi)*beta-4*(1-eta)*xi*beta
	dfunc(2,13) = -4*(1-xi)*xi*beta
	dfunc(3,13) = 4*(1-eta)*(1-xi)*xi
	dfunc(1,12) = -4*(1-eta)*eta*(1-beta)
	dfunc(2,12) = 4*(1-eta)*(1-xi)*(1-beta)-4*eta*(1-xi)*(1-beta)
	dfunc(3,12) = -4*(1-eta)*eta*(1-xi)
	dfunc(1,16) = -4*(1-eta)*eta*beta
	dfunc(2,16) = 4*(1-eta)*(1-xi)*beta-4*eta*(1-xi)*beta
	dfunc(3,16) = 4*(1-eta)*eta*(1-xi)
	dfunc(1,14) = 4*(1-eta)*eta*beta
	dfunc(2,14) = 4*(1-eta)*xi*beta-4*eta*xi*beta
	dfunc(3,14) = 4*(1-eta)*eta*xi
	dfunc(1,10) = 4*(1-eta)*eta*(1-beta)
	dfunc(2,10) = 4*(1-eta)*xi*(1-beta)-4*eta*xi*(1-beta)
	dfunc(3,10) = -4*(1-eta)*eta*xi

	dfunc(1,18) = -4*(1-eta)*(1-beta)*beta
	dfunc(2,18) = -4*(1-xi)*(1-beta)*beta
	dfunc(3,18) = 4*(1-eta)*(1-xi)*(1-beta)-4*(1-eta)*(1-xi)*beta
	dfunc(1,17) = 4*(1-eta)*(1-beta)*beta
	dfunc(2,17) = -4*xi*(1-beta)*beta
	dfunc(3,17) = 4*(1-eta)*xi*(1-beta)-4*(1-eta)*xi*beta
	dfunc(1,19) = 4*eta*(1-beta)*beta
	dfunc(2,19) = 4*xi*(1-beta)*beta
	dfunc(3,19) = 4*eta*xi*(1-beta)-4*eta*xi*beta
	dfunc(1,20) = -4*eta*(1-beta)*beta
	dfunc(2,20) = 4*(1-xi)*(1-beta)*beta
	dfunc(3,20) = 4*eta*(1-xi)*(1-beta)-4*eta*(1-xi)*beta

	return
	
END SUBROUTINE DERBASEFUNC20


SUBROUTINE xyCoor_atGps(ns,xxs,xi,eta,beta,xx)

	implicit none
	
	integer, intent(in) :: ns
	integer :: i,nb
	double precision :: xxs(3,ns),xx(3),xi,eta,beta
	double precision, pointer, dimension(:) :: func
	
	allocate(func(ns))
	
	if (ns==8) call MAPBASEFUNC(func,xi,eta,beta)
	
	if (ns==20) call MAPBASEFUNC20(func,xi,eta,beta)
	
	do i=1,3
		xx(i) = 0.d0
		do nb = 1,ns
			xx(i) = xx(i) + xxs(i,nb)*func(nb)
		end do
	end do
	
	deallocate(func)
	
	return
	
END SUBROUTINE xyCoor_atGps