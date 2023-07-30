!*******************************************************************************
!------------------------------------------------------------------------------
!  write this part out individually by each proc
!  Chunlei Liang  02/17/2011
!-------------------------------------------------------------------------------
	SUBROUTINE tecplotter3d(Np)
	use setup3d	
	IMPLICIT NONE
	INCLUDE 'mpif.h'
!
	integer, intent(in) :: Np
	integer	:: icell, ii,jj,kk, gnid, is,js,ks, i,j, nd,nd1,nd2,iedge,idir
	double precision :: gam1, xsi, eta, bet,v2,dist
	double precision :: Qs(numv),Qp(numv)

	allocate(plotQ(numv,plotnodes))
	! rho u v w p

	plotQ(:,:) = 0.d0
	gam1 = gam - 1.d0
	
	do icell=1,NCELL
	
		do kk=1,Np+1
		do jj=1,Np+1
		do ii=1,Np+1
			xsi = dble(ii-1)/dble(Np)
			eta = dble(jj-1)/dble(Np)
			bet = dble(kk-1)/dble(Np)
			gnid = gnumnode(icell, kk, jj, ii)
			Qs = 0.d0
	
			do ks=1,N
			do js=1,N
			do is=1,N
				Qs(1:numv) = Qs(1:numv) + Q(1:numv,is,js,ks,icell) &
				* hval(N,is,xsi) * hval(N,js,eta) * hval(N,ks,bet)
			end do
			end do
			end do
	
			plotQ(1:numv,gnid) = plotQ(1:numv,gnid) + Qs(1:numv)

		end do
		end do
		end do
		
	end do
	
	
	do nd=1,plotnodes
		plotQ(1:numv,nd) = plotQ(1:numv,nd)/dble(gnode_factor(nd))
		
		plotQ(2,nd) = plotQ(2,nd)/plotQ(1,nd)				! U-velocity
        plotQ(3,nd) = plotQ(3,nd)/plotQ(1,nd)				! V-velocity
        plotQ(4,nd) = plotQ(4,nd)/plotQ(1,nd)				! W-velocity

        v2= plotQ(2,nd)**2+plotQ(3,nd)**2+plotQ(4,nd)**2
        plotQ(5,nd) = gam1*(plotQ(5,nd)-0.5*plotQ(1,nd)*v2)	! pressure	
    end do  


   	contains		
			
	!***************************************
	DOUBLE PRECISION FUNCTION hval(np,i,xval)
	
	use	setup3d
	IMPLICIT NONE
	
	integer, intent(in)	:: np, i
	double precision, intent(in)	:: xval
	double precision	:: hvaln, hvald
	integer	:: s
		
	hvaln = 1.d0
	hvald = 1.d0
	do s=1,np
		if( s/=i) then
			hvaln = hvaln * ( xval - Xs(s) ) 
			hvald = hvald * ( Xs(i) - Xs(s) )
		end if		
	end do
	
	hval = hvaln/hvald

!	write(*,*) i, xval, hval	
			
	END FUNCTION hval		
	
	END SUBROUTINE tecplotter3d	
    	
	