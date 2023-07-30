	SUBROUTINE tecplotter3dsetup(Np)
	
	use setup3d
	
	IMPLICIT NONE
	INCLUDE 'mpif.h'	
	
	integer, intent(in) :: Np
	integer :: icell,k,i,j,gnid,iedge,enode,ii,jj,kk,IC2V(8),kmax
	integer	:: checkmatch, existnode, mm,nn,exmm,exnn,iface	
	double precision :: xsi, eta, bet, rx(3), tolplot
	double precision,allocatable :: xxs(:,:)
	    
!   	write(*,*) "NEDGE : ", NEDGE
!    
    call alloc_memory    
    
    plotnodes = 0
    plotcells = 0
    tolplot = 1d-6*10d7
    gnode_factor(:) = 0
    isedgenod(:) = 0
	isfacenod(:) = 0
	
	if (CURVE_WALL==1) then
		allocate(xxs(3,20))
		kmax = 20
	else
		allocate(xxs(3,8))
		kmax = 8
	end if
    
    ! First set up the 8 vertices of each cell
    ! Basically creating new global nodes
    
    do icell=1,NCELL
    
    	do k=1,8
			IC2V(k) = IVCELL(icell,k)
		end do
		
		! 8 vertices, lets set global node numbers for corr i,j,k
		
		gnumnode(icell,1,1,1) 		= IC2V(1)
		gnumnode(icell, 1, 1, Np+1) 	= IC2V(2)
		gnumnode(icell, 1, Np+1,Np+1) 	= IC2V(3)
		gnumnode(icell, 1,Np+1,1) 	= IC2V(4)
		gnumnode(icell,Np+1,1,1) 	= IC2V(5)
		gnumnode(icell,Np+1,1,Np+1) 	= IC2V(6)
		gnumnode(icell,Np+1,Np+1,Np+1)	= IC2V(7)
		gnumnode(icell, Np+1,Np+1,1) 	= IC2V(8)

		
		do k=1,8
			gnid = IC2V(k)
			if(gnode_factor(gnid) == 0) then
				plotnodes = plotnodes + 1
				plotX(1,gnid) = XV(gnid)
				plotX(2,gnid) = YV(gnid)
				plotX(3,gnid) = ZV(gnid)
			end if
			gnode_factor(gnid) = gnode_factor(gnid) + 1
		end do				
		
	end do	! end of loop over cells
	
	! Just doing a check before we go further
	if(plotnodes/=NVERT) then
		write(*,*) " something wrong in tecplotter3dsetup 1, my friend "
!	
	end if		
    	
    ! We should now set gnumnode for other new nodes within gridcell (icell)
    
	do icell=1,NCELL

		if (CURVE_WALL==1) then
			do j=1,20
				xxs(1,j) = xxf_cell(1,j,icell)
				xxs(2,j) = xxf_cell(2,j,icell)
				xxs(3,j) = xxf_cell(3,j,icell)
			end do
		else
			do j=1,8
				xxs(1,j)=XV(IVCELL(icell,j))
                xxs(2,j)=YV(IVCELL(icell,j))
                xxs(3,j)=ZV(IVCELL(icell,j))
			end do
		end if 	
  	
    	! Lets loop over edges to check if 'nodelized'
    	
    	do k=1,12														
    		iedge=IC2E(icell,k)
    		
    		do enode=1,Np-1
    		
    				select case(k)
    				case(1)
    					ii = enode+1
    					jj = 1
    					kk = 1
    				case(2)
    					ii = Np+1 
    					jj = enode+1
    					kk = 1		
    				case(3)	
    					ii = enode+1
    					jj = Np+1
    					kk = 1
    				case(4)
    					ii = 1 
    					jj = enode+1
    					kk = 1			
    				case(5)
    					ii = enode+1
    					jj = 1
    					kk = Np+1
    				case(6)
    					ii = Np+1 
    					jj = enode+1
    					kk = Np+1
    				case(7)	
    					ii = enode+1
    					jj = Np+1
    					kk = Np+1
    				case(8)
    					ii = 1 
    					jj = enode+1
    					kk = Np+1
    				case(9)
    					ii = 1
    					jj = 1
    					kk = enode+1
    				case(10)
    					ii = Np+1
    					jj = 1
    					kk = enode+1
    				case(11)
    					ii = Np+1
    					jj = Np+1
    					kk = enode+1			
    				case(12)
    					ii = 1
    					jj = Np+1
    					kk = enode+1
    				end select		

    				xsi = dble(ii-1)/dble(Np)
    				eta = dble(jj-1)/dble(Np)
    				bet = dble(kk-1)/dble(Np)

    				if (CURVE_WALL==0.or.CURVE_WALL==1) then
						call xyCoor_atGps(kmax,xxs,xsi,eta,bet,rx)
				 	else
					 	call Coor_transfinite(xxs,xsi,eta,bet,rx)
				 	end if
    				   				
    				if(isedgenod(iedge)==1) then
    					checkmatch = 0
    					do existnode=1,Np-1
    						gnid = gedgenode(iedge,existnode)
    						if( distcalc(rx,gnid) <= tolplot) then
    							! the nodes are identical
    							gnumnode(icell,kk,jj,ii) = gnid
    							gnode_factor(gnid) = gnode_factor(gnid) + 1
    							checkmatch=1
    						end if
    					end do
    					if(checkmatch==0) then
							write(*,*) " Something wrong tecplotter3dsetup 2, edge"
							stop
    					end if			
    				
    				else
						plotnodes = plotnodes + 1
    					plotX(1,plotnodes) = rx(1)
    					plotX(2,plotnodes) = rx(2)
						plotX(3,plotnodes) = rx(3)
    					gnumnode(icell,kk,jj,ii) = plotnodes
    					gedgenode(iedge,enode) = plotnodes
    					gnode_factor(plotnodes) = gnode_factor(plotnodes) + 1
    					
    				end if		
    				
    		end do	! end of loop over nodes on edge
    		
			if(isedgenod(iedge)==0) isedgenod(iedge) = 1
			
		end do	! end of loop over 12 edges of cell
	
	   																						   		
    	! Lets loop over faces to see if nodelized
    	
    	do k=1,6
    		iface = IC2F(icell,k)
    		do mm=1,Np-1
    		do nn=1,Np-1
    			select case(k)
    			case(1)
    				ii = mm+1
    				jj = nn+1
    				kk = 1
    			case(2)
    				ii = mm+1
    				jj = nn+1
    				kk = Np+1
    			case(3) 
    				ii = mm+1
    				jj = 1
    				kk = nn+1	
    			case(4) 
    				ii = mm+1
    				jj = Np+1
    				kk = nn+1
    			case(5)
    				ii = 1
    				jj = mm+1
    				kk = nn+1					
    			case(6)
    				ii = Np+1
    				jj = mm+1
    				kk = nn+1
    			end select
    			
    			xsi = dble(ii-1)/dble(Np)
			    eta = dble(jj-1)/dble(Np)
   			    bet = dble(kk-1)/dble(Np)

   				if (CURVE_WALL==0.or.CURVE_WALL==1) then
				   	call xyCoor_atGps(kmax,xxs,xsi,eta,bet,rx)
				else
					call Coor_transfinite(xxs,xsi,eta,bet,rx)
				end if
	
				if(isfacenod(iface)==1) then
    					
    				checkmatch = 0
    				do exnn=1,Np-1
    				do exmm=1,Np-1
    					gnid = gfacenode(iface,exnn,exmm)
    					if( distcalc(rx,gnid) <= tolplot) then
    						! the nodes are identical
    						gnumnode(icell,kk,jj,ii) = gnid
    						gnode_factor(gnid) = gnode_factor(gnid) + 1
    						checkmatch=1
    					end if
    				end do
    				end do
    				
    				if(checkmatch==0) then
    					write(*,*) " Something wrong tecplotter3dsetup 2, face"
    				 
    				end if			
    			else
    			   	plotnodes = plotnodes + 1
    				plotX(1,plotnodes) = rx(1)
    				plotX(2,plotnodes) = rx(2)
    				plotX(3,plotnodes) = rx(3)
    				gnumnode(icell,kk,jj,ii) = plotnodes
    				gfacenode(iface,nn,mm) = plotnodes
    				gnode_factor(plotnodes) = gnode_factor(plotnodes) + 1
    					
    			end if		
    		end do	! end of loop over nodes over face
    		end do
			if(isfacenod(iface)==0) isfacenod(iface) = 1
		end do	! end of loop over 6 faces of cell
    	! Now to create the new nodes in the interior of cell
    	
    	do kk=2,Np
    	do jj=2,Np
        do ii=2,Np
    	    xsi = dble(ii-1)/dble(Np)
			eta = dble(jj-1)/dble(Np)
			bet = dble(kk-1)/dble(Np)
			
			if (CURVE_WALL==0.or.CURVE_WALL==1) then
				call xyCoor_atGps(kmax,xxs,xsi,eta,bet,rx)
		 	else
			 	call Coor_transfinite(xxs,xsi,eta,bet,rx)
		 	end if

    	    plotnodes = plotnodes + 1

			plotX(1,plotnodes) = rx(1)
			plotX(2,plotnodes) = rx(2)
			plotX(3,plotnodes) = rx(3)
			gnumnode(icell,kk,jj,ii) = plotnodes
			gnode_factor(plotnodes) = gnode_factor(plotnodes) + 1	
		end do	
		end do
		end do
		
		do kk = 1,Np
		do jj = 1,Np
		do ii = 1,Np
			plotcells = plotcells + 1
			gnumcell(icell,kk,jj,ii) = plotcells
			connec_c2n(1,plotcells) = gnumnode(icell,kk,jj,ii)
			connec_c2n(2,plotcells) = gnumnode(icell,kk,jj,ii+1)
			connec_c2n(3,plotcells) = gnumnode(icell,kk,jj+1,ii+1)
			connec_c2n(4,plotcells) = gnumnode(icell,kk,jj+1,ii)
			connec_c2n(5,plotcells) = gnumnode(icell,kk+1,jj,ii)
			connec_c2n(6,plotcells) = gnumnode(icell,kk+1,jj,ii+1)
			connec_c2n(7,plotcells) = gnumnode(icell,kk+1,jj+1,ii+1)
			connec_c2n(8,plotcells) = gnumnode(icell,kk+1,jj+1,ii)
		end do
		end do
		end do
		
	end do	! loop over the cells			
	
!	
!	write(*,*)'rank = ', rank, 'plotnodes,plotcells :', plotnodes, plotcells
	
	contains
	
	DOUBLE PRECISION FUNCTION distcalc(xx,gnid)
	use	setup3d
	IMPLICIT NONE
	
	integer, intent(in) :: gnid
	double precision, intent(in) :: xx(3)
		
	distcalc = sqrt(  (xx(1)-plotX(1,gnid))**2 &
	     	+ (xx(2)-plotX(2,gnid))**2 &
		+ (xx(3)-plotX(3,gnid))**2 )
					
	distcalc = abs(distcalc)				
					
	END FUNCTION distcalc
	
	
	SUBROUTINE alloc_memory
	use setup3d
	IMPLICIT NONE

	allocate(gnumcell(NCELL,NRE,NRE,NRE), connec_c2n(8,NRE**3*NCELL), plotX(3,NRE**3*NVERT) )
    allocate(gnode_factor(NRE**3*NVERT),gnumnode(NCELL,NRE+1,NRE+1,NRE+1))
    allocate(isfacenod(NFACE),gfacenode(NFACE,NRE,NRE))

	allocate(isedgenod(NEDGE),gedgenode(NEDGE,NRE))
	
	allocate(PROC_plotnodes(NPROC))
    allocate(PROC_plotcells(NPROC))


	END SUBROUTINE alloc_memory


    END SUBROUTINE tecplotter3dsetup	

