	SUBROUTINE getedges
	use setup3d	
	IMPLICIT NONE
	INCLUDE 'mpif.h'

	integer	:: icell, n1,n2 ,gedge, en1,en2
	integer	:: IVEDGE(12*ncell,2)

	if (rank.eq.0) print *,'getting edges ...'
	
    allocate(IC2E(NCELL,12))

	NEDGE = 0
	do icell=1,NCELL
		! edge 1 - (v1-v2)
		n1 = IVCELL(icell,1)
		n2 = IVCELL(icell,2)
		gedge = gnumedge(n1,n2)
		IC2E(icell,1) = gedge
				
		! edge 2 - (v2-v3)
		n1 = IVCELL(icell,2)
		n2 = IVCELL(icell,3)
		gedge = gnumedge(n1,n2)
		IC2E(icell,2) = gedge
					
		! edge 3 - (v3-v4)
		n1 = IVCELL(icell,3)
		n2 = IVCELL(icell,4)
		gedge = gnumedge(n1,n2)
		IC2E(icell,3) = gedge
		
		! edge 4 - (v4-v1)
		n1 = IVCELL(icell,4)
		n2 = IVCELL(icell,1)
		gedge = gnumedge(n1,n2)
		IC2E(icell,4) = gedge
		
		! edge 5 - (v5-v6)
		n1 = IVCELL(icell,5)
		n2 = IVCELL(icell,6)
		gedge = gnumedge(n1,n2)
		IC2E(icell,5) = gedge
		
		! edge 6 - (v6-v7)
		n1 = IVCELL(icell,6)
		n2 = IVCELL(icell,7)
		gedge = gnumedge(n1,n2)
		IC2E(icell,6) = gedge
		
		! edge 7 - (v7-v8)
		n1 = IVCELL(icell,7)
		n2 = IVCELL(icell,8)
		gedge = gnumedge(n1,n2)
		IC2E(icell,7) = gedge

		! edge 8 - (v8-v5)
		n1 = IVCELL(icell,8)
		n2 = IVCELL(icell,5)
		gedge = gnumedge(n1,n2)
		IC2E(icell,8) = gedge

		! edge 9 - (v1-v5)
		n1 = IVCELL(icell,1)
		n2 = IVCELL(icell,5)
		gedge = gnumedge(n1,n2)
		IC2E(icell,9) = gedge

		! edge 10 - (v2-v6)
		n1 = IVCELL(icell,2)
		n2 = IVCELL(icell,6)
		gedge = gnumedge(n1,n2)
		IC2E(icell,10) = gedge

		! edge 11 - (v3-v7)
		n1 = IVCELL(icell,3)
		n2 = IVCELL(icell,7)
		gedge = gnumedge(n1,n2)
		IC2E(icell,11) = gedge

		! edge 12 - (v4-v8)
		n1 = IVCELL(icell,4)
		n2 = IVCELL(icell,8)
		gedge = gnumedge(n1,n2)
		IC2E(icell,12) = gedge

	end do	! end of loop over cells
			
	contains
		
	INTEGER FUNCTION gnumedge(n1,n2)
		
		IMPLICIT NONE
		
		integer, intent(in) :: n1, n2
		integer :: checkedge, iedge, en1, en2
										
		! check existing edges if already present
		checkedge = 0
		iedge = 1
		do while( iedge<=NEDGE .and. checkedge==0 )
			en1 = IVEDGE(iedge,1)
			en2 = IVEDGE(iedge,2)
				
			if( (en1==n1 .and. en2==n2) .or. (en1==n2 .and. en2==n1) ) then
				checkedge = 1	
			else
				iedge = iedge + 1
			end if	
		end do	
		
		if(checkedge==1) then
			gnumedge = iedge
		else
			NEDGE = NEDGE + 1
			IVEDGE(NEDGE,1) = n1
			IVEDGE(NEDGE,2) = n2
			gnumedge = NEDGE
		end if
		
	END FUNCTION gnumedge					
					
	END SUBROUTINE getedges		
