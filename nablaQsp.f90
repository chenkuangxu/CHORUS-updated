    SUBROUTINE nablaQsp
    !GET NABLAQ AT THE SOLUTION POINTS

	  use setup3d
	  implicit none
	  include 'mpif.h'

	integer :: is,js,ks,ic,rfp,k,IB
	double precision, dimension(numv) :: dQSxi,dQSeta,dQSbeta
	! double precision :: tempdedxi,tempdedeta,tempdedbeta,enei

	do ic=1,NCELL
	    do ks=1,N
	    do js=1,N
	    do is=1,N
			dQSxi(1:numv) = 0.d0
		    do rfp=1,N+1
				dQSxi(1:numv) = dQSxi(1:numv) + Qvfi(1:numv,rfp,js,ks,ic) &
				           * Mmat(rfp,is)*S1(1,1,rfp,js,ks,ic)
			    !ksi_x
			end do

			dQSeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSeta(1:numv) = dQSeta(1:numv) + Qvfj(1:numv,is,rfp,ks,ic) &
			                * Mmat(rfp,js)*S2(2,1,is,rfp,ks,ic)
				!eta_x
			end do

			dQSbeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSbeta(1:numv) = dQSbeta(1:numv) + Qvfk(1:numv,is,js,rfp,ic) &
			                 * Mmat(rfp,ks)*S3(3,1,is,js,rfp,ic)
				!beta_x
			end do

			nablaQs(1:numv,1,is,js,ks,ic) = (dQSxi(1:numv)+dQSeta(1:numv) &
			+ dQSbeta(1:numv))/Jac(is,js,ks,ic)

			!###################################################################
			dQSxi(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSxi(1:numv) = dQSxi(1:numv) + Qvfi(1:numv,rfp,js,ks,ic) &
						   * Mmat(rfp,is)*S1(1,2,rfp,js,ks,ic)
	            !ksi_y
			end do

			dQSeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSeta(1:numv) = dQSeta(1:numv) + Qvfj(1:numv,is,rfp,ks,ic) &
			                * Mmat(rfp,js)*S2(2,2,is,rfp,ks,ic)
				!eta_y
			end do

			dQSbeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSbeta(1:numv) = dQSbeta(1:numv) + Qvfk(1:numv,is,js,rfp,ic) &
			                * Mmat(rfp,ks)*S3(3,2,is,js,rfp,ic)
				!beta_y
			end do

			nablaQs(1:numv,2,is,js,ks,ic) = (dQSxi(1:numv)+dQSeta(1:numv) &
						+ dQSbeta(1:numv))/Jac(is,js,ks,ic)

			!###################################################################
			dQSxi(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSxi(1:numv) = dQSxi(1:numv) + Qvfi(1:numv,rfp,js,ks,ic)&
			                * Mmat(rfp,is)*S1(1,3,rfp,js,ks,ic)
				!ksi_z
			end do

			dQSeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSeta(1:numv) = dQSeta(1:numv) + Qvfj(1:numv,is,rfp,ks,ic)&
			                * Mmat(rfp,js)*S2(2,3,is,rfp,ks,ic)
				!eta_z
			end do

			dQSbeta(1:numv) = 0.d0
			do rfp=1,N+1
			    dQSbeta(1:numv) = dQSbeta(1:numv) + Qvfk(1:numv,is,js,rfp,ic)&
							 * Mmat(rfp,ks)*S3(3,3,is,js,rfp,ic)
				!eta_z
			end do

			nablaQs(1:numv,3,is,js,ks,ic) = (dQSxi(1:numv)+dQSeta(1:numv) &
										+ dQSbeta(1:numv))/Jac(is,js,ks,ic)
	    end do
	    end do
	    end do
	end do

	! DO IB=1,NINLET
	! 	K = IBFINL(IB)
	! 	ic = IF2C(K,1)
   
	! 	do ks=1,N
	! 	do js=1,N
	! 	do is=1,N
	! 	tempdedxi = 0.d0 
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfi(5,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic)  &
	! 			-0.5*(Qvfi(2,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(3,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(4,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(6,rfp,js,ks,ic)**2+Qvfi(7,rfp,js,ks,ic)**2+&
    !             Qvfi(8,rfp,js,ks,ic)**2+Qvfi(9,rfp,js,ks,ic)**2)/Qvfi(1,rfp,js,ks,ic)
	! 		  	tempdedxi = tempdedxi + enei*Mmat(rfp,is)
	! 	   	end do
	! 	   	dedxi(is,js,ks,ic) = tempdedxi

	! 	   	tempdedeta = 0.d0
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfj(5,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic)  &
	! 			-0.5*(Qvfj(2,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(3,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(4,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(6,is,rfp,ks,ic)**2+Qvfj(7,is,rfp,ks,ic)**2+&
	! 			Qvfj(8,is,rfp,ks,ic)**2+Qvfj(9,is,rfp,ks,ic)**2)/Qvfj(1,is,rfp,ks,ic)
	! 		  	tempdedeta = tempdedeta + enei*Mmat(rfp,js)
	! 	   	end do
	! 	   	dedeta(is,js,ks,ic) = tempdedeta

	! 	   	tempdedbeta = 0.d0
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfk(5,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic)  &
	! 			-0.5*(Qvfk(2,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(3,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(4,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(6,is,js,rfp,ic)**2+Qvfk(7,is,js,rfp,ic)**2+&
	! 			Qvfk(8,is,js,rfp,ic)**2+Qvfk(9,is,js,rfp,ic)**2)/Qvfk(1,is,js,rfp,ic)
	! 		  	tempdedbeta = tempdedbeta + enei*Mmat(rfp,ks)
	! 	   	end do
	! 	   	dedbeta(is,js,ks,ic) = tempdedbeta

	! 	end do
	! 	end do
	! 	end do
	! END DO 


	! DO IB=1,NOUTLET
	! 	K = IBFOUT(IB)
	! 	ic = IF2C(K,1)
   
	! 	do ks=1,N
	! 	do js=1,N
	! 	do is=1,N
	! 	tempdedxi = 0.d0 
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfi(5,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic)  &
	! 			-0.5*(Qvfi(2,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(3,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(4,rfp,js,ks,ic)/Qvfi(1,rfp,js,ks,ic))**2 &
	! 			-0.5*(Qvfi(6,rfp,js,ks,ic)**2+Qvfi(7,rfp,js,ks,ic)**2+&
    !             Qvfi(8,rfp,js,ks,ic)**2+Qvfi(9,rfp,js,ks,ic)**2)/Qvfi(1,rfp,js,ks,ic)
	! 		  	tempdedxi = tempdedxi + enei*Mmat(rfp,is)
	! 	   	end do
	! 	   	dedxi(is,js,ks,ic) = tempdedxi

	! 	   	tempdedeta = 0.d0
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfj(5,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic)  &
	! 			-0.5*(Qvfj(2,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(3,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(4,is,rfp,ks,ic)/Qvfj(1,is,rfp,ks,ic))**2 &
	! 			-0.5*(Qvfj(6,is,rfp,ks,ic)**2+Qvfj(7,is,rfp,ks,ic)**2+&
	! 			Qvfj(8,is,rfp,ks,ic)**2+Qvfj(9,is,rfp,ks,ic)**2)/Qvfj(1,is,rfp,ks,ic)
	! 		  	tempdedeta = tempdedeta + enei*Mmat(rfp,js)
	! 	   	end do
	! 	   	dedeta(is,js,ks,ic) = tempdedeta

	! 	   	tempdedbeta = 0.d0
	! 	   	do rfp=1,N+1
	! 		  	enei =  Qvfk(5,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic)  &
	! 			-0.5*(Qvfk(2,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(3,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(4,is,js,rfp,ic)/Qvfk(1,is,js,rfp,ic))**2 &
	! 			-0.5*(Qvfk(6,is,js,rfp,ic)**2+Qvfk(7,is,js,rfp,ic)**2+&
	! 			Qvfk(8,is,js,rfp,ic)**2+Qvfk(9,is,js,rfp,ic)**2)/Qvfk(1,is,js,rfp,ic)
	! 		  	tempdedbeta = tempdedbeta + enei*Mmat(rfp,ks)
	! 	   	end do
	! 	   	dedbeta(is,js,ks,ic) = tempdedbeta

	! 	end do
	! 	end do
	! 	end do
	! END DO 

    end subroutine nablaQsp