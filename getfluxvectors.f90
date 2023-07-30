        SUBROUTINE getfluxvectors(Qv,Fv,Gv,Hv)
	
	use setup3d
	IMPLICIT NONE
        include 'mpif.h'

	double precision,dimension(5), intent(in)	::	Qv
	double precision,dimension(5), intent(out)	::	Fv,Gv,Hv
	double precision::rho,rhou,rhov,rhow,rhoE,pr

	! the conservative variables	
	rho  = Qv(1)
	rhou = Qv(2)
	rhov = Qv(3)
	rhow = Qv(4)
	rhoE = Qv(5)

	! compute the pressure
	pr = (gam-1)*( rhoE - 0.5 * (rhou**2 + rhov**2+rhow**2)/rho )

	! compute flux vector F
	Fv(1) = rhou
	Fv(2) = (rhou**2)/rho + pr		
	Fv(3) = rhou*rhov / rho
	Fv(4) = rhou*rhow / rho
	Fv(5) = ( rhou/rho )*( rhoE+pr )
				
	! compute flux vector G
	Gv(1) = rhov
	Gv(2) = rhou*rhov / rho		
	Gv(3) = (rhov**2)/rho + pr 
	Gv(4) = rhov*rhow / rho		
	Gv(5) = ( rhov/rho )*( rhoE+pr )	

	Hv(1) = rhow
	Hv(2) = rhou*rhow / rho		
	Hv(3) = rhov*rhow / rho		
	Hv(4) = (rhow**2)/rho + pr 
	Hv(5) = ( rhow/rho )*( rhoE+pr )	

        END SUBROUTINE getfluxvectors	





