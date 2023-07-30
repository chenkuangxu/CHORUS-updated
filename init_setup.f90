    SUBROUTINE init_setup
    use setup3d
    implicit none
    include 'mpif.h'

    integer :: ic,iv,is,js,ks,i,j
    double precision :: xx(3),Xis,Xjs,Xks
    double precision,allocatable :: xxs(:,:)

    allocate(Xs(N),Xf(N+1))
    allocate(Lmat(N+1,N),Mmat(N+1,N))

    CALL setsandpoints(N,Xs,Xf)
    CALL computeLandM

    allocate(XXsolu(3,N,N,N,NCELL))
    allocate(XXfluxi(3,N+1,N,N,NCELL))
    allocate(XXfluxj(3,N,N+1,N,NCELL))
    allocate(XXfluxk(3,N,N,N+1,NCELL))

    allocate(Q(numv,N,N,N,NCELL))
    allocate(resid(numv,N,N,N,NCELL))

    allocate(Jac(N,N,N,NCELL))
    allocate(S1(3,3,N+1,N,N,NCELL))
    allocate(S2(3,3,N,N+1,N,NCELL))
    allocate(S3(3,3,N,N,N+1,NCELL))

    allocate(F1(numv,N+1,N,N,NCELL))
    allocate(G2(numv,N,N+1,N,NCELL))
    allocate(H3(numv,N,N,N+1,NCELL))

    allocate(Fv1(numv,N+1,N,N,NCELL))
    allocate(Gv2(numv,N,N+1,N,NCELL))
    allocate(Hv3(numv,N,N,N+1,NCELL))

    allocate(Gfp2Lfp(10,N,N,NFACE))
    allocate(map1(2,N,N,NPROCINT))
    allocate(map2(2,N,N,NCYCREM))

    if (vismode==1) then
        allocate(nablaQs(numv,3,N,N,N,NCELL))
        allocate(nablaQvfi(numv,3,N+1,N,N,NCELL))
        allocate(nablaQvfj(numv,3,N,N+1,N,NCELL))
        allocate(nablaQvfk(numv,3,N,N,N+1,NCELL))
        allocate(Qvfi(numv,N+1,N,N,NCELL))
        allocate(Qvfj(numv,N,N+1,N,NCELL))
        allocate(Qvfk(numv,N,N,N+1,NCELL))
        !!!
        allocate(dedxi(N,N,N,NCELL))
        allocate(dedeta(N,N,N,NCELL))
        allocate(dedbeta(N,N,N,NCELL))
    end if

    if (CURVE_WALL==1) then
        allocate(xxs(3,20))
    else
        allocate(xxs(3,8))
    end if

    do ic = 1,NCELL
        if (CURVE_WALL==1) then
            do iv = 1,20
                xxs(1,iv) = xxf_cell(1,iv,ic)
                xxs(2,iv) = xxf_cell(2,iv,ic)
                xxs(3,iv) = xxf_cell(3,iv,ic)
            end do
        else
            do iv = 1,8
                xxs(1,iv) = XV(IVCELL(ic,iv))
                xxs(2,iv) = YV(IVCELL(ic,iv))
                xxs(3,iv) = ZV(IVCELL(ic,iv))
            end do
        end if
        
        do ks = 1,N
        do js = 1,N
        do is = 1,N
  
            Xis = Xs(is)
            Xjs = Xs(js)
            Xks = Xs(ks)
  
            if (CURVE_WALL==1) then
                call xyCoor_atGps(20,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==0) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==2) then
                call Coor_transfinite(xxs,Xis,Xjs,Xks,xx)
            end if
  
            XXsolu(1,is,js,ks,ic) = xx(1)
            XXsolu(2,is,js,ks,ic) = xx(2)
            XXsolu(3,is,js,ks,ic) = xx(3)
        end do
        end do
        end do
  
        ! i-direction
        do ks = 1,N
        do js = 1,N
        do is = 1,N+1
  
            Xis = Xf(is)
            Xjs = Xs(js)
            Xks = Xs(ks)
  
            if (CURVE_WALL==1) then
                call xyCoor_atGps(20,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==0) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==2) then
                call Coor_transfinite(xxs,Xis,Xjs,Xks,xx)
            end if
  
            XXfluxi(1,is,js,ks,ic) = xx(1)
            XXfluxi(2,is,js,ks,ic) = xx(2)
            XXfluxi(3,is,js,ks,ic) = xx(3)
        end do
        end do
        end do
        
        ! j-direction
        do ks = 1,N
        do js = 1,N+1
        do is = 1,N
  
            Xis = Xs(is)
            Xjs = Xf(js)
            Xks = Xs(ks)
  
            if (CURVE_WALL==1) then
                call xyCoor_atGps(20,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==0) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==2) then
                call Coor_transfinite(xxs,Xis,Xjs,Xks,xx)
            end if
  
            XXfluxj(1,is,js,ks,ic) = xx(1)
            XXfluxj(2,is,js,ks,ic) = xx(2)
            XXfluxj(3,is,js,ks,ic) = xx(3)
        end do
        end do
        end do
  
        ! k-direction
        do ks = 1,N+1
        do js = 1,N
        do is = 1,N
  
            Xis = Xs(is)
            Xjs = Xs(js)
            Xks = Xf(ks)
  
            if (CURVE_WALL==1) then
                call xyCoor_atGps(20,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==0) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,Xks,xx)
            else if (CURVE_WALL==2) then
                call Coor_transfinite(xxs,Xis,Xjs,Xks,xx)
            end if
  
            XXfluxk(1,is,js,ks,ic) = xx(1)
            XXfluxk(2,is,js,ks,ic) = xx(2)
            XXfluxk(3,is,js,ks,ic) = xx(3)
        end do
        end do
        end do
  
    end do
  
    if (CURVE_WALL==1) then
        allocate(dmdxs(20,3,N,N,N))
        allocate(dmdxf1(20,3,N+1,N,N))
        allocate(dmdxf2(20,3,N,N+1,N))
        allocate(dmdxf3(20,3,N,N,N+1))
        call MAPDER
    else if (CURVE_WALL==0) then
        allocate(dmdxs(8,3,N,N,N))
        allocate(dmdxf1(8,3,N+1,N,N))
        allocate(dmdxf2(8,3,N,N+1,N))
        allocate(dmdxf3(8,3,N,N,N+1))
        call MAPDER
    end if

    call MAPFPOINT

    call MAPPROCINT

    if (CURVE_WALL==0.or.CURVE_WALL==1) call calcjacob
    if (CURVE_WALL==2) call calcjacob_transfinite

    call set_sun_model_variables

    if (restart==0) then
        call SETINITIALCOND
    else
        if (restart_ord==N) then
            CALL READ_ALL_DATA_BINARY
        else 
            CALL READ_ALL_DATA_BINARY_DIFF
        end if
    end if

    !!! Gauss integration
    call compute_Gauss_integral_coefficient
    call calcjacob_gauss_integral

    END SUBROUTINE init_setup

    SUBROUTINE set_sun_model_variables
    use setup3d
    implicit none

    double precision :: c_0,c_1,c_2,alpha_p,XR,radcomp
    integer :: ic,is,js,ks
    double precision :: dTdr,drhodr,dr,flux_t,krad,dedr

    c_0 = 156009750.d0
    c_1 = -45631718.d0
    c_2 = 3337036.8d0
    !!! alpha_p = 849.056886658866d0
    XR = 0.0000000001d0

    !!!!!
    dr = (Routlet - Rinlet) / dble(MAXSM-1)
    dTdr = (TIn(2)-TIn(1)) / dr
    dedr = R_const/(gam-1)*dTdr
    flux_t = constL/(4.d0*pi*RINLET**2.d0*constheati)
    krad = flux_t / (-rhoIn(1)*dedr*gam)
    alpha_p = krad / (c_0 + c_1*Rinlet*XR + c_2*(Rinlet*XR)**2)
    if (rank.eq.0) print *,'alpha_p=',alpha_p
    !!!!!

    allocate(radCfi(N+1,N,N,NCELL))
    allocate(radCfj(N,N+1,N,NCELL))
    allocate(radCfk(N,N,N+1,NCELL))

    allocate(gR(N,N,N,NCELL))

    DO ic = 1,NCELL

      ! i-direction

      DO ks = 1,N
      DO js = 1,N
      DO is = 1,N+1

        radcomp = sqrt(XXfluxi(1,is,js,ks,ic)**2+&
                  XXfluxi(2,is,js,ks,ic)**2+&
                  XXfluxi(3,is,js,ks,ic)**2)

        radCfi(is,js,ks,ic) = alpha_p*(c_0 + c_1*radcomp*XR + c_2*(radcomp*XR)**2)
      
      END DO
      END DO
      END DO

      ! j-direction

      DO ks = 1,N
      DO js = 1,N+1
      DO is = 1,N

        radcomp = sqrt(XXfluxj(1,is,js,ks,ic)**2+&
                  XXfluxj(2,is,js,ks,ic)**2+&
                  XXfluxj(3,is,js,ks,ic)**2)
        
        radCfj(is,js,ks,ic) = alpha_p*(c_0 + c_1*radcomp*XR + c_2*(radcomp*XR)**2)
      
      END DO
      END DO
      END DO

      ! k-direction

      DO ks = 1,N+1
      DO js = 1,N
      DO is = 1,N

        radcomp = sqrt(XXfluxk(1,is,js,ks,ic)**2+&
                  XXfluxk(2,is,js,ks,ic)**2+&
                  XXfluxk(3,is,js,ks,ic)**2)

        radCfk(is,js,ks,ic) = alpha_p*(c_0 + c_1*radcomp*XR + c_2*(radcomp*XR)**2)
      
      END DO
      END DO
      END DO
    
    END DO

    ! calculate gravitational coefficent
    do ic = 1,NCELL

        do ks = 1,N
        do js = 1,N
        do is = 1,N
  
            radcomp = sqrt(XXsolu(1,is,js,ks,ic)**2 + &
                    XXsolu(2,is,js,ks,ic)**2 + &
                    XXsolu(3,is,js,ks,ic)**2)
            gR(is,js,ks,ic) = sunG*sunM/(radcomp*radcomp)
  
        end do
        end do
        end do
  
      end do

    END SUBROUTINE set_sun_model_variables

    SUBROUTINE setsandpoints(np,Xs_tmp,Xf_tmp)
    use setup3d
    IMPLICIT NONE
    
    INTEGER :: i
    INTEGER,INTENT(IN) :: np
    DOUBLE PRECISION,INTENT(OUT) :: Xs_tmp(np),Xf_tmp(np+1)
    
    DO i = 1,np
        Xs_tmp(i) = 0.5D0*(1.0D0-cos(dble(2*i-1)*pi/(2*dble(np))))
    END DO
    
    if(np == 1) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) =  1.d0
    else if(np == 2) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) =  0.5d0
        Xf_tmp(3) =  1.d0
    else if(np == 3) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(3.d0)/3.d0)/2.d0
        Xf_tmp(3) = (1.d0 + dsqrt(3.d0)/3.d0)/2.d0
        Xf_tmp(4) =  1.d0
    else if(np == 4) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(15.d0)/5.d0)/2.d0
        Xf_tmp(3) =  0.5d0
        Xf_tmp(4) = (1.d0 + dsqrt(15.d0)/5.d0)/2.d0
        Xf_tmp(5) =  1.d0
    else if(np == 5) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(525.d0+70.d0*dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(3) = (1.d0 - dsqrt(525.d0-70.d0*dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(4) = (1.d0 + dsqrt(525.d0-70.d0*dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(5) = (1.d0 + dsqrt(525.d0+70.d0*dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(6) =  1.d0
    else if(np == 6) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1. - dsqrt(245.d0+14.d0*dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(3) = (1. - dsqrt(245.d0-14.d0*dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(4) =  0.5d0
        Xf_tmp(5) = (1. + dsqrt(245.d0-14.d0*dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(6) = (1. + dsqrt(245.d0+14.d0*dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(7) =  1.d0
    else if(np == 7) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - 0.932469514203d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.661209386466d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.238619186083d0)/2.d0
        Xf_tmp(5) = (1.d0 + 0.238619186083d0)/2.d0
        Xf_tmp(6) = (1.d0 + 0.661209386466d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.932469514203d0)/2.d0
        Xf_tmp(8) =  1.d0
    else if(np == 8) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - 0.949107912343d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.741531185599d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.405845151377d0)/2.d0
        Xf_tmp(5) =  0.5d0
        Xf_tmp(6) = (1.d0 + 0.405845151377d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.741531185599d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.949107912343d0)/2.d0
        Xf_tmp(9) =  1.d0
    else if(np == 9) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - 0.960289856498d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.796666477414d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.525532409916d0)/2.d0
        Xf_tmp(5) = (1.d0 - 0.183434642496d0)/2.d0
        Xf_tmp(6) = (1.d0 + 0.183434642496d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.525532409916d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.796666477414d0)/2.d0
        Xf_tmp(9) = (1.d0 + 0.960289856498d0)/2.d0
        Xf_tmp(10)=  1.d0
    else if(np == 10) then
        Xf_tmp(1) =  0.d0
        Xf_tmp(2) = (1.d0 - 0.968160239508d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.836031107327d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.613371432701d0)/2.d0
        Xf_tmp(5) = (1.d0 - 0.324253423404d0)/2.d0
        Xf_tmp(6) =  0.5d0
        Xf_tmp(7) = (1.d0 + 0.324253423404d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.613371432701d0)/2.d0
        Xf_tmp(9) = (1.d0 + 0.836031107327d0)/2.d0
        Xf_tmp(10)= (1.d0 + 0.968160239508d0)/2.d0
        Xf_tmp(11)=  1.d0
    end if
    
    END SUBROUTINE setsandpoints
    
    
    SUBROUTINE computeLandM
    
    use setup3d
    IMPLICIT NONE
    
    integer :: np
    integer :: is,ifp,s,k,js
    double precision :: hhval
    double precision :: num, den, xval, sum
    
    np = N
    
    do ifp=1,np+1
        do is=1,np
            Lmat(ifp,is) = hhval(is,Xf(ifp))
        end do
    end do
    
    do ifp=1,np+1
        do is=1,np
       
            xval = Xs(is)
            sum = 0.0d0
            do k=1,np+1
                if(k/=ifp) then
                    num = 1.0d0
                    den = 1.0d0
                    do s=1,np+1
                    if(s/=ifp .and. s/=k) then
                        num = num * (xval - Xf(s))
                    end if
                    if(s/=ifp) then
                        den = den * (Xf(ifp) - Xf(s))
                    end if
       
                    end do
       
                    sum = sum + num/den
                end if
            end do  ! end of loop over k
       
            Mmat(ifp,is) = sum
       
        end do
    end do
    
    ! write(*,*) "Lmat"
    ! do ifp=1,np+1
    !   write(*,*) (Lmat(ifp,is), is=1,np)
    ! end do
    
    ! write(*,*) "Mmat"
    ! do ifp=1,np+1
    !   write(*,*) (Mmat(ifp,is), is=1,np)
    ! end do
    
    ! write (*,*) 'N=', np, 'Lmat and Mmat is generated'
    
    END SUBROUTINE computeLandM
    
    
    FUNCTION hhval(i,xval)
    use setup3d, only: n, Xs
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: i
    INTEGER :: np
    double precision,intent(in) :: xval
    double precision :: hhval,hvaln,hvald
    integer :: s
    
    hvaln = 1.0d0
    hvald = 1.0d0
    np=n
    do s=1,np
        if( s/=i) then
            hvaln = hvaln * ( xval - Xs(s) )
            hvald = hvald * ( Xs(i) - Xs(s) )
        end if
    end do
    
    hhval = hvaln/hvald
    
    END FUNCTION hhval