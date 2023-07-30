    subroutine compute_Gauss_integral_coefficient
    use setup3d,only: G_omega,G_psi,int_num,rank
    IMPLICIT NONE
    include 'mpif.h'
    
    integer :: I

    int_num = 8
    
    allocate(G_omega(int_num))
    allocate(G_psi(int_num))
    
    if(int_num.eq.1) then
        G_omega(1) = 2.0d0
    
        G_psi(1) = 0.0d0
    else if(int_num.eq.2) then
        G_omega(1) = 1.0d0
        G_omega(2) = 1.0d0
     
        G_psi(1) = -sqrt(1.0d0/3.0d0)
        G_psi(2) =  sqrt(1.0d0/3.0d0)
    else if(int_num.eq.3) then
        G_omega(1) = 5.0d0/9.0d0
        G_omega(2) = 8.0d0/9.0d0
        G_omega(3) = 5.0d0/9.0d0
    
        G_psi(1) = -sqrt(3.0d0/5.0d0)
        G_psi(2) =  0.0d0
        G_psi(3) =  sqrt(3.0d0/5.0d0)
    else if(int_num.eq.4) then
        G_omega(1) = (18.0d0-sqrt(30.0d0))/36.0d0
        G_omega(2) = (18.0d0+sqrt(30.0d0))/36.0d0
        G_omega(3) = (18.0d0+sqrt(30.0d0))/36.0d0
        G_omega(4) = (18.0d0-sqrt(30.0d0))/36.0d0
    
        G_psi(1) = -sqrt(3.0d0/7.0d0+2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
        G_psi(2) = -sqrt(3.0d0/7.0d0-2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
        G_psi(3) =  sqrt(3.0d0/7.0d0-2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
        G_psi(4) =  sqrt(3.0d0/7.0d0+2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
    else if(int_num.eq.5) then
        G_omega(1) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        G_omega(2) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        G_omega(3) = 128.0d0/225.0d0
        G_omega(4) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        G_omega(5) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
    
        G_psi(1) = -1.0d0/3.0d0*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
        G_psi(2) = -1.0d0/3.0d0*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        G_psi(3) =  0.0d0
        G_psi(4) =  1.0d0/3.0d0*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        G_psi(5) =  1.0d0/3.0d0*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
    else if(int_num.eq.6) then
        G_omega(1) = 0.1713244923791704d0
        G_omega(2) = 0.3607615730481386d0
        G_omega(3) = 0.4679139345726910d0
        G_omega(4) = 0.4679139345726910d0
        G_omega(5) = 0.3607615730481386d0
        G_omega(6) = 0.1713244923791704d0
    
        G_psi(1) = -0.9324695142031521d0
        G_psi(2) = -0.6612093864662645d0
        G_psi(3) = -0.2386191860831969d0
        G_psi(4) = 0.2386191860831969d0
        G_psi(5) = 0.6612093864662645d0
        G_psi(6) = 0.9324695142031521d0
    else if(int_num.eq.7) then
        G_omega(1) = 0.1294849661688697d0
        G_omega(2) = 0.2797053914892766d0
        G_omega(3) = 0.3818300505051189d0
        G_omega(4) = 0.4179591836734694d0
        G_omega(5) = 0.3818300505051189d0
        G_omega(6) = 0.2797053914892766d0
        G_omega(7) = 0.1294849661688697d0
    
        G_psi(1) = -0.9491079123427585d0
        G_psi(2) = -0.7415311855993945d0
        G_psi(3) = -0.4058451513773972d0
        G_psi(4) = 0.d0
        G_psi(5) = 0.4058451513773972d0
        G_psi(6) = 0.7415311855993945d0
        G_psi(7) = 0.9491079123427585d0
    else if(int_num.eq.8) then
        G_omega(1) = 0.1012285362903763d0
        G_omega(2) = 0.2223810344533745d0
        G_omega(3) = 0.3137066458778873d0
        G_omega(4) = 0.3626837833783620d0
        G_omega(5) = 0.3626837833783620d0
        G_omega(6) = 0.3137066458778873d0
        G_omega(7) = 0.2223810344533745d0
        G_omega(8) = 0.1012285362903763d0
    
        G_psi(1) = -0.9602898564975363d0
        G_psi(2) = -0.7966664774136267d0
        G_psi(3) = -0.5255324099163290d0
        G_psi(4) = -0.1834346424956498d0
        G_psi(5) = 0.1834346424956498d0
        G_psi(6) = 0.5255324099163290d0
        G_psi(7) = 0.7966664774136267d0
        G_psi(8) = 0.9602898564975363d0
    else if(int_num.eq.9) then
        G_omega(1) = 0.0812743883615744d0
        G_omega(2) = 0.1806481606948574d0
        G_omega(3) = 0.2606106964029354d0
        G_omega(4) = 0.3123470770400029d0
        G_omega(5) = 0.3302393550012598d0
        G_omega(6) = 0.3123470770400029d0
        G_omega(7) = 0.2606106964029354d0
        G_omega(8) = 0.1806481606948574d0
        G_omega(9) = 0.0812743883615744d0
    
        G_psi(1) = -0.9681602395076261d0
        G_psi(2) = -0.8360311073266358d0
        G_psi(3) = -0.6133714327005904d0
        G_psi(4) = -0.3242534234038089d0
        G_psi(5) = 0.d0
        G_psi(6) = 0.3242534234038089d0
        G_psi(7) = 0.6133714327005904d0
        G_psi(8) = 0.8360311073266358d0
        G_psi(9) = 0.9681602395076261d0
    else if(int_num.eq.10) then
        G_omega(1) = 0.0666713443086881d0
        G_omega(2) = 0.1494513491505806d0
        G_omega(3) = 0.2190863625159820d0
        G_omega(4) = 0.2692667193099963d0
        G_omega(5) = 0.2955242247147529d0
        G_omega(6) = 0.2955242247147529d0
        G_omega(7) = 0.2692667193099963d0
        G_omega(8) = 0.2190863625159820d0
        G_omega(9) = 0.1494513491505806d0
        G_omega(10)= 0.0666713443086881d0
    
        G_psi(1) = -0.9739065285171717d0
        G_psi(2) = -0.8650633666889845d0
        G_psi(3) = -0.6794095682990244d0
        G_psi(4) = -0.4333953941292472d0
        G_psi(5) = -0.1488743389816312d0
        G_psi(6) = 0.1488743389816312d0
        G_psi(7) = 0.4333953941292472d0
        G_psi(8) = 0.6794095682990244d0
        G_psi(9) = 0.8650633666889845d0
        G_psi(10)= 0.9739065285171717d0
    else
        if (rank.eq.0) print *,'Maximum order of Gauss integration is 10 !!!'
        stop
    endif
    
    do I = 1,int_num
        G_psi(I) = (G_psi(I)+1.d0)/2.d0
        G_omega(I) = G_omega(I)/2.d0
    end do
    
    end subroutine compute_Gauss_integral_coefficient

    SUBROUTINE calcjacob_gauss_integral
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'
    
    integer :: ic,k,is,js,ks
    double precision :: xi,eta,beta,xxs(3,8),coor_der(3,3)
    double precision :: xxi,yxi,zxi,xeta,yeta,zeta,xbeta,ybeta,zbeta

    allocate(JacGa(int_num,int_num,int_num,NCELL))
    
    DO ic = 1,NCELL
    
        DO k = 1,8
            xxs(1,k) = XV(IVCELL(ic,k))
            xxs(2,k) = YV(IVCELL(ic,k))
            xxs(3,k) = ZV(IVCELL(ic,k))
        END DO
    
        DO ks = 1,int_num
        DO js = 1,int_num
        DO is = 1,int_num
    
            xi = G_psi(is)
            eta = G_psi(js)
            beta = G_psi(ks)
    
            CALL Deri_transfinite(xxs,xi,eta,beta,coor_der)
    
            xxi = coor_der(1,1)
            yxi = coor_der(2,1)
            zxi = coor_der(3,1)
    
            xeta = coor_der(1,2)
            yeta = coor_der(2,2)
            zeta = coor_der(3,2)
    
            xbeta = coor_der(1,3)
            ybeta = coor_der(2,3)
            zbeta = coor_der(3,3)
    
            JacGa(is,js,ks,ic) = xxi*yeta*zbeta + xeta*ybeta*zxi+xbeta*yxi*zeta  &
            - xxi*ybeta*zeta - xeta*yxi*zbeta - xbeta*yeta*zxi
    
        END DO
        END DO
        END DO
    
    END DO
    
    END SUBROUTINE calcjacob_gauss_integral

    SUBROUTINE compute_global_kE
    use setup3d
    implicit none
    include 'mpif.h'

    integer :: ig,jg,kg,is,js,ks,ic
    double precision :: xsi,eta,bet,Qs(numv),hhval,kE,v2,rho,volume

    kE = 0.d0
    do ic = 1,NCELL
        do kg = 1,int_num
        do jg = 1,int_num
        do ig = 1,int_num

            Qs = 0.d0
            xsi = G_psi(ig)
            eta = G_psi(jg)
            bet = G_psi(kg)
            do ks=1,N
            do js=1,N
            do is=1,N
                Qs(1:numv) = Qs(1:numv) + Q(1:numv,is,js,ks,ic) &
                * hhval(is,xsi) * hhval(js,eta) * hhval(ks,bet)
            end do
            end do
            end do

            v2 = (Qs(2)/Qs(1))**2+(Qs(3)/Qs(1))**2+(Qs(4)/Qs(1))**2
            rho= Qs(1)
            kE = kE + 0.5d0*rho*v2*&
            G_omega(ig)*G_omega(jg)*G_omega(kg)*JacGa(ig,jg,kg,ic)

        end do
        end do
        end do
    end do

    CALL MPI_ALLREDUCE(kE,kE_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

    volume = 4.d0/3*pi*(Routlet**3.d0-Rinlet**3.d0)

    kE_g = kE_g / volume

    if (rank.eq.0) then
        open(20,file='kE.out',position='append')
        write(20,*) iter,ctime,kE_g
        close(20)
    end if

    END SUBROUTINE compute_global_kE
