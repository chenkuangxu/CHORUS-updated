    function point2hex(xx)
    implicit none
    integer :: point2hex
    double precision :: xx(3),radius,xxn(3),xxf(3)
    double precision :: tol,ratio
    
    tol = 1d-6
    point2hex = 0

    radius = norm2(xx)
    xxn = xx/radius

    if (xxn(1).GE.0.d0) then
        ratio = 1.d0/xxn(1)
        xxf = xxn*ratio
        if ( (xxf(2)+1.d0).GE.tol.AND.(1.d0-xxf(2)).GE.tol&
        .AND.(xxf(3)+1.d0).GE.tol.AND.(1.d0-xxf(3)).GE.tol) then
            point2hex = 1
        end if
    else
        ratio = -1.d0/xxn(1)
        xxf = xxn*ratio
        if ( (xxf(2)+1.d0).GE.tol.AND.(1.d0-xxf(2)).GE.tol&
        .AND.(xxf(3)+1.d0).GE.tol.AND.(1.d0-xxf(3)).GE.tol) then
            point2hex = 3
        end if
    end if
    
    if (xxn(2).GE.0.d0) then
        ratio = 1.d0/xxn(2)
        xxf = xxn*ratio
        if ( (xxf(1)+1.d0).GE.tol.AND.(1.d0-xxf(1)).GE.tol&
        .AND.(xxf(3)+1.d0).GE.tol.AND.(1.d0-xxf(3)).GE.tol) then
            point2hex = 2
        end if
    else
        ratio = -1.d0/xxn(2)
        xxf = xxn*ratio
        if ( (xxf(1)+1.d0).GE.tol.AND.(1.d0-xxf(1)).GE.tol&
        .AND.(xxf(3)+1.d0).GE.tol.AND.(1.d0-xxf(3)).GE.tol) then
            point2hex = 4
        end if
    end if

    if (xxn(3).GE.0.d0) then
        ratio = 1.d0/xxn(3)
        xxf = xxn*ratio
        if ( (xxf(1)+1.d0).GE.tol.AND.(1.d0-xxf(1)).GE.tol&
        .AND.(xxf(2)+1.d0).GE.tol.AND.(1.d0-xxf(2)).GE.tol) then
            point2hex = 5
        end if
    else
        ratio = -1.d0/xxn(3)
        xxf = xxn*ratio
        if ( (xxf(1)+1.d0).GE.tol.AND.(1.d0-xxf(1)).GE.tol&
        .AND.(xxf(2)+1.d0).GE.tol.AND.(1.d0-xxf(2)).GE.tol) then
            point2hex = 6
        end if
    end if    
    
    end function point2hex

    function edge2hex(xx1,xx2)
    implicit none
    integer :: edge2hex,point2hex,num1,num2
    double precision :: xx1(3),xx2(3)

    num1 = point2hex(xx1)
    num2 = point2hex(xx2)

    if (num1.eq.0.and.num2.eq.0) then
        edge2hex = 0
    else if (num1.eq.0) then
        edge2hex = num2
    else if (num2.eq.0) then
        edge2hex = num1
    else if (num1.eq.num2) then
        edge2hex = num1
    else
        print *,'edge cannot find which hex it belongs to !!!'
    end if

    end function edge2hex

    function cell2hex(xxs)
    implicit none
    integer :: cell2hex,edge_num(12),edge2hex,iedge,stat,k,stat_hex(6)
    double precision :: xxs(3,8)
    
    edge_num(1) = edge2hex(xxs(1:3,1),xxs(1:3,2))
    edge_num(2) = edge2hex(xxs(1:3,2),xxs(1:3,3))
    edge_num(3) = edge2hex(xxs(1:3,3),xxs(1:3,4))
    edge_num(4) = edge2hex(xxs(1:3,4),xxs(1:3,1))
    edge_num(5) = edge2hex(xxs(1:3,5),xxs(1:3,6))
    edge_num(6) = edge2hex(xxs(1:3,6),xxs(1:3,7))
    edge_num(7) = edge2hex(xxs(1:3,7),xxs(1:3,8))
    edge_num(8) = edge2hex(xxs(1:3,8),xxs(1:3,5))
    edge_num(9) = edge2hex(xxs(1:3,1),xxs(1:3,5))
    edge_num(10)= edge2hex(xxs(1:3,2),xxs(1:3,6))
    edge_num(11)= edge2hex(xxs(1:3,3),xxs(1:3,7))
    edge_num(12)= edge2hex(xxs(1:3,4),xxs(1:3,8))

    stat = 0
    do iedge = 1,12
        if (edge_num(iedge).ne.0) then
            stat = stat+1
            if (stat.eq.1) then
                cell2hex = edge_num(iedge)
            else if (stat.GT.1) then
                if (cell2hex.ne.edge_num(iedge)) then
                    print *,'cell cannot find which hex it belongs to !!!'
                end if
            end if    
        end if
    end do

    if (stat.eq.0) then
        stat_hex = 0
        do k = 1,8
            if (xxs(1,k).gt.0.d0) stat_hex(1) = stat_hex(1) + 1
            if (xxs(1,k).lt.0.d0) stat_hex(3) = stat_hex(3) + 1
            if (xxs(2,k).gt.0.d0) stat_hex(2) = stat_hex(2) + 1
            if (xxs(2,k).lt.0.d0) stat_hex(4) = stat_hex(4) + 1
            if (xxs(3,k).gt.0.d0) stat_hex(5) = stat_hex(5) + 1
            if (xxs(3,k).lt.0.d0) stat_hex(6) = stat_hex(6) + 1
        end do

        do k = 1,6
            if (stat_hex(k).eq.8) then
                cell2hex = k
            end if
        end do
    end if

    end function cell2hex


    subroutine xyz2surf(xx,num_hex,xsi,bet,radius)
    implicit none

    double precision :: xx(3),xsi,bet,radius
    integer :: num_hex
   
    if (num_hex.eq.1) then
        xsi = datan(xx(2)/xx(1))
        bet = datan(xx(3)/xx(1))
    else if (num_hex.eq.2) then
        xsi = datan(-xx(1)/xx(2))
        bet = datan(xx(3)/xx(2))
    else if (num_hex.eq.3) then
        xsi = datan(xx(2)/xx(1))
        bet = datan(-xx(3)/xx(1))
    else if (num_hex.eq.4) then
        xsi = datan(-xx(1)/xx(2))
        bet = datan(-xx(3)/xx(2))
    else if (num_hex.eq.5) then
        xsi = datan(xx(2)/xx(3))
        bet = datan(-xx(1)/xx(3))
    else if (num_hex.eq.6) then
        xsi = datan(-xx(2)/xx(3))
        bet = datan(-xx(1)/xx(3))
    else
        print *,'error !'
    end if

    radius = norm2(xx)

    end subroutine xyz2surf

    subroutine surf2xyz(xsi,bet,radius,num_hex,xx)
    implicit none

    double precision :: xsi,bet,radius,xx(3),xxn(3)
    integer :: num_hex

    if (num_hex.eq.1) then
        xx(1) = 1.d0
        xx(2) = dtan(xsi)*xx(1)
        xx(3) = dtan(bet)*xx(1)
    else if (num_hex.eq.2) then
        xx(2) = 1.d0
        xx(1) = -dtan(xsi)*xx(2)
        xx(3) = dtan(bet)*xx(2)
    else if (num_hex.eq.3) then
        xx(1) = -1.d0
        xx(2) = dtan(xsi)*xx(1)
        xx(3) = -dtan(bet)*xx(1)
    else if (num_hex.eq.4) then
        xx(2) = -1.d0
        xx(1) = -dtan(xsi)*xx(2)
        xx(3) = -dtan(bet)*xx(2)
    else if (num_hex.eq.5) then
        xx(3) = 1.d0
        xx(2) = dtan(xsi)*xx(3)
        xx(1) = -dtan(bet)*xx(3)
    else if (num_hex.eq.6) then
        xx(3) = -1.d0
        xx(2) = -dtan(xsi)*xx(3)
        xx(1) = -dtan(bet)*xx(3)
    else
        print *,'error !'
    end if

    xxn = xx/norm2(xx)
    xx = xxn*radius

    end subroutine surf2xyz

    subroutine get_vec(vec1,vec2,zeta,vec3,num_hex)
    implicit none
    integer :: num_hex
    double precision,dimension(3),intent(in) :: vec1,vec2
    double precision :: zeta
    double precision,dimension(3),intent(out) :: vec3
    double precision :: tol,rad1,rad2,rad
    double precision :: xsi1,bet1,xsi2,bet2,xsi,bet
    
    call xyz2surf(vec1,num_hex,xsi1,bet1,rad1)
    call xyz2surf(vec2,num_hex,xsi2,bet2,rad2)

    tol = 1d-4

    if (abs(xsi1-xsi2).LT.tol.AND.abs(bet1-bet2).LT.tol) then
        vec3 = (1.d0-zeta)*vec1 + zeta*vec2
    else
        rad = (rad1+rad2)/2.d0
        if (abs(xsi1-xsi2).LT.tol) then
            xsi = (xsi1+xsi2)/2.d0
            bet = (1.d0-zeta)*bet1 + zeta*bet2
        else if (abs(bet1-bet2).LT.tol) then
            bet = (bet1+bet2)/2.d0
            xsi = (1.d0-zeta)*xsi1 + zeta*xsi2
        else
            print *,'something wrong ! num_hex=',num_hex
        end if
        call surf2xyz(xsi,bet,rad,num_hex,vec3)
    end if
    
    end subroutine get_vec

    subroutine get_surface_vec(vec1,vec2,vec3,vec4,zeta_x,zeta_y,vec5,num_hex)
    implicit none
    double precision,dimension(3) :: vec1,vec2,vec3,vec4,vec5
    double precision :: zeta_x,zeta_y
    double precision :: xsi1,bet1,xsi2,bet2,xsi3,bet3,xsi4,bet4,tol,xsi,bet
    double precision :: rad1,rad2,rad3,rad4,rad
    integer :: num_hex
    
    tol = 1d-4
    call xyz2surf(vec1,num_hex,xsi1,bet1,rad1)
    call xyz2surf(vec2,num_hex,xsi2,bet2,rad2)
    call xyz2surf(vec3,num_hex,xsi3,bet3,rad3)
    call xyz2surf(vec4,num_hex,xsi4,bet4,rad4)

    if ((abs(xsi1-xsi2).LT.tol.AND.abs(bet1-bet2).LT.tol).AND.&
    (abs(xsi3-xsi4).LT.tol.AND.abs(bet3-bet4).LT.tol)) then
        rad = (rad1+rad4)/2*(1.d0-zeta_x)+(rad2+rad3)/2*zeta_x
        if (abs(xsi2-xsi3).LT.tol.AND.abs(xsi4-xsi1).LT.tol) then
            xsi = (xsi1+xsi2+xsi3+xsi4)/4
            bet = (bet1+bet2)/2*(1-zeta_y)+(bet3+bet4)/2*zeta_y
        else if (abs(bet2-bet3).LT.tol.AND.abs(bet4-bet1).LT.tol) then
            bet = (bet1+bet2+bet3+bet4)/4
            xsi = (xsi1+xsi2)/2*(1-zeta_y)+(xsi3+xsi4)/2*zeta_y
        end if

    else if ((abs(xsi2-xsi3).LT.tol.AND.abs(bet2-bet3).LT.tol).AND.&
    (abs(xsi1-xsi4).LT.tol.AND.abs(bet1-bet4).LT.tol)) then
        rad = (rad1+rad2)/2*(1-zeta_y)+(rad3+rad4)/2*zeta_y
        if (abs(xsi1-xsi2).LT.tol.AND.abs(xsi3-xsi4).LT.tol) then
            xsi = (xsi1+xsi2+xsi3+xsi4)/4
            bet = (bet1+bet4)/2*(1-zeta_x)+(bet2+bet3)/2*zeta_x
        else if (abs(bet1-bet2).LT.tol.AND.abs(xsi3-xsi4).LT.tol) then
            bet = (bet1+bet2+bet3+bet4)/4
            xsi = (xsi1+xsi4)/2*(1-zeta_x)+(xsi2+xsi3)/2*zeta_x
        end if

    else
        rad = (rad1+rad2+rad3+rad4)/4
        if ((abs(xsi1-xsi2).LT.tol.AND.abs(xsi3-xsi4).LT.tol).AND.&
        (abs(bet1-bet4).LT.tol.AND.abs(bet2-bet3).LT.tol)) then
            xsi = (xsi1+xsi2)/2*(1-zeta_y)+(xsi3+xsi4)/2*zeta_y
            bet = (bet1+bet4)/2*(1-zeta_x)+(bet2+bet3)/2*zeta_x
        else if ((abs(bet1-bet2).LT.tol.AND.abs(bet3-bet4).LT.tol).AND.&
        (abs(xsi1-xsi4).LT.tol.AND.abs(xsi2-xsi3).LT.tol)) then
            xsi = (xsi1+xsi4)/2*(1-zeta_x)+(xsi2+xsi3)/2*zeta_x
            bet = (bet1+bet2)/2*(1-zeta_y)+(bet3+bet4)/2*zeta_y
        else
            print *,'error !'
        end if
        
    end if

    call surf2xyz(xsi,bet,rad,num_hex,vec5)

    end subroutine get_surface_vec


    SUBROUTINE Coor_transfinite(xxs,xi,eta,beta,xx)
    IMPLICIT NONE
    DOUBLE PRECISION:: xxs(3,8),xx(3),xi,eta,beta
    DOUBLE PRECISION,DIMENSION(0:1,0:1,3) :: r_xi,r_eta,r_beta
    DOUBLE PRECISION,DIMENSION(3) :: r_xi0,r_xi1,r_eta0,r_eta1,r_beta0,r_beta1
    DOUBLE PRECISION,DIMENSION(3) :: P1,P2,P3
    DOUBLE PRECISION,DIMENSION(3) :: P1P2,P2P3,P1P3,P1P2P3
    DOUBLE PRECISION,DIMENSION(3) :: vec1,vec2,vec3,vec5
    INTEGER :: num_hex,cell2hex
    
    num_hex = cell2hex(xxs)

    ! r_beta
    ! v1-v5
    vec1 = xxs(:,1)
    vec2 = xxs(:,5)
    CALL get_vec(vec1,vec2,beta,vec3,num_hex)
    r_beta(0,0,:) = vec3(:)

    ! v4-v8
    vec1 = xxs(:,4)
    vec2 = xxs(:,8)
    CALL get_vec(vec1,vec2,beta,vec3,num_hex)
    r_beta(0,1,:) = vec3(:)

    ! v2-v6
    vec1 = xxs(:,2)
    vec2 = xxs(:,6)
    CALL get_vec(vec1,vec2,beta,vec3,num_hex)
    r_beta(1,0,:) = vec3(:)

    ! v3-v7
    vec1 = xxs(:,3)
    vec2 = xxs(:,7)
    CALL get_vec(vec1,vec2,beta,vec3,num_hex)
    r_beta(1,1,:) = vec3(:)

    ! r_xi
    ! v1-v2
    vec1 = xxs(:,1)
    vec2 = xxs(:,2)
    CALL get_vec(vec1,vec2,xi,vec3,num_hex)
    r_xi(0,0,:) = vec3(:)

    ! v5-v6
    vec1 = xxs(:,5)
    vec2 = xxs(:,6)
    CALL get_vec(vec1,vec2,xi,vec3,num_hex)
    r_xi(0,1,:) = vec3(:)

    ! v4-v3
    vec1 = xxs(:,4)
    vec2 = xxs(:,3)
    CALL get_vec(vec1,vec2,xi,vec3,num_hex)
    r_xi(1,0,:) = vec3(:)

    ! v8-v7
    vec1 = xxs(:,8)
    vec2 = xxs(:,7)
    CALL get_vec(vec1,vec2,xi,vec3,num_hex)
    r_xi(1,1,:) = vec3(:)

    ! r_eta
    ! v1-v4
    vec1 = xxs(:,1)
    vec2 = xxs(:,4)
    CALL get_vec(vec1,vec2,eta,vec3,num_hex)
    r_eta(0,0,:) = vec3(:)

    ! v5-v8
    vec1 = xxs(:,5)
    vec2 = xxs(:,8)
    CALL get_vec(vec1,vec2,eta,vec3,num_hex)
    r_eta(0,1,:) = vec3(:)

    ! v2-v3
    vec1 = xxs(:,2)
    vec2 = xxs(:,3)
    CALL get_vec(vec1,vec2,eta,vec3,num_hex)
    r_eta(1,0,:) = vec3(:)

    ! v6-v7
    vec1 = xxs(:,6)
    vec2 = xxs(:,7)
    CALL get_vec(vec1,vec2,eta,vec3,num_hex)
    r_eta(1,1,:) = vec3(:)

    P1P2(:) = (1.d0-xi)*(1.d0-eta)*r_beta(0,0,:)+(1.d0-xi)*eta*r_beta(0,1,:)+&
    xi*(1.d0-eta)*r_beta(1,0,:)+xi*eta*r_beta(1,1,:)

    P2P3(:) = (1.d0-eta)*(1.d0-beta)*r_xi(0,0,:)+(1.d0-eta)*beta*r_xi(0,1,:)+&
    eta*(1.d0-beta)*r_xi(1,0,:)+eta*beta*r_xi(1,1,:)

    P1P3(:) = (1.d0-xi)*(1.d0-beta)*r_eta(0,0,:)+(1.d0-xi)*beta*r_eta(0,1,:)+&
    xi*(1.d0-beta)*r_eta(1,0,:)+xi*beta*r_eta(1,1,:)

    CALL xyCoor_atGps(8,xxs,xi,eta,beta,P1P2P3)

    ! face xi = 0
    call get_surface_vec(xxs(:,5),xxs(:,1),xxs(:,4),xxs(:,8),&
    1-beta,eta,vec5,num_hex)
    r_xi0(1:3) = vec5(1:3)

    ! face xi = 1
    call get_surface_vec(xxs(:,6),xxs(:,2),xxs(:,3),xxs(:,7),&
    1-beta,eta,vec5,num_hex)
    r_xi1(1:3) = vec5(1:3)

    ! face eta = 0
    call get_surface_vec(xxs(:,5),xxs(:,6),xxs(:,2),xxs(:,1),&
    xi,1-beta,vec5,num_hex)
    r_eta0(1:3) = vec5(1:3)

    ! face eta = 1
    CALL get_surface_vec(xxs(:,8),xxs(:,7),xxs(:,3),xxs(:,4),&
    xi,1-beta,vec5,num_hex)
    r_eta1(1:3) = vec5(1:3)

    ! face beta = 0
    call get_surface_vec(xxs(:,1),xxs(:,2),xxs(:,3),xxs(:,4),&
    xi,eta,vec5,num_hex)
    r_beta0(1:3) = vec5(1:3)

    ! face beta = 1
    call get_surface_vec(xxs(:,5),xxs(:,6),xxs(:,7),xxs(:,8),&
    xi,eta,vec5,num_hex)
    r_beta1(1:3) = vec5(1:3)

    P1(:) = (1.d0-xi)*r_xi0(:) + xi*r_xi1(:)
    P2(:) = (1.d0-eta)*r_eta0(:) + eta*r_eta1(:)
    P3(:) = (1.d0-beta)*r_beta0(:) + beta*r_beta1(:)

    xx(:) = P1(:)+P2(:)+P3(:)-P2P3(:)-P1P2(:)-P1P3(:)+P1P2P3(:)

    END SUBROUTINE Coor_transfinite

    SUBROUTINE Deri_transfinite(xxs,xi,eta,beta,coor_der)
    IMPLICIT NONE
    DOUBLE PRECISION:: xxs(3,8),coor_der(3,3),xi,eta,beta,xi_incre,eta_incre,beta_incre
    DOUBLE PRECISION:: xx_add(3),xx_ori(3)
    double precision,parameter::incre = 1d-7
      
    ! mapd(x/y/z,xi/eta/beta)
    ! eg. mapd(2,3) :: ybeta
      
    CALL Coor_transfinite(xxs,xi,eta,beta,xx_ori)
      
    xi_incre = xi + incre
    CALL Coor_transfinite(xxs,xi_incre,eta,beta,xx_add)
    coor_der(1,1) = (xx_add(1) - xx_ori(1)) / incre
    coor_der(2,1) = (xx_add(2) - xx_ori(2)) / incre
    coor_der(3,1) = (xx_add(3) - xx_ori(3)) / incre
      
    eta_incre = eta + incre
    CALL Coor_transfinite(xxs,xi,eta_incre,beta,xx_add)
    coor_der(1,2) = (xx_add(1) - xx_ori(1)) / incre
    coor_der(2,2) = (xx_add(2) - xx_ori(2)) / incre
    coor_der(3,2) = (xx_add(3) - xx_ori(3)) / incre
      
      
    beta_incre = beta + incre
    CALL Coor_transfinite(xxs,xi,eta,beta_incre,xx_add)
    coor_der(1,3) = (xx_add(1) - xx_ori(1)) / incre
    coor_der(2,3) = (xx_add(2) - xx_ori(2)) / incre
    coor_der(3,3) = (xx_add(3) - xx_ori(3)) / incre
      
    END SUBROUTINE Deri_transfinite
