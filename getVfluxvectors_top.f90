      SUBROUTINE getVfluxvectors_top(Qv,nablaQ,Fv,Gv,Hv,krad,nx,ny,nz)

        ! Qv and nablaQ are all smoothed ones
        ! for the interface flux points

        use setup3d
        IMPLICIT NONE

        double precision,dimension(5), intent(in)       ::Qv
        double precision,dimension(:,:), intent(in)     ::nablaQ(5,3)
        double precision, intent(in) :: nx, ny, nz
        double precision,dimension(5), intent(out)      ::Fv,Gv,Hv
        double precision,dimension(:,:) :: T1(3,3), T2(3,3)
        double precision:: rho,rhou,rhov,rhoE,pr,T,krad
        double precision:: u,v,inte,ux,vx,uy,vy,ex,ey,rhox,rhoy
        double precision:: w,rhow,dw(3),rhoz,uz,vz,ez,wx,wy,wz
        double precision :: du(3),dv(3),drho(3),de(3),dKE(3)
        double precision :: diag,tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        double precision :: tauyx,tauzx,tauzy
        double precision :: gam1
        gam1 = gam -1

        rho  = Qv(1)
        rhou = Qv(2)
        rhov = Qv(3)
        rhow = Qv(4)
        rhoE = Qv(5)

        u = rhou/rho
        v = rhov/rho
        w = rhow/rho
        inte = rhoE/rho - 0.5*(u**2+v**2+w**2)

       drho(1:3) = nablaQ(1,1:3)
       du(1:3) = (nablaQ(2,:)-drho(1:3)*u) / rho
       dv(1:3) = (nablaQ(3,:)-drho(1:3)*v) / rho
       dw(1:3) = (nablaQ(4,:)-drho(1:3)*w) / rho

       dKE(1:3) = 0.5*(u**2+v**2+w**2)*drho(1:3)+rho*(u*du(1:3)+v*dv(1:3)+w*dw(1:3))
       de(1:3) = (nablaQ(5,1:3)-dKE(1:3)-drho(1:3)*inte)/rho

      rhox= drho(1)
      rhoy= drho(2)
      rhoz= drho(3)
!
      ux= du(1)
      uy= du(2)
      uz = du(3)
!
      vx= dv(1)
      vy= dv(2)
      vz = dv(3)
!
      wx= dw(1)
      wy= dw(2)
      wz= dw(3)
!
      ex= de(1)
      ey= de(2)
      ez= de(3)
!


      diag = (ux+vy+wz)/3.0
      tauxx = 2.0*mu*rho*(ux-diag)
      tauxy = mu*rho*(uy+vx) 
      tauxz = mu*rho*(wx+uz) 
      tauyy = 2.0*mu*rho*(vy-diag)
      tauyz = mu*rho*(vz+wy)
      tauzz = 2.0*mu*rho*(wz-diag)


        T1(1,1) = tauxx
        T1(1,2) = tauxy
        T1(1,3) = tauxz
        T1(2,1) = T1(1,2)
        T1(2,2) = tauyy
        T1(2,3) = tauyz
        T1(3,1) = T1(1,3)
        T1(3,2) = T1(2,3)
        T1(3,3) = tauzz

! stress-free boundary condition (shear stress vanish on boundary)

        call transform_stress(nx,ny,nz,T1,T2)

        tauxx = T2(1,1)
        tauxy = T2(1,2)
        tauxz = T2(1,3)
        tauyx = T2(2,1)
        tauyy = T2(2,2)
        tauyz = T2(2,3)
        tauzx = T2(3,1)
        tauzy = T2(3,2)
        tauzz = T2(3,3)

! compute flux vector Fv
        Fv(1) = 0.d0
        Fv(2) = tauxx
        Fv(3) = tauxy 
        Fv(4) = tauxz 
       Fv(5) = kappa_s*(rho*ex-gam1*inte*rhox)  + u*tauxx+v*tauxy+w*tauxz + krad*rho*ex*gam
!        Fv(5) = u*tauxx+v*tauxy+w*tauxz+nu/prandt*(rho*ex-gam1*inte*rhox)*0.0
!        Fv(5) = u*tauxx+v*tauxy+w*tauxz
! compute flux vector Gv
        Gv(1) = 0.d0
        Gv(2) = tauyx
        Gv(3) = tauyy 
        Gv(4) = tauyz 
        Gv(5) = kappa_s*(rho*ey-gam1*inte*rhoy)  + u*tauyx+v*tauyy+w*tauyz + krad*rho*ey*gam
!        Gv(5) = u*tauxy+v*tauyy+w*tauyz+nu/prandt*(rho*ey-gam1*inte*rhoy)*0.0
!       Gv(5) = u*tauxy+v*tauyy+w*tauyz
!  compute flux vector Hv
        Hv(1) = 0.d0
        Hv(2) = tauzx
        Hv(3) = tauzy 
        Hv(4) = tauzz 
!  added grad P
        Hv(5) = kappa_s*(rho*ez-gam1*inte*rhoz)  + u*tauzx+v*tauzy+w*tauzz + krad*rho*ez*gam
!        Hv(5) = u*tauxz+v*tauyz+w*tauzz+nu/prandt*(rho*ez-gam1*inte*rhoz)*0.0
!        Hv(5) = u*tauxz+v*tauyz+w*tauzz
!  add rho here in front of luminosity-- 

      END SUBROUTINE getVfluxvectors_top


      SUBROUTINE transform_stress(nx,ny,nz,T1,T2)      

        double precision, dimension(:,:), intent(in) :: T1(3,3)
        double precision, dimension(:,:), intent(out):: T2(3,3)
        double precision, dimension(:,:) :: A(3,3), AT(3,3)
        double precision, dimension(:,:) :: temp_01(3,3), temp_02(3,3)

        double precision :: nx, ny, nz
        double precision :: theta, phi
        integer :: i,j,k

        temp_01(:,:) = 0.d0
        temp_02(:,:) = 0.d0
        T2(:,:) = 0.d0


        phi = datan2(ny,nx)
        theta = dacos(nz)

        A(1,1) = dsin(theta)*dcos(phi)
        A(1,2) = dsin(theta)*dsin(phi)
        A(1,3) = dcos(theta)
        A(2,1) = dcos(theta)*dcos(phi)
        A(2,2) = dcos(theta)*dsin(phi)
        A(2,3) = -dsin(theta)
        A(3,1) = -dsin(phi)
        A(3,2) = dcos(phi)
        A(3,3) = 0.d0


        AT(1,1) = dsin(theta)*dcos(phi)
        AT(1,2) = dcos(theta)*dcos(phi)
        AT(1,3) = -dsin(phi)
        AT(2,1) = dsin(theta)*dsin(phi)
        AT(2,2) = dcos(theta)*dsin(phi)
        AT(2,3) = dcos(phi)
        AT(3,1) = dcos(theta)
        AT(3,2) = -dsin(theta)
        AT(3,3) = 0.d0


        do i=1,3
          do j=1,3
            do k =1,3
                temp_01(i,j) = temp_01(i,j) + A(i,k)*T1(k,j)
            end do
          end do
        end do

        do i=1,3
          do j=1,3
            do k=1,3
                temp_02(i,j) = temp_02(i,j) + temp_01(i,k)*AT(k,j)
            end do
          end do
        end do


! stress-free boundary (shear stress vanish on boundary)

        temp_02(1,2) = 0.d0
        temp_02(1,3) = 0.d0
        temp_02(2,1) = 0.d0
        temp_02(3,1) = 0.d0

! reconstruct viscous stress tensor

        temp_01(:,:) = 0.d0

        do i=1,3
          do j=1,3
            do k =1,3
                temp_01(i,j) = temp_01(i,j) + AT(i,k)*temp_02(k,j)
            end do
          end do
        end do

        do i=1,3
          do j=1,3
            do k=1,3
                T2(i,j) = T2(i,j) + temp_01(i,k)*A(k,j)
            end do
          end do
        end do



      END SUBROUTINE transform_stress
