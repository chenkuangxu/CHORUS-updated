      SUBROUTINE getVfluxvectors(Qv,nablaQ,Fv,Gv,Hv,krad)

        ! Qv and nablaQ are all smoothed ones
        ! for the interface flux points

        use setup3d
        IMPLICIT NONE
        include 'mpif.h'

        double precision,dimension(5), intent(in)       ::Qv
        double precision,dimension(:,:), intent(in)     ::nablaQ(5,3)
        double precision,dimension(5), intent(out)      ::Fv,Gv,Hv
        double precision:: rho,rhou,rhov,rhoE,pr,T,krad
        double precision:: u,v,inte,ux,vx,uy,vy,ex,ey,rhox,rhoy
        double precision:: w,rhow,dw(3),rhoz,uz,vz,ez,wx,wy,wz
        double precision :: du(3),dv(3),drho(3),de(3),dKE(3)
        double precision :: diag,tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
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

! compute flux vector Fv
        Fv(1) = 0.d0
        Fv(2) = tauxx
        Fv(3) = tauxy 
        Fv(4) = tauxz 
        Fv(5) = u*tauxx+v*tauxy+w*tauxz+kappa_s*(rho*ex-gam1*inte*rhox) + krad*rho*cp*(ex*gam/cp) 
        ! Fv(5) = rho*ex
!        Fv(5) = u*tauxx+v*tauxy+w*tauxz+nu/prandt*(rho*ex-gam1*inte*rhox)*0.0
!        Fv(5) = u*tauxx+v*tauxy+w*tauxz
! compute flux vector Gv
        Gv(1) = 0.d0
        Gv(2) = tauxy
        Gv(3) = tauyy 
        Gv(4) = tauyz 
        Gv(5) = u*tauxy+v*tauyy+w*tauyz+kappa_s*(rho*ey-gam1*inte*rhoy) + krad*rho*cp*(ey*gam/cp)
        ! Gv(5) = u*tauxy+v*tauyy+w*tauyz+krad*rho*cp*(ey*gam/cp)
!        Gv(5) = u*tauxy+v*tauyy+w*tauyz+nu/prandt*(rho*ey-gam1*inte*rhoy)*0.0
!       Gv(5) = u*tauxy+v*tauyy+w*tauyz
!  compute flux vector Hv
        Hv(1) = 0.d0
        Hv(2) = tauxz
        Hv(3) = tauyz 
        Hv(4) = tauzz 
!  added grad P
        Hv(5) = u*tauxz+v*tauyz+w*tauzz+kappa_s*(rho*ez-gam1*inte*rhoz) + krad*rho*cp*(ez*gam/cp)
        ! Hv(5) = u*tauxz+v*tauyz+w*tauzz+krad*rho*cp*(ez*gam/cp)
!        Hv(5) = u*tauxz+v*tauyz+w*tauzz+nu/prandt*(rho*ez-gam1*inte*rhoz)*0.0
!        Hv(5) = u*tauxz+v*tauyz+w*tauzz
!  add rho here in front of luminosity-- 

      END SUBROUTINE getVfluxvectors


