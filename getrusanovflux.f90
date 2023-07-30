!
        SUBROUTINE getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
       
        use setup3d 
        IMPLICIT NONE

        include 'mpif.h'

        double precision,dimension(5), intent(in) :: Qfl,Qfr
        double precision,dimension(3), intent(in) :: normf
        double precision,dimension(5), intent(out):: Fnl,Fnr
        double precision, intent(out) :: eigv
        double precision :: rhoi,ul,vl,ur,vr,Vn_l,Vn_r,pl,pr,c_a,magnorm
        double precision :: wr,wl
        integer, intent(in) :: sign_l, sign_r

        ! for the left side
        rhoi= 1./Qfl(1)
        ul      = rhoi * Qfl(2)
        vl      = rhoi * Qfl(3)
        wl      = rhoi * Qfl(4)
        pl      = (gam-1)*( Qfl(5) - 0.5 * Qfl(1) * (ul**2+vl**2+wl**2) )
        Vn_l    = sign_l * (ul*normf(1)+vl*normf(2)+wl*normf(3))

        ! for the right side
        rhoi= 1./Qfr(1)
        ur      = rhoi * Qfr(2)
        vr      = rhoi * Qfr(3)
        wr      = rhoi * Qfr(4)
        pr      = (gam-1)*( Qfr(5) - 0.5 * Qfr(1) * (ur**2+vr**2+wr**2) )
        Vn_r    = sign_l * (ur*normf(1)+vr*normf(2)+wr*normf(3))

        ! Calculate the normal flux component
        Fnl(1) = Qfl(1)*Vn_l + Qfr(1)*Vn_r
        Fnl(2) = Qfl(2)*Vn_l + Qfr(2)*Vn_r + (pl+pr)*normf(1)*sign_l
        Fnl(3) = Qfl(3)*Vn_l + Qfr(3)*Vn_r + (pl+pr)*normf(2)*sign_l
        Fnl(4) = Qfl(4)*Vn_l + Qfr(4)*Vn_r + (pl+pr)*normf(3)*sign_l
        Fnl(5) = ( Qfl(5)+pl )*Vn_l + ( Qfr(5)+pr )*Vn_r

        magnorm= sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
        c_a    = sqrt(gam*(pl+pr)/(Qfl(1)+Qfr(1)))
        eigv   = 0.5*abs(Vn_l+Vn_r) + c_a*magnorm

        Fnl = (Fnl - (Qfr-Qfl)*eigv)*0.5*sign_l
        Fnr = sign_l*sign_r*Fnl

        END SUBROUTINE getrusanovflux




