	! THIS IS THE SUBR TO SET THE BOUNNDARY CONDITIONS AND CAALCULATE FLUXES AT BOUNDARY FACES

	SUBROUTINE BCflux
	use setup3d
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	integer		::	faml,mfp,kfpl
	integer		::	nfp,ifpl,jfpl,ic,iface,sp,ifn,sign_l,sign_r,k
	double precision, dimension(5)	:: Qs,Qfl,Qfr,Fnl,Fnr,Qfl2,Qfr2
	double precision		:: normf(3),VN, eigv,magnorm
	double precision		:: V2, u_exit, v_exit,w_exit, rho_exit, gam1, T_in, p_in, rho_in
	double precision		:: u_l, v_l, w_l,rho_l, p_l, rho_out
	double precision		:: uw, vw,ww, rhow, Tw, pw, energw
	double precision		:: ur, vr, wr,rhor, pr,radist,XFP,YFP,ratio,ZFP
        double precision :: rrx(3), xxf(3,8)
        integer :: VF1,VF2,iv
        gam1 = gam -1.0 

        sign_r = 1


   ! ********************************************************
   ! Wall Faces, now INLET BC with constat Temp flux
   ! ********************************************************

    DO ifn=1,NINLET
        iface = IBFINL(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp=1,N
                faml = Gfp2Lfp(3,mfp,nfp,iface)

                ifpl = Gfp2Lfp(5,mfp,nfp,iface)
                jfpl = Gfp2Lfp(6,mfp,nfp,iface)
                kfpl = Gfp2Lfp(7,mfp,nfp,iface)

                sign_l = 1
                if(faml==1) then
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,sp,jfpl,kfpl,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(ifpl,sp)
                        end do
                        if(ifpl==1) sign_l=-1
                        normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
                else if(faml==2) then
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,ifpl,sp,kfpl,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(jfpl,sp)
                        end do
                        if(jfpl==1) sign_l=-1
                        normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)
                else
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,ifpl,jfpl,sp,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(kfpl,sp)
                        end do
                        if(kfpl==1) sign_l=-1
                        normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)
                end if


                rhow = Qfl(1)
                uw = Qfl(2)/rhow
                vw = Qfl(3)/rhow
                ww = Qfl(4)/rhow

                Qfr(1) = rhow
                magnorm =sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
                VN = (uw*normf(1)+vw*normf(2)+ww*normf(3))/magnorm
                Qfr(2) = Qfr(1) * (uw-VN*normf(1)/magnorm)
                Qfr(3) = Qfr(1) * (vw-VN*normf(2)/magnorm)
                Qfr(4) = Qfr(1) * (ww-VN*normf(3)/magnorm)
                pw = gam1*(Qfl(5) - 0.5*rhow*(uw**2+vw**2+ww**2)) 
                Qfr(5) = pw/gam1+0.5/Qfr(1)*(Qfr(2)**2+Qfr(3)**2+Qfr(4)**2)



                   if(faml==1) then
                      Qvfi(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   else if(faml==2) then
                      Qvfj(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   else
                      Qvfk(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   end if
        

          !!       CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
                 CALL getflux_exact(Qfl,Fnl,normf,sign_l)

                 if(faml==1) then
                        F1(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                else if(faml==2) then
                        G2(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                else
                        H3(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                end if
        end do
        end do

    END DO  !



   ! ********************************************************
   ! Wall Faces, now outlET BC with constat entropy flux
   ! ********************************************************

    DO ifn=1,NOUTLET
        iface = IBFOUT(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp=1,N
                faml = Gfp2Lfp(3,mfp,nfp,iface)

                ifpl = Gfp2Lfp(5,mfp,nfp,iface)
                jfpl = Gfp2Lfp(6,mfp,nfp,iface)
                kfpl = Gfp2Lfp(7,mfp,nfp,iface)

                sign_l = 1
                if(faml==1) then
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,sp,jfpl,kfpl,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(ifpl,sp)
                        end do
                        if(ifpl==1) sign_l=-1
                        normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
                else if(faml==2) then
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,ifpl,sp,kfpl,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(jfpl,sp)
                        end do
                        if(jfpl==1) sign_l=-1
                        normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)
                else
                        Qfl(1:5) = 0.d0
                        do sp=1,N
                                Qs(1:5) = Q(1:5,ifpl,jfpl,sp,ic)
                                Qfl(1:5) = Qfl(1:5) + Qs(1:5)*Lmat(kfpl,sp)
                        end do
                        if(kfpl==1) sign_l=-1
                        normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)
                end if




                rho_l = Qfl(1)
                u_l = Qfl(2)/rho_l
                v_l = Qfl(3)/rho_l
                w_l = Qfl(4)/rho_l
                p_l = (gam-1)*(Qfl(5)-0.5*rho_l*(u_l**2+v_l**2+w_l**2))

! constant temperature wall is used here
                Qfr(1) = rho_l
                energw = pinf/(rinf*(gam-1))
                pr = Qfr(1)*energw*(gam-1)

                magnorm =sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
                VN = (u_l*normf(1)+v_l*normf(2)+w_l*normf(3))/magnorm

                ur = u_l - VN*normf(1)/magnorm
                vr = v_l - VN*normf(2)/magnorm
                wr = w_l - VN*normf(3)/magnorm
                
                Qfr(2) = Qfr(1)*ur
                Qfr(3) = Qfr(1)*vr
                Qfr(4) = Qfr(1)*wr
                Qfr(5) = Qfr(1)*energw + 0.5*Qfr(1)*(ur**2+vr**2+wr**2)
              !  Qfr(5) = p_l/gam1 + 0.5*Qfr(1)*(ur**2+vr**2+wr**2)


                   if(faml==1) then
                      Qvfi(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   else if(faml==2) then
                      Qvfj(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   else
                      Qvfk(1:5,ifpl,jfpl,kfpl,ic)  = Qfr(1:5)
                   end if
        

              !!   CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
                CALL getflux_exact(Qfl,Fnl,normf,sign_l)

                 if(faml==1) then
                        F1(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                else if(faml==2) then
                        G2(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                else
                        H3(1:5,ifpl,jfpl,kfpl,ic)  = Fnl(1:5)
                end if
        end do
        end do

    END DO  !


        END SUBROUTINE BCflux


        SUBROUTINE getflux_exact(Qfl,Fnl,normf,sign_l)
        use setup3d
        IMPLICIT  NONE
        include 'mpif.h'

        double precision, dimension(5), intent(in) :: Qfl
        double precision, dimension(3), intent(in) :: normf
        double precision, dimension(5), intent(out) :: Fnl
        double precision :: rhoi, ul, vl, wl,  pl
        integer, intent(in) :: sign_l

        ! for the left side
        rhoi = 1./Qfl(1)
        ul = rhoi * Qfl(2)
        vl = rhoi * Qfl(3)
        wl = rhoi * Qfl(4)
        pl = (gam-1.0)*(Qfl(5) - 0.5*Qfl(1)*(ul**2+vl**2+wl**2))


        ! Calculate the normal flux component
        Fnl(1) = 0.d0
        Fnl(2) = pl*normf(1)*sign_l
        Fnl(3) = pl*normf(2)*sign_l
        Fnl(4) = pl*normf(3)*sign_l
        Fnl(5) = 0.d0

        Fnl = Fnl*sign_l


        END SUBROUTINE getflux_exact


