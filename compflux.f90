    SUBROUTINE compflux
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: ic,ifp,jfp,kfp,sp
    double precision,dimension(5) :: Qf,Qs,F,G,H,Qv,Fv,Gv,Hv,Qf2,Qf3
    double precision :: nablaQ(5,3),ktemp,magnorm
    double precision :: normf(3),nx,ny,nz
    integer :: mfp,nfp,ifn,sign_r,iface,ifpl,jfpl,kfpl,faml
    
    double precision :: kappai,kappao

    kappai = constL/(4.d0*pi*RINLET**2.d0*constheati)

    do ic = 1,NCELL

        do kfp = 1,N
        do jfp = 1,N
        do ifp = 2,N
            Qf(1:5) = 0.d0
            do sp = 1,N
                Qs(1:5) = Q(1:5,sp,jfp,kfp,ic)
                Qf(1:5) = Qf(1:5) + Qs(1:5)*Lmat(ifp,sp)
            end do
            if (vismode==1) Qvfi(1:5,ifp,jfp,kfp,ic)=Qf(1:5)
            CALL getfluxvectors(Qf,F,G,H)
            F1(1:5,ifp,jfp,kfp,ic) = F(1:5)*S1(1,1,ifp,jfp,kfp,ic) &
                                     + G(1:5)*S1(1,2,ifp,jfp,kfp,ic) &
                                     + H(1:5)*S1(1,3,ifp,jfp,kfp,ic)
        end do
        end do
        end do

        do kfp = 1,N
        do jfp = 2,N
        do ifp = 1,N
            Qf(1:5) = 0.d0
            do sp = 1,N
                Qs(1:5) = Q(1:5,ifp,sp,kfp,ic)
                Qf(1:5) = Qf(1:5) + Qs(1:5)*Lmat(jfp,sp)
            end do
            if (vismode==1) Qvfj(1:5,ifp,jfp,kfp,ic)=Qf(1:5)
            call getfluxvectors(Qf,F,G,H)
            G2(1:5,ifp,jfp,kfp,ic) = F(1:5)*S2(2,1,ifp,jfp,kfp,ic) &
                                     + G(1:5)*S2(2,2,ifp,jfp,kfp,ic) &
                                     + H(1:5)*S2(2,3,ifp,jfp,kfp,ic)
        end do
        end do
        end do

        do kfp = 2,N
        do jfp = 1,N
        do ifp = 1,N
            Qf(1:5) = 0.d0
            do sp=1,N
                Qs(1:5) = Q(1:5,ifp,jfp,sp,ic)
                Qf(1:5) = Qf(1:5) + Qs(1:5)*Lmat(kfp,sp)
            end do
            if (vismode==1) Qvfk(1:5,ifp,jfp,kfp,ic)=Qf(1:5)
            call getfluxvectors(Qf,F,G,H)
            H3(1:5,ifp,jfp,kfp,ic) = F(1:5)*S3(3,1,ifp,jfp,kfp,ic) &
                                     + G(1:5)*S3(3,2,ifp,jfp,kfp,ic) &
                                     + H(1:5)*S3(3,3,ifp,jfp,kfp,ic)
        end do
        end do
        end do

    end do

    call interfaceflux

    call procintflux
    call BCflux

    if(vismode==1) then
        
        call nablaQsp

        do ic=1,NCELL
            ! Family i
            do kfp=1,N
            do jfp=1,N
            do ifp=2,N ! for interior flux points only
                Qf(1:5)  = 0.d0
                Qf2(1:5) = 0.d0
                Qf3(1:5) = 0.d0

                do sp=1,N
                    Qf(1:5) = Qf(1:5) + nablaQs(1:5,1,sp,jfp,kfp,ic)*Lmat(ifp,sp)
                    Qf2(1:5) = Qf2(1:5) + nablaQs(1:5,2,sp,jfp,kfp,ic)*Lmat(ifp,sp)
                    Qf3(1:5) = Qf3(1:5) + nablaQs(1:5,3,sp,jfp,kfp,ic)*Lmat(ifp,sp)
                end do

                nablaQvfi(1:5,1,ifp,jfp,kfp,ic) = Qf(1:5)
                nablaQvfi(1:5,2,ifp,jfp,kfp,ic) = Qf2(1:5)
                nablaQvfi(1:5,3,ifp,jfp,kfp,ic) = Qf3(1:5)

            end do
            end do
            end do

            ! Family j
            do kfp=1,N
            do jfp=2,N
            do ifp=1,N
                Qf(1:5)  = 0.d0
                Qf2(1:5) = 0.d0
                Qf3(1:5) = 0.d0

                do sp=1,N
                    Qf(1:5)  = Qf(1:5) + nablaQs(1:5,1,ifp,sp,kfp,ic)*Lmat(jfp,sp)
                    Qf2(1:5) = Qf2(1:5) + nablaQs(1:5,2,ifp,sp,kfp,ic)*Lmat(jfp,sp)
                    Qf3(1:5) = Qf3(1:5) + nablaQs(1:5,3,ifp,sp,kfp,ic)*Lmat(jfp,sp)
                end do

                nablaQvfj(1:5,1,ifp,jfp,kfp,ic) = Qf(1:5)
                nablaQvfj(1:5,2,ifp,jfp,kfp,ic) = Qf2(1:5)
                nablaQvfj(1:5,3,ifp,jfp,kfp,ic) = Qf3(1:5)

            end do
            end do
            end do

            ! Family k
            do kfp=2,N
            do jfp=1,N
            do ifp=1,N
                Qf(1:5) = 0.d0
                Qf2(1:5) = 0.d0
                Qf3(1:5) = 0.d0
                do sp=1,N
                    Qf(1:5) = Qf(1:5) + nablaQs(1:5,1,ifp,jfp,sp,ic)*Lmat(kfp,sp)
                    Qf2(1:5) = Qf2(1:5) + nablaQs(1:5,2,ifp,jfp,sp,ic)*Lmat(kfp,sp)
                    Qf3(1:5) = Qf3(1:5) + nablaQs(1:5,3,ifp,jfp,sp,ic)*Lmat(kfp,sp)
                end do
                nablaQvfk(1:5,1,ifp,jfp,kfp,ic) = Qf(1:5)
                nablaQvfk(1:5,2,ifp,jfp,kfp,ic) = Qf2(1:5)
                nablaQvfk(1:5,3,ifp,jfp,kfp,ic) = Qf3(1:5)
            end do
            end do
            end do
        end do


        call interfaceVISflux 
        call procintVISflux   
        call BCvisflux

        do ic=1,NCELL
           do kfp=1,N
           do jfp=1,N
           do ifp=1,N+1
              Qv(1:5) = Qvfi(1:5,ifp,jfp,kfp,ic)
              nablaQ(1:5,1) = nablaQvfi(1:5,1,ifp,jfp,kfp,ic)
              nablaQ(1:5,2) = nablaQvfi(1:5,2,ifp,jfp,kfp,ic)
              nablaQ(1:5,3) = nablaQvfi(1:5,3,ifp,jfp,kfp,ic)
              ktemp = radCfi(ifp,jfp,kfp,ic)
              call getVfluxvectors(Qv,nablaQ,Fv,Gv,Hv,ktemp)
              Fv1(1:5,ifp,jfp,kfp,ic) = Fv(1:5)*S1(1,1,ifp,jfp,kfp,ic)&
                                      + Gv(1:5)*S1(1,2,ifp,jfp,kfp,ic) &
                                      + Hv(1:5)*S1(1,3,ifp,jfp,kfp,ic)
           
           end do
           end do
           end do

           do kfp=1,N
           do jfp=1,N+1
           do ifp=1,N
              Qv(1:5) = Qvfj(1:5,ifp,jfp,kfp,ic)
              nablaQ(1:5,1) = nablaQvfj(1:5,1,ifp,jfp,kfp,ic)
              nablaQ(1:5,2) = nablaQvfj(1:5,2,ifp,jfp,kfp,ic)
              nablaQ(1:5,3) = nablaQvfj(1:5,3,ifp,jfp,kfp,ic)
              ktemp = radCfj(ifp,jfp,kfp,ic)
              call getVfluxvectors(Qv,nablaQ,Fv,Gv,Hv,ktemp)
              Gv2(1:5,ifp,jfp,kfp,ic) = Fv(1:5)*S2(2,1,ifp,jfp,kfp,ic)&
                                      + Gv(1:5)*S2(2,2,ifp,jfp,kfp,ic) &
                                      + Hv(1:5)*S2(2,3,ifp,jfp,kfp,ic)
           end do
           end do
           end do

           do kfp=1,N+1
           do jfp=1,N
           do ifp=1,N
              Qv(1:5) = Qvfk(1:5,ifp,jfp,kfp,ic)
              nablaQ(1:5,1) = nablaQvfk(1:5,1,ifp,jfp,kfp,ic)
              nablaQ(1:5,2) = nablaQvfk(1:5,2,ifp,jfp,kfp,ic)
              nablaQ(1:5,3) = nablaQvfk(1:5,3,ifp,jfp,kfp,ic)
              ktemp = radCfk(ifp,jfp,kfp,ic)
              call getVfluxvectors(Qv,nablaQ,Fv,Gv,Hv,ktemp)
              Hv3(1:5,ifp,jfp,kfp,ic) = Fv(1:5)*S3(3,1,ifp,jfp,kfp,ic)&
                                     + Gv(1:5)*S3(3,2,ifp,jfp,kfp,ic) &
                                     + Hv(1:5)*S3(3,3,ifp,jfp,kfp,ic)
           end do
           end do
           end do
        end do
    end if

    DO ifn=1,NINLET
        iface = IBFINL(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp = 1,N  
             faml = Gfp2Lfp(3,mfp,nfp,iface) 

             ifpl = Gfp2Lfp(5,mfp,nfp,iface) 
             jfpl = Gfp2Lfp(6,mfp,nfp,iface) 
	         kfpl = Gfp2Lfp(7,mfp,nfp,iface) 

                sign_r = 1

                if (faml == 1) then
                     if (ifpl == 1) sign_r = -1
                      normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)


                elseif(faml ==2) then
                     if (jfpl == 1) sign_r = -1
                      normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)


                elseif(faml==3) then
                     if (kfpl == 1) sign_r = -1
                      normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)

                end if

                magnorm =sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
!   pointing outward
! then point r direction

                nx = -normf(1)/magnorm*sign_r
                ny = -normf(2)/magnorm*sign_r
                nz = -normf(3)/magnorm*sign_r

                if(faml==1) then
                  Qv(1:5) = Qvfi(1:5,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,1) = nablaQvfi(1:5,1,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,2) = nablaQvfi(1:5,2,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,3) = nablaQvfi(1:5,3,ifpl,jfpl,kfpl,ic)
                  call getVfluxvectors_bot(Qv,nablaQ,Fv,Gv,Hv,kappai,nx,ny,nz,constheati)
                  Fv1(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S1(1,1,ifpl,jfpl,kfpl,ic)&
                                          + Gv(1:5)*S1(1,2,ifpl,jfpl,kfpl,ic) &
                                          + Hv(1:5)*S1(1,3,ifpl,jfpl,kfpl,ic)

                else if(faml==2) then
                  Qv(1:5) = Qvfj(1:5,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,1) = nablaQvfj(1:5,1,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,2) = nablaQvfj(1:5,2,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,3) = nablaQvfj(1:5,3,ifpl,jfpl,kfpl,ic)
                  call getVfluxvectors_bot(Qv,nablaQ,Fv,Gv,Hv,kappai,nx,ny,nz,constheati)
                  Gv2(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S3(3,1,ifpl,jfpl,kfpl,ic)&
                                          + Gv(1:5)*S3(3,2,ifpl,jfpl,kfpl,ic) &
                                          + Hv(1:5)*S3(3,3,ifpl,jfpl,kfpl,ic)
                else
                  Qv(1:5) = Qvfk(1:5,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,1) = nablaQvfk(1:5,1,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,2) = nablaQvfk(1:5,2,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,3) = nablaQvfk(1:5,3,ifpl,jfpl,kfpl,ic)
                  call getVfluxvectors_bot(Qv,nablaQ,Fv,Gv,Hv,kappai,nx,ny,nz,constheati)
                  Hv3(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S3(3,1,ifpl,jfpl,kfpl,ic)&
                                          + Gv(1:5)*S3(3,2,ifpl,jfpl,kfpl,ic) &
                                          + Hv(1:5)*S3(3,3,ifpl,jfpl,kfpl,ic)
                end if
        end do
        end do

       END DO


        DO ifn=1,NOUTLET 
        iface = IBFOUT(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp = 1,N
            faml = Gfp2Lfp(3,mfp,nfp,iface) 
            ifpl = Gfp2Lfp(5,mfp,nfp,iface) 
            jfpl = Gfp2Lfp(6,mfp,nfp,iface) 
            kfpl = Gfp2Lfp(7,mfp,nfp,iface) 

            sign_r = 1

            if (faml == 1) then
                if (ifpl == 1) sign_r = -1
                normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
            elseif(faml ==2) then
                if (jfpl == 1) sign_r = -1
                normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)
            elseif(faml==3) then
                if (kfpl == 1) sign_r = -1
                normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)
            end if

            magnorm =sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
!   pointing outward

            nx = normf(1)/magnorm*sign_r
            ny = normf(2)/magnorm*sign_r
            nz = normf(3)/magnorm*sign_r

            if(faml==1) then
                Qv(1:5) = Qvfi(1:5,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,1) = nablaQvfi(1:5,1,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,2) = nablaQvfi(1:5,2,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,3) = nablaQvfi(1:5,3,ifpl,jfpl,kfpl,ic)
                kappao = radCfi(ifpl,jfpl,kfpl,ic)

                call getVfluxvectors_top(Qv,nablaQ,Fv,GV,Hv,kappao,nx,ny,nz)
                Fv1(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S1(1,1,ifpl,jfpl,kfpl,ic)&
                                             + Gv(1:5)*S1(1,2,ifpl,jfpl,kfpl,ic) &
                                             + Hv(1:5)*S1(1,3,ifpl,jfpl,kfpl,ic)

            else if(faml==2) then
                Qv(1:5) = Qvfj(1:5,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,1) = nablaQvfj(1:5,1,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,2) = nablaQvfj(1:5,2,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,3) = nablaQvfj(1:5,3,ifpl,jfpl,kfpl,ic)
                kappao = radCfj(ifpl,jfpl,kfpl,ic)

                call getVfluxvectors_top(Qv,nablaQ,Fv,GV,Hv,kappao,nx,ny,nz)
                Gv2(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S2(2,1,ifpl,jfpl,kfpl,ic)&
                                             + Gv(1:5)*S2(2,2,ifpl,jfpl,kfpl,ic) &
                                             + Hv(1:5)*S2(2,3,ifpl,jfpl,kfpl,ic)
                else
                  Qv(1:5) = Qvfk(1:5,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,1) = nablaQvfk(1:5,1,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,2) = nablaQvfk(1:5,2,ifpl,jfpl,kfpl,ic)
                  nablaQ(1:5,3) = nablaQvfk(1:5,3,ifpl,jfpl,kfpl,ic)
                  kappao = radCfk(ifpl,jfpl,kfpl,ic)

                  call getVfluxvectors_top(Qv,nablaQ,Fv,Gv,Hv,kappao,nx,ny,nz) 
                  Hv3(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S3(3,1,ifpl,jfpl,kfpl,ic)&
                                               + Gv(1:5)*S3(3,2,ifpl,jfpl,kfpl,ic) &
                                               + Hv(1:5)*S3(3,3,ifpl,jfpl,kfpl,ic)
                end if
        end do
        end do

    END DO


    DO ifn=1,NOUTLET 
        iface = IBFOUT(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp = 1,N
            faml = Gfp2Lfp(3,mfp,nfp,iface) 
            ifpl = Gfp2Lfp(5,mfp,nfp,iface) 
            jfpl = Gfp2Lfp(6,mfp,nfp,iface) 
            kfpl = Gfp2Lfp(7,mfp,nfp,iface) 

            sign_r = 1

            if (faml == 1) then
                if (ifpl == 1) sign_r = -1
                normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
            elseif(faml ==2) then
                if (jfpl == 1) sign_r = -1
                normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)
            elseif(faml==3) then
                if (kfpl == 1) sign_r = -1
                normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)
            end if

            magnorm =sqrt(normf(1)**2+normf(2)**2+normf(3)**2)
!   pointing outward

            nx = normf(1)/magnorm*sign_r
            ny = normf(2)/magnorm*sign_r
            nz = normf(3)/magnorm*sign_r

            if(faml==1) then
                Qv(1:5) = Qvfi(1:5,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,1) = nablaQvfi(1:5,1,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,2) = nablaQvfi(1:5,2,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,3) = nablaQvfi(1:5,3,ifpl,jfpl,kfpl,ic)
                kappao = radCfi(ifpl,jfpl,kfpl,ic)

                call getVfluxvectors_top(Qv,nablaQ,Fv,GV,Hv,kappao,nx,ny,nz)
                Fv1(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S1(1,1,ifpl,jfpl,kfpl,ic)&
                                             + Gv(1:5)*S1(1,2,ifpl,jfpl,kfpl,ic) &
                                             + Hv(1:5)*S1(1,3,ifpl,jfpl,kfpl,ic)

            else if(faml==2) then
                Qv(1:5) = Qvfj(1:5,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,1) = nablaQvfj(1:5,1,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,2) = nablaQvfj(1:5,2,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,3) = nablaQvfj(1:5,3,ifpl,jfpl,kfpl,ic)
                kappao = radCfj(ifpl,jfpl,kfpl,ic)

                call getVfluxvectors_top(Qv,nablaQ,Fv,GV,Hv,kappao,nx,ny,nz)
                Gv2(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S2(2,1,ifpl,jfpl,kfpl,ic)&
                                             + Gv(1:5)*S2(2,2,ifpl,jfpl,kfpl,ic) &
                                             + Hv(1:5)*S2(2,3,ifpl,jfpl,kfpl,ic)
            else
                Qv(1:5) = Qvfk(1:5,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,1) = nablaQvfk(1:5,1,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,2) = nablaQvfk(1:5,2,ifpl,jfpl,kfpl,ic)
                nablaQ(1:5,3) = nablaQvfk(1:5,3,ifpl,jfpl,kfpl,ic)
                kappao = radCfk(ifpl,jfpl,kfpl,ic)

                call getVfluxvectors_top(Qv,nablaQ,Fv,Gv,Hv,kappao,nx,ny,nz) 
                Hv3(1:5,ifpl,jfpl,kfpl,ic) = Fv(1:5)*S3(3,1,ifpl,jfpl,kfpl,ic)&
                                             + Gv(1:5)*S3(3,2,ifpl,jfpl,kfpl,ic) &
                                             + Hv(1:5)*S3(3,3,ifpl,jfpl,kfpl,ic)
            end if
        end do
        end do

    END DO

    END SUBROUTINE compflux
