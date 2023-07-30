    SUBROUTINE interfaceflux

    use setup3d 

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    integer :: icleft,icright,ifacelc,ifacerc,faml,famr
    integer :: ifpl,jfpl,ifpr,jfpr,ic,iface,sp,ifinter,sign_l,sign_r,k
    integer :: kfpl,kfpr,nfp,mfp,count
        
    double precision, dimension(numv)  :: Qs,Qfl,Qfr,Fnl,Fnr
    double precision, dimension(3)  :: normf
    double precision                :: eigv

    ! calculate the interface fluxes at all interior faces

    do ifinter=1,NINTER
        iface   = IBFINTER(ifinter)
        icleft  = IF2C(iface,1)
        icright = IF2C(iface,2)
            
        do nfp=1,N
        do mfp=1,N
            faml = Gfp2Lfp(3,mfp,nfp,iface)
            famr = Gfp2Lfp(4,mfp,nfp,iface)
            ifpl = Gfp2Lfp(5,mfp,nfp,iface)
            jfpl = Gfp2Lfp(6,mfp,nfp,iface)
            kfpl = Gfp2Lfp(7,mfp,nfp,iface)

            ic = icleft

            sign_l = 1
            Qfl(1:numv) = 0.0d0

            if(faml==1) then

                Qfl(1:numv) = 0.0d0
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,sp,jfpl,kfpl,ic)
                    Qfl(1:numv) = Qfl(1:numv) + Qs(1:numv)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
 
            else if(faml==2) then
  
                Qfl(1:numv) = 0.0d0
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpl,sp,kfpl,ic)
                    Qfl(1:numv) = Qfl(1:numv) + Qs(1:numv)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)

            else

                Qfl(1:numv) = 0.0d0
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpl,jfpl,sp,ic)
                    Qfl(1:numv) = Qfl(1:numv) + Qs(1:numv)*Lmat(kfpl,sp)
                end do
                if(kfpl==1) sign_l=-1
                normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)

            end if

            ! lets get the right solution vector at flux point
            ifpr = Gfp2Lfp(8,mfp,nfp,iface)
            jfpr = Gfp2Lfp(9,mfp,nfp,iface)
            kfpr = Gfp2Lfp(10,mfp,nfp,iface)
            ic = icright

            sign_r = 1
            Qfr(1:numv) = 0.0d0
            if(famr==1) then
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,sp,jfpr,kfpr,ic)
                    Qfr(1:numv) = Qfr(1:numv) + Qs(1:numv)*Lmat(ifpr,sp)
                end do
                if(ifpr==N+1) sign_r=-1
            else if(famr==2) then
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpr,sp,kfpr,ic)
                    Qfr(1:numv) = Qfr(1:numv) + Qs(1:numv)*Lmat(jfpr,sp)
                end do
                if(jfpr==N+1) sign_r=-1
            else
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpr,jfpr,sp,ic)
                    Qfr(1:numv) = Qfr(1:numv) + Qs(1:numv)*Lmat(kfpr,sp)
                end do
                if(kfpr==N+1) sign_r=-1
            end if

! *********************************************
            IF (vismode==1) THEN

                if(faml==1) then 
                    Qvfi(1:numv,ifpl,jfpl,kfpl,icleft) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)
                else if(faml==2) then  
                    Qvfj(1:numv,ifpl,jfpl,kfpl,icleft) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)
                else
                    Qvfk(1:numv,ifpl,jfpl,kfpl,icleft) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)   
                end if

                if(famr==1) then       
                    Qvfi(1:numv,ifpr,jfpr,kfpr,icright) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)
                else if(famr==2) then   
                    Qvfj(1:numv,ifpr,jfpr,kfpr,icright) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)
                else
                    Qvfk(1:numv,ifpr,jfpr,kfpr,icright) = &
                    0.5d0*Qfl(1:numv)+0.5d0*Qfr(1:numv)
                end if
                     
            END IF

! *********************************************

            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)

            if(faml==1) then
                F1(1:numv,ifpl,jfpl,kfpl,icleft)  = Fnl(1:numv)
            else if(faml==2) then
                G2(1:numv,ifpl,jfpl,kfpl,icleft)  = Fnl(1:numv)
            else 
                H3(1:numv,ifpl,jfpl,kfpl,icleft)  = Fnl(1:numv)
            end if
             

            if(famr==1) then
                F1(1:numv,ifpr,jfpr,kfpr,icright)  = Fnr(1:numv)
            else if(famr==2) then
                G2(1:numv,ifpr,jfpr,kfpr,icright)  = Fnr(1:numv)
            else 
                H3(1:numv,ifpr,jfpr,kfpr,icright)  = Fnr(1:numv)
            end if


        end do  ! do loop over points on interface
        end do  ! do loop over points on interface

    end do  ! do loop over interior faces


    END SUBROUTINE interfaceflux

    SUBROUTINE interfaceVISflux

    use setup3d
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    integer ::      icleft,icright,faml,famr
    integer ::      ifpl,jfpl,ifpr,jfpr,ic,iface,sp,ifinter
    integer :: kfpl,kfpr,mfp,nfp
    double precision, dimension(numv)  :: Qfl,Qfr,Qfl2,Qfr2
    double precision, dimension(numv)  :: Qfl3,Qfr3

    do ifinter=1,NINTER

        iface   = IBFINTER(ifinter)
        icleft  = IF2C(iface,1)
        icright = IF2C(iface,2)

        do nfp=1,N
        do mfp=1,N
            faml = Gfp2Lfp(3,mfp,nfp,iface) 
            ifpl = Gfp2Lfp(5,mfp,nfp,iface)
            jfpl = Gfp2Lfp(6,mfp,nfp,iface)
            kfpl = Gfp2Lfp(7,mfp,nfp,iface)

            ! lets get the left solution vector at flux point
            ic = icleft
            Qfl(1:numv) = 0.d0
            Qfl2(1:numv) = 0.d0
            Qfl3(1:numv) = 0.d0
            if(faml==1) then

                do sp=1,N
                Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                end do

            else if(faml==2) then

                do sp=1,N
                Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                end do

            else 

                do sp=1,N
                Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                end do
            end if

            famr = Gfp2Lfp(4,mfp,nfp,iface) 
            ifpr = Gfp2Lfp(8,mfp,nfp,iface) 
            jfpr = Gfp2Lfp(9,mfp,nfp,iface) 
            kfpr = Gfp2Lfp(10,mfp,nfp,iface) 
            ic = icright

            Qfr(1:numv) = 0.d0
            Qfr2(1:numv) = 0.d0
            Qfr3(1:numv) = 0.d0

            if(famr==1) then

                do sp=1,N
                Qfr(1:numv) = Qfr(1:numv) + nablaQs(1:numv,1,sp,jfpr,kfpr,ic)*Lmat(ifpr,sp)
                Qfr2(1:numv) = Qfr2(1:numv) + nablaQs(1:numv,2,sp,jfpr,kfpr,ic)*Lmat(ifpr,sp)
                Qfr3(1:numv) = Qfr3(1:numv) + nablaQs(1:numv,3,sp,jfpr,kfpr,ic)*Lmat(ifpr,sp)
                end do

            else if(famr==2) then

                do sp=1,N
                Qfr(1:numv) = Qfr(1:numv) + nablaQs(1:numv,1,ifpr,sp,kfpr,ic)*Lmat(jfpr,sp)
                Qfr2(1:numv) = Qfr2(1:numv) + nablaQs(1:numv,2,ifpr,sp,kfpr,ic)*Lmat(jfpr,sp)
                Qfr3(1:numv) = Qfr3(1:numv) + nablaQs(1:numv,3,ifpr,sp,kfpr,ic)*Lmat(jfpr,sp)
                end do

            else

                do sp=1,N
                Qfr(1:numv) = Qfr(1:numv) + nablaQs(1:numv,1,ifpr,jfpr,sp,ic)*Lmat(kfpr,sp)
                Qfr2(1:numv) = Qfr2(1:numv) + nablaQs(1:numv,2,ifpr,jfpr,sp,ic)*Lmat(kfpr,sp)
                Qfr3(1:numv) = Qfr3(1:numv) + nablaQs(1:numv,3,ifpr,jfpr,sp,ic)*Lmat(kfpr,sp)
                end do

            end if

            if(famr==1) then

                nablaQvfi(1:numv,1,ifpr,jfpr,kfpr,icright) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfi(1:numv,2,ifpr,jfpr,kfpr,icright) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfi(1:numv,3,ifpr,jfpr,kfpr,icright) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else if(famr==2) then

                nablaQvfj(1:numv,1,ifpr,jfpr,kfpr,icright) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfj(1:numv,2,ifpr,jfpr,kfpr,icright) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfj(1:numv,3,ifpr,jfpr,kfpr,icright) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else

                nablaQvfk(1:numv,1,ifpr,jfpr,kfpr,icright) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfk(1:numv,2,ifpr,jfpr,kfpr,icright) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfk(1:numv,3,ifpr,jfpr,kfpr,icright) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            end if

            if(faml==1) then

                nablaQvfi(1:numv,1,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfi(1:numv,2,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfi(1:numv,3,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else if(faml==2) then

                nablaQvfj(1:numv,1,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfj(1:numv,2,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfj(1:numv,3,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else

                nablaQvfk(1:numv,1,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfk(1:numv,2,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfk(1:numv,3,ifpl,jfpl,kfpl,icleft) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            end if

        end do  ! do loop over points on interface, nfp
        end do  ! do loop over points on interface, mfp
    end do        ! do loop over interior faces

    end SUBROUTINE  interfaceVISflux