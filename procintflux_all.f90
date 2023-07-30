    SUBROUTINE PROCINTFLUX 

    use setup3d 

    IMPLICIT NONE
    include 'mpif.h'

    integer :: icright,faml,famr
    integer :: ifpl,jfpl,ifpr,jfpr,ic,sp,ifinter,sign_l,sign_r
    integer :: kfpl,kfpr,nfp,mfp,count
    integer :: k,MLENGTH,mfprm,nfprm,M,IB
    INTEGER :: TAG,DEST,SOURCE,IC2,ivar
    INTEGER :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)
 
    double precision, dimension(numv)  :: Qs,Fnl,Fnr,Qfl,Qfr
    double precision, dimension(3)  :: normf
    double precision                :: eigv
    double precision, dimension(:),allocatable :: temp

    double precision, allocatable :: RBUFC(:,:),Qflc(:,:,:,:),&
    Qfrc(:,:,:,:)

    MLENGTH = N*N*numv

    allocate(temp(MLENGTH),RBUFC(MLENGTH,NPROCINT))
    allocate(Qflc(numv,N,N,NPROCINT),Qfrc(numv,N,N,NPROCINT))

! FIRST POST NON-BLOCKING RECEIVES

    do IB=1,NPROCINT
        K = IBFPROC(IB)
        TAG = K 
        SOURCE=PROCINT2PROC(IB)
        CALL MPI_IRECV(RBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(IB),IERROR)
    end do 

! Then loop over faces and send Qfl

    do IB=1,NPROCINT
        K= IBFPROC(IB)
        ic = IF2C(K,1)

        TAG = PROCINT2F_PROC(IB)
        DEST = PROCINT2PROC(IB)
        M = 1
 
        do nfp=1,N
        do mfp=1,N

            faml = Gfp2Lfp(3,mfp,nfp,K)

            ifpl = Gfp2Lfp(5,mfp,nfp,K)
            jfpl = Gfp2Lfp(6,mfp,nfp,K)
            kfpl = Gfp2Lfp(7,mfp,nfp,K)

            Qflc(1:numv,mfp,nfp,IB) = 0.0d0

            if(faml==1) then

                do sp=1,N
                    Qs(1:numv) = Q(1:numv,sp,jfpl,kfpl,ic)
                    Qflc(1:numv,mfp,nfp,IB) = &
                    Qflc(1:numv,mfp,nfp,IB) + Qs(1:numv)*Lmat(ifpl,sp)
                end do
                  
            else if(faml==2) then
  
                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpl,sp,kfpl,ic)
                    Qflc(1:numv,mfp,nfp,IB) = &
                    Qflc(1:numv,mfp,nfp,IB) + Qs(1:numv)*Lmat(jfpl,sp)
                end do

            else

                do sp=1,N
                    Qs(1:numv) = Q(1:numv,ifpl,jfpl,sp,ic)
                    Qflc(1:numv,mfp,nfp,IB) = &
                    Qflc(1:numv,mfp,nfp,IB) + Qs(1:numv)*Lmat(kfpl,sp)
                end do

            end if

            do ivar = 1,numv
                temp(M+ivar-1) = Qflc(ivar,mfp,nfp,IB)
            end do

            M=M+numv

        end do  ! do loop over points on interface
        end do  ! do loop over points on interface
            
        ! Now send the buffer
        CALL MPI_SEND(temp,MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)

    end do  ! do loop over processor interface faces

! ----------------------------------------------------------
    CALL MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)
! ------------------------------------------------------------

! Now receive the information and store in appropriate Qfr

    do IB=1,NPROCINT
        M = 1
! Loop over remote values of mfp and nfp
        do nfprm=1,N
        do mfprm=1,N

! Local values of mfp and nfp 
            mfp = map1(1,mfprm,nfprm,IB)
            nfp = map1(2,mfprm,nfprm,IB)
            ! Now we have all the left information
            ! Need to receive Qfr 

            do ivar = 1,numv
                Qfrc(ivar,mfp,nfp,IB) = RBUFC(M+ivar-1,IB)
            end do

            M=M+numv 
        end do
        end do
    end do
! -----------------------------------------------------
! Now that we have Qfl and Qfr, we can call getrusanov
! -----------------------------------------------------

    do IB=1,NPROCINT
        K= IBFPROC(IB)
        ic  = IF2C(K,1)
        icright = IF2C(K,2)

        do nfp=1,N
        do mfp=1,N

            faml = Gfp2Lfp(3,mfp,nfp,K)

            ifpl = Gfp2Lfp(5,mfp,nfp,K)
            jfpl = Gfp2Lfp(6,mfp,nfp,K)
            kfpl = Gfp2Lfp(7,mfp,nfp,K)
!
! Qfl has already been calculated so don't need 
! to do anything except calculating S1,S2 and S3
! to do anything except calculating S1,S2 and S3

            sign_l =1

            if(faml==1) then
                normf(1:3) = S1(1,1:3,ifpl,jfpl,kfpl,ic)
                if(ifpl==1) sign_l=-1
            else if(faml==2) then
                normf(1:3) = S2(2,1:3,ifpl,jfpl,kfpl,ic)
                if(jfpl==1) sign_l=-1
            else
                normf(1:3) = S3(3,1:3,ifpl,jfpl,kfpl,ic)
                if(kfpl==1) sign_l=-1
            end if

            ! Now we have all the left information

            if (IF2C(K,4) .eq. 2 .OR. & 
            IF2C(K,4) .eq. 4 .OR. IF2C(K,4) .eq. 6) then
                sign_r = -1
            else
                sign_r = 1
            endif 
            
            Qfl(1:numv) = Qflc(1:numv,mfp,nfp,IB)
            Qfr(1:numv) = Qfrc(1:numv,mfp,nfp,IB)

            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,&
                                normf,sign_l,sign_r,eigv)


            if(faml==1) then
                F1(1:numv,ifpl,jfpl,kfpl,ic)  = Fnl(1:numv)
            else if(faml==2) then
                G2(1:numv,ifpl,jfpl,kfpl,ic)  = Fnl(1:numv)
            else 
                H3(1:numv,ifpl,jfpl,kfpl,ic)  = Fnl(1:numv)
            end if

            if(faml==1) then 
                Qvfi(1:numv,ifpl,jfpl,kfpl,ic) = &
                0.5d0*Qfl(1:numv) + 0.5d0*Qfr(1:numv)
            else if(faml==2) then  
                Qvfj(1:numv,ifpl,jfpl,kfpl,ic) = &
                0.5d0*Qfl(1:numv) + 0.5d0*Qfr(1:numv)
            else
                Qvfk(1:numv,ifpl,jfpl,kfpl,ic) = &
                0.5d0*Qfl(1:numv) + 0.5d0*Qfr(1:numv)   
            end if

        end do  ! do loop over points on interface
        end do  ! do loop over points on interface

    end do  ! do loop over processor interface faces

    deallocate(temp,RBUFC)
    deallocate(Qflc,Qfrc)

    END SUBROUTINE PROCINTFLUX 

    SUBROUTINE PROCINTVISFLUX

    use setup3d

    IMPLICIT NONE
    include 'mpif.h'

    integer :: icright,ifacelc,ifacerc,faml,famr
    integer :: ifpl,jfpl,ifpr,jfpr,ic,iface,sp,ifinter,sign_l,sign_r
    integer :: kfpl,kfpr,nfp,mfp,count
    integer :: K,MLENGTH,mfprm,nfprm,M
    INTEGER :: TAG,DEST,SOURCE,IC2,ivar
    INTEGER :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)

    double precision, allocatable :: RBUF3C(:,:),Qflc(:,:,:,:),&
    Qfrc(:,:,:,:),Qfl2c(:,:,:,:),Qfr2c(:,:,:,:),Qfl3c(:,:,:,:),&
    Qfr3c(:,:,:,:)

    double precision, dimension(numv)  :: Qs,Fnl,Fnr,Qfl,Qfr,Qfl2,Qfr2,Qfl3,Qfr3
    double precision, dimension(3)  :: normf
    double precision                :: eigv
    double precision, dimension(:),allocatable :: temp

    MLENGTH = N*N*numv*3

    allocate(temp(MLENGTH))
    allocate(RBUF3C(MLENGTH,NPROCINT))
    allocate(Qflc(numv,N,N,NPROCINT),Qfrc(numv,N,N,NPROCINT))
    allocate(Qfl2c(numv,N,N,NPROCINT),Qfr2c(numv,N,N,NPROCINT))
    allocate(Qfl3c(numv,N,N,NPROCINT),Qfr3c(numv,N,N,NPROCINT))

! FIRST POST NON-BLOCKING RECEIVES

    do K=1,NPROCINT
        iface = IBFPROC(K)
        TAG = iface
        SOURCE=PROCINT2PROC(K)
        CALL MPI_IRECV(RBUF3C(1,K),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(K),IERROR)
    end do

! Then loop over faces and send Qfl,Qfl2,Qfl3

    do k=1,NPROCINT
        iface   = IBFPROC(k)
        ic  = IF2C(iface,1)

        TAG = PROCINT2F_PROC(k)
        DEST = PROCINT2PROC(k)
        M = 1

        do nfp=1,N
        do mfp=1,N

            faml = Gfp2Lfp(3,mfp,nfp,iface)

            ifpl = Gfp2Lfp(5,mfp,nfp,iface)
            jfpl = Gfp2Lfp(6,mfp,nfp,iface)
            kfpl = Gfp2Lfp(7,mfp,nfp,iface)

            Qflc(1:numv,mfp,nfp,K) = 0.d0
            Qfl2c(1:numv,mfp,nfp,K) = 0.d0
            Qfl3c(1:numv,mfp,nfp,K) = 0.d0

            if(faml==1) then

                do sp=1,N
                    Qflc(1:numv,mfp,nfp,K) = Qflc(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,1,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl2c(1:numv,mfp,nfp,K) = Qfl2c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,2,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl3c(1:numv,mfp,nfp,K) = Qfl3c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,3,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                end do

            else if(faml==2) then

                do sp=1,N
                    Qflc(1:numv,mfp,nfp,K) = Qflc(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,1,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl2c(1:numv,mfp,nfp,K) = Qfl2c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,2,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl3c(1:numv,mfp,nfp,K) = Qfl3c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,3,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                end do

            else

                do sp=1,N
                    Qflc(1:numv,mfp,nfp,K) = Qflc(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,1,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl2c(1:numv,mfp,nfp,K) = Qfl2c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,2,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl3c(1:numv,mfp,nfp,K) = Qfl3c(1:numv,mfp,nfp,K) &
                    + nablaQs(1:numv,3,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                end do

            end if

            do ivar = 1,numv
                temp(M+ivar-1) = Qflc(ivar,mfp,nfp,K)
            end do
            M = M+numv

            do ivar = 1,numv
                temp(M+ivar-1) = Qfl2c(ivar,mfp,nfp,K)
            end do
            M = M+numv

            do ivar = 1,numv
                temp(M+ivar-1) = Qfl3c(ivar,mfp,nfp,K)
            end do
            M = M+numv

        end do  ! do loop over points on interface,mfp
        end do  ! do loop over points on interface,nfp

        ! Now send the buffer
        CALL MPI_SEND(temp,mlength,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)

    end do  ! do loop over processor interface faces

! ----------------------------------------------------------
    CALL MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)
! ------------------------------------------------------------

! Now receive the information and store in appropriate Qfr

    do K=1,NPROCINT
        M = 1
        ! Loop over remote values of mfp and nfp
        do nfprm=1,N
        do mfprm=1,N

            ! Local values of mfp and nfp 
            mfp = map1(1,mfprm,nfprm,K)
            nfp = map1(2,mfprm,nfprm,K)

            ! Now we have all the left information
            ! Need to receive Qfr 

            do ivar = 1,numv
                Qfrc(ivar,mfp,nfp,K) = RBUF3C(M+ivar-1,K)
            end do
            M = M+numv

            do ivar = 1,numv
                Qfr2c(ivar,mfp,nfp,K) = RBUF3C(M+ivar-1,K)
            end do
            M = M+numv

            do ivar = 1,numv
                Qfr3c(ivar,mfp,nfp,K) = RBUF3C(M+ivar-1,K)
            end do
            M = M+numv

        end do
        end do
    end do

! -----------------------------------------------------
! Now that we have Qfl,Qfl2,Qfl3 and Qfr,Qfr2,Qfr3 we can calculate nablaQvf 
! -----------------------------------------------------

    do k=1,NPROCINT
        iface   = IBFPROC(k)
        ic  = IF2C(iface,1)
        icright = IF2C(iface,2)

        do nfp=1,N
        do mfp=1,N

            faml = Gfp2Lfp(3,mfp,nfp,iface)

            ifpl = Gfp2Lfp(5,mfp,nfp,iface)
            jfpl = Gfp2Lfp(6,mfp,nfp,iface)
            kfpl = Gfp2Lfp(7,mfp,nfp,iface)

            Qfl(1:numv)  = Qflc(1:numv,mfp,nfp,K)
            Qfl2(1:numv) = Qfl2c(1:numv,mfp,nfp,K)
            Qfl3(1:numv) = Qfl3c(1:numv,mfp,nfp,K)

            Qfr(1:numv)  = Qfrc(1:numv,mfp,nfp,K)
            Qfr2(1:numv) = Qfr2c(1:numv,mfp,nfp,K)
            Qfr3(1:numv) = Qfr3c(1:numv,mfp,nfp,K)

            if(faml==1) then

                nablaQvfi(1:numv,1,ifpl,jfpl,kfpl,ic) = & 
                0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfi(1:numv,2,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfi(1:numv,3,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else if(faml==2) then

                nablaQvfj(1:numv,1,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfj(1:numv,2,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfj(1:numv,3,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            else

                nablaQvfk(1:numv,1,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl(1:numv)+0.5*Qfr(1:numv)
                nablaQvfk(1:numv,2,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfk(1:numv,3,ifpl,jfpl,kfpl,ic) = &
                0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)

            end if

        end do  ! do loop over points on interface, mfp
        end do  ! do loop over points on interface, nfp

    end do  ! do loop over processor interface faces
    
    deallocate(temp)
    deallocate(RBUF3C,Qflc,Qfrc,Qfl2c,Qfr2c,Qfl3c,Qfr3c)

    END SUBROUTINE PROCINTVISFLUX 
