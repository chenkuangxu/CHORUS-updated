    SUBROUTINE CONNECT_BDRY
    use setup3d
    implicit none
    include 'mpif.h'

    INTEGER :: K,IB,IB2,IVB1,IVB2,IVB3,IVB4,ICSTA1,ICEND1,ICSTA3,ICEND3,ICC
    INTEGER :: K1,K3,ICC1,ICC3,IFA
    INTEGER :: IDIF11,IDIF21,IDIF31,IDIF41,IDIF13,IDIF23,IDIF33,IDIF43
    INTEGER :: IV1,IV2,IV3,IV4,IC1,IC2
    INTEGER,ALLOCATABLE :: IBFPROC_tmp(:),IBFCYCLOC_tmp(:),IBFCYCREM_tmp(:),IVCELL_tmp(:,:)
    DOUBLE PRECISION :: XFIN,YFIN,ZFIN,XFOU,YFOU,ZFOU,epsXY,epsZ,epsXZ,epsY,epsYZ,epsX

    ALLOCATE(BOUNFACEINDEX(NFACE))

    do K = 1,NFACE
        BOUNFACEINDEX(K) = 0
    end do

    ! First we need to find the NINLET,NOUTLET,NSYMP,NWALL

    NINLET=0
    NOUTLET=0
    NSYMP=0
    NWALL=0
    NPRES=0
    NWALM=0
    NCYCLIC=0

    do 10 IB=1,NBOUNDEFINED

        K=0
        IVB1=IVG2IV(IVGBOUN(IB,1))
        IVB2=IVG2IV(IVGBOUN(IB,2))
        IVB3=IVG2IV(IVGBOUN(IB,3))
        IVB4=IVG2IV(IVGBOUN(IB,4))

        if ((IVB1.EQ.0).OR.(IVB2.EQ.0)&
        .OR.(IVB3.EQ.0).OR.IVB4.EQ.0) then
            goto 10
        end if

        ICSTA1=ICVSTA(IVB1)
        ICEND1=ICVSTA(IVB1+1)-1
        ICSTA3=ICVSTA(IVB3)
        ICEND3=ICVSTA(IVB3+1)-1
        ICC=0

        do K1=ICSTA1,ICEND1
            ICC1=ICVERT(K1)
            do K3=ICSTA3,ICEND3
                ICC3=ICVERT(K3)
                if (ICC1.EQ.ICC3) ICC=ICC1
            end do
        end do

        if (ICC.EQ.0) then
            goto 10
        end if

        do IFA=1,6
            IV1=IVCELL(ICC,IVFACE(IFA,1))
            IV2=IVCELL(ICC,IVFACE(IFA,2))
            IV3=IVCELL(ICC,IVFACE(IFA,3))
            IV4=IVCELL(ICC,IVFACE(IFA,4))
            IDIF11=IV1-IVB1
            IDIF21=IV2-IVB1
            IDIF31=IV3-IVB1
            IDIF41=IV4-IVB1
            IDIF13=IV1-IVB3
            IDIF23=IV2-IVB3
            IDIF33=IV3-IVB3
            IDIF43=IV4-IVB3
            if ((IDIF11.EQ.0.OR.IDIF21.EQ.0.OR. &
               IDIF31.EQ.0.OR.IDIF41.EQ.0).AND. &
              (IDIF13.EQ.0.OR.IDIF23.EQ.0.OR. &
               IDIF33.EQ.0.OR.IDIF43.EQ.0)) goto 41
        end do

        write(*,*) 'Should not be here'

41      K=IC2F(ICC,IFA)
        BOUNFACEINDEX(K)=1

        if (K.EQ.0) then
            WRITE(*,*) 'Problem with the attachement of boundaries'
            WRITE(*,*) 'IB =',IB,'ICC =',ICC
            WRITE(*,*) 'Execution paused !'
        end if

        if (IBOUNTYPE(IB).EQ.'INLE') then
            NINLET=NINLET+1
        end if

        if (IBOUNTYPE(IB).EQ.'OUTL') then
            NOUTLET=NOUTLET+1
        end if

        if (IBOUNTYPE(IB).EQ.'SYMP') then
            NSYMP=NSYMP+1
        end if

        if (IBOUNTYPE(IB).EQ.'WALL') then
            NWALL=NWALL+1
        end if

        if (IBOUNTYPE(IB).EQ.'PRES') then
            NPRES=NPRES+1
        end if

        if (IBOUNTYPE(IB).EQ.'WALM') then
            NWALM=NWALM+1
        end if

        if (IBOUNTYPE(IB).EQ.'CYCL') then
            NCYCLIC=NCYCLIC+1
        end if

10  end do

    ALLOCATE(IBFINL(NINLET),IBFOUT(NOUTLET),IBFSYMP(NSYMP),IBFWAL(NWALL),&
    IBFPRES(NPRES),IBFWALM(NWALM),IBFCYC(NCYCLIC))

! =================================================================================================
    NBOUN=0
    NINLET=0
    NOUTLET=0
    NSYMP=0
    NWALL=0
    NPRES=0
    NWALM=0
    NCYCLIC=0

    do 20 IB=1,NBOUNDEFINED

        K=0
        IVB1=IVG2IV(IVGBOUN(IB,1))
        IVB2=IVG2IV(IVGBOUN(IB,2))
        IVB3=IVG2IV(IVGBOUN(IB,3))
        IVB4=IVG2IV(IVGBOUN(IB,4))

        if ((IVB1.EQ.0).OR.(IVB2.EQ.0)&
        .OR.(IVB3.EQ.0).OR.IVB4.EQ.0) then
            goto 20
        end if

        ICSTA1=ICVSTA(IVB1)
        ICEND1=ICVSTA(IVB1+1)-1
        ICSTA3=ICVSTA(IVB3)
        ICEND3=ICVSTA(IVB3+1)-1
        ICC=0

        do K1=ICSTA1,ICEND1
            ICC1=ICVERT(K1)
            do K3=ICSTA3,ICEND3
                ICC3=ICVERT(K3)
                if (ICC1.EQ.ICC3) ICC=ICC1
            end do
        end do

        if (ICC.EQ.0) then
            goto 20
        end if

        do IFA=1,6
            IV1=IVCELL(ICC,IVFACE(IFA,1))
            IV2=IVCELL(ICC,IVFACE(IFA,2))
            IV3=IVCELL(ICC,IVFACE(IFA,3))
            IV4=IVCELL(ICC,IVFACE(IFA,4))
            IDIF11=IV1-IVB1
            IDIF21=IV2-IVB1
            IDIF31=IV3-IVB1
            IDIF41=IV4-IVB1
            IDIF13=IV1-IVB3
            IDIF23=IV2-IVB3
            IDIF33=IV3-IVB3
            IDIF43=IV4-IVB3
            if ((IDIF11.EQ.0.OR.IDIF21.EQ.0.OR. &
               IDIF31.EQ.0.OR.IDIF41.EQ.0).AND. &
              (IDIF13.EQ.0.OR.IDIF23.EQ.0.OR. &
               IDIF33.EQ.0.OR.IDIF43.EQ.0)) goto 21
        end do

        write(*,*) 'Should not be here'

21      K=IC2F(ICC,IFA)
        BOUNFACEINDEX(K)=1

        if (K.EQ.0) then
            WRITE(*,*) 'Problem with the attachement of boundaries'
            WRITE(*,*) 'IB =',IB,'ICC =',ICC
            WRITE(*,*) 'Execution paused !'
        end if

        if (IBOUNTYPE(IB).EQ.'INLE') then
            NINLET=NINLET+1
            NBOUN=NBOUN+1
            IBFINL(NINLET)=K
        end if

        if (IBOUNTYPE(IB).EQ.'OUTL') then
            NOUTLET=NOUTLET+1
            NBOUN=NBOUN+1
            IBFOUT(NOUTLET)=K
        end if

        if (IBOUNTYPE(IB).EQ.'SYMP') then
            NSYMP=NSYMP+1
            NBOUN=NBOUN+1
            IBFSYMP(NSYMP)=K
        end if

        if (IBOUNTYPE(IB).EQ.'WALL') then
            NWALL=NWALL+1
            NBOUN=NBOUN+1
            IBFWAL(NWALL)=K
        end if

        if (IBOUNTYPE(IB).EQ.'PRES') then
            NPRES=NPRES+1
            NBOUN=NBOUN+1
            IBFPRES(NPRES)=K
        end if

        if (IBOUNTYPE(IB).EQ.'WALM') then
            NWALM=NWALM+1
            NBOUN=NBOUN+1
            IBFWALM(NWALM)=K
        end if

        if (IBOUNTYPE(IB).EQ.'CYCL') then
            NCYCLIC=NCYCLIC+1
            NBOUN=NBOUN+1
            IBFCYC(NCYCLIC)=K
        end if

20  end do

    NCELLTOT = NCELL
    !====================================================================================

    allocate(IBFPROC_tmp(NFACE))
    allocate(IF2IBPROC(NFACE))

    IF2IBPROC(:) = 0

    NPROCINT = 0

    do 210 K=1,NFACE
        IC1=IF2C(K,1)
        IC2=IF2C(K,2)

        if (IC2.EQ.0.AND.BOUNFACEINDEX(K).EQ.0) then
            NPROCINT=NPROCINT+1
            IBFPROC_tmp(NPROCINT)=K
            IF2IBPROC(K) = NPROCINT

            NCELLTOT = NCELLTOT + 1
            IF2C(K,2)=NCELLTOT
        end if

210 end do

    allocate(IBFPROC(NPROCINT))
    IBFPROC(1:NPROCINT) = IBFPROC_tmp(1:NPROCINT)
    deallocate(IBFPROC_tmp)

! =========================================================================================

    allocate(IBFCYCLOC_tmp(NFACE),IBFCYCREM_tmp(NFACE))
    allocate(IF2IBCYCREM(NFACE))
    IF2IBCYCREM(:) = 0
    NCYCLOC = 0
    NCYCREM = 0

    do 160 IB=1,NCYCLIC
        K = IBFCYC(IB)
        XFIN = ( XV(IF2V(K,1)) + XV(IF2V(K,2)) &
        + XV(IF2V(K,3)) + XV(IF2V(K,4)) ) / 4
        YFIN = ( YV(IF2V(K,1)) + YV(IF2V(K,2)) &
        + YV(IF2V(K,3)) + YV(IF2V(K,4)) ) / 4
        ZFIN = ( ZV(IF2V(K,1)) + ZV(IF2V(K,2)) &
        + ZV(IF2V(K,3)) + ZV(IF2V(K,4)) ) / 4

        do IB2=1,NCYCLIC
            K = IBFCYC(IB2)
            XFOU = ( XV(IF2V(K,1)) + XV(IF2V(K,2)) &
            + XV(IF2V(K,3)) + XV(IF2V(K,4)) ) / 4
            YFOU = ( YV(IF2V(K,1)) + YV(IF2V(K,2)) &
            + YV(IF2V(K,3)) + YV(IF2V(K,4)) ) / 4
            ZFOU = ( ZV(IF2V(K,1)) + ZV(IF2V(K,2)) &
            + ZV(IF2V(K,3)) + ZV(IF2V(K,4)) ) / 4

            epsXY = dsqrt( (XFIN-XFOU)**2 + (YFIN-YFOU)**2 )
            epsZ  = abs(abs(ZFIN-ZFOU) - DZCYCL)

            epsXZ = dsqrt( (XFIN-XFOU)**2 + (ZFIN-ZFOU)**2 )
            epsY  = abs(abs(YFIN-YFOU) - DYCYCL)

            epsYZ = dsqrt( (YFIN-YFOU)**2 + (ZFIN-ZFOU)**2 )
            epsX  = abs(abs(XFIN-XFOU) - DXCYCL)

            if ( (epsXY.lt.tolCYC.and.epsZ.lt.tolCYC) .or. &
            (epsXZ.lt.tolCYC.and.epsY.lt.tolCYC) .or. &
            (epsYZ.lt.tolCYC.and.epsX.lt.tolCYC) ) then
                NCYCLOC = NCYCLOC + 1
                IBFCYCLOC_tmp(NCYCLOC) = IBFCYC(IB)
                goto 160
            end if
        end do

        NCYCREM = NCYCREM + 1
        IBFCYCREM_tmp(NCYCREM) = IBFCYC(IB)
        IF2IBCYCREM(IBFCYC(IB)) = IB

        NCELLTOT = NCELLTOT + 1
        IF2C(IBFCYC(IB),2) = NCELLTOT

160 continue

    allocate(IBFCYCLOC(NCYCLOC),IBFCYCREM(NCYCREM))
    IBFCYCLOC(1:NCYCLOC) = IBFCYCLOC_tmp(1:NCYCLOC)
    deallocate(IBFCYCLOC_tmp)
    IBFCYCREM(1:NCYCREM) = IBFCYCREM_tmp(1:NCYCREM)
    deallocate(IBFCYCREM_tmp)


    allocate(IVCELL_tmp(NCELL,8))
    IVCELL_tmp = IVCELL
    deallocate(IVCELL)
    allocate(IVCELL(NCELLTOT,8))
    IVCELL = 0
    IVCELL(1:NCELL,:) = IVCELL_tmp(1:NCELL,:)
    deallocate(IVCELL_tmp)

    WRITE(*,*) '--------------------------------'
    WRITE(*,*) 'Info from connectivity'
    WRITE(*,*) '--------------------------------'
    WRITE(*,*) 'RANK =',rank
    WRITE(*,*) 'Number of cells (NCELL) =',NCELL
    WRITE(*,*) 'Total number of cells (NCELLTOT) &
                  (internal+procint)=',NCELLTOT
    WRITE(*,*) 'Total number of faces: ',NFACE
    WRITE(*,*) 'Number of total boundary faces=',NBOUN
    WRITE(*,*) 'Boundary faces attached'
    WRITE(*,*) 'Inlet faces attached=',NINLET
    WRITE(*,*) 'Outlet faces attached=',NOUTLET
    WRITE(*,*) 'Processor interfaces attached=',NPROCINT
    WRITE(*,*) '--------------------------------'

    DEALLOCATE(IBOUNTYPE,BOUNFACEINDEX)

    END SUBROUTINE CONNECT_BDRY
