    SUBROUTINE CONNECTIVITY
    use setup3d
    IMPLICIT NONE

    include 'mpif.h'

    INTEGER :: IC,ICC,IV,K,I
    INTEGER :: IFA,IFFA,L,IFACE,IFACEOLD
    INTEGER :: IVBASE,IVBASE2,IV1,IVBASE4,KSTABASE,KENDBASE
    INTEGER,ALLOCATABLE :: IFILL(:),IDUMMY(:)
    INTEGER :: IVV1,IVV2,IVV3,IVV4,IPROD11,IPROD12,IPROD13,&
    IPROD14,IPROD21,IPROD22,IPROD23,IPROD24
    INTEGER,ALLOCATABLE :: IF2C_tmp(:,:),IF2V_tmp(:,:)
    INTEGER :: MAXFACE
    INTEGER,DIMENSION(4) :: IVtest
    INTEGER,DIMENSION(2) :: ICtest

    MAXFACE = NCELL * 6

    ALLOCATE(IC2F(NCELL,6))
    ALLOCATE(IF2C_tmp(MAXFACE,2))
    ALLOCATE(IF2V_tmp(MAXFACE,4))
    ALLOCATE(IFILL(NVERT),IDUMMY(NVERT))
    ALLOCATE(ICVSTA(NVERT+1))

    DO IV=1,NVERT
      IFILL(IV)=0
    END DO

    DO IC=1,NCELL
      DO K=1,8
        IV = IVCELL(IC,K)
        IFILL(IV)=IFILL(IV)+1
      END DO
    END DO

    K = 1
    DO IV=1,NVERT
      ICVSTA(IV)=K
      IDUMMY(IV)=ICVSTA(IV)
      K=K+IFILL(IV)
    END DO

    ICVSTA(NVERT+1)=ICVSTA(NVERT)+IFILL(NVERT)
    allocate(ICVERT(ICVSTA(NVERT+1)))

    DO IC=1,NCELL
      DO K=1,8
        IV=IVCELL(IC,K)
        L=IDUMMY(IV)
        ICVERT(L)=IC
        IDUMMY(IV)=IDUMMY(IV)+1
      END DO
    END DO

    IVFACE(1,1)=1
    IVFACE(1,2)=2
    IVFACE(1,3)=3
    IVFACE(1,4)=4

    IVFACE(2,1)=5
    IVFACE(2,2)=6
    IVFACE(2,3)=7
    IVFACE(2,4)=8

    IVFACE(3,1)=5
    IVFACE(3,2)=6
    IVFACE(3,3)=2
    IVFACE(3,4)=1

    IVFACE(4,1)=8
    IVFACE(4,2)=7
    IVFACE(4,3)=3
    IVFACE(4,4)=4

    IVFACE(5,1)=5
    IVFACE(5,2)=1
    IVFACE(5,3)=4
    IVFACE(5,4)=8

    IVFACE(6,1)=6
    IVFACE(6,2)=2
    IVFACE(6,3)=3
    IVFACE(6,4)=7

    DO IC=1,NCELL
      DO IFA=1,6
        IC2F(IC,IFA)=0
      END DO
    END DO

    DO IFA=1,MAXFACE
      IF2C_tmp(IFA,1) = 0
      IF2C_tmp(IFA,2) = 0
    END DO

    IFACE = 0

    DO 100 IC=1,NCELL
      DO 110 IFA=1,6
        IFACEOLD=IFACE
        IF (IC2F(IC,IFA).NE.0) GOTO 110
        IVBASE=IVCELL(IC,IVFACE(IFA,1))
        IVBASE2=IVCELL(IC,IVFACE(IFA,2))
        IV1 =IVCELL(IC,IVFACE(IFA,3))
        IVBASE4=IVCELL(IC,IVFACE(IFA,4))

        KSTABASE=ICVSTA(IVBASE)
        KENDBASE=ICVSTA(IVBASE+1)-1

        DO 120 K=KSTABASE,KENDBASE
          ICC=ICVERT(K)
          IF (ICC.EQ.IC) GOTO 120

          DO 130 IFFA=1,6
            IVV1=IVCELL(ICC,IVFACE(IFFA,1))
            IVV2=IVCELL(ICC,IVFACE(IFFA,2))
            IVV3=IVCELL(ICC,IVFACE(IFFA,3))
            IVV4=IVCELL(ICC,IVFACE(IFFA,4))

            IPROD11=(IVBASE-IVV1)
            IPROD12=(IVBASE-IVV2)
            IPROD13=(IVBASE-IVV3)
            IPROD14=(IVBASE-IVV4)

            IPROD21=(IV1-IVV1)
            IPROD22=(IV1-IVV2)
            IPROD23=(IV1-IVV3)
            IPROD24=(IV1-IVV4)

            IF ((IPROD11.EQ.0.OR.IPROD12.EQ.0.OR. &
                 IPROD13.EQ.0.OR.IPROD14.EQ.0).AND. &
                (IPROD21.EQ.0.OR.IPROD22.EQ.0.OR. &
                 IPROD23.EQ.0.OR.IPROD24.EQ.0)) THEN

                 IFACE = IFACE+1
                 IC2F(ICC,IFFA)=IFACE

                 IF2V_tmp(IFACE,1) = IVBASE
                 IF2V_tmp(IFACE,2) = IVBASE2
                 IF2V_tmp(IFACE,3) = IV1
                 IF2V_tmp(IFACE,4) = IVBASE4

                 IF2C_tmp(IFACE,1)=IC
                 IF2C_tmp(IFACE,2)=ICC
  
                 GOTO 115
            END IF
130       CONTINUE
120     CONTINUE

        IF (IF2C_tmp(IFACEOLD+1,1).EQ.0) THEN
          IF2C_tmp(IFACEOLD+1,1)=IC
          IF2C_tmp(IFACEOLD+1,2)=0
          IFACE=IFACEOLD+1

          IF2V_tmp(IFACE,1) = IVBASE
          IF2V_tmp(IFACE,2) = IVBASE2
          IF2V_tmp(IFACE,3) = IV1
          IF2V_tmp(IFACE,4) = IVBASE4
        END IF

115     IC2F(IC,IFA)=IFACE
110   CONTINUE
100 CONTINUE

    NFACE=IFACE

    DEALLOCATE(IFILL,IDUMMY)
    ALLOCATE(IF2C(NFACE,6))
    ALLOCATE(IF2V(NFACE,4))
    IF2C(1:NFACE,1:2)=IF2C_tmp(1:NFACE,1:2)
    IF2V(1:NFACE,:)=IF2V_tmp(1:NFACE,:)
    DEALLOCATE(IF2C_tmp,IF2V_tmp)

    ! ! examine the connectivity of a particular face
    ! IFACE = 2000
    ! do I = 1,4
    !   IVtest(I) = IF2V(IFACE,I)
    ! end do
    ! write(*,*) 'processor',rank,'the four nodal points associated with this face are',&
    ! IVtest(:),'(local)'
    ! write(*,*) 'processor',rank,'the four nodal points associated with this face are',IV2IVG(IVtest(1)),&
    ! IV2IVG(IVtest(2)),IV2IVG(IVtest(3)),IV2IVG(IVtest(4)),'(global)'
    ! do I = 1,2
    !   ICtest(I) = IF2C(IFACE,I)
    ! end do
    ! write(*,*) 'processor',rank,'the two adjacent cells are ',ICtest(:),'(local)'
    ! write(*,*) 'processor',rank,'the two adjacent cells are ',&
    ! IC2ICG(ICtest(1)),IC2ICG(ICtest(2)),'(global)'

    END SUBROUTINE CONNECTIVITY