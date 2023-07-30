    SUBROUTINE MATCHPROC
    use setup3d
    IMPLICIT NONE

    include 'mpif.h'

    integer :: MLENGTH,TAG
    integer :: I,M,K,IB,IFA,KLOC,K1,K3
    integer :: REQHANDLE(NPROC-1),ISTAT(MPI_STATUS_SIZE,NPROC-1)
    integer :: count,match_count,match_count_all
    integer :: IVB1,IVB2,IVB3,IVB4,ICSTA1,ICEND1,ICSTA3,ICEND3
    integer :: IVB1L,IVB2L,IVB3L,IVB4L
    integer :: ICC,ICC1,ICC3,IV1,IV2,IV3,IV4
    integer :: IDIF11,IDIF21,IDIF31,IDIF41,IDIF13,IDIF23,IDIF33,IDIF43

    CALL MPI_ALLREDUCE(NPROCINT,MAXPROCINT,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)

    MAXPROCINT = MAXPROCINT + 1

    allocate(PROCSBUF(MAXPROCINT*5))
    allocate(PROCRBUF(MAXPROCINT*5,NPROC))
    allocate(PROCINT2PROC(NPROCINT),PROCINT2F_PROC(NPROCINT))

    PROCSBUF = 0
    PROCRBUF = 0

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    MLENGTH = MAXPROCINT*5

    count = 0

    DO I=1,NPROC
      IF (I.EQ.RANK+1) THEN
        CYCLE
      ELSE
      count = count + 1
      TAG = I*NPROC+RANK+1
      CALL MPI_IRECV(PROCRBUF(1,I),MLENGTH,MPI_INTEGER,I-1,&
                    TAG,MPI_COMM_WORLD,REQHANDLE(count),IERROR)
      END IF
    END DO

    DO I=1,NPROC
      IF (I.EQ.RANK+1) THEN
        CYCLE
      ELSE
        M = 1
        DO IB=1,NPROCINT
          K = IBFPROC(IB)
          PROCSBUF(M) = K 
          PROCSBUF(M+1) = IV2IVG(IF2V(K,1))
          PROCSBUF(M+2) = IV2IVG(IF2V(K,2))
          PROCSBUF(M+3) = IV2IVG(IF2V(K,3))
          PROCSBUF(M+4) = IV2IVG(IF2V(K,4))
          M=M+5
        END DO
      END IF

      TAG = (rank+1)*NPROC+I
      CALL MPI_SEND(PROCSBUF,MLENGTH,MPI_INTEGER,I-1,TAG,&
                    MPI_COMM_WORLD,IERROR)
    END DO

    CALL MPI_WAITALL(NPROC-1,REQHANDLE,ISTAT,IERROR)

    DO IB=1,NPROCINT
      PROCINT2PROC(IB) = -1
      PROCINT2F_PROC(IB) = -1
    END DO


    match_count_all = 0
    DO I=1,NPROC
      match_count = 0
      if (I.EQ.RANK+1) THEN
        CYCLE
      END IF
      M=1
      DO 40 WHILE (PROCRBUF(M,I).NE.0)

        K = PROCRBUF(M,I)
        IVB1 = IVG2IV(PROCRBUF(M+1,I))
        IVB2 = IVG2IV(PROCRBUF(M+2,I))
        IVB3 = IVG2IV(PROCRBUF(M+3,I))
        IVB4 = IVG2IV(PROCRBUF(M+4,I))

        M = M+5

        if ((IVB1.ne.0).and.(IVB2.ne.0).and.(IVB3.ne.0).and.(IVB4.ne.0)) then

          ICSTA1=ICVSTA(IVB1)
          ICEND1=ICVSTA(IVB1+1)-1
          ICSTA3=ICVSTA(IVB3)
          ICEND3=ICVSTA(IVB3+1)-1
          ICC = 0

          DO K1=ICSTA1,ICEND1
            ICC1=ICVERT(K1)
            DO K3=ICSTA3,ICEND3
              ICC3=ICVERT(K3)
              IF (ICC1.EQ.ICC3) THEN
                ICC=ICC1
                GOTO 15
              END IF
            END DO
          END DO

          IF (ICC.eq.0) THEN
            WRITE(*,*) 'No common cell ---- PROBLEM!!!'
            CYCLE
          END IF 

          
15        DO IFA=1,6

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

            IF ((IDIF11.EQ.0.OR.IDIF21.EQ.0.OR. &
               IDIF31.EQ.0.OR.IDIF41.EQ.0).AND. &
               (IDIF13.EQ.0.OR.IDIF23.EQ.0.OR. &
               IDIF33.EQ.0.OR.IDIF43.EQ.0)) GOTO 41
          END DO

          WRITE(*,*) 'Problem: Could not find common face'
          
41        KLOC=IC2F(ICC,IFA)
          IB = IF2IBPROC(KLOC)
          PROCINT2PROC(IB) = I-1
          PROCINT2F_PROC(IB) = K 
          match_count= match_count + 1

        END IF
      
40    CONTINUE

!       WRITE(*,500) 'processor',rank,'has matched',match_count,&
!       'processor interfaces from rank',I-1
! 500   FORMAT(A,I6,1X,A,I8,1X,A,I6)
      match_count_all = match_count_all + match_count

    END DO

    WRITE(*,200) 'processor',rank,'has matched',match_count_all,&
    'processor interfaces'
200 FORMAT(A,I6,1X,A,I8,1X,A,I6)

    deallocate(PROCSBUF,PROCRBUF)

    END SUBROUTINE MATCHPROC


    SUBROUTINE MAPPROCINT
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: k,mfp,nfp,C2F,IB

    DO IB = 1,NPROCINT
      K = IBFPROC(IB)
      C2F = IF2C(K,4)

      DO mfp = 1,N
      DO nfp = 1,N

        if (C2F.eq.1.OR.C2F.eq.3.OR.C2F.eq.5) then 
          if(IF2C(K,6).eq.1) then
            map1(1,mfp,nfp,IB) = nfp
            ! Meaning that mfp_remote is same point as nfp_local
            map1(2,mfp,nfp,IB) = mfp
            ! Meaning that nfp_remote is same point as mfp_local
          elseif(IF2C(K,6).eq.2) then 
            map1(1,mfp,nfp,IB) = N-mfp+1
            map1(2,mfp,nfp,IB) = nfp
          elseif(IF2C(K,6).eq.3) then 
            map1(1,mfp,nfp,IB) = N-nfp+1
            map1(2,mfp,nfp,IB) = N-mfp+1
          elseif(IF2C(K,6).eq.4) then
            map1(1,mfp,nfp,IB) = mfp
            map1(2,mfp,nfp,IB) = N-nfp+1
          end if 
        elseif (C2F.eq.2.OR.C2F.eq.4.OR.C2F.eq.6) then 
          if(IF2C(K,6).eq.1) then
            map1(1,mfp,nfp,IB) = nfp
            map1(2,mfp,nfp,IB) = mfp
          elseif(IF2C(K,6).eq.2) then 
            map1(1,mfp,nfp,IB) = mfp
            map1(2,mfp,nfp,IB) = N-nfp+1
          elseif(IF2C(K,6).eq.3) then 
            map1(1,mfp,nfp,IB) = N-nfp+1
            map1(2,mfp,nfp,IB) = N-mfp+1
          elseif(IF2C(K,6).eq.4) then
            map1(1,mfp,nfp,IB) = N-mfp+1
            map1(2,mfp,nfp,IB) = nfp
          end if 
        else 
          write(*,*) 'Something WRONG!!!!'
        end if
      
      END DO
      END DO
    
    END DO

    END SUBROUTINE MAPPROCINT

        