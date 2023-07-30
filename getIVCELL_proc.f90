    SUBROUTINE GETIVCELL_PROC

    use setup3d
    IMPLICIT NONE

    include 'mpif.h'
    
    integer :: I,K,IC1,IC2,M,IV
    integer :: TAG,DEST,SOURCE,MLENGTH
    integer :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)

    integer :: VRBUF(8,NPROCINT),VSBUF(8)

    MLENGTH = 8

    DO I=1,NPROCINT

      K = IBFPROC(I)
      TAG = K
      SOURCE = PROCINT2PROC(I)
      CALL MPI_IRECV(VRBUF(1,I),MLENGTH,MPI_INTEGER,SOURCE,TAG,&
                     MPI_COMM_WORLD,REQHANDLE(I),IERROR)
    
    END DO

    DO I=1,NPROCINT

      K = IBFPROC(I)

      IC1=IF2C(K,1)

      TAG = PROCINT2F_PROC(I)

      DEST = PROCINT2PROC(I)

      M=1
      DO IV=1,8
        VSBUF(M) = IV2IVG(IVCELL(IC1,IV))
        M = M+1
      END DO

      CALL MPI_SEND(VSBUF,MLENGTH,MPI_INTEGER,DEST,TAG,&
                    MPI_COMM_WORLD,IERROR)

    END DO

    CALL MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)

    DO I=1,NPROCINT
      K = IBFPROC(I)
      IC2 = IF2C(K,2)

      M=1
      DO IV=1,8
        IVCELL(IC2,IV) = IVG2IV(VRBUF(M,I))
        M = M+1
      END DO
    END DO

    ! if (rank.eq.0) then
    !     I = 1
    !     K = IBFPROC(I)
    !     IC1 = IF2C(K,1)
    !     IC2 = IF2C(K,2)
    !     print *,IC1,IVCELL(IC1,:)
    !     print *,IC2,IVCELL(IC2,:)
    !     PRINT *,'NCELL=',NCELL,'NCELLTOT=',NCELLTOT
    ! end if

    END SUBROUTINE GETIVCELL_PROC