    SUBROUTINE READINPUT

    use setup3d
    IMPLICIT NONE

    include 'mpif.h'

    integer :: I,IOERR

    OPEN(30,FILE='QUAD.INP')
        READ(30,*)
        READ(30,*)
        READ(30,*) NAMECEL
        READ(30,*) NAMEVRT
        READ(30,*) NAMEVRT20
        READ(30,*) NAMEBND
        READ(30,*) METIS_CELL
        READ(30,*)
        READ(30,*) SOLARMODEL
        READ(30,*)
        READ(30,*) N
        READ(30,*) CURVE_WALL
        READ(30,*)
        READ(30,*) k_stage
        READ(30,*) dt
        READ(30,*) MAXITER
        READ(30,*)
        READ(30,*) vismode
        READ(30,*) mu
        READ(30,*) eta0
        READ(30,*) kappa_s
        READ(30,*) 
        READ(30,*) restart
        READ(30,*) restart_ord
        READ(30,*) NAMERESTART
        READ(30,*)
        READ(30,*) NRE
        READ(30,*) nwrite
    CLOSE(30)

    I = 0
    IOERR = 0
    OPEN(31,FILE=SOLARMODEL)
    DO WHILE (IOERR.EQ.0)
      read(31,*,iostat=IOERR)
      I = I+1
    END DO
    CLOSE(31)

    MAXSM = I-1
    ALLOCATE(rIn(MAXSM),rhoIn(MAXSM),TIn(MAXSM),radcIn(MAXSM),dTIn(MAXSM))

    OPEN(31,FILE=SOLARMODEL)
    DO I = 1,MAXSM
      read(31,101) rIn(I),rhoIn(I),TIn(I),radcIn(I),dTIn(I)
    END DO
101 FORMAT(5(D20.12,1X))
    CLOSE(31)

    if (rank.eq.0) WRITE(*,102) 'For setting initial conditions,&
    radius is partitioned into',MAXSM,'intervals'
102 FORMAT(A60,1X,I7,1X,A10)

    rinf = rhoIn(MAXSM)
    pinf = R_const * rinf * TIn(MAXSM)
    if (rank.eq.0) then
      print *,'Top boundaries, rho=',rinf,',pressure =',pinf
    end if

    END SUBROUTINE READINPUT