    PROGRAM CHORUS_MHD

    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    DOUBLE PRECISION:: t_start,t_finish,t_resolution

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    NPROC = size

    if (rank.eq.0) t_start = mpi_wtime()

    CALL READINPUT

    CALL READ_CELL_DATA

    CALL READ_VRT_DATA

    CALL READ_METIS

    CALL GLOB2LOCAL

    if (CURVE_WALL==1) CALL READ_VRT20_DATA

    CALL READ_BND_DATA

    CALL connectivity

    CALL CONNECT_BDRY

    CALL MATCHPROC

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

    CALL GETIVCELL_PROC

    CALL MAPINTERFACE

    CALL getedges

    CALL init_setup

    CALL iterations

    ! CALL tecplotter3dsetup(NRE)
    ! CALL tecplotter3d(NRE)
    ! CALL TECPLOT1

    if (rank.eq.0) then
        t_finish     = mpi_wtime()
        t_resolution = mpi_wtick()
        print *, "Elapsed time ", t_finish - t_start, " seconds, &
                 resolution ", t_resolution
    end if

    CALL MPI_FINALIZE(IERROR)

    END PROGRAM CHORUS_MHD