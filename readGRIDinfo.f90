    SUBROUTINE READ_CELL_DATA

    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: I,IC,K,ICTYPE,ICKEY,IOERR

    !go through the file to get total number of cells

    if(rank.eq.0) then
      WRITE(*,*) 'reading cells ...'    
      I = 0
      IOERR = 0
      open(20,file=NAMECEL)
      do while(IOERR.EQ.0)
        read(20,*,iostat=IOERR)
        I=I+1
      end do
      close(20)

      NCELLGLOB = I-1
      write(*,'(A,I6)') 'number of global cells=', NCELLGLOB
    endif

    call MPI_BCAST(NCELLGLOB,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

    allocate(IVGCELLG(NCELLGLOB,8))

    !go through the file again to read data
    open(20,file=NAMECEL)
    do I=1,NCELLGLOB
      read(20,*) IC,IVGCELLG(I,1:8),ICTYPE,ICKEY
    end do
    close(20)

    END SUBROUTINE READ_CELL_DATA


    ! read XVG YVG ZVG
    SUBROUTINE READ_VRT_DATA

    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: I,IV,IOERR

    !go through the file to get total number of vert

    if(rank.eq.0) then
      
      write(*,*) 'reading vertices ...'

      I = 0
      IOERR = 0
      open(10,file=NAMEVRT)
      do while(IOERR.EQ.0)
            read(10,*,iostat=IOERR)
            I=I+1
      end do
      close(10)

      NVERTGLOB = I-1
      write(*,'(A,I6)') ' number of global vertices = ',NVERTGLOB

    endif

    call MPI_BCAST(NVERTGLOB,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

    allocate(XVG(NVERTGLOB))
    allocate(YVG(NVERTGLOB))
    allocate(ZVG(NVERTGLOB))

    !go throught the file again to get the data

    open(10,file=NAMEVRT)

    do I=1, NVERTGLOB
      read(10,1000) IV,XVG(I),YVG(I),ZVG(I)
1000  FORMAT(I9,6X,3(D20.12,1X))
    end do

    close(10)

    END SUBROUTINE READ_VRT_DATA


    SUBROUTINE READ_VRT20_DATA
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: ICG,K,IC
    double precision,dimension(3) :: xx_tmp

    if (rank.eq.0) then
      print *,'reading quadratic vertices for curved boundaries ...' 
    end if

    ALLOCATE(xxf_cell(3,20,NCELL))

    open(12,FILE=NAMEVRT20)
    do ICG=1,NCELLGLOB
      if (ICG2IC(ICG).NE.0) then
        IC = ICG2IC(ICG)
        do K=1,20
          read(12,*) xxf_cell(1:3,K,IC)
        end do
      else
        do K=1,20
          read(12,*) xx_tmp(1:3)
        end do
      end if
    end do

    ! if (rank.eq.0) then
    !   print *,xxf_cell(1:3,4,100)
    !   print *,'---'
    !   print *,XV(IVCELL(100,4)),YV(IVCELL(100,4)),ZV(IVCELL(100,4))
    ! end if

    END SUBROUTINE READ_VRT20_DATA

    
    SUBROUTINE READ_METIS
    use setup3d
    IMPLICIT NONE

    include 'mpif.h'

    integer :: I
    
    allocate(IC2PROC(NCELLGLOB))
    allocate(PROC_NCELL(NPROC))

    NCELL = 0

    if (NPROC.NE.1) then

      open(13,FILE=METIS_CELL)
      do I = 1,NCELLGLOB
        READ(13,*) IC2PROC(I)
        if (IC2PROC(I).EQ.rank) NCELL = NCELL + 1
      end do
      CLOSE(13)

    else if (NPROC.EQ.1) then

      IC2PROC(:) = 0

      NCELL = NCELLGLOB

    end if

    ! print *,'processor',rank,'has',NCELL,'local cells'

    CALL MPI_ALLGATHER(NCELL,1,MPI_INTEGER,PROC_NCELL,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)

    END SUBROUTINE READ_METIS


    SUBROUTINE GLOB2LOCAL
    use setup3d
    IMPLICIT NONE

    include 'mpif.h'
    integer :: IC,IV,IVNEW,ICNEW,I
    integer :: IVG,IVL
    integer,allocatable :: IV2IVG_tmp(:)
    double precision,allocatable :: XV_tmp(:),YV_tmp(:),ZV_tmp(:)

    allocate(IVCELL(NCELL,8))
    allocate(ICG2IC(NCELLGLOB))
    allocate(IC2ICG(NCELL))
    allocate(IVG2IV(NVERTGLOB))

    allocate(XV_tmp(8*NCELL))
    allocate(YV_tmp(8*NCELL))
    allocate(ZV_tmp(8*NCELL))
    allocate(IV2IVG_tmp(8*NCELL))

    ICG2IC = 0
    IVNEW  = 0
    ICNEW  = 0
    IVG2IV = 0

    do IC = 1,NCELLGLOB
        
      if (IC2PROC(IC).EQ.rank) then
          
        ICNEW = ICNEW + 1
        ICG2IC(IC) = ICNEW
        IC2ICG(ICNEW) = IC

        do IV = 1,8
          IVG = IVGCELLG(IC,IV)
          if (IVG2IV(IVG).EQ.0) then
            IVNEW = IVNEW + 1
            IVL   = IVNEW
            IVG2IV(IVG) = IVL
            IV2IVG_tmp(IVL) = IVG
            XV_tmp(IVL) = XVG(IVG)
            YV_tmp(IVL) = YVG(IVG)
            ZV_tmp(IVL) = ZVG(IVG)
          else
            IVL = IVG2IV(IVG)
          end if

          IVCELL(ICNEW,IV) = IVL

        end do

      end if
    end do

    NVERT = IVNEW
    allocate(XV(NVERT))
    allocate(YV(NVERT))
    allocate(ZV(NVERT))
    allocate(IV2IVG(NVERT))

    XV(1:NVERT) = XV_tmp(1:NVERT)
    YV(1:NVERT) = YV_tmp(1:NVERT)
    ZV(1:NVERT) = ZV_tmp(1:NVERT)
    IV2IVG(1:NVERT) = IV2IVG_tmp(1:NVERT)

    deallocate(XVG)
    deallocate(YVG)
    deallocate(ZVG)
    deallocate(IC2PROC)
    deallocate(IVGCELLG)
    
    ! ! check a specific local vrt
    ! IV = 1786
    ! print *,'this is rank',rank
    ! IVG = IV2IVG(IV)
    ! print *,'global index',IVG
    ! print *,'coordinate',XV(IV),YV(IV),ZV(IV)
    ! print *,IV2IVG
    
    
    print *,'processor',rank,'has',NCELL,'cells'

    END SUBROUTINE GLOB2LOCAL

    SUBROUTINE READ_BND_DATA
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    INTEGER :: I,IOERR,IINLET,IOUTLET,ISYMP,IWALL,IPRES,IWALM,ICYCL
    INTEGER :: IB,IPATCH,IREGION

    if (rank.eq.0) then
      write(*,*) 'reading boundaries ...'
      I = 0
      IOERR = 0
      open(11,FILE=NAMEBND)
      do while (IOERR.EQ.0)
        read(11,*,iostat=IOERR)
        I=I+1
      end do
      close(11)

      NBOUNDEFINED=I-1
      write(*,'(A,I6)') 'number of global boundaries =',NBOUNDEFINED
    end if
    
    call MPI_BCAST(NBOUNDEFINED,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

    ALLOCATE(IVGBOUN(NBOUNDEFINED,4),IBOUNTYPE(NBOUNDEFINED))

    IINLET = 0
    IOUTLET = 0
    ISYMP = 0
    IWALL = 0
    IPRES = 0
    IWALM = 0
    ICYCL = 0

    open(11,FILE=NAMEBND)
    do I = 1,NBOUNDEFINED
      read(11,*) IB,IVGBOUN(I,1:4),IREGION,&
      IPATCH,IBOUNTYPE(I)
      if (IBOUNTYPE(I).EQ.'INLE') IINLET=IINLET+1
      if (IBOUNTYPE(I).EQ.'OUTL') IOUTLET=IOUTLET+1
      if (IBOUNTYPE(I).EQ.'PRES') IPRES=IPRES+1
      if (IBOUNTYPE(I).EQ.'SYMP') ISYMP=ISYMP+1
      if (IBOUNTYPE(I).EQ.'WALL') IWALL=IWALL+1
      if (IBOUNTYPE(I).EQ.'WALM') IWALM=IWALM+1
      if (IBOUNTYPE(I).EQ.'CYCL') ICYCL=ICYCL+1
    end do
    close(11)
    
    if (rank.eq.0) then
      write(*,*) '--------------------------------'
      write(*,*) 'INFO FROM READ BND DATA'
      write(*,*) 'Number of inlet boundary faces=',IINLET
      write(*,*) 'Number of outlet boundary faces=',IOUTLET
      write(*,*) 'Number of pressure inou boundary faces=',IPRES
      write(*,*) 'Number of symmetry boundary faces=',ISYMP
      write(*,*) 'Number of wall boundary faces=',IWALL
      write(*,*) 'Number of WALLM boundary faces=',IWALM
      write(*,*) 'Number of cyclic boundary faces=',ICYCL
		  write(*,*) 'Number of cyclic pairs=',ICYCL/2
		  write(*,*) '--------------------------------------'
    endif
        

    END SUBROUTINE READ_BND_DATA