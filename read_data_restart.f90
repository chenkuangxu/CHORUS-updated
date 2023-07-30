    SUBROUTINE READ_ALL_DATA_BINARY_MPI

    use setup3d
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    integer :: I,is,js,k,ICG,IC,fh
    integer :: status(MPI_STATUS_SIZE)
    double precision,allocatable :: tmp(:,:,:,:)


    !open the restart file...
    call MPI_FILE_OPEN(MPI_COMM_WORLD,NAMERESTART,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)

    !read the header...
    call MPI_FILE_READ(fh,iter,1,MPI_INTEGER,status,ierror)
    call MPI_FILE_READ(fh,ctime,1,MPI_DOUBLE_PRECISION,status,ierror)

    !now read cell no. and solutions...
    allocate(tmp(numv,N,N,N))

    do I=1,NCELLGLOB
        if (rank.eq.0.and.mod(I,500).eq.0) then
            print *,'loading restart data already completed',100.d0*I/NCELLGLOB,'%'
        end if
        call MPI_FILE_READ(fh,ICG,1,MPI_INTEGER,status,ierror)
        call MPI_FILE_READ(fh,tmp,numv*N*N*N,MPI_DOUBLE_PRECISION,status,ierror)

        IC=ICG2IC(ICG)
        if(IC .ne. 0) Q(1:numv,1:N,1:N,1:N,IC)=tmp(1:numv,1:N,1:N,1:N)
    end do

    call MPI_FILE_CLOSE(fh,ierror)

    deallocate(tmp)

    END SUBROUTINE READ_ALL_DATA_BINARY_MPI


    SUBROUTINE READ_ALL_DATA_BINARY

    use setup3d
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    integer :: I,is,js,k,ICG,IC,fh
    integer :: status(MPI_STATUS_SIZE)
    double precision,allocatable :: tmp(:,:,:,:)

    open(newunit=fh,file=NAMERESTART,access='stream')

    read(fh) iter
    read(fh) ctime

    !now read cell no. and solutions...
    allocate(tmp(numv,N,N,N))

    do I=1, NCELLGLOB
   
        if (rank.eq.0.and.mod(I,5000).eq.0) then
            print *,'loading restart data already completed',100.d0*I/NCELLGLOB,'%'
        end if

        read(fh) ICG
        read(fh) tmp(1:numv,1:N,1:N,1:N)

        IC=ICG2IC(ICG)

        if(IC.ne.0) Q(1:numv,1:N,1:N,1:N,IC)=tmp(1:numv,1:N,1:N,1:N)
    end do

    close(fh)

    deallocate(tmp)

    END SUBROUTINE READ_ALL_DATA_BINARY

    SUBROUTINE READ_ALL_DATA_BINARY_DIFF

    use setup3d
    IMPLICIT NONE
    INCLUDE 'mpif.h'
 
    integer :: I,is,js,ks,k,ICG,IC,fh,isr,jsr,ksr
    integer :: status(MPI_STATUS_SIZE)
    double precision :: xsi,eta,bet
    double precision,allocatable :: tmp(:,:,:,:)
    double precision,allocatable :: Qrestart(:,:,:,:,:)

    allocate(Qrestart(numv,restart_ord,restart_ord,restart_ord,NCELL))
 
    open(newunit=fh,file=NAMERESTART,access='stream')
 
    read(fh) iter
    read(fh) ctime
 
    !now read cell no. and solutions...
    allocate(tmp(numv,restart_ord,restart_ord,restart_ord))
 
    do I=1, NCELLGLOB

        if (rank.eq.0.and.mod(I,5000).eq.0) then
            print *,'loading restart data already completed',100.d0*I/NCELLGLOB,'%'
        end if
 
        read(fh) ICG
        read(fh) tmp(1:numv,1:restart_ord,1:restart_ord,1:restart_ord)
 
        IC=ICG2IC(ICG)
 
        if(IC.ne.0) Qrestart(1:numv,1:restart_ord,1:restart_ord,1:restart_ord,&
        IC)=tmp(1:numv,1:restart_ord,1:restart_ord,1:restart_ord)
    end do
 
    close(fh)
 
    deallocate(tmp)

    do IC = 1,NCELL
        do ks = 1,N
        do js = 1,N
        do is = 1,N
            Q(1:numv,is,js,ks,IC) = 0.d0
            xsi = Xs(is)
            eta = Xs(js)
            bet = Xs(ks)
            do ksr = 1,restart_ord
            do jsr = 1,restart_ord
            do isr = 1,restart_ord
                Q(1:numv,is,js,ks,IC) = Q(1:numv,is,js,ks,IC) + &
                Qrestart(1:numv,isr,jsr,ksr,IC)*hval_r(restart_ord,isr,xsi)*&
                hval_r(restart_ord,jsr,eta)*hval_r(restart_ord,ksr,bet)
            end do
            end do
            end do
        end do
        end do
        end do
    end do

    deallocate(Qrestart)

    contains

    DOUBLE PRECISION FUNCTION hval_r(np,i,xval)
	
	use setup3d
	IMPLICIT NONE
	
	integer, intent(in)	:: np, i
	double precision, intent(in)	:: xval
	double precision	:: hvaln, hvald
	integer	:: s
		
	hvaln = 1.d0
	hvald = 1.d0
	do s=1,np
		if( s/=i) then
			hvaln = hvaln * ( xval - Xs_r(s) ) 
			hvald = hvald * ( Xs_r(i) - Xs_r(s) )
		end if		
	end do
	
	hval_r = hvaln/hvald
			
	END FUNCTION hval_r		
 
    END SUBROUTINE READ_ALL_DATA_BINARY_DIFF
