    MODULE setup3d
    IMPLICIT NONE

    !-----------------MPI related variable------------------
    integer :: size,rank,ierror,NPROC

    !-----------------------FILE NAME------------------------
    character(64) :: NAMECEL,NAMEVRT,NAMEBND,METIS_CELL,NAMEVRT20,SOLARMODEL

    !-------------Connectivity and grid variables------------
    ! global
    integer :: NCELLGLOB,NVERTGLOB,NBOUNDEFINED
    integer,allocatable :: IVGCELLG(:,:),IVGBOUN(:,:)
    character(8),allocatable :: IBOUNTYPE(:)
    double precision,allocatable :: XVG(:),YVG(:),ZVG(:)
   
    ! local
    integer :: NCELL,NVERT,NFACE,NCELLTOT,NEDGE,NBOUN
    integer :: NINTER,NINLET,NOUTLET,NSYMP,NWALL,NPROCINT,NPRES,NWALM,NCYCLIC,&
    NCYCLOC,NCYCREM
    integer :: MAXPROCINT,MAXCYCREM
    integer,allocatable :: IBFINTER(:)
    integer,allocatable :: IVCELL(:,:),IC2F(:,:),IF2V(:,:),IF2C(:,:),IC2E(:,:)
    integer,allocatable :: ICVERT(:),ICVSTA(:)
    INTEGER,DIMENSION(6,4) :: IVFACE
    integer,allocatable :: IC2PROC(:)
    integer,allocatable :: ICG2IC(:),IC2ICG(:),IV2IVG(:),IVG2IV(:)
    double precision,allocatable :: XV(:),YV(:),ZV(:)
    double precision,allocatable :: xxf_cell(:,:,:)
    integer,allocatable :: PROCSBUF(:),PROCRBUF(:,:),PROCINT2PROC(:),PROCINT2F_PROC(:),&
    CYCREM2PROC(:),CYCREM2F_PROC(:)
    double precision,allocatable :: CYCRSBUF(:),CYCRRBUF(:,:)
    integer :: NPAIR
    integer,allocatable :: CYC2FPAIR(:,:),mapcl(:,:,:,:)

    !-----------conserved variables and flux variables-------------
    integer,parameter :: numv = 5
    double precision,allocatable,dimension(:,:,:,:,:) :: Q,F1,G2,H3,Fv1,Gv2,Hv3,resid
    double precision,allocatable,dimension(:,:,:,:,:) :: Qvfi,Qvfj,Qvfk
    double precision,allocatable,dimension(:,:,:,:,:,:) :: nablaQs,nablaQvfi,nablaQvfj,nablaQvfk
    double precision,allocatable,dimension(:,:,:,:) :: dedxi,dedeta,dedbeta

    !---------------------coordinate related--------------------
    double precision,allocatable,dimension(:,:,:,:,:)::XXsolu,XXfluxi,XXfluxj,XXfluxk
    double precision,allocatable,dimension(:,:,:,:) :: Jac
    double precision,allocatable,dimension(:,:,:,:,:,:) :: S1,S2,S3
    double precision,allocatable,dimension(:,:,:,:,:) :: dmdxs,dmdxf1,dmdxf2,dmdxf3
    integer :: CURVE_WALL

    !----------------------parallel computing---------------------------
    !! move RBUF,RBUF3C,Qflc,Qfrc,Qfl2c,Qfr2c,Qfl3c,Qfr3c to local subroutines

    !---------------------SD Method related---------------------
    integer :: N
    double precision,allocatable,dimension(:) :: Xs,Xf
    double precision,allocatable,dimension(:,:) :: Lmat,Mmat

    !-----------------------Time Marching-----------------------
    double precision :: ctime,dt
    integer :: k_stage,MAXITER,iter

    !---------------------MAP POINT-----------------------------
    integer,allocatable :: Gfp2Lfp(:,:,:,:),map1(:,:,:,:),map2(:,:,:,:)

    !--------------------CONNECT BOUNDARY-----------------------
    integer,allocatable :: BOUNFACEINDEX(:)
    integer,allocatable :: IBFINL(:),IBFOUT(:),IBFSYMP(:),IBFWAL(:),IBFPROC(:),IF2IBPROC(:),&
    IBFPRES(:),IBFWALM(:),IBFCYC(:),IBFCYCLOC(:),IBFCYCREM(:),IF2IBCYCREM(:)

    !--------------------Gauss Integration----------------------
    double precision,allocatable :: G_omega(:),G_psi(:)
    integer :: int_num
    double precision,allocatable :: JacGa(:,:,:,:)
    double precision :: kE_g

    !------------------------BND KIND---------------------------
    double precision,parameter :: tolCYC = 1d-6

    !--------------------const parameter------------------------
    double precision,parameter :: pi = 3.141592653589793238d0
    double precision,parameter :: lambda = -2.d0/3.d0
    double precision,parameter :: prandt = 0.72d0
    
    !--------------------inflow parameter-----------------------
    double precision :: rinf,pinf,uinf,vinf,winf

    !---------------------viscous related-------------------------
    integer :: vismode
    double precision :: mu, eta0, kappa_s

    !-------------------divergence cleaning--------------------
    double precision :: c_h,alpha,umax,umax_g,eigvmax,eigvmax_g

    !------------------------restart------------------------------
    integer,allocatable :: PROC_NCELL(:)
    integer :: restart,restart_ord
    double precision,allocatable,dimension(:) :: Xs_r,Xf_r
    character(64) :: NAMERESTART

    !----------------------post_processing------------------------
    integer :: NRE,nwrite,plotnodes,plotcells
    integer,allocatable :: gnumcell(:,:,:,:),connec_c2n(:,:),&
    gnode_factor(:),gnumnode(:,:,:,:),isfacenod(:),&
    gfacenode(:,:,:),isedgenod(:),gedgenode(:,:),PROC_plotnodes(:),&
    PROC_plotcells(:)
    double precision,allocatable :: plotX(:,:),plotQ(:,:),mach3d(:)

    !!! does not need cyclic faces in the solar dynamo
    double precision,parameter :: DXCYCL = 10d5
    double precision,parameter :: DYCYCL = 10d5
    double precision,parameter :: DZCYCL = 10d5

    !---------------------INITIAL CONDITION---------------------
    double precision,allocatable :: rIn(:),rhoIn(:),TIn(:),radcIn(:),dTIn(:)
    integer :: MAXSM
    double precision,parameter :: sunM = 1.98891d33
    double precision,parameter :: sunG = 6.67d-8
    double precision,parameter :: omega_rigid = 8.1d-5 
    double precision,parameter :: gam = 1.6666666666666666666666666666666666667d0
    double precision,parameter :: RINLET = 487.d0*1.0d8
    double precision,parameter :: ROUTLET = 661.d0*1.0d8
    double precision,parameter :: constL = -3.846d36
    double precision :: constheati = -1.d0
    double precision,parameter :: R_const = 1.4d8
    double precision,parameter :: Cp = 3.5d8
    double precision,allocatable,dimension(:,:,:,:) :: gR,radCfi,radCfj,radCfk

    END MODULE setup3d


