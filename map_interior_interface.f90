    SUBROUTINE MAPINTERFACE
    
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    INTEGER :: NF,IC1,IC2,IV1,IV2,IV3,IV4
    INTEGER :: IVF2,IVF4
    INTEGER :: IC1V1,IC1V2,IC1V3,IC1V4
    INTEGER :: IC1V5,IC1V6,IC1V7,IC1V8
    INTEGER :: IC2V1,IC2V2,IC2V3,IC2V4
    INTEGER :: IC2V5,IC2V6,IC2V7,IC2V8

    DO NF=1,NFACE
      IF2C(NF,3)=0
      IF2C(NF,4)=0
      IF2C(NF,5)=0
      IF2C(NF,6)=0
    END DO

    DO NF=1,NFACE

      IC1=IF2C(NF,1)
      IC2=IF2C(NF,2)

      IV1= IF2V(NF,1)
      IVF2= IF2V(NF,2)
      IV3= IF2V(NF,3)
      IVF4= IF2V(NF,4)

      IC1V1= IVCELL(IC1,1)
      IC1V2= IVCELL(IC1,2)
      IC1V3= IVCELL(IC1,3)
      IC1V4= IVCELL(IC1,4)

      IC1V5= IVCELL(IC1,5)
      IC1V6= IVCELL(IC1,6)
      IC1V7= IVCELL(IC1,7)
      IC1V8= IVCELL(IC1,8)

      ! Corresponding to face 1, with vertices 1-2-3-4
      IF ((IC1V1.EQ.IV1.AND.IC1V3.EQ.IV3)) THEN 
        IF2C(NF,3)=1
        IF2C(NF,5)=1
        IV2 =IVF2
        IV4 =IVF4
      ! Corresponding to face 2 with vertices 5-6-7-8 
      ELSE IF ((IC1V5.EQ.IV1.AND.IC1V7.EQ.IV3)) THEN
        IF2C(NF,3)=2
        IF2C(NF,5)=1
        IV2 =IVF4
        IV4 =IVF2
      ! Corresponding to face 3 with vertices 5-6-2-1 
      ELSE IF ((IC1V5.EQ.IV1.AND.IC1V2.EQ.IV3)) THEN
        IF2C(NF,3)=3
        IF2C(NF,5)=1
        IV2 =IVF2
        IV4 =IVF4
      ! Corresponding to face 4 with vertices 8-7-3-4 
      ELSE IF ((IC1V8.EQ.IV1.AND.IC1V3.EQ.IV3)) THEN
        IF2C(NF,3)=4
        IF2C(NF,5)=1
        IV2 =IVF4
        IV4 =IVF2
      ! Corresponding to face 5 with vertices 5-1-4-8 
      ELSE IF ((IC1V5.EQ.IV1.AND.IC1V4.EQ.IV3)) THEN
        IF2C(NF,3)=5
        IF2C(NF,5)=1
        IV2 =IVF2
        IV4 =IVF4
      ! Corresponding to face 6 with vertices 6-2-3-7 
      ELSE IF ((IC1V6.EQ.IV1.AND.IC1V3.EQ.IV3)) THEN
        IF2C(NF,3)=6
        IF2C(NF,5)=1
        IV2 =IVF4
        IV4 =IVF2
      ELSE 
        write(*,*) 'Orientation problem left cell no.=', IC1
        write(*,*) '8 v points=', IC1V1,IC1V2,IC1V3,IC1V4,IC1V5,IC1V6,IC1V7,IC1V8 
        write(*,*) '4 v points=', IV1,IV2,IV3,IV4
      END IF

! the right cell each face has four different orientations 
! The reference face indices use the left cells' first face 
! The reference frame uses IV1,IV2,IV3,IV4

      IF (IC2.LE.NCELL+NPROCINT+NCYCREM.AND.IC2.GT.0) THEN
      !  NEEDS TO HAVE IF2C(NF,4) for NPROCINT faces (to be used in procintflux)

        IC2V1= IVCELL(IC2,1)
        IC2V2= IVCELL(IC2,2)
        IC2V3= IVCELL(IC2,3)
        IC2V4= IVCELL(IC2,4)
!
        IC2V5= IVCELL(IC2,5)
        IC2V6= IVCELL(IC2,6)
        IC2V7= IVCELL(IC2,7)
        IC2V8= IVCELL(IC2,8)

!       WRITE (*,*) '8 v points=', IC2V1,IC2V2,IC2V3,IC2V4,IC2V5,IC2V6,IC2V7,IC2V8 
! case 1

        IF (IC2V1.EQ.IV1.AND.IC2V3.EQ.IV3)  THEN
          IF2C(NF,4)=1
          IF2C(NF,6) = 1
        ELSE IF ((IC2V2.EQ.IV1.AND.IC2V4.EQ.IV3)) then
          IF2C(NF,4)=1
          IF2C(NF,6) = 2
        ELSE IF ( (IC2V1.EQ.IV3.AND.IC2V3.EQ.IV1)) then 
          IF2C(NF,4)=1
          IF2C(NF,6) = 3
        ELSE IF ((IC2V4.EQ.IV1.AND.IC2V2.EQ.IV3)) then
          IF2C(NF,4)=1
          IF2C(NF,6) = 4
! case 2
! v 5 and v7 for face 2
        ELSE IF ((IC2V5.EQ.IV1.AND.IC2V7.EQ.IV3)) THEN
          IF2C(NF,4)=2
          IF2C(NF,6)=1
        ELSE IF ((IC2V6.EQ.IV1.AND.IC2V8.EQ.IV3)) THEN 
          IF2C(NF,4)=2
          IF2C(NF,6)=2
        ELSE IF ((IC2V5.EQ.IV3.AND.IC2V7.EQ.IV1)) THEN 
          IF2C(NF,4)=2
          IF2C(NF,6)=3
        ELSE IF ((IC2V6.EQ.IV3.AND.IC2V8.EQ.IV1)) THEN 
          IF2C(NF,4)=2
          IF2C(NF,6)=4
          !write(*,*) 'not case 2'
! Case 3
! v 5 and v2 for face 3
        ELSE IF ((IC2V5.EQ.IV1.AND.IC2V2.EQ.IV3)) THEN
          IF2C(NF,4)=3
          IF2C(NF,6)=1 
        Else IF (IC2V1.EQ.IV3.AND.IC2V6.EQ.IV1) then
          IF2C(NF,4)=3
          IF2C(NF,6)=2 
        Else IF (IC2V5.EQ.IV3.AND.IC2V2.EQ.IV1) then
          IF2C(NF,4)=3
          IF2C(NF,6)=3 
        Else IF (IC2V6.EQ.IV3.AND.IC2V1.EQ.IV1) then
          IF2C(NF,4)=3
          IF2C(NF,6)=4 
! write(*,*) 'not case 3'
! case 4
! v 8 and v3 for face 4
        ELSE IF ((IC2V8.EQ.IV1.AND.IC2V3.EQ.IV3) ) THEN 
          IF2C(NF,4)=4
          IF2C(NF,6)=1
        ELSE IF (IC2V4.EQ.IV3.AND.IC2V7.EQ.IV1) then
          IF2C(NF,4)=4
          IF2C(NF,6)=2
        ELSE IF (IC2V8.EQ.IV3.AND.IC2V3.EQ.IV1) then
          IF2C(NF,4)=4
          IF2C(NF,6)=3
        ELSE IF (IC2V7.EQ.IV3.AND.IC2V4.EQ.IV1) then
          IF2C(NF,4)=4
          IF2C(NF,6)=4
!  write(*,*) 'not case 4'
!  v 5 and v 4 for face5
        ELSE IF (IC2V5.EQ.IV1.AND.IC2V4.EQ.IV3) THEN
          IF2C(NF,4)=5
          IF2C(NF,6)=1
        ELSE IF (IC2V8.EQ.IV3.AND.IC2V1.EQ.IV1) THEN
          IF2C(NF,4)=5
          IF2C(NF,6)=2
        ELSE IF (IC2V5.EQ.IV3.AND.IC2V4.EQ.IV1) THEN
          IF2C(NF,4)=5
          IF2C(NF,6)=3
        ELSE IF (IC2V1.EQ.IV3.AND.IC2V8.EQ.IV1) THEN
          IF2C(NF,4)=5
          IF2C(NF,6)=4
! write(*,*) 'not case 5'
! v 6 and v 3 for face 6
        ELSE IF ((IC2V6.EQ.IV1.AND.IC2V3.EQ.IV3)) THEN
          IF2C(NF,4)=6
          IF2C(NF,6)=1
        ELSE IF (IC2V7.EQ.IV3.AND.IC2V2.EQ.IV1) THEN
          IF2C(NF,4)=6
          IF2C(NF,6)=2
        ELSE IF (IC2V6.EQ.IV3.AND.IC2V3.EQ.IV1) THEN
          IF2C(NF,4)=6
          IF2C(NF,6)=3
        ELSE IF (IC2V2.EQ.IV3.AND.IC2V7.EQ.IV1) THEN
          IF2C(NF,4)=6
          IF2C(NF,6)=4
        ELSE 
          write(*,*) 'LEFT CELL IC2=',IC2
          WRITE(*,*) 'Rank = ',rank
          write (*,*) 'left cell f', IF2C(NF,3)
          write(*,*) '4 v points=', IV1,IVF2,IV3,IVF4
          write(*,*) 'reference 4 v points=', IV1,IV2,IV3,IV4
          write(*,*) '8 v points=', IC2V1,IC2V2,IC2V3,IC2V4,IC2V5,IC2V6,IC2V7,IC2V8 
        END IF
       ! write(*,*) 'not case 6'

      END IF ! IF IC2 .LE.NCELL+NPROCINT+NCYCREM

    ENDDO ! Loop over faces

    END SUBROUTINE MAPINTERFACE


    SUBROUTINE MAPFPOINT
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: mi,ni,k,IFA
    integer :: IC1,IC2,C1F,C2F
    integer,allocatable :: IBFINTER_tmp(:)

    NINTER = 0
    allocate(IBFINTER_tmp(NFACE))

    if(rank.eq.0) WRITE(*,*) '------ START OF MAPFPOINT------'

    DO k = 1,NFACE
      IC2 = IF2C(k,2)

      IF (IC2.LE.NCELL.AND.IC2.GT.0) THEN
        NINTER = NINTER + 1
        IBFINTER_tmp(NINTER) = k
      END IF
    END DO

    allocate(IBFINTER(NINTER))
    IBFINTER(1:NINTER) = IBFINTER_tmp(1:NINTER)
    deallocate(IBFINTER_tmp)

    DO k = 1,NFACE
      IC1 = IF2C(k,1)
      IC2 = IF2C(k,2)

      C1F = IF2C(k,3)
      C2F = IF2C(k,4)

      DO ni = 1,N
      DO mi = 1,N

        !1-2 STORES CELL NO.
        Gfp2Lfp(1,mi,ni,k) =IF2C(k,1)
        Gfp2Lfp(2,mi,ni,k) =IF2C(k,2)

        IF (IC2.LE.NCELL.AND.IC2.GT.0) THEN
!          Here we still need to compute the orientation to be used in procintflux.f90
        if(C2F==1) then
          Gfp2Lfp(4,mi,ni,k)=3
          Gfp2Lfp(10,mi,ni,k) = 1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(8,mi,ni,k) = ni 
            Gfp2Lfp(9,mi,ni,k) = mi 
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(8,mi,ni,k) = N-mi+1 
            Gfp2Lfp(9,mi,ni,k) = ni 
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(8,mi,ni,k) = N-ni+1 
            Gfp2Lfp(9,mi,ni,k) = N-mi+1 
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(8,mi,ni,k) = mi 
            Gfp2Lfp(9,mi,ni,k) = N-ni+1 
          end if
 !         end of face 1 for right cell
        else if(C2F==2) then
          Gfp2Lfp(4,mi,ni,k)=3
          Gfp2Lfp(10,mi,ni,k) = N+1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(8,mi,ni,k) = mi 
            Gfp2Lfp(9,mi,ni,k) = ni 
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(8,mi,ni,k) = N-ni+1 
            Gfp2Lfp(9,mi,ni,k) = mi 
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(8,mi,ni,k) = N-mi+1 
            Gfp2Lfp(9,mi,ni,k) = N-ni +1
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(8,mi,ni,k) = ni 
            Gfp2Lfp(9,mi,ni,k) = N-mi+1 
          end if
!         ! end of face 2 for right cell
        else if(C2F==3) then
          Gfp2Lfp(4,mi,ni,k)=2
          Gfp2Lfp(9,mi,ni,k) = 1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(8,mi,ni,k) = ni 
            Gfp2Lfp(10,mi,ni,k) = N-mi +1
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(8,mi,ni,k) = N-mi+1 
            Gfp2Lfp(10,mi,ni,k) = N-ni +1
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(8,mi,ni,k) = N-ni+1 
            Gfp2Lfp(10,mi,ni,k) =mi 
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(8,mi,ni,k) = mi 
            Gfp2Lfp(10,mi,ni,k) = ni 
          end if
!         ! end of face 3 for right cell

        else if(C2F==4) then
          Gfp2Lfp(4,mi,ni,k) = 2
          Gfp2Lfp(9,mi,ni,k) = N+1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(8,mi,ni,k) = mi 
            Gfp2Lfp(10,mi,ni,k) = N-ni +1
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(8,mi,ni,k) = N-ni+1 
            Gfp2Lfp(10,mi,ni,k) = N-mi +1
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(8,mi,ni,k) = N-mi+1 
            Gfp2Lfp(10,mi,ni,k) =ni 
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(8,mi,ni,k) = ni 
            Gfp2Lfp(10,mi,ni,k) = mi 
          end if

        else if(C2F==5) then
          Gfp2Lfp(4,mi,ni,k) = 1
          Gfp2Lfp(8,mi,ni,k) = 1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(9,mi,ni,k) = mi 
            Gfp2Lfp(10,mi,ni,k) = N-ni +1
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(9,mi,ni,k) = ni 
            Gfp2Lfp(10,mi,ni,k) =mi 
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(9,mi,ni,k) = N-mi+1 
            Gfp2Lfp(10,mi,ni,k) =ni 
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(9,mi,ni,k) = N-ni +1
            Gfp2Lfp(10,mi,ni,k) = N-mi+1 
          end if

        else if(C2F==6) then
          Gfp2Lfp(4,mi,ni,k) = 1
          Gfp2Lfp(8,mi,ni,k) = N+1 
          if(IF2C(k,6)==1)  then
            Gfp2Lfp(9,mi,ni,k) = ni 
            Gfp2Lfp(10,mi,ni,k) = N-mi +1
          else if(IF2C(k,6)==2)  then
            Gfp2Lfp(9,mi,ni,k) = mi 
            Gfp2Lfp(10,mi,ni,k) =ni 
          else if(IF2C(k,6)==3)  then
            Gfp2Lfp(9,mi,ni,k) = N-ni+1 
            Gfp2Lfp(10,mi,ni,k) =mi 
          else if(IF2C(k,6)==4)  then
            Gfp2Lfp(9,mi,ni,k) = N-mi +1
            Gfp2Lfp(10,mi,ni,k) = N-ni+1 
          end if
    
        else 
          write(*,*) 'error'
         
        end if  ! end of face 6 for right cell
        
        END IF

! family xi, eta and beta
        if (C1F==1) then 
          Gfp2Lfp(3,mi,ni,k)=3
          Gfp2Lfp(5,mi,ni,k) = mi 
          Gfp2Lfp(6,mi,ni,k) = ni 
          Gfp2Lfp(7,mi,ni,k) = 1 
        else if (C1F==2) then 
          Gfp2Lfp(3,mi,ni,k)=3
          Gfp2Lfp(5,mi,ni,k) = ni 
          Gfp2Lfp(6,mi,ni,k) = mi 
          Gfp2Lfp(7,mi,ni,k) = N+1 
        else if (C1F==3) then
          Gfp2Lfp(3,mi,ni,k)=2
          Gfp2Lfp(5,mi,ni,k) = mi 
          Gfp2Lfp(6,mi,ni,k) = 1 
          Gfp2Lfp(7,mi,ni,k) = N-ni+1 
        else if (C1F==4) then 
          Gfp2Lfp(3,mi,ni,k)=2 
          Gfp2Lfp(5,mi,ni,k) = ni 
          Gfp2Lfp(6,mi,ni,k) = N+1 
          Gfp2Lfp(7,mi,ni,k) = N-mi +1
        else if (C1F==5) then 
          Gfp2Lfp(3,mi,ni,k)=1
          Gfp2Lfp(5,mi,ni,k) = 1 
          Gfp2Lfp(6,mi,ni,k) = ni 
          Gfp2Lfp(7,mi,ni,k) = N-mi+1 
        else  if (C1F==6)  then
          Gfp2Lfp(3,mi,ni,k)=1
          Gfp2Lfp(5,mi,ni,k) = N+1 
          Gfp2Lfp(6,mi,ni,k) = mi 
          Gfp2Lfp(7,mi,ni,k) = N-ni+1 
        else 
          write(*,*) 'error'
        end if

      END DO
      END DO

    END DO
    
    END SUBROUTINE MAPFPOINT
