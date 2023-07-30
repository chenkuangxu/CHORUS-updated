    SUBROUTINE SETINITIALCOND

    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: ic,is,js,ks,istart,i
    double precision :: radcomp,dist,distm,distp,Ttemp,rhotemp

    if (rank.eq.0) print *,'start setting initial conditions ...'

    istart = 0

    DO ic = 1,NCELL

        ! if (mod(ic,2000).EQ.0.AND.rank.eq.0) print *,'having completed',&
        ! floor(dble(ic)*100/NCELL),'%'

        DO ks = 1,N
        DO js = 1,N
        DO is = 1,N

            radcomp = sqrt(XXsolu(1,is,js,ks,ic)**2+&
            XXsolu(2,is,js,ks,ic)**2+&
            XXsolu(3,is,js,ks,ic)**2)
        
            CALL dichotomy_search(radcomp,istart)

            dist  = rIn(istart+1) - rIn(istart)
            distm = (radcomp - rIn(istart))/dist
            distp =  (rIn(istart+1) - radcomp)/dist

            rhotemp = rhoIn(istart)*distp+rhoIn(istart+1)*distm

            Q(1,is,js,ks,ic) = rhotemp

            Ttemp = TIn(istart)*distp+TIn(istart+1)*distm

            Q(5,is,js,ks,ic) = rhotemp*R_const*Ttemp/(gam-1)

            Q(2:4,is,js,ks,ic) = 0.d0
      
        END DO
        END DO
        END DO
    
    END DO

    if (rank.eq.0) print *,'completed setting initial conditions ...'

    iter = 0
    ctime = 1d-6

    END SUBROUTINE SETINITIALCOND

    SUBROUTINE dichotomy_search(radius,I)
    use setup3d
    IMPLICIT NONE
        
    double precision,intent(in) :: radius
    integer,intent(out) :: I
    integer :: int_a,int_b,int_mid
    
    int_a = 1
    int_b = MAXSM
    
    DO WHILE (int_b-int_a.NE.1)
        int_mid = (int_a+int_b)/2
        if ((radius.LE.rIn(int_mid)).AND.(radius.GT.rIn(int_a))) then
            int_b = int_mid
        else if ((radius.LE.rIn(int_b)).AND.(radius.GT.rIn(int_mid))) then
            int_a = int_mid
        else
            print *,'radius does not fall into any interval,WRONG!!'
        end if
    END DO
    
    I = int_a
    
    RETURN
    END SUBROUTINE dichotomy_search