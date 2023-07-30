    SUBROUTINE calcjacob_transfinite
    use setup3d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: ic,k,is,js,ks,ifp,jfp,kfp
    double precision :: xi,eta,beta,xxs(3,8),coor_der(3,3)
    double precision :: xxi,yxi,zxi,xeta,yeta,zeta,xbeta,ybeta,zbeta

    DO ic = 1,NCELL

      DO k = 1,8
        xxs(1,k) = XV(IVCELL(ic,k))
        xxs(2,k) = YV(IVCELL(ic,k))
        xxs(3,k) = ZV(IVCELL(ic,k))
      END DO

      DO ks = 1,N
      DO js = 1,N
      DO is = 1,N

        xi = Xs(is)
        eta = Xs(js)
        beta = Xs(ks)

        CALL Deri_transfinite(xxs,xi,eta,beta,coor_der)

        xxi = coor_der(1,1)
        yxi = coor_der(2,1)
        zxi = coor_der(3,1)

        xeta = coor_der(1,2)
        yeta = coor_der(2,2)
        zeta = coor_der(3,2)

        xbeta = coor_der(1,3)
        ybeta = coor_der(2,3)
        zbeta = coor_der(3,3)

        Jac(is,js,ks,ic) = xxi*yeta*zbeta + xeta*ybeta*zxi+xbeta*yxi*zeta  &
        - xxi*ybeta*zeta - xeta*yxi*zbeta - xbeta*yeta*zxi

      END DO
      END DO
      END DO

      ! S vector :: elements of inverse Jacobian

      ! loop to find the S vector at the flux points (family 1)
      DO kfp = 1,N
      DO jfp = 1,N
      DO ifp = 1,N+1

        xi = Xf(ifp)
        eta = Xs(jfp)
        beta = Xs(kfp)

        CALL Deri_transfinite(xxs,xi,eta,beta,coor_der)

        xxi = coor_der(1,1)
        yxi = coor_der(2,1)
        zxi = coor_der(3,1)

        xeta = coor_der(1,2)
        yeta = coor_der(2,2)
        zeta = coor_der(3,2)

        xbeta = coor_der(1,3)
        ybeta = coor_der(2,3)
        zbeta = coor_der(3,3)

        S1(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
        S1(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
        S1(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta

        S1(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
        S1(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
        S1(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta

        S1(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
        S1(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
        S1(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi

      END DO
      END DO
      END DO

      ! loop to find the S vector at the flux points (family 2)
      DO kfp = 1,N
      DO jfp = 1,N+1
      DO ifp = 1,N

        xi = Xs(ifp)
        eta = Xf(jfp)
        beta = Xs(kfp)

        CALL Deri_transfinite(xxs,xi,eta,beta,coor_der)

        xxi = coor_der(1,1)
        yxi = coor_der(2,1)
        zxi = coor_der(3,1)

        xeta = coor_der(1,2)
        yeta = coor_der(2,2)
        zeta = coor_der(3,2)

        xbeta = coor_der(1,3)
        ybeta = coor_der(2,3)
        zbeta = coor_der(3,3)

        S2(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
        S2(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
        S2(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta

        S2(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
        S2(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
        S2(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta

        S2(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
        S2(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
        S2(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi

      END DO
      END DO
      END DO

      ! loop to find the S vector at the flux points (family 3)
      DO kfp = 1,N+1
      DO jfp = 1,N
      DO ifp = 1,N

        xi = Xs(ifp)
        eta = Xs(jfp)
        beta = Xf(kfp)

        CALL Deri_transfinite(xxs,xi,eta,beta,coor_der)

        xxi = coor_der(1,1)
        yxi = coor_der(2,1)
        zxi = coor_der(3,1)

        xeta = coor_der(1,2)
        yeta = coor_der(2,2)
        zeta = coor_der(3,2)

        xbeta = coor_der(1,3)
        ybeta = coor_der(2,3)
        zbeta = coor_der(3,3)

        S3(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
        S3(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
        S3(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta

        S3(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
        S3(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
        S3(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta

        S3(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
        S3(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
        S3(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi

      END DO
      END DO
      END DO

    END DO

    END SUBROUTINE calcjacob_transfinite