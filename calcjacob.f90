    SUBROUTINE calcjacob
    use setup3d
    implicit none
    include 'mpif.h'

    integer :: ic,is,js,ks,ifp,jfp,kfp,k,kmax
    double precision :: xi,eta,beta
    double precision,allocatable :: xxs(:,:)
    double precision :: xxi,yxi,zxi,xeta,yeta,zeta,xbeta,ybeta,zbeta

    if (CURVE_WALL==1) then
        allocate(xxs(3,20))
        kmax = 20
    else
        allocate(xxs(3,8))
        kmax = 8
    end if

    do ic = 1,NCELL
        
        if (CURVE_WALL==1) then
            do k = 1,20
                xxs(1,k) = xxf_cell(1,k,ic)
                xxs(2,k) = xxf_cell(2,k,ic)
                xxs(3,k) = xxf_cell(3,k,ic)
            end do
        else
            do k = 1,8
                xxs(1,k) = XV(IVCELL(ic,k))
                xxs(2,k) = YV(IVCELL(ic,k))
                xxs(3,k) = ZV(IVCELL(ic,k))
            end do
        end if

        do ks = 1,N
        do js = 1,N
        do is = 1,N
            xxi =   0.d0
            xeta =  0.d0
            xbeta = 0.d0

            yxi =   0.d0
            yeta =  0.d0
            ybeta = 0.d0

            zxi =   0.d0
            zeta =  0.d0
            zbeta = 0.d0

            do k = 1,kmax
                xxi = xxi + dmdxs(k,1,is,js,ks) * xxs(1,k)
			    yxi = yxi + dmdxs(k,1,is,js,ks) * xxs(2,k)
			    zxi = zxi + dmdxs(k,1,is,js,ks) * xxs(3,k)

			    xeta = xeta + dmdxs(k,2,is,js,ks) * xxs(1,k)
                yeta = yeta + dmdxs(k,2,is,js,ks) * xxs(2,k)
			    zeta = zeta + dmdxs(k,2,is,js,ks) * xxs(3,k)

			    xbeta = xbeta + dmdxs(k,3,is,js,ks) * xxs(1,k)
			    ybeta = ybeta + dmdxs(k,3,is,js,ks) * xxs(2,k)
                zbeta = zbeta + dmdxs(k,3,is,js,ks) * xxs(3,k)
            end do
            
            Jac(is,js,ks,ic) = xxi*yeta*zbeta + xeta*ybeta*zxi+xbeta*yxi*zeta &
            - xxi*ybeta*zeta - xeta*yxi*zbeta - xbeta*yeta*zxi

        end do
        end do
        end do

        ! loop to find the S vector at the flux points (family 1)
        do kfp = 1,N
        do jfp = 1,N
        do ifp = 1,N+1
            xxi =   0.d0
            xeta =  0.d0
            xbeta = 0.d0

            yxi =   0.d0
            yeta =  0.d0
            ybeta = 0.d0

            zxi =   0.d0
            zeta =  0.d0
            zbeta = 0.d0

            do k = 1,kmax
                xxi = xxi + dmdxf1(k,1,ifp,jfp,kfp) * xxs(1,k)
			    yxi = yxi + dmdxf1(k,1,ifp,jfp,kfp) * xxs(2,k)
			    zxi = zxi + dmdxf1(k,1,ifp,jfp,kfp) * xxs(3,k)

			    xeta = xeta + dmdxf1(k,2,ifp,jfp,kfp) * xxs(1,k)
                yeta = yeta + dmdxf1(k,2,ifp,jfp,kfp) * xxs(2,k)
			    zeta = zeta + dmdxf1(k,2,ifp,jfp,kfp) * xxs(3,k)

			    xbeta = xbeta + dmdxf1(k,3,ifp,jfp,kfp) * xxs(1,k)
			    ybeta = ybeta + dmdxf1(k,3,ifp,jfp,kfp) * xxs(2,k)
                zbeta = zbeta + dmdxf1(k,3,ifp,jfp,kfp) * xxs(3,k)
            end do

            S1(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
            S1(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
            S1(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta

            S1(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
            S1(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
            S1(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta

            S1(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
            S1(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
            S1(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi
        end do
        end do
        end do

        ! loop to find the S vector at the flux points (family 2)
        do kfp = 1,N
        do jfp = 1,N+1
        do ifp = 1,N
            xxi =   0.d0
            xeta =  0.d0
            xbeta = 0.d0
    
            yxi =   0.d0
            yeta =  0.d0
            ybeta = 0.d0
    
            zxi =   0.d0
            zeta =  0.d0
            zbeta = 0.d0
    
            do k = 1,kmax
                xxi = xxi + dmdxf2(k,1,ifp,jfp,kfp) * xxs(1,k)
                yxi = yxi + dmdxf2(k,1,ifp,jfp,kfp) * xxs(2,k)
                zxi = zxi + dmdxf2(k,1,ifp,jfp,kfp) * xxs(3,k)
    
                xeta = xeta + dmdxf2(k,2,ifp,jfp,kfp) * xxs(1,k)
                yeta = yeta + dmdxf2(k,2,ifp,jfp,kfp) * xxs(2,k)
                zeta = zeta + dmdxf2(k,2,ifp,jfp,kfp) * xxs(3,k)
    
                xbeta = xbeta + dmdxf2(k,3,ifp,jfp,kfp) * xxs(1,k)
                ybeta = ybeta + dmdxf2(k,3,ifp,jfp,kfp) * xxs(2,k)
                zbeta = zbeta + dmdxf2(k,3,ifp,jfp,kfp) * xxs(3,k)
            end do
    
            S2(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
            S2(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
            S2(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta
    
            S2(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
            S2(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
            S2(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta
    
            S2(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
            S2(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
            S2(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi
        end do
        end do
        end do

        ! loop to find the S vector at the flux points (family 3)
        do kfp = 1,N+1
        do jfp = 1,N
        do ifp = 1,N
            xxi =   0.d0
            xeta =  0.d0
            xbeta = 0.d0
    
            yxi =   0.d0
            yeta =  0.d0
            ybeta = 0.d0
    
            zxi =   0.d0
            zeta =  0.d0
            zbeta = 0.d0
    
            do k = 1,kmax
                xxi = xxi + dmdxf3(k,1,ifp,jfp,kfp) * xxs(1,k)
                yxi = yxi + dmdxf3(k,1,ifp,jfp,kfp) * xxs(2,k)
                zxi = zxi + dmdxf3(k,1,ifp,jfp,kfp) * xxs(3,k)
    
                xeta = xeta + dmdxf3(k,2,ifp,jfp,kfp) * xxs(1,k)
                yeta = yeta + dmdxf3(k,2,ifp,jfp,kfp) * xxs(2,k)
                zeta = zeta + dmdxf3(k,2,ifp,jfp,kfp) * xxs(3,k)
    
                xbeta = xbeta + dmdxf3(k,3,ifp,jfp,kfp) * xxs(1,k)
                ybeta = ybeta + dmdxf3(k,3,ifp,jfp,kfp) * xxs(2,k)
                zbeta = zbeta + dmdxf3(k,3,ifp,jfp,kfp) * xxs(3,k)
            end do

            S3(1,1,ifp,jfp,kfp,ic) = yeta*zbeta - ybeta*zeta
            S3(1,2,ifp,jfp,kfp,ic) = xbeta*zeta - xeta*zbeta
            S3(1,3,ifp,jfp,kfp,ic) = xeta*ybeta - xbeta*yeta

            S3(2,1,ifp,jfp,kfp,ic) = ybeta*zxi - yxi*zbeta
            S3(2,2,ifp,jfp,kfp,ic) = xxi*zbeta - xbeta*zxi
            S3(2,3,ifp,jfp,kfp,ic) = xbeta*yxi - xxi*ybeta

            S3(3,1,ifp,jfp,kfp,ic) = yxi*zeta - yeta*zxi
            S3(3,2,ifp,jfp,kfp,ic) = xeta*zxi- xxi*zeta
            S3(3,3,ifp,jfp,kfp,ic) = xxi*yeta - xeta*yxi
        end do
        end do
        end do

    end do


    END SUBROUTINE calcjacob