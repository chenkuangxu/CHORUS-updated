	SUBROUTINE BCvisflux

	use setup3d
	IMPLICIT NONE
	INCLUDE 'mpif.h'
       
	integer :: faml
	integer :: nfp,ifpl,jfpl,ic,iface,sp,ifn,kfpl,mfp
	double precision, dimension(numv) :: Qfl,Qfr,Qfl2,Qfr2,Qfl3,Qfr3

    DO ifn=1,NINLET
        iface = IBFINL(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp=1,N  
            faml = Gfp2Lfp(3,mfp,nfp,iface) 

            ifpl = Gfp2Lfp(5,mfp,nfp,iface) 
            jfpl = Gfp2Lfp(6,mfp,nfp,iface) 
	        kfpl = Gfp2Lfp(7,mfp,nfp,iface) 

            Qfl(1:numv) = 0.d0
            Qfl2(1:numv) = 0.d0
	        Qfl3(1:numv) = 0.d0

            if(faml==1) then
                do sp=1,N
                    Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                end do   

            else if(faml==2) then
                do sp=1,N
                    Qfl(1:numv)  = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                end do
            
            else if (faml==3) then
                do sp=1,N
                    Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)

                end do
            end if !faml loop

            Qfr(1:numv)  = Qfl(1:numv)
            Qfr2(1:numv) = Qfl2(1:numv)
            Qfr3(1:numv) = Qfl3(1:numv)

            if(faml==1) then
                nablaQvfi(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfi(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfi(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            else if(faml==2) then
                nablaQvfj(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfj(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfj(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            else
                nablaQvfk(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfk(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfk(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            end if
        end do
        end do

    END DO


    DO ifn=1,NOUTLET
        iface = IBFOUT(ifn)
        ic  = IF2C(iface,1)
        do nfp=1,N
        do mfp=1,N  
            faml = Gfp2Lfp(3,mfp,nfp,iface) 

            ifpl = Gfp2Lfp(5,mfp,nfp,iface) 
            jfpl = Gfp2Lfp(6,mfp,nfp,iface) 
	        kfpl = Gfp2Lfp(7,mfp,nfp,iface) 

            Qfl(1:numv) = 0.d0
            Qfl2(1:numv) = 0.d0
	        Qfl3(1:numv) = 0.d0

            if(faml==1) then
                do sp=1,N
                    Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,sp,jfpl,kfpl,ic)*Lmat(ifpl,sp)   
                end do

             else if(faml==2) then
                do sp=1,N
                    Qfl(1:numv)  = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,sp,kfpl,ic)*Lmat(jfpl,sp)
                end do
             
            else if (faml==3) then
                do sp=1,N
                    Qfl(1:numv) = Qfl(1:numv) + nablaQs(1:numv,1,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl2(1:numv) = Qfl2(1:numv) + nablaQs(1:numv,2,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)
                    Qfl3(1:numv) = Qfl3(1:numv) + nablaQs(1:numv,3,ifpl,jfpl,sp,ic)*Lmat(kfpl,sp)

                end do
            end if !faml loop
!
            Qfr(1:numv)  = Qfl(1:numv)
            Qfr2(1:numv) = Qfl2(1:numv)
            Qfr3(1:numv) = Qfl3(1:numv)

            if(faml==1) then
                nablaQvfi(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfi(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfi(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            else if(faml==2) then
                nablaQvfj(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfj(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfj(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            else
                nablaQvfk(1:numv,1,ifpl,jfpl,kfpl,ic) = 0.5*Qfl(1:numv) +0.5*Qfr(1:numv)
                nablaQvfk(1:numv,2,ifpl,jfpl,kfpl,ic) = 0.5*Qfl2(1:numv)+0.5*Qfr2(1:numv)
                nablaQvfk(1:numv,3,ifpl,jfpl,kfpl,ic) = 0.5*Qfl3(1:numv)+0.5*Qfr3(1:numv)
            end if
        end do
        end do

    END DO

    end SUBROUTINE BCvisflux


