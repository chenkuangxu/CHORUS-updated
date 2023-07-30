        SUBROUTINE compresid

        use setup3d

        IMPLICIT NONE
        include 'mpif.h'   
        integer :: is,js,ks,ic,rfp
        double precision, dimension(5) :: dFpsi,dGeta,dHbeta
        double precision, dimension(5) :: dFvpsi,dGveta,dHvbeta
        double precision :: radi,gx,gy,gz

!-------------------------------  add gravity term first st ---------------------------    
!  residual  === negative rhs
        do ic = 1,NCELL
           do ks=1,N
           do js=1,N
           do is=1,N
            radi = sqrt(XXsolu(1,is,js,ks,ic)**2+ XXsolu(2,is,js,ks,ic)**2+ &
                        XXsolu(3,is,js,ks,ic)**2)

!  decompose
      gx = -XXsolu(1,is,js,ks,ic)/radi*gR(is,js,ks,ic)
      gy = -XXsolu(2,is,js,ks,ic)/radi*gR(is,js,ks,ic)
      gz = -XXsolu(3,is,js,ks,ic)/radi*gR(is,js,ks,ic)

            resid(1,is,js,ks,ic) = 0.d0
            resid(2,is,js,ks,ic) = Q(1,is,js,ks,ic)*gx 
            resid(3,is,js,ks,ic) = Q(1,is,js,ks,ic)*gy 
            resid(4,is,js,ks,ic) = Q(1,is,js,ks,ic)*gz 

            resid(5,is,js,ks,ic) = Q(2,is,js,ks,ic)*gx + &
                                   Q(3,is,js,ks,ic)*gy + &
                                   Q(4,is,js,ks,ic)*gz 

           end do
           end do
           end do
        end do
!add in Coriolis force
        do ic = 1,NCELL
           do ks=1,N
           do js=1,N
           do is=1,N
            resid(2,is,js,ks,ic) = resid(2,is,js,ks,ic) + &
                                 Q(3,is,js,ks,ic)*2.0*omega_rigid  
           
            resid(3,is,js,ks,ic) = resid(3,is,js,ks,ic) - &
                                 Q(2,is,js,ks,ic)*2.0*omega_rigid  
!   zero for z

!add in Centrifugal force
!            resid(2,is,js,ks,ic) = resid(2,is,js,ks,ic) + &
!                      Q(1,is,js,ks,ic)*omega_rigid**2*XXsolu(1,is,js,ks,ic) 
          
!            resid(3,is,js,ks,ic) = resid(3,is,js,ks,ic) + &
!                      Q(1,is,js,ks,ic)*omega_rigid**2*XXsolu(2,is,js,ks,ic) 
!   zero for z
           end do
           end do
           end do
        end do

        do ic = 1,NCELL
          do ks = 1,N
          do js = 1,N
          do is = 1,N
            resid(1:5,is,js,ks,ic) = Jac(is,js,ks,ic)*resid(1:5,is,js,ks,ic)
          end do
          end do
          end do
        end do

! ------------------------------------------------------------------------------------------------

        call compflux

        
        do ic=1,NCELL

! need to find the flux derivatives at each solution pt

           do ks=1,N
           do js=1,N
           do is=1,N

              ! compute the dFdpsi derivative		
              dFpsi(1:5) = 0.d0
              dFvpsi(1:5) = 0.d0
              do rfp=1,N+1
                 dFpsi(1:5)  = dFpsi(1:5)  + F1(1:5,rfp,js,ks,ic) *Mmat(rfp,is)
                 dFvpsi(1:5) = dFvpsi(1:5) + Fv1(1:5,rfp,js,ks,ic)*Mmat(rfp,is)
              end do

              ! compute the dGdeta derivative
              dGeta(1:5) = 0.d0
              dGveta(1:5) = 0.d0
              do rfp=1,N+1
                 dGeta(1:5)  = dGeta(1:5)  + G2(1:5,is,rfp,ks,ic) *Mmat(rfp,js)
                 dGveta(1:5) = dGveta(1:5) + Gv2(1:5,is,rfp,ks,ic)*Mmat(rfp,js)
              end do

! compute the dHdbeta derivative
              dHbeta(1:5) = 0.d0
              dHvbeta(1:5) = 0.d0
              do rfp=1,N+1
                  dHbeta(1:5) = dHbeta(1:5) + H3(1:5,is,js,rfp,ic)*Mmat(rfp,ks)
                 dHvbeta(1:5) = dHvbeta(1:5)+ Hv3(1:5,is,js,rfp,ic)*Mmat(rfp,ks)
              end do

              resid(1:5,is,js,ks,ic) = resid(1:5,is,js,ks,ic)  &
                                     -( dFpsi(1:5) + dGeta(1:5) + dHbeta(1:5))


           if (vismode.eq.1) resid(1:5,is,js,ks,ic) = resid(1:5,is,js,ks,ic) &
             +(dFvpsi(1:5)+dGveta(1:5)+dHvbeta(1:5))

             resid(1:5,is,js,ks,ic) = resid(1:5,is,js,ks,ic)/Jac(is,js,ks,ic)

           end do
           end do
           end do

        end do

        END SUBROUTINE compresid
