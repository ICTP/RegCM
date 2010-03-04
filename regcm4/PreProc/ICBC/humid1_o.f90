      subroutine humid1_o(t,q,ps,sigma,ptop,im,jm,km)
      use mod_constants , only : tzero , ep2 , lh0 , lh1 , lsvp1 , lsvp2
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: tr = 1./tzero
!
! Dummy arguments
!
      integer :: im , jm , km
      real :: ptop
      real , dimension(im,jm) :: ps
      real , dimension(im,jm,km) :: q , t
      real , dimension(km) :: sigma
      intent (in) im , jm , km , ps , ptop , sigma , t
      intent (inout) q
!
! Local variables
!
      real :: hl , p , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!     DATA ON SIGMA LEVELS
!
      do k = 1 , km
        do j = 1 , jm
          do i = 1 , im
            p = sigma(k)*(ps(i,j)-ptop) + ptop
            hl = lh0 - lh1*(t(i,j,k)-tzero)       ! LATENT HEAT OF EVAP.
            satvp = lsvp1*exp(lsvp2*hl*(tr-1./t(i,j,k)))
                                                      ! SATURATION VAP PRESS.
            qs = ep2*satvp/(p-satvp)                 ! SAT. MIXING RATIO
            q(i,j,k) = amax1(q(i,j,k)/qs,0.0)
          end do
        end do
      end do
      end subroutine humid1_o
