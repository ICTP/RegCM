      subroutine humid1fv(t,q,p3d,ni,nj,nk)
      use mod_constants , only : tzero , ep2 , lh0 , lh1 , lsvp1 , lsvp2
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: tr = 1./tzero
!
! Dummy arguments
!
      integer :: ni , nj , nk
      real , dimension(ni,nj,nk) :: p3d , q , t
      intent (in) ni , nj , nk , p3d , t
      intent (inout) q
!
! Local variables
!
      real :: hl , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
      do i = 1 , ni
        do j = 1 , nj
          do k = 1 , nk
            if ( p3d(i,j,k)>-9990. ) then
              hl = lh0 - lh1*(t(i,j,k)-tzero)  ! LATENT HEAT OF EVAP.
              satvp = lsvp1*exp(lsvp2*hl*(tr-1./t(i,j,k)))
                                                   ! SATURATION VAP PRESS.
              qs = ep2*satvp/(p3d(i,j,k)-satvp)   ! SAT. MIXING RATIO
              q(i,j,k) = amax1(q(i,j,k)/qs,0.0)    !ALREADY MIXING RATIO
            else
              q(i,j,k) = -9999.
            end if
          end do
        end do
      end do
      end subroutine humid1fv
