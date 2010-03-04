      subroutine humid1(t,q,ps,pt,sigma,ni,nj,nk)
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
      real :: ps , pt
      real , dimension(ni,nj,nk) :: q , t
      real , dimension(nk) :: sigma
      intent (in) ni , nj , nk , ps , pt , sigma , t
      intent (inout) q
!
! Local variables
!
      real :: lh , p , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
      do i = 1 , ni
        do j = 1 , nj
          do k = 1 , nk
            p = (pt+sigma(k)*ps)*10.        ! PRESSURE AT LEVEL K
            lh = lh0 - lh1*(t(i,j,k)-tzero)
            satvp = lsvp1*dexp(lsvp2*lh*(1./tzero-1./t(i,j,k)))
            qs = ep2*satvp/(p-satvp)        ! SAT. MIXING RATIO
            q(i,j,k) = amax1(q(i,j,k)/qs,0.0)
          end do
        end do
      end do
      end subroutine humid1
