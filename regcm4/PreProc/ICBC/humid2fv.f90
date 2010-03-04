      subroutine humid2fv(t,q,ps,pt,sigma,ni,nj,nk)
      use mod_constants , only : tzero , lh0 , lh1 , lsvp1 , lsvp2 , ep2
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: tr = 1./tzero
!
! Dummy arguments
!
      integer :: ni , nj , nk
      real :: pt
      real , dimension(ni,nj) :: ps
      real , dimension(ni,nj,nk) :: q , t
      real , dimension(nk) :: sigma
      intent (in) ni , nj , nk , ps , pt , sigma , t
      intent (inout) q
!
! Local variables
!
      real :: hl , p , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
!
      do i = 1 , ni
        do j = 1 , nj
          do k = 1 , nk
            p = (pt+sigma(k)*ps(i,j))*10.
            hl = lh0 - lh1*(t(i,j,k)-tzero)
            satvp = lsvp1*exp(lsvp2*hl*(tr-1./t(i,j,k)))
            qs = ep2*satvp/(p-satvp)
            q(i,j,k) = amax1(q(i,j,k)*qs,0.0)
          end do
        end do
      end do
!
      end subroutine humid2fv
