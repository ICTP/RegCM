      subroutine p1p2(pd,px,ni,nj)
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj
      real , dimension(ni,nj) :: pd , px
      intent (in) ni , nj , px
      intent (out) pd
!
! Local variables
!
      integer :: i , j , ni1 , nj1
!
!     THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
!     ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
!     SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE
!     I'S.
!
      ni1 = ni - 1
      nj1 = nj - 1
!
      do j = 2 , nj1
        do i = 2 , ni1
          pd(i,j) = 0.25*(px(i,j)+px(i-1,j)+px(i,j-1)+px(i-1,j-1))
        end do
      end do
!
      do i = 2 , ni1
        pd(i,1) = 0.5*(px(i,1)+px(i-1,1))
        pd(i,nj) = 0.5*(px(i,nj1)+px(i-1,nj1))
      end do
!
      do j = 2 , nj1
        pd(1,j) = 0.5*(px(1,j)+px(1,j-1))
        pd(ni,j) = 0.5*(px(ni1,j)+px(ni1,j-1))
      end do
!
      pd(1,1) = px(1,1)
      pd(1,nj) = px(1,nj1)
      pd(ni,1) = px(ni1,1)
      pd(ni,nj) = px(ni1,nj1)
!
      end subroutine p1p2
