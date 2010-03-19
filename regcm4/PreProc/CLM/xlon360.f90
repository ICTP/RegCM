      subroutine xlon360(nx,ny,xlon)
      implicit none
!
! Dummy arguments
!
      integer :: nx , ny
      real , dimension(ny,nx) :: xlon
      intent (in) nx , ny
      intent (inout) xlon
!
! Local variables
!
      integer :: i , j
      real , dimension(ny,nx) :: xlon1
!
      do i = 1 , nx
        do j = 1 , ny
          xlon1(j,i) = xlon(j,i)
        end do
      end do
 
      do i = 1 , nx
        do j = 1 , ny
          if ( xlon1(j,i)<0. ) then
            xlon(j,i) = xlon1(j,i) + 360.
          else
            xlon(j,i) = xlon1(j,i)
          end if
        end do
      end do
 
      end subroutine xlon360
