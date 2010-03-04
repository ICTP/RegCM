      subroutine mxmn3d(var,cvar,jx,iy,np)
      implicit none
!
! Dummy arguments
!
      character(2) :: cvar
      integer :: iy , jx , np
      real , dimension(jx,iy,np) :: var
      intent (in) cvar , iy , jx , np , var
!
! Local variables
!
      integer :: i , j , k
      real :: smax , smin
!
      do k = 1 , np
        smax = -1.E8
        smin = 1.E8
        do j = 1 , iy
          do i = 1 , jx
            if ( smax<var(i,j,k) ) smax = var(i,j,k)
            if ( smin>var(i,j,k) ) smin = var(i,j,k)
          end do
        end do
        write (*,*) cvar , k , smax , smin
      end do
      end subroutine mxmn3d
