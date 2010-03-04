      subroutine top2btm(x,nlon1,nlat1,nlev1)
      implicit none
!
! Dummy arguments
!
      integer :: nlat1 , nlev1 , nlon1
      real , dimension(nlon1,nlat1,nlev1) :: x
      intent (in) nlat1 , nlev1 , nlon1
      intent (inout) x
!
! Local variables
!
      integer :: i , j , k , kr
      real , dimension(30) :: work
!
      do i = 1 , nlon1
        do j = 1 , nlat1
          do k = 1 , nlev1
            work(k) = x(i,j,k)
          end do
          do k = 1 , nlev1
            kr = nlev1 - k + 1
            x(i,j,k) = work(kr)
          end do
        end do
      end do
      end subroutine top2btm
