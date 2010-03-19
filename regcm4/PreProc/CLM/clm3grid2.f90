      subroutine clm3grid2(nlon,nlat,nfld,glon,glat,istart,icount,ifld, &
                         & zlon,zlat,zlev)
 
      implicit none
!
! Dummy arguments
!
      integer :: ifld , nfld , nlat , nlon
      real , dimension(nlat) :: glat
      real , dimension(nlon) :: glon
      integer , dimension(4) :: icount , istart
      real , dimension(icount(2)) :: zlat
      real , dimension(icount(3)) :: zlev
      real , dimension(icount(1)) :: zlon
      intent (in) glat , glon , icount , istart , nlat , nlon
      intent (out) zlat , zlev , zlon
!
! Local variables
!
      integer :: i , j , k
! 
      do i = 1 , icount(1)
        zlon(i) = glon(i+istart(1)-1)
      end do
      do j = 1 , icount(2)
        zlat(j) = glat(j+istart(2)-1)
      end do
      do k = 1 , icount(3)
        zlev(k) = icount(3) - k + 1
      end do
 
      end subroutine clm3grid2
