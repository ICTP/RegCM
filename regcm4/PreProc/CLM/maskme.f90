      subroutine maskme(landmask,vals,vmisdat,nlon,nlat,nlev,ntim)
      implicit none
!
! Dummy arguments
!
      integer :: nlat , nlev , nlon , ntim
      real :: vmisdat
      real , dimension(nlon,nlat) :: landmask
      real , dimension(nlon,nlat,nlev,ntim) :: vals
      intent (in) landmask , nlat , nlev , nlon , ntim , vmisdat
      intent (inout) vals
!
! Local variables
!
      integer :: i , j , k , l
!
!     ** Variables Passed in
!     ** Local variables
 
      do l = 1 , ntim
        do k = 1 , nlev
          do j = 1 , nlat
            do i = 1 , nlon
              if ( landmask(i,j)>0.5 ) then
                vals(i,j,k,l) = vals(i,j,k,l)
              else
                vals(i,j,k,l) = vmisdat
              end if
            end do
          end do
        end do
      end do
 
      end subroutine maskme
