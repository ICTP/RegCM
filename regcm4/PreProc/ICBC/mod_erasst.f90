      module mod_erasst
      use mod_domain , only : iy , jx
      implicit none
      integer , parameter :: ilon = 240 , jlat = 121
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      integer :: nnnend , nstart
      real(8) :: xadd , xscale
      real(4) , dimension(ilon,jlat) :: sst
      real , dimension(iy,jx) :: lu , sstmm , xlat , xlon
      end module mod_erasst
