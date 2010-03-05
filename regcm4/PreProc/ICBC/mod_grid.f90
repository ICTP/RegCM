      module mod_grid

      use mod_param , only : jx , iy , kz

      implicit none

      real , dimension(jx,iy) :: coriol , dlat , dlon , msfx , snowcv , &
                               & topogm , toposdgm , xlandu , xlat ,    &
                               & xlon , zs , zssd
      real , dimension(kz) :: dsigma , sigma2
      real , dimension(kz+1) :: sigmaf
      real :: delx , grdfac
      character(6) :: lgtype
      integer :: i0 , i1 , j0 , j1
      real :: lat0 , lat1 , lon0 , lon1

      end module mod_grid
