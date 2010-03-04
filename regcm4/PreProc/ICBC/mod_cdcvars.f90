      module mod_cdcvars

      use mod_domain, only : jx , iy , kz

      implicit none

      real , dimension(jx,iy) :: b3pd
      real , dimension(jx,iy,kz) :: h4 , q4 , t4 , u4 , v4
      real , dimension(jx,iy) :: ps4 , ts4
      real , dimension(jx,iy) :: coriol , dlat , dlon , msfx , snowcv , &
                               & topogm , toposdgm , xlandu , xlat ,    &
                               & xlon
      real , dimension(kz) :: dsigma , sigma2
      real , dimension(kz+1) :: sigmaf

      end module mod_cdcvars
