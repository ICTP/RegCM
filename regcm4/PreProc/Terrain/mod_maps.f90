      module mod_maps
      use mod_domain , only : iy , jx , nveg , ntex
      implicit none
      real(4) , dimension(iy,jx) :: claya , clayb , coriol , dlat ,     &
                                  & dlon , dmap , htgrid , htsdgrid ,   &
                                  & lndout , mask , sanda , sandb ,     &
                                  & snowam , texout , xlat , xlon , xmap
      real(4) , dimension(iy,jx,nveg) :: frac_lnd
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      integer , dimension(iy,jx) :: intext , lnduse
      end module mod_maps
