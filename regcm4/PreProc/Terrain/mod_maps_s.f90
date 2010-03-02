      module mod_maps_s
      use mod_domain , only : iy , jx , nsg , nveg , ntex
      implicit none
      real(4) , dimension(iy*nsg,jx*nsg) :: claya_s , clayb_s ,         &
           & coriol_s , dlat_s , dlon_s , dmap_s , htgrid_s ,           &
           & htsdgrid_s ,  lndout_s , mask_s , sanda_s , sandb_s ,      &
           & snowam_s , texout_s , xlat_s , xlon_s , xmap_s
      real(4) , dimension(iy*nsg,jx*nsg,nveg) :: frac_lnd_s
      real(4) , dimension(iy*nsg,jx*nsg,ntex) :: frac_tex_s
      integer , dimension(iy*nsg,jx*nsg) :: intext_s , lnduse_s
      end module mod_maps_s
