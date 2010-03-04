!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
