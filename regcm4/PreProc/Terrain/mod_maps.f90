!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_maps

      use mod_regcm_param , only : iy , jx , nveg , iysg , jxsg
      use mod_preproc_param , only : ntex

      implicit none

      real(4) , dimension(iy,jx) :: claya , clayb , coriol , dlat ,     &
                                  & dlon , dmap , htgrid , htsdgrid ,   &
                                  & lndout , mask , sanda , sandb ,     &
                                  & snowam , texout , xlat , xlon , xmap
      real(4) , dimension(iy,jx,nveg) :: frac_lnd
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      integer , dimension(iy,jx) :: intext , lnduse

      real(4) , dimension(iysg,jxsg) :: claya_s , clayb_s , coriol_s ,  &
           & dlat_s , dlon_s , dmap_s , htgrid_s , htsdgrid_s ,         &
           & lndout_s , mask_s , sanda_s , sandb_s , snowam_s ,         &
           & texout_s , xlat_s , xlon_s , xmap_s
      real(4) , dimension(iysg,jxsg,nveg) :: frac_lnd_s
      real(4) , dimension(iysg,jxsg,ntex) :: frac_tex_s
      integer , dimension(iysg,jxsg) :: intext_s , lnduse_s

      end module mod_maps
