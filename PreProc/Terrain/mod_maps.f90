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

  use mod_intkinds
  use mod_realkinds
  use mod_memutil

  real(rk4) , pointer , dimension(:,:) :: coriol , dlat , dlon ,     &
                   dmap , htgrid , lndout , mask , dpth , snowam ,  &
                   texout , xlat , xlon , xmap
  real(rk4) , pointer , dimension(:,:,:) :: frac_tex

  real(rk4) , pointer , dimension(:,:) :: coriol_s , dlat_s ,   &
                      dlon_s , dmap_s , htgrid_s , lndout_s ,  &
                      mask_s , dpth_s , snowam_s , texout_s ,  &
                      xlat_s , xlon_s , xmap_s
  real(rk4) , pointer , dimension(:,:,:) :: frac_tex_s

  real(rk8) , pointer , dimension(:) :: sigma

  contains

  subroutine prepare_grid(iy,jx,kz,ntex)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , kz , ntex
    call getmem1d(sigma,1,kz+1,'maps:sigma')
    call getmem2d(coriol,1,iy,1,jx,'maps:coriol')
    call getmem2d(xlat,1,iy,1,jx,'maps:xlat')
    call getmem2d(xlon,1,iy,1,jx,'maps:xlon')
    call getmem2d(dlat,1,iy,1,jx,'maps:dlat')
    call getmem2d(dlon,1,iy,1,jx,'maps:dlon')
    call getmem2d(dmap,1,iy,1,jx,'maps:dmap')
    call getmem2d(xmap,1,iy,1,jx,'maps:xmap')
    call getmem2d(htgrid,1,iy,1,jx,'maps:htgrid')
    call getmem2d(dpth,1,iy,1,jx,'maps:dpth')
    call getmem2d(lndout,1,iy,1,jx,'maps:lndout')
    call getmem2d(mask,1,iy,1,jx,'maps:mask')
    call getmem2d(snowam,1,iy,1,jx,'maps:snowam')
    call getmem2d(texout,1,iy,1,jx,'maps:texout')
    call getmem3d(frac_tex,1,iy,1,jx,1,ntex,'maps:frac_tex')
  end subroutine prepare_grid

  subroutine prepare_subgrid(iysg,jxsg,ntex)
    implicit none
    integer(ik4) , intent(in) :: iysg , jxsg , ntex
    call getmem2d(coriol_s,1,iysg,1,jxsg,'maps:coriol_s')
    call getmem2d(xlat_s,1,iysg,1,jxsg,'maps:xlat_s')
    call getmem2d(xlon_s,1,iysg,1,jxsg,'maps:xlon_s')
    call getmem2d(dlat_s,1,iysg,1,jxsg,'maps:dlat_s')
    call getmem2d(dlon_s,1,iysg,1,jxsg,'maps:dlon_s')
    call getmem2d(dmap_s,1,iysg,1,jxsg,'maps:dmap_s')
    call getmem2d(xmap_s,1,iysg,1,jxsg,'maps:xmap_s')
    call getmem2d(htgrid_s,1,iysg,1,jxsg,'maps:htgrid_s')
    call getmem2d(dpth_s,1,iysg,1,jxsg,'maps:dpth_s')
    call getmem2d(lndout_s,1,iysg,1,jxsg,'maps:lndout_s')
    call getmem2d(mask_s,1,iysg,1,jxsg,'maps:mask_s')
    call getmem2d(snowam_s,1,iysg,1,jxsg,'maps:snowam_s')
    call getmem2d(texout_s,1,iysg,1,jxsg,'maps:texout_s')
    call getmem3d(frac_tex_s,1,iysg,1,jxsg,1,ntex,'maps:frac_tex_s')
  end subroutine prepare_subgrid

end module mod_maps
