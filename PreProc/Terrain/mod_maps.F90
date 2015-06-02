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

  real(rk8) , pointer , dimension(:,:) :: coriol , dlat , dlon ,   &
                   dmap , htgrid , lndout , mask , dpth , snowam , &
                   smoist , texout , xlat , xlon , xmap , ps0
  real(rk8) , pointer , dimension(:,:,:) :: frac_tex , rmoist
  real(rk8) , pointer , dimension(:,:,:) :: pr0 , t0 , rho0

  real(rk8) , pointer , dimension(:,:) :: coriol_s , dlat_s , &
                      dlon_s , dmap_s , htgrid_s , lndout_s , &
                      mask_s , dpth_s , snowam_s , smoist_s , &
                      texout_s , xlat_s , xlon_s , xmap_s , ps0_s
  real(rk8) , pointer , dimension(:,:,:) :: frac_tex_s , rmoist_s
  real(rk8) , pointer , dimension(:,:,:) :: pr0_s , t0_s , rho0_s

  real(rk8) , pointer , dimension(:) :: sigma

  contains

  subroutine prepare_grid(jx,iy,kz,ntex,nsoil,idyn)
    implicit none
    integer(ik4) , intent(in) :: jx , iy , kz , ntex , nsoil , idyn
    call getmem1d(sigma,1,kz+1,'maps:sigma')
    call getmem2d(coriol,1,jx,1,iy,'maps:coriol')
    call getmem2d(xlat,1,jx,1,iy,'maps:xlat')
    call getmem2d(xlon,1,jx,1,iy,'maps:xlon')
    call getmem2d(dlat,1,jx,1,iy,'maps:dlat')
    call getmem2d(dlon,1,jx,1,iy,'maps:dlon')
    call getmem2d(dmap,1,jx,1,iy,'maps:dmap')
    call getmem2d(xmap,1,jx,1,iy,'maps:xmap')
    call getmem2d(htgrid,1,jx,1,iy,'maps:htgrid')
    call getmem2d(dpth,1,jx,1,iy,'maps:dpth')
    call getmem2d(lndout,1,jx,1,iy,'maps:lndout')
    call getmem2d(mask,1,jx,1,iy,'maps:mask')
    call getmem2d(snowam,1,jx,1,iy,'maps:snowam')
    call getmem2d(smoist,1,jx,1,iy,'maps:smoist')
    call getmem2d(texout,1,jx,1,iy,'maps:texout')
    call getmem3d(frac_tex,1,jx,1,iy,1,ntex,'maps:frac_tex')
    call getmem3d(rmoist,1,jx,1,iy,1,nsoil,'maps:rmoist')
    if ( idyn == 2 ) then
      call getmem2d(ps0,1,jx,1,iy,'maps:ps0')
      call getmem3d(pr0,1,jx,1,iy,1,kz+1,'maps:pr0')
      call getmem3d(t0,1,jx,1,iy,1,kz+1,'maps:t0')
      call getmem3d(rho0,1,jx,1,iy,1,kz+1,'maps:rho0')
    end if
  end subroutine prepare_grid

  subroutine prepare_subgrid(jxsg,iysg,kz,ntex,nsoil,idyn)
    implicit none
    integer(ik4) , intent(in) :: jxsg , iysg , kz , ntex , nsoil , idyn
    call getmem2d(coriol_s,1,jxsg,1,iysg,'maps:coriol_s')
    call getmem2d(xlat_s,1,jxsg,1,iysg,'maps:xlat_s')
    call getmem2d(xlon_s,1,jxsg,1,iysg,'maps:xlon_s')
    call getmem2d(dlat_s,1,jxsg,1,iysg,'maps:dlat_s')
    call getmem2d(dlon_s,1,jxsg,1,iysg,'maps:dlon_s')
    call getmem2d(dmap_s,1,jxsg,1,iysg,'maps:dmap_s')
    call getmem2d(xmap_s,1,jxsg,1,iysg,'maps:xmap_s')
    call getmem2d(htgrid_s,1,jxsg,1,iysg,'maps:htgrid_s')
    call getmem2d(dpth_s,1,jxsg,1,iysg,'maps:dpth_s')
    call getmem2d(lndout_s,1,jxsg,1,iysg,'maps:lndout_s')
    call getmem2d(mask_s,1,jxsg,1,iysg,'maps:mask_s')
    call getmem2d(snowam_s,1,jxsg,1,iysg,'maps:snowam_s')
    call getmem2d(smoist_s,1,jxsg,1,iysg,'maps:smoist_s')
    call getmem2d(texout_s,1,jxsg,1,iysg,'maps:texout_s')
    call getmem3d(frac_tex_s,1,jxsg,1,iysg,1,ntex,'maps:frac_tex_s')
    call getmem3d(rmoist_s,1,jxsg,1,iysg,1,nsoil,'maps:rmoist_s')
    if ( idyn == 2 ) then
      call getmem2d(ps0_s,1,jxsg,1,iysg,'maps:ps0_s')
      call getmem3d(pr0_s,1,jxsg,1,iysg,1,kz+1,'maps:pr0_s')
      call getmem3d(t0_s,1,jxsg,1,iysg,1,kz+1,'maps:t0_s')
      call getmem3d(rho0_s,1,jxsg,1,iysg,1,kz+1,'maps:rho0_s')
    end if
  end subroutine prepare_subgrid

end module mod_maps
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
