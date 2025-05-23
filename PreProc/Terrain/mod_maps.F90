!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_maps

  use mod_intkinds
  use mod_realkinds
  use mod_memutil

  real(rkx), pointer, contiguous, dimension(:,:) :: coriol, dlat, dlon,   &
                   dmap, htgrid, lndout, mask, dpth, snowam, &
                   smoist, texout, xlat, xlon, xmap, ps0,    &
                   ulat, ulon, vlat, vlon, umap, vmap
  real(rkx), pointer, contiguous, dimension(:,:,:) :: frac_tex, rmoist, rts
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pr0, t0, rho0, z0
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zeta, fmz

  real(rkx), pointer, contiguous, dimension(:,:) :: coriol_s, dlat_s, &
                      dlon_s, dmap_s, htgrid_s, lndout_s, &
                      mask_s, dpth_s, snowam_s, smoist_s, &
                      texout_s, xlat_s, xlon_s, xmap_s,   &
                      ulat_s, ulon_s, vlat_s, vlon_s,     &
                      umap_s, vmap_s, ps0_s
  real(rkx), pointer, contiguous, dimension(:,:,:) :: frac_tex_s, rmoist_s, rts_s
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pr0_s, t0_s, rho0_s, z0_s
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zeta_s, fmz_s

  real(rkx), pointer, contiguous, dimension(:) :: sigma
  real(rkx), pointer, contiguous, dimension(:) :: zita
  real(rkx), pointer, contiguous, dimension(:) :: ak
  real(rkx), pointer, contiguous, dimension(:) :: bk

  contains

  subroutine prepare_grid(jx,iy,kz,ntex,nsoil,idyn)
    implicit none
    integer(ik4), intent(in) :: jx, iy, kz, ntex, nsoil, idyn
    call getmem1d(sigma,1,kz+1,'maps:sigma')
    call getmem2d(coriol,1,jx,1,iy,'maps:coriol')
    call getmem2d(xlat,1,jx,1,iy,'maps:xlat')
    call getmem2d(xlon,1,jx,1,iy,'maps:xlon')
    call getmem2d(dlat,1,jx,1,iy,'maps:dlat')
    call getmem2d(dlon,1,jx,1,iy,'maps:dlon')
    call getmem2d(htgrid,1,jx,1,iy,'maps:htgrid')
    call getmem2d(dpth,1,jx,1,iy,'maps:dpth')
    call getmem2d(lndout,1,jx,1,iy,'maps:lndout')
    call getmem2d(mask,1,jx,1,iy,'maps:mask')
    call getmem2d(snowam,1,jx,1,iy,'maps:snowam')
    call getmem2d(smoist,1,jx,1,iy,'maps:smoist')
    call getmem2d(texout,1,jx,1,iy,'maps:texout')
    call getmem3d(frac_tex,1,jx,1,iy,1,ntex,'maps:frac_tex')
    call getmem3d(rmoist,1,jx,1,iy,1,nsoil,'maps:rmoist')
    call getmem3d(rts,1,jx,1,iy,1,nsoil,'maps:rts')
    if ( idyn == 2 ) then
      call getmem2d(ps0,1,jx,1,iy,'maps:ps0')
      call getmem3d(pr0,1,jx,1,iy,1,kz+1,'maps:pr0')
      call getmem3d(t0,1,jx,1,iy,1,kz+1,'maps:t0')
      call getmem3d(rho0,1,jx,1,iy,1,kz+1,'maps:rho0')
      call getmem3d(z0,1,jx,1,iy,1,kz+1,'maps:z0')
    end if
    if ( idyn == 3 ) then
      call getmem1d(zita,1,kz+1,'maps:zita')
      call getmem1d(ak,1,kz+1,'maps:ak')
      call getmem1d(bk,1,kz+1,'maps:bk')
      call getmem3d(zeta,1,jx,1,iy,1,kz+1,'maps:zeta')
      call getmem3d(fmz,1,jx,1,iy,1,kz+1,'maps:fmz')
      call getmem2d(ulat,1,jx,1,iy,'maps:ulat')
      call getmem2d(ulon,1,jx,1,iy,'maps:ulon')
      call getmem2d(vlat,1,jx,1,iy,'maps:vlat')
      call getmem2d(vlon,1,jx,1,iy,'maps:vlon')
      call getmem2d(xmap,1,jx,1,iy,'maps:xmap')
      call getmem2d(umap,1,jx,1,iy,'maps:umap')
      call getmem2d(vmap,1,jx,1,iy,'maps:vmap')
    else
      call getmem2d(dmap,1,jx,1,iy,'maps:dmap')
      call getmem2d(xmap,1,jx,1,iy,'maps:xmap')
    end if
  end subroutine prepare_grid

  subroutine prepare_subgrid(jxsg,iysg,kz,ntex,nsoil,idyn)
    implicit none
    integer(ik4), intent(in) :: jxsg, iysg, kz, ntex, nsoil, idyn
    call getmem2d(coriol_s,1,jxsg,1,iysg,'maps:coriol_s')
    call getmem2d(xlat_s,1,jxsg,1,iysg,'maps:xlat_s')
    call getmem2d(xlon_s,1,jxsg,1,iysg,'maps:xlon_s')
    call getmem2d(dlat_s,1,jxsg,1,iysg,'maps:dlat_s')
    call getmem2d(dlon_s,1,jxsg,1,iysg,'maps:dlon_s')
    call getmem2d(htgrid_s,1,jxsg,1,iysg,'maps:htgrid_s')
    call getmem2d(dpth_s,1,jxsg,1,iysg,'maps:dpth_s')
    call getmem2d(lndout_s,1,jxsg,1,iysg,'maps:lndout_s')
    call getmem2d(mask_s,1,jxsg,1,iysg,'maps:mask_s')
    call getmem2d(snowam_s,1,jxsg,1,iysg,'maps:snowam_s')
    call getmem2d(smoist_s,1,jxsg,1,iysg,'maps:smoist_s')
    call getmem2d(texout_s,1,jxsg,1,iysg,'maps:texout_s')
    call getmem3d(frac_tex_s,1,jxsg,1,iysg,1,ntex,'maps:frac_tex_s')
    call getmem3d(rmoist_s,1,jxsg,1,iysg,1,nsoil,'maps:rmoist_s')
    call getmem3d(rts_s,1,jxsg,1,iysg,1,nsoil,'maps:rts_s')
    if ( idyn == 2 ) then
      call getmem2d(ps0_s,1,jxsg,1,iysg,'maps:ps0_s')
      call getmem3d(pr0_s,1,jxsg,1,iysg,1,kz+1,'maps:pr0_s')
      call getmem3d(t0_s,1,jxsg,1,iysg,1,kz+1,'maps:t0_s')
      call getmem3d(rho0_s,1,jxsg,1,iysg,1,kz+1,'maps:rho0_s')
      call getmem3d(z0_s,1,jxsg,1,iysg,1,kz+1,'maps:z0_s')
    end if
    if ( idyn == 3 ) then
      call getmem3d(zeta_s,1,jxsg,1,iysg,1,kz+1,'maps:zeta_s')
      call getmem3d(fmz_s,1,jxsg,1,iysg,1,kz+1,'maps:fmz_s')
      call getmem2d(ulat_s,1,jxsg,1,iysg,'maps:ulat_s')
      call getmem2d(ulon_s,1,jxsg,1,iysg,'maps:ulon_s')
      call getmem2d(vlat_s,1,jxsg,1,iysg,'maps:vlat_s')
      call getmem2d(vlon_s,1,jxsg,1,iysg,'maps:vlon_s')
      call getmem2d(xmap_s,1,jxsg,1,iysg,'maps:xmap_s')
      call getmem2d(umap_s,1,jxsg,1,iysg,'maps:umap_s')
      call getmem2d(vmap_s,1,jxsg,1,iysg,'maps:vmap_s')
    else
      call getmem2d(dmap_s,1,jxsg,1,iysg,'maps:dmap_s')
      call getmem2d(xmap_s,1,jxsg,1,iysg,'maps:xmap_s')
    end if
  end subroutine prepare_subgrid

end module mod_maps
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
