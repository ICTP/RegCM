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

module mod_grid

  use mod_realkinds
  use mod_dynparam
  use mod_memutil

  private

  real(rk4), public :: clatx, clonx

  real(rk4), public, pointer, contiguous, dimension(:,:) :: xlat, xlon, xmask
  real(rk4), public, pointer, contiguous, dimension(:) :: xlat1d
  real(rk4), public, pointer, contiguous, dimension(:) :: xlon1d
  real(rk4), public, pointer, contiguous, dimension(:) :: sigx
  real(rkx), public, pointer, contiguous, dimension(:) :: zita
  real(rk4), public, pointer, contiguous, dimension(:) :: ax, bx

  public :: init_domain

  contains

  subroutine init_domain
    implicit none
    call getmem2d(xlat,1,jx,1,iy,'mod_read_domain:xlat')
    call getmem2d(xlon,1,jx,1,iy,'mod_read_domain:xlon')
    call getmem2d(xmask,1,jx,1,iy,'mod_read_domain:xmask')
    call getmem1d(xlon1d,1,jx,'mod_read_domain:xlon1d')
    call getmem1d(xlat1d,1,iy,'mod_read_domain:xlat1d')
    call getmem1d(sigx,1,kzp1,'mod_read_domain:sigx')
    call getmem1d(zita,1,kzp1,'mod_read_domain:zita')
    call getmem1d(ax,1,kzp1,'mod_read_domain:ax')
    call getmem1d(bx,1,kzp1,'mod_read_domain:bx')
  end subroutine init_domain

end module mod_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
