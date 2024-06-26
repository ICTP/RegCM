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

module mod_grid

  use mod_realkinds
  use mod_dynparam
  use mod_memutil

  private

  real(rk4) , public :: clatx , clonx

  real(rk4) , public , pointer , dimension(:,:) :: xlat, xlon , xmask
  real(rk4) , public , pointer , dimension(:) :: xlat1d
  real(rk4) , public , pointer , dimension(:) :: xlon1d
  real(rk4) , public , pointer , dimension(:) :: sigx
  real(rkx) , public , pointer , dimension(:) :: zita
  real(rk4) , public , pointer , dimension(:) :: ax , bx

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
