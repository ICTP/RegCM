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

  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_stdio
  use mod_message
  use mod_nchelper
  use mod_domain

  private

  real(rk8) , public , pointer , dimension(:,:) :: xlat , xlon , dlat , dlon
  real(rk8) , public , pointer , dimension(:,:) :: topogm , mask , landuse
  real(rk8) , public , pointer , dimension(:,:) :: pa , tlayer , za
  real(rk8) , public , pointer , dimension(:,:) :: b3pd
  real(rk8) , public , pointer , dimension(:) :: sigma2
  real(rk8) , public , pointer , dimension(:) :: sigmaf
  real(rk8) , public :: delx
  integer(ik4) , public :: i0 , i1 , j0 , j1
  real(rk8) , public :: lat0 , lat1 , lon0 , lon1

  public :: init_grid

  contains

  subroutine init_grid(nx,ny,nz)
    implicit none
    integer(ik4) , intent(in) :: nx , ny , nz
    call getmem2d(xlat,1,nx,1,ny,'mod_grid:xlat')
    call getmem2d(xlon,1,nx,1,ny,'mod_grid:xlon')
    call getmem2d(dlat,1,nx,1,ny,'mod_grid:dlat')
    call getmem2d(dlon,1,nx,1,ny,'mod_grid:dlon')
    call getmem2d(topogm,1,nx,1,ny,'mod_grid:topogm')
    call getmem2d(mask,1,nx,1,ny,'mod_grid:mask')
    call getmem2d(landuse,1,nx,1,ny,'mod_grid:landuse')
    call getmem2d(pa,1,nx,1,ny,'mod_grid:pa')
    call getmem2d(tlayer,1,nx,1,ny,'mod_grid:tlayer')
    call getmem2d(za,1,nx,1,ny,'mod_grid:za')
    call getmem2d(b3pd,1,nx,1,ny,'mod_grid:b3pd')
    call getmem1d(sigma2,1,nz,'mod_grid:sigma2')
    call getmem1d(sigmaf,1,nz+1,'mod_grid:sigmaf')
    call read_domain_info
  end subroutine init_grid

  subroutine read_domain_info
    use mod_dynparam
    implicit none
    integer(ik4) :: incin
    character(len=256) :: fname
    integer(ik4) :: k
    fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    call openfile_withname(fname,incin)
    call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,topogm,mask,landuse)
    call closefile(incin)
    do k = 1 , kz
      sigma2(k) = d_half*(sigmaf(k+1)+sigmaf(k))
    end do
  end subroutine read_domain_info
!
end module mod_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
