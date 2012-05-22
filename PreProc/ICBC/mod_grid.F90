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

  use netcdf
  use mod_memutil
  use mod_realkinds
  use mod_stdio
  use mod_message
  use mod_nchelper
  use mod_domain

  private

  real(sp) , public , pointer , dimension(:,:) :: xlat , xlon , dlat , dlon
  real(sp) , public , pointer , dimension(:,:) :: topogm , mask , landuse
  real(sp) , public , pointer , dimension(:,:) :: pa , tlayer , za
  real(sp) , public , pointer , dimension(:,:) :: b3pd
  real(sp) , public , pointer , dimension(:) :: dsigma , sigma2
  real(sp) , public , pointer , dimension(:) :: sigmaf
  real(dp) , public :: delx
  integer , public :: i0 , i1 , j0 , j1
  real(dp) , public :: lat0 , lat1 , lon0 , lon1

  public :: init_grid

  contains

  subroutine init_grid(nx,ny,nz)
    implicit none
    integer , intent(in) :: nx , ny , nz
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
    call getmem1d(dsigma,1,nz,'mod_grid:dsigma')
    call getmem1d(sigma2,1,nz,'mod_grid:sigma2')
    call getmem1d(sigmaf,1,nz+1,'mod_grid:sigmaf')
    call read_domain_info
  end subroutine init_grid

  subroutine read_domain_info
    use mod_dynparam
    implicit none
    integer :: istatus
    integer :: incin
    character(256) :: fname
    integer :: k
    fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    call openfile_withname(fname,incin)
    istatus = nf90_get_att(incin,nf90_global,'grid_factor',xcone)
    call checkncerr(istatus,__FILE__,__LINE__,'Attribute grid_factor missing')
    call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,topogm,mask,landuse)
    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error closing file '//trim(fname))
    do k = 1 , kz
      sigma2(k) = 0.5*(sigmaf(k+1)+sigmaf(k))
      dsigma(k) = sigmaf(k+1) - sigmaf(k)
    end do
  end subroutine read_domain_info
!
end module mod_grid
