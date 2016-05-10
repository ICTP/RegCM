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
  use mod_nhinterp

  private

  real(rkx) , public , pointer , dimension(:,:) :: xlat , xlon , dlat , dlon
  real(rkx) , public , pointer , dimension(:,:) :: topogm , mask , landuse
  real(rkx) , public , pointer , dimension(:,:) :: msfx , msfd
  real(rkx) , public , pointer , dimension(:,:) :: pa , tlayer , za
  real(rkx) , public , pointer , dimension(:,:) :: pd4
  real(rkx) , public , pointer , dimension(:) :: sigmah
  real(rkx) , public , pointer , dimension(:) :: sigmaf
  real(rkx) , public , pointer , dimension(:) :: dsigma
  real(rkx) , public , pointer , dimension(:,:,:) :: pr0, t0, rho0
  real(rkx) , public , pointer , dimension(:,:) :: ps0

  real(rkx) , public :: delx
  integer(ik4) , public :: i0 , i1 , j0 , j1
  real(rkx) , public :: lat0 , lat1 , lon0 , lon1

  real(rkx) , public :: logp_lrate = 50.0_rkx

  public :: init_grid

  contains

  subroutine init_grid(nx,ny,nz)
    use mod_dynparam
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
    call getmem2d(pd4,1,nx,1,ny,'mod_grid:pd4')
    call getmem1d(sigmah,1,nz,'mod_grid:sigmah')
    call getmem1d(sigmaf,1,nz+1,'mod_grid:sigmaf')
    call getmem1d(dsigma,1,nz,'mod_write:dsigma')
    if ( idynamic == 2 ) then
      call getmem2d(msfx,1,nx,1,ny,'mod_write:msfx')
      call getmem2d(msfd,1,nx,1,ny,'mod_write:msfd')
      call getmem2d(ps0,1,nx,1,ny,'mod_write:ps0')
      call getmem3d(pr0,1,nx,1,ny,1,nz,'mod_write:pr0')
      call getmem3d(rho0,1,nx,1,ny,1,nz,'mod_write:rho0')
      call getmem3d(t0,1,nx,1,ny,1,nz,'mod_write:t0')
    end if
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
    if ( idynamic == 2 ) then
      call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,topogm, &
                       mask,landuse,msfx,msfd)
      call nhsetup(ptop,stdp,stdt,logp_lrate)
    else
      call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,topogm, &
                       mask,landuse)
    end if
    call closefile(incin)
    do k = 1 , kz
      sigmah(k) = d_half*(sigmaf(k+1)+sigmaf(k))
      dsigma(k) = (sigmaf(k+1)-sigmaf(k))
    end do
    if ( idynamic == 2 ) then
      call nhbase(1,iy,1,jx,kz,sigmah,topogm,ps0,pr0,t0,rho0)
    end if
  end subroutine read_domain_info

end module mod_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
