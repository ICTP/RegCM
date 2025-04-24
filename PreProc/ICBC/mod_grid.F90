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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_stdio
  use mod_message
  use mod_nchelper
  use mod_domain
  use mod_nhinterp
  use mod_zita
  use mod_dynparam, only : idynamic, base_state_pressure, logp_lrate
  use mod_dynparam, only : mo_ztop
  use mod_projections

  private

  real(rkx), public, pointer, contiguous, dimension(:,:) :: xlat, xlon, dlat, dlon
  real(rkx), public, pointer, contiguous, dimension(:,:) :: ulat, ulon, vlat, vlon
  real(rkx), public, pointer, contiguous, dimension(:,:) :: topogm, mask, landuse
  real(rkx), public, pointer, contiguous, dimension(:,:) :: msfx, msfd
  real(rkx), public, pointer, contiguous, dimension(:,:) :: pa, tlayer, za
  real(rkx), public, pointer, contiguous, dimension(:) :: sigmah
  real(rkx), public, pointer, contiguous, dimension(:) :: sigmaf
  real(rkx), public, pointer, contiguous, dimension(:) :: dsigma
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: pr0, t0, rho0, z0
  real(rkx), public, pointer, contiguous, dimension(:,:) :: ps0

  real(rkx), public :: delx
  integer(ik4), public :: i0, i1, j0, j1
  real(rkx), public :: lat0, lat1, lon0, lon1, ts0

  type(regcm_projection), public :: pju, pjv, pjd

  public :: init_grid, init_hgrid

  contains

  subroutine init_hgrid(nx,ny,nz)
    use mod_dynparam
    implicit none
    integer(ik4), intent(in) :: nx, ny, nz
    call getmem2d(xlat,1,nx,1,ny,'mod_grid:xlat')
    call getmem2d(xlon,1,nx,1,ny,'mod_grid:xlon')
    call getmem2d(topogm,1,nx,1,ny,'mod_grid:topogm')
    call getmem2d(mask,1,nx,1,ny,'mod_grid:mask')
    call getmem1d(sigmaf,1,nz+1,'mod_grid:sigmaf')
    call read_domain_hinfo
  end subroutine init_hgrid

  subroutine init_grid(nx,ny,nz)
    use mod_dynparam
    implicit none
    integer(ik4), intent(in) :: nx, ny, nz
    call getmem2d(xlat,1,nx,1,ny,'mod_grid:xlat')
    call getmem2d(xlon,1,nx,1,ny,'mod_grid:xlon')
    if ( idynamic == 3 ) then
      call getmem2d(ulat,1,nx,1,ny,'mod_grid:ulat')
      call getmem2d(ulon,1,nx,1,ny,'mod_grid:ulon')
      call getmem2d(vlat,1,nx,1,ny,'mod_grid:vlat')
      call getmem2d(vlon,1,nx,1,ny,'mod_grid:vlon')
    else
      call getmem2d(dlat,1,nx,1,ny,'mod_grid:dlat')
      call getmem2d(dlon,1,nx,1,ny,'mod_grid:dlon')
    end if
    call getmem2d(topogm,1,nx,1,ny,'mod_grid:topogm')
    call getmem2d(mask,1,nx,1,ny,'mod_grid:mask')
    call getmem2d(landuse,1,nx,1,ny,'mod_grid:landuse')
    call getmem2d(pa,1,nx,1,ny,'mod_grid:pa')
    call getmem2d(tlayer,1,nx,1,ny,'mod_grid:tlayer')
    call getmem2d(za,1,nx,1,ny,'mod_grid:za')
    call getmem1d(sigmah,1,nz,'mod_grid:sigmah')
    call getmem1d(sigmaf,1,nz+1,'mod_grid:sigmaf')
    call getmem1d(dsigma,1,nz,'mod_write:dsigma')
    if ( idynamic == 2 ) then
      call getmem2d(msfx,1,nx,1,ny,'mod_write:msfx')
      call getmem2d(msfd,1,nx,1,ny,'mod_write:msfd')
      call getmem2d(ps0,1,nx,1,ny,'mod_write:ps0')
      call getmem3d(pr0,1,nx,1,ny,1,nz,'mod_write:pr0')
      call getmem3d(rho0,1,nx,1,ny,1,nz,'mod_write:rho0')
      call getmem3d(z0,1,nx,1,ny,1,nz,'mod_write:z0')
      call getmem3d(t0,1,nx,1,ny,1,nz,'mod_write:t0')
    else if ( idynamic == 3 ) then
      call getmem3d(z0,1,nx,1,ny,1,nz,'mod_write:z0')
    end if
    call read_domain_info
  end subroutine init_grid

  subroutine read_domain_hinfo
    use mod_dynparam
    implicit none
    integer(ik4) :: incin
    character(len=256) :: fname
    fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    call openfile_withname(fname,incin)
    call read_domain(incin,sigmaf,xlat,xlon,ht=topogm,mask=mask)
    call closefile(incin)
  end subroutine read_domain_hinfo

  subroutine read_domain_info
    use mod_dynparam
    implicit none
    integer(ik4) :: incin
    character(len=256) :: fname
    integer(ik4) :: i, j, k
    real(rkx), dimension(kz) :: zitah
    fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    call openfile_withname(fname,incin)
    if ( idynamic == 2 ) then
      call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,ht=topogm, &
                       mask=mask,lndcat=landuse,msfx=msfx,msfd=msfd)
      call read_reference_surface_temp(incin,ts0)
      call nhsetup(ptop,base_state_pressure,logp_lrate,ts0)
      do k = 1, kz
        sigmah(k) = d_half*(sigmaf(k+1)+sigmaf(k))
        dsigma(k) = (sigmaf(k+1)-sigmaf(k))
      end do
    else if ( idynamic == 3 ) then
      call read_domain(incin,sigmaf,xlat,xlon,ulat=ulat,ulon=ulon, &
                       vlat=vlat,vlon=vlon,ht=topogm,mask=mask,    &
                       lndcat=landuse)
      call model_zitah(zitah,mo_ztop)
      sigmah = sigmazita(zitah,mo_ztop)
      do k = 1, kz
        do i = 1, iy
          do j = 1, jx
            z0(j,i,k) = md_zeta_h(zitah(k),topogm(j,i),mo_ztop,mo_h,mo_a0)
          end do
        end do
        dsigma(k) = (sigmaf(k+1)-sigmaf(k))
      end do
    else
      call read_domain(incin,sigmaf,xlat,xlon,dlat,dlon,ht=topogm, &
                       mask=mask,lndcat=landuse)
      do k = 1, kz
        sigmah(k) = d_half*(sigmaf(k+1)+sigmaf(k))
        dsigma(k) = (sigmaf(k+1)-sigmaf(k))
      end do
    end if
    call closefile(incin)
    if ( idynamic == 2 ) then
      call nhbase(1,iy,1,jx,kz,sigmah,topogm,ps0,pr0,t0,rho0,z0)
    end if
  end subroutine read_domain_info

end module mod_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
