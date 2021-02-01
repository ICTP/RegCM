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

module mod_pgw

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_vertint
  use mod_earth
  use mod_hgt
  use mod_humid
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use mod_kdinterp
  use netcdf

  private

  public :: init_pgw , get_pgw , conclude_pgw

  integer :: ncid = -1
  integer :: ilon , jlat , klev

  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , dimension(:) :: gtemp
  real(rkx) , pointer , dimension(:) :: plevs
  real(rkx) , pointer , dimension(:) :: sigmar

  real(rkx) , pointer , dimension(:,:,:) :: u2 , v2 , q2 , z2 , t2

  real(rkx) , pointer , dimension(:,:,:) :: u3 , v3 , q3 , z3 , t3
  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v
  real(rkx) , pointer , dimension(:,:) :: ps3
  real(rkx) , pointer , dimension(:,:) :: topou , topov

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  real(rkx) :: pss , pst

  character(len=256) :: pgwfile

  contains

  subroutine init_pgw(filename)
    implicit none
    character(len=*) , intent(in) :: filename
    integer :: istatus , idimid , ivarid
    integer :: k

    pgwfile = filename
    istatus = nf90_open(pgwfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lat dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat dimelen in file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lon dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon dimelen in file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'plev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing plev dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading plev dimelen in file '//trim(pgwfile))

    call getmem1d(plevs,1,klev,'mod_era5:plevs')
    call getmem1d(sigmar,1,klev,'mod_era5:sigmar')
    call getmem1d(glat,1,jlat,'mod_era5:glat')
    call getmem1d(glon,1,ilon,'mod_era5:glon')
    call getmem1d(gtemp,1,max(ilon,jlat),'mod_era5:gtemp')

    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lat variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat variable in file '//trim(pgwfile))
    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lon variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon variable in file '//trim(pgwfile))
    istatus = nf90_inq_varid(ncid,'plev',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing plev variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading plev variable in file '//trim(pgwfile))

    do k = 1 , klev
      sigmar(k) = (plevs(klev-k+1)-plevs(1))/(plevs(klev)-plevs(1))
    end do
    pss = (plevs(klev)-plevs(1))/1000.0_rkx ! Pa -> cb
    pst = plevs(1)/1000.0_rkx ! Pa -> cb

    call get_window(glat,glon,xlat,xlon,i_band,gdomain)

    gtemp(1:jlat) = glat
    jlat = gdomain%nj
    call getmem1d(glat,1,jlat,'mod_era5:glat')
    glat = gtemp(gdomain%jgstart:gdomain%jgstop)
    gtemp(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_era5:glon')
    glon(1:gdomain%ni(1)) = gtemp(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = gtemp(gdomain%igstart(2):gdomain%igstop(2))
    end if
    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    call getmem3d(u2,1,ilon,1,jlat,1,klev,'u2')
    call getmem3d(v2,1,ilon,1,jlat,1,klev,'v2')
    call getmem3d(t2,1,ilon,1,jlat,1,klev,'t2')
    call getmem3d(q2,1,ilon,1,jlat,1,klev,'q2')
    call getmem3d(z2,1,ilon,1,jlat,1,klev,'z2')

    call getmem3d(u3,1,jx,1,iy,1,klev,'u3')
    call getmem3d(v3,1,jx,1,iy,1,klev,'v3')
    call getmem3d(t3,1,jx,1,iy,1,klev,'t3')
    call getmem3d(q3,1,jx,1,iy,1,klev,'q3')
    call getmem3d(z3,1,jx,1,iy,1,klev,'z3')
    call getmem2d(ps3,1,jx,1,iy,'ps3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev,'mod_era5:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev,'mod_era5:d3v')
      call getmem2d(topou,1,jx,1,iy,'mod_era5:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_era5:topov')
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    end if
  end subroutine init_pgw

  subroutine get_pgw(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    

  end subroutine get_pgw

  subroutine conclude_pgw
    implicit none
    integer :: istatus
    if ( ncid >= 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file '//pgwfile)
    end if
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_pgw

end module mod_pgw
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
