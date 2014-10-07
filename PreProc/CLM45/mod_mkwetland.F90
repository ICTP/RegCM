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
module mod_mkwetland
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_getwindow
  use mod_bilinear
  use mod_nchelper
  use mod_memutil
  use netcdf

  implicit none

  private

  public :: mkwetland

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: latvar = 'lat'
  character(len=16) , parameter :: lonvar = 'lon'
  character(len=16) , parameter :: varname1 = 'PCT_WETLAND'
  character(len=16) , parameter :: varname2 = 'PCT_LAKE'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rk8) :: vmin = 0.0D0
  real(rk8) :: vmisdat = -9999.0D0

  contains

  subroutine mkwetland(wetfile,wetland,lake)
    implicit none
    character(len=*) , intent(in) :: wetfile
    real(rk8) , dimension(:,:) , intent(out) :: wetland , lake
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarid1 , ivarid2 , ivarmask , ncid
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) :: istatus , i , li , lo
    real(rk8) , dimension(:,:) , allocatable :: rvar1 , rvar2 , rmask
    real(rk8) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//wetfile
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,latdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lat in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lat in file '//trim(inpfile))

    allocate(glat(nlat))

    istatus = nf90_inq_varid(ncid,latvar,ivarid1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid1,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid1,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lon in file '//trim(inpfile))

    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0D0 )
      glon = glon - 360.0D0
    end where
    where ( glon < -180.0D0 )
      glon = glon + 360.0D0
    end where

    call get_window(glat,glon,domain)

    allocate(rvar1(sum(domain%ni),domain%nj))
    allocate(rvar2(sum(domain%ni),domain%nj))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname1,ivarid1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable wetland in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname2,ivarid2)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable wetland in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarid1,rvar1(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable wetland from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarid2,rvar2(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable lake from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable landmask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar1,rmask,rlon,rlat,wetland,xlon,xlat,vmin,vmisdat)
    call bilinear(rvar2,rmask,rlon,rlat,lake,xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar1,rvar2,rmask)
  end subroutine mkwetland

end module mod_mkwetland
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
