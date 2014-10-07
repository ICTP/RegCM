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
module mod_mkvocef
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

  public :: mkvocef

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: latvar = 'lat'
  character(len=16) , parameter :: lonvar = 'lon'
  character(len=16) , parameter :: varname1 = 'ef_btr'
  character(len=16) , parameter :: varname2 = 'ef_crp'
  character(len=16) , parameter :: varname3 = 'ef_fdt'
  character(len=16) , parameter :: varname4 = 'ef_fet'
  character(len=16) , parameter :: varname5 = 'ef_grs'
  character(len=16) , parameter :: varname6 = 'ef_shr'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rk8) :: vmin = 0.0D0
  real(rk8) :: vmisdat = -9999.0D0

  contains

  subroutine mkvocef(vocfile,vocef)
    implicit none
    character(len=*) , intent(in) :: vocfile
    real(rk8) , dimension(:,:,:) , intent(out) :: vocef
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarmask , ncid
    integer(ik4) , dimension(6) :: ivarid
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) :: istatus , i , j , li , lo
    real(rk8) , dimension(:,:,:) , allocatable :: rvar
    real(rk8) , dimension(:,:) , allocatable :: rmask
    real(rk8) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//vocfile
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

    istatus = nf90_inq_varid(ncid,latvar,ivarid(1))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid(1),glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid(1))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid(1),glon)
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

    allocate(rvar(sum(domain%ni),domain%nj,6))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname1,ivarid(1))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 1 in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname2,ivarid(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 2 in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname3,ivarid(3))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 3 in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname4,ivarid(4))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 4 in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname5,ivarid(5))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 5 in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname6,ivarid(6))
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable vocef 6 in file '//trim(inpfile))
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
      do j = 1 , 6
        istatus = nf90_get_var(ncid,ivarid(j),rvar(li:lo,:,j),istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot read variable vocef from file '//trim(inpfile))
      end do
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar,rmask,rlon,rlat,vocef,xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkvocef

end module mod_mkvocef
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
