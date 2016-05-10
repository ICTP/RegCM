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
module mod_mklaisai
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

  public :: mklaisai

  character(len=16) , parameter :: timedim = 'time'
  character(len=16) , parameter :: pftdim = 'pft'
  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter :: varname1 = 'MONTHLY_LAI'
  character(len=16) , parameter :: varname2 = 'MONTHLY_SAI'
  character(len=32) , parameter :: varname3 = 'MONTHLY_HEIGHT_BOT'
  character(len=32) , parameter :: varname4 = 'MONTHLY_HEIGHT_TOP'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rkx) :: vmin = 0.0_rkx
  real(rkx) :: vmisdat = -9999.0_rkx

  contains

  subroutine mklaisai(laisaifile, &
                  monthly_lai,monthly_sai,monthly_top,monthly_bot)
    implicit none
    character(len=*) , intent(in) :: laisaifile
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_sai , monthly_lai
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_top , monthly_bot
    integer(ik4) :: nlat , nlon , npft , nmon
    integer(ik4) :: idimid , ivarid1 , ivarid2 , ivarid3 , ivarid4
    integer(ik4) :: ivarmask , ncid
    integer(ik4) , dimension(4) :: istart , icount
    integer(ik4) :: istatus , i , li , lo
    real(rkx) , dimension(:,:,:,:) , allocatable :: rvar1 , rvar2
    real(rkx) , dimension(:,:,:,:) , allocatable :: rvar3 , rvar4
    real(rkx) , dimension(:,:) , allocatable :: rmask
    real(rkx) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//laisaifile
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,timedim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension time in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nmon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension time in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,pftdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension pft in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=npft)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension pft in file '//trim(inpfile))

    if ( nmon > size(monthly_sai,4) .or. &
         npft > size(monthly_sai,3) .or. &
         any(shape(monthly_sai) /= shape(monthly_lai)) ) then
      call die(__FILE__, &
              'Size too small/invalid for lai/sai in mklaisai',__LINE__)
    end if

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
    where ( glon >  180.0_rkx )
      glon = glon - 360.0_rkx
    end where
    where ( glon < -180.0_rkx )
      glon = glon + 360.0_rkx
    end where

    call get_window(glat,glon,domain)

    allocate(rvar1(sum(domain%ni),domain%nj,npft,nmon))
    allocate(rvar2(sum(domain%ni),domain%nj,npft,nmon))
    allocate(rvar3(sum(domain%ni),domain%nj,npft,nmon))
    allocate(rvar4(sum(domain%ni),domain%nj,npft,nmon))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname1,ivarid1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable MONTHLY_LAI in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname2,ivarid2)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable MONTHLY_SAI in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname3,ivarid3)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable MONTHLY_HEIGHT_TOP in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname4,ivarid4)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable MONTHLY_HEIGHT_BOT in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      istart(3) = 1
      icount(3) = npft
      istart(4) = 1
      icount(4) = nmon
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarid1,rvar1(li:lo,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable laisai from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarid2,rvar2(li:lo,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable lake from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarid3,rvar3(li:lo,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable lake from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarid4,rvar4(li:lo,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable lake from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:), &
              istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable landmask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar1,rmask,rlon,rlat,monthly_lai(:,:,1:npft,:),&
            xlon,xlat,vmin,vmisdat)
    call bilinear(rvar2,rmask,rlon,rlat,monthly_sai(:,:,1:npft,:),&
            xlon,xlat,vmin,vmisdat)
    call bilinear(rvar3,rmask,rlon,rlat,monthly_bot(:,:,1:npft,:),&
            xlon,xlat,vmin,vmisdat)
    call bilinear(rvar4,rmask,rlon,rlat,monthly_top(:,:,1:npft,:),&
            xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar1,rvar2,rvar3,rvar4,rmask)
  end subroutine mklaisai

end module mod_mklaisai
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
