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
module mod_mkurban
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

  public :: mkurban_base , mkurban_param

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: levurbdim = 'nlevurb'
  character(len=16) , parameter :: solardim = 'numsolar'
  character(len=16) , parameter :: raddim = 'numrad'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter :: varname = 'PCT_URBAN'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  integer , parameter :: nparam = 24
  character(len=16) , dimension(nparam) , parameter :: parmname = &
    (/'ALB_IMPROAD     ', 'ALB_PERROAD     ', 'ALB_ROOF        ', &
      'ALB_WALL        ', 'CV_IMPROAD      ', 'CV_ROOF         ', &
      'CV_WALL         ', 'TK_IMPROAD      ', 'TK_ROOF         ', &
      'TK_WALL         ', 'CANYON_HWR      ', 'EM_IMPROAD      ', &
      'EM_PERROAD      ', 'EM_ROOF         ', 'EM_WALL         ', &
      'HT_ROOF         ', 'NLEV_IMPROAD    ', 'THICK_ROOF      ', &
      'THICK_WALL      ', 'T_BUILDING_MAX  ', 'T_BUILDING_MIN  ', &
      'WIND_HGT_CANYON ', 'WTLUNIT_ROOF    ', 'WTROAD_PERV     '/)

  integer(ik4) , public :: npu2d = 14
  integer(ik4) , public :: npu3d = 6
  integer(ik4) , public :: npu4d = 4

  integer(ik4) , dimension(nparam) , parameter :: parmdim = &
    (/4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
  real(rk8) :: vmin = 0.0D0
  real(rk8) :: vmisdat = -9999.0D0

  contains

  subroutine mkurban_base(urbanfile,urban)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rk8) , dimension(:,:) , intent(out) :: urban
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarid , ivarmask , ncid
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) :: istatus , i , li , lo
    real(rk8) , dimension(:,:) , allocatable :: rvar , rmask
    real(rk8) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
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

    istatus = nf90_inq_varid(ncid,latvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
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

    allocate(rvar(sum(domain%ni),domain%nj))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable urban in file '//trim(inpfile))
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
      istatus = nf90_get_var(ncid,ivarid,rvar(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable urban from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable urban from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar,rmask,rlon,rlat,urban,xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkurban_base

  subroutine mkurban_param(urbanfile,urban2d,urban3d,urban4d)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rk8) , dimension(:,:,:) , intent(out) :: urban2d
    real(rk8) , dimension(:,:,:,:) , intent(out) :: urban3d
    real(rk8) , dimension(:,:,:,:,:) , intent(out) :: urban4d
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarid , ivarmask , ncid
    integer(ik4) , dimension(4) :: istart , icount
    integer(ik4) :: istatus , i , n , li , lo
    integer(ik4) :: nurb , nrad , nsol
    real(rk8) , dimension(:,:,:) , allocatable :: rvar3d
    real(rk8) , dimension(:,:,:,:) , allocatable :: rvar4d
    integer(ik4) :: i2 , i3 , i4
    real(rk8) , dimension(:,:) , allocatable :: rvar , rmask
    real(rk8) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
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

    istatus = nf90_inq_varid(ncid,latvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lon in file '//trim(inpfile))

    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0D0 )
      glon = glon - 360.0D0
    end where
    where ( glon < -180.0D0 )
      glon = glon + 360.0D0
    end where

    istatus = nf90_inq_dimid(ncid,levurbdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension nlevurb in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nurb)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension nlevurb in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,solardim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension numsolar in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nsol)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension numsolar in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,raddim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension numrad in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nrad)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension numrad in file '//trim(inpfile))

    call get_window(glat,glon,domain)

    allocate(rvar(sum(domain%ni),domain%nj))
    allocate(rvar3d(sum(domain%ni),domain%nj,nurb))
    allocate(rvar4d(sum(domain%ni),domain%nj,nrad,nsol))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

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
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:), &
                             istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    i2 = lbound(urban2d,3)
    i3 = lbound(urban3d,4)
    i4 = lbound(urban4d,5)
    do n = 1 , nparam
      istatus = nf90_inq_varid(ncid,parmname(n),ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find variable '//trim(parmname(n))// &
        ' in file '//trim(inpfile))
      select case (parmdim(n))
        case (2)
          li = 1
          do i = 1 , domain%ntiles
            istart(1) = domain%igstart(i)
            icount(1) = domain%ni(i)
            istart(2) = domain%jgstart
            icount(2) = domain%nj
            lo = li+domain%ni(i)-1
            istatus = nf90_get_var(ncid,ivarid,rvar(li:lo,:), &
                                   istart(1:2),icount(1:2))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Cannot read variable '//trim(parmname(n))// &
              ' from file '//trim(inpfile))
            li = li + domain%ni(i)
          end do
          call bilinear(rvar,rmask,rlon,rlat,urban2d(:,:,i2), &
                   xlon,xlat,vmin,vmisdat)
          i2 = i2 + 1
        case (3)
          li = 1
          do i = 1 , domain%ntiles
            istart(1) = domain%igstart(i)
            icount(1) = domain%ni(i)
            istart(2) = domain%jgstart
            icount(2) = domain%nj
            lo = li+domain%ni(i)-1
            istart(3) = 1
            icount(3) = nurb
            istatus = nf90_get_var(ncid,ivarid,rvar3d(li:lo,:,:), &
                                   istart(1:3),icount(1:3))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Cannot read variable '//trim(parmname(n))// &
              ' from file '//trim(inpfile))
            li = li + domain%ni(i)
          end do
          call bilinear(rvar3d,rmask,rlon,rlat,urban3d(:,:,1:nurb,i3), &
                   xlon,xlat,vmin,vmisdat)
          i3 = i3 + 1
        case (4)
          li = 1
          do i = 1 , domain%ntiles
            istart(1) = domain%igstart(i)
            icount(1) = domain%ni(i)
            istart(2) = domain%jgstart
            icount(2) = domain%nj
            lo = li+domain%ni(i)-1
            istart(3) = 1
            icount(3) = nrad
            istart(4) = 1
            icount(4) = nsol
            istatus = nf90_get_var(ncid,ivarid,rvar4d(li:lo,:,:,:), &
                                   istart(1:4),icount(1:4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Cannot read variable '//trim(parmname(n))// &
              ' from file '//trim(inpfile))
            li = li + domain%ni(i)
          end do
          call bilinear(rvar4d,rmask,rlon,rlat,urban4d(:,:,1:nsol,1:nrad,i4), &
                   xlon,xlat,vmin,vmisdat)
          i4 = i4 +1
        case default
          call die(__FILE__,'Variable dimension not implemented',__LINE__)
      end select
    end do

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkurban_param

end module mod_mkurban
