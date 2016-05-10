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
module mod_mkdynpft
#ifdef CN
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

  public :: mkdynpft

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: pftdim = 'pft'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter :: varname = 'PCT_PFT'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rkx) :: vmin = 0.0_rkx
  real(rkx) :: vmisdat = -9999.0_rkx

  contains

  subroutine mkdynpft(dynpft,year)
    implicit none
    real(rkx) , dimension(:,:,:) , intent(out) :: dynpft
    integer(ik4) , intent(in) :: year
    integer(ik4) :: nlat , nlon , npft
    integer(ik4) :: idimid , ivarid , ivarmask , ncid
    integer(ik4) , dimension(3) :: istart , icount
    integer(ik4) :: istatus , i , li , lo
    real(rkx) , dimension(:,:,:) , allocatable :: rvar
    real(rkx) , dimension(:,:) , allocatable :: rmask
    real(rkx) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain
    character(len=32) :: p1 , p2
    character(len=4) :: cy

    character(len=256) :: inpfile

!    if ( year > 2100 ) year = 2100
    p1 = 'dynamic'
    p2 = '.'
    write(cy,'(i0.4)') year

    if ( year > 2005 ) then
      select case (dattyp(4:5))
        case ('RF')
          continue
        case ('26')
          p2 = 'SCENARIO'//pthsep//'RCP2.6'
        case ('45')
          p2 = 'SCENARIO'//pthsep//'RCP4.5'
        case ('60')
          p2 = 'SCENARIO'//pthsep//'RCP6.0'
        case ('85', '15')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          if ( dattyp /= "EIN15" .and. &
               dattyp(1:4) /= "NNRP" .and. &
               dattyp /= "JRA55" ) then
            call die(__FILE__, &
              'Dynamic landuse only supported for CMIP5',__LINE__)
          end if
      end select
    end if

    inpfile = trim(inpglob)//pthsep//'CLM45'//pthsep//'surface'// &
            pthsep//trim(p1)//pthsep//trim(p2)//pthsep//&
            'mksrf_landuse_'//cy//'.nc'
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,pftdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension pft in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=npft)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension pft in file '//trim(inpfile))

    if ( npft > size(dynpft,3) ) then
      call die(__FILE__,'Size too small for dynpft in mkdynpft',__LINE__)
    end if

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
    where ( glon >  180.0_rkx )
      glon = glon - 360.0_rkx
    end where
    where ( glon < -180.0_rkx )
      glon = glon + 360.0_rkx
    end where

    call get_window(glat,glon,domain)

    allocate(rvar(sum(domain%ni),domain%nj,npft))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable dynpft in file '//trim(inpfile))
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
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarid,rvar(li:lo,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable dynpft from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask, &
              rmask(li:lo,:),istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable landmask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    dynpft = 0.0_rkx
    call bilinear(rvar,rmask,rlon,rlat,dynpft(:,:,1:npft), &
            xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkdynpft
#endif

end module mod_mkdynpft
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
