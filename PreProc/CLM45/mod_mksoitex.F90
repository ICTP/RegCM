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
module mod_mksoitex
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

  public :: mksoitex

  character(len=16) , parameter :: levdim = 'number_of_layers'
  character(len=32) , parameter :: mpudim = 'max_value_mapunit'
  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter :: varname1 = 'MAPUNITS'
  character(len=16) , parameter :: varname2 = 'PCT_SAND'
  character(len=16) , parameter :: varname3 = 'PCT_CLAY'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real :: vmin = 0.0
  real :: vmisdat = -9999.0

  contains

  subroutine mksoitex(soitexfile,sand,clay)
    implicit none
    character(len=*) , intent(in) :: soitexfile
    real(rk4) , dimension(:,:,:) , intent(out) :: sand , clay
    integer(ik4) :: nlat , nlon , nlev , nmpu
    integer(ik4) :: idimid , ivarid1 , ivarid2 , ivarid3 , ivarmask , ncid
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) :: istatus , i , j , l , li , lo , number_of_mapunits
    real(rk4) , dimension(:,:,:) , allocatable :: rvar1 , rvar2
    real(rk4) , dimension(:,:) , allocatable :: temp1 , temp2
    real(rk4) , dimension(:,:) , allocatable :: rmask , mpu
    real(rk4) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//soitexfile
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,levdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension level in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension level in file '//trim(inpfile))

    if ( nlev > size(sand,3) .or. &
         nlev > size(clay,3) ) then
      call die(__FILE__,'Size too small for clay/sand in mksoitex',__LINE__)
    end if

    istatus = nf90_inq_dimid(ncid,mpudim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension mpu in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nmpu)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension mpu in file '//trim(inpfile))

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

    istatus = nf90_inq_dimid(ncid,'number_of_mapunits',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension number_of_mapunits in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=number_of_mapunits)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension number_of_mapunits in file '//trim(inpfile))

    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0 )
      glon = glon - 360.0
    end where
    where ( glon < -180.0 )
      glon = glon + 360.0
    end where

    call get_window(glat,glon,domain)

    allocate(rvar1(sum(domain%ni),domain%nj,nlev))
    allocate(rvar2(sum(domain%ni),domain%nj,nlev))
    allocate(temp1(nmpu,nlev))
    allocate(temp2(nmpu,nlev))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(mpu(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname1,ivarid1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable soitex in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname2,ivarid2)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable soitex in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,varname3,ivarid3)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable soitex in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    istatus = nf90_get_var(ncid,ivarid2,temp1)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable sand from file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid3,temp2)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable clay from file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarid1,mpu(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable soitex from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable landmask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    do j = 1 , domain%nj
      do i = 1 , sum(domain%ni)
        if ( mpu(i,j) > 0 .and. mpu(i,j) <= number_of_mapunits ) then
          rvar1(i,j,:) = temp1(int(mpu(i,j)),:)
          rvar2(i,j,:) = temp2(int(mpu(i,j)),:)
          if ( any(rvar1(i,j,:) > 100.0) ) then
            write(stderr,*) 'CLAY > 100.0'
          end if
          if ( any(rvar2(i,j,:) > 100.0) ) then
            write(stderr,*) 'SAND > 100.0'
          end if
        else
          rvar1(i,j,:) = vmisdat
          rvar2(i,j,:) = vmisdat
        end if
      end do
    end do

    call bilinear(rvar1,rmask,rlon,rlat,sand,xlon,xlat,vmin,vmisdat)
    call bilinear(rvar2,rmask,rlon,rlat,clay,xlon,xlat,vmin,vmisdat)

    do l = 1  , size(clay,3)
      do i = 1 , size(clay,2)
        do j = 1 , size(clay,1)
          if ( clay(j,i,l) + sand(j,i,l) /= 100.0 ) then
            if ( clay(j,i,l) > sand(j,i,l) ) then
              sand(j,i,l) = 100.0 - clay(j,i,l)
            else
              clay(j,i,l) = 100.0 - sand(j,i,l)
            end if
          end if
        end do
      end do
    end do

    deallocate(glat,glon,rlat,rlon,rvar1,rvar2,temp1,temp2,rmask)
  end subroutine mksoitex

end module mod_mksoitex
