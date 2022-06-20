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

module mod_cmip6_normm

  use mod_intkinds
  use mod_realkinds
  use mod_message
  use mod_date
  use mod_stdio
  use mod_memutil
  use mod_cmip6_helper
  use netcdf

  implicit none

  private

  character(len=*) , parameter :: normm_version = 'v20191108'
  character(len=*) , parameter :: normm_version1 = 'v20210319'
  character(len=*) , parameter :: normm_version2 = 'v20200702'
  character(len=*) , parameter :: normm_version3 = 'v20200218'
  character(len=*) , parameter :: normm_version4 = 'v20210319'

  public :: read_3d_normm , read_2d_normm , read_fx_normm

  contains

    subroutine read_hcoord_normm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_nor:lon')
      call getmem1d(lat,1,nlat,'cmip6_nor:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_normm

    subroutine read_vcoord_normm(ncid,a,b,p0)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: a , b
      real(rkx) , intent(out) :: p0
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(a,1,nlev,'cmip6_nor:a')
      call getmem1d(b,1,nlev,'cmip6_nor:b')
      istatus = nf90_inq_varid(ncid,'a',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,a)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(ncid,'p0',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(ncid,ivarid,p0)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read p0 var')
    end subroutine read_vcoord_normm

    recursive subroutine read_3d_normm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( idate < 2010010106 ) then
          y = (year / 10) * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version3,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        else if ( idate < 2015010106 ) then
          y = 2010
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version3,v%vname)) // &
            '201001010600-201501010000.nc'
        else if ( idate < 2022010106 ) then
          y = 2015
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)) // &
            '201501010600-202101010000.nc'
        else
          y = (year-2021)/10*10 + 2021
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_normm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_normm(v%ncid,v%vcoord%ak,v%vcoord%bk,v%vcoord%p0)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_normm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_normm

    recursive subroutine read_2d_normm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( idate < 2010010106 ) then
          y = (year / 10) * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version1,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        else if ( idate < 2015010106 ) then
          y = 2010
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version1,v%vname)) // &
            '201001010600-201501010000.nc'
        else if ( idate < 2021010106 ) then
          y = 2015
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)) // &
            '201501010600-202101010000.nc'
        else
          y = (year-2021)/10*10 + 2021
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_normm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_normm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_normm

    recursive subroutine read_fx_normm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(normm_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_normm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_normm

end module mod_cmip6_normm

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
