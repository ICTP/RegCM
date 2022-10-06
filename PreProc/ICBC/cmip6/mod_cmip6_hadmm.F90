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

module mod_cmip6_hadmm

  use mod_intkinds
  use mod_realkinds
  use mod_message
  use mod_date
  use mod_stdio
  use mod_memutil
  use mod_kdinterp
  use mod_cmip6_helper
  use netcdf

  implicit none

  private

  character(len=*) , parameter :: hadmm_version = 'v20191207'
  character(len=*) , parameter :: hadmm_version1 = 'v20210331'
  character(len=*) , parameter :: hadmm_version2 = 'v20201019'
  character(len=*) , parameter :: hadmm_version3 = 'v20201113'
  character(len=*) , parameter :: hadmm_version4 = 'v20200515'

  public :: read_3d_hadmm , read_2d_hadmm , read_fx_hadmm , read_sst_hadmm

  contains

    subroutine read_hcoord_hadmm(ncid,lon,lat)
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
      call getmem1d(lon,1,nlon,'cmip6_mpi:lon')
      call getmem1d(lat,1,nlat,'cmip6_mpi:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_hadmm

    subroutine read_hcoord_sst_hadmm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'i',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'j',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem2d(lon,1,nlon,1,nlat,'cmip6_had:lon')
      call getmem2d(lat,1,nlon,1,nlat,'cmip6_had:lat')
      istatus = nf90_inq_varid(ncid,'longitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find longitude var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read longitude var')
      istatus = nf90_inq_varid(ncid,'latitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find latitude var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read latitude var')
    end subroutine read_hcoord_sst_hadmm

    subroutine read_vcoord_hadmm(ncid,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(ap,1,nlev,'cmip6_mpi:lev')
      call getmem1d(b,1,nlev,'cmip6_mpi:b')
      istatus = nf90_inq_varid(ncid,'lev',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lev var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_vcoord_hadmm

    recursive subroutine read_3d_hadmm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y , yp , m , mp
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = year
        m = month / 3 * 3 + 1
        if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 1
          m = 10
        end if
        yp = y
        mp = m + 3
        if ( mp > 12 ) then
          yp = yp + 1
          mp = 1
        end if
        if ( y > 2015 ) then
          ver = hadmm_version2
        else
          ver = hadmm_version3
        end if
        write(v%filename,'(a,i4,i0.2,a,i4,i0.2,a)') &
          trim(cmip6_path(y,'6hrLev',ver,v%vname)), &
          y, m, '010600-', yp, mp, '010000.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_hadmm(v%ncid,v%vcoord%ak,v%vcoord%bk)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:hadmm:'//trim(v%vname))
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
        call read_3d_hadmm(idate,v)
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
    end subroutine read_3d_hadmm

    recursive subroutine read_2d_hadmm(idate,v,lonlyc)
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
        if ( year >= 2015 .and. year < 2020 ) then
          if ( year == 2015 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2010,'6hrLev',hadmm_version2,v%vname)), &
              '201001010600-', '201501010000.nc'
          else
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2015,'6hrLev',hadmm_version3,v%vname)), &
              '201501010600-', '202001010000.nc'
          end if
        else if ( year == 2020 .and. month == 1 .and. &
                 day == 1 .and. hour == 0 ) then
          write(v%filename,'(a,a,a)') &
            trim(cmip6_path(2015,'6hrLev',hadmm_version3,v%vname)), &
            '201501010600-', '202001010000.nc'
        else if ( year >= 2010 .and. year < 2015 ) then
          if ( year == 2010 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2000,'6hrLev',hadmm_version2,v%vname)), &
              '200001010600-', '201001010000.nc'
          else
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2010,'6hrLev',hadmm_version2,v%vname)), &
              '201001010600-', '201501010000.nc'
          end if
        else
          y = (year / 10) * 10
          if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
            y = y - 10
          end if
          if ( y > 2015 ) then
            write(v%filename,'(a,i4,a,i4,a)') &
              trim(cmip6_path(y,'6hrLev',hadmm_version3,v%vname)), &
              y, '01010600-', y+10, '01010000.nc'
          else
            write(v%filename,'(a,i4,a,i4,a)') &
              trim(cmip6_path(y,'6hrLev',hadmm_version2,v%vname)), &
              y, '01010600-', y+10, '01010000.nc'
          end if
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
        call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:hadmm:'//trim(v%vname))
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
        call read_2d_hadmm(idate,v)
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
    end subroutine read_2d_hadmm

    recursive subroutine read_fx_hadmm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(hadmm_version1,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:hadmm:'//trim(v%vname))
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
    end subroutine read_fx_hadmm

    recursive subroutine read_sst_hadmm(idate,v,lat,lon)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , intent(inout) :: v
      real(rkx) , pointer , dimension(:,:) , intent(in) :: lat , lon
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = (year / 5) * 5
        write(v%filename,'(a,i4,a,i4,a)') &
          trim(cmip6_path(y,'Oday',hadmm_version4,v%vname)), &
          y, '0101-', y+4, '1230.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_sst_hadmm(v%ncid,v%hcoord%lon2d,v%hcoord%lat2d)
        call h_interpolator_create(v%hint(1), &
                                   v%hcoord%lat2d,v%hcoord%lon2d,lat,lon)
        call getmem2d(v%var,1,size(v%hcoord%lon2d,1), &
                            1,size(v%hcoord%lat2d,2), &
                            'cmip6_had:'//trim(v%vname))
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,trim(v%vname),v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)// &
          ' var in file '//trim(v%filename)//'.')
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
      irec = nint((tohours(tdif)+12)/24) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
        call read_sst_hadmm(idate,v,lat,lon)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//trim(v%vname)// &
          ' from '//trim(v%filename)//'.')
      where ( v%var > 0.9E+20_rkx )
        v%var = -9999.0_rkx
      else where
        v%var = v%var + 273.15_rkx
      end where
    end subroutine read_sst_hadmm

end module mod_cmip6_hadmm

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
