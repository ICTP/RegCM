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

module mod_cmip6_cesm

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

  character(len=*), parameter :: cesm_version = 'v20190514'
  character(len=*), parameter :: cesm_version1 = 'v20200528'
  character(len=*), parameter :: cesm_version2 = 'v20200210'

  public :: read_3d_cesm, read_2d_cesm, read_fx_cesm, read_sst_cesm

  contains

    subroutine read_hcoord_cesm(ncid,lon,lat)
      implicit none
      integer(ik4), intent(in) :: ncid
      real(rkx), pointer, contiguous, dimension(:), intent(inout) :: lon, lat
      integer(ik4) :: istatus, idimid, ivarid
      integer(ik4) :: nlon, nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_cesm:lon')
      call getmem1d(lat,1,nlat,'cmip6_cesm:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_cesm

    subroutine read_hcoord_sst_cesm(ncid,lon,lat)
      implicit none
      integer(ik4), intent(in) :: ncid
      real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: lon, lat
      integer(ik4) :: istatus, idimid, ivarid
      integer(ik4) :: nlon, nlat
      istatus = nf90_inq_dimid(ncid,'nlon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon nlon')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire nlon dim')
      istatus = nf90_inq_dimid(ncid,'nlat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find nlat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire nlat dim')
      call getmem2d(lon,1,nlon,1,nlat,'cmip6_cesm:lon')
      call getmem2d(lat,1,nlon,1,nlat,'cmip6_cesm:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_sst_cesm

    subroutine read_vcoord_cesm(ncid,a,b,p0)
      implicit none
      integer(ik4), intent(in) :: ncid
      real(rkx), pointer, contiguous, dimension(:), intent(inout) :: a, b
      real(rkx), intent(out) :: p0
      integer(ik4) :: istatus, idimid, ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(a,1,nlev,'cmip6_cesm:a')
      call getmem1d(b,1,nlev,'cmip6_cesm:b')
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
    end subroutine read_vcoord_cesm

    recursive subroutine read_3d_cesm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date), intent(in) :: idate
      type(cmip6_3d_var), pointer, intent(inout) :: v
      logical, optional, intent(in) :: lonlyc
      integer(ik4) :: istatus, idimid, it, irec
      integer(ik4) :: year, month, day, hour, y1, y2
      character(len=32) :: timecal, timeunit
      integer(ik4), dimension(4) :: istart, icount
      real(rk8), dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year < 2010 ) then
          if ( v%vname /= 'ta' ) then
            ver = cesm_version
          else
            ver = cesm_version2
          end if
          y1 = year/10*10
          y2 = y1+9
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',ver,v%vname)), &
            y1, '01010000-', y2, '12311800.nc'
        else if ( year >= 2015 .and. year < 2095 ) then
          y1 = (year-2015)/10*10 + 2015
          y2 = y1+9
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
            y1, '01010000-', y2, '12311800.nc'
        else
          if ( idate < 2015010100 ) then
            if ( v%vname /= 'ta' ) then
              ver = cesm_version
            else
              ver = cesm_version2
            end if
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',ver,v%vname)), &
              '201001010000-201501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
              '209501010000-210101010000.nc'
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
        call read_hcoord_cesm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_cesm(v%ncid,v%vcoord%ak,v%vcoord%bk,v%vcoord%p0)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:cesm:'//trim(v%vname))
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
        call read_3d_cesm(idate,v)
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
    end subroutine read_3d_cesm

    recursive subroutine read_2d_cesm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date), intent(in) :: idate
      type(cmip6_2d_var), pointer, intent(inout) :: v
      logical, optional, intent(in) :: lonlyc
      integer(ik4) :: istatus, idimid, it, irec
      integer(ik4) :: y1
      integer(ik4) :: year, month, day, hour
      character(len=32) :: timecal, timeunit
      integer(ik4), dimension(3) :: istart, icount
      real(rk8), dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year < 2010 ) then
          y1 = year/10 * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version,v%vname)), &
            y1,'01010000-',y1+9,'12311800.nc'
        else if ( year >= 2015 .and. year < 2095 ) then
          y1 = (year-2015)/10*10 + 2015
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
            y1,'01010000-',y1+9,'12311800.nc'
        else
          if ( idate < 2015010100 ) then
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version,v%vname)), &
              '201001010000-201501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
              '209501010000-210101010000.nc'
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
        call read_hcoord_cesm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cesm:'//trim(v%vname))
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
        call read_2d_cesm(idate,v)
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
    end subroutine read_2d_cesm

    recursive subroutine read_fx_cesm(v)
      implicit none
      type(cmip6_2d_var), pointer, intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(cesm_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_cesm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cesm:'//trim(v%vname))
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
    end subroutine read_fx_cesm

    recursive subroutine read_sst_cesm(idate,v,lat,lon)
      implicit none
      type(rcm_time_and_date), intent(in) :: idate
      type(cmip6_2d_var), intent(inout) :: v
      real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: lon, lat
      integer(ik4) :: istatus, idimid, it, irec
      integer(ik4) :: year, month, day, hour
      character(len=32) :: timecal, timeunit
      integer(ik4), dimension(3) :: istart, icount
      real(rk8), dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( idate < 1900010200 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version,v%vname)), &
            '18500102-19000101.nc'
        else if ( idate < 1950010200 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version,v%vname)), &
            '19000102-19500101.nc'
        else if ( idate < 2000010200 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version,v%vname)), &
            '19500102-20000101.nc'
        else if ( idate < 2015010200 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version,v%vname)), &
            '20000102-20150101.nc'
        else if ( idate < 2065010200 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version1,v%vname)), &
            '20150102-20650101.nc'
        else
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'Oday',cesm_version1,v%vname)), &
            '20650102-21010101.nc'
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(trim(v%filename),nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_sst_cesm(v%ncid,v%hcoord%lon2d,v%hcoord%lat2d)
        call h_interpolator_create(v%hint(1), &
                                   v%hcoord%lat2d,v%hcoord%lon2d,lat,lon)
        call getmem2d(v%var,1,size(v%hcoord%lon2d,1), &
                            1,size(v%hcoord%lat2d,2), &
                            'cmip6_cesm:'//trim(v%vname))
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
      irec = int((tohours(tdif)+12)/24) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
        call read_sst_cesm(idate,v,lat,lon)
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
    end subroutine read_sst_cesm

end module mod_cmip6_cesm

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
