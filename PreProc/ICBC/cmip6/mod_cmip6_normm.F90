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
  use mod_kdinterp
  use mod_cmip6_helper
  use netcdf

  implicit none

  private

  integer(ik4) , parameter :: nhistory = 7
  integer(ik4) , parameter :: nscenario = 9
  integer(ik4) , parameter :: nfiles = nhistory + nscenario

  integer(ik8) , parameter , dimension(nfiles) :: nrmm_start = [ &
   195001010300, 196001010300, 197001010300, 198001010300,  &
   199001010300, 200001010300, 201001010300, 201501010600,  &
   202101010600, 203101010600, 204101010600, 205101010600,  &
   206101010600, 207101010600, 208101010600, 209101010600 ]
  integer(ik8) , parameter , dimension(nfiles) :: nrmm_end = [   &
   195912312100, 196912312100, 197912312100, 198912312100,  &
   199912312100, 200912312100, 201412312100, 202101010000,  &
   203101010000, 204101010000, 205101010000, 206101010000,  &
   207101010000, 208101010000, 209101010000, 210101010000 ]

  character(len=*) , parameter :: normm_version = 'v20191108'
  character(len=*) , parameter :: normm_version1 = 'v20200218'
  character(len=*) , parameter :: normm_version2 = 'v20210319'

  public :: read_3d_normm , read_2d_normm , read_fx_normm , read_sst_normm

  contains

    integer(ik4) function isequence(icode) result(ic)
      implicit none
      integer(ik8) , intent(in) :: icode
      integer(ik4) :: i
      do i = 1 , nfiles
        ic = i
        if ( icode >= nrmm_start(ic) .and. icode < nrmm_end(ic) ) exit
      end do
    end function isequence

    character(len=1024) function fname(vname,idate,offset)
      implicit none
      character(len=*) , intent(in) :: vname
      type(rcm_time_and_date) , intent(in) :: idate
      type(rcm_time_interval) , intent(in) , optional :: offset
      type(rcm_time_and_date) :: hdate
      integer(ik4) :: year , month , day , hour
      integer(ik4) :: iseq
      integer(ik8) :: icode

      if ( present(offset) ) then
        hdate = idate + offset
      else
        hdate = idate
      end if
      call split_idate(hdate, year, month, day, hour)
      icode = int(year,8)*100000000 + int(month,8)*1000000 + &
              int(day,8)*10000 + int(hour,8)*100
      iseq = isequence(icode)
      if ( idate < 2015010112 ) then
        write(fname,'(a,i12,a,i12,a)') &
              trim(cmip6_path(int(nrmm_start(iseq)/100000000),&
                              '6hrLev',normm_version,vname)), &
              nrmm_start(iseq), '-', nrmm_end(iseq), '.nc'
      else
        !write(fname,'(a,i12,a,i12,a)') &
        !      trim(cmip6_path(int(nrmm_start(iseq)/100000000),&
        !                      '6hrLev',normm_version1,vname)), &
        !      nrmm_start(iseq), '-', nrmm_end(iseq), '.nc'
        write(fname,'(a,i12,a,i12,a)') &
              trim(cmip6_path(int(nrmm_start(iseq)/100000000),&
                              '6hrLev',normm_version2,vname)), &
              nrmm_start(iseq), '-', nrmm_end(iseq), '.nc'
      end if
    end function fname

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

    subroutine read_hcoord_sst_normm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'i',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon i')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire i dim')
      istatus = nf90_inq_dimid(ncid,'j',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find j dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire j dim')
      call getmem2d(lon,1,nlon,1,nlat,'cmip6_norm:lon')
      call getmem2d(lat,1,nlon,1,nlat,'cmip6_norm:lat')
      istatus = nf90_inq_varid(ncid,'longitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find longitude var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read longitude var')
      istatus = nf90_inq_varid(ncid,'latitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find latitude var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read latitude var')
    end subroutine read_hcoord_sst_normm

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
      character(len=1024) :: seqfile
      integer :: ncid_next , ivar_next , icloseme
      real(rkx) , dimension(:,:,:) , allocatable :: v1

      tdif = rcm_time_interval(3,uhrs)
      if ( v%ncid == -1 ) then
        if ( idate < 2015010112 ) then
          v%filename = fname(v%vname,idate,rcm_time_interval(3,uhrs))
        else
          v%filename = fname(v%vname,idate)
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
          'Error searching '//trim(v%vname)//' var in file '// &
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

      if ( idate > 2015010106 ) then
        tdif = idate - v%first_date
        irec = nint(tohours(tdif)/6.0) + 1
      else
        allocate(v1,source=v%var)
        tdif = idate - v%first_date
        irec = int(tohours(tdif)/6.0) + 1
      end if

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
      if ( year < 2015 ) then
        if ( irec+1 <= v%nrec ) then
          ncid_next = v%ncid
          ivar_next = v%ivar
          istart(4) = irec + 1
          seqfile = v%filename
          icloseme = 0
        else
          seqfile = fname(v%vname,idate,rcm_time_interval(9,uhrs))
#ifdef DEBUG
          write(stderr,*) 'Opening ',trim(seqfile)
#endif
          istatus = nf90_open(seqfile,nf90_nowrite,ncid_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error opening file '//trim(seqfile)//'.')
          istatus = nf90_inq_varid(ncid_next,v%vname,ivar_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error searching '//trim(v%vname)//' var in file '// &
            trim(v%filename)//'.')
          istart(4) = 1
          icloseme = 1
        end if
        istatus = nf90_get_var(ncid_next,ivar_next,v1,istart,icount)
        call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error read variable '//v%vname//' from '//trim(seqfile)//'.')
        if ( idate /= 2015010100 ) then
          v%var = 0.5 * (v%var + v1)
        else
          v%var = 0.4 * v%var + 0.6 * v1
        end if
        if ( icloseme > 0 ) then
          istatus = nf90_close(ncid_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error close file '//trim(seqfile)//'.')
        end if
        deallocate(v1)
      end if
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
      character(len=1024) :: seqfile
      integer :: ncid_next , ivar_next , icloseme
      real(rkx) , dimension(:,:) , allocatable :: v1

      tdif = rcm_time_interval(3,uhrs)
      if ( v%ncid == -1 ) then
        if ( idate < 2015010112 ) then
          v%filename = fname(v%vname,idate,rcm_time_interval(3,uhrs))
        else
          v%filename = fname(v%vname,idate)
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
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching '//trim(v%vname)//' var in file '// &
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

      if ( idate > 2015010106 ) then
        tdif = idate - v%first_date
        irec = nint(tohours(tdif)/6.0) + 1
      else
        allocate(v1,source=v%var)
        tdif = idate - v%first_date
        irec = int(tohours(tdif)/6.0) + 1
      end if

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
      if ( year < 2015 ) then
        if ( irec+1 <= v%nrec ) then
          ncid_next = v%ncid
          ivar_next = v%ivar
          istart(3) = irec + 1
          seqfile = v%filename
          icloseme = 0
        else
          seqfile = fname(v%vname,idate,rcm_time_interval(9,uhrs))
#ifdef DEBUG
          write(stderr,*) 'Opening ',trim(seqfile)
#endif
          istatus = nf90_open(seqfile,nf90_nowrite,ncid_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error opening file '//trim(seqfile)//'.')
          istatus = nf90_inq_varid(ncid_next,v%vname,ivar_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error searching '//trim(v%vname)//' var in file '// &
            trim(v%filename)//'.')
          istart(3) = 1
          icloseme = 1
        end if
        istatus = nf90_get_var(ncid_next,ivar_next,v1,istart,icount)
        call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error read variable '//v%vname//' from '//trim(seqfile)//'.')
        if ( idate /= 2015010100 ) then
          v%var = 0.5 * (v%var + v1)
        else
          v%var = 0.4 * v%var + 0.6 * v1
        end if
        if ( icloseme > 0 ) then
          istatus = nf90_close(ncid_next)
          call cmip6_error(istatus,__FILE__,__LINE__, &
            'Error close file '//trim(seqfile)//'.')
        end if
        deallocate(v1)
      end if
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
        'Error searching '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_normm

    recursive subroutine read_sst_normm(idate,v,lat,lon)
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
        if ( year < 2010 ) then
          y = (year / 10) * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'Oday',normm_version,v%vname)), &
            y, '0101-', y+9, '1231.nc'
        else if ( year >= 2010 .and. year < 2015 ) then
          write(v%filename,'(a,a,a)') &
            trim(cmip6_path(year,'Oday',normm_version,v%vname)), &
            '20100101-', '20141231.nc'
        else if ( year >= 2015 .and. year < 2021 ) then
          write(v%filename,'(a,a,a)') &
            trim(cmip6_path(year,'Oday',normm_version,v%vname)), &
            '20150101-', '20201231.nc'
        else
          y = ((year-1) / 10) * 10 + 1
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'Oday',normm_version,v%vname)), &
            y, '0101-', y+9, '1231.nc'
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
        call read_hcoord_sst_normm(v%ncid,v%hcoord%lon2d,v%hcoord%lat2d)
        call h_interpolator_create(v%hint(1), &
                                   v%hcoord%lat2d,v%hcoord%lon2d,lat,lon)
        call getmem2d(v%var,1,size(v%hcoord%lon2d,1), &
                            1,size(v%hcoord%lat2d,2), &
                            'cmip6_normm:'//trim(v%vname))
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,trim(v%vname),v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching '//trim(v%vname)// &
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
        call read_sst_normm(idate,v,lat,lon)
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
    end subroutine read_sst_normm

end module mod_cmip6_normm

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
