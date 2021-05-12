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

module mod_cmip6

  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_message
  use netcdf
  use mod_dynparam
  use mod_memutil

  implicit none

  integer(ik4) , parameter :: nmodels = 1

  character(len=*) , parameter , dimension(nmodels) :: supported_models = &
    [ 'MPI-ESM1-2-HR' ]

  interface read_coord
    module procedure read_coord_mpihr
  end interface read_coord

  public :: cmip6_readrecord

  contains

    subroutine cmip6_readrecord(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
    end subroutine cmip6_readrecord

    subroutine read_coord_mpihr(ncid,lon,lat,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat , ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat , nlev
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(lon,1,nlon,'cmip6_mpi:lon')
      call getmem1d(lat,1,nlat,'cmip6_mpi:lat')
      call getmem1d(ap,1,nlev,'cmip6_mpi:ap')
      call getmem1d(b,1,nlev,'cmip6_mpi:b')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lat var')
      istatus = nf90_inq_varid(ncid,'ap',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find ap var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read ap var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_coord_mpihr

    recursive subroutine read_ps_mpihr(idate,ps)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      real(rkx) , pointer , dimension(:,:) :: ps
      integer(ik4) , save :: ncid = -1
      integer(ik4) , save :: ivar = -1
      integer(ik4) , save :: irec = -1
      integer(ik4) , save :: nrec = -1
      character(len=1024) , save :: psfilename
      real(rk8) , save :: freq
      type(rcm_time_and_date) , dimension(2) , save :: dates
      integer(ik4) :: istatus , idimid , it
      integer(ik4) :: year , month , day , hour , y
      character(len=6) :: sspname
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = (year / 5) * 5
        if ( month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 5
        end if
        if ( year >= 2015 ) then
          write(psfilename,'(a,i4,a,i4,a)') &
            trim(inpglob)//pthsep//'CMIP6'//pthsep//'MPI-ESM1-2-HR'// &
            pthsep//'ssp'//pthsep//'ps'//pthsep// &
            'ps_6hrLev_MPI-ESM1-2-HR_historical_'// &
             trim(ssp_variant)//'_'//trim(ssp_grid)//'_', &
             y, '01010600-', y+5, '01010000.nc'
        else
          write(psfilename,'(a,i4,a,i4,a)') &
            trim(inpglob)//pthsep//'CMIP6'//pthsep//'MPI-ESM1-2-HR'// &
            pthsep//trim(ssp_code)//pthsep//'ps'//pthsep// &
            'ps_6hrLev_MPI-ESM1-2-HR_'//trim(ssp_code)//'_'// &
            trim(ssp_variant)//'_'//trim(ssp_grid)//'_', &
            y, '01010600-', y+5, '01010000.nc'
        end if
        istatus = nf90_open(psfilename,nf90_nowrite,ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(psfilename)//'.')
      end if
      if ( ivar == -1 ) then
        istatus = nf90_inq_varid(ncid,'ps',ivar)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error searchong ps var in file '//trim(psfilename)//'.')
      end if
      if ( nrec == -1 ) then
        istatus = nf90_inq_dimid(ncid,'time',idimid)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(ncid,idimid,len=nrec)
        call checkncerr(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(ncid,'time',it)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(psfilename)//'.')
        istatus = nf90_get_att(ncid,it,"calendar",timecal)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(psfilename)//'.')
        istatus = nf90_get_att(ncid,it,"units",timeunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(psfilename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(ncid,it,times,istart(1:1),icount(1:1))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(psfilename)//'.')
        dates(1) = timeval2date(times(1),timeunit,timecal)
        dates(2) = timeval2date(times(2),timeunit,timecal)
        tdif = dates(2) - dates(1)
        freq = tohours(tdif)
      end if

      tdif = idate - dates(1)
      irec = nint(tohours(tdif)/freq) + 1

      if ( irec > nrec ) then
        istatus = nf90_close(ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(psfilename)//'.')
        ncid = -1
        ivar = -1
        nrec = -1
        irec = 1
        call read_ps_mpihr(idate,ps)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(ps,1)
      icount(2) = size(ps,2)
      icount(3) = 1
      istatus = nf90_get_var(ncid,ivar,ps,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
          'Error read variable ps from '//trim(psfilename)//'.')
    end subroutine read_ps_mpihr

    subroutine checkncerr(ival,filename,line,arg)
      implicit none
      integer(ik4) , intent(in) :: ival , line
      character(len=8) :: cline
      character(*) , intent(in) :: filename , arg
      if ( ival /= nf90_noerr ) then
        write (cline,'(i8)') line
        write (stderr,*) nf90_strerror(ival)
        call die(filename,trim(cline)//':'//arg,ival)
      end if
    end subroutine checkncerr

end module mod_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
