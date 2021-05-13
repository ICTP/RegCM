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

  public :: cmip6_read_record

  type cmip6_horizontal_coordinates
    real(rkx) , pointer , dimension(:) :: lon1d => null( )
    real(rkx) , pointer , dimension(:) :: lat1d => null( )
    real(rkx) , pointer , dimension(:,:) :: lon2d => null( )
    real(rkx) , pointer , dimension(:,:) :: lat2d => null( )
  end type cmip6_horizontal_coordinates

  type cmip6_vertical_coordinate
    real(rkx) , pointer , dimension(:) :: plev => null( )
    real(rkx) , pointer , dimension(:) :: ak => null( )
    real(rkx) , pointer , dimension(:) :: bk => null( )
    real(rkx) , pointer , dimension(:,:) :: topo => null( )
  end type cmip6_vertical_coordinate

  type cmip6_file
    character(len=1024) :: filename
    integer(ik4) :: ncid = -1
    integer(ik4) :: ivar = -1
    integer(ik4) :: nrec = -1
    real(rk8) :: freq
    type(rcm_time_and_date) , dimension(2) :: dates
  end type cmip6_file

  type, extends(cmip6_file) :: cmip6_2d_var
    real(rkx) , pointer , dimension(:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
  end type cmip6_2d_var

  type, extends(cmip6_file) :: cmip6_3d_var
    real(rkx) , pointer , dimension(:,:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
    type(cmip6_vertical_coordinate) , pointer :: vcoord => null( )
  end type cmip6_3d_var

  type(cmip6_2d_var) , pointer :: sst => null( )
  type(cmip6_2d_var) , pointer :: ps => null( )
  type(cmip6_3d_var) , pointer :: ta => null( )
  type(cmip6_3d_var) , pointer :: qa => null( )
  type(cmip6_3d_var) , pointer :: ua => null( )
  type(cmip6_3d_var) , pointer :: va => null( )
  type(cmip6_3d_var) , pointer :: zg => null( )

  real(rkx) , pointer , dimension(:,:,:) :: pa

  contains

    subroutine cmip6_read_record(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      select case (cmip6_model)
        case ('MPI-ESM1-2-HR')
          call read_ps_mpihr(idate)
        case default
          call die(__FILE__,'__LINE__ : Unsupported cmip6 model.',-1)
      end select
    end subroutine cmip6_read_record

    subroutine read_hcoord_mpihr(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_mpi:lon')
      call getmem1d(lat,1,nlat,'cmip6_mpi:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_mpihr

    subroutine read_vcoord_mpihr(ncid,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(ap,1,nlev,'cmip6_mpi:ap')
      call getmem1d(b,1,nlev,'cmip6_mpi:b')
      istatus = nf90_inq_varid(ncid,'ap',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find ap var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read ap var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_vcoord_mpihr

    recursive subroutine read_ps_mpihr(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=6) :: sspname
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( .not. associated(ps) ) then
        allocate(ps)
      end if
      if ( ps%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = (year / 5) * 5
        if ( month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 5
        end if
        if ( year >= 2015 ) then
          write(ps%filename,'(a,i4,a,i4,a)') &
            trim(inpglob)//pthsep//'CMIP6'//pthsep//trim(cmip6_model)// &
            pthsep//'historical'//pthsep//'ps'//pthsep// &
            'ps_6hrLev_'//trim(cmip6_model)//'_historical_'// &
             trim(ssp_variant)//'_'//trim(ssp_grid)//'_', &
             y, '01010600-', y+5, '01010000.nc'
        else
          write(ps%filename,'(a,i4,a,i4,a)') &
            trim(inpglob)//pthsep//'CMIP6'//pthsep//trim(cmip6_model)// &
            pthsep//trim(ssp_code)//pthsep//'ps'//pthsep// &
            'ps_6hrLev_'//trim(cmip6_model)//'_'//trim(ssp_code)//'_'// &
            trim(ssp_variant)//'_'//trim(ssp_grid)//'_', &
            y, '01010600-', y+5, '01010000.nc'
        end if
        istatus = nf90_open(ps%filename,nf90_nowrite,ps%ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(ps%filename)//'.')
      end if
      if ( .not. associated(ps%hcoord) ) then
        allocate(ps%hcoord)
        call read_hcoord_mpihr(ps%ncid,ps%hcoord%lon1d,ps%hcoord%lat1d)
        call getmem2d(ps%var,1,size(ps%hcoord%lon1d), &
                             1,size(ps%hcoord%lat1d),'cmip6_mpi:ps')
      end if
      if ( ps%ivar == -1 ) then
        istatus = nf90_inq_varid(ps%ncid,'ps',ps%ivar)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error searchong ps var in file '//trim(ps%filename)//'.')
      end if
      if ( ps%nrec == -1 ) then
        istatus = nf90_inq_dimid(ps%ncid,'time',idimid)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(ps%ncid,idimid,len=ps%nrec)
        call checkncerr(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(ps%ncid,'time',it)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(ps%filename)//'.')
        istatus = nf90_get_att(ps%ncid,it,"calendar",timecal)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(ps%filename)//'.')
        istatus = nf90_get_att(ps%ncid,it,"units",timeunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(ps%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(ps%ncid,it,times,istart(1:1),icount(1:1))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(ps%filename)//'.')
        ps%dates(1) = timeval2date(times(1),timeunit,timecal)
        ps%dates(2) = timeval2date(times(2),timeunit,timecal)
        tdif = ps%dates(2) - ps%dates(1)
        ps%freq = tohours(tdif)
      end if

      tdif = idate - ps%dates(1)
      irec = nint(tohours(tdif)/ps%freq) + 1

      if ( irec > ps%nrec ) then
        istatus = nf90_close(ps%ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(ps%filename)//'.')
        ps%ncid = -1
        ps%ivar = -1
        ps%nrec = -1
        ps%freq = -1
        irec = 1
        call read_ps_mpihr(idate)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(ps%var,1)
      icount(2) = size(ps%var,2)
      icount(3) = 1
      istatus = nf90_get_var(ps%ncid,ps%ivar,ps%var,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
          'Error read variable ps from '//trim(ps%filename)//'.')
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
