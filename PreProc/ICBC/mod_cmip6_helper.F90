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

module mod_cmip6_helper

  use mod_intkinds
  use mod_realkinds
  use mod_cmip6
  use mod_message
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_kdinterp
  use netcdf

  implicit none

  private

  public :: init_cmip6 , get_cmip6 , conclude_cmip6

  type(cmip6_2d_var) , pointer :: ps => null( )
  type(cmip6_3d_var) , pointer :: ta => null( )
  type(cmip6_3d_var) , pointer :: qa => null( )
  type(cmip6_3d_var) , pointer :: ua => null( )
  type(cmip6_3d_var) , pointer :: va => null( )
  type(cmip6_3d_var) , pointer :: zg => null( )

  contains

    subroutine init_cmip6(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      select case (cmip6_model)
        case ('MPI-ESM1-2-HR')
          allocate(ua,va,ta,qa,zg)
          call read_3d_mpihr(idate,ua)
          if ( idynamic == 3 ) then
            allocate(ua%hint(3))
            call h_interpolator_create(ua%hint(1),ua%hcoord%lat1d, &
              ua%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(ua%hint(2),ua%hcoord%lat1d, &
              ua%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(ua%hint(3),ua%hcoord%lat1d, &
              ua%hcoord%lon1d, vlat, vlon)
          else
            allocate(ua%hint(2))
            call h_interpolator_create(ua%hint(1),ua%hcoord%lat1d, &
              ua%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(ua%hint(2),ua%hcoord%lat1d, &
              ua%hcoord%lon1d, dlat, dlon)
          end if
          va%hint => ua%hint
          ta%hint => ua%hint
          qa%hint => ua%hint
          zg%hint => ua%hint
        case default
          call die(__FILE__,'__LINE__ : Unsupported cmip6 model.',-1)
      end select
    end subroutine init_cmip6

    subroutine get_cmip6(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      select case (cmip6_model)
        case ('MPI-ESM1-2-HR')
          call read_3d_mpihr(idate,ua)
          call read_3d_mpihr(idate,va)
          call read_3d_mpihr(idate,ta)
          call read_3d_mpihr(idate,qa)
          call read_3d_mpihr(idate,zg)
        case default
          call die(__FILE__,'__LINE__ : Unsupported cmip6 model.',-1)
      end select
    end subroutine get_cmip6

    subroutine conclude_cmip6( )
      implicit none
      if ( idynamic == 3 ) then
        call h_interpolator_destroy(ua%hint(1))
        call h_interpolator_destroy(ua%hint(2))
        call h_interpolator_destroy(ua%hint(3))
      else
        call h_interpolator_destroy(ua%hint(1))
        call h_interpolator_destroy(ua%hint(2))
      end if
      deallocate(ua%hint)
    end subroutine conclude_cmip6

    subroutine read_hcoord_mpihr(ncid,lon,lat)
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
    end subroutine read_hcoord_mpihr

    subroutine read_vcoord_mpihr(ncid,plev)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: plev
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'plev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find plev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire plev dim')
      call getmem1d(plev,1,nlev,'cmip6_mpi:plev')
      istatus = nf90_inq_varid(ncid,'plev',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find plev var')
      istatus = nf90_get_var(ncid,ivarid,plev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read plev var')
    end subroutine read_vcoord_mpihr

    recursive subroutine read_3d_mpihr(idate,v)
      implicit none
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      type(rcm_time_and_date) , intent(in) :: idate
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = (year / 5) * 5
        if ( month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 5
        end if
        write(v%filename,'(a,i4,a,i4,a)') &
          trim(cmip6_path(y,'6hrPlevPt',cmip6_version,v%vname)), &
          y, '01010600-', y+5, '01010000.nc'
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        allocate(v%vcoord)
        call read_hcoord_mpihr(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
        call read_vcoord_mpihr(v%ncid,v%vcoord%plev)
        call getmem3d(v%var,1,size(v%hcoord%lon1d), &
                            1,size(v%hcoord%lat1d), &
                            1,size(v%vcoord%plev), &
                            'cmip6_mpi:'//trim(v%vname))
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
        call read_3d_mpihr(idate,v)
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
    end subroutine read_3d_mpihr

end module mod_cmip6_helper

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
