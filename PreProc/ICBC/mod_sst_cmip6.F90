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

module mod_sst_cmip6

  use mod_intkinds
  use mod_realkinds
  use mod_cmip6_helper
  use mod_message
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_kdinterp
  use mod_date
  use mod_stdio
  use netcdf

  implicit none

  private

  public :: cmip6_sst

  type(cmip6_2d_var) :: sst

  character(len=*) , parameter , public :: mpihr_version = 'v20190710'

  abstract interface
    subroutine read_cmip6_sst(id,var)
      import
      type(rcm_time_and_date) , intent(in) :: id
      type(cmip6_2d_var) , intent(inout) :: var
    end subroutine read_cmip6_sst
  end interface

  contains

    subroutine cmip6_sst
      implicit none
      type(rcm_time_and_date) :: idate , idatef , idateo
      type(rcm_time_interval) :: tdif , step
      procedure(read_cmip6_sst) , pointer :: read_func
      integer :: nsteps , n

      idateo = globidate1
      idatef = globidate2
      tdif = idatef-idateo

      select case (cmip6_model)
        case ('MPI-ESM1-2-HR')
          read_func => read_sst_mpihr
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case default
      end select

      write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
      write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
      write (stdout,*) 'NSTEPS = ', nsteps

      call open_sstfile(idateo)

      allocate(sst%hint(1))

      idate = idateo
      do n = 1 , nsteps
        call read_func(idate,sst)
        call h_interpolate_cont(sst%hint(1),sst%var,sstmm)
        call writerec(idate)
        write (stdout,*) 'WRITEN OUT SST DATA : ' , tochar(idate)
        idate = idate + step
      end do

      call h_interpolator_destroy(sst%hint(1))
      deallocate(sst%hint)
    end subroutine cmip6_sst

    subroutine read_hcoord_sst_mpihr(ncid,lon,lat)
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
      call getmem2d(lon,1,nlon,1,nlat,'cmip6_mpi:lon')
      call getmem2d(lat,1,nlon,1,nlat,'cmip6_mpi:lat')
      istatus = nf90_inq_varid(ncid,'longitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find longitude var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read longitude var')
      istatus = nf90_inq_varid(ncid,'latitude',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find latitude var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read latitude var')
    end subroutine read_hcoord_sst_mpihr

    recursive subroutine read_sst_mpihr(idate,v)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , intent(inout) :: v
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
          trim(cmip6_path(y,'Oday',mpihr_version,v%vname)), &
          y, '0101-', y+4, '1231.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_sst_mpihr(v%ncid,v%hcoord%lon2d,v%hcoord%lat2d)
        call h_interpolator_create(v%hint(1), &
                                   v%hcoord%lat2d,v%hcoord%lon2d,xlat,xlon)
        call getmem2d(v%var,1,size(v%hcoord%lon2d,1), &
                            1,size(v%hcoord%lat2d,2), &
                            'cmip6_mpi:'//trim(v%vname))
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
        call read_sst_mpihr(idate,v)
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
    end subroutine read_sst_mpihr

end module mod_sst_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
