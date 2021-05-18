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
  use mod_dynparam
  use mod_kdinterp
  use netcdf

  implicit none

  private

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
    type(rcm_time_and_date) :: first_date
    type(h_interpolator) , pointer , dimension(:) :: hint
  end type cmip6_file

  type, extends(cmip6_file) :: cmip6_2d_var
    character(len=8) :: vname
    real(rkx) , pointer , dimension(:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
  end type cmip6_2d_var

  type, extends(cmip6_file) :: cmip6_3d_var
    character(len=8) :: vname
    real(rkx) , pointer , dimension(:,:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
    type(cmip6_vertical_coordinate) , pointer :: vcoord => null( )
  end type cmip6_3d_var

  public :: cmip6_2d_var , cmip6_3d_var
  public :: cmip6_path , cmip6_error

  contains

    character(len=1024) function cmip6_path(year,freq,ver,var) result(fpath)
      implicit none
      character(len=*) , intent(in) :: var , freq , ver
      integer(ik4) , intent(in) :: year
      fpath = trim(inpglob)//pthsep//'cmip6'//pthsep
      if ( year < 2015 ) then
        fpath = trim(fpath)//'CMIP'//pthsep
      else
        fpath = trim(fpath)//'ScenarioMIP'//pthsep
      end if
      select case ( cmip6_model )
        case ( 'MPI-ESM1-2-HR' )
          if ( year < 2015 ) then
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-HR'//pthsep
          else
            fpath = trim(fpath)//'DKRZ'//pthsep//'MPI-ESM1-2-HR'//pthsep
          end if
        case default
          call die(__FILE__, &
            '__LINE__ : Unsupported cmip6 model: '//trim(cmip6_model),-1)
      end select
      if ( year < 2015 ) then
        fpath = trim(fpath)//'historical'//pthsep
      else
        fpath = trim(fpath)//trim(cmip6_ssp)//pthsep
      end if
      fpath = trim(fpath)//trim(cmip6_variant)//pthsep//trim(freq)//pthsep// &
        trim(var)//pthsep//trim(cmip6_grid)//pthsep//trim(ver)// &
        pthsep//trim(var)//'_'//trim(freq)//'_'//trim(cmip6_model)//'_'
      if ( year < 2015 ) then
        fpath = trim(fpath)//'historical'//'_'
      else
        fpath = trim(fpath)//trim(cmip6_ssp)//'_'
      end if
      fpath = trim(fpath)//trim(cmip6_variant)//'_'//trim(cmip6_grid)//'_'
    end function cmip6_path

    subroutine cmip6_error(ival,filename,line,arg)
      implicit none
      integer(ik4) , intent(in) :: ival , line
      character(len=8) :: cline
      character(*) , intent(in) :: filename , arg
      if ( ival /= nf90_noerr ) then
        write (cline,'(i8)') line
        write (stderr,*) nf90_strerror(ival)
        call die(filename,trim(cline)//':'//arg,ival)
      end if
    end subroutine cmip6_error

end module mod_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
