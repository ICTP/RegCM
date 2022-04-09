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
  use mod_date
  use mod_message
  use mod_dynparam
  use mod_kdinterp
  use mod_stdio
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
    real(rkx) , pointer , dimension(:) :: sigmar => null( )
    real(rkx) :: pss , pst , p0
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
    integer(ik4) :: ni , nj
    real(rkx) , pointer , dimension(:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
  end type cmip6_2d_var

  type, extends(cmip6_file) :: cmip6_3d_var
    character(len=8) :: vname
    integer(ik4) :: ni , nj , nk
    real(rkx) , pointer , dimension(:,:,:) :: var
    type(cmip6_horizontal_coordinates) , pointer :: hcoord => null( )
    type(cmip6_vertical_coordinate) , pointer :: vcoord => null( )
  end type cmip6_3d_var

  public :: cmip6_2d_var , cmip6_3d_var
  public :: cmip6_fxpath , cmip6_path , cmip6_error

  contains

    character(len=1024) function cmip6_fxpath(ver,var) result(fpath)
      implicit none
      character(len=*) , intent(in) :: var , ver
      character(len=12) :: fx_variant , fx_experiment
      select case ( cmip6_model )
        case ( 'MPI-ESM1-2-HR' )
          fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep//'CMIP'//pthsep
          fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-HR'//pthsep
          fpath = trim(fpath)//'historical'//pthsep
          fx_variant = cmip6_variant
          fx_experiment = '_historical_'
        case ( 'HadGEM3-GC31-MM' )
          fpath = trim(cmip6_inp)//pthsep//'esg_cmip6'// &
            pthsep//'CMIP6'//pthsep//'HighResMIP'//pthsep
          fpath = trim(fpath)//'MOHC'//pthsep//'HadGEM3-GC31-MM'//pthsep
          fpath = trim(fpath)//'hist-1950'//pthsep
          fx_variant = 'r1i1p1f1'
          fx_experiment = '_hist-1950_'
        case ( 'NorESM2-MM' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep//'cmor'// &
            pthsep//'CMIP6'//pthsep//'CMIP/NCC'//pthsep//'NorESM2-MM'//pthsep
          fpath = trim(fpath)//'historical'//pthsep
          fx_variant = 'r1i1p1f1'
          fx_experiment = '_historical_'
        case ( 'CNRM-ESM2-1' )
          fpath = trim(cmip6_inp)//pthsep//'CMIP6_CNRM'//pthsep//'CMIP'// &
            pthsep//'CNRM-CERFACS'//pthsep//'CNRM-ESM2-1'//pthsep
          fpath = trim(fpath)//'amip'//pthsep
          fx_variant = 'r1i1p1f2'
          fx_experiment = '_amip_'
        case ( 'EC-Earth3' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot1'// &
            pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep//'CMIP'//&
            pthsep//'EC-Earth-Consortium'// &
            pthsep//'EC-Earth3'//pthsep//'historical'//pthsep
          fx_variant = 'r1i1p1f1'
          fx_experiment = '_historical_'
        case ( 'CESM2' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
            pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'NCAR'// &
            pthsep//'CESM2'//pthsep//'historical'//pthsep
          fx_variant = 'r11i1p1f1'
          fx_experiment = '_historical_'
        case ( 'GFDL-ESM4' )
          fpath = trim(cmip6_inp)//pthsep//'gfdl_dataroot4'// &
            pthsep//'AerChemMIP'//pthsep
          fpath = trim(fpath)//'NOAA-GFDL'//pthsep//'GFDL-ESM4'//pthsep
          fpath = trim(fpath)//'ssp370-lowNTCFCH4'//pthsep
          fx_variant = 'r1i1p1f1'
          fx_experiment = '_ssp370-lowNTCFCH4_'
        case default
          call die(__FILE__, &
            'Unsupported cmip6 model: '//trim(cmip6_model),-1)
      end select
      fpath = trim(fpath)//trim(fx_variant)// &
        pthsep//'fx'//pthsep//trim(var)//pthsep//trim(cmip6_grid)// &
        pthsep//trim(ver)//pthsep//trim(var)//'_fx_'//trim(cmip6_model)// &
        trim(fx_experiment)//trim(fx_variant)//'_'//trim(cmip6_grid)//'.nc'
    end function cmip6_fxpath

    character(len=1024) function cmip6_path(year,freq,ver,var) result(fpath)
      implicit none
      character(len=*) , intent(in) :: var , freq , ver
      integer(ik4) , intent(in) :: year
      character(len=12) :: experiment
      character(len=4) :: grid
      select case ( cmip6_model )
        case ( 'MPI-ESM1-2-HR' )
          fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-HR'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            fpath = trim(fpath)//'DKRZ'//pthsep//'MPI-ESM1-2-HR'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          grid = cmip6_grid
        case ( 'HadGEM3-GC31-MM' )
          fpath = trim(cmip6_inp)//pthsep//'esg_cmip6'//pthsep//'CMIP6'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'MOHC'//pthsep//'HadGEM3-GC31-MM'//pthsep
          grid = cmip6_grid
        case ( 'NorESM2-MM' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep//'cmor'// &
            pthsep//'CMIP6'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'NCC'//pthsep//'NorESM2-MM'//pthsep
          grid = cmip6_grid
        case ( 'CNRM-ESM2-1' )
          fpath = trim(cmip6_inp)//pthsep//'CMIP6_CNRM'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'CNRM-CERFACS'//pthsep//'CNRM-ESM2-1'//pthsep
          if ( var == 'tos' ) then
            grid = 'gn'
          else
            grid = cmip6_grid
          end if
        case ( 'EC-Earth3' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot1'// &
            pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'EC-Earth-Consortium'//pthsep// &
            'EC-Earth3'//pthsep
          if ( var == 'tos' ) then
            grid = 'gn'
          else
            grid = cmip6_grid
          end if
        case ( 'CESM2' )
          fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep// &
            'CMIP6'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'historical'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'NCAR'//pthsep//'CESM2'//pthsep
          grid = cmip6_grid
        case ( 'GFDL-ESM4' )
          fpath = trim(cmip6_inp)//pthsep//'gfdl_dataroot4'//pthsep
          if ( year < 2015 ) then
            fpath = trim(fpath)//'CMIP'//pthsep
            experiment = 'esm-hist'
          else
            fpath = trim(fpath)//'ScenarioMIP'//pthsep
            experiment = trim(cmip6_ssp)
          end if
          fpath = trim(fpath)//'NOAA-GFDL'//pthsep//'GFDL-ESM4'//pthsep
          if ( var == 'tos' ) then
            grid = 'gn'
          else
            grid = cmip6_grid
          end if
        case default
          call die(__FILE__, &
            'Unsupported cmip6 model: '//trim(cmip6_model),-1)
      end select
      fpath = trim(fpath)//trim(experiment)//pthsep
      fpath = trim(fpath)//trim(cmip6_variant)//pthsep//trim(freq)//pthsep// &
        trim(var)//pthsep//trim(grid)//pthsep//trim(ver)// &
        pthsep//trim(var)//'_'//trim(freq)//'_'//trim(cmip6_model)//'_'
      fpath = trim(fpath)//trim(experiment)//'_'
      fpath = trim(fpath)//trim(cmip6_variant)//'_'//trim(grid)//'_'
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

end module mod_cmip6_helper

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
