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

module mod_cmip7_helper

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

  type cmip7_horizontal_coordinates
    real(rkx), pointer, contiguous, dimension(:) :: lon1d => null( )
    real(rkx), pointer, contiguous, dimension(:) :: lat1d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lon2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lat2d => null( )
  end type cmip7_horizontal_coordinates

  type cmip7_vertical_coordinate
    real(rkx), pointer, contiguous, dimension(:) :: plev => null( )
    real(rkx), pointer, contiguous, dimension(:) :: sigmar => null( )
    real(rkx) :: pss, pst, p0
    real(rkx), pointer, contiguous, dimension(:) :: ak => null( )
    real(rkx), pointer, contiguous, dimension(:) :: bk => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: topo => null( )
  end type cmip7_vertical_coordinate

  type cmip7_file
    character(len=1024) :: filename
    integer(ik4) :: ncid = -1
    integer(ik4) :: ivar = -1
    integer(ik4) :: nrec = -1
    type(rcm_time_and_date) :: first_date
    type(h_interpolator), pointer, contiguous, dimension(:) :: hint => null( )
  end type cmip7_file

  type, extends(cmip7_file) :: cmip7_2d_var
    character(len=8) :: vname
    integer(ik4) :: ni, nj
    real(rkx), pointer, contiguous, dimension(:,:) :: var => null( )
    type(cmip7_horizontal_coordinates), pointer :: hcoord => null( )
  end type cmip7_2d_var

  type, extends(cmip7_file) :: cmip7_3d_var
    character(len=8) :: vname
    integer(ik4) :: ni, nj, nk
    real(rkx), pointer, contiguous, dimension(:,:,:) :: var => null( )
    type(cmip7_horizontal_coordinates), pointer :: hcoord => null( )
    type(cmip7_vertical_coordinate), pointer :: vcoord => null( )
  end type cmip7_3d_var

  public :: cmip7_2d_var, cmip7_3d_var
  public :: cmip7_fxpath, cmip7_path
  public :: cmip7_error

  contains

    character(len=1024) function cmip7_fxpath(ver,var) result(fpath)
      implicit none
      character(len=*), intent(in) :: var, ver
      character(len=24) :: fx_variant, fx_experiment, fx_model, fx_label
      character(len=24) :: fx_freq, fx_grid_label
      if ( dattyp == 'CMIP7' ) then
        select case ( cmip7_model )
          case ( 'EC-Earth3-ESM-1-1' )
            fpath = trim(cmip7_inp)//pthsep//'cmip7'//pthsep//'CMIP'//pthsep
            fpath = trim(fpath)//'EC-Earth-Consortium'//pthsep// &
                    cmip7_model//pthsep//'esm-hist'//pthsep
            fx_variant = cmip7_variant
            fx_freq = 'fx_'
            fx_label = '_ti-u-hxy-u_'
            fx_experiment = '_esm-hist_'
            fx_model = cmip7_model
            fx_grid_label = cmip7_atmo_grid_label
          case default
            call die(__FILE__, &
              'Unsupported cmip7 model: '//trim(cmip7_model),-1)
        end select
        fpath = trim(fpath)//trim(fx_variant)//pthsep//trim(cmip7_region)// &
          pthsep//'fx'//pthsep//trim(var)//pthsep//trim(fx_grid_label)// &
          pthsep//trim(fx_label)//pthsep//trim(fx_grid_label)//pthsep// &
          trim(ver)//pthsep//trim(var)//trim(fx_label)//trim(fx_freq)// &
          trim(fx_label)//trim(cmip7_region)//'_'//trim(fx_grid_label)// &
          '_'//trim(fx_model)//trim(fx_experiment)//trim(fx_variant)//'.nc'
      end if
    end function cmip7_fxpath

    character(len=1024) function cmip7_path(year,freq,ver,var) result(fpath)
      implicit none
      integer(ik4), intent(in) :: year
      character(len=*), intent(in) :: var, freq, ver
      character(len=12) :: experiment, grid
      character(len=16) :: vlabel

      if ( dattyp == 'CMIP7' ) then
        select case ( cmip7_model )
          case ( 'EC-Earth3-ESM-1-1' )
            fpath = trim(cmip7_inp)//pthsep//'cmip7'//pthsep
            if ( var == 'tos' ) then
              vlabel = 'tavg-u-hxy-sea'
              grid = cmip7_ocn_grid_label
            else
              vlabel = 'tpt-al-hxy-u'
              grid = cmip7_atmo_grid_label
            end if
            if ( year < 2015 ) then
              fpath = trim(fpath)//pthsep//'CMIP'//pthsep
              fpath = trim(fpath)//'EC-Earth-Consortium'//pthsep// &
                      trim(cmip7_model)//pthsep
              experiment = 'esm-hist'
            else
              fpath = trim(fpath)//pthsep//'ScenarioMIP'//pthsep
              fpath = trim(fpath)//'EC-Earth-Consortium'//pthsep// &
                      trim(cmip7_model)//pthsep
              experiment = trim(cmip7_experiment)
            end if
          case default
            call die(__FILE__, &
              'Unsupported cmip7 model: '//trim(cmip7_model),-1)
        end select
        fpath = trim(fpath)//trim(experiment)//pthsep
        fpath = trim(fpath)//trim(cmip7_variant)//pthsep// &
          trim(cmip7_region)//pthsep//trim(freq)//pthsep// &
          trim(var)//pthsep//trim(vlabel)//pthsep// &
          trim(grid)//pthsep//trim(ver)//pthsep// &
          trim(var)//'_'//trim(vlabel)//'_'//trim(freq)//'_'// &
          trim(cmip7_region)//'_'//trim(grid)//'_'// &
          trim(cmip7_model)//'_'//trim(experiment)//'_'// &
          trim(cmip7_variant)//'_'
      end if
    end function cmip7_path

    subroutine cmip7_error(ival,filename,line,arg)
      implicit none
      integer(ik4), intent(in) :: ival, line
      character(len=8) :: cline
      character(*), intent(in) :: filename, arg
      if ( ival /= nf90_noerr ) then
        write (cline,'(i8)') line
        write (stderr,*) nf90_strerror(ival)
        call die(filename,trim(cline)//':'//arg,ival)
      end if
    end subroutine cmip7_error

end module mod_cmip7_helper

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
