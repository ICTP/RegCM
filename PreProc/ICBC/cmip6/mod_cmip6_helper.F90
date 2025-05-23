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
    real(rkx), pointer, contiguous, dimension(:) :: lon1d => null( )
    real(rkx), pointer, contiguous, dimension(:) :: lat1d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lon2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lat2d => null( )
  end type cmip6_horizontal_coordinates

  type cmip6_vertical_coordinate
    real(rkx), pointer, contiguous, dimension(:) :: plev => null( )
    real(rkx), pointer, contiguous, dimension(:) :: sigmar => null( )
    real(rkx) :: pss, pst, p0
    real(rkx), pointer, contiguous, dimension(:) :: ak => null( )
    real(rkx), pointer, contiguous, dimension(:) :: bk => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: topo => null( )
  end type cmip6_vertical_coordinate

  type cmip6_file
    character(len=1024) :: filename
    integer(ik4) :: ncid = -1
    integer(ik4) :: ivar = -1
    integer(ik4) :: nrec = -1
    type(rcm_time_and_date) :: first_date
    type(h_interpolator), pointer, contiguous, dimension(:) :: hint => null( )
  end type cmip6_file

  type, extends(cmip6_file) :: cmip6_2d_var
    character(len=8) :: vname
    integer(ik4) :: ni, nj
    real(rkx), pointer, contiguous, dimension(:,:) :: var => null( )
    type(cmip6_horizontal_coordinates), pointer :: hcoord => null( )
  end type cmip6_2d_var

  type, extends(cmip6_file) :: cmip6_3d_var
    character(len=8) :: vname
    integer(ik4) :: ni, nj, nk
    real(rkx), pointer, contiguous, dimension(:,:,:) :: var => null( )
    type(cmip6_horizontal_coordinates), pointer :: hcoord => null( )
    type(cmip6_vertical_coordinate), pointer :: vcoord => null( )
  end type cmip6_3d_var

  public :: cmip6_2d_var, cmip6_3d_var
  public :: cmip6_fxpath, cmip6_path
  public :: cmip6_error

  contains

    character(len=1024) function cmip6_fxpath(ver,var) result(fpath)
      implicit none
      character(len=*), intent(in) :: var, ver
      character(len=24) :: fx_variant, fx_experiment, fx_model
      if ( dattyp == 'CMIP6' ) then
        select case ( cmip6_model )
          case ( 'MPI-ESM1-2-HR' )
            fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep//'CMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-HR'//pthsep
            fpath = trim(fpath)//'historical'//pthsep
            fx_variant = cmip6_variant
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'HadGEM3-GC31-MM' )
            if ( cmip6_inp(1:8) == 'https://' ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_cmip6'// &
                pthsep//'CMIP6'//pthsep//'HighResMIP'//pthsep
            else
              fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep// &
                'HighResMIP'//pthsep
            end if
            fpath = trim(fpath)//'MOHC'//pthsep//'HadGEM3-GC31-MM'//pthsep
            fpath = trim(fpath)//'hist-1950'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_hist-1950_'
            fx_model = cmip6_model
          case ( 'NorESM2-MM' )
            if ( index(cmip6_inp,'noresg.nird.sigma2.no') > 0 ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep//'cmor'//&
                pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'NCC'//pthsep// &
                'NorESM2-MM'//pthsep
            else
              ! This should work for esgf3
              fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep// &
                'CMIP'//pthsep//'NCC'//pthsep//'NorESM2-MM'//pthsep
            end if
            fpath = trim(fpath)//'historical'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'CNRM-ESM2-1' )
            if ( index(cmip6_inp,'umr-cnrm.fr') > 0 ) then
              fpath = trim(cmip6_inp)//pthsep//'CMIP6_CNRM'//pthsep//'CMIP'// &
                pthsep//'CNRM-CERFACS'//pthsep//'CNRM-ESM2-1'//pthsep
              fpath = trim(fpath)//'amip'//pthsep
            else
              fpath = trim(cmip6_inp)//pthsep//'CMIP6'//pthsep//'CMIP'// &
                pthsep//'CNRM-CERFACS'//pthsep//'CNRM-ESM2-1'//pthsep
              fpath = trim(fpath)//'amip'//pthsep
            end if
            fx_variant = 'r1i1p1f2'
            fx_experiment = '_amip_'
            fx_model = cmip6_model
          case ( 'EC-Earth3-Veg' )
            if ( cmip6_inp(1:8) == 'https://' ) then
              fpath = 'https://esg-dn2.nsc.liu.se/thredds/dodsC'//pthsep// &
                'esg_dataroot1'//pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep//&
                'CMIP'//pthsep//'EC-Earth-Consortium'// &
                pthsep//'EC-Earth3'//pthsep//'historical'//pthsep
            else
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot1'//pthsep// &
                'cmip6data'//pthsep//'CMIP6'//pthsep// &
                'CMIP'//pthsep//'EC-Earth-Consortium'// &
                pthsep//'EC-Earth3'//pthsep//'historical'//pthsep
            end if
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_historical_'
            fx_model = 'EC-Earth3'
          case ( 'CESM2' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
              pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'NCAR'// &
              pthsep//'CESM2'//pthsep//'historical'//pthsep
            fx_variant = 'r11i1p1f1'
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'CMCC-ESM2' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
              pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'CMCC'// &
              pthsep//'CMCC-ESM2'//pthsep//'historical'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'GFDL-ESM4' )
            fpath = trim(cmip6_inp)//pthsep//'gfdl_dataroot4'// &
              pthsep//'AerChemMIP'//pthsep
            fpath = trim(fpath)//'NOAA-GFDL'//pthsep//'GFDL-ESM4'//pthsep
            fpath = trim(fpath)//'ssp370-lowNTCFCH4'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_ssp370-lowNTCFCH4_'
            fx_model = cmip6_model
          case ( 'CanESM5' )
            fpath = trim(cmip6_inp)//pthsep//'esgA_dataroot'// &
              pthsep//'AR6'//pthsep//'CMIP6'//pthsep//'CMIP'// &
              pthsep//'CCCma'//pthsep//'CanESM5'//pthsep//'historical'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'MIROC6' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
              pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'MIROC'// &
              pthsep//'MIROC6'//pthsep//'historical'//pthsep
            fx_variant = 'r1i1p1f1'
            fx_experiment = '_historical_'
            fx_model = cmip6_model
          case ( 'MIROC-ES2L')
            fx_variant = cmip6_variant
            fx_model = cmip6_model
            if ( cmip6_experiment(1:3) == 'ssp' ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
                pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'MIROC'// &
                pthsep//'MIROC-ES2L'//pthsep//'historical'//pthsep
              fx_experiment = '_historical_'
            else
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'// &
                pthsep//'CMIP6'//pthsep//'CMIP'//pthsep//'MIROC'// &
                pthsep//'MIROC-ES2L'//pthsep//trim(cmip6_experiment)//pthsep
              fx_experiment = '_'//trim(cmip6_experiment)//'_'
            end if
          case ( 'MPI-ESM1-2-LR' )
            fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep//'CMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-LR'//pthsep
            fpath = trim(fpath)//trim(cmip6_experiment)//pthsep
            fx_variant = cmip6_variant
            fx_experiment = '_'//trim(cmip6_experiment)//'_'
            fx_model = cmip6_model
          case default
            call die(__FILE__, &
              'Unsupported cmip6 model: '//trim(cmip6_model),-1)
        end select
        fpath = trim(fpath)//trim(fx_variant)// &
          pthsep//'fx'//pthsep//trim(var)//pthsep//trim(cmip6_grid)// &
          pthsep//trim(ver)//pthsep//trim(var)//'_fx_'//trim(fx_model)// &
          trim(fx_experiment)//trim(fx_variant)//'_'//trim(cmip6_grid)//'.nc'
      else
        select case ( pmip4_model )
          case ( 'MPI-ESM1-2-LR' )
            fpath = trim(pmip4_inp)//pthsep//'cmip6'//pthsep//'PMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-LR'//pthsep
            fpath = trim(fpath)//trim(pmip4_experiment)//pthsep
            fx_variant = pmip4_variant
            fx_experiment = '_'//trim(pmip4_experiment)//'_'
            fx_model = pmip4_model
          case ( 'IPSL-CM6A-LR' )
            fpath = trim(pmip4_inp)//pthsep//'cmip6'//pthsep//'PMIP'//pthsep
            fpath = trim(fpath)//'IPSL'//pthsep//'IPSL-CM6A-LR'//pthsep
            fpath = trim(fpath)//trim(pmip4_experiment)//pthsep
            fx_variant = pmip4_variant
            fx_experiment = '_'//trim(pmip4_experiment)//'_'
            fx_model = pmip4_model
          case default
            call die(__FILE__, &
              'Unsupported pmip4 model: '//trim(pmip4_model),-1)
        end select
        fpath = trim(fpath)//trim(fx_variant)// &
          pthsep//'fx'//pthsep//trim(var)//pthsep//trim(pmip4_grid)// &
          pthsep//trim(ver)//pthsep//trim(var)//'_fx_'//trim(fx_model)// &
          trim(fx_experiment)//trim(fx_variant)//'_'//trim(pmip4_grid)//'.nc'
      end if
    end function cmip6_fxpath

    character(len=1024) function cmip6_path(year,freq,ver,var) result(fpath)
      implicit none
      integer(ik4), intent(in) :: year
      character(len=*), intent(in) :: var, freq, ver
      character(len=12) :: experiment
      character(len=12) :: grid

      if ( dattyp == 'CMIP6' ) then
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
              experiment = trim(cmip6_experiment)
            end if
            grid = cmip6_grid
          case ( 'HadGEM3-GC31-MM' )
            if ( cmip6_inp(1:8) == 'https://' ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_cmip6'// &
                pthsep//'CMIP6'//pthsep
            else
              fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep
            end if
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'MOHC'//pthsep//'HadGEM3-GC31-MM'//pthsep
            grid = cmip6_grid
          case ( 'NorESM2-MM' )
            if ( index(cmip6_inp,'noresg.nird.sigma2.no') > 0 ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep// &
                'cmor'//pthsep//'CMIP6'//pthsep
            else
              ! This should work for esgf3
              fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep
            end if
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'NCC'//pthsep//'NorESM2-MM'//pthsep
            grid = cmip6_grid
          case ( 'CNRM-ESM2-1' )
            if ( index(cmip6_inp,'umr-cnrm.fr') > 0 ) then
              fpath = trim(cmip6_inp)//pthsep//'CMIP6_CNRM'//pthsep
            else
              fpath = trim(cmip6_inp)//pthsep//'CMIP6'//pthsep
            end if
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'CNRM-CERFACS'//pthsep//'CNRM-ESM2-1'//pthsep
            if ( var == 'tos' ) then
              grid = 'gn'
            else
              grid = cmip6_grid
            end if
          case ( 'EC-Earth3-Veg' )
            if ( var == 'tos' ) then
              fpath = trim(cmip6_inp)//pthsep//'esg_dataroot6'// &
                pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep
            else if ( var == 'ps' ) then
              if ( year < 2015 ) then
                fpath = trim(cmip6_inp)//pthsep//'esg_dataroot6'// &
                  pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep
              else
                if ( cmip6_inp(1:8) == 'https://' ) then
                  fpath = &
                    'https://esg-dn3.nsc.liu.se/thredds/dodsC/esg_dataroot2'// &
                    pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep
                else
                  fpath = trim(cmip6_inp)//pthsep//'cmip6data'// &
                      pthsep//'CMIP6'//pthsep
                end if
              end if
            else
              if ( cmip6_inp(1:8) == 'https://' ) then
                fpath = &
                  'https://esg-dn3.nsc.liu.se/thredds/dodsC/esg_dataroot2'// &
                  pthsep//'cmip6data'//pthsep//'CMIP6'//pthsep
              else
                fpath = trim(cmip6_inp)//pthsep//'cmip6data'// &
                    pthsep//'CMIP6'//pthsep
              end if
            end if
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'EC-Earth-Consortium'//pthsep// &
              'EC-Earth3-Veg'//pthsep
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
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'NCAR'//pthsep//'CESM2'//pthsep
            grid = cmip6_grid
          case ( 'CMCC-ESM2' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep// &
              'CMIP6'//pthsep
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'CMCC'//pthsep//'CMCC-ESM2'//pthsep
            grid = cmip6_grid
          case ( 'GFDL-ESM4' )
            fpath = trim(cmip6_inp)//pthsep//'gfdl_dataroot4'//pthsep
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'esm-hist'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'NOAA-GFDL'//pthsep//'GFDL-ESM4'//pthsep
            if ( var == 'tos' ) then
              grid = 'gn'
            else
              grid = cmip6_grid
            end if
          case ( 'CanESM5' )
            if ( year < 2015 ) then
              if ( var == 'tos' ) then
                fpath = trim(cmip6_inp)//pthsep//'esgC_dataroot'//pthsep// &
                  'AR6'//pthsep//'CMIP6'//pthsep
              else
                fpath = trim(cmip6_inp)//pthsep//'esgA_dataroot'//pthsep// &
                  'AR6'//pthsep//'CMIP6'//pthsep
              end if
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              if ( var == 'tos' ) then
                fpath = trim(cmip6_inp)//pthsep//'esgD_dataroot'//pthsep// &
                  'AR6'//pthsep//'CMIP6'//pthsep
              else
                fpath = trim(cmip6_inp)//pthsep//'esgF_dataroot'//pthsep// &
                  'AR6'//pthsep//'CMIP6'//pthsep
              end if
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'CCCma'//pthsep//'CanESM5'//pthsep
            grid = cmip6_grid
          case ( 'MIROC6' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep// &
              'CMIP6'//pthsep
            if ( year < 2015 ) then
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = 'historical'
            else
              fpath = trim(fpath)//'ScenarioMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'MIROC'//pthsep//'MIROC6'//pthsep
            grid = cmip6_grid
          case ( 'MIROC-ES2L' )
            fpath = trim(cmip6_inp)//pthsep//'esg_dataroot'//pthsep// &
              'CMIP6'//pthsep
            if ( cmip6_experiment(1:3) == 'ssp' ) then
              if ( year < 2015 ) then
                fpath = trim(fpath)//'CMIP'//pthsep
                experiment = 'historical'
              else
                fpath = trim(fpath)//'ScenarioMIP'//pthsep
                experiment = trim(cmip6_experiment)
              end if
            else
              fpath = trim(fpath)//'CMIP'//pthsep
              experiment = trim(cmip6_experiment)
            end if
            fpath = trim(fpath)//'MIROC'//pthsep//'MIROC-ES2L'//pthsep
            grid = cmip6_grid
          case ( 'MPI-ESM1-2-LR' )
            fpath = trim(cmip6_inp)//pthsep//'cmip6'//pthsep
            fpath = trim(fpath)//'CMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-LR'//pthsep
            experiment = cmip6_experiment
            grid = cmip6_grid
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
      else
        select case ( pmip4_model )
          case ( 'MPI-ESM1-2-LR' )
            fpath = trim(pmip4_inp)//pthsep//'cmip6'//pthsep
            fpath = trim(fpath)//'PMIP'//pthsep
            fpath = trim(fpath)//'MPI-M'//pthsep//'MPI-ESM1-2-LR'//pthsep
            experiment = pmip4_experiment
            grid = pmip4_grid
          case ( 'IPSL-CM6A-LR' )
            fpath = trim(pmip4_inp)//pthsep//'cmip6'//pthsep
            fpath = trim(fpath)//'PMIP'//pthsep
            fpath = trim(fpath)//'IPSL'//pthsep//'IPSL-CM6A-LR'//pthsep
            experiment = pmip4_experiment
            if ( var == 'tos' ) then
              grid = 'gn'
            else
              grid = pmip4_grid
            end if
          case default
            call die(__FILE__, &
              'Unsupported pmip4 model: '//trim(pmip4_model),-1)
        end select
        fpath = trim(fpath)//trim(experiment)//pthsep
        fpath = trim(fpath)//trim(pmip4_variant)//pthsep//trim(freq)//pthsep// &
          trim(var)//pthsep//trim(grid)//pthsep//trim(ver)// &
          pthsep//trim(var)//'_'//trim(freq)//'_'//trim(pmip4_model)//'_'
        fpath = trim(fpath)//trim(experiment)//'_'
        fpath = trim(fpath)//trim(pmip4_variant)//'_'//trim(grid)//'_'
      end if
    end function cmip6_path

    subroutine cmip6_error(ival,filename,line,arg)
      implicit none
      integer(ik4), intent(in) :: ival, line
      character(len=8) :: cline
      character(*), intent(in) :: filename, arg
      if ( ival /= nf90_noerr ) then
        write (cline,'(i8)') line
        write (stderr,*) nf90_strerror(ival)
        call die(filename,trim(cline)//':'//arg,ival)
      end if
    end subroutine cmip6_error

end module mod_cmip6_helper

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
