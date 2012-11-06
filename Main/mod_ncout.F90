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
!
module mod_ncout
!
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_ncstream_types
  use mod_ncstream
!
  integer(ik4) , parameter :: nbasevars = 5
  type(ncvariable2d_real) , target , dimension(nbasevars) :: v2dvar_base

  type regcm_stream
    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    integer(ik4) :: nvar
    type(nc_varlist) :: ncvars
  end type regcm_stream

  integer(ik4) , parameter :: natm2dvars = 2
  integer(ik4) , parameter :: natm3dvars = 10
  integer(ik4) , public , parameter :: natmvars = natm2dvars+natm3dvars
  logical , public , dimension(natmvars) :: enable_atm_vars
  type(ncvariable2d_real) , target , dimension(natm2dvars) :: v2dvar_atm
  type(ncvariable3d_real) , target , dimension(natm3dvars) :: v3dvar_atm

  integer(ik4) , parameter :: nsrf2dvars = 18
  integer(ik4) , parameter :: nsrf3dvars = 5
  integer(ik4) , public , parameter :: nsrfvars = nsrf2dvars+nsrf3dvars
  logical , public , dimension(nsrfvars) :: enable_srf_vars
  type(ncvariable2d_real) , target , dimension(nsrf2dvars) :: v2dvar_srf
  type(ncvariable3d_real) , target , dimension(nsrf3dvars) :: v3dvar_srf

  integer(ik4) , parameter :: nsts2dvars = 6
  integer(ik4) , parameter :: nsts3dvars = 4
  integer(ik4) , public , parameter :: nstsvars = nsts2dvars+nsts3dvars
  logical , public , dimension(nstsvars) :: enable_sts_vars
  type(ncvariable2d_real) , target , dimension(nsts2dvars) :: v2dvar_sts
  type(ncvariable2d_real) , target , dimension(nsts3dvars) :: v3dvar_sts

  integer(ik4) , parameter :: nsub2dvars = 10
  integer(ik4) , parameter :: nsub3dvars = 4
  integer(ik4) , public , parameter :: nsubvars = nsub2dvars+nsub3dvars
  logical , public , dimension(nsubvars) :: enable_sub_vars
  type(ncvariable2d_real) , target , dimension(nsub2dvars) :: v2dvar_sub
  type(ncvariable2d_real) , target , dimension(nsub3dvars) :: v3dvar_sub

  integer(ik4) , parameter :: nlak2dvars = 13
  integer(ik4) , parameter :: nlak3dvars = 1
  integer(ik4) , public , parameter :: nlakvars = nlak2dvars+nlak3dvars
  logical , public , dimension(nlakvars) :: enable_lak_vars
  type(ncvariable2d_real) , target , dimension(nlak2dvars) :: v2dvar_lak
  type(ncvariable2d_real) , target , dimension(nlak3dvars) :: v3dvar_lak

  integer(ik4) , parameter :: nrad2dvars = 12
  integer(ik4) , parameter :: nrad3dvars = 4
  integer(ik4) , public , parameter :: nradvars = nrad2dvars+nrad3dvars
  logical , public , dimension(nradvars) :: enable_rad_vars
  type(ncvariable2d_real) , target , dimension(nrad2dvars) :: v2dvar_rad
  type(ncvariable2d_real) , target , dimension(nrad3dvars) :: v3dvar_rad

  integer(ik4) , parameter :: maxstreams = 100

  integer(ik4) , public , parameter :: atm_stream = 1
  integer(ik4) , public , parameter :: srf_stream = 2
  integer(ik4) , public , parameter :: sub_stream = 3
  integer(ik4) , public , parameter :: rad_stream = 4
  integer(ik4) , public , parameter :: lak_stream = 5
  integer(ik4) , public , parameter :: sts_stream = 6
!
  type(regcm_stream) , dimension(maxstreams) :: outstream
  logical , public , dimension(maxstreams) :: enable_flag

  public :: init_output_streams

contains

  subroutine init_output_streams(lparallel,enable_flag,stream_number)
    implicit none
    logical , intent(in) :: lparallel
    logical , intent(in) , dimension(:) :: enable_flag
    integer(ik4) , intent(in) :: stream_number
    integer :: nstream , i , vcount
    integer :: imax1 , imax2 , imax3

    if ( .not. lparallel ) then
      if ( myid /= iocpu ) then
        return
      end if
    end if

    v2dvar_base(1)%vname = 'xlon'
    v2dvar_base(1)%vunit = 'degrees_east'
    v2dvar_base(1)%long_name = 'Longitude on Cross Points'
    v2dvar_base(1)%standard_name = 'longitude'
    v2dvar_base(2)%vname = 'xlat'
    v2dvar_base(2)%vunit = 'degrees_north'
    v2dvar_base(2)%long_name = 'Latitude on Cross Points'
    v2dvar_base(2)%standard_name = 'latitude'
    v2dvar_base(3)%vname = 'mask'
    v2dvar_base(3)%vunit = '1'
    v2dvar_base(3)%long_name = 'Land Mask'
    v2dvar_base(3)%standard_name = 'land_binary_mask'
    v2dvar_base(4)%vname = 'topo'
    v2dvar_base(4)%vunit = 'm'
    v2dvar_base(4)%long_name = 'Surface Model Elevation'
    v2dvar_base(4)%standard_name = 'surface_altitude'
    v2dvar_base(5)%vname = 'ps'
    v2dvar_base(5)%vunit = 'hPa'
    v2dvar_base(5)%long_name = 'Surface Pressure'
    v2dvar_base(5)%standard_name = 'surface_air_pressure'
    v2dvar_base(5)%lrecords = .true.

    do nstream = 1 , stream_number
      if ( .not. enable_flag(nstream) ) return

      if ( nstream == atm_stream ) then

        ! The enable flag respect the below order

        v2dvar_atm(1)%vname = 'tpr'
        v2dvar_atm(1)%vunit = 'kg m-2 day-1'
        v2dvar_atm(1)%long_name = 'Total rain precipitation flux'
        v2dvar_atm(1)%standard_name = 'precipitation_flux'
        v2dvar_atm(1)%lrecords = .true.
        v2dvar_atm(1)%cell_method = 'time: mean'
        v2dvar_atm(2)%vname = 'tgb'
        v2dvar_atm(2)%vunit = 'K'
        v2dvar_atm(2)%long_name = 'Lower groud temperature'
        v2dvar_atm(2)%standard_name = 'soil_temperature'
        v2dvar_atm(2)%lrecords = .true.

        v3dvar_atm(1)%vname = 'u'
        v3dvar_atm(1)%vunit = 'm s-1'
        v3dvar_atm(1)%long_name = 'Zonal component of wind (westerly)'
        v3dvar_atm(1)%standard_name = 'eastward_wind'
        v3dvar_atm(1)%lrecords = .true.
        v3dvar_atm(2)%vname = 'v'
        v3dvar_atm(2)%vunit = 'm s-1'
        v3dvar_atm(2)%long_name = 'Meridional component of wind (southerly)'
        v3dvar_atm(2)%standard_name = 'northward_wind'
        v3dvar_atm(2)%lrecords = .true.
        v3dvar_atm(3)%vname = 't'
        v3dvar_atm(3)%vunit = 'K'
        v3dvar_atm(3)%long_name = 'Air Temperature'
        v3dvar_atm(3)%standard_name = 'air_temperature'
        v3dvar_atm(3)%lrecords = .true.
        v3dvar_atm(4)%vname = 'omega'
        v3dvar_atm(4)%vunit = 'hPa s-1'
        v3dvar_atm(4)%long_name = 'Pressure velocity'
        v3dvar_atm(4)%standard_name = 'lagrangian_tendency_of_air_pressure'
        v3dvar_atm(4)%lrecords = .true.
        v3dvar_atm(5)%vname = 'qv'
        v3dvar_atm(5)%vunit = 'kg kg-1'
        v3dvar_atm(5)%long_name = 'Water vapor mixing ratio'
        v3dvar_atm(5)%standard_name = 'humidity_mixing_ratio'
        v3dvar_atm(5)%lrecords = .true.
        v3dvar_atm(6)%vname = 'qc'
        v3dvar_atm(6)%vunit = 'kg kg-1'
        v3dvar_atm(6)%long_name = 'Cloud liquid water mixing ratio'
        v3dvar_atm(6)%standard_name = 'cloud_liquid_water_mixing_ratio'
        v3dvar_atm(6)%lrecords = .true.
        v3dvar_atm(7)%vname = 'tke'
        v3dvar_atm(7)%vunit = 'm2 s2'
        v3dvar_atm(7)%long_name = 'Turbulent Kinetic Energy'
        v3dvar_atm(7)%standard_name = 'turbulent_kinetic_energy'
        v3dvar_atm(7)%lrecords = .true.
        v3dvar_atm(8)%vname = 'kth'
        v3dvar_atm(8)%vunit = 'm2 s-1'
        v3dvar_atm(8)%long_name = 'Vertical Turbulent Viscosity'
        v3dvar_atm(8)%standard_name = 'vertical_momentum_diffusivity'
        v3dvar_atm(8)%lrecords = .true.
        v3dvar_atm(9)%vname = 'kzm'
        v3dvar_atm(9)%vunit = 'm2 s-1'
        v3dvar_atm(9)%long_name = 'Vertical Turbulent Diffusivity'
        v3dvar_atm(9)%standard_name = 'vertical_scalar_diffusivity'
        v3dvar_atm(9)%lrecords = .true.
        v3dvar_atm(10)%vname = 'swt'
        v3dvar_atm(10)%vunit = 'kg m-2'
        v3dvar_atm(10)%long_name = 'Total soil water'
        v3dvar_atm(10)%standard_name = 'moisture_content_of_soil_layer'
        v3dvar_atm(10)%lrecords = .true.
        v3dvar_atm(10)%axis = 'xys'

        outstream(atm_stream)%nvar = countvars(enable_atm_vars,natmvars)
        allocate(outstream(atm_stream)%ncvars%vlist(outstream(atm_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(atm_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , natm2dvars
          if ( enable_atm_vars(i) ) then
            outstream(atm_stream)%ncvars%vlist(vcount)%vp => v2dvar_atm(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , natm3dvars
          if ( enable_atm_vars(i+natm2dvars) ) then
            outstream(atm_stream)%ncvars%vlist(vcount)%vp => v3dvar_atm(i)
            vcount = vcount + 1
          end if
        end do
      end if

      if ( nstream == srf_stream ) then

        ! The enable flag respect the below order

        v2dvar_srf(1)%vname = 'uvdrag'
        v2dvar_srf(1)%vunit = '1'
        v2dvar_srf(1)%long_name = 'Surface drag stress coefficient in air'
        v2dvar_srf(1)%standard_name = 'surface_drag_coefficient_in_air'
        v2dvar_srf(1)%lrecords = .true.
        v2dvar_srf(2)%vname = 'tg'
        v2dvar_srf(2)%vunit = 'K'
        v2dvar_srf(2)%long_name = 'Ground surface temperature'
        v2dvar_srf(2)%standard_name = 'surface_temperature'
        v2dvar_srf(2)%lrecords = .true.
        v2dvar_srf(3)%vname = 'tlef'
        v2dvar_srf(3)%vunit = 'K'
        v2dvar_srf(3)%long_name = 'Foliage canopy temperature'
        v2dvar_srf(3)%standard_name = 'canopy_temperature'
        v2dvar_srf(3)%lrecords = .true.
        v2dvar_srf(4)%vname = 'tpr'
        v2dvar_srf(4)%vunit = 'kg m-2 day-1'
        v2dvar_srf(4)%long_name = 'Total precipitation flux'
        v2dvar_srf(4)%standard_name = 'precipitation_flux'
        v2dvar_srf(4)%lrecords = .true.
        v2dvar_srf(4)%cell_method = 'time: mean'
        v2dvar_srf(5)%vname = 'evp'
        v2dvar_srf(5)%vunit = 'kg m-2 day-1'
        v2dvar_srf(5)%long_name = 'Total evapotranspiration flux'
        v2dvar_srf(5)%standard_name = 'water_evaporation_flux'
        v2dvar_srf(5)%lrecords = .true.
        v2dvar_srf(5)%cell_method = 'time: mean'
        v2dvar_srf(6)%vname = 'runoff'
        v2dvar_srf(6)%vunit = 'kg m-2 day-1'
        v2dvar_srf(6)%long_name = 'Surface runoff flux'
        v2dvar_srf(6)%standard_name = 'surface_runoff_flux'
        v2dvar_srf(6)%lrecords = .true.
        v2dvar_srf(6)%cell_method = 'time: mean'
        v2dvar_srf(7)%vname = 'scv'
        v2dvar_srf(7)%vunit = 'kg m-2'
        v2dvar_srf(7)%long_name = 'Liquid water equivalent of snow thickness'
        v2dvar_srf(7)%standard_name = 'lwe_thickness_of_surface_snow_amount'
        v2dvar_srf(7)%lrecords = .true.
        v2dvar_srf(7)%cell_method = 'time: mean'
        v2dvar_srf(8)%vname = 'sena'
        v2dvar_srf(8)%vunit = 'W m-2'
        v2dvar_srf(8)%long_name = 'Sensible heat flux'
        v2dvar_srf(8)%standard_name = 'surface_downward_sensible_heat_flux'
        v2dvar_srf(8)%lrecords = .true.
        v2dvar_srf(8)%cell_method = 'time: mean'
        v2dvar_srf(9)%vname = 'flw'
        v2dvar_srf(9)%vunit = 'W m-2'
        v2dvar_srf(9)%long_name = 'Net upward longwave energy flux'
        v2dvar_srf(9)%standard_name = 'net_upward_longwave_flux_in_air'
        v2dvar_srf(9)%lrecords = .true.
        v2dvar_srf(9)%cell_method = 'time: mean'
        v2dvar_srf(10)%vname = 'fsw'
        v2dvar_srf(10)%vunit = 'W m-2'
        v2dvar_srf(10)%long_name = 'Net downward shortwave energy flux'
        v2dvar_srf(10)%standard_name = 'net_downward_shortwave_flux_in_air'
        v2dvar_srf(10)%lrecords = .true.
        v2dvar_srf(10)%cell_method = 'time: mean'
        v2dvar_srf(11)%vname = 'fld'
        v2dvar_srf(11)%vunit = 'W m-2'
        v2dvar_srf(11)%long_name = 'SUrface downward longwave flux in air'
        v2dvar_srf(11)%standard_name = &
          'surface_downwelling_longwave_flux_in_air'
        v2dvar_srf(11)%lrecords = .true.
        v2dvar_srf(11)%cell_method = 'time: mean'
        v2dvar_srf(12)%vname = 'sina'
        v2dvar_srf(12)%vunit = 'W m-2'
        v2dvar_srf(12)%long_name = 'Surface downward shortwave flux in air'
        v2dvar_srf(12)%standard_name = &
          'surface_downwelling_shortwave_flux_in_air'
        v2dvar_srf(12)%lrecords = .true.
        v2dvar_srf(12)%cell_method = 'time: mean'
        v2dvar_srf(13)%vname = 'prcv'
        v2dvar_srf(13)%vunit = 'kg m-2 day-1'
        v2dvar_srf(13)%long_name = 'Convective precipitation flux'
        v2dvar_srf(13)%standard_name = 'convective_rainfall_flux'
        v2dvar_srf(13)%lrecords = .true.
        v2dvar_srf(13)%cell_method = 'time: mean'
        v2dvar_srf(14)%vname = 'zpbl'
        v2dvar_srf(14)%vunit = 'm'
        v2dvar_srf(14)%long_name = 'PBL thickness'
        v2dvar_srf(14)%standard_name = 'atmosphere_boundary_layer_thickness'
        v2dvar_srf(14)%lrecords = .true.
        v2dvar_srf(15)%vname = 'aldirs'
        v2dvar_srf(15)%vunit = '1'
        v2dvar_srf(15)%long_name = &
          'Surface albedo to direct shortwave radiation'
        v2dvar_srf(15)%standard_name = 'surface_albedo_short_wave_direct'
        v2dvar_srf(15)%lrecords = .true.
        v2dvar_srf(16)%vname = 'aldifs'
        v2dvar_srf(16)%vunit = '1'
        v2dvar_srf(16)%long_name = &
          'Surface albedo to diffuse shortwave radiation'
        v2dvar_srf(16)%standard_name = 'surface_albedo_short_wave_diffuse'
        v2dvar_srf(16)%lrecords = .true.
        v2dvar_srf(17)%vname = 'sund'
        v2dvar_srf(17)%vunit = 's'
        v2dvar_srf(17)%long_name = 'Duration of sunshine'
        v2dvar_srf(17)%standard_name = 'duration_of_sunshine'
        v2dvar_srf(17)%lrecords = .true.
        v2dvar_srf(17)%cell_method = 'time: sum'
        v2dvar_srf(18)%vname = 'seaice'
        v2dvar_srf(18)%vunit = '1'
        v2dvar_srf(18)%long_name = 'Sea ice mask'
        v2dvar_srf(18)%standard_name = 'seaice_binary_mask'
        v2dvar_srf(18)%lrecords = .true.

        v3dvar_srf(1)%vname = 'u10m'
        v3dvar_srf(1)%vunit = 'm s-1'
        v3dvar_srf(1)%long_name = '10 meter zonal (westerly) wind component'
        v3dvar_srf(1)%standard_name = 'eastward_wind'
        v3dvar_srf(1)%lrecords = .true.
        v3dvar_srf(1)%axis = 'xyw'
        v3dvar_srf(2)%vname = 'v10m'
        v3dvar_srf(2)%vunit = 'm s-1'
        v3dvar_srf(2)%long_name = &
          '10 meter meridional (southerly) wind component'
        v3dvar_srf(2)%standard_name = 'northward_wind'
        v3dvar_srf(2)%lrecords = .true.
        v3dvar_srf(2)%axis = 'xyw'
        v3dvar_srf(3)%vname = 't2m'
        v3dvar_srf(3)%vunit = 'K'
        v3dvar_srf(3)%long_name = '2 meter air temperature'
        v3dvar_srf(3)%standard_name = 'air_temperature'
        v3dvar_srf(3)%lrecords = .true.
        v3dvar_srf(3)%axis = 'xy2'
        v3dvar_srf(4)%vname = 'q2m'
        v3dvar_srf(4)%vunit = '1'
        v3dvar_srf(4)%long_name = '2 meter air specific humidity'
        v3dvar_srf(4)%standard_name = 'specific_humidity'
        v3dvar_srf(4)%lrecords = .true.
        v3dvar_srf(4)%axis = 'xy2'
        v3dvar_srf(5)%vname = 'smw'
        v3dvar_srf(5)%vunit = 'kg m-2'
        v3dvar_srf(5)%long_name = 'Moisture content of the soil layers'
        v3dvar_srf(5)%standard_name = 'soil_moisture_content'
        v3dvar_srf(5)%lrecords = .true.
        v3dvar_srf(5)%axis = 'xys'

        outstream(srf_stream)%nvar = countvars(enable_srf_vars,nsrfvars)
        allocate(outstream(srf_stream)%ncvars%vlist(outstream(srf_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(srf_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , nsrf2dvars
          if ( enable_srf_vars(i) ) then
            outstream(srf_stream)%ncvars%vlist(vcount)%vp => v2dvar_srf(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsrf3dvars
          if ( enable_srf_vars(i+nsrf2dvars) ) then
            outstream(srf_stream)%ncvars%vlist(vcount)%vp => v3dvar_srf(i)
            vcount = vcount + 1
          end if
        end do
      end if

      if ( nstream == sts_stream ) then

        ! The enable flag respect the below order

        v2dvar_sts(1)%vname = 'tgmax'
        v2dvar_sts(1)%vunit = 'K'
        v2dvar_sts(1)%long_name = 'Maximum surface temperature'
        v2dvar_sts(1)%standard_name = 'surface_temperature'
        v2dvar_sts(1)%lrecords = .true.
        v2dvar_sts(1)%cell_method = 'time: maximum'
        v2dvar_sts(2)%vname = 'tgmin'
        v2dvar_sts(2)%vunit = 'K'
        v2dvar_sts(2)%long_name = 'Minimum surface temperature'
        v2dvar_sts(2)%standard_name = 'surface_temperature'
        v2dvar_sts(2)%lrecords = .true.
        v2dvar_sts(2)%cell_method = 'time: minimum'
        v2dvar_sts(3)%vname = 'pcpmax'
        v2dvar_sts(3)%vunit = 'kg m-2 s-1'
        v2dvar_sts(3)%long_name = 'Maximum total precipitation flux'
        v2dvar_sts(3)%standard_name = 'precipitation_flux'
        v2dvar_sts(3)%lrecords = .true.
        v2dvar_sts(3)%cell_method = 'time: maximum'
        v2dvar_sts(4)%vname = 'pcpavg'
        v2dvar_sts(4)%vunit = 'kg m-2 s-1'
        v2dvar_sts(4)%long_name = 'Mean total precipitation flux'
        v2dvar_sts(4)%standard_name = 'precipitation_flux'
        v2dvar_sts(4)%lrecords = .true.
        v2dvar_sts(4)%cell_method = 'time: mean'
        v2dvar_sts(5)%vname = 'sund'
        v2dvar_sts(5)%vunit = 's'
        v2dvar_sts(5)%long_name = 'Duration of sunshine'
        v2dvar_sts(5)%standard_name = 'duration_of_sunshine'
        v2dvar_sts(5)%lrecords = .true.
        v2dvar_sts(5)%cell_method = 'time: sum'
        v2dvar_sts(6)%vname = 'ps_min'
        v2dvar_sts(6)%vunit = 'hPa'
        v2dvar_sts(6)%long_name = 'Minimum of surface pressure'
        v2dvar_sts(6)%standard_name = 'air_pressure'
        v2dvar_sts(6)%lrecords = .true.
        v2dvar_sts(6)%cell_method = 'time: minimum'

        v3dvar_sts(1)%vname = 't2max'
        v3dvar_sts(1)%vunit = 'K'
        v3dvar_sts(1)%long_name = 'Maximum 2 meter temperature'
        v3dvar_sts(1)%standard_name = 'air_temperature'
        v3dvar_sts(1)%lrecords = .true.
        v3dvar_sts(1)%cell_method = 'time: maximum'
        v3dvar_sts(1)%axis = 'xy2'
        v3dvar_sts(2)%vname = 't2min'
        v3dvar_sts(2)%vunit = 'K'
        v3dvar_sts(2)%long_name = 'Minimum 2 meter temperature'
        v3dvar_sts(2)%standard_name = 'air_temperature'
        v3dvar_sts(2)%lrecords = .true.
        v3dvar_sts(2)%cell_method = 'time: minimum'
        v3dvar_sts(2)%axis = 'xy2'
        v3dvar_sts(3)%vname = 't2avg'
        v3dvar_sts(3)%vunit = 'K'
        v3dvar_sts(3)%long_name = 'Mean 2 meter temperature'
        v3dvar_sts(3)%standard_name = 'air_temperature'
        v3dvar_sts(3)%lrecords = .true.
        v3dvar_sts(3)%cell_method = 'time: mean'
        v3dvar_sts(3)%axis = 'xy2'
        v3dvar_sts(4)%vname = 'w10max'
        v3dvar_sts(4)%vunit = 'm s-1'
        v3dvar_sts(4)%long_name = 'Maximum speed of 10m wind'
        v3dvar_sts(4)%standard_name = 'wind_speed'
        v3dvar_sts(4)%lrecords = .true.
        v3dvar_sts(4)%cell_method = 'time: maximum'
        v3dvar_sts(4)%axis = 'xyw'

        outstream(sts_stream)%nvar = countvars(enable_sts_vars,nstsvars)
        allocate(outstream(sts_stream)%ncvars%vlist(outstream(sts_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(sts_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , nsts2dvars
          if ( enable_sts_vars(i) ) then
            outstream(sts_stream)%ncvars%vlist(vcount)%vp => v2dvar_sts(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsts3dvars
          if ( enable_sts_vars(i+nsts2dvars) ) then
            outstream(sts_stream)%ncvars%vlist(vcount)%vp => v3dvar_sts(i)
            vcount = vcount + 1
          end if
        end do
      end if

      if ( nstream == sub_stream ) then

        ! The enable flag respect the below order

        v2dvar_sub(1)%vname = 'uvdrag'
        v2dvar_sub(1)%vunit = '1'
        v2dvar_sub(1)%long_name = 'Surface drag stress coefficient'
        v2dvar_sub(1)%standard_name = 'surface_drag_coefficient_in_air'
        v2dvar_sub(1)%lrecords = .true.
        v2dvar_sub(2)%vname = 'tg'
        v2dvar_sub(2)%vunit = 'K'
        v2dvar_sub(2)%long_name = 'Ground temperature'
        v2dvar_sub(2)%standard_name = 'surface_temperature'
        v2dvar_sub(2)%lrecords = .true.
        v2dvar_sub(3)%vname = 'tlef'
        v2dvar_sub(3)%vunit = 'K'
        v2dvar_sub(3)%long_name = 'Foliage canopy temperature'
        v2dvar_sub(3)%standard_name = 'canopy_temperature'
        v2dvar_sub(3)%lrecords = .true.
        v2dvar_sub(4)%vname = 'tpr'
        v2dvar_sub(4)%vunit = 'kg m-2 day-1'
        v2dvar_sub(4)%long_name = 'Total precipitation flux'
        v2dvar_sub(4)%standard_name = 'precipitation_flux'
        v2dvar_sub(4)%lrecords = .true.
        v2dvar_sub(4)%cell_method = 'time: mean'
        v2dvar_sub(5)%vname = 'evp'
        v2dvar_sub(5)%vunit = 'kg m-2 day-1'
        v2dvar_sub(5)%long_name = 'Total evapotranspiration flux'
        v2dvar_sub(5)%standard_name = 'water_evaporation_flux'
        v2dvar_sub(5)%cell_method = 'time: mean'
        v2dvar_sub(5)%lrecords = .true.
        v2dvar_sub(6)%vname = 'runoff'
        v2dvar_sub(6)%vunit = 'kg m-2 day-1'
        v2dvar_sub(6)%long_name = 'Surface runoff flux'
        v2dvar_sub(6)%standard_name = 'surface_runoff_flux'
        v2dvar_sub(6)%lrecords = .true.
        v2dvar_sub(6)%cell_method = 'time: mean'
        v2dvar_sub(7)%vname = 'scv'
        v2dvar_sub(7)%vunit = 'kg m-2'
        v2dvar_sub(7)%long_name = 'Liquid water equivalent of snow thickness'
        v2dvar_sub(7)%standard_name = 'lwe_thickness_of_surface_snow_amount'
        v2dvar_sub(7)%lrecords = .true.
        v2dvar_sub(7)%cell_method = 'time: mean'
        v2dvar_sub(8)%vname = 'sena'
        v2dvar_sub(8)%vunit = 'W m-2'
        v2dvar_sub(8)%long_name = 'Sensible heat flux'
        v2dvar_sub(8)%standard_name = 'surface_downward_sensible_heat_flux'
        v2dvar_sub(8)%lrecords = .true.
        v2dvar_sub(8)%cell_method = 'time: mean'
        v2dvar_sub(9)%vname = 'prcv'
        v2dvar_sub(9)%vunit = 'kg m-2 day-1'
        v2dvar_sub(9)%long_name = 'Convective precipitation flux'
        v2dvar_sub(9)%standard_name = 'convective_rainfall_flux'
        v2dvar_sub(9)%lrecords = .true.
        v2dvar_sub(9)%cell_method = 'time: mean'
        v2dvar_sub(10)%vname = 'tlake'
        v2dvar_sub(10)%vunit = 'K'
        v2dvar_sub(10)%long_name = 'Lake water surface temperature'
        v2dvar_sub(10)%standard_name = 'water_temperature'
        v2dvar_sub(10)%lrecords = .true.

        v3dvar_sub(1)%vname = 'u10m'
        v3dvar_sub(1)%vunit = 'm s-1'
        v3dvar_sub(1)%long_name = '10 meter zonal wind component (westerly)'
        v3dvar_sub(1)%standard_name = 'eastward_wind'
        v3dvar_sub(1)%lrecords = .true.
        v3dvar_sub(1)%axis = 'xyw'
        v3dvar_sub(2)%vname = 'v10m'
        v3dvar_sub(2)%vunit = 'm s-1'
        v3dvar_sub(2)%long_name = &
          '10 meter meridional wind component (southerly)'
        v3dvar_sub(2)%standard_name = 'northward_wind'
        v3dvar_sub(2)%lrecords = .true.
        v3dvar_sub(2)%axis = 'xyw'
        v3dvar_sub(3)%vname = 't2m'
        v3dvar_sub(3)%vunit = 'K'
        v3dvar_sub(3)%long_name = '2 meter air temperature'
        v3dvar_sub(3)%standard_name = 'air_temperature'
        v3dvar_sub(3)%lrecords = .true.
        v3dvar_sub(3)%axis = 'xy2'
        v3dvar_sub(3)%vname = 'q2m'
        v3dvar_sub(3)%vunit = '1'
        v3dvar_sub(3)%long_name = '2 meter air specific humidity'
        v3dvar_sub(3)%standard_name = 'specific_humidity'
        v3dvar_sub(3)%lrecords = .true.
        v3dvar_sub(3)%axis = 'xy2'
        v3dvar_sub(4)%vname = 'smw'
        v3dvar_sub(4)%vunit = 'kg kg-1'
        v3dvar_sub(4)%long_name = 'Soil moisture content'
        v3dvar_sub(4)%standard_name = 'soil_moisture_content'
        v3dvar_sub(4)%lrecords = .true.
        v3dvar_sub(4)%axis = 'xys'

        outstream(sub_stream)%nvar = countvars(enable_sub_vars,nsubvars)
        allocate(outstream(sub_stream)%ncvars%vlist(outstream(sub_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(sub_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , nsub2dvars
          if ( enable_sub_vars(i) ) then
            outstream(sub_stream)%ncvars%vlist(vcount)%vp => v2dvar_sub(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsub3dvars
          if ( enable_sub_vars(i+nsub2dvars) ) then
            outstream(sub_stream)%ncvars%vlist(vcount)%vp => v3dvar_sub(i)
            vcount = vcount + 1
          end if
        end do
      end if

      if ( nstream == rad_stream ) then

        ! The enable flag respect the below order

        v2dvar_rad(1)%vname = 'frsa'
        v2dvar_rad(1)%vunit = 'W m-2'
        v2dvar_rad(1)%long_name = 'Surface net downward shortwave flux'
        v2dvar_rad(1)%standard_name = 'surface_net_downward_shortwave_flux'
        v2dvar_rad(1)%lrecords = .true.
        v2dvar_rad(2)%vname = 'frla'
        v2dvar_rad(2)%vunit = 'W m-2'
        v2dvar_rad(2)%long_name = 'Surface net upward longwave flux'
        v2dvar_rad(2)%standard_name = 'surface_net_upward_longwave_flux'
        v2dvar_rad(2)%lrecords = .true.
        v2dvar_rad(3)%vname = 'clrst'
        v2dvar_rad(3)%vunit = 'W m-2'
        v2dvar_rad(3)%long_name = &
          'Clearsky top of atmosphere net downward shortwave flux'
        v2dvar_rad(3)%standard_name =  &
          'toa_net_downward_shortwave_flux_assuming_clear_sky'
        v2dvar_rad(3)%lrecords = .true.
        v2dvar_rad(4)%vname = 'clrss'
        v2dvar_rad(4)%vunit = 'W m-2'
        v2dvar_rad(4)%long_name = 'Clearsky surface net downward shortwave flux'
        v2dvar_rad(4)%standard_name =  &
          'surface_net_downward_shortwave_flux_assuming_clear_sky'
        v2dvar_rad(4)%lrecords = .true.
        v2dvar_rad(5)%vname = 'clrls'
        v2dvar_rad(5)%vunit = 'W m-2'
        v2dvar_rad(5)%long_name = 'Clearsky net upward longwave flux'
        v2dvar_rad(5)%standard_name =  &
          'surface_net_upward_longwave_flux_assuming_clear_sky'
        v2dvar_rad(5)%lrecords = .true.
        v2dvar_rad(6)%vname = 'clrlt'
        v2dvar_rad(6)%vunit = 'W m-2'
        v2dvar_rad(6)%long_name = &
          'Clearsky top of atmosphere net upward longwave flux'
        v2dvar_rad(6)%standard_name =  &
          'toa_net_upward_longwave_flux_assuming_clear_sky'
        v2dvar_rad(6)%lrecords = .true.
        v2dvar_rad(7)%vname = 'solin'
        v2dvar_rad(7)%vunit = 'W m-2'
        v2dvar_rad(7)%long_name = 'Top of atmosphere incoming shortwave flux'
        v2dvar_rad(7)%standard_name = 'toa_incoming_shortwave_flux'
        v2dvar_rad(7)%lrecords = .true.
        v2dvar_rad(8)%vname = 'sabtp'
        v2dvar_rad(8)%vunit = 'W m-2'
        v2dvar_rad(8)%long_name = 'Net top of atmosphere upward shortwave flux'
        v2dvar_rad(8)%standard_name = 'toa_net_upward_shortwave_flux'
        v2dvar_rad(8)%lrecords = .true.
        v2dvar_rad(9)%vname = 'totcf'
        v2dvar_rad(9)%vunit = '1'
        v2dvar_rad(9)%long_name = 'Total cloud fraction'
        v2dvar_rad(9)%standard_name = 'cloud_area_fraction'
        v2dvar_rad(9)%lrecords = .true.
        v2dvar_rad(10)%vname = 'totcl'
        v2dvar_rad(10)%vunit = 'kg m-2'
        v2dvar_rad(10)%long_name = 'Total columnar liquid water content'
        v2dvar_rad(10)%standard_name = &
          'atmosphere_cloud_condensed_water_content'
        v2dvar_rad(10)%lrecords = .true.
        v2dvar_rad(11)%vname = 'totci'
        v2dvar_rad(11)%vunit = 'kg m-2'
        v2dvar_rad(11)%long_name = 'Total columnar ice water content'
        v2dvar_rad(11)%standard_name = 'atmosphere_ice_condensed_water_content'
        v2dvar_rad(11)%lrecords = .true.
        v2dvar_rad(12)%vname = 'firtp'
        v2dvar_rad(12)%vunit = 'W m-2'
        v2dvar_rad(12)%long_name = 'Top of atmosphere net upward longwave flux'
        v2dvar_rad(12)%standard_name = 'toa_net_upward_longwave_flux'
        v2dvar_rad(12)%lrecords = .true.

        v3dvar_rad(1)%vname = 'cld'
        v3dvar_rad(1)%vunit = '1'
        v3dvar_rad(1)%long_name = 'Cloud fractional cover'
        v3dvar_rad(1)%standard_name = 'cloud_area_fraction_in_atmosphere_layer'
        v3dvar_rad(1)%lrecords = .true.
        v3dvar_rad(2)%vname = 'clwp'
        v3dvar_rad(2)%vunit = 'g m-2'
        v3dvar_rad(2)%long_name = 'Cloud liquid water path'
        v3dvar_rad(2)%standard_name = 'thickness_of_liquid_water_cloud'
        v3dvar_rad(2)%lrecords = .true.
        v3dvar_rad(3)%vname = 'qrs'
        v3dvar_rad(3)%vunit = 'K s-1'
        v3dvar_rad(3)%long_name = 'Shortwave radiation heating rate'
        v3dvar_rad(3)%standard_name = &
          'tendency_of_air_temperature_due_to_shortwave_heating'
        v3dvar_rad(3)%lrecords = .true.
        v3dvar_rad(4)%vname = 'qrl'
        v3dvar_rad(4)%vunit = 'K s-1'
        v3dvar_rad(4)%long_name = 'Longwave radiation heating rate'
        v3dvar_rad(4)%standard_name = &
          'tendency_of_air_temperature_due_to_longwave_heating'
        v3dvar_rad(4)%lrecords = .true.

        outstream(rad_stream)%nvar = countvars(enable_rad_vars,nradvars)
        allocate(outstream(rad_stream)%ncvars%vlist(outstream(rad_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(rad_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , nrad2dvars
          if ( enable_rad_vars(i) ) then
            outstream(rad_stream)%ncvars%vlist(vcount)%vp => v2dvar_rad(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nrad3dvars
          if ( enable_rad_vars(i+nrad2dvars) ) then
            outstream(rad_stream)%ncvars%vlist(vcount)%vp => v3dvar_rad(i)
            vcount = vcount + 1
          end if
        end do
      end if

      if ( nstream == lak_stream ) then

        ! The enable flag respect the below order

        v2dvar_lak(1)%vname = 'tg'
        v2dvar_lak(1)%vunit = 'K'
        v2dvar_lak(1)%long_name = 'Ground temperature'
        v2dvar_lak(1)%standard_name = 'surface_temperature'
        v2dvar_lak(1)%lrecords = .true.
        v2dvar_lak(2)%vname = 'tpr'
        v2dvar_lak(2)%vunit = 'kg m-2 day-1'
        v2dvar_lak(2)%long_name = 'Total precipitation flux'
        v2dvar_lak(2)%standard_name = 'precipitation_flux'
        v2dvar_lak(2)%lrecords = .true.
        v2dvar_lak(3)%vname = 'scv'
        v2dvar_lak(3)%vunit = 'kg m-2'
        v2dvar_lak(3)%long_name = 'Liquid water equivalent snow depth'
        v2dvar_lak(3)%standard_name = 'lwe_thickness_of_surface_snow_amount'
        v2dvar_lak(3)%lrecords = .true.
        v2dvar_lak(4)%vname = 'sena'
        v2dvar_lak(4)%vunit = 'W m-2'
        v2dvar_lak(4)%long_name = 'Sensible heat flux'
        v2dvar_lak(4)%standard_name = 'surface_downward_sensible_heat_flux'
        v2dvar_lak(4)%lrecords = .true.
        v2dvar_lak(5)%vname = 'flw'
        v2dvar_lak(5)%vunit = 'W m-2'
        v2dvar_lak(5)%long_name = 'Net longwave energy flux'
        v2dvar_lak(5)%standard_name = 'net_upward_longwave_flux_in_air'
        v2dvar_lak(5)%lrecords = .true.
        v2dvar_lak(6)%vname = 'fsw'
        v2dvar_lak(6)%vunit = 'W m-2'
        v2dvar_lak(6)%long_name = 'Net shortwave absorbed energy flux'
        v2dvar_lak(6)%standard_name = 'net_downward_shortwave_flux_in_air'
        v2dvar_lak(6)%lrecords = .true.
        v2dvar_lak(7)%vname = 'fld'
        v2dvar_lak(7)%vunit = 'W m-2'
        v2dvar_lak(7)%long_name = 'Downward longwave flux at surface in air'
        v2dvar_lak(7)%standard_name = 'surface_downwelling_longwave_flux_in_air'
        v2dvar_lak(7)%lrecords = .true.
        v2dvar_lak(8)%vname = 'sina'
        v2dvar_lak(8)%vunit = 'W m-2'
        v2dvar_lak(8)%long_name = 'Downward shortwave flux at surface in air'
        v2dvar_lak(8)%standard_name = &
          'surface_downwelling_shortwave_flux_in_air'
        v2dvar_lak(8)%lrecords = .true.
        v2dvar_lak(9)%vname = 'aldirs'
        v2dvar_lak(9)%vunit = '1'
        v2dvar_lak(9)%long_name = 'Surface albedo to direct shortwave radiation'
        v2dvar_lak(9)%standard_name = 'surface_albedo_short_wave_direct'
        v2dvar_lak(9)%lrecords = .true.
        v2dvar_lak(10)%vname = 'aldifs'
        v2dvar_lak(10)%vunit = '1'
        v2dvar_lak(10)%long_name = &
          'Surface albedo to diffuse shortwave radiation'
        v2dvar_lak(10)%standard_name = 'surface_albedo_short_wave_diffuse'
        v2dvar_lak(10)%lrecords = .true.
        v2dvar_lak(11)%vname = 'evp'
        v2dvar_lak(11)%vunit = 'mm s-1'
        v2dvar_lak(11)%long_name = 'Water evaporation flux'
        v2dvar_lak(11)%standard_name = 'water_evaporation_flux_where_sea_ice'
        v2dvar_lak(11)%lrecords = .true.
        v2dvar_lak(12)%vname = 'aveice'
        v2dvar_lak(12)%vunit = 'mm'
        v2dvar_lak(12)%long_name = 'Floating ice thickness'
        v2dvar_lak(12)%standard_name = 'floating_ice_thickness'
        v2dvar_lak(12)%lrecords = .true.
        v2dvar_lak(13)%vname = 'hsnow'
        v2dvar_lak(13)%vunit = 'mm'
        v2dvar_lak(13)%long_name = 'Floating snow thickness'
        v2dvar_lak(13)%standard_name = 'surface_snow_thickness_where_sea_ice'
        v2dvar_lak(13)%lrecords = .true.

        v3dvar_lak(1)%vname = 'tlake'
        v3dvar_lak(1)%vunit = 'K'
        v3dvar_lak(1)%long_name = 'Lake water temperature'
        v3dvar_lak(1)%standard_name = 'water_temperature'
        v3dvar_lak(1)%lrecords = .true.
        v3dvar_lak(1)%axis = 'xyd'

        outstream(lak_stream)%nvar = countvars(enable_lak_vars,nlakvars)
        allocate(outstream(lak_stream)%ncvars%vlist(outstream(lak_stream)%nvar))

        vcount = 1
        do i = 1 , nbasevars
          outstream(lak_stream)%ncvars%vlist(vcount)%vp => v2dvar_base(i)
          vcount = vcount + 1
        end do
        do i = 1 , nlak2dvars
          if ( enable_lak_vars(i) ) then
            outstream(lak_stream)%ncvars%vlist(vcount)%vp => v2dvar_lak(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nlak3dvars
          if ( enable_lak_vars(i+nlak2dvars) ) then
            outstream(lak_stream)%ncvars%vlist(vcount)%vp => v3dvar_lak(i)
            vcount = vcount + 1
          end if
        end do
      end if

      outstream(nstream)%opar%pname = 'RegCM Model'
      outstream(nstream)%opar%zero_date = idate1
    end do

  end subroutine init_output_streams

  integer function countvars(eflags,ntot)
    implicit none
    logical , dimension(ntot) , intent(in) :: eflags
    integer , intent(in) :: ntot
    integer :: i
    countvars = nbasevars
    do i = 1 , ntot
      if ( eflags(i) ) countvars = countvars + 1
    end do
  end function countvars

  subroutine newoutfile(idate,itype)
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=16) :: fbname
    integer , intent(in) :: itype

    select case (itype)
      case (atm_stream)
        write (fbname,'(a,i10)') 'ATM.', toint10(idate)
      case (srf_stream)
        write (fbname,'(a,i10)') 'SRF.', toint10(idate)
      case (sub_stream)
        write (fbname,'(a,i10)') 'SUB.', toint10(idate)
      case (rad_stream)
        write (fbname,'(a,i10)') 'RAD.', toint10(idate)
      case (lak_stream)
        write (fbname,'(a,i10)') 'LAK.', toint10(idate)
      case (sts_stream)
        write (fbname,'(a,i10)') 'STS.', toint10(idate)
      case default
        write(stderr,*) 'Undefined output stream. Skipping it.'
        return
    end select
    outstream(itype)%opar%fname = &
      trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
    call outstream_setup(outstream(itype)%ncout,outstream(itype)%opar)

    ! Buffer Zone Control relaxation + diffusion term params

    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('boundary_nspgx',nspgx))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('boundary_nspgd',nspgd))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('boundary_high_nudge',high_nudge))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('boundary_medium_nudge',medium_nudge))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('boundary_low_nudge',low_nudge))

    ! Perturbation control for ensembles

    if ( ensemble_run ) then
      if ( lperturb_topo ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_topo_percent',perturb_frac_topo))
      end if
      if ( lperturb_ts ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_ts_percent',perturb_frac_ts))
      end if
      if ( lperturb_ps ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_ps_percent',perturb_frac_ps))
      end if
      if ( lperturb_t ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_t_percent',perturb_frac_t))
      end if
      if ( lperturb_u ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_u_percent',perturb_frac_u))
      end if
      if ( lperturb_v ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_v_percent',perturb_frac_v))
      end if
      if ( lperturb_q ) then
        call outstream_addatt(outstream(itype)%ncout, &
          ncattribute_real8('perturbation_q_percent',perturb_frac_q))
      end if
    end if

    ! Model start/restart control

    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_logical('model_is_restarted',ifrest))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_initial_start',tochar(idate0)))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_start',tochar(idate1)))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_end',tochar(idate2)))

    ! Model timing parameters

    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('dynamic_time_step_in_seconds',dt))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('surface_interaction_time_step_in_seconds',dtsrf))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('radiation_scheme_time_step_in_minuts',dtrad))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_real8('absorption_emission_time_step_in_hours',dtabem))

    ! Model Physics

    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('lateral_boundary_condition_scheme',iboudy))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('boundary_layer_scheme',ibltyp))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('cumulus_convection_scheme',icup))
    if ( icup == 2 .or. icup == 98 .or. icup == 99 ) then
      call outstream_addatt(outstream(itype)%ncout, &
        ncattribute_integer('grell_scheme_closure',igcc))
    end if
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('moisture_scheme',ipptls))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('ocean_flux_scheme',iocnflx))
    if ( iocnflx == 2 ) then
      call outstream_addatt(outstream(itype)%ncout, &
        ncattribute_integer('zeng_ocean_roughness_formula',iocnrough))
    end if
    if ( iocncpl == 1 ) then
      call outstream_addatt(outstream(itype)%ncout, &
        ncattribute_integer('coupled_ocean_run',iocncpl))
    end if
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('pressure_gradient_scheme',ipgf))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('surface_emissivity_factor_computed',iemiss))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('lake_model_activated',lakemod))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('chemical_aerosol_scheme_activated',ichem))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('ipcc_scenario_code',scenario))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('diurnal_cycle_sst_scheme',idcsst))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('simple_sea_ice_scheme',iseaice))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('seasonal_desert_albedo',idesseas))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('convective_lwp_as_large_scale',iconvlwp))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('rrtm_radiation_scheme_activated',irrtm))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('climatic_ozone_input_dataset',iclimao3))
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_integer('static_solar_constant_used',isolconst))

  end subroutine newoutfile

end module mod_ncout
