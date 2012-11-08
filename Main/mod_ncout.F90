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
  use mod_atm_interface
  use mod_lm_interface
  use mod_cu_interface
  use mod_rad_interface
  use mod_precip

  use mpi
  use netcdf
!
  public :: init_output_streams
  public :: dispose_output_streams
  public :: write_record_output_stream
!
  type varspan
    integer(ik4) :: j1 , j2 , i1 , i2 , k1 , k2
  end type varspan

  type regcm_stream
    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    integer(ik4) :: nvar = 0
    type(nc_varlist) :: ncvars
    integer(ik4) :: jg1 , jg2 , ig1 , ig2
    integer(ik4) :: jl1 , jl2 , il1 , il2
  end type regcm_stream

  integer(ik4) , parameter :: nbase = 5

  integer(ik4) , parameter :: natm2dvars = 2 + nbase
  integer(ik4) , parameter :: natm3dvars = 10
  integer(ik4) , parameter :: natmvars = natm2dvars+natm3dvars

  integer(ik4) , parameter :: nsrf2dvars = 17 + nbase
  integer(ik4) , parameter :: nsrf3dvars = 6
  integer(ik4) , parameter :: nsrfvars = nsrf2dvars+nsrf3dvars

  integer(ik4) , parameter :: nsts2dvars = 6 + nbase
  integer(ik4) , parameter :: nsts3dvars = 4
  integer(ik4) , parameter :: nstsvars = nsts2dvars+nsts3dvars

  integer(ik4) , parameter :: nsub2dvars = 7 + nbase
  integer(ik4) , parameter :: nsub3dvars = 6
  integer(ik4) , parameter :: nsubvars = nsub2dvars+nsub3dvars

  integer(ik4) , parameter :: nlak2dvars = 13 + nbase
  integer(ik4) , parameter :: nlak3dvars = 1
  integer(ik4) , parameter :: nlakvars = nlak2dvars+nlak3dvars

  integer(ik4) , parameter :: nrad2dvars = 12 + nbase
  integer(ik4) , parameter :: nrad3dvars = 4
  integer(ik4) , parameter :: nradvars = nrad2dvars+nrad3dvars

  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_atm => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_atm => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_srf => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_srf => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_sts => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_sts => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_sub => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_sub => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_lak => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_lak => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_rad => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_rad => null()

  integer(ik4) , parameter :: maxstreams = 6

  integer(ik4) , public , parameter :: atm_stream = 1
  integer(ik4) , public , parameter :: srf_stream = 2
  integer(ik4) , public , parameter :: sub_stream = 3
  integer(ik4) , public , parameter :: rad_stream = 4
  integer(ik4) , public , parameter :: lak_stream = 5
  integer(ik4) , public , parameter :: sts_stream = 6
!
  type(regcm_stream) , save , dimension(maxstreams) :: outstream

  logical , public , dimension(maxstreams) :: enable_flag
  logical , public , dimension(natmvars) :: enable_atm_vars
  logical , public , dimension(nsrfvars) :: enable_srf_vars
  logical , public , dimension(nstsvars) :: enable_sts_vars
  logical , public , dimension(nsubvars) :: enable_sub_vars
  logical , public , dimension(nlakvars) :: enable_lak_vars
  logical , public , dimension(nradvars) :: enable_rad_vars

  integer(ik4), parameter :: atm_xlon  = 1
  integer(ik4), parameter :: atm_xlat  = 2
  integer(ik4), parameter :: atm_mask  = 3
  integer(ik4), parameter :: atm_topo  = 4
  integer(ik4), parameter :: atm_ps    = 5
  integer(ik4), parameter :: atm_tpr   = 6
  integer(ik4), parameter :: atm_tgb   = 7

  integer(ik4), parameter :: atm_u     = 1
  integer(ik4), parameter :: atm_v     = 2
  integer(ik4), parameter :: atm_t     = 3
  integer(ik4), parameter :: atm_omega = 4
  integer(ik4), parameter :: atm_qv    = 5
  integer(ik4), parameter :: atm_qc    = 6
  integer(ik4), parameter :: atm_tke   = 7
  integer(ik4), parameter :: atm_kth   = 8
  integer(ik4), parameter :: atm_kzm   = 9
  integer(ik4), parameter :: atm_swt   = 10

  integer(ik4), parameter :: srf_xlon   = 1
  integer(ik4), parameter :: srf_xlat   = 2
  integer(ik4), parameter :: srf_mask   = 3
  integer(ik4), parameter :: srf_topo   = 4
  integer(ik4), parameter :: srf_ps     = 5
  integer(ik4), parameter :: srf_uvdrag = 6
  integer(ik4), parameter :: srf_tg     = 7
  integer(ik4), parameter :: srf_tlef   = 8
  integer(ik4), parameter :: srf_tpr    = 9
  integer(ik4), parameter :: srf_evp    = 10
  integer(ik4), parameter :: srf_scv    = 11
  integer(ik4), parameter :: srf_sena   = 12
  integer(ik4), parameter :: srf_flw    = 13
  integer(ik4), parameter :: srf_fsw    = 14
  integer(ik4), parameter :: srf_fld    = 15
  integer(ik4), parameter :: srf_sina   = 16
  integer(ik4), parameter :: srf_prcv   = 17
  integer(ik4), parameter :: srf_zpbl   = 18
  integer(ik4), parameter :: srf_aldirs = 19
  integer(ik4), parameter :: srf_aldifs = 20
  integer(ik4), parameter :: srf_sund   = 21
  integer(ik4), parameter :: srf_seaice = 22

  integer(ik4), parameter :: srf_u10m   = 1
  integer(ik4), parameter :: srf_v10m   = 2
  integer(ik4), parameter :: srf_t2m    = 3
  integer(ik4), parameter :: srf_q2m    = 4
  integer(ik4), parameter :: srf_smw    = 5
  integer(ik4), parameter :: srf_runoff = 6

  integer(ik4), parameter :: sts_xlon   = 1
  integer(ik4), parameter :: sts_xlat   = 2
  integer(ik4), parameter :: sts_mask   = 3
  integer(ik4), parameter :: sts_topo   = 4
  integer(ik4), parameter :: sts_ps     = 5
  integer(ik4), parameter :: sts_tgmax  = 6
  integer(ik4), parameter :: sts_tgmin  = 7
  integer(ik4), parameter :: sts_pcpmax = 8
  integer(ik4), parameter :: sts_pcpavg = 9
  integer(ik4), parameter :: sts_sund   = 10
  integer(ik4), parameter :: sts_psmin  = 11

  integer(ik4), parameter :: sts_t2max  = 1
  integer(ik4), parameter :: sts_t2min  = 2
  integer(ik4), parameter :: sts_t2avg  = 3
  integer(ik4), parameter :: sts_w10max = 4

  integer(ik4), parameter :: sub_xlon   = 1
  integer(ik4), parameter :: sub_xlat   = 2
  integer(ik4), parameter :: sub_mask   = 3
  integer(ik4), parameter :: sub_topo   = 4
  integer(ik4), parameter :: sub_ps     = 5
  integer(ik4), parameter :: sub_uvdrag = 6
  integer(ik4), parameter :: sub_tg     = 7
  integer(ik4), parameter :: sub_tlef   = 8
  integer(ik4), parameter :: sub_evp    = 9
  integer(ik4), parameter :: sub_scv    = 10
  integer(ik4), parameter :: sub_sena   = 11
  integer(ik4), parameter :: sub_tlake  = 12

  integer(ik4), parameter :: sub_u10m   = 1
  integer(ik4), parameter :: sub_v10m   = 2
  integer(ik4), parameter :: sub_t2m    = 3
  integer(ik4), parameter :: sub_q2m    = 4
  integer(ik4), parameter :: sub_smw    = 5
  integer(ik4), parameter :: sub_runoff = 6

  integer(ik4), parameter :: rad_xlon   = 1
  integer(ik4), parameter :: rad_xlat   = 2
  integer(ik4), parameter :: rad_mask   = 3
  integer(ik4), parameter :: rad_topo   = 4
  integer(ik4), parameter :: rad_ps     = 5
  integer(ik4), parameter :: rad_frsa   = 6
  integer(ik4), parameter :: rad_frla   = 7
  integer(ik4), parameter :: rad_clrst  = 8
  integer(ik4), parameter :: rad_clrss  = 9
  integer(ik4), parameter :: rad_clrlt  = 10
  integer(ik4), parameter :: rad_clrls  = 11
  integer(ik4), parameter :: rad_solin  = 12
  integer(ik4), parameter :: rad_sabtp  = 13
  integer(ik4), parameter :: rad_totcf  = 14
  integer(ik4), parameter :: rad_totcl  = 15
  integer(ik4), parameter :: rad_totci  = 16
  integer(ik4), parameter :: rad_firtp  = 17

  integer(ik4), parameter :: rad_cld    = 1
  integer(ik4), parameter :: rad_clwp   = 2
  integer(ik4), parameter :: rad_qrs    = 3
  integer(ik4), parameter :: rad_qrl    = 4

  integer(ik4), parameter :: lak_xlon   = 1
  integer(ik4), parameter :: lak_xlat   = 2
  integer(ik4), parameter :: lak_mask   = 3
  integer(ik4), parameter :: lak_topo   = 4
  integer(ik4), parameter :: lak_ps     = 5
  integer(ik4), parameter :: lak_tg     = 6
  integer(ik4), parameter :: lak_tpr    = 7
  integer(ik4), parameter :: lak_scv    = 8
  integer(ik4), parameter :: lak_sena   = 9
  integer(ik4), parameter :: lak_flw    = 10
  integer(ik4), parameter :: lak_fsw    = 11
  integer(ik4), parameter :: lak_fld    = 12
  integer(ik4), parameter :: lak_sina   = 13
  integer(ik4), parameter :: lak_aldirs = 14
  integer(ik4), parameter :: lak_aldifs = 15
  integer(ik4), parameter :: lak_evp    = 16
  integer(ik4), parameter :: lak_aveice = 17
  integer(ik4), parameter :: lak_hsnow  = 18

  integer(ik4), parameter :: lak_tlake  = 1

  real(rk8) , pointer , dimension(:,:) :: io2d , io2dsg
  real(rk8) , pointer , dimension(:,:,:) :: io3d

  contains

  subroutine init_output_streams(lparallel)
    implicit none
    logical , intent(in) :: lparallel
    integer(ik4) :: nstream , i , vcount
    integer(ik4) :: kkz
    type(varspan) :: vsize
    logical , dimension(natm2dvars) :: enable_atm2d_vars
    logical , dimension(natm3dvars) :: enable_atm3d_vars
    logical , dimension(nsrf2dvars) :: enable_srf2d_vars
    logical , dimension(nsrf3dvars) :: enable_srf3d_vars
    logical , dimension(nsts2dvars) :: enable_sts2d_vars
    logical , dimension(nsts3dvars) :: enable_sts3d_vars
    logical , dimension(nsub2dvars) :: enable_sub2d_vars
    logical , dimension(nsub3dvars) :: enable_sub3d_vars
    logical , dimension(nlak2dvars) :: enable_lak2d_vars
    logical , dimension(nlak3dvars) :: enable_lak3d_vars
    logical , dimension(nrad2dvars) :: enable_rad2d_vars
    logical , dimension(nrad3dvars) :: enable_rad3d_vars

    do nstream = 1 , maxstreams

      if ( .not. enable_flag(nstream) ) cycle

      vsize%j1 = jci1
      vsize%j2 = jci2
      vsize%i1 = ici1
      vsize%i2 = ici2
      vsize%k1 = 1

      if ( nstream == atm_stream ) then

        allocate(v2dvar_atm(natm2dvars))
        allocate(v3dvar_atm(natm3dvars))
        enable_atm2d_vars = enable_atm_vars(1:natm2dvars)
        enable_atm3d_vars = enable_atm_vars(natm2dvars+1:natmvars)

        ! This variables are always present

        call setup_var(v2dvar_atm(atm_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_atm(atm_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_atm(atm_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_atm(atm_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_atm(atm_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_atm3d_vars(atm_tpr) ) &
          call setup_var(v2dvar_atm(atm_tpr),vsize,'tpr','kg m-2 day-1', &
            'Total rain precipitation flux','precipitation_flux',.true., &
            'time: mean')
        if ( enable_atm3d_vars(atm_tgb) ) &
          call setup_var(v2dvar_atm(atm_tgb),vsize,'tgb','K', &
            'Lower groud temperature','soil_temperature',.true.)

        vsize%k2 = kz
        if ( enable_atm3d_vars(atm_u) ) &
          call setup_var(v3dvar_atm(atm_u),vsize,'u','m s-1', &
            'Zonal component of wind (westerly)','eastward_wind',.true.)
        if ( enable_atm3d_vars(atm_v) ) &
          call setup_var(v3dvar_atm(atm_v),vsize,'v','m s-1', &
            'Meridional component of wind (southerly)','northward_wind',.true.)
        if ( enable_atm3d_vars(atm_t) ) &
          call setup_var(v3dvar_atm(atm_t),vsize,'t','K', &
            'Air Temperature','air_temperature',.true.)
        if ( enable_atm3d_vars(atm_omega) ) &
          call setup_var(v3dvar_atm(atm_omega),vsize,'omega','hPa s-1', &
            'Pressure velocity','lagrangian_tendency_of_air_pressure',.true.)
        if ( enable_atm3d_vars(atm_qv) ) &
          call setup_var(v3dvar_atm(atm_qv),vsize,'qv','kg kg-1', &
            'Water vapor mixing ratio','humidity_mixing_ratio',.true.)
        if ( enable_atm3d_vars(atm_qc) ) &
          call setup_var(v3dvar_atm(atm_qc),vsize,'qc','kg kg-1', &
            'Cloud liquid water mixing ratio', &
            'cloud_liquid_water_mixing_ratio',.true.)

        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          if ( enable_atm3d_vars(atm_tke) ) &
            call setup_var(v3dvar_atm(atm_tke),vsize,'tke','m2 s2', &
              'Turbulent Kinetic Energy','turbulent_kinetic_energy', .true.)
          if ( enable_atm3d_vars(atm_kth) ) &
            call setup_var(v3dvar_atm(atm_kth),vsize,'kth','m2 s-1', &
              'Vertical Turbulent Viscosity','vertical_momentum_diffusivity', &
              .true.)
          if ( enable_atm3d_vars(atm_kzm) ) &
            call setup_var(v3dvar_atm(atm_kzm),vsize,'kzm','m2 s-1', &
              'Vertical Turbulent Diffusivity','vertical_scalar_diffusivity', &
              .true.)
        else
          enable_atm3d_vars(atm_tke) = .false.
          enable_atm3d_vars(atm_kth) = .false.
          enable_atm3d_vars(atm_kzm) = .false.
        end if

        vsize%k2 = 2
        v3dvar_atm(atm_swt)%axis = 'xys'
        if ( enable_atm3d_vars(atm_swt) ) &
          call setup_var(v3dvar_atm(atm_swt),vsize,'swt','kg m-2', &
            'Total soil water','moisture_content_of_soil_layer',.true.)

        enable_atm_vars(1:natm2dvars) = enable_atm2d_vars
        enable_atm_vars(natm2dvars+1:natmvars) = enable_atm3d_vars
        outstream(atm_stream)%nvar = countvars(enable_atm_vars,natmvars)
        allocate(outstream(atm_stream)%ncvars%vlist(outstream(atm_stream)%nvar))

        vcount = 1
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
        outstream(atm_stream)%jl1 = vsize%j1
        outstream(atm_stream)%jl2 = vsize%j2
        outstream(atm_stream)%il1 = vsize%i1
        outstream(atm_stream)%il2 = vsize%i2
        outstream(atm_stream)%jg1 = jout1
        outstream(atm_stream)%jg2 = jout2
        outstream(atm_stream)%ig1 = iout1
        outstream(atm_stream)%ig2 = iout2
      end if

      if ( nstream == srf_stream ) then

        allocate(v2dvar_srf(nsrf2dvars))
        allocate(v3dvar_srf(nsrf3dvars))
        enable_srf2d_vars = enable_srf_vars(1:nsrf2dvars)
        enable_srf3d_vars = enable_srf_vars(nsrf2dvars+1:nsrfvars)

        ! This variables are always present

        call setup_var(v2dvar_srf(srf_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_srf(srf_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_srf(srf_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_srf(srf_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_srf(srf_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_srf2d_vars(srf_uvdrag) ) &
          call setup_var(v2dvar_srf(srf_uvdrag),vsize,'uvdrag','1', &
            'Surface drag stress coefficient in air', &
            'surface_drag_coefficient_in_air',.true.)
        if ( enable_srf2d_vars(srf_tg) ) &
          call setup_var(v2dvar_srf(srf_tg),vsize,'tg','K', &
            'Ground surface temperature','surface_temperature',.true.)
        if ( enable_srf2d_vars(srf_tlef) ) &
          call setup_var(v2dvar_srf(srf_tlef),vsize,'tlef','K', &
            'Foliage canopy temperature','canopy_temperature',.true.)
        if ( enable_srf2d_vars(srf_tpr) ) &
          call setup_var(v2dvar_srf(srf_tpr),vsize,'tpr','kg m-2 day-1', &
            'Total precipitation flux','precipitation_flux',.true., &
            'time: mean')
        if ( enable_srf2d_vars(srf_evp) ) &
          call setup_var(v2dvar_srf(srf_evp),vsize,'evp','kg m-2 day-1', &
            'Total evapotranspiration flux','water_evaporation_flux',.true., &
            'time: mean')
        if ( enable_srf2d_vars(srf_scv) ) &
          call setup_var(v2dvar_srf(srf_scv),vsize,'scv','kg m-2 day-1', &
            'Liquid water equivalent of snow thickness', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean')
        if ( enable_srf2d_vars(srf_sena) ) &
          call setup_var(v2dvar_srf(srf_sena),vsize,'sena','W m-2', &
            'Sensible heat flux','surface_downward_sensible_heat_flux', &
            .true.,'time: mean')
        if ( enable_srf2d_vars(srf_flw) ) &
          call setup_var(v2dvar_srf(srf_flw),vsize,'flw','W m-2', &
            'Net upward longwave energy flux', &
            'net_upward_longwave_flux_in_air',.true.,'time: mean')
        if ( enable_srf2d_vars(srf_fsw) ) &
          call setup_var(v2dvar_srf(srf_fsw),vsize,'fsw','W m-2', &
            'Net downward shortwave energy flux', &
            'net_downward_shortwave_flux_in_air',.true.,'time: mean')
        if ( enable_srf2d_vars(srf_fld) ) &
          call setup_var(v2dvar_srf(srf_fld),vsize,'fld','W m-2', &
            'Surface downward longwave flux in air', &
            'surface_downwelling_longwave_flux_in_air',.true.,'time: mean')
        if ( enable_srf2d_vars(srf_sina) ) &
          call setup_var(v2dvar_srf(srf_sina),vsize,'sina','W m-2', &
            'Surface downward shortwave flux in air', &
            'surface_downwelling_shortwave_flux_in_air',.true.,'time: mean')
        if ( enable_srf2d_vars(srf_prcv) ) &
          call setup_var(v2dvar_srf(srf_prcv),vsize,'prcv','kg m-2 day-1', &
            'Convective precipitation flux','convective_rainfall_flux', &
            .true.,'time: mean')
        if ( enable_srf2d_vars(srf_zpbl) ) &
          call setup_var(v2dvar_srf(srf_zpbl),vsize,'zpbl','m', &
            'PBL thickness','atmosphere_boundary_layer_thickness',.true.)
        if ( enable_srf2d_vars(srf_aldirs) ) &
          call setup_var(v2dvar_srf(srf_aldirs),vsize,'aldirs','1', &
            'Surface albedo to direct shortwave radiation', &
            'surface_albedo_short_wave_direct',.true.)
        if ( enable_srf2d_vars(srf_aldifs) ) &
          call setup_var(v2dvar_srf(srf_aldifs),vsize,'aldifs','1', &
            'Surface albedo to diffuse shortwave radiation', &
            'surface_albedo_short_wave_diffuse',.true.)
        if ( enable_srf2d_vars(srf_sund) ) &
          call setup_var(v2dvar_srf(srf_sund),vsize,'sund','s', &
            'Duration of sunshine','duration_of_sunshine',.true.,'time: sum')

        if ( iseaice == 1 ) then
          if ( enable_srf2d_vars(srf_seaice) ) &
            call setup_var(v2dvar_srf(srf_seaice),vsize,'seaice','1', &
              'Sea ice mask','seaice_binary_mask',.true.)
        else
          enable_srf2d_vars(srf_seaice) = .false.
        end if

        vsize%k2 = 1
        v3dvar_srf(srf_u10m)%axis = 'xyw'
        v3dvar_srf(srf_v10m)%axis = 'xyw'
        v3dvar_srf(srf_t2m)%axis = 'xy2'
        v3dvar_srf(srf_q2m)%axis = 'xy2'

        if ( enable_srf3d_vars(srf_u10m) ) &
          call setup_var(v3dvar_srf(srf_u10m),vsize,'u10m','m s-1', &
            '10 meter zonal (westerly) wind component', &
            'eastward_wind',.true.)
        if ( enable_srf3d_vars(srf_v10m) ) &
          call setup_var(v3dvar_srf(srf_v10m),vsize,'v10m','m s-1', &
            '10 meter meridional (southerly) wind component' , &
            'northward_wind',.true.)
        if ( enable_srf3d_vars(srf_t2m) ) &
          call setup_var(v3dvar_srf(srf_t2m),vsize,'t2m','K', &
            '2 meter air temperature','air_temperature',.true.)
        if ( enable_srf3d_vars(srf_q2m) ) &
          call setup_var(v3dvar_srf(srf_q2m),vsize,'q2m','1', &
            '2 meter air specific humidity','specific_humidity',.true.)

        vsize%k2 = 2
        v3dvar_srf(srf_smw)%axis = 'xys'
        v3dvar_srf(srf_runoff)%axis = 'xys'
        if ( enable_srf3d_vars(srf_smw) ) &
          call setup_var(v3dvar_srf(srf_smw),vsize,'smw','kg m-2', &
            'Moisture content of the soil layers', &
            'soil_moisture_content',.true.)
        if ( enable_srf3d_vars(srf_runoff) ) &
          call setup_var(v3dvar_srf(srf_runoff),vsize,'runoff','kg m-2 day-1', &
            'Runoff flux','runoff_flux',.true.,'time: mean')

        enable_srf_vars(1:nsrf2dvars) = enable_srf2d_vars
        enable_srf_vars(nsrf2dvars+1:nsrfvars) = enable_srf3d_vars
        outstream(srf_stream)%nvar = countvars(enable_srf_vars,nsrfvars)
        allocate(outstream(srf_stream)%ncvars%vlist(outstream(srf_stream)%nvar))

        vcount = 1
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
        outstream(srf_stream)%jl1 = vsize%j1
        outstream(srf_stream)%jl2 = vsize%j2
        outstream(srf_stream)%il1 = vsize%i1
        outstream(srf_stream)%il2 = vsize%i2
        outstream(srf_stream)%jg1 = jout1
        outstream(srf_stream)%jg2 = jout2
        outstream(srf_stream)%ig1 = iout1
        outstream(srf_stream)%ig2 = iout2
      end if

      if ( nstream == sts_stream ) then

        allocate(v2dvar_sts(nsts2dvars))
        allocate(v3dvar_sts(nsts3dvars))
        enable_sts2d_vars = enable_sts_vars(1:nsts2dvars)
        enable_sts3d_vars = enable_sts_vars(nsts2dvars+1:nstsvars)

        ! This variables are always present

        call setup_var(v2dvar_sts(sts_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_sts(sts_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_sts(sts_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_sts(sts_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_sts(sts_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_sts2d_vars(sts_tgmax) ) &
          call setup_var(v2dvar_sts(sts_tgmax),vsize,'tgmax','K', &
            'Maximum surface temperature','surface_temperature', &
            .true.,'time: maximum')
        if ( enable_sts2d_vars(sts_tgmin) ) &
          call setup_var(v2dvar_sts(sts_tgmin),vsize,'tgmin','K', &
            'Minimum surface temperature','surface_temperature', &
            .true.,'time: minimum')
        if ( enable_sts2d_vars(sts_pcpmax) ) &
          call setup_var(v2dvar_sts(sts_pcpmax),vsize,'pcpmax','kg m-2 s-1', &
            'Maximum total precipitation flux','precipitation_flux', &
            .true.,'time: maximum')
        if ( enable_sts2d_vars(sts_pcpavg) ) &
          call setup_var(v2dvar_sts(sts_pcpavg),vsize,'pcpavg','kg m-2 s-1', &
            'Mean total precipitation flux','precipitation_flux', &
            .true.,'time: mean')
        if ( enable_sts2d_vars(sts_sund) ) &
          call setup_var(v2dvar_sts(sts_sund),vsize,'sund','s', &
            'Duration of sunshine','duration_of_sunshine',.true.,'time: sum')
        if ( enable_sts2d_vars(sts_psmin) ) &
          call setup_var(v2dvar_sts(sts_psmin),vsize,'psmin','hPa', &
            'Minimum of surface pressure','air_pressure',.true.,'time: minimum')

        vsize%k2 = 1
        v3dvar_sts(sts_t2max)%axis = 'xy2'
        v3dvar_sts(sts_t2min)%axis = 'xy2'
        v3dvar_sts(sts_t2avg)%axis = 'xy2'
        v3dvar_sts(sts_w10max)%axis = 'xyw'
        if ( enable_sts3d_vars(sts_t2max) ) &
          call setup_var(v3dvar_sts(sts_t2max),vsize,'t2max','K', &
            'Maximum 2 meter temperature','air_temperature',.true., &
            'time: maximum')
        if ( enable_sts3d_vars(sts_t2min) ) &
          call setup_var(v3dvar_sts(sts_t2min),vsize,'t2min','K', &
            'Minimum 2 meter temperature','air_temperature',.true., &
            'time: minimum')
        if ( enable_sts3d_vars(sts_t2avg) ) &
          call setup_var(v3dvar_sts(sts_t2avg),vsize,'t2avg','K', &
            'Mean 2 meter temperature','air_temperature',.true., &
            'time: mean')
        if ( enable_sts3d_vars(sts_w10max) ) &
          call setup_var(v3dvar_sts(sts_w10max),vsize,'w10max','m s-1', &
            'Maximum speed of 10m wind','wind_speed',.true.,'time: maximum')

        enable_sts_vars(1:nsts2dvars) = enable_sts2d_vars
        enable_sts_vars(nsts2dvars+1:nstsvars) = enable_sts3d_vars
        outstream(sts_stream)%nvar = countvars(enable_sts_vars,nstsvars)
        allocate(outstream(sts_stream)%ncvars%vlist(outstream(sts_stream)%nvar))

        vcount = 1
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
        outstream(sts_stream)%jl1 = vsize%j1
        outstream(sts_stream)%jl2 = vsize%j2
        outstream(sts_stream)%il1 = vsize%i1
        outstream(sts_stream)%il2 = vsize%i2
        outstream(sts_stream)%jg1 = jout1
        outstream(sts_stream)%jg2 = jout2
        outstream(sts_stream)%ig1 = iout1
        outstream(sts_stream)%ig2 = iout2
      end if

      if ( nstream == sub_stream ) then

        allocate(v2dvar_sub(nsub2dvars))
        allocate(v3dvar_sub(nsub3dvars))
        enable_sub2d_vars = enable_sub_vars(1:nsub2dvars)
        enable_sub3d_vars = enable_sub_vars(nsub2dvars+1:nsubvars)

        vsize%j1 = (vsize%j1-1)*nsg
        vsize%j2 = vsize%j2*nsg
        vsize%i1 = (vsize%i1-1)*nsg
        vsize%i2 = vsize%i2*nsg

        ! This variables are always present

        call setup_var(v2dvar_sub(sub_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_sub(sub_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_sub(sub_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_sub(sub_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_sub(sub_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_sub2d_vars(sub_uvdrag) ) &
          call setup_var(v2dvar_sub(sub_uvdrag),vsize,'uvdrag','1', &
            'Surface drag stress coefficient in air', &
            'surface_drag_coefficient_in_air',.true.)
        if ( enable_sub2d_vars(sub_tg) ) &
          call setup_var(v2dvar_sub(sub_tg),vsize,'tg','K', &
            'Ground temperature','surface_temperature',.true.)
        if ( enable_sub2d_vars(sub_tlef) ) &
          call setup_var(v2dvar_sub(sub_tlef),vsize,'tlef','K', &
            'Foliage canopy temperature','canopy_temperature',.true.)
        if ( enable_sub2d_vars(sub_evp) ) &
          call setup_var(v2dvar_sub(sub_evp),vsize,'evp','kg m-2 day-1', &
            'Total evapotranspiration flux','water_evaporation_flux',.true., &
            'time: mean')
        if ( enable_sub2d_vars(sub_scv) ) &
          call setup_var(v2dvar_sub(sub_scv),vsize,'scv','kg m-2', &
            'Liquid water equivalent of snow thickness', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean')
        if ( enable_sub2d_vars(sub_sena) ) &
          call setup_var(v2dvar_sub(sub_sena),vsize,'sena','W m-2', &
            'Sensible heat flux','surface_downward_sensible_heat_flux', &
            .true.,'time: mean')
        if ( lakemod == 1 ) then
          if ( enable_sub2d_vars(sub_tlake) ) &
            call setup_var(v2dvar_sub(sub_tlake),vsize,'tlake','K', &
              'Lake water surface temperature','water_temperature',.true.)
        else
          enable_sub2d_vars(sub_tlake) = .false.
        end if

        vsize%k2 = 1
        v3dvar_sub(sub_u10m)%axis = 'xyw'
        v3dvar_sub(sub_v10m)%axis = 'xyw'
        v3dvar_sub(sub_t2m)%axis = 'xy2'
        v3dvar_sub(sub_q2m)%axis = 'xy2'
        if ( enable_sub3d_vars(sub_u10m) ) &
          call setup_var(v3dvar_sts(sub_u10m),vsize,'u10m','m s-1', &
            '10 meter zonal wind component (westerly)', &
            'eastward_wind',.true.)
        if ( enable_sub3d_vars(sub_v10m) ) &
          call setup_var(v3dvar_sts(sub_v10m),vsize,'v10m','m s-1', &
            '10 meter meridional wind component (southerly)', &
            'northward_wind',.true.)
        if ( enable_sub3d_vars(sub_t2m) ) &
          call setup_var(v3dvar_sts(sub_t2m),vsize,'t2m','K', &
            '2 meter air temperature','air_temperature',.true.)
        if ( enable_sub3d_vars(sub_q2m) ) &
          call setup_var(v3dvar_sts(sub_q2m),vsize,'q2m','1', &
            '2 meter air specific humidity','specific_humidity',.true.)
        vsize%k2 = 2
        v3dvar_sub(sub_smw)%axis = 'xys'
        v3dvar_sub(sub_runoff)%axis = 'xys'
        if ( enable_sub3d_vars(sub_smw) ) &
          call setup_var(v3dvar_sts(sub_smw),vsize,'smw','kg kg-1', &
            'Soil moisture content','soil_moisture_content',.true.)
        if ( enable_sub3d_vars(sub_runoff) ) &
          call setup_var(v3dvar_sts(sub_runoff),vsize,'runoff','kg m-2 day-1', &
            'Runoff flux','runoff_flux',.true.,'time: mean')

        enable_sub_vars(1:nsub2dvars) = enable_sub2d_vars
        enable_sub_vars(nsub2dvars+1:nsubvars) = enable_sub3d_vars
        outstream(sub_stream)%nvar = countvars(enable_sub_vars,nsubvars)
        allocate(outstream(sub_stream)%ncvars%vlist(outstream(sub_stream)%nvar))

        vcount = 1
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
        outstream(sub_stream)%jl1 = vsize%j1
        outstream(sub_stream)%jl2 = vsize%j2
        outstream(sub_stream)%il1 = vsize%i1
        outstream(sub_stream)%il2 = vsize%i2
        outstream(sub_stream)%jg1 = joutsg1
        outstream(sub_stream)%jg2 = joutsg2
        outstream(sub_stream)%ig1 = ioutsg1
        outstream(sub_stream)%ig2 = ioutsg2
      end if

      if ( nstream == rad_stream ) then

        allocate(v2dvar_rad(nrad2dvars))
        allocate(v3dvar_rad(nrad3dvars))
        enable_rad2d_vars = enable_rad_vars(1:nrad2dvars)
        enable_rad3d_vars = enable_rad_vars(nrad2dvars+1:nradvars)

        ! This variables are always present

        call setup_var(v2dvar_rad(rad_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_rad(rad_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_rad(rad_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_rad(rad_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_rad(rad_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_rad2d_vars(rad_frsa) ) &
          call setup_var(v2dvar_rad(rad_frsa),vsize,'frsa','W m-2', &
            'Surface net downward shortwave flux', &
            'surface_net_downward_shortwave_flux', .true.)
        if ( enable_rad2d_vars(rad_frla) ) &
          call setup_var(v2dvar_rad(rad_frla),vsize,'frla','W m-2', &
            'Surface net upward longwave flux', &
            'surface_net_upward_longwave_flux', .true.)
        if ( enable_rad2d_vars(rad_clrst) ) &
          call setup_var(v2dvar_rad(rad_clrst),vsize,'clrst','W m-2', &
            'Clearsky top of atmosphere net downward shortwave flux', &
            'toa_net_downward_shortwave_flux_assuming_clear_sky',.true.)
        if ( enable_rad2d_vars(rad_clrss) ) &
          call setup_var(v2dvar_rad(rad_clrss),vsize,'clrss','W m-2', &
            'Clearsky surface net downward shortwave flux', &
            'surface_net_downward_shortwave_flux_assuming_clear_sky',.true.)
        if ( enable_rad2d_vars(rad_clrlt) ) &
          call setup_var(v2dvar_rad(rad_clrlt),vsize,'clrlt','W m-2', &
            'Clearsky top of atmosphere net upward longwave flux', &
            'toa_net_upward_longwave_flux_assuming_clear_sky',.true.)
        if ( enable_rad2d_vars(rad_clrls) ) &
          call setup_var(v2dvar_rad(rad_clrls),vsize,'clrls','W m-2', &
            'Clearsky net upward longwave flux', &
            'surface_net_upward_longwave_flux_assuming_clear_sky',.true.)
        if ( enable_rad2d_vars(rad_solin) ) &
          call setup_var(v2dvar_rad(rad_solin),vsize,'solin','W m-2', &
            'Top of atmosphere incoming shortwave flux', &
            'toa_incoming_shortwave_flux',.true.)
        if ( enable_rad2d_vars(rad_sabtp) ) &
          call setup_var(v2dvar_rad(rad_sabtp),vsize,'sabtp','W m-2', &
            'Net top of atmosphere upward shortwave flux', &
            'toa_net_upward_shortwave_flux',.true.)
        if ( enable_rad2d_vars(rad_totcf) ) &
          call setup_var(v2dvar_rad(rad_totcf),vsize,'totcf','1', &
            'Total cloud fraction','cloud_area_fraction',.true.)
        if ( enable_rad2d_vars(rad_totcl) ) &
          call setup_var(v2dvar_rad(rad_totcl),vsize,'totcl','kg m-2', &
            'Total columnar liquid water content', &
            'atmosphere_cloud_condensed_water_content',.true.)
        if ( enable_rad2d_vars(rad_totci) ) &
          call setup_var(v2dvar_rad(rad_totci),vsize,'totci','kg m-2', &
            'Total columnar ice water content', &
            'atmosphere_ice_condensed_water_content',.true.)
        if ( enable_rad2d_vars(rad_firtp) ) &
          call setup_var(v2dvar_rad(rad_firtp),vsize,'firtp','W m-2', &
            'Top of atmosphere net upward longwave flux', &
            'toa_net_upward_longwave_flux',.true.)

        vsize%k2 = kz
        if ( enable_rad3d_vars(rad_cld) ) &
          call setup_var(v3dvar_rad(rad_cld),vsize,'cld','1', &
            'Cloud fractional cover', &
            'cloud_area_fraction_in_atmosphere_layer',.true.)
        if ( enable_rad3d_vars(rad_clwp) ) &
          call setup_var(v3dvar_rad(rad_clwp),vsize,'clwp','g m-2', &
            'Cloud liquid water path','thickness_of_liquid_water_cloud',.true.)
        if ( enable_rad3d_vars(rad_qrs) ) &
          call setup_var(v3dvar_rad(rad_qrs),vsize,'qrs','K s-1', &
            'Shortwave radiation heating rate', &
            'tendency_of_air_temperature_due_to_shortwave_heating',.true.)
        if ( enable_rad3d_vars(rad_qrl) ) &
          call setup_var(v3dvar_rad(rad_qrl),vsize,'qrl','K s-1', &
            'Longwave radiation heating rate', &
            'tendency_of_air_temperature_due_to_longwave_heating',.true.)

        enable_rad_vars(1:nrad2dvars) = enable_rad2d_vars
        enable_rad_vars(nrad2dvars+1:nradvars) = enable_rad3d_vars
        outstream(rad_stream)%nvar = countvars(enable_rad_vars,nradvars)
        allocate(outstream(rad_stream)%ncvars%vlist(outstream(rad_stream)%nvar))

        vcount = 1
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
        outstream(rad_stream)%jl1 = vsize%j1
        outstream(rad_stream)%jl2 = vsize%j2
        outstream(rad_stream)%il1 = vsize%i1
        outstream(rad_stream)%il2 = vsize%i2
        outstream(rad_stream)%jg1 = jout1
        outstream(rad_stream)%jg2 = jout2
        outstream(rad_stream)%ig1 = iout1
        outstream(rad_stream)%ig2 = iout2
      end if

      if ( nstream == lak_stream ) then

        allocate(v2dvar_lak(nlak2dvars))
        allocate(v3dvar_lak(nlak3dvars))
        enable_lak2d_vars = enable_lak_vars(1:nlak2dvars)
        enable_lak3d_vars = enable_lak_vars(nlak2dvars+1:nlakvars)

        ! This variables are always present

        call setup_var(v2dvar_lak(lak_xlon),vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_lak(lak_xlat),vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_lak(lak_mask),vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_lak(lak_topo),vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_lak(lak_ps),vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.)

        ! The following may be enabled/disabled

        if ( enable_lak2d_vars(lak_tg) ) &
          call setup_var(v2dvar_lak(lak_tg),vsize,'tg','K', &
            'Ground temperature','surface_temperature',.true.)
        if ( enable_lak2d_vars(lak_tpr) ) &
          call setup_var(v2dvar_lak(lak_tpr),vsize,'tpr','kg m-2 day-1', &
            'Total precipitation flux','precipitation_flux',.true.,'time: mean')
        if ( enable_lak2d_vars(lak_scv) ) &
          call setup_var(v2dvar_lak(lak_scv),vsize,'scv','kg m-2', &
            'Liquid water equivalent snow depth', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean')
        if ( enable_lak2d_vars(lak_sena) ) &
          call setup_var(v2dvar_lak(lak_sena),vsize,'sena','W m-2', &
            'Sensible heat flux','surface_downward_sensible_heat_flux',.true.)
        if ( enable_lak2d_vars(lak_flw) ) &
          call setup_var(v2dvar_lak(lak_flw),vsize,'flw','W m-2', &
            'Net longwave energy flux','net_upward_longwave_flux_in_air',.true.)
        if ( enable_lak2d_vars(lak_fsw) ) &
          call setup_var(v2dvar_lak(lak_fsw),vsize,'fsw','W m-2', &
            'Net shortwave energy flux', 'net_downward_shortwave_flux_in_air', &
            .true.)
        if ( enable_lak2d_vars(lak_fld) ) &
          call setup_var(v2dvar_lak(lak_fld),vsize,'fld','W m-2', &
            'Downward longwave flux at surface in air', &
            'surface_downwelling_longwave_flux_in_air',.true.)
        if ( enable_lak2d_vars(lak_sina) ) &
          call setup_var(v2dvar_lak(lak_sina),vsize,'sina','W m-2', &
            'Downward shortwave flux at surface in air', &
            'surface_downwelling_shortwave_flux_in_air',.true.)
        if ( enable_lak2d_vars(lak_aldirs) ) &
          call setup_var(v2dvar_lak(lak_aldirs),vsize,'aldirs','1', &
            'Surface albedo to direct shortwave radiation', &
            'surface_albedo_short_wave_direct',.true.)
        if ( enable_lak2d_vars(lak_aldifs) ) &
          call setup_var(v2dvar_lak(lak_aldifs),vsize,'aldifs','1', &
            'Surface albedo to diffuse shortwave radiation', &
            'surface_albedo_short_wave_diffuse',.true.)
        if ( enable_lak2d_vars(lak_evp) ) &
          call setup_var(v2dvar_lak(lak_evp),vsize,'evp','mm s-1', &
            'Water evaporation flux', &
            'water_evaporation_flux_where_sea_ice',.true.)
        if ( enable_lak2d_vars(lak_aveice) ) &
          call setup_var(v2dvar_lak(lak_aveice),vsize,'aveice','mm', &
            'Floating ice thickness','floating_ice_thickness',.true.)
        if ( enable_lak2d_vars(lak_hsnow) ) &
          call setup_var(v2dvar_lak(lak_hsnow),vsize,'hsnow','mm', &
            'Floating snow thickness','surface_snow_thickness_where_sea_ice', &
            .true.)

        vsize%k2 = ndpmax
        v3dvar_lak(lak_tlake)%axis = 'xyd'
        if ( enable_lak3d_vars(lak_tlake) ) &
          call setup_var(v3dvar_lak(lak_tlake),vsize,'tlake','K', &
            'Lake water temperature','water_temperature',.true.)

        enable_lak_vars(1:nlak2dvars) = enable_lak2d_vars
        enable_lak_vars(nlak2dvars+1:nlakvars) = enable_lak3d_vars
        outstream(lak_stream)%nvar = countvars(enable_lak_vars,nlakvars)
        allocate(outstream(lak_stream)%ncvars%vlist(outstream(lak_stream)%nvar))

        vcount = 1
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
        outstream(lak_stream)%jl1 = vsize%j1
        outstream(lak_stream)%jl2 = vsize%j2
        outstream(lak_stream)%il1 = vsize%i1
        outstream(lak_stream)%il2 = vsize%i2
        outstream(lak_stream)%jg1 = jout1
        outstream(lak_stream)%jg2 = jout2
        outstream(lak_stream)%ig1 = iout1
        outstream(lak_stream)%ig2 = iout2
      end if

      outstream(nstream)%opar%pname = 'RegCM Model'
      outstream(nstream)%opar%zero_date = idate1
      outstream(nstream)%opar%l_sync = lsync

      if ( lparallel ) then
        outstream(nstream)%opar%mpi_comm = mycomm
        outstream(nstream)%opar%mpi_info = mpi_info_null
#ifdef NETCDF4_HDF5
        outstream(nstream)%opar%mpi_iotype = nf90_mpiio
#endif
        ! The "global" indexes in the output stream refer to the INTERNAL
        ! CROSS grid, i.e. for processor 0 this is (2,2) => (1,1) so we must
        ! subtract 1 line/column to rebase on pixel (2,2) of the internal model
        ! cross points grid to point (1,1).
        outstream(nstream)%opar%global_jstart = global_cross_jstart - 1
        outstream(nstream)%opar%global_jend   = global_cross_jend   - 1
        outstream(nstream)%opar%global_istart = global_cross_istart - 1
        outstream(nstream)%opar%global_iend   = global_cross_iend   - 1
        if ( nstream == sub_stream ) then
          outstream(nstream)%opar%global_jstart = &
            (outstream(nstream)%opar%global_jstart-1)*nsg+1
          outstream(nstream)%opar%global_jend =   &
            outstream(nstream)%opar%global_jend*nsg
          outstream(nstream)%opar%global_istart = &
            (outstream(nstream)%opar%global_istart-1)*nsg+1
          outstream(nstream)%opar%global_iend =   &
            outstream(nstream)%opar%global_iend*nsg
        end if
        outstream(nstream)%opar%l_subgrid = .true.
      else
        ! Allocate space to collect all from all CPUs
        if ( myid == iocpu ) then
          kkz = 2
          if ( enable_flag(atm_stream) .or. enable_flag(rad_stream) ) then
            kkz = max(kz,kkz)
          end if
          if ( lakemod == 1 .and. enable_flag(lak_stream) ) then
            kkz = max(ndpmax,kkz)
          end if
          call getmem2d(io2d,jout1,jout2,iout1,iout2,'ncout:io2d')
          if ( kkz > 0 ) then
            call getmem3d(io3d,jout1,jout2,iout1,iout2,1,kkz,'ncout:io3d')
          end if
          if ( enable_flag(sub_stream) ) then
            call getmem2d(io2dsg,joutsg1,joutsg2,ioutsg1,ioutsg2,'ncout:io2dsg')
          end if
        end if
      end if
    end do
  end subroutine init_output_streams

  integer(ik4) function countvars(eflags,ntot)
    implicit none
    integer(ik4) , intent(in) :: ntot
    logical , dimension(ntot) , intent(in) :: eflags
    integer(ik4) :: i
    countvars = nbase
    do i = nbase , ntot
      if ( eflags(i) ) countvars = countvars + 1
    end do
  end function countvars

  subroutine newoutfile(idate,itype)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=16) :: fbname
    character(len=36) :: cdate
    integer(ik4) , intent(in) :: itype

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
    cdate = tochar(idate0)
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_initial_start',cdate))
    cdate = tochar(idate1)
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_start',cdate))
    cdate = tochar(idate1)
    call outstream_addatt(outstream(itype)%ncout, &
      ncattribute_string('model_simulation_end',cdate))

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

  subroutine setup_var(var,vsize,vname,vunit,long_name,standard_name, &
                       l_rec,cell_method,l_fill,rmissval)
    implicit none
    class(ncvariable_standard) , intent(inout) :: var
    type(varspan) , intent(in) :: vsize
    character(len=*) , intent(in) :: vname , vunit , long_name , standard_name
    character(len=*) , intent(in) , optional :: cell_method
    logical , intent(in) , optional :: l_rec , l_fill
    real(rk4) , intent(in) , optional :: rmissval
    var%vname = vname
    var%vunit = vunit
    var%long_name = long_name
    var%standard_name = standard_name
    if ( present(cell_method) ) var%cell_method = cell_method
    if ( present(rmissval) .or. present(l_fill) ) then
      var%lfillvalue = .true.
      if ( present(rmissval) ) var%rmissval = rmissval
    end if
    select type(var)
      type is (ncvariable2d_real)
        if ( present(l_rec) ) var%lrecords = l_rec
        call getmem2d(var%rval,vsize%j1,vsize%j2, &
          vsize%i1,vsize%i2,'ncout:setup_var:'//trim(var%vname))
        var%j1 = vsize%j1
        var%j2 = vsize%j2
        var%i1 = vsize%i1
        var%i2 = vsize%i2
      type is (ncvariable3d_real)
        if ( present(l_rec) ) var%lrecords = l_rec
        call getmem3d(var%rval,vsize%j1,vsize%j2,vsize%i1,vsize%i2, &
          vsize%k1,vsize%k2,'ncout:setup_var:'//trim(var%vname))
        var%j1 = vsize%j1
        var%j2 = vsize%j2
        var%i1 = vsize%i1
        var%i2 = vsize%i2
        var%k1 = vsize%k1
        var%k2 = vsize%k2
    end select
  end subroutine setup_var

  subroutine dispose_output_streams
    implicit none
    integer(ik4) :: nstream
    if ( associated(v2dvar_atm) ) deallocate(v2dvar_atm)
    if ( associated(v3dvar_atm) ) deallocate(v3dvar_atm)
    if ( associated(v2dvar_srf) ) deallocate(v2dvar_srf)
    if ( associated(v3dvar_srf) ) deallocate(v3dvar_srf)
    if ( associated(v2dvar_sts) ) deallocate(v2dvar_sts)
    if ( associated(v3dvar_sts) ) deallocate(v3dvar_sts)
    if ( associated(v2dvar_rad) ) deallocate(v2dvar_rad)
    if ( associated(v3dvar_rad) ) deallocate(v3dvar_rad)
    if ( associated(v2dvar_sub) ) deallocate(v2dvar_sub)
    if ( associated(v3dvar_sub) ) deallocate(v3dvar_sub)
    if ( associated(v2dvar_lak) ) deallocate(v2dvar_lak)
    if ( associated(v3dvar_lak) ) deallocate(v3dvar_lak)
    do nstream = 1 , maxstreams
      call outstream_dispose(outstream(nstream)%ncout)
      if ( associated(outstream(nstream)%ncvars%vlist) ) then
        deallocate(outstream(nstream)%ncvars%vlist)
      end if
    end do
  end subroutine dispose_output_streams

  subroutine write_record_output_stream(lparallel,istream,idate)
    implicit none
    logical , intent(in) :: lparallel
    integer(ik4) , intent(in) :: istream
    type(rcm_time_and_date) , intent(in) :: idate
    real(rk8) , pointer , dimension(:,:) :: tmp2d
    real(rk8) , pointer , dimension(:,:,:) :: tmp3d
    class(ncvariable_standard) , pointer :: vp
    integer(ik4) :: ivar

    call outstream_addrec(outstream(istream)%ncout,idate)
    do ivar = 1 , outstream(istream)%nvar
      vp => outstream(istream)%ncvars%vlist(ivar)%vp
      if ( .not. lparallel ) then
        select type(vp)
          type is (ncvariable2d_real)
            call grid_collect(io2d,vp%rval,vp%j1,vp%j2,vp%i1,vp%i2)
            vp%j1 = outstream(istream)%jg1
            vp%j2 = outstream(istream)%jg2
            vp%i1 = outstream(istream)%ig1
            vp%i2 = outstream(istream)%ig2
            tmp2d => vp%rval
            vp%rval => io2d
          type is (ncvariable3d_real)
            call grid_collect(io3d,vp%rval,vp%j1,vp%j2,vp%i1,vp%i2, &
              1,size(vp%rval,3))
            vp%j1 = outstream(istream)%jg1
            vp%j2 = outstream(istream)%jg2
            vp%i1 = outstream(istream)%ig1
            vp%i2 = outstream(istream)%ig2
            tmp3d => vp%rval
            vp%rval => io3d
        end select
      end if
      call outstream_writevar(outstream(istream)%ncout,vp)
      if ( .not. lparallel ) then
        select type(vp)
          type is (ncvariable2d_real)
            vp%rval => tmp2d
            vp%j1 = outstream(istream)%jl1
            vp%j2 = outstream(istream)%jl2
            vp%i1 = outstream(istream)%il1
            vp%i2 = outstream(istream)%il2
          type is (ncvariable3d_real)
            vp%rval => tmp3d
            vp%j1 = outstream(istream)%jl1
            vp%j2 = outstream(istream)%jl2
            vp%i1 = outstream(istream)%il1
            vp%i2 = outstream(istream)%il2
        end select
      end if
    end do
  end subroutine write_record_output_stream

end module mod_ncout
