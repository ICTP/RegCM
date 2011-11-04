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

module mod_runparams

  use mod_constants
  use mod_dynparam
  use mod_mpmessage
  use mod_service 
  use mod_memutil

  implicit none
 
  type(rcm_time_and_date) , save :: idate0 , idate1 , idate2

  type(rcm_time_and_date) , save :: idatex
  integer :: xyear , xmonth , xday , xhour

  type(rcm_time_and_date) , save :: bdydate1 , bdydate2

  type(rcm_time_interval) , save :: intmdl
  type(rcm_time_interval) , save :: intbdy

  real(8) :: declin , deltmx
  real(8) :: xbctime
  real(8) :: calday , twodt

  ! Step counter. Is zero at idate0, always increasing, never reset.
  integer(8) :: ktau
  ! Final number of step for THIS run
  integer(8) :: mtau
  ! How many steps for an hour (updates date fields Y m d H)
  integer(8) :: khour
  ! Output k values for I/O operations.
  integer(8) :: katm , krad , kche , ksav , kdbg , kbdy , ksrf
  ! Seconds counter in between boundary conditions read
  integer(8) :: nbdytime
  ! Step counters to activate surface and radiation schemes
  integer(8) :: ntsrf , ntrad
  ! Model timestep in seconds (real and integer)
  integer(8) :: ntsec
  real(8) :: dtsec
  ! Internal count for how many SRF outputs every LAK output
  integer :: klak
!
  real(8) :: dt , dt2 , dtbdys
  real(8) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(8) :: c200 , rdxsq , dtsrf , dtabem , dtrad, cpldt
  real(8) :: fnudge , gnudge
  real(8) :: xkhmax , xkhz

  integer :: iboudy , ichem , ipgf , ipptls, cplexvars, cplinterp

  logical :: ifrest , rfstrt , doing_restart, cplbdysmooth 

  integer :: ispgd , ispgx , kchi , kclo , kcmd, cpldbglevel
!
  real(8) :: akht1 , akht2

  real(8) , pointer , dimension(:) :: dtau
  real(8) , pointer , dimension(:) :: a , anudg , dsigma , qcon
  real(8) , pointer , dimension(:) :: sigma
  real(8) , pointer , dimension(:,:) :: twt
  real(8) , pointer , dimension(:) :: wgtd
  real(8) , pointer , dimension(:) :: wgtx

  character(len=3) :: scenario

  integer , parameter :: n_atmvar = 15
  integer , parameter :: n_srfvar = 31
  integer , parameter :: n_subvar = 16
  integer , parameter :: n_radvar = 15
  integer , parameter :: n_chevar = 17
  integer , parameter :: n_lakvar = 16

  integer, private  :: ierr 
  real(8) , private :: total_allocation_size

  type output_variable
    character(len=8) :: vname
    character(len=128) :: vstd_name
    character(len=128) :: vdesc
    character(len=16) :: vunit
    logical :: enabled
  end type

  type(output_variable) , dimension(n_atmvar) :: atm_variables
  type(output_variable) , dimension(n_srfvar) :: srf_variables
  type(output_variable) , dimension(n_subvar) :: sub_variables
  type(output_variable) , dimension(n_radvar) :: rad_variables
  type(output_variable) , dimension(n_chevar) :: che_variables
  type(output_variable) , dimension(n_lakvar) :: lak_variables

  data total_allocation_size /d_zero/
  data doing_restart /.false./

  data atm_variables / &
    output_variable('time','','','',.true.),                                  &
    output_variable('ps','','','',.true.),                                    &
    output_variable('u','eastward_wind','U component (westerly) of wind',     &
                    'm s-1',.true.),                                          &
    output_variable('v','northward_wind','V component (southerly) of wind',   &
                    'm s-1',.true.),                                          &
    output_variable('omega','lagrangian_tendency_of_air_pressure',            &
                    'Pressure velocity','hPa s-1',.true.),                    &
    output_variable('t','air_temperature','Temperature','K',.true.),          &
    output_variable('qv','humidity_mixing_ratio',                             &
                    'Water vapor mixing ratio','kg kg-1',.true.),             &
    output_variable('qc','cloud_liquid_water_mixing_ratio',                   &
                    'Cloud water mixing ratio','kg kg-1',.true.),             &
    output_variable('tke','turbulent_kinetic_energy',                         &
                    'Turbulent Kinetic Energy','m2 s2',.true.),               &
    output_variable('kth','vertical_momentum_diffusivity',                    &
                    'Vertical Turbulent Viscosity','m2 s-1',.true.),          &
    output_variable('kzm','vertical_scalar_diffusivity',                      &
                    'Vertical Turbulent Diffusivity','m2 s-1',.true.),        &
    output_variable('tpr','precipitation_flux',                               &
                    'Total daily precipitation rate','kg m-2 day-1',.true.),  &
    output_variable('tgb','soil_temperature',                                 &
                    'Lower groud temperature','K',.true.),                    &
    output_variable('swt','moisture_content_of_soil_layer',                   &
                    'Total soil water','kg m-2',.true.),                      &
    output_variable('rno','runoff_flux',                                      &
                    'Runoff accumulated infiltration','kg m-2 day-1',.true.) /

  data srf_variables / &
    output_variable('time','','','',.true.),                                      &
    output_variable('tbnds','','','',.true.),                                     &
    output_variable('ps','','','',.true.),                                        &
    output_variable('u10m','eastward_wind',                                       &
                    '10 meters U component (westerly) of wind','m s-1',.true.),   &
    output_variable('v10m','northward_wind',                                      &
                    '10 meters V component (southerly) of wind','m s-1',.true.),  &
    output_variable('uvdrag','surface_drag_coefficient_in_air',                   &
                    'Surface drag stress','1',.true.),                            &
    output_variable('tg','surface_temperature',                                   &
                    'Ground temperature','K',.true.),                             &
    output_variable('tlef','canopy_temperature',                                  &
                    'Foliage temperature','K',.true.),                            &
    output_variable('t2m','air_temperature','2 meters temperature','K',.true.),   &
    output_variable('q2m','humidity_mixing_ratio',                                &
                    '2 meters vapour mixing ratio','kg kg-1',.true.),             &
    output_variable('smw','soil_moisture_content',                                &
                    'Moisture content','kg kg-1',.true.),                         &
    output_variable('tpr','precipitation_flux',                                   &
                    'Total precipitation','kg m-2 day-1',.true.),                 &
    output_variable('evp','water_evaporation_flux',                               &
                    'Total evapotranspiration','kg m-2 day-1',.true.),            &
    output_variable('runoff','surface_runoff_flux',                               &
                    'Surface runoff','kg m-2 day-1',.true.),                      &
    output_variable('scv','snowfall_flux',                                        &
                    'Snow precipitation','kg m-2 day-1',.true.),                  &
    output_variable('sena','surface_downward_sensible_heat_flux',                 &
                    'Sensible heat flux','W m-2',.true.),                         &
    output_variable('flw','net_upward_longwave_flux_in_air',                      &
                    'Net infrared energy flux','W m-2',.true.),                   &
    output_variable('fsw','net_downward_shortwave_flux_in_air',                   &
                    'Net solar absorbed energy flux','W m-2',.true.),             &
    output_variable('fld','surface_downwelling_longwave_flux_in_air',             &
                    'Downward LW flux','W m-2',.true.),                           &
    output_variable('sina','surface_downwelling_shortwave_flux_in_air',           &
                    'Incident solar energy flux','W m-2',.true.),                 &
    output_variable('prcv','convective_rainfall_flux',                            &
                    'Convective precipitation','kg m-2 day-1',.true.),            &
    output_variable('zpbl','atmosphere_boundary_layer_thickness',                 &
                    'PBL layer thickness','m',.true.),                            &
    output_variable('tgmax','surface_temperature',                                &
                    'Maximum surface temperature','K',.true.),                    &
    output_variable('tgmin','surface_temperature',                                &
                    'Minimum surface temperature','K',.true.),                    &
    output_variable('t2max','air_temperature',                                    &
                    'Maximum 2 meters temperature','K',.true.),                   &
    output_variable('t2min','air_temperature',                                    &
                    'Minimum 2 meters temperature','K',.true.),                   &
    output_variable('w10max','wind_speed',                                        &
                    'Maximum speed of 10m wind','m s-1',.true.),                  &
    output_variable('ps_min','air_pressure',                                      &
                    'Minimum of surface pressure','hPa',.true.),                  &
    output_variable('aldirs','surface_albedo_short_wave_direct',                  &
                    'Surface albedo to direct short wave radiation','1',.true.),  &
    output_variable('aldifs','surface_albedo_short_wave_diffuse',                 &
                    'Surface albedo to diffuse short wave radiation','1',.true.), &
    output_variable('seaice','seaice_binary_mask',                                &
                    'Sea ice mask','1',.false.) /

  data sub_variables / &
    output_variable('time','','','',.true.),                                     &
    output_variable('ps','','','',.true.),                                       &
    output_variable('u10m','eastward_wind',                                      &
                    '10 meters U component (westerly) of wind','m s-1',.true.),  &
    output_variable('v10m','northward_wind',                                     &
                    '10 meters V component (southerly) of wind','m s-1',.true.), &
    output_variable('uvdrag','surface_drag_coefficient_in_air',                  &
                    'Surface drag stress','1',.true.),                           &
    output_variable('tg','surface_temperature',                                  &
                    'Ground temperature','K',.true.),                            &
    output_variable('tlef','canopy_temperature',                                 &
                    'Foliage temperature','K',.true.),                           &
    output_variable('t2m','air_temperature','2 meters temperature','K',.true.),  &
    output_variable('q2m','humidity_mixing_ratio',                               &
                    '2 meters vapour mixing ratio','kg kg-1',.true.),            &
    output_variable('smw','soil_moisture_content',                               &
                    'Moisture content','kg kg-1',.true.),                        &
    output_variable('tpr','precipitation_flux',                                  &
                    'Total precipitation','kg m-2 day-1',.true.),                &
    output_variable('evp','water_evaporation_flux',                              &
                    'Total evapotranspiration','kg m-2 day-1',.true.),           &
    output_variable('runoff','surface_runoff_flux',                              &
                    'Surface runoff','kg m-2 day-1',.true.),                     &
    output_variable('scv','snowfall_flux',                                       &
                    'Snow precipitation','kg m-2 day-1',.true.),                 &
    output_variable('sena','surface_downward_sensible_heat_flux',                &
                    'Sensible heat flux','W m-2',.true.),                        &
    output_variable('prcv','convective_rainfall_flux',                           &
                    'Convective precipitation','kg m-2 day-1',.true.) /

  data rad_variables / &
    output_variable('time','','','',.true.),                                         &
    output_variable('ps','','','',.true.),                                           &
    output_variable('cld','cloud_area_fraction_in_atmosphere_layer',                 &
                    'Cloud fractional cover','1',.true.),                            &
    output_variable('clwp','atmosphere_cloud_liquid_water_content',                  &
                    'Cloud liquid water content','g m-2',.true.),                    &
    output_variable('qrs','tendency_of_air_temperature_due_to_shortwave_heating',    &
                    'Solar heating rate','K s-1',.true.),                            &
    output_variable('qrl','tendency_of_air_temperature_due_to_longwave_heating',     &
                    'Longwave cooling rate','K s-1',.true.),                         &
    output_variable('frsa','surface_downwelling_shortwave_flux_in_air',              &
                    'Surface absorbed solar flux','W m-2',.true.),                   &
    output_variable('frla','downwelling_longwave_flux_in_air',                       &
                    'Longwave cooling of surface flux','W m-2',.true.),              &
    output_variable('clrst','downwelling_shortwave_flux_in_air_assuming_clear_sky',  &
                    'clearsky total column absorbed solar flux','W m-2',.true.),     &
    output_variable('clrss','net_downward_shortwave_flux_in_air_assuming_clear_sky', &
                    'clearsky surface absorbed solar flux','W m-2',.true.),          &
    output_variable('clrlt','toa_net_upward_longwave_flux_assuming_clear_sky',       &
                    'clearsky net upward LW flux at TOA','W m-2',.true.),            &
    output_variable('clrls','net_upward_longwave_flux_in_air_assuming_clear_sky',    &
                    'clearsky LW cooling at surface','W m-2',.true.),                &
    output_variable('solin','toa_instantaneous_shortwave_forcing',                   &
                    'Instantaneous incident solar','W m-2',.true.),                  &
    output_variable('sabtp','atmosphere_net_rate_of_absorption_of_shortwave_energy', &
                    'Total column absorbed solar flux','W m-2',.true.),              &
    output_variable('firtp','atmosphere_net_rate_of_absorption_of_longwave_energy',  &
                    'net upward LW flux at TOA','W m-2',.true.) /

  data che_variables / &
    output_variable('time','','','',.true.),                                    &
    output_variable('ps','','','',.true.),                                      &
    output_variable('trac','atmosphere_mixing_ratio_of_tracer',                 &
                    'Tracers mixing ratios','kg kg-1',.true.),                  &
    output_variable('aext8','aerosol_optical_depth',                            &
                    'aer mix. aod.','1',.true.),                                &
    output_variable('assa8','aerosol_single_scattering_albedo',                 &
                    'aer mix. sin. scat. alb','1',.true.),                      &
    output_variable('agfu8','aerosol_asymmetry_parameter',                      &
                    'aer mix. sin. scat. asy','1',.true.),                      &
    output_variable('colb','instantaneous_column_burden',                       &
                    'columnburden inst','mg m-2',.true.),                       &
    output_variable('wdlsc',                                                    &
       'tendency_of_wet_deposition_of_tracer_due_to_large_scale_precipitation', &
                    'wet dep lgscale','mg m-2 day-1',.true.),                   &
    output_variable('wdcvc',                                                    &
       'tendency_of_wet_deposition_of_tracer_due_to_convective_precipitation',  &
                    'wet dep convect','mg m-2 day-1',.true.),                   &
    output_variable('sdrdp','tendency_of_dry_deposition_of_tracer',             &
                    'surf dry depos','mg m-2 day-1',.true.),                    &
    output_variable('xgasc','tendency_of_gas_conversion_of_tracer',             &
                    'chem gas conv','mg m-2 day-1',.true.),                     &
    output_variable('xaquc','tendency_of_aqueous_conversion_of_tracer',         &
                    'chem aqu conv','mg m-2 day-1',.true.),                     &
    output_variable('emiss','tendency_of_surface_emission_of_tracer',           &
                    'surf emission','mg m-2 day-1',.true.),                     &
    output_variable('acstoarf','toa_instantaneous_shortwave_radiative_forcing', &
                    'TOArad SW forcing av.','W m-2',.true.),                    &
    output_variable('acstsrrf','surface_shortwave_radiative_forcing',           &
                    'SRFrad SW forcing av.','W m-2',.true.),                    &
    output_variable('acstalrf','toa_longwave_radiative_forcing',                &
                    'TOArad LW forcing av.','W m-2',.true.),                    &
    output_variable('acssrlrf','surface_longwave_radiative_forcing',            &
                    'SRFrad LW forcing av.','W m-2',.true.) /

  data lak_variables / &
    output_variable('time','','','',.true.),                                      &
    output_variable('ps','','','',.true.),                                        &
    output_variable('tg','surface_temperature',                                   &
                    'Ground temperature','K',.true.),                             &
    output_variable('tpr','precipitation_flux',                                   &
                    'Total precipitation','kg m-2 day-1',.true.),                 &
    output_variable('scv','snowfall_flux',                                        &
                    'Snow precipitation','kg m-2 day-1',.true.),                  &
    output_variable('sena','surface_downward_sensible_heat_flux',                 &
                    'Sensible heat flux','W m-2',.true.),                         &
    output_variable('flw','net_upward_longwave_flux_in_air',                      &
                    'Net infrared energy flux','W m-2',.true.),                   &
    output_variable('fsw','net_downward_shortwave_flux_in_air',                   &
                    'Net solar absorbed energy flux','W m-2',.true.),             &
    output_variable('fld','surface_downwelling_longwave_flux_in_air',             &
                    'Downward LW flux','W m-2',.true.),                           &
    output_variable('sina','surface_downwelling_shortwave_flux_in_air',           &
                    'Incident solar energy flux','W m-2',.true.),                 &
    output_variable('aldirs','surface_albedo_short_wave_direct',                  &
                    'Surface albedo to direct short wave radiation','1',.true.),  &
    output_variable('aldifs','surface_albedo_short_wave_diffuse',                 &
                    'Surface albedo to diffuse short wave radiation','1',.true.), &
    output_variable('evl','water_evaporation_flux_where_sea_ice',                 &
                    'Water evaporation','mm sec-1',.true.),                       &
    output_variable('aveice','floating_ice_thickness',                            &
                    'Floating ice thickness','mm',.true.),                        &
    output_variable('hsnow','surface_snow_thickness_where_sea_ice',               &
                    'Floating snow thickness','mm',.true.),                       &
    output_variable('tlake','water_temperature',                                  &
                    'Lake water temperature','K',.true.) /

  contains

  subroutine allocate_mod_runparams
    implicit none
    call getmem1d(a,1,kz,'mod_runparams:a')
    call getmem1d(anudg,1,kz,'mod_runparams:anudg')
    call getmem1d(dsigma,1,kz,'mod_runparams:dsigma')
    call getmem1d(qcon,1,kz,'mod_runparams:qcon')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
    call getmem2d(twt,1,kz,1,2,'mod_runparams:twt')
    call getmem1d(wgtd,1,nspgd,'mod_runparams:wgtd')
    call getmem1d(wgtx,1,nspgx,'mod_runparams:wgtx')
    call getmem1d(dtau,1,nsplit,'mod_runparams:nsplit')
  end subroutine allocate_mod_runparams
!
  logical function iswater(a)
    real(8) , intent(in) :: a
    iswater = .false.
    if (a > 13.5D0 .and. a < 15.5D0) iswater = .true.
  end function
!
  logical function isocean(a)
    real(8) , intent(in) :: a
    isocean = .false.
    if (a > 14.5D0 .and. a < 15.5D0) isocean = .true.
  end function
!
  logical function islake(a)
    real(8) , intent(in) :: a
    islake = .false.
    if (a > 13.5D0 .and. a < 14.5D0) islake = .true.
  end function
!
end module mod_runparams
