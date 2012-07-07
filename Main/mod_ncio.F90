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
module mod_ncio
!
  use netcdf
  use mod_runparams
  use mod_mpmessage
  use mod_memutil
  use mod_nchelper
  use mod_domain
  use mod_runparams , only : iqc , iqv
!
  integer , parameter :: n_atmvar = 14
  integer , parameter :: n_srfvar = 26
  integer , parameter :: n_subvar = 16
  integer , parameter :: n_radvar = 18
  integer , parameter :: n_lakvar = 16
  integer , parameter :: n_stsvar = 13

  type output_variable
    character(len=8) :: vname
    character(len=128) :: vstd_name
    character(len=128) :: vdesc
    character(len=16) :: vunit
    character(len=16) :: time_meth
    logical :: enabled
  end type

  type(output_variable) , dimension(n_atmvar) :: atm_variables
  type(output_variable) , dimension(n_srfvar) :: srf_variables
  type(output_variable) , dimension(n_subvar) :: sub_variables
  type(output_variable) , dimension(n_radvar) :: rad_variables
  type(output_variable) , dimension(n_lakvar) :: lak_variables
  type(output_variable) , dimension(n_stsvar) :: sts_variables

  data atm_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('u','eastward_wind','U component (westerly) of wind',     &
            'm s-1','point',.true.),                                          &
    output_variable('v','northward_wind','V component (southerly) of wind',   &
            'm s-1','point',.true.),                                          &
    output_variable('omega','lagrangian_tendency_of_air_pressure',            &
            'Pressure velocity','hPa s-1','point',.true.),                    &
    output_variable('t','air_temperature',                                    &
            'Temperature','K','point',.true.),                                &
    output_variable('qv','humidity_mixing_ratio',                             &
            'Water vapor mixing ratio','kg kg-1','point',.true.),             &
    output_variable('qc','cloud_liquid_water_mixing_ratio',                   &
            'Cloud water mixing ratio','kg kg-1','point',.true.),             &
    output_variable('tke','turbulent_kinetic_energy',                         &
            'Turbulent Kinetic Energy','m2 s2','point',.true.),               &
    output_variable('kth','vertical_momentum_diffusivity',                    &
            'Vertical Turbulent Viscosity','m2 s-1','point',.true.),          &
    output_variable('kzm','vertical_scalar_diffusivity',                      &
            'Vertical Turbulent Diffusivity','m2 s-1','point',.true.),        &
    output_variable('tpr','precipitation_flux',                               &
            'Total daily precipitation rate','kg m-2 day-1','point',.true.),  &
    output_variable('tgb','soil_temperature',                                 &
            'Lower groud temperature','K','point',.true.),                    &
    output_variable('swt','moisture_content_of_soil_layer',                   &
            'Total soil water','kg m-2','point',.true.) /

  data srf_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('tbnds','','','','',.true.),                              &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('u10m','eastward_wind',                                   &
          '10 meters U component (westerly) of wind','m s-1','point',.true.), &
    output_variable('v10m','northward_wind',                                  &
          '10 meters V component (southerly) of wind','m s-1','point',.true.),&
    output_variable('uvdrag','surface_drag_coefficient_in_air',               &
          'Surface drag stress','1','point',.true.),                          &
    output_variable('tg','surface_temperature',                               &
          'Ground temperature','K','point',.true.),                           &
    output_variable('tlef','canopy_temperature',                              &
          'Foliage temperature','K','point',.true.),                          &
    output_variable('t2m','air_temperature',                                  &
          '2 meters temperature','K','point',.true.),                         &
    output_variable('q2m','specific_humidity',                                &
          '2 meters Specific humidity','1','point',.true.),                   &
    output_variable('smw','soil_moisture_content',                            &
          'Moisture content','kg m-2','point',.true.),                        &
    output_variable('tpr','precipitation_flux',                               &
          'Total precipitation','kg m-2 day-1','mean',.true.),                &
    output_variable('evp','water_evaporation_flux',                           &
          'Total evapotranspiration','kg m-2 day-1','mean',.true.),           &
    output_variable('runoff','surface_runoff_flux',                           &
          'Surface runoff','kg m-2 day-1','mean',.true.),                     &
    output_variable('scv','lwe_thickness_of_surface_snow_amount',             &
          'Snow thickness','kg m-2','mean',.true.),                           &
    output_variable('sena','surface_downward_sensible_heat_flux',             &
          'Sensible heat flux','W m-2','mean',.true.),                        &
    output_variable('flw','net_upward_longwave_flux_in_air',                  &
          'Net infrared energy flux','W m-2','mean',.true.),                  &
    output_variable('fsw','net_downward_shortwave_flux_in_air',               &
          'Net solar absorbed energy flux','W m-2','mean',.true.),            &
    output_variable('fld','surface_downwelling_longwave_flux_in_air',         &
          'Downward LW flux','W m-2','mean',.true.),                          &
    output_variable('sina','surface_downwelling_shortwave_flux_in_air',       &
          'Incident visible solar energy flux','W m-2','mean',.true.),        &
    output_variable('prcv','convective_rainfall_flux',                        &
          'Convective precipitation','kg m-2 day-1','mean',.true.),           &
    output_variable('zpbl','atmosphere_boundary_layer_thickness',             &
          'PBL layer thickness','m','point',.true.),                          &
    output_variable('aldirs','surface_albedo_short_wave_direct',              &
          'Surface albedo to direct short wave radiation','1','point',.true.),&
    output_variable('aldifs','surface_albedo_short_wave_diffuse',             &
         'Surface albedo to diffuse short wave radiation','1','point',.true.),&
    output_variable('sund','duration_of_sunshine',                            &
         'Duration of sunshine','s','sum',.true.),                            &
    output_variable('seaice','seaice_binary_mask',                            &
          'Sea ice mask','1','point',.false.) /

  data sts_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('tbnds','','','','',.true.),                              &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('tgmax','surface_temperature',                            &
          'Maximum surface temperature','K','maximum',.true.),                &
    output_variable('tgmin','surface_temperature',                            &
          'Minimum surface temperature','K','minimum',.true.),                &
    output_variable('t2max','air_temperature',                                &
          'Maximum 2 meters temperature','K','maximum',.true.),               &
    output_variable('t2min','air_temperature',                                &
          'Minimum 2 meters temperature','K','minimum',.true.),               &
    output_variable('t2avg','air_temperature',                                &
          'Average 2 meters temperature','K','mean',.true.),                  &
    output_variable('w10max','wind_speed',                                    &
          'Maximum speed of 10m wind','m s-1','maximum',.true.),              &
    output_variable('pcpmax','precipitation_flux',                            &
          'Maximum precipitation flux','kg m-2 s-1','maximum',.true.),        &
    output_variable('pcpavg','precipitation_flux',                            &
          'Average precipitation flux','kg m-2 s-1','mean',.true.),           &
    output_variable('sund','duration_of_sunshine',                            &
          'Duration of sunshine','s','sum',.true.),                           &
    output_variable('ps_min','air_pressure',                                  &
          'Minimum of surface pressure','hPa','minimum',.true.) /

  data sub_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('u10m','eastward_wind',                                   &
          '10 meters U component (westerly) of wind','m s-1','point',.true.), &
    output_variable('v10m','northward_wind',                                  &
          '10 meters V component (southerly) of wind','m s-1','point',.true.),&
    output_variable('uvdrag','surface_drag_coefficient_in_air',               &
          'Surface drag stress','1','point',.true.),                          &
    output_variable('tg','surface_temperature',                               &
          'Ground temperature','K','point',.true.),                           &
    output_variable('tlef','canopy_temperature',                              &
          'Foliage temperature','K','point',.true.),                          &
    output_variable('t2m','air_temperature',                                  &
          '2 meters temperature','K','point',.true.),                         &
    output_variable('q2m','specific_humidity',                                &
          '2 meters Specific humidity','1','point',.true.),                   &
    output_variable('smw','soil_moisture_content',                            &
          'Soil moisture content','kg kg-1','point',.true.),                  &
    output_variable('tpr','precipitation_flux',                               &
          'Total precipitation','kg m-2 day-1','mean',.true.),                &
    output_variable('evp','water_evaporation_flux',                           &
          'Total evapotranspiration','kg m-2 day-1','point',.true.),          &
    output_variable('runoff','surface_runoff_flux',                           &
          'Surface runoff','kg m-2 day-1','mean',.true.),                     &
    output_variable('scv','lwe_thickness_of_surface_snow_amount',             &
          'Snow thickness','kg m-2','mean',.true.),                           &
    output_variable('sena','surface_downward_sensible_heat_flux',             &
          'Sensible heat flux','W m-2','mean',.true.),                        &
    output_variable('prcv','convective_rainfall_flux',                        &
          'Convective precipitation','kg m-2 day-1','mean',.true.) /

  data rad_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('cld','cloud_area_fraction_in_atmosphere_layer',          &
       'Cloud fractional cover','1','point',.true.),                          &
    output_variable('clwp','thickness_of_liquid_water_cloud',                 &
       'Cloud liquid water path','g m-2','point',.true.),                     &
    output_variable('qrs',                                                    &
       'tendency_of_air_temperature_due_to_shortwave_heating',                &
       'Solar heating rate','K s-1','point',.true.),                          &
    output_variable('qrl',                                                    &
       'tendency_of_air_temperature_due_to_longwave_heating',                 &
       'Longwave cooling rate','K s-1','point',.true.),                       &
    output_variable('frsa','surface_net_downward_shortwave_flux',             &
       'Surface net downward shortwave flux','W m-2','point',.true.),         &
    output_variable('frla','surface_net_upward_longwave_flux',                &
       'Surface net upward longwave flux','W m-2','point',.true.),            &
    output_variable('clrst',                                                  &
       'toa_net_downward_shortwave_flux_assuming_clear_sky',                  &
       'Clearsky TOA net downward shortwave flux','W m-2','point',.true.),    &
    output_variable('clrss',                                                  &
       'surface_net_downward_shortwave_flux_assuming_clear_sky',              &
       'Surface net downward shortwave flux','W m-2','point',.true.),         &
    output_variable('clrlt','toa_net_upward_longwave_flux_assuming_clear_sky',&
       'Clearsky TOA net upward longwave flux','W m-2','point',.true.),       &
    output_variable('clrls',                                                  &
       'surface_net_upward_longwave_flux_assuming_clear_sky',                 &
       'Clearsky net upward longwave flux','W m-2','point',.true.),           &
    output_variable('solin','toa_incoming_shortwave_flux',                    &
       'Incoming solar flux','W m-2','point',.true.),                         &
    output_variable('sabtp','toa_net_upward_shortwave_flux',                  &
       'Net TOA upward shortwave flux','W m-2','point',.true.),               &
    output_variable('totcf','cloud_area_fraction',                            &
       'Total cloud fraction','1','point',.true.),                            &
    output_variable('totcl','atmosphere_cloud_condensed_water_content',       &
       'Total columnar water content','kg m-2','point',.true.),               &
    output_variable('totci','atmosphere_ice_condensed_water_content',         &
       'Total columnar ice content','kg m-2','point',.true.),                 &
    output_variable('firtp','toa_net_upward_longwave_flux',                   &
       'net upward LW flux at TOA','W m-2','point',.true.) /

  data lak_variables / &
    output_variable('time','','','','',.true.),                               &
    output_variable('ps','','','','point',.true.),                            &
    output_variable('tg','surface_temperature',                               &
         'Ground temperature','K','point',.true.),                            &
    output_variable('tpr','precipitation_flux',                               &
         'Total precipitation','kg m-2 day-1','point',.true.),                &
    output_variable('scv','lwe_thickness_of_surface_snow_amount',             &
          'Snow thickness','kg m-2','mean',.true.),                           &
    output_variable('sena','surface_downward_sensible_heat_flux',             &
         'Sensible heat flux','W m-2','point',.true.),                        &
    output_variable('flw','net_upward_longwave_flux_in_air',                  &
         'Net infrared energy flux','W m-2','point',.true.),                  &
    output_variable('fsw','net_downward_shortwave_flux_in_air',               &
         'Net solar absorbed energy flux','W m-2','point',.true.),            &
    output_variable('fld','surface_downwelling_longwave_flux_in_air',         &
         'Downward LW flux','W m-2','point',.true.),                          &
    output_variable('sina','surface_downwelling_shortwave_flux_in_air',       &
         'Incident solar energy flux','W m-2','point',.true.),                &
    output_variable('aldirs','surface_albedo_short_wave_direct',              &
         'Surface albedo to direct short wave radiation','1','point',.true.), &
    output_variable('aldifs','surface_albedo_short_wave_diffuse',             &
         'Surface albedo to diffuse short wave radiation','1','point',.true.),&
    output_variable('evl','water_evaporation_flux_where_sea_ice',             &
         'Water evaporation','mm sec-1','point',.true.),                      &
    output_variable('aveice','floating_ice_thickness',                        &
         'Floating ice thickness','mm','point',.true.),                       &
    output_variable('hsnow','surface_snow_thickness_where_sea_ice',           &
         'Floating snow thickness','mm','point',.true.),                      &
    output_variable('tlake','water_temperature',                              &
         'Lake water temperature','K','point',.true.) /
!
  public :: ivarname_lookup
  public :: init_mod_ncio , release_mod_ncio
  public :: read_domain_info , read_domain_lake,   &
            read_subdomain , read_subdomain_lake,  &
            close_domain
  public :: open_icbc , read_icbc , icbc_search
  public :: prepare_common_out
  public :: writerec_atm , writerec_srf , writerec_sub , &
            writerec_rad , writerec_lak , writerec_sts
!
  integer :: idmin , isdmin , ibcin , ncatm , ncsrf , ncsts , &
             ncsub , ncrad , nclak
  integer :: istatus
  integer :: ibcrec , ibcnrec
  integer :: iatmrec , isrfrec , istsrec , isubrec , iradrec , ilakrec
  integer , dimension(n_atmvar) :: iatmvar
  integer , dimension(n_srfvar) :: isrfvar
  integer , dimension(n_stsvar) :: istsvar
  integer , dimension(n_subvar) :: isubvar
  integer , dimension(n_radvar) :: iradvar
  integer , dimension(n_lakvar) :: ilakvar
  character(256) :: dname , sdname , icbcname
  type(rcm_time_and_date) , dimension(:) , allocatable :: icbc_idate
  integer , dimension(7) :: icbc_ivar
  real(dp) :: tpd
  real(dp) :: xns2d
  real(sp) :: xns2r
  type(rcm_time_and_date) , save :: cordex_refdate

  ! DIM1 is iy ,   DIM2 is jx , DIM3 is time ,       DIM4 is kz
  ! DIM5 is m10 ,  DIM6 is m2 , DIM7 is soil_layer , DIM8 is ntimes
  ! DIM9 is ntr ,  DIM10 is depth for lake
  integer , dimension(10) :: idims

  integer :: o_is
  integer :: o_ie
  integer :: o_js
  integer :: o_je
  integer :: o_ni
  integer :: o_nj
  integer :: o_isg
  integer :: o_ieg
  integer :: o_jsg
  integer :: o_jeg
  integer :: o_nig
  integer :: o_njg
  integer :: o_nz
  logical :: lwrap , lmaskfill

  real(sp) , dimension(:,:) , pointer :: ioxlat
  real(sp) , dimension(:,:) , pointer :: ioxlon
  real(sp) , dimension(:,:) , pointer :: iotopo
  real(sp) , dimension(:,:) , pointer :: iomask
  real(sp) , dimension(:,:) , pointer :: iolnds
  real(sp) , dimension(:,:) , pointer :: ioxlat_s
  real(sp) , dimension(:,:) , pointer :: ioxlon_s
  real(sp) , dimension(:,:) , pointer :: iotopo_s
  real(sp) , dimension(:,:) , pointer :: iomask_s
  real(sp) , dimension(:,:) , pointer :: subio
  real(sp) , dimension(:,:,:) , pointer :: dumio
  real(sp) , dimension(:,:) , pointer :: sp2d
  real(sp) , dimension(:,:) , pointer :: sp2d1
  real(sp) , dimension(:,:,:) , pointer :: atmsrfmask
  real(sp) , dimension(:,:) , pointer :: atmsrfsum

  integer , dimension(numbat) :: lak_fbats

  data lmaskfill /.false./
  data idmin   /-1/
  data isdmin  /-1/
  data ibcin   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/
  data ncatm   /-1/
  data iatmrec / 1/
  data ncsrf   /-1/
  data ncsts   /-1/
  data isrfrec / 1/
  data istsrec / 1/
  data ncsub   /-1/
  data isubrec / 1/
  data ncrad   /-1/
  data iradrec / 1/
  data nclak   /-1/
  data ilakrec / 1/

  data lak_fbats / 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, &
                   1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                   0, 0/

contains

  function ivarname_lookup(ctype,sname)
    implicit none
    integer :: ivarname_lookup
    character(3) , intent(in) :: ctype
    character(len=*) , intent(in) :: sname
    integer :: i

    ivarname_lookup = -1

    if (ctype == 'ATM') then
      do i = 1 , n_atmvar
        if (sname == atm_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'SRF') then
      do i = 1 , n_srfvar
        if (sname == srf_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'SUB') then
      do i = 1 , n_subvar
        if (sname == sub_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'STS') then
      do i = 1 , n_stsvar
        if (sname == sts_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'RAD') then
      do i = 1 , n_radvar
        if (sname == rad_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'LAK') then
      do i = 1 , n_lakvar
        if (sname == lak_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    endif
  end function ivarname_lookup

  subroutine init_mod_ncio(lband)
    implicit none
    logical , intent(in) :: lband
    character(3) :: sbstring
    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    write (sbstring,'(i0.3)') nsg
    sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//sbstring//'.nc'
    icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'//'YYYYMMDDHH.nc'

    xns2r = 1.0/real(nnsg)
    xns2d = 1.0D0/dble(nnsg)
    cordex_refdate = 1949120100
    call setcal(cordex_refdate,ical)

    if (lband) then
      o_is = 2
      o_ie = iym2
      o_js = 1
      o_je = jx
      o_ni = iym3
      o_nj = jx
      o_isg = nsg+1
      o_ieg = iysg-2*nsg
      o_jsg = 1
      o_jeg = jxsg
      o_nig = iym3sg
      o_njg = jxsg
      o_nz = kz
      lwrap = .true.
    else
      o_is = 2
      o_ie = iym2
      o_js = 2
      o_je = jxm2
      o_ni = iym3
      o_nj = jxm3
      o_isg = nsg+1
      o_ieg = iym2sg
      o_jsg = nsg+1
      o_jeg = jxm2sg
      o_nig = iym3sg
      o_njg = jxm3sg
      o_nz = kz
      lwrap = .false.
    end if
    call getmem2d(ioxlat,1,o_nj,1,o_ni,'ncio:ioxlat')
    call getmem2d(ioxlon,1,o_nj,1,o_ni,'ncio:ioxlon')
    call getmem2d(iotopo,1,o_nj,1,o_ni,'ncio:iotopo')
    call getmem2d(iomask,1,o_nj,1,o_ni,'ncio:iomask')
    call getmem2d(iolnds,1,o_nj,1,o_ni,'ncio:iolnds')
    call getmem3d(dumio,1,o_nj,1,o_ni,1,o_nz,'ncio:dumio')
    call getmem2d(sp2d,1,jx,1,iy,'ncio:sp2d')
    call getmem3d(atmsrfmask,1,nnsg,1,o_nj,1,o_ni,'ncio:atmsrfmask')
    call getmem2d(atmsrfsum,1,o_nj,1,o_ni,'ncio:atmsrfsum')
    if (nsg > 1) then
      call getmem2d(ioxlat_s,1,o_njg,1,o_nig,'ncio:ioxlat_s')
      call getmem2d(ioxlon_s,1,o_njg,1,o_nig,'ncio:ioxlon_s')
      call getmem2d(iotopo_s,1,o_njg,1,o_nig,'ncio:iotopo_s')
      call getmem2d(iomask_s,1,o_njg,1,o_nig,'ncio:iomask_s')
      call getmem2d(subio,1,o_njg,1,o_nig,'ncio:subio')
      call getmem2d(sp2d1,1,jxsg,1,iysg,'ncio:sp2d1')
    end if
  end subroutine init_mod_ncio

  subroutine read_domain_info
    implicit none
    write (aline,*) 'open_domain: READING HEADER FILE:', dname
    call say
    call openfile_withname(dname,idmin)
    call read_domain(idmin)
    tpd = houpd/atmfrq
    ioxlat(:,:) = real(mddom_io%xlat(o_js:o_je,o_is:o_ie))
    ioxlon(:,:) = real(mddom_io%xlon(o_js:o_je,o_is:o_ie))
    iotopo(:,:) = real(mddom_io%ht(o_js:o_je,o_is:o_ie))
    iomask(:,:) = real(mddom_io%mask(o_js:o_je,o_is:o_ie))
    iolnds(:,:) = real(mddom_io%lndcat(o_js:o_je,o_is:o_ie))
  end subroutine read_domain_info

  subroutine read_domain_lake(hlake)
    implicit none

    real(dp) , pointer , dimension(:,:,:) , intent(out) :: hlake
    integer :: n , ivarid

    if (idmin < 0) then
      write (6,*) 'Error : Domain file not in open state'
      call fatal(__FILE__,__LINE__,'DOMAIN FILE')
    end if

    istatus = nf90_inq_varid(idmin, 'dhlake', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dhlake miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable dhlake read error','DOMAIN FILE')
    do n = 1 , nnsg
      hlake(n,:,:) = dble(sp2d)
    end do
  end subroutine read_domain_lake

  subroutine read_subdomain(ht1,lnd1,xlat1,xlon1)
    implicit none

    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ht1
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: lnd1
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: xlat1
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: xlon1

    integer :: ivarid
    
    if ( nsg > 1 ) then
      write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
      call say
      call openfile_withname(sdname,isdmin)
    end if
    istatus = nf90_inq_varid(isdmin, 'topo', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable topo miss', 'SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable topo read error','SUBDOMAIN FILE')
    call reorder_2_3(sp2d1,ht1)
    iotopo_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'landuse', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable landuse miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable landuse read error', &
                  'SUBDOMAIN FILE')
    call reorder_2_3(sp2d1,lnd1)
    istatus = nf90_inq_varid(isdmin, 'xlat', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlat miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlat read error','SUBDOMAIN FILE')
    call reorder_2_3(sp2d1,xlat1)
    ioxlat_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'xlon', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlon miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlon read error','SUBDOMAIN FILE')
    call reorder_2_3(sp2d1,xlon1)
    ioxlon_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable mask read error','SUBDOMAIN FILE')
    iomask_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
  end subroutine read_subdomain

  subroutine read_subdomain_lake(hlake1)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: hlake1
    integer :: ivarid
    if (isdmin < 0) then
      write (6,*) 'Error : Subdom file not in open state'
      call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE')
    end if
    istatus = nf90_inq_varid(isdmin, 'dhlake', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dhlake miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable dhlake read error', &
                  'SUBDOMAIN FILE')
    call reorder_2_3(sp2d1,hlake1)
  end subroutine read_subdomain_lake

  subroutine close_domain
    implicit none
    if (idmin >= 0) then
      istatus = nf90_close(idmin)
      call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
      idmin = -1
    end if
    if ( nsg>1 .and. isdmin >=0 ) then
      istatus = nf90_close(isdmin)
      call check_ok(__FILE__,__LINE__,'SubDomain file close error', &
                   'SUBDOMAIN FILE')
      isdmin = -1
    end if
  end subroutine close_domain

  integer function icbc_search(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1, appdat2
    if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
      icbc_search = -1
    else
      tdif = idate-icbc_idate(1)
      ibcrec = (idnint(tohours(tdif))/ibdyfrq)+1
      if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
        appdat1 = tochar(idate)
        write (6,*) 'Record is not found in ICBC file for ',appdat1
        appdat1 = tochar(icbc_idate(1))
        appdat2 = tochar(icbc_idate(ibcnrec))
        write (6,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      icbc_search = ibcrec
    end if 
  end function icbc_search

  subroutine open_icbc(idate)
    type(rcm_time_and_date) , intent(in) :: idate
    character(10) :: ctime
    integer :: idimid , itvar , i , chkdiff
    real(dp) , dimension(:) , allocatable :: icbc_nctime
    character(64) :: icbc_timeunits , icbc_timecal

    call close_icbc
    write (ctime, '(i10)') toint10(idate)
    icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'//ctime//'.nc'
    call openfile_withname(icbcname,ibcin)
    ibcrec = 1
    ibcnrec = 0
    call check_domain(ibcin,.true.)
    istatus = nf90_inq_dimid(ibcin, 'time', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
    call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
    if ( ibcnrec < 1 ) then
      write (6,*) 'Time var in ICBC has zero dim.'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_inq_varid(ibcin, 'time', itvar)
    call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
    call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'calendar', icbc_timecal)
    call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
    allocate(icbc_nctime(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(6,*) 'Memory allocation error in ICBC for time real values'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    allocate(icbc_idate(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(6,*) 'Memory allocation error in ICBC for time array'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_get_var(ibcin, itvar, icbc_nctime)
    call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
    do i = 1 , ibcnrec
      icbc_idate(i) = timeval2date(icbc_nctime(i), icbc_timeunits, icbc_timecal)
    end do
    if ( ibcnrec > 1 ) then
      chkdiff = idnint(icbc_nctime(2) - icbc_nctime(1))
      if (chkdiff /= ibdyfrq) then
        write (6,*) 'Time var in ICBC inconsistency.'
        write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
        write (6,*) 'Found     ibdyfrq = ', chkdiff
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
    end if
    deallocate(icbc_nctime)
    istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable ps miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
    call check_ok(__FILE__,__LINE__,'variable ts miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
    call check_ok(__FILE__,__LINE__,'variable u miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
    call check_ok(__FILE__,__LINE__,'variable v miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
    call check_ok(__FILE__,__LINE__,'variable t miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
    call check_ok(__FILE__,__LINE__,'variable qv miss', 'ICBC FILE')
  end subroutine open_icbc

  subroutine read_icbc(ps,ts,u,v,t,qv)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: u
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: v
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: t
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: qv
    real(dp) , pointer , dimension(:,:) , intent(out) :: ps
    real(dp) , pointer , dimension(:,:) , intent(out) :: ts

    integer , dimension(4) :: istart , icount

    istart(3) = ibcrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin,icbc_ivar(1),ps,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(2),ts,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    istart(4) = ibcrec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin,icbc_ivar(3),u,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(4),v,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(5),t,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(6),qv,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
  end subroutine read_icbc

  subroutine close_icbc
    implicit none
    if (ibcin >= 0) then
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__, &
           'Error Close ICBC file '//trim(icbcname),'ICBC FILE')
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
  end subroutine close_icbc

  subroutine close_common(ncid, ctype)
    implicit none
    integer , intent(inout) :: ncid
    character(3) , intent(in) :: ctype
    if (ncid >= 0) then
      istatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Close '//ctype//' file', ctype//' FILE')
      ncid = -1
    end if
  end subroutine close_common

  subroutine prepare_common_out(idate,ctype)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(3) , intent(in) :: ctype
    character(32) :: fbname , ctime
    character(16) :: fterr
    character(256) :: ofname
    real(dp) :: hptop
    real(dp) :: rdum1
    real(dp) , dimension(2) :: rdum2
    real(dp) , dimension(iysg) :: yiy
    real(dp) , dimension(jxsg) :: xjx
    real(sp) , dimension(ndpmax) :: depth
    integer :: ncid
    integer , dimension(3) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(4) :: isrvvar
    integer , dimension(5) :: illtpvar
    integer :: itvar , imapvar , i , j

    integer , dimension(5) :: tyx
    integer , dimension(5) :: tzyx
    integer , dimension(5) :: t10yx
    integer , dimension(5) :: t2yx
    integer , dimension(5) :: tlyx
    integer , dimension(5) :: tcyx
    integer , dimension(5) :: tczyx
    integer , dimension(5) :: tdyx
    character(len=128) :: cdum

    if (ctype == 'ATM') then
      ncid = ncatm
      iatmrec = 1
    else if (ctype == 'SRF') then
      ncid = ncsrf
      isrfrec = 1
    else if (ctype == 'SUB') then
      ncid = ncsub
      isubrec = 1
    else if (ctype == 'STS') then
      ncid = ncsts
      istsrec = 1
    else if (ctype == 'RAD') then
      ncid = ncrad
      iradrec = 1
    else if (ctype == 'LAK') then
      ncid = nclak
      ilakrec = 1
    else
      write (aline,*) 'UNKNOWN IO TYPE : ', ctype
      call say
      write (aline,*) 'NOTHING TO DO'
      call say
      return
    end if

    call close_common(ncid, ctype)

    write (fterr, '(a3,a)') ctype, ' FILE'
    write (fbname,'(a,a,i10)') trim(ctype), '.', toint10(idate)
    ofname = trim(dirout)//pthsep//trim(domname)// &
             '_'//trim(fbname)//'.nc'
    ctime = tochar(cordex_refdate)
    write (aline, *) 'Opening new output file ', trim(ofname)
    call say

    call createfile_withname(ofname,ncid)
    call add_common_global_params(ncid,'Model ('//ctype//')')
!
!         ADD RUN PARAMETERS
!
    istatus = nf90_put_att(ncid, nf90_global, 'model_IPCC_scenario', scenario)
    call check_ok(__FILE__,__LINE__,'Error add scenario', fterr)
    call cdumlbcs(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_boundary_conditions' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add lbcs', fterr)
    call cdumcums(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_cumulous_convection_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add icup', fterr)
    if (icup == 2 .or. icup == 99 .or. icup == 98) then
      call cdumcumcl(cdum)
      istatus = nf90_put_att(ncid, nf90_global,  &
            'model_convective_closure_assumption' , trim(cdum))
      call check_ok(__FILE__,__LINE__,'Error add igcc', fterr)
    end if
    call cdumpbl(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_boundary_layer_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ibltyp', fterr)
    call cdummoist(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_moist_physics_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ipptls', fterr)
    call cdumocnflx(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_ocean_flux_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add iocnflx', fterr)
    call cdumpgfs(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_pressure_gradient_force_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ipgf', fterr)
    call cdumemiss(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_use_emission_factor' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add iemiss', fterr)
    call cdumlakes(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_use_lake_model' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add lakemod', fterr)
    call cdumchems(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_chemistry' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ichem', fterr)
    call cdumdcsst(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_diurnal_cycle_sst' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add dcsst', fterr)
    call cdumseaice(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_seaice_effect' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add seaice', fterr)
    call cdumdesseas(cdum)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_seasonal_desert_albedo_effect' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add desseas', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_initial_start' , tochar(globidate1))
    call check_ok(__FILE__,__LINE__,'Error add globidate1', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_start' , tochar(idate1))
    call check_ok(__FILE__,__LINE__,'Error add idate1', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_expected_end' , tochar(idate2))
    call check_ok(__FILE__,__LINE__,'Error add idate2', fterr)
    if (ifrest) then
      cdum = 'Yes'
    else
      cdum = 'No'
    end if
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_is_a_restart' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ifrest', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_seconds' , dt)
    call check_ok(__FILE__,__LINE__,'Error add dt', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_minutes_solar_rad_calc' , dtrad)
    call check_ok(__FILE__,__LINE__,'Error add dtrad', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_seconds_bats_calc' , dtsrf)
    call check_ok(__FILE__,__LINE__,'Error add dtsrf', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_hours_radiation_calc' , dtabem)
    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_hours_boundary_input' , ibdyfrq)
    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
!
!         ADD DIMENSIONS
!
    if (ctype == 'SUB') then
      istatus = nf90_def_dim(ncid, 'iy', o_nig, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_njg, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
    else
      istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
    end if
    istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
    call check_ok(__FILE__,__LINE__,'Error create dim time', fterr)
    istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
    call check_ok(__FILE__,__LINE__,'Error create dim kz', fterr)
!
!         OUT TYPE DEPENDENT DIMENSIONS
!
    if (ctype == 'SRF' .or. ctype == 'SUB' .or. ctype == 'STS') then
      istatus = nf90_def_dim(ncid, 'm10', 1, idims(5))
      call check_ok(__FILE__,__LINE__,'Error create dim m10', fterr)
      istatus = nf90_def_dim(ncid, 'm2', 1, idims(6))
      call check_ok(__FILE__,__LINE__,'Error create dim m2', fterr)
      istatus = nf90_def_dim(ncid, 'soil_layer', 2, idims(7))
      call check_ok(__FILE__,__LINE__,'Error create dim soil_layer', fterr)
      istatus = nf90_def_var(ncid, 'm10', nf90_double, idims(5), isrvvar(1))
      call check_ok(__FILE__,__LINE__,'Error add var m10', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'standard_name', 'height')
      call check_ok(__FILE__,__LINE__,'Error add m10 standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'long_name','10 m height level')
      call check_ok(__FILE__,__LINE__,'Error add m10 long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'positive', 'up')
      call check_ok(__FILE__,__LINE__,'Error add m10 positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'units', 'm')
      call check_ok(__FILE__,__LINE__,'Error add m10 units', fterr)
      istatus = nf90_def_var(ncid, 'm2', nf90_double, idims(6), isrvvar(2))
      call check_ok(__FILE__,__LINE__,'Error add var m2', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'standard_name', 'height')
      call check_ok(__FILE__,__LINE__,'Error add m2 standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'long_name','2 m height level')
      call check_ok(__FILE__,__LINE__,'Error add m2 long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'positive', 'up')
      call check_ok(__FILE__,__LINE__, 'Error add m2 positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'units', 'm')
      call check_ok(__FILE__,__LINE__, 'Error add m2 units', fterr)
      istatus = nf90_def_var(ncid, 'layer', nf90_double, idims(7), isrvvar(3))
      call check_ok(__FILE__,__LINE__,'Error add var layer', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'standard_name', &
                         'model_level_number')
      call check_ok(__FILE__,__LINE__,'Error add layer standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'long_name', &
                         'Surface and root zone')
      call check_ok(__FILE__,__LINE__,'Error add layer long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'positive', 'down')
      call check_ok(__FILE__,__LINE__,'Error add layer positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'units', '1')
      call check_ok(__FILE__,__LINE__,'Error add layer units', fterr)
    end if
    if (ctype == 'STS' .or. ctype == 'SRF') then
      istatus = nf90_def_dim(ncid, 'ntimes', 2, idims(8))
      call check_ok(__FILE__,__LINE__,'Error create dim ntimes', fterr)
    end  if
    if (ctype == 'LAK') then
      istatus = nf90_def_dim(ncid, 'depth', ndpmax, idims(10))
      call check_ok(__FILE__,__LINE__,'Error create dim depth', fterr)
    end if
    istatus = nf90_def_var(ncid, 'rcm_map', nf90_int, varid=imapvar)
    call check_ok(__FILE__,__LINE__,'Error add var rcm_map', fterr)
    if (iproj == 'LAMCON') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'lambert_conformal_conic')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'POLSTR') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'stereographic')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'NORMER') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'mercator')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, imapvar, &
            'grid_mapping_name', 'rotated_latitude_longitude')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    end if

    istatus = nf90_def_var(ncid, 'sigma', nf90_double, idims(4), izvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var sigma', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                         'atmosphere_sigma_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add sigma standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                         'Sigma at model layers')
    call check_ok(__FILE__,__LINE__,'Error add sigma long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add sigma units', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
    call check_ok(__FILE__,__LINE__,'Error add sigma axis', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
    call check_ok(__FILE__,__LINE__,'Error add sigma positive', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                         'sigma: sigma ps: ps ptop: ptop')
    call check_ok(__FILE__,__LINE__,'Error add sigma formula_terms', fterr)
    if (ctype == 'LAK') then
      istatus = nf90_def_var(ncid, 'depth', nf90_float, idims(10), izvar(3))
      call check_ok(__FILE__,__LINE__,'Error add var depth', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'standard_name', 'depth')
      call check_ok(__FILE__,__LINE__,'Error add depth standard_name', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'long_name', &
                           'Depth below surface')
      call check_ok(__FILE__,__LINE__,'Error add depth long_name', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'units', 'm')
      call check_ok(__FILE__,__LINE__,'Error add depth units', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'axis', 'Z')
      call check_ok(__FILE__,__LINE__,'Error add depth axis', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'positive', 'down')
      call check_ok(__FILE__,__LINE__,'Error add depth positive', fterr)
    end if
    istatus = nf90_def_var(ncid, 'ptop', nf90_double, varid=izvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var ptop', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'standard_name', 'air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ptop standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'long_name', 'Pressure at model top')
    call check_ok(__FILE__,__LINE__,'Error add ptop long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ptop units', fterr)
    istatus = nf90_def_var(ncid, 'iy', nf90_double, idims(2), ivvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var iy', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'standard_name', &
                           'projection_y_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add iy standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                           'y-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add iy long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add iy units', fterr)
    istatus = nf90_def_var(ncid, 'jx', nf90_double, idims(1), ivvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var jx', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                         'projection_x_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add jx standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                         'x-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add jx long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add jx units', fterr)
    istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), illtpvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var xlat', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', 'latitude')
    call check_ok(__FILE__,__LINE__,'Error add xlat standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                         'Latitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlat long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'units', 'degrees_north')
    call check_ok(__FILE__,__LINE__,'Error add xlat units', fterr)
    istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), illtpvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var xlon', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', 'longitude')
    call check_ok(__FILE__,__LINE__,'Error add xlon standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                         'Longitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlon long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'units', 'degrees_east')
    call check_ok(__FILE__,__LINE__,'Error add xlon units', fterr)
    istatus = nf90_def_var(ncid, 'topo', nf90_float, idims(1:2), illtpvar(3))
    call check_ok(__FILE__,__LINE__,'Error add var topo', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                         'surface_altitude')
    call check_ok(__FILE__,__LINE__,'Error add topo standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                         'Domain surface elevation')
    call check_ok(__FILE__,__LINE__,'Error add topo long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
    call check_ok(__FILE__,__LINE__,'Error add topo units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add topo coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add topo grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'mask', nf90_float, idims(1:2), illtpvar(4))
    call check_ok(__FILE__,__LINE__,'Error add var mask', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', 'landmask')
    call check_ok(__FILE__,__LINE__,'Error add mask standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                         'Domain land/ocean mask')
    call check_ok(__FILE__,__LINE__,'Error add mask long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add mask units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add mask coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add mask grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), itvar)
    call check_ok(__FILE__,__LINE__,'Error add var time', fterr)
    istatus = nf90_put_att(ncid, itvar, 'standard_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time standard_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'long_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time long_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'calendar', calstr(idate%calendar))
    call check_ok(__FILE__,__LINE__,'Error add time calendar', fterr)
    istatus = nf90_put_att(ncid, itvar, 'units', 'hours since '//ctime)
    call check_ok(__FILE__,__LINE__,'Error add time units', fterr)
    if (ctype == 'STS' .or. ctype == 'SRF') then
      istatus = nf90_put_att(ncid, itvar, 'bounds', 'time_bnds')
      call check_ok(__FILE__,__LINE__,'Error add time bounds', fterr)
    end if
    istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3), illtpvar(5))
    call check_ok(__FILE__,__LINE__,'Error add var ps', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'standard_name', &
                         'surface_air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'long_name', 'Surface pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ps units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add ps coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add ps grid_mapping', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'cell_methods', 'time: point')
    call check_ok(__FILE__,__LINE__,'Error add ps cell_methods', fterr)

    tyx = (/idims(1),idims(2),idims(3),-1,-1/)
    tzyx = (/idims(1),idims(2),idims(4),idims(3),-1/)
    t10yx = (/idims(1),idims(2),idims(5),idims(3),-1/)
    t2yx = (/idims(1),idims(2),idims(6),idims(3),-1/)
    tlyx = (/idims(1),idims(2),idims(7),idims(3),-1/)
    tcyx = (/idims(1),idims(2),idims(9),idims(3),-1/)
    tczyx = (/idims(1),idims(2),idims(4),idims(9),idims(3)/)
    tdyx = (/idims(1),idims(2),idims(10),idims(3),-1/)

    if (ctype == 'ATM') then
      iatmvar = -1
      iatmvar(1) = itvar
      iatmvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tzyx,.false.,3)
      call addvara(ncid,ctype,tzyx,.false.,4)
      call addvara(ncid,ctype,tzyx,.false.,5)
      call addvara(ncid,ctype,tzyx,.false.,6)
      call addvara(ncid,ctype,tzyx,.false.,7)
      call addvara(ncid,ctype,tzyx,.false.,8)
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call addvara(ncid,ctype,tzyx,.false.,9)
        call addvara(ncid,ctype,tzyx,.false.,10)
        call addvara(ncid,ctype,tzyx,.false.,11)
      end if
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
    else if (ctype == 'SRF') then
      isrfvar = -1
      isrfvar(1) = itvar
      istatus = nf90_def_var(ncid, 'time_bnds', nf90_double, &
                             (/idims(8),idims(3)/), isrfvar(2))
      call check_ok(__FILE__,__LINE__,'Error add var time_bnds', fterr)
      istatus = nf90_put_att(ncid, isrfvar(2), &
                             'calendar', calstr(idate%calendar))
      call check_ok(__FILE__,__LINE__,'Error add time_bnds calendar', fterr)
      istatus = nf90_put_att(ncid, isrfvar(2), 'units', 'hours since '//ctime)
      call check_ok(__FILE__,__LINE__,'Error add time_bnds units', fterr)
      isrfvar(3) = illtpvar(5)
      if ( iseaice /= 1 .and. lakemod /= 1 ) then
        srf_variables(26)%enabled = .false.
      end if
      call addvara(ncid,ctype,t10yx,.false.,4)
      call addvara(ncid,ctype,t10yx,.false.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.true.,8)
      call addvara(ncid,ctype,t2yx,.false.,9)
      call addvara(ncid,ctype,t2yx,.false.,10)
      call addvara(ncid,ctype,tlyx,.true.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.true.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
      call addvara(ncid,ctype,tyx,.false.,17)
      call addvara(ncid,ctype,tyx,.false.,18)
      call addvara(ncid,ctype,tyx,.false.,19)
      call addvara(ncid,ctype,tyx,.false.,20)
      call addvara(ncid,ctype,tyx,.false.,21)
      call addvara(ncid,ctype,tyx,.false.,22)
      call addvara(ncid,ctype,tyx,.false.,23)
      call addvara(ncid,ctype,tyx,.false.,24)
      call addvara(ncid,ctype,tyx,.false.,25)
      call addvara(ncid,ctype,tyx,.false.,26)
    else if (ctype == 'STS') then
      istsvar = -1
      istsvar(1) = itvar
      istatus = nf90_def_var(ncid, 'time_bnds', nf90_double, &
                             (/idims(8),idims(3)/), istsvar(2))
      call check_ok(__FILE__,__LINE__,'Error add var time_bnds', fterr)
      istatus = nf90_put_att(ncid, istsvar(2), &
                             'calendar', calstr(idate%calendar))
      call check_ok(__FILE__,__LINE__,'Error add time_bnds calendar', fterr)
      istatus = nf90_put_att(ncid, istsvar(2), 'units', 'hours since '//ctime)
      call check_ok(__FILE__,__LINE__,'Error add time_bnds units', fterr)
      istsvar(3) = illtpvar(5)
      call addvara(ncid,ctype,tyx,.false.,4)
      call addvara(ncid,ctype,tyx,.false.,5)
      call addvara(ncid,ctype,t2yx,.false.,6)
      call addvara(ncid,ctype,t2yx,.false.,7)
      call addvara(ncid,ctype,t2yx,.false.,8)
      call addvara(ncid,ctype,t10yx,.false.,9)
      call addvara(ncid,ctype,tyx,.false.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
    else if (ctype == 'SUB') then
      isubvar = -1
      isubvar(1) = itvar
      isubvar(2) = illtpvar(5)
      call addvara(ncid,ctype,t10yx,.false.,3)
      call addvara(ncid,ctype,t10yx,.false.,4)
      call addvara(ncid,ctype,tyx,.false.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.true.,7)
      call addvara(ncid,ctype,t2yx,.false.,8)
      call addvara(ncid,ctype,t2yx,.false.,9)
      call addvara(ncid,ctype,tlyx,.true.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.true.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.false.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
    else if (ctype == 'RAD') then
      iradvar = -1
      iradvar(1) = itvar
      iradvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tzyx,.false.,3)
      call addvara(ncid,ctype,tzyx,.false.,4)
      call addvara(ncid,ctype,tzyx,.false.,5)
      call addvara(ncid,ctype,tzyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.false.,8)
      call addvara(ncid,ctype,tyx,.false.,9)
      call addvara(ncid,ctype,tyx,.false.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.false.,14)
      call addvara(ncid,ctype,tyx,.false.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
      call addvara(ncid,ctype,tyx,.false.,17)
      call addvara(ncid,ctype,tyx,.false.,18)
    else if (ctype == 'LAK') then
      ilakvar = -1
      ilakvar(1) = itvar
      ilakvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tyx,.false.,3)
      call addvara(ncid,ctype,tyx,.false.,4)
      call addvara(ncid,ctype,tyx,.true.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.false.,8)
      call addvara(ncid,ctype,tyx,.false.,9)
      call addvara(ncid,ctype,tyx,.false.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.true.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.true.,15)
      call addvara(ncid,ctype,tdyx,.true.,16)
    end if

    istatus = nf90_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Error End Definitions NetCDF output',fterr)

    istatus = nf90_put_var(ncid, izvar(1), hsigma)
    call check_ok(__FILE__,__LINE__,'Error var sigma write', fterr)
    hptop = ptop*d_10
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call check_ok(__FILE__,__LINE__,'Error var ptop write', fterr)
    if (ctype == 'LAK') then
      do i = 1 , ndpmax
        depth(i) = real(i)
      end do
      istatus = nf90_put_var(ncid, izvar(3), depth)
      call check_ok(__FILE__,__LINE__,'Error var depth write', fterr)
    end if
    if (ctype == 'SUB') then
      yiy(1) = -((dble((o_nig-1)-1)/2.0D0)*ds)
      xjx(1) = -((dble((o_njg-1)-1)/2.0D0)*ds)
      do i = 2 , o_nig
        yiy(i) = yiy(i-1)+ds
      end do
      do j = 2 , o_njg
        xjx(j) = xjx(j-1)+ds
      end do
      istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_nig))
      call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
      istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_njg))
      call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(1), ioxlat_s)
      call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(2), ioxlon_s)
      call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(3), iotopo_s)
      call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(4), iomask_s)
      call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
    else
      yiy(1) = -((dble(o_ni-1)/2.0D0)*ds)
      xjx(1) = -((dble(o_nj-1)/2.0D0)*ds)
      do i = 2 , o_ni
        yiy(i) = yiy(i-1)+ds
      end do
      do j = 2 , o_nj
        xjx(j) = xjx(j-1)+ds
      end do
      istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_ni))
      call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
      istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_nj))
      call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
      call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
      call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
      call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(4), iomask)
      call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
    end if
    if (ctype == 'SRF' .or. ctype == 'SUB') then
      rdum1 = 10.0D0
      istatus = nf90_put_var(ncid, isrvvar(1), rdum1)
      call check_ok(__FILE__,__LINE__,'Error var m10 write', fterr)
      rdum1 = 2.0D0
      istatus = nf90_put_var(ncid, isrvvar(2), rdum1)
      call check_ok(__FILE__,__LINE__,'Error var m2 write', fterr)
      rdum2(1) = 0.0D0
      rdum2(2) = 1.0D0
      istatus = nf90_put_var(ncid, isrvvar(3), rdum2)
      call check_ok(__FILE__,__LINE__,'Error var layer write', fterr)
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncid)
      call check_ok(__FILE__,__LINE__,'Error initial sync', fterr)
    end if

    if (ctype == 'ATM') then
      ncatm = ncid
    else if (ctype == 'SRF') then
      ncsrf = ncid
    else if (ctype == 'SUB') then
      ncsub = ncid
    else if (ctype == 'STS') then
      ncsts = ncid
    else if (ctype == 'RAD') then
      ncrad = ncid
    else if (ctype == 'LAK') then
      nclak = ncid
    end if
  end subroutine prepare_common_out

  subroutine addvara(ncid,ctype,idims,lmiss,nvar)
    implicit none
    integer , intent(in) :: ncid
    character(3) , intent(in) :: ctype
    integer , dimension(5) , intent(in) :: idims
    logical , intent(in) :: lmiss
    integer , intent(in) :: nvar

    character(len=8)   :: vname
    character(len=128) :: vst , vln
    character(len=16)  :: vuni , vmeth
    logical :: lreq
    integer :: ivar
    character(64) :: cmethodpnt , cmethodmax , cmethodmin
    character(64) :: cmethodsum , cmethodmean
    character(len=128) :: cdum

    integer :: i , ndims

    ndims = 0
    do i = 1 , 5
      if (idims(i) > 0) ndims = ndims+1
    end do

    if ( ctype == 'STS' ) then
      write (cmethodpnt,  '(a)') 'time: point'
      write (cmethodmax,  '(a)') 'time: maximum (interval: 1 day)'
      write (cmethodmin,  '(a)') 'time: minimum (interval: 1 day)'
      write (cmethodmean, '(a)') 'time: mean (interval: 1 day)'
      write (cmethodsum,  '(a)') 'time: sum (interval: 1 day)'
    else
      write (cmethodpnt,  '(a)') 'time: point'
      write (cmethodmax,  '(a,i2,a)') 'time: maximum (interval: ', & 
                           idnint(srffrq), ' hours)'
      write (cmethodmin,  '(a,i2,a)') 'time: minimum (interval: ', &
                           idnint(srffrq), ' hours)'
      write (cmethodmean, '(a,i2,a)') 'time: mean (interval: ', &
                           idnint(srffrq), ' hours)'
      write (cmethodsum,  '(a,i2,a)') 'time: sum (interval: ', &
                           idnint(srffrq), ' hours)'
    end if

    select case (ctype)
      case ('ATM')
        vname = atm_variables(nvar)%vname
        vst   = atm_variables(nvar)%vstd_name
        vln   = atm_variables(nvar)%vdesc
        vuni  = atm_variables(nvar)%vunit
        vmeth = atm_variables(nvar)%time_meth
        lreq  = atm_variables(nvar)%enabled
      case ('SRF')
        vname = srf_variables(nvar)%vname
        vst   = srf_variables(nvar)%vstd_name
        vln   = srf_variables(nvar)%vdesc
        vuni  = srf_variables(nvar)%vunit
        vmeth = srf_variables(nvar)%time_meth
        lreq  = srf_variables(nvar)%enabled
      case ('SUB')
        vname = sub_variables(nvar)%vname
        vst   = sub_variables(nvar)%vstd_name
        vln   = sub_variables(nvar)%vdesc
        vuni  = sub_variables(nvar)%vunit
        vmeth = sub_variables(nvar)%time_meth
        lreq  = sub_variables(nvar)%enabled
      case ('STS')
        vname = sts_variables(nvar)%vname
        vst   = sts_variables(nvar)%vstd_name
        vln   = sts_variables(nvar)%vdesc
        vuni  = sts_variables(nvar)%vunit
        vmeth = sts_variables(nvar)%time_meth
        lreq  = sts_variables(nvar)%enabled
      case ('LAK')
        vname = lak_variables(nvar)%vname
        vst   = lak_variables(nvar)%vstd_name
        vln   = lak_variables(nvar)%vdesc
        vuni  = lak_variables(nvar)%vunit
        vmeth = lak_variables(nvar)%time_meth
        lreq  = lak_variables(nvar)%enabled
      case ('RAD')
        vname = rad_variables(nvar)%vname
        vst   = rad_variables(nvar)%vstd_name
        vln   = rad_variables(nvar)%vdesc
        vuni  = rad_variables(nvar)%vunit
        vmeth = rad_variables(nvar)%time_meth
        lreq  = rad_variables(nvar)%enabled
      case default
        call fatal(__FILE__,__LINE__,ctype//': Not defined output type')
    end select

    if (lreq) then
      cdum = vname
      istatus = nf90_def_var(ncid, cdum, nf90_float, idims(1:ndims), ivar)
      call check_ok(__FILE__,__LINE__,'Error add var '//vname, ctype//' FILE')
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
      call check_ok(__FILE__,__LINE__,'Error setting deflate on var '//vname, &
                    ctype//' FILE')
#endif
      cdum = vst
      istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' standard_name', &
                    ctype//' FILE')
      cdum = vln
      istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' long_name', &
                    ctype//' FILE')
      cdum = vuni
      istatus = nf90_put_att(ncid, ivar, 'units', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' units', &
                    ctype//' FILE')
      istatus = nf90_put_att(ncid, ivar, 'coordinates', 'xlat xlon')
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' coord', &
                    ctype//' FILE')
      istatus = nf90_put_att(ncid, ivar, 'grid_mapping', 'rcm_map')
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' grid_mapping', &
                    ctype//' FILE')
      select case (vmeth)
        case ( 'point' )
          istatus = nf90_put_att(ncid, ivar, 'cell_methods', cmethodpnt)
          call check_ok(__FILE__,__LINE__,'Error add '//vname//' cell_methods',&
                        ctype//' FILE')
        case ( 'maximum' )
          istatus = nf90_put_att(ncid, ivar, 'cell_methods', cmethodmax)
          call check_ok(__FILE__,__LINE__,'Error add '//vname//' cell_methods',&
                        ctype//' FILE')
        case ( 'minimum' )
          istatus = nf90_put_att(ncid, ivar, 'cell_methods', cmethodmin)
          call check_ok(__FILE__,__LINE__,'Error add '//vname//' cell_methods',&
                        ctype//' FILE')
        case ( 'mean' )
          istatus = nf90_put_att(ncid, ivar, 'cell_methods', cmethodmean)
          call check_ok(__FILE__,__LINE__,'Error add '//vname//' cell_methods',&
                        ctype//' FILE')
        case ( 'sum' )
          istatus = nf90_put_att(ncid, ivar, 'cell_methods', cmethodsum)
          call check_ok(__FILE__,__LINE__,'Error add '//vname//' cell_methods',&
                        ctype//' FILE')
      end select
      if (lmiss) then
        istatus = nf90_put_att(ncid, ivar, '_FillValue', smissval)
        call check_ok(__FILE__,__LINE__,'Error add '//vname//' coord', &
                      ctype//' FILE')
      end if
      select case (ctype)
        case ('ATM')
          iatmvar(nvar) = ivar
        case ('SRF')
          isrfvar(nvar) = ivar
        case ('SUB')
          isubvar(nvar) = ivar
        case ('STS')
          istsvar(nvar) = ivar
        case ('LAK')
          ilakvar(nvar) = ivar
        case ('RAD')
          iradvar(nvar) = ivar
        case default
          call fatal(__FILE__,__LINE__,ctype//': Not defined output type')
      end select
    end if
  end subroutine addvara

  subroutine writerec_sts(fbat, idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: fbat
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(2) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    ctime = tochar(idate)

    istart(2) = istsrec
    istart(1) = 1
    icount(2) = 1
    icount(1) = 2
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)-12
    istatus = nf90_put_var(ncsts, istsvar(1), nctime(1:1), &
                           istart(2:2), icount(2:2))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'STS FILE')
    nctime(2) = tohours(tdif)
    nctime(1) = nctime(2)-24
    istatus = nf90_put_var(ncsts, istsvar(2), nctime, &
                           istart(1:2), icount(1:2))
    call check_ok(__FILE__,__LINE__,  &
                  'Error writing time_bnds '//ctime, 'STS FILE')

    istart(3) = istsrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncsts, istsvar(3), & 
                           fbat(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//sts_variables(3)%vname// &
                  ' at '//ctime, 'STS FILE')
    ivar = 4
    do n = 25 , numbat
      if ( sts_variables(ivar)%enabled ) then
        if (ivar == ivarname_lookup('STS', 'w10max') .or. &
            ivar == ivarname_lookup('STS', 't2avg')  .or. &
            ivar == ivarname_lookup('STS', 't2max')  .or. &
            ivar == ivarname_lookup('STS', 't2min')) then
          istart(4) = istsrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsts, istsvar(ivar), &
                                 fbat(:,:,n), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sts_variables(ivar)%vname// &
                        ' at '//ctime, 'STS FILE')
        else
          istart(3) = istsrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsts, istsvar(ivar), & 
                                 fbat(:,:,n), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sts_variables(ivar)%vname// &
                        ' at '//ctime, 'STS FILE')
        end if
      end if
      ivar = ivar + 1
    end do

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncsts)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'STS FILE')
    end if
    istsrec = istsrec + 1
  end subroutine writerec_sts

  subroutine writerec_srf(fbat, mask , idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: fbat
    integer , pointer , dimension(:,:) , intent(in) :: mask
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(2) :: nctime
    type(rcm_time_interval) :: tdif
    logical :: lskip
    character(len=36) :: ctime

    ctime = tochar(idate)

    istart(2) = isrfrec
    istart(1) = 1
    icount(2) = 1
    icount(1) = 2
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncsrf, isrfvar(1), nctime(1:1), &
                           istart(2:2), icount(2:2))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'SRF FILE')
    nctime(2) = tohours(tdif)
    nctime(1) = nctime(2)-srffrq
    istatus = nf90_put_var(ncsrf, isrfvar(2), nctime, istart(1:2), icount(1:2))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing time_bnds '//ctime, 'SRF FILE')
    istart(3) = isrfrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncsrf, isrfvar(3), fbat(:,:,1), istart(1:3), &
                           icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'SRF FILE')

    ivar = 4
    lskip = .false.
    do n = 2 , 24
      if (lskip) then
        lskip = .false.
        cycle
      end if
      if ( srf_variables(ivar)%enabled ) then
        if (ivar == ivarname_lookup('SRF', 'u10m')   .or. &
            ivar == ivarname_lookup('SRF', 'v10m')   .or. &
            ivar == ivarname_lookup('SRF', 't2m')    .or. &
            ivar == ivarname_lookup('SRF', 'q2m')) then
          istart(4) = isrfrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), &
                          fbat(:,:,n), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        else if (ivar == ivarname_lookup('SRF', 'smw')) then
          istart(4) = isrfrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                          fbat(:,:,n), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
          istart(3) = 2
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                          fbat(:,:,n+1), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        else
          istart(3) = isrfrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                   fbat(:,:,n), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        end if
      end if
      if (ivar == ivarname_lookup('SRF', 'smw')) then
        lskip = .true.
      end if
      ivar = ivar + 1
    end do

    if ( srf_variables(26)%enabled ) then
      dumio(:,:,1) = real(mask(o_js:o_je,o_is:o_ie))
      istart(3) = isrfrec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      istatus = nf90_put_var(ncsrf, isrfvar(26), & 
               dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//srf_variables(26)%vname// &
                    ' at '//ctime, 'SRF FILE')
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncsrf)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'SRF FILE')
    end if
    isrfrec = isrfrec + 1
  end subroutine writerec_srf

  subroutine writerec_sub(fsub, idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: fsub
    integer :: ivar
    integer :: n , nxb , nyb
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime
    logical :: lskip

    nxb = o_njg / nsg
    nyb = o_nig / nsg

    ctime = tochar(idate)

    istart(1) = isubrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncsub, isubvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'SUB FILE')
    ivar = 2
    lskip = .false.
    do n = 1 , numsub
      if (lskip) then
        lskip = .false.
        cycle
      end if
      if ( sub_variables(ivar)%enabled ) then
        call reorder_3_2(fsub,subio,n)
        if (ivar == ivarname_lookup('SUB', 'u10m')   .or. &
            ivar == ivarname_lookup('SUB', 'v10m')   .or. &
            ivar == ivarname_lookup('SUB', 't2m')    .or. &
            ivar == ivarname_lookup('SUB', 'q2m') ) then
          istart(4) = isubrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        else if (ivar == ivarname_lookup('SUB', 'smw')) then
          istart(4) = isubrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
          istart(3) = 2
          call reorder_3_2(fsub,subio,n+1)
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        else
          istart(3) = isubrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), & 
                                 subio, istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        end if
      end if
      if (ivar == ivarname_lookup('SUB', 'smw')) then
        lskip = .true.
      end if
      ivar = ivar + 1
    end do
    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncsub)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'SUB FILE')
    end if
    isubrec = isubrec + 1
  end subroutine writerec_sub

  subroutine writerec_rad(nrad3d,nrad2d,frad3d,frad2d,ps,idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nrad3d , nrad2d
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: frad3d
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: frad2d
    real(dp) , pointer , dimension(:,:) , intent(in) :: ps
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    ctime = tochar(idate)

    istart(1) = iradrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncrad, iradvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'RAD FILE')

    istart(3) = iradrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    dumio(:,:,1) = real((ps(o_js:o_je,o_is:o_ie)+ptop)*d_10)
    istatus = nf90_put_var(ncrad, iradvar(2), dumio, istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'RAD FILE')

    ivar = 3
    do n = 1 , nrad3d
      if ( rad_variables(ivar)%enabled ) then
        istart(4) = iradrec
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = o_nz
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncrad, iradvar(ivar), &
                               frad3d(:,:,:,n), istart, icount)
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//rad_variables(ivar)%vname// &
                      ' at '//ctime, 'RAD FILE')
      end if
      ivar = ivar + 1
    end do
    do n = 1 , nrad2d
      if ( rad_variables(ivar)%enabled ) then
        istart(3) = iradrec
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncrad, iradvar(ivar), & 
                               frad2d(:,:,n), istart(1:3), icount(1:3))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//rad_variables(ivar)%vname// &
                      ' at '//ctime, 'RAD FILE')
      end if
      ivar = ivar + 1
    end do

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncrad)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'RAD FILE')
    end if
    iradrec = iradrec + 1
  end subroutine writerec_rad

  subroutine writerec_atm(u,v,omega,t,qx,tke,kth,kzm,ps,rc,rnc, &
                          tgb,swt,mask,idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: u
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: v
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: omega
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: t
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: qx
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: tke
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: kth
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: kzm
    real(dp) , pointer , dimension(:,:) , intent(in) :: ps
    real(dp) , pointer , dimension(:,:) , intent(in) :: rc
    real(dp) , pointer , dimension(:,:) , intent(in) :: rnc
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: tgb
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: swt
    integer , pointer , dimension(:,:,:) , intent(in) :: mask
    integer :: i , j , n , ip1 , ip2 , jp1 , jp2 , k
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    ctime = tochar(idate)

    if (.not. lmaskfill) then
      do n = 1 , nnsg
        atmsrfmask(n,:,:) = real(mask(n,o_js:o_je,o_is:o_ie))
      end do
      where ( atmsrfmask > 0.5 )
        atmsrfmask = 1.0
      elsewhere
        atmsrfmask = 0.0
      end where
      atmsrfsum = sum(atmsrfmask,dim=1)
      lmaskfill = .true.
    end if

    istart(1) = iatmrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncatm, iatmvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'ATM FILE')

    dumio(:,:,1) = real((ps(o_js:o_je,o_is:o_ie) + ptop)*d_10)
    istart(3) = iatmrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncatm, iatmvar(2), &
                           dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'ATM FILE')

    istart(4) = iatmrec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = o_nz
    icount(2) = o_ni
    icount(1) = o_nj

    if ( atm_variables(3)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          ip2 = i+2
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
              jp2 = j+2
            else
              jp1 = j
              jp2 = j+1
              if (j == o_nj) jp2 = 1
            end if
            dumio(j,i,k) = real(((u(jp1,ip1,k)+u(jp2,ip1,k) + &
                                  u(jp1,ip2,k)+u(jp2,ip2,k))*d_rfour) / &
                                  ps(jp1,ip1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(3), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//atm_variables(3)%vname// &
                    ' at '//ctime, 'ATM FILE')
    end if

    if ( atm_variables(4)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          ip2 = i+2
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
              jp2 = j+2
            else
              jp1 = j
              jp2 = j+1
              if (j == o_nj) jp2 = 1
            end if
            dumio(j,i,k) = real(((v(jp1,ip1,k)+v(jp2,ip1,k) + &
                                  v(jp1,ip2,k)+v(jp2,ip2,k))*d_rfour) / &
                                  ps(jp1,ip1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(4), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(4)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(5)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            dumio(j,i,k) = real(omega(jp1,ip1,k)*d_10)
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(5), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(5)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(6)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            dumio(j,i,k) = real(t(jp1,ip1,k)/ps(jp1,ip1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(6), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(6)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(7)%enabled ) then
      dumio = 0.0
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            if (qx(jp1,ip1,k,iqv) > dlowval) then
              dumio(j,i,k) = real(qx(jp1,ip1,k,iqv)/ps(jp1,ip1))
            end if
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(7), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(7)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(8)%enabled ) then
      dumio = 0.0
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            if (qx(jp1,ip1,k,iqc) > dlowval) then
              dumio(j,i,k) = real(qx(jp1,ip1,k,iqc)/ps(jp1,ip1))
            end if
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(8), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(8)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      if ( atm_variables(9)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (tke(jp1,ip1,k) > dlowval) then
                dumio(j,i,k) = real(tke(jp1,ip1,k))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(9), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                      'Error writing '//atm_variables(9)%vname//' at '//ctime, &
                      'ATM FILE')
      end if
      if ( atm_variables(10)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (kth(jp1,ip1,k) > dlowval) then
                dumio(j,i,k) = real(kth(jp1,ip1,k))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(10), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                'Error writing '//atm_variables(10)%vname//' at '//ctime, &
                'ATM FILE')
      end if
      if ( atm_variables(11)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (kzm(jp1,ip1,k) > dlowval) then
                dumio(j,i,k) = real(kzm(jp1,ip1,k))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(11), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                'Error writing '//atm_variables(11)%vname//' at '//ctime, &
                'ATM FILE')
      end if
    end if

    istart(3) = iatmrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj

    if ( atm_variables(12)%enabled ) then
      dumio(:,:,1) = 0.0
      where (rc(o_js:o_je,o_is:o_ie) > dlowval)
        dumio(:,:,1) = real(rc(o_js:o_je,o_is:o_ie))
      end where
      where (rnc(o_js:o_je,o_is:o_ie) > dlowval)
        dumio(:,:,1) = dumio(:,:,1) + real(rnc(o_js:o_je,o_is:o_ie))
      end where
      dumio(:,:,1) = dumio(:,:,1)*real(tpd)
      istatus = nf90_put_var(ncatm, iatmvar(12), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(12)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(13)%enabled ) then
      dumio(:,:,1) = real(sum(tgb(:,o_js:o_je,o_is:o_ie), dim=1)*xns2d)
      istatus = nf90_put_var(ncatm, iatmvar(13), & 
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(13)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(14)%enabled ) then
      dumio(:,:,1) = 0.0
      do n = 1 , nnsg
        where (atmsrfmask(n,:,:) > 0)
          dumio(:,:,1) = dumio(:,:,1) + real(swt(n,o_js:o_je,o_is:o_ie))
        end where
      end do
      where (atmsrfsum > 0)
        dumio(:,:,1) = dumio(:,:,1)/amax1(atmsrfsum/2.0,1.0)
      elsewhere
        dumio(:,:,1) = smissval
      end where
      istatus = nf90_put_var(ncatm, iatmvar(14), & 
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//atm_variables(14)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncatm)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'ATM FILE')
    end if
    iatmrec = iatmrec + 1
  end subroutine writerec_atm

  subroutine writerec_lak(fbat,evl,aveice,hsnow,tlake,idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: fbat
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: evl
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: aveice
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: hsnow
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: tlake
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    ctime = tochar(idate)

    istart(1) = ilakrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(nclak, ilakvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'LAK FILE')

    ivar = 2
    do n = 1 , numbat
      if (lak_fbats(n) == 0) cycle
      istart(3) = ilakrec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                      fbat(:,:,n), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//lak_variables(ivar)%vname// &
                    ' at '//ctime, 'LAK FILE')
      ivar = ivar + 1
    end do

    ! Add lake model output
    dumio(:,:,1) =  real(sum(evl(:,o_js:o_je,o_is:o_ie),1)*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                    dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    dumio(:,:,1) = real(sum(aveice(:,o_js:o_je,o_is:o_ie),1)*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                    dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    dumio(:,:,1) = real(sum(hsnow(:,o_js:o_je,o_is:o_ie),1)*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
             dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    do n = 1 , ndpmax
      istart(4) = ilakrec
      istart(3) = n
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      dumio(:,:,1) = real(sum(tlake(:,o_js:o_je,o_is:o_ie,n),1)*xns2d)
      where (iolnds == 14)
        dumio(:,:,1) = dumio(:,:,1) + real(tzero)
      end where
      istatus = nf90_put_var(nclak, ilakvar(ivar), dumio(:,:,1), istart, icount)
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    end do

    if ( debug_level > 2 ) then
      istatus = nf90_sync(nclak)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'LAK FILE')
    end if
    ilakrec = ilakrec + 1
  end subroutine writerec_lak

  subroutine check_ok(f,l,m1,mf)
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer , intent(in) :: l
    if (istatus /= nf90_noerr) then 
      write (6,*) trim(m1)
      write (6,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

  subroutine release_mod_ncio
    implicit none
    call close_domain
    call close_icbc
    call close_common(ncatm,'ATM')
    call close_common(ncsrf,'SRF')
    call close_common(ncsts,'STS')
    call close_common(ncsub,'SUB')
    call close_common(ncrad,'RAD')
    call close_common(nclak,'LAK')
  end subroutine release_mod_ncio

  subroutine reorder_2_3(m1,m2)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: m1
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: m2
    integer :: i , j , ii , jj , n
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        m2(n,jj,ii) = dble(m1(j,i))
      end do
    end do
  end subroutine reorder_2_3

  subroutine reorder_3_2(m2,m1,m)
    implicit none
    integer , intent(in) :: m
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: m2
    real(sp) , pointer , dimension(:,:) , intent(out) :: m1
    integer :: i , ii , j , jj , n
    do j = 1 , o_njg
      do i = 1 , o_nig
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        m1(j,i) = m2(n,jj+1,ii+1,m)
      end do
    end do
  end subroutine reorder_3_2

end module mod_ncio
