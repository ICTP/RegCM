module mod_clm_balancecheck
  !
  ! Water and energy balance check.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_runparams , only : ktau , dtsrf
  use mod_clm_type , only : clm3
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varpar , only : nlevgrnd , nlevsoi , nlevurb
  use mod_clm_varcon , only : icol_roof , icol_sunwall , icol_shadewall
  use mod_clm_varcon , only : icol_road_perv , icol_road_imperv
  use mod_clm_varcon , only : isturb , spval , istdlak
  use mod_clm_varcon , only : istslak , istsoil , istcrop , istwet
  use mod_clm_subgridave

  implicit none

  private

  save

  public :: BeginWaterBalance  ! Initialize water balance check
  public :: BalanceCheck       ! Water and energy balance check

  contains
  !
  ! Initialize column-level water balance at beginning of time step
  !
  subroutine BeginWaterBalance(lbc,ubc,num_nolakec,filter_nolakec, &
                  num_lakec,filter_lakec,num_hydrologyc,filter_hydrologyc)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! column-index bounds
    ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: num_nolakec
    ! column filter for non-lake points
    integer(ik4) , intent(in) , dimension(ubc-lbc+1) :: filter_nolakec
    ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: num_lakec
    ! column filter for non-lake points
    integer(ik4) , intent(in) , dimension(ubc-lbc+1) :: filter_lakec
    ! number of column soil points in column filter
    integer(ik4) , intent(in) :: num_hydrologyc
    ! column filter for soil points
    integer(ik4) , intent(in) , dimension(ubc-lbc+1) :: filter_hydrologyc
    real(rk8) , pointer , dimension(:) :: h2osfc       ! surface water (mm)
    real(rk8) , pointer , dimension(:) :: londeg       ! longitude
    real(rk8) , pointer , dimension(:) :: latdeg       ! latitude
    integer(ik4) , pointer , dimension(:) :: cgridcell ! column's gridcell index
    integer(ik4) , pointer , dimension(:) :: clandunit ! column's landunit
    integer(ik4) , pointer , dimension(:) :: ltype     ! landunit type
    real(rk8) , pointer , dimension(:) :: h2osno       ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice ! ice lens (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq ! liquid water (kg/m2)
    ! canopy water (mm H2O) (pft-level)
    real(rk8) , pointer , dimension(:) :: h2ocan_pft
    ! water in the unconfined aquifer (mm)
    real(rk8) , pointer , dimension(:) :: wa
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    real(rk8) , pointer , dimension(:) :: zwt ! water table depth (m)
    ! interface level below a "z" level (m)
    real(rk8) , pointer , dimension(:,:) :: zi
    ! canopy water (mm H2O) (column level)
    real(rk8) , pointer , dimension(:) :: h2ocan_col
    ! water mass begining of the time step
    real(rk8) , pointer , dimension(:) :: begwb
    integer(ik4) :: c , f , j ! indices
    real(rk8) , pointer , dimension(:,:) :: dz , watsat

    ! Assign local pointers to derived type members (column-level)

    h2osfc             => clm3%g%l%c%cws%h2osfc
    londeg             => clm3%g%londeg
    latdeg             => clm3%g%latdeg
    cgridcell          => clm3%g%l%c%gridcell
    clandunit          => clm3%g%l%c%landunit
    ltype              => clm3%g%l%itype
    dz                 => clm3%g%l%c%cps%dz
    watsat             => clm3%g%l%c%cps%watsat
    h2osno             => clm3%g%l%c%cws%h2osno
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    begwb              => clm3%g%l%c%cwbal%begwb
    h2ocan_col         => clm3%g%l%c%cws%pws_a%h2ocan
    wa                 => clm3%g%l%c%cws%wa
    ctype              => clm3%g%l%c%itype
    zwt                => clm3%g%l%c%cws%zwt
    zi                 => clm3%g%l%c%cps%zi

    ! Assign local pointers to derived type members (pft-level)

    h2ocan_pft         => clm3%g%l%c%p%pws%h2ocan

    ! Determine beginning water balance for time step
    ! pft-level canopy water averaged to column

    call p2c(num_nolakec,filter_nolakec,h2ocan_pft,h2ocan_col)

    do f = 1 , num_hydrologyc
      c = filter_hydrologyc(f)
      if ( zwt(c) <= zi(c,nlevsoi) ) then
        wa(c) = 5000.D0
      end if
    end do

    do f = 1 , num_nolakec
      c = filter_nolakec(f)
      if ( ctype(c) == icol_roof .or.      &
           ctype(c) == icol_sunwall .or.   &
           ctype(c) == icol_shadewall .or. &
           ctype(c) == icol_road_imperv ) then
        begwb(c) = h2ocan_col(c) + h2osno(c)
      else
        begwb(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
      end if
    end do
    do j = 1 , nlevgrnd
      do f = 1 , num_nolakec
        c = filter_nolakec(f)
        if ( (ctype(c) == icol_sunwall .or.   &
              ctype(c) == icol_shadewall .or. &
              ctype(c) == icol_roof) .and. j > nlevurb ) then
          continue
        else
          begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
        end if
      end do
    end do

    do f = 1 , num_lakec
      c = filter_lakec(f)
      begwb(c) = h2osno(c)
    end do
  end subroutine BeginWaterBalance
  !
  ! This subroutine accumulates the numerical truncation errors of the water
  ! and energy balance calculation. It is helpful to see the performance of
  ! the process of integration.
  !
  ! The error for energy balance:
  !
  ! error = abs(Net radiation - change of internal energy - Sensible heat
  !             - Latent heat)
  !
  ! The error for water balance:
  !
  ! error = abs(precipitation - change of water storage - evaporation - runoff)
  !
  subroutine BalanceCheck(lbp,ubp,lbc,ubc,lbl,ubl,lbg,ubg)
    implicit none
    integer(ik4) :: lbp , ubp ! pft-index bounds
    integer(ik4) :: lbc , ubc ! column-index bounds
    integer(ik4) :: lbl , ubl ! landunit-index bounds
    integer(ik4) :: lbg , ubg ! grid-index bounds
    real(rk8) , pointer , dimension(:) :: tws  !total water storage (mm H2O)
    real(rk8) , pointer , dimension(:) :: volr !river water storage (m3)
    real(rk8) , pointer , dimension(:) :: area !gridcell area (km2)
    ! true => do snow capping
    logical , pointer , dimension(:) :: do_capsnow
    ! rain on ground after interception (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_rain_grnd_col
    ! snow on ground after interception (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_snow_grnd_col
    ! snow falling on surface water (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_snow_h2osfc
    ! effective snow fraction
    real(rk8) , pointer , dimension(:) :: frac_sno_eff
    ! conversion of h2osfc to ice
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_to_ice
    ! snow melt (net)
    real(rk8) , pointer , dimension(:) :: qflx_snow_melt
    ! fraction of ground covered by snow (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_sno
    ! sub-surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_drain_perched
    ! total runoff due to flooding
    real(rk8) , pointer , dimension(:) :: qflx_floodc
    ! soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_evap_soi
    !surface water runoff (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_surf
    ! solar radiation absorbed by soil (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg_soil
    ! solar radiation absorbed by snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg_snow
    ! sum of soil/snow using current fsno, for balance check
    real(rk8) , pointer , dimension(:) :: sabg_chk
    ! pft's column index
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: cactive
    ! pft's gridcell index
    integer(ik4) , pointer , dimension(:) :: pgridcell
    ! pft's landunit index
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! column's gridcell index
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! column's landunit index
    integer(ik4) , pointer , dimension(:) :: clandunit
    integer(ik4) , pointer , dimension(:) :: ltype  ! landunit type
    integer(ik4) , pointer , dimension(:) :: ctype  ! column type
    real(rk8) , pointer , dimension(:) :: forc_rain ! rain rate [mm/s]
    real(rk8) , pointer , dimension(:) :: forc_snow ! snow rate [mm/s]
    ! downward infrared (longwave) radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: forc_lwrad
    ! water mass end of the time step
    real(rk8) , pointer , dimension(:) :: endwb
    ! water mass begining of the time step
    real(rk8) , pointer , dimension(:) :: begwb
    ! solar radiation absorbed (total) (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsa
    ! solar radiation reflected (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr
    ! emitted infrared (longwave) radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_out
    ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_net
    ! solar radiation absorbed by vegetation (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabv
    ! solar radiation absorbed by ground (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg
    ! total sensible heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_tot
    ! total sensible heat flux at grid level (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_totg
    ! energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_dynbal
    ! total latent heat flux (W/m8*2)  [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_tot
    ! soil heat flux (W/m**2) [+ = into soil]
    real(rk8) , pointer , dimension(:) :: eflx_soil_grnd
    ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8) , pointer , dimension(:) :: qflx_evap_tot
    ! irrigation flux (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_irrig
    ! surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_surf
    ! qflx_surf at glaciers, wetlands, lakes
    real(rk8) , pointer , dimension(:) :: qflx_qrgwl
    ! sub-surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_drain
    ! total runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_runoff
    ! total runoff at gridcell level inc land cover change flux (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_runoffg
    ! liq runoff due to dynamic land cover change (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_liq_dynbal
    ! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_ice
    ! excess snowfall due to snow cap inc land cover change flux (mm H20/s)
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_iceg
    ! ice runoff due to dynamic land cover change (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_ice_dynbal
    ! direct beam radiation (vis=forc_sols , nir=forc_soll )
    real(rk8) , pointer , dimension(:,:) :: forc_solad
    ! diffuse radiation     (vis=forc_solsd, nir=forc_solld)
    real(rk8) , pointer , dimension(:,:) :: forc_solai
    ! traffic sensible heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_traffic_pft
    ! sensible heat flux from urban heating/cooling sources
    !   of waste heat (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_wasteheat_pft
    ! ratio of building height to street width
    real(rk8) , pointer , dimension(:) :: canyon_hwr
    ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_heat_from_ac_pft
    ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    ! snow water (mm H2O) at previous time step
    real(rk8) , pointer , dimension(:) :: h2osno_old
    ! surface dew added to snow pack (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_dew_snow
    ! sublimation rate from snow pack (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_sub_snow
    ! net water input into soil from top (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_top_soil
    ! ground surface dew formation (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_dew_grnd
    ! ground surface evaporation rate (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_evap_grnd
    ! water onto ground including canopy runoff [kg/(m2 s)]
    real(rk8) , pointer , dimension(:) :: qflx_prec_grnd
    ! excess liquid water due to snow capping (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_liq
    ! liquid water + ice from layer above soil to top soil layer or
    !  sent to qflx_qrgwl (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_sl_top_soil
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! water conservation error (mm H2O)
    real(rk8) , pointer , dimension(:) :: errh2o
    ! solar radiation conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errsol
    ! longwave radiation conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errlon
    ! surface energy conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errseb
    ! net radiation (positive downward) (W/m**2)
    real(rk8) , pointer , dimension(:) :: netrad
    ! column-level soil/lake energy conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errsoi_col
    ! snow sources (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: snow_sources
    ! snow sinks (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: snow_sinks
    ! error in h2osno (kg m-2)
    real(rk8) , pointer , dimension(:) :: errh2osno
    integer(ik4) :: p , c , l , g ! indices
    logical :: found                       ! flag in search loop
    ! index of first found in search loop
    integer(ik4) :: indexp , indexc , indexg
    ! column level rain rate [mm/s]
    real(rk8) , dimension(lbc:ubc) :: forc_rain_col
    ! column level snow rate [mm/s]
    real(rk8) , dimension(lbc:ubc) :: forc_snow_col

    ! Assign local pointers to derived type scalar members (gridcell-level)

    tws                 => clm3%g%tws
    area                => clm3%g%area
    volr                => clm_a2l%volr
    do_capsnow          => clm3%g%l%c%cps%do_capsnow
    qflx_rain_grnd_col  => clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd
    qflx_snow_grnd_col  => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd
    qflx_snow_h2osfc    => clm3%g%l%c%cwf%qflx_snow_h2osfc
    frac_sno_eff        => clm3%g%l%c%cps%frac_sno_eff
    qflx_h2osfc_to_ice  => clm3%g%l%c%cwf%qflx_h2osfc_to_ice
    frac_sno            => clm3%g%l%c%cps%frac_sno
    qflx_drain_perched  => clm3%g%l%c%cwf%qflx_drain_perched
    qflx_floodc         => clm3%g%l%c%cwf%qflx_floodc
    qflx_evap_soi       => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    qflx_h2osfc_surf    => clm3%g%l%c%cwf%qflx_h2osfc_surf
    qflx_snow_melt      => clm3%g%l%c%cwf%qflx_snow_melt
    sabg_soil           => clm3%g%l%c%p%pef%sabg_soil
    sabg_snow           => clm3%g%l%c%p%pef%sabg_snow
    sabg_chk            => clm3%g%l%c%p%pef%sabg_chk
    pcolumn             => clm3%g%l%c%p%column
    forc_rain           => clm_a2l%forc_rain
    forc_snow           => clm_a2l%forc_snow
    forc_lwrad          => clm_a2l%forc_lwrad
    forc_solad          => clm_a2l%forc_solad
    forc_solai          => clm_a2l%forc_solai

    ! Assign local pointers to derived type scalar members (landunit-level)

    ltype               => clm3%g%l%itype
    canyon_hwr          => clm3%g%l%canyon_hwr

    ! Assign local pointers to derived type scalar members (column-level)

    cactive             => clm3%g%l%c%active
    ctype               => clm3%g%l%c%itype
    cgridcell           => clm3%g%l%c%gridcell
    clandunit           => clm3%g%l%c%landunit
    endwb               => clm3%g%l%c%cwbal%endwb
    begwb               => clm3%g%l%c%cwbal%begwb
    qflx_irrig          => clm3%g%l%c%cwf%qflx_irrig
    qflx_surf           => clm3%g%l%c%cwf%qflx_surf
    qflx_qrgwl          => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_drain          => clm3%g%l%c%cwf%qflx_drain
    qflx_runoff         => clm3%g%l%c%cwf%qflx_runoff
    qflx_snwcp_ice      => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    qflx_evap_tot       => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    errh2o              => clm3%g%l%c%cwbal%errh2o
    errsoi_col          => clm3%g%l%c%cebal%errsoi
    h2osno              => clm3%g%l%c%cws%h2osno
    h2osno_old          => clm3%g%l%c%cws%h2osno_old
    qflx_dew_snow       => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow       => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_top_soil       => clm3%g%l%c%cwf%qflx_top_soil
    qflx_evap_grnd      => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    qflx_dew_grnd       => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_prec_grnd      => clm3%g%l%c%cwf%pwf_a%qflx_prec_grnd
    qflx_snwcp_liq      => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_liq
    qflx_sl_top_soil    => clm3%g%l%c%cwf%qflx_sl_top_soil
    snow_sources        => clm3%g%l%c%cws%snow_sources
    snow_sinks          => clm3%g%l%c%cws%snow_sinks
    errh2osno           => clm3%g%l%c%cws%errh2osno
    snl                 => clm3%g%l%c%cps%snl

    ! Assign local pointers to derived type scalar members (pft-level)

    pactive             => clm3%g%l%c%p%active
    pgridcell           => clm3%g%l%c%p%gridcell
    plandunit           => clm3%g%l%c%p%landunit
    fsa                 => clm3%g%l%c%p%pef%fsa
    fsr                 => clm3%g%l%c%p%pef%fsr
    eflx_lwrad_out      => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net      => clm3%g%l%c%p%pef%eflx_lwrad_net
    sabv                => clm3%g%l%c%p%pef%sabv
    sabg                => clm3%g%l%c%p%pef%sabg
    eflx_sh_tot         => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lh_tot         => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_soil_grnd      => clm3%g%l%c%p%pef%eflx_soil_grnd
    errsol              => clm3%g%l%c%p%pebal%errsol
    errseb              => clm3%g%l%c%p%pebal%errseb
    errlon              => clm3%g%l%c%p%pebal%errlon
    netrad              => clm3%g%l%c%p%pef%netrad
    eflx_wasteheat_pft  => clm3%g%l%c%p%pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => clm3%g%l%c%p%pef%eflx_heat_from_ac_pft
    eflx_traffic_pft    => clm3%g%l%c%p%pef%eflx_traffic_pft

    ! Assign local pointers to derived type scalar members (gridcell-level)

    qflx_runoffg       => clm3%g%gwf%qflx_runoffg
    qflx_liq_dynbal    => clm3%g%gwf%qflx_liq_dynbal
    qflx_snwcp_iceg    => clm3%g%gwf%qflx_snwcp_iceg
    qflx_ice_dynbal    => clm3%g%gwf%qflx_ice_dynbal
    eflx_sh_totg       => clm3%g%gef%eflx_sh_totg
    eflx_dynbal        => clm3%g%gef%eflx_dynbal

    ! Determine column level incoming snow and rain
    ! Assume no incident precipitation on urban wall columns
    ! (as in Hydrology1Mod.F90).

    do c = lbc , ubc
      g = cgridcell(c)
      if ( ctype(c) == icol_sunwall .or. &
           ctype(c) == icol_shadewall ) then
        forc_rain_col(c) = 0.0D0
        forc_snow_col(c) = 0.0D0
      else
        forc_rain_col(c) = forc_rain(g)
        forc_snow_col(c) = forc_snow(g)
      end if
    end do

    ! Water balance check

    do c = lbc , ubc
      g = cgridcell(c)
      l = clandunit(c)

      ! add qflx_drain_perched and qflx_flood
      if ( cactive(c) ) then
        errh2o(c) = endwb(c) - begwb(c) - &
                (forc_rain_col(c) + forc_snow_col(c) +                    &
                 qflx_floodc(c) + qflx_irrig(c) -                         &
                 qflx_evap_tot(c) - qflx_surf(c)  - qflx_h2osfc_surf(c) - &
                 qflx_qrgwl(c) - qflx_drain(c) - qflx_drain_perched(c) -  &
                 qflx_snwcp_ice(c)) * dtsrf
      else
        errh2o(c) = 0.0D0
      end if
    end do

    found = .false.
    do c = lbc , ubc
      if ( abs(errh2o(c)) > 1D-7 ) then
        found = .true.
        indexc = c
      end if
    end do

    if ( found ) then
      write(stderr,*) 'WARNING:  ktau         : ', ktau
      write(stderr,*) 'WARNING:  indexc       : ', indexc
      write(stderr,*) 'WARNING:  ctype(indexc): ', ctype(indexc)
      write(stderr,*) 'WARNING:  water balance error ', errh2o(indexc)
      if ( (ctype(indexc) == icol_roof .or.        &
            ctype(indexc) == icol_road_imperv .or. &
            ctype(indexc) == icol_road_perv) .and. &
            abs(errh2o(indexc)) > 0.10D0 .and. (ktau > 2) ) then
        write(stderr,*) &
              'clm urban model is stopping - error is greater than .10'
        write(stderr,*) &
              'ktau = ',ktau,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
        write(stderr,*)'ctype(indexc): ',ctype(indexc)
        write(stderr,*)'forc_rain    = ',forc_rain_col(indexc)
        write(stderr,*)'forc_snow    = ',forc_snow_col(indexc)
        write(stderr,*)'endwb        = ',endwb(indexc)
        write(stderr,*)'begwb        = ',begwb(indexc)
        write(stderr,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
        write(stderr,*)'qflx_irrig   = ',qflx_irrig(indexc)
        write(stderr,*)'qflx_surf    = ',qflx_surf(indexc)
        write(stderr,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
        write(stderr,*)'qflx_drain   = ',qflx_drain(indexc)
        write(stderr,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
        call fatal(__FILE__,__LINE__,'clm model is stopping')
      else if ( abs(errh2o(indexc)) > 0.10D0 .and. (ktau > 2) ) then
        write(stderr,*)'clm model is stopping - error is greater than .10'
        write(stderr,*) &
             'ktau = ',ktau,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
        write(stderr,*)'ctype(indexc): ',ctype(indexc)
        write(stderr,*)'forc_rain    = ',forc_rain_col(indexc)
        write(stderr,*)'forc_snow    = ',forc_snow_col(indexc)
        write(stderr,*)'endwb        = ',endwb(indexc)
        write(stderr,*)'begwb        = ',begwb(indexc)
        write(stderr,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
        write(stderr,*)'qflx_irrig   = ',qflx_irrig(indexc)
        write(stderr,*)'qflx_surf    = ',qflx_surf(indexc)
        write(stderr,*)'qflx_h2osfc_surf    = ',qflx_h2osfc_surf(indexc)
        write(stderr,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
        write(stderr,*)'qflx_drain   = ',qflx_drain(indexc)
        write(stderr,*)'qflx_drain_perched   = ',qflx_drain_perched(indexc)
        write(stderr,*)'qflx_flood   = ',qflx_floodc(indexc)
        write(stderr,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
        call fatal(__FILE__,__LINE__,'clm model is stopping')
      end if
    end if

    ! Snow balance check

    do c = lbc , ubc
      g = cgridcell(c)
      l = clandunit(c)
      ! As defined here, snow_sources - snow_sinks will equal the change
      ! in h2osno at any given time step but only if there is at least one
      ! snow layer.  h2osno also includes snow that is part of the soil
      ! column (an initial snow layer is only created if h2osno > 10mm).
      if ( snl(c) < 0 ) then
        snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + &
                qflx_dew_grnd(c)
        snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) + &
                qflx_snow_melt(c) + qflx_snwcp_ice(c) +          &
                qflx_snwcp_liq(c) + qflx_sl_top_soil(c)
        if ( ltype(l) == istdlak ) then
          if ( do_capsnow(c) ) then
            snow_sources(c) = qflx_snow_grnd_col(c) + &
                 frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) )
            snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + &
                    qflx_evap_grnd(c) ) + (qflx_snwcp_ice(c) +    &
                    qflx_snwcp_liq(c) - qflx_prec_grnd(c)) +      &
                    qflx_snow_melt(c)  + qflx_sl_top_soil(c)
          else
            snow_sources(c) = qflx_snow_grnd_col(c) +          &
                    frac_sno_eff(c) * (qflx_rain_grnd_col(c) + &
                    qflx_dew_snow(c) + qflx_dew_grnd(c) )
            snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + &
                    qflx_evap_grnd(c) ) + qflx_snow_melt(c) +     &
                    qflx_sl_top_soil(c)
          end if
        end if

        if ( ltype(l) == istsoil .or. &
             ltype(l) == istcrop .or. &
             ltype(l) == istwet ) then
          if ( do_capsnow(c) ) then
            snow_sources(c) = frac_sno_eff(c) * (qflx_dew_snow(c) + &
                    qflx_dew_grnd(c) ) + qflx_h2osfc_to_ice(c) +    &
                    qflx_prec_grnd(c)
            snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) +  &
                    qflx_evap_grnd(c)) + qflx_snwcp_ice(c) +       &
                    qflx_snwcp_liq(c) + qflx_snow_melt(c) + qflx_sl_top_soil(c)
          else
            snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) + &
                    frac_sno_eff(c) * (qflx_rain_grnd_col(c) +                 &
                    qflx_dew_snow(c) + qflx_dew_grnd(c) ) +                    &
                    qflx_h2osfc_to_ice(c)
            snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + &
                    qflx_evap_grnd(c)) + qflx_snow_melt(c) + qflx_sl_top_soil(c)
          end if
        end if

        errh2osno(c) = (h2osno(c) - h2osno_old(c)) - &
                (snow_sources(c) - snow_sinks(c)) * dtsrf
      else
        snow_sources(c) = 0.D0
        snow_sinks(c) = 0.D0
        errh2osno(c) = 0.D0
      end if
    end do

    found = .false.
    do c = lbc , ubc
      if ( cactive(c) .and. abs(errh2osno(c)) > 1.0D-7 ) then
        found = .true.
        indexc = c
      end if
    end do
    if ( found ) then
      write(stderr,*) 'WARNING:  snow balance error '
      write(stderr,*) ' ktau      = ',ktau
      write(stderr,*) ' indexc    = ',indexc
      write(stderr,*) ' ctype     = ',ctype(indexc)
      write(stderr,*) ' ltype     = ',ltype(clandunit(indexc))
      write(stderr,*) ' errh2osno = ',errh2osno(indexc)
      if ( abs(errh2osno(indexc)) > 0.1D0 .and. (ktau > 2) ) then
        write(stderr,*)'clm model is stopping - error is greater than .10'
        write(stderr,*)'snl: ',snl(indexc)
        write(stderr,*)'h2osno: ',h2osno(indexc)
        write(stderr,*)'h2osno_old: ',h2osno_old(indexc)
        write(stderr,*)'snow_sources: ', snow_sources(indexc)
        write(stderr,*)'snow_sinks: ', snow_sinks(indexc)
        write(stderr,*)'qflx_prec_grnd: ',qflx_prec_grnd(indexc)*dtsrf
        write(stderr,*)'qflx_sub_snow: ',qflx_sub_snow(indexc)*dtsrf
        write(stderr,*)'qflx_evap_grnd: ',qflx_evap_grnd(indexc)*dtsrf
        write(stderr,*)'qflx_top_soil: ',qflx_top_soil(indexc)*dtsrf
        write(stderr,*)'qflx_dew_snow: ',qflx_dew_snow(indexc)*dtsrf
        write(stderr,*)'qflx_dew_grnd: ',qflx_dew_grnd(indexc)*dtsrf
        write(stderr,*)'qflx_snwcp_ice: ',qflx_snwcp_ice(indexc)*dtsrf
        write(stderr,*)'qflx_snwcp_liq: ',qflx_snwcp_liq(indexc)*dtsrf
        write(stderr,*)'qflx_sl_top_soil: ',qflx_sl_top_soil(indexc)*dtsrf
        call fatal(__FILE__,__LINE__,'clm model is stopping')
      end if
    end if

    ! Energy balance checks

    do p = lbp , ubp
      if ( pactive(p) ) then
        l = plandunit(p)
        g = pgridcell(p)

        ! Solar radiation energy balance
        ! Do not do this check for an urban pft since it will not
        ! balance on a per-column level because of interactions between
        ! columns and since a separate check is done in the urban radiation
        ! module
        if ( ltype(l) /= isturb ) then
          errsol(p) = fsa(p) + fsr(p) -                &
                  (forc_solad(g,1) + forc_solad(g,2) + &
                   forc_solai(g,1) + forc_solai(g,2))
        else
          errsol(p) = spval
        end if

        ! Longwave radiation energy balance
        ! Do not do this check for an urban pft since it will not
        ! balance on a per-column level because of interactions between
        ! columns and since a separate check is done in the urban radiation
        ! module
        if ( ltype(l) /= isturb ) then
          errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(g)
        else
          errlon(p) = spval
        end if

        ! Surface energy balance
        ! Changed to using (eflx_lwrad_net) here instead of
        ! (forc_lwrad - eflx_lwrad_out) because
        ! there are longwave interactions between urban columns
        ! (and therefore pfts).
        ! For surfaces other than urban, (eflx_lwrad_net) equals
        ! (forc_lwrad - eflx_lwrad_out),
        ! and a separate check is done above for these terms.

        if ( ltype(l) /= isturb ) then
          c = pcolumn(p)
          errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(g) - &
                  eflx_lwrad_out(p) - eflx_sh_tot(p) -        &
                  eflx_lh_tot(p) - eflx_soil_grnd(p)
        else
          errseb(p) = sabv(p) + sabg(p) - eflx_lwrad_net(p) -           &
                  eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) + &
                  eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) +    &
                  eflx_traffic_pft(p)
        end if
        netrad(p) = fsa(p) - eflx_lwrad_net(p)
      end if
    end do

    ! Solar radiation energy balance check

    found = .false.
    indexg = 0
    indexp = 0
    do p = lbp , ubp
      if (pactive(p)) then
        if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > .10D0) ) then
          found = .true.
          indexp = p
          indexg = pgridcell(p)
        end if
      end if
    end do
    if ( found  .and. (ktau > 2) ) then
      write(stderr,100) &
         'BalanceCheck: solar radiation balance error', ktau, indexp, &
         errsol(indexp)
      write(stderr,*)'fsa          = ',fsa(indexp)
      write(stderr,*)'fsr          = ',fsr(indexp)
      write(stderr,*)'forc_solad(1)= ',forc_solad(indexg,1)
      write(stderr,*)'forc_solad(2)= ',forc_solad(indexg,2)
      write(stderr,*)'forc_solai(1)= ',forc_solai(indexg,1)
      write(stderr,*)'forc_solai(2)= ',forc_solai(indexg,2)
      write(stderr,*)'forc_tot     = ',forc_solad(indexg,1) + &
                                       forc_solad(indexg,2) + &
                                       forc_solai(indexg,1) + &
                                       forc_solai(indexg,2)
      call fatal(__FILE__,__LINE__,'clm model is stopping')
    end if

    ! Longwave radiation energy balance check

    found = .false.
    do p = lbp , ubp
      if ( pactive(p) ) then
        if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > .10D0) ) then
          found = .true.
          indexp = p
        end if
      end if
    end do
    if ( found  .and. (ktau > 2) ) then
      write(stderr,100) &
          'BalanceCheck: longwave enery balance error',ktau,indexp, &
          errlon(indexp)
      call fatal(__FILE__,__LINE__,'clm model is stopping')
    end if

    ! Surface energy balance check

    found = .false.
    do p = lbp , ubp
      if ( pactive(p) ) then
        if ( abs(errseb(p)) > 0.10D0 ) then
          found = .true.
          indexp = p
        end if
      end if
    end do
    if ( found  .and. (ktau > 2) ) then
      write(stderr,100) &
         'BalanceCheck: surface flux energy balance error',ktau,indexp, &
         errseb(indexp)
      write(stderr,*)' sabv           = ',sabv(indexp)
      c = pcolumn(indexp)
      write(stderr,*)' column      = ',c
      write(stderr,*)' sabg           = ',sabg(indexp), &
               ((1.D0- frac_sno(c))*sabg_soil(indexp) + &
                  frac_sno(c)*sabg_snow(indexp)),sabg_chk(indexp)
      write(stderr,*)' eflx_lwrad_net = ',eflx_lwrad_net(indexp)
      write(stderr,*)' eflx_sh_tot    = ',eflx_sh_tot(indexp)
      write(stderr,*)' eflx_lh_tot    = ',eflx_lh_tot(indexp)
      write(stderr,*)' eflx_soil_grnd = ',eflx_soil_grnd(indexp)
      call fatal(__FILE__,__LINE__,'clm model is stopping')
    end if

    ! Soil energy balance check

    found = .false.
    do c = lbc , ubc
      if ( abs(errsoi_col(c)) > 1.0D-7 ) then
        found = .true.
        indexc = c
      end if
    end do
    if ( found ) then
      if (abs(errsoi_col(indexc)) > .10D0 .and. (ktau > 2) ) then
        write(stderr,100) &
          'BalanceCheck: soil balance error',ktau,indexc,errsoi_col(indexc)
        write(stderr,*)'ktau = ',ktau,' indexc= ',indexc, &
                ' errsoi_col= ',errsoi_col(indexc)
        call fatal(__FILE__,__LINE__,'clm model is stopping')
      end if
    end if

    ! Update SH and RUNOFF for dynamic land cover change energy and water fluxes
    call c2g( lbc,ubc,lbl,ubl,lbg,ubg,                    &
              qflx_runoff(lbc:ubc),qflx_runoffg(lbg:ubg), &
              c2l_scale_type= 'urbanf',l2g_scale_type='unity' )
    do g = lbg , ubg
      qflx_runoffg(g) = qflx_runoffg(g) - qflx_liq_dynbal(g)
    enddo

    call c2g( lbc,ubc,lbl,ubl,lbg,ubg,                          &
              qflx_snwcp_ice(lbc:ubc),qflx_snwcp_iceg(lbg:ubg), &
              c2l_scale_type= 'urbanf',l2g_scale_type='unity' )
    do g = lbg , ubg
      qflx_snwcp_iceg(g) = qflx_snwcp_iceg(g) - qflx_ice_dynbal(g)
    enddo

    call p2g( lbp,ubp,lbc,ubc,lbl,ubl,lbg,ubg,            &
              eflx_sh_tot(lbp:ubp),eflx_sh_totg(lbg:ubg), &
              p2c_scale_type='unity',                     &
              c2l_scale_type='urbanf',                    &
              l2g_scale_type='unity' )
    do g = lbg , ubg
      eflx_sh_totg(g) =  eflx_sh_totg(g) - eflx_dynbal(g)
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb

    call c2g( lbc,ubc,lbl,ubl,lbg,ubg,                &
              endwb(lbc:ubc),tws(lbg:ubg),            &
              c2l_scale_type= 'urbanf',l2g_scale_type='unity' )

    ! second add river storage as gridcell average depth
    ! 1.e-3 converts [m3/km2] to [mm]
    do g = lbg , ubg
      tws(g) = tws(g) + volr(g) / area(g) * 1.D-3
    enddo

100 format (1x,a,' ktau =',i10,' point =',i6,' imbalance =',f12.6,' W/m2')
! 200 format (1x,a,' ktau =',i10,' point =',i6,' imbalance =',f12.6,' mm')

  end subroutine BalanceCheck

end module mod_clm_balancecheck
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
