module mod_clm_hydrology2
  !
  ! Calculation of soil/snow hydrology.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams

  implicit none

  private

  save

  public :: Hydrology2        ! Calculates soil/snow hydrology

  contains
  !
  ! This is the main subroutine to execute the calculation of soil/snow
  ! hydrology
  ! Calling sequence is:
  !  Hydrology2:             surface hydrology driver
  !    -> SnowWater:         change of snow mass and snow water onto soil
  !    -> SurfaceRunoff:     surface runoff
  !    -> Infiltration:      infiltration into surface soil layer
  !    -> SoilWater:         soil water movement between layers
  !          -> Tridiagonal  tridiagonal matrix solution
  !    -> Drainage:          subsurface runoff
  !    -> SnowCompaction:    compaction of snow layers
  !    -> CombineSnowLayers: combine snow layers that are thinner than minimum
  !    -> DivideSnowLayers:  subdivide snow layers that are thicker than maximum
  !
  subroutine Hydrology2(lbc, ubc, &
                        num_nolakec, filter_nolakec, &
                        num_hydrologyc, filter_hydrologyc, &
                        num_urbanc, filter_urbanc, &
                        num_snowc, filter_snowc, &
                        num_nosnowc, filter_nosnowc)
    use mod_clm_type
    use mod_clm_atmlnd, only : clm_a2l
    use mod_clm_varcon, only : denh2o, denice, istice, istwet, istsoil, &
            isturb, spval, icol_roof, icol_road_imperv, icol_road_perv,  &
            icol_sunwall, icol_shadewall, istdlak, tfrz, hfus, grav
    use mod_clm_varcon, only : istcrop
    use mod_clm_varpar, only : nlevgrnd, nlevsno, nlevsoi, nlevurb
    use mod_clm_snowhydrology, only : SnowCompaction, CombineSnowLayers, &
            DivideSnowLayers, SnowWater, BuildSnowFilter
    use mod_clm_soilhydrology, only : Infiltration, SoilWater, Drainage, &
            SurfaceRunoff
#if (defined VICHYDRO)
    use mod_clm_vicmap, only : CLMVICMap
#endif
    implicit none
    integer(ik4), intent(in) :: lbc, ubc    ! column bounds
    ! number of column non-lake points in column filter
    integer(ik4), intent(in) :: num_nolakec
    ! column filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakec(ubc-lbc+1)
    ! number of column soil points in column filter
    integer(ik4), intent(in) :: num_hydrologyc
    ! column filter for soil points
    integer(ik4), intent(in) :: filter_hydrologyc(ubc-lbc+1)
    ! number of column urban points in column filter
    integer(ik4), intent(in) :: num_urbanc
    ! column filter for urban points
    integer(ik4), intent(in) :: filter_urbanc(ubc-lbc+1)
    ! number of column snow points
    integer(ik4)  :: num_snowc
    ! column filter for snow points
    integer(ik4)  :: filter_snowc(ubc-lbc+1)
    ! number of column non-snow points
    integer(ik4)  :: num_nosnowc
    ! column filter for non-snow points
    integer(ik4)  :: filter_nosnowc(ubc-lbc+1)

    real(rk8), pointer, contiguous :: snow_depth(:)   ! snow height of snow covered area (m)
    real(rk8), pointer, contiguous :: snowdp(:)       ! gridcell averaged snow height (m)
    !eff.  snow cover fraction (col) [frc]
    real(rk8), pointer, contiguous :: frac_sno_eff(:)
    real(rk8), pointer, contiguous :: qflx_evap_soi(:) ! soil evaporation
    real(rk8), pointer, contiguous :: h2osfc(:)        ! surface water (mm)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer, contiguous :: frac_h2osfc(:)
    real(rk8), pointer, contiguous :: t_h2osfc(:)      ! surface water temperature
    ! sub-surface runoff from perched zwt (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_drain_perched(:)
    ! gridcell flux of flood water from RTM
    real(rk8), pointer, contiguous :: qflx_floodg(:)
    ! surface water runoff (mm/s)
    real(rk8), pointer, contiguous :: qflx_h2osfc_surf(:)
    ! true=>do computations on this column (see reweightMod for details)
    logical, pointer, contiguous :: cactive(:)
    integer(ik4), pointer, contiguous :: cgridcell(:)  ! column's gridcell
    integer(ik4), pointer, contiguous :: clandunit(:)  ! column's landunit
    integer(ik4), pointer, contiguous :: ityplun(:)    ! landunit type
    integer(ik4), pointer, contiguous :: ctype(:)      ! column type
    integer(ik4), pointer, contiguous :: snl(:)        ! number of snow layers
    real(rk8), pointer, contiguous :: h2ocan(:)         ! canopy water (mm H2O)
    real(rk8), pointer, contiguous :: h2osno(:)         ! snow water (mm H2O)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer, contiguous :: watsat(:,:)
    real(rk8), pointer, contiguous :: sucsat(:,:)    ! minimum soil suction (mm)
    real(rk8), pointer, contiguous :: bsw(:,:)       ! Clapp and Hornberger "b"
    real(rk8), pointer, contiguous :: z(:,:)         ! layer depth  (m)
    real(rk8), pointer, contiguous :: forc_rain(:)   ! rain rate [mm/s]
    real(rk8), pointer, contiguous :: forc_snow(:)   ! snow rate [mm/s]
    real(rk8), pointer, contiguous :: begwb(:)       ! water mass begining of the time step
    ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8), pointer, contiguous :: qflx_evap_tot(:)
    ! restriction for min of soil potential (mm)
    real(rk8), pointer, contiguous :: smpmin(:)

    real(rk8), pointer, contiguous :: dz(:,:)  ! layer thickness depth (m)
    real(rk8), pointer, contiguous :: zi(:,:)  ! interface depth (m)
    real(rk8), pointer, contiguous :: zwt(:)   ! water table depth (m)
    real(rk8), pointer, contiguous :: fcov(:)  ! fractional impermeable area
    real(rk8), pointer, contiguous :: fsat(:)  ! fractional area with water table at surface
    real(rk8), pointer, contiguous :: wa(:)    ! water in the unconfined aquifer (mm)
    real(rk8), pointer, contiguous :: qcharge(:)   ! aquifer recharge rate (mm/s)
    real(rk8), pointer, contiguous :: smp_l(:,:)   ! soil matrix potential [mm]
    real(rk8), pointer, contiguous :: hk_l(:,:)    ! hydraulic conductivity (mm/s)
    real(rk8), pointer, contiguous :: qflx_rsub_sat(:) ! soil saturation excess [mm h2o/s]

    real(rk8), pointer, contiguous :: endwb(:)  ! water mass end of the time step
    ! soil water as frac. of whc for top 0.05 m  !F. Li and S. Levis
    real(rk8), pointer, contiguous :: wf(:)
    ! soil water as frac. of whc for top 0.17 m added by F. Li and S. Levis
    real(rk8), pointer, contiguous :: wf2(:)
    real(rk8), pointer, contiguous :: snowice(:)       ! average snow ice lens
    real(rk8), pointer, contiguous :: snowliq(:)       ! average snow liquid water
    real(rk8), pointer, contiguous :: t_grnd(:)        ! ground temperature (Kelvin)
    real(rk8), pointer, contiguous :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(rk8), pointer, contiguous :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(rk8), pointer, contiguous :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    ! soil temperature in top 10cm of soil (Kelvin)
    real(rk8), pointer, contiguous :: t_soi_10cm(:)
    ! soil temperature in top 17cm of soil (Kelvin) added by F. Li and S. Levis
    real(rk8), pointer, contiguous :: tsoi17(:)
    ! liquid water + ice lens in top 10cm of soil (kg/m2)
    real(rk8), pointer, contiguous :: h2osoi_liqice_10cm(:)
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8), pointer, contiguous :: h2osoi_vol(:,:)
    real(rk8), pointer, contiguous :: qflx_drain(:) ! sub-surface runoff (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_surf(:)  ! surface runoff (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_infl(:)  ! infiltration (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_qrgwl(:) ! qflx_surf at glaciers, wetlands, lakes
    real(rk8), pointer, contiguous :: qflx_irrig(:) ! irrigation flux (mm H2O /s)
    ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_runoff(:)
    ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_runoff_u(:)
    ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8), pointer, contiguous :: qflx_runoff_r(:)
    real(rk8), pointer, contiguous :: t_grnd_u(:)   ! Urban ground temperature (Kelvin)
    real(rk8), pointer, contiguous :: t_grnd_r(:)   ! Rural ground temperature (Kelvin)
    ! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(rk8), pointer, contiguous :: qflx_snwcp_ice(:)
    ! soil water potential in each soil layer (MPa)
    real(rk8), pointer, contiguous :: soilpsi(:,:)

    real(rk8), pointer, contiguous :: snot_top(:)  ! snow temperature in top layer (col) [K]
    ! temperature gradient in top layer (col) [K m-1]
    real(rk8), pointer, contiguous :: dTdz_top(:)
    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(rk8), pointer, contiguous :: snw_rds(:,:)
    ! effective snow grain size, top layer(col) [microns]
    real(rk8), pointer, contiguous :: snw_rds_top(:)
    ! liquid water fraction in top snow layer (col) [frc]
    real(rk8), pointer, contiguous :: sno_liq_top(:)
    real(rk8), pointer, contiguous :: frac_sno(:)    ! snow cover fraction (col) [frc]
    real(rk8), pointer, contiguous :: h2osno_top(:)  ! mass of snow in top layer (col) [kg]

    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_bcpho(:,:)
    ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_bcphi(:,:)
    ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_bctot(:,:)
    ! total mass of BC in snow column (col) [kg]
    real(rk8), pointer, contiguous :: mss_bc_col(:)
    ! total mass of BC in top snow layer (col) [kg]
    real(rk8), pointer, contiguous :: mss_bc_top(:)
    ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_bcphi(:,:)
    ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_bcpho(:,:)
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_ocpho(:,:)
    ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_ocphi(:,:)
    ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_octot(:,:)
    ! total mass of OC in snow column (col) [kg]
    real(rk8), pointer, contiguous :: mss_oc_col(:)
    ! total mass of OC in top snow layer (col) [kg]
    real(rk8), pointer, contiguous :: mss_oc_top(:)
    ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_ocphi(:,:)
    ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_ocpho(:,:)

    ! mass of dust species 1 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_dst1(:,:)
    ! mass of dust species 2 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_dst2(:,:)
    ! mass of dust species 3 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_dst3(:,:)
    ! mass of dust species 4 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_dst4(:,:)
    ! total mass of dust in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous :: mss_dsttot(:,:)
    ! total mass of dust in snow column (col) [kg]
    real(rk8), pointer, contiguous :: mss_dst_col(:)
    ! total mass of dust in top snow layer (col) [kg]
    real(rk8), pointer, contiguous :: mss_dst_top(:)
    ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_dst1(:,:)
    ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_dst2(:,:)
    ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_dst3(:,:)
    ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous :: mss_cnc_dst4(:,:)
    ! true => do snow capping
    logical, pointer, contiguous :: do_capsnow(:)

    integer(ik4)  :: g,l,c,j,fc    ! indices
    ! partial volume of liquid water in layer
    real(rk8) :: vol_liq(lbc:ubc,1:nlevgrnd)
    ! change in soil water
    real(rk8) :: dwat(lbc:ubc,1:nlevgrnd)
    ! hydraulic conductivity (mm h2o/s)
    real(rk8) :: hk(lbc:ubc,1:nlevgrnd)
    real(rk8) :: dhkdw(lbc:ubc,1:nlevgrnd)  ! d(hk)/d(vol_liq)
    ! temporary variables for soilpsi calculation
#if (defined CN)
    real(rk8) :: psi,vwc,fsattmp
    real(rk8) :: watdry       ! temporary
    ! soil water wgted by depth to maximum depth of 0.5 m
    real(rk8) :: rwat(lbc:ubc)
    ! same as rwat but at saturation
    real(rk8) :: swat(lbc:ubc)
    ! thickness of soil layers contributing to rwat (m)
    real(rk8) :: rz(lbc:ubc)
    real(rk8) :: tsw       ! volumetric soil water to 0.5 m
    real(rk8) :: stsw      ! volumetric soil water to 0.5 m at saturation
#endif
    real(rk8) :: snowmass  ! liquid+ice snow mass in a layer [kg/m2]
    ! temporary factor used to correct for snow capping
    real(rk8) :: snowcap_scl_fct
    ! fraction of soil layer contributing to 10cm total soil water
    real(rk8) :: fracl
    real(rk8) :: s_node ! soil wetness (-)
    real(rk8) :: icefrac(lbc:ubc,1:nlevsoi)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_rain => clm_a2l%forc_rain
    forc_snow => clm_a2l%forc_snow

    ! Assign local pointers to derived subtypes components (landunit-level)

    ityplun => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    snow_depth        => clm3%g%l%c%cps%snow_depth
    snowdp            => clm3%g%l%c%cps%snowdp
    frac_sno_eff      => clm3%g%l%c%cps%frac_sno_eff
    qflx_evap_soi     => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    h2osfc            => clm3%g%l%c%cws%h2osfc
    frac_h2osfc       => clm3%g%l%c%cps%frac_h2osfc
    t_h2osfc          => clm3%g%l%c%ces%t_h2osfc
    qflx_drain_perched=> clm3%g%l%c%cwf%qflx_drain_perched
    qflx_floodg       => clm_a2l%forc_flood
    qflx_h2osfc_surf  => clm3%g%l%c%cwf%qflx_h2osfc_surf
    cactive           => clm3%g%l%c%active
    cgridcell         => clm3%g%l%c%gridcell
    clandunit         => clm3%g%l%c%landunit
    ctype             => clm3%g%l%c%itype
    snl               => clm3%g%l%c%cps%snl
    t_grnd            => clm3%g%l%c%ces%t_grnd
    h2ocan            => clm3%g%l%c%cws%pws_a%h2ocan
    h2osno            => clm3%g%l%c%cws%h2osno
    wf                => clm3%g%l%c%cps%wf
    wf2               => clm3%g%l%c%cps%wf2
    snowice           => clm3%g%l%c%cws%snowice
    snowliq           => clm3%g%l%c%cws%snowliq
    zwt               => clm3%g%l%c%cws%zwt
    fcov              => clm3%g%l%c%cws%fcov
    fsat              => clm3%g%l%c%cws%fsat
    wa                => clm3%g%l%c%cws%wa
    qcharge           => clm3%g%l%c%cws%qcharge
    watsat            => clm3%g%l%c%cps%watsat
    sucsat            => clm3%g%l%c%cps%sucsat
    bsw               => clm3%g%l%c%cps%bsw
    z                 => clm3%g%l%c%cps%z
    dz                => clm3%g%l%c%cps%dz
    zi                => clm3%g%l%c%cps%zi
    t_soisno          => clm3%g%l%c%ces%t_soisno
    h2osoi_ice        => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq        => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol        => clm3%g%l%c%cws%h2osoi_vol
    t_soi_10cm         => clm3%g%l%c%ces%t_soi_10cm
    tsoi17             => clm3%g%l%c%ces%tsoi17
    h2osoi_liqice_10cm => clm3%g%l%c%cws%h2osoi_liqice_10cm
    qflx_evap_tot     => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    qflx_drain        => clm3%g%l%c%cwf%qflx_drain
    qflx_surf         => clm3%g%l%c%cwf%qflx_surf
    qflx_infl         => clm3%g%l%c%cwf%qflx_infl
    qflx_qrgwl        => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_irrig        => clm3%g%l%c%cwf%qflx_irrig
    endwb             => clm3%g%l%c%cwbal%endwb
    begwb             => clm3%g%l%c%cwbal%begwb
    soilpsi           => clm3%g%l%c%cps%soilpsi
    smp_l             => clm3%g%l%c%cws%smp_l
    hk_l              => clm3%g%l%c%cws%hk_l
    qflx_rsub_sat     => clm3%g%l%c%cwf%qflx_rsub_sat
    qflx_runoff       => clm3%g%l%c%cwf%qflx_runoff
    qflx_runoff_u     => clm3%g%l%c%cwf%qflx_runoff_u
    qflx_runoff_r     => clm3%g%l%c%cwf%qflx_runoff_r
    t_grnd_u          => clm3%g%l%c%ces%t_grnd_u
    t_grnd_r          => clm3%g%l%c%ces%t_grnd_r
    snot_top          => clm3%g%l%c%cps%snot_top
    dTdz_top          => clm3%g%l%c%cps%dTdz_top
    snw_rds           => clm3%g%l%c%cps%snw_rds
    snw_rds_top       => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top       => clm3%g%l%c%cps%sno_liq_top
    frac_sno          => clm3%g%l%c%cps%frac_sno
    h2osno_top        => clm3%g%l%c%cps%h2osno_top
    mss_bcpho         => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi         => clm3%g%l%c%cps%mss_bcphi
    mss_bctot         => clm3%g%l%c%cps%mss_bctot
    mss_bc_col        => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top        => clm3%g%l%c%cps%mss_bc_top
    mss_cnc_bcphi     => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho     => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_ocpho         => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi         => clm3%g%l%c%cps%mss_ocphi
    mss_octot         => clm3%g%l%c%cps%mss_octot
    mss_oc_col        => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top        => clm3%g%l%c%cps%mss_oc_top
    mss_cnc_ocphi     => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho     => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_dst1          => clm3%g%l%c%cps%mss_dst1
    mss_dst2          => clm3%g%l%c%cps%mss_dst2
    mss_dst3          => clm3%g%l%c%cps%mss_dst3
    mss_dst4          => clm3%g%l%c%cps%mss_dst4
    mss_dsttot        => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col       => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top       => clm3%g%l%c%cps%mss_dst_top
    mss_cnc_dst1      => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2      => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3      => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4      => clm3%g%l%c%cps%mss_cnc_dst4
    do_capsnow        => clm3%g%l%c%cps%do_capsnow
    qflx_snwcp_ice    => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    smpmin            => clm3%g%l%c%cps%smpmin

    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(lbc,ubc,num_snowc,filter_snowc,num_nosnowc,filter_nosnowc)

    ! Determine soil hydrology
#if (defined VICHYDRO)
    ! mapping soilmoist from CLM to VIC layers for runoff calculations
    call CLMVICMap(lbc, ubc, num_hydrologyc, filter_hydrologyc)
#endif

    ! moved vol_liq from SurfaceRunoff to Infiltration
    call SurfaceRunoff(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                       num_urbanc, filter_urbanc, icefrac )

    call Infiltration(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                      num_urbanc, filter_urbanc, vol_liq)

    call SoilWater(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                   dwat, hk, dhkdw)

#if (defined VICHYDRO)
    ! mapping soilmoist from CLM to VIC layers for runoff calculations
    call CLMVICMap(lbc, ubc, num_hydrologyc, filter_hydrologyc)
#endif

    call Drainage(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                  num_urbanc, filter_urbanc, icefrac)

    ! Natural compaction and metamorphosis.

    call SnowCompaction(lbc, ubc, num_snowc, filter_snowc)

    ! Combine thin snow elements

    call CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)

    ! Divide thick snow elements

    call DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)

    ! Set empty snow layers to zero

    do j = -nlevsno+1,0
      do fc = 1, num_snowc
        c = filter_snowc(fc)
        if (j <= snl(c) .and. snl(c) > -nlevsno) then
          h2osoi_ice(c,j) = 0._rk8
          h2osoi_liq(c,j) = 0._rk8
          t_soisno(c,j) = 0._rk8
          dz(c,j) = 0._rk8
          z(c,j) = 0._rk8
          zi(c,j-1) = 0._rk8
        end if
      end do
    end do

    ! Build new snow filter

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      snowice(c) = 0._rk8
      snowliq(c) = 0._rk8
    end do

    do j = -nlevsno+1, 0
      do fc = 1, num_snowc
        c = filter_snowc(fc)
        if (j >= snl(c)+1) then
          snowice(c) = snowice(c) + h2osoi_ice(c,j)
          snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
        end if
      end do
    end do

    ! Calculate column average snow depth
    do c = lbc,ubc
      snowdp(c) = snow_depth(c) * frac_sno_eff(c)
    end do

    ! Determine ground temperature, ending water balance and volumetric
    ! soil water
    ! Calculate soil temperature and total water (liq+ice) in top 10cm of soil
    ! Calculate soil temperature and total water (liq+ice) in top 17cm of soil
    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)
      if (ityplun(l) /= isturb) then
        t_soi_10cm(c) = 0._rk8
        tsoi17(c) = 0._rk8
        h2osoi_liqice_10cm(c) = 0._rk8
      end if
    end do
    do j = 1, nlevsoi
      do fc = 1, num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if (ityplun(l) /= isturb) then
          ! soil T at top 17 cm added by F. Li and S. Levis
          if (zi(c,j) <= 0.17_rk8) then
            fracl = 1._rk8
            tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
          else
            if (zi(c,j) > 0.17_rk8 .and. zi(c,j-1) < 0.17_rk8) then
              fracl = (0.17_rk8 - zi(c,j-1))/dz(c,j)
              tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
            end if
          end if

          if (zi(c,j) <= 0.1_rk8) then
            fracl = 1._rk8
            t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
            h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                  (h2osoi_liq(c,j)+h2osoi_ice(c,j))*fracl
          else
            if (zi(c,j) > 0.1_rk8 .and. zi(c,j-1) < 0.1_rk8) then
              fracl = (0.1_rk8 - zi(c,j-1))/dz(c,j)
              t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
              h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                 (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                                          fracl
            end if
          end if
        end if
      end do
    end do

    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)

      ! t_grnd is weighted average of exposed soil and snow
      if (snl(c) < 0) then
        t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
               + (1.0_rk8 - frac_sno_eff(c)- frac_h2osfc(c)) * t_soisno(c,1) &
               + frac_h2osfc(c) * t_h2osfc(c)
      else
        t_grnd(c) = (1.0_rk8 - frac_h2osfc(c)) * t_soisno(c,1) + &
                frac_h2osfc(c) * t_h2osfc(c)
      end if

      if (ityplun(l)==isturb) then
        t_grnd_u(c) = t_soisno(c,snl(c)+1)
      else
        t_soi_10cm(c) = t_soi_10cm(c)/0.1_rk8
        tsoi17(c) =  tsoi17(c)/0.17_rk8         ! F. Li and S. Levis
      end if
      if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
        t_grnd_r(c) = t_soisno(c,snl(c)+1)
      end if
      if (ctype(c) == icol_roof .or.      &
          ctype(c) == icol_sunwall .or.   &
          ctype(c) == icol_shadewall .or. &
          ctype(c) == icol_road_imperv) then
        endwb(c) = h2ocan(c) + h2osno(c)
      else
        ! add h2osfc to water balance
        endwb(c) = h2ocan(c) + h2osno(c) + h2osfc(c) + wa(c)
      end if
    end do

    do j = 1, nlevgrnd
      do fc = 1, num_nolakec
        c = filter_nolakec(fc)
        if ( ( ctype(c) == icol_sunwall .or.   &
               ctype(c) == icol_shadewall .or. &
               ctype(c) == icol_roof) .and. j > nlevurb) then
          continue
        else
          endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + &
                            h2osoi_ice(c,j)/(dz(c,j)*denice)
        end if
      end do
    end do

    ! Determine wetland and land ice hydrology (must be placed here
    ! since need snow updated from CombineSnowLayers)

    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)
      g = cgridcell(c)
      if ( ityplun(l) == istwet .or. ityplun(l) == istice ) then
        qflx_drain(c)         = 0._rk8
        qflx_drain_perched(c) = 0._rk8
        qflx_h2osfc_surf(c)   = 0._rk8
        qflx_irrig(c)         = 0._rk8
        qflx_surf(c)          = 0._rk8
        qflx_infl(c)          = 0._rk8
        ! add flood water flux to runoff for wetlands/glaciers
        qflx_qrgwl(c) = forc_rain(g) + forc_snow(g) + &
                qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
                (endwb(c)-begwb(c))/dtsrf
        fcov(c)       = spval
        fsat(c)       = spval
        qcharge(c)    = spval
        qflx_rsub_sat(c) = spval
      else if (ityplun(l) == isturb .and. ctype(c) /= icol_road_perv) then
        fcov(c)               = spval
        fsat(c)               = spval
        qflx_drain_perched(c) = 0._rk8
        qflx_h2osfc_surf(c)   = 0._rk8
        qcharge(c)            = spval
        qflx_rsub_sat(c)      = spval
      end if

      qflx_runoff(c) = qflx_drain(c) + qflx_surf(c) + &
              qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

      if ((ityplun(l)==istsoil .or. ityplun(l)==istcrop) &
           .and. cactive(c)) then
        qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)
      end if
      if (ityplun(l)==isturb) then
        qflx_runoff_u(c) = qflx_runoff(c)
      else if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
        qflx_runoff_r(c) = qflx_runoff(c)
      end if
    end do

#if (defined CN)
    ! Update soilpsi.
    ! ZMS: Note this could be merged with the following loop updating
    ! smp_l in the future.
    do j = 1, nlevgrnd
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)

        if (h2osoi_liq(c,j) > 0._rk8) then

          vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)

          ! the following limit set to catch very small values of
          ! fractional saturation that can crash the calculation of psi

          ! use the same contants used in the supercool so that psi for
          ! frozen soils is consistent
          fsattmp = max(vwc/watsat(c,j), 0.001_rk8)
          psi = sucsat(c,j) * (-9.8e-6_rk8) * (fsattmp)**(-bsw(c,j))  ! Mpa
          soilpsi(c,j) = min(max(psi,-15.0_rk8),0._rk8)

        else
          soilpsi(c,j) = -15.0_rk8
        end if
      end do
    end do
#endif

    ! Update smp_l for history and for ch4Mod.
    ! ZMS: Note, this form, which seems to be the same as used in
    ! SoilWater, DOES NOT distinguish between ice and water volume,
    ! in contrast to the soilpsi calculation above.
    ! It won't be used in ch4Mod if t_soisno <= tfrz, though.
    do j = 1, nlevgrnd
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)

        s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_rk8)
        s_node = min(1.0_rk8, s_node)

        smp_l(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
        smp_l(c,j) = max(smpmin(c), smp_l(c,j))
      end do
    end do

#if (defined CN)
    ! Available soil water up to a depth of 0.05 m.
    ! Potentially available soil water (=whc) up to a depth of 0.05 m.
    ! Water content as fraction of whc up to a depth of 0.05 m.

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      rwat(c) = 0._rk8
      swat(c) = 0._rk8
      rz(c)   = 0._rk8
    end do

    do j = 1, nlevgrnd
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        !if (z(c,j)+0.5_rk8*dz(c,j) <= 0.5_rk8) then
        if (z(c,j)+0.5_rk8*dz(c,j) <= 0.05_rk8) then
          watdry = watsat(c,j) * (316230._rk8/sucsat(c,j)) ** (-1._rk8/bsw(c,j))
          rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
          swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
          rz(c) = rz(c) + dz(c,j)
        end if
      end do
    end do

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      if (rz(c) /= 0._rk8) then
        tsw  = rwat(c)/rz(c)
        stsw = swat(c)/rz(c)
      else
        watdry = watsat(c,1) * (316230._rk8/sucsat(c,1)) ** (-1._rk8/bsw(c,1))
        tsw = h2osoi_vol(c,1) - watdry
        stsw = watsat(c,1) - watdry
      end if
      wf(c) = tsw/stsw
    end do

    do j = 1, nlevgrnd
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        if (z(c,j)+0.5_rk8*dz(c,j) <= 0.17_rk8) then
          watdry = watsat(c,j) * (316230._rk8/sucsat(c,j)) ** (-1._rk8/bsw(c,j))
          rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
          swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
          rz(c) = rz(c) + dz(c,j)
        end if
      end do
    end do

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      if (rz(c) /= 0._rk8) then
        tsw  = rwat(c)/rz(c)
        stsw = swat(c)/rz(c)
      else
        watdry = watsat(c,1) * (316230._rk8/sucsat(c,1)) ** (-1._rk8/bsw(c,1))
        tsw = h2osoi_vol(c,1) - watdry
        stsw = watsat(c,1) - watdry
      end if
      wf2(c) = tsw/stsw
    end do
#endif

    !  Calculate column-integrated aerosol masses, and
    !  mass concentrations for radiative calculations and output
    !  (based on new snow level state, after SnowFilter is rebuilt.
    !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
    !  can be zero snow layers but an active column in filter)

    do fc = 1, num_snowc
      c = filter_snowc(fc)

      ! Zero column-integrated aerosol mass before summation
      mss_bc_col(c)  = 0._rk8
      mss_oc_col(c)  = 0._rk8
      mss_dst_col(c) = 0._rk8

      do j = -nlevsno+1, 0

        ! layer mass of snow:
        snowmass = h2osoi_ice(c,j)+h2osoi_liq(c,j)

        ! Correct the top layer aerosol mass to account for snow capping.
        ! This approach conserves the aerosol mass concentration
        ! (but not the aerosol amss) when snow-capping is invoked

        if (j == snl(c)+1) then
          if (do_capsnow(c)) then
            snowcap_scl_fct = snowmass / (snowmass+(qflx_snwcp_ice(c)*dtsrf))

            mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
            mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
            mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
            mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct

            mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
            mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
            mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
            mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct
          end if
        end if

        if (j >= snl(c)+1) then
          mss_bctot(c,j)     = mss_bcpho(c,j) + mss_bcphi(c,j)
          mss_bc_col(c)      = mss_bc_col(c)  + mss_bctot(c,j)
          mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
          mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

          mss_octot(c,j)     = mss_ocpho(c,j) + mss_ocphi(c,j)
          mss_oc_col(c)      = mss_oc_col(c)  + mss_octot(c,j)
          mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
          mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass

          mss_dsttot(c,j)    = mss_dst1(c,j)  + &
                  mss_dst2(c,j) + mss_dst3(c,j) + mss_dst4(c,j)
          mss_dst_col(c)     = mss_dst_col(c) + mss_dsttot(c,j)
          mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
          mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
          mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
          mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass

        else
          !set variables of empty snow layers to zero
          snw_rds(c,j)       = 0._rk8

          mss_bcpho(c,j)     = 0._rk8
          mss_bcphi(c,j)     = 0._rk8
          mss_bctot(c,j)     = 0._rk8
          mss_cnc_bcphi(c,j) = 0._rk8
          mss_cnc_bcpho(c,j) = 0._rk8

          mss_ocpho(c,j)     = 0._rk8
          mss_ocphi(c,j)     = 0._rk8
          mss_octot(c,j)     = 0._rk8
          mss_cnc_ocphi(c,j) = 0._rk8
          mss_cnc_ocpho(c,j) = 0._rk8

          mss_dst1(c,j)      = 0._rk8
          mss_dst2(c,j)      = 0._rk8
          mss_dst3(c,j)      = 0._rk8
          mss_dst4(c,j)      = 0._rk8
          mss_dsttot(c,j)    = 0._rk8
          mss_cnc_dst1(c,j)  = 0._rk8
          mss_cnc_dst2(c,j)  = 0._rk8
          mss_cnc_dst3(c,j)  = 0._rk8
          mss_cnc_dst4(c,j)  = 0._rk8
        end if
      end do

      ! top-layer diagnostics
      h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
      mss_bc_top(c)  = mss_bctot(c,snl(c)+1)
      mss_oc_top(c)  = mss_octot(c,snl(c)+1)
      mss_dst_top(c) = mss_dsttot(c,snl(c)+1)
    end do

    ! Zero mass variables in columns without snow
    do fc = 1, num_nosnowc
      c = filter_nosnowc(fc)

      h2osno_top(c)      = 0._rk8
      snw_rds(c,:)       = 0._rk8

      mss_bc_top(c)      = 0._rk8
      mss_bc_col(c)      = 0._rk8
      mss_bcpho(c,:)     = 0._rk8
      mss_bcphi(c,:)     = 0._rk8
      mss_bctot(c,:)     = 0._rk8
      mss_cnc_bcphi(c,:) = 0._rk8
      mss_cnc_bcpho(c,:) = 0._rk8

      mss_oc_top(c)      = 0._rk8
      mss_oc_col(c)      = 0._rk8
      mss_ocpho(c,:)     = 0._rk8
      mss_ocphi(c,:)     = 0._rk8
      mss_octot(c,:)     = 0._rk8
      mss_cnc_ocphi(c,:) = 0._rk8
      mss_cnc_ocpho(c,:) = 0._rk8

      mss_dst_top(c)     = 0._rk8
      mss_dst_col(c)     = 0._rk8
      mss_dst1(c,:)      = 0._rk8
      mss_dst2(c,:)      = 0._rk8
      mss_dst3(c,:)      = 0._rk8
      mss_dst4(c,:)      = 0._rk8
      mss_dsttot(c,:)    = 0._rk8
      mss_cnc_dst1(c,:)  = 0._rk8
      mss_cnc_dst2(c,:)  = 0._rk8
      mss_cnc_dst3(c,:)  = 0._rk8
      mss_cnc_dst4(c,:)  = 0._rk8

      ! top-layer diagnostics
      ! (spval is not averaged when computing history fields)
      snot_top(c)        = spval
      dTdz_top(c)        = spval
      snw_rds_top(c)     = spval
      sno_liq_top(c)     = spval
    end do
  end subroutine Hydrology2

end module mod_clm_hydrology2
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
