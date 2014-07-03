module mod_clm_soiltemperature
  !
  ! Calculates snow and soil temperatures including phase change
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams , only : dtsrf
  use mod_clm_type
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varcon , only : sb, capr, cnfac, hvap, isturb,  &
                    icol_roof, icol_sunwall, icol_shadewall,  &
                    icol_road_perv, icol_road_imperv, istwet, &
                    denh2o, denice, cpice, cpliq, hfus, tfrz, &
                    istice, istsoil, istcrop, grav, tkwat,    &
                    tkice, tkair, thk_bedrock
  use mod_clm_varpar , only : nlevsno, nlevgrnd, nlevurb, nlevsoi, &
                    max_pft_per_col
  use mod_clm_banddiagonal , only : BandDiagonal

  implicit none

  private

  save

  ! Snow and soil temperatures including phase change
  public :: SoilTemperature

  ! Set therm conduct. and heat cap of snow/soil layers
  private :: SoilThermProp
  ! When surface water freezes move ice to bottom snow layer
  private :: PhaseChangeH2osfc
  ! Calculation of the phase change within snow and soil layers
  private :: PhaseChange_beta

  contains
  !
  ! Snow and soil temperatures including phase change
  ! o The volumetric heat capacity is calculated as a linear combination
  !   in terms of the volumetric fraction of the constituent phases.
  ! o The thermal conductivity of soil is computed from
  !   the algorithm of Johansen (as reported by Farouki 1981), and the
  !   conductivity of snow is from the formulation used in
  !   SNTHERM (Jordan 1991).
  ! o Boundary conditions:
  !   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
  ! o Soil / snow temperature is predicted from heat conduction
  !   in 10 soil layers and up to 5 snow layers.
  !   The thermal conductivities at the interfaces between two
  !   neighboring layers (j, j+1) are derived from an assumption that
  !   the flux across the interface is equal to that from the node j
  !   to the interface and the flux from the interface to the node j+1.
  !   The equation is solved using the Crank-Nicholson method and
  !   results in a tridiagonal system equation.
  !
  subroutine SoilTemperature(lbl, ubl, lbc, ubc, num_urbanl, filter_urbanl, &
                             num_nolakec, filter_nolakec, xmf, fact,        &
                             c_h2osfc, xmf_h2osfc)
    implicit none
    ! column bounds
    integer(ik4) , intent(in) :: lbc , ubc
    ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: num_nolakec
    ! column filter for non-lake points
    integer(ik4) , intent(in) :: filter_nolakec(ubc-lbc+1)
    ! landunit-index bounds
    integer(ik4) , intent(in) :: lbl, ubl
    ! number of urban landunits
    integer(ik4) , intent(in) :: num_urbanl
    ! urban landunit filter
    integer(ik4) , intent(in)  :: filter_urbanl(ubl-lbl+1)
    ! total latent heat of phase change of ground water
    real(rk8), intent(out) :: xmf(lbc:ubc)
    ! used in computing tridiagonal matrix
    real(rk8), intent(out) :: fact(lbc:ubc, -nlevsno+1:nlevgrnd)
    !latent heat of phase change of surface water
    real(rk8), intent(out) :: xmf_h2osfc(lbc:ubc)
    !heat capacity of surface water
    real(rk8), intent(out) :: c_h2osfc(lbc:ubc)

    real(rk8), pointer :: snow_depth(:)    ! snow height (m)
    ! eff. fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno_eff(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno(:)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer :: frac_h2osfc(:)
    ! surface water (mm)
    real(rk8), pointer :: h2osfc(:)
    ! surface water temperature
    real(rk8), pointer :: t_h2osfc(:)
    ! saved surface water temperature
    real(rk8), pointer :: t_h2osfc_bef(:)
    ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_snow(:)
    ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_soil(:)
    ! sensible heat flux from surface water (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_h2osfc(:)
    ! evaporation flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_snow(:)
    ! evaporation flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_soil(:)
    ! evaporation flux from h2osfc (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_h2osfc(:)
    ! solar radiation absorbed by soil (W/m**2)
    real(rk8), pointer :: sabg_soil(:)
    ! solar radiation absorbed by snow (W/m**2)
    real(rk8), pointer :: sabg_snow(:)
    ! sum of soil/snow using current fsno, for balance check
    real(rk8), pointer :: sabg_chk(:)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: pgridcell(:) ! pft's gridcell index
    integer(ik4) , pointer :: plandunit(:) ! pft's landunit index
    integer(ik4) , pointer :: clandunit(:) ! column's landunit
    integer(ik4) , pointer :: ltype(:)     ! landunit type
    integer(ik4) , pointer :: ctype(:)     ! column type
    integer(ik4) , pointer :: npfts(:)     ! column's number of pfts
    integer(ik4) , pointer :: pfti(:)      ! column's beginning pft index
    real(rk8), pointer :: pwtcol(:)        ! weight of pft relative to column
    ! downward infrared (longwave) radiation (W/m**2)
    real(rk8), pointer :: forc_lwrad(:)
    integer(ik4) , pointer :: snl(:)    ! number of snow layers
    ! latent heat of vapor of water (or sublimation) [j/kg]
    real(rk8), pointer :: htvp(:)
    real(rk8), pointer :: emg(:)     ! ground emissivity
    ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rk8), pointer :: cgrnd(:)
    ! downward longwave radiation blow the canopy [W/m2]
    real(rk8), pointer :: dlrad(:)
    ! solar radiation absorbed by ground (W/m**2)
    real(rk8), pointer :: sabg(:)
    ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
    integer(ik4) , pointer :: frac_veg_nosno(:)
    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_grnd(:)
    ! soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer :: qflx_evap_soi(:)
    ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8), pointer :: qflx_tran_veg(:)
    ! interface level below a "z" level (m)
    real(rk8), pointer :: zi(:,:)
    real(rk8), pointer :: dz(:,:)       ! layer depth (m)
    real(rk8), pointer :: z(:,:)        ! layer thickness (m)
    real(rk8), pointer :: t_soisno(:,:) ! soil temperature (Kelvin)
    ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8), pointer :: eflx_lwrad_net(:)
    ! temperature at previous time step [K]
    real(rk8), pointer :: tssbef(:,:)
    ! internal building temperature (K)
    real(rk8), pointer :: t_building(:)
    ! maximum internal building temperature (K)
    real(rk8), pointer :: t_building_max(:)
    ! minimum internal building temperature (K)
    real(rk8), pointer :: t_building_min(:)
    ! soil heat content (MJ/m2)
    real(rk8), pointer :: hc_soi(:)
    ! soil plus snow plus lake heat content (MJ/m2)
    real(rk8), pointer :: hc_soisno(:)
    ! heat flux between soil layer 1 and 2 (W/m2)
    real(rk8), pointer :: eflx_fgr12(:)
    ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
    real(rk8), pointer :: eflx_fgr(:,:)
    ! traffic sensible heat flux (W/m**2)
    real(rk8), pointer :: eflx_traffic(:)
    ! sensible heat flux from urban heating/cooling sources
    ! of waste heat (W/m**2)
    real(rk8), pointer :: eflx_wasteheat(:)
    ! sensible heat flux from urban heating/cooling sources
    ! of waste heat (W/m**2)
    real(rk8), pointer :: eflx_wasteheat_pft(:)
    !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(rk8), pointer :: eflx_heat_from_ac(:)
    !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(rk8), pointer :: eflx_heat_from_ac_pft(:)
    ! traffic sensible heat flux (W/m**2)
    real(rk8), pointer :: eflx_traffic_pft(:)
    ! total anthropogenic heat flux (W/m**2)
    real(rk8), pointer :: eflx_anthro(:)
    ! urban canyon height to width ratio
    real(rk8), pointer :: canyon_hwr(:)
    ! weight of roof with respect to landunit
    real(rk8), pointer :: wtlunit_roof(:)
    ! heat flux from beneath column (W/m**2) [+ = upward]
    real(rk8), pointer :: eflx_bot(:)

    real(rk8), pointer :: t_grnd(:)   ! ground surface temperature [K]

    ! net ground heat flux into the surface (W/m**2)
    real(rk8), pointer :: eflx_gnet(:)
    ! temperature derivative of ground net heat flux
    real(rk8), pointer :: dgnetdT(:)
    ! heat flux from urban building interior to walls, roof (W/m**2)
    real(rk8), pointer :: eflx_building_heat(:)

    ! variables needed for SNICAR
    ! absorbed solar radiation (pft,lyr) [W/m2]
    real(rk8), pointer :: sabg_lyr(:,:)
    real(rk8), pointer :: h2osno(:)        ! total snow water (col) [kg/m2]
    real(rk8), pointer :: h2osoi_liq(:,:)  ! liquid water (col,lyr) [kg/m2]
    real(rk8), pointer :: h2osoi_ice(:,:)  ! ice content (col,lyr) [kg/m2]

    ! Urban building HAC fluxes
    ! urban air conditioning flux (W/m**2)
    real(rk8), pointer :: eflx_urban_ac(:)
    ! urban heating flux (W/m**2)
    real(rk8), pointer :: eflx_urban_heat(:)

    integer(ik4)  :: j,c,p,l,g,pi   !  indices
    integer(ik4)  :: fc             ! lake filtered column indices
    integer(ik4)  :: fl             ! urban filtered landunit indices
    integer(ik4)  :: jtop(lbc:ubc)  ! top level at each column
    ! "a" vector for tridiagonal matrix
    real(rk8) :: at (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! "b" vector for tridiagonal matrix
    real(rk8) :: bt (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! "c" vector for tridiagonal matrix
    real(rk8) :: ct (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! "r" vector for tridiagonal solution
    real(rk8) :: rt (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! heat capacity [J/(m2 K)]
    real(rk8) :: cv (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! thermal conductivity [W/(m K)]
    real(rk8) :: tk (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! heat diffusion through the layer interface [W/m2]
    real(rk8) :: fn (lbc:ubc,-nlevsno+1:nlevgrnd)
    ! heat diffusion through the layer interface [W/m2]
    real(rk8) :: fn1(lbc:ubc,-nlevsno+1:nlevgrnd)
    real(rk8) :: dzm   ! used in computing tridiagonal matrix
    real(rk8) :: dzp   ! used in computing tridiagonal matrix
    real(rk8) :: hs(lbc:ubc)  ! net energy flux into the surface (w/m2)
    real(rk8), pointer :: dhsdT(:) ! d(hs)/dT
    real(rk8) :: lwrad_emit(lbc:ubc) ! emitted longwave radiation
    ! time derivative of emitted longwave radiation
    real(rk8) :: dlwrad_emit(lbc:ubc)
    ! index of top layer of snowpack (-4 to 0) [idx]
    integer(ik4)  :: lyr_top
    ! absorbed solar radiation (col,lyr) [W/m2]
    real(rk8) :: sabg_lyr_col(lbc:ubc,-nlevsno+1:1)
    ! net energy flux into surface layer, pft-level [W/m2]
    real(rk8) :: eflx_gnet_top
    ! net energy flux into surface layer (col) [W/m2]
    real(rk8) :: hs_top(lbc:ubc)
    ! is urban air conditioning on?
    logical  :: cool_on(lbl:ubl)
    ! is urban heating on?
    logical  :: heat_on(lbl:ubl)
    real(rk8), pointer :: tk_h2osfc(:)
    real(rk8) :: fn_h2osfc(lbc:ubc)
    real(rk8) :: fn1_h2osfc(lbc:ubc)
    real(rk8) :: dz_h2osfc(lbc:ubc)
    integer(ik4), parameter :: nband=5
    real(rk8),pointer :: bmatrix(:,:,:)
    real(rk8),pointer :: tvector(:,:)
    real(rk8),pointer :: rvector(:,:)
    real(rk8),dimension(lbc:ubc) :: hs_snow
    real(rk8),dimension(lbc:ubc) :: hs_soil
    real(rk8),dimension(lbc:ubc) :: hs_top_snow
    real(rk8),dimension(lbc:ubc) :: hs_top_soil
    real(rk8),dimension(lbc:ubc) :: hs_h2osfc
    real(rk8),dimension(lbc:ubc)  :: lwrad_emit_snow
    real(rk8),dimension(lbc:ubc)  :: lwrad_emit_soil
    real(rk8),dimension(lbc:ubc)  :: lwrad_emit_h2osfc
    real(rk8) :: eflx_gnet_snow
    real(rk8) :: eflx_gnet_soil
    real(rk8) :: eflx_gnet_h2osfc
    integer(ik4)  :: jbot(lbc:ubc) ! bottom level at each column

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_lwrad     => clm_a2l%forc_lwrad

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype          => clm3%g%l%itype
    t_building     => clm3%g%l%lps%t_building
    t_building_max => clm3%g%l%lps%t_building_max
    t_building_min => clm3%g%l%lps%t_building_min
    eflx_traffic   => clm3%g%l%lef%eflx_traffic
    canyon_hwr     => clm3%g%l%canyon_hwr
    eflx_wasteheat => clm3%g%l%lef%eflx_wasteheat
    eflx_heat_from_ac => clm3%g%l%lef%eflx_heat_from_ac
    wtlunit_roof   => clm3%g%l%wtlunit_roof

    ! Assign local pointers to derived subtypes components (column-level)

    frac_sno_eff   => clm3%g%l%c%cps%frac_sno_eff
    frac_sno       => clm3%g%l%c%cps%frac_sno
    frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc
    h2osfc         => clm3%g%l%c%cws%h2osfc
    t_h2osfc       => clm3%g%l%c%ces%t_h2osfc
    t_h2osfc_bef   => clm3%g%l%c%ces%t_h2osfc_bef
    snow_depth     => clm3%g%l%c%cps%snow_depth
    eflx_sh_snow   => clm3%g%l%c%p%pef%eflx_sh_snow
    eflx_sh_soil   => clm3%g%l%c%p%pef%eflx_sh_soil
    eflx_sh_h2osfc => clm3%g%l%c%p%pef%eflx_sh_h2osfc
    qflx_ev_snow   => clm3%g%l%c%p%pwf%qflx_ev_snow
    qflx_ev_soil   => clm3%g%l%c%p%pwf%qflx_ev_soil
    qflx_ev_h2osfc => clm3%g%l%c%p%pwf%qflx_ev_h2osfc
    sabg_soil      => clm3%g%l%c%p%pef%sabg_soil
    sabg_snow      => clm3%g%l%c%p%pef%sabg_snow
    sabg_chk       => clm3%g%l%c%p%pef%sabg_chk
    ctype          => clm3%g%l%c%itype
    clandunit      => clm3%g%l%c%landunit
    npfts          => clm3%g%l%c%npfts
    pfti           => clm3%g%l%c%pfti
    snl            => clm3%g%l%c%cps%snl
    htvp           => clm3%g%l%c%cps%htvp
    emg            => clm3%g%l%c%cps%emg
    t_grnd         => clm3%g%l%c%ces%t_grnd
    hc_soi         => clm3%g%l%c%ces%hc_soi
    hc_soisno      => clm3%g%l%c%ces%hc_soisno
    eflx_fgr12     => clm3%g%l%c%cef%eflx_fgr12
    eflx_fgr       => clm3%g%l%c%cef%eflx_fgr
    zi             => clm3%g%l%c%cps%zi
    dz             => clm3%g%l%c%cps%dz
    z              => clm3%g%l%c%cps%z
    t_soisno       => clm3%g%l%c%ces%t_soisno
    eflx_building_heat => clm3%g%l%c%cef%eflx_building_heat
    tssbef             => clm3%g%l%c%ces%tssbef
    eflx_urban_ac      => clm3%g%l%c%cef%eflx_urban_ac
    eflx_urban_heat    => clm3%g%l%c%cef%eflx_urban_heat
    eflx_bot           => clm3%g%l%c%cef%eflx_bot

    ! Assign local pointers to derived subtypes components (pft-level)

    pactive        => clm3%g%l%c%p%active
    pgridcell      => clm3%g%l%c%p%gridcell
    plandunit      => clm3%g%l%c%p%landunit
    pwtcol         => clm3%g%l%c%p%wtcol
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    cgrnd          => clm3%g%l%c%p%pef%cgrnd
    dlrad          => clm3%g%l%c%p%pef%dlrad
    sabg           => clm3%g%l%c%p%pef%sabg
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
    eflx_gnet      => clm3%g%l%c%p%pef%eflx_gnet
    dgnetdT        => clm3%g%l%c%p%pef%dgnetdT
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net
    eflx_wasteheat_pft => clm3%g%l%c%p%pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => clm3%g%l%c%p%pef%eflx_heat_from_ac_pft
    eflx_traffic_pft => clm3%g%l%c%p%pef%eflx_traffic_pft
    eflx_anthro => clm3%g%l%c%p%pef%eflx_anthro

    sabg_lyr       => clm3%g%l%c%p%pef%sabg_lyr
    h2osno         => clm3%g%l%c%cws%h2osno
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice

    ! Compute ground surface and soil temperatures

    ! Thermal conductivity and Heat capacity

    allocate( tk_h2osfc(lbc:ubc) )
    tk_h2osfc(:) = nan
    allocate( dhsdT(lbc:ubc) )

    call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv, tk_h2osfc)

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      lwrad_emit(c)  = emg(c) * sb * t_grnd(c)**4
      dlwrad_emit(c) = 4.D0*emg(c) * sb * t_grnd(c)**3
      ! fractionate lwrad_emit; balanced in CanopyFluxes & Biogeophysics2
      lwrad_emit_snow(c)    = emg(c) * sb * t_soisno(c,snl(c)+1)**4
      lwrad_emit_soil(c)    = emg(c) * sb * t_soisno(c,1)**4
      lwrad_emit_h2osfc(c)  = emg(c) * sb * t_h2osfc(c)**4
    end do

    hs_snow(lbc:ubc)   = 0.D0
    hs_soil(lbc:ubc)   = 0.D0
    hs_h2osfc(lbc:ubc) = 0.D0
    hs(lbc:ubc)        = 0.D0
    dhsdT(lbc:ubc)     = 0.D0
    do pi = 1,max_pft_per_col
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if ( pi <= npfts(c) ) then
          p = pfti(c) + pi - 1
          l = plandunit(p)
          g = pgridcell(p)

          if (pactive(p)) then
            if (ltype(l) /= isturb) then
              eflx_gnet(p) = sabg(p) + dlrad(p) + &
                     (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                     lwrad_emit(c) - &
                     (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
              ! save sabg for balancecheck, in case frac_sno is set
              ! to zero later
              sabg_chk(p) = frac_sno_eff(c) * sabg_snow(p) + &
                (1.D0 - frac_sno_eff(c) ) * sabg_soil(p)

              eflx_gnet_snow = sabg_snow(p) + dlrad(p) + &
                      dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                      lwrad_emit_snow(c) - &
                      (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

              eflx_gnet_soil = sabg_soil(p) + dlrad(p) + &
                      dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                      lwrad_emit_soil(c) - &
                      (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

              eflx_gnet_h2osfc = sabg_soil(p) + dlrad(p) + &
                      dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                      lwrad_emit_h2osfc(c) - &
                      (eflx_sh_h2osfc(p)+qflx_ev_h2osfc(p)*htvp(c))
            else
              ! For urban columns we use the net longwave radiation
              ! (eflx_lwrad_net) because of interactions between urban columns.
              ! All wasteheat and traffic flux goes into canyon floor
              if (ctype(c) == icol_road_perv .or. &
                  ctype(c) == icol_road_imperv) then
                eflx_wasteheat_pft(p) = eflx_wasteheat(l) / &
                                     (1.D0-wtlunit_roof(l))
                eflx_heat_from_ac_pft(p) = eflx_heat_from_ac(l) / &
                                     (1.D0-wtlunit_roof(l))
                eflx_traffic_pft(p) = eflx_traffic(l)/(1.D0-wtlunit_roof(l))
              else
                eflx_wasteheat_pft(p) = 0.D0
                eflx_heat_from_ac_pft(p) = 0.D0
                eflx_traffic_pft(p) = 0.D0
              end if
              ! Include transpiration term because needed for previous road
              ! and include wasteheat and traffic flux
              eflx_gnet(p) = sabg(p) + dlrad(p)  - eflx_lwrad_net(p) - &
                             (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + &
                             qflx_tran_veg(p)*hvap) + &
                             eflx_wasteheat_pft(p) +  &
                             eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
              eflx_anthro(p) = eflx_wasteheat_pft(p) + eflx_traffic_pft(p)
              eflx_gnet_snow = eflx_gnet(p)
              eflx_gnet_soil = eflx_gnet(p)
              eflx_gnet_h2osfc = eflx_gnet(p)
            end if
            dgnetdT(p) = - cgrnd(p) - dlwrad_emit(c)
            hs(c) = hs(c) + eflx_gnet(p) * pwtcol(p)
            dhsdT(c) = dhsdT(c) + dgnetdT(p) * pwtcol(p)
            ! separate surface fluxes for soil/snow
            hs_snow(c) = hs_snow(c) + eflx_gnet_snow * pwtcol(p)
            hs_soil(c) = hs_soil(c) + eflx_gnet_soil * pwtcol(p)
            hs_h2osfc(c) = hs_h2osfc(c) + eflx_gnet_h2osfc * pwtcol(p)
          end if
        end if
      end do
    end do

    ! Additional calculations with SNICAR:
    ! Set up tridiagonal matrix in a new manner. There is now
    ! absorbed solar radiation in each snow layer, instead of
    ! only the surface. Following the current implementation,
    ! absorbed solar flux should be: S + ((delS/delT)*dT),
    ! where S is absorbed radiation, and T is temperature. Now,
    ! assume delS/delT is zero, then it is OK to just add S
    ! to each layer

    ! Initialize:
    sabg_lyr_col(lbc:ubc,-nlevsno+1:1) = 0.D0
    hs_top(lbc:ubc)      = 0.D0
    hs_top_snow(lbc:ubc) = 0.D0
    hs_top_soil(lbc:ubc) = 0.D0

    do pi = 1 , max_pft_per_col
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        lyr_top = snl(c) + 1
        if ( pi <= npfts(c) ) then
          p = pfti(c) + pi - 1
          if (pactive(p)) then
            g = pgridcell(p)
            l = plandunit(p)
            if (ltype(l) /= isturb ) then
              eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) + &
                dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                lwrad_emit(c) - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

              hs_top(c) = hs_top(c) + eflx_gnet_top*pwtcol(p)

              eflx_gnet_snow = sabg_lyr(p,lyr_top) + dlrad(p) + &
                dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                lwrad_emit_snow(c) - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

              eflx_gnet_soil = sabg_lyr(p,lyr_top) + dlrad(p) + &
                dble(1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) - &
                lwrad_emit_soil(c) - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

              hs_top_snow(c) = hs_top_snow(c) + eflx_gnet_snow*pwtcol(p)
              hs_top_soil(c) = hs_top_soil(c) + eflx_gnet_soil*pwtcol(p)

              do j = lyr_top,1,1
                sabg_lyr_col(c,j) = sabg_lyr_col(c,j) + &
                  sabg_lyr(p,j) * pwtcol(p)
              end do
            else
              hs_top(c)      = hs_top(c) + eflx_gnet(p)*pwtcol(p)
              hs_top_snow(c) = hs_top_snow(c) + eflx_gnet(p)*pwtcol(p)
              hs_top_soil(c) = hs_top_soil(c) + eflx_gnet(p)*pwtcol(p)
              sabg_lyr_col(c,lyr_top) = sabg_lyr_col(c,lyr_top) + &
                sabg(p) * pwtcol(p)
            end if
          end if
        end if
      end do
    end do

    ! Restrict internal building temperature to between min and max
    ! and determine if heating or air conditioning is on

    do fl = 1,num_urbanl
      l = filter_urbanl(fl)
      if (ltype(l) == isturb) then
        cool_on(l) = .false.
        heat_on(l) = .false.
        if (t_building(l) > t_building_max(l)) then
          t_building(l) = t_building_max(l)
          cool_on(l) = .true.
          heat_on(l) = .false.
        else if (t_building(l) < t_building_min(l)) then
          t_building(l) = t_building_min(l)
          cool_on(l) = .false.
          heat_on(l) = .true.
        end if
      end if
    end do

    ! Determine heat diffusion through the layer interface and factor
    ! used in computing tridiagonal matrix and set up vector r and
    ! vectors a, b, c that define tridiagonal matrix and solve system

    do j = -nlevsno+1 , nlevgrnd
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if ((ctype(c) == icol_sunwall .or. &
             ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof) .and. j <= nlevurb) then
          if (j >= snl(c)+1) then
            if (j == snl(c)+1) then
              fact(c,j) = dtsrf/cv(c,j)
              fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j)) / &
                (z(c,j+1)-z(c,j))
            else if (j <= nlevurb-1) then
              fact(c,j) = dtsrf/cv(c,j)
              fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j)) / &
                (z(c,j+1)-z(c,j))
              dzm     = (z(c,j)-z(c,j-1))
            else if (j == nlevurb) then
              fact(c,j) = dtsrf/cv(c,j)
              ! For urban sunwall, shadewall, and roof columns,
              ! there is a non-zero heat flux across the bottom
              ! "soil" layer and the equations are derived assuming
              ! a prescribed internal building temperature.
              ! (See Oleson urban notes of 6/18/03).
              fn(c,j) = tk(c,j) * &
                (t_building(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
            end if
          end if
        else if (ctype(c) /= icol_sunwall .and. &
                 ctype(c) /= icol_shadewall .and. &
                 ctype(c) /= icol_roof) then
          if (j >= snl(c)+1) then
            if (j == snl(c)+1) then
              fact(c,j) = dtsrf/cv(c,j) * dz(c,j) / &
                (0.5D0*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
              fn(c,j) = tk(c,j) * &
                (t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
            else if (j <= nlevgrnd-1) then
              fact(c,j) = dtsrf/cv(c,j)
              fn(c,j) = tk(c,j) * &
                (t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
              dzm = (z(c,j)-z(c,j-1))
            else if (j == nlevgrnd) then
              fact(c,j) = dtsrf/cv(c,j)
              fn(c,j) = eflx_bot(c)
            end if
          end if
        end if
      end do
    end do
    !
    ! urban aboveground columns
    !
    do j = -nlevsno+1 , nlevurb
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if ((ctype(c) == icol_sunwall .or. &
             ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof)) then
          if (j >= snl(c)+1) then
            if (j == snl(c)+1) then
              dzp = z(c,j+1)-z(c,j)
              at(c,j) = 0.D0
              bt(c,j) = 1.0D0+(1.D0-cnfac)*fact(c,j)*tk(c,j) / &
                         dzp-fact(c,j)*dhsdT(c)
              ct(c,j) = -(1.D0-cnfac)*fact(c,j)*tk(c,j)/dzp
              ! changed hs to hs_top
              rt(c,j) = t_soisno(c,j) + &
                fact(c,j)*(hs_top(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
            else if (j <= nlevurb-1) then
              dzm = (z(c,j)-z(c,j-1))
              dzp = (z(c,j+1)-z(c,j))
              at(c,j) = -(1.D0-cnfac)*fact(c,j)* tk(c,j-1)/dzm
              bt(c,j) = 1.D0+(1.D0-cnfac)*fact(c,j)*(tk(c,j)/dzp+tk(c,j-1)/dzm)
              ct(c,j) = -(1.D0-cnfac)*fact(c,j)* tk(c,j)/dzp

              ! if this is a snow layer or the top soil layer,
              ! add absorbed solar flux to factor 'rt'
              if (j <= 1) then
                rt(c,j) = t_soisno(c,j) + &
                  cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) ) + &
                       (fact(c,j)*sabg_lyr_col(c,j))
              else
                rt(c,j) = t_soisno(c,j) + &
                  cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
              end if
            else if (j == nlevurb) then
              ! For urban sunwall, shadewall, and roof columns,
              ! there is a non-zero heat flux across the bottom
              ! "soil" layer and the equations are derived assuming
              ! a prescribed internal building temperature.
              ! (See Oleson urban notes of 6/18/03).
              dzm = ( z(c,j)-z(c,j-1))
              dzp = (zi(c,j)-z(c,j))
              at(c,j) = - (1.D0-cnfac)*fact(c,j)*(tk(c,j-1)/dzm)
              bt(c,j) = 1.D0+(1.D0-cnfac)*fact(c,j)*(tk(c,j-1)/dzm+tk(c,j)/dzp)
              ct(c,j) = 0.D0
              rt(c,j) = t_soisno(c,j) + fact(c,j)*( fn(c,j) - cnfac*fn(c,j-1) )
            end if
          end if
        end if
      end do
    end do
    !
    ! non-urban landunits + urban road columns
    !
    do j = -nlevsno+1 , nlevgrnd
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        if (ctype(c) /= icol_sunwall .and. &
            ctype(c) /= icol_shadewall .and. &
            ctype(c) /= icol_roof) then
          if (j >= snl(c)+1) then
            if (j == snl(c)+1) then
              dzp = z(c,j+1)-z(c,j)
              at(c,j) = 0.D0
              bt(c,j) = 1.0D0+(1.D0-cnfac)*fact(c,j)*tk(c,j) / &
                dzp-fact(c,j)*dhsdT(c)
              ct(c,j) = -(1.D0-cnfac)*fact(c,j)*tk(c,j)/dzp
              rt(c,j) = t_soisno(c,j) + fact(c,j)*( hs_top_snow(c) - &
                        dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
            else if (j == 1) then
              ! this is the snow/soil interface layer
              dzm = (z(c,j)-z(c,j-1))
              dzp = (z(c,j+1)-z(c,j))
              at(c,j) =  -frac_sno_eff(c) * (1.D0-cnfac) * fact(c,j) * &
                         tk(c,j-1)/dzm
              bt(c,j) = 1.D0 + (1.D0-cnfac)*fact(c,j)*(tk(c,j)/dzp + &
                        frac_sno_eff(c) * tk(c,j-1)/dzm) - &
                        (1.D0 - frac_sno_eff(c))*fact(c,j)*dhsdT(c)
              ct(c,j) = -(1.D0-cnfac)*fact(c,j)*tk(c,j)/dzp
              rt(c,j) = t_soisno(c,j) + fact(c,j) * &
                    ((1.D0-frac_sno_eff(c)) * &
                    (hs_soil(c) - dhsdT(c)*t_soisno(c,j)) + &
                     cnfac*(fn(c,j) - frac_sno_eff(c) * fn(c,j-1)))
              rt(c,j) = rt(c,j) + frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(c,j)
            else if (j <= nlevgrnd-1) then
              dzm = (z(c,j)-z(c,j-1))
              dzp = (z(c,j+1)-z(c,j))
              at(c,j) = -(1.D0-cnfac)*fact(c,j)* tk(c,j-1)/dzm
              bt(c,j) = 1.D0+(1.D0-cnfac) * &
                fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
              ct(c,j) =  -(1.D0-cnfac)*fact(c,j)* tk(c,j)/dzp
              rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
              if (j < 1) rt(c,j) = rt(c,j) + fact(c,j)*sabg_lyr_col(c,j)
            else if (j == nlevgrnd) then
              dzm = (z(c,j)-z(c,j-1))
              at(c,j) = -(1.D0-cnfac)*fact(c,j)*tk(c,j-1)/dzm
              bt(c,j) = 1.D0+(1.D0-cnfac)*fact(c,j)*tk(c,j-1)/dzm
              ct(c,j) = 0.D0
              rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j) * &
                fn(c,j-1) + fact(c,j)*fn(c,j)
            end if
          end if
        end if
      end do
    end do

    ! compute thermal properties of h2osfc
    do fc = 1 , num_nolakec
      c = filter_nolakec(fc)
      dz_h2osfc(c) = max(1.0D-6,1.0e-3*h2osfc(c))
      !"areametric" heat capacity [J/K/m^2]
      c_h2osfc(c) = cpliq*denh2o*dz_h2osfc(c)
    end do

    ! set up compact matrix for band diagonal solver, requires additional
    !     sub/super diagonals (1 each), and one additional row for t_h2osfc
    jtop = -9999
    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      jtop(c) = snl(c)
      ! compute jbot
      if ((ctype(c) == icol_sunwall .or. &
           ctype(c) == icol_shadewall .or. &
           ctype(c) == icol_roof) ) then
        jbot(c) = nlevurb
      else
        jbot(c) = nlevgrnd
      end if
    end do

    ! allocate matrices for BandDiagonal
    allocate(bmatrix(lbc:ubc,nband,-nlevsno:nlevgrnd))
    bmatrix(:,:,:) = 0.0
    allocate(tvector(lbc:ubc,-nlevsno:nlevgrnd))
    tvector(:,:) = nan
    allocate(rvector(lbc:ubc,-nlevsno:nlevgrnd))
    rvector(:,:) = nan

    ! the solution will be organized as (snow:h2osfc:soil) to minimize
    !     bandwidth; this requires a 5-element band instead of 3
    do fc = 1 , num_nolakec
      c = filter_nolakec(fc)
      do j = snl(c)+1, 0
        ! snow layers; bottom layer will have one offset coefficient
        bmatrix(c,2,j-1) = ct(c,j)
        bmatrix(c,3,j-1) = bt(c,j)
        bmatrix(c,4,j-1) = at(c,j)
        rvector(c,j-1) = rt(c,j)
        tvector(c,j-1) = t_soisno(c,j)
      end do
      ! bottom snow layer has super coef shifted to 2nd super diagonal
      bmatrix(c,2,-1) = 0.0D0
      bmatrix(c,1,-1) = ct(c,0) !flux to top soil layer
      ! surface water layer has two coefficients
      dzm = (0.5D0*dz_h2osfc(c)+z(c,1))
      !flux to top soil layer
      bmatrix(c,2,0)= -(1.D0-cnfac)*(dtsrf/c_h2osfc(c))*tk_h2osfc(c)/dzm
      !interaction from atm
      bmatrix(c,3,0)= 1.0D0+(1.D0-cnfac)*(dtsrf/c_h2osfc(c)) * &
            tk_h2osfc(c)/dzm -(dtsrf/c_h2osfc(c))*dhsdT(c)
      fn_h2osfc(c) = tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzm
      !rhs for h2osfc
      rvector(c,0) = t_h2osfc(c) +  (dtsrf/c_h2osfc(c)) * &
            ( hs_h2osfc(c) - dhsdT(c)*t_h2osfc(c) + cnfac*fn_h2osfc(c) )
      tvector(c,0) = t_h2osfc(c)
      ! soil layers; top layer will have one offset and one extra coefficient
      bmatrix(c,2,1:nlevgrnd) = ct(c,1:nlevgrnd)
      bmatrix(c,3,1:nlevgrnd) = bt(c,1:nlevgrnd)
      bmatrix(c,4,1:nlevgrnd) = at(c,1:nlevgrnd)
      ! top soil layer has sub coef shifted to 2nd super diagonal
      if ( frac_h2osfc(c) /= 0.0D0 )then
        !flux from h2osfc
        bmatrix(c,4,1)=  -frac_h2osfc(c) * (1.D0-cnfac) * fact(c,1) * &
                tk_h2osfc(c)/dzm
      else
        bmatrix(c,4,1) = 0.0D0
      end if
      bmatrix(c,5,1) = at(c,1)
      ! diagonal element correction for presence of h2osfc
      if ( frac_h2osfc(c) /= 0.0D0 )then
        bmatrix(c,3,1) = bmatrix(c,3,1)+ frac_h2osfc(c) * &
              ((1.D0-cnfac)*fact(c,1)*tk_h2osfc(c)/dzm + fact(c,1)*dhsdT(c))
      end if
      rvector(c,1:nlevgrnd) = rt(c,1:nlevgrnd)
      ! rhs correction for h2osfc
      if ( frac_h2osfc(c) /= 0.0D0 )then
        rvector(c,1) = rvector(c,1) - frac_h2osfc(c)*fact(c,1) * &
          ((hs_soil(c) - dhsdT(c)*t_soisno(c,1))+cnfac*fn_h2osfc(c))
      end if
      tvector(c,1:nlevgrnd) = t_soisno(c,1:nlevgrnd)
    end do

    call BandDiagonal(lbc, ubc, -nlevsno, nlevgrnd, jtop, jbot, num_nolakec, &
         filter_nolakec, nband, bmatrix, rvector, tvector)

    ! return temperatures to original array
    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      do j = snl(c)+1, 0
        t_soisno(c,j) = tvector(c,j-1) !snow layers
      end do
      t_soisno(c,1:nlevgrnd) = tvector(c,1:nlevgrnd)  !soil layers
      if ( frac_h2osfc(c) == 0.D0 ) then
        t_h2osfc(c) = t_soisno(c,1)
      else
        t_h2osfc(c) = tvector(c,0)           !surface water
      end if
    end do
    deallocate(bmatrix)
    deallocate(tvector)
    deallocate(rvector)

    ! Melting or Freezing

    do j = -nlevsno+1,nlevgrnd
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if ((ctype(c) == icol_sunwall .or. &
             ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof) .and. j <= nlevurb) then
          if (j >= snl(c)+1) then
            if (j <= nlevurb-1) then
              fn1(c,j) = tk(c,j) * &
                (t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
            else if (j == nlevurb) then
              ! For urban sunwall, shadewall, and roof columns,
              ! there is a non-zero heat flux across the bottom
              ! "soil" layer and the equations are derived assuming
              ! a prescribed internal building temperature.
              ! (See Oleson urban notes of 6/18/03).
              ! Note new formulation for fn, this will be used below
              ! in net energey flux computations
              fn1(c,j) = tk(c,j) * (t_building(l)-t_soisno(c,j)) / &
                                   (zi(c,j) - z(c,j))
              fn(c,j) = tk(c,j) * (t_building(l) - tssbef(c,j)) / &
                                   (zi(c,j) - z(c,j))
            end if
          end if
        else if (ctype(c) /= icol_sunwall .and. &
                 ctype(c) /= icol_shadewall .and. &
                 ctype(c) /= icol_roof) then
          if (j >= snl(c)+1) then
            if (j <= nlevgrnd-1) then
              fn1(c,j) = tk(c,j) * &
                (t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
            else if (j == nlevgrnd) then
              fn1(c,j) = 0.D0
            end if
          end if
        end if
      end do
    end do

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)
      if (ltype(l) == isturb) then
        eflx_building_heat(c) = cnfac*fn(c,nlevurb) + &
           (1.0D0-cnfac)*fn1(c,nlevurb)
        if (cool_on(l)) then
          eflx_urban_ac(c) = abs(eflx_building_heat(c))
          eflx_urban_heat(c) = 0.D0
        else if (heat_on(l)) then
          eflx_urban_ac(c) = 0.D0
          eflx_urban_heat(c) = abs(eflx_building_heat(c))
        else
          eflx_urban_ac(c) = 0.D0
          eflx_urban_heat(c) = 0.D0
        end if
      end if
    end do

    ! compute terms needed for phase change of h2osfc
    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      dzp = (0.5D0*dz_h2osfc(c)+z(c,1))
      fn1_h2osfc(c)=tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzp
    end do

    xmf_h2osfc = 0.0D0
    ! compute phase change of h2osfc
    call PhaseChangeH2osfc (lbc, ubc, num_nolakec, &
      filter_nolakec, fact, dhsdT, c_h2osfc, xmf_h2osfc)

    call Phasechange_beta (lbc, ubc, num_nolakec, &
      filter_nolakec, fact, dhsdT, xmf)

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      ! this expression will (should) work whether there is snow or not
      if (snl(c) < 0) then
        if ( frac_h2osfc(c) /= 0.D0 ) then
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) + &
               (1.0D0 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) + &
                frac_h2osfc(c) * t_h2osfc(c)
        else
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) + &
               (1.0D0 - frac_sno_eff(c)) * t_soisno(c,1)
        end if
      else
        if ( frac_h2osfc(c) /= 0.D0 ) then
          t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + &
            frac_h2osfc(c) * t_h2osfc(c)
        else
          t_grnd(c) = t_soisno(c,1)
        end if
      end if
    end do

    ! Initialize soil heat content
    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)
      if (ltype(l) /= isturb) then
        hc_soisno(c) = 0.D0
        hc_soi(c)    = 0.D0
      end if
      eflx_fgr12(c)= 0.D0
    end do

    ! Calculate soil heat content and soil plus snow heat content
    do j = -nlevsno+1,nlevgrnd
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if (j == 1) then ! this only needs to be done once
          eflx_fgr12(c) = -cnfac*fn(c,1) - (1.D0-cnfac)*fn1(c,1)
        end if
        if (j > 0 .and. j < nlevgrnd .and. &
            (ltype(l) == istsoil .or. ltype(l) == istcrop)) then
          eflx_fgr(c,j) = -cnfac*fn(c,j) - (1.D0-cnfac)*fn1(c,j)
        else if (j == nlevgrnd .and. &
            (ltype(l) == istsoil .or. ltype(l) == istcrop)) then
          eflx_fgr(c,j) = 0.D0
        end if
        if (ltype(l) /= isturb) then
          if (j >= snl(c)+1) then
            hc_soisno(c) = hc_soisno(c) + cv(c,j)*t_soisno(c,j) / 1.D6
          end if
          if (j >= 1) then
            hc_soi(c) = hc_soi(c) + cv(c,j)*t_soisno(c,j) / 1.D6
          end if
        end if
      end do
    end do
    deallocate( tk_h2osfc )
    deallocate( dhsdT )
  end subroutine SoilTemperature
  !
  ! Calculation of thermal conductivities and heat capacities of
  ! snow/soil layers
  ! (1) The volumetric heat capacity is calculated as a linear combination
  !     in terms of the volumetric fraction of the constituent phases.
  !
  ! (2) The thermal conductivity of soil is computed from the algorithm of
  !     Johansen (as reported by Farouki 1981), and of snow is from the
  !     formulation used in SNTHERM (Jordan 1991).
  ! The thermal conductivities at the interfaces between two neighboring
  ! layers (j, j+1) are derived from an assumption that the flux across
  ! the interface is equal to that from the node j to the interface and the
  ! flux from the interface to the node j+1.
  !
  subroutine SoilThermProp (lbc, ubc,  num_nolakec, filter_nolakec, &
       tk, cv,tk_h2osfc)
    implicit none
    integer(ik4) , intent(in)  :: lbc, ubc                       ! column bounds
    integer(ik4) , intent(in)  :: num_nolakec                    ! number of column non-lake points in column filter
    integer(ik4) , intent(in)  :: filter_nolakec(ubc-lbc+1)      ! column filter for non-lake points
    real(rk8), intent(out) :: cv(lbc:ubc,-nlevsno+1:nlevgrnd)! heat capacity [J/(m2 K)]
    real(rk8), intent(out) :: tk(lbc:ubc,-nlevsno+1:nlevgrnd)! thermal conductivity [W/(m K)]
    real(rk8), intent(out) :: tk_h2osfc(lbc:ubc)             ! thermal conductivity of h2osfc [W/(m K)]

    real(rk8), pointer :: h2osfc(:)        ! surface (mm H2O)
    real(rk8), pointer :: frac_sno(:)      ! fractional snow covered area
    integer(ik4) , pointer :: ctype(:)         ! column type
    integer(ik4) , pointer :: clandunit(:)     ! column's landunit
    integer(ik4) , pointer :: ltype(:)         ! landunit type
    integer(ik4) , pointer :: snl(:)           ! number of snow layers
    real(rk8), pointer :: h2osno(:)        ! snow water (mm H2O)

    real(rk8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: tksatu(:,:)      ! thermal conductivity, saturated soil [W/m-K]
    real(rk8), pointer :: tkmg(:,:)        ! thermal conductivity, soil minerals  [W/m-K]
    real(rk8), pointer :: tkdry(:,:)       ! thermal conductivity, dry soil (W/m/Kelvin)
    real(rk8), pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(rk8), pointer :: dz(:,:)          ! layer depth (m)
    real(rk8), pointer :: zi(:,:)          ! interface level below a "z" level (m)
    real(rk8), pointer :: z(:,:)           ! layer thickness (m)
    real(rk8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(rk8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(rk8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(rk8), pointer :: tk_wall(:,:)     ! thermal conductivity of urban wall
    real(rk8), pointer :: tk_roof(:,:)     ! thermal conductivity of urban roof
    real(rk8), pointer :: tk_improad(:,:)  ! thermal conductivity of urban impervious road
    real(rk8), pointer :: cv_wall(:,:)     ! thermal conductivity of urban wall
    real(rk8), pointer :: cv_roof(:,:)     ! thermal conductivity of urban roof
    real(rk8), pointer :: cv_improad(:,:)  ! thermal conductivity of urban impervious road
    integer(ik4),  pointer :: nlev_improad(:)  ! number of impervious road layers

    integer(ik4)  :: l,c,j                     ! indices
    integer(ik4)  :: fc                        ! lake filtered column indices
    real(rk8) :: bw                        ! partial density of water (ice + liquid)
    real(rk8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(rk8) :: dke                       ! kersten number
    real(rk8) :: fl                        ! volume fraction of liquid or unfrozen water to total water
    real(rk8) :: satw                      ! relative total water content of soil.
    real(rk8) :: thk(lbc:ubc,-nlevsno+1:nlevgrnd) ! thermal conductivity of layer
    real(rk8) :: zh2osfc

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype    => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    h2osfc     => clm3%g%l%c%cws%h2osfc
    frac_sno   => clm3%g%l%c%cps%frac_sno_eff
    ctype      => clm3%g%l%c%itype
    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    h2osno     => clm3%g%l%c%cws%h2osno
    watsat     => clm3%g%l%c%cps%watsat
    tksatu     => clm3%g%l%c%cps%tksatu
    tkmg       => clm3%g%l%c%cps%tkmg
    tkdry      => clm3%g%l%c%cps%tkdry
    csol       => clm3%g%l%c%cps%csol
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    tk_wall    => clm3%g%l%lps%tk_wall
    tk_roof    => clm3%g%l%lps%tk_roof
    tk_improad => clm3%g%l%lps%tk_improad
    cv_wall    => clm3%g%l%lps%cv_wall
    cv_roof    => clm3%g%l%lps%cv_roof
    cv_improad => clm3%g%l%lps%cv_improad
    nlev_improad => clm3%g%l%lps%nlev_improad

    ! Thermal conductivity of soil from Farouki (1981)
    ! Urban values are from Masson et al. 2002, Evaluation of the
    ! Town Energy Balance (TEB) scheme with direct measurements from
    ! dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = -nlevsno+1 , nlevgrnd
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)

        ! Only examine levels from 1->nlevgrnd
        if (j >= 1) then
          l = clandunit(c)
          if ((ctype(c) == icol_sunwall .or. &
               ctype(c) == icol_shadewall) .and. j <= nlevurb) then
            thk(c,j) = tk_wall(l,j)
          else if (ctype(c) == icol_roof .and. j <= nlevurb) then
            thk(c,j) = tk_roof(l,j)
          else if (ctype(c) == icol_road_imperv .and. &
                   j >= 1 .and. j <= nlev_improad(l)) then
            thk(c,j) = tk_improad(l,j)
          else if (ltype(l) /= istwet .AND.  &
                   ltype(l) /= istice  .AND. &
                   ctype(c) /= icol_sunwall .AND. &
                   ctype(c) /= icol_shadewall .AND. &
                   ctype(c) /= icol_roof) then
            satw = (h2osoi_liq(c,j)/denh2o + &
              h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
            satw = min(1.D0, satw)
            if (satw > .1D-6) then
              if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                dke = max(0.D0, log10(satw) + 1.0D0)
              else                               ! Frozen soil
                dke = satw
              end if
              fl = (h2osoi_liq(c,j)/(denh2o*dz(c,j))) / &
                    (h2osoi_liq(c,j)/(denh2o*dz(c,j)) + &
                     h2osoi_ice(c,j)/(denice*dz(c,j)))
              dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))* &
                      tkice**((1.D0-fl)*watsat(c,j))
              thk(c,j) = dke*dksat + (1.D0-dke)*tkdry(c,j)
            else
              thk(c,j) = tkdry(c,j)
            end if
            if (j > nlevsoi) thk(c,j) = thk_bedrock
          else if ( ltype(l) == istice ) then
            thk(c,j) = tkwat
            if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
          else if (ltype(l) == istwet) then
            if (j > nlevsoi) then
              thk(c,j) = thk_bedrock
            else
              thk(c,j) = tkwat
              if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
            end if
          end if
        end if

        ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
        ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
        if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
          bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/(frac_sno(c)*dz(c,j))
          thk(c,j) = tkair + (7.75D-5 *bw + 1.105D-6*bw*bw)*(tkice-tkair)
        end if
      end do
    end do

    ! Thermal conductivity at the layer interface

    do j = -nlevsno+1,nlevgrnd
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if ((ctype(c) == icol_sunwall .or. &
             ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof) .and. j <= nlevurb) then
          if (j >= snl(c)+1 .AND. j <= nlevurb-1) then
            tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) / &
                   (thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
          else if (j == nlevurb) then
            ! For urban sunwall, shadewall, and roof columns,
            ! there is a non-zero heat flux across the bottom
            ! "soil" layer and the equations are derived assuming
            ! a prescribed internal building temperature.
            ! (See Oleson urban notes of 6/18/03).
            tk(c,j) = thk(c,j)
          end if
        else if (ctype(c) /= icol_sunwall .and. &
                 ctype(c) /= icol_shadewall .and. &
                 ctype(c) /= icol_roof) then
          if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
            tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) / &
                    (thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
          else if (j == nlevgrnd) then
            tk(c,j) = 0.D0
          end if
        end if
      end do
    end do

    ! calculate thermal conductivity of h2osfc
    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      zh2osfc = 1.0D-3*(0.5D0*h2osfc(c)) !convert to [m] from [mm]
      tk_h2osfc(c)= tkwat*thk(c,1)*(z(c,1)+zh2osfc) / &
                    (tkwat*z(c,1)+thk(c,1)*zh2osfc)
    end do

    ! Soil heat capacity, from de Vires (1963)
    ! Urban values are from Masson et al. 2002, Evaluation of
    ! the Town Energy Balance (TEB) scheme with direct measurements
    ! from dry districts in two cities, J. Appl. Meteorol., 41, 1011-1026.

    do j = 1, nlevgrnd
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        if ((ctype(c) == icol_sunwall .OR. &
             ctype(c) == icol_shadewall) .and. j <= nlevurb) then
          cv(c,j) = cv_wall(l,j) * dz(c,j)
        else if (ctype(c) == icol_roof .and. j <= nlevurb) then
          cv(c,j) = cv_roof(l,j) * dz(c,j)
        else if (ctype(c) == icol_road_imperv .and. &
                 j >= 1 .and. j <= nlev_improad(l)) then
          cv(c,j) = cv_improad(l,j) * dz(c,j)
        else if (ltype(l) /= istwet .AND. &
                 ltype(l) /= istice .AND. &
                 ctype(c) /= icol_sunwall .AND. &
                 ctype(c) /= icol_shadewall .AND. &
                 ctype(c) /= icol_roof) then
          cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) + &
            (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
        else if (ltype(l) == istwet) then
          cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          if (j > nlevsoi) cv(c,j) = csol(c,j)*dz(c,j)
        else if ( ltype(l) == istice ) then
          cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
        end if
        if (j == 1) then
          if (snl(c)+1 == 1 .AND. h2osno(c) > 0.D0) then
            cv(c,j) = cv(c,j) + cpice*h2osno(c)
          end if
        end if
      end do
    end do

    ! Snow heat capacity

    do j = -nlevsno+1,0
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
          cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
        end if
      end do
    end do
  end subroutine SoilThermProp
  !
  ! Only freezing is considered.
  ! When water freezes, move ice to bottom snow layer.
  !
  subroutine PhaseChangeH2osfc (lbc, ubc, num_nolakec, filter_nolakec, fact, &
                                dhsdT,c_h2osfc,xmf_h2osfc)
    implicit none
    integer(ik4) , intent(in) :: lbc, ubc                             ! column bounds
    integer(ik4) , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: filter_nolakec(ubc-lbc+1)            ! column filter for non-lake points
    real(rk8), intent(inout) :: fact  (lbc:ubc, -nlevsno+1:nlevgrnd) ! temporary
    real(rk8), intent(in) :: dhsdT (lbc:ubc)                      ! temperature derivative of "hs"
    real(rk8), intent(in) :: c_h2osfc(lbc:ubc)                    ! heat capacity of surface water
    real(rk8), intent(out) :: xmf_h2osfc(lbc:ubc)                 ! latent heat of phase change of surface water

    real(rk8), pointer :: qflx_h2osfc_to_ice(:) ! conversion of h2osfc to ice
    real(rk8), pointer :: frac_sno(:)      ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_h2osfc(:)   ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer :: t_h2osfc(:)     ! surface water temperature
    real(rk8), pointer :: t_h2osfc_bef(:)  ! saved surface water temperature
    real(rk8), pointer :: h2osfc(:)        ! surface water (mm)
    real(rk8), pointer :: int_snow(:)      ! integrated snowfall [mm]
    integer(ik4) , pointer :: snl(:)           ! number of snow layers
    real(rk8), pointer :: h2osno(:)        ! snow water (mm H2O)

    real(rk8), pointer :: snow_depth(:)        !snow height (m)

    real(rk8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2) (new)
    real(rk8), pointer :: tssbef(:,:)      !temperature at previous time step [K]
    real(rk8), pointer :: dz(:,:)          !layer thickness (m)

    real(rk8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)

    integer(ik4)  :: c   !do loop index
    integer(ik4)  :: fc  !lake filtered column indices
    real(rk8) :: temp1   !temporary variables [kg/m2]
    real(rk8) :: hm(lbc:ubc) !energy residual [W/m2]
    real(rk8) :: xm(lbc:ubc) !melting or freezing within a time step [kg/m2]
    real(rk8) :: tinc        !t(n+1)-t(n) (K)
    real(rk8) :: rho_avg
    real(rk8) :: z_avg
    real(rk8) :: dcv(lbc:ubc)
    real(rk8) :: t_h2osfc_new
    real(rk8) :: c1
    real(rk8) :: c2
    real(rk8) :: c_h2osfc_ice

    ! Assign local pointers to derived subtypes components (column-level)

    frac_sno     => clm3%g%l%c%cps%frac_sno_eff
    frac_h2osfc  => clm3%g%l%c%cps%frac_h2osfc
    t_h2osfc     => clm3%g%l%c%ces%t_h2osfc
    t_h2osfc_bef => clm3%g%l%c%ces%t_h2osfc_bef
    h2osfc       => clm3%g%l%c%cws%h2osfc
    int_snow     => clm3%g%l%c%cws%int_snow
    qflx_h2osfc_to_ice       => clm3%g%l%c%cwf%qflx_h2osfc_to_ice
    snl          => clm3%g%l%c%cps%snl
    h2osno       => clm3%g%l%c%cws%h2osno
    snow_depth       => clm3%g%l%c%cps%snow_depth
    h2osoi_ice   => clm3%g%l%c%cws%h2osoi_ice
    t_soisno     => clm3%g%l%c%ces%t_soisno
    tssbef       => clm3%g%l%c%ces%tssbef
    dz           => clm3%g%l%c%cps%dz

    ! Initialization

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      xmf_h2osfc(c) = 0.D0
      hm(c) = 0.D0
      xm(c) = 0.D0
      qflx_h2osfc_to_ice(c) = 0.D0
    end do

    ! Freezing identification
    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      ! If liquid exists below melt point, freeze some to ice.
      if ( frac_h2osfc(c) > 0.D0 .AND. t_h2osfc(c) <= tfrz) then
        tinc = t_h2osfc(c)-tfrz
        t_h2osfc(c) = tfrz
        ! energy absorbed beyond freezing temperature
        hm(c) = dhsdT(c)*tinc - tinc*c_h2osfc(c)/dtsrf

        ! mass of water converted from liquid to ice
        xm(c) = hm(c)*dtsrf/hfus
        temp1 = h2osfc(c) - xm(c)

        ! compute change in cv due to additional ice
        dcv(c) = cpice*min(xm(c),h2osfc(c))

        z_avg = frac_sno(c)*snow_depth(c)
        if (z_avg > 0.D0) then
          rho_avg = min(800.D0,h2osno(c)/z_avg)
        else
          rho_avg = 200.D0
        end if
        ! xm < h2osfc
        if ( temp1 >= 0.D0 ) then ! add some frozen water to snow column
          ! add ice to snow column
          h2osno(c) = h2osno(c) + xm(c)
          int_snow(c) = int_snow(c) + xm(c)

          if ( snl(c) < 0 ) h2osoi_ice(c,0) = h2osoi_ice(c,0) + xm(c)

          ! remove ice from h2osfc
          h2osfc(c) = h2osfc(c) - xm(c)

          xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)

          qflx_h2osfc_to_ice(c) = xm(c)/dtsrf

          ! update snow depth
          if (frac_sno(c) > 0 .and. snl(c) < 0) then
            snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
          else
            snow_depth(c)=h2osno(c)/denice
          end if
        ! xm > h2osfc
        else
          !all h2osfc converted to ice, apply residual heat to top soil layer
          rho_avg = (h2osno(c)*rho_avg + h2osfc(c)*denice) / &
            (h2osno(c) + h2osfc(c))
          h2osno(c) = h2osno(c) + h2osfc(c)
          int_snow(c) = int_snow(c) + h2osfc(c)

          qflx_h2osfc_to_ice(c) = h2osfc(c)/dtsrf

          ! excess energy is used to cool ice layer
          if ( snl(c) < 0 ) h2osoi_ice(c,0) = h2osoi_ice(c,0) + h2osfc(c)

          ! compute heat capacity of frozen h2osfc layer
          c_h2osfc_ice = cpice*denice*(1.0D-3*h2osfc(c)) !h2osfc in [m]

          ! cool frozen h2osfc layer with extra heat
          t_h2osfc_new = t_h2osfc(c) - &
            temp1*hfus/(dtsrf*dhsdT(c) - c_h2osfc_ice)

          ! next, determine equilibrium temperature of combined ice/snow layer
          xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)
          if ( snl(c) == 0 ) then
            t_soisno(c,0) = t_h2osfc_new
          else if ( snl(c) == -1 ) then
            c1 = frac_sno(c)/fact(c,0) - dhsdT(c)*dtsrf
            if ( frac_h2osfc(c) /= 0.0D0 ) then
              c2 = frac_h2osfc(c)*(c_h2osfc_ice/dtsrf)
            else
              c2 = 0.0D0
            end if
            ! account for the change in t_soisno(c,0) via xmf_h2osfc(c)
            xmf_h2osfc(c) = xmf_h2osfc(c) + frac_sno(c)*t_soisno(c,0)/fact(c,0)
            t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new)/(c1 + c2)
            xmf_h2osfc(c) = xmf_h2osfc(c) - frac_sno(c)*t_soisno(c,0)/fact(c,0)
          else
            c1 = frac_sno(c)/fact(c,0)
            if ( frac_h2osfc(c) /= 0.0D0 ) then
              c2 = frac_h2osfc(c)*(c_h2osfc_ice/dtsrf)
            else
              c2 = 0.0D0
            end if
            xmf_h2osfc(c) = xmf_h2osfc(c) + c1*t_soisno(c,0)
            t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new)/(c1 + c2)
            xmf_h2osfc(c) = xmf_h2osfc(c) - c1*t_soisno(c,0)
          end if

          ! set h2osfc to zero (all liquid converted to ice)
          h2osfc(c) = 0.D0

          ! update snow depth
          if ( frac_sno(c) > 0 .and. snl(c) < 0 ) then
            snow_depth(c) = h2osno(c)/(rho_avg*frac_sno(c))
          else
            snow_depth(c) = h2osno(c)/denice
          end if
        end if
      end if
    end do
  end subroutine PhaseChangeH2osfc
  !
  ! Calculation of the phase change within snow and soil layers:
  ! (1) Check the conditions for which the phase change may take place,
  !     i.e., the layer temperature is great than the freezing point
  !     and the ice mass is not equal to zero (i.e. melting),
  !     or the layer temperature is less than the freezing point
  !     and the liquid water mass is greater than the allowable supercooled
  !     liquid water calculated from freezing point depression (i.e. freezing).
  ! (2) Assess the rate of phase change from the energy excess (or deficit)
  !     after setting the layer temperature to freezing point.
  ! (3) Re-adjust the ice and liquid mass, and the layer temperature
  !
  subroutine Phasechange_beta (lbc, ubc, num_nolakec, filter_nolakec, fact, &
                          dhsdT, xmf)
    implicit none
    integer(ik4) , intent(in) :: lbc, ubc                             ! column bounds
    integer(ik4) , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: filter_nolakec(ubc-lbc+1)            ! column filter for non-lake points
    real(rk8), intent(in) :: fact  (lbc:ubc, -nlevsno+1:nlevgrnd) ! temporary
    real(rk8), intent(in) :: dhsdT (lbc:ubc)                      ! temperature derivative of "hs"
    real(rk8), intent(out):: xmf   (lbc:ubc)                      ! total latent heat of phase change

    real(rk8), pointer :: qflx_snow_melt(:)! net snow melt
    real(rk8), pointer :: frac_sno_eff(:)  ! eff. fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno(:)      ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_h2osfc(:)   ! fraction of ground covered by surface water (0 to 1)
    integer(ik4) , pointer :: ctype(:)         !column type
    integer(ik4) , pointer :: ltype(:)         ! landunit type
    integer(ik4) , pointer :: clandunit(:)     ! column's landunit
    integer(ik4) , pointer :: snl(:)           ! number of snow layers
    real(rk8), pointer :: h2osno(:)        ! snow water (mm H2O)

    real(rk8), pointer :: snow_depth(:)        ! snow height (m)

    real(rk8), pointer :: qflx_snomelt(:)  ! snow melt (mm H2O /s)
    real(rk8), pointer :: eflx_snomelt(:)  !snow melt heat flux (W/m**2)
    real(rk8), pointer :: eflx_snomelt_u(:)!urban snow melt heat flux (W/m**2)
    real(rk8), pointer :: eflx_snomelt_r(:)!rural snow melt heat flux (W/m**2)
    real(rk8), pointer :: qflx_snofrz_lyr(:,:)  !snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(rk8), pointer :: qflx_snofrz_col(:)  !column-integrated snow freezing rate (positive definite) [kg m-2 s-1]

    real(rk8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2) (new)
    real(rk8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2) (new)
    real(rk8), pointer :: tssbef(:,:)      !temperature at previous time step [K]
    real(rk8), pointer :: sucsat(:,:)      !minimum soil suction (mm)
    real(rk8), pointer :: watsat(:,:)      !volumetric soil water at saturation (porosity)
    real(rk8), pointer :: bsw(:,:)         !Clapp and Hornberger "b"
    real(rk8), pointer :: dz(:,:)          !layer thickness (m)

    real(rk8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)

    integer(ik4), pointer :: imelt(:,:)        !flag for melting (=1), freezing (=2), Not=0 (new)
!   real(rk8), pointer :: eflx_snomelt_u(:)!urban snow melt heat flux (W/m**2)

    integer(ik4)  :: j,c,l     !do loop index
    integer(ik4)  :: fc        !lake filtered column indices
    real(rk8) :: heatr         !energy residual or loss after melting or freezing
    real(rk8) :: temp1                              !temporary variables [kg/m2]
    real(rk8) :: hm(lbc:ubc,-nlevsno+1:nlevgrnd)    !energy residual [W/m2]
    real(rk8) :: xm(lbc:ubc,-nlevsno+1:nlevgrnd)    !melting or freezing within a time step [kg/m2]
    real(rk8) :: wmass0(lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of ice and liquid (kg/m2)
    real(rk8) :: wice0 (lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of ice (kg/m2)
    real(rk8) :: wliq0 (lbc:ubc,-nlevsno+1:nlevgrnd)!initial mass of liquid (kg/m2)
    real(rk8) :: supercool(lbc:ubc,nlevgrnd)        !supercooled water in soil (kg/m2)
    real(rk8) :: propor                             !proportionality constant (-)
    real(rk8) :: tinc(lbc:ubc,-nlevsno+1:nlevgrnd)  !t(n+1)-t(n) (K)
    real(rk8) :: smp                                !frozen water potential (mm)

    ! Assign local pointers to derived subtypes components (column-level)

    qflx_snow_melt => clm3%g%l%c%cwf%qflx_snow_melt
    frac_sno_eff => clm3%g%l%c%cps%frac_sno_eff
    frac_sno     => clm3%g%l%c%cps%frac_sno
    frac_h2osfc  => clm3%g%l%c%cps%frac_h2osfc
    clandunit    => clm3%g%l%c%landunit
    ltype        => clm3%g%l%itype
    ctype        => clm3%g%l%c%itype
    snl          => clm3%g%l%c%cps%snl
    h2osno       => clm3%g%l%c%cws%h2osno
    snow_depth   => clm3%g%l%c%cps%snow_depth
    qflx_snomelt => clm3%g%l%c%cwf%qflx_snomelt
    eflx_snomelt => clm3%g%l%c%cef%eflx_snomelt
    eflx_snomelt_u => clm3%g%l%c%cef%eflx_snomelt_u
    eflx_snomelt_r => clm3%g%l%c%cef%eflx_snomelt_r
    h2osoi_liq   => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice   => clm3%g%l%c%cws%h2osoi_ice
    imelt        => clm3%g%l%c%cps%imelt
    t_soisno     => clm3%g%l%c%ces%t_soisno
    tssbef       => clm3%g%l%c%ces%tssbef
    bsw          => clm3%g%l%c%cps%bsw
    sucsat       => clm3%g%l%c%cps%sucsat
    watsat       => clm3%g%l%c%cps%watsat
    dz           => clm3%g%l%c%cps%dz
    qflx_snofrz_lyr => clm3%g%l%c%cwf%qflx_snofrz_lyr
    qflx_snofrz_col => clm3%g%l%c%cwf%qflx_snofrz_col

    ! Initialization

    do fc = 1 , num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)

      qflx_snomelt(c) = 0.D0
      xmf(c) = 0.D0
      qflx_snofrz_lyr(c,-nlevsno+1:0) = 0.D0
      qflx_snofrz_col(c) = 0.D0
      qflx_snow_melt(c) = 0.D0
    end do

    do j = -nlevsno+1 , nlevgrnd       ! all layers
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        if ( j >= snl(c)+1 ) then
          ! Initialization
          imelt(c,j) = 0
          hm(c,j) = 0.D0
          xm(c,j) = 0.D0
          wice0(c,j) = h2osoi_ice(c,j)
          wliq0(c,j) = h2osoi_liq(c,j)
          wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
        end if   ! end of snow layer if-block
      end do   ! end of column-loop
    end do   ! end of level-loop

    ! snow layers
    do j = -nlevsno+1 , 0
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        if ( j >= snl(c)+1 ) then
          ! Melting identification
          ! If ice exists above melt point, melt some to liquid.
          if (h2osoi_ice(c,j) > 0.D0 .AND. t_soisno(c,j) > tfrz) then
            imelt(c,j) = 1
            ! tinc(c,j) = t_soisno(c,j) - tfrz
            tinc(c,j) = tfrz - t_soisno(c,j)
            t_soisno(c,j) = tfrz
          end if

          ! Freezing identification
          ! If liquid exists below melt point, freeze some to ice.
          if ( h2osoi_liq(c,j) > 0.D0 .AND. t_soisno(c,j) < tfrz ) then
            imelt(c,j) = 2
            ! tinc(c,j) = t_soisno(c,j) - tfrz
            tinc(c,j) = tfrz - t_soisno(c,j)
            t_soisno(c,j) = tfrz
          end if
        end if   ! end of snow layer if-block
      end do   ! end of column-loop
    end do   ! end of level-loop

    ! soil layers
    do j = 1 , nlevgrnd
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)
        l = clandunit(c)
        supercool(c,j) = 0.0D0
        ! add in urban condition if-block
        if ((ctype(c) /= icol_sunwall .and. &
             ctype(c) /= icol_shadewall .and. &
             ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then
          if ( h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz ) then
            imelt(c,j) = 1
            ! tinc(c,j) = t_soisno(c,j) - tfrz
            tinc(c,j) = tfrz - t_soisno(c,j)
            t_soisno(c,j) = tfrz
          end if

          ! from Zhao (1997) and Koren (1999)
          supercool(c,j) = 0.0D0
          if (ltype(l) == istsoil .or. &
              ltype(l) == istcrop .or. &
              ctype(c) == icol_road_perv) then
            if ( t_soisno(c,j) < tfrz ) then
              smp = hfus*(tfrz-t_soisno(c,j)) / &
                   (grav*t_soisno(c,j)) * 1000.D0  !(mm)
              supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1.D0/bsw(c,j))
              supercool(c,j) = supercool(c,j)*dz(c,j)*1000.D0       ! (mm)
            end if
          end if

          if ( h2osoi_liq(c,j) > supercool(c,j) .AND. &
               t_soisno(c,j) < tfrz) then
            imelt(c,j) = 2
            ! tinc(c,j) = t_soisno(c,j) - tfrz
            tinc(c,j) = tfrz - t_soisno(c,j)
            t_soisno(c,j) = tfrz
          end if

          ! If snow exists, but its thickness is less than
          ! the critical value (0.01 m)
          if ( snl(c)+1 == 1 .AND. h2osno(c) > 0.D0 .AND. j == 1 ) then
            if ( t_soisno(c,j) > tfrz ) then
              imelt(c,j) = 1
              ! tincc,j) = t_soisno(c,j) - tfrz
              tinc(c,j) = tfrz - t_soisno(c,j)
              t_soisno(c,j) = tfrz
            end if
          end if
        end if
      end do
    end do

    do j = -nlevsno+1 , nlevgrnd       ! all layers
      do fc = 1 , num_nolakec
        c = filter_nolakec(fc)

        if ((ctype(c) /= icol_sunwall .and. &
             ctype(c) /= icol_shadewall .and. &
             ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then
          if ( j >= snl(c)+1 ) then

            ! Calculate the energy surplus and loss for melting and freezing
            if ( imelt(c,j) > 0 ) then
              ! added unique cases for this calculation,
              ! to account for absorbed solar radiation in each layer
              !======================================================
              if ( j == snl(c)+1 ) then ! top layer
                hm(c,j) = dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)

                if ( j==1 .and. frac_h2osfc(c) /= 0.0D0 ) then
                  hm(c,j) = hm(c,j) - frac_h2osfc(c)*(dhsdT(c)*tinc(c,j))
                end if
              else if (j == 1) then
                hm(c,j) = (1.0D0 - frac_sno_eff(c) - frac_h2osfc(c)) * &
                           dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)
              else ! non-interfacial snow/soil layers
                hm(c,j) = - tinc(c,j)/fact(c,j)
              end if
            end if

            ! These two errors were checked carefully (Y. Dai).
            ! They result from the computed error of "Tridiagonal-Matrix"
            ! in subroutine "thermal".
            if (imelt(c,j) == 1 .AND. hm(c,j) < 0.D0) then
              hm(c,j) = 0.D0
              imelt(c,j) = 0
            end if
            if (imelt(c,j) == 2 .AND. hm(c,j) > 0.D0) then
              hm(c,j) = 0.D0
              imelt(c,j) = 0
            end if

            ! The rate of melting and freezing

            if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0.D0) then
              xm(c,j) = hm(c,j)*dtsrf/hfus ! kg/m2

              ! If snow exists, but its thickness is less than
              ! the critical value (1 cm). Note: more work is needed
              ! to determine how to tune the snow depth for this case
              if (j == 1) then
                if (snl(c)+1 == 1 .AND. h2osno(c) > 0.D0 .AND. &
                    xm(c,j) > 0.D0) then
                  temp1 = h2osno(c)   ! kg/m2
                  h2osno(c) = max(0.D0,temp1-xm(c,j))
                  propor = h2osno(c)/temp1
                  snow_depth(c) = propor * snow_depth(c)
                  heatr = hm(c,j) - hfus*(temp1-h2osno(c))/dtsrf   ! W/m2
                  if (heatr > 0.D0) then
                    xm(c,j) = heatr*dtsrf/hfus  ! kg/m2
                    hm(c,j) = heatr             ! W/m2
                  else
                    xm(c,j) = 0.D0
                    hm(c,j) = 0.D0
                  end if
                  ! kg/(m2 s)
                  qflx_snomelt(c) = max(0.D0,(temp1-h2osno(c)))/dtsrf
                  xmf(c) = hfus*qflx_snomelt(c)
                  qflx_snow_melt(c) = qflx_snomelt(c)
                end if
              end if

              heatr = 0.D0
              if (xm(c,j) > 0.D0) then
                h2osoi_ice(c,j) = max(0.D0, wice0(c,j)-xm(c,j))
                heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtsrf
              else if (xm(c,j) < 0.D0) then
                if (j <= 0) then
                  h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                else
                  if (wmass0(c,j) < supercool(c,j)) then
                    h2osoi_ice(c,j) = 0.D0
                  else
                    h2osoi_ice(c,j) = min(wmass0(c,j) - &
                               supercool(c,j),wice0(c,j)-xm(c,j))
                  end if
                end if
                heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtsrf
              end if

              h2osoi_liq(c,j) = max(0.D0,wmass0(c,j)-h2osoi_ice(c,j))

              if ( abs(heatr) > 0.D0 ) then
                if ( j == snl(c)+1 ) then
                  if ( j == 1 ) then
                    t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr / &
                            (1.D0-(1.0D0 - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                  else
                    t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr / &
                            (1.D0-fact(c,j)*dhsdT(c))
                  end if
                else if (j == 1) then
                  t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr / &
                            (1.D0-(1.0D0 - frac_sno_eff(c) - &
                            frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                else
                  t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                end if

                if ( j <= 0 ) then    ! snow
                  if (h2osoi_liq(c,j)*h2osoi_ice(c,j) > 0.D0) then
                    t_soisno(c,j) = tfrz
                  end if
                end if
              end if  ! end of heatr > 0 if-block

              if ( j >= 1 ) then
                xmf(c) = xmf(c) + hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtsrf
              else
                xmf(c) = xmf(c) + &
                  frac_sno_eff(c)*hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtsrf
              end if

              if (imelt(c,j) == 1 .AND. j < 1) then
                qflx_snomelt(c) = qflx_snomelt(c) + &
                  max(0.D0,(wice0(c,j)-h2osoi_ice(c,j)))/dtsrf
              end if

              ! layer freezing mass flux (positive):
              if (imelt(c,j) == 2 .AND. j < 1) then
                qflx_snofrz_lyr(c,j) = max(0.D0, &
                  (h2osoi_ice(c,j)-wice0(c,j)))/dtsrf
              end if
            end if
          end if   ! end of snow layer if-block
        end if
      end do   ! end of column-loop
    end do   ! end of level-loop

    ! Needed for history file output

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      eflx_snomelt(c) = qflx_snomelt(c) * hfus
      l = clandunit(c)
      if (ltype(l) == isturb) then
        eflx_snomelt_u(c) = eflx_snomelt(c)
      else if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
        eflx_snomelt_r(c) = eflx_snomelt(c)
      end if
    end do
    do j = -nlevsno+1,0
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
      end do
    end do
  end subroutine Phasechange_beta

end module mod_clm_soiltemperature
