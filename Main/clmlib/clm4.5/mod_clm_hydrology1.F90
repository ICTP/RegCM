module mod_clm_hydrology1
  !
  ! Calculation of
  ! (1) water storage of intercepted precipitation
  ! (2) direct throughfall and canopy drainage of precipitation
  ! (3) the fraction of foliage covered by water and the fraction
  !     of foliage that is dry and transpiring.
  ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_stdio

  implicit none

  private

  save

  public :: Hydrology1_readnl ! Read namelist
  public :: Hydrology1        ! Run

  ! use old fsno parameterization (N&Y07)
  integer(ik4) , public :: oldfflag = 0

  contains
  !
  ! Read the namelist for Hydrology1
  !
  subroutine Hydrology1_readnl( NLFilename )
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    integer(ik4) :: ierr          ! error code
    integer(ik4) :: unitn         ! unit for namelist file
    character(len=32) :: subname = 'Hydrology1_readnl'  ! subroutine name

    namelist / clm_hydrology1_inparm / oldfflag

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( myid == iocpu ) then
      unitn = file_getUnit()
      write(stdout,*) 'Read in clm_hydrology1_inparm namelist'
      open(unitn,file=NLFilename,status='old',action='read',iostat=ierr)
      if (ierr /= 0) then
        call fatal(__FILE__,__LINE__, &
           subname // ':: ERROR open namelist file '//NLFilename)
      end if
      read(unitn, clm_hydrology1_inparm, iostat=ierr)
      if (ierr /= 0) then
        call fatal(__FILE__,__LINE__, &
             subname // ':: ERROR reading clm_hydrology1_inparm namelist')
      end if
      call file_freeUnit( unitn )
    end if
    ! Broadcast namelist variables read in
    call bcast(oldfflag)
  end subroutine Hydrology1_readnl
  !
  ! Calculation of
  ! (1) water storage of intercepted precipitation
  ! (2) direct throughfall and canopy drainage of precipitation
  ! (3) the fraction of foliage covered by water and the fraction
  !     of foliage that is dry and transpiring.
  ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
  ! Note:  The evaporation loss is taken off after the calculation of leaf
  ! temperature in the subroutine clm\_leaftem.f90, not in this subroutine.
  !
  subroutine Hydrology1(lbc, ubc, lbp, ubp, num_nolakec, filter_nolakec, &
                        num_nolakep, filter_nolakep)
    use mod_clm_type
    use mod_clm_atmlnd , only : clm_a2l
    use mod_clm_varcon , only : tfrz, istice, istwet, istsoil, isturb,    &
           istcrop, icol_roof, icol_sunwall, icol_shadewall, hfus,denice, &
           zlnd,rpi,spval
    use mod_clm_varctl , only : subgridflag
    use mod_clm_varpar , only : nlevsoi,nlevsno
    use mod_clm_h2osfc , only : FracH2oSfc
    use mod_clm_fracwet , only : FracWet
    use mod_clm_subgridave , only : p2c
    use mod_clm_snicar , only : snw_rds_min
    implicit none
    integer(ik4), intent(in) :: lbp, ubp    ! pft bounds
    integer(ik4), intent(in) :: lbc, ubc    ! column bounds
    ! number of column non-lake points in column filter
    integer(ik4), intent(in) :: num_nolakec
    ! column filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakec(ubc-lbc+1)
    ! number of pft non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakep(ubp-lbp+1)

    real(rkx), pointer :: swe_old(:,:)   ! snow water before update
    ! eff. fraction of ground covered by snow (0 to 1)
    real(rkx), pointer :: frac_sno_eff(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rkx), pointer :: frac_sno(:)
    real(rkx), pointer :: h2osfc(:)         ! surface water (mm)
    ! fraction of ground covered by surface water (0 to 1)
    real(rkx), pointer :: frac_h2osfc(:)
    !snow falling on surface water (mm/s)
    real(rkx), pointer :: qflx_snow_h2osfc(:)
    ! integrated snowfall [mm]
    real(rkx), pointer :: int_snow(:)
    ! gridcell flux of flood water from RTM
    real(rkx), pointer :: qflx_floodg(:)
    real(rkx), pointer :: qflx_floodc(:) ! column flux of flood water from RTM
    real(rkx), pointer :: qflx_snow_melt(:) ! snow melt from previous time step
    real(rkx), pointer :: n_melt(:)         ! SCA shape parameter
    integer(ik4) , pointer :: cgridcell(:)  ! columns's gridcell
    integer(ik4) , pointer :: clandunit(:)  ! columns's landunit
    integer(ik4) , pointer :: pgridcell(:)  ! pft's gridcell
    integer(ik4) , pointer :: plandunit(:)  ! pft's landunit
    integer(ik4) , pointer :: pcolumn(:)    ! pft's column
    integer(ik4) , pointer :: npfts(:)      ! number of pfts in column
    integer(ik4) , pointer :: pfti(:)       ! column's beginning pft index
    integer(ik4) , pointer :: ltype(:)      ! landunit type
    integer(ik4) , pointer :: ctype(:)      ! column type
    real(rkx), pointer :: forc_rain(:)      ! rain rate [mm/s]
    real(rkx), pointer :: forc_snow(:)      ! snow rate [mm/s]
    real(rkx), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    logical , pointer :: do_capsnow(:)      ! true => do snow capping
    real(rkx), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
    real(rkx), pointer :: dewmx(:)          ! Maximum allowed dew [mm]
    ! fraction of veg not covered by snow (0/1 now) [-]
    integer(ik4) , pointer :: frac_veg_nosno(:)
    ! one-sided leaf area index with burying by snow
    real(rkx), pointer :: elai(:)
    ! one-sided stem area index with burying by snow
    real(rkx), pointer :: esai(:)
    ! canopy water mass balance term (column)
    real(rkx), pointer :: h2ocan_loss(:)
    ! current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]
    real(rkx), pointer :: irrig_rate(:)

    integer(ik4) , pointer :: snl(:)        ! number of snow layers
    real(rkx), pointer :: snow_depth(:)     ! snow height (m)
    real(rkx), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(rkx), pointer :: h2ocan(:)         ! total canopy water (mm H2O)
    real(rkx), pointer :: qflx_irrig(:)     ! irrigation amount (mm/s)
    ! number of time steps for which we still need to irrigate today
    integer(ik4), pointer  :: n_irrig_steps_left(:)

    ! interception of precipitation [mm/s]
    real(rkx), pointer :: qflx_prec_intr(:)
    ! water onto ground including canopy runoff [kg/(m2 s)]
    real(rkx), pointer :: qflx_prec_grnd(:)
    ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(rkx), pointer :: qflx_snwcp_liq(:)
    ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(rkx), pointer :: qflx_snwcp_ice(:)
    ! snow on ground after interception (mm H2O/s) [+]
    real(rkx), pointer :: qflx_snow_grnd_pft(:)
    ! snow on ground after interception (mm H2O/s) [+]
    real(rkx), pointer :: qflx_snow_grnd_col(:)
    ! rain on ground after interception (mm H2O/s) [+]
    real(rkx), pointer :: qflx_rain_grnd(:)
    ! fraction of canopy that is wet (0 to 1)
    real(rkx), pointer :: fwet(:)
    ! fraction of foliage that is green and dry [-] (new)
    real(rkx), pointer :: fdry(:)
    ! interface level below a "z" level (m)
    real(rkx), pointer :: zi(:,:)
    real(rkx), pointer :: dz(:,:)   ! layer depth (m)
    real(rkx), pointer :: z(:,:)    ! layer thickness (m)
    real(rkx), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)
    real(rkx), pointer :: h2osoi_ice(:,:)    ! ice lens (kg/m2)
    real(rkx), pointer :: h2osoi_liq(:,:)    ! liquid water (kg/m2)
    ! fraction of ice relative to the tot water
    real(rkx), pointer :: frac_iceold(:,:)
    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(rkx), pointer :: snw_rds(:,:)
    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_bcpho(:,:)
    ! mass of hydrophilic BC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_bcphi(:,:)
    ! total mass of BC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_bctot(:,:)
    ! total column mass of BC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_bc_col(:)
    ! total top-layer mass of BC (col,lyr) [kg]
    real(rkx), pointer :: mss_bc_top(:)
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_ocpho(:,:)
    ! mass of hydrophilic OC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_ocphi(:,:)
    ! total mass of OC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_octot(:,:)
    ! total column mass of OC in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_oc_col(:)
    ! total top-layer mass of OC (col,lyr) [kg]
    real(rkx), pointer :: mss_oc_top(:)
    ! mass of dust species 1 in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst1(:,:)
    ! mass of dust species 2 in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst2(:,:)
    ! mass of dust species 3 in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst3(:,:)
    ! mass of dust species 4 in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst4(:,:)
    ! total mass of dust in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dsttot(:,:)
    ! total column mass of dust in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst_col(:)
    ! total top-layer mass of dust in snow (col,lyr) [kg]
    real(rkx), pointer :: mss_dst_top(:)

    integer(ik4)  :: f     ! filter index
    integer(ik4)  :: p     ! pft index
    integer(ik4)  :: c     ! column index
    integer(ik4)  :: l     ! landunit index
    integer(ik4)  :: g     ! gridcell index
    integer(ik4)  :: newnode ! flag when new snow node is set, (1=yes, 0=no)
    real(rkx) :: h2ocanmx  ! maximum allowed water on canopy [mm]
    real(rkx) :: fpi       ! coefficient of interception
    real(rkx) :: xrun      ! excess water that exceeds the leaf capacity [mm/s]
    ! layer thickness rate change due to precipitation [mm/s]
    real(rkx) :: dz_snowf
    ! bulk density of newly fallen dry snow [kg/m3]
    real(rkx) :: bifall
    real(rkx) :: fracsnow(lbp:ubp)     ! frac of precipitation that is snow
    real(rkx) :: fracrain(lbp:ubp)     ! frac of precipitation that is rain
    ! rate of canopy runoff and snow falling off canopy [mm/s]
    real(rkx) :: qflx_candrip(lbp:ubp)
    real(rkx) :: qflx_through_rain(lbp:ubp)   ! direct rain throughfall [mm/s]
    real(rkx) :: qflx_through_snow(lbp:ubp)   ! direct snow throughfall [mm/s]
    ! snow precipitation incident on ground [mm/s]
    real(rkx) :: qflx_prec_grnd_snow(lbp:ubp)
    ! rain precipitation incident on ground [mm/s]
    real(rkx) :: qflx_prec_grnd_rain(lbp:ubp)
    real(rkx) :: z_avg                        ! grid cell average snow depth
    real(rkx) :: temp_snow_depth,temp_intsnow ! temporary variables
    real(rkx) :: fmelt
    real(rkx) :: smr
    real(rkx) :: fsno_new
    real(rkx) :: accum_factor
    real(rkx) :: newsnow(lbc:ubc)
    real(rkx) :: snowmelt(lbc:ubc)
    integer(ik4)  :: j

    ! Assign local pointers to derived type members (gridcell-level)

    pgridcell          => clm3%g%l%c%p%gridcell
    forc_rain          => clm_a2l%forc_rain
    forc_snow          => clm_a2l%forc_snow

    ! Assign local pointers to derived type members (landunit-level)

    ltype              => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    swe_old            => clm3%g%l%c%cws%swe_old
    frac_sno_eff       => clm3%g%l%c%cps%frac_sno_eff
    frac_sno           => clm3%g%l%c%cps%frac_sno
    frac_h2osfc        => clm3%g%l%c%cps%frac_h2osfc
    h2osfc             => clm3%g%l%c%cws%h2osfc
    qflx_snow_h2osfc   => clm3%g%l%c%cwf%qflx_snow_h2osfc
    int_snow           => clm3%g%l%c%cws%int_snow
    qflx_floodg        => clm_a2l%forc_flood
    qflx_floodc        => clm3%g%l%c%cwf%qflx_floodc
    qflx_snow_melt     => clm3%g%l%c%cwf%qflx_snow_melt
    n_melt             => clm3%g%l%c%cps%n_melt

    cgridcell          => clm3%g%l%c%gridcell
    clandunit          => clm3%g%l%c%landunit
    ctype              => clm3%g%l%c%itype
    pfti               => clm3%g%l%c%pfti
    npfts              => clm3%g%l%c%npfts
    do_capsnow         => clm3%g%l%c%cps%do_capsnow
    forc_t             => clm3%g%l%c%ces%forc_t
    t_grnd             => clm3%g%l%c%ces%t_grnd
    snl                => clm3%g%l%c%cps%snl
    snow_depth             => clm3%g%l%c%cps%snow_depth
    h2osno             => clm3%g%l%c%cws%h2osno
    zi                 => clm3%g%l%c%cps%zi
    dz                 => clm3%g%l%c%cps%dz
    z                  => clm3%g%l%c%cps%z
    frac_iceold        => clm3%g%l%c%cps%frac_iceold
    t_soisno           => clm3%g%l%c%ces%t_soisno
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    qflx_snow_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd
    h2ocan_loss        => clm3%g%l%c%cwf%h2ocan_loss
    snw_rds            => clm3%g%l%c%cps%snw_rds
    mss_bcpho          => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi          => clm3%g%l%c%cps%mss_bcphi
    mss_bctot          => clm3%g%l%c%cps%mss_bctot
    mss_bc_col         => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top         => clm3%g%l%c%cps%mss_bc_top
    mss_ocpho          => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi          => clm3%g%l%c%cps%mss_ocphi
    mss_octot          => clm3%g%l%c%cps%mss_octot
    mss_oc_col         => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top         => clm3%g%l%c%cps%mss_oc_top
    mss_dst1           => clm3%g%l%c%cps%mss_dst1
    mss_dst2           => clm3%g%l%c%cps%mss_dst2
    mss_dst3           => clm3%g%l%c%cps%mss_dst3
    mss_dst4           => clm3%g%l%c%cps%mss_dst4
    mss_dsttot         => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col        => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top        => clm3%g%l%c%cps%mss_dst_top

    ! Assign local pointers to derived type members (pft-level)

    plandunit          => clm3%g%l%c%p%landunit
    pcolumn            => clm3%g%l%c%p%column
    dewmx              => clm3%g%l%c%p%pps%dewmx
    frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno
    elai               => clm3%g%l%c%p%pps%elai
    esai               => clm3%g%l%c%p%pps%esai
    h2ocan             => clm3%g%l%c%p%pws%h2ocan
    qflx_prec_intr     => clm3%g%l%c%p%pwf%qflx_prec_intr
    qflx_prec_grnd     => clm3%g%l%c%p%pwf%qflx_prec_grnd
    qflx_snwcp_liq     => clm3%g%l%c%p%pwf%qflx_snwcp_liq
    qflx_snwcp_ice     => clm3%g%l%c%p%pwf%qflx_snwcp_ice
    qflx_snow_grnd_pft => clm3%g%l%c%p%pwf%qflx_snow_grnd
    qflx_rain_grnd     => clm3%g%l%c%p%pwf%qflx_rain_grnd
    fwet               => clm3%g%l%c%p%pps%fwet
    fdry               => clm3%g%l%c%p%pps%fdry
    irrig_rate         => clm3%g%l%c%cps%irrig_rate
    n_irrig_steps_left => clm3%g%l%c%cps%n_irrig_steps_left
    qflx_irrig         => clm3%g%l%c%cwf%qflx_irrig

    ! Start pft loop

    do f = 1, num_nolakep
      p = filter_nolakep(f)
      g = pgridcell(p)
      l = plandunit(p)
      c = pcolumn(p)

      ! Canopy interception and precipitation onto ground surface
      ! Add precipitation to leaf water

      if ( ltype(l) == istsoil .or. &
           ltype(l) == istwet .or.  &
           ltype(l) == isturb .or.  &
           ltype(l) == istcrop) then
        qflx_candrip(p) = 0._rkx      ! rate of canopy runoff
        qflx_through_snow(p) = 0._rkx ! rain precipitation direct through canopy
        qflx_through_rain(p) = 0._rkx ! snow precipitation direct through canopy
        qflx_prec_intr(p) = 0._rkx    ! total intercepted precipitation
        fracsnow(p) = 0._rkx          ! fraction of input precip that is snow
        fracrain(p) = 0._rkx          ! fraction of input precip that is rain

        if ( ctype(c) /= icol_sunwall .and. &
             ctype(c) /= icol_shadewall ) then
          if ( frac_veg_nosno(p) == 1 .and. &
               (forc_rain(g) + forc_snow(g)) > 0._rkx ) then

            ! determine fraction of input precipitation that is snow and rain

            fracsnow(p) = forc_snow(g)/(forc_snow(g) + forc_rain(g))
            fracrain(p) = forc_rain(g)/(forc_snow(g) + forc_rain(g))

            ! The leaf water capacities for solid and liquid are different,
            ! generally double for snow, but these are of somewhat less
            ! significance for the water budget because of lower evap. rate at
            ! lower temperature.  Hence, it is reasonable to assume that
            ! vegetation storage of solid water is the same as liquid water.
            h2ocanmx = dewmx(p) * (elai(p) + esai(p))

            ! Coefficient of interception
            ! set fraction of potential interception to max 0.25
            fpi = 0.25_rkx*(1._rkx - exp(-0.5_rkx*(elai(p) + esai(p))))

            ! Direct throughfall
            qflx_through_snow(p) = forc_snow(g) * (1._rkx-fpi)
            qflx_through_rain(p) = forc_rain(g) * (1._rkx-fpi)

            ! Intercepted precipitation [mm/s]
            qflx_prec_intr(p) = (forc_snow(g) + forc_rain(g)) * fpi

            ! Water storage of intercepted precipitation and dew
            h2ocan(p) = max(0._rkx, h2ocan(p) + dtsrf*qflx_prec_intr(p))
            if ( h2ocan(p) < 1.0e-10_rkx ) h2ocan(p) = 0.0_rkx

            ! Initialize rate of canopy runoff and snow falling off canopy
            qflx_candrip(p) = 0._rkx

            ! Excess water that exceeds the leaf capacity
            xrun = (h2ocan(p) - h2ocanmx)/dtsrf

            ! Test on maximum dew on leaf
            ! Note if xrun > 0 then h2ocan must be at least h2ocanmx
            if ( xrun > 0._rkx ) then
              qflx_candrip(p) = xrun
              h2ocan(p) = h2ocanmx
            end if
          end if
        end if
      else if ( ltype(l) == istice ) then
        h2ocan(p)            = 0._rkx
        qflx_candrip(p)      = 0._rkx
        qflx_through_snow(p) = 0._rkx
        qflx_through_rain(p) = 0._rkx
        qflx_prec_intr(p)    = 0._rkx
        fracsnow(p)          = 0._rkx
        fracrain(p)          = 0._rkx
      end if

      ! Precipitation onto ground (kg/(m2 s))
      ! PET, 1/18/2005: Added new terms for mass balance correction
      ! due to dynamic pft weight shifting (column-level h2ocan_loss)
      ! Because the fractionation between rain and snow is indeterminate if
      ! rain + snow = 0, I am adding this very small flux only to the rain
      ! components.

      if ( ctype(c) /= icol_sunwall .and. &
           ctype(c) /= icol_shadewall ) then
        if ( frac_veg_nosno(p) == 0 ) then
          qflx_prec_grnd_snow(p) = forc_snow(g)
          qflx_prec_grnd_rain(p) = forc_rain(g) + h2ocan_loss(c)
        else
          qflx_prec_grnd_snow(p) = qflx_through_snow(p) + &
                  (qflx_candrip(p) * fracsnow(p))
          qflx_prec_grnd_rain(p) = qflx_through_rain(p) + &
                  (qflx_candrip(p) * fracrain(p)) + h2ocan_loss(c)
        end if
      ! Urban sunwall and shadewall have no intercepted precipitation
      else
        qflx_prec_grnd_snow(p) = 0.0_rkx
        qflx_prec_grnd_rain(p) = 0.0_rkx
      end if

      ! Determine whether we're irrigating here; set qflx_irrig appropriately
      if ( n_irrig_steps_left(c) > 0 ) then
        qflx_irrig(c)         = irrig_rate(c)
        n_irrig_steps_left(c) = n_irrig_steps_left(c) - 1
      else
        qflx_irrig(c) = 0._rkx
      end if

      ! Add irrigation water directly onto ground (bypassing canopy
      ! interception)
      ! Note that it's still possible that (some of) this irrigation water
      ! will runoff (as runoff is computed later)
      qflx_prec_grnd_rain(p) = qflx_prec_grnd_rain(p) + qflx_irrig(c)

      ! Done irrigation

      qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

      if ( do_capsnow(c) ) then
        qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
        qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)

        qflx_snow_grnd_pft(p) = 0._rkx
        qflx_rain_grnd(p) = 0._rkx
      else
        qflx_snwcp_liq(p) = 0._rkx
        qflx_snwcp_ice(p) = 0._rkx
        ! ice onto ground (mm/s)
        qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)
        ! liquid water onto ground (mm/s)
        qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)
      end if
    end do ! (end pft loop)

    ! Determine the fraction of foliage covered by water and the
    ! fraction of foliage that is dry and transpiring.

    call FracWet(num_nolakep, filter_nolakep)

    ! Update column level state variables for snow.

    call p2c(num_nolakec,filter_nolakec,qflx_snow_grnd_pft,qflx_snow_grnd_col)

    ! apply gridcell flood water flux to non-lake columns
    do f = 1 , num_nolakec
      c = filter_nolakec(f)
      g = cgridcell(c)
      if ( ctype(c) /= icol_sunwall .and. &
           ctype(c) /= icol_shadewall ) then
        qflx_floodc(c) = qflx_floodg(g)
      else
        qflx_floodc(c) = 0._rkx
      end if
    end do

    ! Determine snow height and snow water

    do f = 1 , num_nolakec
      c = filter_nolakec(f)
      l = clandunit(c)
      g = cgridcell(c)

      ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
      ! U.S.Department of Agriculture Forest Service, Project F,
      ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

      qflx_snow_h2osfc(c) = 0._rkx
      ! set temporary variables prior to updating
      temp_snow_depth=snow_depth(c)
      ! save initial snow content
      do j = -nlevsno+1 , snl(c)
        swe_old(c,j) = 0.0_rkx
      end do
      do j = snl(c)+1 , 0
        swe_old(c,j) = h2osoi_liq(c,j)+h2osoi_ice(c,j)
      end do

      if ( do_capsnow(c) ) then
        dz_snowf = 0._rkx
        newsnow(c) = (1._rkx - frac_h2osfc(c)) * qflx_snow_grnd_col(c) * dtsrf
        frac_sno(c) = 1._rkx
        int_snow(c) = 5.e2_rkx
      else
        if ( forc_t(c) > tfrz + 2._rkx ) then
          bifall = 50._rkx + 1.7_rkx*(17.0_rkx)**1.5_rkx
        else if ( forc_t(c) > tfrz - 15._rkx ) then
          bifall = 50._rkx + 1.7_rkx*(forc_t(c) - tfrz + 15._rkx)**1.5_rkx
        else
          bifall = 50._rkx
        end if

        ! newsnow is all snow that doesn't fall on h2osfc
        newsnow(c) = (1._rkx - frac_h2osfc(c)) * qflx_snow_grnd_col(c) * dtsrf

        ! update int_snow
        !h2osno could be larger due to frost
        int_snow(c) = max(int_snow(c),h2osno(c))

        ! snowmelt from previous time step * dtsrf
        snowmelt(c) = qflx_snow_melt(c) * dtsrf

        ! set shape factor for accumulation of snow
        accum_factor = 0.1_rkx

        if ( h2osno(c) > 0.0 ) then

          !==================  FSCA PARAMETERIZATIONS  ===================
          ! fsca parameterization based on *changes* in swe
          ! first compute change from melt during previous time step
          if ( snowmelt(c) > 0._rkx ) then
            smr = min(1._rkx,(h2osno(c))/(int_snow(c)))
            frac_sno(c) = 1.0_rkx - &
                    (acos(min(1._rkx,(2.*smr - 1._rkx)))/rpi)**(n_melt(c))
          end if

          ! update fsca by new snow event, add to previous fsca
          if ( newsnow(c) > 0._rkx ) then
            fsno_new = 1._rkx - &
                    (1._rkx - tanh(accum_factor*newsnow(c)))*(1._rkx - frac_sno(c))
            frac_sno(c) = fsno_new

            ! reset int_snow after accumulation events
            temp_intsnow = (h2osno(c) + newsnow(c)) / &
                   (0.5_rkx*(cos(rpi*(1._rkx-max(frac_sno(c),1e-6_rkx))** &
                   (1.0_rkx/n_melt(c)))+1._rkx))
            int_snow(c) = min(1.e8_rkx,temp_intsnow)
          end if

          !====================================================================

          ! for subgrid fluxes
          if ( subgridflag == 1 .and. ltype(l) /= isturb ) then
            if ( frac_sno(c) > 0._rkx )then
              snow_depth(c)=snow_depth(c) + newsnow(c)/(bifall * frac_sno(c))
            else
              snow_depth(c) = 0._rkx
            end if
          else
            ! for uniform snow cover
            snow_depth(c) = snow_depth(c)+newsnow(c)/bifall
          end if

          ! use original fsca formulation (n&y 07)
          if ( oldfflag == 1 ) then
            ! snow cover fraction in Niu et al. 2007
            if ( snow_depth(c) > 0.0_rkx )  then
              frac_sno(c) = tanh(snow_depth(c)/(2.5_rkx*zlnd* &
                  (min(800._rkx,(h2osno(c) + newsnow(c)) / &
                                   snow_depth(c))/100._rkx)**1._rkx) )
            end if
            if ( h2osno(c) < 1.0_rkx ) then
              frac_sno(c) = min(frac_sno(c),h2osno(c))
            end if
          end if
        else ! h2osno == 0
          ! initialize frac_sno and snow_depth when no snow present initially
          if ( newsnow(c) > 0._rkx ) then
            z_avg = newsnow(c)/bifall
            fmelt = newsnow(c)
            frac_sno(c) = tanh(accum_factor*newsnow(c))

            ! make int_snow consistent w/ new fsno, h2osno
            int_snow(c) = 0.0_rkx !reset prior to adding newsnow below
            temp_intsnow= (h2osno(c) + newsnow(c)) / &
                     (0.5_rkx*(cos(rpi*(1._rkx-max(frac_sno(c),1e-6_rkx))** &
                     (1.0_rkx/n_melt(c)))+1._rkx))
            int_snow(c) = min(1.e8_rkx,temp_intsnow)

            ! update snow_depth and h2osno to be consistent with frac_sno, z_avg
            if ( subgridflag == 1 .and. ltype(l) /= isturb ) then
              snow_depth(c) = z_avg/frac_sno(c)
            else
              snow_depth(c) = newsnow(c)/bifall
            end if
            ! use n&y07 formulation
            if ( oldfflag == 1 ) then
              ! snow cover fraction in Niu et al. 2007
              if ( snow_depth(c) > 0.0_rkx ) then
                frac_sno(c) = tanh(snow_depth(c)/(2.5_rkx*zlnd* &
                  (min(800._rkx,newsnow(c)/snow_depth(c))/100._rkx)**1._rkx) )
              end if
            end if
          else
            z_avg = 0._rkx
            snow_depth(c) = 0._rkx
            frac_sno(c) = 0._rkx
          end if
        end if ! end of h2osno > 0

        ! snow directly falling on surface water melts, increases h2osfc
        qflx_snow_h2osfc(c) = frac_h2osfc(c)*qflx_snow_grnd_col(c)

        ! update h2osno for new snow
        h2osno(c) = h2osno(c) + newsnow(c)
        int_snow(c) = int_snow(c) + newsnow(c)

        ! update change in snow depth
        dz_snowf = (snow_depth(c) - temp_snow_depth) / dtsrf

      end if !end of do_capsnow construct

      ! set frac_sno_eff variable
      if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
        if ( subgridflag == 1 ) then
          frac_sno_eff(c) = frac_sno(c)
        else
          frac_sno_eff(c) = 1._rkx
        end if
      else
        frac_sno_eff(c) = 1._rkx
      end if

      if ( ltype(l) == istwet .and. t_grnd(c) > tfrz ) then
        h2osno(c) = 0._rkx
        snow_depth(c) = 0._rkx
      end if

      ! When the snow accumulation exceeds 10 mm, initialize snow layer
      ! Currently, the water temperature for the precipitation is simply set
      ! as the surface air temperature

      newnode = 0    ! flag for when snow node will be initialized
      if ( snl(c) == 0 .and. &
           qflx_snow_grnd_col(c) > 0.0_rkx .and. &
           frac_sno(c)*snow_depth(c) >= 0.01_rkx ) then
        newnode = 1
        snl(c) = -1
        dz(c,0) = snow_depth(c)                       ! meter
        z(c,0) = -0.5_rkx*dz(c,0)
        zi(c,-1) = -dz(c,0)
        t_soisno(c,0) = min(tfrz, forc_t(c))      ! K
        h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
        h2osoi_liq(c,0) = 0._rkx                   ! kg/m2
        frac_iceold(c,0) = 1._rkx

        ! intitialize SNICAR variables for fresh snow:
        snw_rds(c,0)    = snw_rds_min

        mss_bcpho(c,:)  = 0._rkx
        mss_bcphi(c,:)  = 0._rkx
        mss_bctot(c,:)  = 0._rkx
        mss_bc_col(c)   = 0._rkx
        mss_bc_top(c)   = 0._rkx

        mss_ocpho(c,:)  = 0._rkx
        mss_ocphi(c,:)  = 0._rkx
        mss_octot(c,:)  = 0._rkx
        mss_oc_col(c)   = 0._rkx
        mss_oc_top(c)   = 0._rkx

        mss_dst1(c,:)   = 0._rkx
        mss_dst2(c,:)   = 0._rkx
        mss_dst3(c,:)   = 0._rkx
        mss_dst4(c,:)   = 0._rkx
        mss_dsttot(c,:) = 0._rkx
        mss_dst_col(c)  = 0._rkx
        mss_dst_top(c)  = 0._rkx
      end if

      ! The change of ice partial density of surface node due to precipitation.
      ! Only ice part of snowfall is added here, the liquid part will be added
      ! later.

      if ( snl(c) < 0 .and. newnode == 0 ) then
        h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+newsnow(c)
        dz(c,snl(c)+1) = dz(c,snl(c)+1) + dz_snowf*dtsrf
      end if
    end do

    ! update surface water fraction (this may modify frac_sno)
    call FracH2oSfc(lbc, ubc, num_nolakec, filter_nolakec,frac_h2osfc)

  end subroutine Hydrology1

end module mod_clm_hydrology1
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
