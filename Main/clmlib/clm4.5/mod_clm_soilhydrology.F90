module mod_clm_soilhydrology
  !
  ! Calculate soil hydrology
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_dynparam
  use mod_mpmessage
  use mod_runparams , only : dtsrf
  use mod_clm_type
  use mod_clm_varcon , only : e_ice , denh2o , denice , rpi , pc , mu
  use mod_clm_varcon , only : wimp , pondmx_urban , roverg ,  &
           icol_roof, icol_sunwall, icol_shadewall, istsoil , &
           icol_road_imperv, icol_road_perv, isturb , tfrz , istcrop
  use mod_clm_varcon , only : grav , hfus , pondmx , watmin
  use mod_clm_varpar , only : maxpatch_pft
#if (defined VICHYDRO)
  use mod_clm_varpar , only : nlayer , nlayert
  use mod_clm_varcon , only : secspday , nlvic
  use mod_clm_vicmap , only : CLMVICMap
#endif
  use mod_clm_varpar , only : nlevsoi , max_pft_per_col , nlevgrnd
  use mod_clm_h2osfc , only : FracH2oSfc
  use mod_clm_tridiagonal , only : Tridiagonal

  implicit none

  private

  save

  public :: SoilHydrology_readnl ! Initialization for Soil Hydrology

  public :: SurfaceRunoff    ! Calculate surface runoff
  public :: Infiltration     ! Calculate infiltration into surface soil layer
  public :: SoilWater        ! Calculate soil hydrology
  public :: Drainage         ! Calculate subsurface drainage

  !If surface water is active or not
  integer(ik4) , public :: h2osfcflag = 1

  !use control soil hydraulic properties
  integer(ik4) :: origflag = 0

  contains
  !
  ! Read namelist for SoilHydrology
  !
  subroutine SoilHydrology_readnl( NLFilename )
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    integer(ik4) :: ierr                 ! error code
    integer(ik4) :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'SoilHydrology_readnl'  ! subroutine name

    namelist / clm_soilhydrology_inparm / h2osfcflag, origflag

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( myid == iocpu )then
      unitn = file_getUnit( )
      write(stdout,*) 'Read in clm_soilhydrology_inparm namelist'
      open(unitn,file=NLFilename,status='old',action='read',iostat=ierr)
      if (ierr /= 0) then
        call fatal(__FILE__,__LINE__, &
           subname // ':: ERROR open namelist file '//NLFilename)
      end if
      read(unitn, clm_soilhydrology_inparm, iostat=ierr)
      if (ierr /= 0) then
        call fatal(__FILE__,__LINE__,&
            subname // ':: ERROR reading clm_soilhydrology_inparm namelist')
      end if
      call file_freeUnit( unitn )
    end if
    ! Broadcast namelist variables read in
    call bcast(h2osfcflag)
    call bcast(origflag)
  end subroutine SoilHydrology_readnl
  !
  ! Calculate surface runoff
  !
  subroutine SurfaceRunoff (lbc, ubc, num_hydrologyc, &
                  filter_hydrologyc, num_urbanc, filter_urbanc, icefrac)
    implicit none
    integer(ik4) , intent(in)  :: lbc, ubc   ! column bounds
    ! number of column soil points in column filter
    integer(ik4) , intent(in)  :: num_hydrologyc
    ! column filter for soil points
    integer(ik4) , intent(in)  :: filter_hydrologyc(ubc-lbc+1)
    ! number of column urban points in column filter
    integer(ik4) , intent(in)  :: num_urbanc
    ! column filter for urban points
    integer(ik4) , intent(in)  :: filter_urbanc(ubc-lbc+1)
    ! fraction of ice in layer (-)
    real(rk8), intent(out) :: icefrac(lbc:ubc,1:nlevsoi)
    ! snow falling on surface water (mm/s)
    real(rk8), pointer :: qflx_snow_h2osfc(:)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer :: frac_h2osfc(:)
    ! frost table depth (m)
    real(rk8), pointer :: frost_table(:)
    ! perched water table depth (m)
    real(rk8), pointer :: zwt_perched(:)
    ! column flux of flood water from RTM
    real(rk8), pointer :: qflx_floodc(:)
    integer(ik4) , pointer :: clandunit(:)  ! column's landunit
    integer(ik4) , pointer :: ltype(:)      ! landunit type
    integer(ik4) , pointer :: cgridcell(:)  ! gridcell index for each column
    integer(ik4) , pointer :: ctype(:)      ! column type index
    ! net water input into soil from top (mm/s)
    real(rk8), pointer :: qflx_top_soil(:)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: watsat(:,:)
    ! decay factor (m)
    real(rk8), pointer :: hkdepth(:)
    ! water table depth (m)
    real(rk8), pointer :: zwt(:)
    ! fractional impermeable area
    real(rk8), pointer :: fcov(:)
    ! fractional area with water table at surface
    real(rk8), pointer :: fsat(:)
    ! layer depth (m)
    real(rk8), pointer :: dz(:,:)
    ! ice lens (kg/m2)
    real(rk8), pointer :: h2osoi_ice(:,:)
    ! liquid water (kg/m2)
    real(rk8), pointer :: h2osoi_liq(:,:)
    ! maximum saturated fraction for a gridcell
    real(rk8), pointer :: wtfact(:)
    ! hydraulic conductivity at saturation (mm H2O /s)
    real(rk8), pointer :: hksat(:,:)
    ! Clapp and Hornberger "b"
    real(rk8), pointer :: bsw(:,:)
    ! minimum soil suction (mm)
    real(rk8), pointer :: sucsat(:,:)
    ! minus number of snow layers
    integer(ik4) , pointer :: snl(:)
    ! ground surface evaporation rate (mm H2O/s) [+]
    real(rk8), pointer :: qflx_evap_grnd(:)
    ! interface level below a "z" level (m)
    real(rk8), pointer :: zi(:,:)
    ! surface runoff (mm H2O /s)
    real(rk8), pointer :: qflx_surf(:)
    ! effective porosity = porosity - vol_ice
    real(rk8), pointer :: eff_porosity(:,:)
    !fractional impermeability (-)
    real(rk8), pointer :: fracice(:,:)
#if (defined VICHYDRO)
    !VIC b infiltration parameter
    real(rk8), pointer :: b_infil(:)
    !maximum soil moisture (ice + liq, mm)
    real(rk8), pointer :: max_moist(:,:)
    !soil moisture in each VIC layers (liq, mm)
    real(rk8), pointer :: moist(:,:)
    !ice len in each VIC layers(ice, mm)
    real(rk8), pointer :: ice(:,:)
    !maximum infiltration capacity in VIC (mm)
    real(rk8), pointer :: max_infil(:)
    !column average soil moisture in top VIC layers (mm)
    real(rk8), pointer :: i_0(:)
#endif
    integer(ik4)  :: c,j,fc    !indices
    ! excess soil water above urban ponding limit
    real(rk8) :: xs(lbc:ubc)
    !partial volume of ice lens in layer
    real(rk8) :: vol_ice(lbc:ubc,1:nlevsoi)
    !decay factor (m-1)
    real(rk8) :: fff(lbc:ubc)
#if (defined VICHYDRO)
    !fraction of the saturated area
    real(rk8) :: A(lbc:ubc)
    !temporary variable (exponent)
    real(rk8) :: ex(lbc:ubc)
    !temporary, soil moisture in top VIC layers
    real(rk8) :: top_moist(lbc:ubc)
    !temporary, maximum soil moisture in top VIC layers
    real(rk8) :: top_max_moist(lbc:ubc)
    !temporary, ice len in top VIC layers
    real(rk8) :: top_ice(lbc:ubc)
    !temporary, ice fraction in top VIC layers
    real(rk8) :: top_icefrac
    !temporary, fraction covered by ice for runoff calculations
    real(rk8) :: top_fracice
    character(len=32) :: subname = 'SurfaceRunoff'  ! subroutine name
#endif

    ! Assign local pointers to derived subtype components (column-level)

    qflx_snow_h2osfc  => clm3%g%l%c%cwf%qflx_snow_h2osfc
    frac_h2osfc       => clm3%g%l%c%cps%frac_h2osfc
    frost_table       => clm3%g%l%c%cws%frost_table
    zwt_perched       => clm3%g%l%c%cws%zwt_perched
    qflx_floodc       => clm3%g%l%c%cwf%qflx_floodc
    cgridcell         => clm3%g%l%c%gridcell
    clandunit         => clm3%g%l%c%landunit
    ltype             => clm3%g%l%itype
    ctype             => clm3%g%l%c%itype
    qflx_top_soil     => clm3%g%l%c%cwf%qflx_top_soil
    qflx_surf         => clm3%g%l%c%cwf%qflx_surf
    watsat            => clm3%g%l%c%cps%watsat
    hkdepth           => clm3%g%l%c%cps%hkdepth
    dz                => clm3%g%l%c%cps%dz
    h2osoi_ice        => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq        => clm3%g%l%c%cws%h2osoi_liq
    fcov              => clm3%g%l%c%cws%fcov
    fsat              => clm3%g%l%c%cws%fsat
    eff_porosity      => clm3%g%l%c%cps%eff_porosity
    wtfact            => clm3%g%l%c%cps%wtfact
    zwt               => clm3%g%l%c%cws%zwt
    fracice           => clm3%g%l%c%cps%fracice
    hksat             => clm3%g%l%c%cps%hksat
    bsw               => clm3%g%l%c%cps%bsw
    sucsat            => clm3%g%l%c%cps%sucsat
    snl               => clm3%g%l%c%cps%snl
    qflx_evap_grnd    => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    zi                => clm3%g%l%c%cps%zi
#if (defined VICHYDRO)
    b_infil           => clm3%g%l%c%cps%b_infil
    max_moist         => clm3%g%l%c%cps%max_moist
    moist             => clm3%g%l%c%cws%moist
    ice               => clm3%g%l%c%cws%ice
    max_infil         => clm3%g%l%c%cws%max_infil
    i_0               => clm3%g%l%c%cws%i_0
#endif

    do j = 1 , nlevsoi
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)

        ! Porosity of soil, partial volume of ice and liquid,
        ! fraction of ice in each layer, fractional impermeability

        vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
        if ( origflag == 1 ) then
          icefrac(c,j) = min(1._rk8,h2osoi_ice(c,j)/(h2osoi_ice(c,j) + &
                                  h2osoi_liq(c,j)))
        else
          icefrac(c,j) = min(1._rk8,vol_ice(c,j)/watsat(c,j))
        end if

        fracice(c,j) = max(0._rk8,exp(-3._rk8*(1._rk8-icefrac(c,j))) - &
                exp(-3._rk8))/(1.0_rk8-exp(-3._rk8))
      end do
    end do

    ! Saturated fraction

    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)
      fff(c) = 0.5_rk8
#if (defined VICHYDRO)
      top_moist(c) = 0._rk8
      top_ice(c) = 0._rk8
      top_max_moist(c) = 0._rk8
      do j = 1 , nlayer - 1
        top_ice(c) = top_ice(c) + ice(c,j)
        top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
        top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
      end do
      if ( top_moist(c)> top_max_moist(c) ) top_moist(c)= top_max_moist(c)
      top_ice(c)     = max(0._rk8,top_ice(c))
      max_infil(c)   = (1._rk8+b_infil(c)) * top_max_moist(c)
      ex(c)          = b_infil(c) / (1._rk8 + b_infil(c))
      A(c)           = 1._rk8 - (1._rk8 - top_moist(c) / top_max_moist(c))**ex(c)
      i_0(c)         = max_infil(c) * (1._rk8 - (1._rk8 - A(c))**(1._rk8/b_infil(c)))
      fsat(c)        = A(c)  !for output
#else
      fsat(c) = wtfact(c) * exp(-0.5_rk8*fff(c)*zwt(c))
#endif

      ! use perched water table to determine fsat (if present)
      if ( frost_table(c) > zwt(c) ) then
#if (defined VICHYDRO)
        fsat(c) =  A(c)
#else
        fsat(c) = wtfact(c) * exp(-0.5_rk8*fff(c)*zwt(c))
#endif
      else
        if ( frost_table(c) > zwt_perched(c)) then
          fsat(c) = wtfact(c) * exp(-0.5_rk8*fff(c)*zwt_perched(c))
          !*( frost_table(c) - zwt_perched(c))/4.0
        end if
      end if
      if ( origflag == 1 ) then
#if (defined VICHYDRO)
        call fatal(__FILE__,__LINE__,&
             subname // ':: VICHYDRO is not available for origflag = 1')
#else
        fcov(c) = (1._rk8 - fracice(c,1)) * fsat(c) + fracice(c,1)
#endif
      else
        fcov(c) = fsat(c)
      end if
    end do

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)

      ! assume qinmax large relative to qflx_top_soil in control
      if (origflag == 1) then
        qflx_surf(c) = fcov(c) * qflx_top_soil(c)
      else
        ! only send fast runoff directly to streams
        qflx_surf(c) = fsat(c) * qflx_top_soil(c)
      end if
    end do

    ! Determine water in excess of ponding limit for urban roof and
    ! impervious road. Excess goes to surface runoff.
    ! No surface runoff for sunwall and shadewall.

    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      if ( ctype(c) == icol_roof .or. &
           ctype(c) == icol_road_imperv ) then
        ! If there are snow layers then all qflx_top_soil goes to surface runoff
        if (snl(c) < 0) then
          qflx_surf(c) = max(0._rk8,qflx_top_soil(c))
        else
          xs(c) = max(0._rk8, &
            h2osoi_liq(c,1)/dtsrf + qflx_top_soil(c) - &
                                    qflx_evap_grnd(c) - pondmx_urban/dtsrf)
          if (xs(c) > 0.0_rk8 ) then
            h2osoi_liq(c,1) = pondmx_urban
          else
            h2osoi_liq(c,1) = max(0._rk8,h2osoi_liq(c,1)+ &
                     (qflx_top_soil(c)-qflx_evap_grnd(c))*dtsrf)
          end if
          qflx_surf(c) = xs(c)
        end if
      else if ( ctype(c) == icol_sunwall .or. &
                ctype(c) == icol_shadewall ) then
        qflx_surf(c) = 0._rk8
      end if
      ! send flood water flux to runoff for all urban columns
      qflx_surf(c) = qflx_surf(c) + qflx_floodc(c)
    end do

    ! remove stormflow and snow on h2osfc from qflx_top_soil
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      ! add flood water flux to qflx_top_soil
      qflx_top_soil(c) = qflx_top_soil(c) + qflx_snow_h2osfc(c) + qflx_floodc(c)
    end do
  end subroutine SurfaceRunoff
  !
  ! Calculate infiltration into surface soil layer (minus the evaporation)
  !
  subroutine Infiltration(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                          num_urbanc, filter_urbanc, vol_liq)
    implicit none
    integer(ik4), intent(in) :: lbc, ubc    ! column bounds
    ! number of column soil points in column filter
    integer(ik4), intent(in) :: num_hydrologyc
    ! column filter for soil points
    integer(ik4), intent(in) :: filter_hydrologyc(ubc-lbc+1)
    ! number of column urban points in column filter
    integer(ik4), intent(in) :: num_urbanc
    ! column filter for urban points
    integer(ik4), intent(in) :: filter_urbanc(ubc-lbc+1)
    ! partial volume of liquid water in layer
    real(rk8), intent(out) :: vol_liq(lbc:ubc,1:nlevsoi)

    real(rk8), pointer :: frost_table(:)    ! frost table depth (m)
    real(rk8), pointer :: zwt_perched(:)    ! perched water table depth (m)
    ! fractional area with water table at surface
    real(rk8), pointer :: fsat(:)
    integer(ik4) , pointer :: clandunit(:)  ! column's landunit
    integer(ik4) , pointer :: ltype(:)      ! landunit type
    real(rk8), pointer :: h2osfc_thresh(:)  ! level at which h2osfc "percolates"
    ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno(:)
    ! ground surface evaporation rate (mm H2O/s) [+]
    real(rk8), pointer :: qflx_evap_soi(:)
    real(rk8), pointer :: qflx_h2osfc_surf(:) ! surface water runoff (mm/s)
    real(rk8), pointer :: h2osfc(:)           ! surface water (mm)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer :: frac_h2osfc(:)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer :: frac_h2osfc_temp(:)
    real(rk8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
    real(rk8), pointer :: h2osoi_ice(:,:) ! ice lens (kg/m2)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: watsat(:,:)
    real(rk8), pointer :: sucsat(:,:)    ! minimum soil suction (mm)
    real(rk8), pointer :: bsw(:,:)       ! Clapp and Hornberger "b"
    real(rk8), pointer :: t_soisno(:,:)  ! soil temperature (Kelvin)
    real(rk8), pointer :: t_h2osfc(:)    ! soil temperature (Kelvin)
    ! restriction for min of soil potential (mm)
    real(rk8), pointer :: smpmin(:)
    real(rk8), pointer :: dz(:,:)    ! layer depth (m)
    ! hydraulic conductivity at saturation (mm H2O /s)
    real(rk8), pointer :: hksat(:,:)
    real(rk8), pointer :: hksat_min(:,:)    ! mineral hksat
    ! fractional area with water table at surface
    real(rk8), pointer :: fcov(:)
    ! effective porosity = porosity - vol_ice
    real(rk8), pointer :: eff_porosity(:,:)
    real(rk8), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(rk8), pointer :: snow_depth(:)     ! snow height (m)
    real(rk8), pointer :: topo_slope(:)     ! topographic slope
    ! evaporation flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_snow(:)
    ! evaporation flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_soil(:)
    ! evaporation flux from h2osfc (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_h2osfc(:)
    real(rk8), pointer :: zwt(:)           ! water table depth (m)
    integer(ik4) , pointer :: ctype(:)     ! column type index
    integer(ik4) , pointer :: snl(:)       ! minus number of snow layers
    ! net water input into soil from top (mm/s)
    real(rk8), pointer :: qflx_top_soil(:)
    real(rk8), pointer :: qflx_surf(:)      ! surface runoff (mm H2O /s)
    ! ground surface evaporation rate (mm H2O/s) [+]
    real(rk8), pointer :: qflx_evap_grnd(:)

#if (defined VICHYDRO)
    real(rk8), pointer :: b_infil(:)   !VIC b infiltration parameter
    !maximum soil moisture (ice + liq, mm)
    real(rk8), pointer :: max_moist(:,:)
    !soil moisture in each VIC layers (liq, mm)
    real(rk8), pointer :: moist(:,:)
    !ice len in each VIC layers(ice, mm)
    real(rk8), pointer :: ice(:,:)
    !maximum infiltration capacity in VIC (mm)
    real(rk8), pointer :: max_infil(:)
    !column average soil moisture in top VIC layers (mm)
    real(rk8), pointer :: i_0(:)
#endif
    real(rk8), pointer :: qflx_infl(:)  !infiltration (mm H2O /s)

    integer(ik4) :: c,j,fc  ! indices
    real(rk8) :: qinmax     ! maximum infiltration capacity (mm/s)
    ! partial volume of ice lens in layer
    real(rk8) :: vol_ice(lbc:ubc,1:nlevsoi)
    real(rk8) :: qflx_evap(lbc:ubc)         ! local evaporation array
    real(rk8) :: qflx_h2osfc_drain(lbc:ubc) ! bottom drainage from h2osfc
    real(rk8) :: qflx_in_h2osfc(lbc:ubc)    ! surface input to h2osfc
    real(rk8) :: qflx_in_soil(lbc:ubc)      ! surface input to soil
    ! infiltration excess runoff -> h2osfc
    real(rk8) :: qflx_infl_excess(lbc:ubc)
    real(rk8) :: frac_infclust  ! fraction of submerged area that is connected
    real(rk8) :: fsno     ! copy of frac_sno
    real(rk8) :: k_wet    ! linear reservoir coefficient for h2osfc
    real(rk8) :: icefrac(lbc:ubc,1:nlevsoi) !
#if (defined VICHYDRO)
    ! temporary, variable soil moisture holding capacity
    ! in top VIC layers for runoff calculation
    real(rk8) :: basis
    ! temp VIC surface runoff
    real(rk8) :: rsurf_vic
    ! temporary, soil moisture in top VIC layers
    real(rk8) :: top_moist(lbc:ubc)
    ! temporary, maximum soil moisture in top VIC layers
    real(rk8) :: top_max_moist(lbc:ubc)
    ! temporary, ice len in top VIC layers
    real(rk8) :: top_ice(lbc:ubc)
    ! temporary, ice fraction in top VIC layers
    real(rk8) :: top_icefrac
    ! temporary, fraction covered by ice for runoff calculations
    real(rk8) :: top_fracice
#endif

    ! Assign local pointers to derived type members (column-level)

    frost_table    => clm3%g%l%c%cws%frost_table
    zwt_perched    => clm3%g%l%c%cws%zwt_perched
    fsat           => clm3%g%l%c%cws%fsat
    h2osfc_thresh  => clm3%g%l%c%cps%h2osfc_thresh
    frac_sno       => clm3%g%l%c%cps%frac_sno_eff
    qflx_evap_soi  => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    qflx_h2osfc_surf  => clm3%g%l%c%cwf%qflx_h2osfc_surf
    frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc
    frac_h2osfc_temp => clm3%g%l%c%cps%frac_h2osfc_temp
    h2osfc         => clm3%g%l%c%cws%h2osfc
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    sucsat         => clm3%g%l%c%cps%sucsat
    watsat         => clm3%g%l%c%cps%watsat
    bsw            => clm3%g%l%c%cps%bsw
    t_soisno       => clm3%g%l%c%ces%t_soisno
    smpmin         => clm3%g%l%c%cps%smpmin
    fcov           => clm3%g%l%c%cws%fcov
    eff_porosity   => clm3%g%l%c%cps%eff_porosity
    hksat          => clm3%g%l%c%cps%hksat
    hksat_min      => clm3%g%l%c%cps%hksat_min
    dz             => clm3%g%l%c%cps%dz
    h2osno         => clm3%g%l%c%cws%h2osno
    snow_depth     => clm3%g%l%c%cps%snow_depth
    t_h2osfc       => clm3%g%l%c%ces%t_h2osfc
    clandunit      => clm3%g%l%c%landunit
    ltype          => clm3%g%l%itype
    topo_slope     => clm3%g%l%c%cps%topo_slope
    qflx_ev_snow   => clm3%g%l%c%cwf%pwf_a%qflx_ev_snow
    qflx_ev_soil   => clm3%g%l%c%cwf%pwf_a%qflx_ev_soil
    qflx_ev_h2osfc => clm3%g%l%c%cwf%pwf_a%qflx_ev_h2osfc
    zwt            => clm3%g%l%c%cws%zwt
    ctype          => clm3%g%l%c%itype
    snl            => clm3%g%l%c%cps%snl
    qflx_top_soil  => clm3%g%l%c%cwf%qflx_top_soil
    qflx_surf      => clm3%g%l%c%cwf%qflx_surf
    qflx_infl      => clm3%g%l%c%cwf%qflx_infl
    qflx_evap_grnd => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
#if (defined VICHYDRO)
    b_infil        => clm3%g%l%c%cps%b_infil
    max_moist      => clm3%g%l%c%cps%max_moist
    moist          => clm3%g%l%c%cws%moist
    ice            => clm3%g%l%c%cws%ice
    max_infil      => clm3%g%l%c%cws%max_infil
    i_0            => clm3%g%l%c%cws%i_0
#endif

    ! Infiltration into surface soil layer (minus the evaporation)
    do j = 1,nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        ! Porosity of soil, partial volume of ice and liquid
        vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
        eff_porosity(c,j) = max(0.01_rk8,watsat(c,j)-vol_ice(c,j))
        vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(dz(c,j)*denh2o))
        icefrac(c,j) = min(1._rk8,vol_ice(c,j)/watsat(c,j))
      end do
    end do

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      ! partition moisture fluxes between soil and h2osfc
      if ( ltype(clandunit(c)) == istsoil .or. &
           ltype(clandunit(c)) == istcrop ) then

        ! explicitly use frac_sno=0 if snl=0
        if (snl(c) >= 0) then
          fsno = 0._rk8
          ! if no snow layers, sublimation is removed from h2osoi_ice
          ! in drainage
          qflx_evap(c) = qflx_evap_grnd(c)
        else
          fsno = frac_sno(c)
          qflx_evap(c) = qflx_ev_soil(c)
        end if

        !1. partition surface inputs between soil and h2osfc
        qflx_in_soil(c) = (1._rk8 - frac_h2osfc(c)) * &
                (qflx_top_soil(c)  - qflx_surf(c))
        qflx_in_h2osfc(c) = frac_h2osfc(c) * &
                (qflx_top_soil(c)  - qflx_surf(c))

        !2. remove evaporation (snow treated in SnowHydrology)
        qflx_in_soil(c) = qflx_in_soil(c) - &
                (1.0_rk8 - fsno - frac_h2osfc(c))*qflx_evap(c)
        qflx_in_h2osfc(c) = qflx_in_h2osfc(c) - &
                frac_h2osfc(c) * qflx_ev_h2osfc(c)

        !3. determine maximum infiltration rate
#if (defined VICHYDRO)
        top_moist(c) = 0._rk8
        top_ice(c) = 0._rk8
        top_max_moist(c) = 0._rk8
        do j = 1, nlayer - 1
          top_ice(c) = top_ice(c) + ice(c,j)
          top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
          top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
        end do
        top_icefrac = min(1._rk8,top_ice(c)/top_max_moist(c))
        if ( qflx_in_soil(c) <= 0._rk8 ) then
          rsurf_vic = 0._rk8
        else if ( max_infil(c) <= 0._rk8 ) then
          rsurf_vic = qflx_in_soil(c)
        else if ( (i_0(c) + qflx_in_soil(c)*dtsrf) > max_infil(c) ) then
          !(Eq.(3a) Wood et al. 1992)
          rsurf_vic = (qflx_in_soil(c)*dtsrf - &
                  top_max_moist(c) + top_moist(c))/dtsrf
        else
          !(Eq.(3b) Wood et al. 1992)
          basis = 1._rk8 - (i_0(c) + qflx_in_soil(c)*dtsrf)/max_infil(c)
          rsurf_vic = (qflx_in_soil(c)*dtsrf - &
                  top_max_moist(c) + top_moist(c)    &
                + top_max_moist(c) * basis**(1._rk8 + b_infil(c)))/dtsrf
        end if
        rsurf_vic = max(0._rk8, rsurf_vic)
        qinmax = (1._rk8 - fsat(c)) * &
                10._rk8**(-e_ice*top_icefrac)*(qflx_in_soil(c) - rsurf_vic)
#else
        qinmax = (1._rk8 - fsat(c)) * &
                minval(10._rk8**(-e_ice*(icefrac(c,1:3)))*hksat(c,1:3))
#endif
        qflx_infl_excess(c) = max(0._rk8,qflx_in_soil(c) - &
                (1.0_rk8 - frac_h2osfc(c))*qinmax)

        !4. soil infiltration and h2osfc "run-on"
        qflx_infl(c) = qflx_in_soil(c) - qflx_infl_excess(c)
        qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) + qflx_infl_excess(c)

        !5. surface runoff from h2osfc
        if ( h2osfcflag == 1 ) then
          ! calculate runoff from h2osfc
          if ( frac_h2osfc(c) <= pc ) then
            frac_infclust = 0.0_rk8
            ! there is a potential conflict between frac_sno and frac_h2osfc
            ! calculate temporary surface water fraction to enable runoff
            ! when frac_sno is large
            if ( h2osfc(c) >= h2osfc_thresh(c) ) then
              call FracH2oSfc(lbc, ubc, num_hydrologyc, &
                      filter_hydrologyc,frac_h2osfc_temp,1)
              frac_infclust = (frac_h2osfc_temp(c)-pc)**mu
            end if
          else
            frac_infclust = (frac_h2osfc(c)-pc)**mu
          end if
        end if

        ! limit runoff to value of storage above S(pc)
        if ( h2osfc(c) >= h2osfc_thresh(c) .and. &
             h2osfcflag /= 0 ) then
          ! spatially variable k_wet
          k_wet = 1.0_rk8 * sin((rpi/180.0_rk8) * topo_slope(c))
          qflx_h2osfc_surf(c) = k_wet * frac_infclust * &
                  (h2osfc(c) - h2osfc_thresh(c))

          qflx_h2osfc_surf(c) = min(qflx_h2osfc_surf(c),(h2osfc(c) - &
                  h2osfc_thresh(c))/dtsrf)
        else
          qflx_h2osfc_surf(c) = 0._rk8
        end if

        ! cutoff lower limit
        if ( qflx_h2osfc_surf(c) < 1.0e-8_rk8 ) qflx_h2osfc_surf(c) = 0._rk8

        ! use this for non-h2osfc code
        if ( h2osfcflag == 0 ) then
          qflx_h2osfc_surf(c) = 0._rk8
          ! shift infiltration excess from h2osfc input to surface runoff
          qflx_in_h2osfc(c) = qflx_in_h2osfc(c) - qflx_infl_excess(c)
          qflx_surf(c) = qflx_surf(c) + qflx_infl_excess(c)
          qflx_infl_excess(c) = 0._rk8
        end if

        qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) - qflx_h2osfc_surf(c)

        !6. update h2osfc prior to calculating bottom drainage from h2osfc
        h2osfc(c) = h2osfc(c) + qflx_in_h2osfc(c) * dtsrf
        !--  if all water evaporates, there will be no bottom drainage
        if ( h2osfc(c) < 0.0 ) then
          qflx_infl(c) = qflx_infl(c) + h2osfc(c)/dtsrf
          h2osfc(c) = 0.0
          qflx_h2osfc_drain(c) = 0._rk8
        else
          qflx_h2osfc_drain(c) = min(frac_h2osfc(c)*qinmax,h2osfc(c)/dtsrf)
        end if

        if ( h2osfcflag == 0 ) then
          qflx_h2osfc_drain(c) = max(0._rk8,h2osfc(c)/dtsrf) !ensure no h2osfc
        end if

        !7. remove drainage from h2osfc and add to qflx_infl
        h2osfc(c) = h2osfc(c) - qflx_h2osfc_drain(c) * dtsrf
        qflx_infl(c) = qflx_infl(c) + qflx_h2osfc_drain(c)

        !#######################################################
      else
        ! non-vegetated landunits (i.e. urban) use original CLM4 code
        if ( snl(c) >= 0 ) then
          ! when no snow present, sublimation is removed in Drainage
          qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) - qflx_evap_grnd(c)
        else
          qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) &
               - (1.0_rk8 - frac_sno(c)) * qflx_ev_soil(c)
        end if
        qflx_h2osfc_surf(c) = 0._rk8
      end if
    end do

    ! No infiltration for impervious urban surfaces

    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      if ( ctype(c) /= icol_road_perv ) then
        qflx_infl(c) = 0._rk8
      end if
    end do
  end subroutine Infiltration
  !
  ! Soil hydrology
  ! Soil moisture is predicted from a 10-layer model (as with soil
  ! temperature), in which the vertical soil moisture transport is governed
  ! by infiltration, runoff, gradient diffusion, gravity, and root
  ! extraction through canopy transpiration.  The net water applied to the
  ! surface layer is the snowmelt plus precipitation plus the throughfall
  ! of canopy dew minus surface runoff and evaporation.
  ! CLM3.5 uses a zero-flow bottom boundary condition.
  !
  ! The vertical water flow in an unsaturated porous media is described by
  ! Darcy's law, and the hydraulic conductivity and the soil negative
  ! potential vary with soil water content and soil texture based on the work
  ! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
  ! integrated over the layer thickness, in which the time rate of change in
  ! water mass must equal the net flow across the bounding interface, plus the
  ! rate of internal source or sink. The terms of water flow across the layer
  ! interfaces are linearly expanded by using first-order Taylor expansion.
  ! The equations result in a tridiagonal system equation.
  !
  ! Note: length units here are all millimeter
  ! (in temperature subroutine uses same soil layer
  ! structure required but lengths are m)
  !
  ! Richards equation:
  !
  ! d wat      d     d wat d psi
  ! ----- = - -- [ k(----- ----- - 1) ] + S
  !   dt      dz       dz  d wat
  !
  ! where: wat = volume of water per volume of soil (mm**3/mm**3)
  ! psi = soil matrix potential (mm)
  ! dt  = time step (s)
  ! z   = depth (mm)
  ! dz  = thickness (mm)
  ! qin = inflow at top (mm h2o /s)
  ! qout= outflow at bottom (mm h2o /s)
  ! s   = source/sink flux (mm h2o /s)
  ! k   = hydraulic conductivity (mm h2o /s)
  !
  !                       d qin                  d qin
  ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
  !                       d wat(j-1)             d wat(j)
  !                ==================|=================
  !                                  < qin
  !
  !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
  !
  !                                  > qout
  !                ==================|=================
  !                        d qout               d qout
  ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
  !                        d wat(j)             d wat(j+1)
  !
  !
  ! Solution: linearize k and psi about d wat and use tridiagonal
  ! system of equations to solve for d wat,
  ! where for layer j
  !
  !
  ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
  !
  subroutine SoilWater(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                       dwat, hk, dhkdw)
    implicit none
    integer(ik4) , intent(in)  :: lbc, ubc      ! column bounds
    ! number of column soil points in column filter
    integer(ik4) , intent(in)  :: num_hydrologyc
    ! column filter for soil points
    integer(ik4) , intent(in)  :: filter_hydrologyc(ubc-lbc+1)
    ! change of soil water [m3/m3]
    real(rk8), intent(out) :: dwat(lbc:ubc,1:nlevsoi)
    ! hydraulic conductivity [mm h2o/s]
    real(rk8), intent(out) :: hk(lbc:ubc,1:nlevsoi)
    ! d(hk)/d(vol_liq)
    real(rk8), intent(out) :: dhkdw(lbc:ubc,1:nlevsoi)

    real(rk8), pointer :: h2osoi_ice(:,:)  ! ice water (kg/m2)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: ctype(:)  ! column type index
    integer(ik4) , pointer :: npfts(:)  ! column's number of pfts - ADD
    real(rk8), pointer :: pwtcol(:)     ! weight relative to column for each pft
    real(rk8), pointer :: z(:,:)        ! layer depth (m)
    real(rk8), pointer :: dz(:,:)       ! layer thickness (m)
    ! restriction for min of soil potential (mm)
    real(rk8), pointer :: smpmin(:)
    real(rk8), pointer :: qflx_infl(:)  ! infiltration (mm H2O /s)
    ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8), pointer :: qflx_tran_veg_pft(:)
    ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8), pointer :: qflx_tran_veg_col(:)
    ! effective porosity = porosity - vol_ice
    real(rk8), pointer :: eff_porosity(:,:)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: watsat(:,:)
    ! hydraulic conductivity at saturation (mm H2O /s)
    real(rk8), pointer :: hksat(:,:)
    real(rk8), pointer :: bsw(:,:)     ! Clapp and Hornberger "b"
    real(rk8), pointer :: sucsat(:,:)  ! minimum soil suction (mm)
    real(rk8), pointer :: t_soisno(:,:) ! soil temperature (Kelvin)
    ! effective fraction of roots in each soil layer
    real(rk8), pointer :: rootr_pft(:,:)
    ! beginning pft index for each column
    integer(ik4) , pointer :: pfti(:)
    ! fractional impermeability (-)
    real(rk8), pointer :: fracice(:,:)
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8), pointer :: h2osoi_vol(:,:)
    real(rk8), pointer :: qcharge(:)  ! aquifer recharge rate (mm/s)
    real(rk8), pointer :: hkdepth(:)  ! decay factor (m)
    real(rk8), pointer :: zwt(:)      ! water table depth (m)
    real(rk8), pointer :: zi(:,:)     ! interface level below a "z" level (m)

    real(rk8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
    ! effective fraction of roots in each soil layer
    real(rk8), pointer :: rootr_col(:,:)
    real(rk8), pointer :: smp_l(:,:)   ! soil matrix potential [mm]
    real(rk8), pointer :: hk_l(:,:)    ! hydraulic conductivity (mm/s)

    integer(ik4)  :: p,c,fc,j          ! do loop indices
    integer(ik4)  :: jtop(lbc:ubc)     ! top level at each column
    ! "a" left off diagonal of tridiagonal matrix
    real(rk8) :: amx(lbc:ubc,1:nlevsoi+1)
    ! "b" diagonal column for tridiagonal matrix
    real(rk8) :: bmx(lbc:ubc,1:nlevsoi+1)
    ! "c" right off diagonal tridiagonal matrix
    real(rk8) :: cmx(lbc:ubc,1:nlevsoi+1)
    ! "r" forcing term of tridiagonal matrix
    real(rk8) :: rmx(lbc:ubc,1:nlevsoi+1)
    real(rk8) :: zmm(lbc:ubc,1:nlevsoi+1)  ! layer depth [mm]
    real(rk8) :: dzmm(lbc:ubc,1:nlevsoi+1) ! layer thickness [mm]
    real(rk8) :: den                       ! used in calculating qin, qout
    real(rk8) :: dqidw0(lbc:ubc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i-1))
    real(rk8) :: dqidw1(lbc:ubc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i))
    real(rk8) :: dqodw1(lbc:ubc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i))
    real(rk8) :: dqodw2(lbc:ubc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i+1))
    real(rk8) :: dsmpdw(lbc:ubc,1:nlevsoi+1) ! d(smp)/d(vol_liq)
    real(rk8) :: num                         ! used in calculating qin, qout
    ! flux of water into soil layer [mm h2o/s]
    real(rk8) :: qin(lbc:ubc,1:nlevsoi+1)
    ! flux of water out of soil layer [mm h2o/s]
    real(rk8) :: qout(lbc:ubc,1:nlevsoi+1)
    real(rk8) :: s_node                    ! soil wetness
    real(rk8) :: s1                        ! "s" at interface of layer
    real(rk8) :: s2                        ! k*s**(2b+2)
    real(rk8) :: smp(lbc:ubc,1:nlevsoi)    ! soil matrix potential [mm]
    real(rk8) :: sdamp      ! extrapolates soiwat dependence of evaporation
    integer(ik4)  :: pi     ! pft index
    real(rk8) :: temp(lbc:ubc)  ! accumulator for rootr weighting
    ! index of the soil layer right above the water table (-)
    integer(ik4)  :: jwt(lbc:ubc)
    real(rk8) :: smp1,dsmpdw1,wh,wh_zwt,ka
    real(rk8) :: dwat2(lbc:ubc,1:nlevsoi+1)
    ! used in calculating qin, qout (difference in equilbirium matric potential)
    real(rk8) :: dzq
    real(rk8) :: zimm(lbc:ubc,0:nlevsoi)  ! layer interface depth [mm]
    ! equilibrium matric potential for each layer [mm]
    real(rk8) :: zq(lbc:ubc,1:nlevsoi+1)
    ! equilibrium volumetric water content
    real(rk8) :: vol_eq(lbc:ubc,1:nlevsoi+1)
    real(rk8) :: tempi    ! temp variable for calculating vol_eq
    real(rk8) :: temp0    ! temp variable for calculating vol_eq
    real(rk8) :: voleq1   ! temp variable for calculating vol_eq
    real(rk8) :: zwtmm(lbc:ubc)     ! water table depth [mm]
    real(rk8) :: imped(lbc:ubc,1:nlevsoi)
    real(rk8) :: vol_ice(lbc:ubc,1:nlevsoi)
    real(rk8) :: icefrac(lbc:ubc,1:nlevsoi)
    real(rk8) :: vwc_zwt(lbc:ubc)
    real(rk8) :: vwc_liq(lbc:ubc,1:nlevsoi+1) ! liquid volumetric water content

    ! Assign local pointers to derived type members (column-level)

    h2osoi_ice        => clm3%g%l%c%cws%h2osoi_ice
    qcharge           => clm3%g%l%c%cws%qcharge
    hkdepth           => clm3%g%l%c%cps%hkdepth
    zi                => clm3%g%l%c%cps%zi
    zwt               => clm3%g%l%c%cws%zwt
    ctype             => clm3%g%l%c%itype
    npfts             => clm3%g%l%c%npfts
    z                 => clm3%g%l%c%cps%z
    dz                => clm3%g%l%c%cps%dz
    smpmin            => clm3%g%l%c%cps%smpmin
    watsat            => clm3%g%l%c%cps%watsat
    hksat             => clm3%g%l%c%cps%hksat
    bsw               => clm3%g%l%c%cps%bsw
    sucsat            => clm3%g%l%c%cps%sucsat
    eff_porosity      => clm3%g%l%c%cps%eff_porosity
    rootr_col         => clm3%g%l%c%cps%rootr_column
    t_soisno          => clm3%g%l%c%ces%t_soisno
    h2osoi_liq        => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol        => clm3%g%l%c%cws%h2osoi_vol
    qflx_infl         => clm3%g%l%c%cwf%qflx_infl
    fracice           => clm3%g%l%c%cps%fracice
    qflx_tran_veg_col => clm3%g%l%c%cwf%pwf_a%qflx_tran_veg
    pfti              => clm3%g%l%c%pfti
    smp_l             => clm3%g%l%c%cws%smp_l
    hk_l              => clm3%g%l%c%cws%hk_l

    ! Assign local pointers to derived type members (pft-level)

    pactive           => clm3%g%l%c%p%active
    qflx_tran_veg_pft => clm3%g%l%c%p%pwf%qflx_tran_veg
    rootr_pft         => clm3%g%l%c%p%pps%rootr
    pwtcol            => clm3%g%l%c%p%wtcol

    ! Because the depths in this routine are in mm, use local
    ! variable arrays instead of pointers

    amx = 0.0_rk8
    bmx = 0.0_rk8
    cmx = 0.0_rk8
    rmx = 0.0_rk8

    do j = 1, nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        zmm(c,j) = z(c,j)*1.e3_rk8
        dzmm(c,j) = dz(c,j)*1.e3_rk8
        zimm(c,j) = zi(c,j)*1.e3_rk8
        ! calculate icefrac up here
        vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
        icefrac(c,j) = min(1._rk8,vol_ice(c,j)/watsat(c,j))
        vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_rk8)/(dz(c,j)*denh2o)
      end do
    end do

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      zimm(c,0) = 0.0_rk8
      zwtmm(c)  = zwt(c)*1.e3_rk8
    end do

    ! First step is to calculate the column-level effective rooting
    ! fraction in each soil layer. This is done outside the usual
    ! PFT-to-column averaging routines because it is not a simple
    ! weighted average of the PFT level rootr arrays. Instead, the
    ! weighting depends on both the per-unit-area transpiration
    ! of the PFT and the PFTs area relative to all PFTs.

    temp(:) = 0._rk8

    do j = 1 , nlevsoi
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        rootr_col(c,j) = 0._rk8
      end do
    end do

    do pi = 1 , max_pft_per_col
      do j = 1 , nlevsoi
        do fc = 1 , num_hydrologyc
          c = filter_hydrologyc(fc)
          if (pi <= npfts(c)) then
            p = pfti(c) + pi - 1
            if (pactive(p)) then
              rootr_col(c,j) = rootr_col(c,j) + &
                      rootr_pft(p,j) * qflx_tran_veg_pft(p) * pwtcol(p)
            end if
          end if
        end do
      end do
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        if (pi <= npfts(c)) then
          p = pfti(c) + pi - 1
          if (pactive(p)) then
            temp(c) = temp(c) + qflx_tran_veg_pft(p) * pwtcol(p)
          end if
        end if
      end do
    end do

    do j = 1, nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        if ( temp(c) /= 0._rk8 ) then
          rootr_col(c,j) = rootr_col(c,j)/temp(c)
        end if
      end do
    end do

    !compute jwt index
    ! The layer index of the first unsaturated layer, i.e., the
    ! layer right above the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       ! allow jwt to equal zero when zwt is in top layer
       do j = 1,nlevsoi
        if ( zwt(c) <= zi(c,j) ) then
          jwt(c) = j-1
          exit
        end if
      end do

      ! compute vwc at water table depth (mainly for case when t < tfrz)
      !     this will only be used when zwt is below the soil column
      vwc_zwt(c) = watsat(c,nlevsoi)
      if ( t_soisno(c,jwt(c)+1) < tfrz ) then
        vwc_zwt(c) = vwc_liq(c,nlevsoi)
        do j = nlevsoi,nlevgrnd
          if ( zwt(c) <= zi(c,j) ) then
            smp1 = hfus*(tfrz-t_soisno(c,j)) / &
                    (grav*t_soisno(c,j)) * 1000._rk8  !(mm)
            !smp1 = max(0._rk8,smp1)
            smp1 = max(sucsat(c,nlevsoi),smp1)
            vwc_zwt(c) = watsat(c,nlevsoi) * &
                    (smp1/sucsat(c,nlevsoi))**(-1._rk8/bsw(c,nlevsoi))
            ! for temperatures close to tfrz, limit vwc to total water content
            vwc_zwt(c) = min(vwc_zwt(c), &
                    0.5_rk8*(watsat(c,nlevsoi) + h2osoi_vol(c,nlevsoi)) )
            exit
          end if
        end do
      end if
    end do

    ! calculate the equilibrium water content based on the water table depth

    do j = 1 , nlevsoi
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        if ( (zwtmm(c) <= zimm(c,j-1)) ) then
          vol_eq(c,j) = watsat(c,j)

        ! use the weighted average from the saturated part (depth > wtd)
        ! and the equilibrium solution for the rest of the layer

        else if ( (zwtmm(c) < zimm(c,j)) .and. &
                  (zwtmm(c) > zimm(c,j-1)) ) then
          tempi = 1.0_rk8
          temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1)) / &
                  sucsat(c,j)))**(1._rk8-1._rk8/bsw(c,j))
          voleq1 = -sucsat(c,j)*watsat(c,j)/(1._rk8-1._rk8/bsw(c,j)) / &
                  (zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
          vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + &
                  watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
          vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
          vol_eq(c,j) = max(vol_eq(c,j),0.0_rk8)
        else
          tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j)) / &
                     sucsat(c,j)))**(1._rk8-1._rk8/bsw(c,j))
          temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1)) / &
                     sucsat(c,j)))**(1._rk8-1._rk8/bsw(c,j))
          vol_eq(c,j) = -sucsat(c,j)*watsat(c,j) / &
                  (1._rk8-1._rk8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
          vol_eq(c,j) = max(vol_eq(c,j),0.0_rk8)
          vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
        end if
        zq(c,j) = -sucsat(c,j) * &
                (max(vol_eq(c,j)/watsat(c,j),0.01_rk8))**(-bsw(c,j))
        zq(c,j) = max(smpmin(c), zq(c,j))
      end do
    end do

    ! If water table is below soil column calculate zq for the 11th layer
    j = nlevsoi
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      if ( jwt(c) == nlevsoi ) then
        tempi = 1._rk8
        temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j)) / &
                   sucsat(c,j)))**(1._rk8-1._rk8/bsw(c,j))
        vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j) / &
                (1._rk8-1._rk8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
        vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_rk8)
        vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
        zq(c,j+1) = -sucsat(c,j) * &
                (max(vol_eq(c,j+1)/watsat(c,j),0.01_rk8))**(-bsw(c,j))
        zq(c,j+1) = max(smpmin(c), zq(c,j+1))
      end if
    end do

    ! Hydraulic conductivity and soil matric potential and their derivatives

    sdamp = 0._rk8
    do j = 1, nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        ! compute hydraulic conductivity based on liquid water content only

        if ( origflag == 1 ) then
          s1 = 0.5_rk8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
              (0.5_rk8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
        else
          s1 = 0.5_rk8*(vwc_liq(c,j) + vwc_liq(c,min(nlevsoi, j+1))) / &
              (0.5_rk8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
        end if
        s1 = min(1._rk8, s1)
        s2 = hksat(c,j)*s1**(2._rk8*bsw(c,j)+2._rk8)

        ! replace fracice with impedance factor, as in zhao 97,99
        if ( origflag == 1 ) then
          imped(c,j) = (1._rk8 - 0.5_rk8*(fracice(c,j) + &
            fracice(c,min(nlevsoi, j+1))))
        else
          imped(c,j) = 10._rk8**(-e_ice*(0.5_rk8*(icefrac(c,j) + &
                  icefrac(c,min(nlevsoi, j+1)))))
        end if
        hk(c,j) = imped(c,j)*s1*s2
        dhkdw(c,j) = imped(c,j)*(2._rk8*bsw(c,j)+3._rk8)*s2* &
                       (1._rk8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))

        ! compute matric potential and derivative based on
        ! liquid water content only
        if ( origflag == 1 ) then
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_rk8)
        else
          s_node = max(vwc_liq(c,j)/watsat(c,j), 0.01_rk8)
        end if
        s_node = min(1.0_rk8, s_node)

        smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
        smp(c,j) = max(smpmin(c), smp(c,j))

        if ( origflag == 1 ) then
          dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))
        else
          dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/vwc_liq(c,j)
        end if

        smp_l(c,j) = smp(c,j)
        hk_l(c,j) = hk(c,j)
      end do
    end do

    ! aquifer (11th) layer
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      zmm(c,nlevsoi+1) = 0.5_rk8*(1.e3_rk8*zwt(c) + zmm(c,nlevsoi))
      if ( jwt(c) < nlevsoi ) then
        dzmm(c,nlevsoi+1) = dzmm(c,nlevsoi)
      else
        dzmm(c,nlevsoi+1) = (1.e3_rk8*zwt(c) - zmm(c,nlevsoi))
      end if
    end do

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      qin(c,j)    = qflx_infl(c)
      den    = (zmm(c,j+1)-zmm(c,j))
      dzq    = (zq(c,j+1)-zq(c,j))
      num    = (smp(c,j+1)-smp(c,j)) - dzq
      qout(c,j)   = -hk(c,j)*num/den
      dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
      dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
      rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c) * rootr_col(c,j)
      amx(c,j) =  0._rk8
      bmx(c,j) =  dzmm(c,j)*(sdamp+1._rk8/dtsrf) + dqodw1(c,j)
      cmx(c,j) =  dqodw2(c,j)
    end do

    ! Nodes j=2 to j=nlevsoi-1

    do j = 2 , nlevsoi - 1
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        den    = (zmm(c,j) - zmm(c,j-1))
        dzq    = (zq(c,j)-zq(c,j-1))
        num    = (smp(c,j)-smp(c,j-1)) - dzq
        qin(c,j)    = -hk(c,j-1)*num/den
        dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
        dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
        den    = (zmm(c,j+1)-zmm(c,j))
        dzq    = (zq(c,j+1)-zq(c,j))
        num    = (smp(c,j+1)-smp(c,j)) - dzq
        qout(c,j)   = -hk(c,j)*num/den
        dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
        dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
        rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
        amx(c,j) = -dqidw0(c,j)
        bmx(c,j) =  dzmm(c,j)/dtsrf - dqidw1(c,j) + dqodw1(c,j)
        cmx(c,j) =  dqodw2(c,j)
      end do
    end do

    ! Node j=nlevsoi (bottom)

    j = nlevsoi
    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)
      if ( j > jwt(c) ) then !water table is in soil column
        den    = (zmm(c,j) - zmm(c,j-1))
        dzq    = (zq(c,j)-zq(c,j-1))
        num    = (smp(c,j)-smp(c,j-1)) - dzq
        qin(c,j)    = -hk(c,j-1)*num/den
        dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
        dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
        qout(c,j)   =  0._rk8
        dqodw1(c,j) =  0._rk8
        rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
        amx(c,j) = -dqidw0(c,j)
        bmx(c,j) =  dzmm(c,j)/dtsrf - dqidw1(c,j) + dqodw1(c,j)
        cmx(c,j) =  0._rk8

        ! next set up aquifer layer; hydrologically inactive
        rmx(c,j+1) = 0._rk8
        amx(c,j+1) = 0._rk8
        bmx(c,j+1) = dzmm(c,j+1)/dtsrf
        cmx(c,j+1) = 0._rk8
      else ! water table is below soil column

        ! compute aquifer soil moisture as average of layer 10 and saturation
        if ( origflag == 1 ) then
          s_node = max(0.5_rk8*(1.0_rk8+h2osoi_vol(c,j)/watsat(c,j)), 0.01_rk8)
        else
          s_node = max(0.5_rk8*((vwc_zwt(c)+vwc_liq(c,j))/watsat(c,j)), 0.01_rk8)
        end if
        s_node = min(1.0_rk8, s_node)

        ! compute smp for aquifer layer
        smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
        smp1 = max(smpmin(c), smp1)

        ! compute dsmpdw for aquifer layer
        dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

        ! first set up bottom layer of soil column
        den    = (zmm(c,j) - zmm(c,j-1))
        dzq    = (zq(c,j)-zq(c,j-1))
        num    = (smp(c,j)-smp(c,j-1)) - dzq
        qin(c,j)    = -hk(c,j-1)*num/den
        dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
        dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
        den    = (zmm(c,j+1)-zmm(c,j))
        dzq    = (zq(c,j+1)-zq(c,j))
        num    = (smp1-smp(c,j)) - dzq
        qout(c,j)   = -hk(c,j)*num/den
        dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
        dqodw2(c,j) = -( hk(c,j)*dsmpdw1 + num*dhkdw(c,j))/den

        rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
        amx(c,j) = -dqidw0(c,j)
        bmx(c,j) =  dzmm(c,j)/dtsrf - dqidw1(c,j) + dqodw1(c,j)
        cmx(c,j) =  dqodw2(c,j)

        ! next set up aquifer layer; den/num unchanged, qin=qout
        qin(c,j+1)    = qout(c,j)
        dqidw0(c,j+1) = -(-hk(c,j)*dsmpdw(c,j) + num*dhkdw(c,j))/den
        dqidw1(c,j+1) = -( hk(c,j)*dsmpdw1   + num*dhkdw(c,j))/den
        qout(c,j+1)   =  0._rk8  ! zero-flow bottom boundary condition
        dqodw1(c,j+1) =  0._rk8  ! zero-flow bottom boundary condition
        rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
        amx(c,j+1) = -dqidw0(c,j+1)
        bmx(c,j+1) =  dzmm(c,j+1)/dtsrf - dqidw1(c,j+1) + dqodw1(c,j+1)
        cmx(c,j+1) =  0._rk8
      end if
    end do

    ! Solve for dwat

    jtop(:) = 1
    where ( abs(amx) < 1.e-20_rk8 )
      amx = sign(1.0e-20_rk8,amx)
    end where
    where ( abs(bmx) < 1.e-20_rk8 )
      bmx = sign(1.0e-20_rk8,bmx)
    end where
    where ( abs(cmx) < 1.e-20_rk8 )
      cmx = sign(1.0e-20_rk8,cmx)
    end where
    where ( abs(rmx) < 1.e-20_rk8 )
      rmx = sign(1.0e-20_rk8,rmx)
    end where
    call Tridiagonal(lbc, ubc, 1, nlevsoi+1, jtop,      &
                     num_hydrologyc, filter_hydrologyc, &
                     amx, bmx, cmx, rmx, dwat2 )
    ! set dwat
    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)
      do j = 1 , nlevsoi
        dwat(c,j) = dwat2(c,j)
      end do
    end do

    ! Renew the mass of liquid water
    ! also compute qcharge from dwat in aquifer layer
    ! update in drainage for case jwt < nlevsoi

    do fc = 1,num_hydrologyc
      c = filter_hydrologyc(fc)
      do j = 1 , nlevsoi
        h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat2(c,j)*dzmm(c,j)
      end do

      ! calculate qcharge for case jwt < nlevsoi
      if ( jwt(c) < nlevsoi ) then
        !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat
        wh_zwt = 0._rk8

        ! Recharge rate qcharge to groundwater (positive to aquifer)
        s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_rk8)
        s1 = min(1._rk8, s_node)

        !scs: this is the expression for unsaturated hk
        ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
               *s1**(2._rk8*bsw(c,jwt(c)+1)+3._rk8)

        ! Recharge rate qcharge to groundwater (positive to aquifer)
        smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
        wh      = smp1 - zq(c,max(1,jwt(c)))

        !scs: original formulation
        if ( jwt(c) == 0 ) then
          qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._rk8)
        else
          ! qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._rk8)
          !scs: 1/2, assuming flux is at zwt interface, saturation
          ! deeper than zwt
          qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._rk8*2.0)
        end if

        ! To limit qcharge  (for the first several timesteps)
        qcharge(c) = max(-10.0_rk8/dtsrf,qcharge(c))
        qcharge(c) = min( 10.0_rk8/dtsrf,qcharge(c))
      else
        ! if water table is below soil column, compute qcharge from dwat2(11)
        qcharge(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtsrf
      end if
    end do
  end subroutine SoilWater
  !
  ! Calculate subsurface drainage
  !
  subroutine Drainage(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                      num_urbanc, filter_urbanc, icefrac)
    implicit none
    integer(ik4) , intent(in) :: lbc, ubc                     ! column bounds
    integer(ik4) , intent(in) :: num_hydrologyc               ! number of column soil points in column filter
    integer(ik4) , intent(in) :: num_urbanc                   ! number of column urban points in column filter
    integer(ik4) , intent(in) :: filter_urbanc(ubc-lbc+1)     ! column filter for urban points
    integer(ik4) , intent(in) :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    real(rk8), intent(in) :: icefrac(lbc:ubc,1:nlevsoi)   ! fraction of ice in layer

    real(rk8), pointer :: h2osfc(:)         ! surface water (mm)
    real(rk8), pointer :: frac_h2osfc(:)    !
    real(rk8), pointer :: topo_ndx(:)       ! topographic index
    real(rk8), pointer :: topo_slope(:)     ! topographic slope
    real(rk8), pointer :: frost_table(:)    ! frost table depth (m)
    real(rk8), pointer :: zwt_perched(:)    ! perched water table depth (m)
    ! perched wt sub-surface runoff (mm H2O /s)
    real(rk8), pointer :: qflx_drain_perched(:)
    integer(ik4) , pointer :: ltype(:)      ! landunit type
    integer(ik4) , pointer :: clandunit(:)  ! column's landunit
    integer(ik4) , pointer :: ctype(:)      ! column type index
    integer(ik4) , pointer :: snl(:)        ! number of snow layers
    ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(rk8), pointer :: qflx_snwcp_liq(:)
    ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(rk8), pointer :: qflx_snwcp_ice(:)
    ! ground surface dew formation (mm H2O /s) [+]
    real(rk8), pointer :: qflx_dew_grnd(:)
    ! surface dew added to snow pack (mm H2O /s) [+]
    real(rk8), pointer :: qflx_dew_snow(:)
    ! sublimation rate from snow pack (mm H2O /s) [+]
    real(rk8), pointer :: qflx_sub_snow(:)
    real(rk8), pointer :: dz(:,:)           ! layer depth (m)
    real(rk8), pointer :: bsw(:,:)          ! Clapp and Hornberger "b"
    ! effective porosity = porosity - vol_ice
    real(rk8), pointer :: eff_porosity(:,:)
    real(rk8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)
    ! hydraulic conductivity at saturation (mm H2O /s)
    real(rk8), pointer :: hksat(:,:)
    real(rk8), pointer :: sucsat(:,:)       ! minimum soil suction (mm)
    real(rk8), pointer :: z(:,:)            ! layer depth (m)
    ! interface level below a "z" level (m)
    real(rk8), pointer :: zi(:,:)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: watsat(:,:)
    real(rk8), pointer :: hkdepth(:)    ! decay factor (m)
    real(rk8), pointer :: zwt(:)        ! water table depth (m)
    real(rk8), pointer :: wa(:)         ! water in the unconfined aquifer (mm)
    real(rk8), pointer :: qcharge(:)    ! aquifer recharge rate (mm/s)
#if (defined VICHYDRO)
    real(rk8), pointer :: moist(:,:)    !soil layer moisture (mm)
    real(rk8), pointer :: ice(:,:)      !soil layer moisture (mm)
    !fracton of Dsmax where non-linear baseflow begins
    real(rk8), pointer :: Ds(:)
    real(rk8), pointer :: Dsmax(:)      !max. velocity of baseflow (mm/day)
    !fraction of maximum soil moisutre where non-liear base flow occurs
    real(rk8), pointer :: Wsvic(:)
    real(rk8), pointer :: c_param(:)     !baseflow exponent (Qb)
    real(rk8), pointer :: max_moist(:,:) !maximum soil moisture (ice + liq)
    real(rk8), pointer :: depth(:,:)     !VIC soil depth
    real(rk8), pointer :: hk_l(:,:)      !hydraulic conductivity (mm/s)
    character(len=32) :: subname = 'Drainage'  ! subroutine name
#endif

    real(rk8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(rk8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)

    ! sub-surface runoff (mm H2O /s)
    real(rk8), pointer :: qflx_drain(:)
    ! irrigation flux (mm H2O /s)
    real(rk8), pointer :: qflx_irrig(:)
    ! qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
    real(rk8), pointer :: qflx_qrgwl(:)
    ! implicit evaporation for soil temperature equation
    real(rk8), pointer :: eflx_impsoil(:)
    ! soil saturation excess [mm h2o/s]
    real(rk8), pointer :: qflx_rsub_sat(:)

    integer(ik4)  :: c,j,fc,i   ! indices
    ! water needed to bring soil moisture to watmin (mm)
    real(rk8) :: xs(lbc:ubc)
    ! layer thickness (mm)
    real(rk8) :: dzmm(lbc:ubc,1:nlevsoi)
    ! index of the soil layer right above the water table (-)
    integer(ik4)  :: jwt(lbc:ubc)
    ! subsurface runoff - bottom drainage (mm/s)
    real(rk8) :: rsub_bot(lbc:ubc)
    ! subsurface runoff - topographic control (mm/s)
    real(rk8) :: rsub_top(lbc:ubc)
    ! decay factor (m-1)
    real(rk8) :: fff(lbc:ubc)
    ! excess soil water above saturation at layer i (mm)
    real(rk8) :: xsi(lbc:ubc)
    ! excess soil water above saturation at layer 1 (mm)
    real(rk8) :: xs1(lbc:ubc)
    ! summation of hk*dzmm for layers below water table (mm**2/s)
    real(rk8) :: wtsub
    real(rk8) :: rous      ! aquifer yield (-)
    ! summation of dzmm of layers below water table (mm)
    real(rk8) :: dzsum
    ! summation of icefrac*dzmm of layers below water table (-)
    real(rk8) :: icefracsum
    ! fractional impermeability of soil layers (-)
    real(rk8) :: fracice_rsub(lbc:ubc)
    ! available soil liquid water in a layer
    real(rk8) :: available_h2osoi_liq
    real(rk8) :: rsub_top_max
    real(rk8) :: h2osoi_vol
    real(rk8) :: imped
    real(rk8) :: rsub_top_tot
    real(rk8) :: rsub_top_layer
    real(rk8) :: qcharge_tot
    real(rk8) :: qcharge_layer
    real(rk8) :: s_y
    integer(ik4)  :: k,k_frz,k_perch
    real(rk8) :: sat_lev
    real(rk8) :: s1
    real(rk8) :: s2
    real(rk8) :: m
    real(rk8) :: b
    real(rk8) :: q_perch
    real(rk8) :: q_perch_max
#if (defined VICHYDRO)
    ! temporary variable for ARNO subsurface runoff calculation
    real(rk8) :: dsmax_tmp(lbc:ubc)
    ! temporary variable for ARNO subsurface runoff calculation
    real(rk8) :: rsub_tmp
    ! temporary variable for ARNO subsurface runoff calculation
    real(rk8) :: frac
    ! relative moisture, temporary variable
    real(rk8) :: rel_moist
    ! summation of hk*dzmm for layers in the third VIC layer
    real(rk8) :: wtsub_vic
#endif
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    h2osfc         => clm3%g%l%c%cws%h2osfc
    frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc
    topo_ndx       => clm3%g%l%c%cps%topo_ndx
    topo_slope     => clm3%g%l%c%cps%topo_slope
    frost_table    => clm3%g%l%c%cws%frost_table
    zwt_perched    => clm3%g%l%c%cws%zwt_perched
    qflx_drain_perched    => clm3%g%l%c%cwf%qflx_drain_perched
    clandunit      => clm3%g%l%c%landunit
    ltype          => clm3%g%l%itype
    ctype          => clm3%g%l%c%itype
    snl           => clm3%g%l%c%cps%snl
    dz            => clm3%g%l%c%cps%dz
    bsw           => clm3%g%l%c%cps%bsw
    t_soisno      => clm3%g%l%c%ces%t_soisno
    hksat         => clm3%g%l%c%cps%hksat
    sucsat        => clm3%g%l%c%cps%sucsat
    z             => clm3%g%l%c%cps%z
    zi            => clm3%g%l%c%cps%zi
    watsat        => clm3%g%l%c%cps%watsat
    hkdepth       => clm3%g%l%c%cps%hkdepth
    zwt           => clm3%g%l%c%cws%zwt
    wa            => clm3%g%l%c%cws%wa
    qcharge       => clm3%g%l%c%cws%qcharge
    eff_porosity  => clm3%g%l%c%cps%eff_porosity
    qflx_snwcp_liq => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_liq
    qflx_snwcp_ice => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    qflx_dew_grnd => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_dew_snow => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_drain    => clm3%g%l%c%cwf%qflx_drain
    qflx_irrig    => clm3%g%l%c%cwf%qflx_irrig
    qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_rsub_sat => clm3%g%l%c%cwf%qflx_rsub_sat
    eflx_impsoil  => clm3%g%l%c%cef%eflx_impsoil
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
#if (defined VICHYDRO)
    Dsmax          => clm3%g%l%c%cps%dsmax
    Ds             => clm3%g%l%c%cps%ds
    Wsvic          => clm3%g%l%c%cps%Wsvic
    c_param        => clm3%g%l%c%cps%c_param
    max_moist      => clm3%g%l%c%cps%max_moist
    depth          => clm3%g%l%c%cps%depth
    moist          => clm3%g%l%c%cws%moist
    ice            => clm3%g%l%c%cws%ice
    hk_l          => clm3%g%l%c%cws%hk_l
#endif

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        dzmm(c,j) = dz(c,j)*1.e3_rk8
      end do
    end do

    ! Initial set

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      qflx_drain(c)    = 0._rk8
      rsub_bot(c)      = 0._rk8
      qflx_rsub_sat(c) = 0._rk8
      rsub_top(c)      = 0._rk8
      fracice_rsub(c)  = 0._rk8
    end do

    ! The layer index of the first unsaturated layer, i.e., the layer
    ! right above the water table

    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      jwt(c) = nlevsoi
      ! allow jwt to equal zero when zwt is in top layer
      do j = 1 , nlevsoi
        if ( zwt(c) <= zi(c,j) ) then
          jwt(c) = j-1
          exit
        end if
      end do
    end do

    rous = 0.2_rk8

    ! QCHARGE =========================================
    ! Water table changes due to qcharge
    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)

      ! use analytical expression for aquifer specific yield
      rous = watsat(c,nlevsoi) * ( 1.0_rk8 - &
              (1.0_rk8+1.e3_rk8*zwt(c)/sucsat(c,nlevsoi))**(-1.0_rk8/bsw(c,nlevsoi)))
      rous = max(rous,0.02_rk8)

       ! water table is below the soil column
      if ( jwt(c) == nlevsoi ) then
        wa(c)  = wa(c) + qcharge(c)  * dtsrf
        zwt(c) = zwt(c) - (qcharge(c)  * dtsrf)/1000._rk8/rous
      else
        ! water table within soil layers 1-9
        ! try to raise water table to account for qcharge
        qcharge_tot = qcharge(c) * dtsrf
        if ( qcharge_tot > 0.0_rk8 ) then !rising water table
          do j = jwt(c)+1 , 1 , -1
            ! use analytical expression for specific yield
            s_y = watsat(c,j) * ( 1.0_rk8 - &
                    (1.0_rk8+1.e3_rk8*zwt(c)/sucsat(c,j))**(-1.0_rk8/bsw(c,j)))
            s_y = max(s_y,0.02_rk8)
            qcharge_layer = min(qcharge_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3_rk8))
            qcharge_layer = max(qcharge_layer,0._rk8)

            if ( s_y > 0._rk8 ) zwt(c) = zwt(c) - qcharge_layer/s_y/1000._rk8

            qcharge_tot = qcharge_tot - qcharge_layer
            if ( qcharge_tot <= 0.0_rk8 ) exit
          end do
        else ! deepening water table (negative qcharge)
          do j = jwt(c)+1 , nlevsoi
            ! use analytical expression for specific yield
            s_y = watsat(c,j) * ( 1.0_rk8 - &
                    (1.0_rk8 + 1.e3_rk8*zwt(c)/sucsat(c,j))**(-1.0_rk8/bsw(c,j)))
            s_y = max(s_y,0.02_rk8)
            qcharge_layer = max(qcharge_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3_rk8))
            qcharge_layer = min(qcharge_layer,0._rk8)
            qcharge_tot = qcharge_tot - qcharge_layer

            if ( qcharge_tot >= 0.0_rk8 ) then
              zwt(c) = zwt(c) - qcharge_layer/s_y/1000._rk8
              exit
            else
              zwt(c) = zi(c,j)
            end if
          end do
          if ( qcharge_tot > 0.0_rk8 ) zwt(c) = zwt(c) - qcharge_tot/1000._rk8/rous
        end if

        ! recompute jwt for following calculations
        ! allow jwt to equal zero when zwt is in top layer
        jwt(c) = nlevsoi
        do j = 1 ,nlevsoi
          if ( zwt(c) <= zi(c,j) ) then
            jwt(c) = j-1
            exit
          end if
        end do
      end if
    end do

    !==  BASEFLOW ==================================================
    ! perched water table code
    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)

      !  specify maximum drainage rate
      q_perch_max = 1.e-5_rk8 * sin(topo_slope(c) * (rpi/180._rk8))

      ! if layer containing water table is frozen, compute the following:
      ! frost table, perched water table, and drainage from perched
      ! saturated layer

      ! define frost table as first frozen layer with unfrozen layer above it
      if (t_soisno(c,1) > tfrz ) then
        k_frz = nlevsoi
      else
        k_frz=1
      end if

      do k = 2 , nlevsoi
        if ( t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz ) then
          k_frz = k
          exit
        end if
      end do

      frost_table(c)=z(c,k_frz)

      ! initialize perched water table to frost table, and
      ! qflx_drain_perched(c) to zero
      zwt_perched(c)=frost_table(c)
      qflx_drain_perched(c) = 0._rk8

      ! water table above frost table  =============================
      ! if water table is above frost table, do not use topmodel
      ! baseflow formulation
      if ( zwt(c) < frost_table(c) .and. &
           t_soisno(c,k_frz) <= tfrz .and. origflag == 0 ) then
        ! compute drainage from perched saturated region
        wtsub = 0._rk8
        q_perch = 0._rk8
        do k = jwt(c)+1 , k_frz
          imped = 10._rk8**(-e_ice*(0.5_rk8*(icefrac(c,k) + &
                  icefrac(c,min(nlevsoi, k+1)))))
          q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
          wtsub = wtsub + dzmm(c,k)
        end do
        if ( wtsub > 0._rk8 ) q_perch = q_perch/wtsub

        qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

        ! remove drainage from perched saturated layers
        rsub_top_tot = - qflx_drain_perched(c) * dtsrf
        do k = jwt(c)+1 , k_frz
          rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
          rsub_top_layer=min(rsub_top_layer,0._rk8)
          rsub_top_tot = rsub_top_tot - rsub_top_layer

          h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

          if ( rsub_top_tot >= 0.0_rk8 ) then
            zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._rk8
            exit
          else
            zwt(c) = zi(c,k)
          end if
        end do

        ! if rsub_top_tot is greater than available water (above frost table),
        ! then decrease qflx_drain_perched by residual amount for water balance
        qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtsrf

        ! recompute jwt
        ! allow jwt to equal zero when zwt is in top layer
        jwt(c) = nlevsoi
        do j = 1 , nlevsoi
          if ( zwt(c) <= zi(c,j) ) then
            jwt(c) = j-1
            exit
          end if
        end do
      else
        ! water table below frost table  =============================
        ! compute possible perched water table *and* groundwater table
        ! afterwards locate perched water table from bottom up starting
        ! at frost table sat_lev is an arbitrary saturation level used
        ! to determine perched water table
        sat_lev=0.9

        k_perch = 1
        do k = k_frz,1,-1
          h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) + &
                       h2osoi_ice(c,k)/(dz(c,k)*denice)
          if ( h2osoi_vol/watsat(c,k) <= sat_lev ) then
            k_perch = k
            exit
          end if
        end do

        ! if frost_table = nlevsoi, only compute perched water table if frozen
        if ( t_soisno(c,k_frz) > tfrz ) k_perch=k_frz

        ! if perched water table exists
        if ( k_frz > k_perch ) then
          ! interpolate between k_perch and k_perch+1 to find perched
          ! water table height
          s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) + &
                h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice)) / &
                watsat(c,k_perch)
          s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) + &
                h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice)) / &
                watsat(c,k_perch+1)

          m = (z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
          b = z(c,k_perch+1)-m*s2
          zwt_perched(c)=max(0._rk8,m*sat_lev+b)

          ! compute drainage from perched saturated region
          wtsub = 0._rk8
          q_perch = 0._rk8
          do k = k_perch, k_frz
            imped = 10._rk8**(-e_ice*(0.5_rk8*(icefrac(c,k) + &
                                 icefrac(c,min(nlevsoi, k+1)))))
            q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
            wtsub = wtsub + dzmm(c,k)
          end do
          if ( wtsub > 0._rk8 ) q_perch = q_perch/wtsub

          qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt_perched(c))

          ! no perched water table drainage if using original formulation
          if ( origflag == 1 ) qflx_drain_perched(c) = 0._rk8

          ! remove drainage from perched saturated layers
          rsub_top_tot = -  qflx_drain_perched(c) * dtsrf
          do k = k_perch+1 , k_frz
            rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
            rsub_top_layer=min(rsub_top_layer,0._rk8)
            rsub_top_tot = rsub_top_tot - rsub_top_layer
            h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

            if ( rsub_top_tot >= 0.0_rk8 ) then
              zwt_perched(c) = zwt_perched(c) - &
                      rsub_top_layer/eff_porosity(c,k)/1000._rk8
              exit
            else
              zwt_perched(c) = zi(c,k)
            end if
          end do

          ! if rsub_top_tot is greater than available water (above frost
          ! table), then decrease qflx_drain_perched by residual amount
          ! for water balance
          qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtsrf

        else
          qflx_drain_perched(c) = 0._rk8
        end if !k_frz > k_perch

        ! Topographic runoff
        fff(c)         = 1._rk8/ hkdepth(c)
        dzsum = 0._rk8
        icefracsum = 0._rk8
        do j = max(jwt(c),1), nlevsoi
          dzsum  = dzsum + dzmm(c,j)
          icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
        end do
        ! add ice impedance factor to baseflow
        if ( origflag == 1 ) then
#if (defined VICHYDRO)
          call fatal(__FILE__,__LINE__, &
                  subname // ':: VICHYDRO is not available for origflag = 1')
#else
          fracice_rsub(c) = max(0._rk8,exp(-3._rk8*(1._rk8-(icefracsum/dzsum))) &
               - exp(-3._rk8))/(1.0_rk8-exp(-3._rk8))
          imped = (1._rk8 - fracice_rsub(c))
          rsub_top_max = 5.5e-3_rk8
#endif
        else
#if (defined VICHYDRO)
          imped = 10._rk8**(-e_ice*min(1.0_rk8,ice(c,nlayer)/max_moist(c,nlayer)))
          dsmax_tmp(c) = Dsmax(c) * dtsrf/ secspday !mm/day->mm/dtsrf
          rsub_top_max = dsmax_tmp(c)
#else
          imped = 10._rk8**(-e_ice*(icefracsum/dzsum))
          rsub_top_max = 10._rk8 * sin((rpi/180.) * topo_slope(c))
#endif
        end if
#if (defined VICHYDRO)
        ! ARNO model for the bottom soil layer (based on bottom soil layer
        ! moisture from previous time step
        ! use watmin instead for resid_moist to be consistent with
        ! default hydrology
        rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin)
        frac = (Ds(c) * rsub_top_max )/Wsvic(c)
        rsub_tmp = (frac * rel_moist)/dtsrf
        if ( rel_moist > Wsvic(c) ) then
          frac = (rel_moist - Wsvic(c))/(1.0_rk8 - Wsvic(c))
          rsub_tmp = rsub_tmp + (rsub_top_max * &
            (1.0_rk8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtsrf
        end if
        rsub_top(c) = imped * rsub_tmp
        ! make sure baseflow isn't negative
        rsub_top(c) = max(0._rk8, rsub_top(c))
#else
        rsub_top(c) = imped * rsub_top_max* exp(-fff(c)*zwt(c))
#endif
        ! use analytical expression for aquifer specific yield
        rous = watsat(c,nlevsoi) * ( 1.0_rk8 - &
          (1.0_rk8+1.e3_rk8*zwt(c)/sucsat(c,nlevsoi))**(-1.0_rk8/bsw(c,nlevsoi)))
        rous = max(rous,0.02_rk8)

        ! water table is below the soil column
        if ( jwt(c) == nlevsoi ) then
          wa(c)  = wa(c) - rsub_top(c) * dtsrf
          zwt(c) = zwt(c) + (rsub_top(c) * dtsrf)/1000._rk8/rous
          h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) + &
                  max(0._rk8,(wa(c)-5000._rk8))
          wa(c)  = min(wa(c), 5000._rk8)
        else
          ! water table within soil layers 1-9
          !--  Now remove water via rsub_top
          rsub_top_tot = - rsub_top(c) * dtsrf
          !should never be positive... but include for completeness
          if ( rsub_top_tot > 0.) then !rising water table
            write(stderr,*) 'RSUB_TOP IS POSITIVE in Drainage!'
            write(stderr,*)'clm model is stopping'
            call fatal(__FILE__,__LINE__,'clm now stopping')


          else ! deepening water table
#if (defined VICHYDRO)
            wtsub_vic = 0._rk8
            do j = (nlvic(1)+nlvic(2)+1), nlevsoi
              wtsub_vic = wtsub_vic + hk_l(c,j)*dzmm(c,j)
            end do

            do j = (nlvic(1)+nlvic(2)+1) , nlevsoi
              rsub_top_layer = max(rsub_top_tot, &
                      rsub_top_tot*hk_l(c,j)*dzmm(c,j)/wtsub_vic)
              rsub_top_layer = min(rsub_top_layer,0._rk8)
              h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
              rsub_top_tot = rsub_top_tot - rsub_top_layer
            end do
#else
            do j = jwt(c)+1 , nlevsoi
              ! use analytical expression for specific yield
              s_y = watsat(c,j) * ( 1.0_rk8 - &
                (1.0_rk8+1.e3_rk8*zwt(c)/sucsat(c,j))**(-1.0_rk8/bsw(c,j)))
              s_y=max(s_y,0.02_rk8)

              rsub_top_layer = max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3_rk8))
              rsub_top_layer = min(rsub_top_layer,0._rk8)
              h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer

              rsub_top_tot = rsub_top_tot - rsub_top_layer

              if ( rsub_top_tot >= 0.) then
                zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._rk8
                exit
              else
                zwt(c) = zi(c,j)
              end if
            end do
#endif
            ! remove residual rsub_top
            zwt(c) = zwt(c) - rsub_top_tot/1000._rk8/rous
            wa(c) = wa(c) + rsub_top_tot
          end if

          ! recompute jwt
          ! allow jwt to equal zero when zwt is in top layer
          jwt(c) = nlevsoi
          do j = 1,nlevsoi
            if ( zwt(c) <= zi(c,j) ) then
              jwt(c) = j-1
              exit
            end if
          end do
        end if ! end of jwt if construct

        zwt(c) = max(0.0_rk8,zwt(c))
        zwt(c) = min(80._rk8,zwt(c))
      end if
    end do


    ! excessive water above saturation added to the above unsaturated
    ! layer like a bucket if column fully saturated, excess water goes
    ! to runoff

    do j = nlevsoi , 2 , -1
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        xsi(c) = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._rk8)
        h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
        h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
      end do
    end do

    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)
      xs1(c) = max(max(h2osoi_liq(c,1),0._rk8) - &
                   max(0._rk8,(pondmx+watsat(c,1)*dzmm(c,1) - &
                       h2osoi_ice(c,1))),0._rk8)
      h2osoi_liq(c,1) = min(max(0._rk8,pondmx+watsat(c,1)*dzmm(c,1) - &
              h2osoi_ice(c,1)), h2osoi_liq(c,1))
      if ( ltype(clandunit(c)) == isturb ) then
        qflx_rsub_sat(c)     = xs1(c) / dtsrf
      else
        if ( h2osfcflag == 1 ) then
          ! send this water up to h2osfc rather than sending to drainage
          h2osfc(c) = h2osfc(c) + xs1(c)
          qflx_rsub_sat(c)     = 0._rk8
        else
          ! use original code to send water to drainage (non-h2osfc case)
          qflx_rsub_sat(c)     = xs1(c) / dtsrf
        end if
      end if
      ! add in ice check
      xs1(c) = max(max(h2osoi_ice(c,1),0._rk8) - &
                   max(0._rk8,(pondmx+watsat(c,1)*dzmm(c,1) - &
                   h2osoi_liq(c,1))),0._rk8)
      h2osoi_ice(c,1) = min(max(0._rk8,pondmx+watsat(c,1)*dzmm(c,1) - &
                               h2osoi_liq(c,1)), h2osoi_ice(c,1))
      qflx_snwcp_ice(c) = qflx_snwcp_ice(c) + xs1(c) / dtsrf
    end do

    ! Limit h2osoi_liq to be greater than or equal to watmin.
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water

    do j = 1 , nlevsoi-1
      do fc = 1 , num_hydrologyc
        c = filter_hydrologyc(fc)
        if ( h2osoi_liq(c,j) < watmin ) then
          xs(c) = watmin - h2osoi_liq(c,j)
          ! deepen water table if water is passed from below zwt layer
          if ( j == jwt(c) ) then
            zwt(c) = zwt(c) + xs(c)/eff_porosity(c,j)/1000._rk8
          endif
        else
          xs(c) = 0._rk8
        end if
        h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
        h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
      end do
    end do

    ! Get water for bottom layer from layers above if possible
    j = nlevsoi
    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)
      if ( h2osoi_liq(c,j) < watmin ) then
        xs(c) = watmin - h2osoi_liq(c,j)
        searchforwater: &
        do i = nlevsoi-1, 1, -1
          available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._rk8)
          if ( available_h2osoi_liq >= xs(c) ) then
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
            h2osoi_liq(c,i) = h2osoi_liq(c,i) - xs(c)
            xs(c) = 0._rk8
            exit searchforwater
          else
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + available_h2osoi_liq
            h2osoi_liq(c,i) = h2osoi_liq(c,i) - available_h2osoi_liq
            xs(c) = xs(c) - available_h2osoi_liq
          end if
        end do searchforwater
      else
        xs(c) = 0._rk8
      end if
      ! Needed in case there is no water to be found
      h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
      ! Instead of removing water from aquifer where it eventually
      ! shows up as excess drainage to the ocean, take it back out of
      ! drainage
      rsub_top(c) = rsub_top(c) - xs(c)/dtsrf
    end do

    do fc = 1 , num_hydrologyc
      c = filter_hydrologyc(fc)

      ! Sub-surface runoff and drainage

      qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

      ! Set imbalance for snow capping

      qflx_qrgwl(c) = qflx_snwcp_liq(c)

      ! Implicit evaporation term is now zero

      eflx_impsoil(c) = 0._rk8

      ! Renew the ice and liquid mass due to condensation

      if ( snl(c)+1 >= 1 ) then
        ! make consistent with how evap_grnd removed in infiltration
        h2osoi_liq(c,1) = h2osoi_liq(c,1) + &
                (1._rk8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtsrf
        h2osoi_ice(c,1) = h2osoi_ice(c,1) + &
                (1._rk8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtsrf
        if (qflx_sub_snow(c)*dtsrf > h2osoi_ice(c,1)) then
          qflx_sub_snow(c) = h2osoi_ice(c,1)/dtsrf
          h2osoi_ice(c,1) = 0._rk8
        else
          h2osoi_ice(c,1) = h2osoi_ice(c,1) - &
                  (1._rk8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtsrf
        end if
      end if
    end do

    ! No drainage for urban columns
    ! (except for pervious road as computed above)

    do fc = 1 , num_urbanc
      c = filter_urbanc(fc)
      if ( ctype(c) /= icol_road_perv ) then
        qflx_drain(c) = 0._rk8
        qflx_irrig(c) = 0._rk8
        ! This must be done for roofs and impervious road (walls will be zero)
        qflx_qrgwl(c) = qflx_snwcp_liq(c)
        eflx_impsoil(c) = 0._rk8
      end if

      ! Renew the ice and liquid mass due to condensation for urban
      ! roof and impervious road

      if ( ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv ) then
        if ( snl(c)+1 >= 1 ) then
          h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtsrf
          h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtsrf)
          if (qflx_sub_snow(c)*dtsrf > h2osoi_ice(c,1)) then
            qflx_sub_snow(c) = h2osoi_ice(c,1)/dtsrf
            h2osoi_ice(c,1) = 0._rk8
          else
            h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtsrf)
          end if
        end if
      end if
    end do
  end subroutine Drainage

end module mod_clm_soilhydrology
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
