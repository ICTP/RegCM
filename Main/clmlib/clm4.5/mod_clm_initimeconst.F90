module mod_clm_initimeconst
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_runparams
  use mod_dynparam
  use mod_mpmessage
  use mod_mppparam
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_nchelper
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varpar , only : nlevsoi , nlevgrnd , nlevlak , numpft
  use mod_clm_varpar , only : numrad , nlevurb , mach_eps , toplev_equalspace
  use mod_clm_varpar , only : nlev_equalspace , more_vertlayers , nlevsoifl
  use mod_clm_varcon , only : albsat , albdry , spval , secspday , rpi
  use mod_clm_varcon , only : istice , istdlak , istwet , isturb , istsoil
  use mod_clm_varcon , only : istcrop , icol_roof , icol_sunwall
  use mod_clm_varcon , only : icol_shadewall , icol_road_perv
  use mod_clm_varcon , only : icol_road_imperv, zlak , dzlak
  use mod_clm_varcon , only : zsoi , dzsoi , zisoi , dzsoi_decomp
  use mod_clm_varcon , only : pc , mu
  use mod_clm_varctl , only : fsurdat , ialblawr
  use mod_clm_varctl , only : fsnowoptics, fsnowaging
  use mod_clm_varsur , only : pctspec
  use mod_clm_pftvarcon
  use mod_clm_organicfile , only : organicrd
#ifdef CN
  use mod_clm_cninispecial , only : CNiniSpecial
#endif
#if (defined VICHYDRO)
  use mod_clm_vicmap , only : initCLMVICMap
  use mod_clm_initsoilparvic , only : initSoilParVIC
  use mod_clm_varpar , only : nlayer , nlayert
  use mod_clm_varcon , only : nlvic
#endif
  use mod_clm_snicar , only : SnowAge_init, SnowOptics_init
#if (defined LCH4)
  use mod_clm_ch4varcon , only : usephfact, fin_use_fsat
#endif
#ifdef CN
#ifndef CENTURY_DECOMP
  use mod_clm_cndecompcascadebgc , only : init_decompcascade
#else
  use mod_clm_cndecompcascadecentury , only : init_decompcascade
#endif
#endif
  use mod_clm_soilhydrology , only : h2osfcflag

  implicit none

  private

  save

  public :: iniTimeConst

  contains
  !
  ! Initialize time invariant clm variables
  ! 1) removed references to shallow lake - since it is not used
  ! 2) ***Make c%z, c%zi and c%dz allocatable depending on if you
  !    have lake or soil
  ! 3) rootfr only initialized for soil points
  !
  subroutine iniTimeConst
    implicit none
    ! gridcell elevation standard deviation
    real(rk8), pointer :: topo_std(:)
    ! gridcell topographic slope
    real(rk8), pointer :: topo_slope(:)
    ! microtopography pdf sigma (m)
    real(rk8), pointer :: micro_sigma(:)
    ! level at which h2osfc "percolates"
    real(rk8), pointer :: h2osfc_thresh(:)
    ! mineral hksat
    real(rk8), pointer :: hksat_min(:,:)
    ! SCA shape parameter
    real(rk8), pointer :: n_melt(:)
    integer(ik4) , pointer :: ivt(:)       !  vegetation type index
    integer(ik4) , pointer :: pcolumn(:)   ! column index of corresponding pft
    integer(ik4) , pointer :: pgridcell(:) ! gridcell index of corresponding pft
    integer(ik4) , pointer :: clandunit(:) ! landunit index of column
    integer(ik4) , pointer :: cgridcell(:) ! gridcell index of column
    integer(ik4) , pointer :: ctype(:)     ! column type index
    integer(ik4) , pointer :: ltype(:)     ! landunit type index
    real(rk8), pointer :: thick_wall(:)    ! total thickness of urban wall
    real(rk8), pointer :: thick_roof(:)    ! total thickness of urban roof
    real(rk8), pointer :: lat(:)           ! gridcell latitude (radians)

    real(rk8), pointer :: z(:,:)      ! layer depth (m)
    real(rk8), pointer :: zi(:,:)     ! interface level below a "z" level (m)
    real(rk8), pointer :: dz(:,:)     ! layer thickness depth (m)
    real(rk8), pointer :: rootfr(:,:) ! fraction of roots in each soil layer
    ! fraction of roots in each soil layer for urban pervious road
    real(rk8), pointer :: rootfr_road_perv(:,:)
    !root resistance by layer (0-1)  (nlevgrnd)
    real(rk8), pointer :: rresis(:,:)
    real(rk8), pointer :: dewmx(:)   ! maximum allowed dew [mm]
    real(rk8), pointer :: bsw(:,:)   ! Clapp and Hornberger "b" (nlevgrnd)
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8), pointer :: watsat(:,:)
    ! volumetric soil water at field capacity (nlevsoi)
    real(rk8), pointer :: watfc(:,:)
    real(rk8), pointer :: watdry(:,:)  ! btran parameter for btran=0
    real(rk8), pointer :: watopt(:,:)  ! btran parameter for btran = 1
    ! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(rk8), pointer :: hksat(:,:)
    ! minimum soil suction (mm) (nlevgrnd)
    real(rk8), pointer :: sucsat(:,:)
    ! heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)
    real(rk8), pointer :: csol(:,:)
    ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd)
    real(rk8), pointer :: tkmg(:,:)
    ! thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
    real(rk8), pointer :: tkdry(:,:)
    ! thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
    real(rk8), pointer :: tksatu(:,:)
    ! maximum saturated fraction for a gridcell
    real(rk8), pointer :: wtfact(:)
#ifdef LCH4
    ! coefficient for determining finundated (m)
    !real(rk8), pointer :: zwt0(:)
    ! maximum inundated fractional area for gridcell
    !real(rk8), pointer :: f0(:)
    ! coefficient for determining finundated (m)
    !real(rk8), pointer :: p3(:)
    real(rk8), pointer :: k(:)
    real(rk8), pointer :: q(:)
    real(rk8), pointer :: v(:)
    real(rk8), pointer :: maxf(:)
#endif
    ! bulk density of dry soil material [kg/m^3]
    real(rk8), pointer :: bd(:,:)
    ! added by Lei Meng
#ifdef LCH4
    real(rk8), pointer :: pH(:)      ! pH values for methane code
#endif
    ! restriction for min of soil potential (mm) (new)
    real(rk8), pointer :: smpmin(:)
    ! decay factor (m)
    real(rk8), pointer :: hkdepth(:)
    integer(ik4) , pointer :: isoicol(:)  ! soil color class
#ifdef CN
    real(rk8) , pointer :: q10(:)
    real(rk8) , pointer :: ndep(:)
#endif

    ! added by F. Li and S. Levis
    real(rk8), pointer :: gdp_lf(:)      ! global gdp data
    real(rk8), pointer :: peatf_lf(:)    ! global peatf data
    integer(ik4), pointer :: abm_lf(:)   ! global abm data

    ! threshold soil moisture based on clay content
    real(rk8), pointer :: gwc_thr(:)
    ! [frc] Mass fraction clay limited to 0.20
    real(rk8), pointer :: mss_frc_cly_vld(:)
    ! emission factors for isoprene (ug isoprene m-2 h-1)
    real(rk8), pointer :: efisop(:,:)
    real(rk8), pointer :: max_dayl(:) ! maximum daylength (s)
    real(rk8), pointer :: sandfrac(:)
    real(rk8), pointer :: clayfrac(:)
#if (defined VICHYDRO)
    !CLM column level sand fraction for calculating VIC parameters
    real(rk8), pointer :: sandcol(:,:)
    !CLM column level clay fraction for calculating VIC parameters
    real(rk8), pointer :: claycol(:,:)
    !CLM column level organic matter fraction for calculating VIC parameters
    real(rk8), pointer :: om_fraccol(:,:)
    !b infiltration parameter
    real(rk8), pointer :: b_infil(:)
    real(rk8), pointer :: dsmax(:)   !maximum baseflow rate
    !fracton of Dsmax where non-linear baseflow begins
    real(rk8), pointer :: ds(:)
    !fraction of maximum soil moisutre where non-liear base flow occurs
    real(rk8), pointer :: Wsvic(:)
    !pore-size distribution related paramter(Q12)
    real(rk8), pointer :: expt(:,:)
    !Saturated hydrologic conductivity
    real(rk8), pointer :: ksat(:,:)
    !soil moisture dissusion parameter
    real(rk8), pointer :: phi_s(:,:)
    !layer depth of upper layer(m)
    real(rk8), pointer :: depth(:,:)
    real(rk8), pointer :: porosity(:,:)  !soil porosity
    real(rk8), pointer :: max_moist(:,:) !maximum soil moisture (ice + liq)
#endif
    ! For lakes
    real(rk8), pointer :: cellsand(:,:)    ! column 3D sand
    real(rk8), pointer :: cellclay(:,:)    ! column 3D clay
    real(rk8), pointer :: cellorg(:,:)     ! column 3D org content
    real(rk8), pointer :: lakedepth(:)     ! variable lake depth
    ! extinction coefficient from surface data (1/m)
    real(rk8), pointer :: etal(:)
    real(rk8), pointer :: lakefetch(:)  ! lake fetch from surface data (m)

    type(clm_filetype)  :: ncid   ! netcdf id
    integer(ik4)  :: j , ib , lev ! indices
    ! integer(ik4)  :: bottom
    integer(ik4)  :: g , l , c , p             ! indices
    integer(ik4)  :: m                         ! vegetation type index
    real(rk8) :: tkm     ! mineral conductivity
    real(rk8) :: xksat   ! maximum hydraulic conductivity of soil [mm/s]
    real(rk8) :: thick_equal = 0.2_rk8
    real(rk8), pointer :: zsoifl(:)   ! original soil midpoint
    real(rk8), pointer :: zisoifl(:)  ! original soil interface depth
    real(rk8), pointer :: dzsoifl(:)  ! original soil thickness
    real(rk8) :: clay,sand        ! temporaries
    ! real(rk8) :: slope,intercept  ! temporary, for rooting distribution
    real(rk8) :: temp, max_decl   ! temporary, for calculation of max_dayl
    integer(ik4)  :: begp, endp  ! per-proc beginning and ending pft indices
    integer(ik4)  :: begc, endc  ! per-proc beginning and ending column indices
    integer(ik4)  :: begl, endl  ! per-proc beginning and ending ldunit indices
    integer(ik4)  :: begg, endg  ! per-proc gridcell ending gridcell indices
    integer(ik4)  :: numg    ! total number of gridcells across all processors
    integer(ik4)  :: numl    ! total number of landunits across all processors
    integer(ik4)  :: numc    ! total number of columns across all processors
    integer(ik4)  :: nump    ! total number of pfts across all processors
#if (defined VICHYDRO)
    integer(ik4)  :: ivic , ivicstrt , ivicend  ! indices
    real(rk8) , pointer :: b2d(:)         ! read in - VIC b
    real(rk8) , pointer :: ds2d(:)        ! read in - VIC Ds
    real(rk8) , pointer :: dsmax2d(:)     ! read in - VIC Dsmax
    real(rk8) , pointer :: ws2d(:)        ! read in - VIC Ws
#endif

    real(rk8) , pointer :: temp_ef(:)     ! read in - temporary EFs
    real(rk8) , pointer :: efisop2d(:,:)  ! read in - isoprene emission factors

    integer(ik4) ,pointer :: soic2d(:)    ! read in - soil color

    ! added by F. Li and S. Levis
    real(rk8) , pointer :: gdp(:)         ! global gdp data
    real(rk8) , pointer :: peatf(:)       ! global peatf data
    integer(ik4) , pointer :: abm(:)      ! global abm data

    real(rk8) , pointer :: sand3d(:,:)    ! read in - soil texture: percent sand
    real(rk8) , pointer :: clay3d(:,:)    ! read in - soil texture: percent clay
    real(rk8) , pointer :: organic3d(:,:) ! read in - organic matter: kg/m3
    real(rk8) , pointer :: gti(:)         ! read in - fmax
#ifdef LCH4
    !real(rk8) , pointer :: zwt0_in(:)     ! read in - zwt0
    !real(rk8) , pointer :: f0_in(:)       ! read in - f0
    !real(rk8) , pointer :: p3_in(:)       ! read in - p3
    real(rk8) , pointer :: k_in(:)
    real(rk8) , pointer :: q_in(:)
    real(rk8) , pointer :: v_in(:)
    real(rk8) , pointer :: maxf_in(:)
    real(rk8) , pointer :: pH_in(:)       ! read in - pH
#endif

#ifdef CN
    real(rk8) , pointer :: q10_in(:)
    real(rk8) , pointer :: ndep_in(:)
#endif

    real(rk8) , pointer :: lakedepth_in(:) ! read in - lakedepth
    real(rk8) , pointer :: etal_in(:)     ! read in - etal
    real(rk8) , pointer :: lakefetch_in(:) ! read in - lakefetch
    real(rk8) :: om_frac                 ! organic matter fraction
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(rk8) :: om_tkm = 0.25_rk8
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(rk8) :: om_csol = 2.5_rk8
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(rk8) :: om_tkd = 0.05_rk8
    real(rk8) :: om_watsat  ! porosity of organic soil
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(rk8) :: om_hksat
    ! saturated suction for organic matter (mm)(Letts, 2000)
    real(rk8) :: om_sucsat
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(rk8) :: om_b
    ! depth (m) that organic matter takes on characteristics of sapric peat
    real(rk8) :: zsapric = 0.5_rk8
    ! organic matter (kg/m3) where soil is assumed to act like peat
    real(rk8) :: organic_max = 130._rk8
    ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(rk8) :: csol_bedrock = 2.0e6_rk8
    real(rk8) :: pcalpha      = 0.5_rk8   ! percolation threshold
    real(rk8) :: pcbeta       = 0.139_rk8 ! percolation exponent
    real(rk8) :: perc_frac   ! "percolating" fraction of organic soil
    real(rk8) :: perc_norm   ! normalize to 1 when 100% organic soil
    real(rk8) :: uncon_hksat ! series conductivity of mineral/organic soil
    real(rk8) :: uncon_frac  ! fraction of "unconnected" soil
    integer(ik4)  :: ier                           ! error status
    character(len= 32) :: subname = 'iniTimeConst' ! subroutine name
    integer(ik4) :: mxsoil_color       ! maximum number of soil color classes
    real(rk8), allocatable :: zurb_wall(:,:)  ! wall (layer node depth)
    real(rk8), allocatable :: zurb_roof(:,:)  ! roof (layer node depth)
    real(rk8), allocatable :: dzurb_wall(:,:) ! wall (layer thickness)
    real(rk8), allocatable :: dzurb_roof(:,:) ! roof (layer thickness)
    real(rk8), allocatable :: ziurb_wall(:,:) ! wall (layer interface)
    real(rk8), allocatable :: ziurb_roof(:,:) ! roof (layer interface)

    integer(ik4) :: nzero_slope  ! Number of points to zero out slope

    real(rk8) , pointer :: std(:)     ! read in - topo_std
    real(rk8) , pointer :: tslope(:)  ! read in - topo_slope
    real(rk8) :: maxslope , slopemax , minslope , d , fd
    real(rk8) :: dfdd , slope0 , slopebeta

#ifdef __PGI
    !real(rk8) , external :: erf
#endif
    if ( myid == italk ) then
      write(stdout,*) 'Attempting to initialize time invariant variables'
    end if

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! added by F. Li and S. Levis
    allocate(gdp(begg:endg))
    allocate(peatf(begg:endg))
    allocate(abm(begg:endg))

    allocate(soic2d(begg:endg), gti(begg:endg))
#ifdef LCH4
    !allocate(zwt0_in(begg:endg))
    !allocate(f0_in(begg:endg))
    !allocate(p3_in(begg:endg))
    allocate(k_in(begg:endg))
    allocate(q_in(begg:endg))
    allocate(v_in(begg:endg))
    allocate(maxf_in(begg:endg))

    if ( usephfact ) allocate(ph_in(begg:endg))
#endif
    allocate(lakedepth_in(begg:endg))
    allocate(etal_in(begg:endg))
    allocate(lakefetch_in(begg:endg))
#ifdef CN
    allocate(q10_in(begg:endg))
    allocate(ndep_in(begg:endg))
#endif

    allocate(temp_ef(begg:endg),efisop2d(6,begg:endg))
#if (defined VICHYDRO)
    allocate(b2d(begg:endg), ds2d(begg:endg), &
             dsmax2d(begg:endg),ws2d(begg:endg))
    allocate(sandcol(begc:endc,1:nlevgrnd), &
             claycol(begc:endc,1:nlevgrnd), &
             om_fraccol(begc:endc,1:nlevgrnd)) ! allocation for local variables
#endif

    efisop          => clm3%g%gve%efisop

    ! Assign local pointers to derived subtypes components (gridcell-level)
    lat             => clm3%g%lat

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype               => clm3%g%l%itype
    thick_wall          => clm3%g%l%lps%thick_wall
    thick_roof          => clm3%g%l%lps%thick_roof

    ! Assign local pointers to derived subtypes components (column-level)

    topo_std        => clm3%g%l%c%cps%topo_std
    topo_slope      => clm3%g%l%c%cps%topo_slope
    micro_sigma     => clm3%g%l%c%cps%micro_sigma
    h2osfc_thresh   => clm3%g%l%c%cps%h2osfc_thresh
    hksat_min       => clm3%g%l%c%cps%hksat_min
    n_melt          => clm3%g%l%c%cps%n_melt
    ctype           => clm3%g%l%c%itype
    clandunit       => clm3%g%l%c%landunit
    cgridcell       => clm3%g%l%c%gridcell
    z               => clm3%g%l%c%cps%z
    dz              => clm3%g%l%c%cps%dz
    zi              => clm3%g%l%c%cps%zi
    bsw             => clm3%g%l%c%cps%bsw
    watsat          => clm3%g%l%c%cps%watsat
    watfc           => clm3%g%l%c%cps%watfc
    watdry          => clm3%g%l%c%cps%watdry
    watopt          => clm3%g%l%c%cps%watopt
    rootfr_road_perv => clm3%g%l%c%cps%rootfr_road_perv
    hksat           => clm3%g%l%c%cps%hksat
    sucsat          => clm3%g%l%c%cps%sucsat
    tkmg            => clm3%g%l%c%cps%tkmg
    tksatu          => clm3%g%l%c%cps%tksatu
    tkdry           => clm3%g%l%c%cps%tkdry
    csol            => clm3%g%l%c%cps%csol
    smpmin          => clm3%g%l%c%cps%smpmin
    hkdepth         => clm3%g%l%c%cps%hkdepth
    wtfact          => clm3%g%l%c%cps%wtfact
    bd              => clm3%g%l%c%cps%bd
#ifdef LCH4
    !zwt0            => clm3%g%l%c%cps%zwt0
    !f0              => clm3%g%l%c%cps%f0
    !p3              => clm3%g%l%c%cps%p3
    k               => clm3%g%l%c%cps%k
    q               => clm3%g%l%c%cps%q
    v               => clm3%g%l%c%cps%v
    maxf            => clm3%g%l%c%cps%maxf
    pH              => clm3%g%l%c%cps%pH
#endif
    isoicol         => clm3%g%l%c%cps%isoicol

#ifdef CN
    q10            => clm3%g%l%c%cps%q10
    ndep           => clm3%g%l%c%cps%ndep
#endif

    ! added by F. Li and S. Levis
    gdp_lf          => clm3%g%l%c%cps%gdp_lf
    peatf_lf        => clm3%g%l%c%cps%peatf_lf
    abm_lf          => clm3%g%l%c%cps%abm_lf

    gwc_thr         => clm3%g%l%c%cps%gwc_thr
    mss_frc_cly_vld => clm3%g%l%c%cps%mss_frc_cly_vld
    max_dayl        => clm3%g%l%c%cps%max_dayl
    cellsand        => clm3%g%l%c%cps%cellsand
    cellclay        => clm3%g%l%c%cps%cellclay
    cellorg         => clm3%g%l%c%cps%cellorg
    lakedepth       => clm3%g%l%c%cps%lakedepth
    etal            => clm3%g%l%c%cps%etal
    lakefetch       => clm3%g%l%c%cps%lakefetch
#if (defined VICHYDRO)
    b_infil        => clm3%g%l%c%cps%b_infil
    dsmax          => clm3%g%l%c%cps%dsmax
    ds             => clm3%g%l%c%cps%ds
    Wsvic          => clm3%g%l%c%cps%Wsvic
    expt           => clm3%g%l%c%cps%expt
    ksat           => clm3%g%l%c%cps%ksat
    phi_s          => clm3%g%l%c%cps%phi_s
    depth          => clm3%g%l%c%cps%depth
    porosity       => clm3%g%l%c%cps%porosity
    max_moist      => clm3%g%l%c%cps%max_moist
#endif

    ! Assign local pointers to derived subtypes components (pft-level)

    ivt             => clm3%g%l%c%p%itype
    pgridcell       => clm3%g%l%c%p%gridcell
    pcolumn         => clm3%g%l%c%p%column
    dewmx           => clm3%g%l%c%p%pps%dewmx
    rootfr          => clm3%g%l%c%p%pps%rootfr
    rresis          => clm3%g%l%c%p%pps%rresis
    sandfrac        => clm3%g%l%c%p%pps%sandfrac
    clayfrac        => clm3%g%l%c%p%pps%clayfrac

    if (nlevurb > 0) then
      allocate(zurb_wall(begl:endl,nlevurb),    &
               zurb_roof(begl:endl,nlevurb),    &
               dzurb_wall(begl:endl,nlevurb),   &
               dzurb_roof(begl:endl,nlevurb),   &
               ziurb_wall(begl:endl,0:nlevurb), &
               ziurb_roof(begl:endl,0:nlevurb), &
               stat=ier)
      if (ier /= 0) then
        call fatal(__FILE__,__LINE__, 'iniTimeConst: allocation error')
      end if
    end if

    ! --------------------------------------------------------------------
    ! Read soil color, sand and clay from surface dataset
    ! --------------------------------------------------------------------

    if (myid == italk) then
      write(stdout,*) 'Attempting to read soil color, sand &
              &and clay boundary data .....'
    end if

    call clm_openfile(fsurdat,ncid)

    call clm_inqdim(ncid,'nlevsoi',nlevsoifl)
    if ( .not. more_vertlayers ) then
       if ( nlevsoifl /= nlevsoi ) then
         call fatal(__FILE__,__LINE__, &
                 trim(subname)//' ERROR: Number of soil layers on &
                 &file does NOT match the number being used' )
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if
    allocate(sand3d(begg:endg,nlevsoifl), clay3d(begg:endg,nlevsoifl))
    allocate(organic3d(begg:endg,nlevsoifl))

    ! Determine number of soil color classes
    ! if number of soil color classes is not on input dataset set it to 8

    if ( .not. clm_check_var(ncid,'mxsoil_color') ) then
      mxsoil_color = 8
    else
      call clm_readvar(ncid,'mxsoil_color',mxsoil_color)
    end if

    ! Methane code parameters for finundated
#ifdef LCH4
    if ( .not. fin_use_fsat ) then
      !call clm_readvar(ncid,'ZWT0',zwt0_in,gcomm_gridcell)
      !call clm_readvar(ncid,'F0',f0_in,gcomm_gridcell)
      !call clm_readvar(ncid,'P3',p3_in,gcomm_gridcell)
      call clm_readvar(ncid,'K_PAR',k_in,gcomm_gridcell)
      call clm_readvar(ncid,'XM_PAR',q_in,gcomm_gridcell)
      call clm_readvar(ncid,'V_PAR',v_in,gcomm_gridcell)
      call clm_readvar(ncid,'MAXF',maxf_in,gcomm_gridcell)
    end if
    ! pH factor for methane model
    if ( usephfact ) then
      call clm_readvar(ncid,'PH',ph_in,gcomm_gridcell)
    end if
#endif
    ! def CH4

#ifdef CN
    call clm_readvar(ncid,'Q10',q10_in,gcomm_gridcell)
    call clm_readvar(ncid,'NDEP',ndep_in,gcomm_gridcell)
#endif

    ! Read lakedepth
    if ( .not. clm_check_var(ncid,'LAKEDEPTH') ) then
      if ( myid == italk ) then
        write(stdout,*) 'WARNING:: LAKEDEPTH not found on surface data set.'
        write(stdout,*) 'All lake columns will have lake depth', &
                        ' set equal to default value.'
      end if
      lakedepth_in(:) = 10.0_rk8
    else
      call clm_readvar(ncid,'LAKEDEPTH',lakedepth_in,gcomm_gridcell)
    end if

    ! Read lake eta
    if ( .not. clm_check_var(ncid,'ETALAKE') ) then
      if ( myid == italk ) then
        write(stdout,*) 'WARNING:: ETALAKE not found on surface data set.'
        write(stdout,*) 'All lake columns will have eta', &
                        ' set equal to default value'
      end if
      etal_in(:) = -1.0_rk8
    else
      call clm_readvar(ncid,'ETALAKE',etal_in,gcomm_gridcell)
    end if

    ! lake fetch
    if ( .not. clm_check_var(ncid,'LAKEFETCH') ) then
      if ( myid == italk ) then
        write(stdout,*) 'WARNING:: LAKEFETCH not found on surface data set.'
        write(stdout,*) 'All lake columns will have fetch', &
                        ' set equal to default value'
      end if
      lakefetch_in(:) = -1.0_rk8
    else
      call clm_readvar(ncid,'LAKEFETCH',lakefetch_in,gcomm_gridcell)
    end if

    ! Read in topographic index and slope
    allocate(tslope(begg:endg))
    allocate(std(begg:endg))

    call clm_readvar(ncid,'SLOPE',tslope,gcomm_gridcell)
    call clm_readvar(ncid,'STD_ELEV',std,gcomm_gridcell)

    ! Read fmax

    call clm_readvar(ncid,'FMAX',gti,gcomm_gridcell)
#if (defined VICHYDRO)
    call clm_readvar(ncid,'binfl',b2d,gcomm_gridcell)
    call clm_readvar(ncid,'Ds',ds2d,gcomm_gridcell)
    call clm_readvar(ncid,'Dsmax',dsmax2d,gcomm_gridcell)
    call clm_readvar(ncid,'Ws',ws2d,gcomm_gridcell)
#endif

    ! Read in soil color, sand and clay fraction
    call clm_readvar(ncid,'SOIL_COLOR',soic2d,gcomm_gridcell)

    ! Read in GDP data added by F. Li and S. Levis
    call clm_readvar(ncid,'gdp',gdp,gcomm_gridcell)

    ! Read in peatf data added by F. Li and S. Levis
    call clm_readvar(ncid,'peatf',peatf,gcomm_gridcell)

    ! Read in ABM data added by F. Li and S. Levis
    call clm_readvar(ncid,'abm',abm,gcomm_gridcell)

    ! Read in emission factors
    call clm_readvar(ncid,'EF1_BTR',temp_ef,gcomm_gridcell)
    efisop2d(1,:)=temp_ef(:)

    call clm_readvar(ncid,'EF1_FET',temp_ef,gcomm_gridcell)
    efisop2d(2,:)=temp_ef(:)

    call clm_readvar(ncid,'EF1_FDT',temp_ef,gcomm_gridcell)
    efisop2d(3,:)=temp_ef(:)

    call clm_readvar(ncid,'EF1_SHR',temp_ef,gcomm_gridcell)
    efisop2d(4,:)=temp_ef(:)

    call clm_readvar(ncid,'EF1_GRS',temp_ef,gcomm_gridcell)
    efisop2d(5,:)=temp_ef(:)

    call clm_readvar(ncid,'EF1_CRP',temp_ef,gcomm_gridcell)
    efisop2d(6,:)=temp_ef(:)

    call clm_readvar(ncid,'PCT_SAND',sand3d,gcomm_gridcell)

    call clm_readvar(ncid,'PCT_CLAY',clay3d,gcomm_gridcell)

    call clm_closefile(ncid)

    if (myid == italk) then
      write(stdout,*) &
              'Successfully read fmax, soil color, sand and clay boundary data'
      write(stdout,*)
    end if

    ! Determine saturated and dry soil albedos for n color classes and
    ! numrad wavebands (1=vis, 2=nir)

    allocate(albsat(mxsoil_color,numrad), &
             albdry(mxsoil_color,numrad), stat=ier)
    if (ier /= 0) then
      write(stderr,*) 'iniTimeConst: allocation error for albsat, albdry'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    if ( mxsoil_color == 8 ) then
      albsat(1:8,1) = &
              [0.12_rk8,0.11_rk8,0.10_rk8,0.09_rk8,0.08_rk8,0.07_rk8,0.06_rk8,0.05_rk8]
      albsat(1:8,2) = &
              [0.24_rk8,0.22_rk8,0.20_rk8,0.18_rk8,0.16_rk8,0.14_rk8,0.12_rk8,0.10_rk8]
      albdry(1:8,1) = &
              [0.24_rk8,0.22_rk8,0.20_rk8,0.18_rk8,0.16_rk8,0.14_rk8,0.12_rk8,0.10_rk8]
      albdry(1:8,2) = &
              [0.48_rk8,0.44_rk8,0.40_rk8,0.36_rk8,0.32_rk8,0.28_rk8,0.24_rk8,0.20_rk8]
    else if ( mxsoil_color == 20 ) then
      if ( ialblawr == 1 ) then
        albsat(1:20,1) = &
                [0.26_rk8,0.24_rk8,0.22_rk8,0.20_rk8,0.19_rk8,0.18_rk8,0.17_rk8,0.16_rk8,&
                  0.15_rk8,0.14_rk8,0.13_rk8,0.12_rk8,0.11_rk8,0.10_rk8,0.09_rk8,0.08_rk8,&
                  0.07_rk8,0.06_rk8,0.05_rk8,0.04_rk8]
        albsat(1:20,2) = &
                [0.52_rk8,0.48_rk8,0.44_rk8,0.40_rk8,0.38_rk8,0.36_rk8,0.34_rk8,0.32_rk8,&
                  0.30_rk8,0.28_rk8,0.26_rk8,0.24_rk8,0.22_rk8,0.20_rk8,0.18_rk8,0.16_rk8,&
                  0.14_rk8,0.12_rk8,0.10_rk8,0.08_rk8]
        albdry(1:20,1) = &
                [0.37_rk8,0.35_rk8,0.33_rk8,0.31_rk8,0.30_rk8,0.29_rk8,0.28_rk8,0.27_rk8,&
                  0.26_rk8,0.25_rk8,0.24_rk8,0.23_rk8,0.22_rk8,0.21_rk8,0.20_rk8,0.19_rk8,&
                  0.18_rk8,0.17_rk8,0.16_rk8,0.15_rk8]
        albdry(1:20,2) = &
                [0.63_rk8,0.59_rk8,0.55_rk8,0.51_rk8,0.49_rk8,0.47_rk8,0.45_rk8,0.43_rk8,&
                  0.41_rk8,0.39_rk8,0.37_rk8,0.35_rk8,0.33_rk8,0.31_rk8,0.29_rk8,0.27_rk8,&
                  0.25_rk8,0.23_rk8,0.21_rk8,0.19_rk8]
      else
        albsat(1:20,1) = &
                [0.25_rk8,0.23_rk8,0.21_rk8,0.20_rk8,0.19_rk8,0.18_rk8,0.17_rk8,0.16_rk8,&
                  0.15_rk8,0.14_rk8,0.13_rk8,0.12_rk8,0.11_rk8,0.10_rk8,0.09_rk8,0.08_rk8,&
                  0.07_rk8,0.06_rk8,0.05_rk8,0.04_rk8]
        albsat(1:20,2) = &
                [0.50_rk8,0.46_rk8,0.42_rk8,0.40_rk8,0.38_rk8,0.36_rk8,0.34_rk8,0.32_rk8,&
                  0.30_rk8,0.28_rk8,0.26_rk8,0.24_rk8,0.22_rk8,0.20_rk8,0.18_rk8,0.16_rk8,&
                  0.14_rk8,0.12_rk8,0.10_rk8,0.08_rk8]
        albdry(1:20,1) = &
                [0.36_rk8,0.34_rk8,0.32_rk8,0.31_rk8,0.30_rk8,0.29_rk8,0.28_rk8,0.27_rk8,&
                  0.26_rk8,0.25_rk8,0.24_rk8,0.23_rk8,0.22_rk8,0.20_rk8,0.18_rk8,0.16_rk8,&
                  0.14_rk8,0.12_rk8,0.10_rk8,0.08_rk8]
        albdry(1:20,2) = &
                [0.61_rk8,0.57_rk8,0.53_rk8,0.51_rk8,0.49_rk8,0.48_rk8,0.45_rk8,0.43_rk8,&
                  0.41_rk8,0.39_rk8,0.37_rk8,0.35_rk8,0.33_rk8,0.31_rk8,0.29_rk8,0.27_rk8,&
                  0.25_rk8,0.23_rk8,0.21_rk8,0.16_rk8]
      end if
    else
      write(stderr,*)'maximum color class = ',mxsoil_color,' is not supported'
      call fatal(__FILE__,__LINE__,'clm_now stopping')
    end if

    do p = begp , endp
      g = pgridcell(p)
      if ( sand3d(g,1)+clay3d(g,1) == 0.0_rk8 )then
        if ( any( sand3d(g,:)+clay3d(g,:) /= 0.0_rk8 ) )then
          call fatal(__FILE__,__LINE__, &
              'found depth points that do NOT sum to zero when surface does' )
        end if
        sand3d(g,:) = 1.0_rk8
        clay3d(g,:) = 1.0_rk8
      end if
      if ( any( sand3d(g,:)+clay3d(g,:) == 0.0_rk8 ) )then
        call fatal(__FILE__,__LINE__, &
                'after setting, found points sum to zero' )
      end if
      sandfrac(p) = sand3d(g,1)/100.0_rk8
      clayfrac(p) = clay3d(g,1)/100.0_rk8
    end do

    ! --------------------------------------------------------------------
    ! If a organic matter dataset has been specified, read it
    ! --------------------------------------------------------------------

    call organicrd(organic3d)

    ! --------------------------------------------------------------------
    ! Initialize time constant arrays of ecophysiological constants and
    ! arrays of dgvm ecophysiological constants
    ! --------------------------------------------------------------------

    do m = 0 , numpft
      if (m <= ntree) then
        pftcon%tree(m) = 1
      else
        pftcon%tree(m) = 0
      end if
      pftcon%z0mr(m) = z0mr(m)
      pftcon%displar(m) = displar(m)
      pftcon%dleaf(m) = dleaf(m)
      pftcon%xl(m) = xl(m)
      do ib = 1 , numrad
        pftcon%rhol(m,ib) = rhol(m,ib)
        pftcon%rhos(m,ib) = rhos(m,ib)
        pftcon%taul(m,ib) = taul(m,ib)
        pftcon%taus(m,ib) = taus(m,ib)
      end do
      pftcon%c3psn(m) = c3psn(m)
      pftcon%slatop(m) = slatop(m)
      pftcon%dsladlai(m) = dsladlai(m)
      pftcon%leafcn(m) = leafcn(m)
      pftcon%flnr(m) = flnr(m)
      pftcon%smpso(m) = smpso(m)
      pftcon%smpsc(m) = smpsc(m)
      pftcon%fnitr(m) = fnitr(m)
      pftcon%woody(m) = woody(m)
      pftcon%lflitcn(m) = lflitcn(m)
      pftcon%frootcn(m) = frootcn(m)
      pftcon%livewdcn(m) = livewdcn(m)
      pftcon%deadwdcn(m) = deadwdcn(m)
      pftcon%graincn(m) = graincn(m)
      pftcon%froot_leaf(m) = froot_leaf(m)
      pftcon%stem_leaf(m) = stem_leaf(m)
      pftcon%croot_stem(m) = croot_stem(m)
      pftcon%flivewd(m) = flivewd(m)
      pftcon%fcur(m) = fcur(m)
      pftcon%lf_flab(m) = lf_flab(m)
      pftcon%lf_fcel(m) = lf_fcel(m)
      pftcon%lf_flig(m) = lf_flig(m)
      pftcon%fr_flab(m) = fr_flab(m)
      pftcon%fr_fcel(m) = fr_fcel(m)
      pftcon%fr_flig(m) = fr_flig(m)
      pftcon%leaf_long(m) = leaf_long(m)
      pftcon%evergreen(m) = evergreen(m)
      pftcon%stress_decid(m) = stress_decid(m)
      pftcon%season_decid(m) = season_decid(m)
      pftcon%dwood(m) = dwood
      pftcon%fertnitro(m) = fertnitro(m)
      pftcon%fleafcn(m)   = fleafcn(m)
      pftcon%ffrootcn(m)  = ffrootcn(m)
      pftcon%fstemcn(m)   = fstemcn(m)
    end do

#ifdef CNDV
    do m = 0,numpft
      dgv_pftcon%crownarea_max(m) = pftpar20(m)
      dgv_pftcon%tcmin(m) = pftpar28(m)
      dgv_pftcon%tcmax(m) = pftpar29(m)
      dgv_pftcon%gddmin(m) = pftpar30(m)
      dgv_pftcon%twmax(m) = pftpar31(m)
      dgv_pftcon%reinickerp(m) = reinickerp
      dgv_pftcon%allom1(m) = allom1
      dgv_pftcon%allom2(m) = allom2
      dgv_pftcon%allom3(m) = allom3
      ! modification for shrubs by X.D.Z
      if (m > ntree .and. m <= nbrdlf_dcd_brl_shrub ) then
        dgv_pftcon%allom1(m) = allom1s
        dgv_pftcon%allom2(m) = allom2s
      end if
    end do
#endif

    ! --------------------------------------------------------------------
    ! Define layer structure for soil, lakes, urban walls and roof
    ! Vertical profile of snow is not initialized here
    ! --------------------------------------------------------------------

    ! Soil layers and interfaces (assumed same for all non-lake patches)
    ! "0" refers to soil surface and "nlevsoi" refers to the bottom of
    ! model soil

    if ( more_vertlayers )then
      ! replace standard exponential grid with a grid that starts out
      ! exponential, then has several evenly spaced layers, then finishes
      ! off exponential.
      ! this allows the upper soil to behave as standard, but then
      ! continues with higher resolution to a deeper depth, so that,
      ! for example, permafrost dynamics are not lost due to an inability
      ! to resolve temperature, moisture, and biogeochemical dynamics
      ! at the base of the active layer
      do j = 1 , toplev_equalspace
        zsoi(j) = scalez*(exp(0.5_rk8*(j-0.5_rk8))-1._rk8)    !node depths
      end do

      do j = toplev_equalspace+1 , toplev_equalspace + nlev_equalspace
        zsoi(j) = zsoi(j-1) + thick_equal
      end do

      do j = toplev_equalspace + nlev_equalspace +1 , nlevgrnd
        zsoi(j) = scalez*(exp(0.5_rk8*((j - nlev_equalspace)-0.5_rk8))-1._rk8) + &
                nlev_equalspace * thick_equal
      end do
    else
      do j = 1 , nlevgrnd
        zsoi(j) = scalez*(exp(0.5_rk8*(j-0.5_rk8))-1._rk8)    !node depths
      end do
    end if

    dzsoi(1) = 0.5_rk8*(zsoi(1)+zsoi(2))  !thickness b/n two interfaces
    do j = 2 , nlevgrnd-1
      dzsoi(j)= 0.5_rk8*(zsoi(j+1)-zsoi(j-1))
    end do
    dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)

    zisoi(0) = 0._rk8
    do j = 1, nlevgrnd-1
      zisoi(j) = 0.5_rk8*(zsoi(j)+zsoi(j+1))  !interface depths
    end do
    zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_rk8*dzsoi(nlevgrnd)

    if (myid == italk) call vprntv(zsoi,size(zsoi),'zsoi')
    if (myid == italk) call vprntv(zisoi,size(zisoi),'zisoi')
    if (myid == italk) call vprntv(dzsoi,size(dzsoi),'dzsoi')

#if (defined VICHYDRO)
    !define the depth of VIC soil layers here
    nlvic(1) = 3
    nlvic(2) = 3
    nlvic(3) = nlevsoi-(nlvic(1)+nlvic(2))
#endif

    ! define a vertical grid spacing such that it is the normal dzsoi
    ! if nlevdecomp =nlevgrnd, or else 1 meter
#ifdef VERTSOILC
    !thickness b/n two interfaces
    dzsoi_decomp(1) = 0.5_rk8*(zsoi(1)+zsoi(2))
    do j = 2 , nlevgrnd-1
      dzsoi_decomp(j)= 0.5_rk8*(zsoi(j+1)-zsoi(j-1))
    end do
    dzsoi_decomp(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)
#else
    dzsoi_decomp(1) = 1.0_rk8
#endif
    if (myid == italk) write(stdout, *) 'dzsoi_decomp', dzsoi_decomp(:)

    ! get original soil depths to be used in interpolation of sand and clay
    allocate(zsoifl(1:nlevsoifl),zisoifl(0:nlevsoifl),dzsoifl(1:nlevsoifl))
    do j = 1 , nlevsoifl
      zsoifl(j) = 0.025*(exp(0.5_rk8*(j-0.5_rk8))-1._rk8)    !node depths
    end do

    dzsoifl(1) = 0.5_rk8*(zsoifl(1)+zsoifl(2)) !thickness b/n two interfaces
    do j = 2 , nlevsoifl-1
      dzsoifl(j)= 0.5_rk8*(zsoifl(j+1)-zsoifl(j-1))
    end do
    dzsoifl(nlevsoifl) = zsoifl(nlevsoifl)-zsoifl(nlevsoifl-1)

    zisoifl(0) = 0._rk8
    do j = 1 , nlevsoifl-1
      zisoifl(j) = 0.5_rk8*(zsoifl(j)+zsoifl(j+1))         !interface depths
    end do
    zisoifl(nlevsoifl) = zsoifl(nlevsoifl) + 0.5_rk8*dzsoifl(nlevsoifl)

    ! Column level initialization for urban wall and roof layers and interfaces
    do l = begl , endl
      ! "0" refers to urban wall/roof surface and "nlevsoi" refers
      ! to urban wall/roof bottom
      if ( ltype(l) == isturb ) then
#if (defined VANCOUVER)
        zurb_wall(l,1) = 0.010_rk8/2._rk8
        zurb_wall(l,2) = zurb_wall(l,1) + 0.010_rk8/2._rk8 + 0.020_rk8/2._rk8
        zurb_wall(l,3) = zurb_wall(l,2) + 0.020_rk8/2._rk8 + 0.070_rk8/2._rk8
        zurb_wall(l,4) = zurb_wall(l,3) + 0.070_rk8/2._rk8 + 0.070_rk8/2._rk8
        zurb_wall(l,5) = zurb_wall(l,4) + 0.070_rk8/2._rk8 + 0.030_rk8/2._rk8

        zurb_roof(l,1) = 0.010_rk8/2._rk8
        zurb_roof(l,2) = zurb_roof(l,1) + 0.010_rk8/2._rk8 + 0.010_rk8/2._rk8
        zurb_roof(l,3) = zurb_roof(l,2) + 0.010_rk8/2._rk8 + 0.010_rk8/2._rk8
        zurb_roof(l,4) = zurb_roof(l,3) + 0.010_rk8/2._rk8 + 0.010_rk8/2._rk8
        zurb_roof(l,5) = zurb_roof(l,4) + 0.010_rk8/2._rk8 + 0.030_rk8/2._rk8

        dzurb_wall(l,1) = 0.010_rk8
        dzurb_wall(l,2) = 0.020_rk8
        dzurb_wall(l,3) = 0.070_rk8
        dzurb_wall(l,4) = 0.070_rk8
        dzurb_wall(l,5) = 0.030_rk8
        write(stdout,*)'Total thickness of wall: ',sum(dzurb_wall(l,:))
        write(stdout,*)'Wall layer thicknesses: ',dzurb_wall(l,:)

        dzurb_roof(l,1) = 0.010_rk8
        dzurb_roof(l,2) = 0.010_rk8
        dzurb_roof(l,3) = 0.010_rk8
        dzurb_roof(l,4) = 0.010_rk8
        dzurb_roof(l,5) = 0.030_rk8
        write(stdout,*)'Total thickness of roof: ',sum(dzurb_roof(l,:))
        write(stdout,*)'Roof layer thicknesses: ',dzurb_roof(l,:)

        ziurb_wall(l,0) = 0.
        ziurb_wall(l,1) = dzurb_wall(l,1)
        do j = 2 , nlevurb
          ziurb_wall(l,j) = sum(dzurb_wall(l,1:j))
        end do
        write(stdout,*)'Wall layer interface depths: ',ziurb_wall(l,:)

        ziurb_roof(l,0) = 0.
        ziurb_roof(l,1) = dzurb_roof(l,1)
        do j = 2 , nlevurb
          ziurb_roof(l,j) = sum(dzurb_roof(l,1:j))
        end do
        write(stdout,*)'Roof layer interface depths: ',ziurb_roof(l,:)
#elif (defined MEXICOCITY)
        zurb_wall(l,1) = 0.015_rk8/2._rk8
        zurb_wall(l,2) = zurb_wall(l,1) + 0.015_rk8/2._rk8 + 0.120_rk8/2._rk8
        zurb_wall(l,3) = zurb_wall(l,2) + 0.120_rk8/2._rk8 + 0.150_rk8/2._rk8
        zurb_wall(l,4) = zurb_wall(l,3) + 0.150_rk8/2._rk8 + 0.150_rk8/2._rk8
        zurb_wall(l,5) = zurb_wall(l,4) + 0.150_rk8/2._rk8 + 0.015_rk8/2._rk8

        zurb_roof(l,1) = 0.010_rk8/2._rk8
        zurb_roof(l,2) = zurb_roof(l,1) + 0.010_rk8/2._rk8 + 0.050_rk8/2._rk8
        zurb_roof(l,3) = zurb_roof(l,2) + 0.050_rk8/2._rk8 + 0.050_rk8/2._rk8
        zurb_roof(l,4) = zurb_roof(l,3) + 0.050_rk8/2._rk8 + 0.050_rk8/2._rk8
        zurb_roof(l,5) = zurb_roof(l,4) + 0.050_rk8/2._rk8 + 0.025_rk8/2._rk8

        dzurb_wall(l,1) = 0.015_rk8
        dzurb_wall(l,2) = 0.120_rk8
        dzurb_wall(l,3) = 0.150_rk8
        dzurb_wall(l,4) = 0.150_rk8
        dzurb_wall(l,5) = 0.015_rk8
        write(stdout,*)'Total thickness of wall: ',sum(dzurb_wall(l,:))
        write(stdout,*)'Wall layer thicknesses: ',dzurb_wall(l,:)

        dzurb_roof(l,1) = 0.010_rk8
        dzurb_roof(l,2) = 0.050_rk8
        dzurb_roof(l,3) = 0.050_rk8
        dzurb_roof(l,4) = 0.050_rk8
        dzurb_roof(l,5) = 0.025_rk8
        write(stdout,*)'Total thickness of roof: ',sum(dzurb_roof(l,:))
        write(stdout,*)'Roof layer thicknesses: ',dzurb_roof(l,:)

        ziurb_wall(l,0) = 0.
        ziurb_wall(l,1) = dzurb_wall(l,1)
        do j = 2 , nlevurb
          ziurb_wall(l,j) = sum(dzurb_wall(l,1:j))
        end do
        write(stdout,*)'Wall layer interface depths: ',ziurb_wall(l,:)

        ziurb_roof(l,0) = 0.
        ziurb_roof(l,1) = dzurb_roof(l,1)
        do j = 2 , nlevurb
          ziurb_roof(l,j) = sum(dzurb_roof(l,1:j))
        end do
        write(stdout,*)'Roof layer interface depths: ',ziurb_roof(l,:)
#else
        do j = 1 , nlevurb
          zurb_wall(l,j) = (j-0.5)*(thick_wall(l)/float(nlevurb))  !node depths
        end do
        do j = 1 , nlevurb
          zurb_roof(l,j) = (j-0.5)*(thick_roof(l)/float(nlevurb))  !node depths
        end do

        !thickness b/n two interfaces
        dzurb_roof(l,1) = 0.5*(zurb_roof(l,1)+zurb_roof(l,2))
        do j = 2 , nlevurb-1
          dzurb_roof(l,j)= 0.5*(zurb_roof(l,j+1)-zurb_roof(l,j-1))
        end do
        dzurb_roof(l,nlevurb) = zurb_roof(l,nlevurb)-zurb_roof(l,nlevurb-1)

        !thickness b/n two interfaces
        dzurb_wall(l,1) = 0.5*(zurb_wall(l,1)+zurb_wall(l,2))
        do j = 2 , nlevurb-1
          dzurb_wall(l,j)= 0.5*(zurb_wall(l,j+1)-zurb_wall(l,j-1))
        end do
        dzurb_wall(l,nlevurb) = zurb_wall(l,nlevurb)-zurb_wall(l,nlevurb-1)

        ziurb_wall(l,0) = 0.
        do j = 1 , nlevurb-1
          !interface depths
          ziurb_wall(l,j) = 0.5*(zurb_wall(l,j)+zurb_wall(l,j+1))
        end do
        ziurb_wall(l,nlevurb) = zurb_wall(l,nlevurb) + &
                0.5*dzurb_wall(l,nlevurb)

        ziurb_roof(l,0) = 0.
        do j = 1 , nlevurb-1
          !interface depths
          ziurb_roof(l,j) = 0.5*(zurb_roof(l,j)+zurb_roof(l,j+1))
        end do
        ziurb_roof(l,nlevurb) = zurb_roof(l,nlevurb) + &
                0.5*dzurb_roof(l,nlevurb)
#endif
      end if
    end do

    ! Grid level initialization
    do g = begg , endg
      ! VOC emission factors
      ! Set gridcell and landunit indices
      efisop(:,g)=efisop2d(:,g)
    end do

    ! --------------------------------------------------------------------
    ! Initialize soil and lake levels
    ! Initialize soil color, thermal and hydraulic properties
    ! --------------------------------------------------------------------

    nzero_slope = 0
    ! Column level initialization
    do c = begc , endc

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)

      ! initialize maximum daylength, based on latitude and maximum declination
      ! maximum declination hardwired for present-day orbital parameters,
      ! +/- 23.4667 degrees = +/- 0.409571 radians, use negative value
      ! for S. Hem
      max_decl = 0.409571
      if (lat(g) < 0._rk8) max_decl = -max_decl
      temp = -(sin(lat(g))*sin(max_decl))/(cos(lat(g)) * cos(max_decl))
      temp = min(1._rk8,max(-1._rk8,temp))
      max_dayl(c) = 2.0_rk8 * 13750.9871_rk8 * acos(temp)

      ! Initialize restriction for min of soil potential (mm)
      smpmin(c) = -1.e8_rk8

      ! Decay factor (m)
      hkdepth(c) = 1._rk8/2.5_rk8

      ! Maximum saturated fraction
      wtfact(c) = gti(g)
#if (defined VICHYDRO)
      b_infil(c) = b2d(g)
      ds(c)      = ds2d(g)
      dsmax(c)   = dsmax2d(g)
      Wsvic(c)   = ws2d(g)
#endif

      ! GDP data added by F. Li and S. Levis
      gdp_lf(c) = gdp(g)

      ! peatf data added by F. Li and S. Levis
      peatf_lf(c) = peatf(g)

      ! abm data added by F. Li and S. Levis
      abm_lf(c) = abm(g)

      ! Parameters for calculation of finundated
#ifdef LCH4
      if ( .not. fin_use_fsat ) then
        !zwt0(c) = zwt0_in(g)
        !f0(c)   = f0_in(g)
        !p3(c)   = p3_in(g)
        k(c) = k_in(g)
        q(c) = q_in(g)
        v(c) = v_in(g)
        maxf(c) = maxf_in(g)
      end if
      ! Methane pH factor
      if ( usephfact ) pH(c) = pH_in(g)
#endif

#ifdef CN
       q10(c) = min(max(q10_in(g),1.0_rk8),3.0_rk8)
       ndep(c) = ndep_in(g)
#endif

      ! Lake data
      lakedepth(c) = lakedepth_in(g)
      etal(c) = etal_in(g)
      lakefetch(c) = lakefetch_in(g)

      ! Topographic variables
      topo_std(c) = std(g)
      if ( pctspec(g) >= 100.0_rk8-mach_eps )then
        ! Zero out slope over ALL special land-units
        topo_slope(c) = 0.0_rk8
        nzero_slope   = nzero_slope + 1
      else
        ! check for near zero slopes, set minimum value
        topo_slope(c) = max(tslope(g),0.2_rk8)
      end if

      n_melt(c) = 200.0/max(10.0_rk8,topo_std(c))

      ! microtopographic parameter, units are meters
      minslope=0.05
      slopemax=0.4_rk8
      maxslope=(slopemax - minslope)/(slopemax)

      ! try smooth function of slope
      slopebeta=3._rk8
      slopemax=0.4_rk8
      slope0=slopemax**(-1._rk8/slopebeta)
      micro_sigma(c) = (topo_slope(c) + slope0)**(-slopebeta)

      ! determine h2osfc threshold ("fill & spill" concept)
      if ( micro_sigma(c) > 1.e-6_rk8 ) then
        d = 0.0_rk8
        do p=1,4
          fd = 0.5_rk8*(1.0_rk8+erf(d/(micro_sigma(c)*sqrt(2.0_rk8)))) - pc
          dfdd = exp(-d**2/(2.0_rk8*micro_sigma(c)**2)) / &
                  (micro_sigma(c)*sqrt(2.0_rk8*rpi))
          d = d - fd/dfdd
        end do
        h2osfc_thresh(c) = 0.5_rk8*d*(1.0_rk8 + &
                erf(d/(micro_sigma(c)*sqrt(2.0_rk8)))) + &
                micro_sigma(c)/sqrt(2.0_rk8*rpi) * &
                exp(-d**2/(2.0_rk8*micro_sigma(c)**2))
        h2osfc_thresh(c) = 1.e3_rk8 * h2osfc_thresh(c) !convert to mm from meters
      else
        micro_sigma(c) = 9.99e-7_rk8
        h2osfc_thresh(c) = 0._rk8
      end if

      if ( h2osfcflag == 0 ) then
        slopemax = 0.05_rk8
        micro_sigma(c) = slopemax
        ! set to zero for no h2osfc (w/frac_infclust =large)
        h2osfc_thresh(c) = 0._rk8
      end if

      ! Soil color
      isoicol(c) = soic2d(g)

      ! Soil hydraulic and thermal properties
      ! Note that urban roof, sunwall and shadewall thermal properties used to
      ! derive thermal conductivity and heat capacity are set to special
      ! value because thermal conductivity and heat capacity for urban
      ! roof, sunwall and shadewall are prescribed in SoilThermProp.F90 in
      ! SoilTemperatureMod.F90
      ! Lakes will be set in initSLake. This could also be done here to
      ! facilitate changing soil properties everywhere, but there may
      ! also be reasons to keep lakes separate (e.g. if lake-specific
      ! soil data became available, if thermokarst lakes were treated, etc.)
      if ( ltype(l)==istwet .or. ltype(l)==istice ) then
        do lev = 1 , nlevgrnd
          bsw(c,lev)    = spval
          watsat(c,lev) = spval
          watfc(c,lev)  = spval
          hksat(c,lev)  = spval
          sucsat(c,lev) = spval
          tkmg(c,lev)   = spval
          tksatu(c,lev) = spval
          tkdry(c,lev)  = spval
          if ( ltype(l) == istwet .and. lev > nlevsoi ) then
            csol(c,lev) = csol_bedrock
          else
            csol(c,lev)= spval
          end if
          watdry(c,lev) = spval
          watopt(c,lev) = spval
          bd(c,lev) = spval
          if ( lev <= nlevsoi ) then
            cellsand(c,lev) = spval
            cellclay(c,lev) = spval
            cellorg(c,lev)  = spval
          end if
        end do
#if (defined VICHYDRO)
        do lev = 1 , nlayer
          sandcol(c,lev)   = spval
          claycol(c,lev)   = spval
          om_fraccol(c,lev) = spval
          porosity(c,lev)  = spval
          max_moist(c,lev) = spval
          expt(c,lev)      = spval
          ksat(c,lev)      = spval
          phi_s(c,lev)     = spval
          depth(c,lev)     = spval
        end do
#endif
      else if (ltype(l)==isturb .and. &
              (ctype(c) /= icol_road_perv) .and. &
              (ctype(c) /= icol_road_imperv) ) then
        ! Urban Roof, sunwall, shadewall properties set to special value
        do lev = 1 , nlevgrnd
          watsat(c,lev) = spval
          watfc(c,lev)  = spval
          bsw(c,lev)    = spval
          hksat(c,lev)  = spval
          sucsat(c,lev) = spval
          tkmg(c,lev)   = spval
          tksatu(c,lev) = spval
          tkdry(c,lev)  = spval
          csol(c,lev)   = spval
          watdry(c,lev) = spval
          watopt(c,lev) = spval
          bd(c,lev) = spval
          if ( lev <= nlevsoi ) then
            cellsand(c,lev) = spval
            cellclay(c,lev) = spval
            cellorg(c,lev)  = spval
          end if
        end do
#if (defined VICHYDRO)
        do lev = 1 , nlayer
          sandcol(c,lev)   = spval
          claycol(c,lev)   = spval
          om_fraccol(c,lev) = spval
          porosity(c,lev)  = spval
          max_moist(c,lev) = spval
          expt(c,lev)      = spval
          ksat(c,lev)      = spval
          phi_s(c,lev)     = spval
          depth(c,lev)     = spval
        end do
#endif
      !else if (ltype(l) /= istdlak) then
      ! soil columns of both urban and non-urban types
      else
        do lev = 1 , nlevgrnd
          ! duplicate clay and sand values from last soil layer
          if ( more_vertlayers ) then
            if ( lev == 1 ) then
              clay    = clay3d(g,1)
              sand    = sand3d(g,1)
              om_frac = organic3d(g,1)/organic_max
            else if ( lev <= nlevsoi ) then
              do j = 1 , nlevsoifl-1
                if (zisoi(lev) >= zisoifl(j) .AND. &
                    zisoi(lev) < zisoifl(j+1)) then
                  clay    = clay3d(g,j+1)
                  sand    = sand3d(g,j+1)
                  om_frac = organic3d(g,j+1)/organic_max
                end if
              end do
            else
              clay    = clay3d(g,nlevsoifl)
              sand    = sand3d(g,nlevsoifl)
              om_frac = 0._rk8
            end if
          else
            ! duplicate clay and sand values from 10th soil layer
            if (lev <= nlevsoi) then
              clay    = clay3d(g,lev)
              sand    = sand3d(g,lev)
              om_frac = (organic3d(g,lev)/organic_max)**2._rk8
            else
              clay    = clay3d(g,nlevsoi)
              sand    = sand3d(g,nlevsoi)
              om_frac = 0._rk8
            end if
          end if
          if (ltype(l) == istdlak) then
            if (lev <= nlevsoi) then
              cellsand(c,lev) = sand
              cellclay(c,lev) = clay
              cellorg(c,lev)  = om_frac*organic_max
              watsat(c,lev) = 1.0_rk8
            end if
          else if (ltype(l) /= istdlak) then
            ! soil columns of both urban and non-urban types
            ! No organic matter for urban
            if (ltype(l)==isturb) then
              om_frac = 0._rk8
            end if
            if (lev <= nlevsoi) then
              cellsand(c,lev) = sand
              cellclay(c,lev) = clay
              cellorg(c,lev)  = om_frac*organic_max
            end if
#if (defined VICHYDRO)
            claycol(c,lev)    = clay
            sandcol(c,lev)    = sand
            om_fraccol(c,lev) = om_frac
#endif
            ! Note that the following properties are overwritten for urban
            ! impervious road layers that are not soil in SoilThermProp.F90
            ! within SoilTemperatureMod.F90
            watsat(c,lev) = 0.489_rk8 - 0.00126_rk8*sand
            bsw(c,lev) = 2.91 + 0.159*clay
            sucsat(c,lev) = 10._rk8 * ( 10._rk8**(1.88_rk8-0.0131_rk8*sand) )
            om_watsat = max(0.93_rk8 - 0.1_rk8*(zsoi(lev)/zsapric), 0.83_rk8)
            om_b = min(2.7_rk8 + 9.3_rk8*(zsoi(lev)/zsapric), 12.0_rk8)
            om_sucsat = min(10.3_rk8 - 0.2_rk8*(zsoi(lev)/zsapric), 10.1_rk8)
            om_hksat = max(0.28_rk8 - 0.2799_rk8*(zsoi(lev)/zsapric), 0.0001_rk8)

            bd = (1._rk8-watsat(c,lev))*2.7e3_rk8
            watsat(c,lev) = (1._rk8 - om_frac)*watsat(c,lev) + om_watsat*om_frac
            tkm = (1._rk8-om_frac)*(8.80_rk8*sand+2.92_rk8*clay) / &
                    (sand+clay)+om_tkm*om_frac ! W/(m K)
            bsw(c,lev) = (1._rk8-om_frac)*(2.91_rk8 + 0.159_rk8*clay) + om_frac*om_b
            sucsat(c,lev) = (1._rk8-om_frac)*sucsat(c,lev) + om_sucsat*om_frac
            xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s
            hksat_min(c,lev) = xksat

            ! perc_frac is zero unless perf_frac greater than percolation
            ! threshold
            if (om_frac > pcalpha) then
              perc_norm = (1._rk8 - pcalpha)**(-pcbeta)
              perc_frac = perc_norm*(om_frac - pcalpha)**pcbeta
            else
              perc_frac = 0._rk8
            end if
            ! uncon_frac is fraction of mineral soil plus fraction of
            ! "nonpercolating" organic soil
            uncon_frac=(1._rk8-om_frac)+(1._rk8-perc_frac)*om_frac
            ! uncon_hksat is series addition of mineral/organic conductivites
            if (om_frac < 1._rk8) then
              uncon_hksat=uncon_frac/((1._rk8-om_frac)/xksat &
                   +((1._rk8-perc_frac)*om_frac)/om_hksat)
            else
              uncon_hksat = 0._rk8
            end if
            hksat(c,lev) = uncon_frac*uncon_hksat + &
                    (perc_frac*om_frac)*om_hksat
            tkmg(c,lev) = tkm ** (1._rk8- watsat(c,lev))
            tksatu(c,lev) = tkmg(c,lev)*0.57_rk8**watsat(c,lev)
            tkdry(c,lev) = ((0.135_rk8*bd(c,lev) + 64.7_rk8) / &
                     (2.7e3_rk8 - 0.947_rk8*bd(c,lev)))*(1._rk8-om_frac) + &
                      om_tkd*om_frac
            csol(c,lev) = ((1._rk8-om_frac)*(2.128_rk8*sand+2.385_rk8*clay) / &
                    (sand+clay) + om_csol*om_frac)*1.e6_rk8  ! J/(m3 K)
            if (lev > nlevsoi) then
              csol(c,lev) = csol_bedrock
            end if
            watdry(c,lev) = watsat(c,lev) * &
                    (316230._rk8/sucsat(c,lev)) ** (-1._rk8/bsw(c,lev))
            watopt(c,lev) = watsat(c,lev) * &
                    (158490._rk8/sucsat(c,lev)) ** (-1._rk8/bsw(c,lev))
            !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
            ! water content at field capacity, defined as hk = 0.1 mm/day
            ! used eqn (7.70) in CLM3 technote with
            ! k = 0.1 (mm/day) / secspday (day/sec)
            watfc(c,lev) = watsat(c,lev) * (0.1_rk8 / &
                    (hksat(c,lev)*secspday))**(1._rk8/(2._rk8*bsw(c,lev)+3._rk8))
          end if
        end do
        !
        ! Urban pervious and impervious road
        !
        ! Impervious road layers -- same as above except set watdry
        ! and watopt as missing
        if (ctype(c) == icol_road_imperv) then
          do lev = 1,nlevgrnd
            watdry(c,lev) = spval
            watopt(c,lev) = spval
          end do
          ! pervious road layers -- same as above except also set
          ! rootfr_road_perv
          ! Currently, pervious road has same properties as soil
        else if (ctype(c) == icol_road_perv) then
          do lev = 1, nlevgrnd
            rootfr_road_perv(c,lev) = 0._rk8
          end do
          do lev = 1,nlevsoi
            rootfr_road_perv(c,lev) = 0.1_rk8  ! uniform profile
          end do
        end if
      end if

      ! Define non-lake levels, layers and interfaces
      ! Lakes will be set in initSLake
      if (ltype(l) == isturb) then
        if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall) then
          z(c,1:nlevurb)  = zurb_wall(l,1:nlevurb)
          zi(c,0:nlevurb) = ziurb_wall(l,0:nlevurb)
          dz(c,1:nlevurb) = dzurb_wall(l,1:nlevurb)
          if (nlevurb < nlevgrnd) then
            z(c,nlevurb+1:nlevgrnd)  = spval
            zi(c,nlevurb+1:nlevgrnd) = spval
            dz(c,nlevurb+1:nlevgrnd) = spval
          end if
        else if (ctype(c)==icol_roof) then
          z(c,1:nlevurb)  = zurb_roof(l,1:nlevurb)
          zi(c,0:nlevurb) = ziurb_roof(l,0:nlevurb)
          dz(c,1:nlevurb) = dzurb_roof(l,1:nlevurb)
          if (nlevurb < nlevgrnd) then
            z(c,nlevurb+1:nlevgrnd)  = spval
            zi(c,nlevurb+1:nlevgrnd) = spval
            dz(c,nlevurb+1:nlevgrnd) = spval
          end if
        else
          z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
#if (defined VICHYDRO)
          depth(c,:) = 0._rk8
          ivicstrt = 1
          do ivic = 1,nlayer
            ivicend = ivicstrt+nlvic(ivic)-1
            do j = ivicstrt,ivicend
              depth(c,ivic) = depth(c,ivic)+dz(c,j)
            end do
            ivicstrt = ivicend+1
          end do
          depth(c, nlayer+1:nlayert) = dz(c, nlevsoi+1:nlevgrnd)
          ! Column level initialization
          ! create weights to map soil moisture profiles (10 layer) to
          ! 3 layers for VIC hydrology, M.Huang
          call initCLMVICMap(c)
          call initSoilParVIC(c, claycol, sandcol, om_fraccol)
#endif
        end if
      else if (ltype(l) /= istdlak) then
        z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
        zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
        dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
#if (defined VICHYDRO)
        depth(c,:) = 0._rk8
        ivicstrt = 1
        do ivic = 1,nlayer
          ivicend = ivicstrt+nlvic(ivic)-1
          do j = ivicstrt,ivicend
            depth(c,ivic) = depth(c,ivic)+dz(c,j)
          end do
          ivicstrt = ivicend+1
        end do
        depth(c, nlayer+1:nlayert) = dz(c, nlevsoi+1:nlevgrnd)
        ! Column level initialization
        ! create weights to map soil moisture profiles (10 layer) to
        ! 3 layers for VIC hydrology, M.Huang
        call initCLMVICMap(c)
        call initSoilParVIC(c, claycol, sandcol, om_fraccol)
#endif
      end if

      ! Initialize terms needed for dust model
      clay = clay3d(g,1)
      gwc_thr(c) = 0.17_rk8 + 0.14_rk8*clay*0.01_rk8
      mss_frc_cly_vld(c) = min(clay*0.01_rk8, 0.20_rk8)
    end do

    if ( nzero_slope > 0 )then
      write(stdout,'(A,I6,A)') " Set", nzero_slope, &
                             " 100% special land-units points to zero slope"
    end if

    ! pft level initialization
    do p = begp, endp

      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).

      c = pcolumn(p)
      if (ivt(p) /= noveg) then
        do lev = 1, nlevgrnd
          rootfr(p,lev) = 0._rk8
        end do
        do lev = 1, nlevsoi-1
          rootfr(p,lev) = .5_rk8*( exp(-roota_par(ivt(p)) * zi(c,lev-1))  &
                               + exp(-rootb_par(ivt(p)) * zi(c,lev-1))  &
                               - exp(-roota_par(ivt(p)) * zi(c,lev  ))  &
                               - exp(-rootb_par(ivt(p)) * zi(c,lev  )) )
        end do
        rootfr(p,nlevsoi) = .5_rk8*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
                                + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )
        rootfr(p,nlevsoi+1:nlevgrnd) =  0.0_rk8

!#if (defined CN)
!        ! replacing the exponential rooting distribution
!        ! with a linear decrease, going to zero at the bottom of the lowest
!        ! soil layer for woody pfts, but going to zero at the bottom of
!        ! layer 8 for non-woody pfts.  This corresponds to 3.43 m for woody
!        ! bottom, vs 1.38 m for non-woody bottom.
!        if ( abs(woody(ivt(p))-1._rk8) < epsilon(1.0) ) then
!           bottom = nlevsoi
!           slope = -2._rk8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._rk8/zi(c,bottom)
!           do lev = 1, bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_rk8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._rk8
!              end do
!           end if
!        else
!           bottom = 8
!           slope = -2._rk8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._rk8/zi(c,bottom)
!           do lev=1,bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_rk8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._rk8
!              end do
!           end if
!        end if
!#endif
      else
        rootfr(p,1:nlevsoi) = 0._rk8
      end if

      ! Initialize maximum allowed dew

      dewmx(p) = 0.1_rk8

      ! initialize rresis, for use in ecosystemdyn
      do lev = 1,nlevgrnd
        rresis(p,lev) = 0._rk8
      end do
    end do ! end pft level initialization

#ifdef CN
    ! ----------------------------
    ! Initialize time-constant arrays of decomposition constants
    ! -----------------------------
    if (myid == italk) then
      write(stdout,*) 'Initializing decomposition pools and transitions ...'
    end if
    call init_decompcascade(begc, endc)

    ! initialize the CN variables for special landunits, including lake points
    call CNiniSpecial()

#endif
    deallocate(gdp,peatf,abm) ! F. Li and S. Levis
    deallocate(soic2d,sand3d,clay3d,gti,organic3d)
    deallocate(zisoifl,zsoifl,dzsoifl)
    deallocate(temp_ef,efisop2d)
#ifdef LCH4
    !deallocate(zwt0_in, f0_in, p3_in)
     deallocate(k_in, q_in, v_in, maxf_in)
    if (usephfact) deallocate(pH_in)
#endif

#ifdef CN
    deallocate(q10_in)
    deallocate(ndep_in)
#endif

    deallocate(lakedepth_in)
    deallocate(etal_in)
    deallocate(lakefetch_in)
    if (nlevurb > 0) then
      deallocate(zurb_wall, zurb_roof, dzurb_wall, &
                 dzurb_roof, ziurb_wall, ziurb_roof)
    end if

    deallocate(std)
    deallocate(tslope)
#if (defined VICHYDRO)
    deallocate(b2d, ds2d, dsmax2d,ws2d)
    deallocate(sandcol, claycol, om_fraccol)
#endif

    ! Initialize SNICAR optical and aging parameters:
    call SnowOptics_init( )

    call SnowAge_init( )

    if (myid == italk) then
      write(stdout,*) 'Successfully initialized time invariant variables'
    end if
  end subroutine iniTimeConst

end module mod_clm_initimeconst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
