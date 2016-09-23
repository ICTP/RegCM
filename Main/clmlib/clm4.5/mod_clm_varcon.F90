module mod_clm_varcon
  !
  ! Module containing various model constants
  !
  use mod_realkinds
  use mod_intkinds
  use mod_constants , only : cpd , pdbratio , secpd , earthrad , tzero , &
    rgasmol , rhoice , rhoh2o , wlhf , wlhs , wlhv , spcpice , spcpfw ,  &
    egrav , rwat , vonkar , sigm , mathpi , regrav , rdry ,              &
    lnd_sfcemiss , ocn_sfcemiss
  use mod_clm_varpar , only : numrad , nlevgrnd , nlevlak , &
     nlevdecomp_full , ngases
#if (defined VICHYDRO)
  use mod_clm_varpar   , only: nlayer
#endif

  implicit none

  private

  save

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(rk8) , public , parameter :: rpi    = mathpi

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  ! Fsca shape parameter
  real(rk8) , public , parameter :: n_melt = 0.7_rk8
  ! Soil ice impedance factor
  real(rk8) , public , parameter :: e_ice = 6.0_rk8
  ! Threshold probability
  real(rk8) , public , parameter :: pc = 0.4_rk8
  ! Connectivity exponent
  real(rk8) , public , parameter :: mu = 0.13889_rk8
  ! Gravity constant [m/s2]
  real(rk8) , public , parameter :: grav   = egrav
  ! Stefan-Boltzmann constant  [W/m2/K4]
  real(rk8) , public , parameter :: sb     = sigm
  ! Von Karman constant [-]
  real(rk8) , public , parameter :: vkc    = vonkar
  ! Gas constant for dry air [J/kg/K]
  real(rk8) , public , parameter :: rair   = rdry
  ! Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(rk8) , public , parameter :: roverg = rwat*regrav*1000._rk8
  ! Specific heat of water [J/kg-K]
  real(rk8) , public , parameter :: cpliq  = spcpfw
  ! Specific heat of ice [J/kg-K]
  real(rk8) , public , parameter :: cpice  = spcpice
  ! Specific heat of dry air [J/kg/K]
  real(rk8) , public , parameter :: cpair  = cpd
  ! Latent heat of evap for water [J/kg]
  real(rk8) , public , parameter :: hvap   = wlhv
  ! Latent heat of sublimation    [J/kg]
  real(rk8) , public , parameter :: hsub   = wlhs
  ! Latent heat of fusion for ice [J/kg]
  real(rk8) , public , parameter :: hfus   = wlhf
  ! Density of liquid water [kg/m3]
  real(rk8) , public , parameter :: denh2o = rhoh2o
  ! Density of ice [kg/m3]
  real(rk8) , public , parameter :: denice = rhoice
  ! Universal gas constant [J/K/kmole]
  real(rk8) , public , parameter :: rgas   = rgasmol*1000.0_rk8
  ! Thermal conductivity of air   [W/m/K]
  real(rk8) , public , parameter :: tkair  = 0.023_rk8
  ! Thermal conductivity of ice   [W/m/K]
  real(rk8) , public , parameter :: tkice  = 2.290_rk8
  ! Thermal conductivity of water [W/m/K]
  real(rk8) , public , parameter :: tkwat  = 0.57_rk8
  ! Freezing temperature [K]
  real(rk8) , public , parameter :: tfrz   = tzero
  ! Constant atmospheric O2 molar ratio (mol/mol)
  real(rk8) , public , parameter :: o2_molar_const = 0.209_rk8

  ! Bulk density snow (kg/m**3)
  real(rk8) , public , parameter :: bdsno = 250._rk8
  ! Constant for aerodynamic parameter weighting
  real(rk8) , public , parameter :: alpha_aero = 1.0_rk8
  ! Critical value of elai+esai for which aerodynamic parameters are maximum
  real(rk8) , public , parameter :: tlsai_crit = 2.0_rk8
  ! Minimum soil moisture (mm)
  real(rk8) , public , parameter :: watmin = 0.1_rk8

  ! Radius of earth (km)
  real(rk8) , public , parameter :: re = earthrad*0.001_rk8

  ! Degree's earth rotates per second
  real(rk8) , public , parameter :: degpsec = 15._rk8/3600.0_rk8

  ! Seconds per day
  real(rk8) , public , parameter ::  secspday = secpd
  ! Integer seconds per day
  integer(ik4) , public , parameter :: isecspday = 86400
  ! Special value for real data
  real(rk8) , public , parameter ::  spval = 1.e20_rk8
  real(rk4) , public , parameter ::  rspval = 1.e20_rk4
  ! Special value for int data
  integer(ik4) , public , parameter :: ispval = -9999

  ! These are tunable constants from clm2_3

  ! Roughness length for soil [m]
  real(rk8) , public , parameter :: zlnd = 0.01_rk8
  ! Roughness length for snow [m]
  real(rk8) , public , parameter :: zsno = 0.0024_rk8
  ! Drag coefficient for soil under canopy [-]
  real(rk8) , public , parameter :: csoilc = 0.004_rk8
  ! Tuning factor to turn first layer T into surface T
  real(rk8) , public , parameter :: capr   = 0.34_rk8
  ! Crank Nicholson factor between 0 and 1
  real(rk8) , public , parameter :: cnfac  = 0.5_rk8
  ! Irreducible water saturation of snow
  real(rk8) , public , parameter :: ssi    = 0.033_rk8
  ! Water impremeable if porosity less than wimp
  real(rk8) , public , parameter :: wimp   = 0.05_rk8
  ! Ponding depth (mm)
  real(rk8) , public , parameter :: pondmx = 0.0_rk8
  ! Ponding depth for urban roof and impervious road (mm)
  real(rk8) , public , parameter :: pondmx_urban = 1.0_rk8

  ! Thermal conductivity of 'typical' saturated granitic rock
  ! (Clauser and Huenges, 1995)(W/m/K)
  real(rk8) , public , parameter :: thk_bedrock = 3.0_rk8

  !!! C13
  ! preindustrial value for atmospheric del13C
  real(rk8) , public , parameter :: preind_atm_del13c = -6.0_rk8
  real(rk8) , public , parameter :: preind_atm_ratio = pdbratio + &
          (preind_atm_del13c * pdbratio)/1000.0_rk8  ! 13C/12C
  real(rk8) , public , parameter :: c13ratio = preind_atm_ratio / &
          (1.0_rk8+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  !!! C14
  real(rk8) , public , parameter :: c14ratio = 1.e-12_rk8 ! 1._rk8

  ! Note that the wasteheat factors are currently set to zero until a better
  ! parameterization can be developed
  ! The prior parameterization appeared to be significantly overestimating
  ! wasteheat
  ! Wasteheat factor for urban heating (-)
  real(rk8) , public , parameter :: ht_wasteheat_factor = 0.0_rk8
  ! Wasteheat factor for urban air conditioning (-)
  real(rk8) , public , parameter :: ac_wasteheat_factor = 0.0_rk8
  ! Limit on wasteheat (W/m2)
  real(rk8) , public , parameter :: wasteheat_limit = 100._rk8

  ! Max allowed snow thickness (mm H2O)
  real(rk8) , public , parameter :: h2osno_max = 1000._rk8
  integer(ik4) :: i  ! loop index


#ifdef NITRIF_DENITRIF
  ! Fraction of N lost as N2O in nitrification (Parton et al., 2001)
  ! real(rk8) , public , parameter :: nitrif_n2o_loss_frac = 0.02_rk8
  ! Fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(rk8) , public , parameter :: nitrif_n2o_loss_frac = 6.e-4_rk8
  ! Fraction of N mineralized that is dieverted to the nitrification stream
  ! (Parton et al., 2001)
  real(rk8) , public , parameter :: frac_minrlztn_to_no3 = 0.2_rk8
#endif

  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "land unit " types
  !   1     soil (includes vegetated landunits)
  !   2     land ice (glacier)
  !   3     deep lake
  !  (DEPRECATED: New lake model has variable depth) 4     shallow lake
  !   5     wetland (swamp, marsh, etc.)
  !   6     urban
  !   7     land ice (glacier) with multiple elevation classes
  !   8     crop

  ! Soil landunit type
  integer(ik4) , public , parameter :: istsoil    = 1
  ! Land ice landunit type
  integer(ik4) , public , parameter :: istice     = 2
  ! Deep lake landunit type
  ! Not used; now 3 is used for all lakes, which have variable depth.
  integer(ik4) , public , parameter :: istdlak    = 3
  ! Shallow lake landunit type
  integer(ik4) , public , parameter :: istslak    = 4
  ! Wetland landunit type
  integer(ik4) , public , parameter :: istwet     = 5
  ! Urban landunit type
  integer(ik4) , public , parameter :: isturb     = 6
  ! Land ice (multiple elevation classes) landunit type
  integer(ik4) , public , parameter :: istice_mec = 7
  ! Crop landunit type
  integer(ik4) , public , parameter :: istcrop    = 8
  ! Maximum value that clm3%g%l%itype can have
  ! (i.e., largest value in the above list)
  integer(ik4) , public , parameter :: max_lunit  = 8

  ! urban column types

  integer(ik4) , public , parameter :: icol_roof        = 61
  integer(ik4) , public , parameter :: icol_sunwall     = 62
  integer(ik4) , public , parameter :: icol_shadewall   = 63
  integer(ik4) , public , parameter :: icol_road_imperv = 64
  integer(ik4) , public , parameter :: icol_road_perv   = 65

  ! urban density types

  integer(ik4) , public , parameter :: udens_base     = 600
  integer(ik4) , public , parameter :: udens_tbd      = 601
  integer(ik4) , public , parameter :: udens_hd       = 602
  integer(ik4) , public , parameter :: udens_md       = 603

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------

  ! Wet soil albedo by color class and waveband (1=vis,2=nir)
  real(rk8) , public , allocatable :: albsat(:,:)
  ! Dry soil albedo by color class and waveband (1=vis,2=nir)
  real(rk8) , public , allocatable :: albdry(:,:)
  ! Two-stream parameter betad for snow
  real(rk8) , public , parameter :: betads  = 0.5_rk8
  ! Two-stream parameter betai for snow
  real(rk8) , public , parameter :: betais  = 0.5_rk8
  ! Two-stream parameter omega for snow by band
  real(rk8) , public :: omegas(numrad)
  data (omegas(i),i=1,numrad) /0.8_rk8, 0.4_rk8/

  ! Lake Model Constants will be defined in SLakeCon.

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  ! Lake z (layers)
  real(rk8) , public , pointer , dimension(:) :: zlak
  ! Lake dz (thickness)
  real(rk8) , public , pointer , dimension(:) :: dzlak
  ! Soil z (layers)
  real(rk8) , public , pointer , dimension(:) :: zsoi
  ! Soil dz (thickness)
  real(rk8) , public , pointer , dimension(:) :: dzsoi
  ! Soil zi (interfaces)
  real(rk8) , public , pointer , dimension(:) :: zisoi
  ! Soil dz (thickness)
  real(rk8) , public , pointer , dimension(:) :: dzsoi_decomp

#if (defined VICHYDRO)
  ! Number of CLM layers in each VIC layer (#)
  integer(ik4) , public , pointer  , dimension(:) :: nlvic
#endif

  !------------------------------------------------------------------
  ! (Non-tunable) Constants for the CH4 submodel
  !       (Tuneable constants in ch4varcon)
  !------------------------------------------------------------------
  ! Note some of these constants are also used in CNNitrifDenitrifMod

  ! Molar mass of C atoms (g/mol)
  real(rk8) , public , parameter :: catomw = 12.011_rk8

  ! Schmidt # calculation constants (spp, #)
  real(rk8) , public :: s_con(ngases,4)
  data (s_con(1,i),i=1,4) /1898_rk8, -110.1_rk8, 2.8340_rk8, -0.027910_rk8/ ! CH4
  data (s_con(2,i),i=1,4) /1801_rk8, -120.1_rk8, 3.7818_rk8, -0.047608_rk8/ ! O2
  data (s_con(3,i),i=1,4) /1911_rk8, -113.7_rk8, 2.9670_rk8, -0.029430_rk8/ ! CO2

  ! Water diffusivity constants (spp, #)  (mult. by 10^-4)
  real(rk8) , public :: d_con_w(ngases,3)
  data (d_con_w(1,i),i=1,3) /0.9798_rk8, 0.02986_rk8, 0.0004381_rk8/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.1720_rk8, 0.03443_rk8, 0.0005048_rk8/ ! O2
  data (d_con_w(3,i),i=1,3) /0.9390_rk8, 0.02671_rk8, 0.0004095_rk8/ ! CO2

  ! Gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  real(rk8) , public :: d_con_g(ngases,2)
  data (d_con_g(1,i),i=1,2) /0.1875_rk8, 0.00130_rk8/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_rk8, 0.00117_rk8/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_rk8, 0.00090_rk8/ ! CO2

  ! Vonstant (K) for Henry's law (4.12, Wania)
  real(rk8) , public :: c_h_inv(ngases)
  data c_h_inv(1:3) /1600._rk8, 1500._rk8, 2400._rk8/ ! CH4, O2, CO2
  ! Henry's constant (L.atm/mol) at standard temperature (298K)
  real(rk8) , public :: kh_theta(ngases)
  data kh_theta(1:3) /714.29_rk8, 769.23_rk8, 29.4_rk8/ ! CH4, O2, CO2
  ! Base temperature for calculation of Henry's constant (K)
  real(rk8) , public , parameter :: kh_tbase = 298._rk8

  ! Initialze constants that need to be initialized
  public :: clm_varcon_init

  contains
  !
  ! This subroutine initializes constants in clm_varcon. MUST be called
  ! after the clm_varpar_init.
  !
  subroutine clm_varcon_init()
    implicit none
    allocate( zlak(1:nlevlak) )
    allocate( dzlak(1:nlevlak) )
    allocate( zsoi(1:nlevgrnd) )
    allocate( dzsoi(1:nlevgrnd) )
    allocate( zisoi(0:nlevgrnd) )
    allocate( dzsoi_decomp(1:nlevdecomp_full) )
#if (defined VICHYDRO)
    allocate( nlvic(1:nlayer))
#endif
  end subroutine clm_varcon_init

end module mod_clm_varcon
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
