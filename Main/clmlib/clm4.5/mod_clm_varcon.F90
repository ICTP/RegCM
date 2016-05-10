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

  real(rkx) , public , parameter :: rpi    = mathpi

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  ! Fsca shape parameter
  real(rkx) , public , parameter :: n_melt = 0.7_rkx
  ! Soil ice impedance factor
  real(rkx) , public , parameter :: e_ice = 6.0_rkx
  ! Threshold probability
  real(rkx) , public , parameter :: pc = 0.4_rkx
  ! Connectivity exponent
  real(rkx) , public , parameter :: mu = 0.13889_rkx
  ! Gravity constant [m/s2]
  real(rkx) , public , parameter :: grav   = egrav
  ! Stefan-Boltzmann constant  [W/m2/K4]
  real(rkx) , public , parameter :: sb     = sigm
  ! Von Karman constant [-]
  real(rkx) , public , parameter :: vkc    = vonkar
  ! Gas constant for dry air [J/kg/K]
  real(rkx) , public , parameter :: rair   = rdry
  ! Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(rkx) , public , parameter :: roverg = rwat*regrav*1000._rkx
  ! Specific heat of water [J/kg-K]
  real(rkx) , public , parameter :: cpliq  = spcpfw
  ! Specific heat of ice [J/kg-K]
  real(rkx) , public , parameter :: cpice  = spcpice
  ! Specific heat of dry air [J/kg/K]
  real(rkx) , public , parameter :: cpair  = cpd
  ! Latent heat of evap for water [J/kg]
  real(rkx) , public , parameter :: hvap   = wlhv
  ! Latent heat of sublimation    [J/kg]
  real(rkx) , public , parameter :: hsub   = wlhs
  ! Latent heat of fusion for ice [J/kg]
  real(rkx) , public , parameter :: hfus   = wlhf
  ! Density of liquid water [kg/m3]
  real(rkx) , public , parameter :: denh2o = rhoh2o
  ! Density of ice [kg/m3]
  real(rkx) , public , parameter :: denice = rhoice
  ! Universal gas constant [J/K/kmole]
  real(rkx) , public , parameter :: rgas   = rgasmol*1000.0_rkx
  ! Thermal conductivity of air   [W/m/K]
  real(rkx) , public , parameter :: tkair  = 0.023_rkx
  ! Thermal conductivity of ice   [W/m/K]
  real(rkx) , public , parameter :: tkice  = 2.290_rkx
  ! Thermal conductivity of water [W/m/K]
  real(rkx) , public , parameter :: tkwat  = 0.57_rkx
  ! Freezing temperature [K]
  real(rkx) , public , parameter :: tfrz   = tzero
  ! Critical temperature to determine rain or snow
  real(rkx) , public , parameter :: tcrit  = 2.5_rkx
  ! Constant atmospheric O2 molar ratio (mol/mol)
  real(rkx) , public , parameter :: o2_molar_const = 0.209_rkx

  ! Bulk density snow (kg/m**3)
  real(rkx) , public , parameter :: bdsno = 250._rkx
  ! Constant for aerodynamic parameter weighting
  real(rkx) , public , parameter :: alpha_aero = 1.0_rkx
  ! Critical value of elai+esai for which aerodynamic parameters are maximum
  real(rkx) , public , parameter :: tlsai_crit = 2.0_rkx
  ! Minimum soil moisture (mm)
  real(rkx) , public , parameter :: watmin = 0.01_rkx

  ! Radius of earth (km)
  real(rkx) , public , parameter :: re = earthrad*0.001_rkx

  ! Degree's earth rotates per second
  real(rkx) , public , parameter :: degpsec = 15._rkx/3600.0_rkx

  ! Seconds per day
  real(rkx) , public , parameter ::  secspday = secpd
  ! Integer seconds per day
  integer(ik4) , public , parameter :: isecspday = 86400
  ! Special value for real data
  real(rkx) , public , parameter ::  spval = 1.e36_rkx
  real(rk4) , public , parameter ::  rspval = 1.e36_rkx
  ! Special value for int data
  integer(ik4) , public , parameter :: ispval = -9999

  ! These are tunable constants from clm2_3

  ! Roughness length for soil [m]
  real(rkx) , public , parameter :: zlnd = 0.01_rkx
  ! Roughness length for snow [m]
  real(rkx) , public , parameter :: zsno = 0.0024_rkx
  ! Drag coefficient for soil under canopy [-]
  real(rkx) , public , parameter :: csoilc = 0.004_rkx
  ! Tuning factor to turn first layer T into surface T
  real(rkx) , public , parameter :: capr   = 0.34_rkx
  ! Crank Nicholson factor between 0 and 1
  real(rkx) , public , parameter :: cnfac  = 0.5_rkx
  ! Irreducible water saturation of snow
  real(rkx) , public , parameter :: ssi    = 0.033_rkx
  ! Water impremeable if porosity less than wimp
  real(rkx) , public , parameter :: wimp   = 0.05_rkx
  ! Ponding depth (mm)
  real(rkx) , public , parameter :: pondmx = 0.0_rkx
  ! Ponding depth for urban roof and impervious road (mm)
  real(rkx) , public , parameter :: pondmx_urban = 1.0_rkx

  ! Thermal conductivity of 'typical' saturated granitic rock
  ! (Clauser and Huenges, 1995)(W/m/K)
  real(rkx) , public , parameter :: thk_bedrock = 3.0_rkx

  !!! C13
  ! preindustrial value for atmospheric del13C
  real(rkx) , public , parameter :: preind_atm_del13c = -6.0_rkx
  real(rkx) , public , parameter :: preind_atm_ratio = pdbratio + &
          (preind_atm_del13c * pdbratio)/1000.0_rkx  ! 13C/12C
  real(rkx) , public , parameter :: c13ratio = preind_atm_ratio / &
          (1.0_rkx+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  !!! C14
  real(rkx) , public , parameter :: c14ratio = 1.e-12_rkx ! 1._rkx

  ! Note that the wasteheat factors are currently set to zero until a better
  ! parameterization can be developed
  ! The prior parameterization appeared to be significantly overestimating
  ! wasteheat
  ! Wasteheat factor for urban heating (-)
  real(rkx) , public , parameter :: ht_wasteheat_factor = 0.0_rkx
  ! Wasteheat factor for urban air conditioning (-)
  real(rkx) , public , parameter :: ac_wasteheat_factor = 0.0_rkx
  ! Limit on wasteheat (W/m2)
  real(rkx) , public , parameter :: wasteheat_limit = 100._rkx

  ! Max allowed snow thickness (mm H2O)
  real(rkx) , public , parameter :: h2osno_max = 1000._rkx
  integer(ik4) :: i  ! loop index


#ifdef NITRIF_DENITRIF
  ! Fraction of N lost as N2O in nitrification (Parton et al., 2001)
  ! real(rkx) , public , parameter :: nitrif_n2o_loss_frac = 0.02_rkx
  ! Fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(rkx) , public , parameter :: nitrif_n2o_loss_frac = 6.e-4_rkx
  ! Fraction of N mineralized that is dieverted to the nitrification stream
  ! (Parton et al., 2001)
  real(rkx) , public , parameter :: frac_minrlztn_to_no3 = 0.2_rkx
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
  real(rkx) , public , allocatable :: albsat(:,:)
  ! Dry soil albedo by color class and waveband (1=vis,2=nir)
  real(rkx) , public , allocatable :: albdry(:,:)
  ! Two-stream parameter betad for snow
  real(rkx) , public , parameter :: betads  = 0.5_rkx
  ! Two-stream parameter betai for snow
  real(rkx) , public , parameter :: betais  = 0.5_rkx
  ! Two-stream parameter omega for snow by band
  real(rkx) , public :: omegas(numrad)
  data (omegas(i),i=1,numrad) /0.8_rkx, 0.4_rkx/

  ! Lake Model Constants will be defined in SLakeCon.

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  ! Lake z (layers)
  real(rkx) , public , pointer , dimension(:) :: zlak
  ! Lake dz (thickness)
  real(rkx) , public , pointer , dimension(:) :: dzlak
  ! Soil z (layers)
  real(rkx) , public , pointer , dimension(:) :: zsoi
  ! Soil dz (thickness)
  real(rkx) , public , pointer , dimension(:) :: dzsoi
  ! Soil zi (interfaces)
  real(rkx) , public , pointer , dimension(:) :: zisoi
  ! Soil dz (thickness)
  real(rkx) , public , pointer , dimension(:) :: dzsoi_decomp

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
  real(rkx) , public , parameter :: catomw = 12.011_rkx

  ! Schmidt # calculation constants (spp, #)
  real(rkx) , public :: s_con(ngases,4)
  data (s_con(1,i),i=1,4) /1898_rkx, -110.1_rkx, 2.8340_rkx, -0.027910_rkx/ ! CH4
  data (s_con(2,i),i=1,4) /1801_rkx, -120.1_rkx, 3.7818_rkx, -0.047608_rkx/ ! O2
  data (s_con(3,i),i=1,4) /1911_rkx, -113.7_rkx, 2.9670_rkx, -0.029430_rkx/ ! CO2

  ! Water diffusivity constants (spp, #)  (mult. by 10^-4)
  real(rkx) , public :: d_con_w(ngases,3)
  data (d_con_w(1,i),i=1,3) /0.9798_rkx, 0.02986_rkx, 0.0004381_rkx/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.1720_rkx, 0.03443_rkx, 0.0005048_rkx/ ! O2
  data (d_con_w(3,i),i=1,3) /0.9390_rkx, 0.02671_rkx, 0.0004095_rkx/ ! CO2

  ! Gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  real(rkx) , public :: d_con_g(ngases,2)
  data (d_con_g(1,i),i=1,2) /0.1875_rkx, 0.00130_rkx/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_rkx, 0.00117_rkx/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_rkx, 0.00090_rkx/ ! CO2

  ! Vonstant (K) for Henry's law (4.12, Wania)
  real(rkx) , public :: c_h_inv(ngases)
  data c_h_inv(1:3) /1600._rkx, 1500._rkx, 2400._rkx/ ! CH4, O2, CO2
  ! Henry's constant (L.atm/mol) at standard temperature (298K)
  real(rkx) , public :: kh_theta(ngases)
  data kh_theta(1:3) /714.29_rkx, 769.23_rkx, 29.4_rkx/ ! CH4, O2, CO2
  ! Base temperature for calculation of Henry's constant (K)
  real(rkx) , public , parameter :: kh_tbase = 298._rkx

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
