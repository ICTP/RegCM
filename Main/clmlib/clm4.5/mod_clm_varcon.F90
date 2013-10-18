module mod_clm_varcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varcon
!
! !DESCRIPTION:
! Module containing various model constants
!
! !USES:
  use mod_realkinds 
  use mod_constants , only : cpd , pdbratio , secpd , earthrad , tzero , &
    rgasmol , rhoice , rhoh2o , wlhf , wlhs , wlhv , spcpice , spcpfw ,  &
    egrav , rwat , vonkar , sigm , mathpi , rdry
  use mod_clm_varpar , only : numrad , nlevgrnd , nlevlak , &
     nlevdecomp_full , nsompools , ngases
#if (defined VICHYDRO)
  use mod_clm_varpar   , only: nlayer
#endif

!
! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
!
!EOP
!-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(rk8) :: rpi    = mathpi

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(rk8), parameter :: n_melt=0.7     !fsca shape parameter
  real(rk8), parameter :: e_ice=6.0      !soil ice impedance factor
  real(rk8), parameter :: pc = 0.4       !threshold probability
  real(rk8), parameter :: mu = 0.13889   !connectivity exponent 
  real(rk8) :: grav   = egrav            !gravity constant [m/s2]
  real(rk8) :: sb     = sigm             !stefan-boltzmann constant  [W/m2/K4]
  real(rk8) :: vkc    = vonkar           !von Karman constant [-]
  real(rk8) :: rair   = rdry             !gas constant for dry air [J/kg/K]
  real(rk8) :: roverg = rwat/egrav*1000.D0 !Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(rk8) :: cpliq  = spcpfw           !Specific heat of water [J/kg-K]
  real(rk8) :: cpice  = spcpice          !Specific heat of ice [J/kg-K]
  real(rk8) :: cpair  = cpd              !specific heat of dry air [J/kg/K]
  real(rk8) :: hvap   = wlhv             !Latent heat of evap for water [J/kg]
  real(rk8) :: hsub   = wlhs             !Latent heat of sublimation    [J/kg]
  real(rk8) :: hfus   = wlhf             !Latent heat of fusion for ice [J/kg]
  real(rk8) :: denh2o = rhoh2o           !density of liquid water [kg/m3]
  real(rk8) :: denice = rhoice           !density of ice [kg/m3]
  real(rk8) :: rgas   = rgasmol*1000.0D0 !universal gas constant [J/K/kmole]
  real(rk8) :: tkair  = 0.023D0          !thermal conductivity of air   [W/m/K]
  real(rk8) :: tkice  = 2.290D0          !thermal conductivity of ice   [W/m/K]
  real(rk8) :: tkwat  = 0.57D0           !thermal conductivity of water [W/m/K]
  real(rk8) :: tfrz   = tzero            !freezing temperature [K]
  real(rk8) :: tcrit  = 2.5D0       !critical temperature to determine rain or snow
  real(rk8) :: o2_molar_const = 0.209D0   !constant atmospheric O2 molar ratio (mol/mol)

  real(rk8) :: bdsno = 250.D0       !bulk density snow (kg/m**3)
  real(rk8) :: alpha_aero = 1.0D0   !constant for aerodynamic parameter weighting
  real(rk8) :: tlsai_crit = 2.0D0   !critical value of elai+esai for which aerodynamic parameters are maximum
  real(rk8) :: watmin = 0.01D0      !minimum soil moisture (mm)

  real(rk8) :: re = earthrad*0.001D0 !radius of earth (km)

  real(rk8), public, parameter :: degpsec = 15.D0/3600.0D0 ! Degree's earth rotates per second

  real(rk8), public, parameter ::  secspday= secpd    ! Seconds per day
  integer,  public, parameter :: isecspday= secspday  ! Integer seconds per day
  real(rk8), public, parameter ::  spval = 1.D36  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data

  ! These are tunable constants from clm2_3

  real(rk8) :: zlnd = 0.01D0      !Roughness length for soil [m]
  real(rk8) :: zsno = 0.0024D0    !Roughness length for snow [m]
  real(rk8) :: csoilc = 0.004D0   !Drag coefficient for soil under canopy [-]
  real(rk8) :: capr   = 0.34D0    !Tuning factor to turn first layer T into surface T
  real(rk8) :: cnfac  = 0.5D0     !Crank Nicholson factor between 0 and 1
  real(rk8) :: ssi    = 0.033D0   !Irreducible water saturation of snow
  real(rk8) :: wimp   = 0.05D0    !Water impremeable if porosity less than wimp
  real(rk8) :: pondmx = 0.0D0     !Ponding depth (mm)
  real(rk8) :: pondmx_urban = 1.0D0  !Ponding depth for urban roof and impervious road (mm)

  real(rk8) :: thk_bedrock = 3.0D0      ! thermal conductivity of 'typical' saturated granitic rock 
                                        ! (Clauser and Huenges, 1995)(W/m/K)

  !!! C13
  real(rk8), parameter :: preind_atm_del13c = -6.0   ! preindustrial value for atmospheric del13C
  real(rk8), parameter :: preind_atm_ratio = pdbratio + (preind_atm_del13c * pdbratio)/1000.0  ! 13C/12C
  real(rk8) :: c13ratio = preind_atm_ratio/(1.0+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  !!! C14
  real(rk8) :: c14ratio = 1.D-12
  ! real(rk8) :: c14ratio = 1.D0  ! debug lets set to 1 to try to avoid numerical errors

  ! Note that the wasteheat factors are currently set to zero until a better parameterization can be developed
  ! The prior parameterization appeared to be significantly overestimating wasteheat
  real(rk8) :: ht_wasteheat_factor = 0.0D0  !wasteheat factor for urban heating (-)
  real(rk8) :: ac_wasteheat_factor = 0.0D0  !wasteheat factor for urban air conditioning (-)
  real(rk8) :: wasteheat_limit = 100.D0  !limit on wasteheat (W/m2)

  real(rk8), parameter :: h2osno_max = 1000.D0    ! max allowed snow thickness (mm H2O)
  real(rk8), parameter :: lapse_glcmec = 0.006D0  ! surface temperature lapse rate (deg m-1)
                                                  ! Pritchard et al. (GRL, 35, 2008) use 0.006  
  integer, private :: i  ! loop index


#ifdef NITRIF_DENITRIF
 !  real(rk8), parameter :: nitrif_n2o_loss_frac = 0.02D0   !fraction of N lost as N2O in nitrification (Parton et al., 2001)
  real(rk8), parameter :: nitrif_n2o_loss_frac = 6.D-4   !fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(rk8), parameter :: frac_minrlztn_to_no3 = 0.2D0   !fraction of N mineralized that is dieverted to the nitrification stream (Parton et al., 2001)
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

  integer, parameter :: istsoil    = 1  !soil         landunit type
  integer, parameter :: istice     = 2  !land ice     landunit type
  integer, parameter :: istdlak    = 3  !deep lake    landunit type
  ! Not used; now 3 is used for all lakes, which have variable depth.
  integer, parameter :: istslak    = 4  !shallow lake landunit type
  integer, parameter :: istwet     = 5  !wetland      landunit type
  integer, parameter :: isturb     = 6  !urban        landunit type
  integer, parameter :: istice_mec = 7  !land ice (multiple elevation classes) landunit type
  integer, parameter :: istcrop    = 8  !crop         landunit type
  integer, parameter :: max_lunit  = 8  !maximum value that clm3%g%l%itype can have
                             !(i.e., largest value in the above list)

  ! urban column types

  integer, parameter :: icol_roof        = 61
  integer, parameter :: icol_sunwall     = 62
  integer, parameter :: icol_shadewall   = 63
  integer, parameter :: icol_road_imperv = 64
  integer, parameter :: icol_road_perv   = 65

  ! urban density types

  integer, parameter :: udens_base     = 600
  integer, parameter :: udens_tbd      = 601
  integer, parameter :: udens_hd       = 602
  integer, parameter :: udens_md       = 603

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------


  real(rk8), allocatable :: albsat(:,:) ! wet soil albedo by color class and waveband (1=vis,2=nir)
  real(rk8), allocatable :: albdry(:,:) ! dry soil albedo by color class and waveband (1=vis,2=nir)
  real(rk8) :: betads  = 0.5D0            ! two-stream parameter betad for snow
  real(rk8) :: betais  = 0.5D0            ! two-stream parameter betai for snow
  real(rk8) :: omegas(numrad)           ! two-stream parameter omega for snow by band
  data (omegas(i),i=1,numrad) /0.8D0, 0.4D0/

  ! Lake Model Constants will be defined in SLakeCon.

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(rk8), pointer :: zlak(:)         !lake z  (layers)
  real(rk8), pointer :: dzlak(:)        !lake dz (thickness)
  real(rk8), pointer :: zsoi(:)         !soil z  (layers)
  real(rk8), pointer :: dzsoi(:)        !soil dz (thickness)
  real(rk8), pointer :: zisoi(:)        !soil zi (interfaces)
  real(rk8), pointer :: dzsoi_decomp(:) !soil dz (thickness)
#if (defined VICHYDRO)
  integer, pointer  :: nlvic(:)        !number of CLM layers in each VIC layer (#)
#endif

  !------------------------------------------------------------------
  ! (Non-tunable) Constants for the CH4 submodel (Tuneable constants in ch4varcon)
  !------------------------------------------------------------------
  ! Note some of these constants are also used in CNNitrifDenitrifMod

  real(rk8), parameter :: catomw = 12.011D0 ! molar mass of C atoms (g/mol)

  real(rk8) :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
  data (s_con(1,i),i=1,4) /1898D0, -110.1D0, 2.834D0, -0.02791D0/ ! CH4
  data (s_con(2,i),i=1,4) /1801D0, -120.1D0, 3.7818D0, -0.047608D0/ ! O2
  data (s_con(3,i),i=1,4) /1911D0, -113.7D0, 2.967D0, -0.02943D0/ ! CO2

  real(rk8) :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
  data (d_con_w(1,i),i=1,3) /0.9798D0, 0.02986D0, 0.0004381D0/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.172D0, 0.03443D0, 0.0005048D0/ ! O2
  data (d_con_w(3,i),i=1,3) /0.939D0, 0.02671D0, 0.0004095D0/ ! CO2

  real(rk8) :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  data (d_con_g(1,i),i=1,2) /0.1875D0, 0.0013D0/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759D0, 0.00117D0/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325D0, 0.0009D0/ ! CO2

  real(rk8) :: c_h_inv(ngases)    ! constant (K) for Henry's law (4.12, Wania)
  data c_h_inv(1:3) /1600.D0, 1500.D0, 2400.D0/ ! CH4, O2, CO2
  real(rk8) :: kh_theta(ngases)    ! Henry's constant (L.atm/mol) at standard temperature (298K)
  data kh_theta(1:3) /714.29D0, 769.23D0, 29.4D0/ ! CH4, O2, CO2
  real(rk8) :: kh_tbase = 298.D0 ! base temperature for calculation of Henry's constant (K)

! !PUBLIC MEMBER FUNCTIONS:
  public clm_varcon_init          ! Initialze constants that need to be initialized

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!-----------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varcon_init
!
! !INTERFACE:
  subroutine clm_varcon_init()
!
! !DESCRIPTION:
! This subroutine initializes constants in clm_varcon. MUST be called 
! after the clm_varpar_init.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by E. Kluzek
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------
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
