module mod_clm_drydep

  !========================================================================
  ! Module for handling dry depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !
  ! !REVISION HISTORY:
  !     2008-Nov-12 - F. Vitt - creation.
  !     2009-Feb-19 - E. Kluzek - merge shr_drydep_tables module in.
  !     2009-Feb-20 - E. Kluzek - use shr_ coding standards, and check for namelist file.
  !     2009-Feb-20 - E. Kluzek - Put D0 on all constants, remove namelist read out.
  !     2009-Mar-23 - F. Vitt - Some corrections/cleanup and addition of drydep_method.
  !     2009-Mar-27 - E. Kluzek - Get description and units from J.F. Lamarque.
  !========================================================================

  ! !USES:

  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_stdio
  use mod_mpmessage

  implicit none 
  save 

  private

  save

  ! !PUBLIC MEMBER FUNCTIONS

  public :: seq_drydep_read         ! Read namelist
  public :: seq_drydep_init         ! Initialization of drydep data
  public :: seq_drydep_setHCoeff    ! Calculate Henry's law coefficients

  ! !PRIVATE ARRAY SIZES

  integer, private, parameter :: maxspc = 100              ! Maximum number of species
  integer, public,  parameter :: n_species_table = 55      ! Number of species to work with
  integer, private, parameter :: NSeas = 5                 ! Number of seasons
  integer, private, parameter :: NLUse = 11                ! Number of land-use types

  ! !PUBLIC DATA MEMBERS:

  ! method specification
  character(16),public,parameter :: DD_XATM = 'xactive_atm'! dry-dep atmosphere
  character(16),public,parameter :: DD_XLND = 'xactive_lnd'! dry-dep land
  character(16),public,parameter :: DD_TABL = 'table'      ! dry-dep table (atm and lnd)
  character(16),public :: drydep_method = DD_XLND          ! Which option choosen

  real(rk8), public, parameter :: ph     = 1.D-5         ! measure of the acidity (dimensionless)

  logical, public  :: lnd_drydep                           ! If dry-dep fields passed
  integer, public  :: n_drydep = 0                         ! Number in drypdep list
  character(len=32), public, dimension(maxspc) :: drydep_list = ''   ! List of dry-dep species

  character(len=80), public :: drydep_fields_token = ''   ! First drydep fields token

  real(rk8), public, allocatable, dimension(:) :: foxd      ! reactivity factor for oxidation (dimensioness)
  real(rk8), public, allocatable, dimension(:) :: drat      ! ratio of molecular diffusivity (D_H2O/D_species; dimensionless)
  integer,  public, allocatable, dimension(:) :: mapping   ! mapping to species table
  ! --- Indices for each species ---
  integer,  public :: h2_ndx, ch4_ndx, co_ndx, pan_ndx, mpan_ndx, so2_ndx, o3_ndx, o3a_ndx, xpan_ndx

  !---------------------------------------------------------------------------
  ! Table 1 from Wesely, Atmos. Environment, 1989, p1293
  ! Table 2 from Sheih, microfiche PB86-218104 and Walcek, Atmos.  Environment, 1986, p949
  ! Table 3-5 compiled by P. Hess
  !
  ! index #1 : season
  !           1 -> midsummer with lush vegetation
  !           2 -> autumn with unharvested cropland
  !           3 -> late autumn after frost, no snow
  !           4 -> winter, snow on ground, and subfreezing
  !           5 -> transitional spring with partially green short annuals
  !
  ! index #2 : landuse type
  !           1 -> urban land
  !           2 -> agricultural land
  !           3 -> range land
  !           4 -> deciduous forest
  !           5 -> coniferous forest
  !           6 -> mixed forest including wetland
  !           7 -> water, both salt and fresh
  !           8 -> barren land, mostly desert
  !           9 -> nonforested wetland
  !           10 -> mixed agricultural and range land
  !           11 -> rocky open areas with low growing shrubs
  !
  ! JFL August 2000
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! table to parameterize the impact of soil moisture on the deposition of H2 and
  ! CO on soils (from Sanderson et al., J. Atmos. Chem., 46, 15-28, 2003).
  !---------------------------------------------------------------------------

  !--- deposition of h2 and CO on soils ---
  real(rk8), parameter, public :: h2_a(NLUse) = &
                (/  0.000D0,  0.000D0, 0.270D0,  0.000D0,  0.000D0,  &
                    0.000D0,  0.000D0, 0.000D0,  0.000D0,  0.000D0, 0.000D0/)
  !--- deposition of h2 and CO on soils ---
  real(rk8), parameter, public :: h2_b(NLUse) = &
                (/  0.000D0,-41.390D0, -0.472D0,-41.900D0,-41.900D0,  &
                  -41.900D0,  0.000D0,  0.000D0,  0.000D0,-41.390D0,  0.000D0/)
  !--- deposition of h2 and CO on soils ---
  real(rk8), parameter, public :: h2_c(NLUse) = &
                (/  0.000D0, 16.850D0, 1.235D0, 19.700D0, 19.700D0, &
                   19.700D0,  0.000D0, 0.000D0,  0.000D0, 17.700D0, 1.000D0/)

  !--- deposition of h2 and CO on soils
  !
  !--- ri:   Richardson number                      (dimensionless)
  !--- rlu:  Resistance of leaves in upper canopy   (s.m-1)
  !--- rac:  Aerodynamic resistance to lower canopy (s.m-1)
  !--- rgss: Ground surface resistance for SO2      (s.m-1)
  !--- rgso: Ground surface resistance for O3       (s.m-1)
  !--- rcls: Lower canopy resistance for SO2        (s.m-1)
  !--- rclo: Lower canopy resistance for O3         (s.m-1)
  !
  real(rk8), public, dimension(NSeas,NLUse) :: ri, rlu, rac, rgss, rgso, rcls, rclo

  data ri  (1,1:NLUse) &
       /1.D36,  60.D0, 120.D0,  70.D0, 130.D0, 100.D0,1.D36,1.D36,  80.D0, 100.D0, 150.D0/
  data rlu (1,1:NLUse) &
       /1.D36,2000.D0,2000.D0,2000.D0,2000.D0,2000.D0,1.D36,1.D36,2500.D0,2000.D0,4000.D0/
  data rac (1,1:NLUse) &
       / 100.D0, 200.D0, 100.D0,2000.D0,2000.D0,2000.D0,   0.D0,   0.D0, 300.D0, 150.D0, 200.D0/
  data rgss(1,1:NLUse) &
       / 400.D0, 150.D0, 350.D0, 500.D0, 500.D0, 100.D0,   0.D0,1000.D0,  0.D0, 220.D0, 400.D0/
  data rgso(1,1:NLUse) &
       / 300.D0, 150.D0, 200.D0, 200.D0, 200.D0, 300.D0,2000.D0, 400.D0,1000.D0, 180.D0, 200.D0/
  data rcls(1,1:NLUse) &
       /1.D36,2000.D0,2000.D0,2000.D0,2000.D0,2000.D0,1.D36,1.D36,2500.D0,2000.D0,4000.D0/
  data rclo(1,1:NLUse) &
       /1.D36,1000.D0,1000.D0,1000.D0,1000.D0,1000.D0,1.D36,1.D36,1000.D0,1000.D0,1000.D0/

  data ri  (2,1:NLUse) &
       /1.D36,1.D36,1.D36,1.D36, 250.D0, 500.D0,1.D36,1.D36,1.D36,1.D36,1.D36/
  data rlu (2,1:NLUse) &
       /1.D36,9000.D0,9000.D0,9000.D0,4000.D0,8000.D0,1.D36,1.D36,9000.D0,9000.D0,9000.D0/
  data rac (2,1:NLUse) &
       / 100.D0, 150.D0, 100.D0,1500.D0,2000.D0,1700.D0,   0.D0,   0.D0, 200.D0, 120.D0, 140.D0/
  data rgss(2,1:NLUse) &
       / 400.D0, 200.D0, 350.D0, 500.D0, 500.D0, 100.D0,   0.D0,1000.D0,   0.D0, 300.D0, 400.D0/
  data rgso(2,1:NLUse) &
       / 300.D0, 150.D0, 200.D0, 200.D0, 200.D0, 300.D0,2000.D0, 400.D0, 800.D0, 180.D0, 200.D0/
  data rcls(2,1:NLUse) &
       /1.D36,9000.D0,9000.D0,9000.D0,2000.D0,4000.D0,1.D36,1.D36,9000.D0,9000.D0,9000.D0/
  data rclo(2,1:NLUse) &
       /1.D36, 400.D0, 400.D0, 400.D0,1000.D0, 600.D0,1.D36,1.D36, 400.D0, 400.D0, 400.D0/

  data ri  (3,1:NLUse) &
       /1.D36,1.D36,1.D36,1.D36, 250.D0, 500.D0,1.D36,1.D36,1.D36,1.D36,1.D36/
  data rlu (3,1:NLUse) &
       /1.D36,1.D36,9000.D0,9000.D0,4000.D0,8000.D0,1.D36,1.D36,9000.D0,9000.D0,9000.D0/
  data rac (3,1:NLUse) &
       / 100.D0,  10.D0, 100.D0,1000.D0,2000.D0,1500.D0,   0.D0,   0.D0, 100.D0, 50.D0, 120.D0/
  data rgss(3,1:NLUse) &
       / 400.D0, 150.D0, 350.D0, 500.D0, 500.D0, 200.D0,   0.D0,1000.D0,   0.D0, 200.D0, 400.D0/
  data rgso(3,1:NLUse) &
       / 300.D0, 150.D0, 200.D0, 200.D0, 200.D0, 300.D0,2000.D0, 400.D0,1000.D0, 180.D0, 200.D0/
  data rcls(3,1:NLUse) &
       /1.D36,1.D36,9000.D0,9000.D0,3000.D0,6000.D0,1.D36,1.D36,9000.D0,9000.D0,9000.D0/
  data rclo(3,1:NLUse) &
       /1.D36,1000.D0, 400.D0, 400.D0,1000.D0, 600.D0,1.D36,1.D36, 800.D0, 600.D0, 600.D0/

  data ri  (4,1:NLUse) &
       /1.D36,1.D36,1.D36,1.D36, 400.D0, 800.D0,1.D36,1.D36,1.D36,1.D36,1.D36/
  data rlu (4,1:NLUse) &
       /1.D36,1.D36,1.D36,1.D36,6000.D0,9000.D0,1.D36,1.D36,9000.D0,9000.D0,9000.D0/
  data rac (4,1:NLUse) &
       / 100.D0,  10.D0,  10.D0,1000.D0,2000.D0,1500.D0,   0.D0,   0.D0,  50.D0,  10.D0,  50.D0/
  data rgss(4,1:NLUse) &
       / 100.D0, 100.D0, 100.D0, 100.D0, 100.D0, 100.D0,   0.D0,1000.D0, 100.D0, 100.D0,  50.D0/
  data rgso(4,1:NLUse) &
       / 600.D0,3500.D0,3500.D0,3500.D0,3500.D0,3500.D0,2000.D0, 400.D0,3500.D0,3500.D0,3500.D0/
  data rcls(4,1:NLUse) &
       /1.D36,1.D36,1.D36,9000.D0, 200.D0, 400.D0,1.D36,1.D36,9000.D0,1.D36,9000.D0/
  data rclo(4,1:NLUse) &
       /1.D36,1000.D0,1000.D0, 400.D0,1500.D0, 600.D0,1.D36,1.D36, 800.D0,1000.D0, 800.D0/

  data ri  (5,1:NLUse) &
       /1.D36, 120.D0, 240.D0, 140.D0, 250.D0, 190.D0,1.D36,1.D36, 160.D0, 200.D0, 300.D0/
  data rlu (5,1:NLUse) &
       /1.D36,4000.D0,4000.D0,4000.D0,2000.D0,3000.D0,1.D36,1.D36,4000.D0,4000.D0,8000.D0/
  data rac (5,1:NLUse) &
       / 100.D0,  50.D0,  80.D0,1200.D0,2000.D0,1500.D0,   0.D0,   0.D0, 200.D0, 60.D0, 120.D0/
  data rgss(5,1:NLUse) &
       / 500.D0, 150.D0, 350.D0, 500.D0, 500.D0, 200.D0,   0.D0,1000.D0,   0.D0, 250.D0, 400.D0/
  data rgso(5,1:NLUse) &
       / 300.D0, 150.D0, 200.D0, 200.D0, 200.D0, 300.D0,2000.D0, 400.D0,1000.D0, 180.D0, 200.D0/
  data rcls(5,1:NLUse) &
       /1.D36,4000.D0,4000.D0,4000.D0,2000.D0,3000.D0,1.D36,1.D36,4000.D0,4000.D0,8000.D0/
  data rclo(5,1:NLUse) &
       /1.D36,1000.D0, 500.D0, 500.D0,1500.D0, 700.D0,1.D36,1.D36, 600.D0, 800.D0, 800.D0/

  !---------------------------------------------------------------------------
  !         ... roughness length
  !---------------------------------------------------------------------------
  real(rk8), public, dimension(NSeas,NLUse) :: z0

  data z0  (1,1:NLUse) &
       /1.000D0,0.250D0,0.050D0,1.000D0,1.000D0,1.000D0,0.0006D0,0.002D0,0.150D0,0.100D0,0.100D0/
  data z0  (2,1:NLUse) &
       /1.000D0,0.100D0,0.050D0,1.000D0,1.000D0,1.000D0,0.0006D0,0.002D0,0.100D0,0.080D0,0.080D0/
  data z0  (3,1:NLUse) &
       /1.000D0,0.005D0,0.050D0,1.000D0,1.000D0,1.000D0,0.0006D0,0.002D0,0.100D0,0.020D0,0.060D0/
  data z0  (4,1:NLUse) &
       /1.000D0,0.001D0,0.001D0,1.000D0,1.000D0,1.000D0,0.0006D0,0.002D0,0.001D0,0.001D0,0.040D0/
  data z0  (5,1:NLUse) &
       /1.000D0,0.030D0,0.020D0,1.000D0,1.000D0,1.000D0,0.0006D0,0.002D0,0.010D0,0.030D0,0.060D0/

  !real(rk8), private, dimension(11,5), parameter :: z0xxx = reshape ( &
  ! (/   1.000,0.250,0.050,1.000,1.000,1.000,0.0006,0.002,0.150,0.100,0.100 ,  &
  !      1.000,0.100,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.080,0.080 ,  &
  !      1.000,0.005,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.020,0.060 ,  &
  !      1.000,0.001,0.001,1.000,1.000,1.000,0.0006,0.002,0.001,0.001,0.040 ,  &
  !      1.000,0.030,0.020,1.000,1.000,1.000,0.0006,0.002,0.010,0.030,0.060  /), (/11,5/) )

  !---------------------------------------------------------------------------
  ! public chemical data
  !---------------------------------------------------------------------------

  !--- data for foxd (reactivity factor for oxidation) ----
  real(rk8), public, parameter :: dfoxd(n_species_table) = &
          (/  1.D0     &
             ,1.D0     &
             ,1.D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 &
             ,1.D0     &
             ,.1D0     &
             ,1.D-36 &
             ,0.D0     &
             ,0.D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 &
             ,1.D-36 &
             ,.1D0     &
             ,1.D0     &
             ,1.D-36 &
             ,.1D0     &
             ,1.D0     &
             ,1.D-36 &
             ,.1D0     &
             ,.1D0     &
             ,.1D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 &
             ,.1D0     &
             ,1.D-36 &
             ,.1D0     &
             ,1.D-36 &
             ,.1D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 &
             ,1.D-36 &
             ,1.D-36 &
             ,.1D0     &
             ,1.D-36 &
             ,.1D0     &
             ,1.D-36 &
             ,.1D0     &
             ,.1D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 &  
             ,1.D-36 &  
             ,1.D-36 &  
             ,1.D-36 &
             ,.1D0     &
             ,.1D0     &
             ,.1D0     &
             ,1.D-36 &
             ,1.D-36 & ! HCN
             ,1.D-36 & ! CH3CN
            /)
!EOP

! PRIVATE DATA:

  Interface seq_drydep_setHCoeff                      ! overload subroutine
     Module Procedure set_hcoeff_scalar
     Module Procedure set_hcoeff_vector
  End Interface

  real(rk8), private, parameter :: small_value = 1.D-36          !--- smallest value to use ---

  !---------------------------------------------------------------------------
  ! private chemical data
  !---------------------------------------------------------------------------

  !--- Names of species that can work with ---
  character(len=20), public, parameter :: species_name_table(n_species_table) = &
                         (/ 'OX      '                       & 
                           ,'H2O2    '                       & 
                           ,'OH      '                       & 
                           ,'HO2     '                       & 
                           ,'CO      '                       & 
                           ,'CH4     '                       & 
                           ,'CH3O2   '                       & 
                           ,'CH3OOH  '                       & 
                           ,'CH2O    '                       & 
                           ,'CHOOH   '                       & 
                           ,'NO      '                       & 
                           ,'NO2     '                       & 
                           ,'HNO3    '                       & 
                           ,'CO2     '                       & 
                           ,'NH3     '                       & 
                           ,'N2O5    '                       & 
                           ,'NO3     '                       & 
                           ,'CH3OH   '                       & 
                           ,'HO2NO2  '                       & 
                           ,'O1D     '                       & 
                           ,'C2H6    '                       & 
                           ,'C2H5O2  '                       & 
                           ,'PO2     '                       & 
                           ,'MACRO2  '                       & 
                           ,'ISOPO2  '                       & 
                           ,'C4H10   '                       & 
                           ,'CH3CHO  '                       & 
                           ,'C2H5OOH '                       & 
                           ,'C3H6    '                       & 
                           ,'POOH    '                       & 
                           ,'C2H4    '                       & 
                           ,'PAN     '                       & 
                           ,'CH3COOOH'                       & 
                           ,'C10H16  '                       & 
                           ,'CHOCHO  '                       & 
                           ,'CH3COCHO'                       & 
                           ,'GLYALD  '                       & 
                           ,'CH3CO3  '                       & 
                           ,'C3H8    '                       & 
                           ,'C3H7O2  '                       & 
                           ,'CH3COCH3'                       & 
                           ,'C3H7OOH '                       & 
                           ,'RO2     '                       & 
                           ,'ROOH    '                       & 
                           ,'Rn      '                       & 
                           ,'ISOP    '                       &
                           ,'MVK     '                       &
                           ,'MACR    '                       &
                           ,'C2H5OH  '                       &
                           ,'ONITR   '                       & 
                           ,'ONIT    '                       & 
                           ,'ISOPNO3 '                       & 
                           ,'HYDRALD '                       &
                           ,'HCN     '                       &
                           ,'CH3CN   '                       &
                           /)

  !--- data for effective Henry's Law coefficient ---
  real(rk8), public, parameter :: dheff(n_species_table*6) = &
            (/1.15D-02, 2560.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,8.33D+04, 7379.D0,2.2D-12,-3730.D0,0.D0     ,    0.D0  &
             ,3.00D+01,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,2.00D+03, 6600.D0,3.5D-05,    0.D0,0.D0     ,    0.D0  &
             ,1.00D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.11D+02, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,6.30D+03, 6425.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,5.53D+03, 5700.D0,1.8D-04,-1510.D0,0.D0     ,    0.D0  &
             ,1.90D-03, 1480.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,6.40D-03, 2500.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,0.D0      ,    0.D0,2.6D+06, 8700.D0,0.D0     ,    0.D0  &
             ,3.40D-02, 2420.D0,4.5D-07,-1000.D0,3.6D-11,-1760.D0  &
             ,7.40D+01, 3400.D0,1.7D-05, -450.D0,1.0D-14,-6716.D0  &
             ,2.14D+00, 3362.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,0.65D+00,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,2.20D+02, 4934.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,0.D0      ,    0.D0,3.2D+01,    0.D0,0.D0     ,    0.D0  &
             ,1.00D-16,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &  
             ,1.14D+01, 6267.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.36D+02, 5995.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,2.20D+02, 5653.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,5.00D+00,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,8.37D+02, 5308.D0,1.8D-04,-1510.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.00D+05,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.71D+03, 7541.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,4.14D+04, 4630.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.45D-03, 2700.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.00D+06,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,2.70D+01, 5300.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.36D+02, 5995.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &  
             ,7.47D+00, 5241.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,3.36D+02, 5995.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,0.00D+00,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.70D-03,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,2.00D+02, 6500.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.51D+03, 6485.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.00D+03, 6000.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.00D+01,    0.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,7.00D+01, 6000.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,1.20D+01, 5000.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
             ,5.00D+01, 4000.D0,0.D0     ,    0.D0,0.D0     ,    0.D0  &
            /)

  real(rk8), private, parameter :: wh2o = amw
  real(rk8), private, parameter :: mol_wgts(n_species_table) = &
       (/ 47.9981995D0, 34.0135994D0, 17.0067997D0, 33.0061989D0, 28.0104008D0, &
          16.0405998D0, 47.0320015D0, 48.0393982D0, 30.0251999D0, 46.0246010D0, &
          30.0061398D0, 46.0055389D0, 63.0123405D0, 44.0098000D0, 17.0289402D0, &
          108.010483D0, 62.0049400D0, 32.0400009D0, 79.0117416D0, 15.9994001D0, &
          30.0664005D0, 61.0578003D0, 91.0830002D0, 119.093399D0, 117.119797D0, &
          58.1180000D0, 44.0509987D0, 62.0652008D0, 42.0774002D0, 92.0904007D0, &
          28.0515995D0, 121.047943D0, 76.0497971D0, 136.228394D0, 58.0355988D0, &
          72.0614014D0, 60.0503998D0, 75.0423965D0, 44.0922012D0, 75.0836029D0, &
          58.0768013D0, 76.0910034D0, 31.9988003D0, 33.0061989D0, 222.000000D0, &
          68.1141968D0, 70.0877991D0, 70.0877991D0, 46.0657997D0, 147.125946D0, &
          119.074341D0, 162.117935D0, 100.112999D0, 27.0256D0   , 41.0524D0  /)

!===============================================================================
CONTAINS
!===============================================================================

!====================================================================================

  subroutine seq_drydep_read(NLFilename, seq_drydep_fields)

    !========================================================================
    ! reads drydep_inparm namelist and sets up CCSM driver list of fields for 
    ! land-atmosphere communications.
    ! 
    ! !REVISION HISTORY:
    !  2009-Feb-20 - E. Kluzek - Separate out as subroutine from previous input_init
    !========================================================================

    implicit none

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    character(len=*), intent(out) :: seq_drydep_fields  

    !----- local -----
    integer :: i                ! Indices
    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    character(len=8) :: token   ! dry dep field name to add

    !----- formats -----
    character(*),parameter :: subName = '(seq_drydep_read) ' 
    character(*),parameter :: F00   = "('(seq_drydep_read) ',8a)" 
    character(*),parameter :: FI1   = "('(seq_drydep_init) ',a,I2)" 

    namelist /drydep_inparm/ drydep_list, drydep_method

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the drydep field list to pass
    ! First check if file exists and if not, n_drydep will be zero
    !-----------------------------------------------------------------------------

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0  )then
       call fatal(__FILE__,__LINE__,subName//'ERROR: nlfilename not set' )
    end if
    inquire( file=trim(NLFileName), exist=exists)
    if ( exists ) then
       unitn = file_getUnit()
       open( unitn, file=trim(NLFilename), status='old' )
       if ( debug_level > 0 ) write(stdout,F00) &
            'Read in drydep_inparm namelist from: ', trim(NLFilename)
       ierr = 1
       do while ( ierr /= 0 )
          read(unitn, drydep_inparm, iostat=ierr)
          if (ierr < 0) then
             call fatal(__FILE__,__LINE__, &
               subName//'ERROR: encountered end-of-file on namelist read' )
          endif
       end do
       close( unitn )
       call file_freeUnit( unitn )
    end if

    n_drydep = 0

    !--- Loop over species to fill list of fields to communicate for drydep ---
    seq_drydep_fields = ' '  
    do i=1,maxspc
       if ( len_trim(drydep_list(i))==0 ) exit
       write(token,333) i
       seq_drydep_fields = trim(seq_drydep_fields)//':'//trim(token)                 
       if ( i == 1 ) then
          seq_drydep_fields = trim(token)                 
          drydep_fields_token = trim(token)
       endif
       n_drydep = n_drydep+1
    enddo

    !--- Make sure method is valid and determine if land is passing drydep fields ---
    lnd_drydep = n_drydep>0 .and. drydep_method == DD_XLND

    if ( debug_level > 0 ) then
       write(stdout,*) 'seq_drydep_read: drydep_method: ', trim(drydep_method)
       if ( n_drydep == 0 )then
          write(stdout,F00) 'No dry deposition fields will be transfered'
       else
          write(stdout,FI1) 'Number of dry deposition fields transfered is ', &
                                n_drydep
       end if
    end if

    if ( trim(drydep_method)/=trim(DD_XATM) .and. &
         trim(drydep_method)/=trim(DD_XLND) .and. &
         trim(drydep_method)/=trim(DD_TABL) ) then
       if ( debug_level > 0 ) then
          write(stdout,*) 'seq_drydep_read: drydep_method : ', trim(drydep_method)
          write(stdout,*) 'seq_drydep_read: drydep_method must be set to : ', &
               DD_XATM,', ', DD_XLND,', or ', DD_TABL
       end if
       call fatal(__FILE__,__LINE__, &
       'seq_drydep_read: incorrect dry deposition method specification')
    endif

    ! Need to explicitly add Sl_ based on naming convention
333 format ('Sl_dd',i3.3)

  end subroutine seq_drydep_read

!====================================================================================

  subroutine seq_drydep_init( )

    !========================================================================
    ! Initialization of dry deposition fields
    ! reads drydep_inparm namelist and sets up CCSM driver list of fields for 
    ! land-atmosphere communications.
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !  2009-Feb-20 - E. Kluzek - Check for existance of file if not return, set n_drydep=0
    !  2009-Feb-20 - E. Kluzek - Move namelist read to separate subroutine
    !========================================================================

    implicit none

    !----- local -----
    integer :: i, l                      ! Indices
    character(len=32) :: test_name       ! field test name
    !----- formats -----
    character(*),parameter :: subName = '(seq_drydep_init) ' 
    character(*),parameter :: F00   = "('(seq_drydep_init) ',8a)" 

    !-----------------------------------------------------------------------------
    ! Allocate and fill foxd, drat and mapping as well as species indices
    !-----------------------------------------------------------------------------

    if ( n_drydep > 0 ) then

       allocate( foxd(n_drydep) )
       allocate( drat(n_drydep) )
       allocate( mapping(n_drydep) )

       ! This initializes these variables to infinity.
       foxd = inf
       drat = inf

       mapping(:) = 0

    end if

    h2_ndx=-1; ch4_ndx=-1; co_ndx=-1; mpan_ndx = -1; pan_ndx = -1; so2_ndx=-1; o3_ndx=-1; xpan_ndx=-1

    !--- Loop over drydep species that need to be worked with ---
    do i=1,n_drydep
       if ( len_trim(drydep_list(i))==0 ) exit

       test_name = drydep_list(i)

       if( trim(test_name) == 'O3' ) then
          test_name = 'OX'
       end if

       !--- Figure out if species maps to a species in the species table ---
       do l = 1,n_species_table
          if(  trim( test_name ) == trim( species_name_table(l) ) ) then
             mapping(i)  = l
             exit
          end if
       end do

       !--- If it doesn't map to a species in the species table find species close enough ---
       if( mapping(i) < 1 ) then
          select case( trim(test_name) )
          case( 'H2' )
             test_name = 'CO'
          case( 'HYAC', 'CH3COOH', 'EOOH' )
             test_name = 'CH2O'
          case( 'O3S', 'O3INERT', 'MPAN' )
             test_name = 'OX'
          case( 'ISOPOOH', 'MACROOH', 'Pb', 'XOOH', 'H2SO4' )
             test_name = 'HNO3'
          case( 'ALKOOH', 'MEKOOH', 'TOLOOH', 'BENOOH', 'XYLOOH', 'TERPOOH','SOGM','SOGI','SOGT','SOGB','SOGX' )
                test_name = 'CH3OOH'
          case( 'SOA', 'SO2', 'SO4', 'CB1', 'CB2', 'OC1', 'OC2', 'NH3', 'NH4', 'SA1', 'SA2', 'SA3', 'SA4','HCN','CH3CN','HCOOH' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'O3A', 'XMPAN' )
             test_name = 'OX'
          case( 'XPAN' )
             test_name = 'PAN'
          case( 'XNO' )
             test_name = 'NO'
          case( 'XNO2' )
             test_name = 'NO2'
          case( 'XHNO3' )
             test_name = 'HNO3'
          case( 'XONIT' )
             test_name = 'ONIT'
          case( 'XONITR' )
             test_name = 'ONITR'
          case( 'XHO2NO2')
             test_name = 'HO2NO2'
          case( 'XNH4NO3' )
             test_name = 'HNO3'

          case( 'COhc','COme')
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd
          case( 'CO01','CO02','CO03','CO04','CO05','CO06','CO07','CO08','CO09','CO10' )
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd
          case( 'CO11','CO12','CO13','CO14','CO15','CO16','CO17','CO18','CO19','CO20' )
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd
          case( 'CO21','CO22','CO23','CO24','CO25','CO26','CO27','CO28','CO29','CO30' )
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd
          case( 'CO31','CO32','CO33','CO34','CO35','CO36','CO37','CO38','CO39','CO40' )
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd
          case( 'CO41','CO42','CO43','CO44','CO45','CO46','CO47','CO48','CO49','CO50' )
                test_name = 'CO'  ! this is just a place holder. values are set in drydep_fromlnd

          case( 'NH4NO3' )
             test_name = 'HNO3'
          case default
             test_name = 'blank'
          end select

          !--- If found a match check the species table again ---
          if( trim(test_name) /= 'blank' ) then
             do l = 1,n_species_table
                if( trim( test_name ) == trim( species_name_table(l) ) ) then
                   mapping(i)  = l
                   exit
                end if
             end do
          else
             if ( debug_level > 0 ) write(stdout,F00) trim(drydep_list(i)), &
                                 ' not in tables; will have dep vel = 0'
             call fatal(__FILE__,__LINE__, &
              subName//': '//trim(drydep_list(i))//' is not in tables' )
          end if
       end if

       !--- Figure out the specific species indices ---
       if ( trim(drydep_list(i)) == 'H2' )   h2_ndx   = i
       if ( trim(drydep_list(i)) == 'CO' )   co_ndx   = i
       if ( trim(drydep_list(i)) == 'CH4' )  ch4_ndx  = i
       if ( trim(drydep_list(i)) == 'MPAN' ) mpan_ndx = i
       if ( trim(drydep_list(i)) == 'PAN' )  pan_ndx  = i
       if ( trim(drydep_list(i)) == 'SO2' )  so2_ndx  = i
       if ( trim(drydep_list(i)) == 'OX' .or. trim(drydep_list(i)) == 'O3' ) o3_ndx  = i 
       if ( trim(drydep_list(i)) == 'O3A' ) o3a_ndx  = i 
       if ( trim(drydep_list(i)) == 'XPAN' ) xpan_ndx = i

       if( mapping(i) > 0) then
         l = mapping(i)
         foxd(i) = dfoxd(l)
         drat(i) = sqrt(mol_wgts(l)/wh2o)
       endif

    enddo

    where( rgss < 1.D0 )
       rgss = 1.D0
    endwhere

    where( rac < small_value)
       rac = small_value
    endwhere

  end subroutine seq_drydep_init

!====================================================================================

  subroutine set_hcoeff_scalar( sfc_temp, heff )

    !========================================================================
    ! Interface to seq_drydep_setHCoeff when input is scalar
    ! wrapper routine used when surface temperature is a scalar (single column) rather
    ! than an array (multiple columns).
    !
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !========================================================================

    implicit none

    real(rk8), intent(in)     :: sfc_temp         ! Input surface temperature
    real(rk8), intent(out)    :: heff(n_drydep)   ! Output Henry's law coefficients

    !----- local -----
    real(rk8) :: sfc_temp_tmp(1)    ! surface temp

    sfc_temp_tmp(:) = sfc_temp
    call set_hcoeff_vector( 1, sfc_temp_tmp, heff(:n_drydep) )

  end subroutine set_hcoeff_scalar

!====================================================================================

  subroutine set_hcoeff_vector( ncol, sfc_temp, heff )

    !========================================================================
    ! Interface to seq_drydep_setHCoeff when input is vector
    ! sets dry depositions coefficients -- used by both land and atmosphere models
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !========================================================================

    implicit none

    integer, intent(in)      :: ncol                  ! Input size of surface-temp vector
    real(rk8), intent(in)     :: sfc_temp(ncol)        ! Surface temperature
    real(rk8), intent(out)    :: heff(ncol,n_drydep)   ! Henry's law coefficients

    !----- local -----
    real(rk8), parameter :: t0     = 298.D0    ! Standard Temperature
    real(rk8), parameter :: ph_inv = 1.D0/ph   ! Inverse of PH
    integer  :: m, l, id       ! indices
    real(rk8) :: e298           ! Henry's law coefficient @ standard temperature (298K)
    real(rk8) :: dhr            ! temperature dependence of Henry's law coefficient
    real(rk8) :: dk1s(ncol)     ! DK Work array 1
    real(rk8) :: dk2s(ncol)     ! DK Work array 2
    real(rk8) :: wrk(ncol)      ! Work array

    !----- formats -----
    character(*),parameter :: subName = '(seq_drydep_set_hcoeff) ' 
    character(*),parameter :: F00   = "('(seq_drydep_set_hcoeff) ',8a)" 

    !-------------------------------------------------------------------------------
    ! notes: 
    !-------------------------------------------------------------------------------

    wrk(:) = (t0 - sfc_temp(:))/(t0*sfc_temp(:))
    do m = 1,n_drydep
       l    = mapping(m)
       id   = 6*(l - 1)
       e298 = dheff(id+1)
       dhr  = dheff(id+2)
       heff(:,m) = e298*exp( dhr*wrk(:) )
       !--- Calculate coefficients based on the drydep tables ---
       if( dheff(id+3) /= 0.D0 .and. dheff(id+5) == 0.D0 ) then
          e298 = dheff(id+3)
          dhr  = dheff(id+4)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          where( heff(:,m) /= 0.D0 )
             heff(:,m) = heff(:,m)*(1.D0 + dk1s(:)*ph_inv)
          elsewhere
             heff(:,m) = dk1s(:)*ph_inv
          endwhere
       end if
       !--- For coefficients that are non-zero AND CO2 or NH3 handle things this way ---
       if( dheff(id+5) /= 0.D0 ) then
          if( trim( drydep_list(m) ) == 'CO2' .or. trim( drydep_list(m) ) == 'NH3' ) then
             e298 = dheff(id+3)
             dhr  = dheff(id+4)
             dk1s(:) = e298*exp( dhr*wrk(:) )
             e298 = dheff(id+5)
             dhr  = dheff(id+6)
             dk2s(:) = e298*exp( dhr*wrk(:) )
             !--- For Carbon dioxide ---
             if( trim(drydep_list(m)) == 'CO2' ) then
                heff(:,m) = heff(:,m)*(1.D0 + dk1s(:)*ph_inv)*(1.D0 + dk2s(:)*ph_inv)
             !--- For NH3 ---
             else if( trim( drydep_list(m) ) == 'NH3' ) then
                heff(:,m) = heff(:,m)*(1.D0 + dk1s(:)*ph/dk2s(:))
             !--- This can't happen ---
             else
                write(stdout,F00) 'Bad species ',drydep_list(m)
                call fatal(__FILE__,__LINE__, &
                 subName//'ERROR: in assigning coefficients' )
             end if
          end if
       end if
    end do

  end subroutine set_hcoeff_vector

!===============================================================================

end module mod_clm_drydep
