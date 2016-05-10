module mod_clm_drydep
  !
  ! Module for handling dry depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_stdio
  use mod_mpmessage

  implicit none

  private

  save

  public :: seq_drydep_read         ! Read namelist
  public :: seq_drydep_init         ! Initialization of drydep data
  public :: seq_drydep_setHCoeff    ! Calculate Henry's law coefficients

  ! Maximum number of species
  integer(ik4), private, parameter :: maxspc = 100
  ! Number of species to work with
  integer(ik4), public,  parameter :: n_species_table = 55
  integer(ik4), private, parameter :: NSeas = 5  ! Number of seasons
  integer(ik4), private, parameter :: NLUse = 11 ! Number of land-use types

  ! method specification
  ! dry-dep atmosphere
  character(16) , public , parameter :: DD_XATM = 'xactive_atm'
  ! dry-dep land
  character(16) , public , parameter :: DD_XLND = 'xactive_lnd'
  ! dry-dep table (atm and lnd)
  character(16) , public , parameter :: DD_TABL = 'table'

  ! Which option choosen
  character(16) , public :: drydep_method = DD_XLND

  ! measure of the acidity (dimensionless)
  real(rkx), public, parameter :: ph     = 1.e-5_rkx

  logical, public  :: lnd_drydep  ! If dry-dep fields passed
  integer(ik4), public  :: n_drydep = 0 ! Number in drypdep list

  ! List of dry-dep species
  character(len=32), public, dimension(maxspc) :: drydep_list = ''

  ! First drydep fields token
  character(len=80), public :: drydep_fields_token = ''

  ! reactivity factor for oxidation (dimensioness)
  real(rkx), public, allocatable, dimension(:) :: foxd
  ! ratio of molecular diffusivity (D_H2O/D_species; dimensionless)
  real(rkx), public, allocatable, dimension(:) :: drat
  ! mapping to species table
  integer(ik4),  public, allocatable, dimension(:) :: mapping
  ! --- Indices for each species ---
  integer(ik4),  public :: h2_ndx, ch4_ndx, co_ndx, pan_ndx, &
          mpan_ndx, so2_ndx, o3_ndx, o3a_ndx, xpan_ndx

  !---------------------------------------------------------------------------
  ! Table 1 from Wesely, Atmos. Environment, 1989, p1293
  ! Table 2 from Sheih, microfiche PB86-218104 and
  !              Walcek, Atmos.  Environment, 1986, p949
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
  ! table to parameterize the impact of soil moisture on the deposition
  ! of H2 and CO on soils
  ! (from Sanderson et al., J. Atmos. Chem., 46, 15-28, 2003).
  !---------------------------------------------------------------------------

  !--- deposition of h2 and CO on soils ---
  real(rkx), parameter, public :: h2_a(NLUse) = &
                (/  0.000_rkx,  0.000_rkx, 0.270_rkx,  0.000_rkx,  0.000_rkx,  &
                    0.000_rkx,  0.000_rkx, 0.000_rkx,  0.000_rkx,  0.000_rkx, 0.000_rkx/)
  !--- deposition of h2 and CO on soils ---
  real(rkx), parameter, public :: h2_b(NLUse) = &
                (/  0.000_rkx,-41.390_rkx, -0.472_rkx,-41.900_rkx,-41.900_rkx,  &
                  -41.900_rkx,  0.000_rkx,  0.000_rkx,  0.000_rkx,-41.390_rkx,  0.000_rkx/)
  !--- deposition of h2 and CO on soils ---
  real(rkx), parameter, public :: h2_c(NLUse) = &
                (/  0.000_rkx, 16.850_rkx, 1.235_rkx, 19.700_rkx, 19.700_rkx, &
                   19.700_rkx,  0.000_rkx, 0.000_rkx,  0.000_rkx, 17.700_rkx, 1.000_rkx/)

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
  real(rkx), public, dimension(NSeas,NLUse) :: ri, rlu, rac, &
          rgss, rgso, rcls, rclo

  data ri  (1,1:NLUse) &
       /1.e36_rkx,  60._rkx, 120._rkx,  70._rkx, 130._rkx, 100._rkx,1.e36_rkx,1.e36_rkx,  &
        80._rkx, 100._rkx, 150._rkx/
  data rlu (1,1:NLUse) &
       /1.e36_rkx,2000._rkx,2000._rkx,2000._rkx,2000._rkx,2000._rkx,1.e36_rkx,1.e36_rkx,&
        2500._rkx,2000._rkx,4000._rkx/
  data rac (1,1:NLUse) &
       / 100._rkx, 200._rkx, 100._rkx,2000._rkx,2000._rkx,2000._rkx,   0._rkx,   &
         0._rkx, 300._rkx, 150._rkx, 200._rkx/
  data rgss(1,1:NLUse) &
       / 400._rkx, 150._rkx, 350._rkx, 500._rkx, 500._rkx, 100._rkx,   0._rkx, &
         1000._rkx,  0._rkx, 220._rkx, 400._rkx/
  data rgso(1,1:NLUse) &
       / 300._rkx, 150._rkx, 200._rkx, 200._rkx, 200._rkx, 300._rkx,2000._rkx, &
         400._rkx,1000._rkx, 180._rkx, 200._rkx/
  data rcls(1,1:NLUse) &
       /1.e36_rkx,2000._rkx,2000._rkx,2000._rkx,2000._rkx,2000._rkx,1.e36_rkx,1.e36_rkx, &
        2500._rkx,2000._rkx,4000._rkx/
  data rclo(1,1:NLUse) &
       /1.e36_rkx,1000._rkx,1000._rkx,1000._rkx,1000._rkx,1000._rkx,1.e36_rkx,1.e36_rkx,&
        1000._rkx,1000._rkx,1000._rkx/

  data ri  (2,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx, 250._rkx, 500._rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx/
  data rlu (2,1:NLUse) &
       /1.e36_rkx,9000._rkx,9000._rkx,9000._rkx,4000._rkx,8000._rkx,1.e36_rkx,1.e36_rkx,&
       9000._rkx,9000._rkx,9000._rkx/
  data rac (2,1:NLUse) &
       / 100._rkx, 150._rkx, 100._rkx,1500._rkx,2000._rkx,1700._rkx,   0._rkx,  &
         0._rkx, 200._rkx, 120._rkx, 140._rkx/
  data rgss(2,1:NLUse) &
       / 400._rkx, 200._rkx, 350._rkx, 500._rkx, 500._rkx, 100._rkx,   0._rkx, &
        1000._rkx,   0._rkx, 300._rkx, 400._rkx/
  data rgso(2,1:NLUse) &
       / 300._rkx, 150._rkx, 200._rkx, 200._rkx, 200._rkx, 300._rkx,2000._rkx, &
         400._rkx, 800._rkx, 180._rkx, 200._rkx/
  data rcls(2,1:NLUse) &
       /1.e36_rkx,9000._rkx,9000._rkx,9000._rkx,2000._rkx,4000._rkx,1.e36_rkx,1.e36_rkx, &
        9000._rkx,9000._rkx,9000._rkx/
  data rclo(2,1:NLUse) &
       /1.e36_rkx, 400._rkx, 400._rkx, 400._rkx,1000._rkx, 600._rkx,1.e36_rkx,1.e36_rkx, &
       400._rkx, 400._rkx, 400._rkx/

  data ri  (3,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx, 250._rkx, 500._rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx/
  data rlu (3,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,9000._rkx,9000._rkx,4000._rkx,8000._rkx,1.e36_rkx,1.e36_rkx, &
       9000._rkx,9000._rkx,9000._rkx/
  data rac (3,1:NLUse) &
       / 100._rkx,  10._rkx, 100._rkx,1000._rkx,2000._rkx,1500._rkx,   0._rkx, &
         0._rkx, 100._rkx, 50._rkx, 120._rkx/
  data rgss(3,1:NLUse) &
       / 400._rkx, 150._rkx, 350._rkx, 500._rkx, 500._rkx, 200._rkx,   0._rkx, &
        1000._rkx,   0._rkx, 200._rkx, 400._rkx/
  data rgso(3,1:NLUse) &
       / 300._rkx, 150._rkx, 200._rkx, 200._rkx, 200._rkx, 300._rkx,2000._rkx, &
         400._rkx,1000._rkx, 180._rkx, 200._rkx/
  data rcls(3,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,9000._rkx,9000._rkx,3000._rkx,6000._rkx,1.e36_rkx,1.e36_rkx, &
        9000._rkx,9000._rkx,9000._rkx/
  data rclo(3,1:NLUse) &
       /1.e36_rkx,1000._rkx, 400._rkx, 400._rkx,1000._rkx, 600._rkx,1.e36_rkx,1.e36_rkx, &
        800._rkx, 600._rkx, 600._rkx/

  data ri  (4,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx, 400._rkx, 800._rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx/
  data rlu (4,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,1.e36_rkx,1.e36_rkx,6000._rkx,9000._rkx,1.e36_rkx,1.e36_rkx, &
        9000._rkx,9000._rkx,9000._rkx/
  data rac (4,1:NLUse) &
       / 100._rkx,  10._rkx,  10._rkx,1000._rkx,2000._rkx,1500._rkx,   0._rkx, &
          0._rkx,  50._rkx,  10._rkx,  50._rkx/
  data rgss(4,1:NLUse) &
       / 100._rkx, 100._rkx, 100._rkx, 100._rkx, 100._rkx, 100._rkx,   0._rkx, &
        1000._rkx, 100._rkx, 100._rkx,  50._rkx/
  data rgso(4,1:NLUse) &
       / 600._rkx,3500._rkx,3500._rkx,3500._rkx,3500._rkx,3500._rkx,2000._rkx, &
         400._rkx,3500._rkx,3500._rkx,3500._rkx/
  data rcls(4,1:NLUse) &
       /1.e36_rkx,1.e36_rkx,1.e36_rkx,9000._rkx, 200._rkx, 400._rkx,1.e36_rkx,1.e36_rkx, &
        9000._rkx,1.e36_rkx,9000._rkx/
  data rclo(4,1:NLUse) &
       /1.e36_rkx,1000._rkx,1000._rkx, 400._rkx,1500._rkx, 600._rkx,1.e36_rkx,1.e36_rkx, &
        800._rkx,1000._rkx, 800._rkx/

  data ri  (5,1:NLUse) &
       /1.e36_rkx, 120._rkx, 240._rkx, 140._rkx, 250._rkx, 190._rkx,1.e36_rkx,1.e36_rkx, &
       160._rkx, 200._rkx, 300._rkx/
  data rlu (5,1:NLUse) &
       /1.e36_rkx,4000._rkx,4000._rkx,4000._rkx,2000._rkx,3000._rkx,1.e36_rkx,1.e36_rkx, &
        4000._rkx,4000._rkx,8000._rkx/
  data rac (5,1:NLUse) &
       / 100._rkx,  50._rkx,  80._rkx,1200._rkx,2000._rkx,1500._rkx,   0._rkx,  &
         0._rkx, 200._rkx, 60._rkx, 120._rkx/
  data rgss(5,1:NLUse) &
       / 500._rkx, 150._rkx, 350._rkx, 500._rkx, 500._rkx, 200._rkx,   0._rkx, &
        1000._rkx,   0._rkx, 250._rkx, 400._rkx/
  data rgso(5,1:NLUse) &
       / 300._rkx, 150._rkx, 200._rkx, 200._rkx, 200._rkx, 300._rkx,2000._rkx, &
         400._rkx,1000._rkx, 180._rkx, 200._rkx/
  data rcls(5,1:NLUse) &
       /1.e36_rkx,4000._rkx,4000._rkx,4000._rkx,2000._rkx,3000._rkx,1.e36_rkx,1.e36_rkx, &
         4000._rkx,4000._rkx,8000._rkx/
  data rclo(5,1:NLUse) &
       /1.e36_rkx,1000._rkx, 500._rkx, 500._rkx,1500._rkx, 700._rkx,1.e36_rkx,1.e36_rkx, &
         600._rkx, 800._rkx, 800._rkx/

  !---------------------------------------------------------------------------
  !         ... roughness length
  !---------------------------------------------------------------------------
  real(rkx), public, dimension(NSeas,NLUse) :: z0

  data z0  (1,1:NLUse) &
       /1.000_rkx,0.250_rkx,0.050_rkx,1.000_rkx,1.000_rkx,1.000_rkx, &
        0.0006_rkx,0.002_rkx,0.150_rkx,0.100_rkx,0.100_rkx/
  data z0  (2,1:NLUse) &
       /1.000_rkx,0.100_rkx,0.050_rkx,1.000_rkx,1.000_rkx,1.000_rkx, &
        0.0006_rkx,0.002_rkx,0.100_rkx,0.080_rkx,0.080_rkx/
  data z0  (3,1:NLUse) &
       /1.000_rkx,0.005_rkx,0.050_rkx,1.000_rkx,1.000_rkx,1.000_rkx, &
        0.0006_rkx,0.002_rkx,0.100_rkx,0.020_rkx,0.060_rkx/
  data z0  (4,1:NLUse) &
       /1.000_rkx,0.001_rkx,0.001_rkx,1.000_rkx,1.000_rkx,1.000_rkx, &
        0.0006_rkx,0.002_rkx,0.001_rkx,0.001_rkx,0.040_rkx/
  data z0  (5,1:NLUse) &
       /1.000_rkx,0.030_rkx,0.020_rkx,1.000_rkx,1.000_rkx,1.000_rkx, &
        0.0006_rkx,0.002_rkx,0.010_rkx,0.030_rkx,0.060_rkx/

  !real(rkx), private, dimension(11,5), parameter :: z0xxx = reshape ( &
  ! (/   1.000,0.250,0.050,1.000,1.000,1.000,0.0006,0.002,0.150,0.100,0.100 ,  &
  !      1.000,0.100,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.080,0.080 ,  &
  !      1.000,0.005,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.020,0.060 ,  &
  !      1.000,0.001,0.001,1.000,1.000,1.000,0.0006,0.002,0.001,0.001,0.040 ,  &
  !      1.000,0.030,0.020,1.000,1.000,1.000,0.0006,0.002,0.010,0.030,0.060  /), (/11,5/) )

  !---------------------------------------------------------------------------
  ! public chemical data
  !---------------------------------------------------------------------------

  !--- data for foxd (reactivity factor for oxidation) ----
  real(rkx), public, parameter :: dfoxd(n_species_table) = &
          (/  1._rkx     &
             ,1._rkx     &
             ,1._rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1._rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,0._rkx     &
             ,0._rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1._rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1._rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,.1_rkx     &
             ,.1_rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,.1_rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,1.e-36_rkx &
             ,.1_rkx     &
             ,.1_rkx     &
             ,.1_rkx     &
             ,1.e-36_rkx &
             ,1.e-36_rkx & ! HCN
             ,1.e-36_rkx & ! CH3CN
            /)

  Interface seq_drydep_setHCoeff                      ! overload subroutine
     Module Procedure set_hcoeff_scalar
     Module Procedure set_hcoeff_vector
  End Interface

  !--- smallest value to use ---
  real(rkx), private, parameter :: small_value = 1.e-36_rkx

  !---------------------------------------------------------------------------
  ! private chemical data
  !---------------------------------------------------------------------------

  !--- Names of species that can work with ---
  character(len=20) , public , &
          parameter :: species_name_table(n_species_table) = &
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
  real(rkx), public, parameter :: dheff(n_species_table*6) = &
            (/1.15e-2_rkx, 2560._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,8.33e4_rkx, 7379._rkx,2.2e-12_rkx,-3730._rkx,0._rkx     ,    0._rkx  &
             ,3.00e1_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,2.00e3_rkx, 6600._rkx,3.5e-5_rkx,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.00e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.11e2_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,6.30e3_rkx, 6425._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,5.53e3_rkx, 5700._rkx,1.8e-4_rkx,-1510._rkx,0._rkx     ,    0._rkx  &
             ,1.90e-3_rkx, 1480._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,6.40e-3_rkx, 2500._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,0._rkx      ,    0._rkx,2.6e6_rkx, 8700._rkx,0._rkx     ,    0._rkx  &
             ,3.40e-2_rkx, 2420._rkx,4.5e-7_rkx,-1000._rkx,3.6e-11_rkx,-1760._rkx  &
             ,7.40e1_rkx, 3400._rkx,1.7e-5_rkx, -450._rkx,1.0e-14_rkx,-6716._rkx  &
             ,2.14e0_rkx, 3362._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,0.65e0_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,2.20e2_rkx, 4934._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,0._rkx      ,    0._rkx,3.2e1_rkx,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.00e-16_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.14e1_rkx, 6267._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.36e2_rkx, 5995._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,2.20e2_rkx, 5653._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,5.00e0_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,8.37e2_rkx, 5308._rkx,1.8e-4_rkx,-1510._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.00e5_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.71e3_rkx, 7541._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,4.14e4_rkx, 4630._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.45e-3_rkx, 2700._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.00e6_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,2.70e1_rkx, 5300._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.36e2_rkx, 5995._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.47e0_rkx, 5241._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,3.36e2_rkx, 5995._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,0.00e0_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.70e-3_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,2.00e2_rkx, 6500._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.51e3_rkx, 6485._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.00e3_rkx, 6000._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.00e1_rkx,    0._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,7.00e1_rkx, 6000._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,1.20e1_rkx, 5000._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
             ,5.00e1_rkx, 4000._rkx,0._rkx     ,    0._rkx,0._rkx     ,    0._rkx  &
            /)

  real(rkx), private, parameter :: wh2o = amw
  real(rkx), private, parameter :: mol_wgts(n_species_table) = &
     (/ 47.9981995_rkx, 34.0135994_rkx, 17.0067997_rkx, 33.0061989_rkx, 28.0104008_rkx, &
        16.0405998_rkx, 47.0320015_rkx, 48.0393982_rkx, 30.0251999_rkx, 46.0246010_rkx, &
        30.0061398_rkx, 46.0055389_rkx, 63.0123405_rkx, 44.0098000_rkx, 17.0289402_rkx, &
        108.010483_rkx, 62.0049400_rkx, 32.0400009_rkx, 79.0117416_rkx, 15.9994001_rkx, &
        30.0664005_rkx, 61.0578003_rkx, 91.0830002_rkx, 119.093399_rkx, 117.119797_rkx, &
        58.1180000_rkx, 44.0509987_rkx, 62.0652008_rkx, 42.0774002_rkx, 92.0904007_rkx, &
        28.0515995_rkx, 121.047943_rkx, 76.0497971_rkx, 136.228394_rkx, 58.0355988_rkx, &
        72.0614014_rkx, 60.0503998_rkx, 75.0423965_rkx, 44.0922012_rkx, 75.0836029_rkx, &
        58.0768013_rkx, 76.0910034_rkx, 31.9988003_rkx, 33.0061989_rkx, 222.000000_rkx, &
        68.1141968_rkx, 70.0877991_rkx, 70.0877991_rkx, 46.0657997_rkx, 147.125946_rkx, &
        119.074341_rkx, 162.117935_rkx, 100.112999_rkx, 27.0256_rkx   , 41.0524_rkx  /)

  contains
  !
  ! reads drydep_inparm namelist and sets up CCSM driver list of fields for
  ! land-atmosphere communications.
  !
  subroutine seq_drydep_read(NLFilename, seq_drydep_fields)
    implicit none

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    character(len=*), intent(out) :: seq_drydep_fields

    integer(ik4) :: i                ! Indices
    integer(ik4) :: unitn            ! namelist unit number
    integer(ik4) :: ierr             ! error code
    logical :: exists           ! if file exists or not
    character(len=8) :: token   ! dry dep field name to add

    character(*),parameter :: subName = '(seq_drydep_read) '
    character(*),parameter :: F00   = "('(seq_drydep_read) ',8a)"
    character(*),parameter :: FI1   = "('(seq_drydep_init) ',a,I2)"

    namelist /drydep_inparm/ drydep_list, drydep_method

    !-------------------------------------------------------------------------
    ! Read namelist and figure out the drydep field list to pass
    ! First check if file exists and if not, n_drydep will be zero
    !-------------------------------------------------------------------------

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
        end if
      end do
      close( unitn )
      call file_freeUnit( unitn )
    end if

    n_drydep = 0

    !--- Loop over species to fill list of fields to communicate for drydep ---
    seq_drydep_fields = ' '
    do i = 1 , maxspc
      if ( len_trim(drydep_list(i))==0 ) exit
      write(token,333) i
      seq_drydep_fields = trim(seq_drydep_fields)//':'//trim(token)
      if ( i == 1 ) then
        seq_drydep_fields = trim(token)
        drydep_fields_token = trim(token)
      end if
      n_drydep = n_drydep+1
    end do

    !--- Make sure method is valid and determine if land is
    !--- passing drydep fields ---
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
    end if

    ! Need to explicitly add Sl_ based on naming convention
333 format ('Sl_dd',i3.3)

  end subroutine seq_drydep_read
  !
  ! Initialization of dry deposition fields
  ! reads drydep_inparm namelist and sets up CCSM driver list of fields for
  ! land-atmosphere communications.
  !
  subroutine seq_drydep_init( )
    implicit none

    integer(ik4) :: i, l                      ! Indices
    character(len=32) :: test_name       ! field test name

    character(*),parameter :: subName = '(seq_drydep_init) '
    character(*),parameter :: F00   = "('(seq_drydep_init) ',8a)"

    !-------------------------------------------------------------------------
    ! Allocate and fill foxd, drat and mapping as well as species indices
    !-------------------------------------------------------------------------

    if ( n_drydep > 0 ) then

      allocate( foxd(n_drydep) )
      allocate( drat(n_drydep) )
      allocate( mapping(n_drydep) )

      ! This initializes these variables to infinity.
      foxd = inf
      drat = inf

      mapping(:) = 0

    end if

    h2_ndx   = -1
    ch4_ndx  = -1
    co_ndx   = -1
    mpan_ndx = -1
    pan_ndx  = -1
    so2_ndx  = -1
    o3_ndx   = -1
    xpan_ndx = -1

    !--- Loop over drydep species that need to be worked with ---
    do i = 1 , n_drydep
      if ( len_trim(drydep_list(i))==0 ) exit

      test_name = drydep_list(i)

      if ( trim(test_name) == 'O3' ) then
        test_name = 'OX'
      end if

      !--- Figure out if species maps to a species in the species table ---
      do l = 1 , n_species_table
        if ( trim( test_name ) == trim( species_name_table(l) ) ) then
          mapping(i)  = l
          exit
        end if
      end do

      !--- If it doesn't map to a species in the species table
      !--- find species close enough ---
      if ( mapping(i) < 1 ) then
        select case( trim(test_name) )
          case( 'H2' )
            test_name = 'CO'
          case( 'HYAC', 'CH3COOH', 'EOOH' )
            test_name = 'CH2O'
          case( 'O3S', 'O3INERT', 'MPAN' )
            test_name = 'OX'
          case( 'ISOPOOH', 'MACROOH', 'Pb', 'XOOH', 'H2SO4' )
            test_name = 'HNO3'
          case( 'ALKOOH', 'MEKOOH', 'TOLOOH', 'BENOOH', &
                'XYLOOH', 'TERPOOH','SOGM','SOGI','SOGT','SOGB','SOGX' )
            test_name = 'CH3OOH'
          case( 'SOA', 'SO2', 'SO4', 'CB1', 'CB2', 'OC1', &
                'OC2', 'NH3', 'NH4', 'SA1', 'SA2', 'SA3', &
                'SA4','HCN','CH3CN','HCOOH' )
            ! this is just a place holder. values are explicitly set below
            test_name = 'OX'
          case( 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
            ! this is just a place holder. values are explicitly set below
            test_name = 'OX'
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
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'CO01','CO02','CO03','CO04','CO05','CO06', &
                'CO07','CO08','CO09','CO10' )
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'CO11','CO12','CO13','CO14','CO15','CO16', &
                'CO17','CO18','CO19','CO20' )
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'CO21','CO22','CO23','CO24','CO25','CO26', &
                'CO27','CO28','CO29','CO30' )
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'CO31','CO32','CO33','CO34','CO35','CO36', &
                'CO37','CO38','CO39','CO40' )
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'CO41','CO42','CO43','CO44','CO45','CO46', &
                'CO47','CO48','CO49','CO50' )
            ! this is just a place holder. values are set in drydep_fromlnd
            test_name = 'CO'
          case( 'NH4NO3' )
            test_name = 'HNO3'
          case default
            test_name = 'blank'
        end select

        !--- If found a match check the species table again ---
        if ( trim(test_name) /= 'blank' ) then
          do l = 1 , n_species_table
            if ( trim( test_name ) == trim( species_name_table(l) ) ) then
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
      if ( trim(drydep_list(i)) == 'OX' .or. &
           trim(drydep_list(i)) == 'O3' ) o3_ndx  = i
      if ( trim(drydep_list(i)) == 'O3A' ) o3a_ndx  = i
      if ( trim(drydep_list(i)) == 'XPAN' ) xpan_ndx = i

      if ( mapping(i) > 0 ) then
        l = mapping(i)
        foxd(i) = dfoxd(l)
        drat(i) = sqrt(mol_wgts(l)/wh2o)
      end if
    end do

    where( rgss < 1._rkx )
      rgss = 1._rkx
    endwhere

    where( rac < small_value)
      rac = small_value
    endwhere
  end subroutine seq_drydep_init
  !
  ! Interface to seq_drydep_setHCoeff when input is scalar
  ! wrapper routine used when surface temperature is a scalar (single column)
  ! rather than an array (multiple columns).
  !
  subroutine set_hcoeff_scalar( sfc_temp, heff )
    implicit none
    real(rkx), intent(in) :: sfc_temp        ! Input surface temperature
    real(rkx), intent(out) :: heff(n_drydep) ! Output Henry's law coefficients

    real(rkx) :: sfc_temp_tmp(1)    ! surface temp

    sfc_temp_tmp(:) = sfc_temp
    call set_hcoeff_vector( 1, sfc_temp_tmp, heff(:n_drydep) )

  end subroutine set_hcoeff_scalar
  !
  ! Interface to seq_drydep_setHCoeff when input is vector
  ! sets dry depositions coefficients -- used by both land and atmosphere models
  !
  subroutine set_hcoeff_vector( ncol, sfc_temp, heff )
    implicit none
    integer(ik4), intent(in) :: ncol ! Input size of surface-temp vector
    real(rkx), intent(in) :: sfc_temp(ncol)        ! Surface temperature
    real(rkx), intent(out) :: heff(ncol,n_drydep)  ! Henry's law coefficients

    real(rkx), parameter :: t0     = 298._rkx    ! Standard Temperature
    real(rkx), parameter :: ph_inv = 1._rkx/ph   ! Inverse of PH
    integer(ik4)  :: m, l, id   ! indices
    real(rkx) :: e298  ! Henry's law coefficient @ standard temperature (298K)
    real(rkx) :: dhr   ! temperature dependence of Henry's law coefficient
    real(rkx) :: dk1s(ncol)  ! DK Work array 1
    real(rkx) :: dk2s(ncol)  ! DK Work array 2
    real(rkx) :: wrk(ncol)   ! Work array

    character(*),parameter :: subName = '(seq_drydep_set_hcoeff) '
    character(*),parameter :: F00   = "('(seq_drydep_set_hcoeff) ',8a)"

    wrk(:) = (t0 - sfc_temp(:))/(t0*sfc_temp(:))
    do m = 1 , n_drydep
      l    = mapping(m)
      id   = 6*(l - 1)
      e298 = dheff(id+1)
      dhr  = dheff(id+2)
      heff(:,m) = e298*exp( dhr*wrk(:) )
      !--- Calculate coefficients based on the drydep tables ---
      if ( dheff(id+3) /= 0._rkx .and. dheff(id+5) == 0._rkx ) then
        e298 = dheff(id+3)
        dhr  = dheff(id+4)
        dk1s(:) = e298*exp( dhr*wrk(:) )
        where( heff(:,m) /= 0._rkx )
          heff(:,m) = heff(:,m)*(1._rkx + dk1s(:)*ph_inv)
        elsewhere
          heff(:,m) = dk1s(:)*ph_inv
        end where
      end if
      !--- For coefficients that are non-zero AND CO2 or NH3 handle
      !--- things this way ---
      if ( dheff(id+5) /= 0._rkx ) then
        if ( trim( drydep_list(m) ) == 'CO2' .or. &
             trim( drydep_list(m) ) == 'NH3' ) then
          e298 = dheff(id+3)
          dhr  = dheff(id+4)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          e298 = dheff(id+5)
          dhr  = dheff(id+6)
          dk2s(:) = e298*exp( dhr*wrk(:) )
          !--- For Carbon dioxide ---
          if ( trim(drydep_list(m)) == 'CO2' ) then
            heff(:,m) = heff(:,m)*(1._rkx + dk1s(:)*ph_inv) * &
                    (1._rkx + dk2s(:)*ph_inv)
            !--- For NH3 ---
          else if ( trim( drydep_list(m) ) == 'NH3' ) then
            heff(:,m) = heff(:,m)*(1._rkx + dk1s(:)*ph/dk2s(:))
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

end module mod_clm_drydep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
