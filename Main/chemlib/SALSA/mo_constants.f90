MODULE mo_constants

!- Description:
!
!  This module contains basic constants and derived constants
!
!- Author:
!
!  M. Giorgetta, MPI, April 1999
!  I. Kirchner, MPI, December 2000, time control
!  L. Kornblueh, MPI, January 2001, cleanup

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

! Universal constants
  REAL(dp), PARAMETER :: api   = 3.14159265358979323846_dp ! pi
  REAL(dp), PARAMETER :: asqrt2= 1.41421356237309504880_dp
! REAL(dp), PARAMETER :: ar    = 8.314e3_dp       ! universal gas constant in J/K/kmol
  REAL(dp), PARAMETER :: argas = 8.314409_dp      ! universal gas constant (R) in J/K/mol
  REAL(dp), PARAMETER :: avo   = 6.02214e23_dp    ! Avogadro constant in 1/mol
  REAL(dp), PARAMETER :: ak    = 1.380662e-23_dp  ! Boltzmann constant in J/K
  REAL(dp), PARAMETER :: stbo  = 5.67E-8_dp       ! Stephan-Boltzmann constant in W/m2/K4
  REAL(dp), PARAMETER :: ap0   = 101325._dp       ! Standard pressure in Pa

! Molar weights in g/mol
  REAL(dp), PARAMETER :: amco2 = 44.011_dp        ! molecular weight of carbon dioxide
  REAL(dp), PARAMETER :: amch4 = 16.043_dp        ! molecular weight of methane
  REAL(dp), PARAMETER :: amo3  = 47.9982_dp       ! molecular weight of ozone
  REAL(dp), PARAMETER :: amn2o = 44.013_dp        ! molecular weight of N2O
  REAL(dp), PARAMETER :: amc11 =137.3686_dp       ! molecular weight of CFC11
  REAL(dp), PARAMETER :: amc12 =120.9140_dp       ! molecular weight of CFC12
  REAL(dp), PARAMETER :: amo2  = 31.9988_dp       ! molecular weight of molecular oxygen
  REAL(dp), PARAMETER :: amw   = 18.0154_dp       ! molecular weight of water vapor
  REAL(dp), PARAMETER :: amd   = 28.970_dp        ! molecular weight of dry air

! Dry air and water vapour thermodynamic constants
  REAL(dp), PARAMETER :: cpd   = 1005.46_dp       ! specific heat of dry air at constant
                                                  ! pressure in J/K/kg
  REAL(dp), PARAMETER :: cpv   = 1869.46_dp       ! specific heat of water vapour at
                                                  ! constant pressure in J/K/kg
  REAL(dp), PARAMETER :: rd    = 287.05_dp        ! gas constant for dry air in J/K/kg
  REAL(dp), PARAMETER :: rv    = 461.51_dp        ! gas constant for water vapour
                                                  ! in J/K/kg
  REAL(dp), PARAMETER :: rcpd  = 1.0_dp/cpd       ! auxiliary constant in K*kg/J
  REAL(dp)            :: vtmpc1= rv/rd-1.0_dp     ! dimensionless auxiliary constant
  REAL(dp)            :: vtmpc2= cpv/cpd-1.0_dp   ! dimensionless auxiliary constant

! H2O related constants  (liquid, ice, snow), phase change constants
  REAL(dp), PARAMETER :: rhoh2o = 1000.0_dp        ! density of liquid water in kg/m3
  REAL(dp), PARAMETER :: rhosea = 1025.0_dp        ! density of sea water in kg/m3
  REAL(dp), PARAMETER :: rhoice = 917.0_dp         ! density of ice in kg/m3
  REAL(dp), PARAMETER :: rhosno = 300.0_dp         ! density of snow in kg/m3
  REAL(dp), PARAMETER :: rhoiw  = rhoice/rhoh2o    ! density ratio (ice/water)
  REAL(dp), PARAMETER :: alv    = 2.5008e6_dp      ! latent heat for vaporisation in J/kg
  REAL(dp), PARAMETER :: als    = 2.8345e6_dp      ! latent heat for sublimation in J/kg
  REAL(dp), PARAMETER :: alf    = als-alv          ! latent heat for fusion in J/kg
  REAL(dp), PARAMETER :: cpliq  = 4218._dp         ! specific heat for liquid water J/K/kg
  REAL(dp), PARAMETER :: cpsea  = 3994._dp         ! specific heat for sea water J/K/kg
  REAL(dp), PARAMETER :: cpice  = 2106._dp         ! specific heat for ice J/K/kg
  REAL(dp), PARAMETER :: cpsno  = 2090._dp         ! specific heat for snow J/K/kg
  REAL(dp), PARAMETER :: alice  = 2.1656_dp        ! thermal conductivity of ice in W/K/m
  REAL(dp), PARAMETER :: alsno  = 0.31_dp          ! thermal conductivity of snow in W/K/m
  REAL(dp), PARAMETER :: tmelt  = 273.15_dp        ! melting temperature of ice/snow in K

! Earth and earth orbit parameters
  REAL(dp), PARAMETER :: a     = 6371000.0_dp     ! radius of the earth in m
  REAL(dp), PARAMETER :: rae   = 0.1277E-2_dp     ! ratio of atmosphere to earth radius
  REAL(dp), PARAMETER :: omega = .7292E-4_dp      ! solid rotation velocity of the earth
                                                  ! in 1/s
  REAL(dp), PARAMETER :: secperday = 86400._dp    ! seconds per day
  REAL(dp), PARAMETER :: g     = 9.80665_dp       ! gravity acceleration in m/s2
  REAL(dp), PARAMETER :: qg    = 1.0_dp/g         ! inverse of gravity acceleration
 
! Constants used for computation of saturation mixing ratio
! over liquid water (*c_les*) or ice(*c_ies*)
  REAL(dp), PARAMETER :: c1es  = 610.78_dp           !
  REAL(dp), PARAMETER :: c2es  = c1es*rd/rv          !
  REAL(dp), PARAMETER :: c3les = 17.269_dp           !
  REAL(dp), PARAMETER :: c3ies = 21.875_dp           !
  REAL(dp), PARAMETER :: c4les = 35.86_dp            !
  REAL(dp), PARAMETER :: c4ies = 7.66_dp             !
  REAL(dp), PARAMETER :: c5les = c3les*(tmelt-c4les) !
  REAL(dp), PARAMETER :: c5ies = c3ies*(tmelt-c4ies) !
  REAL(dp), PARAMETER :: c5alvcp = c5les*alv/cpd     !
  REAL(dp), PARAMETER :: c5alscp = c5ies*als/cpd     !
  REAL(dp), PARAMETER :: alvdcp  = alv/cpd           !
  REAL(dp), PARAMETER :: alsdcp  = als/cpd           !

! Specifications, thresholds, and derived constants for the following subroutines:
! s_lake, s_licetemp, s_sicetemp, meltpond, meltpond_ice, update_albedo_ice_meltpond

  REAL(dp), PARAMETER :: dmix     = 10.0_dp   ! mixed-layer depth of lakes in m
  REAL(dp), PARAMETER :: dmixsea  = 50.0_dp   ! mixed-layer depth of ocean in m
  REAL(dp), PARAMETER :: dice     = 0.05_dp   ! minimum ice thickness in m
  REAL(dp), PARAMETER :: dicepond = 0.01_dp   ! minimum ice thickness of pond ice in m
  REAL(dp), PARAMETER :: dicelim  = 0.10_dp   ! threshold ice thickness for pond closing in m
  REAL(dp), PARAMETER :: dpondmin = 0.01_dp   ! minimum pond depth for pond fraction in m
  REAL(dp), PARAMETER :: albpondi = 0.30_dp   ! albedo of pond ice

  REAL(dp), PARAMETER :: snicecond = alice/alsno * rhoh2o/rhosno
  REAL(dp), PARAMETER :: hcapmix   = rhoh2o*cpliq*dmix     ! heat capacity of lake mixed layer in J/K/m2
  REAL(dp), PARAMETER :: hcapice   = rhoice*cpice*dice     ! heat capacity of upper ice layer
  REAL(dp), PARAMETER :: hcapicep  = rhoice*cpice*dicepond ! heat capacity of upper pond ice layer
  REAL(dp), PARAMETER :: rhoilf    = rhoice*alf            ! [J/m3]
  REAL(dp), PARAMETER :: rhowlf    = rhoh2o*alf            ! [J/m3]
  REAL(dp), PARAMETER :: hcaprilf  = hcapmix/rhoilf        ! [m/K]
  REAL(dp), PARAMETER :: rilfhcap  = rhoilf/hcapmix        ! [K/m]
  REAL(dp), PARAMETER :: tfreez    = dice*rilfhcap         ! cooling below tmelt required to form dice
 
END MODULE mo_constants
