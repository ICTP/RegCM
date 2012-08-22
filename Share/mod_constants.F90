!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_constants

  use mod_realkinds

  implicit none

  ! numbers
  real(dp) , parameter :: d_zero = 0.0D+00
  real(dp) , parameter :: d_one = 1.0D+00
  real(dp) , parameter :: d_two = 2.0D+00
  real(dp) , parameter :: d_three = 3.0D+00
  real(dp) , parameter :: d_four = 4.0D+00
  real(dp) , parameter :: d_five = 5.0D+00
  real(dp) , parameter :: d_six = 6.0D+00
  real(dp) , parameter :: d_nine = 6.0D+00
  real(dp) , parameter :: d_half = 0.50D+00
  real(dp) , parameter :: d_rfour = 0.250D+00
  real(dp) , parameter :: d_twelve = 12.0D+00
  real(dp) , parameter :: d_60 = 60.0D+00
  real(dp) , parameter :: d_10 = 1.0D+01
  real(dp) , parameter :: d_r10 = 1.0D-01
  real(dp) , parameter :: d_100 = 1.0D+02
  real(dp) , parameter :: d_r100 = 1.0D-02
  real(dp) , parameter :: d_1000 = 1.0D+03
  real(dp) , parameter :: d_r1000 = 1.0D-03
  real(dp) , parameter :: onet = d_one/d_three
  real(dp) , parameter :: twot = d_two/d_three
  real(dp) , parameter :: fourt = d_four/d_three

  ! Angles degrees
  real(dp) , parameter :: deg00  = 0.0D+00
  real(dp) , parameter :: deg45  = 45.0D+00
  real(dp) , parameter :: deg90  = 90.0D+00
  real(dp) , parameter :: deg180 = 180.0D+00
  real(dp) , parameter :: deg360 = 360.0D+00

  real(dp) , parameter :: minqx   = 1.0D-08
  real(dp) , parameter :: dlowval = 1.0D-30
  real(sp) , parameter :: slowval = 1.0E-30
  real(dp) , parameter :: dmissval = 1.0D+20
  real(sp) , parameter :: smissval = 1.0E+20

  ! time conversion
  real(dp) , parameter :: secpm = 60.0D+00
  real(dp) , parameter :: secph = 3600.0D+00
  real(dp) , parameter :: secpd = 86400.0D+00
  real(dp) , parameter :: minph = 60.0D+00
  real(dp) , parameter :: minpd = 1440.0D+00
  real(dp) , parameter :: houpd = 24.0D+00

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(dp) , parameter :: egrav = 9.80665D+00

  ! Gas constant for dry air in Joules/kg/K
  real(dp) , parameter :: rgas = 287.0058D+00
  ! Gas constant for water vapor in Joules/kg/K
  real(dp) , parameter :: rwat = 461.90D+00
  ! Ratio of the two above
  real(dp) , parameter :: rgow = rgas/rwat
  ! Reverse of the above
  real(dp) , parameter :: rgowi = rwat/rgas
  ! Helper value to ease calculations
  real(dp) , parameter :: retv = rwat/rgas - 1.0D0

  ! Specific heat at constant pressure for dry air J/kg/K
  real(dp) , parameter :: cpd = 1005.46D+00
  ! Specific heat at constant pressure for moist air J/kg/K
  real(dp) , parameter :: cpv = 1869.46D+00
  ! Specific heat of water at 15 Celsius J/kg/K
  real(dp) , parameter :: cpw = 4186.95D+00
  ! Specific heat of water at 0 Celsius J/kg/K
  real(dp) , parameter :: cpw0 = 4218.0D+00

  ! Specific heats per m**3  (joules/m**3/k)
  real(dp) , parameter :: ch2o = 4.186D+06
  real(dp) , parameter :: cice = 0.45D+00*ch2o
  real(dp) , parameter :: cwi = d_one/0.45D+00
  real(dp) , parameter :: csnw = 0.49D+00*ch2o
  real(dp) , parameter :: cws = d_one/0.49D+00

  ! Latent heats (Joules/kg)
  real(dp) , parameter :: wlhf = 0.3336D+06
  real(dp) , parameter :: wlhv = 2.51040D+06
  real(dp) , parameter :: wlhs = wlhv + wlhf

  ! Various utility terms used in calculations
  real(dp) , parameter :: regrav = d_one/egrav
  real(dp) , parameter :: rcpd = d_one/cpd
  real(dp) , parameter :: rovcp = rgas*rcpd
  real(dp) , parameter :: rovg  = rgas/egrav
  real(dp) , parameter :: govr  = egrav/rgas
  real(dp) , parameter :: vtmpc1 = rwat/rgas - d_one
  real(dp) , parameter :: vtmpc2 = cpv*rcpd - d_one
  real(dp) , parameter :: rhoh2o = 1000.0D+00
  real(dp) , parameter :: rhosnow = 330.0D+00
  real(dp) , parameter :: rhoice = 917.0D+00
  real(dp) , parameter :: tzero = 273.15D+00
  real(dp) , parameter :: rtzero = d_one/tzero
  real(dp) , parameter :: wattp = 273.16D+00
  real(dp) , parameter :: tboil = 373.1339D+00
  real(dp) , parameter :: c1es = 610.78D+00
  real(dp) , parameter :: c2es = c1es*rgas/rwat
  real(dp) , parameter :: c3les = 17.2693882D+00
  real(dp) , parameter :: c3ies = 21.875D+00
  real(dp) , parameter :: c4les = 35.86D+00
  real(dp) , parameter :: c4ies = 7.66D+00
  real(dp) , parameter :: c5les = c3les*(tzero-c4les)
  real(dp) , parameter :: c5ies = c3ies*(tzero-c4ies)
  real(dp) , parameter :: c5alvcp = c5les*wlhv*rcpd
  real(dp) , parameter :: c5alscp = c5ies*wlhs*rcpd
  real(dp) , parameter :: wlhvocp = wlhv*rcpd
  real(dp) , parameter :: wlhsocp = wlhs*rcpd
  real(dp) , parameter :: wlhfocp = wlhf*rcpd
  real(dp) , parameter :: cpowlhv = cpd/wlhv
  real(dp) , parameter :: cpowlhs = cpd/wlhs
  real(dp) , parameter :: cpowlhf = cpd/wlhf
  real(dp) , parameter :: rtice = tzero-23.0D+00
  real(dp) , parameter :: rtwat = tzero
  real(dp) , parameter :: mpcrt = rtice+(rtwat-rtice)/sqrt(2.0D+00)
  real(dp) , parameter :: rtwat_rtice_r = 1.0D+00/(rtwat-rtice)
  real(dp) , parameter :: pq0 = 379.90516D+00
  ! value used for the latent heat term in the exponent for
  ! calculating equivalent potential temperature
  real(dp) , parameter :: eliwv = 2.72D+06

  ! Standard atmosphere ICAO 1993
  real(dp) , parameter :: stdp = 1.013250D+05
  real(dp) , parameter :: stdpmb = 1013.250D+00
  real(dp) , parameter :: stdt = 288.15D+00
  real(dp) , parameter :: lrate = 0.00649D+00 ! K/km from MSL up to 11 km
  real(dp) , parameter :: bltop = 0.960D+00
 
  ! Stefan-Boltzmann  constant CODATA 2007
  real(dp) , parameter :: sigm = 5.670400D-08
  ! Boltzman Constant k CODATA 2007
  real(dp) , parameter :: boltzk = 1.3806504D-23
  ! Avogadro Constant
  real(dp) , parameter :: navgdr = 6.02214084D23
  ! Same in lit/cm3
  real(dp) , parameter :: avogadrl = navgdr*d_1000

  ! Atmos. surface pressure mol/cm3
  real(dp) , parameter :: atmos = 2.247D19
  ! Conversion parameter for Henry L-atm/mol-K
  real(dp) , parameter :: rtcon = 8.314D-02
  ! RU g-cm2/s2-mol-K
  real(dp) , parameter :: rumolec = 8.314D7
  ! Droplet diffusion coefficient Lelieveld+Crutzen 1991 cm2/sec
  real(dp) , parameter :: dropdif = 2.0D-05
  ! Gas-phase diffusion coeff. Lelieveld and Crutzen, 1991 cm2/s
  real(dp) , parameter :: difgas = 0.1D0

  ! Fixed emissivity of water
  real(dp) , parameter :: emsw = 0.97D+00

  ! Trigonometric constants. 
  real(dp) , parameter :: mathpi =                                   &
                      &   3.1415926535897932384626433832795029D+00
  real(dp) , parameter :: invpi = d_one/mathpi
  real(dp) , parameter :: halfpi = mathpi*d_half
  real(dp) , parameter :: twopi = mathpi*d_two
  real(dp) , parameter :: degrad = mathpi/180.0D+00
  real(dp) , parameter :: raddeg = 180.0D+00/mathpi

  ! Maximum stomatl resistance (s/m)
  real(dp) , parameter :: rmax0 = 2.0D+04

  ! Maximum allowed dew(mm) and inverse (dewmaxi)
  real(dp) , parameter :: dewmax = 0.1D+00
  real(dp) , parameter :: dewmaxi = d_one/dewmax

  ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
  real(dp) , parameter :: trsmx0 = 2.D-04

  ! Drainage out of 10m layer bottom (mm/s)
  ! drain is set fairly large to prevent swamping the soil
  real(dp) , parameter :: drain = 1.0D-04

  ! Earth radius in meters
  real(dp) , parameter :: earthrad = 6.371229D+06
  real(dp) , parameter :: erkm = earthrad/d_1000
  ! Angular velocity of rotation of Earth
  real(dp) , parameter :: eomeg = 7.2921159D-05
  real(dp) , parameter :: eomeg2 = d_two*eomeg
  ! Solar Constant in W/m**2
  real(dp) , parameter :: solcon = 1367.0D+00
  ! Solar Constant in erg/cm**2/sec
  real(dp) , parameter :: scon = solcon*d_1000

  ! Soil roughness length
  real(dp) , parameter :: zlnd = 0.01D+00
  ! Ocean roughness length
  real(dp) , parameter :: zoce = 0.00040D+00
  ! Snow roughness length
  real(dp) , parameter :: zsno = 0.00040D+00
  ! Von Karman constant
  real(dp) , parameter :: vonkar = 0.4D+00
  ! Drag coefficient for soil under canopy
  real(dp) , parameter :: csoilc = 4.0D-03
  ! Turbulent wind for stable conditions (m/sec)
  real(dp) , parameter :: wtur = 0.1D+00

  ! Constant used in computing virture temperature.
  real(dp) , parameter :: ep1 = 0.608D+00
  ! Constant used in computing saturation mixing ratio.
  ! Ratio of mean molecular weight of water (18.016 g/mole)
  ! to that of dry air (28.966 g/mole)
  real(dp) , parameter :: ep2 = 0.62197D+00
  ! Constants used in computing saturation vapor pressure.
  real(dp) , parameter :: svp1 = 0.6112D+00
  real(dp) , parameter :: svp2 = 17.67D+00
  real(dp) , parameter :: svp3 = 29.65D+00
  real(dp) , parameter :: svp4 = 0.611D+00
  real(dp) , parameter :: svp5 = 22.514D+00
  real(dp) , parameter :: svp6 = 6.15D+03
  ! Constants used in computing evaporation latent heat.
  real(dp) , parameter :: lh0 = 597.3D+00
  real(dp) , parameter :: lh1 = 0.566D+00
  ! Constants from latent heat and temperature to saturation vapor p.
  real(dp) , parameter :: lsvp1 = 6.11D+00
  real(dp) , parameter :: lsvp2 = 9.045D+00

  ! Terminal velocity constants
  real(dp) , parameter :: avt = 841.99667D+00
  real(dp) , parameter :: bvt = 0.8D+00
  real(dp) , parameter :: g4pb = 17.837825D+00
  real(dp) , parameter :: g3pb = g4pb/(3.0D+00+bvt)
  real(dp) , parameter :: g5pb = 1.8273D+00
  real(dp) , parameter :: vtc = avt*g4pb/6.0D+00
  real(dp) , parameter :: trel = 3000.0D+00

  ! Dynamic parameters
  ! alpha = .2495 in brown-campana; = 0. in split explicit
  real(dp) , parameter :: alpha = 0.0D+00
  real(dp) , parameter :: beta = d_one - d_two*alpha
  real(dp) , parameter :: gnu = 0.10D+00
  real(dp) , parameter :: omu = d_one - d_two*gnu
  real(dp) , parameter :: gnuhf = d_half*gnu
  real(dp) , parameter :: omuhf = d_one - d_two*gnuhf

  ! Cumulous parameters
  real(dp) , parameter :: tauht = 7200.0D+00

  ! Aerosol densities
  ! now defined in chemistry modules since they are not constant

  ! Constants used in Betts Miller
  real(dp) , parameter :: aliq = 613.3D+00
  real(dp) , parameter :: bliq = 17.502D+00
  real(dp) , parameter :: cliq = 4780.8D+00
  real(dp) , parameter :: dliq = 32.19D+00
  real(dp) , parameter :: aice = 613.2D+00
  real(dp) , parameter :: bice = 22.452D+00
  real(dp) , parameter :: cice1 = 6133.0D+00
  real(dp) , parameter :: dice = 0.61D+00
  real(dp) , parameter :: xls0 = 2.905D+06
  real(dp) , parameter :: xls1 = 259.532D+00

  ! GTS system constants
  real(dp) , parameter :: egravgts = egrav*d_100
  real(dp) , parameter :: regravgts = d_one/egravgts
  real(dp) , parameter :: cpdgts = cpd*1.0D+04
  real(dp) , parameter :: gocp = egravgts/cpdgts
  real(dp) , parameter :: sslp = stdp*d_10 ! dynes/cm^2
  real(dp) , parameter :: rsslp = d_one/sslp
  real(dp) , parameter :: stebol = sigm*d_1000
  real(dp) , parameter :: rgsslp = d_half/(egravgts*sslp)
  ! Effective molecular weight of dry air (g/mol)
  real(dp) , parameter :: amd = 28.9644D+00
  ! Effective molecular weight of dry air (kg/mol)
  real(dp) , parameter :: amdk = amd*d_r1000
  ! Molecular weight of ozone (g/mol)
  real(dp) , parameter :: amo = 48.0D+00
  ! Molecular weight of co2 (g/mol)
  real(dp) , parameter :: amco2 = 44.0D+00

  ! Radiation constants
  real(dp) , parameter :: dpfco2 = 5.0D-03
  real(dp) , parameter :: dpfo3 = 2.5D-03

  ! Pressure gradient force calculations (Why not standard atmosphere?)
  real(dp) , parameter :: t00pg = 287.0D+00       ! stdt ?
  real(dp) , parameter :: p00pg = 101.325D+00     ! stdp ?
  real(dp) , parameter :: alam  = 6.5D-03         ! Lapse rate ?
  real(dp) , parameter :: pgfaa1 = alam*rgas*regrav ! Utility constant

  ! Molecular heat diffusion coefficient in water
  real(dp) , parameter :: hdmw = 1.3889D-07  ! m^2/s

  ! Seaice temperature from ICBC trigger value
  real(dp) , parameter :: icetemp = 271.38D0

end module mod_constants
