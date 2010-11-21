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

      use m_realkinds

      implicit none

      ! Standard Gravity (m/sec**2) 3rd CGPM
      real(DP) , parameter :: gti = 9.80665D+00

      ! Gas constant for dry air in Joules/kg/K
      real(DP) , parameter :: rgas = 287.0058D+00
      ! Gas constant for water vapor in Joules/kg/K
      real(DP) , parameter :: rwat = 461.90D+00

      ! Specific heat at constant pressure for dry air J/kg/K
      real(DP) , parameter :: cpd = 1005.46D+00
      ! Specific heat at constant pressure for moist air J/kg/K
      real(DP) , parameter :: cpv = 1869.46D+00
      ! Specific heat of water at 15 Celsius J/kg/K
      real(DP) , parameter :: cpw = 4186.95D+00
      ! Specific heat of water at 0 Celsius J/kg/K
      real(DP) , parameter :: cpw0 = 4218.0D+00

      ! Specific heats per m**3  (joules/m**3/k)
      real(DP) , parameter :: ch2o = 4.186D+06
      real(DP) , parameter :: cice = 0.45D+00*ch2o
      real(DP) , parameter :: cwi = 1.0D+00/0.45D+00
      real(DP) , parameter :: csnw = 0.49D+00*ch2o
      real(DP) , parameter :: cws = 1.0D+00/0.49D+00

      ! Latent heats (Joules/kg)
      real(DP) , parameter :: wlhf = 0.3336D+06
      real(DP) , parameter :: wlhv = 2.51040D+06
      real(DP) , parameter :: wlhs = wlhv + wlhf

      ! Various utility terms used in calculations
      real(DP) , parameter :: rgti = 1.0D+00/gti
      real(DP) , parameter :: rcpd = 1.0D+00/cpd
      real(DP) , parameter :: rovcp = rgas*rcpd
      real(DP) , parameter :: rovg  = rgas/gti
      real(DP) , parameter :: govr  = gti/rgas
      real(DP) , parameter :: vtmpc1 = rwat/rgas - 1.0D+00
      real(DP) , parameter :: vtmpc2 = cpv*rcpd - 1.0D+00
      real(DP) , parameter :: rhoh2o = 1000.0D+00
      real(DP) , parameter :: rhos = 330.0D+00
      real(DP) , parameter :: rhoi = 917.0D+00
      real(DP) , parameter :: tzero = 273.15D+00
      real(DP) , parameter :: rtzero = 1.0D+00/tzero
      real(DP) , parameter :: wattp = 273.16D+00
      real(DP) , parameter :: tboil = 373.1339D+00
      real(DP) , parameter :: c1es = 610.78D+00
      real(DP) , parameter :: c2es = c1es*rgas/rwat
      real(DP) , parameter :: c3les = 17.2693882D+00
      real(DP) , parameter :: c3ies = 21.875D+00
      real(DP) , parameter :: c4les = 35.86D+00
      real(DP) , parameter :: c4ies = 7.66D+00
      real(DP) , parameter :: c5les = c3les*(tzero-c4les)
      real(DP) , parameter :: c5ies = c3ies*(tzero-c4ies)
      real(DP) , parameter :: c5alvcp = c5les*wlhv*rcpd
      real(DP) , parameter :: c5alscp = c5ies*wlhs*rcpd
      real(DP) , parameter :: wlhvocp = wlhv*rcpd
      real(DP) , parameter :: wlhsocp = wlhs*rcpd
      real(DP) , parameter :: pq0 = 379.90516D+00
      ! value used for the latent heat term in the exponent for
      ! calculating equivalent potential temperature
      real(DP) , parameter :: eliwv = 2.72D+06

      ! Standard atmosphere ICAO 1993
      real(DP) , parameter :: stdp = 1.013250D+05
      real(DP) , parameter :: stdpmb = 1013.250D+00
      real(DP) , parameter :: stdt = 288.15D+00
      real(DP) , parameter :: lrate = 0.00649D+00 ! K/km from MSL up to 11 km
      real(DP) , parameter :: bltop = 0.960D+00
 
      ! Stefan-Boltzmann  constant CODATA 2007
      real(DP) , parameter :: sigm = 5.670400D-08
      ! Boltzman Constant k CODATA 2007
      real(DP) , parameter :: boltzk = 1.3806504D-23
      ! Avogadro Constant
      real(DP) , parameter :: navgdr = 6.02214084D23

      ! Fixed emissivity of water
      real(DP) , parameter :: emsw = 0.97D+00

      ! Trigonometric constants. 
      real(DP) , parameter :: mathpi =                                   &
                          &   3.1415926535897932384626433832795029D+00
      real(DP) , parameter :: invpi = 1.0D+00/mathpi
      real(DP) , parameter :: halfpi = mathpi*0.5D+00
      real(DP) , parameter :: twopi = mathpi*2.0D+00
      real(DP) , parameter :: degrad = mathpi/180.0D+00
      real(DP) , parameter :: raddeg = 180.0D+00/mathpi

      ! Maximum stomatl resistance (s/m)
      real(DP) , parameter :: rmax0 = 2.0D+04

      ! Maximum allowed dew(mm) and inverse (dewmxi)
      real(DP) , parameter :: dewmx = 0.1D+00
      real(DP) , parameter :: dewmxi = 1.0D+00/dewmx

      ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
      real(DP) , parameter :: trsmx0 = 2.D-04

      ! Drainage out of 10m layer bottom (mm/s)
      ! drain is set fairly large to prevent swamping the soil
      real(DP) , parameter :: drain = 1.D-04

      ! Earth radius in meters
      real(DP) , parameter :: earthrad = 6.371229D+06
      real(DP) , parameter :: erkm = earthrad/1000.0D+00
      ! Length of day in seconds
      real(DP) , parameter :: tau1 = 8.64D+04
      ! Days per year
      real(DP) , parameter :: dayspy = 365.2422D+00
      ! Degrees per day
      real(DP) , parameter :: dpd = 360.0/dayspy
      ! Angular velocity of rotation of Earth
      real(DP) , parameter :: eomeg = 7.2921159D-05
      real(DP) , parameter :: eomeg2 = 2.0D0*eomeg
      ! Solar Constant in W/m**2
      real(DP) , parameter :: solcon = 1367.0D+00
      ! Solar Constant in erg/cm**2/sec
      real(DP) , parameter :: scon = solcon*1000.0D+00

      ! Soil roughness length
      real(DP) , parameter :: zlnd = 0.01D+00
      ! Ocean roughness length
      real(DP) , parameter :: zoce = 0.00040D+00
      ! Snow roughness length
      real(DP) , parameter :: zsno = 0.00040D+00
      ! Von Karman constant
      real(DP) , parameter :: vonkar = 0.4D+00
      ! Drag coefficient for soil under canopy
      real(DP) , parameter :: csoilc = 4.0D-03
      ! Turbulent wind for stable conditions (m/sec)
      real(DP) , parameter :: wtur = 0.1D+00

      ! Constant used in computing virture temperature.
      real(DP) , parameter :: ep1 = 0.608D+00
      ! Constant used in computing saturation mixing ratio.
      ! Ratio of mean molecular weight of water (18.016 g/mole)
      ! to that of dry air (28.966 g/mole)
      real(DP) , parameter :: ep2 = 0.62197D+00
      ! Constants used in computing saturation vapor pressure.
      real(DP) , parameter :: svp1 = 0.6112D+00
      real(DP) , parameter :: svp2 = 17.67D+00
      real(DP) , parameter :: svp3 = 29.65D+00
      real(DP) , parameter :: svp4 = 0.611D+00
      real(DP) , parameter :: svp5 = 22.514D+00
      real(DP) , parameter :: svp6 = 6.15D+03
      ! Constants used in computing evaporation latent heat.
      real(DP) , parameter :: lh0 = 597.3D+00
      real(DP) , parameter :: lh1 = 0.566D+00
      ! Constants from latent heat and temperature to saturation vapor p.
      real(DP) , parameter :: lsvp1 = 6.11D+00
      real(DP) , parameter :: lsvp2 = 9.045D+00

      ! Terminal velocity constants
      real(DP) , parameter :: avt = 841.99667D+00
      real(DP) , parameter :: bvt = 0.8D+00
      real(DP) , parameter :: g4pb = 17.837825D+00
      real(DP) , parameter :: g3pb = g4pb/(3.0D+00+bvt)
      real(DP) , parameter :: g5pb = 1.8273D+00
      real(DP) , parameter :: vtc = avt*g4pb/6.0D+00
      real(DP) , parameter :: trel = 3000.0D+00

      ! Dynamic parameters
      ! alpha = .2495 in brown-campana; = 0. in split explicit
      real(DP) , parameter :: alpha = 0.0D+00
      real(DP) , parameter :: beta = 1.0D+00 - 2.0D+00*alpha
      real(DP) , parameter :: gnu = 0.10D+00
      real(DP) , parameter :: omu = 1.0D+00 - 2.0D+00*gnu
      real(DP) , parameter :: gnuhf = 0.5D+00*gnu
      real(DP) , parameter :: omuhf = 1.0D+00 - 2.0D+00*gnuhf

      ! Cumulous parameters
      real(DP) , parameter :: tauht = 7200.0D+00

      ! Aerosol densities
      real(DP) , parameter :: rhoso4 = 1.76D+00
      real(DP) , parameter :: rhobc = 1.D+00
      real(DP) , parameter :: rhooc = 1.D+00
      real(DP) , parameter :: rhodust = 2.5D+00

      ! Constants used in Betts Miller
      real(DP) , parameter :: aliq = 613.3D+00
      real(DP) , parameter :: bliq = 17.502D+00
      real(DP) , parameter :: cliq = 4780.8D+00
      real(DP) , parameter :: dliq = 32.19D+00
      real(DP) , parameter :: aice = 613.2D+00
      real(DP) , parameter :: bice = 22.452D+00
      real(DP) , parameter :: cice1 = 6133.0D+00
      real(DP) , parameter :: dice = 0.61D+00
      real(DP) , parameter :: xls0 = 2.905D+06
      real(DP) , parameter :: xls1 = 259.532D+00

      ! GTS system constants
      real(DP) , parameter :: gtigts = gti*100.0D+00
      real(DP) , parameter :: rga = 1.0D+00/gtigts
      real(DP) , parameter :: cpdgts = cpd*1.0D+04
      real(DP) , parameter :: gocp = gtigts/cpdgts
      real(DP) , parameter :: sslp = stdp*10.0D+00 ! dynes/cm^2
      real(DP) , parameter :: rsslp = 1.0D+00/sslp
      real(DP) , parameter :: stebol = sigm*1.0D+03
      real(DP) , parameter :: rgsslp = 0.5D+00/(gtigts*sslp)
      ! Effective molecular weight of dry air (g/mol)
      real(DP) , parameter :: amd = 28.9644D+00
      ! Molecular weight of ozone (g/mol)
      real(DP) , parameter :: amo = 48.0D+00
      ! Molecular weight of co2 (g/mol)
      real(DP) , parameter :: amco2 = 44.0D+00

      ! Radiation constants
      real(DP) , parameter :: dpfco2 = 5.0D-03
      real(DP) , parameter :: dpfo3 = 2.5D-03

      ! Pressure gradient force calculations (Why not standard atmosphere?)
      real(DP) , parameter :: t00pg = 287.0D+00       ! stdt ?
      real(DP) , parameter :: p00pg = 101.325D+00     ! stdp ?
      real(DP) , parameter :: alam  = 6.5D-03         ! Lapse rate ?
      real(DP) , parameter :: pgfaa1 = alam*rgas*rgti ! Utility constant

      ! Molecular heat diffusion coefficient in water
      real(DP) , parameter :: hdmw = 1.3889D-07  ! m^2/s

      end module mod_constants
