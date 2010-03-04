!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_constants

      implicit none

      ! Standard Gravity (m/sec**2) 3rd CGPM
      real(8) , parameter :: gti = 9.80665D+00

      ! Gas constant for dry air in Joules/kg/K
      real(8) , parameter :: rgas = 287.0058D+00
      ! Gas constant for water vapor in Joules/kg/K
      real(8) , parameter :: rwat = 461.90D+00

      ! Specific heat at constant pressure for dry air J/kg/K
      real(8) , parameter :: cpd = 1005.46D+00
      ! Specific heat at constant pressure for moist air J/kg/K
      real(8) , parameter :: cpv = 1869.46D+00
      ! Specific heat of water at 15 Celsius J/kg/K
      real(8) , parameter :: cpw = 4186.95D+00

      ! Various utility terms used in calculations
      real(8) , parameter :: rgti = 1.0D+00/gti
      real(8) , parameter :: rcpd = 1.0D+00/cpd
      real(8) , parameter :: rovcp = rgas*rcpd
      real(8) , parameter :: rovg  = rgas/gti
      real(8) , parameter :: govr  = gti/rgas
      real(8) , parameter :: vtmpc1 = rwat/rgas - 1.0D+00
      real(8) , parameter :: vtmpc2 = cpv*rcpd - 1.0D+00
      real(8) , parameter :: rhoh2o = 1000.0D+00
      real(8) , parameter :: alv = 2.5008D+06
      real(8) , parameter :: als = 2.8345D+06
      real(8) , parameter :: alf = als - alv
      real(8) , parameter :: tzero = 273.15D+00
      real(8) , parameter :: wattp = 273.16D+00
      real(8) , parameter :: tboil = 373.1339D+00
      real(8) , parameter :: c1es = 610.78D+00
      real(8) , parameter :: c2es = c1es*rgas/rwat
      real(8) , parameter :: c3les = 17.2693882D+00
      real(8) , parameter :: c3ies = 21.875D+00
      real(8) , parameter :: c4les = 35.86D+00
      real(8) , parameter :: c4ies = 7.66D+00
      real(8) , parameter :: c5les = c3les*(tzero-c4les)
      real(8) , parameter :: c5ies = c3ies*(tzero-c4ies)
      real(8) , parameter :: c5alvcp = c5les*alv*rcpd
      real(8) , parameter :: c5alscp = c5ies*als*rcpd
      real(8) , parameter :: alvdcp = alv*rcpd
      real(8) , parameter :: alsdcp = als*rcpd
      real(8) , parameter :: pq0 = 379.90516D+00
      ! value used for the latent heat term in the exponent for
      ! calculating equivalent potential temperature
      real(8) , parameter :: eliwv = 2.72D+06

      ! Standard atmosphere ICAO 1993
      real(8) , parameter :: stdp = 1.013250D+05
      real(8) , parameter :: stdpmb = 1013.250D+00
      real(8) , parameter :: stdt = 288.15D+00
      real(8) , parameter :: lrate = 0.00649D+00 ! K/km from MSL up to 11 km
 
      ! Stefan-Boltzmann  constant CODATA 2007
      real(8) , parameter :: sigm = 5.670400D-08

      ! Trigonometric constants. 
      real(8) , parameter :: mathpi =                                   &
                          &   3.1415926535897932384626433832795029D+00
      real(8) , parameter :: invpi = 1.0D+00/mathpi
      real(8) , parameter :: halfpi = mathpi*0.5D+00
      real(8) , parameter :: twopi = mathpi*2.0D+00
      real(8) , parameter :: degrad = mathpi/180.0D+00
      real(8) , parameter :: raddeg = 180.0D+00/mathpi

      ! Specific heats per m**3  (joules/m**3/k)
      real(8) , parameter :: ch2o = 4.186D+06
      real(8) , parameter :: cice = 0.45D+00*ch2o
      real(8) , parameter :: cwi = 1.0D+00/0.45D+00
      real(8) , parameter :: csnw = 0.49D+00*ch2o
      real(8) , parameter :: cws = 1.0D+00/0.49D+00

      ! Latent heats (Joules/kg)
      real(8) , parameter :: wlhf = 0.3336D+06
      real(8) , parameter :: wlhv = 2.51040D+06
      real(8) , parameter :: wlhs = wlhv + wlhf
      real(8) , parameter :: wlhvocp = wlhv*rcpd

      ! Maximum stomatl resistance (s/m)
      real(8) , parameter :: rmax0 = 2.0D+04

      ! Maximum allowed dew(mm) and inverse (dewmxi)
      real(8) , parameter :: dewmx = 0.1D+00
      real(8) , parameter :: dewmxi = 1.0D+00/dewmx

      ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
      real(8) , parameter :: trsmx0 = 2.D-04

      ! Drainage out of 10m layer bottom (mm/s)
      ! drain is set fairly large to prevent swamping the soil
      real(8) , parameter :: drain = 1.D-04

      ! Earth radius
      real(8) , parameter :: earthrad = 6.371229D+06
      ! Length of day in seconds
      real(8) , parameter :: tau1 = 8.64D+04
      ! Days per year
      real(8) , parameter :: dayspy = 365.24D+00
      ! Degrees per day
      real(8) , parameter :: dpd = 360.0/dayspy
      ! Angular velocity of rotation of Earth
      real(8) , parameter :: eomeg = 7.2921159D-05
      real(8) , parameter :: eomeg2 = 2.0D0*eomeg
      ! Solar Constant in W/m**2
      real(8) , parameter :: solcon = 1.3956D+03

      ! Soil roughness length
      real(8) , parameter :: zlnd = 0.01D+00
      ! Ocean roughness length
      real(8) , parameter :: zoce = 0.00040D+00
      ! Snow roughness length
      real(8) , parameter :: zsno = 0.00040D+00
      ! Von Karman constant
      real(8) , parameter :: vonkar = 0.378D+00
      ! Drag coefficient for soil under canopy
      real(8) , parameter :: csoilc = 4.0D-03
      ! Turbulent wind for stable conditions (m/sec)
      real(8) , parameter :: wtur = 0.1D+00

      ! Constant used in computing virture temperature.
      real(8) , parameter :: ep1 = 0.608D+00
      ! Constant used in computing saturation mixing ratio.
      ! Ratio of mean molecular weight of water (18.016 g/mole)
      ! to that of dry air (28.966 g/mole)
      real(8) , parameter :: ep2 = 0.62197D+00
      ! Constants used in computing saturation vapor pressure.
      real(8) , parameter :: svp1 = 0.6112D+00
      real(8) , parameter :: svp2 = 17.67D+00
      real(8) , parameter :: svp3 = 29.65D+00
      real(8) , parameter :: svp4 = 0.611D+00
      real(8) , parameter :: svp5 = 22.514D+00
      real(8) , parameter :: svp6 = 6.15D+03
      ! Constants used in computing evaporation latent heat.
      real(8) , parameter :: lh0 = 597.3D+00
      real(8) , parameter :: lh1 = 0.566D+00
      ! Constants from latent heat and temperature to saturation vapor p.
      real(8) , parameter :: lsvp1 = 6.11D+00
      real(8) , parameter :: lsvp2 = 9.045D+00

      ! Terminal velocity constants
      real(8) , parameter :: avt = 841.99667D+00
      real(8) , parameter :: bvt = 0.8D+00
      real(8) , parameter :: g4pb = 17.837825D+00
      real(8) , parameter :: g3pb = g4pb/(3.0D+00+bvt)
      real(8) , parameter :: g5pb = 1.8273D+00
      real(8) , parameter :: vtc = avt*g4pb/6.0D+00
      real(8) , parameter :: trel = 3000.0D+00
      ! Marshall Palmer
      real(8) , parameter :: n0r = 8.0D+06
      real(8) , parameter :: ppi = 1.0D+00/(mathpi*n0r)
      real(8) , parameter :: prac = mathpi*n0r*avt*g3pb*0.25D+00
      real(8) , parameter :: prec1 = 2.0D+00*mathpi*n0r*0.78D+00
      real(8) , parameter :: prec2 = 2.0D+00*mathpi*n0r*0.32D+00*       &
                                   & avt**0.5D+00*g5pb

      ! Dynamic parameters
      ! alpha = .2495 in brown-campana; = 0. in split explicit
      real(8) , parameter :: alpha = 0.0D+00
      real(8) , parameter :: beta = 1.0D+00 - 2.0D+00*alpha
      real(8) , parameter :: gnu = 0.10D+00
      real(8) , parameter :: omu = 1.0D+00 - 2.0D+00*gnu
      real(8) , parameter :: gnuhf = 0.5D+00*gnu
      real(8) , parameter :: omuhf = 1.0D+00 - 2.0D+00*gnuhf

      ! Cumulous parameters
      real(8) , parameter :: tauht = 7200.0D+00

      ! Aerosol densities
      real(8) , parameter :: rhoso4 = 1.76D+00
      real(8) , parameter :: rhobc = 1.D+00
      real(8) , parameter :: rhooc = 1.D+00
      real(8) , parameter :: rhodust = 2.5D+00

      ! Constants used in Betts Miller
      real(8) , parameter :: aliq = 613.3D+00
      real(8) , parameter :: bliq = 17.502D+00
      real(8) , parameter :: cliq = 4780.8D+00
      real(8) , parameter :: dliq = 32.19D+00
      real(8) , parameter :: aice = 613.2D+00
      real(8) , parameter :: bice = 22.452D+00
      real(8) , parameter :: cice1 = 6133.0D+00
      real(8) , parameter :: dice = 0.61D+00
      real(8) , parameter :: xls0 = 2.905D+06
      real(8) , parameter :: xls1 = 259.532D+00

      ! GTS system constants
      real(8) , parameter :: gtigts = gti*100.0D+00
      real(8) , parameter :: rga = 1.0D+00/gtigts
      real(8) , parameter :: cpdgts = cpd*1.0D+04
      real(8) , parameter :: gocp = gtigts/cpdgts
      real(8) , parameter :: sslp = stdp*10.0D+00
      real(8) , parameter :: rsslp = 1.0D+00/sslp
      real(8) , parameter :: stebol = sigm*1.0D+03
      real(8) , parameter :: rgsslp = 0.5D+00/(gtigts*sslp)

      ! Radiation constants
      real(8) , parameter :: dpfco2 = 5.0D-03
      real(8) , parameter :: dpfo3 = 2.5D-03

      ! Pressure gradient force calculations (Why not standard atmosphere?)
      real(8) , parameter :: t00pg = 287.0D+00       ! stdt ?
      real(8) , parameter :: p00pg = 101.325D+00     ! stdp ?
      real(8) , parameter :: alam  = 6.5D-03         ! Lapse rate ?
      real(8) , parameter :: pgfaa1 = alam*rgas*rgti ! Utility constant

      end module mod_constants
