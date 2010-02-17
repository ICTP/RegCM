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

      ! Gravity (m/sec**2)
      real(8) , parameter :: gti = 9.80665D0
      real(8) , parameter :: rgti = 1.0D0/gti

      ! Specific heat at constant pressure for dry air
      real(8) , parameter :: cpd = 1005.46D0
      real(8) , parameter :: rcpd = 1.D0/cpd
      ! Specific heat at constant volume for dry air
      real(8) , parameter :: cpv = 1869.46D0
      ! Specific heat of water (salt)
      real(8) , parameter :: cpw = 4179.98D0

      ! Gas constant for dry air in Joules/kg/K
      real(8) , parameter :: rgas = 287.04D0
      ! Gas constant for water vapor in Joules/kg/K
      real(8) , parameter :: rwat = 461.90D0

      real(8) , parameter :: rovcp = rgas*rcpd
      real(8) , parameter :: rovg  = rgas/gti

      ! Various utility terms used in calculations
      real(8) , parameter :: vtmpc1 = rwat/rgas - 1.0D0
      real(8) , parameter :: vtmpc2 = cpv*rcpd - 1.0D0
      real(8) , parameter :: rhoh2o = 1000.0D0
      real(8) , parameter :: alv = 2.5008D6
      real(8) , parameter :: als = 2.8345D6
      real(8) , parameter :: alf = als - alv
      real(8) , parameter :: tmelt = 273.16D0
      real(8) , parameter :: tboil = 373.16D0
      real(8) , parameter :: c1es = 610.78D0
      real(8) , parameter :: c2es = c1es*rgas/rwat
      real(8) , parameter :: c3les = 17.2693882D0
      real(8) , parameter :: c3ies = 21.875D0
      real(8) , parameter :: c4les = 35.86D0
      real(8) , parameter :: c4ies = 7.66D0
      real(8) , parameter :: c5les = c3les*(tmelt-c4les)
      real(8) , parameter :: c5ies = c3ies*(tmelt-c4ies)
      real(8) , parameter :: c5alvcp = c5les*alv*rcpd
      real(8) , parameter :: c5alscp = c5ies*als*rcpd
      real(8) , parameter :: alvdcp = alv*rcpd
      real(8) , parameter :: alsdcp = als*rcpd
      real(8) , parameter :: pq0 = 379.90516D0
      ! value used for the latent heat term in the exponent fori
      ! calculating equivalent potential temperature
      real(8) , parameter :: eliwv = 2.72D6

      ! Standard atmosphere
      real(8) , parameter :: stdp = 1.013250D5
      real(8) , parameter :: stdpmb = 1013.250D5
      real(8) , parameter :: stdt = 288.15D0
      real(8) , parameter :: lrate = 0.0065D0
 
      ! Stefan-Boltzmann  constant
      real(8) , parameter :: sigm = 5.76383D-8

      ! Trigonometric constants
      real(8) , parameter :: mathpi = 3.14159265358979323846D0
      real(8) , parameter :: twopi = mathpi*2.0D0
      real(8) , parameter :: degrad = mathpi/180.0

      ! Specific heats per m**3  (joules/m**3/k)
      real(8) , parameter :: ch2o = 4.186D6
      real(8) , parameter :: cice = 0.45D0*ch2o
      real(8) , parameter :: cwi = 1.0D0/0.45D0
      real(8) , parameter :: csnw = 0.49D0*ch2o
      real(8) , parameter :: cws = 1.0D0/0.49D0

      ! Latent heats (Joules/kg)
      real(8) , parameter :: wlhf = 0.3336D6
      real(8) , parameter :: wlhv = 2.51040D6
      real(8) , parameter :: wlhs = wlhv + wlhf
      real(8) , parameter :: wlhvocp = wlhv*rcpd

      ! Maximum stomatl resistance (s/m)
      real(8) , parameter :: rmax0 = 2.0D4

      ! Maximum allowed dew(mm) and inverse (dewmxi)
      real(8) , parameter :: dewmx = 0.1D0
      real(8) , parameter :: dewmxi = 1.0D0/dewmx

      ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
      real(8) , parameter :: trsmx0 = 2.D-4

      ! Drainage out of 10m layer bottom (mm/s)
      ! drain is set fairly large to prevent swamping the soil
      real(8) , parameter :: drain = 1.D-4

      ! Length of day in seconds
      real(8) , parameter :: tau1 = 8.64D4
      ! Days per year
      real(8) , parameter :: dayspy = 365.24
      ! Degrees per day
      real(8) , parameter :: dpd = 360.0/dayspy
      ! Angular velocity of rotation of Earth
      real(8) , parameter :: eomeg = 7.292D-5
      ! Solar Constant in W/m**2
      real(8) , parameter :: solcon = 1.3956D3

      ! Soil roughness length
      real(8) , parameter :: zlnd = 0.01D0
      ! Ocean roughness length
      real(8) , parameter :: zoce = 0.00040D0
      ! Snow roughness length
      real(8) , parameter :: zsno = 0.00040D0
      ! Von Karman constant
      real(8) , parameter :: vonkar = 0.378D0
      ! Drag coefficient for soil under canopy
      real(8) , parameter :: csoilc = 4.0D-3
      ! Turbulent wind for stable conditions (m/sec)
      real(8) , parameter :: wtur = 0.1D0

      ! Constant used in computing virture temperature.
      real(8) , parameter :: ep1 = 0.608D0
      ! Constant used in computing saturation mixing ratio.
      ! Ratio of mean molecular weight of water (18.016 g/mole)
      ! to that of dry air (28.966 g/mole)
      real(8) , parameter :: ep2 = 0.62197D0
      ! Constants used in computing saturation vapor pressure.
      real(8) , parameter :: svp1 = 0.6112D0
      real(8) , parameter :: svp2 = 17.67D0
      real(8) , parameter :: svp3 = 29.65D0

      ! Terminal velocity constants
      real(8) , parameter :: avt = 841.99667D0
      real(8) , parameter :: bvt = 0.8D0
      real(8) , parameter :: g4pb = 17.837825D0
      real(8) , parameter :: g3pb = g4pb/(3.0D0+bvt)
      real(8) , parameter :: g5pb = 1.8273D0
      real(8) , parameter :: vtc = avt*g4pb/6.0D0
      real(8) , parameter :: trel = 3000.0D0
      ! Marshall Palmer
      real(8) , parameter :: n0r = 8.0D6
      real(8) , parameter :: ppi = 1.0D0/(mathpi*n0r)
      real(8) , parameter :: prac = mathpi*n0r*avt*g3pb*0.25D0
      real(8) , parameter :: prec1 = 2.0D0*mathpi*n0r*0.78D0
      real(8) , parameter :: prec2 = 2.0D0*mathpi*n0r*0.32D0*           &
                                   & avt**0.5D0*g5pb

      ! Dynamic parameters
      ! alpha = .2495 in brown-campana; = 0. in split explicit
      real(8) , parameter :: alpha = 0.0D0
      real(8) , parameter :: beta = 1.0D0 - 2.0D0*alpha
      real(8) , parameter :: gnu = 0.10D0
      real(8) , parameter :: omu = 1.0D0 - 2.0D0*gnu
      real(8) , parameter :: gnuhf = 0.5D0*gnu
      real(8) , parameter :: omuhf = 1.0D0 - 2.0D0*gnuhf

      ! Cumulous parameters
      real(8) , parameter :: tauht = 7200.0D0

      ! Aerosol densities
      real(8) , parameter :: rhoso4 = 1.76D0
      real(8) , parameter :: rhobc = 1.D0
      real(8) , parameter :: rhooc = 1.D0
      real(8) , parameter :: rhodust = 2.5D0

      ! Constants used in Betts Miller
      real(8) , parameter :: aliq = 613.3D0
      real(8) , parameter :: bliq = 17.502D0
      real(8) , parameter :: cliq = 4780.8D0
      real(8) , parameter :: dliq = 32.19D0
      real(8) , parameter :: aice = 613.2D0
      real(8) , parameter :: bice = 22.452D0
      real(8) , parameter :: cice1 = 6133.0D0
      real(8) , parameter :: dice = 0.61D0
      real(8) , parameter :: xls0 = 2.905D6
      real(8) , parameter :: xls1 = 259.532D0

      ! GTS system constants
      real(8) , parameter :: gtigts = gti*100.0D0
      real(8) , parameter :: rga = 1.0D0/gtigts
      real(8) , parameter :: cpdgts = cpd*1.0D4
      real(8) , parameter :: gocp = gtigts/cpdgts
      real(8) , parameter :: sslp = stdp*10.0D0
      real(8) , parameter :: rsslp = 1.0D0/sslp
      real(8) , parameter :: stebol = sigm*1.0D3
      real(8) , parameter :: rgsslp = 0.5/(gtigts*sslp)

      ! Radiation constants
      real(8) , parameter :: dpfco2 = 5.0D-3
      real(8) , parameter :: dpfo3 = 2.5D-3

      end module mod_constants
