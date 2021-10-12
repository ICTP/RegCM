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
  real(rkx) , parameter :: d_zero = 0.0_rkx
  real(rkx) , parameter :: d_one = 1.0_rkx
  real(rkx) , parameter :: d_two = 2.0_rkx
  real(rkx) , parameter :: d_three = 3.0_rkx
  real(rkx) , parameter :: d_four = 4.0_rkx
  real(rkx) , parameter :: d_five = 5.0_rkx
  real(rkx) , parameter :: d_six = 6.0_rkx
  real(rkx) , parameter :: d_nine = 9.0_rkx
  real(rkx) , parameter :: d_half = 0.50_rkx
  real(rkx) , parameter :: d_rfour = 0.250_rkx
  real(rkx) , parameter :: d_twelve = 12.0_rkx
  real(rkx) , parameter :: d_60 = 60.0_rkx
  real(rkx) , parameter :: d_10 = 10.0_rkx
  real(rkx) , parameter :: d_r10 = 0.1_rkx
  real(rkx) , parameter :: d_100 = 100.0_rkx
  real(rkx) , parameter :: d_r100 = 0.01_rkx
  real(rkx) , parameter :: d_1000 = 1000.0_rkx
  real(rkx) , parameter :: d_r1000 = 0.001_rkx
  real(rkx) , parameter :: onet = d_one/d_three
  real(rkx) , parameter :: twot = d_two/d_three
  real(rkx) , parameter :: fourt = d_four/d_three

  ! Angles degrees : Operations must be always performed double precision
  real(rk8) , parameter :: deg00  = 0.0_rk8
  real(rk8) , parameter :: deg45  = 45.0_rk8
  real(rk8) , parameter :: deg90  = 90.0_rk8
  real(rk8) , parameter :: deg180 = 180.0_rk8
  real(rk8) , parameter :: deg360 = 360.0_rk8

  ! minimum values for uncoupled/coupled variables which require them
  real(rkx) , parameter :: minqq   = 1.0e-8_rkx
  real(rkx) , parameter :: minqc   = 1.0e-10_rkx
  real(rkx) , parameter :: minqv   = minqq * 100.0_rkx
#ifdef SINGLE_PRECISION_REAL
  real(rkx) , parameter :: mintr   = 1.0e-20_rkx
#else
  real(rkx) , parameter :: mintr   = 1.0e-30_rkx
#endif
  real(rkx) , parameter :: minqx   = 1.0e-16_rkx

  ! Low/Hi values
  real(rkx) , parameter :: dlowval = 1.0e-20_rkx
  real(rkx) , parameter :: dhival  = 1.0e+20_rkx
  real(rk4) , parameter :: slowval = 1.0e-20_rk4
  real(rk4) , parameter :: shival  = 1.0e+20_rk4
  real(rkx) , parameter :: dmissval = 1.0e+20_rkx
  real(rk4) , parameter :: smissval = 1.0e+20_rk4

  ! time conversion
  real(rkx) , parameter :: secpm = 60.0_rkx
  real(rkx) , parameter :: secph = 3600.0_rkx
  real(rkx) , parameter :: secpd = 86400.0_rkx
  real(rkx) , parameter :: rsecpd = 1.0_rkx/86400.0_rkx
  real(rkx) , parameter :: minph = 60.0_rkx
  real(rkx) , parameter :: minpd = 1440.0_rkx
  real(rkx) , parameter :: houpd = 24.0_rkx

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(rkx) , parameter :: egrav = 9.80665_rkx

  real(rkx) , parameter :: speedoflight = 299792458.0_rkx
  real(rkx) , parameter :: plankconstant = 6.62607550e-34_rkx
  ! Stefan-Boltzmann  constant CODATA 2007
  real(rkx) , parameter :: sigm = 5.670400e-8_rkx
  ! Boltzman Constant k CODATA 2007
  real(rkx) , parameter :: boltzk = 1.3806504e-23_rkx
  ! Avogadro Constant
  real(rkx) , parameter :: navgdr = 6.02214129e23_rkx
  ! CRC Handbook of Chemistry and Physics, 1997
  ! Effective molecular weight of dry air (g/mol)
  ! 78.084       % N2
  ! 20.9476      % O2
  !  0.00934     % Ar
  !  0.000314    % CO2
  !  0.00001818  % Ne
  !  0.000002    % CH4
  !  0.00000524  % He
  !  0.00000114  % Kr
  !  0.0000005   % H2
  !  0.000000087 % Xe
  real(rkx) , parameter :: amd   = 28.96454_rkx
  ! Effective molecular weight of water (g/mol)
  real(rkx) , parameter :: amw   = 18.01528_rkx

  ! From http://pubchem.ncbi.nlm.nih.gov
  ! Molecular weight of oxygen molecule (g/mol)
  real(rkx) , parameter :: amo2   = 31.9988_rkx
  ! Molecular weight of ozone (g/mol)
  real(rkx) , parameter :: amo3   = 47.99820_rkx
  ! Molecular weight of carbon dioxide (g/mol)
  real(rkx) , parameter :: amco2 = 44.00950_rkx
  ! Molecular weight of nitrous oxide (g/mol)
  real(rkx) , parameter :: amn2o = 44.0128_rkx
  ! Molecular weight of methane (g/mol)
  real(rkx) , parameter :: amch4 = 16.04246_rkx
  ! Molecular weight of cfc11 (g/mol)
  real(rkx) , parameter :: amcfc11 = 137.368103_rkx
  ! Molecular weight of cfc12 (g/mol)
  real(rkx) , parameter :: amcfc12 = 120.913506_rkx

  ! Ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)
  real(rkx) , parameter :: pdbratio = 0.0112372_rkx

  real(rkx) , parameter :: rgasmol = navgdr*boltzk ! 8.3144717808
  ! Gas constant for dry air
  real(rkx) , parameter :: c287 = rgasmol/amd      ! 0.2870569248
  ! Gas constant for dry air in Joules/kg/K
  real(rkx) , parameter :: rgas = c287*1000.0_rkx  ! 287.0569248
  real(rkx) , parameter :: rdry = rgas
  ! Gas constant for water vapor in Joules/kg/K
  real(rkx) , parameter :: rwat = (rgasmol/amw)*1000.0_rkx ! 461.5233169
  ! Ratio of the two above
  real(rkx) , parameter :: rgow = rgas/rwat
  ! Reverse of the above
  real(rkx) , parameter :: rgowi = rwat/rgas

  ! Specific heat at constant pressure for dry air J/kg/K
  real(rkx) , parameter :: cpd = 3.5_rkx*rgas  ! 1004.6992368
  real(rkx) , parameter :: cvd = 2.5_rkx*rgas  !  717.6423120
  ! Specific heat at constant pressure for water vapor J/kg/K
  real(rkx) , parameter :: cpv = 4.0_rkx*rwat  ! 1846.0932676000
  ! Specific heat of water at 15 Celsius J/kg/K
  real(rkx) , parameter :: cpw = 4186.95_rkx
  ! Specific heat of ice at 0 Celsius J/kg/K
  real(rkx) , parameter :: cpi = 2117.27_rkx
  ! Specific heat of water at 0 Celsius J/kg/K
  real(rkx) , parameter :: cpw0 = 4218.0_rkx
  ! Specific heat of water vapour
  real(rkx) , parameter :: cp_h2o = cpd * (4.0_rkx*amd) / (3.5_rkx*amw)

  ! Specific heats per m**3  (joules/m**3/k)
  real(rkx) , parameter :: ch2o = 4.18695e6_rkx
  real(rkx) , parameter :: shice = 0.45_rkx*ch2o
  real(rkx) , parameter :: cwi = d_one/0.45_rkx
  real(rkx) , parameter :: csnw = 0.49_rkx*ch2o
  real(rkx) , parameter :: cws = d_one/0.49_rkx

  ! Specific heats per kg (J/kg/K)
  real(rkx) , parameter :: spcpfw = 4.188e3_rkx    ! fresh h2o
  real(rkx) , parameter :: spcpsw = 3.996e3_rkx    ! Sea water
  real(rkx) , parameter :: spcpice = 2.11727e3_rkx ! fresh ice

  ! Latent heats (Joules/kg)
  ! Water vaporization latent heat at T 0 Celsius
  real(rkx) , parameter :: wlhv = 2.50080e6_rkx
  ! Water fusion latent heat at T 0 Celsius
  real(rkx) , parameter :: wlhf = 0.33355e6_rkx
  ! Water sublimation latent heat at T 0 Celsius
  real(rkx) , parameter :: wlhs = wlhv + wlhf
  ! Reverse helpers
  real(rkx) , parameter :: rwlhv = d_one/wlhv
  real(rkx) , parameter :: rwlhf = d_one/wlhf
  real(rkx) , parameter :: rwlhs = d_one/wlhs

  ! Various utility terms used in calculations
  real(rkx) , parameter :: regrav = d_one/egrav
  real(rkx) , parameter :: rcpd = d_one/cpd
  real(rkx) , parameter :: rovcp = rgas*rcpd
  real(rkx) , parameter :: rdrcv = rgas/cvd
  real(rkx) , parameter :: cpovr = cpd/rgas
  real(rkx) , parameter :: rovg  = rgas/egrav
  real(rkx) , parameter :: govr  = egrav/rgas
  real(rkx) , parameter :: gdry  = -egrav/cpd
  real(rkx) , parameter :: hcratio = cpv*rcpd
  real(rkx) , parameter :: hcrm1 = hcratio - d_one
  real(rkx) , parameter :: rhoh2o = 1000.0_rkx
  real(rkx) , parameter :: rhosea = 1026.0_rkx
  real(rkx) , parameter :: rhosnow = 100.0_rkx
  real(rkx) , parameter :: rhosnowp = 330.0_rkx
  real(rkx) , parameter :: rhoice = 917.0_rkx
  real(rkx) , parameter :: tzero = 273.15_rkx
  real(rkx) , parameter :: tiso =  216.65_rkx
  real(rkx) , parameter :: rtzero = d_one/tzero
  real(rkx) , parameter :: wattp = 273.16_rkx
  real(rkx) , parameter :: tboil = 373.1339_rkx
  real(rkx) , parameter :: c1es = 610.78_rkx
  real(rkx) , parameter :: c2es = c1es*amw/amd
  real(rkx) , parameter :: c3les = 17.2693882_rkx
  real(rkx) , parameter :: c3ies = 21.875_rkx
  real(rkx) , parameter :: c4les = 35.86_rkx
  real(rkx) , parameter :: c4ies = 7.66_rkx
  real(rkx) , parameter :: c5les = c3les*(tzero-c4les)
  real(rkx) , parameter :: c5ies = c3ies*(tzero-c4ies)
  real(rkx) , parameter :: c5alvcp = c5les*wlhv*rcpd
  real(rkx) , parameter :: c5alscp = c5ies*wlhs*rcpd
  real(rkx) , parameter :: wlhvocp = wlhv*rcpd
  real(rkx) , parameter :: wlhsocp = wlhs*rcpd
  real(rkx) , parameter :: wlhfocp = wlhf*rcpd
  real(rkx) , parameter :: cpowlhv = cpd*rwlhv
  real(rkx) , parameter :: cpowlhs = cpd*rwlhs
  real(rkx) , parameter :: cpowlhf = cpd*rwlhf
  real(rkx) , parameter :: rtber = tzero-5.0_rkx
  real(rkx) , parameter :: rtice = tzero-23.0_rkx
  real(rkx) , parameter :: rtwat = tzero
  real(rkx) , parameter :: mpcrt = rtice+(rtwat-rtice)/sqrt(2.0_rkx)
  real(rkx) , parameter :: rtwat_rtice_r = 1.0_rkx/(rtwat-rtice)
  real(rkx) , parameter :: pq0 = 379.90516_rkx
  ! value used for the latent heat term in the exponent for
  ! calculating equivalent potential temperature
  real(rkx) , parameter :: eliwv = 2.72e6_rkx

  ! Standard atmosphere ICAO 1993
  real(rkx) , parameter :: p00 = 1.000000e5_rkx
  real(rkx) , parameter :: stdp = 1.013250e5_rkx
  real(rkx) , parameter :: stdpmb = 1013.250_rkx
  real(rkx) , parameter :: stdpcb = 101.3250_rkx
  real(rkx) , parameter :: stdt = 288.15_rkx
  real(rkx) , parameter :: stdrho = 1.28_rkx
  real(rkx) , parameter :: lrate = 0.00649_rkx ! K/m from MSL up to 11 km

  ! Atmos. surface pressure mol/cm3
  real(rkx) , parameter :: atmos = 2.247e19_rkx
  ! Conversion parameter for Henry L-atm/mol-K
  real(rkx) , parameter :: rtcon = 8.314e-2_rkx
  ! RU g-cm2/s2-mol-K
  real(rkx) , parameter :: rumolec = 8.314e7_rkx
  ! Droplet diffusion coefficient Lelieveld+Crutzen 1991 cm2/sec
  real(rkx) , parameter :: dropdif = 2.0e-5_rkx
  ! Gas-phase diffusion coeff. Lelieveld and Crutzen, 1991 cm2/s
  real(rkx) , parameter :: difgas = 0.1_rkx

  ! Fixed emissivity of water
  real(rkx) , parameter :: emsw = 0.97_rkx

  real(rk8) , parameter :: m_euler = 0.577215664901532860606512090082402431_rk8

  ! Trigonometric constants.
  real(rk8) , parameter :: mathpi =                                   &
                      &   3.1415926535897932384626433832795029_rk8
  real(rk8) , parameter :: invpi = d_one/mathpi
  real(rk8) , parameter :: halfpi = mathpi*d_half
  real(rk8) , parameter :: quartpi = halfpi*d_half
  real(rk8) , parameter :: twopi = mathpi*d_two
  real(rk8) , parameter :: pisqr = mathpi*mathpi
  real(rk8) , parameter :: degrad = mathpi/180.0_rk8
  real(rk8) , parameter :: raddeg = 180.0_rk8/mathpi

  ! Maximum stomatl resistance (s/m)
  real(rkx) , parameter :: rmax0 = 2.0e4_rkx

  ! Maximum allowed dew(mm) and inverse (dewmaxi)
  real(rkx) , parameter :: dewmax = 0.1_rkx
  real(rkx) , parameter :: dewmaxi = d_one/dewmax

  ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
  real(rkx) , parameter :: trsmx0 = 2.e-4_rkx

  ! Drainage out of 10m layer bottom (mm/s)
  ! drain is set fairly large to prevent swamping the soil
  real(rkx) , parameter :: drain = 4.0e-4_rkx

  ! Minimum ratio between potential and actual water content
  real(rkx) , parameter :: minwrat = 1.0e-4_rkx

  ! Earth radius in meters
  real(rk8) , parameter :: earthrad = 6.371229e6_rk8
  real(rk8) , parameter :: erkm = earthrad/1000.0_rk8
  real(rk8) , parameter :: rearthrad = d_one/earthrad
  ! Angular velocity of rotation of Earth
  real(rk8) , parameter :: eomeg = 7.2921159e-5_rk8
  real(rk8) , parameter :: eomeg2 = d_two*eomeg

  ! Soil roughness length
  real(rkx) , parameter :: zlnd = 0.01_rkx
  ! Ocean roughness length
  real(rkx) , parameter :: zoce = 0.00023_rkx
  ! Snow roughness length
  real(rkx) , parameter :: zsno = 0.00040_rkx
  ! Von Karman constant
  real(rkx) , parameter :: vonkar = 0.4_rkx
  ! Drag coefficient for soil under canopy
  real(rkx) , parameter :: csoilc = 4.0e-3_rkx
  ! Turbulent wind for stable conditions (m/sec)
  real(rkx) , parameter :: wtur = 0.1_rkx

  ! Constant used in computing virtual temperature.
  real(rkx) , parameter :: ep1 = amd/amw - d_one ! 0.6077762876
  ! Constant used in computing saturation mixing ratio.
  ! Ratio of mean molecular weight of water to that of dry air
  real(rkx) , parameter :: ep2 = amw/amd   ! 0.6219770795
  real(rkx) , parameter :: rep2 = amd/amw  ! 1.6077762876

  ! Terminal velocity constants
  real(rkx) , parameter :: avt = 841.99667_rkx
  real(rkx) , parameter :: bvt = 0.8_rkx
  real(rkx) , parameter :: g4pb = 17.837825_rkx
  real(rkx) , parameter :: g3pb = g4pb/(3.0_rkx+bvt)
  real(rkx) , parameter :: g5pb = 1.8273_rkx
  real(rkx) , parameter :: vtc = avt*g4pb/6.0_rkx

  ! Dynamic parameters
  ! alpha = .2495 in brown-campana; = 0. in split explicit
  real(rkx) , parameter :: alpha_hyd = 0.0_rkx
  real(rkx) , parameter :: beta_hyd = d_one - d_two*alpha_hyd

  ! Constant surface Long Wave emissivity
  real(rkx) , parameter :: lnd_sfcemiss = 0.985_rkx
  real(rkx) , parameter :: ocn_sfcemiss = 0.984_rkx

  ! Constants used in Kain-Fritsch and WSM5
  real(rkx) , parameter :: aliq = 613.3_rkx
  real(rkx) , parameter :: bliq = 17.502_rkx
  real(rkx) , parameter :: cliq = 4780.8_rkx
  real(rkx) , parameter :: dliq = 32.19_rkx
  real(rkx) , parameter :: aice = 613.20_rkx
  real(rkx) , parameter :: bice = 22.452_rkx
  real(rkx) , parameter :: cice = 6133.0_rkx
  real(rkx) , parameter :: dice = 0.61_rkx

  ! GTS system constants
  real(rkx) , parameter :: egravgts = egrav*d_100
  real(rkx) , parameter :: regravgts = d_one/egravgts
  real(rkx) , parameter :: cpdgts = cpd*1.0e4_rkx
  real(rkx) , parameter :: gocp = egravgts/cpdgts
  real(rkx) , parameter :: sslp = stdp*d_10 ! dynes/cm^2
  real(rkx) , parameter :: rsslp = d_one/sslp
  real(rkx) , parameter :: stebol = sigm*d_1000
  real(rkx) , parameter :: rgsslp = d_half/(egravgts*sslp)
  ! Effective molecular weight of dry air (kg/mol)
  real(rkx) , parameter :: amdk = amd*d_r1000
  ! Avogadro Constant in lit/cm3
  real(rkx) , parameter :: avogadrl = navgdr*d_1000

  ! Radiation constants
  real(rkx) , parameter :: dpfco2 = 5.0e-3_rkx
  real(rkx) , parameter :: dpfo3 = 2.5e-3_rkx

  ! Pressure gradient force calculations (Why not standard atmosphere?)
  real(rkx) , parameter :: t00pg = 287.0_rkx       ! stdt ?
  real(rkx) , parameter :: p00pg = 101.325_rkx     ! stdpcb ?
  real(rkx) , parameter :: alam  = 6.5e-3_rkx      ! Lapse rate ?
  real(rkx) , parameter :: pgfaa1 = alam*rgas*regrav ! Utility constant

  ! Molecular heat diffusion coefficient in water
  real(rkx) , parameter :: hdmw = 1.3889e-7_rkx  ! m^2/s

  ! Seaice minimum depth value
  real(rkx) , parameter :: iceminh = 0.01_rkx

  ! Allowed range for cloud fraction
  real(rkx) , parameter :: lowcld = 0.0001_rkx
  real(rkx) , parameter :: hicld  = 0.9999_rkx

end module mod_constants
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
