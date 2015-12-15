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
  real(rk8) , parameter :: d_zero = 0.0D+00
  real(rk8) , parameter :: d_one = 1.0D+00
  real(rk8) , parameter :: d_two = 2.0D+00
  real(rk8) , parameter :: d_three = 3.0D+00
  real(rk8) , parameter :: d_four = 4.0D+00
  real(rk8) , parameter :: d_five = 5.0D+00
  real(rk8) , parameter :: d_six = 6.0D+00
  real(rk8) , parameter :: d_nine = 6.0D+00
  real(rk8) , parameter :: d_half = 0.50D+00
  real(rk8) , parameter :: d_rfour = 0.250D+00
  real(rk8) , parameter :: d_twelve = 12.0D+00
  real(rk8) , parameter :: d_60 = 60.0D+00
  real(rk8) , parameter :: d_10 = 1.0D+01
  real(rk8) , parameter :: d_r10 = 1.0D-01
  real(rk8) , parameter :: d_100 = 1.0D+02
  real(rk8) , parameter :: d_r100 = 1.0D-02
  real(rk8) , parameter :: d_1000 = 1.0D+03
  real(rk8) , parameter :: d_r1000 = 1.0D-03
  real(rk8) , parameter :: onet = d_one/d_three
  real(rk8) , parameter :: twot = d_two/d_three
  real(rk8) , parameter :: fourt = d_four/d_three

  ! Angles degrees
  real(rk8) , parameter :: deg00  = 0.0D+00
  real(rk8) , parameter :: deg45  = 45.0D+00
  real(rk8) , parameter :: deg90  = 90.0D+00
  real(rk8) , parameter :: deg180 = 180.0D+00
  real(rk8) , parameter :: deg360 = 360.0D+00

  ! minimum values for uncoupled/coupled variables which require them
  real(rk8) , parameter :: minqq   = 1.0D-8
  real(rk8) , parameter :: minqv   = 1.0D-10
  real(rk8) , parameter :: minqx   = 1.0D-12
  real(rk8) , parameter :: minww   = 1.0D-20
  real(rk8) , parameter :: mintr   = 1.0D-50

  ! Low/Hi values
  real(rk8) , parameter :: dlowval = 1.0D-60
  real(rk8) , parameter :: dhival  = 1.0D+60
  real(rk4) , parameter :: slowval = 1.0E-30
  real(rk4) , parameter :: shival  = 1.0E+30
  real(rk8) , parameter :: dmissval = 1.0D+20
  real(rk4) , parameter :: smissval = 1.0E+20

  ! time conversion
  real(rk8) , parameter :: secpm = 60.0D+00
  real(rk8) , parameter :: secph = 3600.0D+00
  real(rk8) , parameter :: secpd = 86400.0D+00
  real(rk8) , parameter :: rsecpd = 1.0D0/86400.0D+00
  real(rk8) , parameter :: minph = 60.0D+00
  real(rk8) , parameter :: minpd = 1440.0D+00
  real(rk8) , parameter :: houpd = 24.0D+00

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(rk8) , parameter :: egrav = 9.80665D+00

  real(rk8) , parameter :: speedoflight = 299792458.0D+00
  real(rk8) , parameter :: plankconstant = 6.62607550D-34
  ! Stefan-Boltzmann  constant CODATA 2007
  real(rk8) , parameter :: sigm = 5.670400D-08
  ! Boltzman Constant k CODATA 2007
  real(rk8) , parameter :: boltzk = 1.3806504D-23
  ! Avogadro Constant
  real(rk8) , parameter :: navgdr = 6.02214129D23
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
  real(rk8) , parameter :: amd   = 28.96454D+00
  ! Effective molecular weight of water (g/mol)
  real(rk8) , parameter :: amw   = 18.01528D+00

  ! From http://pubchem.ncbi.nlm.nih.gov
  ! Molecular weight of oxygen molecule (g/mol)
  real(rk8) , parameter :: amo2   = 31.9988D+00
  ! Molecular weight of ozone (g/mol)
  real(rk8) , parameter :: amo3   = 47.99820D+00
  ! Molecular weight of carbon dioxide (g/mol)
  real(rk8) , parameter :: amco2 = 44.00950D+00
  ! Molecular weight of nitrous oxide (g/mol)
  real(rk8) , parameter :: amn2o = 44.0128D+00
  ! Molecular weight of methane (g/mol)
  real(rk8) , parameter :: amch4 = 16.04246D+00
  ! Molecular weight of cfc11 (g/mol)
  real(rk8) , parameter :: amcfc11 = 137.368103D+00
  ! Molecular weight of cfc12 (g/mol)
  real(rk8) , parameter :: amcfc12 = 120.913506D+00

  ! Ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)
  real(rk8) , parameter :: pdbratio = 0.0112372D+00

  real(rk8) , parameter :: rgasmol = navgdr*boltzk ! 8.3144717808
  ! Gas constant for dry air
  real(rk8) , parameter :: c287 = rgasmol/amd
  ! Gas constant for dry air in Joules/kg/K
  real(rk8) , parameter :: rgas = c287*1000.0D+00  ! 287.0569248
  real(rk8) , parameter :: rdry = rgas
  ! Gas constant for water vapor in Joules/kg/K
  real(rk8) , parameter :: rwat = (rgasmol/amw)*1000.0D+00 ! 461.5233169
  ! Ratio of the two above
  real(rk8) , parameter :: rgow = rgas/rwat
  ! Reverse of the above
  real(rk8) , parameter :: rgowi = rwat/rgas

  ! Specific heat at constant pressure for dry air J/kg/K
  real(rk8) , parameter :: cpd = 3.5D+00*rgas  ! 1004.6992368000
  ! Specific heat at constant pressure for moist air J/kg/K
  real(rk8) , parameter :: cpv = 4.0D+00*rwat  ! 1846.0932676000
  ! Specific heat of water at 15 Celsius J/kg/K
  real(rk8) , parameter :: cpw = 4186.95D+00
  ! Specific heat of water at 0 Celsius J/kg/K
  real(rk8) , parameter :: cpw0 = 4218.0D+00
  ! Specific heat of water vapour
  real(rk8) , parameter :: cp_h2o = cpd * (4.0D0*amd) / (3.5D0*amw)

  ! Specific heats per m**3  (joules/m**3/k)
  real(rk8) , parameter :: ch2o = 4.18695D+06
  real(rk8) , parameter :: shice = 0.45D+00*ch2o
  real(rk8) , parameter :: cwi = d_one/0.45D+00
  real(rk8) , parameter :: csnw = 0.49D+00*ch2o
  real(rk8) , parameter :: cws = d_one/0.49D+00

  ! Specific heats per kg (J/kg/K)
  real(rk8) , parameter :: spcpfw = 4.188D+03    ! fresh h2o
  real(rk8) , parameter :: spcpsw = 3.996D+03    ! Sea water
  real(rk8) , parameter :: spcpice = 2.11727D+03 ! fresh ice

  ! Latent heats (Joules/kg)
  real(rk8) , parameter :: wlhf = 0.33355D+06
  real(rk8) , parameter :: wlhv = 2.50080D+06
  real(rk8) , parameter :: wlhs = wlhv + wlhf
  real(rk8) , parameter :: rwlhv = d_one/wlhv
  real(rk8) , parameter :: rwlhf = d_one/wlhf
  real(rk8) , parameter :: rwlhs = d_one/wlhs

  ! Various utility terms used in calculations
  real(rk8) , parameter :: regrav = d_one/egrav
  real(rk8) , parameter :: rcpd = d_one/cpd
  real(rk8) , parameter :: rovcp = rgas*rcpd
  real(rk8) , parameter :: rovg  = rgas/egrav
  real(rk8) , parameter :: govr  = egrav/rgas
  real(rk8) , parameter :: gdry  = -egrav/cpd
  real(rk8) , parameter :: hcratio = cpv*rcpd
  real(rk8) , parameter :: hcrm1 = hcratio - d_one
  real(rk8) , parameter :: rhoh2o = 1000.0D+00
  real(rk8) , parameter :: rhosea = 1026.0D+00
  real(rk8) , parameter :: rhosnow = 330.0D+00
  real(rk8) , parameter :: rhoice = 917.0D+00
  real(rk8) , parameter :: tzero = 273.15D+00
  real(rk8) , parameter :: tiso = 216.65D+00
  real(rk8) , parameter :: rtzero = d_one/tzero
  real(rk8) , parameter :: wattp = 273.16D+00
  real(rk8) , parameter :: tboil = 373.1339D+00
  real(rk8) , parameter :: c1es = 610.78D+00
  real(rk8) , parameter :: c2es = c1es*rgas/rwat
  real(rk8) , parameter :: c3les = 17.2693882D+00
  real(rk8) , parameter :: c3ies = 21.875D+00
  real(rk8) , parameter :: c4les = 35.86D+00
  real(rk8) , parameter :: c4ies = 7.66D+00
  real(rk8) , parameter :: c5les = c3les*(tzero-c4les)
  real(rk8) , parameter :: c5ies = c3ies*(tzero-c4ies)
  real(rk8) , parameter :: c5alvcp = c5les*wlhv*rcpd
  real(rk8) , parameter :: c5alscp = c5ies*wlhs*rcpd
  real(rk8) , parameter :: wlhvocp = wlhv*rcpd
  real(rk8) , parameter :: wlhsocp = wlhs*rcpd
  real(rk8) , parameter :: wlhfocp = wlhf*rcpd
  real(rk8) , parameter :: cpowlhv = cpd*rwlhv
  real(rk8) , parameter :: cpowlhs = cpd*rwlhs
  real(rk8) , parameter :: cpowlhf = cpd*rwlhf
  real(rk8) , parameter :: rtber = tzero-5.0D0
  real(rk8) , parameter :: rtice = tzero-23.0D+00
  real(rk8) , parameter :: rtwat = tzero
  real(rk8) , parameter :: mpcrt = rtice+(rtwat-rtice)/sqrt(2.0D+00)
  real(rk8) , parameter :: rtwat_rtice_r = 1.0D+00/(rtwat-rtice)
  real(rk8) , parameter :: pq0 = 379.90516D+00
  ! value used for the latent heat term in the exponent for
  ! calculating equivalent potential temperature
  real(rk8) , parameter :: eliwv = 2.72D+06

  ! Standard atmosphere ICAO 1993
  real(rk8) , parameter :: stdp = 1.013250D+05
  real(rk8) , parameter :: stdpmb = 1013.250D+00
  real(rk8) , parameter :: stdt = 288.15D+00
  real(rk8) , parameter :: lrate = 0.00649D+00 ! K/m from MSL up to 11 km
  real(rk8) , parameter :: bltop = 0.960D+00

  ! Atmos. surface pressure mol/cm3
  real(rk8) , parameter :: atmos = 2.247D19
  ! Conversion parameter for Henry L-atm/mol-K
  real(rk8) , parameter :: rtcon = 8.314D-02
  ! RU g-cm2/s2-mol-K
  real(rk8) , parameter :: rumolec = 8.314D7
  ! Droplet diffusion coefficient Lelieveld+Crutzen 1991 cm2/sec
  real(rk8) , parameter :: dropdif = 2.0D-05
  ! Gas-phase diffusion coeff. Lelieveld and Crutzen, 1991 cm2/s
  real(rk8) , parameter :: difgas = 0.1D0

  ! Fixed emissivity of water
  real(rk8) , parameter :: emsw = 0.97D+00

  ! Trigonometric constants.
  real(rk8) , parameter :: mathpi =                                   &
                      &   3.1415926535897932384626433832795029D+00
  real(rk8) , parameter :: invpi = d_one/mathpi
  real(rk8) , parameter :: halfpi = mathpi*d_half
  real(rk8) , parameter :: twopi = mathpi*d_two
  real(rk8) , parameter :: degrad = mathpi/180.0D+00
  real(rk8) , parameter :: raddeg = 180.0D+00/mathpi

  ! Maximum stomatl resistance (s/m)
  real(rk8) , parameter :: rmax0 = 2.0D+04

  ! Maximum allowed dew(mm) and inverse (dewmaxi)
  real(rk8) , parameter :: dewmax = 0.1D+00
  real(rk8) , parameter :: dewmaxi = d_one/dewmax

  ! Maximum rate of transpiration with saturated soil (kg/m**2/s)
  real(rk8) , parameter :: trsmx0 = 2.D-04

  ! Drainage out of 10m layer bottom (mm/s)
  ! drain is set fairly large to prevent swamping the soil
  real(rk8) , parameter :: drain = 4.0D-04

  ! Minimum ratio between potential and actual water content
  real(rk8) , parameter :: minwrat = 1.0D-04

  ! Earth radius in meters
  real(rk8) , parameter :: earthrad = 6.371229D+06
  real(rk8) , parameter :: erkm = earthrad/d_1000
  real(rk8) , parameter :: rearthrad = d_one/earthrad
  ! Angular velocity of rotation of Earth
  real(rk8) , parameter :: eomeg = 7.2921159D-05
  real(rk8) , parameter :: eomeg2 = d_two*eomeg

  ! Soil roughness length
  real(rk8) , parameter :: zlnd = 0.01D+00
  ! Ocean roughness length
  real(rk8) , parameter :: zoce = 0.00023D+00
  ! Snow roughness length
  real(rk8) , parameter :: zsno = 0.00040D+00
  ! Von Karman constant
  real(rk8) , parameter :: vonkar = 0.4D+00
  ! Drag coefficient for soil under canopy
  real(rk8) , parameter :: csoilc = 4.0D-03
  ! Turbulent wind for stable conditions (m/sec)
  real(rk8) , parameter :: wtur = 0.1D+00

  ! Constant used in computing virtual temperature.
  real(rk8) , parameter :: ep1 = amd/amw - d_one ! 0.6077762876
  ! Constant used in computing saturation mixing ratio.
  ! Ratio of mean molecular weight of water to that of dry air
  real(rk8) , parameter :: ep2 = amw/amd   ! 0.6219770795
  real(rk8) , parameter :: rep2 = amd/amw  ! 1.6077762876

  ! Terminal velocity constants
  real(rk8) , parameter :: avt = 841.99667D+00
  real(rk8) , parameter :: bvt = 0.8D+00
  real(rk8) , parameter :: g4pb = 17.837825D+00
  real(rk8) , parameter :: g3pb = g4pb/(3.0D+00+bvt)
  real(rk8) , parameter :: g5pb = 1.8273D+00
  real(rk8) , parameter :: vtc = avt*g4pb/6.0D+00
  real(rk8) , parameter :: trel = 3000.0D+00

  ! Dynamic parameters
  ! alpha = .2495 in brown-campana; = 0. in split explicit
  real(rk8) , parameter :: alpha_hyd = 0.0D+00
  real(rk8) , parameter :: beta_hyd = d_one - d_two*alpha_hyd

  real(rk8) , parameter :: gnu = 0.10D+00
  real(rk8) , parameter :: omu = d_one - d_two*gnu
  real(rk8) , parameter :: gnuhf = d_half*gnu
  real(rk8) , parameter :: omuhf = d_one - d_two*gnuhf

  ! Constant surface Long Wave emissivity
  real(rk8) , parameter :: lnd_sfcemiss = 0.985D0
  real(rk8) , parameter :: ocn_sfcemiss = 0.984D0

  ! Constants used in Kain-Fritsch
  real(rk8) , parameter :: aliq = 613.3D0
  real(rk8) , parameter :: bliq = 17.502D0
  real(rk8) , parameter :: cliq = 4780.8D0
  real(rk8) , parameter :: dliq = 32.19D0
  real(rk8) , parameter :: aice = 613.20D0
  real(rk8) , parameter :: bice = 22.452D0
  real(rk8) , parameter :: cice = 6133.0D0
  real(rk8) , parameter :: dice = 0.61D0
  real(rk8) , parameter :: xlv0 = 3.15D6
  real(rk8) , parameter :: xlv1 = 2370.0D0
  real(rk8) , parameter :: xls0 = 2.905D+06
  real(rk8) , parameter :: xls1 = 259.532D+00

  ! GTS system constants
  real(rk8) , parameter :: egravgts = egrav*d_100
  real(rk8) , parameter :: regravgts = d_one/egravgts
  real(rk8) , parameter :: cpdgts = cpd*1.0D+04
  real(rk8) , parameter :: gocp = egravgts/cpdgts
  real(rk8) , parameter :: sslp = stdp*d_10 ! dynes/cm^2
  real(rk8) , parameter :: rsslp = d_one/sslp
  real(rk8) , parameter :: stebol = sigm*d_1000
  real(rk8) , parameter :: rgsslp = d_half/(egravgts*sslp)
  ! Effective molecular weight of dry air (kg/mol)
  real(rk8) , parameter :: amdk = amd*d_r1000
  ! Avogadro Constant in lit/cm3
  real(rk8) , parameter :: avogadrl = navgdr*d_1000

  ! Radiation constants
  real(rk8) , parameter :: dpfco2 = 5.0D-03
  real(rk8) , parameter :: dpfo3 = 2.5D-03

  ! Pressure gradient force calculations (Why not standard atmosphere?)
  real(rk8) , parameter :: t00pg = 287.0D+00       ! stdt ?
  real(rk8) , parameter :: p00pg = 101.325D+00     ! stdp ?
  real(rk8) , parameter :: alam  = 6.5D-03         ! Lapse rate ?
  real(rk8) , parameter :: pgfaa1 = alam*rgas*regrav ! Utility constant

  ! Molecular heat diffusion coefficient in water
  real(rk8) , parameter :: hdmw = 1.3889D-07  ! m^2/s

  ! Seaice temperature from ICBC trigger value
  real(rk8) , parameter :: icetemp = 271.38D0
  real(rk8) , parameter :: iceminh = 0.01D0

  ! Allowed range for cloud fraction
  real(rk8) , parameter :: lowcld = 0.0001D0
  real(rk8) , parameter :: hicld  = 0.9999D0

  contains

  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
  pure real(rk8) function pfesat(t) result(es)
    implicit none
    real(rk8) , intent(in) :: t     ! Temperature (K)

    real(rk8) :: td , t_limit
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  6.11213476D0
    real(rk8) , parameter :: a1 =  0.444007856D0
    real(rk8) , parameter :: a2 =  0.143064234D-01
    real(rk8) , parameter :: a3 =  0.264461437D-03
    real(rk8) , parameter :: a4 =  0.305903558D-05
    real(rk8) , parameter :: a5 =  0.196237241D-07
    real(rk8) , parameter :: a6 =  0.892344772D-10
    real(rk8) , parameter :: a7 = -0.373208410D-12
    real(rk8) , parameter :: a8 =  0.209339997D-15
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  6.11123516D0
    real(rk8) , parameter :: c1 =  0.503109514D0
    real(rk8) , parameter :: c2 =  0.188369801D-01
    real(rk8) , parameter :: c3 =  0.420547422D-03
    real(rk8) , parameter :: c4 =  0.614396778D-05
    real(rk8) , parameter :: c5 =  0.602780717D-07
    real(rk8) , parameter :: c6 =  0.387940929D-09
    real(rk8) , parameter :: c7 =  0.149436277D-11
    real(rk8) , parameter :: c8 =  0.262655803D-14

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      es = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      es = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es   = es * d_100 ! pa
  end function pfesat
  !
  ! Computes derivative in temperature of saturation pressure
  !
  pure real(rk8) function pfesdt(t) result(esdt)
    implicit none
    real(rk8), intent(in)  :: t     ! Temperature (K)

    real(rk8) :: td , t_limit
    !
    ! For derivative:water vapor
    !
    real(rk8), parameter :: b0 =  0.444017302D0
    real(rk8), parameter :: b1 =  0.286064092D-01
    real(rk8), parameter :: b2 =  0.794683137D-03
    real(rk8), parameter :: b3 =  0.121211669D-04
    real(rk8), parameter :: b4 =  0.103354611D-06
    real(rk8), parameter :: b5 =  0.404125005D-09
    real(rk8), parameter :: b6 = -0.788037859D-12
    real(rk8), parameter :: b7 = -0.114596802D-13
    real(rk8), parameter :: b8 =  0.381294516D-16
    !
    ! For derivative:ice
    !
    real(rk8), parameter :: d0 =  0.503277922D0
    real(rk8), parameter :: d1 =  0.377289173D-01
    real(rk8), parameter :: d2 =  0.126801703D-02
    real(rk8), parameter :: d3 =  0.249468427D-04
    real(rk8), parameter :: d4 =  0.313703411D-06
    real(rk8), parameter :: d5 =  0.257180651D-08
    real(rk8), parameter :: d6 =  0.133268878D-10
    real(rk8), parameter :: d7 =  0.394116744D-13
    real(rk8), parameter :: d8 =  0.498070196D-16

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      esdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      esdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if
    esdt = esdt * 100.D0 ! pa/K
  end function pfesdt

  pure real(rk8) function pfqsat(t,p,e) result(qs)
    implicit none
    real(rk8) , intent(in) :: t             ! Temperature (K)
    real(rk8) , intent(in) :: p             ! Pressure (Pa)
    real(rk8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rk8) :: es , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ! Bolton 1980
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qs = max(es * vp1, minqq)  ! kg/kg
  end function pfqsat

  pure real(rk8) function pfqsdt(t,p,e,dedt) result(qsdt)
    implicit none
    real(rk8) , intent(in) :: t             ! Temperature (K)
    real(rk8) , intent(in) :: p             ! Pressure (Pa)
    real(rk8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rk8) , intent(in) , optional :: dedt ! derivative of e in dt (Pa/K)
    real(rk8) :: es , esdt , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    if ( present(dedt) ) then
      esdt = dedt
    else
      esdt = pfesdt(t)
    end if
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qsdt = esdt * vp2 * p ! 1 / K
  end function pfqsdt

  pure real(rk8) function sig2p(ps,sigma,ptop) result(p)
    implicit none
    real(rk8) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
    p = (sigma * ( ps - ptop ) + ptop) * d_1000 ! Pressure in Pa
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
  end function sig2p

  pure real(rk8) function msig2p(ps,sigma,ptop) result(p)
    implicit none
    real(rk8) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!!!!! Assume input is in cbar !!!!!!!!!!!!!
    p = (sigma * ps + ptop) * d_1000 ! Pressure in Pa !
    !!!!!!!!! In the model, ps is ps-ptop !!!!!!!!!!!!!
    !!!!!!!!!!!!! Assume input is in cbar !!!!!!!!!!!!!
  end function msig2p

  ! Latent heat for condensation/sublimation of water in J/kg
  pure real(rk8) function wlh(t)
    implicit none
    real(rk8) , intent(in) :: t
    real(rk8) :: tc
    tc = t - tzero
    if ( tc > d_zero ) then
      wlh = 2500.79D0 - 2.36418D0*tc + 1.58927D-3*tc*tc - 6.14342D-5*tc*tc*tc
    else
      wlh = 2834.1D0 - 0.29D0*tc - 0.004D0*tc*tc
    end if
    wlh = wlh*1.0D3
  end function wlh

end module mod_constants
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
