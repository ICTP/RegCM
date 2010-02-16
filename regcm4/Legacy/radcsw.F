C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine radcsw(pint    ,h2ommr  ,o3mmr   ,
     $                  cld     ,clwp    ,rel     ,rei     ,fice    ,
     $                  eccf             ,asdir   ,asdif   ,aldir   ,
     $                  aldif   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $                  fsds    ,fsnsc   ,fsntc   ,sols    ,soll    ,
     $                  solsd   ,solld   ,fsnirt  ,fsnrtc  ,fsnirtsq,
     $                  tauxar_mix, tauasc_mix, gtota_mix,ftota_mix,
     $      tauxar_mix_cs, tauasc_mix_cs, gtota_mix_cs,ftota_mix_cs,
     $                   aeradfo,aeradfos)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cqian 30/06/99,  csm new scheme: hygroscopic growth effect of
cqian            sulfate has been included
cqian            main changed codes: radcsw,

C-----------------------------------------------------------------------
C
C Solar radiation code
C
C Basic method is Delta-Eddington as described in:
C
C    Briegleb, Bruce P., 1992: Delta-Eddington
C    Approximation for Solar Radiation in the NCAR Community Climate Model,
C    Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C Two changes to the basic method described above are: (1) the distinction
C between liquid and ice particle clouds, and (2) the addition of an
C aerosol with sulfate radiative properties.
C
C Divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
C solar flux fractions specified for each interval. allows for
C seasonally and diurnally varying solar input.  Includes molecular,
C cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud,
C and surface absorption. Computes delta-eddington reflections and
C transmissions assuming homogeneously mixed layers. Adds the layers
C assuming scattering between layers to be isotropic, and distinguishes
C direct solar beam from scattered radiation.
C
C Longitude loops are broken into 1 or 2 sections, so that only daylight
C (i.e. coszrs > 0) computations are done.
C
C Note that an extra layer above the model top layer is added.
C
C cgs units are used.
C
C Special diagnostic calculation of the clear sky surface and total column
C absorbed flux is also done for cloud forcing diagnostics.
C
C
C---------------------------Code history--------------------------------
C
C Modified March 1995 to add aerosols
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C 19th Band Added:   W. Collins March 1997
C Merge Wong optics: T. Schneider, Mar 1998
C
C-----------------------------------------------------------------------
c
c $Id: radcsw.F,v 1.1.2.1.4.2 1998/07/10 20:51:26 eaton Exp $
c
      implicit none
CEES
C------------------------------Parameters-------------------------------
#     include "parameter.inc"
#     include "include/parrad.cb"
#     include "include/param2.cb"
#     include "include/bats2.cb"
#     include "include/bats.cb"
#     include "include/tracer.cb"
#     include "include/aeroppt.cb"

c
C-----------------------------------------------------------------------
      real(kind=8) scon                ! Solar constant (in erg/cm**2/sec)
      parameter (scon = 1.367e6)
CEES added scon here: CSM vers. passes scon into radcsw

C-----------------------------------------------------------------------
      integer nspint           ! Num of spctrl intervals across solar spectrum

      parameter ( nspint = 19 )
C-----------------------Constants for new bands--------------------------
      real(kind=8)
     $     V_RAYTAU_35,
     $     V_RAYTAU_64,
     $     V_ABO3_35,
     $     V_ABO3_64
      parameter(
     $     V_RAYTAU_35 = 0.155208,
     $     V_RAYTAU_64 = 0.0392,
     $     V_ABO3_35 = 2.4058030e+01,
     $     V_ABO3_64 = 2.210e+01
     $     )
C------------------------------Commons----------------------------------
c
C Radiation constants
C
#     include "include/crdcon.cb"
C
C Input arguments
C
      real(kind=8)  pint(plond,plevp)   ! Interface pressure
      real(kind=8)  h2ommr(plond,plev)  ! Specific humidity (h2o mass mix ratio)
      real(kind=8)  o3mmr(plond,plev)   ! Ozone mass mixing ratio
C
      real(kind=8)  cld(plond,plevp)    ! Fractional cloud cover
      real(kind=8)  clwp(plond,plev)    ! Layer liquid water path
      real(kind=8)  rel(plond,plev)     ! Liquid effective drop size (microns)
      real(kind=8)  rei(plond,plev)     ! Ice effective drop size (microns)
      real(kind=8)  fice(plond,plev)    ! Fractional ice content within cloud
C
      real(kind=8)  eccf                ! Eccentricity factor (1./earth-sun dist ** 2)
CEES  real(kind=8)  coszrs(plond)       ! Cosine solar zenith angle
      real(kind=8)  asdir(plond)        ! 0.2-0.7 micro-meter srfc alb to direct rad
      real(kind=8)  aldir(plond)        ! 0.7-5.0 micro-meter srfc alb to direct rad
      real(kind=8)  asdif(plond)        ! 0.2-0.7 micro-meter srfc alb to diffuse rad
      real(kind=8)  aldif(plond)        ! 0.7-5.0 micro-meter srfc alb to diffuse rad

CEES parameter set in radcsw      real scon                ! solar constant
C
C Output arguments
C
CEES these differ, need to check
      real(kind=8)  solin(plond)    ! Incident solar flux
      real(kind=8)  qrs(plond,plev) ! Solar heating rate
      real(kind=8)  fsns(plond)     ! Surface absorbed solar flux
      real(kind=8)  fsnt(plond)     ! Total column absorbed solar flux
      real(kind=8)  fsds(plond)     ! Flux Shortwave Downwelling Surface
C
      real(kind=8)  fsnsc(plond)    ! Clear sky surface absorbed solar flux
      real(kind=8)  fsntc(plond)    ! Clear sky total column absorbed solar flx
      real(kind=8)  sols(plond)     ! Direct solar rad incident on surface (< 0.7)
      real(kind=8)  soll(plond)     ! Direct solar rad incident on surface (>= 0.7)
      real(kind=8)  solsd(plond)    ! Diffuse solar rad incident on surface (< 0.7)
      real(kind=8)  solld(plond)    ! Diffuse solar rad incident on surface (>= 0.7)
      real(kind=8)  fsnirt(plond)   ! Near-IR flux absorbed at toa
      real(kind=8)  fsnrtc(plond)   ! Clear sky near-IR flux absorbed at toa
      real(kind=8)  fsnirtsq(plond) ! Near-IR flux absorbed at toa >= 0.7 microns
C
C------------------------------Externals--------------------------------
C
      integer   isrchfgt      ! Search for first array element > 0
      integer   isrchfle      ! Search for first array element < 0
      external  isrchfgt,isrchfle
C
C---------------------------Local variables-----------------------------
C
      integer ns              ! Spectral loop index
      integer i               ! Longitude loop index
      integer k               ! Level loop index
      integer n               ! Loop index for daylight
      integer nloop           ! Number of daylight loops
      integer is(2)           ! Daytime start indices
      integer ie(2)           ! Daytime end indices
      integer indxsl          ! Index for cloud particle properties
C
C A. Slingo's data for cloud particle radiative properties (from 'A GCM
C Parameterization for the Shortwave Properties of Water Clouds' JAS
C vol. 46 may 1989 pp 1419-1427)
C
      real(kind=8)  abarl(4)        ! A coefficient for extinction optical depth
      real(kind=8)  bbarl(4)        ! B coefficient for extinction optical depth
      real(kind=8)  cbarl(4)        ! C coefficient for single particle scat albedo
      real(kind=8)  dbarl(4)        ! D coefficient for single particle scat albedo
      real(kind=8)  ebarl(4)        ! E coefficient for asymmetry parameter
      real(kind=8)  fbarl(4)        ! F coefficient for asymmetry parameter

      save abarl, bbarl, cbarl, dbarl, ebarl, fbarl

      data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
      data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
      data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
      data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
      data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
      data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/

      real(kind=8)  abarli          ! A coefficient for current spectral interval
      real(kind=8)  bbarli          ! B coefficient for current spectral interval
      real(kind=8)  cbarli          ! C coefficient for current spectral interval
      real(kind=8)  dbarli          ! D coefficient for current spectral interval
      real(kind=8)  ebarli          ! E coefficient for current spectral interval
      real(kind=8)  fbarli          ! F coefficient for current spectral interval
C
C Caution... A. Slingo recommends no less than 4.0 micro-meters nor
C greater than 20 micro-meters
C
c ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
c
      real(kind=8)  abari(4)        ! a coefficient for extinction optical depth
      real(kind=8)  bbari(4)        ! b coefficient for extinction optical depth
      real(kind=8)  cbari(4)        ! c coefficient for single particle scat albedo
      real(kind=8)  dbari(4)        ! d coefficient for single particle scat albedo
      real(kind=8)  ebari(4)        ! e coefficient for asymmetry parameter
      real(kind=8)  fbari(4)        ! f coefficient for asymmetry parameter

      save abari, bbari, cbari, dbari, ebari, fbari

      data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
      data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
      data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
      data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
      data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
      data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/

      real(kind=8)  abarii          ! A coefficient for current spectral interval
      real(kind=8)  bbarii          ! B coefficient for current spectral interval
      real(kind=8)  cbarii          ! C coefficient for current spectral interval
      real(kind=8)  dbarii          ! D coefficient for current spectral interval
      real(kind=8)  ebarii          ! E coefficient for current spectral interval
      real(kind=8)  fbarii          ! F coefficient for current spectral interval
C
      real(kind=8)  delta           ! Pressure (atmospheres) for stratos. h2o limit
      real(kind=8)  o2mmr           ! O2 mass mixing ratio:

      save delta, o2mmr

C CO2 info:
C
c      real(kind=8)  mmwair,         ! Mean molecular weight of air
c    $        mmwco2          ! Mean molecular weight of co2
c     save    mmwair, mmwco2
c     data    mmwair / 28.9644 /
c     data    mmwco2 / 44.0000 /

      data    delta /  1.70e-3 /
      data    o2mmr / .23143 /
C
C Next series depends on spectral interval
C
      real(kind=8)  frcsol(nspint) ! Fraction of solar flux in each spectral interval
      real(kind=8)  wavmin(nspint) ! Min wavelength (micro-meters) of interval
      real(kind=8)  wavmax(nspint) ! Max wavelength (micro-meters) of interval
      real(kind=8)  raytau(nspint) ! Rayleigh scattering optical depth
      real(kind=8)  abh2o(nspint)  ! Absorption coefficiant for h2o (cm2/g)
      real(kind=8)  abo3 (nspint)  ! Absorption coefficiant for o3  (cm2/g)
      real(kind=8)  abco2(nspint)  ! Absorption coefficiant for co2 (cm2/g)
      real(kind=8)  abo2 (nspint)  ! Absorption coefficiant for o2  (cm2/g)
      real(kind=8)  ph2o(nspint)   ! Weight of h2o in spectral interval
      real(kind=8)  pco2(nspint)   ! Weight of co2 in spectral interval
      real(kind=8)  po2 (nspint)   ! Weight of o2  in spectral interval
      real(kind=8)  nirwgt(nspint) ! Weight for intervals to simulate satellite filter
      real(kind=8)  wgtint         ! Weight for specific spectral interval

      save frcsol ,wavmin ,wavmax ,raytau ,abh2o ,abo3 ,
     $     abco2  ,abo2   ,ph2o   ,pco2   ,po2   ,nirwgt

      data frcsol / .001488, .001389, .001290, .001686, .002877,
     $              .003869, .026336, .360739, .065392, .526861,
     $              .526861, .526861, .526861, .526861, .526861,
     $              .526861, .006239, .001834, .001834/
C
C weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
C
      data nirwgt /  0.0,   0.0,   0.0,      0.0,   0.0,
     $               0.0,   0.0,   0.0, 0.320518,   1.0,  1.0,
     $               1.0,   1.0,   1.0,      1.0,   1.0,
     $               1.0,   1.0,   1.0 /

      data wavmin / .200,  .245,  .265,  .275,  .285,
     $              .295,  .305,  .350,  .640,  .700,  .701,
     $              .701,  .701,  .701,  .702,  .702,
     $             2.630, 4.160, 4.160/

      data wavmax / .245,  .265,  .275,  .285,  .295,
     $              .305,  .350,  .640,  .700, 5.000, 5.000,
     $             5.000, 5.000, 5.000, 5.000, 5.000,
     $             2.860, 4.550, 4.550/

      data raytau / 4.020, 2.180, 1.700, 1.450, 1.250,
     $              1.085, 0.730, V_RAYTAU_35, V_RAYTAU_64, 0.020,
     $              .0001, .0001, .0001, .0001, .0001, .0001,
     $              .0001, .0001, .0001/
C
C Absorption coefficients
C
      data abh2o /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .002,
     $                .035,     .377,    1.950,   9.400,  44.600,
     $             190.000,     .000,     .000,    .000/

      data abo3  /
     $ 5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04,
     $ 4.441e+03,  1.775e+02, V_ABO3_35, V_ABO3_64,      .000,
     $      .000,   .000    ,   .000   ,   .000   ,      .000,
     $      .000,   .000    ,   .000   ,   .000    /

      data abco2  /    .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,     .000,    .000,    .000,
     $                 .000,     .094,     .196,   1.963/

      data abo2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,1.11e-05,6.69e-05,
     $                .000,     .000,     .000,    .000,    .000,
     $                .000,     .000,    .000,    .000/
C
C Spectral interval weights
C
      data ph2o  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .505,
     $                .210,     .120,     .070,    .048,    .029,
     $                .018,     .000,     .000,    .000/

      data pco2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,     .000,    .000,    .000,
     $                .000,    1.000,     .640,    .360/

      data po2   /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,     .000,   1.000,   1.000,
     $                .000,     .000,     .000,    .000,    .000,
     $                .000,     .000,     .000,    .000/
C
C Diagnostic and accumulation arrays; note that sfltot, fswup, and
C fswdn are not used in the computation,but are retained for future use.
C
      real(kind=8)  solflx(plond)        ! Solar flux in current interval
      real(kind=8)  sfltot(plond)        ! Spectrally summed total solar flux
      real(kind=8)  totfld(plond,0:plev) ! Spectrally summed flux divergence
      real(kind=8)  fswup(plond,0:plevp) ! Spectrally summed up flux
      real(kind=8)  fswdn(plond,0:plevp) ! Spectrally summed down flux
C
C Cloud radiative property arrays
C
      real(kind=8)  tauxcl(plond,0:plev) ! water cloud extinction optical depth
      real(kind=8)  tauxci(plond,0:plev) ! ice cloud extinction optical depth
      real(kind=8)  wcl(plond,0:plev)    ! liquid cloud single scattering albedo
      real(kind=8)  gcl(plond,0:plev)    ! liquid cloud asymmetry parameter
      real(kind=8)  fcl(plond,0:plev)    ! liquid cloud forward scattered fraction
      real(kind=8)  wci(plond,0:plev)    ! ice cloud single scattering albedo
      real(kind=8)  gci(plond,0:plev)    ! ice cloud asymmetry parameter
      real(kind=8)  fci(plond,0:plev)    ! ice cloud forward scattered fraction
C 
C Various arrays and other constants:
C
      real(kind=8)  pflx(plond,0:plevp)  ! Interface press, including extra layer
      real(kind=8)  zenfac(plond)        ! Square root of cos solar zenith angle
      real(kind=8)  sqrco2               ! Square root of the co2 mass mixg ratio
      real(kind=8)  tmp1                 ! Temporary constant array
      real(kind=8)  tmp2                 ! Temporary constant array
      real(kind=8)  pdel                 ! Pressure difference across layer
      real(kind=8)  path                 ! Mass path of layer
      real(kind=8)  ptop                 ! Lower interface pressure of extra layer
      real(kind=8)  ptho2                ! Used to compute mass path of o2
      real(kind=8)  ptho3                ! Used to compute mass path of o3
      real(kind=8)  pthco2               ! Used to compute mass path of co2
      real(kind=8)  pthh2o               ! Used to compute mass path of h2o
      real(kind=8)  h2ostr               ! Inverse square root h2o mass mixing ratio
      real(kind=8)  wavmid               ! Spectral interval middle wavelength
      real(kind=8)  trayoslp             ! Rayleigh optical depth/standard pressure
      real(kind=8)  tmp1l                ! Temporary constant array
      real(kind=8)  tmp2l                ! Temporary constant array
      real(kind=8)  tmp3l                ! Temporary constant array
      real(kind=8)  tmp1i                ! Temporary constant array
      real(kind=8)  tmp2i                ! Temporary constant array
      real(kind=8)  tmp3i                ! Temporary constant array
      real(kind=8)  rdenom               ! Multiple scattering term
      real(kind=8)  psf                  ! Frac of solar flux in spect interval
      real(kind=8)  gocp                 ! Gravity/cp
C
C Layer absorber amounts; note that 0 refers to the extra layer added
C above the top model layer
C
      real(kind=8)  uh2o(plond,0:plev)   ! Layer absorber amount of h2o
      real(kind=8)  uo3(plond,0:plev)    ! Layer absorber amount of  o3
      real(kind=8)  uco2(plond,0:plev)   ! Layer absorber amount of co2
      real(kind=8)  uo2(plond,0:plev)    ! Layer absorber amount of  o2
      real(kind=8)  uaer(plond,0:plev)   ! Layer aerosol amount
C
C Total column absorber amounts:
C
      real(kind=8)  uth2o(plond)         ! Total column  absorber amount of h2o
      real(kind=8)  uto3(plond)          ! Total column  absorber amount of  o3
      real(kind=8)  utco2(plond)         ! Total column  absorber amount of co2
      real(kind=8)  uto2(plond)          ! Total column  absorber amount of  o2
C
C These arrays are defined for plev model layers; 0 refers to the extra
C layer on top:
C
      real(kind=8)  rdir(plond,0:plev)    ! Layer reflectivity to direct rad
      real(kind=8)  rdif(plond,0:plev)    ! Layer reflectivity to diffuse rad
      real(kind=8)  tdir(plond,0:plev)    ! Layer transmission to direct rad
      real(kind=8)  tdif(plond,0:plev)    ! Layer transmission to diffuse rad
      real(kind=8)  explay(plond,0:plev)  ! Solar beam exp transmission for layer
      real(kind=8)  flxdiv(plond,0:plev)  ! Flux divergence for layer
C
C These arrays are defined at model interfaces; 0 is the top of the
C extra layer above the model top; plevp is the earth surface:
C
      real(kind=8)  rupdir(plond,0:plevp) ! Ref to dir rad for layers below
      real(kind=8)  rupdif(plond,0:plevp) ! Ref to dif rad for layers below
      real(kind=8)  rdndif(plond,0:plevp) ! Ref to dif rad for layers above
      real(kind=8)  exptdn(plond,0:plevp) ! Solar beam exp down transm from top
      real(kind=8)  tottrn(plond,0:plevp) ! Total transmission for layers above
      real(kind=8)  fluxup(plond,0:plevp) ! Up   flux at model interface
      real(kind=8)  fluxdn(plond,0:plevp) ! Down flux at model interface
C

      real(kind=8)  tauxar_mix(plond,0:plev,nspi) ! aerosol extinction optical depth
      real(kind=8)  tauasc_mix(plond,0:plev,nspi)
      real(kind=8)  gtota_mix(plond,0:plev,nspi)
      real(kind=8)  ftota_mix(plond,0:plev,nspi)
      real(kind=8)  tauxar_mix_cs(plond,nspi)  ! aerosol extinction optical depth
      real(kind=8)  tauasc_mix_cs(plond,nspi)
      real(kind=8)  gtota_mix_cs(plond,nspi)
      real(kind=8)  ftota_mix_cs(plond,nspi)

      real(kind=8)  wkaer(plond,0:plev,4) !works table 

c variable for radiative forcing calculation 
  
      real(kind=8)  zero(plond,4) 
      real(kind=8)  x0fsntc(plond), x0fsnsc(plond),x0fsnrtc(plond)  
      real(kind=8)  aeradfo(plond),aeradfos(plond)  ! spectrally integrated aerosol radiative forcing ( TOA) 
C-----------------------------------------------------------------------
C
          
C Initialize output fields:
C
      do i=1, plon
        fsds(i)  = 0.0
        fsnirt(i) = 0.0
        fsnrtc(i) = 0.0
        fsnirtsq(i) = 0.0
        fsnt(i)  = 0.0
        fsns(i)  = 0.0
        solin(i) = 0.0
        fsnsc(i) = 0.0
        fsntc(i) = 0.0
        sols(i) = 0.0
        soll(i) = 0.0
        solsd(i) = 0.0
        solld(i) = 0.0
        sabveg(i) = 0.0
        solis(i) = 0.0
        solvs(i) = 0.0
        solvd(i) = 0.0
c
        aeradfo(i)=0.0
        aeradfos(i) =0.0
        x0fsntc(i)=0.0
        x0fsnsc(i)=0.0
        x0fsnrtc(i)=0.0
      end do
      do k=1, plev
        do i=1, plon
          qrs(i,k) = 0.0
        end do
      end do
      do k=1, plev
        do i=1, plon
          pdel   = pint(i,k+1) - pint(i,k)
          path   = pdel / gravit
        end do
      end do
C
C Compute starting, ending daytime loop indices:
C
      nloop = 0
      is(1) = isrchfgt(plon,coszrs,1,0.d0)
C
C If night everywhere, return:
C
      if(is(1).gt.plon) return
      ie(1) = isrchfle(plon-is(1),coszrs(is(1)+1),1,0.d0) + is(1) - 1
      nloop = 1
C
C Possibly 2 daytime loops needed:
C
      if (ie(1).ne.plon) then
        is(2) = isrchfgt(plon-ie(1),coszrs(ie(1)+1),1,0.d0) + ie(1)
        if(is(2).le.plon) then
          nloop = 2
          ie(2) = plon
        end if
      end if
C
C Define solar incident radiation and interface pressures:
C
      do n=1,nloop
        do i=is(n),ie(n)
          solin(i) = scon*eccf*coszrs(i)
          pflx(i,0) = 0.
        end do
      end do
      do k=1,plevp
        do n=1,nloop
          do i=is(n),ie(n)
            pflx(i,k) = pint(i,k)
          end do
        end do
      end do
C
C Compute optical paths:
c  CO2, use old scheme(as constant)
C
      tmp1   = 0.5/(gravit*sslp)
c     co2mmr = co2vmr*(mmwco2/mmwair)

      sqrco2 = sqrt(co2mmr)
      do n=1,nloop
        do i=is(n),ie(n)
          ptop      = pflx(i,1)
          ptho2     = o2mmr * ptop / gravit
          ptho3     = o3mmr(i,1) * ptop / gravit
          pthco2    = sqrco2 * (ptop / gravit)
          h2ostr    = sqrt( 1. / h2ommr(i,1) )
          zenfac(i) = sqrt(coszrs(i))
          pthh2o    = ptop**2*tmp1 + (ptop*rga)*(h2ostr*zenfac(i)*delta)
          uh2o(i,0) = h2ommr(i,1)*pthh2o
          uco2(i,0) = zenfac(i)*pthco2
          uo2 (i,0) = zenfac(i)*ptho2
          uo3 (i,0) = ptho3
          uaer(i,0) = 0.0
        end do
      end do
C
      tmp2 = delta/gravit
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            pdel   = pflx(i,k+1) - pflx(i,k)
            path   = pdel / gravit
            ptho2  = o2mmr * path
            ptho3  = o3mmr(i,k) * path
            pthco2 = sqrco2 * path
            h2ostr = sqrt(1.0/h2ommr(i,k))
            pthh2o = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 +
     $              pdel*h2ostr*zenfac(i)*tmp2
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2 (i,k) = zenfac(i)*ptho2
            uo3 (i,k) = ptho3
          end do
        end do
      end do
C
C Compute column absorber amounts for the clear sky computation:
C
      do n=1,nloop
        do i=is(n),ie(n)
          uth2o(i) = 0.0
          uto3(i)  = 0.0
          utco2(i) = 0.0
          uto2(i)  = 0.0
        end do
      end do
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            uth2o(i) = uth2o(i) + uh2o(i,k)
            uto3(i)  = uto3(i)  + uo3(i,k)
            utco2(i) = utco2(i) + uco2(i,k)
            uto2(i)  = uto2(i)  + uo2(i,k)
          end do
        end do
      end do
C
C Initialize spectrally integrated totals:
C
      do k=0,plev
        do i=1,plon
          totfld(i,k) = 0.0
          fswup (i,k) = 0.0
          fswdn (i,k) = 0.0
        end do
      end do
      do i=1,plon
        sfltot(i)       = 0.0
        fswup (i,plevp) = 0.0
        fswdn (i,plevp) = 0.0
      end do
C
C Set cloud properties for top (0) layer; so long as tauxcl is zero,
C there is no cloud above top of model; the other cloud properties
C are arbitrary:
C
      do n=1,nloop
        do i=is(n),ie(n)
          tauxcl(i,0) = 0.
          wcl(i,0)     = 0.999999
          gcl(i,0)     = 0.85
          fcl(i,0)     = 0.725
          tauxci(i,0) = 0.
          wci(i,0)     = 0.999999
          gci(i,0)     = 0.85
          fci(i,0)     = 0.725
        end do
      end do
C
C Begin spectral loop
C
      do 100 ns=1,nspint
        wgtint = nirwgt(ns)
C
C Set index for cloud particle properties based on the wavelength,
C according to A. Slingo (1989) equations 1-3:
C Use index 1 (0.25 to 0.69 micrometers) for visible
C Use index 2 (0.69 - 1.19 micrometers) for near-infrared
C Use index 3 (1.19 to 2.38 micrometers) for near-infrared
C Use index 4 (2.38 to 4.00 micrometers) for near-infrared
C
C Note that the minimum wavelength is encoded (with .001, .002, .003)
C in order to specify the index appropriate for the near-infrared
C cloud absorption properties
C
        if(wavmax(ns) .le. 0.7) then
          indxsl = 1
        else if(wavmin(ns) .eq. 0.700) then
          indxsl = 2
        else if(wavmin(ns) .eq. 0.701) then
          indxsl = 3
        else if(wavmin(ns) .eq. 0.702 .or. wavmin(ns) .gt. 2.38) then
          indxsl = 4
        end if
C
C Set cloud extinction optical depth, single scatter albedo,
C asymmetry parameter, and forward scattered fraction:
C
        abarli = abarl(indxsl)
        bbarli = bbarl(indxsl)
        cbarli = cbarl(indxsl)
        dbarli = dbarl(indxsl)
        ebarli = ebarl(indxsl)
        fbarli = fbarl(indxsl)
c
        abarii = abari(indxsl)
        bbarii = bbari(indxsl)
        cbarii = cbari(indxsl)
        dbarii = dbari(indxsl)
        ebarii = ebari(indxsl)
        fbarii = fbari(indxsl)
c
        do k=1,plev
          do n=1,nloop
            do i=is(n),ie(n)
c
c liquid
c
              tmp1l = abarli + bbarli/rel(i,k)
              tmp2l = 1. - cbarli - dbarli*rel(i,k)
              tmp3l = fbarli*rel(i,k)
c
c ice
c
              tmp1i = abarii + bbarii/rei(i,k)
              tmp2i = 1. - cbarii - dbarii*rei(i,k)
              tmp3i = fbarii*rei(i,k)
C
C Cloud fraction incorporated into cloud extinction optical depth
C
Cfound
c April 12 2000, Filippo found the different scheme here:

Cscheme 1
Cccm3.6.6
c             tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
c    $                     *cld(i,k)*sqrt(cld(i,k))
c             tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
c    $                     *cld(i,k)*sqrt(cld(i,k))
c
Cscheme 2
cKN
              tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
     $   *cld(i,k)/(1.+(1.-0.85)*(1.-cld(i,k))*clwp(i,k)*tmp1l
     $                                        *(1.-fice(i,k)))
              tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
     $   *cld(i,k)/(1.+(1.-0.78)*(1.-cld(i,k))*clwp(i,k)*tmp1i
     $                                        *fice(i,k))

Cscheme 3
CEES below replaced
c             tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
c    $                     *cld(i,k)**0.85
c             tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
c    $                     *cld(i,k)**0.85
Cfound_

C
C Do not let single scatter albedo be 1; delta-eddington solution
C for non-conservative case:
C
cqian   30/06/99        wcl(i,k) = dmin1(tmp2l,.999999)
              wcl(i,k) = dmin1(tmp2l,.999999d0)
              gcl(i,k) = ebarli + tmp3l
              fcl(i,k) = gcl(i,k)*gcl(i,k)
C
              wci(i,k) = dmin1(tmp2i,.999999d0)
              gci(i,k) = ebarii + tmp3i
              fci(i,k) = gci(i,k)*gci(i,k)
C             
            end do
          end do
        end do
c
C Set reflectivities for surface based on mid-point wavelength
C
        wavmid = 0.5*(wavmin(ns) + wavmax(ns))
C
C Wavelength less  than 0.7 micro-meter
C
        if (wavmid .lt. 0.7 ) then
          do n=1,nloop
            do i=is(n),ie(n)
              albdir(i) = asdir(i)
              albdif(i) = asdif(i)
            end do
          end do
C
C Wavelength greater than 0.7 micro-meter
C
        else
          do n=1,nloop
            do i=is(n),ie(n)
              albdir(i) = aldir(i)
              albdif(i) = aldif(i)
            end do
          end do
        end if
        trayoslp = raytau(ns)/sslp
C
C Layer input properties now completely specified; compute the
C delta-Eddington solution reflectivities and transmissivities
C for each layer, starting from the top and working downwards:

c options for aerosol: no climatic feedback if idirect .eq. 1 
        if (idirect .eq. 2) then 
          do k=0,plev
          do i=1,plond
            wkaer(i,k,1) = tauxar_mix(i,k,ns)
            wkaer(i,k,2) = tauasc_mix(i,k,ns)
            wkaer(i,k,3) = gtota_mix(i,k,ns)
            wkaer(i,k,4) = ftota_mix(i,k,ns)
          end do
          end do
        else if (idirect .eq.1) then
          do k=0,plev
          do i=1,plond
            wkaer(i,k,1) = 0.
            wkaer(i,k,2) = 0.
            wkaer(i,k,3) = 0.
            wkaer(i,k,4) = 0.
          end do
          end do
        end if
            
        call radded(coszrs   ,trayoslp,pflx   ,abh2o(ns),abo3(ns),
     $              abco2(ns),abo2(ns),uh2o   ,uo3      ,uco2    ,
     $              uo2      ,tauxcl  ,wcl    ,gcl      ,fcl     ,
     $              tauxci   ,wci     ,gci    ,fci      ,
     &              wkaer(1,0,1),wkaer(1,0,2),wkaer(1,0,3),wkaer(1,0,4),
     &              nloop    ,is      , 
     $              ie       ,rdir    ,rdif   ,tdir     ,tdif    ,
     $              explay   ,exptdn  ,rdndif ,tottrn   )
C
C Compute reflectivity to direct and diffuse radiation for layers below
C by adding succesive layers starting from the surface and working
C upwards:
C
        do n=1,nloop
          do i=is(n),ie(n)
            rupdir(i,plevp) = albdir(i)
            rupdif(i,plevp) = albdif(i)
          end do
        end do
        do k=plev,0,-1
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) +
     $                      rupdif(i,k+1)*tdif(i,k)**2*rdenom
            end do
          end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,plevp
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $               exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
C
C Compute flux divergence in each layer using the interface up and down
C fluxes:
C
        do k=0,plev
          do n=1,nloop
            do i=is(n),ie(n)
              flxdiv(i,k) = (fluxdn(i,k  ) - fluxdn(i,k+1)) +
     $                      (fluxup(i,k+1) - fluxup(i,k  ))
            end do
          end do
        end do
C
C Monochromatic computation completed; accumulate in totals; adjust
C fraction within spectral interval to allow for the possibility of
C sub-divisions within a particular interval:
C
        psf = 1.0
        if(ph2o(ns).ne.0.) psf = psf*ph2o(ns)
        if(pco2(ns).ne.0.) psf = psf*pco2(ns)
        if(po2 (ns).ne.0.) psf = psf*po2 (ns)
        do n=1,nloop
          do i=is(n),ie(n)
            solflx(i)  = solin(i)*frcsol(ns)*psf
            fsnt(i) = fsnt(i) + solflx(i)*(fluxdn(i,1) - fluxup(i,1))
           
            fsns(i) = fsns(i) + solflx(i)*
     $               (fluxdn(i,plevp) - fluxup(i,plevp))
 
            sfltot(i)  = sfltot(i) + solflx(i)
            fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
            fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
C
C Down spectral fluxes need to be in mks; thus the .001 conversion factors
C
            if (wavmid .lt. 0.7) then
              sols(i) = sols(i) + exptdn(i,plevp)*solflx(i)*0.001
              solsd(i) = solsd(i) + (fluxdn(i,plevp) -
     $                  exptdn(i,plevp)) * solflx(i)*0.001
cKN   added below
                 sabveg(i)=sabveg(i) + (solflx(i)*(fluxdn(i,plevp)-
     $                fluxup(i,plevp)))*(1.-albvs(i))/(1.-albdir(i))
     $                *0.001
cKN   added above
            else
              soll(i) = soll(i) + exptdn(i,plevp)*solflx(i)*0.001
              solld(i) = solld(i) + (fluxdn(i,plevp) -
     $                  exptdn(i,plevp)) * solflx(i)*0.001
              fsnirtsq(i) = fsnirtsq(i) +
     $                      solflx(i)*(fluxdn(i,0) - fluxup(i,0))
cKN   added below
              sabveg(i)=sabveg(i) + (solflx(i)*(fluxdn(i,plevp)-
     $             fluxup(i,plevp)))*(1.-albvl(i))/(1.-albdir(i))
     $             *0.001
cKN   added above
            end if
            fsnirt(i) = fsnirt(i) +
     $                  wgtint * solflx(i)*
     $                  (fluxdn(i,0) - fluxup(i,0))

C
          end do
        end do
        do k=0,plev
          do n=1,nloop
            do i=is(n),ie(n)
              totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv(i,k)
              fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
              fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
            end do
          end do
        end do

c    solis is incident visible solar radiation
        if (ns.eq.8) then
C -trapuv
C        do i=1,plond
C         solis(i)=solflx(i)*0.001*fluxdn(i,plevp)
C        end do
C -trapuv_
          do n=1,nloop
            do i=is(n),ie(n)
              solvs(i)=exptdn(i,plevp)*solflx(i)*0.001
              solvd(i)=(fluxdn(i,plevp)-exptdn(i,plevp))*solflx(i)*0.001
              solis(i)=solflx(i)*0.001*fluxdn(i,plevp)
            end do
          end do
        end if
CEES apr 20

CFAB
C CLEAR SKY CALCULATION PLUS AUTOMATCI CALCULATEION OF AEROSOL FORCING
c RAD CLR is called 2 times , one with O aerosol OP , and one with actual 
c aerosol. DIFFERENCE  in net TOA SW for the two case is saved as one more variable 
c in the rad file. The outputed TOASW ( fsntc, clrst) is accounting for aerosol. 
        if ( idirect .ge. 1 ) then 

          do i=1,plond
            zero(i,1) =0.
            zero(i,2) =0.
            zero(i,3) =0.
            zero(i,4) =0.
          end do

C Following code is the diagnostic clear sky computation:
C
C Compute delta-Eddington solution reflectivities and transmissivities
C for the entire column; note, for convenience, we use the same
C reflectivity and transmissivity arrays as for the full calculation
C above, where 0 for layer quantities refers to the entire atmospheric
C column, and where 0 for interface quantities refers to top of atmos-
C phere, while 1 refers to the surface:
C
          call radclr(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns) ,
     $              abco2(ns),abo2(ns),uth2o   ,uto3     ,utco2    ,
     $              uto2 ,zero(1,1), zero(1,2), zero(1,3), zero(1,4) ,
     $              nloop    ,is      ,ie      ,rdir     ,rdif     ,
     $              tdir     ,tdif    ,explay  ,exptdn   ,rdndif   ,
     $              tottrn   )

       
C
C Compute reflectivity to direct and diffuse radiation for entire
C column; 0,1 on layer quantities refers to two effective layers
C overlying surface; 0 on interface quantities refers to top of column;
C 2 on interface quantities refers to the surface:
C
          do n=1,nloop
            do i=is(n),ie(n)
              rupdir(i,2) = albdir(i)
              rupdif(i,2) = albdif(i)
            end do
          end do
C
          do k=1,0,-1
            do n=1,nloop
              do i=is(n),ie(n)
                rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
                rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
                rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
              end do
            end do
          end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
          do k=0,2
            do n=1,nloop
              do i=is(n),ie(n)
                rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
                fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
                fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
              end do
            end do
          end do
C
          do n=1,nloop
            do i=is(n),ie(n)
             x0fsntc(i)=x0fsntc(i)+solflx(i)*(fluxdn(i,0)-fluxup(i,0))
             x0fsnsc(i)=x0fsnsc(i)+solflx(i)*(fluxdn(i,2)-fluxup(i,2))
             x0fsnrtc(i) = x0fsnrtc(i) +
     $                  wgtint * solflx(i) *
     $                  (fluxdn(i,0) - fluxup(i,0))

c SAVE the ref net TOA flux ( and put back the cumul variables to 0.)  

            end do
          end do
c         if(ns==8)  print *,'DANS radscw FUPav',
c     &    maxval(fluxup(:,0)* solflx(:)*1.E-3)
C
C End of clear sky calculation with O aerosol OP
C
        end if
C
C Following code is the diagnostic clear sky computation:
C
C Compute delta-Eddington solution reflectivities and transmissivities
C for the entire column; note, for convenience, we use the same
C reflectivity and transmissivity arrays as for the full calculation
C above, where 0 for layer quantities refers to the entire atmospheric
C column, and where 0 for interface quantities refers to top of atmos-
C phere, while 1 refers to the surface:
C
        call radclr(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns) ,
     $              abco2(ns),abo2(ns),uth2o   ,uto3     ,utco2    ,
     $              uto2 , tauxar_mix_cs(1,ns), tauasc_mix_cs(1,ns), 
     $              gtota_mix_cs(1,ns), ftota_mix_cs(1,ns)       ,
     $              nloop    ,is      ,ie      ,rdir     ,rdif     ,
     $              tdir     ,tdif    ,explay  ,exptdn   ,rdndif   ,
     $              tottrn   )
        
C
C Compute reflectivity to direct and diffuse radiation for entire
C column; 0,1 on layer quantities refers to two effective layers
C overlying surface; 0 on interface quantities refers to top of column;
C 2 on interface quantities refers to the surface:
C
   
        do n=1,nloop
          do i=is(n),ie(n)

            rupdir(i,2) = albdir(i)
            rupdif(i,2) = albdif(i)

          end do
        end do
C
        do k=1,0,-1
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
            end do
          end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,2
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
C
        do n=1,nloop
          do i=is(n),ie(n)
            fsntc(i) = fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
            fsnsc(i) = fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
            fsnrtc(i) = fsnrtc(i) +
     $                  wgtint * solflx(i) *
     $                  (fluxdn(i,0) - fluxup(i,0))

          end do
        end do

c       if(ns==8)  print *,'DANS radscw FUPav',
c     &    maxval(fluxup(:,0)* solflx(:)*1.E-3)

C
C End of clear sky calculation
C
  100 continue                  ! End of spectral interval loop

C FAB calculation of TOA aerosol radiative forcing 
      if (idirect .ge. 1) then
        do n=1,nloop
          do i=is(n),ie(n)
c test
            aeradfo(i) = -(x0fsntc(i) - fsntc(i)) 
            aeradfos(i)= -(x0fsnsc(i) - fsnsc(i)) 
          end do
        end do
cTEST      
      end if
C
C Compute solar heating rate (k/s)
C
      gocp = gravit/cpair
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k) - pint(i,k+1))
          end do
        end do
      end do
c
c Set the downwelling flux at the surface
c
      do i=1,plon
        fsds(i) = fswdn(i,plevp)
      end do
C
      return
      end subroutine radcsw
