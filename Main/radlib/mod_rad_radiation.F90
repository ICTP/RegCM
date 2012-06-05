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
 
module mod_rad_radiation

  use mod_dynparam
  use mod_mpmessage
  use mod_service
  use mod_runparams , only : iemiss , idirect
  use mod_memutil

! Used by this module only

  use mod_rad_common
  use mod_rad_tracer
  use mod_rad_scenarios
  use mod_rad_aerosol

  private

! Maximum total cloud fraction for radiation model
  real(dp) , public :: cftotmax

  public :: allocate_mod_rad_radiation , radini , radctl

  integer :: npoints

  real(dp) , pointer , dimension(:) :: co2plk , dtx , dty
  real(dp) , pointer , dimension(:) :: xptrop , dlat
  real(dp) , pointer , dimension(:) :: tco2 , th2o , to3 , xsum
  real(dp) , pointer , dimension(:,:) :: co2ems , emstrc , h2oems , o3ems
  real(dp) , pointer , dimension(:,:) :: dbvtit , pnmsq , term6 , term9
  real(dp) , pointer , dimension(:) :: abstrc , dw , pnew , to3co2 , ux
  real(dp) , pointer , dimension(:,:) :: emplnk , pinpl , uinpl , winpl , tbar
  real(dp) , pointer , dimension(:,:,:) :: bplnk
  real(dp) , pointer , dimension(:) :: diralb , difalb
  real(dp) , pointer , dimension(:,:,:) :: absnxt , abstot , xuinpl
  real(dp) , pointer , dimension(:,:) :: emstot
!
! Trace gas variables
!
! bch4     - pressure factor for ch4
  real(dp) , pointer , dimension(:,:) :: bch4
! bn2o0   - pressure factor for n2o
  real(dp) , pointer , dimension(:,:) :: bn2o0
! bn2o1   - pressure factor for n2o
  real(dp) , pointer , dimension(:,:) :: bn2o1
! co2em   - Layer co2 normalized planck funct. derivative
  real(dp) , pointer , dimension(:,:) :: co2em
! co2t    - Prs wghted temperature path
  real(dp) , pointer , dimension(:,:) :: co2t
! h2otr   - H2o trnmsn for o3 overlap
  real(dp) , pointer , dimension(:,:) :: h2otr
! ucfc11  - CFC11 path length
  real(dp) , pointer , dimension(:,:) :: ucfc11
! ucfc12  - CFC12 path length
  real(dp) , pointer , dimension(:,:) :: ucfc12
! un2o0   - N2O path length
  real(dp) , pointer , dimension(:,:) :: un2o0
! un2o1   - N2O path length (hot band)
  real(dp) , pointer , dimension(:,:) :: un2o1
! uch4    - CH4 path length
  real(dp) , pointer , dimension(:,:) :: uch4
! uco211  - CO2 9.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco211
! uco212  - CO2 9.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco212
! uco213  - CO2 9.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco213
! uco221  - CO2 10.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco221
! uco222  - CO2 10.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco222
! uco223  - CO2 10.4 micron band path length
  real(dp) , pointer , dimension(:,:) :: uco223
! uptype   - continuum path length
  real(dp) , pointer , dimension(:,:) :: uptype
!
! plol     - Ozone prs wghted path length
! plos     - Ozone path length
  real(dp) , pointer , dimension(:,:) :: plol , plos
! tplnka   - Planck fnctn level temperature
  real(dp) , pointer , dimension(:,:) :: tplnka
! tint    - Interface temperature
! tint4   - Interface temperature**4
  real(dp) , pointer , dimension(:,:) :: tint , tint4
! tlayr   - Level temperature
! tlayr4  - Level temperature**4
  real(dp) , pointer , dimension(:,:) :: tlayr , tlayr4
! w       - H2o path
  real(dp) , pointer , dimension(:,:) :: w
! dbvtly   - Level drvtv plnck fnctn for o3
  real(dp) , pointer , dimension(:,:) :: dbvtly
! s2c     - H2o cont amount
! s2t     - H2o cont temperature
  real(dp) , pointer , dimension(:,:) :: s2c , s2t
! co2eml  - Interface co2 normalized planck funct. deriv.
  real(dp) , pointer , dimension(:,:) :: co2eml
! ful     - Total upwards longwave flux
! fsul    - Clear sky upwards longwave flux
! fdl     - Total downwards longwave flux
! fsdl    - Clear sky downwards longwv flux
! rtclrsf - d_one/tclrsf(n,k)
  real(dp) , pointer , dimension(:,:) :: fdl , fsdl , fsul , ful
  real(dp) , pointer , dimension(:,:) :: fdl0 , fsdl0 , fsul0 , ful0
  real(dp) , pointer , dimension(:,:) :: rtclrsf
  real(dp) , pointer , dimension(:) :: taugab , tauray
! tmp     - Temporary
! delt    - Diff t**4 mid layer to top interface
! delt1   - Diff t**4 lower intrfc to mid layer
! tplnke  - Planck fnctn temperature
  real(dp) , pointer , dimension(:) :: delt , delt1 , tmp , tplnke
! s       - Flx integral sum
  real(dp) , pointer , dimension(:,:,:) :: s , s0
! fclb4   - Sig t**4 for cld bottom interfc
! fclt4   - Sig t**4 for cloud top interfc
  real(dp) , pointer , dimension(:,:) :: fclb4 , fclt4
! klov    - Cloud lowest level index
! khiv    - Cloud highest level index
! khivm   - khiv(n) - 1
  integer , pointer , dimension(:) :: khiv , khivm , klov
  logical , pointer , dimension(:) :: seldo , done , start
!
!     These arrays are defined for kz model layers; 0 refers to the
!     extra layer on top:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer reflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transmission for layer
! flxdiv   - Flux divergence for layer
! totfld   - Spectrally summed flux divergence
!
!     Cloud radiative property arrays
!
! tauxcl   - water cloud extinction optical depth
! tauxci   - ice cloud extinction optical depth
! wcl      - liquid cloud single scattering albedo
! gcl      - liquid cloud asymmetry parameter
! fcl      - liquid cloud forward scattered fraction
! wci      - ice cloud single scattering albedo
! gci      - ice cloud asymmetry parameter
! fci      - ice cloud forward scattered fraction
!
!     Layer absorber amounts
!
! uh2o     - Layer absorber amount of h2o
! uo3      - Layer absorber amount of  o3
! uco2     - Layer absorber amount of co2
! uo2      - Layer absorber amount of  o2
!
  real(dp) , pointer , dimension(:,:) :: rdir , rdif , tdir , tdif ,    &
        explay , flxdiv , totfld , wcl , gcl , fcl , &
        wci , gci , fci , uh2o , uo3 , uco2 , uo2
!
!     These arrays are defined at model interfaces; 0 is the top of the
!     extra layer above the model top; kzp1 is the earth surface:
!
! rupdir   - Ref to dir rad for layers below
! rupdif   - Ref to dif rad for layers below
! rdndif   - Ref to dif rad for layers above
! exptdn   - Solar beam exp down transm from top
! tottrn   - Total transmission for layers above
! fluxup   - Up   flux at model interface
! fluxdn   - Down flux at model interface
! pflx     - Interface press, including extra layer
! fswup    - Spectrally summed up flux
! fswdn    - Spectrally summed down flux
!
  real(dp) , pointer , dimension(:,:) :: rupdir , rupdif , rdndif , &
        exptdn , tottrn , fluxup , fluxdn , pflx , fswup , fswdn
!
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
!
  real(dp) , pointer , dimension(:,:,:) :: abplnk1 , abplnk2
!
! zenfac   - Square root of cos solar zenith angle
! solflx   - Solar flux in current interval
! uth2o    - Total column  absorber amount of h2o
! uto3     - Total column  absorber amount of  o3
! utco2    - Total column  absorber amount of co2
! uto2     - Total column  absorber amount of  o2
!
  real(dp) , pointer , dimension(:) :: solflx , utco2 , uth2o ,   &
         uto2 , uto3 , x0fsnsc ,  x0fsntc , zenfac
!
!   o3mmr    - Ozone mass mixing ratio
!   pbr      - Model mid-level pressures (dynes/cm2)
!   pnm      - Model interface pressures (dynes/cm2)
!   rh       - level relative humidity (fraction)
!   plco2    - Prs weighted CO2 path
!   plh2o    - Prs weighted H2O path
!   tclrsf   - Total clear sky fraction, level to space
!   cfc11    - cfc11 mass mixing ratio
!   cfc12    - cfc12 mass mixing ratio
!   ch4      - methane mass mixing ratio
!   n2o      - nitrous oxide mass mixing ratio
!
  real(dp) , pointer , dimension(:) :: fslwdcs
  real(dp) , pointer , dimension(:,:) :: cfc11 , cfc12 , ch4 , n2o , &
          o3mmr , pbr , rh
  real(dp) , pointer , dimension(:,:) :: plco2 , plh2o , pnm , tclrsf
!
  real(dp) , dimension(2) :: a1 , a2 , b1 , b2 , realk , st
  real(dp) , dimension(4) :: c1 , c2 , c3 , c4 , c5 , c6 , c7
  real(dp) :: c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 ,  &
             c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 ,  &
             c28 , c29 , c30 , c31 , c8 , c9 , cfa1
  real(dp) :: co2vmr
  real(dp) , dimension(3,4) :: coefa , coefc , coefe
  real(dp) , dimension(4,4) :: coefb , coefd
  real(dp) , dimension(6,2) :: coeff , coefi
  real(dp) , dimension(2,4) :: coefg , coefh
  real(dp) , dimension(3,2) :: coefj , coefk
!
  real(dp) , parameter :: verynearone = 0.999999D0
!
! r80257   - Conversion factor for h2o pathlength
  real(dp) , parameter :: r80257 = d_one/8.0257D-04
  real(dp) , parameter :: r293 = d_one/293.0D0
  real(dp) , parameter :: r250 = d_one/250.0D0
! r3205    - Line width factor for o3 (see R&Di)
  real(dp) , parameter :: r3205 = d_one/0.3205D0
  real(dp) , parameter :: r300 = d_one/300.0D0
! r2sslp   - 1/2 of rsslp
  real(dp) , parameter :: r2sslp = d_one/(d_two*sslp)
! r296   - Inverse stand temp for h2o continuum
  real(dp) , parameter :: r296 = d_one/296.0D0
! repsil - Inver ratio mol weight h2o to dry air
  real(dp) , parameter :: repsil = d_one/ep2
!
! Initialize further longwave constants referring to far wing
! correction; R&D refers to:
!
! Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
! Emissivity and Absorptivity Formulation for Water Vapor
! Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
  real(dp) , parameter :: fwcoef = 0.1D0      ! See eq(33) R&D
  real(dp) , parameter :: fwc1 = 0.30D0       ! See eq(33) R&D
  real(dp) , parameter :: fwc2 = 4.5D0        ! See eq(33) and eq(34) in R&D
  real(dp) , parameter :: fc1 = 2.6D0         ! See eq(34) R&D
!
! Initialize ozone data.
!
  real(dp) , parameter :: v0 = 22.4136D0  ! Volume of a gas at stp (m**3/kmol)
  real(dp) , parameter :: p0 = 0.1D0*sslp ! Standard pressure (pascals)
!
! Constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
  real(dp) , parameter :: cplos = v0/(amd*egrav)*d_100
  real(dp) , parameter :: cplol = v0/(amd*egrav*p0)*d_half*d_100
!
! v_raytau_xx - Constants for new bands
! v_abo3_xx   - Constants for new bands
!
  real(dp) , parameter :: v_raytau_35 = 0.155208D0
  real(dp) , parameter :: v_raytau_64 = 0.0392D0
  real(dp) , parameter :: v_abo3_35   = 2.4058030D+01
  real(dp) , parameter :: v_abo3_64   = 2.210D+01
!
! delta    - Pressure (atmospheres) for stratos. h2o limit
! o2mmr    - O2 mass mixing ratio
!
  real(dp) , parameter :: delta = 1.70D-3
  real(dp) , parameter :: o2mmr = 0.23143D0
!
! Minimum total transmission below which no layer computation are done:
!
! trmin   - Minimum total transmission allowed
! wray    - Rayleigh single scatter albedo
! gray    - Rayleigh asymetry parameter
! fray    - Rayleigh forward scattered fraction
!
  real(dp) , parameter :: trmin = 1.0D-3
  real(dp) , parameter :: wray = 0.999999D0
  real(dp) , parameter :: gray = 0.0D0
  real(dp) , parameter :: fray = 0.1D0
! 
! A. Slingo's data for cloud particle radiative properties
! (from 'A GCM Parameterization for the Shortwave Properties of Water
! Clouds' JAS vol. 46 may 1989 pp 1419-1427)
!
! abarl    - A coefficient for extinction optical depth
! bbarl    - B coefficient for extinction optical depth
! cbarl    - C coefficient for single particle scat albedo
! dbarl    - D coefficient for single particle scat albedo
! ebarl    - E coefficient for asymmetry parameter
! fbarl    - F coefficient for asymmetry parameter
!
! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
!
! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
! abari    - a coefficient for extinction optical depth
! bbari    - b coefficient for extinction optical depth
! cbari    - c coefficient for single particle scat albedo
! dbari    - d coefficient for single particle scat albedo
! ebari    - e coefficient for asymmetry parameter
! fbari    - f coefficient for asymmetry parameter
!
  real(dp) , dimension(4) :: abari , abarl , bbari , bbarl , cbari , &
                             cbarl , dbari , dbarl , ebari , ebarl , &
                             fbari , fbarl
!
! Next series depends on spectral interval
!
! frcsol   - Fraction of solar flux in each spectral interval
! wavmin   - Min wavelength (micro-meters) of interval
! wavmax   - Max wavelength (micro-meters) of interval
! raytau   - Rayleigh scattering optical depth
! abh2o    - Absorption coefficiant for h2o (cm2/g)
! abo3     - Absorption coefficiant for o3  (cm2/g)
! abco2    - Absorption coefficiant for co2 (cm2/g)
! abo2     - Absorption coefficiant for o2  (cm2/g)
! ph2o     - Weight of h2o in spectral interval
! pco2     - Weight of co2 in spectral interval
! po2      - Weight of o2  in spectral interval
! nirwgt   - Weight for intervals to simulate satellite filter
!
  real(dp) , dimension(nspi) :: abco2 , abh2o , abo2 , abo3 ,   &
                                frcsol , nirwgt , pco2 , ph2o , &
                                po2 , raytau , wavmax , wavmin
!
! H2O DMISSIVITY AND ABSORTIVITY CODFFICIDNTS
!
  data coefa/1.01400D+00 , 6.41695D-03 , 2.85787D-05 , &
             1.01320D+00 , 6.86400D-03 , 2.96961D-05 , &
             1.02920D+00 , 1.01680D-02 , 5.30226D-05 , &
             1.02743D+00 , 9.85113D-03 , 5.00233D-05/
!
  data coefb/8.85675D+00 , -3.51620D-02 ,  2.38653D-04 , &
            -1.71439D-06 ,  5.73841D+00 , -1.91919D-02 , &
             1.65993D-04 , -1.54665D-06 ,  6.64034D+00 , &
             1.56651D-02 , -9.73357D-05 ,  0.00000D+00 , &
             7.09281D+00 ,  1.40056D-02 , -1.15774D-04 , &
             0.00000D+00/
!
  data coefc/9.90127D-01 , 1.22475D-03 , 4.90135D-06 , &
             9.89753D-01 , 1.97081D-03 , 3.42046D-06 , &
             9.75230D-01 , 1.03341D-03 , 0.00000D+00 , &
             9.77366D-01 , 8.60014D-04 , 0.00000D+00/
!
  data coefd/7.03047D-01 , -2.63501D-03 , -1.57023D-06 , &
             0.00000D+00 ,  5.29269D-01 , -3.14754D-03 , &
             4.39595D-06 ,  0.00000D+00 ,  7.88193D-02 , &
             1.31290D-03 ,  4.25827D-06 , -1.23982D-08 , &
             1.62744D-01 ,  2.22847D-03 ,  2.60102D-06 , &
            -4.30133D-08/
!
  data coefe/3.93137D-02 , -4.34341D-05 , 3.74545D-07 , &
             3.67785D-02 , -3.10794D-05 , 2.94436D-07 , &
             7.42500D-02 ,  3.97397D-05 , 0.00000D+00 , &
             7.52859D-02 ,  4.18073D-05 , 0.00000D+00/
!
  data coeff/2.20370D-01 , 1.39719D-03 , -7.32011D-06 , &
            -1.40262D-08 , 2.13638D-10 , -2.35955D-13 , &
             3.07431D-01 , 8.27225D-04 , -1.30067D-05 , &
             3.49847D-08 , 2.07835D-10 , -1.98937D-12/
!
  data coefg/9.04489D+00 , -9.56499D-03 ,  1.80898D+01 , &
            -1.91300D-02 ,  8.72239D+00 , -9.53359D-03 , &
             1.74448D+01 , -1.90672D-02/
!
  data coefh/5.46557D+01 , -7.30387D-02 ,  1.09311D+02 ,  &
            -1.46077D-01 ,  5.11479D+01 , -6.82615D-02 ,  &
             1.02296D+02 , -1.36523D-01/
!
  data coefi/3.31654D-01 , -2.86103D-04 , -7.87860D-06 , &
             5.88187D-08 , -1.25340D-10 , -1.37731D-12 , &
             3.14365D-01 , -1.33872D-03 , -2.15585D-06 , &
             6.07798D-08 , -3.45612D-10 , -9.34139D-15/
!
  data coefj/2.82096D-02 , 2.47836D-04 , 1.16904D-06 , &
             9.27379D-02 , 8.04454D-04 , 6.88844D-06/
!
  data coefk/2.48852D-01 , 2.09667D-03 , 2.60377D-06 , &
             1.03594D+00 , 6.58620D-03 , 4.04456D-06/
!
! Narrow band data for H2O
! 200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
!
  data realk/0.18967069430426D-04 ,  0.70172244841851D-04/
  data st   /0.31930234492350D-03 ,  0.97907319939060D-03/
  data a1   /0.28775403075736D-01 ,  0.23236701470511D-01/
  data a2  /-0.57966222388131D-04 , -0.95105504388411D-04/
  data b1   /0.29927771523756D-01 ,  0.21737073577293D-01/
  data b2  /-0.86322071248593D-04 , -0.78543550629536D-04/
!
  data abarl / 2.817D-02 ,  2.682D-02 , 2.264D-02 , 1.281D-02/
  data bbarl / 1.305D+00 ,  1.346D+00 , 1.454D+00 , 1.641D+00/
  data cbarl /-5.620D-08 , -6.940D-06 , 4.640D-04 , 0.201D+00/
  data dbarl / 1.630D-07 ,  2.350D-05 , 1.240D-03 , 7.560D-03/
  data ebarl / 0.829D+00 ,  0.794D+00 , 0.754D+00 , 0.826D+00/
  data fbarl / 2.482D-03 ,  4.226D-03 , 6.560D-03 , 4.353D-03/
! 
  data abari / 3.4480D-03 , 3.4480D-03 , 3.4480D-03 , 3.44800D-03/
  data bbari / 2.4310D+00 , 2.4310D+00 , 2.4310D+00 , 2.43100D+00/
  data cbari / 1.0000D-05 , 1.1000D-04 , 1.8610D-02 , 0.46658D+00/
  data dbari / 0.0000D+00 , 1.4050D-05 , 8.3280D-04 , 2.05000D-05/
  data ebari / 0.7661D+00 , 0.7730D+00 , 0.7940D+00 , 0.95950D+00/
  data fbari / 5.8510D-04 , 5.6650D-04 , 7.2670D-04 , 1.07600D-04/

  data frcsol/0.001488D0 , 0.001389D0 , 0.001290D0 , 0.001686D0 , &
              0.002877D0 , 0.003869D0 , 0.026336D0 , 0.360739D0 , &
              0.065392D0 , 0.526861D0 , 0.526861D0 , 0.526861D0 , &
              0.526861D0 , 0.526861D0 , 0.526861D0 , 0.526861D0 , &
              0.006239D0 , 0.001834D0 , 0.001834D0/

! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans

  data nirwgt/0.000000D0 , 0.000000D0 , 0.000000D0 , 0.000000D0 , &
              0.000000D0 , 0.000000D0 , 0.000000D0 , 0.000000D0 , &
              0.320518D0 , 1.000000D0 , 1.000000D0 , 1.000000D0 , &
              1.000000D0 , 1.000000D0 , 1.000000D0 , 1.000000D0 , &
              1.000000D0 , 1.000000D0 , 1.000000D0/
 
  data wavmin/0.200D0 , 0.245D0 , 0.265D0 , 0.275D0 , 0.285D0 ,&
              0.295D0 , 0.305D0 , 0.350D0 , 0.640D0 , 0.700D0 ,&
              0.701D0 , 0.701D0 , 0.701D0 , 0.701D0 , 0.702D0 ,&
              0.702D0 , 2.630D0 , 4.160D0 , 4.160D0/
 
  data wavmax/0.245D0 , 0.265D0 , 0.275D0 , 0.285D0 , 0.295D0 ,&
              0.305D0 , 0.350D0 , 0.640D0 , 0.700D0 , 5.000D0 ,&
              5.000D0 , 5.000D0 , 5.000D0 , 5.000D0 , 5.000D0 ,&
              5.000D0 , 2.860D0 , 4.550D0 , 4.550D0/
 
  data raytau/4.0200D0 , 2.1800D0 , 1.7000D0 , 1.4500D0 , &
              1.2500D0 , 1.0850D0 , 0.7300D0 ,            &
              v_raytau_35 , v_raytau_64 ,                 &
              0.0200D0 , 0.0001D0 , 0.0001D0 , 0.0001D0 , &
              0.0001D0 , 0.0001D0 , 0.0001D0 , 0.0001D0 , &
              0.0001D0 , 0.0001D0/
!
! Absorption coefficients
!
  data abh2o/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,  0.000D0 , &
             0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,  0.002D0 , &
             0.035D0 , 0.377D0 , 1.950D0 , 9.400D0 , 44.600D0 , &
           190.000D0 , 0.000D0 , 0.000D0 , 0.000D0/
 
  data abo3/5.370D+04 , 13.080D+04 , 9.292D+04 , 4.530D+04 , &
            1.616D+04 ,  4.441D+03 , 1.775D+02 , v_abo3_35 , &
            v_abo3_64 ,  0.000D+00 , 0.000D+00 , 0.000D+00 , &
            0.000D+00 ,  0.000D+00 , 0.000D+00 , 0.000D+00 , &
            0.000D+00 ,  0.000D+00 , 0.000D+00/
 
  data abco2/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
             0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
             0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
             0.000D0 , 0.094D0 , 0.196D0 , 1.963D0/
 
  data abo2/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,  0.000D0 ,  &
            0.000D0 , 0.000D0 , 0.000D0 , 1.11D-05 , 6.69D-05 , &
            0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,  0.000D0 ,  &
            0.000D0 , 0.000D0 , 0.000D0 , 0.000D0/
!
! Spectral interval weights
!
  data ph2o/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
            0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.505D0 , &
            0.210D0 , 0.120D0 , 0.070D0 , 0.048D0 , 0.029D0 , &
            0.018D0 , 0.000D0 , 0.000D0 , 0.000D0/
 
  data pco2/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
            0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
            0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
            0.000D0 , 1.000D0 , 0.640D0 , 0.360D0/
 
  data po2/0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
           0.000D0 , 0.000D0 , 0.000D0 , 1.000D0 , 1.000D0 , &
           0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , &
           0.000D0 , 0.000D0 , 0.000D0 , 0.000D0/
!
  contains 

  subroutine allocate_mod_rad_radiation 
    implicit none        
!
    npoints = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem3d(absnxt,1,npoints,1,kz,1,4,'rad:absnxt')
    call getmem3d(xuinpl,1,npoints,1,kz,1,4,'rad:xuinpl')
    call getmem3d(abstot,1,npoints,1,kzp1,1,kzp1,'rad:abstot')
    call getmem2d(emstot,1,npoints,1,kzp1,'rad:emstot')
    call getmem1d(diralb,1,npoints,'rad:diralb')
    call getmem1d(difalb,1,npoints,'rad:difalb')
    call getmem1d(co2plk,1,npoints,'rad:co2plk')
    call getmem1d(dtx,1,npoints,'rad:dtx')
    call getmem1d(dty,1,npoints,'rad:dty')
    call getmem1d(tco2,1,npoints,'rad:tco2')
    call getmem1d(th2o,1,npoints,'rad:th2o')
    call getmem1d(to3,1,npoints,'rad:to3')
    call getmem1d(xsum,1,npoints,'rad:xsum')
    call getmem1d(abstrc,1,npoints,'rad:abstrc')
    call getmem1d(dw,1,npoints,'rad:dw')
    call getmem1d(pnew,1,npoints,'rad:pnew')
    call getmem1d(to3co2,1,npoints,'rad:to3co2')
    call getmem1d(ux,1,npoints,'rad:ux')
    call getmem1d(taugab,1,npoints,'rad:taugab')
    call getmem1d(tauray,1,npoints,'rad:tauray')
    call getmem1d(seldo,1,npoints,'rad:seldo')
    call getmem1d(khiv,1,npoints,'rad:khiv')
    call getmem1d(khivm,1,npoints,'rad:khivm')
    call getmem1d(klov,1,npoints,'rad:klov')
    call getmem1d(delt,1,npoints,'rad:delt')
    call getmem1d(delt1,1,npoints,'rad:delt1')
    call getmem1d(tmp,1,npoints,'rad:tmp')
    call getmem1d(tplnke,1,npoints,'rad:tplnke')
    call getmem1d(done,1,npoints,'rad:done')
    call getmem1d(start,1,npoints,'rad:start')

    call getmem2d(co2ems,1,npoints,1,kzp1,'rad:co2ems')
    call getmem2d(emstrc,1,npoints,1,kzp1,'rad:emstrc')
    call getmem2d(h2oems,1,npoints,1,kzp1,'rad:h2oems')
    call getmem2d(o3ems,1,npoints,1,kzp1,'rad:o3ems')
    call getmem2d(dbvtit,1,npoints,1,kzp1,'rad:dbvtit')
    call getmem2d(pnmsq,1,npoints,1,kzp1,'rad:pnmsq')
    call getmem2d(term6,1,npoints,1,kzp1,'rad:term6')
    call getmem2d(term9,1,npoints,1,kzp1,'rad:term9')
    call getmem2d(bch4,1,npoints,1,kzp1,'rad:bch4')
    call getmem2d(bn2o0,1,npoints,1,kzp1,'rad:bn2o0')
    call getmem2d(bn2o1,1,npoints,1,kzp1,'rad:bn2o1')
    call getmem2d(co2em,1,npoints,1,kzp1,'rad:co2em')
    call getmem2d(co2t,1,npoints,1,kzp1,'rad:co2t')
    call getmem2d(h2otr,1,npoints,1,kzp1,'rad:h2otr')
    call getmem2d(ucfc11,1,npoints,1,kzp1,'rad:ucfc11')
    call getmem2d(ucfc12,1,npoints,1,kzp1,'rad:ucfc12')
    call getmem2d(un2o0,1,npoints,1,kzp1,'rad:un2o0')
    call getmem2d(un2o1,1,npoints,1,kzp1,'rad:un2o1')
    call getmem2d(uch4,1,npoints,1,kzp1,'rad:uch4')
    call getmem2d(uco211,1,npoints,1,kzp1,'rad:uco211')
    call getmem2d(uco212,1,npoints,1,kzp1,'rad:uco212')
    call getmem2d(uco213,1,npoints,1,kzp1,'rad:uco213')
    call getmem2d(uco221,1,npoints,1,kzp1,'rad:uco221')
    call getmem2d(uco222,1,npoints,1,kzp1,'rad:uco222')
    call getmem2d(uco223,1,npoints,1,kzp1,'rad:uco223')
    call getmem2d(uptype,1,npoints,1,kzp1,'rad:uptype')
    call getmem2d(plol,1,npoints,1,kzp1,'rad:plol')
    call getmem2d(plos,1,npoints,1,kzp1,'rad:plos')
    call getmem2d(tplnka,1,npoints,1,kzp1,'rad:tplnka')
    call getmem2d(tint,1,npoints,1,kzp1,'rad:tint')
    call getmem2d(tint4,1,npoints,1,kzp1,'rad:tint4')
    call getmem2d(tlayr,1,npoints,1,kzp1,'rad:tlayr')
    call getmem2d(tlayr4,1,npoints,1,kzp1,'rad:tlayr4')
    call getmem2d(w,1,npoints,1,kzp1,'rad:w')
    call getmem2d(s2c,1,npoints,1,kzp1,'rad:s2c')
    call getmem2d(s2t,1,npoints,1,kzp1,'rad:s2t')
    call getmem2d(co2eml,1,npoints,1,kzp1,'rad:co2eml')
    call getmem2d(fdl,1,npoints,1,kzp1,'rad:fdl')
    call getmem2d(fsdl,1,npoints,1,kzp1,'rad:fsdl')
    call getmem2d(ful,1,npoints,1,kzp1,'rad:ful')
    call getmem2d(fsul,1,npoints,1,kzp1,'rad:fsul')
    call getmem2d(fdl0,1,npoints,1,kzp1,'rad:fdl0')
    call getmem2d(fsdl0,1,npoints,1,kzp1,'rad:fsdl0')
    call getmem2d(ful0,1,npoints,1,kzp1,'rad:ful0')
    call getmem2d(fsul0,1,npoints,1,kzp1,'rad:fsul0')
    call getmem2d(rtclrsf,1,npoints,1,kzp1,'rad:rtclrsf')

    call getmem2d(dbvtly,1,npoints,1,kz,'rad:dbvtly')
    call getmem2d(fclb4,1,npoints,1,kz,'rad:fclb4')
    call getmem2d(fclt4,1,npoints,1,kz,'rad:fclt4')

    call getmem2d(emplnk,1,14,1,npoints,'rad:emplnk')
    call getmem3d(bplnk,1,14,1,npoints,1,4,'rad:bplnk')

    call getmem3d(abplnk1,1,14,1,npoints,1,kzp1,'rad:abplnk1')
    call getmem3d(abplnk2,1,14,1,npoints,1,kzp1,'rad:abplnk2')

    call getmem2d(exptdn,1,npoints,0,kzp1,'rad:exptdn')
    call getmem2d(fluxdn,1,npoints,0,kzp1,'rad:fluxdn')
    call getmem2d(fluxup,1,npoints,0,kzp1,'rad:fluxup')
    call getmem2d(rdndif,1,npoints,0,kzp1,'rad:rdndif')
    call getmem2d(rupdif,1,npoints,0,kzp1,'rad:rupdif')
    call getmem2d(rupdir,1,npoints,0,kzp1,'rad:rupdir')
    call getmem2d(tottrn,1,npoints,0,kzp1,'rad:tottrn')
    call getmem2d(pflx,1,npoints,0,kzp1,'rad:pflx')
    call getmem2d(fswup,1,npoints,0,kzp1,'rad:fswup')
    call getmem2d(fswdn,1,npoints,0,kzp1,'rad:fswdn')

    call getmem2d(rdir,1,npoints,0,kz,'rad:rdir')
    call getmem2d(rdif,1,npoints,0,kz,'rad:rdif')
    call getmem2d(tdir,1,npoints,0,kz,'rad:tdir')
    call getmem2d(tdif,1,npoints,0,kz,'rad:tdif')
    call getmem2d(explay,1,npoints,0,kz,'rad:explay')
    call getmem2d(flxdiv,1,npoints,0,kz,'rad:flxdiv')
    call getmem2d(totfld,1,npoints,0,kz,'rad:totfld')
    call getmem2d(wcl,1,npoints,0,kz,'rad:wcl')
    call getmem2d(gcl,1,npoints,0,kz,'rad:gcl')
    call getmem2d(fcl,1,npoints,0,kz,'rad:fcl')
    call getmem2d(wci,1,npoints,0,kz,'rad:wci')
    call getmem2d(gci,1,npoints,0,kz,'rad:gci')
    call getmem2d(fci,1,npoints,0,kz,'rad:fci')
    call getmem2d(uh2o,1,npoints,0,kz,'rad:uh2o')
    call getmem2d(uo3,1,npoints,0,kz,'rad:uo3')
    call getmem2d(uco2,1,npoints,0,kz,'rad:uco2')
    call getmem2d(uo2,1,npoints,0,kz,'rad:uo2')

    call getmem3d(s,1,npoints,1,kzp1,1,kzp1,'rad:s')
    call getmem3d(s0,1,npoints,1,kzp1,1,kzp1,'rad:s0')
    call getmem2d(pinpl,1,npoints,1,4,'rad:pinpl')
    call getmem2d(uinpl,1,npoints,1,4,'rad:uinpl')
    call getmem2d(winpl,1,npoints,1,4,'rad:winpl')
    call getmem2d(tbar,1,npoints,1,4,'rad:tbar')

    call getmem1d(solflx,1,npoints,'rad:solflx')
    call getmem1d(utco2,1,npoints,'rad:utco2')
    call getmem1d(uth2o,1,npoints,'rad:uth2o')
    call getmem1d(uto2,1,npoints,'rad:uto2')
    call getmem1d(uto3,1,npoints,'rad:uto3')
    call getmem1d(x0fsnsc,1,npoints,'rad:x0fsnsc')
    call getmem1d(x0fsntc,1,npoints,'rad:x0fsntc')
    call getmem1d(zenfac,1,npoints,'rad:zenfac')

    call getmem1d(fslwdcs,1,npoints,'rad:fslwdcs')
    call getmem2d(cfc11,1,npoints,1,kz,'rad:cfc11')
    call getmem2d(cfc12,1,npoints,1,kz,'rad:cfc12')
    call getmem2d(ch4,1,npoints,1,kz,'rad:ch4')
    call getmem2d(n2o,1,npoints,1,kz,'rad:n2o')
    call getmem2d(o3mmr,1,npoints,1,kz,'rad:o3mmr')
    call getmem2d(pbr,1,npoints,1,kz,'rad:pbr')

    call getmem2d(plco2,1,npoints,1,kzp1,'rad:plco2')
    call getmem2d(plh2o,1,npoints,1,kzp1,'rad:plh2o')
    call getmem2d(pnm,1,npoints,1,kzp1,'rad:pnm')
    call getmem2d(tclrsf,1,npoints,1,kzp1,'rad:tclrsf')
!
  end subroutine allocate_mod_rad_radiation 
!
!-----------------------------------------------------------------------
!
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
  subroutine radini(iyear)
!
    implicit none
    integer , intent(in) :: iyear
!
!   iband  - H2O band index
!
    integer :: iband

    character (len=64) :: subroutine_name='radini'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Set general radiation consts; convert to cgs units where
!   appropriate:
!
!IPCC
!1991-1995
!   co2vmr  =  3.55e-4
!   ch40 = 0.55241 * 1.714e-6
!   n2o0 = 1.51913 * 0.311e-6
!1961-1965
!   co2vmr  =  3.10e-4
!   ch40 = 0.55241 * 1.414e-6
!   n2o0 = 1.51913 * 0.287e-6
!  
!   cfc110 = 4.69548 * 0.280e-9
!   cfc120 = 4.14307 * 0.503e-9
!   co2mmr = 1.51913 * co2vmr
    if ( iyear >= 1750 .and. iyear <= 2100 ) then
      co2vmr = cgas(2,iyear)*1.0D-6
      co2mmr = co2vmr*44.0D0/28.9644D0
      ch40 = cgas(3,iyear)*1.0D-9*0.55241D0
      n2o0 = cgas(4,iyear)*1.0D-9*1.51913D0
      cfc110 = cgas(5,iyear)*1.0D-12*4.69548D0
      cfc120 = cgas(6,iyear)*1.0D-12*4.14307D0
    else
      write (aline,*) 'Loading gas scenario for simulation year: ', iyear
      call say
      call fatal(__FILE__,__LINE__,                                   &
            'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)')
    end if
!
!   Coefficients for h2o emissivity and absorptivity.
!
    do iband = 1 , 4
      c1(iband) = coefe(3,iband)/coefe(2,iband)
      c2(iband) = coefb(3,iband)/coefb(2,iband)
      c3(iband) = coefb(4,iband)/coefb(3,iband)
      c4(iband) = coefd(3,iband)/coefd(2,iband)
      c5(iband) = coefd(4,iband)/coefd(3,iband)
      c6(iband) = coefa(3,iband)/coefa(2,iband)
      c7(iband) = coefc(3,iband)/coefc(2,iband)
    end do
    c8 = coeff(3,1)/coeff(2,1)
    c9 = coeff(3,2)/coeff(2,2)
    c10 = coeff(4,1)/coeff(3,1)
    c11 = coeff(4,2)/coeff(3,2)
    c12 = coeff(5,1)/coeff(4,1)
    c13 = coeff(5,2)/coeff(4,2)
    c14 = coeff(6,1)/coeff(5,1)
    c15 = coeff(6,2)/coeff(5,2)
    c16 = coefj(3,1)/coefj(2,1)
    c17 = coefk(3,1)/coefk(2,1)
    c18 = coefi(3,1)/coefi(2,1)
    c19 = coefi(3,2)/coefi(2,2)
    c20 = coefi(4,1)/coefi(3,1)
    c21 = coefi(4,2)/coefi(3,2)
    c22 = coefi(5,1)/coefi(4,1)
    c23 = coefi(5,2)/coefi(4,2)
    c24 = coefi(6,1)/coefi(5,1)
    c25 = coefi(6,2)/coefi(5,2)
    c26 = coefj(3,2)/coefj(2,2)
    c27 = coefk(3,2)/coefk(2,2)
    c28 = d_half
    c29 = 0.002053D0
    c30 = 0.1D0
    c31 = 3.0D-5
    cfa1 = 0.61D0
!
    call time_end(subroutine_name,indx)
  end subroutine radini
!
!-----------------------------------------------------------------------
!
! Driver for radiation computation.
!
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Calling sequence:
!
!     radinp      Converts units of model fields and computes ozone
!                 mixing ratio for solar scheme
!
!     radcsw      Performs solar computation
!       radalb    Computes surface albedos
!       radded    Computes delta-Eddington solution
!       radclr    Computes diagnostic clear sky fluxes
!
!     radclw      Performs longwave computation
!
!       radtpl    Computes path quantities
!       radems    Computes emissivity
!       radabs    Computes absorptivity
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
! Modified:          B. Briegleb, March 1995 to add aerosol
!                    to shortwave code
!
!-----------------------------------------------------------------------
!
  subroutine radctl(n1,n2,dlat,xptrop,ts,pmid,pint,pmln,piln,   &
                    t,h2ommr,rh,cld,effcld,clwp,aermmr,fsns,qrs,qrl,  &
                    flwds,rel,rei,fice,sols,soll,solsd,solld,emsvt,   &
                    fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin,alb, &
                    albc,fsds,fsnirt,fsnrtc,fsnirtsq,totcf,eccf,o3vmr,&
                    czen,czengt0,adirsw,adifsw,adirlw,adiflw,asw,alw, &
                    abv,sol,aeradfo,aeradfos,aerlwfo,aerlwfos,        &
                    absgasnxt,absgastot,emsgastot,tauxcl,tauxci,labsem)
!
    implicit none
!
!     Input arguments
!
!   ts      - Surface (skin) temperature
!   pmid    - Model level pressures
!   pint    - Model interface pressures
!   pmln    - Natural log of pmid
!   rel     - liquid cloud particle effective radius
!   rei     - ice effective drop size (microns)
!   fice    - fractional ice content within cloud
!   piln    - Natural log of pint
!   t       - Model level temperatures
!   h2ommr  - Model level specific humidity
!   cld     - Fractional cloud cover
!   effcld  - Effective fractional cloud cover
!   clwp    - Cloud liquid water path
!
!     Output solar arguments
!
!   fsns    - Surface absorbed solar flux
!   sols    - Downward solar rad onto surface (sw direct)
!   soll    - Downward solar rad onto surface (lw direct)
!   solsd   - Downward solar rad onto surface (sw diffuse)
!   solld   - Downward solar rad onto surface (lw diffuse)
!   qrs     - Solar heating rate
!
!       Output longwave arguments
!
!   qrl     - Longwave cooling rate
!   flwds   - Surface down longwave flux
!
    integer , intent(in) :: n1 , n2
    logical , intent(in) :: labsem
    real(dp) , intent(in) :: eccf
    real(dp) , pointer , dimension(:) :: alb , albc , emsvt , &
            flns , flnsc , flnt , flntc , flwds , fsds , fsnirt , fsnirtsq , &
            fsnrtc , fsns , fsnsc , fsnt , fsntc , solin , soll , solld ,    &
            sols , solsd , ts , totcf , aeradfo , aeradfos , aerlwfo ,       &
            aerlwfos , dlat , xptrop , adirsw , adifsw , adirlw , adiflw ,   &
            asw , alw , abv , sol , czen
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: absgasnxt
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: absgastot
    real(dp) , pointer , dimension(:,:) , intent(out) :: emsgastot
    logical , pointer , dimension(:) , intent(in) :: czengt0
    real(dp) , pointer , dimension(:,:) :: cld , effcld , piln , pint
    real(dp) , pointer , dimension(:,:,:) :: tauxcl , tauxci
    real(dp) , pointer , dimension(:,:) :: clwp , fice , h2ommr , pmid ,  &
            pmln , qrl , qrs , rei , rel , t , rh
    real(dp) , pointer , dimension(:,:,:) :: aermmr
    real(dp) , pointer , dimension(:,:) :: o3vmr
    intent (out) alb , albc , abv , sol , tauxcl , tauxci
    intent (out) flns , flnsc , flnt , flntc , flwds , fsds ,       &
                   fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   fsntc , solin , totcf
!
!---------------------------Local variables-----------------------------
!
!   i        - index
!   solin    - Solar incident flux
!   fsnt     - Net column abs solar flux at model top
!   fsntc    - Clear sky total column abs solar flux
!   fsnsc    - Clear sky surface abs solar flux
!   fsnirt   - Near-IR flux absorbed at toa
!   fsnrtc   - Clear sky near-IR flux absorbed at toa
!   fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
!   fsds     - Flux Shortwave Downwelling Surface
!   flnt     - Net outgoing lw flux at model top
!   flns     - Srf longwave cooling (up-down) flux
!   flntc    - Clear sky lw flux at model top
!   flnsc    - Clear sky lw flux at srf (up-down)
!   o3vmr    - Ozone volume mixing ratio
!   eccf     - Earth/sun distance factor
!
    integer :: n , k
!
    character (len=64) :: subroutine_name='radctl'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Instead of interpolating the o3vmr from the time-interpolated
!   values, we pass compute o3vmr in getdat() and pass it directly
!   into radctl(). o3mmr will be computed in radinp().
!
!   Set latitude dependent radiation input
!
    call radinp(n1,n2,pmid,pint,h2ommr,cld,o3vmr,pbr,pnm, &
                plco2,plh2o,tclrsf,o3mmr)
!
!   Solar radiation computation
!
    if ( dosw ) then
!
!     Specify aerosol mass mixing ratio
!
      call aermix(pnm,n1,n2)
 
      call aeroppt(rh,aermmr,pint,n1,n2)
!
      call radcsw(n1,n2,pnm,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,  &
                  solin,qrs,fsns,fsnt,fsds,fsnsc,fsntc,sols,soll,solsd,solld, &
                  fsnirt,fsnrtc,fsnirtsq,adirsw,adifsw,adirlw,adiflw,asw,alw, &
                  abv,sol,czen,czengt0,aeradfo,aeradfos,tauxcl,tauxci)
!
!     Convert units of shortwave fields needed by rest of model from CGS to MKS
!
      do n = n1 , n2

        solin(n) = solin(n)*1.0D-3
        fsnt(n) = fsnt(n)*1.0D-3
        fsns(n) = fsns(n)*1.0D-3
        fsntc(n) = fsntc(n)*1.0D-3
        fsnsc(n) = fsnsc(n)*1.0D-3
!
!       clear sky column partitioning for surface flux  
!       note : should be generalised to the whole column to be
!              really in energy balance !
!
        totcf(n) = d_one
        do k = 1 , kzp1
          totcf(n) = totcf(n) * (d_one - cld(n,k)) 
        end do
        totcf(n) = d_one - totcf(n)

!       maximum cld cover considered        
!       fsns(n) = fsns(n) * maxval(cld(n,:)) + &
!                 fsnsc(n) * (1-maxval(cld(n,:)))
!       random overlap assumption is tocf(n)
!       Now average btw rand ov and maximum cloud cover as fil suggest
        totcf(n) =  d_half * ( totcf(n) + maxval(cld(n,:)) )

!       Fil suggestion of putting a max on column cloud fraction
!       TAO: implement a user-specified CF maximum (default of 0.75d0,
!       as suggested by Erika)
        if ( totcf(n) > cftotmax ) totcf(n) = cftotmax
        if ( totcf(n) < d_zero ) totcf(n) = d_zero

        fsns(n) = fsns(n) * totcf(n) + fsnsc(n) * (d_one-totcf(n))
        fsds(n) = fsds(n)*1.0D-3
        fsnirt(n) = fsnirt(n)*1.0D-3
        fsnrtc(n) = fsnrtc(n)*1.0D-3
        fsnirtsq(n) = fsnirtsq(n)*1.0D-3
      end do
!
!     Calculate/outfld albedo and clear sky albedo
!
      do n = n1 , n2
        if ( solin(n) > d_zero ) then
          alb(n) = (solin(n)-fsnt(n))/solin(n)
        else
          alb(n) = d_zero
        end if
      end do
!
      do n = n1 , n2
        if ( solin(n) > d_zero ) then
          albc(n) = (solin(n)-fsntc(n))/solin(n)
        else
          albc(n) = d_zero
        end if
      end do
    end if
!
!   Longwave radiation computation
!
    if ( dolw ) then
!
!     Specify trace gas mixing ratios
!
      call trcmix(n1,n2,dlat,xptrop,pmid,n2o,ch4,cfc11,cfc12)
!
      call radclw(n1,n2,ts,t,h2ommr,o3vmr,pbr,pnm,pmln, &
                  piln,n2o,ch4,cfc11,cfc12,effcld,tclrsf,qrl,flns,   &
                  flnt,flnsc,flntc,flwds,fslwdcs,emsvt,aerlwfo,      &
                  aerlwfos,absgasnxt,absgastot,emsgastot,labsem)
!
!     Convert units of longwave fields needed by rest of model from CGS to MKS
!
      do n = n1 , n2
        flnt(n) = flnt(n)*1.0D-3
        flns(n) = flns(n)*1.0D-3
        flntc(n) = flntc(n)*1.0D-3
        flnsc(n) = flnsc(n)*1.0D-3
        flwds(n) = flwds(n)*1.0D-3
!FAB
        fslwdcs(n) = fslwdcs(n)*1.0D-3
!       essai clear sky column
!
!       flwds(n) = flwds(n) * maxval(cld((n,:))) + &
!                  flwds(n) * (1-maxval(cld((n,:))))
!       flwds(n) = flwds(n) * maxval(cld(n,:)) + &
!                  fslwdcs(n)*(d_one-maxval(cld(n,:)))
!       flns(n) = flns(n) * maxval(cld(n,:)) + &
!                 flnsc(n)*(d_one-maxval(cld(n,:))) 
!
!       totcf(n) has been calculated for the SW, dolw is always true 
        flwds(n) = flwds(n) * totcf(n) + &
                   fslwdcs(n) * (d_one - totcf(n))
        flns(n) = flns(n) * totcf(n) + &
                  flnsc(n) * (d_one - totcf(n))
      end do
    end if

    call time_end(subroutine_name,indx)
  end subroutine radctl
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!qian 30/06/99,  csm new scheme: hygroscopic growth effect of
!qian            sulfate has been included
!qian            main changed codes: radcsw,
 
!-----------------------------------------------------------------------
!
! Solar radiation code
!
! Basic method is Delta-Eddington as described in:
!
!    Briegleb, Bruce P., 1992: Delta-Eddington
!    Approximation for Solar Radiation in the NCAR Community Climate Model,
!    Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! Two changes to the basic method described above are: (1) the distinction
! between liquid and ice particle clouds, and (2) the addition of an
! aerosol with sulfate radiative properties.
!
! Divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud,
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers
! assuming scattering between layers to be isotropic, and distinguishes
! direct solar beam from scattered radiation.
!
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. czen > 0) computations are done.
!
! Note that an extra layer above the model top layer is added.
!
! cgs units are used.
!
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
!
!
!---------------------------Code history--------------------------------
!
! Modified March 1995 to add aerosols
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Briegleb, May 1996
! 19th Band Added:   W. Collins March 1997
! Merge Wong optics: T. Schneider, Mar 1998
!
!-----------------------------------------------------------------------
!
!     Input arguments
!
! pint    - Interface pressure
! h2ommr  - Specific humidity (h2o mass mix ratio)
! o3mmr   - Ozone mass mixing ratio
! cld     - Fractional cloud cover
! clwp    - Layer liquid water path
! rel     - Liquid effective drop size (microns)
! rei     - Ice effective drop size (microns)
! fice    - Fractional ice content within cloud
! eccf    - Eccentricity factor (d_one/earth-sun dist ** 2)
!
!     Output arguments
!
! solin    - Incident solar flux
! qrs      - Solar heating rate
! fsns     - Surface absorbed solar flux
! fsnt     - Total column absorbed solar flux
! fsds     - Flux Shortwave Downwelling Surface
! fsnsc    - Clear sky surface absorbed solar flux
! fsntc    - Clear sky total column absorbed solar flx
! sols     - Direct solar rad incident on surface (< 0.7)
! soll     - Direct solar rad incident on surface (>= 0.7)
! solsd    - Diffuse mod_solar rad incident on surface (< 0.7)
! solld    - Diffuse mod_solar rad incident on surface (>= 0.7)
! fsnirt   - Near-IR flux absorbed at toa
! fsnrtc   - Clear sky near-IR flux absorbed at toa
! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
!
  subroutine radcsw(n1,n2,pint,h2ommr,o3mmr,cld,clwp,rel,rei,fice,  &
                    eccf,solin,qrs,fsns,fsnt,fsds,fsnsc,fsntc,sols,soll,  &
                    solsd,solld,fsnirt,fsnrtc,fsnirtsq,adirsw,adifsw,     &
                    adirlw,adiflw,asw,alw,abv,sol,czen,czengt0,aeradfo,   &
                    aeradfos,tauxcl,tauxci)
! 
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) :: eccf
    real(dp) , pointer , dimension(:) :: aeradfo , aeradfos , fsds , fsnirt , &
             fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , fsntc , solin , soll , &
             solld , sols , solsd , adirsw , adifsw , adirlw , adiflw , asw , &
             alw , abv , sol , czen
    logical , pointer , dimension(:) , intent(in) :: czengt0
    real(dp) , pointer , dimension(:,:) :: cld , pint
    real(dp) , pointer , dimension(:,:,:) :: tauxcl , tauxci
    real(dp) , pointer , dimension(:,:) :: clwp , fice , h2ommr , o3mmr , &
             qrs , rei , rel
    intent (in) cld , clwp , eccf , fice , h2ommr , o3mmr , pint , rei , rel , &
           adirsw , adifsw , adirlw , adiflw , asw , alw
    intent (out) aeradfo , aeradfos , fsds , qrs , abv , sol , tauxcl , tauxci
    intent (inout) fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   fsntc , solin , soll , solld , sols , solsd
!
!---------------------------Local variables-----------------------------
!
! ns       - Spectral loop index
! i        - Longitude loop index
! k        - Level loop index
! n        - Loop index for daylight
! indxsl   - Index for cloud particle properties
!
!     A. Slingo's data for cloud particle radiative properties (from 'A
!     GCM Parameterization for the Shortwave Properties of Water
!     Clouds' JAS vol. 46 may 1989 pp 1419-1427)
!
! abarl    - A coefficient for extinction optical depth
! bbarl    - B coefficient for extinction optical depth
! cbarl    - C coefficient for single particle scat albedo
! dbarl    - D coefficient for single particle scat albedo
! ebarl    - E coefficient for asymmetry parameter
! fbarl    - F coefficient for asymmetry parameter
! abarli   - A coefficient for current spectral interval
! bbarli   - B coefficient for current spectral interval
! cbarli   - C coefficient for current spectral interval
! dbarli   - D coefficient for current spectral interval
! ebarli   - E coefficient for current spectral interval
! fbarli   - F coefficient for current spectral interval
!
!     Caution... A. Slingo recommends no less than 4.0 micro-meters nor
!     greater than 20 micro-meters
!
!     ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
! abari    - a coefficient for extinction optical depth
! bbari    - b coefficient for extinction optical depth
! cbari    - c coefficient for single particle scat albedo
! dbari    - d coefficient for single particle scat albedo
! ebari    - e coefficient for asymmetry parameter
! fbari    - f coefficient for asymmetry parameter
! abarii   - A coefficient for current spectral interval
! bbarii   - B coefficient for current spectral interval
! cbarii   - C coefficient for current spectral interval
! dbarii   - D coefficient for current spectral interval
! ebarii   - E coefficient for current spectral interval
! fbarii   - F coefficient for current spectral interval
! wgtint   - Weight for specific spectral interval
!
!     Diagnostic and accumulation : note that sfltot, fswup, and
!     fswdn are not used in the computation,but are retained for future
!     use.
!
! sfltot   - Spectrally summed total solar flux
!
!     Various arrays and other constants:
!
! sqrco2   - Square root of the co2 mass mixg ratio
! tmp1     - Temporary constant array
! tmp2     - Temporary constant array
! pdel     - Pressure difference across layer
! path     - Mass path of layer
! xptop    - Lower interface pressure of extra layer
! ptho2    - Used to compute mass path of o2
! ptho3    - Used to compute mass path of o3
! pthco2   - Used to compute mass path of co2
! pthh2o   - Used to compute mass path of h2o
! h2ostr   - Inverse square root h2o mass mixing ratio
! wavmid   - Spectral interval middle wavelength
! trayoslp - Rayleigh optical depth/standard pressure
! tmp1l    - Temporary constant array
! tmp2l    - Temporary constant array
! tmp3l    - Temporary constant array
! tmp1i    - Temporary constant array
! tmp2i    - Temporary constant array
! tmp3i    - Temporary constant array
! rdenom   - Multiple scattering term
! psf      - Frac of solar flux in spect interval
! aeradfo  - spectrally integrated aerosol radiative forcing ( TOA)
!-----------------------------------------------------------------------
!
    real(dp) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli ,  &
               dbarii , dbarli , ebarii , ebarli , fbarii , fbarli ,  &
               h2ostr , path , pdel , psf , pthco2 , pthh2o , ptho2 , &
               ptho3 , xptop , rdenom , sqrco2 , tmp1 , tmp1i ,       &
               tmp1l , tmp2 , tmp2i , tmp2l , tmp3i , tmp3l ,         &
               trayoslp , wavmid , wgtint , sfltot , x0fsnrtc
    integer :: n , indxsl , k , ns
!
    character (len=64) :: subroutine_name='radcsw'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Initialize output fields:
!
    fsds(:) = d_zero
    fsnirt(:) = d_zero
    fsnrtc(:) = d_zero
    fsnirtsq(:) = d_zero
    fsnt(:) = d_zero
    fsns(:) = d_zero
    solin(:) = d_zero
    fsnsc(:) = d_zero
    fsntc(:) = d_zero
    sols(:) = d_zero
    soll(:) = d_zero
    solsd(:) = d_zero
    solld(:) = d_zero
    abv(:) = d_zero
    sol(:) = d_zero
!
    aeradfo(:) = d_zero
    aeradfos(:) = d_zero
    x0fsntc(:) = d_zero
    x0fsnsc(:) = d_zero
!
    qrs(:,:) = d_zero
!
!   Define solar incident radiation and interface pressures:
!
    do n = n1 , n2
      if ( czengt0(n) ) then
        solin(n) = scon*eccf*czen(n)
        pflx(n,0) = d_zero
      end if
    end do
    do k = 1 , kzp1
      do n = n1 , n2
        if ( czengt0(n) ) then
          pflx(n,k) = pint(n,k)
        end if
      end do
    end do
!
!   Compute optical paths:
!   CO2, use old scheme(as constant)
!
    tmp1 = d_half/(egravgts*sslp)
!   co2mmr = co2vmr*(mmwco2/mmwair)
   
    sqrco2 = dsqrt(co2mmr)
    do n = n1 , n2
      if ( czengt0(n) ) then
        xptop = pflx(n,1)
        ptho2 = o2mmr*xptop*regravgts
        ptho3 = o3mmr(n,1)*xptop*regravgts
        pthco2 = sqrco2*(xptop*regravgts)
        h2ostr = dsqrt(d_one/h2ommr(n,1))
        zenfac(n) = dsqrt(czen(n))
        pthh2o = xptop**d_two*tmp1+(xptop*regravgts)*(h2ostr*zenfac(n)*delta)
        uh2o(n,0) = h2ommr(n,1)*pthh2o
        uco2(n,0) = zenfac(n)*pthco2
        uo2(n,0) = zenfac(n)*ptho2
        uo3(n,0) = ptho3
      end if
    end do
!
    tmp2 = delta*regravgts
    do k = 1 , kz
      do n = n1 , n2
        if ( czengt0(n) ) then
          pdel = pflx(n,k+1) - pflx(n,k)
          path = pdel*regravgts
          ptho2 = o2mmr*path
          ptho3 = o3mmr(n,k)*path
          pthco2 = sqrco2*path
          h2ostr = dsqrt(d_one/h2ommr(n,k))
          pthh2o = (pflx(n,k+1)**d_two-pflx(n,k)**d_two) * &
                    tmp1 + pdel*h2ostr*zenfac(n)*tmp2
          uh2o(n,k) = h2ommr(n,k)*pthh2o
          uco2(n,k) = zenfac(n)*pthco2
          uo2(n,k) = zenfac(n)*ptho2
          uo3(n,k) = ptho3
        end if
      end do
    end do
!
!   Compute column absorber amounts for the clear sky computation:
!
    do n = n1 , n2
      if ( czengt0(n) ) then
        uth2o(n) = d_zero
        uto3(n) = d_zero
        utco2(n) = d_zero
        uto2(n) = d_zero
      end if
    end do
    do k = 1 , kz
      do n = n1 , n2
        if ( czengt0(n) ) then
          uth2o(n) = uth2o(n) + uh2o(n,k)
          uto3(n) = uto3(n) + uo3(n,k)
          utco2(n) = utco2(n) + uco2(n,k)
          uto2(n) = uto2(n) + uo2(n,k)
        end if
      end do
    end do
!
!   Initialize spectrally integrated totals:
!
    do k = 0 , kz
      do n = n1 , n2
        totfld(n,k) = d_zero
        fswup(n,k) = d_zero
        fswdn(n,k) = d_zero
      end do
    end do
    do n = n1 , n2
      fswup(n,kzp1) = d_zero
      fswdn(n,kzp1) = d_zero
    end do
!
!   Set cloud properties for top (0) layer; so long as tauxcl is zero,
!   there is no cloud above top of model; the other cloud properties
!   are arbitrary:
!
    do n = n1 , n2
      if ( czengt0(n) ) then
        tauxcl(n,0,:) = d_zero
        wcl(n,0) = verynearone
        gcl(n,0) = 0.850D0
        fcl(n,0) = 0.725D0
        tauxci(n,0,:) = d_zero
        wci(n,0) = verynearone
        gci(n,0) = 0.850D0
        fci(n,0) = 0.725D0
      end if
    end do
!
!   Begin spectral loop
!
    do ns = 1 , nspi
      wgtint = nirwgt(ns)
!
!     Set index for cloud particle properties based on the wavelength,
!     according to A. Slingo (1989) equations 1-3:
!     Use index 1 (0.25 to 0.69 micrometers) for visible
!     Use index 2 (0.69 - 1.19 micrometers) for near-infrared
!     Use index 3 (1.19 to 2.38 micrometers) for near-infrared
!     Use index 4 (2.38 to 4.00 micrometers) for near-infrared
!
!     Note that the minimum wavelength is encoded (with 0.001, 0.002,
!     0.003) in order to specify the index appropriate for the
!     near-infrared cloud absorption properties
!
      indxsl = 0
      if ( wavmax(ns) <= 0.70D0 ) then
        indxsl = 1
      else if ( dabs(wavmin(ns)-0.700D0) < dlowval ) then
        indxsl = 2
      else if ( dabs(wavmin(ns)-0.701D0) < dlowval ) then
        indxsl = 3
      else if ( dabs(wavmin(ns)-0.702D0) < dlowval .or. &
                     wavmin(ns) > 2.38D0 ) then
        indxsl = 4
      end if
!
!     Set cloud extinction optical depth, single scatter albedo,
!     asymmetry parameter, and forward scattered fraction:
!
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)
!
      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)
!
      do k = 1 , kz
        do n = n1 , n2
          if ( czengt0(n) ) then
!
!           liquid
!
            tmp1l = abarli + bbarli/rel(n,k)
            tmp2l = d_one - cbarli - dbarli*rel(n,k)
            tmp3l = fbarli*rel(n,k)
!
!           ice
!
            tmp1i = abarii + bbarii/rei(n,k)
            tmp2i = d_one - cbarii - dbarii*rei(n,k)
            tmp3i = fbarii*rei(n,k)
!
!           Cloud fraction incorporated into cloud extinction optical depth
!found
!           April 12 2000, Filippo found the different scheme here:
   
!scheme     1
!ccm3.6.6
!           tauxcl(n,k,ns) = clwp(n,k)*tmp1l*(d_one-fice(n,k))*cld(n,k)*dsqrt(cld(n,k))
!           tauxci(n,k,ns) = clwp(n,k)*tmp1i*fice(n,k)*cld(n,k)*dsqrt(cld(n,k))
!
!scheme     2
!KN
            tauxcl(n,k,ns) = clwp(n,k)*tmp1l*(d_one-fice(n,k))*cld(n,k) / &
                          (d_one+(d_one-0.85D0)*(d_one-cld(n,k))*      &
                          clwp(n,k)*tmp1l*(d_one-fice(n,k)))
            tauxci(n,k,ns) = clwp(n,k)*tmp1i*fice(n,k)*cld(n,k) /     &
                          (d_one+(d_one-0.78D0)*(d_one-cld(n,k)) * &
                          clwp(n,k)*tmp1i*fice(n,k))
   
!scheme     3
!EES        below replaced
!           tauxcl(n,k,ns) = clwp(n,k)*tmp1l*(d_one-fice(n,k))*cld(n,k)**0.85
!           tauxci(n,k,ns) = clwp(n,k)*tmp1i*fice(n,k)*cld(n,k)**0.85
!found_
!
!           Do not let single scatter albedo be 1; delta-eddington
!           solution for non-conservative case:
!
!qian       30/06/99        wcl(n,k) = dmin1(tmp2l,0.999999)
            wcl(n,k) = dmin1(tmp2l,verynearone)
            gcl(n,k) = ebarli + tmp3l
            fcl(n,k) = gcl(n,k)*gcl(n,k)
!
            wci(n,k) = dmin1(tmp2i,verynearone)
            gci(n,k) = ebarii + tmp3i
            fci(n,k) = gci(n,k)*gci(n,k)
!
          end if
        end do
      end do
!
!     Set reflectivities for surface based on mid-point wavelength
!
      wavmid = (wavmin(ns)+wavmax(ns))*d_half
!
!     Wavelength less  than 0.7 micro-meter
!
      if ( wavmid < 0.7D0 ) then
        do n = n1 , n2
          if ( czengt0(n) ) then
            diralb(n) = adirsw(n)
            difalb(n) = adifsw(n)
          end if
        end do
!
!       Wavelength greater than 0.7 micro-meter
!
      else
        do n = n1 , n2
          if ( czengt0(n) ) then
            diralb(n) = adirlw(n)
            difalb(n) = adiflw(n)
          end if
        end do
      end if
      trayoslp = raytau(ns)/sslp
!
!     Layer input properties now completely specified; compute the
!     delta-Eddington solution reflectivities and transmissivities
!     for each layer, starting from the top and working downwards:
!
      call radded(n1,n2,trayoslp,czen,czengt0,tauxcl,tauxci,ns)
!
!     Compute reflectivity to direct and diffuse mod_radiation for layers
!     below by adding succesive layers starting from the surface and
!     working upwards:
!
      do n = n1 , n2
        if ( czengt0(n) ) then
          rupdir(n,kzp1) = diralb(n)
          rupdif(n,kzp1) = difalb(n)
        end if
      end do
      do k = kz , 0 , -1
        do n = n1 , n2
          if ( czengt0(n) ) then
            rdenom = d_one/(d_one-rdif(n,k)*rupdif(n,k+1))
            rupdir(n,k) = rdir(n,k) + tdif(n,k) *                    &
                          (rupdir(n,k+1)*explay(n,k)+rupdif(n,k+1) * &
                          (tdir(n,k)-explay(n,k)))*rdenom
            rupdif(n,k) = rdif(n,k) + rupdif(n,k+1)*tdif(n,k)**d_two*rdenom
          end if
        end do
      end do
!
!     Compute up and down fluxes for each interface, using the added
!     atmospheric layer properties at each interface:
!
      do k = 0 , kzp1
        do n = n1 , n2
          if ( czengt0(n) ) then
            rdenom = d_one/(d_one-rdndif(n,k)*rupdif(n,k))
            fluxup(n,k) = (exptdn(n,k)*rupdir(n,k)+                   &
                          (tottrn(n,k)-exptdn(n,k))*rupdif(n,k))*rdenom
            fluxdn(n,k) = exptdn(n,k) +                          &
                          (tottrn(n,k)-exptdn(n,k)+exptdn(n,k) * &
                           rupdir(n,k)*rdndif(n,k))*rdenom
          end if
        end do
      end do
!
!     Compute flux divergence in each layer using the interface up
!     and down fluxes:
!
      do k = 0 , kz
        do n = n1 , n2
          if ( czengt0(n) ) then
            flxdiv(n,k) = (fluxdn(n,k)-fluxdn(n,k+1)) + &
                          (fluxup(n,k+1)-fluxup(n,k))
          end if
        end do
      end do
!
!     Monochromatic computation completed; accumulate in totals;
!     adjust fraction within spectral interval to allow for the
!     possibility of sub-divisions within a particular interval:
!
      psf = d_one
      if ( dabs(ph2o(ns)) > dlowval ) psf = psf*ph2o(ns)
      if ( dabs(pco2(ns)) > dlowval ) psf = psf*pco2(ns)
      if ( dabs(po2(ns)) > dlowval ) psf = psf*po2(ns)
      sfltot = d_zero
      do n = n1 , n2
        if ( czengt0(n) ) then
          solflx(n) = solin(n)*frcsol(ns)*psf
          fsnt(n) = fsnt(n) + solflx(n)*(fluxdn(n,1)-fluxup(n,1))
   
          fsns(n) = fsns(n) + solflx(n)*(fluxdn(n,kzp1)-fluxup(n,kz + 1))
   
          sfltot = sfltot + solflx(n)
          fswup(n,0) = fswup(n,0) + solflx(n)*fluxup(n,0)
          fswdn(n,0) = fswdn(n,0) + solflx(n)*fluxdn(n,0)
!
!           Down spectral fluxes need to be in mks; thus the 0.001
!           conversion factors
          if ( wavmid < 0.7D0 ) then
            sols(n) = sols(n) + exptdn(n,kzp1)*solflx(n)*d_r1000
            solsd(n) = solsd(n) + &
                   (fluxdn(n,kzp1)-exptdn(n,kz + 1))*solflx(n)*d_r1000
!KN         added below
            abv(n) = abv(n) + (solflx(n) *            &
                       (fluxdn(n,kzp1)-fluxup(n,kz + 1)))*  &
                       (d_one-asw(n))/(d_one-diralb(n))*d_r1000
!KN         added above
          else
            soll(n) = soll(n) + exptdn(n,kzp1)*solflx(n)*d_r1000
            solld(n) = solld(n) + &
                   (fluxdn(n,kzp1)-exptdn(n,kz + 1))*solflx(n)*d_r1000
            fsnirtsq(n) = fsnirtsq(n) + solflx(n)*(fluxdn(n,0)-fluxup(n,0))
!KN         added below
            abv(n) = abv(n) + &
                       (solflx(n)*(fluxdn(n,kzp1)-fluxup(n,kz + 1)))* &
                        (d_one-alw(n))/(d_one-diralb(n))*d_r1000
!KN         added above
          end if
          fsnirt(n) = fsnirt(n) + wgtint*solflx(n) * (fluxdn(n,0)-fluxup(n,0))
!
        end if
      end do
      do k = 0 , kz
        do n = n1 , n2
          if ( czengt0(n) ) then
            totfld(n,k) = totfld(n,k) + solflx(n)*flxdiv(n,k)
            fswup(n,k+1) = fswup(n,k+1) + solflx(n)*fluxup(n,k+1)
            fswdn(n,k+1) = fswdn(n,k+1) + solflx(n)*fluxdn(n,k+1)
          end if
        end do
      end do
   
!     solar is incident visible solar radiation
      if ( ns == 8 ) then
        do n = n1 , n2
          if ( czengt0(n) ) then
            sol(n) = solflx(n)*d_r1000*fluxdn(n,kzp1)
          end if
        end do
      end if
   
!FAB
!     CLEAR SKY CALCULATION PLUS AEROSOL
!     FORCING RAD CLR is called 2 times , one with O aerosol OP , and
!     one with actual aerosol. DIFFERENCE  in net TOA SW for the two
!     case is saved as one more variable in the rad file. The
!     outputed TOASW ( fsntc, clrst) is accounting for aerosol.
      if ( lchem .and. idirect >= 1 ) then
   
!       Following code is the diagnostic clear sky computation:
!
!       Compute delta-Eddington solution reflectivities and
!       transmissivities for the entire column; note, for
!       convenience, we use mod_the same reflectivity and transmissivity
!       arrays as for the full calculation above, where 0 for layer
!       quantities refers to the entire atmospheric column, and where
!       0 for interface quantities refers to top of atmos- phere,
!       while 1 refers to the surface:
!
        call radclr(n1,n2,trayoslp,czen,czengt0,ns)
!
!       Compute reflectivity to direct and diffuse mod_radiation for
!       entire column; 0,1 on layer quantities refers to two
!       effective layers overlying surface; 0 on interface quantities
!       refers to top of column; 2 on interface quantities refers to
!       the surface:
        do n = n1 , n2
          if ( czengt0(n) ) then
            rupdir(n,2) = diralb(n)
            rupdif(n,2) = difalb(n)
          end if
        end do
!
        do k = 1 , 0 , -1
          do n = n1 , n2
            if ( czengt0(n) ) then
              rdenom = d_one/(d_one-rdif(n,k)*rupdif(n,k+1))
              rupdir(n,k) = rdir(n,k) + tdif(n,k) *                    &
                            (rupdir(n,k+1)*explay(n,k)+rupdif(n,k+1) * &
                            (tdir(n,k)-explay(n,k)))*rdenom
              rupdif(n,k) = rdif(n,k) + rupdif(n,k+1)*tdif(n,k)**d_two*rdenom
            end if
          end do
        end do
!
!       Compute up and down fluxes for each interface, using the added
!       atmospheric layer properties at each interface:
!
        do k = 0 , 2
          do n = n1 , n2
            if ( czengt0(n) ) then
              rdenom = d_one/(d_one-rdndif(n,k)*rupdif(n,k))
              fluxup(n,k) = (exptdn(n,k)*rupdir(n,k)+(tottrn(n,k) - &
                            exptdn(n,k))*rupdif(n,k))*rdenom
              fluxdn(n,k) = exptdn(n,k) +                           &
                            (tottrn(n,k)-exptdn(n,k)+exptdn(n,k) *  &
                             rupdir(n,k)*rdndif(n,k))*rdenom
            end if
          end do
        end do
!
        x0fsnrtc = d_zero  
        do n = n1 , n2
          if ( czengt0(n) ) then
!           SAVE the ref net TOA flux ( and put back the cumul variables to 0.)
            x0fsntc(n) = x0fsntc(n) + solflx(n)*(fluxdn(n,0)-fluxup(n,0))
            x0fsnsc(n) = x0fsnsc(n) + solflx(n)*(fluxdn(n,2)-fluxup(n,2))
            x0fsnrtc = x0fsnrtc + wgtint*solflx(n)*(fluxdn(n,0)-fluxup(n,0))
          end if
        end do
!
!       End of clear sky calculation with O aerosol OP
!
      end if
!
!     Following code is the diagnostic clear sky computation:
!
!     Compute delta-Eddington solution reflectivities and
!     transmissivities for the entire column; note, for convenience,
!     we use mod_the same reflectivity and transmissivity arrays as for
!     the full calculation above, where 0 for layer quantities refers
!     to the entire atmospheric column, and where 0 for interface
!     quantities refers to top of atmos- phere, while 1 refers to the
!     surface:
!
      call radclr(n1,n2,trayoslp,czen,czengt0,ns)
!
!     Compute reflectivity to direct and diffuse mod_radiation for entire
!     column; 0,1 on layer quantities refers to two effective layers
!     overlying surface; 0 on interface quantities refers to top of
!     column; 2 on interface quantities refers to the surface:
!
      do n = n1 , n2
        if ( czengt0(n) ) then
          rupdir(n,2) = diralb(n)
          rupdif(n,2) = difalb(n)
        end if
      end do
!
      do k = 1 , 0 , -1
        do n = n1 , n2
          if ( czengt0(n) ) then
            rdenom = d_one/(d_one-rdif(n,k)*rupdif(n,k+1))
            rupdir(n,k) = rdir(n,k) + tdif(n,k) *     &
                          (rupdir(n,k+1)*explay(n,k)+ &
                           rupdif(n,k+1)*(tdir(n,k)-explay(n,k)))*rdenom
            rupdif(n,k) = rdif(n,k) + rupdif(n,k+1)*tdif(n,k)**d_two*rdenom
          end if
        end do
      end do
!
!     Compute up and down fluxes for each interface, using the added
!     atmospheric layer properties at each interface:
!
      do k = 0 , 2
        do n = n1 , n2
          if ( czengt0(n) ) then
            rdenom = d_one/(d_one-rdndif(n,k)*rupdif(n,k))
            fluxup(n,k) = (exptdn(n,k)*rupdir(n,k)+(tottrn(n,k) - &
                           exptdn(n,k))*rupdif(n,k))*rdenom
            fluxdn(n,k) = exptdn(n,k) +                           &
                          (tottrn(n,k)-exptdn(n,k)+exptdn(n,k) *  &
                           rupdir(n,k)*rdndif(n,k))*rdenom
          end if
        end do
      end do
!
      do n = n1 , n2
        if ( czengt0(n) ) then
          fsntc(n) = fsntc(n) + solflx(n)*(fluxdn(n,0)-fluxup(n,0))
          fsnsc(n) = fsnsc(n) + solflx(n)*(fluxdn(n,2)-fluxup(n,2))
          fsnrtc(n) = fsnrtc(n) + wgtint*solflx(n) * (fluxdn(n,0)-fluxup(n,0))
        end if
      end do
!
!     End of clear sky calculation
!
    end do  ! End of spectral interval loop
   
!   FAB calculation of TOA aerosol radiative forcing
    if ( lchem .and. idirect >= 1 ) then
      do n = n1 , n2
        if ( czengt0(n) ) then
          aeradfo(n) = -(x0fsntc(n)-fsntc(n))
          aeradfos(n) = -(x0fsnsc(n)-fsnsc(n))
        end if
      end do
    end if
!
!   Compute solar heating rate (k/s)
!
    do k = 1 , kz
      do n = n1 , n2
        if ( czengt0(n) ) then
          qrs(n,k) = -gocp*totfld(n,k)/(pint(n,k)-pint(n,k+1))
        end if
      end do
    end do
!
!   Set the downwelling flux at the surface
!
    do n = n1 , n2
      fsds(n) = fswdn(n,kzp1)
    end do
!
    call time_end(subroutine_name,indx)
  end subroutine radcsw
!
!-----------------------------------------------------------------------
!
! Compute longwave radiation heating rates and boundary fluxes
!
! Uses broad band absorptivity/emissivity method to compute clear sky;
! assumes randomly overlapped clouds with variable cloud emissivity to
! include effects of clouds.
!
! Computes clear sky absorptivity/emissivity at lower frequency (in
! general) than the model radiation frequency; uses previously computed
! and stored values for efficiency
!
! Note: This subroutine contains vertical indexing which proceeds
!       from bottom to top rather than the top to bottom indexing
!       used in the rest of the model.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
!     Input arguments
!
! ts      - Ground (skin) temperature
! emsvt   - Emissivity of surface
!
!     Input arguments which are only passed to other routines
!
! tnm     - Level temperature
! qnm     - Level moisture field
! o3vmr   - ozone volume mixing ratio
! pmid    - Level pressure
! pint    - Model interface pressure
! pmln    - Ln(pmid)
! piln    - Ln(pint)
! plco2   - Path length co2
! plh2o   - Path length h2o
! n2o     - nitrous oxide mass mixing ratio
! ch4     - methane mass mixing ratio
! cfc11   - cfc11 mass mixing ratio
! cfc12   - cfc12 mass mixing ratio
!
!     Input/Output arguments
!
! cld     - Cloud cover
! tclrsf  - Clear sky fraction
!
!     Output arguments
!
! qrl     - Longwave heating rate
! flns    - Surface cooling flux
! flnt    - Net outgoing flux
! flnsc   - Clear sky surface cooing
! flntc   - Net clear sky outgoing flux
! flwds   - Down longwave flux at surface
!
  subroutine radclw(n1,n2,ts,tnm,qnm,o3vmr,pmid,pint,pmln,     &
                    piln,n2o,ch4,cfc11,cfc12,cld,tclrsf,qrl,flns,flnt, &
                    flnsc,flntc,flwds,fslwdcs,emsvt,aerlwfo,aerlwfos,  &
                    absgasnxt,absgastot,emsgastot,labsem)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    logical , intent(in) :: labsem
    real(dp) , pointer , dimension(:,:) :: cfc11 , cfc12 , ch4 , n2o , &
               o3vmr , pmid , pmln , qnm , qrl , tnm
    real(dp) , pointer , dimension(:,:) :: cld , piln , pint , tclrsf
    real(dp) , pointer , dimension(:,:,:) :: absgasnxt , absgastot
    real(dp) , pointer , dimension(:,:) :: emsgastot
    real(dp) , pointer , dimension(:) :: emsvt , flns , flnsc , flnt , &
               flntc , flwds , fslwdcs , ts
    real(dp), pointer , dimension(:) :: aerlwfo , aerlwfos
    intent (in) cld , emsvt
    intent (out) flns , flnsc , flnt , flntc , flwds , qrl , aerlwfo , aerlwfos
    intent (inout) tclrsf
!
!---------------------------Local variables-----------------------------
!
    real(dp) :: absbt , bk1 , bk2 , tmp1
    integer :: n , k , k1 , k2 , k3 , khighest , km , km1 , km2 , &
               km3 , km4 , ns , irad , nradaer
!
    character (len=64) :: subroutine_name='radclw'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
    do n = n1 , n2
      rtclrsf(n,1) = d_one/tclrsf(n,1)
    end do
!
    do k = 1 , kz
      do n = n1 , n2
        fclb4(n,k) = d_zero
        fclt4(n,k) = d_zero
        tclrsf(n,k+1) = tclrsf(n,k)*(d_one-cld(n,k+1))
        rtclrsf(n,k+1) = d_one/tclrsf(n,k+1)
      end do
    end do
!
!   Calculate some temperatures needed to derive absorptivity and
!   emissivity, as well as some h2o path lengths
!
    call radtpl(n1,n2,ts,tnm,pmln,qnm,piln,pint)
   
!   do emissivity and absorptivity calculations
!   only if abs/ems computation
!
    if ( labsem ) then
!
!     Compute ozone path lengths at frequency of a/e calculation.
!
      call radoz2(n1,n2,o3vmr,pint)
!
!     Compute trace gas path lengths
!
      call trcpth(n1,n2,tnm,pint,cfc11,cfc12,n2o,ch4,qnm, &
                  ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212, &
                  uco213,uco221,uco222,uco223,bn2o0,bn2o1,bch4,uptype)
!
!
!     Compute total emissivity:
!
      call radems(n1,n2,pint,emsgastot)
!
!     Compute total absorptivity:
!
      call radabs(n1,n2,pint,pmid,piln,pmln,absgasnxt,absgastot)

    end if
!
!   Find the lowest and highest level cloud for each grid point
!   Note: Vertical indexing here proceeds from bottom to top
!
    do n = n1 , n2
      klov(n) = 0
      done(n) = .false.
    end do
    do k = 1 , kz
      do n = n1 , n2
        if ( .not.done(n) .and. cld(n,kzp2-k) > d_zero ) then
          done(n) = .true.
          klov(n) = k
        end if
      end do
    end do
    where ( klov > 0 )
      seldo = .true.
    elsewhere
      seldo = .false.
    end where
    do n = n1 , n2
      khiv(n) = klov(n)
      done(n) = .false.
    end do
    do k = kz , 1 , -1
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( .not.done(n) .and. cld(n,kzp2-k) > d_zero ) then
          done(n) = .true.
          khiv(n) = k
        end if
      end do
    end do
    do n = n1 , n2
      khivm(n) = khiv(n) - 1
    end do
!
!   Note: Vertical indexing here proceeds from bottom to top
!
    do n = n1 , n2
      if ( .not. seldo(n) ) cycle
      do k = klov(n) , khiv(n)
        fclt4(n,kzp1-k) = stebol*tint4(n,kzp2-k)
        fclb4(n,kzp1-k) = stebol*tint4(n,kzp3-k)
      end do
    end do

!
!   option to calculate LW aerosol radiative forcing
!
!   FAB LW radiative forcing ( rad=1 : avec dust)
    if ( .not. lchem .and. idirect > 0 ) then
      nradaer = 2
      fsul0(:,:) = d_zero
      fsdl0(:,:) = d_zero
    else
      nradaer = 1
    end if

    abstot(:,:,:) = absgastot(:,:,:)
    absnxt(:,:,:) = absgasnxt(:,:,:)
    emstot(:,:)   = emsgastot(:,:)

    do irad = 1 , nradaer

      if ( lchem .and. idirect > 0 .and. irad == 2 ) then
        abstot(:,:,:) = d_one-(d_one-abstot(:,:,:))*aertrlw(:,:,:)
        emstot(:,:) = d_one-(d_one-emstot(:,:))*aertrlw(:,:,1)
        do k = 1 , kz  ! aertrlw defined on plev levels
          do ns = 1 , 4
            absnxt(:,k,ns) = d_one-(d_one-absnxt(:,k,ns)) *   &
                              (aertrlw(:,k,k+1)**xuinpl(:,k,ns))
          end do
        end do
      end if
!
!     Compute sums used in integrals (all longitude points)
!
!     Definition of bk1 & bk2 depends on finite differencing.  for
!     trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
!     layers only.
!
!     delt=t**4 in layer above current sigma level km.
!     delt1=t**4 in layer below current sigma level km.
!
      do n = n1 , n2
        delt(n) = tint4(n,kz) - tlayr4(n,kzp1)
        delt1(n) = tlayr4(n,kzp1) - tint4(n,kzp1)
        s(n,kzp1,kzp1) = stebol*(delt1(n)*absnxt(n,kz,1) + &
                         delt(n)*absnxt(n,kz,4))
        s(n,kz,kzp1) = stebol*(delt(n)*absnxt(n,kz,2) + &
                         delt1(n)*absnxt(n,kz,3))
      end do
      do k = 1 , kz - 1
        do n = n1 , n2
          bk2 = (abstot(n,k,kz)+abstot(n,k,kzp1))*d_half
          bk1 = bk2
          s(n,k,kzp1) = stebol*(bk2*delt(n)+bk1*delt1(n))
        end do
      end do
!
!     All k, km>1
!
      do km = kz , 2 , -1
        do n = n1 , n2
          delt(n) = tint4(n,km-1) - tlayr4(n,km)
          delt1(n) = tlayr4(n,km) - tint4(n,km)
        end do
        do k = kzp1 , 1 , -1
          do n = n1 , n2
            if ( k == km ) then
              bk2 = absnxt(n,km-1,4)
              bk1 = absnxt(n,km-1,1)
            else if ( k == km-1 ) then
              bk2 = absnxt(n,km-1,2)
              bk1 = absnxt(n,km-1,3)
            else
              bk2 = (abstot(n,k,km-1)+ &
                     abstot(n,k,km))*d_half
              bk1 = bk2
            end if
            s(n,k,km) = s(n,k,km+1)+stebol*(bk2*delt(n)+bk1*delt1(n))
          end do
        end do
      end do
!
!     Computation of clear sky fluxes always set first level of fsul
!
      do n = n1 , n2
        if ( iemiss == 1 ) then
          fsul(n,kzp1) = emsvt(n)*(stebol*(ts(n)**d_four))
        else
          fsul(n,kzp1) = stebol*(ts(n)**d_four)
        end if
      end do
!
!     Downward clear sky fluxes store intermediate quantities in down
!     flux Initialize fluxes to clear sky values.
!
      do n = n1 , n2
        tmp(n) = fsul(n,kzp1) - stebol*tint4(n,kzp1)
        fsul(n,1) = fsul(n,kzp1) - abstot(n,1,kzp1)*tmp(n)+s(n,1,2)
        fsdl(n,1) = stebol*(tplnke(n)**d_four)*emstot(n,1)
        ful(n,1) = fsul(n,1)
        fdl(n,1) = fsdl(n,1)
      end do
!
!     fsdl(n,kzp1) assumes isothermal layer
!
      do k = 2 , kz
        do n = n1 , n2
          fsul(n,k) = fsul(n,kzp1) - abstot(n,k,kzp1)*tmp(n)+s(n,k,k+1)
          ful(n,k) = fsul(n,k)
          fsdl(n,k) = stebol*(tplnke(n)**d_four)*emstot(n,k) - &
                              (s(n,k,2)-s(n,k,k+1))
          fdl(n,k) = fsdl(n,k)
        end do
      end do
!
!     Store the downward emission from level 1 = total gas emission *
!     sigma t**4.  fsdl does not yet include all terms
!
      do n = n1 , n2
        ful(n,kzp1) = fsul(n,kzp1)
        absbt = stebol*(tplnke(n)**d_four)*emstot(n,kzp1)
        fsdl(n,kzp1) = absbt - s(n,kzp1,2)
        fdl(n,kzp1) = fsdl(n,kzp1)
      end do

!     FAB radiative forcing sur fsul

      if ( lchem .and. idirect > 0 .and. irad == 1 ) then
        fsul0(:,:) = fsul(:,:)! save fsul0 = no dust
        fsdl0(:,:) = fsdl(:,:)!
        ful0(:,:) = ful(:,:)
        fdl0(:,:) = fdl(:,:)
        s0(:,:,:) = s(:,:,:)
      end if

    end do ! end rad loop

!   FAB after this DO loop fsul account for dust LW effect
!   which is OK in case of idirect=2

    if ( lchem .and. idirect > 0 ) then

      aerlwfo(:) = fsul0(:,1) - fsul(:,1)

!     surface lw net ! fsul(n,plevp) - fsdl(n,plevp)
!     aerlwfos(:)= fsdl0(:,kz)-fsdl(:,kz)
      aerlwfos(:) = (fsul0(:,kzp1)-fsdl0(:,kzp1))- &
                    (fsul(:,kzp1) - fsdl(:,kzp1))
       
!     return to no aerosol LW effect  situation if idirect ==1
      if ( idirect == 1 ) then
        fsul(:,:) = fsul0(:,:)
        fsdl(:,:) = fsdl0(:,:)
        ful(:,:) = ful0(:,:)
        fdl(:,:) = fdl0(:,:)
        s(:,:,:) = s0(:,:,:)
      end if 

    end if ! end aersol rad diagnostic

!   FAB MODIF  : save downward clear sky long wave flux in surface  
    fslwdcs(:) = fsdl(:,kzp1)
!
!   Modifications for clouds
!
!   Further qualify longitude subset for computations.  Select only
!   those locations where there are clouds
!   (total cloud fraction <= 1.e-3 treated as clear)
!
    where ( tclrsf(:,kzp1) < verynearone )
      seldo = .true.
    elsewhere
      seldo = .false.
    end where
!
!   Compute downflux at level 1 for cloudy sky
!
    do n = n1 , n2
      if ( .not. seldo(n) ) cycle
!
!     First clear sky flux plus flux from cloud at level 1
!
      fdl(n,kzp1) = fsdl(n,kzp1)*tclrsf(n,kz) * &
                    rtclrsf(n,kzp1-khiv(n))+fclb4(n,kz-1)*cld(n,kz)
    end do
!
!   Flux emitted by other layers
!   Note: Vertical indexing here proceeds from bottom to top
!
    khighest = khiv(intmax(khiv))
    do km = 3 , khighest
      km1 = kzp1 - km
      km2 = kzp2 - km
      km4 = kzp4 - km
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( km <= khiv(n) ) then
          tmp1 = cld(n,km2)*tclrsf(n,kz)*rtclrsf(n,km2)
          fdl(n,kzp1) = fdl(n,kzp1) + (fclb4(n,km1)-s(n,kzp1,km4))*tmp1
        end if
      end do
    end do
!
!   Note: Vertical indexing here proceeds from bottom to top
!
    do k = 1 , khighest - 1
      k1 = kzp1 - k
      k2 = kzp2 - k
      k3 = kzp3 - k
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( k >= klov(n) .and. k <= khivm(n) ) then
          ful(n,k2) = fsul(n,k2)*(tclrsf(n,kzp1)*rtclrsf(n,k1))
        end if
      end do
      do km = 1 , k
        km1 = kzp1 - km
        km2 = kzp2 - km
        km3 = kzp3 - km
        do n = n1 , n2
          if ( .not. seldo(n) ) cycle
!
          if ( k <= khivm(n) .and. km >= klov(n) .and. km <= khivm(n)) then
            ful(n,k2) = ful(n,k2) + (fclt4(n,km1)+s(n,k2,k3)-s(n,k2,km3)) * &
                        cld(n,km2)*(tclrsf(n,km1)*rtclrsf(n,k1))
          end if
        end do
      end do ! km = 1 , k
    end do   ! k = 1 , khighest-1
!
    do k = 1 , kzp1
      k2 = kzp2 - k
      k3 = kzp3 - k
      do n = n1 , n2
        start(n) = .false.
      end do
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( k >= khiv(n) ) then
          start(n) = .true.
          ful(n,k2) = fsul(n,k2)*tclrsf(n,kzp1)*rtclrsf(n,kzp1-khiv(n))
        end if
      end do
      do km = 1 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km3 = kzp3 - km
        do n = n1 , n2
          if ( .not. seldo(n) ) cycle
          if ( start(n) .and. km >= klov(n) .and. km <= khiv(n) ) then
            ful(n,k2) = ful(n,k2) + (cld(n,km2)*tclrsf(n,km1)* &
                  rtclrsf(n,kzp1-khiv(n)))*(fclt4(n,km1)+s(n,k2,k3)-s(n,k2,km3))
          end if
        end do
      end do  ! km = 1 , khighest
    end do    ! k = 1 , kzp1
!
!   Computation of the downward fluxes
!
    do k = 2 , khighest - 1
      k1 = kzp1 - k
      k2 = kzp2 - k
      k3 = kzp3 - k
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( k <= khivm(n) ) fdl(n,k2) = d_zero
      end do
      do km = k + 1 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km4 = kzp4 - km
        do n = n1 , n2
          if ( .not. seldo(n) ) cycle
!
          if ( k <= khiv(n) .and. &
               km >= max0(k+1,klov(n)) .and. km <= khiv(n) ) then
            fdl(n,k2) = fdl(n,k2)+(cld(n,km2)*tclrsf(n,k1)*rtclrsf(n,km2)) * &
                        (fclb4(n,km1)-s(n,k2,km4)+s(n,k2,k3))
          end if
        end do
      end do ! km = k+1 , khighest
      do n = n1 , n2
        if ( .not. seldo(n) ) cycle
        if ( k <= khivm(n) ) then
          fdl(n,k2) = fdl(n,k2) + &
                  fsdl(n,k2)*(tclrsf(n,k1)*rtclrsf(n,kzp1-khiv(n)))
        end if
      end do
    end do  ! k = 1 , khighest-1
!
!   End cloud modification loops
!
!   All longitudes: store history tape quantities
!
    do n = n1 , n2
!
!     Downward longwave flux
!
      flwds(n) = fdl(n,kzp1)
!
!     Net flux
!
      flns(n) = ful(n,kzp1) - fdl(n,kzp1)
!
!     Clear sky flux at top of atmosphere
!
      flntc(n) = fsul(n,1)
      flnsc(n) = fsul(n,kzp1) - fsdl(n,kzp1)
!
!     Outgoing ir
!
      flnt(n) = ful(n,1) - fdl(n,1)
    end do
!
!   Computation of longwave heating (k per sec)
!
    do k = 1 , kz
      do n = n1 , n2
        qrl(n,k) = (ful(n,k)-fdl(n,k)-ful(n,k+1)+fdl(n,k+1))*gocp / &
                    ((pint(n,k)-pint(n,k+1)))
      end do
    end do
!
    call time_end(subroutine_name,indx)

    contains

      integer function intmax(imax)
        implicit none
        integer , pointer , dimension(:) , intent(in) :: imax
        integer :: i , n , is , ie , mx
        is = lbound(imax,1)
        ie = ubound(imax,1)
        intmax = is
        n = ie-is+1
        if ( n == 1 ) return
        mx = imax(is)
        do i = is+1 , ie
          if ( imax(i) > mx ) then
            mx = imax(i)
            intmax = i
          end if
        end do
      end function intmax

  end subroutine radclw
!
!-----------------------------------------------------------------------
!
! Delta-Eddington solution for special clear sky computation
!
! Computes total reflectivities and transmissivities for two atmospheric
! layers: an overlying purely ozone absorbing layer, and the rest of the
! column below.
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! trayoslp - Tray/sslp
!
  subroutine radclr(n1,n2,trayoslp,czen,czengt0,ns)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , intent(in) :: trayoslp
    real(dp) , pointer , dimension(:) , intent(in) :: czen
    logical , pointer , dimension(:) , intent(in) :: czengt0
    integer , intent(in) :: ns
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! ii       - Longitude index
! nval     - Number of long values satisfying criteria
! indx     - Array of longitude indices
! taugab   - Total column gas absorption optical depth
! tauray   - Column rayleigh optical depth
! tautot   - Total column optical depth
! wtot     - Total column single scatter albedo
! gtot     - Total column asymmetry parameter
! ftot     - Total column forward scatter fraction
! ts       - Column scaled extinction optical depth
! ws       - Column scaled single scattering albedo
! gs       - Column scaled asymmetry parameter
! rdenom   - Mulitiple scattering term
! rdirexp  - Layer direct ref times exp transmission
! tdnmexp  - Total transmission minus exp transmission
!
    real(dp) :: alp , amg , apg , arg , extins , ftot ,   &
               gam , gs , gtot , lm , ne , rdenom , rdirexp , &
               tautot , tdnmexp , ts , ue , ws , wtot
    integer :: n , k
!
    character (len=64) :: subroutine_name='radclw'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!-----------------------------------------------------------------------
!
!   Initialize all total transmimission values to 0, so that nighttime
!   values from previous computations are not used:
!
    tottrn(:,:) = d_zero
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above each interface by starting from the top and adding layers down:
!
!   The top layer is assumed to be a purely absorbing ozone layer, and
!   that the mean diffusivity for diffuse mod_transmission is 1.66:
!
    do n = n1 , n2
      if ( czengt0(n) ) then
!
        taugab(n) = abo3(ns)*uto3(n)
!
!       Limit argument of exponential to 25, in case czen is very small:
        arg = dmin1(taugab(n)/czen(n),25.0D0)
        explay(n,0) = dexp(-arg)
        tdir(n,0) = explay(n,0)
!
!       Same limit for diffuse mod_transmission:
!
        arg = dmin1(1.66D0*taugab(n),25.0D0)
        tdif(n,0) = dexp(-arg)
!
        rdir(n,0) = d_zero
        rdif(n,0) = d_zero
!
!       Initialize top interface of extra layer:
!
        exptdn(n,0) = d_one
        rdndif(n,0) = d_zero
        tottrn(n,0) = d_one
!
        rdndif(n,1) = rdif(n,0)
        tottrn(n,1) = tdir(n,0)
!
      end if
    end do
!
!   Now, complete the rest of the column; if the total transmission
!   through the top ozone layer is less than trmin, then no
!   delta-Eddington computation for the underlying column is done:
!
    do k = 1 , 1
!
!     Initialize current layer properties to zero;only if total
!     transmission to the top interface of the current layer exceeds
!     the minimum, will these values be computed below:
!
      do n = n1 , n2
        if ( czengt0(n) ) then
!
          rdir(n,k) = d_zero
          rdif(n,k) = d_zero
          tdir(n,k) = d_zero
          tdif(n,k) = d_zero
          explay(n,k) = d_zero
!
!         Calculates the solar beam transmission, total transmission,
!         and reflectivity for diffuse mod_radiation from below at the
!         top of the current layer:
!
          exptdn(n,k) = exptdn(n,k-1)*explay(n,k-1)
          rdenom = d_one/(d_one-rdif(n,k-1)*rdndif(n,k-1))
          rdirexp = rdir(n,k-1)*exptdn(n,k-1)
          tdnmexp = tottrn(n,k-1) - exptdn(n,k-1)
          tottrn(n,k) = exptdn(n,k-1)*tdir(n,k-1) + &
                        tdif(n,k-1)*(tdnmexp+rdndif(n,k-1)*rdirexp)*rdenom
          rdndif(n,k) = rdif(n,k-1) + &
                        (rdndif(n,k-1)*tdif(n,k-1))*(tdif(n,k-1)*rdenom)
!
        end if
      end do
!
      do n = n1 , n2
!
!       Compute next layer delta-Eddington solution only if total
!       transmission of radiation to the interface just above the layer
!       exceeds trmin.
!
        if ( tottrn(n,k) > trmin ) then
!
!         Remember, no ozone absorption in this layer:
!
          tauray(n) = trayoslp*pflx(n,kzp1)
          taugab(n) = abh2o(ns)*uth2o(n) + abco2(ns)*utco2(n) + abo2(ns)*uto2(n)
          tautot = tauray(n) + taugab(n) + tauxar(n,ns)
          wtot = (wray*tauray(n)+tauasc(n,ns))/tautot
          gtot = (gray*wray*tauray(n)+gtota(n,ns))/(wtot*tautot)
          ftot = (fray*wray*tauray(n)+ftota(n,ns))/(wtot*tautot)
!
          ts = taus(wtot,ftot,tautot)
          ws = omgs(wtot,ftot)
          gs = asys(gtot,ftot)
          lm = el(ws,gs)
          alp = xalpha(ws,czen(n),gs,lm)
          gam = xgamma(ws,czen(n),gs,lm)
          ue = f_u(ws,gs,lm)
!
!         Limit argument of exponential to 25, in case lm very large:
!
          arg = dmin1(lm*ts,25.0D0)
          extins = dexp(-arg)
          ne = f_n(ue,extins)
!
          rdif(n,k) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
          tdif(n,k) = d_four*ue/ne
!
!         Limit argument of exponential to 25, in case czen is very small:
          arg = dmin1(ts/czen(n),25.0D0)
          explay(n,k) = dexp(-arg)
!
          apg = alp + gam
          amg = alp - gam
          rdir(n,k) = amg*(tdif(n,k)*explay(n,k)-d_one)+apg*rdif(n,k)
          tdir(n,k) = apg*tdif(n,k) + (amg*rdif(n,k)-(apg-d_one))*explay(n,k)
!
!         Under rare conditions, reflectivies and transmissivities
!         can be negative; zero out any negative values
!
          rdir(n,k) = dmax1(rdir(n,k),d_zero)
          tdir(n,k) = dmax1(tdir(n,k),d_zero)
          rdif(n,k) = dmax1(rdif(n,k),d_zero)
          tdif(n,k) = dmax1(tdif(n,k),d_zero)
!
        end if
      end do
!
    end do
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for both layers
!   above the surface:
!
    k = 2
    do n = n1 , n2
      if ( czengt0(n) ) then
        exptdn(n,k) = exptdn(n,k-1)*explay(n,k-1)
        rdenom = d_one/(d_one-rdif(n,k-1)*rdndif(n,k-1))
        rdirexp = rdir(n,k-1)*exptdn(n,k-1)
        tdnmexp = tottrn(n,k-1) - exptdn(n,k-1)
        tottrn(n,k) = exptdn(n,k-1)*tdir(n,k-1) + &
                      tdif(n,k-1)*(tdnmexp+rdndif(n,k-1)*rdirexp)*rdenom
        rdndif(n,k) = rdif(n,k-1) + &
                      (rdndif(n,k-1)*tdif(n,k-1))*(tdif(n,k-1)*rdenom)
      end if
    end do
!
    call time_end(subroutine_name,indx)
  end subroutine radclr
!
!-----------------------------------------------------------------------
!
! Computes layer reflectivities and transmissivities, from the top down
! to the surface using the delta-Eddington solutions for each layer;
! adds layers from top down to surface as well.
!
! If total transmission to the interface above a particular layer is
! less than trmin, then no further delta-Eddington solutions are
! evaluated for layers below
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
  subroutine radded(n1,n2,trayoslp,czen,czengt0,tauxcl,tauxci,ns)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:) , intent(in) :: czen
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: tauxcl , tauxci
    logical , pointer , dimension(:) , intent(in) :: czengt0
    real(dp) :: trayoslp
    integer :: ns
    intent (in) trayoslp , ns
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! ii       - Longitude index
! nval     - Number of long values satisfying criteria
! indx     - Array of longitude indices
! taugab   - Layer total gas absorption optical depth
! tauray   - Layer rayleigh optical depth
! taucsc   - Layer cloud scattering optical depth
! tautot   - Total layer optical depth
! wtot     - Total layer single scatter albedo
! gtot     - Total layer asymmetry parameter
! ftot     - Total layer forward scatter fraction
! wtau     - rayleigh layer scattering optical depth
! wt       - layer total single scattering albedo
! ts       - layer scaled extinction optical depth
! ws       - layer scaled single scattering albedo
! gs       - layer scaled asymmetry parameter
! rdenom   - mulitiple scattering term
! rdirexp  - layer direct ref times exp transmission
! tdnmexp  - total transmission minus exp transmission
!
!---------------------------Statement functions-------------------------
!
!     Statement functions and other local variables
!
! alpha    - Term in direct reflect and transmissivity
! xgamm    - Term in direct reflect and transmissivity
! el       - Term in alpha,xgamm,n,u
! taus     - Scaled extinction optical depth
! omgs     - Scaled single particle scattering albedo
! asys     - Scaled asymmetry parameter
! u        - Term in diffuse reflect and transmissivity
! n        - Term in diffuse reflect and transmissivity
! lm       - Temporary for el
! ne       - Temporary for n
!
!     Intermediate terms for delta-eddington solution
!
! alp      - Temporary for alpha
! gam      - Temporary for xgamm
! ue       - Temporary for u
! arg      - Exponential argument
! extins   - Extinction
! amg      - Alp - gam
! apg      - Alp + gam
!
    real(dp) :: alp , amg , apg , arg , extins , ftot ,        &
               gam , gs , gtot , lm , ne , rdenom , rdirexp , &
               taucsc , tautot , tdnmexp , ts , ue , ws ,     &
               wt , wtau , wtot
    integer :: n , k
!
    character (len=64) :: subroutine_name='radded'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!-----------------------------------------------------------------------
!
!   Initialize all total transmission values to 0, so that nighttime
!   values from previous computations are not used:
!
    tottrn(:,:) = d_zero
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above each interface by starting from the top and adding layers down:
!   For the extra layer above model top:
!
    do n = n1 , n2
      if ( czengt0(n) ) then
!
        tauray(n) = trayoslp*(pflx(n,1)-pflx(n,0))
        taugab(n) = abh2o(ns)*uh2o(n,0) + abo3(ns)*uo3(n,0) + &
                    abco2(ns)*uco2(n,0) + abo2(ns)*uo2(n,0)
        tautot = tauxcl(n,0,ns)+tauxci(n,0,ns) + &
                 tauray(n)+taugab(n)+tauxar3d(n,0,ns)
        taucsc = tauxcl(n,0,ns)*wcl(n,0)+tauxci(n,0,ns)*wci(n,0) + &
                 tauasc3d(n,0,ns)
        wtau = wray*tauray(n)
        wt = wtau + taucsc
        wtot = wt/tautot
        gtot = (wtau*gray+gcl(n,0)*tauxcl(n,0,ns)*wcl(n,0)+gci(n,0) * &
                tauxci(n,0,ns)*wci(n,0)+gtota3d(n,0,ns))/wt
        ftot = (wtau*fray+fcl(n,0)*tauxcl(n,0,ns)*wcl(n,0)+fci(n,0) * &
                tauxci(n,0,ns)*wci(n,0)+ftota3d(n,0,ns))/wt
        ts = taus(wtot,ftot,tautot)
        ws = omgs(wtot,ftot)
        gs = asys(gtot,ftot)
        lm = el(ws,gs)
        alp = xalpha(ws,czen(n),gs,lm)
        gam = xgamma(ws,czen(n),gs,lm)
        ue = f_u(ws,gs,lm)
!
!       Limit argument of exponential to 25, in case lm*ts very large:
!
        arg = dmin1(lm*ts,25.0D0)
        extins = dexp(-arg)
        ne = f_n(ue,extins)
!
        rdif(n,0) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
        tdif(n,0) = d_four*ue/ne
!
!       Limit argument of exponential to 25, in case czen is very small:
        arg = dmin1(ts/czen(n),25.0D0)
        explay(n,0) = dexp(-arg)
!
        apg = alp + gam
        amg = alp - gam
        rdir(n,0) = amg*(tdif(n,0)*explay(n,0)-d_one) + apg*rdif(n,0)
        tdir(n,0) = apg*tdif(n,0) + (amg*rdif(n,0)-(apg-d_one))*explay(n,0)
!
!       Under rare conditions, reflectivies and transmissivities can
!       be negative; zero out any negative values
!
        rdir(n,0) = dmax1(rdir(n,0),d_zero)
        tdir(n,0) = dmax1(tdir(n,0),d_zero)
        rdif(n,0) = dmax1(rdif(n,0),d_zero)
        tdif(n,0) = dmax1(tdif(n,0),d_zero)
!
!       Initialize top interface of extra layer:
!
        exptdn(n,0) = d_one
        rdndif(n,0) = d_zero
        tottrn(n,0) = d_one
!
        rdndif(n,1) = rdif(n,0)
        tottrn(n,1) = tdir(n,0)
!
      end if
    end do
!
!   Now, continue down one layer at a time; if the total transmission
!   to the interface just above a given layer is less than trmin,
!   then no delta-eddington computation for that layer is done:
!
    do k = 1 , kz
!
!     Initialize current layer properties to zero; only if total
!     transmission to the top interface of the current layer exceeds
!     the minimum, will these values be computed below:
!
      do n = n1 , n2
        if ( czengt0(n) ) then
!
          rdir(n,k) = d_zero
          rdif(n,k) = d_zero
          tdir(n,k) = d_zero
          tdif(n,k) = d_zero
          explay(n,k) = d_zero
!
!         Calculates the solar beam transmission, total transmission,
!         and reflectivity for diffuse mod_radiation from below at the
!         top of the current layer:
!
          exptdn(n,k) = exptdn(n,k-1)*explay(n,k-1)
!KN       modified below (for computational stability)
!         rdenom      = d_one/(1. - rdif(n,k-1)*rdndif(n,k-1))
          rdenom = d_one/(d_one-dmin1(rdif(n,k-1)*rdndif(n,k-1),verynearone))
!KN       modified above
          rdirexp = rdir(n,k-1)*exptdn(n,k-1)
          tdnmexp = tottrn(n,k-1) - exptdn(n,k-1)
          tottrn(n,k) = exptdn(n,k-1)*tdir(n,k-1) + tdif(n,k-1) *     &
                        (tdnmexp+rdndif(n,k-1)*rdirexp)*rdenom
          rdndif(n,k) = rdif(n,k-1) + (rdndif(n,k-1)*tdif(n,k-1)) *   &
                        (tdif(n,k-1)*rdenom)
!
        end if
      end do
!
      do n = n1 , n2
!
!       Compute next layer delta-eddington solution only if total
!       transmission of radiation to the interface just above the layer
!       exceeds trmin.
!
        if ( tottrn(n,k) > trmin ) then
!
          tauray(n) = trayoslp*(pflx(n,k+1)-pflx(n,k))
          taugab(n) = abh2o(ns)*uh2o(n,k) + abo3(ns)*uo3(n,k) +  &
                      abco2(ns)*uco2(n,k) + abo2(ns)*uo2(n,k)
!
          tautot = tauxcl(n,k,ns) + tauxci(n,k,ns) + tauray(n) + &
                   taugab(n) + tauxar3d(n,k,ns)
          taucsc = tauxcl(n,k,ns)*wcl(n,k) + tauxci(n,k,ns)*wci(n,k) + &
                   tauasc3d(n,k,ns)
          wtau = wray*tauray(n)
          wt = wtau + taucsc
          wtot = wt/tautot
          gtot = (wtau*gray+gcl(n,k)*wcl(n,k)*tauxcl(n,k,ns)+gci(n,k) *  &
                  wci(n,k)*tauxci(n,k,ns)+gtota3d(n,k,ns))/wt
          ftot = (wtau*fray+fcl(n,k)*wcl(n,k)*tauxcl(n,k,ns)+fci(n,k) *  &
                  wci(n,k)*tauxci(n,k,ns)+ftota3d(n,k,ns))/wt
!
          ts = taus(wtot,ftot,tautot)
          ws = omgs(wtot,ftot)
          gs = asys(gtot,ftot)
          lm = el(ws,gs)
          alp = xalpha(ws,czen(n),gs,lm)
          gam = xgamma(ws,czen(n),gs,lm)
          ue = f_u(ws,gs,lm)
!
!         Limit argument of exponential to 25, in case lm very large:
!
          arg = dmin1(lm*ts,25.0D0)
          extins = dexp(-arg)
          ne = f_n(ue,extins)
!
          rdif(n,k) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
          tdif(n,k) = d_four*ue/ne
!
!         Limit argument of exponential to 25, in case czen is very small:
          arg = dmin1(ts/czen(n),25.0D0)
          explay(n,k) = dexp(-arg)
!
          apg = alp + gam
          amg = alp - gam
          rdir(n,k) = amg*(tdif(n,k)*explay(n,k)-d_one)+apg*rdif(n,k)
          tdir(n,k) = apg*tdif(n,k) + (amg*rdif(n,k)-(apg-d_one))*explay(n,k)
!
!         Under rare conditions, reflectivies and transmissivities
!         can be negative; zero out any negative values
!
          rdir(n,k) = dmax1(rdir(n,k),d_zero)
          tdir(n,k) = dmax1(tdir(n,k),d_zero)
          rdif(n,k) = dmax1(rdif(n,k),d_zero)
          tdif(n,k) = dmax1(tdif(n,k),d_zero)
        end if
      end do
!
    end do
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above the surface:
!
    k = kzp1
    do n = n1 , n2
      if ( czengt0(n) ) then
        exptdn(n,k) = exptdn(n,k-1)*explay(n,k-1)
!KN     modified below (for computational stability)
!       rdenom = d_one/(1. - rdif(n,k-1)*rdndif(n,k-1))
        rdenom = d_one/(d_one-dmin1(rdif(n,k-1)*rdndif(n,k-1),verynearone))
!KN     modified above
        rdirexp = rdir(n,k-1)*exptdn(n,k-1)
        tdnmexp = tottrn(n,k-1) - exptdn(n,k-1)
        tottrn(n,k) = exptdn(n,k-1)*tdir(n,k-1) + tdif(n,k-1) *       &
                      (tdnmexp+rdndif(n,k-1)*rdirexp)*rdenom
        rdndif(n,k) = rdif(n,k-1) + (rdndif(n,k-1)*tdif(n,k-1)) *     &
                      (tdif(n,k-1)*rdenom)
      end if
    end do
!
    call time_end(subroutine_name,indx)
  end subroutine radded
!
!-----------------------------------------------------------------------
!
! Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
!
! h2o  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
! co2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019.
!            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
!            are also included.
!
! o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R.Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of CFC11
!
! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of CFC12
!
!
! Computes individual absorptivities for non-adjacent layers, accounting
! for band overlap, and sums to obtain the total; then, computes the
! nearest layer contribution.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Briegleb, May 1996
!
!-----------------------------------------------------------------------
!
!------------------------------Arguments--------------------------------
!
!     Trace gas variables
!
! abstrc    - total trace gas absorptivity
! bplnk     - Planck functions for sub-divided layers
!
!     Output arguments (radbuf)
!
! Nearest layer absorptivities
! Non-adjacent layer absorptivites
! Total emissivity
!
  subroutine radabs(n1,n2,pint,pmid,piln,pmln,absgasnxt,absgastot)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:,:) , intent(in) :: pint , pmid , &
      piln , pmln
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: absgasnxt
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: absgastot
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! k1       - Level index
! k2       - Level index
! kn       - Nearest level index
! iband    - Band  index
! pnew     - Effective pressure for H2O vapor linewidth
! trline   - Transmission due to H2O lines in window
! ux       - Pressure weighted H2O path length
! tbar     - Mean layer temperature
! emm      - Mean co2 emissivity
! o3emm    - Mean o3 emissivity
! o3bndi   - Ozone band parameter
! temh2o   - Mean layer temperature equivalent to tbar
! k21      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 650-800 cm-1 region (tr1)
! k22      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 500-650 cm-1 region (tr2)
! uc1      - H2o continuum pathlength in 500-800 cm-1
! to3h2o   - H2o trnsmsn for overlap with o3
! pi       - For co2 absorptivity computation
! sqti     - Used to store sqrt of mean temperature
! et       - Co2 hot band factor
! et2      - Co2 hot band factor squared
! et4      - Co2 hot band factor to fourth power
! omet     - Co2 stimulated emission term
! f1co2    - Co2 central band factor
! f2co2    - Co2 weak band factor
! f3co2    - Co2 weak band factor
! t1co2    - Overlap factr weak bands on strong band
! sqwp     - Sqrt of co2 pathlength
! f1sqwp   - Main co2 band factor
! oneme    - Co2 stimulated emission term
! alphat   - Part of the co2 stimulated emission term
! wco2     - Constants used to define co2 pathlength
! posqt    - Effective pressure for co2 line width
! u7       - Co2 hot band path length
! u8       - Co2 hot band path length
! u9       - Co2 hot band path length
! u13      - Co2 hot band path length
! rbeta7   - Inverse of co2 hot band line width par
! rbeta8   - Inverse of co2 hot band line width par
! rbeta9   - Inverse of co2 hot band line width par
! rbeta13  - Inverse of co2 hot band line width par
! tpatha   - For absorptivity computation
! a        - Eq(2) in table A3a of R&D
! abso     - Absorptivity for various gases/bands
! dtp      - Path temp minus 300 K used in h2o rotation band absorptivity
! dtx      - Planck temperature minus 250 K
! dty      - Path temperature minus 250 K
! dtz      - Planck temperature minus 300 K
! term1    - Equation(5) in table A3a of R&D(1986)
! term2    - Delta a(Te) in table A3a of R&D(1986)
! term3    - DB/dT function for rotation and
! term4    - Equation(6) in table A3a of R&D(1986)
! term5    - Delta a(Tp) in table A3a of R&D(1986)
! term6    - DB/dT function for window region
! term7    - Kl_inf(i) in eq(8) of table A3a of R&D
! term8    - Delta kl_inf(i) in eq(8)
! term9    - DB/dT function for 500-800 cm-1 region
! tr1      - Eqn(6) in table A2 of R&D for 650-800
! tr10     - Eqn (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! tr2      - Eqn(6) in table A2 of R&D for 500-650
! tr5      - Eqn(4) in table A2 of R&D for 650-800
! tr6      - Eqn(4) in table A2 of R&D for 500-650
! tr9      - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! uc       - Y + 0.002U in eq(8) of table A2 of R&D
! sqrtu    - Sqrt of pressure weighted h20 pathlength
! fwk      - Equation(33) in R&D far wing correction
! fwku     - GU term in eqs(1) and (6) in table A2
! r2st     - 1/(2*beta) in eq(10) in table A2
! dtyp15   - DeltaTp in eqs(11) & (12) in table A3a
! dtyp15sq - (DeltaTp)^2 in eqs(11) & (12) table A3a
! to3co2   - P weighted temp in ozone band model
! dpnm     - Pressure difference between two levels
! pnmsq    - Pressure squared
! dw       - Amount of h2o between two levels
! uinpl    - Nearest layer subdivision factor
! winpl    - Nearest layer subdivision factor
! zinpl    - Nearest layer subdivision factor
! pinpl    - Nearest layer subdivision factor
! dplh2o   - Difference in press weighted h2o amount
! ds2c     - Y in eq(7) in table A2 of R&D
! a11      - A1 in table A3b for rotation band absorptivity
! a31      - A3 in table A3b for rotation band absorptivity
! a21      - First part in numerator of A2 in table A3b
! a22      - Second part in numerator of A2 in table A3b
! a23      - Denominator of A2 in table A3b (rotation band)
! t1t4     - Eq(3) in table A3a of R&D
! t2t5     - Eq(4) in table A3a of R&D
! rsum     - Eq(1) in table A2 of R&D
! a41      - Numerator in A2 in Vib-rot abstivity(table A3b)
! a51      - Denominator in A2 in Vib-rot (table A3b)
! a61      - A3 factor for Vib-rot band in table A3b
! phi      - Eq(11) in table A3a of R&D
! psi      - Eq(12) in table A3a of R&D
! cf812    - Eq(11) in table A2 of R&D
! ubar     - H2o scaled path see comment for eq(10) table A2
! pbar     - H2o scaled pres see comment for eq(10) table A2
! g4       - Arguement in exp() in eq(10) table A2
! dplos    - Ozone pathlength eq(A2) in R&Di
! dplol    - Presure weighted ozone pathlength
! beta     - Local interface temperature (includes Voigt line correction factor)
! rphat    - Effective pressure for ozone beta
! tcrfac   - Ozone temperature factor table 1 R&Di
! tmp1     - Ozone band factor see eq(A1) in R&Di
! u1       - Effective ozone pathlength eq(A2) in R&Di
! realnu   - 1/beta factor in ozone band model eq(A1)
! tmp2     - Ozone band factor see eq(A1) in R&Di
! u2       - Effective ozone pathlength eq(A2) in R&Di
! rsqti    - Reciprocal of sqrt of path temperature
! tpath    - Path temperature used in co2 band model
! tmp3     - Weak band factor see K&B
! rdpnmsq  - Reciprocal of difference in press^2
! rdpnm    - Reciprocal of difference in press
! p1       - Mean pressure factor
! p2       - Mean pressure factor
! dtym10   - T - 260 used in eq(9) and (10) table A3a
! dplco2   - Co2 pathlength
! corfac   - Correction factors in table A3b
! g2       - Part of arguement in eq(10) in table A2
! te       - A_0 T factor in ozone model table 1 of R&Di
! denom    - Denominator in eq(8) of table A3a of R&D
! trab2    - Transmission terms for H2o  500 -  800 cm-1
! trab4    - Transmission terms for H2o  800 - 1000 cm-1
! trab6    - Transmission terms for H2o 1000 - 1200 cm-1
! absbnd   - Proportional to co2 band absorptance
! dbvtit   - Intrfc drvtv plnck fnctn for o3
! dbvt     - Planck fnctn tmp derivative for o3
!
    real(dp) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,    &
               absbnd , alphat , beta , cf812 , corfac , denom ,      &
               dplco2 , dplol , dplos , ds2c , dtym10 , et , et2 ,    &
               et4 , f1co2 , g2 , g4 , k21 , k22 , o3bndi , omet ,    &
               oneme , p1 , p2 , pbar , phi , pi , posqt , psi ,      &
               rbeta13 , rbeta8 , rbeta9 , rdpnm , rdpnmsq , realnu , &
               rphat , rsqti , rsum , sqwp , t1t4 , t2t5 ,            &
               tcrfac , te , tlocal , tmp1 , tmp2 , tmp3 , tpath ,    &
               tr1 , tr2 , tr5 , tr6 , tr9 , tr10 , u1 , u13 , u2 ,   &
               u8 , u9 , ubar , wco2 , dplh2o , dtp , dtz , sqti ,    &
               dpnm , dtyp15 , dtyp15sq , f1sqwp , f2co2 , f3co2 ,    &
               fwk , fwku , rbeta7 , sqrtu , t1co2 , to3h2o ,         &
               tpatha , trab2 , trab4 , trab6 , u7 , uc1 , uc
    real(dp) , dimension(6) :: abso
    real(dp) , dimension(4) :: emm , o3emm , term1 , term2 , &
               term3 , term4 , term5 , zinpl , temh2o
    real(dp) , dimension(2) :: r2st , term7 , term8 , trline
    integer :: n , iband , k , k1 , k2 , kn , wvl
!
    character (len=64) :: subroutine_name='radabs'
    integer :: indx = 0
!
!-----------------------------------------------------------------------
    call time_begin(subroutine_name,indx)
!
!   Initialize
!
    do k = 1 , kz
      do n = n1 , n2
        dbvtly(n,k) = dbvt(tlayr(n,k+1))
        dbvtit(n,k) = dbvt(tint(n,k))
      end do
    end do
    do n = n1 , n2
      dbvtit(n,kzp1) = dbvt(tint(n,kz + 1))
    end do
!
    r2st(1) = d_one/(d_two*st(1))
    r2st(2) = d_one/(d_two*st(2))
!   bndfct  = 2.0*22.18d0/(dsqrt(196.d0)*300.)
!
!   Non-adjacent layer absorptivity:
!
!   abso(1)     0 -  800 cm-1   h2o rotation band
!   abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
!   abso(3)   800 - 1200 cm-1   h2o window
!   abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
!   abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!   abso(6)   co2 15  micrometer band system
!
    do k = 1 , kzp1
      do n = n1 , n2
        pnmsq(n,k) = pint(n,k)**d_two
        dtx(n) = tplnka(n,k) - 250.0D0
        term6(n,k) = coeff(1,2) + coeff(2,2)*dtx(n) *          &
                      (d_one+c9*dtx(n)*(d_one+c11*dtx(n) *     &
                      (d_one+c13*dtx(n)*(d_one+c15*dtx(n)))))
        term9(n,k) = coefi(1,2) + coefi(2,2)*dtx(n) *          &
                      (d_one+c19*dtx(n)*(d_one+c21*dtx(n) *    &
                      (d_one+c23*dtx(n)*(d_one+c25*dtx(n)))))
      end do
    end do
!
!   Non-nearest layer level loops
!
    do k1 = kzp1 , 1 , -1
      do k2 = kzp1 , 1 , -1
        if ( k1 /= k2 ) then
          do n = n1 , n2
            dplh2o = plh2o(n,k1) - plh2o(n,k2)
            ux(n) = dabs(dplh2o)
            sqrtu = dsqrt(ux(n))
            ds2c = dabs(s2c(n,k1)-s2c(n,k2))
            dw(n) = dabs(w(n,k1)-w(n,k2))
            uc1 = (ds2c+1.7D-3*ux(n))*(d_one+d_two*ds2c)/(d_one+15.0D0*ds2c)
            uc = ds2c + 2.0D-3*ux(n)
            pnew(n) = ux(n)/dw(n)
            tpatha = (s2t(n,k1)-s2t(n,k2))/dplh2o
            dtx(n) = tplnka(n,k2) - 250.0D0
            dty(n) = tpatha - 250.0D0
            dtyp15 = dty(n) + 15.0D0
            dtyp15sq = dtyp15**d_two
            dtz = dtx(n) - 50.0D0
            dtp = dty(n) - 50.0D0
            do iband = 2 , 4 , 2
              term1(iband) = coefe(1,iband) + &
                             coefe(2,iband)*dtx(n)*(d_one+c1(iband)*dtx(n))
              term2(iband) = coefb(1,iband) + &
                             coefb(2,iband)*dtx(n)*(d_one+c2(iband)*dtx(n) * &
                             (d_one+c3(iband)*dtx(n)))
              term3(iband) = coefd(1,iband) + &
                             coefd(2,iband)*dtx(n)*(d_one+c4(iband)*dtx(n) * &
                             (d_one+c5(iband)*dtx(n)))
              term4(iband) = coefa(1,iband) + &
                             coefa(2,iband)*dty(n)*(d_one+c6(iband)*dty(n))
              term5(iband) = coefc(1,iband) + &
                             coefc(2,iband)*dty(n)*(d_one+c7(iband)*dty(n))
            end do
!
!           abso(1)     0 -  800 cm-1   h2o rotation band
!
            a11 = 0.44D0 + 3.380D-4*dtz - 1.520D-6*dtz*dtz
            a31 = 1.05D0 - 6.000D-3*dtp + 3.000D-6*dtp*dtp
            a21 = 1.00D0 + 1.717D-3*dtz - 1.133D-5*dtz*dtz
            a22 = 1.00D0 + 4.443D-3*dtp + 2.750D-5*dtp*dtp
            a23 = 1.00D0 + 3.600D0*sqrtu
            corfac = a31*(a11+((d_two*a21*a22)/a23))
            t1t4 = term1(2)*term4(2)
            t2t5 = term2(2)*term5(2)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
            fwk = fwcoef + fwc1/(d_one+fwc2*ux(n))
            fwku = fwk*ux(n)
            rsum = dexp(-a*(sqrtu+fwku))
            abso(1) = (d_one-rsum)*term3(2)
!           trab1(n)  = rsum
!
!           abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
            a41 = 1.75D0 - 3.960D-03*dtz
            a51 = 1.00D0 + 1.3D0*sqrtu
            a61 = 1.00D0 + 1.250D-03*dtp + 6.250D-05*dtp*dtp
            corfac = 0.29D0*(d_one+a41/a51)*a61
            t1t4 = term1(4)*term4(4)
            t2t5 = term2(4)*term5(4)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
            rsum = dexp(-a*(sqrtu+fwku))
            abso(2) = (d_one-rsum)*term3(4)
!           trab7(n)  = rsum
!
!         Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
            do k = 1 , 2
              phi = dexp(a1(k)*dtyp15+a2(k)*dtyp15sq)
              psi = dexp(b1(k)*dtyp15+b2(k)*dtyp15sq)
              ubar = dw(n)*phi*1.66D0*r80257
              pbar = pnew(n)*(psi/phi)
              cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
              g2 = d_one + ubar*d_four*st(k)*cf812/pbar
              g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-d_one)
              trline(k) = dexp(-g4)
            end do
            term7(1) = coefj(1,1)+coefj(2,1)*dty(n)*(d_one+c16*dty(n))
            term8(1) = coefk(1,1)+coefk(2,1)*dty(n)*(d_one+c17*dty(n))
            term7(2) = coefj(1,2)+coefj(2,2)*dty(n)*(d_one+c26*dty(n))
            term8(2) = coefk(1,2)+coefk(2,2)*dty(n)*(d_one+c27*dty(n))
!
!           abso(3)   800 - 1200 cm-1   h2o window
!           abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
            k21 = term7(1) + term8(1)/(d_one+(c30+c31*(dty(n)-d_10)* &
                  (dty(n)-d_10))*sqrtu)
            k22 = term7(2) + term8(2)/(d_one+(c28+c29*(dty(n)-d_10))*sqrtu)
            tr1 = dexp(-(k21*(sqrtu+fc1*fwku)))
            tr2 = dexp(-(k22*(sqrtu+fc1*fwku)))
            tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(n))*uc1))
            tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(n))*uc1))
            tr9 = tr1*tr5
            tr10 = tr2*tr6
            th2o(n) = tr10
            trab2 = 0.65D0*tr9 + 0.35D0*tr10
            trab4 = dexp(-(coefg(1,3)+coefg(2,3)*dtx(n))*uc)
            trab6 = dexp(-(coefg(1,4)+coefg(2,4)*dtx(n))*uc)
            abso(3) = term6(n,k2)*(d_one-trab4*d_half*trline(2)- &
                        trab6*d_half*trline(1))
            abso(4) = term9(n,k2)*d_half*(tr1-tr9+tr2-tr10)
            if ( k2 < k1 ) then
              to3h2o = h2otr(n,k1)/h2otr(n,k2)
            else
              to3h2o = h2otr(n,k2)/h2otr(n,k1)
            end if
!
!           abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
            dpnm = pint(n,k1) - pint(n,k2)
            to3co2(n) = (pint(n,k1)*co2t(n,k1)-pint(n,k2)*co2t(n,k2))/dpnm
            te = (to3co2(n)*r293)**0.7D0
            dplos = plos(n,k1) - plos(n,k2)
            dplol = plol(n,k1) - plol(n,k2)
            u1 = 18.29D0*dabs(dplos)/te
            u2 = 0.5649D0*dabs(dplos)/te
            rphat = dplol/dplos
            tlocal = tint(n,k2)
            tcrfac = dsqrt(tlocal*r250)*te
            beta = r3205*(rphat+dpfo3*tcrfac)
            realnu = te/beta
            tmp1 = u1/dsqrt(d_four+u1*(d_one+realnu))
            tmp2 = u2/dsqrt(d_four+u2*(d_one+realnu))
            o3bndi = 74.0D0*te*dlog(d_one+tmp1+tmp2)
            abso(5) = o3bndi*to3h2o*dbvtit(n,k2)
            to3(n) = d_one/(d_one+0.1D0*tmp1+0.1D0*tmp2)
!           trab5(n)  = d_one-(o3bndi/(1060-980.))
!
!           abso(6)      co2 15  micrometer band system
!
            sqwp = dsqrt(dabs(plco2(n,k1)-plco2(n,k2)))
            et = dexp(-480.0D0/to3co2(n))
            sqti = dsqrt(to3co2(n))
            rsqti = d_one/sqti
            et2 = et*et
            et4 = et2*et2
            omet = d_one - 1.5D0*et2
            f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
            f1sqwp = f1co2*sqwp
            t1co2 = d_one/(d_one+(245.18D0*omet*sqwp*rsqti))
            oneme = d_one - et2
            alphat = oneme**3.0D0*rsqti
            pi = dabs(dpnm)
            wco2 = 2.5221D0*co2vmr*pi*regravgts
            u7 = 4.9411D4*alphat*et2*wco2
            u8 = 3.9744D4*alphat*et4*wco2
            u9 = 1.0447D5*alphat*et4*et2*wco2
            u13 = 2.8388D3*alphat*et4*wco2
            tpath = to3co2(n)
            tlocal = tint(n,k2)
            tcrfac = dsqrt(tlocal*r250*tpath*r300)
            posqt = ((pint(n,k2)+pint(n,k1))*r2sslp+dpfco2*tcrfac)*rsqti
            rbeta7 = d_one/(5.3228D0*posqt)
            rbeta8 = d_one/(10.6576D0*posqt)
            rbeta9 = rbeta7
            rbeta13 = rbeta9
            f2co2 = (u7/dsqrt(d_four+u7*(d_one+rbeta7))) + &
                    (u8/dsqrt(d_four+u8*(d_one+rbeta8))) + &
                    (u9/dsqrt(d_four+u9*(d_one+rbeta9)))
            f3co2 = u13/dsqrt(d_four+u13*(d_one+rbeta13))
            if ( k2 >= k1 ) then
              sqti = dsqrt(tlayr(n,k2))
            end if
!
            tmp1 = dlog(d_one+f1sqwp)
            tmp2 = dlog(d_one+f2co2)
            tmp3 = dlog(d_one+f3co2)
            absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
            abso(6) = trab2*co2em(n,k2)*absbnd
            tco2(n) = d_one/(d_one+d_10*(u7/dsqrt(d_four+u7*(d_one+rbeta7))))
!           trab3(n)  = 1. - bndfct*absbnd
            absgastot(n,k1,k2) = abso(1) + abso(2) + abso(3) + &
                                 abso(4) + abso(5) + abso(6)
          end do
!
!         Calculate absorptivity due to trace gases
!
          call trcab(n1,n2,k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4,    &
                     uco211,uco212,uco213,uco221,uco222,uco223,bn2o0,     &
                     bn2o1,bch4,to3co2,pint,dw,pnew,s2c,uptype,ux,abplnk1, &
                     tco2,th2o,to3,abstrc)
!
!         Sum total absorptivity
!
          do n = n1 , n2
            absgastot(n,k1,k2) = absgastot(n,k1,k2) + abstrc(n)
          end do
        end if
      end do
    end do  ! End of non-nearest layer level loops
!
!   Non-adjacent layer absorptivity:
!
!   abso(1)     0 -  800 cm-1   h2o rotation band
!   abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
!   abso(3)   800 - 1200 cm-1   h2o window
!   abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
!   abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!   abso(6)   co2 15  micrometer band system
!
!   Nearest layer level loop
!
    do k2 = kz , 1 , -1
      do n = n1 , n2
        tbar(n,1) = (tint(n,k2+1)+tlayr(n,k2+1))*d_half
        emm(1) = (co2em(n,k2+1)+co2eml(n,k2))*d_half
        tbar(n,2) = (tlayr(n,k2+1)+tint(n,k2))*d_half
        emm(2) = (co2em(n,k2)+co2eml(n,k2))*d_half
        tbar(n,3) = (tbar(n,2)+tbar(n,1))*d_half
        emm(3) = emm(1)
        tbar(n,4) = tbar(n,3)
        emm(4) = emm(2)
        o3emm(1) = (dbvtit(n,k2+1)+dbvtly(n,k2))*d_half
        o3emm(2) = (dbvtit(n,k2)+dbvtly(n,k2))*d_half
        o3emm(3) = o3emm(1)
        o3emm(4) = o3emm(2)
        temh2o(1) = tbar(n,1)
        temh2o(2) = tbar(n,2)
        temh2o(3) = tbar(n,1)
        temh2o(4) = tbar(n,2)
        dpnm = pint(n,k2+1) - pint(n,k2)
!
!     Weighted Planck functions for trace gases
!
        do wvl = 1 , 14
          bplnk(wvl,n,1) = (abplnk1(wvl,n,k2+1)+abplnk2(wvl,n,k2))*d_half
          bplnk(wvl,n,2) = (abplnk1(wvl,n,k2)+abplnk2(wvl,n,k2))*d_half
          bplnk(wvl,n,3) = bplnk(wvl,n,1)
          bplnk(wvl,n,4) = bplnk(wvl,n,2)
        end do
!
        rdpnmsq = d_one/(pnmsq(n,k2+1)-pnmsq(n,k2))
        rdpnm = d_one/dpnm
        p1 = (pmid(n,k2)+pint(n,k2+1))*d_half
        p2 = (pmid(n,k2)+pint(n,k2))*d_half
        uinpl(n,1) = (pnmsq(n,k2+1)-p1**d_two)*rdpnmsq
        uinpl(n,2) = -(pnmsq(n,k2)-p2**d_two)*rdpnmsq
        uinpl(n,3) = -(pnmsq(n,k2)-p1**d_two)*rdpnmsq
        uinpl(n,4) = (pnmsq(n,k2+1)-p2**d_two)*rdpnmsq
        winpl(n,1) = ((pint(n,k2+1)-pmid(n,k2))*d_half)*rdpnm
        winpl(n,2) = ((-pint(n,k2)+pmid(n,k2))*d_half)*rdpnm
        winpl(n,3) = ((pint(n,k2+1)+pmid(n,k2))*d_half-pint(n,k2))*rdpnm
        winpl(n,4) = ((-pint(n,k2)-pmid(n,k2))*d_half+pint(n,k2+1))*rdpnm
        tmp1 = d_one/(piln(n,k2+1)-piln(n,k2))
        tmp2 = piln(n,k2+1) - pmln(n,k2)
        tmp3 = piln(n,k2) - pmln(n,k2)
        zinpl(1) = (tmp2*d_half)*tmp1
        zinpl(2) = (-tmp3*d_half)*tmp1
        zinpl(3) = (tmp2*d_half-tmp3)*tmp1
        zinpl(4) = (tmp2-tmp3*d_half)*tmp1
        pinpl(n,1) = (p1+pint(n,k2+1))*d_half
        pinpl(n,2) = (p2+pint(n,k2))*d_half
        pinpl(n,3) = (p1+pint(n,k2))*d_half
        pinpl(n,4) = (p2+pint(n,k2+1))*d_half

!       FAB AER SAVE uinpl  for aerosl LW forcing calculation

        if ( lchem .and. idirect > 0 ) then
          do kn = 1 , 4
            xuinpl (n,k2,kn) =  uinpl(n,kn)
          end do
        end if
!       FAB AER SAVE uinpl  for aerosl LW forcing calculation

        do kn = 1 , 4
          ux(n) = uinpl(n,kn)*dabs(plh2o(n,k2)-plh2o(n,k2+1))
          sqrtu = dsqrt(ux(n))
          dw(n) = dabs(w(n,k2)-w(n,k2+1))
          pnew(n) = ux(n)/(winpl(n,kn)*dw(n))
          ds2c = dabs(s2c(n,k2)-s2c(n,k2+1))
          uc1 = uinpl(n,kn)*ds2c
          uc1 = (uc1+1.7D-3*ux(n))*(d_one+d_two*uc1)/(d_one+15.0D0*uc1)
          uc = uinpl(n,kn)*ds2c + 2.0D-3*ux(n)
          dtx(n) = temh2o(kn) - 250.0D0
          dty(n) = tbar(n,kn) - 250.0D0
          dtyp15 = dty(n) + 15.0D0
          dtyp15sq = dtyp15**d_two
          dtz = dtx(n) - 50.0D0
          dtp = dty(n) - 50.0D0
          do iband = 2 , 4 , 2
            term1(iband) = coefe(1,iband) + coefe(2,iband)*dtx(n) * &
                           (d_one+c1(iband)*dtx(n))
            term2(iband) = coefb(1,iband) + coefb(2,iband)*dtx(n) * &
                           (d_one+c2(iband)*dtx(n)                * &
                           (d_one+c3(iband)*dtx(n)))
            term3(iband) = coefd(1,iband) + coefd(2,iband)*dtx(n) * &
                           (d_one+c4(iband)*dtx(n)                * &
                           (d_one+c5(iband)*dtx(n)))
            term4(iband) = coefa(1,iband) + coefa(2,iband)*dty(n) * &
                           (d_one+c6(iband)*dty(n))
            term5(iband) = coefc(1,iband) + coefc(2,iband)*dty(n) * &
                           (d_one+c7(iband)*dty(n))
          end do
!
!         abso(1)     0 -  800 cm-1   h2o rotation band
!
          a11 = 0.44D0 + 3.380D-4*dtz - 1.520D-6*dtz*dtz
          a31 = 1.05D0 - 6.000D-3*dtp + 3.000D-6*dtp*dtp
          a21 = 1.00D0 + 1.717D-3*dtz - 1.133D-5*dtz*dtz
          a22 = 1.00D0 + 4.443D-3*dtp + 2.750D-5*dtp*dtp
          a23 = 1.00D0 + 3.600D0*sqrtu
          corfac = a31*(a11+((d_two*a21*a22)/a23))
          t1t4 = term1(2)*term4(2)
          t2t5 = term2(2)*term5(2)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
          fwk = fwcoef + fwc1/(d_one+fwc2*ux(n))
          fwku = fwk*ux(n)
          rsum = dexp(-a*(sqrtu+fwku))
          abso(1) = (d_one-rsum)*term3(2)
!         trab1(n) = rsum
!
!         abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
          a41 = 1.75D0 - 3.960D-03*dtz
          a51 = 1.00D0 + 1.3D0*sqrtu
          a61 = 1.00D0 + 1.250D-03*dtp + 6.250D-05*dtp*dtp
          corfac = 0.29D0*(d_one+a41/a51)*a61
          t1t4 = term1(4)*term4(4)
          t2t5 = term2(4)*term5(4)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
          rsum = dexp(-a*(sqrtu+fwku))
          abso(2) = (d_one-rsum)*term3(4)
!         trab7(n) = rsum
!
!       Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
          do k = 1 , 2
            phi = dexp(a1(k)*dtyp15+a2(k)*dtyp15sq)
            psi = dexp(b1(k)*dtyp15+b2(k)*dtyp15sq)
            ubar = dw(n)*phi*winpl(n,kn)*1.66D0*r80257
            pbar = pnew(n)*(psi/phi)
            cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
            g2 = d_one + ubar*d_four*st(k)*cf812/pbar
            g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-d_one)
            trline(k) = dexp(-g4)
          end do
          term7(1) = coefj(1,1)+coefj(2,1)*dty(n)*(d_one+c16*dty(n))
          term8(1) = coefk(1,1)+coefk(2,1)*dty(n)*(d_one+c17*dty(n))
          term7(2) = coefj(1,2)+coefj(2,2)*dty(n)*(d_one+c26*dty(n))
          term8(2) = coefk(1,2)+coefk(2,2)*dty(n)*(d_one+c27*dty(n))
!
!         abso(3)   800 - 1200 cm-1   h2o window
!         abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
!
          dtym10 = dty(n) - d_10
          denom = d_one + (c30+c31*dtym10*dtym10)*sqrtu
          k21 = term7(1) + term8(1)/denom
          denom = d_one + (c28+c29*dtym10)*sqrtu
          k22 = term7(2) + term8(2)/denom
          term9(n,2) = coefi(1,2) + coefi(2,2)*dtx(n) *        &
                       (d_one+c19*dtx(n)*(d_one+c21*dtx(n) *   &
                       (d_one+c23*dtx(n)*(d_one+c25*dtx(n)))))
          tr1 = dexp(-(k21*(sqrtu+fc1*fwku)))
          tr2 = dexp(-(k22*(sqrtu+fc1*fwku)))
          tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(n))*uc1))
          tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(n))*uc1))
          tr9 = tr1*tr5
          tr10 = tr2*tr6
          trab2 = 0.65D0*tr9 + 0.35D0*tr10
          th2o(n) = tr10
          trab4 = dexp(-(coefg(1,3)+coefg(2,3)*dtx(n))*uc)
          trab6 = dexp(-(coefg(1,4)+coefg(2,4)*dtx(n))*uc)
          term6(n,2) = coeff(1,2) + coeff(2,2)*dtx(n) *     &
                       (d_one+c9*dtx(n)*(d_one+c11*dtx(n) * &
                       (d_one+c13*dtx(n)*(d_one+c15*dtx(n)))))
          abso(3) = term6(n,2)*(d_one-trab4*d_half*trline(2) - &
                                      trab6*d_half*trline(1))
          abso(4) = term9(n,2)*d_half*(tr1-tr9+tr2-tr10)
!
!         abso(5)  o3  9.6 micrometer (nu3 and nu1 bands)
!
          te = (tbar(n,kn)*r293)**0.7D0
          dplos = dabs(plos(n,k2+1)-plos(n,k2))
          u1 = zinpl(kn)*18.29D0*dplos/te
          u2 = zinpl(kn)*0.5649D0*dplos/te
          tlocal = tbar(n,kn)
          tcrfac = dsqrt(tlocal*r250)*te
          beta = r3205*(pinpl(n,kn)*rsslp+dpfo3*tcrfac)
          realnu = te/beta
          tmp1 = u1/dsqrt(d_four+u1*(d_one+realnu))
          tmp2 = u2/dsqrt(d_four+u2*(d_one+realnu))
          o3bndi = 74.0D0*te*dlog(d_one+tmp1+tmp2)
          abso(5) = o3bndi*o3emm(kn)*(h2otr(n,k2+1)/h2otr(n,k2))
          to3(n) = d_one/(d_one+0.1D0*tmp1+0.1D0*tmp2)
!         trab5(n) = d_one-(o3bndi/(1060-980.))
!
!         abso(6)   co2 15  micrometer band system
!
          dplco2 = plco2(n,k2+1) - plco2(n,k2)
          sqwp = dsqrt(uinpl(n,kn)*dplco2)
          et = dexp(-480.0D0/tbar(n,kn))
          sqti = dsqrt(tbar(n,kn))
          rsqti = d_one/sqti
          et2 = et*et
          et4 = et2*et2
          omet = (d_one-1.5D0*et2)
          f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
          f1sqwp = f1co2*sqwp
          t1co2 = d_one/(d_one+(245.18D0*omet*sqwp*rsqti))
          oneme = d_one - et2
          alphat = oneme**3.0D0*rsqti
          pi = dabs(dpnm)*winpl(n,kn)
          wco2 = 2.5221D0*co2vmr*pi*regravgts
          u7 = 4.9411D4*alphat*et2*wco2
          u8 = 3.9744D4*alphat*et4*wco2
          u9 = 1.0447D5*alphat*et4*et2*wco2
          u13 = 2.8388D3*alphat*et4*wco2
          tpath = tbar(n,kn)
          tlocal = tbar(n,kn)
          tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
          posqt = (pinpl(n,kn)*rsslp+dpfco2*tcrfac)*rsqti
          rbeta7 = d_one/(5.3228D0*posqt)
          rbeta8 = d_one/(10.6576D0*posqt)
          rbeta9 = rbeta7
          rbeta13 = rbeta9
          f2co2 = u7/dsqrt(d_four+u7*(d_one+rbeta7)) + &
                  u8/dsqrt(d_four+u8*(d_one+rbeta8)) + &
                  u9/dsqrt(d_four+u9*(d_one+rbeta9))
          f3co2 = u13/dsqrt(d_four+u13*(d_one+rbeta13))
          tmp1 = dlog(d_one+f1sqwp)
          tmp2 = dlog(d_one+f2co2)
          tmp3 = dlog(d_one+f3co2)
          absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
          abso(6) = trab2*emm(kn)*absbnd
          tco2(n) = d_one/(d_one+d_10*u7/dsqrt(d_four+u7*(d_one+rbeta7)))
!         trab3(n) = 1. - bndfct*absbnd
          absgasnxt(n,k2,kn) = abso(1) + abso(2) + abso(3) + &
                               abso(4) + abso(5) + abso(6) 
        end do
      end do
!
!       Calculate trace gas absorptivity for nearest layer
!
      do kn = 1 , 4
        call trcabn(n1,n2,k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4, &
                    uco211,uco212,uco213,uco221,uco222,uco223,tbar,   &
                    bplnk,winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,ux, &
                    pnew,abstrc,uinpl)
!
!       Total next layer absorptivity:
!
        do n = n1 , n2
          absgasnxt(n,k2,kn) = absgasnxt(n,k2,kn) + abstrc(n)
        end do
      end do
    end do  !  end of nearest layer level loop
!
    call time_end(subroutine_name,indx)
  end subroutine radabs
!
!-----------------------------------------------------------------------
!
! Compute emissivity for H2O, CO2, O3
!
! H2O  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
!
!
! CO2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019
!
! O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R. Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! Computes individual emissivities, accounting for band overlap, and
! sums to obtain the total.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
!------------------------------Arguments--------------------------------
!
!     Output arguments
!
! emplnk  - emissivity Planck factor
! emstrc  - total trace gas emissivity
!
  subroutine radems(n1,n2,pint,emsgastot)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:,:) , intent(in) :: pint
    real(dp) , pointer , dimension(:,:) , intent(out) :: emsgastot
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! iband   - H2o band index
!
!     Local variables for H2O:
!
! h2oems  - H2o emissivity
! tpathe  - Used to compute h2o emissivity
! a       - Eq(2) in table A3a of R&D
! corfac  - Correction factors in table A3b rotation band absorptivity
! dtp     - Path temperature minus 300 K used in
! dtx     - Planck temperature minus 250 K
! dty     - Path temperature minus 250 K
! dtz     - Planck temperature minus 300 K
! emis    - Total emissivity (h2o+co2+o3)
! rsum    - Eq(1) in table A2 of R&D
! term1   - Equation(5) in table A3a of R&D(1986)
! term2   - Delta a(Te) in table A3a of R&D(1986)
! term3   - B(T) function for rotation and vibration-rotation band emissivity
! term4   - Equation(6) in table A3a of R&D(1986)
! term5   - Delta a(Tp) in table A3a of R&D(1986)
! xterm6  - B(T) function for window region
! term7   - Kl_inf(i) in eq(8) of table A3a of R&D
! term8   - Delta kl_inf(i) in eq(8)
! xterm9  - B(T) function for 500-800 cm-1 region
! tr1     - Equation(6) in table A2 for 650-800
! tr2     - Equation(6) in table A2 for 500-650
! tr3     - Equation(4) in table A2 for 650-800
! tr4     - Equation(4),table A2 of R&D for 500-650 
! tr7     - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! tr8     - Equation (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! uc      - Y + 0.002U in eq(8) of table A2 of R&D
! xpnew   - Effective pressure for h2o linewidth
! trline  - Transmission due to H2O lines in window
! k21     - Exponential coefficient used to calc rot band transmissivity
!           in the 650-800 cm-1 region (tr1)
! k22     - Exponential coefficient used to calc rot band transmissivity
!           in the 500-650 cm-1 region (tr2)
! u       - Pressure weighted H2O path length
! uc1     - H2o continuum pathlength 500-800 cm-1
! a11     - A1 in table A3b for rotation band emiss
! a31     - A3 in table A3b for rotation band emiss
! a21     - First part in numerator of A2 table A3b
! a22     - Second part in numertor of A2 table A3b
! a23     - Denominator of A2 table A3b (rot band)
! t1t4    - Eq(3) in table A3a of R&D
! t2t5    - Eq(4) in table A3a of R&D
! fwk     - Equation(33) in R&D far wing correction
! a41     - Numerator in A2 in Vib-rot (table A3b)
! a51     - Denominator in A2 in Vib-rot(table A3b)
! a61     - A3 factor for Vib-rot band in table A3b
! phi     - Eq(11) in table A3a of R&D
! psi     - Eq(12) in table A3a of R&D
! ubar    - H2o scaled path comment eq(10) table A2
! g1      - Part of eq(10) table A2
! pbar    - H2o scaled pres comment eq(10) table A2
! g3      - Part of eq(10) table A2
! g2      - Part of arguement in eq(10) in table A2
! g4      - Arguement in exp() in eq(10) table A2
! cf812   - Eq(11) in table A2 of R&D
! troco2  - H2o overlap factor for co2 absorption
!
!     Local variables for CO2:
!
! co2ems  - Co2 emissivity
! co2plk  - Used to compute co2 emissivity
! xsum    - Used to calculate path temperature
! t1i     - Co2 hot band temperature factor
! sqti    - Sqrt of temperature
! pi      - Pressure used in co2 mean line width
! et      - Co2 hot band factor
! et2     - Co2 hot band factor
! et4     - Co2 hot band factor
! omet    - Co2 stimulated emission term
! ex      - Part of co2 planck function
! f1co2   - Co2 weak band factor
! f2co2   - Co2 weak band factor
! f3co2   - Co2 weak band factor
! t1co2   - Overlap factor weak bands strong band
! sqwp    - Sqrt of co2 pathlength
! f1sqwp  - Main co2 band factor
! oneme   - Co2 stimulated emission term
! alphat  - Part of the co2 stimulated emiss term
! wco2    - Consts used to define co2 pathlength
! posqt   - Effective pressure for co2 line width
! rbeta7  - Inverse of co2 hot band line width par
! rbeta8  - Inverse of co2 hot band line width par
! rbeta9  - Inverse of co2 hot band line width par
! rbeta13 - Inverse of co2 hot band line width par
! tpath   - Path temp used in co2 band model
! tmp1    - Co2 band factor
! tmp2    - Co2 band factor
! tmp3    - Co2 band factor
! tlayr5  - Temperature factor in co2 Planck func
! rsqti   - Reciprocal of sqrt of temperature
! exm1sq  - Part of co2 Planck function
! u7      - Absorber amount for various co2 band systems
! u8      - Absorber amount for various co2 band systems
! u9      - Absorber amount for various co2 band systems
! u13     - Absorber amount for various co2 band systems
!
!     Local variables for O3:
!
! o3ems   - Ozone emissivity
! dbvtt   - Tmp drvtv of planck fctn for tplnke
! te      - Temperature factor
! u1      - Path length factor
! u2      - Path length factor
! phat    - Effecitive path length pressure
! tlocal  - Local planck function temperature
! tcrfac  - Scaled temperature factor
! beta    - Absorption funct factor voigt effect
! realnu  - Absorption function factor
! o3bndi  - Band absorption factor
!
!     Transmission terms for various spectral intervals:
!
! trem4   - H2o   800 - 1000 cm-1
! trem6   - H2o  1000 - 1200 cm-1
! absbnd  - Proportional to co2 band absorptance
! tco2    - co2 overlap factor
! th2o    - h2o overlap factor
! to3     - o3 overlap factor
!
    real(dp) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,    &
               absbnd , alphat , beta , cf812 , et , et2 , et4 , ex , &
               exm1sq , f1co2 , f1sqwp , f2co2 , f3co2 , fwk , g1 ,   &
               g2 , g3 , g4 , o3bndi , omet , oneme , pbar , phat ,   &
               phi , pi , posqt , psi , k21 , k22 , trem4 , trem6 ,   &
               rbeta13 , rbeta7 , rbeta8 , rbeta9 , realnu , rsqti ,  &
               sqti , sqwp , t1co2 , t1i , t1t4 , t2t5 , tpathe ,     &
               tcrfac , te , tlayr5 , tlocal , tmp1 , tmp2 , tmp3 ,   &
               tpath , u1 , u13 , u2 , u7 , u8 , u9 , ubar , wco2 ,   &
               tr1 , tr2 , tr3 , tr4 , tr7 , tr8 , corfac , dbvtt ,   &
               dtp , dtz , xpnew , rsum , u , uc , uc1 , troco2
    real(dp) , dimension(4) :: term1 , term2 , term3 , term4 , term5
    real(dp) , dimension(4) :: emis
    real(dp) :: xterm6 , xterm9
    real(dp) , dimension(2) :: term7 , term8 , trline
    integer :: n , k , kk , iband , l
!
    character (len=64) :: subroutine_name='radems'
    integer :: indx = 0
!
!-----------------------------------------------------------------------
    call time_begin(subroutine_name,indx)
!
!   Planck function for co2
!
    do n = n1 , n2
      ex = dexp(960.0D0/tplnke(n))
      co2plk(n) = 5.0D8/((tplnke(n)**d_four)*(ex-d_one))
      co2t(n,1) = tplnke(n)
      xsum(n) = co2t(n,1)*pint(n,1)
    end do
    kk = 1
    do k = kzp1 , 2 , -1
      kk = kk + 1
      do n = n1 , n2
        xsum(n) = xsum(n) + tlayr(n,kk)*(pint(n,kk)-pint(n,kk-1))
        ex = dexp(960.0D0/tlayr(n,kk))
        tlayr5 = tlayr(n,kk)*tlayr4(n,kk)
        co2eml(n,kk-1) = 1.2D11*ex/(tlayr5*(ex-d_one)**d_two)
        co2t(n,kk) = xsum(n)/pint(n,kk)
      end do
    end do
!
!   bndfct = 2.d0*22.18/(dsqrt(196.d0)*300.)
!
!   Calculate trace gas Planck functions
!
    call trcplk(n1,n2,tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
!   Interface loop
!
    do k = 1 , kzp1
      do n = n1 , n2
!
!       H2O emissivity
!
!       emis(1)     0 -  800 cm-1   rotation band
!       emis(2)  1200 - 2200 cm-1   vibration-rotation band
!       emis(3)   800 - 1200 cm-1   window
!       emis(4)   500 -  800 cm-1   rotation band overlap with co2
!
!       For the p type continuum
!
        uc = s2c(n,k) + 2.0D-3*plh2o(n,k)
        u = plh2o(n,k)
!
!       Apply scaling factor for 500-800 continuum
!
        uc1 = (s2c(n,k)+1.7D-3*plh2o(n,k)) * &
              (d_one+d_two*s2c(n,k))/(d_one+15.0D0*s2c(n,k))
        tpathe = s2t(n,k)/plh2o(n,k)
        dtx(n) = tplnke(n) - 250.0D0
        dty(n) = tpathe - 250.0D0
!
!       emis(1)     0 -  800 cm-1   rotation band
!
        do iband = 1 , 3 , 2
          term1(iband) = coefe(1,iband) + coefe(2,iband)*dtx(n) *   &
                         (d_one+c1(iband)*dtx(n))
          term2(iband) = coefb(1,iband) + coefb(2,iband)*dtx(n) *   &
                         (d_one+c2(iband)*dtx(n)*(d_one+c3(iband)*dtx(n)))
          term3(iband) = coefd(1,iband) + coefd(2,iband)*dtx(n) *   &
                         (d_one+c4(iband)*dtx(n)*(d_one+c5(iband)*dtx(n)))
          term4(iband) = coefa(1,iband) + coefa(2,iband)*dty(n) *   &
                         (d_one+c6(iband)*dty(n))
          term5(iband) = coefc(1,iband) + coefc(2,iband)*dty(n) *   &
                         (d_one+c7(iband)*dty(n))
        end do
        dtp = dty(n) - 50.0D0
        dtz = dtx(n) - 50.0D0
        a11 = 0.37D0 - 3.33D-5*dtz + 3.33D-6*dtz*dtz
        a31 = 1.07D0 - 1.00D-3*dtp + 1.475D-5*dtp*dtp
        a21 = 1.3870D0 + 3.80D-3*dtz - 7.8D-6*dtz*dtz
        a22 = d_one - 1.21D-3*dtp - 5.33D-6*dtp*dtp
        a23 = 0.9D0 + 2.62D0*dsqrt(u)
        corfac = a31*(a11+((a21*a22)/a23))
        t1t4 = term1(1)*term4(1)
        t2t5 = term2(1)*term5(1)
        a = t1t4 + t2t5/(d_one+t2t5*dsqrt(u)*corfac)
        fwk = fwcoef + fwc1/(d_one+fwc2*u)
        rsum = dexp(-a*(dsqrt(u)+fwk*u))
        emis(1) = (d_one-rsum)*term3(1)
!       trem1  = rsum
!
!       emis(2)  1200 - 2200 cm-1   vibration-rotation band
!
        a41 = 1.75D0 - 3.96D-3*dtz
        a51 = 1.00D0 + 1.3D0*dsqrt(u)
        a61 = 1.00D0 + 1.25D-3*dtp + 6.25D-5*dtp*dtp
        corfac = 0.3D0*(d_one+(a41)/(a51))*a61
        t1t4 = term1(3)*term4(3)
        t2t5 = term2(3)*term5(3)
        a = t1t4 + t2t5/(d_one+t2t5*dsqrt(u)*corfac)
        fwk = fwcoef + fwc1/(d_one+fwc2*u)
        rsum = dexp(-a*(dsqrt(u)+fwk*u))
        emis(2) = (d_one-rsum)*term3(3)
!       trem7 = rsum
!
!       Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
!       emis(3)   800 - 1200 cm-1   window
!
        do l = 1 , 2
          phi = a1(l)*(dty(n)+15.0D0)+a2(l)*(dty(n)+15.0D0)**d_two
          psi = b1(l)*(dty(n)+15.0D0)+b2(l)*(dty(n)+15.0D0)**d_two
          phi = dexp(phi)
          psi = dexp(psi)
          ubar = w(n,k)*phi
          ubar = (ubar*1.66D0)*r80257
          xpnew = u/w(n,k)
          pbar = xpnew*(psi/phi)
          cf812 = cfa1 + ((d_one-cfa1)/(d_one+ubar*pbar*d_10))
          g1 = (realk(l)*pbar)/(d_two*st(l))
          g2 = d_one + (ubar*d_four*st(l)*cf812)/pbar
          g3 = dsqrt(g2) - d_one
          g4 = g1*g3
          trline(l) = dexp(-g4)
        end do
        xterm6 = coeff(1,1) + coeff(2,1)*dtx(n) *        &
                 (d_one+c8*dtx(n)*(d_one+c10*dtx(n) *    &
                 (d_one+c12*dtx(n)*(d_one+c14*dtx(n)))))
        term7(1) = coefj(1,1)+coefj(2,1)*dty(n)*(d_one+c16*dty(n))
        term8(1) = coefk(1,1)+coefk(2,1)*dty(n)*(d_one+c17*dty(n))
        term7(2) = coefj(1,2)+coefj(2,2)*dty(n)*(d_one+c26*dty(n))
        term8(2) = coefk(1,2)+coefk(2,2)*dty(n)*(d_one+c27*dty(n))
        trem4 = dexp(-(coefg(1,1)+coefg(2,1)*dtx(n))*uc)*trline(2)
        trem6 = dexp(-(coefg(1,2)+coefg(2,2)*dtx(n))*uc)*trline(1)
        emis(3) = xterm6*(d_one-trem4*d_half-trem6*d_half)
!
!       emis(4)   500 -  800 cm-1   rotation band overlap with co2
!
        k21 = term7(1) + term8(1)/(d_one+(c30+c31*(dty(n)-d_10) * &
                 (dty(n)-d_10))*dsqrt(u))
        k22 = term7(2) + term8(2)/(d_one+(c28+c29*(dty(n)-d_10))*dsqrt(u))
        xterm9 = coefi(1,1) + coefi(2,1)*dtx(n) *        &
                (d_one+c18*dtx(n)*(d_one+c20*dtx(n) *   &
                (d_one+c22*dtx(n)*(d_one+c24*dtx(n)))))
        fwk = fwcoef + fwc1/(d_one+fwc2*u)
        tr1 = dexp(-(k21*(dsqrt(u)+fc1*fwk*u)))
        tr2 = dexp(-(k22*(dsqrt(u)+fc1*fwk*u)))
        tr3 = dexp(-((coefh(1,1)+coefh(2,1)*dtx(n))*uc1))
        tr4 = dexp(-((coefh(1,2)+coefh(2,2)*dtx(n))*uc1))
        tr7 = tr1*tr3
        tr8 = tr2*tr4
        emis(4) = xterm9*d_half*(tr1-tr7+tr2-tr8)
        h2oems(n,k) = emis(1) + emis(2) + emis(3) + emis(4)
        troco2 = 0.65D0*tr7 + 0.35D0*tr8
        th2o(n) = tr8
!       trem2(n) = troco2
!
!       CO2 emissivity for 15 micron band system
!
        t1i = dexp(-480.0D0/co2t(n,k))
        sqti = dsqrt(co2t(n,k))
        rsqti = d_one/sqti
        et = t1i
        et2 = et*et
        et4 = et2*et2
        omet = d_one - 1.5D0*et2
        f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
        sqwp = dsqrt(plco2(n,k))
        f1sqwp = f1co2*sqwp
        t1co2 = d_one/(d_one+245.18D0*omet*sqwp*rsqti)
        oneme = d_one - et2
        alphat = oneme**3.0D0*rsqti
        wco2 = 2.5221D0*co2vmr*pint(n,k)*regravgts
        u7 = 4.9411D4*alphat*et2*wco2
        u8 = 3.9744D4*alphat*et4*wco2
        u9 = 1.0447D5*alphat*et4*et2*wco2
        u13 = 2.8388D3*alphat*et4*wco2
!
        tpath = co2t(n,k)
        tlocal = tplnke(n)
        tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
        pi = pint(n,k)*rsslp + d_two*dpfco2*tcrfac
        posqt = pi/(d_two*sqti)
        rbeta7 = d_one/(5.3288D0*posqt)
        rbeta8 = d_one/(10.6576D0*posqt)
        rbeta9 = rbeta7
        rbeta13 = rbeta9
        f2co2 = (u7/dsqrt(d_four+u7*(d_one+rbeta7))) + &
                (u8/dsqrt(d_four+u8*(d_one+rbeta8))) + &
                (u9/dsqrt(d_four+u9*(d_one+rbeta9)))
        f3co2 = u13/dsqrt(d_four+u13*(d_one+rbeta13))
        tmp1 = dlog(d_one+f1sqwp)
        tmp2 = dlog(d_one+f2co2)
        tmp3 = dlog(d_one+f3co2)
        absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
        tco2(n) = d_one/(d_one+d_10*(u7/dsqrt(d_four+u7*(d_one+rbeta7))))
        co2ems(n,k) = troco2*absbnd*co2plk(n)
        ex = dexp(960.0D0/tint(n,k))
        exm1sq = (ex-d_one)**d_two
        co2em(n,k) = 1.2D11*ex/(tint(n,k)*tint4(n,k)*exm1sq)
!       trem3(n) = 1. - bndfct*absbnd
!
!       O3 emissivity
!
        h2otr(n,k) = dexp(-12.0D0*s2c(n,k))
        te = (co2t(n,k)/293.0D0)**0.7D0
        u1 = 18.29D0*plos(n,k)/te
        u2 = 0.5649D0*plos(n,k)/te
        phat = plos(n,k)/plol(n,k)
        tlocal = tplnke(n)
        tcrfac = dsqrt(tlocal*r250)*te
        beta = (d_one/0.3205D0)*((d_one/phat)+(dpfo3*tcrfac))
        realnu = (d_one/beta)*te
        o3bndi = 74.0D0*te*(tplnke(n)/375.0D0)*dlog(d_one+fo3(u1,realnu) + &
                 fo3(u2,realnu))
        dbvtt = dbvt(tplnke(n))
        o3ems(n,k) = dbvtt*h2otr(n,k)*o3bndi
        to3(n) = d_one/(d_one+0.1D0*fo3(u1,realnu)+0.1D0*fo3(u2,realnu))
!       trem5(n)    = d_one-(o3bndi/(1060-980.))
      end do
!
!     Calculate trace gas emissivities
!
      call trcems(n1,n2,k,co2t,pint,ucfc11,ucfc12,un2o0,un2o1,  &
                  bn2o0,bn2o1,uch4,bch4,uco211,uco212,uco213,uco221,  &
                  uco222,uco223,uptype,w,s2c,ux,emplnk,th2o,tco2,to3, &
                  emstrc)
!
!     Total emissivity:
!
      do n = n1 , n2
        emsgastot(n,k) = h2oems(n,k)+co2ems(n,k)+o3ems(n,k)+emstrc(n,k)
      end do
    end do  ! End of interface loop
!
    call time_end(subroutine_name,indx)
  end subroutine radems
!
!-----------------------------------------------------------------------
!
! Computes the path length integrals to the model interfaces given the
! ozone volume mixing ratio
!
!---------------------------Code history--------------------------------
!
! Original version:     CCM1
! Standardized:         J. Rosinski, June 1992
! Reviewed:             J. Kiehl, B. Briegleb, August 1992
! Mixing ratio version: Bruce Biegleb, September 1992
!
!-----------------------------------------------------------------------
!
!------------------------------Input arguments--------------------------
!
! o3vmr   - ozone volume mixing ratio
! pint    - Model interface pressures
!
!----------------------------Output arguments---------------------------
!
!
  subroutine radoz2(n1,n2,o3vmr,pint)
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:,:) :: o3vmr
    real(dp) , pointer , dimension(:,:) :: pint
    intent (in) o3vmr , pint
!
!-----------------------------------------------------------------------
!
    integer :: n , k
    character (len=64) :: subroutine_name='radoz2'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Evaluate the ozone path length integrals to interfaces;
!   factors of 0.1 and 0.01 to convert pressures from cgs to mks:
!
!   Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
!
    do n = n1 , n2
      plos(n,1) = 0.1D0*cplos*o3vmr(n,1)*pint(n,1)
      plol(n,1) = 0.01D0*cplol*o3vmr(n,1)*pint(n,1)*pint(n,1)
    end do
    do k = 2 , kzp1
      do n = n1 , n2
        plos(n,k) = plos(n,k-1) + &
             0.1D0*cplos*o3vmr(n,k-1)*(pint(n,k)-pint(n,k-1))
        plol(n,k) = plol(n,k-1) + 0.01D0*cplol*o3vmr(n,k-1) * &
                    (pint(n,k)*pint(n,k)-pint(n,k-1)*pint(n,k-1))
      end do
    end do
    call time_end(subroutine_name,indx)
  end subroutine radoz2
!
!-----------------------------------------------------------------------
!
! Compute temperatures and path lengths for longwave radiation
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
  subroutine radtpl(n1,n2,ts,tnm,pmln,qnm,piln,pint)
!
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:) , intent(in) :: ts
    real(dp) , pointer , dimension(:,:) , intent(in) :: tnm , pmln , &
      qnm , piln , pint
!
!---------------------------Local variables-----------------------------
!
! i      - Longitude index
! k      - Level index
! dy     - Thickness of layer for tmp interp
! dpnm   - Pressure thickness of layer
! dpnmsq - Prs squared difference across layer
! rtnm   - Inverse level temperature
!
!-----------------------------------------------------------------------
!
    real(dp) :: dpnm , dpnmsq , dy , rtnm
    integer :: n , k
    character (len=64) :: subroutine_name='radtpl'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Set the top and bottom intermediate level temperatures,
!   top level planck temperature and top layer temp**4.
!
!   Tint is lower interface temperature
!   (not available for bottom layer, so use ground temperature)
!
    do n = n1 , n2
      tint(n,kzp1) = ts(n)
      tint4(n,kzp1) = tint(n,kz + 1)**d_four
      tplnka(n,1) = tnm(n,1)
      tint(n,1) = tplnka(n,1)
      tlayr4(n,1) = tplnka(n,1)**d_four
      tint4(n,1) = tlayr4(n,1)
    end do
!
!   Intermediate level temperatures are computed using temperature
!   at the full level below less dy*delta t,between the full level
!
    do k = 2 , kz
      do n = n1 , n2
        dy = (piln(n,k)-pmln(n,k))/(pmln(n,k-1)-pmln(n,k))
        tint(n,k) = tnm(n,k) - dy*(tnm(n,k)-tnm(n,k-1))
        tint4(n,k) = tint(n,k)**d_four
      end do
    end do
!
!   Now set the layer temp=full level temperatures and establish a
!   planck temperature for absorption (tplnka) which is the average
!   the intermediate level temperatures.  Note that tplnka is not
!   equal to the full level temperatures.
!
    do k = 2 , kzp1
      do n = n1 , n2
        tlayr(n,k) = tnm(n,k-1)
        tlayr4(n,k) = tlayr(n,k)**d_four
        tplnka(n,k) = (tint(n,k)+tint(n,k-1))*d_half
      end do
    end do
!
!   Calculate tplank for emissivity calculation.
!   Assume isothermal tplnke i.e. all levels=ttop.
!
    do n = n1 , n2
      tplnke(n) = tplnka(n,1)
      tlayr(n,1) = tint(n,1)
    end do
!
!   Now compute h2o path fields:
!
    do n = n1 , n2
      s2t(n,1) = plh2o(n,1)*tnm(n,1)
!     ccm3.2
!     w(n,1)   = (plh2o(n,1)*2.) / pint(n,1)
!     s2c(n,1) = plh2o(n,1) * qnm(n,1) * repsil
   
!     ccm3.6.6
      w(n,1) = sslp*(plh2o(n,1)*d_two)/pint(n,1)
      rtnm = d_one/tnm(n,1)
      s2c(n,1) = plh2o(n,1)*dexp(1800.0D0*(rtnm-r296))*qnm(n,1)*repsil
    end do
    do k = 1 , kz
      do n = n1 , n2
        dpnm = pint(n,k+1) - pint(n,k)
        dpnmsq = pint(n,k+1)**d_two - pint(n,k)**d_two
        rtnm = d_one/tnm(n,k)
        s2t(n,k+1) = s2t(n,k) + rgsslp*dpnmsq*qnm(n,k)*tnm(n,k)
        w(n,k+1) = w(n,k) + regravgts*qnm(n,k)*dpnm
        s2c(n,k+1) = s2c(n,k) + rgsslp*dpnmsq*qnm(n,k) * &
                     dexp(1800.0D0*(rtnm-r296))*qnm(n,k)*repsil
      end do
    end do
!
    call time_end(subroutine_name,indx)
  end subroutine radtpl
!
!-----------------------------------------------------------------------
!
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
!
! Convert model pressures to cgs, compute path length arrays needed for the
! longwave radiation, and compute ozone mixing ratio, needed for the solar
! radiation.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
  subroutine radinp(n1,n2,pmid,pint,h2ommr,cld,o3vmr,pmidrd, &
                    pintrd,plco2,plh2o,tclrsf,o3mmr)
    implicit none
!
    integer , intent(in) :: n1 , n2
    real(dp) , pointer , dimension(:,:) :: cld , pint , pintrd , plco2 ,   &
           plh2o , tclrsf
    real(dp) , pointer , dimension(:,:) :: h2ommr , o3mmr , o3vmr , pmid , &
           pmidrd
    intent (in) cld , h2ommr , o3vmr , pint , pmid
    intent (out) o3mmr , plco2 , pmidrd
    intent (inout) pintrd , plh2o , tclrsf
!
    real(dp) :: cpwpl , vmmr
    integer :: n , k
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pmid    - Pressure at model mid-levels (pascals)
! pint    - Pressure at model interfaces (pascals)
! h2ommr  - H2o mass mixing ratio
! cld     - Fractional cloud cover
! o3vmr   - ozone volume mixing ratio
!
!     Output arguments
!
! pmidrd  - Pressure at mid-levels (dynes/cm*2)
! pintrd  - Pressure at interfaces (dynes/cm*2)
! plco2   - Vert. pth lngth of co2 (prs-weighted)
! plh2o   - Vert. pth lngth h2o vap.(prs-weighted)
! tclrsf  - Product of clr-sky fractions from top of atmosphere to level.
! o3mmr   - Ozone mass mixing ratio
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude loop index
! k       - Vertical loop index
! cpwpl   - Const in co2 mixing ratio to path length conversn
! vmmr    - Ozone volume mixing ratio
!
!-----------------------------------------------------------------------
!
!   Compute solar distance factor and cosine solar zenith angle usi
!   day value where a round day (such as 213.0) refers to 0z at
!   Greenwich longitude.
!
!   Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
!   Processes in Meterology and Climatology, Elsevier Scientific
!   Publishing Company, New York  p. 57, p. 62,63.
!
    character (len=64) :: subroutine_name='radinp'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
    !
    ! Convert pressure from pascals to dynes/cm2
    !
    do k = 1 , kz
      do n = n1 , n2
        pmidrd(n,k) = pmid(n,k)*d_10
        pintrd(n,k) = pint(n,k)*d_10
      end do
    end do
    do n = n1 , n2
      pintrd(n,kzp1) = pint(n,kzp1)*d_10
    end do
    !
    ! Compute path quantities used in the longwave radiation:
    !
    vmmr = amco2/amd
    cpwpl = vmmr*d_half/(egravgts*sslp)
    do n = n1 , n2
      plh2o(n,1) = rgsslp*h2ommr(n,1)*pintrd(n,1)*pintrd(n,1)
      plco2(n,1) = co2vmr*cpwpl*pintrd(n,1)*pintrd(n,1)
      tclrsf(n,1) = d_one
    end do
    do k = 1 , kz
      do n = n1 , n2
        plh2o(n,k+1) = plh2o(n,k) + rgsslp*(pintrd(n,k+1)**d_two - &
                       pintrd(n,k)**d_two) * h2ommr(n,k)
        plco2(n,k+1) = co2vmr*cpwpl*pintrd(n,k+1)**d_two
        tclrsf(n,k+1) = tclrsf(n,k)*(d_one-cld(n,k+1))
      end do
    end do
    !
    ! Convert ozone volume mixing ratio to mass mixing ratio:
    !
    vmmr = amo/amd
    do k = 1 , kz
      do n = n1 , n2
        o3mmr(n,k) = vmmr*o3vmr(n,k)
      end do
    end do
    call time_end(subroutine_name,indx)
  end subroutine radinp
!
  real(dp) function xalpha(w,uu,g,e)
    implicit none
    real(dp) , intent(in) :: w , uu , g , e
    xalpha = 0.75D0*w*uu*((d_one+g*(d_one-w))/(d_one-e*e*uu*uu))
  end function xalpha
  real(dp) function xgamma(w,uu,g,e)
    implicit none
    real(dp) , intent(in) :: w , uu , g , e
    xgamma = (w*d_half)*((3.0D0*g*(d_one-w)*uu*uu+d_one)/(d_one-e*e*uu*uu))
  end function xgamma
  real(dp) function el(w,g)
    implicit none
    real(dp) , intent(in) :: w , g
    el = dsqrt(3.0D0*(d_one-w)*(d_one-w*g))
  end function el
  real(dp) function taus(w,f,t)
    implicit none
    real(dp) , intent(in) :: w , f , t
    taus = (d_one-w*f)*t
  end function taus
  real(dp) function omgs(w,f)
    implicit none
    real(dp) , intent(in) :: w , f
    omgs = (d_one-f)*w/(d_one-w*f)
  end function omgs
  real(dp) function asys(g,f)
    implicit none
    real(dp) , intent(in) :: g , f
    asys = (g-f)/(d_one-f)
  end function asys
  real(dp) function f_u(w,g,e)
    implicit none
    real(dp) , intent(in) :: w , g , e
    f_u = 1.50D0*(d_one-w*g)/e
  end function f_u
  real(dp) function f_n(uu,et)
    implicit none
    real(dp) , intent(in) :: uu , et
    f_n = ((uu+d_one)*(uu+d_one)/et)-((uu-d_one)*(uu-d_one)*et)
  end function f_n
  real(dp) function dbvt(t)
!   Derivative of planck function at 9.6 micro-meter wavelength
    implicit none
    real(dp) , intent(in) :: t
    dbvt = (-2.8911366682D-4 + (2.3771251896D-6+1.1305188929D-10*t)*t) /  &
            (d_one+(-6.1364820707D-3+1.5550319767D-5*t)*t)
  end function dbvt
  real(dp) function fo3(ux,vx)
!   an absorption function factor
    implicit none
    real(dp) , intent(in) :: ux , vx
    fo3 = ux/dsqrt(d_four+ux*(d_one+vx))
  end function fo3
!
end module mod_rad_radiation
