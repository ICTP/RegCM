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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_mpmessage
  use mod_service
  use mod_runparams , only : idirect , ichem , iclimaaer , rcmtimer
  use mod_runparams , only : scon , cftotmax , lsrfhack , scenario , mincld
  use mod_mppparam , only : italk
  use mod_memutil
  use mod_ipcc_scenario

  ! Used by this module only

  use mod_rad_common
  use mod_rad_tracer
  use mod_rad_aerosol

  implicit none

  private

  public :: allocate_mod_rad_radiation , radini , radctl , radtype

  type radtype
    integer(ik4) :: n1 , n2
    real(rkx) :: eccf
    logical :: labsem
    integer(ik4) , dimension(:) , pointer :: ioro
    real(rkx) , dimension(:) , pointer :: dlat
    real(rkx) , dimension(:) , pointer :: xptrop
    real(rkx) , dimension(:) , pointer :: ts
    real(rkx) , dimension(:) , pointer :: ps
    real(rkx) , dimension(:) , pointer :: totcl
    real(rkx) , dimension(:) , pointer :: totci
    real(rkx) , dimension(:) , pointer :: totwv
    real(rkx) , dimension(:) , pointer :: fsns
    real(rkx) , dimension(:) , pointer :: flwds
    real(rkx) , dimension(:) , pointer :: sols
    real(rkx) , dimension(:) , pointer :: soll
    real(rkx) , dimension(:) , pointer :: solsd
    real(rkx) , dimension(:) , pointer :: solld
    real(rkx) , dimension(:) , pointer :: emiss
    real(rkx) , dimension(:) , pointer :: fsnt
    real(rkx) , dimension(:) , pointer :: fsntc
    real(rkx) , dimension(:) , pointer :: fsnsc
    real(rkx) , dimension(:) , pointer :: flnt
    real(rkx) , dimension(:) , pointer :: lwout
    real(rkx) , dimension(:) , pointer :: lwin
    real(rkx) , dimension(:) , pointer :: flns
    real(rkx) , dimension(:) , pointer :: flntc
    real(rkx) , dimension(:) , pointer :: flnsc
    real(rkx) , dimension(:) , pointer :: solin
    real(rkx) , dimension(:) , pointer :: solout
    real(rkx) , dimension(:) , pointer :: alb
    real(rkx) , dimension(:) , pointer :: albc
    real(rkx) , dimension(:) , pointer :: fsds
    real(rkx) , dimension(:) , pointer :: fsnirt
    real(rkx) , dimension(:) , pointer :: fsnrtc
    real(rkx) , dimension(:) , pointer :: fsnirtsq
    real(rkx) , dimension(:) , pointer :: totcf
    real(rkx) , dimension(:) , pointer :: czen
    real(rkx) , dimension(:) , pointer :: adirsw
    real(rkx) , dimension(:) , pointer :: adifsw
    real(rkx) , dimension(:) , pointer :: adirlw
    real(rkx) , dimension(:) , pointer :: adiflw
    real(rkx) , dimension(:) , pointer :: asw
    real(rkx) , dimension(:) , pointer :: alw
    real(rkx) , dimension(:) , pointer :: abv
    real(rkx) , dimension(:) , pointer :: sol
    real(rkx) , dimension(:) , pointer :: aeradfo
    real(rkx) , dimension(:) , pointer :: aeradfos
    real(rkx) , dimension(:) , pointer :: aerlwfo
    real(rkx) , dimension(:) , pointer :: aerlwfos
    real(rkx) , dimension(:) , pointer :: asaeradfo
    real(rkx) , dimension(:) , pointer :: asaeradfos
    real(rkx) , dimension(:) , pointer :: asaerlwfo
    real(rkx) , dimension(:) , pointer :: asaerlwfos
    logical , dimension(:) , pointer :: czengt0
    real(rkx) , dimension(:,:) , pointer :: pmid
    real(rkx) , dimension(:,:) , pointer :: pint
    real(rkx) , dimension(:,:) , pointer :: pmln
    real(rkx) , dimension(:,:) , pointer :: piln
    real(rkx) , dimension(:,:) , pointer :: t
    real(rkx) , dimension(:,:) , pointer :: q
    real(rkx) , dimension(:,:) , pointer :: ql
    real(rkx) , dimension(:,:) , pointer :: qi
    real(rkx) , dimension(:,:) , pointer :: dz
    real(rkx) , dimension(:,:) , pointer :: rh
    real(rkx) , dimension(:,:) , pointer :: rho
    real(rkx) , dimension(:,:) , pointer :: cld
    real(rkx) , dimension(:,:) , pointer :: effcld
    real(rkx) , dimension(:,:) , pointer :: clwp
    real(rkx) , dimension(:,:) , pointer :: qrs
    real(rkx) , dimension(:,:) , pointer :: qrl
    real(rkx) , dimension(:,:) , pointer :: rel
    real(rkx) , dimension(:,:) , pointer :: rei
    real(rkx) , dimension(:,:) , pointer :: fice
    real(rkx) , dimension(:,:) , pointer :: o3vmr
    real(rkx) , dimension(:,:) , pointer :: emsgastot
    real(rkx) , dimension(:,:,:) , pointer :: absgasnxt
    real(rkx) , dimension(:,:,:) , pointer :: absgastot
    real(rkx) , dimension(:,:,:) , pointer :: tauxcl
    real(rkx) , dimension(:,:,:) , pointer :: tauxci
    real(rkx) , dimension(:,:,:) , pointer :: outtaucl
    real(rkx) , dimension(:,:,:) , pointer :: outtauci
  end type radtype

  integer(ik4) , parameter :: nlwspi = 14

  logical , save :: linteract = .false.
  logical , save :: lzero = .false.

  integer(ik4) :: npoints

  real(rkx) , pointer , dimension(:) :: diralb , difalb
  real(rkx) , pointer , dimension(:) :: cfc11mmr , cfc12mmr
  real(rkx) , pointer , dimension(:) :: ch4mmr , n2ommr
  real(rkx) , pointer , dimension(:) :: co2mmr , co2vmr
  real(rkx) , pointer , dimension(:) :: tauaer , tauasc , gtota , ftota
  real(rkx) , pointer , dimension(:,:) :: emplnk , emstot
  real(rkx) , pointer , dimension(:,:,:) :: absnxt , abstot , xuinpl
  !
  ! Background aerosol mass mixing ratio
  !
  real(rkx) , pointer , dimension(:,:) :: aermmb
  !
  ! Trace gas variables
  !
  ! bch4     - pressure factor for ch4
  real(rkx) , pointer , dimension(:,:) :: bch4
  ! bn2o0   - pressure factor for n2o
  real(rkx) , pointer , dimension(:,:) :: bn2o0
  ! bn2o1   - pressure factor for n2o
  real(rkx) , pointer , dimension(:,:) :: bn2o1
  ! co2em   - Layer co2 normalized planck funct. derivative
  real(rkx) , pointer , dimension(:,:) :: co2em
  ! co2eml  - Interface co2 normalized planck funct. deriv.
  real(rkx) , pointer , dimension(:,:) :: co2eml
  ! co2t    - Prs wghted temperature path
  real(rkx) , pointer , dimension(:,:) :: co2t
  ! h2otr   - H2o trnmsn for o3 overlap
  real(rkx) , pointer , dimension(:,:) :: h2otr
  ! ucfc11  - CFC11 path length
  real(rkx) , pointer , dimension(:,:) :: ucfc11
  ! ucfc12  - CFC12 path length
  real(rkx) , pointer , dimension(:,:) :: ucfc12
  ! un2o0   - N2O path length
  real(rkx) , pointer , dimension(:,:) :: un2o0
  ! un2o1   - N2O path length (hot band)
  real(rkx) , pointer , dimension(:,:) :: un2o1
  ! uch4    - CH4 path length
  real(rkx) , pointer , dimension(:,:) :: uch4
  ! uco211  - CO2 9.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco211
  ! uco212  - CO2 9.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco212
  ! uco213  - CO2 9.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco213
  ! uco221  - CO2 10.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco221
  ! uco222  - CO2 10.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco222
  ! uco223  - CO2 10.4 micron band path length
  real(rkx) , pointer , dimension(:,:) :: uco223
  ! uptype   - continuum path length
  real(rkx) , pointer , dimension(:,:) :: uptype
  ! plol     - Ozone prs wghted path length
  real(rkx) , pointer , dimension(:,:) :: plol
  ! plos     - Ozone path length
  real(rkx) , pointer , dimension(:,:) :: plos
  ! tplnka   - Planck fnctn level temperature
  real(rkx) , pointer , dimension(:,:) :: tplnka
  ! tint    - Interface temperature
  real(rkx) , pointer , dimension(:,:) :: tint
  ! tint4   - Interface temperature**4
  real(rkx) , pointer , dimension(:,:) :: tint4
  ! tlayr   - Level temperature
  real(rkx) , pointer , dimension(:,:) :: tlayr
  ! tlayr4  - Level temperature**4
  real(rkx) , pointer , dimension(:,:) :: tlayr4
  ! w       - H2o path
  real(rkx) , pointer , dimension(:,:) :: wh2op
  ! s2c     - H2o cont amount
  real(rkx) , pointer , dimension(:,:) :: s2c
  ! s2t     - H2o cont temperature
  real(rkx) , pointer , dimension(:,:) :: s2t
  ! ful     - Total upwards longwave flux
  real(rkx) , pointer , dimension(:,:) :: ful , ful0
  ! fsul    - Clear sky upwards longwave flux
  real(rkx) , pointer , dimension(:,:) :: fsul , fsul0
  ! fdl     - Total downwards longwave flux
  real(rkx) , pointer , dimension(:,:) :: fdl , fdl0
  ! fsdl    - Clear sky downwards longwv flux
  real(rkx) , pointer , dimension(:,:) :: fsdl , fsdl0
  ! fis     - Flx integral sum
  real(rkx) , pointer , dimension(:,:,:) :: fis , fis0
  ! rtclrsf - d_one/tclrsf(k,n)
  real(rkx) , pointer , dimension(:,:) :: rtclrsf
  ! tplnke  - Planck fnctn temperature
  real(rkx) , pointer , dimension(:) :: tplnke
  ! fclb4   - Sig t**4 for cld bottom interfc
  real(rkx) , pointer , dimension(:,:) :: fclb4
  ! fclt4   - Sig t**4 for cloud top interfc
  real(rkx) , pointer , dimension(:,:) :: fclt4
  ! klov    - Cloud lowest level index
  integer(ik4) , pointer , dimension(:) :: klov
  ! khiv    - Cloud highest level index
  integer(ik4) , pointer , dimension(:) :: khiv
  ! khivm   - khiv(n) - 1
  integer(ik4) , pointer , dimension(:) :: khivm
  ! Control logicals
  logical , pointer , dimension(:) :: skip , done
  !
  ! These arrays are defined for kz model layers; 0 refers to the
  ! extra layer on top:
  !
  ! rdir     - Layer reflectivity to direct rad
  ! rdif     - Layer reflectivity to diffuse rad
  ! tdir     - Layer transmission to direct rad
  ! tdif     - Layer transmission to diffuse rad
  ! explay   - Solar beam exp transmission for layer
  ! flxdiv   - Flux divergence for layer
  ! totfld   - Spectrally summed flux divergence
  !
  real(rkx) , pointer , dimension(:,:) :: rdir , rdif , tdir , tdif
  real(rkx) , pointer , dimension(:,:) :: explay , flxdiv , totfld
  !
  ! Cloud radiative property arrays
  !
  ! tauxcl   - water cloud extinction optical depth   ( from interface )
  ! tauxci   - ice cloud extinction optical depth     ( from interface )
  ! wcl      - liquid cloud single scattering albedo
  ! gcl      - liquid cloud asymmetry parameter
  ! fcl      - liquid cloud forward scattered fraction
  ! wci      - ice cloud single scattering albedo
  ! gci      - ice cloud asymmetry parameter
  ! fci      - ice cloud forward scattered fraction
  !
  real(rkx) , pointer , dimension(:,:) :: wcl , gcl , fcl , wci , gci , fci
  !
  ! Layer absorber amounts
  !
  ! uh2o     - Layer absorber amount of h2o
  ! uo3      - Layer absorber amount of  o3
  ! uco2     - Layer absorber amount of co2
  ! uo2      - Layer absorber amount of  o2
  !
  real(rkx) , pointer , dimension(:,:) :: uh2o , uo3 , uco2 , uo2
  !
  ! These arrays are defined at model interfaces; 0 is the top of the
  ! extra layer above the model top; kzp1 is the earth surface:
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
  real(rkx) , pointer , dimension(:,:) :: rupdir , rupdif , rdndif , &
        exptdn , tottrn , fluxup , fluxdn , pflx , fswup , fswdn
  !
  ! abplnk1 - non-nearest layer Plack factor
  ! abplnk2 - nearest layer factor
  !
  real(rkx) , pointer , dimension(:,:,:) :: abplnk1 , abplnk2
  !
  ! solflx   - Solar flux in current interval
  ! uth2o    - Total column  absorber amount of h2o
  ! uto3     - Total column  absorber amount of  o3
  ! utco2    - Total column  absorber amount of co2
  ! uto2     - Total column  absorber amount of  o2
  !
  real(rkx) , pointer , dimension(:) :: solflx
  real(rkx) , pointer , dimension(:) :: utco2 , uth2o , uto2 , uto3
  ! ref net TOA flux
  real(rkx) , pointer , dimension(:) :: toafsnsc ,  toafsntc
  !
  ! o3mmr    - Ozone mass mixing ratio
  ! cfc11    - cfc11 mass mixing ratio
  ! cfc12    - cfc12 mass mixing ratio
  ! ch4      - methane mass mixing ratio
  ! n2o      - nitrous oxide mass mixing ratio
  !
  real(rkx) , pointer , dimension(:,:) :: o3mmr , cfc11 , cfc12 , ch4 , n2o
  !
  ! pbr      - Model mid-level pressures (dynes/cm2)
  ! pnm      - Model interface pressures (dynes/cm2)
  !
  real(rkx) , pointer , dimension(:,:) :: pbr , pnm
  !
  ! plco2    - Prs weighted CO2 path
  ! plh2o    - Prs weighted H2O path
  ! tclrsf   - Total clear sky fraction, level to space
  !
  real(rkx) , pointer , dimension(:,:) :: plco2 , plh2o , tclrsf
  !
  ! fslwdcs  - Downward clear sky long wave flux at surface
  !
  real(rkx) , pointer , dimension(:) :: fslwdcs

  real(rkx) , parameter :: verynearone = 0.999999_rkx
  ! r80257   - Conversion factor for h2o pathlength
  real(rkx) , parameter :: r80257 = d_one/8.0257e-4_rkx
  real(rkx) , parameter :: r293 = d_one/293.0_rkx
  real(rkx) , parameter :: r250 = d_one/250.0_rkx
  ! r3205    - Line width factor for o3 (see R&Di)
  real(rkx) , parameter :: r3205 = d_one/0.3205_rkx
  real(rkx) , parameter :: r300 = d_one/300.0_rkx
  ! r2sslp   - 1/2 of rsslp
  real(rkx) , parameter :: r2sslp = d_one/(d_two*sslp)
  ! r296   - Inverse stand temp for h2o continuum
  real(rkx) , parameter :: r296 = d_one/296.0_rkx
  ! repsil - Inver ratio mol weight h2o to dry air
  real(rkx) , parameter :: repsil = d_one/ep2
  !
  ! Initialize further longwave constants referring to far wing
  ! correction; R&D refers to:
  !
  ! Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
  ! Emissivity and Absorptivity Formulation for Water Vapor
  ! Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
  !
  real(rkx) , parameter :: fwcoef = 0.1_rkx      ! See eq(33) R&D
  real(rkx) , parameter :: fwc1 = 0.30_rkx       ! See eq(33) R&D
  real(rkx) , parameter :: fwc2 = 4.5_rkx        ! See eq(33) and eq(34) in R&D
  real(rkx) , parameter :: fc1 = 2.6_rkx         ! See eq(34) R&D
  !
  ! Initialize ozone data.
  !
  real(rkx) , parameter :: v0 = 22.4136_rkx ! Volume of a gas at stp (m**3/kmol)
  real(rkx) , parameter :: p0 = 0.1_rkx*sslp ! Standard pressure (pascals)
  !
  ! Constants for ozone path integrals (multiplication by 100 for unit
  ! conversion to cgs from mks):
  !
  real(rkx) , parameter :: cplos = v0/(amd*egrav)*d_100
  real(rkx) , parameter :: cplol = v0/(amd*egrav*p0)*d_half*d_100
  !
  ! delta    - Pressure (atmospheres) for stratos. h2o limit
  ! o2mmr    - O2 mass mixing ratio
  !
  real(rkx) , parameter :: delta = 1.70e-3_rkx
  real(rkx) , parameter :: o2mmr = 0.23143_rkx
  !
  ! Minimum total transmission below which no layer computation are done:
  !
  ! trmin   - Minimum total transmission allowed
  ! wray    - Rayleigh single scatter albedo
  ! gray    - Rayleigh asymetry parameter
  ! fray    - Rayleigh forward scattered fraction
  !
  real(rkx) , parameter :: trmin = 1.0e-3_rkx
  real(rkx) , parameter :: wray = 0.999999_rkx
  real(rkx) , parameter :: gray = 0.0_rkx
  real(rkx) , parameter :: fray = 0.1_rkx
  !
  ! H2O DMISSIVITY AND ABSORTIVITY CODFFICIDNTS
  !
  real(rkx) , dimension(3,4) , parameter :: coefa = reshape([ &
    1.01400e+0_rkx , 6.41695e-3_rkx , 2.85787e-5_rkx , &
    1.01320e+0_rkx , 6.86400e-3_rkx , 2.96961e-5_rkx , &
    1.02920e+0_rkx , 1.01680e-2_rkx , 5.30226e-5_rkx , &
    1.02743e+0_rkx , 9.85113e-3_rkx , 5.00233e-5_rkx ], [3,4])

  real(rkx) , dimension(4,4) , parameter :: coefb = reshape([ &
    8.85675e+0_rkx , -3.51620e-2_rkx ,  2.38653e-4_rkx , -1.71439e-6_rkx , &
    5.73841e+0_rkx , -1.91919e-2_rkx ,  1.65993e-4_rkx , -1.54665e-6_rkx , &
    6.64034e+0_rkx ,  1.56651e-2_rkx , -9.73357e-5_rkx ,  0.00000e+0_rkx , &
    7.09281e+0_rkx ,  1.40056e-2_rkx , -1.15774e-4_rkx ,  0.00000e+0_rkx], &
   [4,4])

  real(rkx) , dimension(3,4) , parameter :: coefc = reshape([ &
    9.90127e-1_rkx , 1.22475e-3_rkx , 4.90135e-6_rkx , &
    9.89753e-1_rkx , 1.97081e-3_rkx , 3.42046e-6_rkx , &
    9.75230e-1_rkx , 1.03341e-3_rkx , 0.00000e+0_rkx , &
    9.77366e-1_rkx , 8.60014e-4_rkx , 0.00000e+0_rkx],[3,4])

  real(rkx) , dimension(4,4) , parameter :: coefd = reshape([ &
    7.03047e-1_rkx , -2.63501e-3_rkx , -1.57023e-6_rkx ,  0.00000e+0_rkx , &
    5.29269e-1_rkx , -3.14754e-3_rkx ,  4.39595e-6_rkx ,  0.00000e+0_rkx , &
    7.88193e-2_rkx ,  1.31290e-3_rkx ,  4.25827e-6_rkx , -1.23982e-8_rkx , &
    1.62744e-1_rkx ,  2.22847e-3_rkx ,  2.60102e-6_rkx , -4.30133e-8_rkx], &
   [4,4])

  real(rkx) , dimension(3,4) , parameter :: coefe = reshape([ &
    3.93137e-2_rkx , -4.34341e-5_rkx , 3.74545e-8_rkx , &
    3.67785e-2_rkx , -3.10794e-5_rkx , 2.94436e-8_rkx , &
    7.42500e-2_rkx ,  3.97397e-5_rkx , 0.00000e+0_rkx , &
    7.52859e-2_rkx ,  4.18073e-5_rkx , 0.00000e+0_rkx], [3,4])

  real(rkx) , dimension(6,2) , parameter :: coeff = reshape([ &
    2.20370e-1_rkx , 1.39719e-3_rkx , -7.32011e-6_rkx ,   &
   -1.40262e-8_rkx , 2.13638e-10_rkx, -2.35955e-13_rkx ,  &
    3.07431e-1_rkx , 8.27225e-4_rkx , -1.30067e-5_rkx ,   &
    3.49847e-8_rkx , 2.07835e-10_rkx, -1.98937e-12_rkx], [6,2])

  real(rkx) , dimension(2,4) , parameter :: coefg = reshape([ &
    9.04489e+0_rkx , -9.56499e-3_rkx ,  1.80898e+1_rkx , &
   -1.91300e-2_rkx ,  8.72239e+0_rkx , -9.53359e-3_rkx , &
    1.74448e+1_rkx , -1.90672e-2_rkx],[2,4])

  real(rkx) , dimension(2,4) , parameter :: coefh = reshape([ &
    5.46557e+1_rkx , -7.30387e-2_rkx ,  1.09311e+2_rkx ,  &
   -1.46077e-1_rkx ,  5.11479e+1_rkx , -6.82615e-2_rkx ,  &
    1.02296e+2_rkx , -1.36523e-1_rkx],[2,4])

  real(rkx) , dimension(6,2) , parameter :: coefi = reshape([ &
    3.31654e-1_rkx , -2.86103e-4_rkx , -7.87860e-6_rkx ,   &
    5.88187e-8_rkx , -1.25340e-10_rkx , -1.37731e-12_rkx , &
    3.14365e-1_rkx , -1.33872e-3_rkx , -2.15585e-6_rkx ,   &
    6.07798e-8_rkx , -3.45612e-10_rkx , -9.34139e-15_rkx],[6,2])

  real(rkx) , dimension(3,2) , parameter :: coefj = reshape([ &
    2.82096e-2_rkx , 2.47836e-4_rkx , 1.16904e-6_rkx , &
    9.27379e-2_rkx , 8.04454e-4_rkx , 6.88844e-6_rkx],[3,2])

  real(rkx) , dimension(3,2) , parameter :: coefk = reshape([ &
    2.48852e-1_rkx , 2.09667e-3_rkx , 2.60377e-6_rkx , &
    1.03594e+0_rkx , 6.58620e-3_rkx , 4.04456e-6_rkx],[3,2])
  !
  ! Narrow band data for H2O
  ! 200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
  !
  real(rkx) , dimension(2) , parameter :: realk = [ &
       0.18967069430426e-4_rkx ,  0.70172244841851e-4_rkx ]
  real(rkx) , dimension(2) , parameter :: st = [ &
       0.31930234492350e-3_rkx ,  0.97907319939060e-3_rkx ]
  real(rkx) , dimension(2) , parameter :: a1 = [ &
       0.28775403075736e-1_rkx ,  0.23236701470511e-1_rkx ]
  real(rkx) , dimension(2) , parameter :: a2 = [ &
      -0.57966222388131e-4_rkx , -0.95105504388411e-4_rkx ]
  real(rkx) , dimension(2) , parameter :: b1 = [ &
       0.29927771523756e-1_rkx ,  0.21737073577293e-1_rkx ]
  real(rkx) , dimension(2) , parameter :: b2 = [ &
      -0.86322071248593e-4_rkx , -0.78543550629536e-4_rkx ]
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
  real(rkx) , dimension(4) , parameter :: abarl = [ &
      2.817e-2_rkx ,  2.682e-2_rkx , 2.264e-2_rkx , 1.281e-2_rkx ]
  real(rkx) , dimension(4) , parameter :: bbarl = [ &
      1.305e+0_rkx ,  1.346e+0_rkx , 1.454e+0_rkx , 1.641e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: cbarl = [ &
     -5.620e-8_rkx , -6.940e-6_rkx , 4.640e-4_rkx , 0.201e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: dbarl = [ &
      1.630e-8_rkx ,  2.350e-5_rkx , 1.240e-3_rkx , 7.560e-3_rkx ]
  real(rkx) , dimension(4) , parameter :: ebarl = [ &
      0.829e+0_rkx ,  0.794e+0_rkx , 0.754e+0_rkx , 0.826e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: fbarl = [ &
      2.482e-3_rkx ,  4.226e-3_rkx , 6.560e-3_rkx , 4.353e-3_rkx ]
  real(rkx) , dimension(4) , parameter :: abari = [ &
      3.4480e-3_rkx , 3.4480e-3_rkx , 3.4480e-3_rkx , 3.44800e-3_rkx ]
  real(rkx) , dimension(4) , parameter :: bbari = [ &
      2.4310e+0_rkx , 2.4310e+0_rkx , 2.4310e+0_rkx , 2.43100e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: cbari = [ &
      1.0000e-5_rkx , 1.1000e-4_rkx , 1.8610e-2_rkx , 0.46658e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: dbari = [ &
      0.0000e+0_rkx , 1.4050e-5_rkx , 8.3280e-4_rkx , 2.05000e-5_rkx ]
  real(rkx) , dimension(4) , parameter :: ebari = [ &
      0.7661e+0_rkx , 0.7730e+0_rkx , 0.7940e+0_rkx , 0.95950e+0_rkx ]
  real(rkx) , dimension(4) , parameter :: fbari = [ &
      5.8510e-4_rkx , 5.6650e-4_rkx , 7.2670e-4_rkx , 1.07600e-4_rkx ]
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
  real(rkx) , dimension(nspi) , parameter :: frcsol = [           &
      0.001488_rkx , 0.001389_rkx , 0.001290_rkx , 0.001686_rkx , &
      0.002877_rkx , 0.003869_rkx , 0.026336_rkx , 0.360739_rkx , &
      0.065392_rkx , 0.526861_rkx , 0.526861_rkx , 0.526861_rkx , &
      0.526861_rkx , 0.526861_rkx , 0.526861_rkx , 0.526861_rkx , &
      0.006239_rkx , 0.001834_rkx , 0.001834_rkx ]
  !
  ! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
  !
  real(rkx) , dimension(nspi) , parameter :: nirwgt = [           &
      0.000000_rkx , 0.000000_rkx , 0.000000_rkx , 0.000000_rkx , &
      0.000000_rkx , 0.000000_rkx , 0.000000_rkx , 0.000000_rkx , &
      0.320518_rkx , 1.000000_rkx , 1.000000_rkx , 1.000000_rkx , &
      1.000000_rkx , 1.000000_rkx , 1.000000_rkx , 1.000000_rkx , &
      1.000000_rkx , 1.000000_rkx , 1.000000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: raytau = [     &
      4.0200_rkx , 2.1800_rkx , 1.7000_rkx , 1.4500_rkx ,   &
      1.2500_rkx , 1.0850_rkx , 0.7300_rkx , 0.155208_rkx , &
      0.0392_rkx , 0.0200_rkx , 0.0001_rkx , 0.0001_rkx ,   &
      0.0001_rkx , 0.0001_rkx , 0.0001_rkx , 0.0001_rkx ,   &
      0.0001_rkx , 0.0001_rkx , 0.0001_rkx ]
  !
  ! Absorption coefficients
  !
  real(rkx) , dimension(nspi) , parameter :: abh2o = [             &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ,  0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ,  0.002_rkx , &
      0.035_rkx , 0.377_rkx , 1.950_rkx , 9.400_rkx , 44.600_rkx , &
    190.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: abo3 = [                  &
      5.370e+4_rkx , 13.080e+4_rkx , 9.292e+4_rkx , 4.530e+4_rkx ,     &
      1.616e+4_rkx ,  4.441e+3_rkx , 1.775e+2_rkx , 2.4058030e+1_rkx , &
      2.210e+1_rkx ,  0.000e+0_rkx , 0.000e+0_rkx , 0.000e+0_rkx ,     &
      0.000e+0_rkx ,  0.000e+0_rkx , 0.000e+0_rkx , 0.000e+0_rkx ,     &
      0.000e+0_rkx ,  0.000e+0_rkx , 0.000e+0_rkx ]

  real(rkx) , dimension(nspi) , parameter :: abco2 = [            &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.094_rkx , 0.196_rkx , 1.963_rkx ]

  real(rkx) , dimension(nspi) , parameter :: abo2 = [                 &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ,  0.000_rkx ,    &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 1.11e-5_rkx , 6.69e-5_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ,  0.000_rkx ,    &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ]
  !
  ! Spectral interval weights
  !
  real(rkx) , dimension(nspi) , parameter :: ph2o = [             &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.505_rkx , &
      0.210_rkx , 0.120_rkx , 0.070_rkx , 0.048_rkx , 0.029_rkx , &
      0.018_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: pco2 = [             &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 1.000_rkx , 0.640_rkx , 0.360_rkx ]

  real(rkx) , dimension(nspi) , parameter :: po2 = [             &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 1.000_rkx , 1.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx , &
      0.000_rkx , 0.000_rkx , 0.000_rkx , 0.000_rkx ]
  !
  ! Planck function factors - used in trcplk
  !
  real(rkx) , dimension(nlwspi) , parameter :: f1 = &
     [ 5.85713e8_rkx , 7.94950e8_rkx , 1.47009e9_rkx , 1.40031e9_rkx , &
       1.34853e8_rkx , 1.05158e9_rkx , 3.35370e8_rkx , 3.99601e8_rkx , &
       5.35994e8_rkx , 8.42955e8_rkx , 4.63682e8_rkx , 5.18944e8_rkx , &
       8.83202e8_rkx , 1.03279e9_rkx ]
  real(rkx) , dimension(nlwspi) , parameter :: f2 = &
     [ 2.02493e11_rkx , 3.04286e11_rkx , 6.90698e11_rkx , &
       6.47333e11_rkx , 2.85744e10_rkx , 4.41862e11_rkx , &
       9.62780e10_rkx , 1.21618e11_rkx , 1.79905e11_rkx , &
       3.29029e11_rkx , 1.48294e11_rkx , 1.72315e11_rkx , &
       3.50140e11_rkx , 4.31364e11_rkx ]
  real(rkx) , dimension(nlwspi) , parameter :: f3 = &
     [ 1383.0_rkx , 1531.0_rkx , 1879.0_rkx , 1849.0_rkx ,  848.0_rkx , &
       1681.0_rkx , 1148.0_rkx , 1217.0_rkx , 1343.0_rkx , 1561.0_rkx , &
       1279.0_rkx , 1328.0_rkx , 1586.0_rkx , 1671.0_rkx ]
  !
  ! Coefficients for h2o emissivity and absorptivity.
  !
  ! c1(iband) = coefe(3,iband)/coefe(2,iband)
  real(rkx) , dimension(4) , parameter :: c1 = &
    [ coefe(3,1)/coefe(2,1) , coefe(3,2)/coefe(2,2) , &
      coefe(3,3)/coefe(2,3) , coefe(3,4)/coefe(2,4) ]
  ! c2(iband) = coefb(3,iband)/coefb(2,iband)
  real(rkx) , dimension(4) , parameter :: c2 = &
    [ coefb(3,1)/coefb(2,1) , coefb(3,2)/coefb(2,2) , &
      coefb(3,3)/coefb(2,3) , coefb(3,4)/coefb(2,4) ]
  ! c3(iband) = coefb(4,iband)/coefb(3,iband)
  real(rkx) , dimension(4) , parameter :: c3 = &
    [ coefb(4,1)/coefb(3,1) , coefb(4,2)/coefb(3,2) , &
      coefb(4,3)/coefb(3,3) , coefb(4,4)/coefb(3,4) ]
  ! c4(iband) = coefd(3,iband)/coefd(2,iband)
  real(rkx) , dimension(4) , parameter :: c4 = &
    [ coefd(3,1)/coefd(2,1) , coefd(3,2)/coefd(2,2) , &
      coefd(3,3)/coefd(2,3) , coefd(3,4)/coefd(2,4) ]
  ! c5(iband) = coefd(4,iband)/coefd(3,iband)
  real(rkx) , dimension(4) , parameter :: c5 = &
    [ coefd(4,1)/coefd(3,1) , coefd(4,2)/coefd(3,2) , &
      coefd(4,3)/coefd(3,3) , coefd(4,4)/coefd(3,4) ]
  ! c6(iband) = coefa(3,iband)/coefa(2,iband)
  real(rkx) , dimension(4) , parameter :: c6 = &
    [ coefa(3,1)/coefa(2,1) , coefa(3,2)/coefa(2,2) , &
      coefa(3,3)/coefa(2,3) , coefa(3,4)/coefa(2,4) ]
  ! c7(iband) = coefc(3,iband)/coefc(2,iband)
  real(rkx) , dimension(4) , parameter :: c7 = &
    [ coefc(3,1)/coefc(2,1) , coefc(3,2)/coefc(2,2) , &
      coefc(3,3)/coefc(2,3) , coefc(3,4)/coefc(2,4) ]
  real(rkx) , parameter :: c8 = coeff(3,1)/coeff(2,1)
  real(rkx) , parameter :: c9 = coeff(3,2)/coeff(2,2)
  real(rkx) , parameter :: c10 = coeff(4,1)/coeff(3,1)
  real(rkx) , parameter :: c11 = coeff(4,2)/coeff(3,2)
  real(rkx) , parameter :: c12 = coeff(5,1)/coeff(4,1)
  real(rkx) , parameter :: c13 = coeff(5,2)/coeff(4,2)
  real(rkx) , parameter :: c14 = coeff(6,1)/coeff(5,1)
  real(rkx) , parameter :: c15 = coeff(6,2)/coeff(5,2)
  real(rkx) , parameter :: c16 = coefj(3,1)/coefj(2,1)
  real(rkx) , parameter :: c17 = coefk(3,1)/coefk(2,1)
  real(rkx) , parameter :: c18 = coefi(3,1)/coefi(2,1)
  real(rkx) , parameter :: c19 = coefi(3,2)/coefi(2,2)
  real(rkx) , parameter :: c20 = coefi(4,1)/coefi(3,1)
  real(rkx) , parameter :: c21 = coefi(4,2)/coefi(3,2)
  real(rkx) , parameter :: c22 = coefi(5,1)/coefi(4,1)
  real(rkx) , parameter :: c23 = coefi(5,2)/coefi(4,2)
  real(rkx) , parameter :: c24 = coefi(6,1)/coefi(5,1)
  real(rkx) , parameter :: c25 = coefi(6,2)/coefi(5,2)
  real(rkx) , parameter :: c26 = coefj(3,2)/coefj(2,2)
  real(rkx) , parameter :: c27 = coefk(3,2)/coefk(2,2)
  real(rkx) , parameter :: c28 = d_half
  real(rkx) , parameter :: c29 = 0.002053_rkx
  real(rkx) , parameter :: c30 = 0.1_rkx
  real(rkx) , parameter :: c31 = 3.0e-5_rkx
  real(rkx) , parameter :: cfa1 = 0.61_rkx

  logical :: luse_max_rnovl = .true.

#ifdef SINGLE_PRECISION_REAL
  real(rkx) , parameter :: mxarg = 16.0_rkx
#else
  real(rkx) , parameter :: mxarg = 25.0_rkx
#endif

  contains

  subroutine allocate_mod_rad_radiation
    implicit none

    npoints = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem3d(absnxt,1,kz,1,4,1,npoints,'rad:absnxt')
    call getmem3d(xuinpl,1,kz,1,4,1,npoints,'rad:xuinpl')
    call getmem3d(abstot,1,kzp1,1,kzp1,1,npoints,'rad:abstot')
    call getmem2d(emstot,1,kzp1,1,npoints,'rad:emstot')
    call getmem1d(diralb,1,npoints,'rad:diralb')
    call getmem1d(difalb,1,npoints,'rad:difalb')
    call getmem1d(tauaer,1,npoints,'rad:tauaer')
    call getmem1d(tauasc,1,npoints,'rad:tauasc')
    call getmem1d(gtota,1,npoints,'rad:gtota')
    call getmem1d(ftota,1,npoints,'rad:ftota')
    call getmem1d(skip,1,npoints,'rad:skip')
    call getmem1d(khiv,1,npoints,'rad:khiv')
    call getmem1d(khivm,1,npoints,'rad:khivm')
    call getmem1d(klov,1,npoints,'rad:klov')
    call getmem1d(tplnke,1,npoints,'rad:tplnke')
    call getmem1d(done,1,npoints,'rad:done')
    call getmem1d(co2mmr,1,npoints,'rad:co2mmr')
    call getmem1d(co2vmr,1,npoints,'rad:co2vmr')
    call getmem1d(n2ommr,1,npoints,'rad:n2ommr')
    call getmem1d(ch4mmr,1,npoints,'rad:ch4mmr')
    call getmem1d(cfc11mmr,1,npoints,'rad:cfc11mmr')
    call getmem1d(cfc12mmr,1,npoints,'rad:cfc12mmr')
    call getmem2d(bch4,1,kzp1,1,npoints,'rad:bch4')
    call getmem2d(bn2o0,1,kzp1,1,npoints,'rad:bn2o0')
    call getmem2d(bn2o1,1,kzp1,1,npoints,'rad:bn2o1')
    call getmem2d(co2em,1,kzp1,1,npoints,'rad:co2em')
    call getmem2d(co2eml,1,kzp1,1,npoints,'rad:co2eml')
    call getmem2d(co2t,1,kzp1,1,npoints,'rad:co2t')
    call getmem2d(h2otr,1,kzp1,1,npoints,'rad:h2otr')
    call getmem2d(ucfc11,1,kzp1,1,npoints,'rad:ucfc11')
    call getmem2d(ucfc12,1,kzp1,1,npoints,'rad:ucfc12')
    call getmem2d(un2o0,1,kzp1,1,npoints,'rad:un2o0')
    call getmem2d(un2o1,1,kzp1,1,npoints,'rad:un2o1')
    call getmem2d(uch4,1,kzp1,1,npoints,'rad:uch4')
    call getmem2d(uco211,1,kzp1,1,npoints,'rad:uco211')
    call getmem2d(uco212,1,kzp1,1,npoints,'rad:uco212')
    call getmem2d(uco213,1,kzp1,1,npoints,'rad:uco213')
    call getmem2d(uco221,1,kzp1,1,npoints,'rad:uco221')
    call getmem2d(uco222,1,kzp1,1,npoints,'rad:uco222')
    call getmem2d(uco223,1,kzp1,1,npoints,'rad:uco223')
    call getmem2d(uptype,1,kzp1,1,npoints,'rad:uptype')
    call getmem2d(plol,1,kzp1,1,npoints,'rad:plol')
    call getmem2d(plos,1,kzp1,1,npoints,'rad:plos')
    call getmem2d(tplnka,1,kzp1,1,npoints,'rad:tplnka')
    call getmem2d(tint,1,kzp1,1,npoints,'rad:tint')
    call getmem2d(tint4,1,kzp1,1,npoints,'rad:tint4')
    call getmem2d(tlayr,1,kzp1,1,npoints,'rad:tlayr')
    call getmem2d(tlayr4,1,kzp1,1,npoints,'rad:tlayr4')
    call getmem2d(wh2op,1,kzp1,1,npoints,'rad:w')
    call getmem2d(s2c,1,kzp1,1,npoints,'rad:s2c')
    call getmem2d(s2t,1,kzp1,1,npoints,'rad:s2t')
    call getmem2d(fdl,1,kzp1,1,npoints,'rad:fdl')
    call getmem2d(fsdl,1,kzp1,1,npoints,'rad:fsdl')
    call getmem2d(ful,1,kzp1,1,npoints,'rad:ful')
    call getmem2d(fsul,1,kzp1,1,npoints,'rad:fsul')
    call getmem2d(fdl0,1,kzp1,1,npoints,'rad:fdl0')
    call getmem2d(fsdl0,1,kzp1,1,npoints,'rad:fsdl0')
    call getmem2d(ful0,1,kzp1,1,npoints,'rad:ful0')
    call getmem2d(fsul0,1,kzp1,1,npoints,'rad:fsul0')
    call getmem2d(rtclrsf,1,kzp1,1,npoints,'rad:rtclrsf')
    call getmem2d(aermmb,1,kz,1,npoints,'rad:aermmb')

    call getmem2d(fclb4,1,kz,1,npoints,'rad:fclb4')
    call getmem2d(fclt4,1,kz,1,npoints,'rad:fclt4')

    call getmem2d(emplnk,1,nlwspi,1,npoints,'rad:emplnk')

    call getmem3d(abplnk1,1,nlwspi,1,kzp1,1,npoints,'rad:abplnk1')
    call getmem3d(abplnk2,1,nlwspi,1,kzp1,1,npoints,'rad:abplnk2')

    call getmem2d(exptdn,0,kzp1,1,npoints,'rad:exptdn')
    call getmem2d(fluxdn,0,kzp1,1,npoints,'rad:fluxdn')
    call getmem2d(fluxup,0,kzp1,1,npoints,'rad:fluxup')
    call getmem2d(rdndif,0,kzp1,1,npoints,'rad:rdndif')
    call getmem2d(rupdif,0,kzp1,1,npoints,'rad:rupdif')
    call getmem2d(rupdir,0,kzp1,1,npoints,'rad:rupdir')
    call getmem2d(tottrn,0,kzp1,1,npoints,'rad:tottrn')
    call getmem2d(pflx,0,kzp1,1,npoints,'rad:pflx')
    call getmem2d(fswup,0,kzp1,1,npoints,'rad:fswup')
    call getmem2d(fswdn,0,kzp1,1,npoints,'rad:fswdn')
    call getmem2d(plco2,1,kzp1,1,npoints,'rad:plco2')
    call getmem2d(plh2o,1,kzp1,1,npoints,'rad:plh2o')
    call getmem2d(tclrsf,1,kzp1,1,npoints,'rad:tclrsf')

    call getmem2d(rdir,0,kz,1,npoints,'rad:rdir')
    call getmem2d(rdif,0,kz,1,npoints,'rad:rdif')
    call getmem2d(tdir,0,kz,1,npoints,'rad:tdir')
    call getmem2d(tdif,0,kz,1,npoints,'rad:tdif')
    call getmem2d(explay,0,kz,1,npoints,'rad:explay')
    call getmem2d(flxdiv,0,kz,1,npoints,'rad:flxdiv')
    call getmem2d(totfld,0,kz,1,npoints,'rad:totfld')
    call getmem2d(wcl,0,kz,1,npoints,'rad:wcl')
    call getmem2d(gcl,0,kz,1,npoints,'rad:gcl')
    call getmem2d(fcl,0,kz,1,npoints,'rad:fcl')
    call getmem2d(wci,0,kz,1,npoints,'rad:wci')
    call getmem2d(gci,0,kz,1,npoints,'rad:gci')
    call getmem2d(fci,0,kz,1,npoints,'rad:fci')
    call getmem2d(uh2o,0,kz,1,npoints,'rad:uh2o')
    call getmem2d(uo3,0,kz,1,npoints,'rad:uo3')
    call getmem2d(uco2,0,kz,1,npoints,'rad:uco2')
    call getmem2d(uo2,0,kz,1,npoints,'rad:uo2')

    call getmem3d(fis,1,kzp1,1,kzp1,1,npoints,'rad:fis')
    call getmem3d(fis0,1,kzp1,1,kzp1,1,npoints,'rad:fis0')

    call getmem1d(solflx,1,npoints,'rad:solflx')
    call getmem1d(utco2,1,npoints,'rad:utco2')
    call getmem1d(uth2o,1,npoints,'rad:uth2o')
    call getmem1d(uto2,1,npoints,'rad:uto2')
    call getmem1d(uto3,1,npoints,'rad:uto3')
    call getmem1d(toafsnsc,1,npoints,'rad:toafsnsc')
    call getmem1d(toafsntc,1,npoints,'rad:toafsntc')
    call getmem1d(fslwdcs,1,npoints,'rad:fslwdcs')

    call getmem2d(cfc11,1,kz,1,npoints,'rad:cfc11')
    call getmem2d(cfc12,1,kz,1,npoints,'rad:cfc12')
    call getmem2d(ch4,1,kz,1,npoints,'rad:ch4')
    call getmem2d(n2o,1,kz,1,npoints,'rad:n2o')
    call getmem2d(pbr,1,kz,1,npoints,'rad:pbr')
    call getmem2d(pnm,1,kzp1,1,npoints,'rad:pnm')
    call getmem2d(o3mmr,1,kz,1,npoints,'rad:o3mmr')

  end subroutine allocate_mod_rad_radiation
  !
  !-----------------------------------------------------------------------
  !
  ! Initialize various constants for radiation scheme; note that
  ! the radiation scheme uses cgs units.
  !
  !-----------------------------------------------------------------------
  !
  ! RegCM : Most constants moved in mod_constants, GHG gases from CMIP
  !
  subroutine radini(n1,n2,iyear,imonth,lat, &
                    co2vmr,co2mmr,ch4mmr,n2ommr,cfc11mmr,cfc12mmr)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2 , iyear , imonth
    real(rkx) , dimension(n1:n2) , intent(in) :: lat
    real(rkx) , dimension(n1:n2) , intent(out) :: co2vmr , co2mmr
    real(rkx) , dimension(n1:n2) , intent(out) :: ch4mmr , n2ommr
    real(rkx) , dimension(n1:n2) , intent(out) :: cfc11mmr , cfc12mmr
    integer(ik4) :: n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radini'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Set general radiation consts; convert to cgs units where
    ! appropriate.
    !
    ! Cannot be done in parallel because it requires I/O from file
    !
    do n = n1 , n2
      co2vmr(n) = ghgval(igh_co2,iyear,imonth,lat(n))
      co2mmr(n) = co2vmr(n)*(amco2/amd)
      ch4mmr(n) = ghgval(igh_ch4,iyear,imonth,lat(n))*(amch4/amd)
      n2ommr(n) = ghgval(igh_n2o,iyear,imonth,lat(n))*(amn2o/amd)
      cfc11mmr(n) = ghgval(igh_cfc11,iyear,imonth,lat(n))*(amcfc11/amd)
      cfc12mmr(n) = ghgval(igh_cfc12,iyear,imonth,lat(n))*(amcfc12/amd)
    end do
    !
    ! Set execution flag for aerosol and their interaction with radiation
    !
    lzero = .true.
    linteract = ( (ichem == 1 .and. idirect > 0) .or. iclimaaer > 0 )
    if ( ichem == 1 ) then
      if ( idirect == 2 ) lzero = .false.
    end if
    if ( iclimaaer > 0 ) then
      lzero = .false.
    end if
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  !-----------------------------------------------------------------------
  !
  subroutine radctl(rt,iyear,imonth)
    implicit none
    !
    ! Input arguments
    !
    ! ts      - Surface (skin) temperature
    ! pmid    - Model level pressures
    ! pint    - Model interface pressures
    ! pmln    - Natural log of pmid
    ! rel     - liquid cloud particle effective radius
    ! rei     - ice effective drop size (microns)
    ! fice    - fractional ice content within cloud
    ! piln    - Natural log of pint
    ! t       - Model level temperatures
    ! q       - Model level specific humidity
    ! cld     - Fractional cloud cover
    ! effcld  - Effective fractional cloud cover
    ! clwp    - Cloud liquid water path
    !
    ! Output solar arguments
    !
    ! fsns    - Surface absorbed solar flux
    ! sols    - Downward solar rad onto surface (sw direct)
    ! soll    - Downward solar rad onto surface (lw direct)
    ! solsd   - Downward solar rad onto surface (sw diffuse)
    ! solld   - Downward solar rad onto surface (lw diffuse)
    ! qrs     - Solar heating rate
    !
    ! qrl     - Longwave cooling rate
    ! flwds   - Surface down longwave flux
    ! solin    - Solar incident flux
    ! solout   - Solar outgoing flux
    ! fsnt     - Net column abs solar flux at model top
    ! fsntc    - Clear sky total column abs solar flux
    ! fsnsc    - Clear sky surface abs solar flux
    ! fsnirt   - Near-IR flux absorbed at toa
    ! fsnrtc   - Clear sky near-IR flux absorbed at toa
    ! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
    ! fsds     - Flux Shortwave Downwelling Surface
    ! flnt     - Net outgoing lw flux at model top
    ! lwout    - outgoing lw flux at model top
    ! lwin     - incoming lw flux at model top
    ! flns     - Srf longwave cooling (up-down) flux
    ! flntc    - Clear sky lw flux at model top
    ! flnsc    - Clear sky lw flux at srf (up-down)
    ! o3vmr    - Ozone volume mixing ratio
    ! eccf     - Earth/sun distance factor
    !
    integer(ik4) , intent(in) :: iyear , imonth
    type(radtype) , intent(inout) :: rt
    integer(ik4) :: n
    integer(ik4) :: k
#ifndef RCEMIP
    real(rkx) :: alat
#endif
    real(rkx) :: pratio , xcfc11 , xcfc12 , xch4 , xn2o , betafac

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radctl'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Set latitude dependent radiation input: radini sets many
    ! radiation parameters
    !
    call radini(rt%n1,rt%n2,iyear,imonth,rt%dlat, &
                co2vmr,co2mmr,ch4mmr,n2ommr,cfc11mmr,cfc12mmr)
    !
    ! Compute on interface levels and put in vertical profiles
    !
    call radinp(rt%n1,rt%n2,rt%pmid,rt%pint,rt%q,rt%cld,rt%o3vmr, &
                pbr,pnm,plco2,plh2o,tclrsf,o3mmr)
    !
    ! Solar radiation computation
    !
    if ( dosw ) then
      !
      ! Specify aerosol mass mixing ratio
      !
      call aermix(rt%n1,rt%n2,rt%pint,aermmb)

      call aeroppt(rt%rh,rt%pint,rt%n1,rt%n2)

      call radcsw(rt%n1,rt%n2,rt%eccf,pnm,rt%q,o3mmr,aermmb,rt%cld,   &
                  rt%clwp,rt%rel,rt%rei,rt%fice,rt%czen,rt%czengt0,   &
                  rt%adirsw,rt%adifsw,rt%adirlw,rt%adiflw,rt%asw,     &
                  rt%alw,rt%solin,rt%solout,rt%qrs,rt%fsns,rt%fsnt,   &
                  rt%fsds,rt%fsnsc,rt%fsntc,rt%sols,rt%soll,rt%solsd, &
                  rt%solld,rt%fsnirt,rt%fsnrtc,rt%fsnirtsq,rt%abv,    &
                  rt%sol,rt%aeradfo,rt%aeradfos,rt%tauxcl,rt%tauxci,  &
                  rt%outtaucl,rt%outtauci)
      !
      ! Convert units of shortwave fields needed by rest of model
      ! from CGS to MKS
      !
#ifdef STDPAR
      do concurrent ( n = rt%n1:rt%n2 ) local(k,betafac)
#else
      do n = rt%n1 , rt%n2
#endif
        rt%solin(n) = rt%solin(n)*1.0e-3_rkx
        rt%solout(n) = rt%solout(n)*1.0e-3_rkx
        rt%fsnt(n) = rt%fsnt(n)*1.0e-3_rkx
        rt%fsns(n) = rt%fsns(n)*1.0e-3_rkx
        rt%fsntc(n) = rt%fsntc(n)*1.0e-3_rkx
        rt%fsnsc(n) = rt%fsnsc(n)*1.0e-3_rkx
        !
        ! clear sky column partitioning for surface flux
        ! note : should be generalised to the whole column to be
        !        really in energy balance !
        !
        rt%totcf(n) = d_one
        if ( luse_max_rnovl ) then
          do k = 2 , kzp1
            rt%totcf(n) = rt%totcf(n) * &
                   (1.0001_rkx - max(rt%cld(k-1,n),rt%cld(k,n)))/ &
                   (1.0001_rkx - rt%cld(k-1,n))
          end do
        else
          do k = 1 , kzp1
            rt%totcf(n) = rt%totcf(n) * (d_one - rt%cld(k,n))
          end do
        end if
        rt%totcf(n) = d_one - rt%totcf(n)
        !
        ! maximum cld cover considered
        ! rt%fsns(n) = rt%fsns(n) * maxval(rt%cld(:,n)) + &
        !           rt%fsnsc(n) * (1-maxval(rt%cld(:,n)))
        ! random overlap assumption is tocf(n)
        ! Now average btw rand ov and maximum cloud cover as fil suggest
        ! rt%totcf(n) =  d_half * ( rt%totcf(n) + maxval(rt%cld(:,n)) )
        ! abv is proportional to fsns in radcsw : Calculate the factor
        if ( rt%fsns(n) > d_zero ) then
          betafac = rt%abv(n) / rt%fsns(n)
        else
          betafac = d_zero
        end if
        ! Fil suggestion of putting a max on column cloud fraction
        ! TAO: implement a user-specified CF maximum (default of 1.0)
        if ( lsrfhack ) then
          if ( rt%totcf(n) > cftotmax ) rt%totcf(n) = cftotmax
          if ( rt%totcf(n) < d_zero ) rt%totcf(n) = d_zero
          rt%fsns(n) = rt%fsns(n) * rt%totcf(n) + &
                       rt%fsnsc(n) * (d_one-rt%totcf(n))
        end if
        ! Apply the clear-sky / cloudy-sky also to abv using the beta factor
        rt%abv(n) = betafac * rt%fsns(n)
        rt%fsds(n) = rt%fsds(n)*1.0e-3_rkx
        rt%fsnirt(n) = rt%fsnirt(n)*1.0e-3_rkx
        rt%fsnrtc(n) = rt%fsnrtc(n)*1.0e-3_rkx
        rt%fsnirtsq(n) = rt%fsnirtsq(n)*1.0e-3_rkx
        !
        ! Calculate/outfld albedo and clear sky albedo
        !
        if ( rt%solin(n) > d_zero ) then
          rt%alb(n) = (rt%solin(n)-rt%fsnt(n))/rt%solin(n)
          rt%albc(n) = (rt%solin(n)-rt%fsntc(n))/rt%solin(n)
        else
          rt%alb(n) = d_zero
          rt%albc(n) = d_zero
        end if
      end do
    end if
    !
    ! Longwave radiation computation
    !
    if ( dolw ) then
      !
      ! Specify trace gas mixing ratios
      !
#ifdef RCEMIP
      xn2o = 0.3478_rkx
      xch4 = 0.2353_rkx
      xcfc11 = 0.7273_rkx
      xcfc12 = 0.4000_rkx
#endif
#ifdef STDPAR
      do concurrent ( n = rt%n1:rt%n2 ) local(k)
#else
      do n = rt%n1 , rt%n2
#endif
#ifndef RCEMIP
        alat = abs(rt%dlat(n))
        if ( alat <= 45.0_rkx ) then
          xn2o = 0.3478_rkx + 0.00116_rkx*alat
          xch4 = 0.2353_rkx
          xcfc11 = 0.7273_rkx + 0.00606_rkx*alat
          xcfc12 = 0.4000_rkx + 0.00222_rkx*alat
        else
          xn2o = 0.4000_rkx + 0.013333_rkx*(alat-45.0_rkx)
          xch4 = 0.2353_rkx + 0.0225489_rkx*(alat-45.0_rkx)
          xcfc11 = 1.00_rkx + 0.013333_rkx*(alat-45.0_rkx)
          xcfc12 = 0.50_rkx + 0.024444_rkx*(alat-45.0_rkx)
        end if
#endif
        do k = 1 , kz
          !  set stratospheric scale height factor for gases
          if ( rt%pmid(k,n) >= rt%xptrop(n) ) then
            ch4(k,n) = ch4mmr(n)
            n2o(k,n) = n2ommr(n)
            cfc11(k,n) = cfc11mmr(n)
            cfc12(k,n) = cfc12mmr(n)
          else
            pratio = rt%pmid(k,n)/rt%xptrop(n)
            ch4(k,n) = ch4mmr(n)*(pratio**xch4)
            n2o(k,n) = n2ommr(n)*(pratio**xn2o)
            cfc11(k,n) = cfc11mmr(n)*(pratio**xcfc11)
            cfc12(k,n) = cfc12mmr(n)*(pratio**xcfc12)
          end if
        end do
      end do
      call radclw(rt%n1,rt%n2,rt%labsem,rt%ts,rt%emiss,rt%t,rt%q,   &
                  rt%o3vmr,pbr,pnm,rt%pmln,rt%piln,n2o,ch4,cfc11,   &
                  cfc12,rt%effcld,tclrsf,rt%flns,rt%flnt,rt%lwout,  &
                  rt%lwin,rt%flnsc,rt%flntc,rt%flwds,fslwdcs,       &
                  rt%aerlwfo,rt%aerlwfos,rt%absgasnxt,rt%absgastot, &
                  rt%emsgastot,rt%qrl)
      !
      ! Convert units of longwave fields needed by rest of model from CGS to MKS
      !
      do concurrent ( n = rt%n1:rt%n2 )
        rt%flnt(n) = rt%flnt(n)*1.0e-3_rkx
        rt%lwout(n) = rt%lwout(n)*1.0e-3_rkx
        rt%lwin(n) = rt%lwin(n)*1.0e-3_rkx
        rt%flns(n) = rt%flns(n)*1.0e-3_rkx
        rt%flntc(n) = rt%flntc(n)*1.0e-3_rkx
        rt%flnsc(n) = rt%flnsc(n)*1.0e-3_rkx
        rt%flwds(n) = rt%flwds(n)*1.0e-3_rkx
        fslwdcs(n) = fslwdcs(n)*1.0e-3_rkx
        !
        ! essai clear sky column
        !
        ! rt%flwds(n) = rt%flwds(n) * maxval(rt%cld((:,n))) + &
        !            rt%flwds(n) * (1-maxval(rt%cld((:,n))))
        ! rt%flwds(n) = rt%flwds(n) * maxval(rt%cld(:,n)) + &
        !            fslwdcs(n)*(d_one-maxval(rt%cld(:,n)))
        ! rt%flns(n) = rt%flns(n) * maxval(rt%cld(:,n)) + &
        !           rt%flnsc(n)*(d_one-maxval(rt%cld(:,n)))
        !
        ! rt%totcf(n) has been calculated for the SW, dolw is always true
        !
        if ( lsrfhack ) then
          rt%flwds(n) = rt%flwds(n) * rt%totcf(n) + &
                        fslwdcs(n) * (d_one - rt%totcf(n))
          rt%flns(n)  = rt%flns(n) * rt%totcf(n)  + &
                        rt%flnsc(n) * (d_one - rt%totcf(n))
        end if
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  ! Note that an extra layer above the model top layer is added.
  !
  ! cgs units are used.
  !
  ! Special diagnostic calculation of the clear sky surface and total column
  ! absorbed flux is also done for cloud forcing diagnostics.
  !
  ! Input arguments
  !
  ! pnm     - Interface pressure (dynes/cm2)
  ! h2ommr  - Specific humidity (h2o mass mix ratio)
  ! o3mmr   - Ozone mass mixing ratio
  ! cld     - Fractional cloud cover
  ! clwp    - Layer liquid water path
  ! rel     - Liquid effective drop size (microns)
  ! rei     - Ice effective drop size (microns)
  ! fice    - Fractional ice content within cloud
  ! eccf    - Eccentricity factor (d_one/earth-sun dist ** 2)
  !
  ! Output arguments
  !
  ! solin    - Incident solar flux
  ! solout   - Outgoing solar flux
  ! qrs      - Solar heating rate
  ! fsns     - Surface absorbed solar flux
  ! fsnt     - Total column absorbed solar flux
  ! fsds     - Flux Shortwave Downwelling Surface
  ! fsnsc    - Clear sky surface absorbed solar flux
  ! fsntc    - Clear sky total column absorbed solar flx
  ! sols     - Direct solar rad incident on surface (< 0.7)
  ! soll     - Direct solar rad incident on surface (>= 0.7)
  ! solsd    - Diffuse solar rad incident on surface (< 0.7)
  ! solld    - Diffuse solar rad incident on surface (>= 0.7)
  ! fsnirt   - Near-IR flux absorbed at toa
  ! fsnrtc   - Clear sky near-IR flux absorbed at toa
  ! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
  !
  subroutine radcsw(n1,n2,eccf,pnm,h2ommr,o3mmr,aermmb,cld,clwp,rel,rei,   &
                    fice,czen,czengt0,adirsw,adifsw,adirlw,adiflw,asw,alw, &
                    solin,solout,qrs,fsns,fsnt,fsds,fsnsc,fsntc,sols,soll, &
                    solsd,solld,fsnirt,fsnrtc,fsnirtsq,abv,sol,aeradfo,    &
                    aeradfos,tauxcl,tauxci,outtaucl,outtauci)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , intent(in) :: eccf
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pnm
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: o3mmr , h2ommr , aermmb
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: cld
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: clwp , fice , rel , rei
    real(rkx) , dimension(n1:n2) , intent(in) :: czen
    logical , dimension(n1:n2) , intent(in) :: czengt0
    real(rkx) , dimension(n1:n2) , intent(in) :: adirsw
    real(rkx) , dimension(n1:n2) , intent(in) :: adifsw
    real(rkx) , dimension(n1:n2) , intent(in) :: adirlw
    real(rkx) , dimension(n1:n2) , intent(in) :: adiflw
    real(rkx) , dimension(n1:n2) , intent(in) :: asw
    real(rkx) , dimension(n1:n2) , intent(in) :: alw
    real(rkx) , dimension(n1:n2) , intent(out) :: aeradfo
    real(rkx) , dimension(n1:n2) , intent(out) :: aeradfos
    real(rkx) , dimension(n1:n2) , intent(out) :: fsds
    real(rkx) , dimension(n1:n2) , intent(out) :: fsnirt
    real(rkx) , dimension(n1:n2) , intent(out) :: fsnirtsq
    real(rkx) , dimension(n1:n2) , intent(out) :: fsnrtc
    real(rkx) , dimension(n1:n2) , intent(out) :: fsns
    real(rkx) , dimension(n1:n2) , intent(out) :: fsnsc
    real(rkx) , dimension(n1:n2) , intent(out) :: fsnt
    real(rkx) , dimension(n1:n2) , intent(out) :: fsntc
    real(rkx) , dimension(n1:n2) , intent(out) :: solin
    real(rkx) , dimension(n1:n2) , intent(out) :: solout
    real(rkx) , dimension(n1:n2) , intent(out) :: soll
    real(rkx) , dimension(n1:n2) , intent(out) :: solld
    real(rkx) , dimension(n1:n2) , intent(out) :: sols
    real(rkx) , dimension(n1:n2) , intent(out) :: solsd
    real(rkx) , dimension(n1:n2) , intent(out) :: abv
    real(rkx) , dimension(n1:n2) , intent(out) :: sol
    real(rkx) , dimension(kzp1,4,n1:n2) , intent(out) ::  outtaucl , outtauci
    real(rkx) , dimension(0:kz,n1:n2,nspi) , intent(out) :: tauxcl , tauxci
    real(rkx) , dimension(kz,n1:n2) , intent(out) :: qrs
    !
    ! indxsl   - Index for cloud particle properties
    !
    ! A. Slingo's data for cloud particle radiative properties (from 'A
    ! GCM Parameterization for the Shortwave Properties of Water
    ! Clouds' JAS vol. 46 may 1989 pp 1419-1427)
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
    ! abarii   - A coefficient for current spectral interval
    ! bbarii   - B coefficient for current spectral interval
    ! cbarii   - C coefficient for current spectral interval
    ! dbarii   - D coefficient for current spectral interval
    ! ebarii   - E coefficient for current spectral interval
    ! fbarii   - F coefficient for current spectral interval
    ! wgtint   - Weight for specific spectral interval
    !
    ! Diagnostic and accumulation : note that sfltot, fswup, and
    ! fswdn are not used in the computation,but are retained for future
    ! use.
    !
    ! sfltot   - Spectrally summed total solar flux
    !
    ! Various arrays and other constants:
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
    ! zenfac   - Square root of cos solar zenith angle
    ! aeradfo  - spectrally integrated aerosol radiative forcing ( TOA)
    !-----------------------------------------------------------------------
    !
    real(rkx) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli , &
                 dbarii , dbarli , ebarii , ebarli , fbarii , fbarli , &
                 psf , trayoslp , wavmid , wgtint
    real(rkx) , dimension(4) :: ww
    integer(ik4) :: n , k , indxsl , ns , is
    real(rkx) , parameter :: tmp1 = d_half/(egravgts*sslp)
    real(rkx) , parameter :: tmp2 = delta*regravgts
    real(rkx) :: sqrco2 , xptop , pdel , path
    real(rkx) :: ptho2 , ptho3 , pthco2 , pthh2o , h2ostr
    real(rkx) :: tmp1l , tmp2l , tmp3l , tmp1i , tmp2i , tmp3i
    real(rkx) :: rdenom , zenfac
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radcsw'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Initialize output fields:
    !
    fsds(:) = d_zero
    fsnirt(:) = d_zero
    fsnrtc(:) = d_zero
    fsnirtsq(:) = d_zero
    fsnt(:) = d_zero
    fsns(:) = d_zero
    solin(:) = d_zero
    solout(:) = d_zero
    fsnsc(:) = d_zero
    fsntc(:) = d_zero
    sols(:) = d_zero
    soll(:) = d_zero
    solsd(:) = d_zero
    solld(:) = d_zero
    abv(:) = d_zero
    sol(:) = d_zero
    aeradfo(:) = d_zero
    aeradfos(:) = d_zero
    toafsntc(:) = d_zero
    toafsnsc(:) = d_zero
    outtaucl(:,:,:) = d_zero
    outtauci(:,:,:) = d_zero
    tauxcl(:,:,:) = d_zero
    tauxci(:,:,:) = d_zero
    ww(:) = d_zero
    qrs(:,:) = d_zero
    !
    ! Define solar incident radiation and interface pressures:
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) &
      local(sqrco2,xptop,pdel,path,ptho2,ptho3,pthco2,pthh2o, &
      h2ostr,zenfac,k)
#else
    do n = n1 , n2
#endif
      !
      ! Initialize spectrally integrated totals:
      !
      fswup(kzp1,n) = d_zero
      fswdn(kzp1,n) = d_zero
      if ( czengt0(n) ) then
        solin(n) = scon*eccf*czen(n)
        pflx(0,n) = d_zero
        do k = 1 , kzp1
          pflx(k,n) = pnm(k,n)
        end do
        !
        ! Compute optical paths:
        ! CO2, use old scheme(as constant)
        !
        ! co2mmr = co2vmr*(mmwco2/mmwair)
        sqrco2 = sqrt(co2mmr(n))
        xptop = pflx(1,n)
        ptho2 = o2mmr*xptop*regravgts
        ptho3 = o3mmr(1,n)*xptop*regravgts
        pthco2 = sqrco2*(xptop*regravgts)
        h2ostr = sqrt(d_one/h2ommr(1,n))
        zenfac = sqrt(czen(n))
        pthh2o = (xptop**2)*tmp1 + &
          (xptop*regravgts) * (h2ostr*zenfac*delta)
        uh2o(0,n) = h2ommr(1,n)*pthh2o
        uco2(0,n) = zenfac*pthco2
        uo2(0,n) = zenfac*ptho2
        uo3(0,n) = ptho3
        do k = 1 , kz
          pdel = pflx(k+1,n) - pflx(k,n)
          path = pdel*regravgts
          ptho2 = o2mmr*path
          ptho3 = o3mmr(k,n)*path
          pthco2 = sqrco2*path
          h2ostr = sqrt(d_one/h2ommr(k,n))
          pthh2o = (pflx(k+1,n)**2-pflx(k,n)**2) * &
                    tmp1 + pdel*h2ostr*zenfac*tmp2
          uh2o(k,n) = h2ommr(k,n)*pthh2o
          uco2(k,n) = zenfac*pthco2
          uo2(k,n) = zenfac*ptho2
          uo3(k,n) = ptho3
        end do
        !
        ! Compute column absorber amounts for the clear sky computation:
        !
        uth2o(n) = d_zero
        uto3(n) = d_zero
        utco2(n) = d_zero
        uto2(n) = d_zero
        do k = 1 , kz
          uth2o(n) = uth2o(n) + uh2o(k,n)
          uto3(n) = uto3(n) + uo3(k,n)
          utco2(n) = utco2(n) + uco2(k,n)
          uto2(n) = uto2(n) + uo2(k,n)
        end do
        do k = 0 , kz
          totfld(k,n) = d_zero
          fswup(k,n) = d_zero
          fswdn(k,n) = d_zero
        end do
        !
        ! Set cloud properties for top (0) layer; so long as tauxcl is zero,
        ! there is no cloud above top of model; the other cloud properties
        ! are arbitrary:
        !
        wcl(0,n) = verynearone
        gcl(0,n) = 0.850_rkx
        fcl(0,n) = 0.725_rkx
        wci(0,n) = verynearone
        gci(0,n) = 0.850_rkx
        fci(0,n) = 0.725_rkx
      end if
    end do
    !
    ! Begin spectral loop
    !
    do ns = 1 , nspi
      !
      ! Begin spectral loop
      !
      wgtint = nirwgt(ns)
      !
      ! Set index for cloud particle properties based on the wavelength,
      ! according to A. Slingo (1989) equations 1-3:
      ! Use index 1 (0.25 to 0.69 micrometers) for visible
      ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
      ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
      ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
      !
      ! Note that the minimum wavelength is encoded (with 0.001, 0.002,
      ! 0.003) in order to specify the index appropriate for the
      ! near-infrared cloud absorption properties
      !
      indxsl = 0
      if ( wavmax(ns) <= 0.70_rkx ) then
        indxsl = 1
      else if ( abs(wavmin(ns)-0.700_rkx) < dlowval ) then
        indxsl = 2
      else if ( abs(wavmin(ns)-0.701_rkx) < dlowval ) then
        indxsl = 3
      else if ( abs(wavmin(ns)-0.702_rkx) < dlowval .or. &
                    wavmin(ns) > 2.38_rkx ) then
        indxsl = 4
      end if
      !
      ! Set cloud extinction optical depth, single scatter albedo,
      ! asymmetry parameter, and forward scattered fraction:
      !
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)

      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)

      ww(indxsl) = ww(indxsl) + d_one
      !
      ! Set reflectivities for surface based on mid-point wavelength
      !
      wavmid = (wavmin(ns)+wavmax(ns))*d_half
#ifdef STDPAR
      do concurrent ( n = n1:n2 ) &
        local(tmp1l,tmp2l,tmp3l,tmp1i,tmp2i,tmp3i,k)
#else
      do n = n1 , n2
#endif
        if ( czengt0(n) ) then
          do k = 1 , kz
            !
            ! liquid
            !
            tmp1l = abarli + bbarli/rel(k,n)
            tmp2l = d_one - cbarli - dbarli*rel(k,n)
            tmp3l = fbarli*rel(k,n)
            !
            ! ice
            !
            tmp1i = abarii + bbarii/rei(k,n)
            tmp2i = d_one - cbarii - dbarii*rei(k,n)
            tmp3i = fbarii*rei(k,n)
            !
            !  Cloud fraction incorporated into cloud extinction optical depth
            !  found April 12 2000, Filippo found the different scheme here:
            !
            ! Scheme     1
            ! The one in ccm3.6.6
            !tauxcl(k,n,ns) = clwp(k,n) * tmp1l * &
            !          (d_one-fice(k,n)) * cld(k,n) * sqrt(cld(k,n))
            !tauxci(k,n,ns) = clwp(k,n) * tmp1i * &
            !           fice(k,n) * cld(k,n) * sqrt(cld(k,n))
            !
            ! Scheme     2
            ! unknown origin (?????)
            !tauxcl(k,n,ns) = ((clwp(k,n)*cld(k,n))* &
            !              (d_one-fice(k,n))*tmp1l) / &
            !              (d_one+(d_one-0.85_rkx)*((d_one-cld(k,n))*      &
            !              (clwp(k,n)*tmp1l*(d_one-fice(k,n)))))
            !tauxci(k,n,ns) = (clwp(k,n)*cld(k,n)*fice(k,n)*tmp1i) /  &
            !              (d_one+(d_one-0.78_rkx)*((d_one-cld(k,n)) * &
            !              (clwp(k,n)*tmp1i*fice(k,n))))
            !
            tauxcl(k,n,ns) = ((clwp(k,n)*cld(k,n)) * &
              (d_one-fice(k,n))*tmp1l) / &
              (d_one+(d_one-0.85_rkx)*((d_one-cld(k,n))*      &
              (clwp(k,n)*tmp1l*(d_one-fice(k,n)))))
            tauxci(k,n,ns) = (clwp(k,n)*cld(k,n)*fice(k,n)*tmp1i) /  &
                          (d_one+(d_one-0.78_rkx)*((d_one-cld(k,n)) * &
                          (clwp(k,n)*tmp1i*fice(k,n))))
            outtaucl(k,indxsl,n) = outtaucl(k,indxsl,n) + tauxcl(k,n,ns)
            outtauci(k,indxsl,n) = outtauci(k,indxsl,n) + tauxci(k,n,ns)
            !
            !scheme     3
            ! tauxcl(k,n,ns) = clwp(k,n)*tmp1l* &
            !           (d_one-fice(k,n))*cld(k,n)**0.85
            ! tauxci(k,n,ns) = clwp(k,n)*tmp1i*fice(k,n)*cld(k,n)**0.85
            !
            ! Do not let single scatter albedo be 1; delta-eddington
            ! solution for non-conservative case:
            !
            wcl(k,n) = min(tmp2l,verynearone)
            gcl(k,n) = ebarli + tmp3l
            fcl(k,n) = gcl(k,n)*gcl(k,n)

            wci(k,n) = min(tmp2i,verynearone)
            gci(k,n) = ebarii + tmp3i
            fci(k,n) = gci(k,n)*gci(k,n)
          end do
          if ( wavmid < 0.7_rkx ) then
            !
            ! Wavelength less  than 0.7 micro-meter
            !
            diralb(n) = adirsw(n)
            difalb(n) = adifsw(n)
          else
            !
            ! Wavelength greater than 0.7 micro-meter
            !
            diralb(n) = adirlw(n)
            difalb(n) = adiflw(n)
          end if
        end if
      end do

      trayoslp = raytau(ns)/sslp
      !
      ! Layer input properties now completely specified; compute the
      ! delta-Eddington solution reflectivities and transmissivities
      ! for each layer, starting from the top and working downwards:
      !
      ! options for aerosol: no climatic feedback if idirect == 1
      ! should be consistent with aeroppt routine
      !
      call radded(n1,n2,trayoslp,czen,czengt0,pflx,       &
                  abh2o(ns),abo3(ns),abco2(ns),abo2(ns),  &
                  uh2o,uo3,uco2,uo2,tauxcl(:,:,ns),       &
                  wcl,gcl,fcl,tauxci(:,:,ns),wci,gci,fci, &
                  tauxar3d(:,:,ns),tauasc3d(:,:,ns),      &
                  gtota3d(:,:,ns),ftota3d(:,:,ns),        &
                  tottrn,exptdn,rdndif,rdif,tdif,rdir,tdir,explay)
#ifdef STDPAR
      do concurrent ( n = n1:n2 ) local(rdenom,k)
#else
      do n = n1 , n2
#endif
        if ( czengt0(n) ) then
          rupdir(kzp1,n) = diralb(n)
          rupdif(kzp1,n) = difalb(n)
          !
          ! Compute reflectivity to direct and diffuse radiation for layers
          ! below by adding succesive layers starting from the surface and
          ! working upwards:
          !
          do k = kz , 0 , -1
            rdenom = d_one/(d_one-(rdif(k,n)*rupdif(k+1,n)))
            rupdir(k,n) = rdir(k,n) + tdif(k,n) *      &
                          (rupdir(k+1,n)*explay(k,n) + &
                           rupdif(k+1,n)*(tdir(k,n)-explay(k,n)))*rdenom
            rupdif(k,n) = rdif(k,n) + rupdif(k+1,n)*(tdif(k,n)**2)*rdenom
          end do
          !
          ! Compute up and down fluxes for each interface, using the added
          ! atmospheric layer properties at each interface:
          !
          do k = 0 , kzp1
            rdenom = d_one/(d_one-(rdndif(k,n)*rupdif(k,n)))
            fluxup(k,n) = (exptdn(k,n)*rupdir(k,n)+   &
                          (tottrn(k,n)-exptdn(k,n))*rupdif(k,n))*rdenom
            fluxdn(k,n) = exptdn(k,n) +                              &
                          (tottrn(k,n) - exptdn(k,n) + exptdn(k,n) * &
                          (rupdir(k,n)*rdndif(k,n)))*rdenom
          end do
          !
          ! Compute flux divergence in each layer using the interface up
          ! and down fluxes:
          !
          do k = 0 , kz
            flxdiv(k,n) = (fluxdn(k,n) - fluxdn(k+1,n)) + &
                          (fluxup(k+1,n) - fluxup(k,n))
          end do
        end if
      end do
      !
      ! Monochromatic computation completed; accumulate in totals;
      ! adjust fraction within spectral interval to allow for the
      ! possibility of sub-divisions within a particular interval:
      !
      psf = d_one
      if ( abs(ph2o(ns)) > dlowval ) psf = psf*ph2o(ns)
      if ( abs(pco2(ns)) > dlowval ) psf = psf*pco2(ns)
      if ( abs(po2(ns)) > dlowval ) psf = psf*po2(ns)
#ifdef STDPAR
      do concurrent ( n = n1:n2 ) local(k)
#else
      do n = n1 , n2
#endif
        if ( czengt0(n) ) then
          solflx(n) = solin(n)*frcsol(ns)*psf
          fsnt(n) = fsnt(n) + solflx(n)*(fluxdn(1,n)    - fluxup(1,n))
          fsns(n) = fsns(n) + solflx(n)*(fluxdn(kzp1,n) - fluxup(kzp1,n))
          solout(n) = solout(n) + solflx(n)*fluxup(0,n)
          fswup(0,n) = fswup(0,n) + solflx(n)*fluxup(0,n)
          fswdn(0,n) = fswdn(0,n) + solflx(n)*fluxdn(0,n)
          !
          ! Down spectral fluxes need to be in mks; thus the 0.001
          ! conversion factors
          !
          if ( wavmid < 0.7_rkx ) then
            sols(n) = sols(n) + (exptdn(kzp1,n)*solflx(n))*d_r1000
            solsd(n) = solsd(n) + &
                       ((fluxdn(kzp1,n)-exptdn(kzp1,n))*solflx(n))*d_r1000
            abv(n) = abv(n) + ((solflx(n) *               &
                       (fluxdn(kzp1,n)-fluxup(kzp1,n)))*  &
                       (d_one-asw(n))/(d_one-diralb(n)))*d_r1000
          else
            soll(n) = soll(n) + (exptdn(kzp1,n)*solflx(n))*d_r1000
            solld(n) = solld(n) + &
                   ((fluxdn(kzp1,n)-exptdn(kzp1,n))*solflx(n))*d_r1000
            fsnirtsq(n) = fsnirtsq(n) + solflx(n)*(fluxdn(0,n)-fluxup(0,n))
            abv(n) = abv(n) + &
                       ((solflx(n)*(fluxdn(kzp1,n)-fluxup(kzp1,n)))* &
                       (d_one-alw(n))/(d_one-diralb(n)))*d_r1000
          end if
          fsnirt(n) = fsnirt(n)+wgtint*solflx(n)*(fluxdn(0,n)-fluxup(0,n))
          do k = 0 , kz
            totfld(k,n) = totfld(k,n) + solflx(n)*flxdiv(k,n)
            fswup(k+1,n) = fswup(k+1,n) + solflx(n)*fluxup(k+1,n)
            fswdn(k+1,n) = fswdn(k+1,n) + solflx(n)*fluxdn(k+1,n)
          end do
          ! solar is incident visible solar radiation
          if ( ns == 8 ) then
            sol(n) = (solflx(n)*fluxdn(kzp1,n))*d_r1000
          end if
        end if
      end do

      !sfltot = d_zero
      !do n = n1 , n2
      !  if ( czengt0(n) ) then
      !    sfltot = sfltot + solflx(n)
      !  end if
      !end do

#ifdef STDPAR
      do concurrent ( n = n1:n2 ) local(k)
#else
      do n = n1 , n2
#endif
        tauaer(n) = tauxar3d(1,n,ns)
        tauasc(n) = tauasc3d(1,n,ns)
        ftota(n) =  ftota3d(1,n,ns)
        gtota(n) =  gtota3d(1,n,ns)
        do k = 2 , kz
          tauaer(n) = tauaer(n) + tauxar3d(k,n,ns)
          tauasc(n) = tauasc(n) + tauasc3d(k,n,ns)
          ftota(n) =  ftota(n)  + ftota3d(k,n,ns)
          gtota(n) =  gtota(n)  + gtota3d(k,n,ns)
        end do
      end do
      !FAB
      ! CLEAR SKY CALCULATION PLUS AEROSOL
      ! FORCING RAD CLR is called 2 times , one with O aerosol OP , and
      ! one with actual aerosol. DIFFERENCE  in net TOA SW for the two
      ! case is saved as one more variable in the rad file. The
      ! outputed TOASW ( fsntc, clrst) is accounting for aerosol.
      !FAB
      if ( linteract ) then
        !
        ! Following code is the diagnostic clear sky computation:
        !
        ! Compute delta-Eddington solution reflectivities and
        ! transmissivities for the entire column; note, for
        ! convenience, we use mod_the same reflectivity and transmissivity
        ! arrays as for the full calculation above, where 0 for layer
        ! quantities refers to the entire atmospheric column, and where
        ! 0 for interface quantities refers to top of atmos- phere,
        ! while 1 refers to the surface:
        !
        call radclr(n1,n2,trayoslp,czen,czengt0,.true.,pflx,  &
                    abh2o(ns),abco2(ns),abo2(ns),abo3(ns),    &
                    uth2o,uto3,utco2,uto2,tauaer,tauasc,gtota,&
                    ftota,tottrn,exptdn,rdndif,rdif,tdif,rdir,tdir,explay)
        !
        ! Compute reflectivity to direct and diffuse radiation for
        ! entire column; 0,1 on layer quantities refers to two
        ! effective layers overlying surface; 0 on interface quantities
        ! refers to top of column; 2 on interface quantities refers to
        ! the surface:
        !
#ifdef STDPAR
        do concurrent ( n = n1:n2 ) local(k,rdenom)
#else
        do n = n1 , n2
#endif
          if ( czengt0(n) ) then
            rupdir(2,n) = diralb(n)
            rupdif(2,n) = difalb(n)
            do k = 1 , 0 , -1
              rdenom = d_one/(d_one-rdif(k,n)*rupdif(k+1,n))
              rupdir(k,n) = rdir(k,n) + tdif(k,n) *                    &
                            (rupdir(k+1,n)*explay(k,n)+rupdif(k+1,n) * &
                            (tdir(k,n)-explay(k,n)))*rdenom
              rupdif(k,n) = rdif(k,n) + rupdif(k+1,n)*(tdif(k,n)**2)*rdenom
            end do
            !
            ! Compute up and down fluxes for each interface, using the added
            ! atmospheric layer properties at each interface:
            !
            do k = 0 , 2
              rdenom = d_one/(d_one-rdndif(k,n)*rupdif(k,n))
              fluxup(k,n) = (exptdn(k,n)*rupdir(k,n)+(tottrn(k,n) - &
                            exptdn(k,n))*rupdif(k,n))*rdenom
              fluxdn(k,n) = exptdn(k,n) +                           &
                            (tottrn(k,n)-exptdn(k,n)+exptdn(k,n) *  &
                             rupdir(k,n)*rdndif(k,n))*rdenom
            end do
            ! SAVE the ref net TOA flux
            ! ( and put back the cumul variables to 0.)
            toafsntc(n) = toafsntc(n) + solflx(n)*(fluxdn(0,n)-fluxup(0,n))
            toafsnsc(n) = toafsnsc(n) + solflx(n)*(fluxdn(2,n)-fluxup(2,n))
          end if
        end do

        !toafsnrtc = d_zero
        !do n = n1 , n2
        !  if ( czengt0(n) ) then
        !    toafsnrtc = toafsnrtc + wgtint*solflx(n)*(fluxdn(0,n)-fluxup(0,n))
        !  end if
        !end do
        !
        ! End of clear sky calculation with O aerosol OP
        !
      end if
      !
      ! Following code is the diagnostic clear sky computation:
      !
      ! Compute delta-Eddington solution reflectivities and
      ! transmissivities for the entire column; note, for convenience,
      ! we use mod_the same reflectivity and transmissivity arrays as for
      ! the full calculation above, where 0 for layer quantities refers
      ! to the entire atmospheric column, and where 0 for interface
      ! quantities refers to top of atmos- phere, while 1 refers to the
      ! surface:
      !
      call radclr(n1,n2,trayoslp,czen,czengt0,.false.,pflx, &
                  abh2o(ns),abco2(ns),abo2(ns),abo3(ns),    &
                  uth2o,uto3,utco2,uto2,tauaer,tauasc,gtota,&
                  ftota,tottrn,exptdn,rdndif,rdif,tdif,rdir,tdir,explay)
      !
      ! Compute reflectivity to direct and diffuse radiation for entire
      ! column; 0,1 on layer quantities refers to two effective layers
      ! overlying surface; 0 on interface quantities refers to top of
      ! column; 2 on interface quantities refers to the surface:
      !
#ifdef STDPAR
      do concurrent ( n = n1:n2 ) local(rdenom,k)
#else
      do n = n1 , n2
#endif
        if ( czengt0(n) ) then
          rupdir(2,n) = diralb(n)
          rupdif(2,n) = difalb(n)
          do k = 1 , 0 , -1
            rdenom = d_one/(d_one-rdif(k,n)*rupdif(k+1,n))
            rupdir(k,n) = rdir(k,n) + tdif(k,n) *     &
                          (rupdir(k+1,n)*explay(k,n)+ &
                           rupdif(k+1,n)*(tdir(k,n)-explay(k,n)))*rdenom
            rupdif(k,n) = rdif(k,n) + rupdif(k+1,n)*tdif(k,n)**2*rdenom
          end do
          !
          ! Compute up and down fluxes for each interface, using the added
          ! atmospheric layer properties at each interface:
          !
          do k = 0 , 2
            rdenom = d_one/(d_one-rdndif(k,n)*rupdif(k,n))
            fluxup(k,n) = (exptdn(k,n)*rupdir(k,n)+(tottrn(k,n) - &
                           exptdn(k,n))*rupdif(k,n))*rdenom
            fluxdn(k,n) = exptdn(k,n) +                           &
                          (tottrn(k,n)-exptdn(k,n)+exptdn(k,n) *  &
                           rupdir(k,n)*rdndif(k,n))*rdenom
          end do
          fsntc(n) = fsntc(n) + solflx(n)*(fluxdn(0,n)-fluxup(0,n))
          fsnsc(n) = fsnsc(n) + solflx(n)*(fluxdn(2,n)-fluxup(2,n))
          fsnrtc(n) = fsnrtc(n)+wgtint*solflx(n)*(fluxdn(0,n)-fluxup(0,n))
        end if
      end do
      !
      ! End of clear sky calculation
      !
    end do ! End of spectral interval loop

    do is = 1 , 4
      do k = 1 , kzp1
        do n = n1, n2
          outtaucl(k,is,n) = outtaucl(k,is,n) / ww(is)
          outtauci(k,is,n) = outtauci(k,is,n) / ww(is)
        end do
      end do
    end do

    ! FAB calculation of TOA aerosol radiative forcing
    ! convert from cgs to MKS
    if ( linteract ) then
      do concurrent ( n = n1:n2 )
        if ( czengt0(n) ) then
          aeradfo(n) = -(toafsntc(n)-fsntc(n)) * d_r1000
          aeradfos(n) = -(toafsnsc(n)-fsnsc(n)) * d_r1000
        end if
      end do
    end if
    !
    ! Compute solar heating rate (k/s)
    !
    do concurrent( n = n1:n2, k = 1:kz )
      if ( czengt0(n) ) then
        qrs(k,n) = -(gocp*totfld(k,n))/(pnm(k,n)-pnm(k+1,n))
      end if
    end do
    !
    ! Set the downwelling flux at the surface
    !
    do concurrent ( n = n1:n2 )
      fsds(n) = fswdn(kzp1,n)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  ! Input arguments
  !
  ! ts      - Ground (skin) temperature
  ! emiss   - Emissivity of surface
  !
  ! Input arguments which are only passed to other routines
  !
  ! tnm     - Level temperature
  ! qnm     - Level moisture field
  ! o3vmr   - ozone volume mixing ratio
  ! pbr     - Level pressure (dynes/cm*2)
  ! pnm     - Model interface pressure (dynes/cm*2)
  ! pmln    - Ln(pmid)
  ! piln    - Ln(pint)
  ! plco2   - Path length co2
  ! plh2o   - Path length h2o
  ! n2o     - nitrous oxide mass mixing ratio
  ! ch4     - methane mass mixing ratio
  ! cfc11   - cfc11 mass mixing ratio
  ! cfc12   - cfc12 mass mixing ratio
  ! tclrsf  - Clear sky fraction
  !
  ! Input/Output arguments
  !
  ! cld      - Cloud cover
  !
  ! Output arguments
  !
  ! qrl     - Longwave heating rate
  ! flns    - Surface cooling flux
  ! flnt    - Net outgoing flux
  ! lwout   - outgoing LW flux
  ! lwin    - incoming LW flux
  ! flnsc   - Clear sky surface cooing
  ! flntc   - Net clear sky outgoing flux
  ! flwds   - Down longwave flux at surface
  !
  !  cs is clearsky
  !  Aerosol longwave added
  !
  subroutine radclw(n1,n2,labsem,ts,emiss,tnm,qnm,o3vmr,pbr,pnm,    &
                    pmln,piln,n2o,ch4,cfc11,cfc12,cld,tclrsf,       &
                    flns,flnt,lwout,lwin,flnsc,flntc,flwds,fslwdcs, &
                    aerlwfo,aerlwfos,absgasnxt,absgastot,emsgastot,qrl)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    logical , intent(in) :: labsem
    real(rkx) , dimension(n1:n2) , intent(in) :: ts , emiss
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: qnm , tnm
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: pbr , pmln , o3vmr
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: piln , pnm
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: cfc11 , cfc12
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: ch4 , n2o
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: cld
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: tclrsf
    real(rkx) , dimension(n1:n2) , intent(out) :: flns , flnsc , flnt
    real(rkx) , dimension(n1:n2) , intent(out) :: flntc , flwds , fslwdcs
    real(rkx) , dimension(n1:n2) , intent(out) :: lwout , lwin
    real(rkx) , dimension(n1:n2) , intent(out) :: aerlwfo , aerlwfos
    real(rkx) , dimension(kz,4,n1:n2) , intent(out) :: absgasnxt
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: emsgastot
    real(rkx) , dimension(kzp1,kzp1,n1:n2) , intent(out) :: absgastot
    real(rkx) , dimension(kz,n1:n2) , intent(out) :: qrl
    logical :: lstart
    integer(ik4) :: n , khighest , irad , nradaer
    integer(ik4) :: k , km , k1 , k2 , k3
    real(rkx) :: bk1 , bk2 , absbt , tmp , tmp1 , delt , delt1
    integer(ik4) :: km1 , km2 , km3 , km4
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radclw'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k)
#else
    do n = n1 , n2
#endif
      rtclrsf(1,n) = d_one/tclrsf(1,n)
      do k = 1 , kz
        fclb4(k,n) = d_zero
        fclt4(k,n) = d_zero
        tclrsf(k+1,n) = tclrsf(k,n)*(d_one-cld(k+1,n))
        rtclrsf(k+1,n) = d_one/tclrsf(k+1,n)
      end do
    end do
    !
    ! Calculate some temperatures needed to derive absorptivity and
    ! emissivity, as well as some h2o path lengths
    !
    call radtpl(n1,n2,ts,tnm,pnm,qnm,pmln,piln,plh2o, &
                tint,tint4,tlayr,tlayr4,tplnka,s2t,s2c,wh2op,tplnke)
    !
    ! do emissivity and absorptivity calculations
    ! only if abs/ems computation
    !
    if ( labsem ) then
      !
      ! Compute ozone path lengths at frequency of a/e calculation.
      !
      call radoz2(n1,n2,o3vmr,pnm,plos,plol)
      !
      ! Compute trace gas path lengths
      !
      call trcpth(n1,n2,tnm,pnm,qnm,cfc11,cfc12,n2o,ch4,co2mmr,  &
                  ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212,  &
                  uco213,uco221,uco222,uco223,bn2o0,bn2o1,bch4,  &
                  uptype)
      !
      ! Calculate trace gas Planck functions
      !
      call trcplk(n1,n2,tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
      !
      ! Compute total emissivity:
      !
      call radems(n1,n2,pnm,tint,tint4,tlayr,tlayr4,tplnke,plos,plol,    &
                  plh2o,ucfc11,ucfc12,un2o0,un2o1,bn2o0,bn2o1,uch4,bch4, &
                  uco211,uco212,uco213,uco221,uco222,uco223,uptype,      &
                  wh2op,s2c,s2t,emplnk,co2t,co2em,co2eml,h2otr,emsgastot)
      !
      ! Compute total absorptivity:
      !
      call radabs(n1,n2,pnm,pbr,piln,pmln,tint,tlayr,co2em,co2eml,   &
                  tplnka,s2c,s2t,wh2op,h2otr,co2t,plco2,plh2o,plol,  &
                  plos,abplnk1,abplnk2,ucfc11,ucfc12,un2o0,un2o1,    &
                  bn2o0,bn2o1,uch4,bch4,uco211,uco212,uco213,uco221, &
                  uco222,uco223,uptype,absgasnxt,absgastot)
    end if
    !
    ! Find the lowest and highest level cloud for each grid point
    ! Note: Vertical indexing here proceeds from bottom to top
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k)
#else
    do n = n1 , n2
#endif
      klov(n) = 0
      done(n) = .false.
      do k = 1 , kz
        if ( .not. done(n) .and. cld(kzp2-k,n) > 0.0_rkx ) then
          done(n) = .true.
          klov(n) = k
        end if
      end do
      if ( klov(n) > 0 ) then
        skip(n) = .false.
      else
        skip(n) = .true.
      end if
      khiv(n) = klov(n)
      done(n) = .false.
      do k = kz , 1 , -1
        if ( skip(n) ) cycle
        if ( .not. done(n) .and. cld(kzp2-k,n) > 0.0_rkx ) then
          done(n) = .true.
          khiv(n) = k
        end if
      end do
      khivm(n) = khiv(n) - 1
      if ( skip(n) ) cycle
      !
      ! Note: Vertical indexing here proceeds from bottom to top
      !
      do k = klov(n) , khiv(n)
        fclt4(kzp1-k,n) = stebol*tint4(kzp2-k,n)
        fclb4(kzp1-k,n) = stebol*tint4(kzp3-k,n)
      end do
    end do
    !
    ! option to calculate LW aerosol radiative forcing
    !
    ! FAB LW radiative forcing ( rad=1 : avec dust)
    !
    nradaer = 1
    if ( linteract ) then
      nradaer = 2
    end if

    do irad = 1 , nradaer
      !
      ! Compute sums used in integrals (all longitude points)
      !
      ! Definition of bk1 & bk2 depends on finite differencing.  for
      ! trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
      ! layers only.
      !
      ! delt=t**4 in layer above current sigma level km.
      ! delt1=t**4 in layer below current sigma level km.
      !
#ifdef STDPAR
      do concurrent ( n = n1:n2 ) local(bk1,bk2,absbt,delt,delt1,k,km)
#else
      do n = n1 , n2
#endif
        do km = 1 , 4
          do k = 1 , kz
            absnxt(k,km,n) = absgasnxt(k,km,n)
          end do
        end do
        do k = 1 , kzp1
          emstot(k,n)   = emsgastot(k,n)
        end do
        do km = 1 , kzp1
          do k = 1 , kzp1
            abstot(k,km,n) = absgastot(k,km,n)
          end do
        end do
        if  ( linteract .and. irad == 2 ) then
          do km = 1 , 4
            do k = 1 , kz
              absnxt(k,km,n) = d_one-(d_one-absgasnxt(k,km,n)) * &
                              (aertrlw(k,k+1,n)**xuinpl(k,km,n))
            end do
          end do
          do k = 1 , kzp1
            emstot(k,n) = d_one-(d_one-emsgastot(k,n)) * aertrlw(k,1,n)
          end do
          do km = 1 , kzp1
            do k = 1 , kzp1
              abstot(k,km,n) = d_one-(d_one-absgastot(k,km,n)) * &
                               aertrlw(k,km,n)
            end do
          end do
        end if
        delt = tint4(kz,n) - tlayr4(kzp1,n)
        delt1 = tlayr4(kzp1,n) - tint4(kzp1,n)
        fis(kzp1,kzp1,n) = stebol*(delt1*absnxt(kz,1,n) + &
                         delt*absnxt(kz,4,n))
        fis(kz,kzp1,n) = stebol*(delt*absnxt(kz,2,n) + &
                         delt1*absnxt(kz,3,n))
        do k = 1 , kz - 1
          bk2 = (abstot(k,kz,n)+abstot(k,kzp1,n))*d_half
          bk1 = bk2
          fis(k,kzp1,n) = stebol*(bk2*delt+bk1*delt1)
        end do
        do km = kz , 2 , -1
          delt = tint4(km-1,n) - tlayr4(km,n)
          delt1 = tlayr4(km,n) - tint4(km,n)
          !
          ! All k, km>1
          !
          do k = kzp1 , 1 , -1
            if ( k == km ) then
              bk2 = absnxt(km-1,4,n)
              bk1 = absnxt(km-1,1,n)
            else if ( k == km-1 ) then
              bk2 = absnxt(km-1,2,n)
              bk1 = absnxt(km-1,3,n)
            else
              bk2 = d_half * (abstot(k,km-1,n) + abstot(k,km,n))
              bk1 = bk2
            end if
            fis(k,km,n) = fis(k,km+1,n) + stebol*(bk2*delt + bk1*delt1)
          end do
        end do
        !
        ! Computation of clear sky fluxes always set first level of fsul
        !
        fsul(kzp1,n) = emiss(n) * stebol * ts(n)**4
        !
        ! Downward clear sky fluxes store intermediate quantities in down
        ! flux Initialize fluxes to clear sky values.
        !
        tmp = fsul(kzp1,n) - stebol*tint4(kzp1,n)
        fsul(1,n) = fsul(kzp1,n) - abstot(1,kzp1,n) * tmp + fis(1,2,n)
        fsdl(1,n) = emstot(1,n) * stebol * tplnke(n)**4
        ful(1,n) = fsul(1,n)
        fdl(1,n) = fsdl(1,n)
        do k = 2 , kz
          fsul(k,n) = fsul(kzp1,n) - abstot(k,kzp1,n)*tmp + fis(k,k+1,n)
          ful(k,n) = fsul(k,n)
          fsdl(k,n) = stebol*(tplnke(n)**4) * emstot(k,n) - &
                              (fis(k,2,n)-fis(k,k+1,n))
          fdl(k,n) = fsdl(k,n)
        end do
        !
        ! fsdl(kzp1,n) assumes isothermal layer
        !
        !
        ! Store the downward emission from level 1 = total gas emission *
        ! sigma t**4.  fsdl does not yet include all terms
        !
        ful(kzp1,n) = fsul(kzp1,n)
        absbt = emstot(kzp1,n) * stebol * tplnke(n)**4
        fsdl(kzp1,n) = absbt - fis(kzp1,2,n)
        fdl(kzp1,n) = fsdl(kzp1,n)
      end do
      !
      ! FAB radiative forcing sur fsul
      !
      if ( linteract .and. irad == 1 ) then
#ifdef STDPAR
        do concurrent ( n = n1:n2 ) local(k1,k2)
#else
        do n = n1 , n2
#endif
          do k1 = 1 , kzp1
            fsul0(k1,n) = fsul(k1,n) ! save fsul0 = no dust
            fsdl0(k1,n) = fsdl(k1,n) !
            ful0(k1,n) = ful(k1,n)
            fdl0(k1,n) = fdl(k1,n)
            do k2 = 1 , kzp1
              fis0(k2,k1,n) = fis(k2,k1,n)
            end do
          end do
        end do
      end if

    end do ! end rad loop
    !
    ! FAB after this DO loop fsul account for dust LW effect
    ! which is OK in case of idirect=2
    ! also convert rad for to MKSA
    if ( linteract ) then
      do concurrent ( n = n1:n2 )
        aerlwfo(n) = (fsul0(1,n) - fsul(1,n) ) * d_r1000
        aerlwfos(n) = ( (fsul0(kzp1,n) - fsdl0(kzp1,n)) - &
                        (fsul(kzp1,n)  - fsdl(kzp1,n) ) ) * d_r1000
      end do
      ! return to no aerosol LW effect situation if idirect == 1
      if ( lzero ) then
#ifdef STDPAR
        do concurrent ( n = n1:n2 ) local(k1,k2)
#else
        do n = n1 , n2
#endif
          do k1 = 1 , kzp1
            fsul(k1,n) = fsul0(k1,n)
            fsdl(k1,n) = fsdl0(k1,n)
            ful(k1,n) = ful0(k1,n)
            fdl(k1,n) = fdl0(k1,n)
            do k2 = 1 , kzp1
              fis(k2,k1,n) = fis0(k2,k1,n)
            end do
          end do
        end do
      end if
    end if ! end aersol rad diagnostic
    !
    ! FAB MODIF  : save downward clear sky long wave flux in surface
    !
    do concurrent ( n = n1:n2 )
      fslwdcs(n) = fsdl(kzp1,n)
    end do
    !
    ! Modifications for clouds
    !
    ! Further qualify longitude subset for computations.  Select only
    ! those locations where there are clouds
    ! (total cloud fraction <= 1.e-3 treated as clear)
    !
    do concurrent ( n = n1:n2 )
      if ( tclrsf(kzp1,n) < 0.999_rkx ) then
        skip(n) = .false.
      else
        skip(n) = .true.
      end if
    end do
    !
    ! Compute downflux at level 1 for cloudy sky
    !
    do n = n1 , n2
      if ( skip(n) ) cycle
      !
      ! First clear sky flux plus flux from cloud at level 1
      !
      fdl(kzp1,n) = fsdl(kzp1,n)*tclrsf(kz,n) * &
                    rtclrsf(kzp1-khiv(n),n)+fclb4(kz-1,n)*cld(kz,n)
    end do
    !
    ! Flux emitted by other layers
    ! Note: Vertical indexing here proceeds from bottom to top
    !
    khighest = khiv(intmax(khiv))
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) &
      local(tmp1,lstart,km,km1,km2,km3,km4,k,k1,k2,k3)
#else
    do n = n1 , n2
#endif
      if ( skip(n) ) cycle
      do km = 3 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km4 = kzp4 - km
        if ( km <= khiv(n) ) then
          tmp1 = cld(km2,n)*tclrsf(kz,n)*rtclrsf(km2,n)
          fdl(kzp1,n) = fdl(kzp1,n) + (fclb4(km1,n)-fis(kzp1,km4,n))*tmp1
        end if
      end do
      !
      ! Note: Vertical indexing here proceeds from bottom to top
      !
      do k = 1 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        if ( k >= klov(n) .and. k <= khivm(n) ) then
          ful(k2,n) = fsul(k2,n)*(tclrsf(kzp1,n)*rtclrsf(k1,n))
        end if
        do km = 1 , k
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          if ( k <= khivm(n) .and. km >= klov(n) .and. km <= khivm(n)) then
            ful(k2,n) = ful(k2,n) + &
              (fclt4(km1,n)+fis(k2,k3,n)-fis(k2,km3,n)) * &
               cld(km2,n)*(tclrsf(km1,n)*rtclrsf(k1,n))
          end if
        end do ! km = 1 , k
      end do   ! k = 1 , khighest-1
      lstart = .false.
      do k = 1 , kzp1
        k2 = kzp2 - k
        k3 = kzp3 - k
        if ( k >= khiv(n) ) then
          lstart = .true.
          ful(k2,n) = fsul(k2,n)*tclrsf(kzp1,n)*rtclrsf(kzp1-khiv(n),n)
        end if
        do km = 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          if ( lstart .and. km >= klov(n) .and. km <= khiv(n) ) then
            ful(k2,n) = ful(k2,n) + (cld(km2,n)*tclrsf(km1,n)* &
              rtclrsf(kzp1-khiv(n),n))* &
              (fclt4(km1,n)+fis(k2,k3,n)-fis(k2,km3,n))
          end if
        end do  ! km = 1 , khighest
      end do    ! k = 1 , kzp1
      !
      ! Computation of the downward fluxes
      !
      do k = 2 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        if ( k <= khivm(n) ) fdl(k2,n) = d_zero
        do km = k + 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km4 = kzp4 - km
          if ( k <= khiv(n) .and. &
               km >= max0(k+1,klov(n)) .and. km <= khiv(n) ) then
            fdl(k2,n) = fdl(k2,n)+(cld(km2,n)*tclrsf(k1,n)*rtclrsf(km2,n)) * &
                    (fclb4(km1,n)-fis(k2,km4,n)+fis(k2,k3,n))
          end if
        end do ! km = k+1 , khighest
        if ( k <= khivm(n) ) then
           fdl(k2,n) = fdl(k2,n) + &
                fsdl(k2,n)*(tclrsf(k1,n)*rtclrsf(kzp1-khiv(n),n))
        end if
      end do  ! k = 1 , khighest-1
    end do
    !
    ! End cloud modification loops
    !
    !
    ! Downward longwave flux
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k)
#else
    do n = n1 , n2
#endif
      flwds(n) = fdl(kzp1,n)
      !
      ! Net flux
      !
      flns(n) = ful(kzp1,n) - fdl(kzp1,n)
      !
      ! Clear sky flux at top of atmosphere
      !
      flntc(n) = fsul(1,n)
      flnsc(n) = fsul(kzp1,n) - fsdl(kzp1,n)
      !
      ! Outgoing ir
      !
      flnt(n) = ful(1,n) - fdl(1,n)
      lwout(n) = ful(1,n)
      lwin(n) = fdl(1,n)
      !
      ! Computation of longwave heating (k per sec)
      !
      do k = 1 , kz
        qrl(k,n) = (ful(k,n)-fdl(k,n)-ful(k+1,n)+fdl(k+1,n))*gocp / &
                  ((pnm(k,n)-pnm(k+1,n)))
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  ! Input arguments
  !
  !   trayoslp  - Tray/sslp
  !   czen      - Cosine zenith angle
  !   czengt0   - logical cosine zenith angle greater than zero
  !   lcls      - Add or remove aerosol effect
  !   pflx      - Interface pressure
  !   abh2o     - Absorption coefficiant for h2o
  !   abo3      - Absorption coefficiant for o3
  !   abco2     - Absorption coefficiant for co2
  !   abo2      - Absorption coefficiant for o2
  !   uth2o     - Total column absorber amount of h2o
  !   uto3      - Total column absorber amount of  o3
  !   utco2     - Total column absorber amount of co2
  !   uto2      - Total column absorber amount of  o2
  !   tauaer    - Total column aerosol extinction
  !   tauasc    - Aerosol single scattering albedo * extinction
  !   gtota     - Aerosol asymmetry parameter * SSA * EXT
  !   ftota     - Aerosol *forward scattering fraction * SSA * EXT
  !
  ! Output
  !
  !   tottrn    - Total transmission for layers above
  !   exptdn    - Solar beam exp down transmn from top
  !   rdndif    - Added dif ref for layers above
  !   rdif      - Layer refflectivity to diffuse rad
  !   tdif      - Layer transmission to diffuse rad
  !   rdir      - Layer reflectivity to direct rad
  !   tdir      - Layer transmission to direct rad
  !   explay    - Solar beam exp transmn for layer
  !
  subroutine radclr(n1,n2,trayoslp,czen,czengt0,lcls,pflx,  &
                    abh2o,abco2,abo2,abo3,uth2o,uto3,utco2, &
                    uto2,tauaer,tauasc,gtota,ftota,tottrn,  &
                    exptdn,rdndif,rdif,tdif,rdir,tdir,explay)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    logical , intent(in) :: lcls
    real(rkx) , intent(in) :: trayoslp
    real(rkx) , dimension(n1:n2) , intent(in) :: czen
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(in) :: pflx
    real(rkx) , intent(in) :: abh2o , abco2 , abo2 , abo3
    real(rkx) , dimension(n1:n2) , intent(in) :: uth2o
    real(rkx) , dimension(n1:n2) , intent(in) :: uto3
    real(rkx) , dimension(n1:n2) , intent(in) :: utco2
    real(rkx) , dimension(n1:n2) , intent(in) :: uto2
    real(rkx) , dimension(n1:n2) , intent(in) :: tauaer
    real(rkx) , dimension(n1:n2) , intent(in) :: tauasc ! waer * tauaer
    real(rkx) , dimension(n1:n2) , intent(in) :: gtota  ! gaer * waer * tauaer
    real(rkx) , dimension(n1:n2) , intent(in) :: ftota  ! faer * waer * tauaer
    logical , dimension(n1:n2) , intent(in) :: czengt0
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(out) :: tottrn
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(out) :: exptdn
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(out) :: rdndif
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: explay
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: rdir , rdif
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: tdir , tdif
    !
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
    integer(ik4) :: n
    real(rkx) :: arg , rdenom , rdirexp , tdnmexp
    real(rkx) :: tautot , wtot , gtot , ftot , extins
    real(rkx) :: ts , ws , gs , lm , alp , gam , ue , ne
    real(rkx) :: apg , amg , taugab , tauray
    integer(ik4) :: k

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radclr'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Compute total direct beam transmission, total transmission, and
    ! reflectivity for diffuse radiation (from below) for all layers
    ! above each interface by starting from the top and adding layers down:
    !
    ! The top layer is assumed to be a purely absorbing ozone layer, and
    ! that the mean diffusivity for diffuse mod_transmission is 1.66:
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) &
      local(arg,rdenom,rdirexp,tdnmexp,tautot,wtot,gtot,ftot,extins, &
            ts,ws,gs,lm,alp,gam,ue,ne,apg,amg,taugab,tauray,k)
#else
    do n = n1 , n2
#endif
      !-------------------------------------------------------------------
      !
      ! Initialize all total transmimission values to 0, so that nighttime
      ! values from previous computations are not used:
      !
      do k = 1 , kzp1
        tottrn(k,n) = d_zero
      end do
      if ( czengt0(n) ) then
        taugab = abo3*uto3(n)
        ! Limit argument of exponential, in case czen is very small:
        arg = min(taugab/czen(n),mxarg)
        explay(0,n) = exp(-arg)
        tdir(0,n) = explay(0,n)
        !
        ! Same limit for diffuse mod_transmission:
        !
        arg = min(1.66_rkx*taugab,mxarg)
        tdif(0,n) = exp(-arg)
        rdir(0,n) = d_zero
        rdif(0,n) = d_zero
        !
        ! Initialize top interface of extra layer:
        !
        exptdn(0,n) = d_one
        rdndif(0,n) = d_zero
        tottrn(0,n) = d_one
        rdndif(1,n) = rdif(0,n)
        tottrn(1,n) = tdir(0,n)
        !
        ! Now, complete the rest of the column; if the total transmission
        ! through the top ozone layer is less than trmin, then no
        ! delta-Eddington computation for the underlying column is done:
        !
        do k = 1 , 1
          !
          ! Initialize current layer properties to zero;only if total
          ! transmission to the top interface of the current layer exceeds
          ! the minimum, will these values be computed below:
          !
          rdir(k,n) = d_zero
          rdif(k,n) = d_zero
          tdir(k,n) = d_zero
          tdif(k,n) = d_zero
          explay(k,n) = d_zero
          !
          ! Calculates the solar beam transmission, total transmission,
          ! and reflectivity for diffuse radiation from below at the
          ! top of the current layer:
          !
          exptdn(k,n) = exptdn(k-1,n)*explay(k-1,n)
          rdenom = d_one/(d_one-rdif(k-1,n)*rdndif(k-1,n))
          rdirexp = rdir(k-1,n)*exptdn(k-1,n)
          tdnmexp = tottrn(k-1,n) - exptdn(k-1,n)
          tottrn(k,n) = exptdn(k-1,n)*tdir(k-1,n) + &
                        tdif(k-1,n)*(tdnmexp+rdndif(k-1,n)*rdirexp)*rdenom
          rdndif(k,n) = rdif(k-1,n) + &
                        (rdndif(k-1,n)*tdif(k-1,n))*(tdif(k-1,n)*rdenom)
          !
          ! Compute next layer delta-Eddington solution only if total
          ! transmission of radiation to the interface just above the layer
          ! exceeds trmin.
          !
          if ( tottrn(k,n) > trmin ) then
            !
            ! Remember, no ozone absorption in this layer:
            !
            tauray = trayoslp*pflx(kzp1,n)
            taugab = abh2o*uth2o(n) + abco2*utco2(n) + abo2*uto2(n)
            if ( lcls ) then
              tautot = tauray + taugab
              wtot = (wray*tauray)/tautot
              gtot = (gray*wray*tauray)/(wtot*tautot)
              ftot = (fray*wray*tauray/(wtot*tautot))
            else
              tautot = tauray + taugab + tauaer(n)
              wtot = (wray*tauray + tauasc(n))/tautot
              gtot = (gray*wray*tauray + gtota(n))/(wtot*tautot)
              ftot = (fray*wray*tauray + ftota(n))/(wtot*tautot)
            end if
            ts = taus(wtot,ftot,tautot)
            ws = omgs(wtot,ftot)
            gs = asys(gtot,ftot)
            lm = el(ws,gs)
            alp = xalpha(ws,czen(n),gs,lm)
            gam = xgamma(ws,czen(n),gs,lm)
            ue = f_u(ws,gs,lm)
            !
            ! Limit argument of exponential, in case lm very large:
            !
            arg = min(lm*ts,mxarg)
            extins = exp(-arg)
            ne = f_n(ue,extins)
            rdif(k,n) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
            tdif(k,n) = d_four*ue/ne
            ! Limit argument of exponential, in case czen is very small:
            arg = min(ts/czen(n),mxarg)
            explay(k,n) = exp(-arg)
            apg = alp + gam
            amg = alp - gam
            rdir(k,n) = amg*(tdif(k,n)*explay(k,n)-d_one)+apg*rdif(k,n)
            tdir(k,n) = apg*tdif(k,n) + &
                        (amg*rdif(k,n)-(apg-d_one))*explay(k,n)
            !
            ! Under rare conditions, reflectivies and transmissivities
            ! can be negative; zero out any negative values
            !
            rdir(k,n) = max(rdir(k,n),d_zero)
            tdir(k,n) = max(tdir(k,n),d_zero)
            rdif(k,n) = max(rdif(k,n),d_zero)
            tdif(k,n) = max(tdif(k,n),d_zero)
          end if
        end do
        k = 2
        exptdn(k,n) = exptdn(k-1,n)*explay(k-1,n)
        rdenom = d_one/(d_one-rdif(k-1,n)*rdndif(k-1,n))
        rdirexp = rdir(k-1,n)*exptdn(k-1,n)
        tdnmexp = tottrn(k-1,n) - exptdn(k-1,n)
        tottrn(k,n) = exptdn(k-1,n)*tdir(k-1,n) + &
                      tdif(k-1,n)*(tdnmexp+rdndif(k-1,n)*rdirexp)*rdenom
        rdndif(k,n) = rdif(k-1,n) + &
                      (rdndif(k-1,n)*tdif(k-1,n))*(tdif(k-1,n)*rdenom)
      end if
    end do
    !
    ! Compute total direct beam transmission, total transmission, and
    ! reflectivity for diffuse radiation (from below) for both layers
    ! above the surface:
    !
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  ! Input Arguments
  !   trayoslp - Tray/sslp
  !   czen     - Cosine zenith angle
  !   czengt0  - Logical czen greater zero
  !   pflx     - Interface pressure
  !   abh2o    - Absorption coefficiant for h2o
  !   abo3     - Absorption coefficiant for o3
  !   abco2    - Absorption coefficiant for co2
  !   abo2     - Absorption coefficiant for o2
  !   uh2o     - Layer absorber amount of h2o
  !   uo3      - Layer absorber amount of  o3
  !   uco2     - Layer absorber amount of co2
  !   uo2      - Layer absorber amount of  o2
  !   tauxcl   - Cloud extinction optical depth (liquid)
  !   tauxci   - Cloud extinction optical depth (ice)
  !   tauaer   - Aerosol extinction optical depth
  !   wcl      - Cloud single scattering albedo (liquid)
  !   gcl      - Cloud asymmetry parameter (liquid)
  !   fcl      - Cloud forward scattered fraction (liquid)
  !   wci      - Cloud single scattering albedo (ice)
  !   gci      - Cloud asymmetry parameter (ice)
  !   fci      - Cloud forward scattered fraction (ice)
  !   tauasc   - Aerosol single scattering albedo * extinction
  !   gtota    - Aerosol asymmetry parameter * SSA * extinction
  !   ftota    - Aerosol forward scattered fraction * SSA * extinction
  !
  ! Output Arguments
  !   tottrn   - Total transmission for layers above
  !   exptdn   - Solar beam exp down transm from top
  !   rdndif   - Added dif ref for layers above
  !   rdif     - Layer refflectivity to diffuse rad
  !   tdif     - Layer transmission to diffuse rad
  !   rdir     - Layer reflectivity to direct rad
  !   tdir     - Layer transmission to direct rad
  !   explay   - Solar beam exp transm for layer
  !
  !-----------------------------------------------------------------------
  !
  subroutine radded(n1,n2,trayoslp,czen,czengt0,pflx, &
      abh2o,abo3,abco2,abo2,uh2o,uo3,uco2,uo2,        &
      tauxcl,wcl,gcl,fcl,tauxci,wci,gci,fci,          &
      tauaer,tauasc,gtota,ftota,                      &
      tottrn,exptdn,rdndif,rdif,tdif,rdir,tdir,explay)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(n1:n2) , intent(in) :: czen
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(in) :: pflx
    real(rkx) , intent(in) :: abh2o , abo3 , abco2 , abo2
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: uh2o , uo3
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: uco2 , uo2
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: tauxcl , tauxci , tauaer
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: wcl , gcl , fcl
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: wci , gci , fci
    real(rkx) , dimension(0:kz,n1:n2) , intent(in) :: tauasc , gtota , ftota
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(out) :: tottrn
    real(rkx) , dimension(0:kzp1,n1:n2) , intent(out) :: exptdn , rdndif
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: rdif , tdif
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: rdir , tdir
    real(rkx) , dimension(0:kz,n1:n2) , intent(out) :: explay
    logical , dimension(n1:n2) , intent(in) :: czengt0
    real(rkx) , intent(in) :: trayoslp
    !
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
    ! Intermediate terms for delta-eddington solution
    !
    ! alp      - Temporary for xalpha
    ! gam      - Temporary for xgamma
    ! lm       - Temporary for el
    ! ne       - Temporary for f_n
    ! ue       - Temporary for f_u
    ! arg      - Exponential argument
    ! extins   - Extinction
    ! amg      - Alp - gam
    ! apg      - Alp + gam
    !
    integer(ik4) :: n
    integer(ik4) :: k
    real(rkx) :: tautot , taucsc , wtau , wt , wtot , gtot , ftot
    real(rkx) :: ws , gs , ts , lm , alp , gam , ne , ue
    real(rkx) :: apg , amg , extins , rdenom , rdirexp , tdnmexp
    real(rkx) :: taugab , tauray
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radded'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Compute total direct beam transmission, total transmission, and
    ! reflectivity for diffuse radiation (from below) for all layers
    ! above each interface by starting from the top and adding layers down:
    ! For the extra layer above model top:
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) &
      local(tautot,taucsc,wtau,wt,wtot,gtot,ftot,ws,gs,ts,lm,alp,gam, &
            ne,ue,apg,amg,extins,rdenom,rdirexp,tdnmexp,tauray,taugab,k)
#else
    do n = n1 , n2
#endif
      !-----------------------------------------------------------------
      !
      ! Initialize all total transmission values to 0, so that nighttime
      ! values from previous computations are not used:
      !
      do k = 1 , kzp1
        tottrn(k,n) = d_zero
      end do
      if ( czengt0(n) ) then
        tauray = trayoslp * (pflx(1,n)-pflx(0,n))
        taugab = abh2o*uh2o(0,n) + abo3*uo3(0,n) + &
                 abco2*uco2(0,n) + abo2*uo2(0,n)
        if ( lzero ) then
          tautot = tauray + taugab + tauxcl(0,n) + tauxci(0,n)
          taucsc = tauxcl(0,n)*wcl(0,n) + &
                   tauxci(0,n)*wci(0,n)
          wtau = wray*tauray
          wt = wtau + taucsc
          wtot = wt/tautot
          gtot = (wtau*gray + gcl(0,n)*tauxcl(0,n)*wcl(0,n)+ &
                              gci(0,n)*tauxci(0,n)*wci(0,n)) / wt
          ftot = (wtau*fray + fcl(0,n)*tauxcl(0,n)*wcl(0,n) + &
                              fci(0,n)*tauxci(0,n)*wci(0,n)) / wt
        else
          tautot = tauray + taugab + &
                   tauxcl(0,n) + tauxci(0,n) + tauaer(0,n)
          taucsc = tauxcl(0,n)*wcl(0,n) + &
                   tauxci(0,n)*wci(0,n) + &
                   tauasc(0,n)
          wtau = wray*tauray
          wt = wtau + taucsc
          wtot = wt/tautot
          gtot = (wtau*gray + gcl(0,n)*tauxcl(0,n)*wcl(0,n) + &
                              gci(0,n)*tauxci(0,n)*wci(0,n) + &
                              gtota(0,n)) / wt
          ftot = (wtau*fray + fcl(0,n)*tauxcl(0,n)*wcl(0,n) + &
                              fci(0,n)*tauxci(0,n)*wci(0,n) + &
                              ftota(0,n)) / wt
        end if
        ts = taus(wtot,ftot,tautot)
        ws = omgs(wtot,ftot)
        gs = asys(gtot,ftot)
        lm = el(ws,gs)
        alp = xalpha(ws,czen(n),gs,lm)
        gam = xgamma(ws,czen(n),gs,lm)
        ue = f_u(ws,gs,lm)
        !
        ! Limit argument of exponential, in case lm*ts very large:
        !
        extins = exp(-min(lm*ts,mxarg))
        ne = f_n(ue,extins)
        rdif(0,n) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
        tdif(0,n) = d_four*ue/ne
        ! Limit argument of exponential, in case czen is very small:
        explay(0,n) = exp(-min(ts/czen(n),mxarg))
        apg = alp + gam
        amg = alp - gam
        rdir(0,n) = amg*(tdif(0,n)*explay(0,n)-d_one) + apg*rdif(0,n)
        tdir(0,n) = apg*tdif(0,n) + (amg*rdif(0,n)-(apg-d_one))*explay(0,n)
        !
        ! Under rare conditions, reflectivies and transmissivities can
        ! be negative; zero out any negative values
        !
        rdir(0,n) = max(rdir(0,n),d_zero)
        tdir(0,n) = max(tdir(0,n),d_zero)
        rdif(0,n) = max(rdif(0,n),d_zero)
        tdif(0,n) = max(tdif(0,n),d_zero)
        !
        ! Initialize top interface of extra layer:
        !
        exptdn(0,n) = d_one
        rdndif(0,n) = d_zero
        tottrn(0,n) = d_one
        rdndif(1,n) = rdif(0,n)
        tottrn(1,n) = tdir(0,n)
        !
        ! Now, continue down one layer at a time; if the total transmission
        ! to the interface just above a given layer is less than trmin,
        ! then no delta-eddington computation for that layer is done:
        !
        do k = 1 , kz
          !
          ! Initialize current layer properties to zero; only if total
          ! transmission to the top interface of the current layer exceeds
          ! the minimum, will these values be computed below:
          !
          rdir(k,n) = d_zero
          rdif(k,n) = d_zero
          tdir(k,n) = d_zero
          tdif(k,n) = d_zero
          explay(k,n) = d_zero
          !
          ! Calculates the solar beam transmission, total transmission,
          ! and reflectivity for diffuse radiation from below at the
          ! top of the current layer:
          !
          exptdn(k,n) = exptdn(k-1,n)*explay(k-1,n)
          rdenom = d_one/(d_one - rdif(k-1,n)*rdndif(k-1,n))
          rdirexp = rdir(k-1,n)*exptdn(k-1,n)
          tdnmexp = tottrn(k-1,n) - exptdn(k-1,n)
          tottrn(k,n) = exptdn(k-1,n)*tdir(k-1,n) + tdif(k-1,n) *     &
                      (tdnmexp+rdndif(k-1,n)*rdirexp)*rdenom
          rdndif(k,n) = rdif(k-1,n) + (rdndif(k-1,n)*tdif(k-1,n)) *   &
                      (tdif(k-1,n)*rdenom)
          !
          ! Compute next layer delta-eddington solution only if total
          ! transmission of radiation to the interface just above the layer
          ! exceeds trmin.
          !
          if ( tottrn(k,n) > trmin ) then
            tauray = trayoslp*(pflx(k+1,n)-pflx(k,n))
            taugab = abh2o*uh2o(k,n) + abo3*uo3(k,n) +  &
                     abco2*uco2(k,n) + abo2*uo2(k,n)

            if ( lzero ) then
              tautot = tauray + taugab + tauxcl(k,n) + tauxci(k,n)
              taucsc = tauxcl(k,n)*wcl(k,n) + tauxci(k,n)*wci(k,n)
              wtau = wray*tauray
              wt = wtau + taucsc
              wtot = wt/tautot
              gtot = (wtau*gray + gcl(k,n)*wcl(k,n)*tauxcl(k,n) + &
                                  gci(k,n)*wci(k,n)*tauxci(k,n))/wt
              ftot = (wtau*fray + fcl(k,n)*wcl(k,n)*tauxcl(k,n) + &
                                  fci(k,n)*wci(k,n)*tauxci(k,n))/wt
            else
              tautot = tauray+taugab + tauxcl(k,n)+tauxci(k,n)+tauaer(k,n)
              taucsc = tauxcl(k,n)*wcl(k,n) + &
                       tauxci(k,n)*wci(k,n) + &
                       tauasc(k,n)
              wtau = wray*tauray
              wt = wtau + taucsc
              wtot = wt/tautot
              gtot = (wtau*gray + gcl(k,n)*wcl(k,n)*tauxcl(k,n) + &
                                  gci(k,n)*wci(k,n)*tauxci(k,n) + &
                                  gtota(k,n))/wt
              ftot = (wtau*fray + fcl(k,n)*wcl(k,n)*tauxcl(k,n) + &
                                  fci(k,n)*wci(k,n)*tauxci(k,n) + &
                                  ftota(k,n))/wt
            end if
            ts = taus(wtot,ftot,tautot)
            ws = omgs(wtot,ftot)
            gs = asys(gtot,ftot)
            lm = el(ws,gs)
            alp = xalpha(ws,czen(n),gs,lm)
            gam = xgamma(ws,czen(n),gs,lm)
            ue = f_u(ws,gs,lm)
            !
            ! Limit argument of exponential, in case lm very large:
            !
            extins = exp(-min(lm*ts,mxarg))
            ne = f_n(ue,extins)
            rdif(k,n) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
            tdif(k,n) = d_four*ue/ne
            ! Limit argument of exponential, in case czen is very small:
            explay(k,n) = exp(-min(ts/czen(n),mxarg))
            apg = alp + gam
            amg = alp - gam
            rdir(k,n) = amg*(tdif(k,n)*explay(k,n)-d_one)+apg*rdif(k,n)
            tdir(k,n) = apg*tdif(k,n)+(amg*rdif(k,n)-(apg-d_one))*explay(k,n)
            !
            ! Under rare conditions, reflectivies and transmissivities
            ! can be negative; zero out any negative values
            !
            rdir(k,n) = max(rdir(k,n),d_zero)
            tdir(k,n) = max(tdir(k,n),d_zero)
            rdif(k,n) = max(rdif(k,n),d_zero)
            tdif(k,n) = max(tdif(k,n),d_zero)
          end if
        end do
        !
        ! Compute total direct beam transmission, total transmission, and
        ! reflectivity for diffuse radiation (from below) for all layers
        ! above the surface:
        !
        exptdn(kzp1,n) = exptdn(kz,n)*explay(kz,n)
        rdenom = d_one/(d_one-rdif(kz,n)*rdndif(kz,n))
        rdirexp = rdir(kz,n)*exptdn(kz,n)
        tdnmexp = tottrn(kz,n) - exptdn(kz,n)
        tottrn(kzp1,n) = exptdn(kz,n)*tdir(kz,n) + tdif(kz,n) *   &
                    (tdnmexp+rdndif(kz,n)*rdirexp)*rdenom
        rdndif(kzp1,n) = rdif(kz,n) + (rdndif(kz,n)*tdif(kz,n)) * &
                    (tdif(kz,n)*rdenom)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  ! n2o  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
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
  ! Input
  !   pbr     - Prssr at mid-levels (dynes/cm2)
  !   pnm     - Prssr at interfaces (dynes/cm2)
  !   co2em   - Co2 emissivity function
  !   co2eml  - Co2 emissivity function
  !   tplnka  - Planck fnctn level temperature
  !   s2c     - H2o continuum path length
  !   s2t     - H2o tmp and prs wghted path
  !   w       - H2o prs wghted path
  !   h2otr   - H2o trnsmssn fnct for o3 overlap
  !   plco2   - Co2 prs wghted path length
  !   plh2o   - H2o prs wfhted path length
  !   co2t    - Tmp and prs wghted path length
  !   tint    - Interface temperatures
  !   tlayr   - K-1 level temperatures
  !   plol    - Ozone prs wghted path length
  !   plos    - Ozone path length
  !   pmln    - Ln(pmidm1)
  !   piln    - Ln(pintm1)
  !   ucfc11  - CFC11 path length
  !   ucfc12  - CFC12 path length
  !   un2o0   - N2O path length
  !   un2o1   - N2O path length (hot band)
  !   uch4    - CH4 path length
  !   uco211  - CO2 9.4 micron band path length
  !   uco212  - CO2 9.4 micron band path length
  !   uco213  - CO2 9.4 micron band path length
  !   uco221  - CO2 10.4 micron band path length
  !   uco222  - CO2 10.4 micron band path length
  !   uco223  - CO2 10.4 micron band path length
  !   uptype  - continuum path length
  !   bn2o0   - pressure factor for n2o
  !   bn2o1   - pressure factor for n2o
  !   bch4    - pressure factor for ch4
  !   abplnk1 - non-nearest layer Planck factor
  !   abplnk2 - nearest layer factor
  !
  ! Output
  !   absgastot - Total absorptivity
  !   absgasnxt - Total nearest layer absorptivity
  !
  ! Nearest layer absorptivities
  ! Non-adjacent layer absorptivites
  ! Total emissivity
  !
  subroutine radabs(n1,n2,pnm,pbr,piln,pmln,tint,tlayr,co2em,co2eml,   &
                    tplnka,s2c,s2t,wh2op,h2otr,co2t,plco2,plh2o,plol,  &
                    plos,abplnk1,abplnk2,ucfc11,ucfc12,un2o0,un2o1,    &
                    bn2o0,bn2o1,uch4,bch4,uco211,uco212,uco213,uco221, &
                    uco222,uco223,uptype,absgasnxt,absgastot)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: tint , tlayr
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pnm , piln
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: pbr , pmln
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: co2em , co2eml
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: tplnka
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: s2c , s2t , wh2op
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: h2otr , co2t
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: plco2 , plh2o
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: plol , plos
    real(rkx) , dimension(nlwspi,kzp1,n1:n2) , intent(in) :: abplnk1
    real(rkx) , dimension(nlwspi,kzp1,n1:n2) , intent(in) :: abplnk2
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: ucfc11 , ucfc12
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: un2o0 , un2o1
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: bn2o0 , bn2o1
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uch4 , bch4
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco211 , uco212
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco213 , uco221
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco222 , uco223
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uptype
    real(rkx) , dimension(kzp1,kzp1,n1:n2) , intent(out) :: absgastot
    real(rkx) , dimension(kz,nlwspi,n1:n2) , intent(out) :: absgasnxt
    !
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
    ! tr10     - Eqn(6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
    ! tr2      - Eqn(6) in table A2 of R&D for 500-650
    ! tr5      - Eqn(4) in table A2 of R&D for 650-800
    ! tr6      - Eqn(4) in table A2 of R&D for 500-650
    ! tr9      - Equ(6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
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
    ! beta     - Local interface temperature
    !            (includes Voigt line correction factor)
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
    ! dbvtly   - Level drvtv plnck fnctn for o3
    !
    integer(ik4) :: n
    real(rkx) , dimension(2) :: r2st
    integer(ik4) :: k , k1 , k2 , iband , kn , wvl
    real(rkx) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,   &
      absbnd , alphat , beta , cf812 , corfac , denom , dplco2 ,       &
      dplol , dplos , ds2c , dtym10 , et , et2 , et4 , f1co2 , g2 ,    &
      g4 , k21 , k22 , o3bndi , omet , oneme , p1 , p2 , pbar , phi ,  &
      pi , posqt , psi , rbeta13 , rbeta8 , rbeta9 , rdpnm , rdpnmsq , &
      realnu , rphat , rsqti , rsum , sqwp , t1t4 , t2t5 , tcrfac ,    &
      te , tlocal , tmp1 , tmp2 , tmp3 , tpath , tr1 , tr2 , tr5 ,     &
      tr6 , tr9 , tr10 , u1 , u13 , u2 , u8 , u9 , ubar , wco2 ,       &
      dplh2o , dtp , dtz , sqti , dpnm , dtyp15 , dtyp15sq , f1sqwp ,  &
      f2co2 , f3co2 , fwk , fwku , rbeta7 , sqrtu , t1co2 , to3h2o ,   &
      tpatha , trab2 , trab4 , trab6 , u7 , uc1 , uc , ux , tco2 ,     &
      to3 , dw , abstrc , th2o , pnew , dtx , dty , to3co2
    real(rkx) :: duptyp , du1 , du2 , duch4 , dbetac , du01 , du11 ,   &
      dbeta01 , dbeta11 , duco11 , duco12 , duco13 , duco21 , duco22 , &
      duco23 , tpnm
    real(rkx) , dimension(6) :: abso
    real(rkx) , dimension(4) :: emm , o3emm , term1 , term2 , &
                      term3 , term4 , term5 , zinpl , temh2o
    real(rkx) , dimension(2) :: term7 , term8 , trline
    real(rkx) , dimension(kzp1) :: dbvtit
    real(rkx) , dimension(kzp1) :: term6
    real(rkx) , dimension(kzp1) :: term9
    real(rkx) , dimension(kzp1) :: pnmsq
    real(rkx) , dimension(kz) :: dbvtly
    real(rkx) , dimension(4) :: tbar , pinpl , uinpl , winpl
    real(rkx) , dimension(nlwspi,4) :: bplnk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radabs'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Initialize
    !
    r2st(1) = d_one/(d_two*st(1))
    r2st(2) = d_one/(d_two*st(2))
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k,k1,k2,iband,kn,wvl,  &
      a,a11,a21,a22,a23,a31,a41,a51,a61,absbnd,alphat,beta,  &
      cf812,corfac,denom,dplco2,dplol,dplos,ds2c,dtym10,et,  &
      et2,et4,f1co2,g2,g4,k21,k22,o3bndi,omet,oneme,p1,p2,   &
      pbar,phi,pi,posqt,psi,rbeta13,rbeta8,rbeta9,rdpnm,     &
      rdpnmsq,realnu,rphat,rsqti,rsum,sqwp,t1t4,t2t5,te,     &
      tcrfac,tlocal,tmp1,tmp2,tmp3,tpath,tr1,tr2,tr5,tr6,    &
      tr9,tr10,u1,u13,u2,u8,u9,ubar,wco2,dplh2o,dtp,dtz,     &
      sqti,dpnm,dtyp15,dtyp15sq,f1sqwp,f2co2,f3co2,fwk,fwku, &
      rbeta7,sqrtu,t1co2,to3h2o,tpatha,trab2,trab4,trab6,u7, &
      uc1,uc,abso,emm,o3emm,term1,term2,term3,term4,term5,   &
      term6,term9,zinpl,temh2o,term7,term8,trline,abstrc,    &
      ux,dw,to3,th2o,pnew,tco2,dtx,dty,to3co2,dbvtit,dbvtly, &
      pnmsq,tbar,pinpl,uinpl,winpl,bplnk,duptyp,du1,du2,     &
      duch4,dbetac,du01,du11,dbeta01,dbeta11,duco11,duco12,  &
      duco13,duco21,duco22,duco23,tpnm)
#else
    do n = n1 , n2
#endif
      dbvtit(kzp1) = dbvt(tint(kzp1,n))
      do k = 1 , kz
        dbvtly(k) = dbvt(tlayr(k+1,n))
        dbvtit(k) = dbvt(tint(k,n))
      end do
      !
      ! bndfct  = 2.0*22.18d0/(sqrt(196.d0)*300.)
      !
      ! Non-adjacent layer absorptivity:
      !
      ! abso(1)     0 -  800 cm-1   h2o rotation band
      ! abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
      ! abso(3)   800 - 1200 cm-1   h2o window
      ! abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
      ! abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
      ! abso(6)   co2 15  micrometer band system
      !
      do k = 1 , kzp1
        pnmsq(k) = pnm(k,n)**2
        dtx = tplnka(k,n) - 250.0_rkx
        term6(k) = coeff(1,2) + coeff(2,2)*dtx *    &
                   (d_one+c9*dtx*(d_one+c11*dtx *   &
                   (d_one+c13*dtx*(d_one+c15*dtx))))
        term9(k) = coefi(1,2) + coefi(2,2)*dtx *       &
                    (d_one+c19*dtx*(d_one+c21*dtx *    &
                    (d_one+c23*dtx*(d_one+c25*dtx))))
      end do
      !
      ! Non-nearest layer level loops
      !
      do k1 = kzp1 , 1 , -1
        do k2 = kzp1 , 1 , -1
          if ( k1 /= k2 ) then
            dplh2o = plh2o(k1,n) - plh2o(k2,n)
            ux = abs(dplh2o)
            sqrtu = sqrt(ux)
            ds2c = abs(s2c(k1,n)-s2c(k2,n))
            dw = abs(wh2op(k1,n)-wh2op(k2,n))
            uc1 = (ds2c+1.7e-3_rkx*ux) * &
                 (d_one+d_two*ds2c)/(d_one+15.0_rkx*ds2c)
            uc = ds2c + 2.0e-3_rkx*ux
            pnew = ux/dw
            tpatha = (s2t(k1,n)-s2t(k2,n))/dplh2o
            dtx = tplnka(k2,n) - 250.0_rkx
            dty = tpatha - 250.0_rkx
            dtyp15 = dty + 15.0_rkx
            dtyp15sq = dtyp15**2
            dtz = dtx - 50.0_rkx
            dtp = dty - 50.0_rkx
            do iband = 2 , 4 , 2
              term1(iband) = coefe(1,iband) + &
                         coefe(2,iband)*dtx*(d_one+c1(iband)*dtx)
              term2(iband) = coefb(1,iband) + &
                         coefb(2,iband)*dtx*(d_one+c2(iband)*dtx * &
                         (d_one+c3(iband)*dtx))
              term3(iband) = coefd(1,iband) + &
                         coefd(2,iband)*dtx*(d_one+c4(iband)*dtx * &
                         (d_one+c5(iband)*dtx))
              term4(iband) = coefa(1,iband) + &
                         coefa(2,iband)*dty*(d_one+c6(iband)*dty)
              term5(iband) = coefc(1,iband) + &
                         coefc(2,iband)*dty*(d_one+c7(iband)*dty)
            end do
            !
            ! abso(1)     0 -  800 cm-1   h2o rotation band
            !
            a11 = 0.44_rkx + 3.380e-4_rkx*dtz - 1.520e-6_rkx*dtz*dtz
            a31 = 1.05_rkx - 6.000e-3_rkx*dtp + 3.000e-6_rkx*dtp*dtp
            a21 = 1.00_rkx + 1.717e-3_rkx*dtz - 1.133e-5_rkx*dtz*dtz
            a22 = 1.00_rkx + 4.443e-3_rkx*dtp + 2.750e-5_rkx*dtp*dtp
            a23 = 1.00_rkx + 3.600_rkx*sqrtu
            corfac = a31*(a11+((d_two*a21*a22)/a23))
            t1t4 = term1(2)*term4(2)
            t2t5 = term2(2)*term5(2)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
            fwk = fwcoef + fwc1/(d_one+fwc2*ux)
            fwku = fwk*ux
            rsum = exp(-a*(sqrtu+fwku))
            abso(1) = (d_one-rsum)*term3(2)
            ! trab1(n)  = rsum
            !
            ! abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
            !
            a41 = 1.75_rkx - 3.960e-3_rkx*dtz
            a51 = 1.00_rkx + 1.3_rkx*sqrtu
            a61 = 1.00_rkx + 1.250e-3_rkx*dtp + 6.250e-5_rkx*dtp*dtp
            corfac = 0.29_rkx*(d_one+a41/a51)*a61
            t1t4 = term1(4)*term4(4)
            t2t5 = term2(4)*term5(4)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
            rsum = exp(-a*(sqrtu+fwku))
            abso(2) = (d_one-rsum)*term3(4)
            ! trab7(n)  = rsum
            !
            ! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
            !
            do k = 1 , 2
              phi = exp(a1(k)*dtyp15+a2(k)*dtyp15sq)
              psi = exp(b1(k)*dtyp15+b2(k)*dtyp15sq)
              ubar = dw*phi*1.66_rkx*r80257
              pbar = pnew*(psi/phi)
              cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
              g2 = d_one + ubar*d_four*st(k)*cf812/pbar
              g4 = realk(k)*pbar*r2st(k)*(sqrt(g2)-d_one)
              trline(k) = exp(-g4)
            end do
            term7(1) = coefj(1,1)+coefj(2,1)*dty*(d_one+c16*dty)
            term8(1) = coefk(1,1)+coefk(2,1)*dty*(d_one+c17*dty)
            term7(2) = coefj(1,2)+coefj(2,2)*dty*(d_one+c26*dty)
            term8(2) = coefk(1,2)+coefk(2,2)*dty*(d_one+c27*dty)
            !
            ! abso(3)   800 - 1200 cm-1   h2o window
            ! abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
            k21 = term7(1) + term8(1)/(d_one+(c30+c31*(dty-d_10)* &
                  (dty-d_10))*sqrtu)
            k22 = term7(2) + term8(2)/(d_one+(c28+c29*(dty-d_10))*sqrtu)
            tr1 = exp(-(k21*(sqrtu+fc1*fwku)))
            tr2 = exp(-(k22*(sqrtu+fc1*fwku)))
            tr5 = exp(-((coefh(1,3)+coefh(2,3)*dtx)*uc1))
            tr6 = exp(-((coefh(1,4)+coefh(2,4)*dtx)*uc1))
            tr9 = tr1*tr5
            tr10 = tr2*tr6
            th2o = tr10
            trab2 = 0.65_rkx*tr9 + 0.35_rkx*tr10
            trab4 = exp(-(coefg(1,3)+coefg(2,3)*dtx)*uc)
            trab6 = exp(-(coefg(1,4)+coefg(2,4)*dtx)*uc)
            abso(3) = term6(k2)*(d_one-trab4*d_half*trline(2)- &
                      trab6*d_half*trline(1))
            abso(4) = term9(k2)*d_half*(tr1-tr9+tr2-tr10)
            if ( k2 < k1 ) then
              to3h2o = h2otr(k1,n)/h2otr(k2,n)
            else
              to3h2o = h2otr(k2,n)/h2otr(k1,n)
            end if
            !
            ! abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
            !
            dpnm = pnm(k1,n) - pnm(k2,n)
            to3co2 = (pnm(k1,n)*co2t(k1,n)-pnm(k2,n)*co2t(k2,n))/dpnm
            te = (to3co2*r293)**0.7_rkx
            dplos = plos(k1,n) - plos(k2,n)
            dplol = plol(k1,n) - plol(k2,n)
            u1 = 18.29_rkx*abs(dplos)/te
            u2 = 0.5649_rkx*abs(dplos)/te
            rphat = dplol/dplos
            tlocal = tint(k2,n)
            tcrfac = sqrt(tlocal*r250)*te
            beta = r3205*(rphat+dpfo3*tcrfac)
            realnu = te/beta
            tmp1 = u1/sqrt(d_four+u1*(d_one+realnu))
            tmp2 = u2/sqrt(d_four+u2*(d_one+realnu))
            o3bndi = 74.0_rkx*te*log(d_one+tmp1+tmp2)
            abso(5) = o3bndi*to3h2o*dbvtit(k2)
            to3 = d_one/(d_one+0.1_rkx*tmp1+0.1_rkx*tmp2)
            ! trab5(n)  = d_one-(o3bndi/(1060-980.))
            !
            ! abso(6)      co2 15  micrometer band system
            !
            sqwp = sqrt(abs(plco2(k1,n)-plco2(k2,n)))
            et = exp(-480.0_rkx/to3co2)
            sqti = sqrt(to3co2)
            rsqti = d_one/sqti
            et2 = et*et
            et4 = et2*et2
            omet = d_one - 1.5_rkx*et2
            f1co2 = 899.70_rkx*omet*rsqti* &
              (d_one+1.94774_rkx*et+4.73486_rkx*et2)
            f1sqwp = f1co2*sqwp
            t1co2 = d_one/(d_one+(245.18_rkx*omet*sqwp*rsqti))
            oneme = d_one - et2
            alphat = oneme**3*rsqti
            pi = abs(dpnm)
            wco2 = 2.5221_rkx*co2vmr(n)*pi*regravgts
            u7 = 4.9411e4_rkx*alphat*et2*wco2
            u8 = 3.9744e4_rkx*alphat*et4*wco2
            u9 = 1.0447e5_rkx*alphat*et4*et2*wco2
            u13 = 2.8388e3_rkx*alphat*et4*wco2
            tpath = to3co2
            tlocal = tint(k2,n)
            tcrfac = sqrt(tlocal*r250*tpath*r300)
            posqt = ((pnm(k2,n)+pnm(k1,n))*r2sslp+dpfco2*tcrfac)*rsqti
            rbeta7 = d_one/(5.3228_rkx*posqt)
            rbeta8 = d_one/(10.6576_rkx*posqt)
            rbeta9 = rbeta7
            rbeta13 = rbeta9
            f2co2 = (u7/sqrt(d_four+u7*(d_one+rbeta7))) + &
                    (u8/sqrt(d_four+u8*(d_one+rbeta8))) + &
                    (u9/sqrt(d_four+u9*(d_one+rbeta9)))
            f3co2 = u13/sqrt(d_four+u13*(d_one+rbeta13))
            if ( k2 >= k1 ) then
              sqti = sqrt(tlayr(k2,n))
            end if

            tmp1 = log(d_one+f1sqwp)
            tmp2 = log(d_one+f2co2)
            tmp3 = log(d_one+f3co2)
            absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
            abso(6) = trab2*co2em(k2,n)*absbnd
            tco2 = d_one/(d_one+d_10*(u7/sqrt(d_four+u7*(d_one+rbeta7))))
            ! trab3(n)  = 1. - bndfct*absbnd
            !
            ! Calculate absorptivity due to trace gases
            !
            tpnm    = abs(pnm(k1,n)+pnm(k2,n))
            duptyp  = abs(uptype(k1,n)-uptype(k2,n))
            du1     = abs(ucfc11(k1,n)-ucfc11(k2,n))
            du2     = abs(ucfc12(k1,n)-ucfc12(k2,n))
            duch4   = abs(uch4(k1,n)-uch4(k2,n))
            dbetac  = abs(bch4(k1,n)-bch4(k2,n))/duch4
            du01    = abs(un2o0(k1,n)-un2o0(k2,n))
            du11    = abs(un2o1(k1,n)-un2o1(k2,n))
            dbeta01 = abs(bn2o0(k1,n)-bn2o0(k2,n))/du01
            dbeta11 = abs(bn2o1(k1,n)-bn2o1(k2,n))/du11
            duco11  = abs(uco211(k1,n)-uco211(k2,n))
            duco12  = abs(uco212(k1,n)-uco212(k2,n))
            duco13  = abs(uco213(k1,n)-uco213(k2,n))
            duco21  = abs(uco221(k1,n)-uco221(k2,n))
            duco22  = abs(uco222(k1,n)-uco222(k2,n))
            duco23  = abs(uco223(k1,n)-uco223(k2,n))
            abstrc = trcab(tpnm,ds2c,duptyp,du1,du2,duch4,dbetac,  &
                           du01,du11,dbeta01,dbeta11,duco11,duco12, &
                           duco13,duco21,duco22,duco23,dw,pnew,     &
                           to3co2,ux,tco2,th2o,to3,abplnk1(:,k2,n))
            !
            ! Sum total absorptivity
            !
            absgastot(k1,k2,n) = abso(1) + abso(2) + abso(3) + &
                                 abso(4) + abso(5) + abso(6) + abstrc
          end if
        end do
      end do  ! End of non-nearest layer level loops
      !
      ! Non-adjacent layer absorptivity:
      !
      ! abso(1)     0 -  800 cm-1   h2o rotation band
      ! abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
      ! abso(3)   800 - 1200 cm-1   h2o window
      ! abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
      ! abso(5)   o3  9.6 micrometer band (nu3 and nu1 bands)
      ! abso(6)   co2 15  micrometer band system
      !
      ! Nearest layer level loop
      !
      do k2 = kz , 1 , -1
        tbar(1) = (tint(k2+1,n)+tlayr(k2+1,n))*d_half
        emm(1) = (co2em(k2+1,n)+co2eml(k2,n))*d_half
        tbar(2) = (tlayr(k2+1,n)+tint(k2,n))*d_half
        emm(2) = (co2em(k2,n)+co2eml(k2,n))*d_half
        tbar(3) = (tbar(2)+tbar(1))*d_half
        emm(3) = emm(1)
        tbar(4) = tbar(3)
        emm(4) = emm(2)
        o3emm(1) = (dbvtit(k2+1)+dbvtly(k2))*d_half
        o3emm(2) = (dbvtit(k2)+dbvtly(k2))*d_half
        o3emm(3) = o3emm(1)
        o3emm(4) = o3emm(2)
        temh2o(1) = tbar(1)
        temh2o(2) = tbar(2)
        temh2o(3) = tbar(1)
        temh2o(4) = tbar(2)
        dpnm = pnm(k2+1,n) - pnm(k2,n)
        !
        ! Weighted Planck functions for trace gases
        !
        do wvl = 1 , nlwspi
          bplnk(wvl,1) = (abplnk1(wvl,k2+1,n)+abplnk2(wvl,k2,n))*d_half
          bplnk(wvl,2) = (abplnk1(wvl,k2,  n)+abplnk2(wvl,k2,n))*d_half
          bplnk(wvl,3) = bplnk(wvl,1)
          bplnk(wvl,4) = bplnk(wvl,2)
        end do
        rdpnmsq = d_one/(pnmsq(k2+1)-pnmsq(k2))
        rdpnm = d_one/dpnm
        p1 = (pbr(k2,n)+pnm(k2+1,n))*d_half
        p2 = (pbr(k2,n)+pnm(k2,n))*d_half
        uinpl(1) = (pnmsq(k2+1)-p1**2)*rdpnmsq
        uinpl(2) = -(pnmsq(k2)-p2**2)*rdpnmsq
        uinpl(3) = -(pnmsq(k2)-p1**2)*rdpnmsq
        uinpl(4) = (pnmsq(k2+1)-p2**2)*rdpnmsq
        winpl(1) = ((pnm(k2+1,n)-pbr(k2,n))*d_half)*rdpnm
        winpl(2) = ((-pnm(k2,n)+pbr(k2,n))*d_half)*rdpnm
        winpl(3) = ((pnm(k2+1,n)+pbr(k2,n))*d_half-pnm(k2,n))*rdpnm
        winpl(4) = ((-pnm(k2,n)-pbr(k2,n))*d_half+pnm(k2+1,n))*rdpnm
        tmp1 = d_one/(piln(k2+1,n)-piln(k2,n))
        tmp2 = piln(k2+1,n) - pmln(k2,n)
        tmp3 = piln(k2,n)   - pmln(k2,n)
        zinpl(1) = (tmp2*d_half)*tmp1
        zinpl(2) = (-tmp3*d_half)*tmp1
        zinpl(3) = (tmp2*d_half-tmp3)*tmp1
        zinpl(4) = (tmp2-tmp3*d_half)*tmp1
        pinpl(1) = (p1+pnm(k2+1,n))*d_half
        pinpl(2) = (p2+pnm(k2,n))*d_half
        pinpl(3) = (p1+pnm(k2,n))*d_half
        pinpl(4) = (p2+pnm(k2+1,n))*d_half
        ! FAB AER SAVE uinpl  for aerosl LW forcing calculation
        if ( linteract  ) then
          do kn = 1 , 4
            xuinpl(k2,kn,n) = uinpl(kn)
          end do
        end if
        ! FAB AER SAVE uinpl  for aerosl LW forcing calculation
        do kn = 1 , 4
          ux = abs(uinpl(kn)*(plh2o(k2,n)-plh2o(k2+1,n)))
          sqrtu = sqrt(ux)
          dw = abs(wh2op(k2,n)-wh2op(k2+1,n))
          pnew = ux/(winpl(kn)*dw)
          ds2c = abs(s2c(k2,n)-s2c(k2+1,n))
          uc1 = uinpl(kn)*ds2c
          uc1 = (uc1+1.7e-3_rkx*ux)*(d_one+d_two*uc1)/&
                (d_one+15.0_rkx*uc1)
          uc = uinpl(kn)*ds2c + 2.0e-3_rkx*ux
          dtx = temh2o(kn) - 250.0_rkx
          dty = tbar(kn) - 250.0_rkx
          dtyp15 = dty + 15.0_rkx
          dtyp15sq = dtyp15**2
          dtz = dtx - 50.0_rkx
          dtp = dty - 50.0_rkx
          do iband = 2 , 4 , 2
            term1(iband) = coefe(1,iband) + coefe(2,iband)*dtx * &
                           (d_one+c1(iband)*dtx)
            term2(iband) = coefb(1,iband) + coefb(2,iband)*dtx * &
                           (d_one+c2(iband)*dtx                * &
                           (d_one+c3(iband)*dtx))
            term3(iband) = coefd(1,iband) + coefd(2,iband)*dtx * &
                           (d_one+c4(iband)*dtx                * &
                           (d_one+c5(iband)*dtx))
            term4(iband) = coefa(1,iband) + coefa(2,iband)*dty * &
                           (d_one+c6(iband)*dty)
            term5(iband) = coefc(1,iband) + coefc(2,iband)*dty * &
                           (d_one+c7(iband)*dty)
          end do
          !
          ! abso(1)     0 -  800 cm-1   h2o rotation band
          !
          a11 = 0.44_rkx + 3.380e-4_rkx*dtz - 1.520e-6_rkx*dtz*dtz
          a31 = 1.05_rkx - 6.000e-3_rkx*dtp + 3.000e-6_rkx*dtp*dtp
          a21 = 1.00_rkx + 1.717e-3_rkx*dtz - 1.133e-5_rkx*dtz*dtz
          a22 = 1.00_rkx + 4.443e-3_rkx*dtp + 2.750e-5_rkx*dtp*dtp
          a23 = 1.00_rkx + 3.600_rkx*sqrtu
          corfac = a31*(a11+((d_two*a21*a22)/a23))
          t1t4 = term1(2)*term4(2)
          t2t5 = term2(2)*term5(2)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
          fwk = fwcoef + fwc1/(d_one+fwc2*ux)
          fwku = fwk*ux
          rsum = exp(-a*(sqrtu+fwku))
          abso(1) = (d_one-rsum)*term3(2)
          ! trab1(n) = rsum
          !
          ! abso(2)  1200 - 2200 cm-1   h2o vibration-rotation band
          !
          a41 = 1.75_rkx - 3.960e-3_rkx*dtz
          a51 = 1.00_rkx + 1.3_rkx*sqrtu
          a61 = 1.00_rkx + 1.250e-3_rkx*dtp + 6.250e-5_rkx*dtp*dtp
          corfac = 0.29_rkx*(d_one+a41/a51)*a61
          t1t4 = term1(4)*term4(4)
          t2t5 = term2(4)*term5(4)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu*corfac)
          rsum = exp(-a*(sqrtu+fwku))
          abso(2) = (d_one-rsum)*term3(4)
          ! trab7(n) = rsum
          !
          ! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
          !
          do k = 1 , 2
            phi = exp(a1(k)*dtyp15+a2(k)*dtyp15sq)
            psi = exp(b1(k)*dtyp15+b2(k)*dtyp15sq)
            ubar = dw*phi*winpl(kn)*1.66_rkx*r80257
            pbar = pnew*(psi/phi)
            cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
            g2 = d_one + ubar*d_four*st(k)*cf812/pbar
            g4 = realk(k)*pbar*r2st(k)*(sqrt(g2)-d_one)
            trline(k) = exp(-g4)
          end do
          term7(1) = coefj(1,1)+coefj(2,1)*dty*(d_one+c16*dty)
          term8(1) = coefk(1,1)+coefk(2,1)*dty*(d_one+c17*dty)
          term7(2) = coefj(1,2)+coefj(2,2)*dty*(d_one+c26*dty)
          term8(2) = coefk(1,2)+coefk(2,2)*dty*(d_one+c27*dty)
          !
          ! abso(3)   800 - 1200 cm-1   h2o window
          ! abso(4)   500 -  800 cm-1   h2o rotation band overlap with co2
          !
          dtym10 = dty - d_10
          denom = d_one + (c30+c31*dtym10*dtym10)*sqrtu
          k21 = term7(1) + term8(1)/denom
          denom = d_one + (c28+c29*dtym10)*sqrtu
          k22 = term7(2) + term8(2)/denom
          term9(2) = coefi(1,2) + coefi(2,2)*dtx *     &
                     (d_one+c19*dtx*(d_one+c21*dtx *   &
                     (d_one+c23*dtx*(d_one+c25*dtx))))
          tr1 = exp(-(k21*(sqrtu+fc1*fwku)))
          tr2 = exp(-(k22*(sqrtu+fc1*fwku)))
          tr5 = exp(-((coefh(1,3)+coefh(2,3)*dtx)*uc1))
          tr6 = exp(-((coefh(1,4)+coefh(2,4)*dtx)*uc1))
          tr9 = tr1*tr5
          tr10 = tr2*tr6
          trab2 = 0.65_rkx*tr9 + 0.35_rkx*tr10
          th2o = tr10
          trab4 = exp(-(coefg(1,3)+coefg(2,3)*dtx)*uc)
          trab6 = exp(-(coefg(1,4)+coefg(2,4)*dtx)*uc)
          term6(2) = coeff(1,2) + coeff(2,2)*dtx *  &
                       (d_one+c9*dtx*(d_one+c11*dtx * &
                       (d_one+c13*dtx*(d_one+c15*dtx))))
          abso(3) = term6(2)*(d_one-trab4*d_half*trline(2) - &
                                    trab6*d_half*trline(1))
          abso(4) = term9(2)*d_half*(tr1-tr9+tr2-tr10)
          !
          ! abso(5)  o3  9.6 micrometer (nu3 and nu1 bands)
          !
          te = (tbar(kn)*r293)**0.7_rkx
          dplos = abs(plos(k2+1,n)-plos(k2,n))
          u1 = zinpl(kn)*18.29_rkx*dplos/te
          u2 = zinpl(kn)*0.5649_rkx*dplos/te
          tlocal = tbar(kn)
          tcrfac = sqrt(tlocal*r250)*te
          beta = r3205*(pinpl(kn)*rsslp+dpfo3*tcrfac)
          realnu = te/beta
          tmp1 = u1/sqrt(d_four+u1*(d_one+realnu))
          tmp2 = u2/sqrt(d_four+u2*(d_one+realnu))
          o3bndi = 74.0_rkx*te*log(d_one+tmp1+tmp2)
          abso(5) = o3bndi*o3emm(kn)*(h2otr(k2+1,n)/h2otr(k2,n))
          to3 = d_one/(d_one+0.1_rkx*tmp1+0.1_rkx*tmp2)
          ! trab5(n) = d_one-(o3bndi/(1060-980.))
          !
          ! abso(6)   co2 15  micrometer band system
          !
          dplco2 = plco2(k2+1,n) - plco2(k2,n)
          sqwp = sqrt(uinpl(kn)*dplco2)
          et = exp(-480.0_rkx/tbar(kn))
          sqti = sqrt(tbar(kn))
          rsqti = d_one/sqti
          et2 = et*et
          et4 = et2*et2
          omet = (d_one-1.5_rkx*et2)
          f1co2 = 899.70_rkx*omet*rsqti*(d_one+1.94774_rkx*et+4.73486_rkx*et2)
          f1sqwp = f1co2*sqwp
          t1co2 = d_one/(d_one+(245.18_rkx*omet*sqwp*rsqti))
          oneme = d_one - et2
          alphat = oneme**3*rsqti
          pi = abs(dpnm)*winpl(kn)
          wco2 = 2.5221_rkx*co2vmr(n)*pi*regravgts
          u7 = 4.9411e4_rkx*alphat*et2*wco2
          u8 = 3.9744e4_rkx*alphat*et4*wco2
          u9 = 1.0447e5_rkx*alphat*et4*et2*wco2
          u13 = 2.8388e3_rkx*alphat*et4*wco2
          tpath = tbar(kn)
          tlocal = tbar(kn)
          tcrfac = sqrt((tlocal*r250)*(tpath*r300))
          posqt = (pinpl(kn)*rsslp+dpfco2*tcrfac)*rsqti
          rbeta7 = d_one/(5.3228_rkx*posqt)
          rbeta8 = d_one/(10.6576_rkx*posqt)
          rbeta9 = rbeta7
          rbeta13 = rbeta9
          f2co2 = u7/sqrt(d_four+u7*(d_one+rbeta7)) + &
                  u8/sqrt(d_four+u8*(d_one+rbeta8)) + &
                  u9/sqrt(d_four+u9*(d_one+rbeta9))
          f3co2 = u13/sqrt(d_four+u13*(d_one+rbeta13))
          tmp1 = log(d_one+f1sqwp)
          tmp2 = log(d_one+f2co2)
          tmp3 = log(d_one+f3co2)
          absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
          abso(6) = trab2*emm(kn)*absbnd
          tco2 = d_one/(d_one+d_10*u7/sqrt(d_four+u7*(d_one+rbeta7)))
          ! trab3(n) = 1. - bndfct*absbnd
          !
          ! Calculate trace gas absorptivity for nearest layer
          !
          ds2c   = abs(s2c(k2+1,n)-s2c(k2,n))*uinpl(kn)
          duptyp = abs(uptype(k2+1,n)-uptype(k2,n))*uinpl(kn)
          du1    = abs(ucfc11(k2+1,n)-ucfc11(k2,n))*winpl(kn)
          du2    = abs(ucfc12(k2+1,n)-ucfc12(k2,n))*winpl(kn)
          duch4  = abs(uch4(k2+1,n)-uch4(k2,n))*winpl(kn)
          du01   = abs(un2o0(k2+1,n)-un2o0(k2,n))*winpl(kn)
          du11   = abs(un2o1(k2+1,n)-un2o1(k2,n))*winpl(kn)
          duco11 = abs(uco211(k2+1,n)-uco211(k2,n))*winpl(kn)
          duco12 = abs(uco212(k2+1,n)-uco212(k2,n))*winpl(kn)
          duco13 = abs(uco213(k2+1,n)-uco213(k2,n))*winpl(kn)
          duco21 = abs(uco221(k2+1,n)-uco221(k2,n))*winpl(kn)
          duco22 = abs(uco222(k2+1,n)-uco222(k2,n))*winpl(kn)
          duco23 = abs(uco223(k2+1,n)-uco223(k2,n))*winpl(kn)
          abstrc = trcabn(tbar(kn),dw,pnew,tco2,th2o,to3,ux,pinpl(kn),   &
                          winpl(kn),ds2c,duptyp,du1,du2,duch4,du01,du11, &
                          duco11,duco12,duco13,duco21,duco22,duco23,     &
                          bplnk(:,kn))
          !
          ! Total next layer absorptivity:
          !
          absgasnxt(k2,kn,n) = abso(1) + abso(2) + abso(3) + &
                               abso(4) + abso(5) + abso(6) + abstrc
        end do
      end do  !  end of nearest layer level loop
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine radabs
  !
  !-----------------------------------------------------------------------
  !
  ! Compute emissivity for H2O, CO2, O3, CH4, N2O, CFC11 and CFC12
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
  ! CH4  ....  Uses a broad band model for the 7.7 micron band of methane.
  !
  ! N20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
  !            bands of nitrous oxide
  !
  ! CFC11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
  !            micron bands of CFC11
  !
  ! CFC12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
  !            micron bands of CFC12
  !
  ! Computes individual emissivities, accounting for band overlap, and
  ! sums to obtain the total.
  !
  ! Input arguments
  !
  !   s2c     - H2o continuum path length
  !   s2t     - Tmp and prs wghted h2o path length
  !   w       - H2o path length
  !   tplnke  - Layer planck temperature
  !   plh2o   - H2o prs wghted path length
  !   pnm     - Model interface pressure (dynes/cm*2)
  !   plco2   - Prs wghted path of co2
  !   tint    - Model interface temperatures
  !   tint4   - Tint to the 4th power
  !   tlayr   - K-1 model layer temperature
  !   tlayr4  - Tlayr to the 4th power
  !   plol    - Pressure wghtd ozone path
  !   plos    - Ozone path
  !   ucfc11  - CFC11 path length
  !   ucfc12  - CFC12 path length
  !   un2o0   - N2O path length
  !   un2o1   - N2O path length (hot band)
  !   uch4    - CH4 path length
  !   uco211  - CO2 9.4 micron band path length
  !   uco212  - CO2 9.4 micron band path length
  !   uco213  - CO2 9.4 micron band path length
  !   uco221  - CO2 10.4 micron band path length
  !   uco222  - CO2 10.4 micron band path length
  !   uco223  - CO2 10.4 micron band path length
  !   bn2o0   - pressure factor for n2o
  !   bn2o1   - pressure factor for n2o
  !   bch4    - pressure factor for ch4
  !   uptype  - p-type continuum path length
  !
  ! Output arguments
  !
  !   co2em   - Layer co2 normalzd plnck funct drvtv
  !   co2eml  - Intrfc co2 normalzd plnck func drvtv
  !   co2t    - Tmp and prs weighted path length
  !   h2otr   - H2o transmission over o3 band
  !   emplnk  - emissivity Planck factor
  !   emstrc  - total trace gas emissivity
  !
  subroutine radems(n1,n2,pnm,tint,tint4,tlayr,tlayr4,tplnke,plos,plol,    &
                    plh2o,ucfc11,ucfc12,un2o0,un2o1,bn2o0,bn2o1,uch4,bch4, &
                    uco211,uco212,uco213,uco221,uco222,uco223,uptype,      &
                    wh2op,s2c,s2t,emplnk,co2t,co2em,co2eml,h2otr,emsgastot)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pnm
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: tint , tint4
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: tlayr , tlayr4
    real(rkx) , dimension(n1:n2) , intent(in) :: tplnke
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: plol , plos , plh2o
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: ucfc11 , ucfc12
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: un2o0 , un2o1
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: bn2o0 , bn2o1
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uch4 , bch4
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco211 , uco212
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco213 , uco221
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uco222 , uco223
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: uptype
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: wh2op , s2c , s2t
    real(rkx) , dimension(nlwspi,n1:n2) , intent(in) :: emplnk
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: co2t , co2em , co2eml
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: h2otr
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: emsgastot
    !
    ! iband   - H2o band index
    !
    ! Local variables for H2O:
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
    ! term3   - B(T) function for rotation and vibration-rotation band emiss.
    ! term4   - Equation(6) in table A3a of R&D(1986)
    ! term5   - Delta a(Tp) in table A3a of R&D(1986)
    ! term6   - B(T) function for window region
    ! term7   - Kl_inf(i) in eq(8) of table A3a of R&D
    ! term8   - Delta kl_inf(i) in eq(8)
    ! term9   - B(T) function for 500-800 cm-1 region
    ! tr1     - Equation(6) in table A2 for 650-800
    ! tr2     - Equation(6) in table A2 for 500-650
    ! tr3     - Equation(4) in table A2 for 650-800
    ! tr4     - Equation(4),table A2 of R&D for 500-650
    ! tr7     - Equ. (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
    ! tr8     - Equ. (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
    ! uc      - Y + 0.002U in eq(8) of table A2 of R&D
    ! pnew    - Effective pressure for h2o linewidth
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
    ! Local variables for CO2:
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
    ! Local variables for O3:
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
    ! Transmission terms for various spectral intervals:
    !
    ! trem4   - H2o   800 - 1000 cm-1
    ! trem6   - H2o  1000 - 1200 cm-1
    ! absbnd  - Proportional to co2 band absorptance
    ! tco2    - co2 overlap factor
    ! th2o    - h2o overlap factor
    ! to3     - o3 overlap factor
    !
    integer(ik4) :: n
    real(rkx) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,  &
                 absbnd , alphat , beta , cf812 , et , et2 , et4 , ex , &
                 exm1sq , f1co2 , f1sqwp , f2co2 , f3co2 , fwk , g1 ,   &
                 g2 , g3 , g4 , o3bndi , omet , oneme , pbar , phat ,   &
                 phi , pi , posqt , psi , k21 , k22 , trem4 , trem6 ,   &
                 rbeta13 , rbeta7 , rbeta8 , rbeta9 , realnu , rsqti ,  &
                 sqti , sqwp , t1co2 , t1i , t1t4 , t2t5 , tpathe ,     &
                 tcrfac , te , tlayr5 , tlocal , tmp1 , tmp2 , tmp3 ,   &
                 tpath , u1 , u13 , u2 , u7 , u8 , u9 , ubar , wco2 ,   &
                 tr1 , tr2 , tr3 , tr4 , tr7 , tr8 , corfac , dbvtt ,   &
                 dtp , dtz , pnew , rsum , uc , uc1 , ux , troco2 ,     &
                 tco2 , to3 , th2o , emstrc , h2oems , co2ems , o3ems , &
                 xsum , dtx , dty , co2plk
    real(rkx) , dimension(4) :: term1 , term2 , term3 , term4 , term5
    real(rkx) , dimension(4) :: emis
    real(rkx) :: term6 , term9
    real(rkx) , dimension(2) :: term7 , term8 , trline
    integer(ik4) :: k , kk , iband , l
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radems'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k,kk,iband,l,term6,term9,emis,   &
      a,a11,a21,a22,a23,a31,a41,a51,a61,absbnd,alphat,beta,cf812,et,   &
      et2,et4,ex,exm1sq,f1co2,f1sqwp,f2co2,f3co2,fwk,g1,g2,g3,g4,omet, &
      o3bndi,oneme,pbar,phat,phi,pi,posqt,psi,k21,k22,trem4,trem6,     &
      rbeta13,rbeta7,rbeta8,rbeta9,realnu,rsqti,sqti,sqwp,t1co2,t1i,   &
      t1t4,t2t5,tpathe,tcrfac,te,tlayr5,tlocal,tmp1,tmp2,tmp3,tpath,   &
      u1,u13,u2,u7,u8,u9,ubar,wco2,tr1,tr2,tr3,tr4,tr7,tr8,corfac,     &
      dbvtt,dtp,dtz,pnew,rsum,uc,uc1,troco2,term1,term2,term3,term4,   &
      term5,term7,term8,trline,ux,tco2,th2o,to3,emstrc,h2oems,co2ems,  &
      o3ems,xsum,dtx,dty,co2plk)
#else
    do n = n1 , n2
#endif

      ex = exp(960.0_rkx/tplnke(n))
      co2plk = 5.0e8_rkx/((tplnke(n)**4)*(ex-d_one))
      co2t(1,n) = tplnke(n)
      xsum = co2t(1,n)*pnm(1,n)
      kk = 1
      do k = kzp1 , 2 , -1
        kk = kk + 1
        xsum = xsum + tlayr(kk,n)*(pnm(kk,n)-pnm(kk-1,n))
        ex = exp(960.0_rkx/tlayr(kk,n))
        tlayr5 = tlayr(kk,n)*tlayr4(kk,n)
        co2eml(kk-1,n) = 1.2e11_rkx*ex/(tlayr5*(ex-d_one)**2)
        co2t(kk,n) = xsum/pnm(kk,n)
      end do
      !
      ! bndfct = 2.d0*22.18/(sqrt(196.d0)*300.)
      ! Interface loop
      !
      do k = 1 , kzp1
        !
        ! H2O emissivity
        !
        ! emis(1)     0 -  800 cm-1   rotation band
        ! emis(2)  1200 - 2200 cm-1   vibration-rotation band
        ! emis(3)   800 - 1200 cm-1   window
        ! emis(4)   500 -  800 cm-1   rotation band overlap with co2
        !
        ! For the p type continuum
        !
        ux = plh2o(k,n)
        uc = s2c(k,n) + 2.0e-3_rkx*ux
        pnew = ux/wh2op(k,n)
        !
        ! Apply scaling factor for 500-800 continuum
        !
        uc1 = (s2c(k,n)+1.7e-3_rkx*plh2o(k,n)) * &
              (d_one+d_two*s2c(k,n))/(d_one+15.0_rkx*s2c(k,n))
        tpathe = s2t(k,n)/plh2o(k,n)
        dtx = tplnke(n) - 250.0_rkx
        dty = tpathe - 250.0_rkx
        dtz = dtx - 50.0_rkx
        dtp = dty - 50.0_rkx
        !
        ! emis(1)     0 -  800 cm-1   rotation band
        !
        do iband = 1 , 3 , 2
          term1(iband) = coefe(1,iband) + coefe(2,iband)*dtx *   &
                         (d_one+c1(iband)*dtx)
          term2(iband) = coefb(1,iband) + coefb(2,iband)*dtx *   &
                         (d_one+c2(iband)*dtx*(d_one+c3(iband)*dtx))
          term3(iband) = coefd(1,iband) + coefd(2,iband)*dtx *   &
                         (d_one+c4(iband)*dtx*(d_one+c5(iband)*dtx))
          term4(iband) = coefa(1,iband) + coefa(2,iband)*dty *   &
                         (d_one+c6(iband)*dty)
          term5(iband) = coefc(1,iband) + coefc(2,iband)*dty *   &
                         (d_one+c7(iband)*dty)
        end do
        a11 = 0.37_rkx - 3.33e-5_rkx*dtz + 3.33e-6_rkx*dtz*dtz
        a31 = 1.07_rkx - 1.00e-3_rkx*dtp + 1.475e-5_rkx*dtp*dtp
        a21 = 1.3870_rkx + 3.80e-3_rkx*dtz - 7.8e-6_rkx*dtz*dtz
        a22 = d_one - 1.21e-3_rkx*dtp - 5.33e-6_rkx*dtp*dtp
        a23 = 0.9_rkx + 2.62_rkx*sqrt(ux)
        corfac = a31*(a11+((a21*a22)/a23))
        t1t4 = term1(1)*term4(1)
        t2t5 = term2(1)*term5(1)
        a = t1t4 + t2t5/(d_one+t2t5*sqrt(ux)*corfac)
        fwk = fwcoef + fwc1/(d_one+fwc2*ux)
        rsum = exp(-a*(sqrt(ux)+fwk*ux))
        emis(1) = (d_one-rsum)*term3(1)
        ! trem1  = rsum
        !
        ! emis(2)  1200 - 2200 cm-1   vibration-rotation band
        !
        a41 = 1.75_rkx - 3.96e-3_rkx*dtz
        a51 = 1.00_rkx + 1.3_rkx*sqrt(ux)
        a61 = 1.00_rkx + 1.25e-3_rkx*dtp + 6.25e-5_rkx*dtp*dtp
        corfac = 0.3_rkx*(d_one+(a41)/(a51))*a61
        t1t4 = term1(3)*term4(3)
        t2t5 = term2(3)*term5(3)
        a = t1t4 + t2t5/(d_one+t2t5*sqrt(ux)*corfac)
        fwk = fwcoef + fwc1/(d_one+fwc2*ux)
        rsum = exp(-a*(sqrt(ux)+fwk*ux))
        emis(2) = (d_one-rsum)*term3(3)
        ! trem7 = rsum
        !
        ! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
        !
        ! emis(3)   800 - 1200 cm-1   window
        !
        do l = 1 , 2
          phi = a1(l)*(dty+15.0_rkx)+a2(l)*(dty+15.0_rkx)**2
          psi = b1(l)*(dty+15.0_rkx)+b2(l)*(dty+15.0_rkx)**2
          phi = exp(phi)
          psi = exp(psi)
          ubar = wh2op(k,n)*phi
          ubar = (ubar*1.66_rkx)*r80257
          pbar = pnew*(psi/phi)
          cf812 = cfa1 + ((d_one-cfa1)/(d_one+ubar*pbar*d_10))
          g1 = (realk(l)*pbar)/(d_two*st(l))
          g2 = d_one + (ubar*d_four*st(l)*cf812)/pbar
          g3 = sqrt(g2) - d_one
          g4 = g1*g3
          trline(l) = exp(-g4)
        end do
        term6 = coeff(1,1) + coeff(2,1)*dtx *     &
                (d_one+c8*dtx*(d_one+c10*dtx *    &
                (d_one+c12*dtx*(d_one+c14*dtx))))
        term7(1) = coefj(1,1)+coefj(2,1)*dty*(d_one+c16*dty)
        term8(1) = coefk(1,1)+coefk(2,1)*dty*(d_one+c17*dty)
        term7(2) = coefj(1,2)+coefj(2,2)*dty*(d_one+c26*dty)
        term8(2) = coefk(1,2)+coefk(2,2)*dty*(d_one+c27*dty)
        trem4 = exp(-(coefg(1,1)+coefg(2,1)*dtx)*uc)*trline(2)
        trem6 = exp(-(coefg(1,2)+coefg(2,2)*dtx)*uc)*trline(1)
        emis(3) = term6*(d_one-trem4*d_half-trem6*d_half)
        !
        ! emis(4)   500 -  800 cm-1   rotation band overlap with co2
        !
        k21 = term7(1) + term8(1)/(d_one+(c30+c31*(dty-d_10) * &
                 (dty-d_10))*sqrt(ux))
        k22 = term7(2) + term8(2)/(d_one+(c28+c29*(dty-d_10))*sqrt(ux))
        term9 = coefi(1,1) + coefi(2,1)*dtx *  &
                (d_one+c18*dtx*(d_one+c20*dtx * &
                (d_one+c22*dtx*(d_one+c24*dtx))))
        fwk = fwcoef + fwc1/(d_one+fwc2*ux)
        tr1 = exp(-(k21*(sqrt(ux)+fc1*fwk*ux)))
        tr2 = exp(-(k22*(sqrt(ux)+fc1*fwk*ux)))
        tr3 = exp(-((coefh(1,1)+coefh(2,1)*dtx)*uc1))
        tr4 = exp(-((coefh(1,2)+coefh(2,2)*dtx)*uc1))
        tr7 = tr1*tr3
        tr8 = tr2*tr4
        emis(4) = term9*d_half*(tr1-tr7+tr2-tr8)
        h2oems = emis(1) + emis(2) + emis(3) + emis(4)
        troco2 = 0.65_rkx*tr7 + 0.35_rkx*tr8
        th2o = tr8
        ! trem2(n) = troco2
        !
        ! CO2 emissivity for 15 micron band system
        !
        t1i = exp(-480.0_rkx/co2t(k,n))
        sqti = sqrt(co2t(k,n))
        rsqti = d_one/sqti
        et = t1i
        et2 = et*et
        et4 = et2*et2
        omet = d_one - 1.5_rkx*et2
        f1co2 = 899.70_rkx*omet*(d_one+1.94774_rkx*et+4.73486_rkx*et2)*rsqti
        sqwp = sqrt(plco2(k,n))
        f1sqwp = f1co2*sqwp
        t1co2 = d_one/(d_one+245.18_rkx*omet*sqwp*rsqti)
        oneme = d_one - et2
        alphat = oneme**3*rsqti
        wco2 = 2.5221_rkx*co2vmr(n)*pnm(k,n)*regravgts
        u7 = 4.9411e4_rkx*alphat*et2*wco2
        u8 = 3.9744e4_rkx*alphat*et4*wco2
        u9 = 1.0447e5_rkx*alphat*et4*et2*wco2
        u13 = 2.8388e3_rkx*alphat*et4*wco2

        tpath = co2t(k,n)
        tlocal = tplnke(n)
        tcrfac = sqrt((tlocal*r250)*(tpath*r300))
        pi = pnm(k,n)*rsslp + d_two*dpfco2*tcrfac
        posqt = pi/(d_two*sqti)
        rbeta7 = d_one/(5.3288_rkx*posqt)
        rbeta8 = d_one/(10.6576_rkx*posqt)
        rbeta9 = rbeta7
        rbeta13 = rbeta9
        f2co2 = (u7/sqrt(d_four+u7*(d_one+rbeta7))) + &
                (u8/sqrt(d_four+u8*(d_one+rbeta8))) + &
                (u9/sqrt(d_four+u9*(d_one+rbeta9)))
        f3co2 = u13/sqrt(d_four+u13*(d_one+rbeta13))
        tmp1 = log(d_one+f1sqwp)
        tmp2 = log(d_one+f2co2)
        tmp3 = log(d_one+f3co2)
        absbnd = (tmp1+d_two*t1co2*tmp2+d_two*tmp3)*sqti
        tco2 = d_one/(d_one+d_10*(u7/sqrt(d_four+u7*(d_one+rbeta7))))
        co2ems = troco2*absbnd*co2plk
        ex = exp(960.0_rkx/tint(k,n))
        exm1sq = (ex-d_one)**2
        co2em(k,n) = 1.2e11_rkx*ex/(tint(k,n)*tint4(k,n)*exm1sq)
        ! trem3(n) = 1. - bndfct*absbnd
        !
        ! O3 emissivity
        !
        h2otr(k,n) = exp(-12.0_rkx*s2c(k,n))
        te = (co2t(k,n)/293.0_rkx)**0.7_rkx
        u1 = 18.29_rkx*plos(k,n)/te
        u2 = 0.5649_rkx*plos(k,n)/te
        phat = plos(k,n)/plol(k,n)
        tlocal = tplnke(n)
        tcrfac = sqrt(tlocal*r250)*te
        beta = (d_one/0.3205_rkx)*((d_one/phat)+(dpfo3*tcrfac))
        realnu = (d_one/beta)*te
        o3bndi = 74.0_rkx*te*(tplnke(n)/375.0_rkx)* &
                 log(d_one+fo3(u1,realnu)+fo3(u2,realnu))
        dbvtt = dbvt(tplnke(n))
        o3ems = dbvtt*h2otr(k,n)*o3bndi
        to3 = d_one/(d_one+0.1_rkx*fo3(u1,realnu)+0.1_rkx*fo3(u2,realnu))
        ! trem5(n)    = d_one-(o3bndi/(1060-980.))
        !
        ! Calculate trace gas emissivities
        !
        emstrc = trcems(co2t(k,n),pnm(k,n),ucfc11(k,n),ucfc12(k,n),  &
                        un2o0(k,n),un2o1(k,n),bn2o0(k,n),bn2o1(k,n), &
                        uch4(k,n),bch4(k,n),uco211(k,n),uco212(k,n), &
                        uco213(k,n),uco221(k,n),uco222(k,n),         &
                        uco223(k,n),uptype(k,n),wh2op(k,n),s2c(k,n), &
                        ux,emplnk(:,n),th2o,tco2,to3)
        !
        ! Total emissivity:
        !
        emsgastot(k,n) = h2oems + co2ems + o3ems + emstrc
      end do  ! End of interface loop
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine radems
  !
  !-----------------------------------------------------------------------
  !
  ! Computes the path length integrals to the model interfaces given the
  ! ozone volume mixing ratio
  !
  ! Input
  !   o3vmr - ozone volume mixing ratio
  !   pnm  - Model interface pressures
  ! Output
  !   plol  - Ozone prs weighted path length (cm)
  !   plos  - Ozone path length (cm)
  !
  !-----------------------------------------------------------------------
  !
  subroutine radoz2(n1,n2,o3vmr,pnm,plos,plol)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: o3vmr
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pnm
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: plos , plol
    integer(ik4) :: n
    integer(ik4) :: k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radoz2'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Evaluate the ozone path length integrals to interfaces;
    ! factors of 0.1 and 0.01 to convert pressures from cgs to mks:
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k)
#else
    do n = n1 , n2
#endif
      plos(1,n) = 0.1_rkx*cplos*o3vmr(1,n)*pnm(1,n)
      plol(1,n) = 0.01_rkx*cplol*o3vmr(1,n)*pnm(1,n)*pnm(1,n)
      do k = 2 , kzp1
        plos(k,n) = plos(k-1,n) + &
             0.1_rkx*cplos*o3vmr(k-1,n)*(pnm(k,n)-pnm(k-1,n))
        plol(k,n) = plol(k-1,n) + 0.01_rkx*cplol*o3vmr(k-1,n) * &
                    (pnm(k,n)*pnm(k,n)-pnm(k-1,n)*pnm(k-1,n))
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine radoz2
  !
  !-----------------------------------------------------------------------
  !
  ! Compute temperatures and path lengths for longwave radiation
  !
  ! Input arguments
  !   ts     - Surface radiative temperature
  !   tnm    - Model level temperatures
  !   pnm    - Pressure at model interfaces (dynes/cm2)
  !   qnm    - Model level specific humidity
  !   pmln   - Ln(pmidm1)
  !   piln   - Ln(pintm1)
  !   plh2o  - Pressure weighted h2o path
  ! Output arguments
  !   tint   - Layer interface temperature
  !   tint4  - Tint to the 4th power
  !   tlayr  - K-1 level temperature
  !   tlayr4 - Tlayr to the 4th power
  !   tplnka - Level temperature from interface temperatures
  !   s2t    - H2o tmp and prs wghtd path lengt
  !   s2c    - H2o continuum path length
  !   w      - H2o path length
  !   tplnke - Equal to tplnka
  !
  !-----------------------------------------------------------------------
  !
  subroutine radtpl(n1,n2,ts,tnm,pnm,qnm,pmln,piln,plh2o, &
                    tint,tint4,tlayr,tlayr4,tplnka,s2t,s2c,wh2op,tplnke)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(n1:n2) , intent(in) :: ts
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: tnm , qnm , pmln
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: piln , pnm , plh2o
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: tint , tint4 , tplnka
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: tlayr , tlayr4
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: s2t , s2c , wh2op
    real(rkx) , dimension(n1:n2) , intent(out) :: tplnke
    !
    ! dy     - Thickness of layer for tmp interp
    ! dpnm   - Pressure thickness of layer
    ! dpnmsq - Prs squared difference across layer
    ! rtnm   - Inverse level temperature
    !
    integer(ik4) :: n
    integer(ik4) :: k
    real(rkx) :: dpnm , dpnmsq , dy , rtnm
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radtpl'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Set the top and bottom intermediate level temperatures,
    ! top level planck temperature and top layer temp**4.
    !
    ! Tint is lower interface temperature
    ! (not available for bottom layer, so use ground temperature)
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k,dpnm,dpnmsq,dy,rtnm)
#else
    do n = n1 , n2
#endif
      tint(kzp1,n) = ts(n)
      tint4(kzp1,n) = tint(kzp1,n)**4
      tplnka(1,n) = tnm(1,n)
      tint(1,n) = tplnka(1,n)
      tlayr4(1,n) = tplnka(1,n)**4
      tint4(1,n) = tlayr4(1,n)
      !
      ! Intermediate level temperatures are computed using temperature
      ! at the full level below less dy*delta t,between the full level
      !
      do k = 2 , kz
        dy = (piln(k,n)-pmln(k,n))/(pmln(k-1,n)-pmln(k,n))
        tint(k,n) = tnm(k,n) - dy*(tnm(k,n)-tnm(k-1,n))
        tint4(k,n) = tint(k,n)**4
      end do
      !
      ! Now set the layer temp=full level temperatures and establish a
      ! planck temperature for absorption (tplnka) which is the average
      ! the intermediate level temperatures.  Note that tplnka is not
      ! equal to the full level temperatures.
      !
      do k = 2 , kzp1
        tlayr(k,n) = tnm(k-1,n)
        tlayr4(k,n) = tlayr(k,n)**4
        tplnka(k,n) = (tint(k,n)+tint(k-1,n))*d_half
      end do
      !
      ! Calculate tplank for emissivity calculation.
      ! Assume isothermal tplnke i.e. all levels=ttop.
      !
      tplnke(n) = tplnka(1,n)
      tlayr(1,n) = tint(1,n)
      !
      ! Now compute h2o path fields:
      !
      s2t(1,n) = plh2o(1,n)*tnm(1,n)
      ! ccm3.2
      ! wh2op(1,n)   = (plh2o(1,n)*2.) / pnm(1,n)
      ! s2c(1,n) = plh2o(1,n) * qnm(1,n) * repsil
      ! ccm3.6.6
      wh2op(1,n) = sslp*(plh2o(1,n)*d_two)/pnm(1,n)
      rtnm = d_one/tnm(1,n)
      s2c(1,n) = plh2o(1,n)*exp(1800.0_rkx*(rtnm-r296))*qnm(1,n)*repsil
      do k = 1 , kz
        dpnm = pnm(k+1,n) - pnm(k,n)
        dpnmsq = pnm(k+1,n)**2 - pnm(k,n)**2
        rtnm = d_one/tnm(k,n)
        s2t(k+1,n) = s2t(k,n) + rgsslp*dpnmsq*qnm(k,n)*tnm(k,n)
        wh2op(k+1,n) = wh2op(k,n) + regravgts*qnm(k,n)*dpnm
        s2c(k+1,n) = s2c(k,n) + rgsslp*dpnmsq*qnm(k,n) * &
                     exp(1800.0_rkx*(rtnm-r296))*qnm(k,n)*repsil
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
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
  !-----------------------------------------------------------------------
  !
  ! RegCM - The Sun/Earth geometry is moved elsewhere
  !
  subroutine radinp(n1,n2,pmid,pint,h2ommr,cld,o3vmr, &
                    pbr,pnm,plco2,plh2o,tclrsf,o3mmr)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pint , cld
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: pmid , h2ommr , o3vmr
    real(rkx) , dimension(kz,n1:n2) , intent(out) :: pbr , o3mmr
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: plco2 , plh2o
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: pnm , tclrsf
    !
    ! vmmr - Ozone volume mixing ratio
    !
    real(rkx) , parameter :: vmmr = amo3/amd
    real(rkx) , parameter :: cpwpl = d_half*(amco2/amd)/(egravgts*sslp)
    integer(ik4) :: n , k
    !
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    ! pmid    - Pressure at model mid-levels (pascals)
    ! pint    - Pressure at model interfaces (pascals)
    ! h2ommr  - H2o mass mixing ratio
    ! cld     - Fractional cloud cover
    ! o3vmr   - ozone volume mixing ratio
    !
    ! Output arguments
    !
    ! pbr     - Pressure at interfaces (dynes/cm*2)
    ! pnm     - Pressure at mid-levels (dynes/cm*2)
    ! plco2   - Vert. pth lngth of co2 (prs-weighted)
    ! plh2o   - Vert. pth lngth h2o vap.(prs-weighted)
    ! o3mmr   - Ozone mass mixing ratio
    ! tclrsf  - Product of clr-sky fractions from top of atmosphere to level.
    !
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'radinp'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Convert pressure from pascals to dynes/cm2
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k)
#else
    do n = n1 , n2
#endif
      do k = 1 , kz
        pbr(k,n) = pmid(k,n)*d_10
        pnm(k,n) = pint(k,n)*d_10
      end do
      pnm(kzp1,n) = pint(kzp1,n)*d_10
      !
      ! Compute path quantities used in the longwave radiation:
      !
      plh2o(1,n) = rgsslp*h2ommr(1,n)*pnm(1,n)*pnm(1,n)
      plco2(1,n) = co2vmr(n)*cpwpl*pnm(1,n)*pnm(1,n)
      tclrsf(1,n) = d_one
      do k = 1 , kz
        plh2o(k+1,n) = plh2o(k,n) + rgsslp*(pnm(k+1,n)**2 - &
                       pnm(k,n)**2) * h2ommr(k,n)
        plco2(k+1,n) = co2vmr(n)*cpwpl*pnm(k+1,n)**2
        tclrsf(k+1,n) = tclrsf(k,n)*(d_one-cld(k+1,n))
      end do
      !
      ! Convert ozone volume mixing ratio to mass mixing ratio:
      !
      do k = 1 , kz
        o3mmr(k,n) = vmmr*o3vmr(k,n)
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine radinp

  ! xalpha - Term in direct reflect and transmissivity
  pure real(rkx) function xalpha(wi,uui,gi,ei)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , uui , gi , ei
    real(rk8) :: w , uu , g , e
    w = wi
    uu = uui
    g = gi
    e = ei
    xalpha = real(0.75_rk8*(w*uu)*((1.0_rk8+(g*(1.0_rk8-w))) / &
                  (1.0_rk8-((e*e)*(uu*uu)))),rkx)
  end function xalpha

  ! xgamma - Term in direct reflect and transmissivity
  pure real(rkx) function xgamma(wi,uui,gi,ei)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , uui , gi , ei
    real(rk8) :: w , uu , g , e
    w = wi
    uu = uui
    g = gi
    e = ei
    xgamma = real((w*0.5_rk8)*((3.0_rk8*g*(1.0_rk8-w)*(uu*uu)+1.0_rk8) / &
                               (1.0_rk8-((e*e)*(uu*uu)))),rkx)
  end function xgamma

  ! el - Term in xalpha,xgamma,f_n,f_u
  pure real(rkx) function el(wi,gi)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , gi
    real(rk8) :: w , g
    w = wi
    g = gi
    el = real(sqrt(3.0_rk8*(1.0_rk8-w)*(1.0_rk8-w*g)),rkx)
  end function el

  ! taus - Scaled extinction optical depth
  pure real(rkx) function taus(wi,fi,ti)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , fi , ti
    real(rk8) :: w , f , t
    w = wi
    f = fi
    t = ti
    taus = real((1.0_rk8-w*f)*t,rkx)
  end function taus

  ! omgs - Scaled single particle scattering albedo
  pure real(rkx) function omgs(wi,fi)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , fi
    real(rk8) :: w , f
    w = wi
    f = fi
    omgs = real((1.0_rk8-f)*w/(1.0_rk8-w*f),rkx)
  end function omgs

  ! asys - Scaled asymmetry parameter
  pure real(rkx) function asys(gi,fi)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: gi , fi
    real(rk8) :: g , f
    g = gi
    f = fi
    asys = real((g-f)/(1.0_rk8-f),rkx)
  end function asys

  ! f_u - Term in diffuse reflect and transmissivity
  pure real(rkx) function f_u(wi,gi,ei)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: wi , gi , ei
    real(rk8) :: w , g , e
    w = wi
    g = gi
    e = ei
    f_u = real(1.50_rk8*(1.0_rk8-w*g)/e,rkx)
  end function f_u

  ! f_n - Term in diffuse reflect and transmissivity
  pure real(rkx) function f_n(uui,eti)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: uui , eti
    real(rk8) :: uu , et
    uu = uui
    et = eti
    f_n = real(((uu+1.0_rk8)*(uu+1.0_rk8)/et) - &
               ((uu-1.0_rk8)*(uu-1.0_rk8)*et),rkx)
  end function f_n

  ! dbvt - Planck fnctn tmp derivative for o3
  pure real(rkx) function dbvt(ti)
!$acc routine seq
    ! Derivative of planck function at 9.6 micro-meter wavelength
    implicit none
    real(rkx) , intent(in) :: ti
    real(rk8) :: t
    t = ti
    dbvt = real((-2.8911366682e-4_rk8 + &
           (2.3771251896e-6_rk8+1.1305188929e-10_rk8*t)*t) /  &
           (1.0_rk8+(-6.1364820707e-3_rk8+1.5550319767e-5_rk8*t)*t),rkx)
  end function dbvt

  pure real(rkx) function fo3(uxi,vxi)
!$acc routine seq
    ! an absorption function factor
    implicit none
    real(rkx) , intent(in) :: uxi , vxi
    real(rk8) :: ux , vx
    ux = uxi
    vx = vxi
    fo3 = real(ux/sqrt(4.0_rk8+ux*(1.0_rk8+vx)),rkx)
  end function fo3

  pure integer(ik4) function intmax(imax)
!$acc routine seq
    implicit none
    integer(ik4) , pointer , dimension(:) , intent(in) :: imax
    integer(ik4) :: i , n , is , ie , mx
    is = lbound(imax,1)
    ie = ubound(imax,1)
    intmax = is
    n = ie-is+1
    if ( n > 1 ) then
      mx = imax(is)
      do i = is+1 , ie
        if ( imax(i) > mx ) then
          mx = imax(i)
          intmax = i
        end if
      end do
    end if
  end function intmax
  !
  !----------------------------------------------------------------------
  ! Calculate path lengths and pressure factors for CH4, N2O, CFC11
  ! and CFC12.
  !           Coded by J.T. Kiehl, November 21, 1994.
  !
  !-----------------------------------------------------------------------
  !
  !------------------------------Arguments--------------------------------
  !
  ! Input arguments
  !
  ! tnm    - Model level temperatures
  ! pnm    - Pressure at model interfaces (dynes/cm2)
  ! qmn    - h2o specific humidity
  ! cfc11  - CFC11 mass mixing ratio
  ! cfc12  - CFC12 mass mixing ratio
  ! n2o    - N2O mass mixing ratio
  ! ch4    - CH4 mass mixing ratio
  !
  ! Output arguments
  !
  ! ucfc11 - CFC11 path length
  ! ucfc12 - CFC12 path length
  ! un2o0  - N2O path length
  ! un2o1  - N2O path length (hot band)
  ! uch4   - CH4 path length
  ! uco211 - CO2 9.4 micron band path length
  ! uco212 - CO2 9.4 micron band path length
  ! uco213 - CO2 9.4 micron band path length
  ! uco221 - CO2 10.4 micron band path length
  ! uco222 - CO2 10.4 micron band path length
  ! uco223 - CO2 10.4 micron band path length
  ! bn2o0  - pressure factor for n2o
  ! bn2o1  - pressure factor for n2o
  ! bch4   - pressure factor for ch4
  ! uptype - p-type continuum path length
  !
  !-----------------------------------------------------------------------
  !
  subroutine trcpth(n1,n2,tnm,pnm,qnm,cfc11,cfc12,n2o,ch4,co2mmr, &
                    ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212, &
                    uco213,uco221,uco222,uco223,bn2o0,bn2o1,bch4, &
                    uptype)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(n1:n2) , intent(in) :: co2mmr
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: tnm , qnm
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: n2o , ch4
    real(rkx) , dimension(kz,n1:n2) , intent(in) :: cfc11 , cfc12
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: pnm
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: bch4 , uch4
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: bn2o0 , un2o0
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: bn2o1 , un2o1
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: ucfc11 , ucfc12
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: uco211 , uco212
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: uco213 , uco221
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: uco222 , uco223
    real(rkx) , dimension(kzp1,n1:n2) , intent(out) :: uptype
    !
    !   co2fac - co2 factor
    !   alpha1 - stimulated emission term
    !   alpha2 - stimulated emission term
    !   rt     - reciprocal of local temperature
    !   rsqrt  - reciprocal of sqrt of temp
    !   pbar   - mean pressure
    !   dpnm   - difference in pressure
    !
    real(rkx) , parameter :: diff = 1.66_rkx ! diffusivity factor
    real(rkx) :: alpha1 , alpha2 , dpnm , pbar , rsqrt , rt , co2fac
    integer(ik4) :: n , k

#ifdef STDPAR
    do concurrent ( n = n1:n2 ) &
      local(k,alpha1,alpha2,dpnm,pbar,rsqrt,rt,co2fac)
#else
    do n = n1 , n2
#endif
      !-----------------------------------------------------------------------
      !   Calculate path lengths for the trace gases
      !-----------------------------------------------------------------------
      ucfc11(1,n) = 1.8_rkx*cfc11(1,n)*pnm(1,n)*regravgts
      ucfc12(1,n) = 1.8_rkx*cfc12(1,n)*pnm(1,n)*regravgts
      un2o0(1,n) = diff*1.02346e5_rkx*n2o(1,n)*pnm(1,n)*regravgts/sqrt(tnm(1,n))
      un2o1(1,n) = diff*2.01909_rkx*un2o0(1,n)*exp(-847.36_rkx/tnm(1,n))
      uch4(1,n) = diff*8.60957e4_rkx*ch4(1,n)*pnm(1,n)*regravgts/sqrt(tnm(1,n))
      co2fac = diff*co2mmr(n)*pnm(1,n)*regravgts
      alpha1 = (d_one-exp(-1540.0_rkx/tnm(1,n)))**3/sqrt(tnm(1,n))
      alpha2 = (d_one-exp(-1360.0_rkx/tnm(1,n)))**3/sqrt(tnm(1,n))
      uco211(1,n) = 3.42217e3_rkx*co2fac*alpha1*exp(-1849.7_rkx/tnm(1,n))
      uco212(1,n) = 6.02454e3_rkx*co2fac*alpha1*exp(-2782.1_rkx/tnm(1,n))
      uco213(1,n) = 5.53143e3_rkx*co2fac*alpha1*exp(-3723.2_rkx/tnm(1,n))
      uco221(1,n) = 3.88984e3_rkx*co2fac*alpha2*exp(-1997.6_rkx/tnm(1,n))
      uco222(1,n) = 3.67108e3_rkx*co2fac*alpha2*exp(-3843.8_rkx/tnm(1,n))
      uco223(1,n) = 6.50642e3_rkx*co2fac*alpha2*exp(-2989.7_rkx/tnm(1,n))
      bn2o0(1,n) = diff*19.399_rkx*pnm(1,n)**2*n2o(1,n) * &
                 1.02346e5_rkx*regravgts/(sslp*tnm(1,n))
      bn2o1(1,n) = bn2o0(1,n)*exp(-847.36_rkx/tnm(1,n))*2.06646e5_rkx
      bch4(1,n) = diff*2.94449_rkx*ch4(1,n)*pnm(1,n)**2*regravgts * &
                8.60957e4_rkx/(sslp*tnm(1,n))
      uptype(1,n) = diff*qnm(1,n)*pnm(1,n)**2*exp(1800.0_rkx* &
                  (d_one/tnm(1,n)-r296))*regravgts/sslp
      do k = 1 , kz
        rt = d_one/tnm(k,n)
        rsqrt = sqrt(rt)
        pbar = ((pnm(k+1,n)+pnm(k,n))*d_half)/sslp
        dpnm = (pnm(k+1,n)-pnm(k,n))*regravgts
        alpha1 = diff*rsqrt*(d_one-exp(-1540.0_rkx/tnm(k,n)))**3
        alpha2 = diff*rsqrt*(d_one-exp(-1360.0_rkx/tnm(k,n)))**3
        ucfc11(k+1,n) = ucfc11(k,n) + 1.8_rkx*cfc11(k,n)*dpnm
        ucfc12(k+1,n) = ucfc12(k,n) + 1.8_rkx*cfc12(k,n)*dpnm
        un2o0(k+1,n) = un2o0(k,n) + diff*1.02346e5_rkx*n2o(k,n)*rsqrt*dpnm
        un2o1(k+1,n) = un2o1(k,n) + diff*2.06646e5_rkx*n2o(k,n) * &
                     rsqrt*exp(-847.36_rkx/tnm(k,n))*dpnm
        uch4(k+1,n) = uch4(k,n) + diff*8.60957e4_rkx*ch4(k,n)*rsqrt*dpnm
        uco211(k+1,n) = uco211(k,n) + 1.15_rkx*3.42217e3_rkx*alpha1 * &
                      co2mmr(n)*exp(-1849.7_rkx/tnm(k,n))*dpnm
        uco212(k+1,n) = uco212(k,n) + 1.15_rkx*6.02454e3_rkx*alpha1 * &
                      co2mmr(n)*exp(-2782.1_rkx/tnm(k,n))*dpnm
        uco213(k+1,n) = uco213(k,n) + 1.15_rkx*5.53143e3_rkx*alpha1 * &
                      co2mmr(n)*exp(-3723.2_rkx/tnm(k,n))*dpnm
        uco221(k+1,n) = uco221(k,n) + 1.15_rkx*3.88984e3_rkx*alpha2 * &
                      co2mmr(n)*exp(-1997.6_rkx/tnm(k,n))*dpnm
        uco222(k+1,n) = uco222(k,n) + 1.15_rkx*3.67108e3_rkx*alpha2 * &
                      co2mmr(n)*exp(-3843.8_rkx/tnm(k,n))*dpnm
        uco223(k+1,n) = uco223(k,n) + 1.15_rkx*6.50642e3_rkx*alpha2 * &
                      co2mmr(n)*exp(-2989.7_rkx/tnm(k,n))*dpnm
        bn2o0(k+1,n) = bn2o0(k,n) + diff*19.399_rkx*pbar*rt * &
                     1.02346e5_rkx*n2o(k,n)*dpnm
        bn2o1(k+1,n) = bn2o1(k,n) + diff*19.399_rkx*pbar*rt * &
                     2.06646e5_rkx*exp(-847.36_rkx/tnm(k,n))*n2o(k,n)*dpnm
        bch4(k+1,n) = bch4(k,n) + diff*2.94449_rkx*rt*pbar * &
                  8.60957e4_rkx*ch4(k,n)*dpnm
        uptype(k+1,n) = uptype(k,n) + diff*qnm(k,n)*exp(1800.0_rkx*(d_one / &
                      tnm(k,n)-r296))*pbar*dpnm
      end do
    end do
  end subroutine trcpth
  !
  !-----------------------------------------------------------------------
  !
  ! SUBROUTINE AERMIX
  !
  ! Set global mean tropospheric aerosol
  !
  ! Specify aerosol mixing ratio and compute relative humidity for later
  ! adjustment of aerosol optical properties. Aerosol mass mixing ratio
  ! is specified so that the column visible aerosol optical depth is a
  ! specified global number (tauvis). This means that the actual mixing
  ! ratio depends on pressure thickness of the lowest three atmospheric
  ! layers near the surface.
  !
  ! Optical properties and relative humidity parameterization are from:
  !
  ! J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
  ! and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
  ! 16 April 1993
  !
  ! Visible (vis) here means 0.5-0.7 micro-meters
  ! Forward scattering fraction is taken as asymmetry parameter squared
  !
  ! Input
  !   pnm    - Radiation level interface pressures (dynes/cm2)
  ! Output
  !   aermmb - Aerosol background mass mixing ratio
  !
  !-----------------------------------------------------------------------
  !
  ! RegCM : Removed the dependency on relative humidity
  !
  subroutine aermix(n1,n2,pnm,aermmb)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , intent(in) , dimension(kzp1,n1:n2) :: pnm
    real(rkx) , intent(out) , dimension(kz,n1:n2) :: aermmb
    !
    !-----------------------------------------------------------------------
    !
    ! mxaerl - max nmbr aerosol levels counting up from surface
    ! tauvis - visible optical depth
    ! kaervs - visible extinction coefficiant of aerosol (m2/g)
    ! omgvis - visible omega0
    ! gvis   - visible forward scattering asymmetry parameter
    !
    !-----------------------------------------------------------------------
    !
    integer(ik4) , parameter :: mxaerl = 4
    ! multiplication factor for kaer
    real(rkx) , parameter :: kaervs = 5.3012_rkx
    real(rkx) , parameter :: omgvis = 0.999999_rkx
    real(rkx) , parameter :: gvis = 0.694889_rkx
    ! added for efficiency
    real(rkx) , parameter :: rhfac = 1.6718_rkx
    !
    ! real(rkx) , parameter :: a0 = -9.2906106183_rkx
    ! real(rkx) , parameter :: a1 =  0.52570211505_rkx
    ! real(rkx) , parameter :: a2 = -0.0089285760691_rkx
    ! real(rkx) , parameter :: a4 =  5.0877212432e-05_rkx
    !
    integer(ik4) :: n , k
    !fil  tauvis = 0.01_rkx
    real(rkx) , parameter :: tauvis = 0.14_rkx
    !
    !-----------------------------------------------------------------------
    !
    ! Define background aerosol
    ! Set relative humidity and factor; then aerosol amount:
    !
    ! if ( rh(i,k) > 0.9 ) then
    !   rhfac = 2.8
    ! else if (rh(i,k) < 0.6 ) then
    !   rhfac = 1.0
    ! else
    !   rhpc  = 100.0 * rh(i,k)
    !   rhfac = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
    ! end if
    !
    ! Find constant aerosol mass mixing ratio for specified levels
    ! in the column, converting units where appropriate
    ! for the moment no more used
    !
    do n = n1 , n2
      do k = 1 , kz
        if ( k >= kz + 1 - mxaerl ) then
          aermmb(k,n) = egravgts * tauvis / &
                  (1.0e4_rkx*kaervs*rhfac*(d_one-omgvis*gvis*gvis) * &
                  (pnm(kzp1,n)-pnm(kzp1-mxaerl,n)))
        else
          aermmb(k,n) = d_zero
        end if
      end do
    end do
  end subroutine aermix
  !
  !----------------------------------------------------------------------
  !   Calculate Planck factors for absorptivity and emissivity of
  !   CH4, N2O, CFC11 and CFC12
  !
  !-----------------------------------------------------------------------
  !
  ! Input arguments
  !
  ! tint    - interface temperatures
  ! tlayr   - k-1 level temperatures
  ! tplnke  - Top Layer temperature
  !
  ! output arguments
  !
  ! emplnk  - emissivity Planck factor
  ! abplnk1 - non-nearest layer Plack factor
  ! abplnk2 - nearest layer factor
  !
  subroutine trcplk(n1,n2,tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(n1:n2) , intent(in) :: tplnke
    real(rkx) , dimension(kzp1,n1:n2) , intent(in) :: tint , tlayr
    real(rkx) , dimension(nlwspi,n1:n2) , intent(out) :: emplnk
    real(rkx) , dimension(nlwspi,kzp1,n1:n2) , intent(out) :: abplnk1 , abplnk2
    !
    ! wvl   - wavelength index
    ! f1    - Planck function factor
    ! f2    -       "
    ! f3    -       "
    !
    integer(ik4) :: n , k , wvl
    !
    ! Calculate emissivity Planck factor
    !
#ifdef STDPAR
    do concurrent ( n = n1:n2 ) local(k,wvl)
#else
    do n = n1 , n2
#endif
      do wvl = 1 , nlwspi
        emplnk(wvl,n) = f1(wvl)/(tplnke(n)**4*(exp(f3(wvl)/tplnke(n))-d_one))
      end do
      !
      ! Calculate absorptivity Planck factor for tint and tlayr temperatures
      !
      do  k = 1 , kzp1
        do wvl = 1 , nlwspi
          ! non-nearlest layer function
          abplnk1(wvl,k,n) = (f2(wvl)*exp(f3(wvl)/tint(k,n))) / &
                           (tint(k,n)**5*(exp(f3(wvl)/tint(k,n))-d_one)**2)
          ! nearest layer function
          abplnk2(wvl,k,n) = (f2(wvl)*exp(f3(wvl)/tlayr(k,n))) / &
                           (tlayr(k,n)**5*(exp(f3(wvl)/tlayr(k,n))-d_one)**2)
        end do
      end do
    end do
  end subroutine trcplk

end module mod_rad_radiation
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
