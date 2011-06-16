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
 
module mod_radiation

  use mod_runparams
  use mod_bats
  use mod_aerosol
  use mod_message
  use mod_service
  use mod_memutil
#ifdef CLM
  use mod_clm
#endif

! Used by this module only

  use mod_tracer
  use mod_scenarios

  private

! Maximum total cloud fraction for radiation model
  real(8) :: cftotmax

  public :: allocate_mod_radiation , radini , radctl
  public :: absnxt , abstot , emstot
  public :: cftotmax

! absnxt  - Nearest layer absorptivities
! abstot  - Non-adjacent layer absorptivites
! emstot  - Total emissivity

  real(8) , pointer , dimension(:,:,:,:)  :: absnxt , absnxt0
  real(8) , pointer , dimension(:,:,:,:)  :: abstot , abstot0
  real(8) , pointer , dimension(:,:,:) :: emstot , emstot0
  real(8) , pointer , dimension(:,:,:,:):: xuinpl
!
  real(8) , dimension(2) :: a1 , a2 , b1 , b2 , realk , st
  real(8) , dimension(4) :: c1 , c2 , c3 , c4 , c5 , c6 , c7
  real(8) :: c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 ,  &
             c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 ,  &
             c28 , c29 , c30 , c31 , c8 , c9 , cfa1
  real(8) :: co2vmr
  real(8) , dimension(3,4) :: coefa , coefc , coefe
  real(8) , dimension(4,4) :: coefb , coefd
  real(8) , dimension(6,2) :: coeff , coefi
  real(8) , dimension(2,4) :: coefg , coefh
  real(8) , dimension(3,2) :: coefj , coefk
!
  real(8) , parameter :: verynearone = 0.999999D0
!
! r80257   - Conversion factor for h2o pathlength
  real(8) , parameter :: r80257 = d_one/8.0257D-04
  real(8) , parameter :: r293 = d_one/293.0D0
  real(8) , parameter :: r250 = d_one/250.0D0
! r3205    - Line width factor for o3 (see R&Di)
  real(8) , parameter :: r3205 = d_one/0.3205D0
  real(8) , parameter :: r300 = d_one/300.0D0
! r2sslp   - 1/2 of rsslp
  real(8) , parameter :: r2sslp = d_one/(d_two*sslp)
! r296   - Inverse stand temp for h2o continuum
  real(8) , parameter :: r296 = d_one/296.0D0
! repsil - Inver ratio mol weight h2o to dry air
  real(8) , parameter :: repsil = d_one/ep2
!
! Initialize further longwave constants referring to far wing
! correction; R&D refers to:
!
! Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
! Emissivity and Absorptivity Formulation for Water Vapor
! Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
  real(8) , parameter :: fwcoef = 0.1D0      ! See eq(33) R&D
  real(8) , parameter :: fwc1 = 0.30D0       ! See eq(33) R&D
  real(8) , parameter :: fwc2 = 4.5D0        ! See eq(33) and eq(34) in R&D
  real(8) , parameter :: fc1 = 2.6D0         ! See eq(34) R&D
!
! Initialize ozone data.
!
  real(8) , parameter :: v0 = 22.4136D0  ! Volume of a gas at stp (m**3/kmol)
  real(8) , parameter :: p0 = 0.1D0*sslp ! Standard pressure (pascals)
!
! Constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
  real(8) , parameter :: cplos = v0/(amd*egrav)*d_100
  real(8) , parameter :: cplol = v0/(amd*egrav*p0)*d_half*d_100
!
! v_raytau_xx - Constants for new bands
! v_abo3_xx   - Constants for new bands
!
  real(8) , parameter :: v_raytau_35 = 0.155208D0
  real(8) , parameter :: v_raytau_64 = 0.0392D0
  real(8) , parameter :: v_abo3_35   = 2.4058030D+01
  real(8) , parameter :: v_abo3_64   = 2.210D+01
!
! delta    - Pressure (atmospheres) for stratos. h2o limit
! o2mmr    - O2 mass mixing ratio
!
  real(8) , parameter :: delta = 1.70D-3
  real(8) , parameter :: o2mmr = 0.23143D0
!
! Minimum total transmission below which no layer computation are done:
!
! trmin   - Minimum total transmission allowed
! wray    - Rayleigh single scatter albedo
! gray    - Rayleigh asymetry parameter
! fray    - Rayleigh forward scattered fraction
!
  real(8) , parameter :: trmin = 1.0D-3
  real(8) , parameter :: wray = 0.999999D0
  real(8) , parameter :: gray = 0.0D0
  real(8) , parameter :: fray = 0.1D0
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
  real(8) , dimension(4) :: abari , abarl , bbari , bbarl , cbari , &
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
  real(8) , dimension(nspi) :: abco2 , abh2o , abo2 , abo3 ,   &
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

  subroutine allocate_mod_radiation 
    implicit none        
    character (len=50) :: subroutine_name='allocate_mod_radiation'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)

    call getmem4d(absnxt,1,iym1,1,kz,1,4,1,jxp,'radiation:absnxt')
    call getmem4d(absnxt0,1,iym1,1,kz,1,4,1,jxp,'radiation:absnxt0')
    call getmem4d(abstot,1,iym1,1,kzp1,1,kzp1,1,jxp,'radiation:abstot')
    call getmem4d(abstot0,1,iym1,1,kzp1,1,kzp1,1,jxp,'radiation:abstot0')
    call getmem3d(emstot,1,iym1,1,kzp1,1,jxp,'radiation:emstot')
    call getmem3d(emstot0,1,iym1,1,kzp1,1,jxp,'radiation:emstot0')
    call getmem4d(xuinpl,1,iym1,1,kzp1,1,4,1,jxp,'radiation:xuinpl')

    call time_end(subroutine_name,indx)

  end subroutine allocate_mod_radiation 
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
  subroutine radini
!
    implicit none
!
!   iband  - H2O band index
!
    integer :: iband

    character (len=50) :: subroutine_name='radini'
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
    if ( idatex%year >= 1750 .and. idatex%year <= 2100 ) then
      co2vmr = cgas(2,idatex%year)*1.0D-6
      co2mmr = co2vmr*44.0D0/28.9644D0
      ch40 = cgas(3,idatex%year)*1.0D-9*0.55241D0
      n2o0 = cgas(4,idatex%year)*1.0D-9*1.51913D0
      cfc110 = cgas(5,idatex%year)*1.0D-12*4.69548D0
      cfc120 = cgas(6,idatex%year)*1.0D-12*4.14307D0
    else
      write (aline,*) '  Simulation year:  ' , idatex%year
      call say
      call fatal(__FILE__,__LINE__,                                   &
            'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)')
    end if
!   print*,'IN RADINI (TOP)'
!   print*,'  co2vmr= ',co2vmr
!   print*,'  co2mmr= ',co2mmr
!   print*,'  ch40  = ',ch40
!   print*,'  n2o0  = ',n2o0
!   print*,'  cfc110= ',cfc110
!   print*,'  cfc120= ',cfc120
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
  subroutine radctl(jslc,alat,coslat,ts,pmid,pint,pmln,piln,t,       &
                    h2ommr,cld,effcld,clwp,albs,albsd,albl,albld,    &
                    fsns,qrs,qrl,flwds,rel,rei,fice,sols,soll,solsd, &
                    solld,emiss,fsnt,fsntc,fsnsc,flnt,flns,flntc,    &
                    flnsc,solin,alb,albc,fsds,fsnirt,fsnrtc,         &
                    fsnirtsq,eccf,o3vmr)
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
!   alat    - current latitude(radians)
!   coslat  - cosine latitude
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
    real(8) :: eccf
    integer :: jslc
    real(8) , dimension(iym1) :: alb , albc , albl , albld , albs ,  &
                                 albsd , alat , coslat , emiss ,     &
                                 flns , flnsc , flnt , flntc ,       &
                                 flwds , fsds , fsnirt , fsnirtsq ,  &
                                 fsnrtc , fsns , fsnsc , fsnt ,      &
                                 fsntc , solin , soll , solld ,      &
                                 sols , solsd , ts
    real(8) , dimension(iym1,kzp1) :: cld , effcld , piln , pint
    real(8) , dimension(iym1,kz) :: clwp , fice , h2ommr , pmid ,  &
           pmln , qrl , qrs , rei , rel , t
    real(8) , dimension(iym1,kz) :: o3vmr
    intent (out) alb , albc
    intent (inout) flns , flnsc , flnt , flntc , flwds , fsds ,       &
                   fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   fsntc , solin
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
!   pbr      - Model mid-level pressures (dynes/cm2)
!   pnm      - Model interface pressures (dynes/cm2)
!   o3vmr    - Ozone volume mixing ratio
!   o3mmr    - Ozone mass mixing ratio
!   plco2    - Prs weighted CO2 path
!   plh2o    - Prs weighted H2O path
!   tclrsf   - Total clear sky fraction, level to space
!   eccf     - Earth/sun distance factor
!   n2o      - nitrous oxide mass mixing ratio
!   ch4      - methane mass mixing ratio
!   cfc11    - cfc11 mass mixing ratio
!   cfc12    - cfc12 mass mixing ratio
!   rh       - level relative humidity (fraction)
!
    real(8) , dimension(iym1) :: fslwdcs
    real(8) , dimension(iym1) :: aeradfo , aeradfos
    real(8),  dimension(iym1)::  aerlwfo , aerlwfos
    real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o
    integer :: i , k
    real(8) , dimension(iym1,kz) :: o3mmr , pbr , rh
    real(8) , dimension(iym1,kzp1) :: plco2 , plh2o , pnm , tclrsf
    character (len=50) :: subroutine_name='radctl'
    integer :: indx = 0
    real(8) , dimension(iym1) :: totcf
!
    call time_begin(subroutine_name,indx)
!
!   Instead of interpolating the o3vmr from the time-interpolated
!   values, we pass compute o3vmr in getdat() and pass it directly
!   into radctl(). o3mmr will be computed in radinp().
!
!   Set latitude dependent radiation input
!
    call radinp(pmid,pint,h2ommr,cld,o3vmr,pbr,pnm,plco2,plh2o,tclrsf,eccf,o3mmr)
!
!   Solar radiation computation
!
    if ( dosw ) then
!
!     Specify aerosol mass mixing ratio
!
      call aermix(pnm,rh,jslc,1,iym1,iym1,kz,ntr)
 
      call aeroppt(rh,pint)
!
      call radcsw(pnm,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,albs,   &
                  albsd,albl,albld,solin,qrs,fsns,fsnt,fsds,fsnsc,    &
                  fsntc,sols,soll,solsd,solld,fsnirt,fsnrtc,fsnirtsq, &
                  aeradfo,aeradfos)
!
!     call aerout(jslc,aeradfo,aeradfos)
!
!     Convert units of shortwave fields needed by rest of model from CGS to MKS
!
      do i = 1 , iym1
        solin(i) = solin(i)*1.0D-3
        fsnt(i) = fsnt(i)*1.0D-3
        fsns(i) = fsns(i)*1.0D-3
        fsntc(i) = fsntc(i)*1.0D-3
        fsnsc(i) = fsnsc(i)*1.0D-3
!
!       clear sky column partitioning for surface flux  
!       note : should be generalised to the whole column to be
!              really in energy balance !
!
        totcf(i) = d_one
        do k = 1 , kzp1
          totcf(i) = totcf(i) * (d_one - cld(i,k)) 
        end do
        totcf(i) = d_one - totcf(i)

!       maximum cld cover considered        
!       fsns(i) = fsns(i) * maxval(cld(i,:)) + &
!                 fsnsc(i) * (1-maxval(cld(i,:)))
!       random overlap assumption is tocf(i)
!       Now average btw rand ov and maximum cloud cover as fil suggest
        totcf(i) =  d_half * ( totcf(i) + maxval(cld(i,:)) )

!       Fil suggestion of putting a max on column cloud fraction
!       TAO: implement a user-specified CF maximum (default of 0.75d0,
!       as suggested by Erika)
        if ( totcf(i) > cftotmax ) totcf(i) = cftotmax
        if ( totcf(i) < d_zero ) totcf(i) = d_zero

        fsns(i) = fsns(i) * totcf(i) + fsnsc(i) * (d_one-totcf(i))
        fsds(i) = fsds(i)*1.0D-3
        fsnirt(i) = fsnirt(i)*1.0D-3
        fsnrtc(i) = fsnrtc(i)*1.0D-3
        fsnirtsq(i) = fsnirtsq(i)*1.0D-3
      end do
!
!     Calculate/outfld albedo and clear sky albedo
!
      do i = 1 , iym1
        if ( solin(i) > d_zero ) then
          alb(i) = (solin(i)-fsnt(i))/solin(i)
        else
          alb(i) = d_zero
        end if
      end do
!
      do i = 1 , iym1
        if ( solin(i) > d_zero ) then
          albc(i) = (solin(i)-fsntc(i))/solin(i)
        else
          albc(i) = d_zero
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
      call trcmix(pmid,alat,coslat,n2o,ch4,cfc11,cfc12)
!
      call radclw(jslc,ts,t,h2ommr,o3vmr,pbr,pnm,pmln,piln,plco2,     &
                  plh2o,n2o,ch4,cfc11,cfc12,effcld,tclrsf,qrl,flns,   &
                  flnt,flnsc,flntc,flwds,fslwdcs,emiss,aerlwfo,       &
                  aerlwfos)
!
!     Convert units of longwave fields needed by rest of model from CGS to MKS
!
      do i = 1 , iym1
        flnt(i) = flnt(i)*1.0D-3
        flns(i) = flns(i)*1.0D-3
        flntc(i) = flntc(i)*1.0D-3
        flnsc(i) = flnsc(i)*1.0D-3
        flwds(i) = flwds(i)*1.0D-3
!FAB
        fslwdcs(i) = fslwdcs(i)*1.0D-3
!       essai clear sky column
!
!       flwds(i) = flwds(i) * maxval(cld((i,:))) + &
!                  flwds(i) * (1-maxval(cld((i,:))))
!       flwds(i) = flwds(i) * maxval(cld(i,:)) + &
!                  fslwdcs(i)*(d_one-maxval(cld(i,:)))
!       flns(i) = flns(i) * maxval(cld(i,:)) + &
!                 flnsc(i)*(d_one-maxval(cld(i,:))) 
!
!       totcf(i) has been calculated for the SW, dolw is always true 
        flwds(i) = flwds(i) * totcf(i) + &
                   fslwdcs(i) * (d_one - totcf(i))
        flns(i) = flns(i) * totcf(i) + &
                  flnsc(i) * (d_one - totcf(i))
      end do
    end if

    if ( ichem==1 ) then
      call aerout(jslc,aeradfo,aeradfos,aerlwfo,aerlwfos)
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
! (i.e. coszrs > 0) computations are done.
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
  subroutine radcsw(pint,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,   &
                    asdir,asdif,aldir,aldif,solin,qrs,fsns,fsnt,    &
                    fsds,fsnsc,fsntc,sols,soll,solsd,solld,fsnirt,  &
                    fsnrtc,fsnirtsq,aeradfo,aeradfos)
! 
    implicit none
!
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
! asdir   - 0.2-0.7 micro-meter srfc alb to direct rad
! aldir   - 0.7-5.0 micro-meter srfc alb to direct rad
! asdif   - 0.2-0.7 micro-meter srfc alb to diffuse mod_rad
! aldif   - 0.7-5.0 micro-meter srfc alb to diffuse mod_rad
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
!
    real(8) :: eccf
    real(8) , dimension(iym1) :: aeradfo , aeradfos , aldif , aldir ,&
                                  asdif , asdir , fsds , fsnirt ,     &
                                  fsnirtsq , fsnrtc , fsns , fsnsc ,  &
                                  fsnt , fsntc , solin , soll ,       &
                                  solld , sols , solsd
    real(8) , dimension(iym1,kzp1) :: cld , pint
    real(8) , dimension(iym1,kz) :: clwp , fice , h2ommr , o3mmr , &
           qrs , rei , rel
    intent (in) aldif , aldir , asdif , asdir , cld , clwp , eccf ,   &
                fice , h2ommr , o3mmr , pint , rei , rel
    intent (out) aeradfo , aeradfos , fsds , qrs
    intent (inout) fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   fsntc , solin , soll , solld , sols , solsd
!
!---------------------------Local variables-----------------------------
!
! ns       - Spectral loop index
! i        - Longitude loop index
! k        - Level loop index
! n        - Loop index for daylight
! nloop    - Number of daylight loops
! is       - Daytime start indices
! ie       - Daytime end indices
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
!     Diagnostic and accumulation arrays; note that sfltot, fswup, and
!     fswdn are not used in the computation,but are retained for future
!     use.
!
! solflx   - Solar flux in current interval
! sfltot   - Spectrally summed total solar flux
! totfld   - Spectrally summed flux divergence
! fswup    - Spectrally summed up flux
! fswdn    - Spectrally summed down flux
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
!     Various arrays and other constants:
!
! pflx     - Interface press, including extra layer
! zenfac   - Square root of cos solar zenith angle
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
!
!     Layer absorber amounts; note that 0 refers to the extra layer
!     added above the top model layer
!
! uh2o     - Layer absorber amount of h2o
! uo3      - Layer absorber amount of  o3
! uco2     - Layer absorber amount of co2
! uo2      - Layer absorber amount of  o2
!
!     Total column absorber amounts:
!
! uth2o    - Total column  absorber amount of h2o
! uto3     - Total column  absorber amount of  o3
! utco2    - Total column  absorber amount of co2
! uto2     - Total column  absorber amount of  o2
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
!
! wkaer    - works table
! aeradfo  - spectrally integrated aerosol radiative forcing ( TOA)
!-----------------------------------------------------------------------
!
    real(8) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli ,  &
               dbarii , dbarli , ebarii , ebarli , fbarii ,   &
               fbarli , h2ostr , path , pdel , psf ,          &
               pthco2 , pthh2o , ptho2 , ptho3 , xptop , rdenom ,     &
               sqrco2 , tmp1 , tmp1i , tmp1l , tmp2 , tmp2i , tmp2l , &
               tmp3i , tmp3l , trayoslp , wavmid , wgtint
    real(8) , dimension(iym1,0:kz) :: explay , fci , fcl , flxdiv ,&
           gci , gcl , rdif , rdir , tauxci , tauxcl , tdif , tdir ,  &
           totfld , uco2 , uh2o , uo2 , uo3 , wci , wcl
    real(8) , dimension(iym1,0:kzp1) :: exptdn , fluxdn , fluxup ,  &
           fswdn , fswup , pflx , rdndif , rupdif , rupdir , tottrn
    integer :: i , indxsl , k , n , nloop , ns
    integer , dimension(2) :: ie , is
    real(8) , dimension(iym1) :: sfltot , solflx , utco2 , uth2o ,   &
                                  uto2 , uto3 , x0fsnrtc , x0fsnsc ,  &
                                  x0fsntc , zenfac
    real(8) , dimension(iym1,0:kz,4) :: wkaer
    real(8) , dimension(iym1,4) :: zero
!
    character (len=50) :: subroutine_name='radcsw'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
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
    sabveg(:) = d_zero
    solis(:) = d_zero
    solvs(:) = d_zero
    solvd(:) = d_zero
!
    aeradfo(:) = d_zero
    aeradfos(:) = d_zero
    x0fsntc(:) = d_zero
    x0fsnsc(:) = d_zero
    x0fsnrtc(:) = d_zero
!
    qrs(:,:) = d_zero
!
    do k = 1 , kz
      do i = 1 , iym1
        pdel = pint(i,k+1) - pint(i,k)
        path = pdel*regravgts
      end do
    end do
    nloop = 0
    is = 0
    ie = 0
!
!   Compute starting daytime loop index
!
    is(1) = isrchfgt(iym1,coszrs,1,d_zero)
!
!   If night everywhere, return
!
    if ( is(1) > iym1 ) return
!
!   Compute ending daytime loop index
!
    ie(1) = isrchfle(iym1-is(1),coszrs(is(1)+1:iy),1,d_zero) + is(1)-1
    nloop = 1
!
!   Possibly 2 daytime loops needed
!
    if ( ie(1) /= iym1 ) then
      is(2) = isrchfgt(iym1-ie(1),coszrs(ie(1)+1:iy),1,d_zero) + ie(1)
      if ( is(2) < iym1 ) then
        ie(2) = isrchfle(iym1-is(2),coszrs(is(2)+1:iy),1,d_zero) + is(2)-1
        if ( ie(2) > is(2) ) then
          nloop = 2
        end if
      end if
    end if
!
!   Define solar incident radiation and interface pressures:
!
    do n = 1 , nloop
      do i = is(n) , ie(n)
        solin(i) = scon*eccf*coszrs(i)
        pflx(i,0) = d_zero
      end do
    end do
    do k = 1 , kzp1
      do n = 1 , nloop
        do i = is(n) , ie(n)
          pflx(i,k) = pint(i,k)
        end do
      end do
    end do
!
!   Compute optical paths:
!   CO2, use old scheme(as constant)
!
    tmp1 = d_half/(egravgts*sslp)
!   co2mmr = co2vmr*(mmwco2/mmwair)
   
    sqrco2 = dsqrt(co2mmr)
    do n = 1 , nloop
      do i = is(n) , ie(n)
        xptop = pflx(i,1)
        ptho2 = o2mmr*xptop*regravgts
        ptho3 = o3mmr(i,1)*xptop*regravgts
        pthco2 = sqrco2*(xptop*regravgts)
        h2ostr = dsqrt(d_one/h2ommr(i,1))
        zenfac(i) = dsqrt(coszrs(i))
        pthh2o = xptop**d_two*tmp1+(xptop*regravgts)*(h2ostr*zenfac(i)*delta)
        uh2o(i,0) = h2ommr(i,1)*pthh2o
        uco2(i,0) = zenfac(i)*pthco2
        uo2(i,0) = zenfac(i)*ptho2
        uo3(i,0) = ptho3
      end do
    end do
!
    tmp2 = delta*regravgts
    do k = 1 , kz
      do n = 1 , nloop
        do i = is(n) , ie(n)
          pdel = pflx(i,k+1) - pflx(i,k)
          path = pdel*regravgts
          ptho2 = o2mmr*path
          ptho3 = o3mmr(i,k)*path
          pthco2 = sqrco2*path
          h2ostr = dsqrt(d_one/h2ommr(i,k))
          pthh2o = (pflx(i,k+1)**d_two-pflx(i,k)**d_two) * &
                    tmp1 + pdel*h2ostr*zenfac(i)*tmp2
          uh2o(i,k) = h2ommr(i,k)*pthh2o
          uco2(i,k) = zenfac(i)*pthco2
          uo2(i,k) = zenfac(i)*ptho2
          uo3(i,k) = ptho3
        end do
      end do
    end do
!
!   Compute column absorber amounts for the clear sky computation:
!
    do n = 1 , nloop
      do i = is(n) , ie(n)
        uth2o(i) = d_zero
        uto3(i) = d_zero
        utco2(i) = d_zero
        uto2(i) = d_zero
      end do
    end do
    do k = 1 , kz
      do n = 1 , nloop
        do i = is(n) , ie(n)
          uth2o(i) = uth2o(i) + uh2o(i,k)
          uto3(i) = uto3(i) + uo3(i,k)
          utco2(i) = utco2(i) + uco2(i,k)
          uto2(i) = uto2(i) + uo2(i,k)
        end do
      end do
    end do
!
!   Initialize spectrally integrated totals:
!
    do k = 0 , kz
      do i = 1 , iym1
        totfld(i,k) = d_zero
        fswup(i,k) = d_zero
        fswdn(i,k) = d_zero
      end do
    end do
    do i = 1 , iym1
      sfltot(i) = d_zero
      fswup(i,kzp1) = d_zero
      fswdn(i,kzp1) = d_zero
    end do
!
!   Set cloud properties for top (0) layer; so long as tauxcl is zero,
!   there is no cloud above top of model; the other cloud properties
!   are arbitrary:
!
    do n = 1 , nloop
      do i = is(n) , ie(n)
        tauxcl(i,0) = d_zero
        wcl(i,0) = verynearone
        gcl(i,0) = 0.850D0
        fcl(i,0) = 0.725D0
        tauxci(i,0) = d_zero
        wci(i,0) = verynearone
        gci(i,0) = 0.850D0
        fci(i,0) = 0.725D0
      end do
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
        do n = 1 , nloop
          do i = is(n) , ie(n)
!
!           liquid
!
            tmp1l = abarli + bbarli/rel(i,k)
            tmp2l = d_one - cbarli - dbarli*rel(i,k)
            tmp3l = fbarli*rel(i,k)
!
!           ice
!
            tmp1i = abarii + bbarii/rei(i,k)
            tmp2i = d_one - cbarii - dbarii*rei(i,k)
            tmp3i = fbarii*rei(i,k)
!
!           Cloud fraction incorporated into cloud extinction optical depth
!found
!           April 12 2000, Filippo found the different scheme here:
   
!scheme     1
!ccm3.6.6
!           tauxcl(i,k) = clwp(i,k)*tmp1l*(d_one-fice(i,k))*cld(i,k)*dsqrt(cld(i,k))
!           tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)*cld(i,k)*dsqrt(cld(i,k))
!
!scheme     2
!KN
            tauxcl(i,k) = clwp(i,k)*tmp1l*(d_one-fice(i,k))*cld(i,k) / &
                          (d_one+(d_one-0.85D0)*(d_one-cld(i,k))*      &
                          clwp(i,k)*tmp1l*(d_one-fice(i,k)))
            tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)*cld(i,k) /     &
                          (d_one+(d_one-0.78D0)*(d_one-cld(i,k)) * &
                          clwp(i,k)*tmp1i*fice(i,k))
   
!scheme     3
!EES        below replaced
!           tauxcl(i,k) = clwp(i,k)*tmp1l*(d_one-fice(i,k))*cld(i,k)**0.85
!           tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)*cld(i,k)**0.85
!found_
!
!           Do not let single scatter albedo be 1; delta-eddington
!           solution for non-conservative case:
!
!qian       30/06/99        wcl(i,k) = dmin1(tmp2l,0.999999)
            wcl(i,k) = dmin1(tmp2l,verynearone)
            gcl(i,k) = ebarli + tmp3l
            fcl(i,k) = gcl(i,k)*gcl(i,k)
!
            wci(i,k) = dmin1(tmp2i,verynearone)
            gci(i,k) = ebarii + tmp3i
            fci(i,k) = gci(i,k)*gci(i,k)
!
          end do
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
        do n = 1 , nloop
          do i = is(n) , ie(n)
            albdir(i) = asdir(i)
            albdif(i) = asdif(i)
          end do
        end do
!
!       Wavelength greater than 0.7 micro-meter
!
      else
        do n = 1 , nloop
          do i = is(n) , ie(n)
            albdir(i) = aldir(i)
            albdif(i) = aldif(i)
          end do
        end do
      end if
      trayoslp = raytau(ns)/sslp
!
!     Layer input properties now completely specified; compute the
!     delta-Eddington solution reflectivities and transmissivities
!     for each layer, starting from the top and working downwards:
   
!     options for aerosol: no climatic feedback if idirect == 1
!     should be consistent with aeroppt routine

      if (ichem==1 ) then
        if ( idirect == 2 ) then
          do k = 0 , kz
            do i = 1 , iym1
              wkaer(i,k,1) = tauxar_mix(i,k,ns)
              wkaer(i,k,2) = tauasc_mix(i,k,ns)
              wkaer(i,k,3) = gtota_mix(i,k,ns)
              wkaer(i,k,4) = ftota_mix(i,k,ns)
            end do
          end do
        else if ( idirect == 1 ) then
          do k = 0 , kz
            do i = 1 , iym1
              wkaer(i,k,1) = d_zero
              wkaer(i,k,2) = d_zero
              wkaer(i,k,3) = d_zero
              wkaer(i,k,4) = d_zero
            end do
          end do
        end if
      else
        do k = 0 , kz
          do i = 1 , iym1
            wkaer(i,k,1) = d_zero
            wkaer(i,k,2) = d_zero
            wkaer(i,k,3) = d_zero
            wkaer(i,k,4) = d_zero
          end do
        end do
      end if
   
      call radded(coszrs,trayoslp,pflx,ns,uh2o,uo3,uco2,uo2,tauxcl, &
                  wcl,gcl,fcl,tauxci,wci,gci,fci,wkaer(1,0,1),      &
                  wkaer(1,0,2),wkaer(1,0,3),wkaer(1,0,4),nloop,     &
                  is,ie,rdir,rdif,tdir,tdif,explay,exptdn,rdndif,tottrn)
!
!     Compute reflectivity to direct and diffuse mod_radiation for layers
!     below by adding succesive layers starting from the surface and
!     working upwards:
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          rupdir(i,kzp1) = albdir(i)
          rupdif(i,kzp1) = albdif(i)
        end do
      end do
      do k = kz , 0 , -1
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rdenom = d_one/(d_one-rdif(i,k)*rupdif(i,k+1))
            rupdir(i,k) = rdir(i,k) + tdif(i,k) *                    &
                          (rupdir(i,k+1)*explay(i,k)+rupdif(i,k+1) * &
                          (tdir(i,k)-explay(i,k)))*rdenom
            rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)**d_two*rdenom
          end do
        end do
      end do
!
!     Compute up and down fluxes for each interface, using the added
!     atmospheric layer properties at each interface:
!
      do k = 0 , kzp1
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rdenom = d_one/(d_one-rdndif(i,k)*rupdif(i,k))
            fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+                   &
                          (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
            fluxdn(i,k) = exptdn(i,k) +                          &
                          (tottrn(i,k)-exptdn(i,k)+exptdn(i,k) * &
                           rupdir(i,k)*rdndif(i,k))*rdenom
          end do
        end do
      end do
!
!     Compute flux divergence in each layer using the interface up
!     and down fluxes:
!
      do k = 0 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            flxdiv(i,k) = (fluxdn(i,k)-fluxdn(i,k+1))+(fluxup(i,k+1)-fluxup(i,k))
          end do
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
      do n = 1 , nloop
        do i = is(n) , ie(n)
          solflx(i) = solin(i)*frcsol(ns)*psf
          fsnt(i) = fsnt(i) + solflx(i)*(fluxdn(i,1)-fluxup(i,1))
   
          fsns(i) = fsns(i) + solflx(i)*(fluxdn(i,kzp1)-fluxup(i,kz + 1))
   
          sfltot(i) = sfltot(i) + solflx(i)
          fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
          fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
!
!           Down spectral fluxes need to be in mks; thus the 0.001
!           conversion factors
          if ( wavmid < 0.7D0 ) then
            sols(i) = sols(i) + exptdn(i,kzp1)*solflx(i)*d_r1000
            solsd(i) = solsd(i) + (fluxdn(i,kzp1)-exptdn(i,kz + 1))*solflx(i)*d_r1000
!KN         added below
            sabveg(i) = sabveg(i) + (solflx(i) *            &
                       (fluxdn(i,kzp1)-fluxup(i,kz + 1)))*  &
                       (d_one-albvs(i))/(d_one-albdir(i))*d_r1000
!KN         added above
          else
            soll(i) = soll(i) + exptdn(i,kzp1)*solflx(i)*d_r1000
            solld(i) = solld(i) + (fluxdn(i,kzp1)-exptdn(i,kz + 1))*solflx(i)*d_r1000
            fsnirtsq(i) = fsnirtsq(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
!KN         added below
            sabveg(i) = sabveg(i)+(solflx(i)*(fluxdn(i,kzp1)-fluxup(i,kz + 1)))* &
                        (d_one-albvl(i))/(d_one-albdir(i))*d_r1000
!KN         added above
          end if
          fsnirt(i) = fsnirt(i) + wgtint*solflx(i) * (fluxdn(i,0)-fluxup(i,0))
   
!
        end do
      end do
      do k = 0 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            totfld(i,k) = totfld(i,k) + solflx(i)*flxdiv(i,k)
            fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
            fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
          end do
        end do
      end do
   
!     solis is incident visible solar radiation
      if ( ns == 8 ) then
!       -trapuv
!       do i=1,iym1
!         solis(i) = solflx(i)*0.001D0*fluxdn(i,kzp1)
!       end do
!       -trapuv_
        do n = 1 , nloop
          do i = is(n) , ie(n)
            solvs(i) = exptdn(i,kzp1)*solflx(i)*d_r1000
            solvd(i) = (fluxdn(i,kzp1)-exptdn(i,kz + 1))*solflx(i)*d_r1000
            solis(i) = solflx(i)*d_r1000*fluxdn(i,kzp1)
          end do
        end do
      end if
!EES    apr 20
   
!FAB
!     CLEAR SKY CALCULATION PLUS AEROSOL
!     FORCING RAD CLR is called 2 times , one with O aerosol OP , and
!     one with actual aerosol. DIFFERENCE  in net TOA SW for the two
!     case is saved as one more variable in the rad file. The
!     outputed TOASW ( fsntc, clrst) is accounting for aerosol.
      if ( ichem == 1 .and. idirect >= 1 ) then
   
        do i = 1 , iym1
          zero(i,1) = d_zero
          zero(i,2) = d_zero
          zero(i,3) = d_zero
          zero(i,4) = d_zero
        end do
   
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
        call radclr(coszrs,trayoslp,pflx,ns,uth2o,uto3,utco2,uto2, &
                    zero(1,1),zero(1,2),zero(1,3),zero(1,4),nloop, &
                    is,ie,rdir,rdif,tdir,tdif,explay,exptdn,rdndif,tottrn)
!
!       Compute reflectivity to direct and diffuse mod_radiation for
!       entire column; 0,1 on layer quantities refers to two
!       effective layers overlying surface; 0 on interface quantities
!       refers to top of column; 2 on interface quantities refers to
!       the surface:
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rupdir(i,2) = albdir(i)
            rupdif(i,2) = albdif(i)
          end do
        end do
!
        do k = 1 , 0 , -1
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = d_one/(d_one-rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k) *                    &
                            (rupdir(i,k+1)*explay(i,k)+rupdif(i,k+1) * &
                            (tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)**d_two*rdenom
            end do
          end do
        end do
!
!       Compute up and down fluxes for each interface, using the added
!       atmospheric layer properties at each interface:
!
        do k = 0 , 2
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = d_one/(d_one-rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+(tottrn(i,k) - &
                            exptdn(i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k) = exptdn(i,k) +                           &
                            (tottrn(i,k)-exptdn(i,k)+exptdn(i,k) *  &
                             rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
!
        do n = 1 , nloop
          do i = is(n) , ie(n)
            x0fsntc(i) = x0fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
            x0fsnsc(i) = x0fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
            x0fsnrtc(i) = x0fsnrtc(i) + wgtint*solflx(i)*(fluxdn(i,0)-fluxup(i,0))
   
!           SAVE the ref net TOA flux ( and put back the cumul variables to 0.)
          end do
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
      call radclr(coszrs,trayoslp,pflx,ns,uth2o,uto3,utco2,uto2, &
                  tauxar_mix_cs(:,ns),tauasc_mix_cs(:,ns),       &
                  gtota_mix_cs(:,ns),ftota_mix_cs(:,ns),nloop,   &
                  is,ie,rdir,rdif,tdir,tdif,explay,exptdn,rdndif,&
                  tottrn)
!
!     Compute reflectivity to direct and diffuse mod_radiation for entire
!     column; 0,1 on layer quantities refers to two effective layers
!     overlying surface; 0 on interface quantities refers to top of
!     column; 2 on interface quantities refers to the surface:
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          rupdir(i,2) = albdir(i)
          rupdif(i,2) = albdif(i)
        end do
      end do
!
      do k = 1 , 0 , -1
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rdenom = d_one/(d_one-rdif(i,k)*rupdif(i,k+1))
            rupdir(i,k) = rdir(i,k) + tdif(i,k) *     &
                          (rupdir(i,k+1)*explay(i,k)+ &
                           rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
            rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)**d_two*rdenom
          end do
        end do
      end do
!
!     Compute up and down fluxes for each interface, using the added
!     atmospheric layer properties at each interface:
!
      do k = 0 , 2
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rdenom = d_one/(d_one-rdndif(i,k)*rupdif(i,k))
            fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+(tottrn(i,k) - &
                           exptdn(i,k))*rupdif(i,k))*rdenom
            fluxdn(i,k) = exptdn(i,k) +                           &
                          (tottrn(i,k)-exptdn(i,k)+exptdn(i,k) *  &
                           rupdir(i,k)*rdndif(i,k))*rdenom
          end do
        end do
      end do
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          fsntc(i) = fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
          fsnsc(i) = fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
          fsnrtc(i) = fsnrtc(i) + wgtint*solflx(i) * (fluxdn(i,0)-fluxup(i,0))
   
        end do
      end do
!
!     End of clear sky calculation
!
    end do  ! End of spectral interval loop
   
!   FAB calculation of TOA aerosol radiative forcing
    if ( ichem==1 .and. idirect >= 1 ) then
      do n = 1 , nloop
        do i = is(n) , ie(n)
          aeradfo(i) = -(x0fsntc(i)-fsntc(i))
          aeradfos(i) = -(x0fsnsc(i)-fsnsc(i))
        end do
      end do
    end if
!
!   Compute solar heating rate (k/s)
!
    do k = 1 , kz
      do n = 1 , nloop
        do i = is(n) , ie(n)
          qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k)-pint(i,k+1))
        end do
      end do
    end do
!
!   Set the downwelling flux at the surface
!
    do i = 1 , iym1
      fsds(i) = fswdn(i,kzp1)
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
  subroutine radclw(jslc,ts,tnm,qnm,o3vmr,pmid,pint,pmln,piln,plco2,&
                    plh2o,n2o,ch4,cfc11,cfc12,cld,tclrsf,qrl,flns,  &
                    flnt,flnsc,flntc,flwds,fslwdcs,emiss,aerlwfo,   &
                    aerlwfos)
!
    implicit none
!
!     Input arguments
!
! ts      - Ground (skin) temperature
! emiss   - Emissivity of surface
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
    integer :: jslc
    real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           o3vmr , pmid , pmln , qnm , qrl , tnm
    real(8) , dimension(iym1,kzp1) :: cld , piln , pint , plco2 , plh2o , tclrsf
    real(8) , dimension(iym1) :: emiss , flns , flnsc , flnt ,   &
                                  flntc , flwds , fslwdcs , ts
    real(8), dimension(iym1) :: aerlwfo , aerlwfos

    intent (in) cld , emiss
    intent (out) flns , flnsc , flnt , flntc , flwds , qrl , aerlwfo , aerlwfos

    intent (inout) tclrsf
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! k2      - Level index
! k3      - Level index
! km      - Level index
! km1     - Level index
! km2     - Level index
! km3     - Level index
! km4     - Level index
! tmp     - Temporary
! tmp1    - Temporary 1
! absbt   - Downward emission at model top
! plol    - O3 pressure wghted path length
! plos    - O3 path length
! co2em   - Layer co2 normalized planck funct. derivative
! co2eml  - Interface co2 normalized planck funct. deriv.
! delt    - Diff t**4 mid layer to top interface
! delt1   - Diff t**4 lower intrfc to mid layer
! bk1     - Absrptvty for vertical quadrature
! bk2     - Absrptvty for vertical quadrature
! ful     - Total upwards longwave flux
! fsul    - Clear sky upwards longwave flux
! fdl     - Total downwards longwave flux
! fsdl    - Clear sky downwards longwv flux
! fclb4   - Sig t**4 for cld bottom interfc
! fclt4   - Sig t**4 for cloud top interfc
! s       - Flx integral sum
! tplnka  - Planck fnctn temperature
! s2c     - H2o cont amount
! s2t     - H2o cont temperature
! w       - H2o path
! tplnke  - Planck fnctn temperature
! h2otr   - H2o trnmsn for o3 overlap
! co2t    - Prs wghted temperature path
! tint    - Interface temperature
! tint4   - Interface temperature**4
! tlayr   - Level temperature
! tlayr4  - Level temperature**4
! rtclrsf - d_one/tclrsf(i,k)
! klov    - Cloud lowest level index
! khiv    - Cloud highest level index
! khivm   - khiv(i) - 1
!
!     Trace gas variables
!
! ucfc11  - CFC11 path length
! ucfc12  - CFC12 path length
! un2o0   - N2O path length
! un2o1   - N2O path length (hot band)
! uch4    - CH4 path length
! uco211  - CO2 9.4 micron band path length
! uco212  - CO2 9.4 micron band path length
! uco213  - CO2 9.4 micron band path length
! uco221  - CO2 10.4 micron band path length
! uco222  - CO2 10.4 micron band path length
! uco223  - CO2 10.4 micron band path length
! bn2o0   - pressure factor for n2o
! bn2o1   - pressure factor for n2o
! bch4    - pressure factor for ch4
! uptype  - p-type continuum path length
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
!
    real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
    real(8) , dimension(iym1) :: absbt , bk1 , bk2 , delt , delt1 , tmp , tplnke
    real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,  &
           co2t , fdl , fsdl , fsul , ful , h2otr , plol , plos ,     &
           rtclrsf , s2c , s2t , tint , tint4 , tlayr , tlayr4 ,      &
           tplnka , ucfc11 , ucfc12 , uch4 , uco211 , uco212 ,        &
           uco213 , uco221 , uco222 , uco223 , un2o0 , un2o1 ,        &
           uptype , w , fsul0 , fsdl0 , ful0 , fdl0
    real(8) , dimension(iym1,kz) :: co2eml , fclb4 , fclt4
    logical , dimension(iym1) :: done , start
    real(8) :: tmp1
    integer :: i , ii , k , k1 , k2 , k3 , khighest , km , km1 , km2 , &
               km3 , km4 , iym1c , irad , n , nradaer
    integer , dimension(iym1) :: indx , khiv , khivm , klov
    real(8) , dimension(iym1,kzp1,kzp1) :: s , s0
!
!
    character (len=50) :: subroutine_name='radclw'
    integer :: idindx = 0
!
    call time_begin(subroutine_name,idindx)
!
    do i = 1 , iym1
      rtclrsf(i,1) = d_one/tclrsf(i,1)
    end do
!
    do k = 1 , kz
      do i = 1 , iym1
        fclb4(i,k) = d_zero
        fclt4(i,k) = d_zero
        tclrsf(i,k+1) = tclrsf(i,k)*(d_one-cld(i,k+1))
        rtclrsf(i,k+1) = d_one/tclrsf(i,k+1)
      end do
    end do
!
!   Calculate some temperatures needed to derive absorptivity and
!   emissivity, as well as some h2o path lengths
!
    call radtpl(tnm,ts,qnm,pint,plh2o,tplnka,s2c,s2t,w,tplnke,tint, &
                tint4,tlayr,tlayr4,pmln,piln)
   
!   do emissivity and absorptivity calculations
!   only if abs/ems computation
!
    if ( ktau == 0 .or. (mod(ktau+1,ifrabe) == 0) ) then
 
!
!     Compute ozone path lengths at frequency of a/e calculation.
!
      call radoz2(o3vmr,pint,plol,plos)
!
!     Compute trace gas path lengths
!
      call trcpth(tnm,pint,cfc11,cfc12,n2o,ch4,qnm,ucfc11,ucfc12,     &
                  un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,uco222,&
                  uco223,bn2o0,bn2o1,bch4,uptype)
!
!
!     Compute total emissivity:
!
      call radems(s2c,s2t,w,tplnke,plh2o,pint,plco2,tint,tint4,tlayr, &
                  tlayr4,plol,plos,ucfc11,ucfc12,un2o0,un2o1,uch4,    &
                  uco211,uco212,uco213,uco221,uco222,uco223,uptype,   &
                  bn2o0,bn2o1,bch4,co2em,co2eml,co2t,h2otr,abplnk1,   &
                  abplnk2,jslc)
!
!     Compute total absorptivity:
!
      call radabs(pmid,pint,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,plco2,&
                  plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,ucfc11,   &
                  ucfc12,un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,&
                  uco222,uco223,uptype,bn2o0,bn2o1,bch4,abplnk1,      &
                  abplnk2,jslc)

       if (ichem /= 0 .and. idirect > 0) then
         abstot0(:,:,:,jslc) = abstot(:,:,:,jslc)
         emstot0(:,:,jslc) = emstot(:,:,jslc)
         absnxt0(:,:,:,jslc) = absnxt(:,:,:,jslc)   
      end if
    end if
!
!   Find the lowest and highest level cloud for each grid point
!   Note: Vertical indexing here proceeds from bottom to top
!
    do i = 1 , iym1
      klov(i) = 0
      done(i) = .false.
    end do
    do k = 1 , kz
      do i = 1 , iym1
        if ( .not.done(i) .and. cld(i,kzp2-k) > d_zero ) then
          done(i) = .true.
          klov(i) = k
        end if
      end do
    end do
    call whenne(iym1,klov,1,0,indx,iym1c)
    do i = 1 , iym1
      khiv(i) = klov(i)
      done(i) = .false.
    end do
    do k = kz , 1 , -1
      do ii = 1 , iym1c
        i = indx(ii)
        if ( .not.done(i) .and. cld(i,kzp2-k) > d_zero ) then
          done(i) = .true.
          khiv(i) = k
        end if
      end do
    end do
    do i = 1 , iym1
      khivm(i) = khiv(i) - 1
    end do
!
!   Note: Vertical indexing here proceeds from bottom to top
!
    do ii = 1 , iym1c
      i = indx(ii)
      do k = klov(i) , khiv(i)
        fclt4(i,kzp1-k) = stebol*tint4(i,kzp2-k)
        fclb4(i,kzp1-k) = stebol*tint4(i,kzp3-k)
      end do
    end do

!
!   option to calculate LW aerosol radiative forcing
!
!   FAB LW radiative forcing ( rad=1 : avec dust)
    if (ichem /= 0 .and. idirect > 0) then
      nradaer = 2
      fsul0(:,:) = d_zero
      fsdl0(:,:) = d_zero
      abstot(:,:,:,jslc) = abstot0(:,:,:,jslc)
      emstot(:,:,jslc) = emstot0(:,:,jslc)
      absnxt(:,:,:,jslc) = absnxt0(:,:,:,jslc)
    else
      nradaer = 1
    end if

    do irad = 1 , nradaer

      if (ichem==1 .and. idirect > 0 .and. irad==2 ) then
        abstot(:,:,:,jslc) = d_one-(d_one-abstot0(:,:,:,jslc))*aerlwtr(:,:,:)
        emstot(:,:,jslc) = d_one-(d_one-emstot0(:,:,jslc))*aerlwtr(:,:,1)
        do k = 1 , kz  ! aerlwtr defined on plev levels
          do n = 1 , 4
            absnxt(:,k,n,jslc) = d_one-(d_one-absnxt0(:,k,n,jslc)) *   &
                                (aerlwtr(:,k,k+1)**xuinpl(:,k,n,jslc))
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
      do i = 1 , iym1
        delt(i) = tint4(i,kz) - tlayr4(i,kzp1)
        delt1(i) = tlayr4(i,kzp1) - tint4(i,kzp1)
        s(i,kzp1,kzp1) = stebol*(delt1(i)*absnxt(i,kz,1,jslc) + &
                         delt(i)*absnxt(i,kz,4,jslc))
        s(i,kz,kzp1) = stebol*(delt(i)*absnxt(i,kz,2,jslc) + &
                               delt1(i)*absnxt(i,kz,3,jslc))
      end do
      do k = 1 , kz - 1
        do i = 1 , iym1
          bk2(i) = (abstot(i,k,kz,jslc)+abstot(i,k,kzp1,jslc))*d_half
          bk1(i) = bk2(i)
          s(i,k,kzp1) = stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
        end do
      end do
!
!     All k, km>1
!
      do km = kz , 2 , -1
        do i = 1 , iym1
          delt(i) = tint4(i,km-1) - tlayr4(i,km)
          delt1(i) = tlayr4(i,km) - tint4(i,km)
        end do
        do k = kzp1 , 1 , -1
          if ( k == km ) then
            do i = 1 , iym1
              bk2(i) = absnxt(i,km-1,4,jslc)
              bk1(i) = absnxt(i,km-1,1,jslc)
            end do
          else if ( k == km-1 ) then
            do i = 1 , iym1
              bk2(i) = absnxt(i,km-1,2,jslc)
              bk1(i) = absnxt(i,km-1,3,jslc)
            end do
          else
            do i = 1 , iym1
              bk2(i) = (abstot(i,k,km-1,jslc)+ &
                        abstot(i,k,km,jslc))*d_half
              bk1(i) = bk2(i)
            end do
          end if
          do i = 1 , iym1
            s(i,k,km) = s(i,k,km+1)+stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
          end do
        end do
      end do
!
!     Computation of clear sky fluxes always set first level of fsul
!
      do i = 1 , iym1
        if ( iemiss == 1 ) then
          fsul(i,kzp1) = emiss(i)*(stebol*(ts(i)**d_four))
        else
          fsul(i,kzp1) = stebol*(ts(i)**d_four)
        end if
      end do
!
!     Downward clear sky fluxes store intermediate quantities in down
!     flux Initialize fluxes to clear sky values.
!
      do i = 1 , iym1
        tmp(i) = fsul(i,kzp1) - stebol*tint4(i,kzp1)
        fsul(i,1) = fsul(i,kzp1) - abstot(i,1,kzp1,jslc)*tmp(i)+s(i,1,2)
        fsdl(i,1) = stebol*(tplnke(i)**d_four)*emstot(i,1,jslc)
        ful(i,1) = fsul(i,1)
        fdl(i,1) = fsdl(i,1)
      end do
!
!     fsdl(i,kzp1) assumes isothermal layer
!
      do k = 2 , kz
        do i = 1 , iym1
          fsul(i,k) = fsul(i,kzp1) - abstot(i,k,kzp1,jslc)*tmp(i)+s(i,k,k+1)
          ful(i,k) = fsul(i,k)
          fsdl(i,k) = stebol*(tplnke(i)**d_four)*emstot(i,k,jslc) - &
                              (s(i,k,2)-s(i,k,k+1))
          fdl(i,k) = fsdl(i,k)
        end do
      end do
!
!     Store the downward emission from level 1 = total gas emission *
!     sigma t**4.  fsdl does not yet include all terms
!
      do i = 1 , iym1
        ful(i,kzp1) = fsul(i,kzp1)
        absbt(i) = stebol*(tplnke(i)**d_four)*emstot(i,kzp1,jslc)
        fsdl(i,kzp1) = absbt(i) - s(i,kzp1,2)
        fdl(i,kzp1) = fsdl(i,kzp1)
      end do

!     FAB radiative forcing sur fsul

      if (ichem==1 .and. idirect > 0 .and. irad==1 ) then
        fsul0(:,:) = fsul(:,:)! save fsul0 = no dust
        fsdl0(:,:) = fsdl(:,:)!
        ful0(:,:) = ful(:,:)
        fdl0(:,:) = fdl(:,:)
        s0(:,:,:) = s(:,:,:)
      end if

    end do ! end rad loop

!   FAB after this DO loop fsul account for dust LW effect
!   which is OK in case of idirect=2

    if (ichem ==1 .and. idirect > 0 ) then

      aerlwfo(:) = fsul0(:,1) - fsul(:,1)

!     surface lw net ! fsul(i,plevp) - fsdl(i,plevp)
!     aerlwfos(:)= fsdl0(:,kz)-fsdl(:,kz)
      aerlwfos(:) = (fsul0(:,kzp1)-fsdl0(:,kzp1))-                    &
                    (fsul(:,kzp1) - fsdl(:,kzp1))
       
!     return to no aerosol LW effect  situation if idirect ==1
      if ( idirect==1 ) then
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
    call whenflt(iym1,tclrsf(1,kzp1),1,verynearone,indx,iym1c)
!
!   Compute downflux at level 1 for cloudy sky
!
    do ii = 1 , iym1c
      i = indx(ii)
!
!     First clear sky flux plus flux from cloud at level 1
!
      fdl(i,kzp1) = fsdl(i,kzp1)*tclrsf(i,kz) * &
                    rtclrsf(i,kzp1-khiv(i))+fclb4(i,kz-1)*cld(i,kz)
    end do
!
!   Flux emitted by other layers
!   Note: Vertical indexing here proceeds from bottom to top
!
    khighest = khiv(intmax(iym1,khiv,1))
    do km = 3 , khighest
      km1 = kzp1 - km
      km2 = kzp2 - km
      km4 = kzp4 - km
      do ii = 1 , iym1c
        i = indx(ii)
        if ( km <= khiv(i) ) then
          tmp1 = cld(i,km2)*tclrsf(i,kz)*rtclrsf(i,km2)
          fdl(i,kzp1) = fdl(i,kzp1) + (fclb4(i,km1)-s(i,kzp1,km4))*tmp1
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
      do ii = 1 , iym1c
        i = indx(ii)
        if ( k >= klov(i) .and. k <= khivm(i) ) then
          ful(i,k2) = fsul(i,k2)*(tclrsf(i,kzp1)*rtclrsf(i,k1))
        end if
      end do
      do km = 1 , k
        km1 = kzp1 - km
        km2 = kzp2 - km
        km3 = kzp3 - km
        do ii = 1 , iym1c
          i = indx(ii)
!
          if ( k <= khivm(i) .and. km >= klov(i) .and. km <= khivm(i)) then
            ful(i,k2) = ful(i,k2) + (fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3)) * &
                        cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
          end if
        end do
      end do ! km = 1 , k
    end do   ! k = 1 , khighest-1
!
    do k = 1 , kzp1
      k2 = kzp2 - k
      k3 = kzp3 - k
      do i = 1 , iym1
        start(i) = .false.
      end do
      do ii = 1 , iym1c
        i = indx(ii)
        if ( k >= khiv(i) ) then
          start(i) = .true.
          ful(i,k2) = fsul(i,k2)*tclrsf(i,kzp1)*rtclrsf(i,kzp1-khiv(i))
        end if
      end do
      do km = 1 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km3 = kzp3 - km
        do ii = 1 , iym1c
          i = indx(ii)
          if ( start(i) .and. km >= klov(i) .and. km <= khiv(i) ) then
            ful(i,k2) = ful(i,k2) + (cld(i,km2)*tclrsf(i,km1)* &
                  rtclrsf(i,kzp1-khiv(i)))*(fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3))
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
      do ii = 1 , iym1c
        i = indx(ii)
        if ( k <= khivm(i) ) fdl(i,k2) = d_zero
      end do
      do km = k + 1 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km4 = kzp4 - km
        do ii = 1 , iym1c
          i = indx(ii)
!
          if ( k <= khiv(i) .and. &
               km >= max0(k+1,klov(i)) .and. km <= khiv(i) ) then
            fdl(i,k2) = fdl(i,k2)+(cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2)) * &
                        (fclb4(i,km1)-s(i,k2,km4)+s(i,k2,k3))
          end if
        end do
      end do ! km = k+1 , khighest
      do ii = 1 , iym1c
        i = indx(ii)
        if ( k <= khivm(i) ) then
          fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)*(tclrsf(i,k1)*rtclrsf(i,kzp1-khiv(i)))
        end if
      end do
    end do  ! k = 1 , khighest-1
!
!   End cloud modification loops
!
!   All longitudes: store history tape quantities
!
    do i = 1 , iym1
!
!     Downward longwave flux
!
      flwds(i) = fdl(i,kzp1)
!
!     Net flux
!
      flns(i) = ful(i,kzp1) - fdl(i,kzp1)
!
!     Clear sky flux at top of atmosphere
!
      flntc(i) = fsul(i,1)
      flnsc(i) = fsul(i,kzp1) - fsdl(i,kzp1)
!
!     Outgoing ir
!
      flnt(i) = ful(i,1) - fdl(i,1)
    end do
!
!   Computation of longwave heating (k per sec)
!
    do k = 1 , kz
      do i = 1 , iym1
        qrl(i,k) = (ful(i,k)-fdl(i,k)-ful(i,k+1)+fdl(i,k+1))*gocp / &
                    ((pint(i,k)-pint(i,k+1)))
      end do
    end do
!
    call time_end(subroutine_name,idindx)
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
  subroutine radclr(coszrs,trayoslp,pflx,ns,uth2o,uto3,utco2,uto2, &
                    tauxar_mix_css,tauasc_mix_css,gtota_mix_css,   &
                    ftota_mix_css,nloop,is,ie,rdir,rdif,tdir,tdif, &
                    explay,exptdn,rdndif,tottrn)
!
    implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! coszrs   - Cosine zenith angle
! trayoslp - Tray/sslp
! pflx     - Interface pressure
! uth2o    - Total column absorber amount of h2o
! uto3     - Total column absorber amount of  o3
! utco2    - Total column absorber amount of co2
! uto2     - Total column absorber amount of  o2
! nloop    - Number of loops (1 or 2)
! is       - Starting index for 1 or 2 loops
! ie       - Ending index for 1 or 2 loops
!
!     Input/Output arguments
!
!     Following variables are defined for each layer; note, we use
!     layer 0 to refer to the entire atmospheric column:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer refflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transmn for layer
!
!     Note that the following variables are defined on interfaces, with
!     the index k referring to the top interface of the kth layer:
!     exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
!     transmission to the top interface of the 5th layer.
!
! exptdn   - Solar beam exp down transmn from top
! rdndif   - Added dif ref for layers above
! tottrn   - Total transmission for layers above
!
!
    real(8) :: trayoslp
    integer :: nloop , ns
    real(8) , dimension(iym1) :: coszrs , ftota_mix_css ,            &
                                  gtota_mix_css , tauasc_mix_css ,    &
                                  tauxar_mix_css , utco2 , uth2o ,    &
                                  uto2 , uto3
    real(8) , dimension(iym1,0:kz) :: explay , rdif , rdir , tdif ,&
           tdir
    real(8) , dimension(iym1,0:kzp1) :: exptdn , pflx , rdndif ,    &
           tottrn
    integer , dimension(2) :: ie , is
    intent (in) coszrs , ftota_mix_css ,&
                gtota_mix_css , ie , is , nloop , pflx ,              &
                tauasc_mix_css , tauxar_mix_css , trayoslp , utco2 ,  &
                uth2o , uto2 , uto3 , ns
    intent (inout) explay , exptdn , rdif , rdir , rdndif , tdif ,    &
                   tdir , tottrn
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! nn       - Index of longitude loops (max=nloop)
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
    real(8) :: alp , amg , apg , arg , extins , ftot ,   &
               gam , gs , gtot , lm , ne , rdenom , rdirexp , &
               tautot , tdnmexp , ts , ue , ws , wtot
    integer :: i , ii , k , nn , nval
    integer , dimension(iym1) :: indx
    real(8) , dimension(iym1) :: taugab , tauray
!
    character (len=50) :: subroutine_name='radclw'
    integer :: idindx = 0
!
    call time_begin(subroutine_name,idindx)
!-----------------------------------------------------------------------
!
!   Initialize all total transmimission values to 0, so that nighttime
!   values from previous computations are not used:
!
    tottrn = d_zero
!   print*,'dans radclr', maxval(tottrn)
    do k = 0 , kzp1
      do i = 1 , iym1
        tottrn(i,k) = d_zero
      end do
    end do
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above each interface by starting from the top and adding layers down:
!
!   The top layer is assumed to be a purely absorbing ozone layer, and
!   that the mean diffusivity for diffuse mod_transmission is 1.66:
!
    do nn = 1 , nloop
      do i = is(nn) , ie(nn)
!
        taugab(i) = abo3(ns)*uto3(i)
!
!       Limit argument of exponential to 25, in case coszrs is very small:
        arg = dmin1(taugab(i)/coszrs(i),25.0D0)
        explay(i,0) = dexp(-arg)
        tdir(i,0) = explay(i,0)
!
!       Same limit for diffuse mod_transmission:
!
        arg = dmin1(1.66D0*taugab(i),25.0D0)
        tdif(i,0) = dexp(-arg)
!
        rdir(i,0) = d_zero
        rdif(i,0) = d_zero
!
!       Initialize top interface of extra layer:
!
        exptdn(i,0) = d_one
        rdndif(i,0) = d_zero
        tottrn(i,0) = d_one
!
        rdndif(i,1) = rdif(i,0)
        tottrn(i,1) = tdir(i,0)
!
      end do
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
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
!
          rdir(i,k) = d_zero
          rdif(i,k) = d_zero
          tdir(i,k) = d_zero
          tdif(i,k) = d_zero
          explay(i,k) = d_zero
!
!         Calculates the solar beam transmission, total transmission,
!         and reflectivity for diffuse mod_radiation from below at the
!         top of the current layer:
!
          exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
          rdenom = d_one/(d_one-rdif(i,k-1)*rdndif(i,k-1))
          rdirexp = rdir(i,k-1)*exptdn(i,k-1)
          tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
          tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + &
                        tdif(i,k-1)*(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
          rdndif(i,k) = rdif(i,k-1) + &
                        (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
!
        end do
      end do
!
!     Compute next layer delta-Eddington solution only if total
!     transmission of radiation to the interface just above the layer
!     exceeds trmin.
      call whenfgt(iym1,tottrn(1,k),1,trmin,indx,nval)
      if ( nval > 0 ) then
!CDIR$  IVDEP
        do ii = 1 , nval
          i = indx(ii)
!
!         Remember, no ozone absorption in this layer:
!
          tauray(i) = trayoslp*pflx(i,kzp1)
          taugab(i) = abh2o(ns)*uth2o(i) + abco2(ns)*utco2(i) + abo2(ns)*uto2(i)
!
          tautot = tauray(i) + taugab(i) + tauxar_mix_css(i)
!
          wtot = (wray*tauray(i)+tauasc_mix_css(i))/tautot
!
          gtot = (gray*wray*tauray(i)+gtota_mix_css(i))/(wtot*tautot)
!
          ftot = (fray*wray*tauray(i)+ftota_mix_css(i))/(wtot*tautot)
!
          ts = taus(wtot,ftot,tautot)
          ws = omgs(wtot,ftot)
          gs = asys(gtot,ftot)
          lm = el(ws,gs)
          alp = xalpha(ws,coszrs(i),gs,lm)
          gam = xgamma(ws,coszrs(i),gs,lm)
          ue = u(ws,gs,lm)
!
!         Limit argument of exponential to 25, in case lm very large:
!
          arg = dmin1(lm*ts,25.0D0)
          extins = dexp(-arg)
          ne = n(ue,extins)
!
          rdif(i,k) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
          tdif(i,k) = d_four*ue/ne
!
!         Limit argument of exponential to 25, in case coszrs is very small:
          arg = dmin1(ts/coszrs(i),25.0D0)
          explay(i,k) = dexp(-arg)
!
          apg = alp + gam
          amg = alp - gam
          rdir(i,k) = amg*(tdif(i,k)*explay(i,k)-d_one)+apg*rdif(i,k)
          tdir(i,k) = apg*tdif(i,k) + (amg*rdif(i,k)-(apg-d_one))*explay(i,k)
!
!         Under rare conditions, reflectivies and transmissivities
!         can be negative; zero out any negative values
!
          rdir(i,k) = dmax1(rdir(i,k),d_zero)
          tdir(i,k) = dmax1(tdir(i,k),d_zero)
          rdif(i,k) = dmax1(rdif(i,k),d_zero)
          tdif(i,k) = dmax1(tdif(i,k),d_zero)
!
        end do
      end if
!
    end do
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for both layers
!   above the surface:
!
    k = 2
    do nn = 1 , nloop
      do i = is(nn) , ie(nn)
        exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
        rdenom = d_one/(d_one-rdif(i,k-1)*rdndif(i,k-1))
        rdirexp = rdir(i,k-1)*exptdn(i,k-1)
        tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
        tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + &
                      tdif(i,k-1)*(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
        rdndif(i,k) = rdif(i,k-1) + &
                      (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
      end do
    end do
!
    call time_end(subroutine_name,idindx)
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
  subroutine radded(coszrs,trayoslp,pflx,ns,uh2o,uo3,uco2,uo2,tauxcl, &
                    wcl,gcl,fcl,tauxci,wci,gci,fci,tauxar_mixs,       &
                    tauasc_mixs,gtota_mixs,ftota_mixs,nloop,is,ie,    &
                    rdir,rdif,tdir,tdif,explay,exptdn,rdndif,tottrn)
!
    implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! coszrs   - Cosine zenith angle
! trayoslp - Tray/sslp
! pflx     - Interface pressure
! uh2o     - Layer absorber amount of h2o
! uo3      - Layer absorber amount of  o3
! uco2     - Layer absorber amount of co2
! uo2      - Layer absorber amount of  o2
! tauxcl   - Cloud extinction optical depth
! wcl      - Cloud single scattering albedo
! gcl      - Cloud assymetry parameter
! fcl      - Cloud forward scattered fraction
! tauxci   - Cloud extinction optical depth
! wci      - Cloud single scattering albedo
! gci      - Cloud assymetry parameter
! fci      - Cloud forward scattered fraction
! nloop    - Number of loops (1 or 2)
! is       - Starting index for 1 or 2 loops
! ie       - Ending index for 1 or 2 loops
!
!     Input/Output arguments
!
!     Following variables are defined for each layer; 0 refers to extra
!     layer above top of model:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer refflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transm for layer
!
!     (Note that the following variables are defined on interfaces,
!     with the index k referring to the top interface of the kth layer:
!     exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
!     transmission to the top interface of the 5th layer; kzp1 refers
!     to the earth surface)
!
! rdndif   - Added dif ref for layers above
! exptdn   - Solar beam exp down transm from top
! tottrn   - Total transmission for layers above
!
!
    real(8) :: trayoslp
    integer :: nloop , ns
    real(8) , dimension(iym1) :: coszrs
    real(8) , dimension(iym1,0:kz) :: explay , fci , fcl ,             &
           ftota_mixs , gci , gcl , gtota_mixs , rdif , rdir ,         &
           tauasc_mixs , tauxar_mixs , tauxci , tauxcl , tdif , tdir , &
           uco2 , uh2o , uo2 , uo3 , wci , wcl
    real(8) , dimension(iym1,0:kzp1) :: exptdn , pflx , rdndif , tottrn
    integer , dimension(2) :: ie , is
    intent (in) coszrs , fci , fcl ,    &
                ftota_mixs , gci , gcl , gtota_mixs , ie , is ,       &
                nloop , pflx , tauasc_mixs , tauxar_mixs , tauxci ,   &
                tauxcl , trayoslp , uco2 , uh2o , uo2 , uo3 , wci ,   &
                wcl , ns
    intent (inout) explay , exptdn , rdif , rdir , rdndif , tdif ,    &
                   tdir , tottrn
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! nn       - Index of longitude loops (max=nloop)
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
! w        - Dummy argument for statement function
! uu       - Dummy argument for statement function
! g        - Dummy argument for statement function
! e        - Dummy argument for statement function
! f        - Dummy argument for statement function
! t        - Dummy argument for statement function
! et       - Dummy argument for statement function
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
    real(8) :: alp , amg , apg , arg , extins , ftot ,        &
               gam , gs , gtot , lm , ne , rdenom , rdirexp , &
               taucsc , tautot , tdnmexp , ts , ue , ws ,     &
               wt , wtau , wtot
    integer :: i , ii , k , nn , nval
    integer , dimension(iym1) :: indx
    real(8) , dimension(iym1) :: taugab , tauray
!
    character (len=50) :: subroutine_name='radded'
    integer :: idindx = 0
!
    call time_begin(subroutine_name,idindx)
!-----------------------------------------------------------------------
!
!   Initialize all total transmission values to 0, so that nighttime
!   values from previous computations are not used:
!
    tottrn = d_zero
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above each interface by starting from the top and adding layers down:
!   For the extra layer above model top:
!
    do nn = 1 , nloop
      do i = is(nn) , ie(nn)
!
        tauray(i) = trayoslp*(pflx(i,1)-pflx(i,0))
        taugab(i) = abh2o(ns)*uh2o(i,0) + abo3(ns)*uo3(i,0) + &
                    abco2(ns)*uco2(i,0) + abo2(ns)*uo2(i,0)
!
        tautot = tauxcl(i,0)+tauxci(i,0)+tauray(i)+taugab(i)+tauxar_mixs(i,0)
        taucsc = tauxcl(i,0)*wcl(i,0) + tauxci(i,0)*wci(i,0)+tauasc_mixs(i,0)
        wtau = wray*tauray(i)
        wt = wtau + taucsc
        wtot = wt/tautot
        gtot = (wtau*gray+gcl(i,0)*tauxcl(i,0)*wcl(i,0)+gci(i,0) *    &
                tauxci(i,0)*wci(i,0)+gtota_mixs(i,0))/wt
        ftot = (wtau*fray+fcl(i,0)*tauxcl(i,0)*wcl(i,0)+fci(i,0) *    &
                tauxci(i,0)*wci(i,0)+ftota_mixs(i,0))/wt
!
        ts = taus(wtot,ftot,tautot)
        ws = omgs(wtot,ftot)
        gs = asys(gtot,ftot)
        lm = el(ws,gs)
        alp = xalpha(ws,coszrs(i),gs,lm)
        gam = xgamma(ws,coszrs(i),gs,lm)
        ue = u(ws,gs,lm)
!
!       Limit argument of exponential to 25, in case lm*ts very large:
!
        arg = dmin1(lm*ts,25.0D0)
        extins = dexp(-arg)
        ne = n(ue,extins)
!
        rdif(i,0) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
        tdif(i,0) = d_four*ue/ne
!
!       Limit argument of exponential to 25, in case coszrs is very small:
        arg = dmin1(ts/coszrs(i),25.0D0)
        explay(i,0) = dexp(-arg)
!
        apg = alp + gam
        amg = alp - gam
        rdir(i,0) = amg*(tdif(i,0)*explay(i,0)-d_one) + apg*rdif(i,0)
        tdir(i,0) = apg*tdif(i,0) + (amg*rdif(i,0)-(apg-d_one))*explay(i,0)
!
!       Under rare conditions, reflectivies and transmissivities can
!       be negative; zero out any negative values
!
        rdir(i,0) = dmax1(rdir(i,0),d_zero)
        tdir(i,0) = dmax1(tdir(i,0),d_zero)
        rdif(i,0) = dmax1(rdif(i,0),d_zero)
        tdif(i,0) = dmax1(tdif(i,0),d_zero)
!
!       Initialize top interface of extra layer:
!
        exptdn(i,0) = d_one
        rdndif(i,0) = d_zero
        tottrn(i,0) = d_one
!
        rdndif(i,1) = rdif(i,0)
        tottrn(i,1) = tdir(i,0)
!
      end do
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
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
!
          rdir(i,k) = d_zero
          rdif(i,k) = d_zero
          tdir(i,k) = d_zero
          tdif(i,k) = d_zero
          explay(i,k) = d_zero
!
!         Calculates the solar beam transmission, total transmission,
!         and reflectivity for diffuse mod_radiation from below at the
!         top of the current layer:
!
          exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
!KN       modified below (for computational stability)
!         rdenom      = d_one/(1. - rdif(i,k-1)*rdndif(i,k-1))
          rdenom = d_one/(d_one-dmin1(rdif(i,k-1)*rdndif(i,k-1),verynearone))
!KN       modified above
          rdirexp = rdir(i,k-1)*exptdn(i,k-1)
          tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
          tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1) *     &
                        (tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
          rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1)) *   &
                        (tdif(i,k-1)*rdenom)
!
        end do
      end do
!
!     Compute next layer delta-eddington solution only if total
!     transmission of radiation to the interface just above the layer
!     exceeds trmin.
      call whenfgt(iym1,tottrn(1,k),1,trmin,indx,nval)
      if ( nval > 0 ) then
!CDIR$  IVDEP
        do ii = 1 , nval
          i = indx(ii)
!
          tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
          taugab(i) = abh2o(ns)*uh2o(i,k) + abo3(ns)*uo3(i,k) +  &
                      abco2(ns)*uco2(i,k) + abo2(ns)*uo2(i,k)
!
          tautot = tauxcl(i,k)+tauxci(i,k)+tauray(i)+taugab(i)+tauxar_mixs(i,k)
          taucsc = tauxcl(i,k)*wcl(i,k)+tauxci(i,k)*wci(i,k)+tauasc_mixs(i,k)
          wtau = wray*tauray(i)
          wt = wtau + taucsc
          wtot = wt/tautot
          gtot = (wtau*gray+gcl(i,k)*wcl(i,k)*tauxcl(i,k)+gci(i,k) *  &
                  wci(i,k)*tauxci(i,k)+gtota_mixs(i,k))/wt
          ftot = (wtau*fray+fcl(i,k)*wcl(i,k)*tauxcl(i,k)+fci(i,k) *  &
                  wci(i,k)*tauxci(i,k)+ftota_mixs(i,k))/wt
!
          ts = taus(wtot,ftot,tautot)
          ws = omgs(wtot,ftot)
          gs = asys(gtot,ftot)
          lm = el(ws,gs)
          alp = xalpha(ws,coszrs(i),gs,lm)
          gam = xgamma(ws,coszrs(i),gs,lm)
          ue = u(ws,gs,lm)
!
!         Limit argument of exponential to 25, in case lm very large:
!
          arg = dmin1(lm*ts,25.0D0)
          extins = dexp(-arg)
          ne = n(ue,extins)
!
          rdif(i,k) = (ue+d_one)*(ue-d_one)*(d_one/extins-extins)/ne
          tdif(i,k) = d_four*ue/ne
!
!         Limit argument of exponential to 25, in case coszrs is very small:
          arg = dmin1(ts/coszrs(i),25.0D0)
          explay(i,k) = dexp(-arg)
!
          apg = alp + gam
          amg = alp - gam
          rdir(i,k) = amg*(tdif(i,k)*explay(i,k)-d_one)+apg*rdif(i,k)
          tdir(i,k) = apg*tdif(i,k) + (amg*rdif(i,k)-(apg-d_one))*explay(i,k)
!
!         Under rare conditions, reflectivies and transmissivities
!         can be negative; zero out any negative values
!
          rdir(i,k) = dmax1(rdir(i,k),d_zero)
          tdir(i,k) = dmax1(tdir(i,k),d_zero)
          rdif(i,k) = dmax1(rdif(i,k),d_zero)
          tdif(i,k) = dmax1(tdif(i,k),d_zero)
        end do
      end if
!
    end do
!
!   Compute total direct beam transmission, total transmission, and
!   reflectivity for diffuse mod_radiation (from below) for all layers
!   above the surface:
!
    k = kzp1
    do nn = 1 , nloop
      do i = is(nn) , ie(nn)
        exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
!KN     modified below (for computational stability)
!       rdenom = d_one/(1. - rdif(i,k-1)*rdndif(i,k-1))
        rdenom = d_one/(d_one-dmin1(rdif(i,k-1)*rdndif(i,k-1),verynearone))
!KN     modified above
        rdirexp = rdir(i,k-1)*exptdn(i,k-1)
        tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
        tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1) *       &
                      (tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
        rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1)) *     &
                      (tdif(i,k-1)*rdenom)
      end do
    end do
!
    call time_end(subroutine_name,idindx)
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
  subroutine radabs(pbr,pnm,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,    &
                    plco2,plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,&
                    ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212,   &
                    uco213,uco221,uco222,uco223,uptype,bn2o0,bn2o1, &
                    bch4,abplnk1,abplnk2,jslc)
!
    implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pbr      - Prssr at mid-levels (dynes/cm2)
! pnm      - Prssr at interfaces (dynes/cm2)
! co2em    - Co2 emissivity function
! co2eml   - Co2 emissivity function
! tplnka   - Planck fnctn level temperature
! s2c      - H2o continuum path length
! s2t      - H2o tmp and prs wghted path
! w        - H2o prs wghted path
! h2otr    - H2o trnsmssn fnct for o3 overlap
! plco2    - Co2 prs wghted path length
! plh2o    - H2o prs wfhted path length
! co2t     - Tmp and prs wghted path length
! tint     - Interface temperatures
! tlayr    - K-1 level temperatures
! plol     - Ozone prs wghted path length
! plos     - Ozone path length
! pmln     - Ln(pmidm1)
! piln     - Ln(pintm1)
!
!     Trace gas variables
!
! ucfc11   - CFC11 path length
! ucfc12   - CFC12 path length
! un2o0    - N2O path length
! un2o1    - N2O path length (hot band)
! uch4     - CH4 path length
! uco211   - CO2 9.4 micron band path length
! uco212   - CO2 9.4 micron band path length
! uco213   - CO2 9.4 micron band path length
! uco221   - CO2 10.4 micron band path length
! uco222   - CO2 10.4 micron band path length
! uco223   - CO2 10.4 micron band path length
! uptype   - continuum path length
! bn2o0    - pressure factor for n2o
! bn2o1    - pressure factor for n2o
! bch4     - pressure factor for ch4
! abplnk1   - non-nearest layer Plack factor
! abplnk2   - nearest layer factor
! abstrc    - total trace gas absorptivity
! bplnk     - Planck functions for sub-divided layers
!
!     Output arguments (radbuf)
!
! Nearest layer absorptivities
! Non-adjacent layer absorptivites
! Total emissivity
!
    integer :: jslc
    real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
    real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           co2t , h2otr , piln , plco2 , plh2o , plol , plos , pnm ,  &
           s2c , s2t , tint , tlayr , tplnka , ucfc11 , ucfc12 ,      &
           uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,        &
           uco223 , un2o0 , un2o1 , uptype , w
    real(8) , dimension(iym1,kz) :: co2eml , pbr , pmln
    intent (in) abplnk2 , co2em , co2eml , co2t , h2otr , jslc , pbr ,&
                piln , plco2 , plh2o , plol , plos , pmln , s2t ,     &
                tint , tlayr , tplnka , w
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
! u        - Pressure weighted H2O path length
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
! dbvtly   - Level drvtv plnck fnctn for o3
! dbvt     - Planck fnctn tmp derivative for o3
!
    real(8) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,    &
               absbnd , alphat , beta , cf812 , corfac , denom ,      &
               dplco2 , dplol , dplos , ds2c , dtym10 , et , et2 ,    &
               et4 , f1co2 , g2 , g4 , k21 , k22 , o3bndi , omet ,    &
               oneme , p1 , p2 , pbar , phi , pi , posqt , psi ,      &
               rbeta13 , rbeta8 , rbeta9 , rdpnm , rdpnmsq , realnu , &
               rphat , rsqti , rsum , sqwp , t1t4 , t2t5 ,&
               tcrfac , te , tlocal , tmp1 , tmp2 , tmp3 , tpath ,    &
               tr1 , tr2 , tr5 , tr6 , u1 , u13 , u2 , u8 , u9 ,      &
               ubar , wco2
    real(8) , dimension(iym1,6) :: abso
    real(8) , dimension(iym1) :: abstrc , dplh2o , dpnm , dtp , dtx ,  &
                                  dty , dtyp15 , dtyp15sq , dtz , dw , &
                                  f1sqwp , f2co2 , f3co2 , fwk ,       &
                                  fwku , pnew , rbeta7 , sqrtu ,       &
                                  sqti , t1co2 , tco2 , th2o , to3 ,   &
                                  to3co2 , to3h2o , tpatha , tr10 ,    &
                                  tr9 , trab2 , trab4 , trab6 , u ,    &
                                  u7 , uc , uc1
    real(8) , dimension(14,iym1,4) :: bplnk
    real(8) , dimension(iym1,kzp1) :: dbvtit , pnmsq , term6 , term9
    real(8) , dimension(iym1,kz) :: dbvtly
    real(8) , dimension(iym1,4) :: emm , o3emm , pinpl , tbar ,       &
                                    temh2o , term1 , term2 , term3 ,  &
                                    term4 , term5 , uinpl , winpl ,   &
                                    zinpl
    integer :: i , iband , k , k1 , k2 , kn , wvl
    real(8) , dimension(2) :: r2st
    real(8) , dimension(iym1,2) :: term7 , term8 , trline
!
!
    character (len=50) :: subroutine_name='radabs'
    integer :: indx = 0
!
!-----------------------------------------------------------------------
    call time_begin(subroutine_name,indx)
!
!   Initialize
!
    do k = 1 , kz
      do i = 1 , iym1
        dbvtly(i,k) = dbvt(tlayr(i,k+1))
        dbvtit(i,k) = dbvt(tint(i,k))
      end do
    end do
    do i = 1 , iym1
      dbvtit(i,kzp1) = dbvt(tint(i,kz + 1))
    end do
!
    r2st(1) = d_one/(d_two*st(1))
    r2st(2) = d_one/(d_two*st(2))
!   bndfct  = 2.0*22.18d0/(dsqrt(196.d0)*300.)
!
!   Non-adjacent layer absorptivity:
!
!   abso(i,1)     0 -  800 cm-1   h2o rotation band
!   abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!   abso(i,3)   800 - 1200 cm-1   h2o window
!   abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!   abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!   abso(i,6)   co2 15  micrometer band system
!
    do k = 1 , kzp1
      do i = 1 , iym1
        pnmsq(i,k) = pnm(i,k)**d_two
        dtx(i) = tplnka(i,k) - 250.0D0
        term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i) *          &
                      (d_one+c9*dtx(i)*(d_one+c11*dtx(i) *     &
                      (d_one+c13*dtx(i)*(d_one+c15*dtx(i)))))
        term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i) *          &
                      (d_one+c19*dtx(i)*(d_one+c21*dtx(i) *    &
                      (d_one+c23*dtx(i)*(d_one+c25*dtx(i)))))
      end do
    end do
!
!   Non-nearest layer level loops
!
    do k1 = kzp1 , 1 , -1
      do k2 = kzp1 , 1 , -1
        if ( k1 /= k2 ) then
          do i = 1 , iym1
            dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
            u(i) = dabs(dplh2o(i))
            sqrtu(i) = dsqrt(u(i))
            ds2c = dabs(s2c(i,k1)-s2c(i,k2))
            dw(i) = dabs(w(i,k1)-w(i,k2))
            uc1(i) = (ds2c+1.7D-3*u(i))*(d_one+d_two*ds2c)/(d_one+15.0D0*ds2c)
            uc(i) = ds2c + 2.0D-3*u(i)
          end do
          do i = 1 , iym1
            pnew(i) = u(i)/dw(i)
            tpatha(i) = (s2t(i,k1)-s2t(i,k2))/dplh2o(i)
            dtx(i) = tplnka(i,k2) - 250.0D0
            dty(i) = tpatha(i) - 250.0D0
            dtyp15(i) = dty(i) + 15.0D0
            dtyp15sq(i) = dtyp15(i)**d_two
            dtz(i) = dtx(i) - 50.0D0
            dtp(i) = dty(i) - 50.0D0
          end do
          do iband = 2 , 4 , 2
            do i = 1 , iym1
              term1(i,iband) = coefe(1,iband) + &
                               coefe(2,iband)*dtx(i)*(d_one+c1(iband)*dtx(i))
              term2(i,iband) = coefb(1,iband) + &
                               coefb(2,iband)*dtx(i)*(d_one+c2(iband)*dtx(i) * &
                               (d_one+c3(iband)*dtx(i)))
              term3(i,iband) = coefd(1,iband) + &
                               coefd(2,iband)*dtx(i)*(d_one+c4(iband)*dtx(i) * &
                               (d_one+c5(iband)*dtx(i)))
              term4(i,iband) = coefa(1,iband) + &
                               coefa(2,iband)*dty(i)*(d_one+c6(iband)*dty(i))
              term5(i,iband) = coefc(1,iband) + &
                               coefc(2,iband)*dty(i)*(d_one+c7(iband)*dty(i))
            end do
          end do
!
!         abso(i,1)     0 -  800 cm-1   h2o rotation band
!
          do i = 1 , iym1
            a11 = 0.44D0 + 3.380D-4*dtz(i) - 1.520D-6*dtz(i)*dtz(i)
            a31 = 1.05D0 - 6.000D-3*dtp(i) + 3.000D-6*dtp(i)*dtp(i)
            a21 = 1.00D0 + 1.717D-3*dtz(i) - 1.133D-5*dtz(i)*dtz(i)
            a22 = 1.00D0 + 4.443D-3*dtp(i) + 2.750D-5*dtp(i)*dtp(i)
            a23 = 1.00D0 + 3.600D0*sqrtu(i)
            corfac = a31*(a11+((d_two*a21*a22)/a23))
            t1t4 = term1(i,2)*term4(i,2)
            t2t5 = term2(i,2)*term5(i,2)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu(i)*corfac)
            fwk(i) = fwcoef + fwc1/(d_one+fwc2*u(i))
            fwku(i) = fwk(i)*u(i)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,1) = (d_one-rsum)*term3(i,2)
!           trab1(i)  = rsum
          end do
!
!         abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
          do i = 1 , iym1
            a41 = 1.75D0 - 3.960D-03*dtz(i)
            a51 = 1.00D0 + 1.3D0*sqrtu(i)
            a61 = 1.00D0 + 1.250D-03*dtp(i) + 6.250D-05*dtp(i)*dtp(i)
            corfac = 0.29D0*(d_one+a41/a51)*a61
            t1t4 = term1(i,4)*term4(i,4)
            t2t5 = term2(i,4)*term5(i,4)
            a = t1t4 + t2t5/(d_one+t2t5*sqrtu(i)*corfac)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,2) = (d_one-rsum)*term3(i,4)
!           trab7(i)  = rsum
          end do
!
!         Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
          do k = 1 , 2
            do i = 1 , iym1
              phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
              psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
              ubar = dw(i)*phi*1.66D0*r80257
              pbar = pnew(i)*(psi/phi)
              cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
              g2 = d_one + ubar*d_four*st(k)*cf812/pbar
              g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-d_one)
              trline(i,k) = dexp(-g4)
            end do
          end do
          do i = 1 , iym1
            term7(i,1) = coefj(1,1)+coefj(2,1)*dty(i)*(d_one+c16*dty(i))
            term8(i,1) = coefk(1,1)+coefk(2,1)*dty(i)*(d_one+c17*dty(i))
            term7(i,2) = coefj(1,2)+coefj(2,2)*dty(i)*(d_one+c26*dty(i))
            term8(i,2) = coefk(1,2)+coefk(2,2)*dty(i)*(d_one+c27*dty(i))
          end do
!
!         abso(i,3)   800 - 1200 cm-1   h2o window
!         abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
          do i = 1 , iym1
            k21 = term7(i,1) + term8(i,1)/(d_one+(c30+c31*(dty(i)-d_10)* &
                  (dty(i)-d_10))*sqrtu(i))
            k22 = term7(i,2) + term8(i,2)/(d_one+(c28+c29*(dty(i)-d_10))*sqrtu(i))
            tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
            tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
            tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
            tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i) = tr1*tr5
            tr10(i) = tr2*tr6
            th2o(i) = tr10(i)
            trab2(i) = 0.65D0*tr9(i) + 0.35D0*tr10(i)
            trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
            trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
            abso(i,3) = term6(i,k2)*(d_one-trab4(i)*d_half*trline(i,2)- &
                        trab6(i)*d_half*trline(i,1))
            abso(i,4) = term9(i,k2)*d_half*(tr1-tr9(i)+tr2-tr10(i))
          end do
          if ( k2 < k1 ) then
            do i = 1 , iym1
              to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
            end do
          else
            do i = 1 , iym1
              to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
            end do
          end if
!
!           abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
          do i = 1 , iym1
            dpnm(i) = pnm(i,k1) - pnm(i,k2)
            to3co2(i) = (pnm(i,k1)*co2t(i,k1)-pnm(i,k2)*co2t(i,k2))/dpnm(i)
            te = (to3co2(i)*r293)**0.7D0
            dplos = plos(i,k1) - plos(i,k2)
            dplol = plol(i,k1) - plol(i,k2)
            u1 = 18.29D0*dabs(dplos)/te
            u2 = 0.5649D0*dabs(dplos)/te
            rphat = dplol/dplos
            tlocal = tint(i,k2)
            tcrfac = dsqrt(tlocal*r250)*te
            beta = r3205*(rphat+dpfo3*tcrfac)
            realnu = te/beta
            tmp1 = u1/dsqrt(d_four+u1*(d_one+realnu))
            tmp2 = u2/dsqrt(d_four+u2*(d_one+realnu))
            o3bndi = 74.0D0*te*dlog(d_one+tmp1+tmp2)
            abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
            to3(i) = d_one/(d_one+0.1D0*tmp1+0.1D0*tmp2)
!           trab5(i)  = d_one-(o3bndi/(1060-980.))
          end do
!
!         abso(i,6)      co2 15  micrometer band system
!
          do i = 1 , iym1
            sqwp = dsqrt(dabs(plco2(i,k1)-plco2(i,k2)))
            et = dexp(-480.0D0/to3co2(i))
            sqti(i) = dsqrt(to3co2(i))
            rsqti = d_one/sqti(i)
            et2 = et*et
            et4 = et2*et2
            omet = d_one - 1.5D0*et2
            f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
            f1sqwp(i) = f1co2*sqwp
            t1co2(i) = d_one/(d_one+(245.18D0*omet*sqwp*rsqti))
            oneme = d_one - et2
            alphat = oneme**3.0D0*rsqti
            pi = dabs(dpnm(i))
            wco2 = 2.5221D0*co2vmr*pi*regravgts
            u7(i) = 4.9411D4*alphat*et2*wco2
            u8 = 3.9744D4*alphat*et4*wco2
            u9 = 1.0447D5*alphat*et4*et2*wco2
            u13 = 2.8388D3*alphat*et4*wco2
            tpath = to3co2(i)
            tlocal = tint(i,k2)
            tcrfac = dsqrt(tlocal*r250*tpath*r300)
            posqt = ((pnm(i,k2)+pnm(i,k1))*r2sslp+dpfco2*tcrfac)*rsqti
            rbeta7(i) = d_one/(5.3228D0*posqt)
            rbeta8 = d_one/(10.6576D0*posqt)
            rbeta9 = rbeta7(i)
            rbeta13 = rbeta9
            f2co2(i) = (u7(i)/dsqrt(d_four+u7(i)*(d_one+rbeta7(i)))) + &
                       (u8/dsqrt(d_four+u8*(d_one+rbeta8)))          + &
                       (u9/dsqrt(d_four+u9*(d_one+rbeta9)))
            f3co2(i) = u13/dsqrt(d_four+u13*(d_one+rbeta13))
          end do
          if ( k2 >= k1 ) then
            do i = 1 , iym1
              sqti(i) = dsqrt(tlayr(i,k2))
            end do
          end if
!
          do i = 1 , iym1
            tmp1 = dlog(d_one+f1sqwp(i))
            tmp2 = dlog(d_one+f2co2(i))
            tmp3 = dlog(d_one+f3co2(i))
            absbnd = (tmp1+d_two*t1co2(i)*tmp2+d_two*tmp3)*sqti(i)
            abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
            tco2(i) = d_one/(d_one+d_10*(u7(i)/dsqrt(d_four+u7(i)*(d_one+rbeta7(i)))))
!           trab3(i)  = 1. - bndfct*absbnd
          end do
!
!         Calculate absorptivity due to trace gases
!
          call trcab(k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,     &
                     uco212,uco213,uco221,uco222,uco223,bn2o0,bn2o1,  &
                     bch4,to3co2,pnm,dw,pnew,s2c,uptype,u,abplnk1,    &
                     tco2,th2o,to3,abstrc)
!
!         Sum total absorptivity
!
          do i = 1 , iym1
            abstot(i,k1,k2,jslc) = abso(i,1) + abso(i,2) + abso(i,3) + &
                                   abso(i,4) + abso(i,5) + abso(i,6) + abstrc(i)
          end do
        end if
      end do
    end do  ! End of non-nearest layer level loops
!
!   Non-adjacent layer absorptivity:
!
!   abso(i,1)     0 -  800 cm-1   h2o rotation band
!   abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!   abso(i,3)   800 - 1200 cm-1   h2o window
!   abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!   abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!   abso(i,6)   co2 15  micrometer band system
!
!   Nearest layer level loop
!
    do k2 = kz , 1 , -1
      do i = 1 , iym1
        tbar(i,1) = (tint(i,k2+1)+tlayr(i,k2+1))*d_half
        emm(i,1) = (co2em(i,k2+1)+co2eml(i,k2))*d_half
        tbar(i,2) = (tlayr(i,k2+1)+tint(i,k2))*d_half
        emm(i,2) = (co2em(i,k2)+co2eml(i,k2))*d_half
        tbar(i,3) = (tbar(i,2)+tbar(i,1))*d_half
        emm(i,3) = emm(i,1)
        tbar(i,4) = tbar(i,3)
        emm(i,4) = emm(i,2)
        o3emm(i,1) = (dbvtit(i,k2+1)+dbvtly(i,k2))*d_half
        o3emm(i,2) = (dbvtit(i,k2)+dbvtly(i,k2))*d_half
        o3emm(i,3) = o3emm(i,1)
        o3emm(i,4) = o3emm(i,2)
        temh2o(i,1) = tbar(i,1)
        temh2o(i,2) = tbar(i,2)
        temh2o(i,3) = tbar(i,1)
        temh2o(i,4) = tbar(i,2)
        dpnm(i) = pnm(i,k2+1) - pnm(i,k2)
      end do
!
!     Weighted Planck functions for trace gases
!
      do wvl = 1 , 14
        do i = 1 , iym1
          bplnk(wvl,i,1) = (abplnk1(wvl,i,k2+1)+abplnk2(wvl,i,k2))*d_half
          bplnk(wvl,i,2) = (abplnk1(wvl,i,k2)+abplnk2(wvl,i,k2))*d_half
          bplnk(wvl,i,3) = bplnk(wvl,i,1)
          bplnk(wvl,i,4) = bplnk(wvl,i,2)
        end do
      end do
!
      do i = 1 , iym1
        rdpnmsq = d_one/(pnmsq(i,k2+1)-pnmsq(i,k2))
        rdpnm = d_one/dpnm(i)
        p1 = (pbr(i,k2)+pnm(i,k2+1))*d_half
        p2 = (pbr(i,k2)+pnm(i,k2))*d_half
        uinpl(i,1) = (pnmsq(i,k2+1)-p1**d_two)*rdpnmsq
        uinpl(i,2) = -(pnmsq(i,k2)-p2**d_two)*rdpnmsq
        uinpl(i,3) = -(pnmsq(i,k2)-p1**d_two)*rdpnmsq
        uinpl(i,4) = (pnmsq(i,k2+1)-p2**d_two)*rdpnmsq
        winpl(i,1) = ((pnm(i,k2+1)-pbr(i,k2))*d_half)*rdpnm
        winpl(i,2) = ((-pnm(i,k2)+pbr(i,k2))*d_half)*rdpnm
        winpl(i,3) = ((pnm(i,k2+1)+pbr(i,k2))*d_half-pnm(i,k2))*rdpnm
        winpl(i,4) = ((-pnm(i,k2)-pbr(i,k2))*d_half+pnm(i,k2+1))*rdpnm
        tmp1 = d_one/(piln(i,k2+1)-piln(i,k2))
        tmp2 = piln(i,k2+1) - pmln(i,k2)
        tmp3 = piln(i,k2) - pmln(i,k2)
        zinpl(i,1) = (tmp2*d_half)*tmp1
        zinpl(i,2) = (-tmp3*d_half)*tmp1
        zinpl(i,3) = (tmp2*d_half-tmp3)*tmp1
        zinpl(i,4) = (tmp2-tmp3*d_half)*tmp1
        pinpl(i,1) = (p1+pnm(i,k2+1))*d_half
        pinpl(i,2) = (p2+pnm(i,k2))*d_half
        pinpl(i,3) = (p1+pnm(i,k2))*d_half
        pinpl(i,4) = (p2+pnm(i,k2+1))*d_half
      end do

!     FAB AER SAVE uinpl  for aerosl LW forcing calculation
      xuinpl (:,k2,:,jslc) =  uinpl(:,:)
!     FAB AER SAVE uinpl  for aerosl LW forcing calculation

      do kn = 1 , 4
        do i = 1 , iym1
          u(i) = uinpl(i,kn)*dabs(plh2o(i,k2)-plh2o(i,k2+1))
          sqrtu(i) = dsqrt(u(i))
          dw(i) = dabs(w(i,k2)-w(i,k2+1))
          pnew(i) = u(i)/(winpl(i,kn)*dw(i))
          ds2c = dabs(s2c(i,k2)-s2c(i,k2+1))
          uc1(i) = uinpl(i,kn)*ds2c
          uc1(i) = (uc1(i)+1.7D-3*u(i))*(d_one+d_two*uc1(i))/(d_one+15.0D0*uc1(i))
          uc(i) = uinpl(i,kn)*ds2c + 2.0D-3*u(i)
        end do
        do i = 1 , iym1
          dtx(i) = temh2o(i,kn) - 250.0D0
          dty(i) = tbar(i,kn) - 250.0D0
          dtyp15(i) = dty(i) + 15.0D0
          dtyp15sq(i) = dtyp15(i)**d_two
          dtz(i) = dtx(i) - 50.0D0
          dtp(i) = dty(i) - 50.0D0
        end do
        do iband = 2 , 4 , 2
          do i = 1 , iym1
            term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i) * &
                             (d_one+c1(iband)*dtx(i))
            term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i) * &
                             (d_one+c2(iband)*dtx(i)                * &
                             (d_one+c3(iband)*dtx(i)))
            term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i) * &
                             (d_one+c4(iband)*dtx(i)                * &
                             (d_one+c5(iband)*dtx(i)))
            term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i) * &
                             (d_one+c6(iband)*dty(i))
            term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i) * &
                             (d_one+c7(iband)*dty(i))
          end do
        end do
!
!       abso(i,1)     0 -  800 cm-1   h2o rotation band
!
        do i = 1 , iym1
          a11 = 0.44D0 + 3.380D-4*dtz(i) - 1.520D-6*dtz(i)*dtz(i)
          a31 = 1.05D0 - 6.000D-3*dtp(i) + 3.000D-6*dtp(i)*dtp(i)
          a21 = 1.00D0 + 1.717D-3*dtz(i) - 1.133D-5*dtz(i)*dtz(i)
          a22 = 1.00D0 + 4.443D-3*dtp(i) + 2.750D-5*dtp(i)*dtp(i)
          a23 = 1.00D0 + 3.600D0*sqrtu(i)
          corfac = a31*(a11+((d_two*a21*a22)/a23))
          t1t4 = term1(i,2)*term4(i,2)
          t2t5 = term2(i,2)*term5(i,2)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu(i)*corfac)
          fwk(i) = fwcoef + fwc1/(d_one+fwc2*u(i))
          fwku(i) = fwk(i)*u(i)
          rsum = dexp(-a*(sqrtu(i)+fwku(i)))
          abso(i,1) = (d_one-rsum)*term3(i,2)
!         trab1(i) = rsum
        end do
!
!       abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
        do i = 1 , iym1
          a41 = 1.75D0 - 3.960D-03*dtz(i)
          a51 = 1.00D0 + 1.3D0*sqrtu(i)
          a61 = 1.00D0 + 1.250D-03*dtp(i) + 6.250D-05*dtp(i)*dtp(i)
          corfac = 0.29D0*(d_one+a41/a51)*a61
          t1t4 = term1(i,4)*term4(i,4)
          t2t5 = term2(i,4)*term5(i,4)
          a = t1t4 + t2t5/(d_one+t2t5*sqrtu(i)*corfac)
          rsum = dexp(-a*(sqrtu(i)+fwku(i)))
          abso(i,2) = (d_one-rsum)*term3(i,4)
!         trab7(i) = rsum
        end do
!
!       Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
        do k = 1 , 2
          do i = 1 , iym1
            phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
            psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
            ubar = dw(i)*phi*winpl(i,kn)*1.66D0*r80257
            pbar = pnew(i)*(psi/phi)
            cf812 = cfa1 + (d_one-cfa1)/(d_one+ubar*pbar*d_10)
            g2 = d_one + ubar*d_four*st(k)*cf812/pbar
            g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-d_one)
            trline(i,k) = dexp(-g4)
          end do
        end do
        do i = 1 , iym1
          term7(i,1) = coefj(1,1)+coefj(2,1)*dty(i)*(d_one+c16*dty(i))
          term8(i,1) = coefk(1,1)+coefk(2,1)*dty(i)*(d_one+c17*dty(i))
          term7(i,2) = coefj(1,2)+coefj(2,2)*dty(i)*(d_one+c26*dty(i))
          term8(i,2) = coefk(1,2)+coefk(2,2)*dty(i)*(d_one+c27*dty(i))
        end do
!
!       abso(i,3)   800 - 1200 cm-1   h2o window
!       abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!
        do i = 1 , iym1
          dtym10 = dty(i) - d_10
          denom = d_one + (c30+c31*dtym10*dtym10)*sqrtu(i)
          k21 = term7(i,1) + term8(i,1)/denom
          denom = d_one + (c28+c29*dtym10)*sqrtu(i)
          k22 = term7(i,2) + term8(i,2)/denom
          term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i) *        &
                       (d_one+c19*dtx(i)*(d_one+c21*dtx(i) *   &
                       (d_one+c23*dtx(i)*(d_one+c25*dtx(i)))))
          tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
          tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
          tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
          tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
          tr9(i) = tr1*tr5
          tr10(i) = tr2*tr6
          trab2(i) = 0.65D0*tr9(i) + 0.35D0*tr10(i)
          th2o(i) = tr10(i)
          trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
          trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
          term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i) *            &
                       (d_one+c9*dtx(i)*(d_one+c11*dtx(i) *        &
                       (d_one+c13*dtx(i)*(d_one+c15*dtx(i)))))
          abso(i,3) = term6(i,2) *                                 &
                      (d_one-trab4(i)*d_half*trline(i,2)-          &
                             trab6(i)*d_half*trline(i,1))
          abso(i,4) = term9(i,2)*d_half*(tr1-tr9(i)+tr2-tr10(i))
        end do
!
!       abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
!
        do i = 1 , iym1
          te = (tbar(i,kn)*r293)**0.7D0
          dplos = dabs(plos(i,k2+1)-plos(i,k2))
          u1 = zinpl(i,kn)*18.29D0*dplos/te
          u2 = zinpl(i,kn)*0.5649D0*dplos/te
          tlocal = tbar(i,kn)
          tcrfac = dsqrt(tlocal*r250)*te
          beta = r3205*(pinpl(i,kn)*rsslp+dpfo3*tcrfac)
          realnu = te/beta
          tmp1 = u1/dsqrt(d_four+u1*(d_one+realnu))
          tmp2 = u2/dsqrt(d_four+u2*(d_one+realnu))
          o3bndi = 74.0D0*te*dlog(d_one+tmp1+tmp2)
          abso(i,5) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
          to3(i) = d_one/(d_one+0.1D0*tmp1+0.1D0*tmp2)
!         trab5(i) = d_one-(o3bndi/(1060-980.))
        end do
!
!       abso(i,6)   co2 15  micrometer band system
!
        do i = 1 , iym1
          dplco2 = plco2(i,k2+1) - plco2(i,k2)
          sqwp = dsqrt(uinpl(i,kn)*dplco2)
          et = dexp(-480.0D0/tbar(i,kn))
          sqti(i) = dsqrt(tbar(i,kn))
          rsqti = d_one/sqti(i)
          et2 = et*et
          et4 = et2*et2
          omet = (d_one-1.5D0*et2)
          f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
          f1sqwp(i) = f1co2*sqwp
          t1co2(i) = d_one/(d_one+(245.18D0*omet*sqwp*rsqti))
          oneme = d_one - et2
          alphat = oneme**3.0D0*rsqti
          pi = dabs(dpnm(i))*winpl(i,kn)
          wco2 = 2.5221D0*co2vmr*pi*regravgts
          u7(i) = 4.9411D4*alphat*et2*wco2
          u8 = 3.9744D4*alphat*et4*wco2
          u9 = 1.0447D5*alphat*et4*et2*wco2
          u13 = 2.8388D3*alphat*et4*wco2
          tpath = tbar(i,kn)
          tlocal = tbar(i,kn)
          tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
          posqt = (pinpl(i,kn)*rsslp+dpfco2*tcrfac)*rsqti
          rbeta7(i) = d_one/(5.3228D0*posqt)
          rbeta8 = d_one/(10.6576D0*posqt)
          rbeta9 = rbeta7(i)
          rbeta13 = rbeta9
          f2co2(i) = u7(i)/dsqrt(d_four+u7(i)*(d_one+rbeta7(i))) + &
                     u8/dsqrt(d_four+u8*(d_one+rbeta8)) +          &
                     u9/dsqrt(d_four+u9*(d_one+rbeta9))
          f3co2(i) = u13/dsqrt(d_four+u13*(d_one+rbeta13))
          tmp1 = dlog(d_one+f1sqwp(i))
          tmp2 = dlog(d_one+f2co2(i))
          tmp3 = dlog(d_one+f3co2(i))
          absbnd = (tmp1+d_two*t1co2(i)*tmp2+d_two*tmp3)*sqti(i)
          abso(i,6) = trab2(i)*emm(i,kn)*absbnd
          tco2(i) = d_one/(d_one+d_10*u7(i)/dsqrt(d_four+u7(i)*(d_one+rbeta7(i))))
!         trab3(i) = 1. - bndfct*absbnd
        end do
!
!       Calculate trace gas absorptivity for nearest layer
!
        call trcabn(k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,      &
                    uco212,uco213,uco221,uco222,uco223,tbar,bplnk,    &
                    winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,u,pnew,   &
                    abstrc,uinpl)
!
!       Total next layer absorptivity:
!
        do i = 1 , iym1
          absnxt(i,k2,kn,jslc) = abso(i,1) + abso(i,2) + abso(i,3)    &
                                 + abso(i,4) + abso(i,5) + abso(i,6)  &
                                 + abstrc(i)
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
  subroutine radems(s2c,s2t,w,tplnke,plh2o,pnm,plco2,tint,tint4,    &
                    tlayr,tlayr4,plol,plos,ucfc11,ucfc12,un2o0,     &
                    un2o1,uch4,uco211,uco212,uco213,uco221,uco222,  &
                    uco223,uptype,bn2o0,bn2o1,bch4,co2em,co2eml,    &
                    co2t,h2otr,abplnk1,abplnk2,jslc)
!
    implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! s2c     - H2o continuum path length
! s2t     - Tmp and prs wghted h2o path length
! w       - H2o path length
! tplnke  - Layer planck temperature
! plh2o   - H2o prs wghted path length
! pnm     - Model interface pressure
! plco2   - Prs wghted path of co2
! tint    - Model interface temperatures
! tint4   - Tint to the 4th power
! tlayr   - K-1 model layer temperature
! tlayr4  - Tlayr to the 4th power
! plol    - Pressure wghtd ozone path
! plos    - Ozone path
!
!     Trace gas variables
!
! ucfc11  - CFC11 path length
! ucfc12  - CFC12 path length
! un2o0   - N2O path length
! un2o1   - N2O path length (hot band)
! uch4    - CH4 path length
! uco211  - CO2 9.4 micron band path length
! uco212  - CO2 9.4 micron band path length
! uco213  - CO2 9.4 micron band path length
! uco221  - CO2 10.4 micron band path length
! uco222  - CO2 10.4 micron band path length
! uco223  - CO2 10.4 micron band path length
! bn2o0   - pressure factor for n2o
! bn2o1   - pressure factor for n2o
! bch4    - pressure factor for ch4
! uptype  - p-type continuum path length
!
!     Output arguments
!
! co2em   - Layer co2 normalzd plnck funct drvtv
! co2eml  - Intrfc co2 normalzd plnck func drvtv
! co2t    - Tmp and prs weighted path length
! h2otr   - H2o transmission over o3 band
! emplnk  - emissivity Planck factor
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
! emstrc  - total trace gas emissivity
!
    integer :: jslc
    real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
    real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em , &
           co2t , h2otr , plco2 , plh2o , plol , plos , pnm , s2c ,  &
           s2t , tint , tint4 , tlayr , tlayr4 , ucfc11 , ucfc12 ,   &
           uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,       &
           uco223 , un2o0 , un2o1 , uptype , w
    real(8) , dimension(iym1,kz) :: co2eml
    real(8) , dimension(iym1) :: tplnke
    intent (in) jslc , plco2 , plh2o , plol , plos , s2t , tint4 , tlayr4
    intent (out) co2em , co2eml
    intent (inout) co2t , h2otr
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
! term6   - B(T) function for window region
! term7   - Kl_inf(i) in eq(8) of table A3a of R&D
! term8   - Delta kl_inf(i) in eq(8)
! term9   - B(T) function for 500-800 cm-1 region
! tr1     - Equation(6) in table A2 for 650-800
! tr2     - Equation(6) in table A2 for 500-650
! tr3     - Equation(4) in table A2 for 650-800
! tr4     - Equation(4),table A2 of R&D for 500-650 
! tr7     - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! tr8     - Equation (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
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
    real(8) , dimension(iym1) :: a , co2plk , corfac , dbvtt , dtp ,  &
                                  dtx , dty , dtz , k21 , k22 , pnew ,&
                                  rsum , tco2 , th2o , to3 , tpathe , &
                                  tr1 , tr2 , tr3 , tr4 , tr7 , tr8 , &
                                  trem4 , trem6 , u , uc , uc1 , xsum
    real(8) :: a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,        &
               absbnd , alphat , beta , cf812 , et , et2 , et4 , ex , &
               exm1sq , f1co2 , f1sqwp , f2co2 , f3co2 , fwk , g1 ,   &
               g2 , g3 , g4 , o3bndi , omet , oneme , pbar , phat ,   &
               phi , pi , posqt , psi , &
               rbeta13 , rbeta7 , rbeta8 , rbeta9 , realnu , rsqti ,  &
               sqti , sqwp , t1co2 , t1i , t1t4 , t2t5 ,              &
               tcrfac , te , tlayr5 , tlocal , tmp1 , tmp2 , tmp3 ,   &
               tpath , u1 , u13 , u2 , u7 , u8 , u9 , ubar , wco2
    real(8) , dimension(iym1,kzp1) :: co2ems , emstrc , h2oems ,    &
           o3ems , troco2
    real(8) , dimension(iym1,4) :: emis , term1 , term2 , term3 ,    &
                                    term4 , term5
    real(8) , dimension(14,iym1) :: emplnk
    integer :: i , iband , k , k1
    real(8) , dimension(iym1,2) :: term6 , term7 , term8 , term9 ,   &
                                    trline
!
    character (len=50) :: subroutine_name='radems'
    integer :: indx = 0
!
!-----------------------------------------------------------------------
    call time_begin(subroutine_name,indx)
!
!   Planck function for co2
!
    do i = 1 , iym1
      ex = dexp(960.0D0/tplnke(i))
      co2plk(i) = 5.0D8/((tplnke(i)**d_four)*(ex-d_one))
      co2t(i,1) = tplnke(i)
      xsum(i) = co2t(i,1)*pnm(i,1)
    end do
    k = 1
    do k1 = kzp1 , 2 , -1
      k = k + 1
      do i = 1 , iym1
        xsum(i) = xsum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
        ex = dexp(960.0D0/tlayr(i,k1))
        tlayr5 = tlayr(i,k1)*tlayr4(i,k1)
        co2eml(i,k1-1) = 1.2D11*ex/(tlayr5*(ex-d_one)**d_two)
        co2t(i,k) = xsum(i)/pnm(i,k)
      end do
    end do
!   bndfct = 2.d0*22.18/(dsqrt(196.d0)*300.)
!
!   Initialize planck function derivative for O3
!
    do i = 1 , iym1
      dbvtt(i) = dbvt(tplnke(i))
    end do
!
!   Calculate trace gas Planck functions
!
    call trcplk(tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
!   Interface loop
!
    do k1 = 1 , kzp1
!
!     H2O emissivity
!
!     emis(i,1)     0 -  800 cm-1   rotation band
!     emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!     emis(i,3)   800 - 1200 cm-1   window
!     emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
!     For the p type continuum
!
      do i = 1 , iym1
        uc(i) = s2c(i,k1) + 2.0D-3*plh2o(i,k1)
        u(i) = plh2o(i,k1)
        pnew(i) = u(i)/w(i,k1)
!
!       Apply scaling factor for 500-800 continuum
!
        uc1(i) = (s2c(i,k1)+1.7D-3*plh2o(i,k1)) * &
                 (d_one+d_two*s2c(i,k1))/(d_one+15.0D0*s2c(i,k1))
        tpathe(i) = s2t(i,k1)/plh2o(i,k1)
      end do
      do i = 1 , iym1
        dtx(i) = tplnke(i) - 250.0D0
        dty(i) = tpathe(i) - 250.0D0
        dtz(i) = dtx(i) - 50.0D0
        dtp(i) = dty(i) - 50.0D0
      end do
      do iband = 1 , 3 , 2
        do i = 1 , iym1
          term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i) *   &
                           (d_one+c1(iband)*dtx(i))
          term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i) *   &
                           (d_one+c2(iband)*dtx(i)*(d_one+c3(iband)*dtx(i)))
          term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i) *   &
                           (d_one+c4(iband)*dtx(i)*(d_one+c5(iband)*dtx(i)))
          term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i) *   &
                           (d_one+c6(iband)*dty(i))
          term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i) *   &
                           (d_one+c7(iband)*dty(i))
        end do
      end do
!
!     emis(i,1)     0 -  800 cm-1   rotation band
!
      do i = 1 , iym1
        a11 = 0.37D0 - 3.33D-5*dtz(i) + 3.33D-6*dtz(i)*dtz(i)
        a31 = 1.07D0 - 1.00D-3*dtp(i) + 1.475D-5*dtp(i)*dtp(i)
        a21 = 1.3870D0 + 3.80D-3*dtz(i) - 7.8D-6*dtz(i)*dtz(i)
        a22 = d_one - 1.21D-3*dtp(i) - 5.33D-6*dtp(i)*dtp(i)
        a23 = 0.9D0 + 2.62D0*dsqrt(u(i))
        corfac(i) = a31*(a11+((a21*a22)/a23))
        t1t4 = term1(i,1)*term4(i,1)
        t2t5 = term2(i,1)*term5(i,1)
        a(i) = t1t4 + t2t5/(d_one+t2t5*dsqrt(u(i))*corfac(i))
        fwk = fwcoef + fwc1/(d_one+fwc2*u(i))
        rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
        emis(i,1) = (d_one-rsum(i))*term3(i,1)
!       trem1(i)  = rsum(i)
!
!       emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!
        a41 = 1.75D0 - 3.96D-3*dtz(i)
        a51 = 1.00D0 + 1.3D0*dsqrt(u(i))
        a61 = 1.00D0 + 1.25D-3*dtp(i) + 6.25D-5*dtp(i)*dtp(i)
        corfac(i) = 0.3D0*(d_one+(a41)/(a51))*a61
        t1t4 = term1(i,3)*term4(i,3)
        t2t5 = term2(i,3)*term5(i,3)
        a(i) = t1t4 + t2t5/(d_one+t2t5*dsqrt(u(i))*corfac(i))
        fwk = fwcoef + fwc1/(d_one+fwc2*u(i))
        rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
        emis(i,2) = (d_one-rsum(i))*term3(i,3)
!       trem7(i) = rsum(i)
      end do
!
!     Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
      do k = 1 , 2
        do i = 1 , iym1
          phi = a1(k)*(dty(i)+15.0D0)+a2(k)*(dty(i)+15.0D0)**d_two
          psi = b1(k)*(dty(i)+15.0D0)+b2(k)*(dty(i)+15.0D0)**d_two
          phi = dexp(phi)
          psi = dexp(psi)
          ubar = w(i,k1)*phi
          ubar = (ubar*1.66D0)*r80257
          pbar = pnew(i)*(psi/phi)
          cf812 = cfa1 + ((d_one-cfa1)/(d_one+ubar*pbar*d_10))
          g1 = (realk(k)*pbar)/(d_two*st(k))
          g2 = d_one + (ubar*d_four*st(k)*cf812)/pbar
          g3 = dsqrt(g2) - d_one
          g4 = g1*g3
          trline(i,k) = dexp(-g4)
        end do
      end do
      do i = 1 , iym1
        term7(i,1) = coefj(1,1)+coefj(2,1)*dty(i)*(d_one+c16*dty(i))
        term8(i,1) = coefk(1,1)+coefk(2,1)*dty(i)*(d_one+c17*dty(i))
        term7(i,2) = coefj(1,2)+coefj(2,2)*dty(i)*(d_one+c26*dty(i))
        term8(i,2) = coefk(1,2)+coefk(2,2)*dty(i)*(d_one+c27*dty(i))
      end do
!
!     emis(i,3)   800 - 1200 cm-1   window
!
      do i = 1 , iym1
        term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i) *        &
                     (d_one+c8*dtx(i)*(d_one+c10*dtx(i) *    &
                     (d_one+c12*dtx(i)*(d_one+c14*dtx(i)))))
        trem4(i) = dexp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))*trline(i,2)
        trem6(i) = dexp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))*trline(i,1)
        emis(i,3) = term6(i,1)*(d_one-trem4(i)*d_half-trem6(i)*d_half)
!
!       emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
        k21(i) = term7(i,1) + term8(i,1)/(d_one+(c30+c31*(dty(i)-d_10) * &
                 (dty(i)-d_10))*dsqrt(u(i)))
        k22(i) = term7(i,2) + term8(i,2)/(d_one+(c28+c29*(dty(i)-d_10))*dsqrt(u(i)))
        term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i) *        &
                     (d_one+c18*dtx(i)*(d_one+c20*dtx(i) *   &
                     (d_one+c22*dtx(i)*(d_one+c24*dtx(i)))))
        fwk = fwcoef + fwc1/(d_one+fwc2*u(i))
        tr1(i) = dexp(-(k21(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
        tr2(i) = dexp(-(k22(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
        tr3(i) = dexp(-((coefh(1,1)+coefh(2,1)*dtx(i))*uc1(i)))
        tr4(i) = dexp(-((coefh(1,2)+coefh(2,2)*dtx(i))*uc1(i)))
        tr7(i) = tr1(i)*tr3(i)
        tr8(i) = tr2(i)*tr4(i)
        emis(i,4) = term9(i,1)*d_half*(tr1(i)-tr7(i)+tr2(i)-tr8(i))
        h2oems(i,k1) = emis(i,1) + emis(i,2) + emis(i,3) + emis(i,4)
        troco2(i,k1) = 0.65D0*tr7(i) + 0.35D0*tr8(i)
        th2o(i) = tr8(i)
!       trem2(i)     = troco2(i,k1)
      end do
!
!     CO2 emissivity for 15 micron band system
!
      do i = 1 , iym1
        t1i = dexp(-480.0D0/co2t(i,k1))
        sqti = dsqrt(co2t(i,k1))
        rsqti = d_one/sqti
        et = t1i
        et2 = et*et
        et4 = et2*et2
        omet = d_one - 1.5D0*et2
        f1co2 = 899.70D0*omet*(d_one+1.94774D0*et+4.73486D0*et2)*rsqti
        sqwp = dsqrt(plco2(i,k1))
        f1sqwp = f1co2*sqwp
        t1co2 = d_one/(d_one+245.18D0*omet*sqwp*rsqti)
        oneme = d_one - et2
        alphat = oneme**3.0D0*rsqti
        wco2 = 2.5221D0*co2vmr*pnm(i,k1)*regravgts
        u7 = 4.9411D4*alphat*et2*wco2
        u8 = 3.9744D4*alphat*et4*wco2
        u9 = 1.0447D5*alphat*et4*et2*wco2
        u13 = 2.8388D3*alphat*et4*wco2
!
        tpath = co2t(i,k1)
        tlocal = tplnke(i)
        tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
        pi = pnm(i,k1)*rsslp + d_two*dpfco2*tcrfac
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
        tco2(i) = d_one/(d_one+d_10*(u7/dsqrt(d_four+u7*(d_one+rbeta7))))
        co2ems(i,k1) = troco2(i,k1)*absbnd*co2plk(i)
        ex = dexp(960.0D0/tint(i,k1))
        exm1sq = (ex-d_one)**d_two
        co2em(i,k1) = 1.2D11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
!       trem3(i) = 1. - bndfct*absbnd
      end do
!
!     O3 emissivity
!
      do i = 1 , iym1
        h2otr(i,k1) = dexp(-12.0D0*s2c(i,k1))
        te = (co2t(i,k1)/293.0D0)**0.7D0
        u1 = 18.29D0*plos(i,k1)/te
        u2 = 0.5649D0*plos(i,k1)/te
        phat = plos(i,k1)/plol(i,k1)
        tlocal = tplnke(i)
        tcrfac = dsqrt(tlocal*r250)*te
        beta = (d_one/0.3205D0)*((d_one/phat)+(dpfo3*tcrfac))
        realnu = (d_one/beta)*te
        o3bndi = 74.0D0*te*(tplnke(i)/375.0D0)*dlog(d_one+fo3(u1,realnu) + &
                 fo3(u2,realnu))
        o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
        to3(i) = d_one/(d_one+0.1D0*fo3(u1,realnu)+0.1D0*fo3(u2,realnu))
!       trem5(i)    = d_one-(o3bndi/(1060-980.))
      end do
!
!     Calculate trace gas emissivities
!
      call trcems(k1,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,bn2o0,bn2o1,   &
                  uch4,bch4,uco211,uco212,uco213,uco221,uco222,uco223, &
                  uptype,w,s2c,u,emplnk,th2o,tco2,to3,emstrc)
!
!     Total emissivity:
!
      do i = 1 , iym1
        emstot(i,k1,jslc) = h2oems(i,k1)+co2ems(i,k1)+o3ems(i,k1)+emstrc(i,k1)
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
  subroutine radoz2(o3vmr,pint,plol,plos)
!
    implicit none
!
!------------------------------Input arguments--------------------------
!
! o3vmr   - ozone volume mixing ratio
! pint    - Model interface pressures
!
!----------------------------Output arguments---------------------------
!
! plol    - Ozone prs weighted path length (cm)
! plos    - Ozone path length (cm)
!
!
    real(8) , dimension(iym1,kz) :: o3vmr
    real(8) , dimension(iym1,kzp1) :: pint , plol , plos
    intent (in) o3vmr , pint
    intent (inout) plol , plos
!
!-----------------------------------------------------------------------
!
    integer :: i , k
    character (len=50) :: subroutine_name='radoz2'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
!   Evaluate the ozone path length integrals to interfaces;
!   factors of 0.1 and 0.01 to convert pressures from cgs to mks:
!
!   Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
!
    do i = 1 , iym1
      plos(i,1) = 0.1D0*cplos*o3vmr(i,1)*pint(i,1)
      plol(i,1) = 0.01D0*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
    end do
    do k = 2 , kzp1
      do i = 1 , iym1
        plos(i,k) = plos(i,k-1) + 0.1D0*cplos*o3vmr(i,k-1)*(pint(i,k)-pint(i,k-1))
        plol(i,k) = plol(i,k-1) + 0.01D0*cplol*o3vmr(i,k-1) * &
                    (pint(i,k)*pint(i,k)-pint(i,k-1)*pint(i,k-1))
      end do
    end do
!
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
  subroutine radtpl(tnm,ts,qnm,pnm,plh2o,tplnka,s2c,s2t,w,tplnke,   &
                    tint,tint4,tlayr,tlayr4,pmln,piln)
!
    implicit none
!
!     Input arguments
!
! tnm    - Model level temperatures
! ts     - Surface skin temperature
! qnm    - Model level specific humidity
! pnm    - Pressure at model interfaces (dynes/cm2)
! plh2o  - Pressure weighted h2o path
!
!     Output arguments
!
! tplnka - Level temperature from interface temperatures
! s2c    - H2o continuum path length
! s2t    - H2o tmp and prs wghtd path length
! w      - H2o path length
! tplnke - Equal to tplnka
! tint   - Layer interface temperature
! tint4  - Tint to the 4th power
! tlayr  - K-1 level temperature
! tlayr4 - Tlayr to the 4th power
! pmln   - Ln(pmidm1)
! piln   - Ln(pintm1)
!
!
    real(8) , dimension(iym1,kzp1) :: piln , plh2o , pnm , s2c ,    &
           s2t , tint , tint4 , tlayr , tlayr4 , tplnka , w
    real(8) , dimension(iym1,kz) :: pmln , qnm , tnm
    real(8) , dimension(iym1) :: tplnke , ts
    intent (in) piln , plh2o , pmln , pnm , qnm , tnm , ts
    intent (out) tint4 , tplnke
    intent (inout) s2c , s2t , tint , tlayr , tlayr4 , tplnka , w
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
    real(8) :: dpnm , dpnmsq , dy , rtnm
    integer :: i , k
    character (len=50) :: subroutine_name='radtpl'
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
    do i = 1 , iym1
      tint(i,kzp1) = ts(i)
      tint4(i,kzp1) = tint(i,kz + 1)**d_four
      tplnka(i,1) = tnm(i,1)
      tint(i,1) = tplnka(i,1)
      tlayr4(i,1) = tplnka(i,1)**d_four
      tint4(i,1) = tlayr4(i,1)
    end do
!
!   Intermediate level temperatures are computed using temperature
!   at the full level below less dy*delta t,between the full level
!
    do k = 2 , kz
      do i = 1 , iym1
        dy = (piln(i,k)-pmln(i,k))/(pmln(i,k-1)-pmln(i,k))
        tint(i,k) = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
        tint4(i,k) = tint(i,k)**d_four
      end do
    end do
!
!   Now set the layer temp=full level temperatures and establish a
!   planck temperature for absorption (tplnka) which is the average
!   the intermediate level temperatures.  Note that tplnka is not
!   equal to the full level temperatures.
!
    do k = 2 , kzp1
      do i = 1 , iym1
        tlayr(i,k) = tnm(i,k-1)
        tlayr4(i,k) = tlayr(i,k)**d_four
        tplnka(i,k) = (tint(i,k)+tint(i,k-1))*d_half
      end do
    end do
!
!   Calculate tplank for emissivity calculation.
!   Assume isothermal tplnke i.e. all levels=ttop.
!
    do i = 1 , iym1
      tplnke(i) = tplnka(i,1)
      tlayr(i,1) = tint(i,1)
    end do
!
!   Now compute h2o path fields:
!
    do i = 1 , iym1
      s2t(i,1) = plh2o(i,1)*tnm(i,1)
!     ccm3.2
!     w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
!     s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil
   
!     ccm3.6.6
      w(i,1) = sslp*(plh2o(i,1)*d_two)/pnm(i,1)
      rtnm = d_one/tnm(i,1)
      s2c(i,1) = plh2o(i,1)*dexp(1800.0D0*(rtnm-r296))*qnm(i,1)*repsil
    end do
    do k = 1 , kz
      do i = 1 , iym1
        dpnm = pnm(i,k+1) - pnm(i,k)
        dpnmsq = pnm(i,k+1)**d_two - pnm(i,k)**d_two
        rtnm = d_one/tnm(i,k)
        s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
        w(i,k+1) = w(i,k) + regravgts*qnm(i,k)*dpnm
        s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)                &
                     *dexp(1800.0D0*(rtnm-r296))*qnm(i,k)*repsil
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
  subroutine radinp(pmid,pint,h2ommr,cld,o3vmr,pmidrd,pintrd,plco2, &
                    plh2o,tclrsf,eccf,o3mmr)
!
    implicit none
!
    real(8) :: eccf
    real(8) , dimension(iym1,kzp1) :: cld , pint , pintrd , plco2 ,   &
           plh2o , tclrsf
    real(8) , dimension(iym1,kz) :: h2ommr , o3mmr , o3vmr , pmid ,   &
           pmidrd
    intent (in) cld , h2ommr , o3vmr , pint , pmid
    intent (out) eccf , o3mmr , plco2 , pmidrd
    intent (inout) pintrd , plh2o , tclrsf
!
#ifndef CLM
    real(8) :: theta
#endif
    real(8) :: cpwpl , vmmr
    integer :: i , k
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
! eccf    - Earth-sun distance factor
! o3mmr   - Ozone mass mixing ratio
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude loop index
! k       - Vertical loop index
! theta   - Earth orbit seasonal angle in radians
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
!   Compute eccentricity factor (sun-earth distance factor)
!
    character (len=50) :: subroutine_name='radinp'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
#ifdef CLM
    eccf  = r2ceccf
#else
    theta = twopi*calday/dayspy
    eccf = 1.000110D0 + 0.034221D0*dcos(theta) +  &
           0.001280D0 * dsin(theta) + &
           0.000719D0 * dcos(d_two*theta) + &
           0.000077D0 * dsin(d_two*theta)
#endif
!
!   Convert pressure from pascals to dynes/cm2
!
    do k = 1 , kz
      do i = 1 , iym1
        pmidrd(i,k) = pmid(i,k)*d_10
        pintrd(i,k) = pint(i,k)*d_10
      end do
    end do
    do i = 1 , iym1
      pintrd(i,kzp1) = pint(i,kz + 1)*d_10
    end do
!
!   Compute path quantities used in the longwave radiation:
!
    vmmr = amco2/amd
    cpwpl = vmmr*d_half/(egravgts*sslp)
    do i = 1 , iym1
      plh2o(i,1) = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
      plco2(i,1) = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
      tclrsf(i,1) = d_one
    end do
    do k = 1 , kz
      do i = 1 , iym1
        plh2o(i,k+1) = plh2o(i,k) + rgsslp*(pintrd(i,k+1)**d_two - &
                       pintrd(i,k)**d_two) * h2ommr(i,k)
        plco2(i,k+1) = co2vmr*cpwpl*pintrd(i,k+1)**d_two
        tclrsf(i,k+1) = tclrsf(i,k)*(d_one-cld(i,k+1))
      end do
    end do
!
!   Convert ozone volume mixing ratio to mass mixing ratio:
!
    vmmr = amo/amd
    do k = 1 , kz
      do i = 1 , iym1
        o3mmr(i,k) = vmmr*o3vmr(i,k)
      end do
    end do
!
    call time_end(subroutine_name,indx)
  end subroutine radinp
!
  function isrchfgt(n,array,inc,rtarg)
    implicit none
!
    integer :: inc , n
    real(8) :: rtarg
    real(8) , dimension(*) :: array
    integer :: isrchfgt
    intent (in) array , inc , n , rtarg
!
    integer :: i
!
    if ( n <= 0 ) then
      isrchfgt = 0
      return
    end if
    isrchfgt = 1
    do i = 1 , n , inc
      if ( array(i) > rtarg ) then
        return
      else
        isrchfgt = isrchfgt + inc
      end if
    end do
  end function isrchfgt
 
  function isrchfle(n,array,inc,rtarg)
    implicit none
!
    integer :: inc , n
    real(8) :: rtarg
    real(8) , dimension(*) :: array
    integer :: isrchfle
    intent (in) array , inc , n , rtarg
!
    integer :: i
!
    if ( n <= 0 ) then
      isrchfle = 0
      return
    end if
    isrchfle = 1
    do i = 1 , n , inc
      if ( array(i) <= rtarg ) then
        return
      else
        isrchfle = isrchfle + inc
      end if
    end do
  end function isrchfle
!
  subroutine wheneq(n,iarray,inc,itarg,indx,nval)
 
    implicit none
!
    integer :: inc , itarg , n , nval
    integer , dimension(*) :: iarray , indx
    intent (in) iarray , inc , itarg , n
    intent (out) indx
    intent (inout) nval
!
    integer :: i , ina
!
    ina = 1
    nval = 0
    if ( inc < 0 ) ina = (-inc)*(n-1) + 1
    do i = 1 , n
      if ( iarray(ina) == itarg ) then
        nval = nval + 1
        indx(nval) = i
      end if
      ina = ina + inc
    end do
 
  end subroutine wheneq
!
  subroutine whenne(n,array,inc,itarg,indx,nval)
 
    implicit none
!
    integer :: inc , itarg , n , nval
    integer , dimension(*) :: array , indx
    intent (in) array , inc , itarg , n
    intent (out) indx
    intent (inout) nval
!
    integer :: i , ina
!
    ina = 1
    nval = 0
    if ( inc < 0 ) ina = (-inc)*(n-1) + 1
    do i = 1 , n
      if ( array(ina) /= itarg ) then
        nval = nval + 1
        indx(nval) = i
      end if
      ina = ina + inc
    end do
  end subroutine whenne
!
  subroutine whenfgt(n,array,inc,rtarg,indx,nval)
! 
    implicit none
!
    integer :: inc , n , nval
    real(8) :: rtarg
    real(8) , dimension(*) :: array
    integer , dimension(*) :: indx
    intent (in) array , inc , n , rtarg
    intent (out) indx
    intent (inout) nval
!
    integer :: i , ina
!
    ina = 1
    nval = 0
    if ( inc < 0 ) ina = (-inc)*(n-1) + 1
    do i = 1 , n
      if ( array(ina) > rtarg ) then
        nval = nval + 1
        indx(nval) = i
      end if
      ina = ina + inc
    end do
 
  end subroutine whenfgt
!
  subroutine whenflt(n,array,inc,rtarg,indx,nval)
! 
    implicit none
!
    integer :: inc , n , nval
    real(8) :: rtarg
    real(8) , dimension(*) :: array
    integer , dimension(*) :: indx
    intent (in) array , inc , n , rtarg
    intent (out) indx , nval
!
    integer :: i , ina
!
    ina = 1
    nval = 0
    if ( inc < 0 ) ina = (-inc)*(n-1) + 1
    do i = 1 , n
      if ( array(ina) < rtarg ) then
        nval = nval + 1
        indx(nval) = i
      end if
      ina = ina + inc
    end do
  end subroutine whenflt
!
  function intmax(n,imax,inc)
!
    implicit none
!
    integer :: inc , n
    integer :: intmax
    integer , dimension(*) :: imax
    intent (in) inc , imax , n
!
    integer :: i , mx
!
    mx = imax(1)
    intmax = 1
    do i = 1 + inc , inc*n , inc
      if ( imax(i) > mx ) then
        mx = imax(i)
        intmax = i
      end if
    end do
  end function intmax
!
  function xalpha(w,uu,g,e)
    implicit none
    real(8) :: xalpha
    real(8) , intent(in) :: w , uu , g , e
    xalpha = 0.75D0*w*uu*((d_one+g*(d_one-w))/(d_one-e*e*uu*uu))
  end function xalpha
  function xgamma(w,uu,g,e)
    implicit none
    real(8) :: xgamma
    real(8) , intent(in) :: w , uu , g , e
    xgamma = (w*d_half)*((3.0D0*g*(d_one-w)*uu*uu+d_one)/(d_one-e*e*uu*uu))
  end function xgamma
  function el(w,g)
    implicit none
    real(8) :: el
    real(8) , intent(in) :: w , g
    el = dsqrt(3.0D0*(d_one-w)*(d_one-w*g))
  end function el
  function taus(w,f,t)
    implicit none
    real(8) :: taus
    real(8) , intent(in) :: w , f , t
    taus = (d_one-w*f)*t
  end function taus
  function omgs(w,f)
    implicit none
    real(8) :: omgs
    real(8) , intent(in) :: w , f
    omgs = (d_one-f)*w/(d_one-w*f)
  end function omgs
  function asys(g,f)
    implicit none
    real(8) :: asys
    real(8) , intent(in) :: g , f
    asys = (g-f)/(d_one-f)
  end function asys
  function u(w,g,e)
    implicit none
    real(8) :: u
    real(8) , intent(in) :: w , g , e
    u = 1.50D0*(d_one-w*g)/e
  end function u
  function n(uu,et)
    implicit none
    real(8) :: n
    real(8) , intent(in) :: uu , et
    n = ((uu+d_one)*(uu+d_one)/et)-((uu-d_one)*(uu-d_one)*et)
  end function n
  function dbvt(t)
!   Derivative of planck function at 9.6 micro-meter wavelength
    implicit none
    real(8) :: dbvt
    real(8) , intent(in) :: t
    dbvt = (-2.8911366682D-4 + (2.3771251896D-6+1.1305188929D-10*t)*t) /  &
            (d_one+(-6.1364820707D-3+1.5550319767D-5*t)*t)
  end function dbvt
  function fo3(ux,vx)
!   an absorption function factor
    implicit none
    real(8) :: fo3
    real(8) , intent(in) :: ux , vx
    fo3 = ux/dsqrt(d_four+ux*(d_one+vx))
  end function fo3
!
end module mod_radiation
