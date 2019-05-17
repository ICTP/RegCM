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

module mod_sunorbit
  !
  ! Compute solar orbit parameters
  !
  use mod_realkinds
  use mod_intkinds
  use mod_message
  use mod_stdio
  use mod_constants , only : mathpi , degrad , d_zero , d_one , d_two , d_half
  use mod_dynparam , only : dayspy , vernal_equinox

  implicit none

  private

  real(rk8) , parameter :: ORB_ECCEN_MIN  =   0.0_rk8 ! min value for eccen
  real(rk8) , parameter :: ORB_ECCEN_MAX  =   0.1_rk8 ! max value for eccen
  real(rk8) , parameter :: ORB_OBLIQ_MIN  = -90.0_rk8 ! min value for obliq
  real(rk8) , parameter :: ORB_OBLIQ_MAX  = +90.0_rk8 ! max value for obliq
  real(rk8) , parameter :: ORB_MVELP_MIN  =   0.0_rk8 ! min value for mvelp
  real(rk8) , parameter :: ORB_MVELP_MAX  = 360.0_rk8 ! max value for mvelp
  real(rk8) , parameter :: ORB_UNDEF_REAL = 1.e36_rk8
  integer(ik4) , parameter :: ORB_UNDEF_INT  = 2000000000

  public :: orb_cosz , orb_decl , orb_params

  interface orb_cosz
    module procedure orb_cosz_r8
    module procedure orb_cosz_r4
  end interface orb_cosz

  contains
  !
  ! FUNCTION to return the cosine of the solar zenith angle.
  !
  real(rk8) function orb_cosz_r8(jday,lat,lon,declin)
    implicit none
    real(rk8) , intent(in) :: jday   ! Julian cal day
    real(rk8) , intent(in) :: lat    ! Centered latitude (radians)
    real(rk8) , intent(in) :: lon    ! Centered longitude (radians)
    real(rk8) , intent(in) :: declin ! Solar declination (radians)
    orb_cosz_r8 = sin(lat)*sin(declin) - &
                  cos(lat)*cos(declin)*cos(jday*2.0_rk8*mathpi + lon)
  end function orb_cosz_r8
  !
  real(rk4) function orb_cosz_r4(jday,lat,lon,declin)
    implicit none
    real(rk4) , intent(in) :: jday   ! Julian cal day
    real(rk4) , intent(in) :: lat    ! Centered latitude (radians)
    real(rk4) , intent(in) :: lon    ! Centered longitude (radians)
    real(rk8) , intent(in) :: declin ! Solar declination (radians)
    real(rk8) :: dlat , dlon , djday
    dlat = lat
    dlon = lon
    djday = jday
    orb_cosz_r4 = real(sin(dlat)*sin(declin) - &
                   cos(dlat)*cos(declin)*cos(djday*2.0_rk8*mathpi + dlon),rk4)
  end function orb_cosz_r4
  !
  ! Calculate earths orbital parameters using Dave Threshers formula which
  ! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term
  ! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy
  ! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
  !
  !
  subroutine orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp)
    implicit none
    integer(ik4) , intent(in) :: iyear_AD  ! Year to calculate orbit for
    real(rk8) , intent(inout) :: eccen   ! orbital eccentricity
    real(rk8) , intent(inout) :: obliq   ! obliquity in degrees
    real(rk8) , intent(inout) :: mvelp   ! moving vernal equinox long
    real(rk8) , intent(out) :: obliqr    ! Earths obliquity in rad
    real(rk8) , intent(out) :: lambm0    ! Mean long of perihelion at
                                         ! vernal equinox (radians)
    real(rk8) , intent(out)   :: mvelpp  ! moving vernal equinox long
                                         ! of perihelion plus pi (rad)
    ! # of elements in series wrt obliquity
    integer(ik4) , parameter :: poblen = 47
    ! # of elements in series wrt eccentricity
    integer(ik4) , parameter :: pecclen = 19
    ! # of elements in series wrt vernal equinox
    integer(ik4) , parameter :: pmvelen = 78
    ! arc sec to deg conversion
    real(rk8) , parameter :: psecdeg = 1.0_rk8/3600.0_rk8

    real(rk8) :: yb4_1950AD         ! number of years before 1950 AD

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    ! amplitudes for obliquity cos series
    real(rk8) , parameter :: obamp(poblen) =  &
          [   -2462.2214466_rk8, -857.3232075_rk8, -629.3231835_rk8,   &
                -414.2804924_rk8, -311.7632587_rk8,  308.9408604_rk8,   &
                -162.5533601_rk8, -116.1077911_rk8,  101.1189923_rk8,   &
                 -67.6856209_rk8,   24.9079067_rk8,   22.5811241_rk8,   &
                 -21.1648355_rk8,  -15.6549876_rk8,   15.3936813_rk8,   &
                  14.6660938_rk8,  -11.7273029_rk8,   10.2742696_rk8,   &
                   6.4914588_rk8,    5.8539148_rk8,   -5.4872205_rk8,   &
                  -5.4290191_rk8,    5.1609570_rk8,    5.0786314_rk8,   &
                  -4.0735782_rk8,    3.7227167_rk8,    3.3971932_rk8,   &
                  -2.8347004_rk8,   -2.6550721_rk8,   -2.5717867_rk8,   &
                  -2.4712188_rk8,    2.4625410_rk8,    2.2464112_rk8,   &
                  -2.0755511_rk8,   -1.9713669_rk8,   -1.8813061_rk8,   &
                  -1.8468785_rk8,    1.8186742_rk8,    1.7601888_rk8,   &
                  -1.5428851_rk8,    1.4738838_rk8,   -1.4593669_rk8,   &
                   1.4192259_rk8,   -1.1818980_rk8,    1.1756474_rk8,   &
                  -1.1316126_rk8,    1.0896928_rk8]

    ! rates for obliquity cosine series
    real(rk8) , parameter :: obrate(poblen) = &
            [  31.609974_rk8, 32.620504_rk8, 24.172203_rk8,   &
                31.983787_rk8, 44.828336_rk8, 30.973257_rk8,   &
                43.668246_rk8, 32.246691_rk8, 30.599444_rk8,   &
                42.681324_rk8, 43.836462_rk8, 47.439436_rk8,   &
                63.219948_rk8, 64.230478_rk8,  1.010530_rk8,   &
                 7.437771_rk8, 55.782177_rk8,  0.373813_rk8,   &
                13.218362_rk8, 62.583231_rk8, 63.593761_rk8,   &
                76.438310_rk8, 45.815258_rk8,  8.448301_rk8,   &
                56.792707_rk8, 49.747842_rk8, 12.058272_rk8,   &
                75.278220_rk8, 65.241008_rk8, 64.604291_rk8,   &
                 1.647247_rk8,  7.811584_rk8, 12.207832_rk8,   &
                63.856665_rk8, 56.155990_rk8, 77.448840_rk8,   &
                 6.801054_rk8, 62.209418_rk8, 20.656133_rk8,   &
                48.344406_rk8, 55.145460_rk8, 69.000539_rk8,   &
                11.071350_rk8, 74.291298_rk8, 11.047742_rk8,   &
                 0.636717_rk8, 12.844549_rk8]

    ! phases for obliquity cosine series
    real(rk8) , parameter :: obphas(poblen) = &
          [    251.9025_rk8, 280.8325_rk8, 128.3057_rk8,   &
                292.7252_rk8,  15.3747_rk8, 263.7951_rk8,   &
                308.4258_rk8, 240.0099_rk8, 222.9725_rk8,   &
                268.7809_rk8, 316.7998_rk8, 319.6024_rk8,   &
                143.8050_rk8, 172.7351_rk8,  28.9300_rk8,   &
                123.5968_rk8,  20.2082_rk8,  40.8226_rk8,   &
                123.4722_rk8, 155.6977_rk8, 184.6277_rk8,   &
                267.2772_rk8,  55.0196_rk8, 152.5268_rk8,   &
                 49.1382_rk8, 204.6609_rk8,  56.5233_rk8,   &
                200.3284_rk8, 201.6651_rk8, 213.5577_rk8,   &
                 17.0374_rk8, 164.4194_rk8,  94.5422_rk8,   &
                131.9124_rk8,  61.0309_rk8, 296.2073_rk8,   &
                135.4894_rk8, 114.8750_rk8, 247.0691_rk8,   &
                256.6114_rk8,  32.1008_rk8, 143.6804_rk8,   &
                 16.8784_rk8, 160.6835_rk8,  27.5932_rk8,   &
                348.1074_rk8,  82.6496_rk8]

    ! Cosine/sine series data for computation of eccentricity and fixed vernal
    ! equinox longitude of perihelion (fvelp): amplitude,
    ! rate (arc seconds/year), phase (degrees).

    ! ampl for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecamp (pecclen) = &
          [   0.01860798_rk8,  0.01627522_rk8, -0.01300660_rk8,   &
               0.00988829_rk8, -0.00336700_rk8,  0.00333077_rk8,   &
              -0.00235400_rk8,  0.00140015_rk8,  0.00100700_rk8,   &
               0.00085700_rk8,  0.00064990_rk8,  0.00059900_rk8,   &
               0.00037800_rk8, -0.00033700_rk8,  0.00027600_rk8,   &
               0.00018200_rk8, -0.00017400_rk8, -0.00012400_rk8,   &
               0.00001250_rk8]

    ! rates for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecrate(pecclen) = &
          [    4.2072050_rk8,  7.3460910_rk8, 17.8572630_rk8,  &
               17.2205460_rk8, 16.8467330_rk8,  5.1990790_rk8,  &
               18.2310760_rk8, 26.2167580_rk8,  6.3591690_rk8,  &
               16.2100160_rk8,  3.0651810_rk8, 16.5838290_rk8,  &
               18.4939800_rk8,  6.1909530_rk8, 18.8677930_rk8,  &
               17.4255670_rk8,  6.1860010_rk8, 18.4174410_rk8,  &
                0.6678630_rk8]

    ! phases for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecphas(pecclen) = &
          [    28.620089_rk8, 193.788772_rk8, 308.307024_rk8,  &
               320.199637_rk8, 279.376984_rk8,  87.195000_rk8,  &
               349.129677_rk8, 128.443387_rk8, 154.143880_rk8,  &
               291.269597_rk8, 114.860583_rk8, 332.092251_rk8,  &
               296.414411_rk8, 145.769910_rk8, 337.237063_rk8,  &
               152.092288_rk8, 126.839891_rk8, 210.667199_rk8,  &
                72.108838_rk8]

    ! Sine series data for computation of moving vernal equinox longitude of
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    ! amplitudes for mvelp sine series
    real(rk8) , parameter :: mvamp (pmvelen) = &
          [   7391.0225890_rk8, 2555.1526947_rk8, 2022.7629188_rk8,  &
              -1973.6517951_rk8, 1240.2321818_rk8,  953.8679112_rk8,  &
               -931.7537108_rk8,  872.3795383_rk8,  606.3544732_rk8,  &
               -496.0274038_rk8,  456.9608039_rk8,  346.9462320_rk8,  &
               -305.8412902_rk8,  249.6173246_rk8, -199.1027200_rk8,  &
                191.0560889_rk8, -175.2936572_rk8,  165.9068833_rk8,  &
                161.1285917_rk8,  139.7878093_rk8, -133.5228399_rk8,  &
                117.0673811_rk8,  104.6907281_rk8,   95.3227476_rk8,  &
                 86.7824524_rk8,   86.0857729_rk8,   70.5893698_rk8,  &
                -69.9719343_rk8,  -62.5817473_rk8,   61.5450059_rk8,  &
                -57.9364011_rk8,   57.1899832_rk8,  -57.0236109_rk8,  &
                -54.2119253_rk8,   53.2834147_rk8,   52.1223575_rk8,  &
                -49.0059908_rk8,  -48.3118757_rk8,  -45.4191685_rk8,  &
                -42.2357920_rk8,  -34.7971099_rk8,   34.4623613_rk8,  &
                -33.8356643_rk8,   33.6689362_rk8,  -31.2521586_rk8,  &
                -30.8798701_rk8,   28.4640769_rk8,  -27.1960802_rk8,  &
                 27.0860736_rk8,  -26.3437456_rk8,   24.7253740_rk8,  &
                 24.6732126_rk8,   24.4272733_rk8,   24.0127327_rk8,  &
                 21.7150294_rk8,  -21.5375347_rk8,   18.1148363_rk8,  &
                -16.9603104_rk8,  -16.1765215_rk8,   15.5567653_rk8,  &
                 15.4846529_rk8,   15.2150632_rk8,   14.5047426_rk8,  &
                -14.3873316_rk8,   13.1351419_rk8,   12.8776311_rk8,  &
                 11.9867234_rk8,   11.9385578_rk8,   11.7030822_rk8,  &
                 11.6018181_rk8,  -11.2617293_rk8,  -10.4664199_rk8,  &
                 10.4333970_rk8,  -10.2377466_rk8,   10.1934446_rk8,  &
                -10.1280191_rk8,   10.0289441_rk8,  -10.0034259_rk8]

    ! rates for mvelp sine series
    real(rk8) , parameter :: mvrate(pmvelen) = &
          [    31.609974_rk8, 32.620504_rk8, 24.172203_rk8,   &
                 0.636717_rk8, 31.983787_rk8,  3.138886_rk8,   &
                30.973257_rk8, 44.828336_rk8,  0.991874_rk8,   &
                 0.373813_rk8, 43.668246_rk8, 32.246691_rk8,   &
                30.599444_rk8,  2.147012_rk8, 10.511172_rk8,   &
                42.681324_rk8, 13.650058_rk8,  0.986922_rk8,   &
                 9.874455_rk8, 13.013341_rk8,  0.262904_rk8,   &
                 0.004952_rk8,  1.142024_rk8, 63.219948_rk8,   &
                 0.205021_rk8,  2.151964_rk8, 64.230478_rk8,   &
                43.836462_rk8, 47.439436_rk8,  1.384343_rk8,   &
                 7.437771_rk8, 18.829299_rk8,  9.500642_rk8,   &
                 0.431696_rk8,  1.160090_rk8, 55.782177_rk8,   &
                12.639528_rk8,  1.155138_rk8,  0.168216_rk8,   &
                 1.647247_rk8, 10.884985_rk8,  5.610937_rk8,   &
                12.658184_rk8,  1.010530_rk8,  1.983748_rk8,   &
                14.023871_rk8,  0.560178_rk8,  1.273434_rk8,   &
                12.021467_rk8, 62.583231_rk8, 63.593761_rk8,   &
                76.438310_rk8,  4.280910_rk8, 13.218362_rk8,   &
                17.818769_rk8,  8.359495_rk8, 56.792707_rk8,   &
                 8.448301_rk8,  1.978796_rk8,  8.863925_rk8,   &
                 0.186365_rk8,  8.996212_rk8,  6.771027_rk8,   &
                45.815258_rk8, 12.002811_rk8, 75.278220_rk8,   &
                65.241008_rk8, 18.870667_rk8, 22.009553_rk8,   &
                64.604291_rk8, 11.498094_rk8,  0.578834_rk8,   &
                 9.237738_rk8, 49.747842_rk8,  2.147012_rk8,   &
                 1.196895_rk8,  2.133898_rk8,  0.173168_rk8]

    ! phases for mvelp sine series
    real(rk8) , parameter :: mvphas(pmvelen) = &
          [    251.9025_rk8, 280.8325_rk8, 128.3057_rk8,   &
                348.1074_rk8, 292.7252_rk8, 165.1686_rk8,   &
                263.7951_rk8,  15.3747_rk8,  58.5749_rk8,   &
                 40.8226_rk8, 308.4258_rk8, 240.0099_rk8,   &
                222.9725_rk8, 106.5937_rk8, 114.5182_rk8,   &
                268.7809_rk8, 279.6869_rk8,  39.6448_rk8,   &
                126.4108_rk8, 291.5795_rk8, 307.2848_rk8,   &
                 18.9300_rk8, 273.7596_rk8, 143.8050_rk8,   &
                191.8927_rk8, 125.5237_rk8, 172.7351_rk8,   &
                316.7998_rk8, 319.6024_rk8,  69.7526_rk8,   &
                123.5968_rk8, 217.6432_rk8,  85.5882_rk8,   &
                156.2147_rk8,  66.9489_rk8,  20.2082_rk8,   &
                250.7568_rk8,  48.0188_rk8,   8.3739_rk8,   &
                 17.0374_rk8, 155.3409_rk8,  94.1709_rk8,   &
                221.1120_rk8,  28.9300_rk8, 117.1498_rk8,   &
                320.5095_rk8, 262.3602_rk8, 336.2148_rk8,   &
                233.0046_rk8, 155.6977_rk8, 184.6277_rk8,   &
                267.2772_rk8,  78.9281_rk8, 123.4722_rk8,   &
                188.7132_rk8, 180.1364_rk8,  49.1382_rk8,   &
                152.5268_rk8,  98.2198_rk8,  97.4808_rk8,   &
                221.5376_rk8, 168.2438_rk8, 161.1199_rk8,   &
                 55.0196_rk8, 262.6495_rk8, 200.3284_rk8,   &
                201.6651_rk8, 294.6547_rk8,  99.8233_rk8,   &
                213.5577_rk8, 154.1631_rk8, 232.7153_rk8,   &
                138.3034_rk8, 204.6609_rk8, 106.5938_rk8,   &
                250.4676_rk8, 332.3345_rk8,  27.3039_rk8]

    integer(ik4) :: i       ! Index for series summations
    real(rk8) :: obsum   ! Obliquity series summation
    real(rk8) :: cossum  ! Cos series summation for eccentricity/fvelp
    real(rk8) :: sinsum  ! Sin series summation for eccentricity/fvelp
    real(rk8) :: fvelp   ! Fixed vernal equinox long of perihelion
    real(rk8) :: mvsum   ! mvelp series summation
    real(rk8) :: beta    ! Intermediate argument for lambm0
    real(rk8) :: years   ! Years to time of interest ( pos <=> future)
    real(rk8) :: eccen2  ! eccentricity squared
    real(rk8) :: eccen3  ! eccentricity cubed

    ! Check for flag to use input orbit parameters

    if ( iyear_AD == ORB_UNDEF_INT ) then

      ! Check input obliq, eccen, and mvelp to ensure reasonable

      if ( obliq == ORB_UNDEF_REAL ) then
        write(stdout,*) ' Have to specify orbital parameters:'
        write(stdout,*) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
        write(stdout,*) 'iyear_AD is the year to simulate orbit (ie. 1950): '
        write(stdout,*) 'obliq, eccen, mvelp specify the orbit directly:'
        write(stdout,*) 'The AMIP II settings (for a 1995 orbit) are: '
        write(stdout,*) ' obliq =  23.4441'
        write(stdout,*) ' eccen =   0.016715'
        write(stdout,*) ' mvelp = 102.7'
        call die(__FILE__,'unreasonable obliq')
      end if
      if ( (obliq < ORB_OBLIQ_MIN).or.(obliq > ORB_OBLIQ_MAX) ) then
        write(stdout,*) 'Input obliquity unreasonable: ', obliq
        call die(__FILE__,'unreasonable obliq')
      end if
      if ( (eccen < ORB_ECCEN_MIN).or.(eccen > ORB_ECCEN_MAX) ) then
        write(stdout,*) 'Input eccentricity unreasonable: ', eccen
        call die(__FILE__,'unreasonable eccen')
      end if
      if ( (mvelp < ORB_MVELP_MIN).or.(mvelp > ORB_MVELP_MAX) ) then
        write(stdout,*) 'Input mvelp unreasonable: ' , mvelp
        call die(__FILE__,'unreasonable mvelp')
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

    else  ! Otherwise calculate based on years before present

      yb4_1950AD = 1950.0_rk8 - real(iyear_AD,rk8)
      if ( abs(yb4_1950AD) > 1000000.0_rk8 )then
        write(stdout,*) 'orbit only valid for years+-1000000'
        write(stdout,*) 'Relative to 1950 AD'
        write(stdout,*) '# of years before 1950: ',yb4_1950AD
        write(stdout,*) 'Year to simulate was  : ',iyear_AD
        call die(__FILE__,'unreasonable year')
      end if

      ! The following calculates the earths obliquity, orbital eccentricity
      ! (and various powers of it) and vernal equinox mean longitude of
      ! perihelion for years in the past (future = negative of years past),
      ! using constants (see parameter section) given in the program of:
      !
      ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
      ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
      ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
      !
      ! and formulas given in the paper (where less precise constants are also
      ! given):
      !
      ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
      ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
      !
      ! The algorithm is valid only to 1,000,000 years past or hence.
      ! For a solution valid to 5-10 million years past see the above author.
      ! Algorithm below is better for years closer to present than is the
      ! 5-10 million year solution.
      !
      ! Years to time of interest must be negative of years before present
      ! (1950) in formulas that follow.

      years = - yb4_1950AD

      ! In the summations below, cosine or sine arguments, which end up in
      ! degrees, must be converted to radians via multiplication by degrad.
      !
      ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
      ! degrees. Convert the amplitudes and rates, which are in arc secs, into
      ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
      ! factor).  For obliq, first term is Berger 1978 epsilon star; second
      ! term is series summation in degrees.

      obsum = d_zero
      do i = 1 , poblen
        obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
                obphas(i))*degrad)
      end do
      obliq = 23.320556_rk8 + obsum

      ! Summation of cosine and sine series for computation of eccentricity
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation
      ! of moving vernal equinox longitude of perihelion.  Convert the rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.

      cossum = d_zero
      do i = 1 , pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      sinsum = d_zero
      do i = 1 , pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      ! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen

      ! A series of cases for fvelp, which is in radians.

      if (abs(cossum) <= 1.0e-8_rk8) then
        if (sinsum == d_zero) then
          fvelp = d_zero
        else if (sinsum < d_zero) then
          fvelp = 1.5_rk8*mathpi
        else if (sinsum > d_zero) then
          fvelp = d_half*mathpi
        end if
      else if (cossum < d_zero) then
        fvelp = atan(sinsum/cossum) + mathpi
      else if (cossum > d_zero) then
        if (sinsum < d_zero) then
          fvelp = atan(sinsum/cossum) + 2.0_rk8*mathpi
        else
          fvelp = atan(sinsum/cossum)
        end if
      end if

      ! Summation of sin series for computation of moving vernal equinox long
      ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
      ! first term is fvelp in degrees; second term is Berger 1978 psi bar
      ! times years and in degrees; third term is Berger 1978 zeta; fourth
      ! term is series summation in degrees.  Convert the amplitudes and rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
      ! Series summation plus second and third terms constitute Berger 1978
      ! psi, which is the general precession.

      mvsum = d_zero
      do i = 1 , pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
                mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273_rk8*psecdeg*years + 3.392506_rk8 + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp < d_zero)
        mvelp = mvelp + 360.0_rk8
      end do
      do while (mvelp >= 360.0_rk8)
        mvelp = mvelp - 360.0_rk8
      end do

    end if  ! end of test on whether to calculate or use input orbital params

    ! Orbit needs the obliquity in radians

    obliqr = obliq*degrad

    ! 180 degrees must be added to mvelp since observations are made from the
    ! earth and the sun is considered (wrongly for the algorithm) to go around
    ! the earth. For a more graphic explanation see Appendix B in:
    !
    ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
    ! Periods.  J. of Geophysical Research 98:10,341-10,362.
    !
    ! Additionally, orbit will need this value in radians. So mvelp becomes
    ! mvelpp (mvelp plus pi)

    mvelpp = (mvelp + 180._rk8)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = sqrt(d_one - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = d_two*((d_half*eccen+0.125_rk8*eccen3) * &
            (d_one+beta)*sin(mvelpp) - &
              0.250_rk8*eccen2*(d_half+beta)*sin(d_two*mvelpp) + &
              0.125_rk8*eccen3*(d_one/3.0_rk8 + beta)*sin(3.0_rk8*mvelpp))

  end subroutine orb_params
  !
  ! Compute earth/orbit parameters using formula suggested by
  ! Duane Thresher.
  !
  subroutine orb_decl(calday,eccen,mvelpp,lambm0,obliqr,delta,eccf)
    implicit none
    real(rk8) , intent(in) :: calday ! Calendar day, including fraction
    real(rk8) , intent(in) :: eccen  ! Eccentricity
    real(rk8) , intent(in) :: obliqr ! Earths obliquity in radians
    real(rk8) , intent(in) :: lambm0 ! Mean long of perihelion at the
                                     ! vernal equinox (radians)
    real(rk8) , intent(in) :: mvelpp ! moving vernal equinox longitude
                                     ! of perihelion plus pi (radians)
    real(rk8) , intent(out) :: delta ! Solar declination angle in rad
    real(rk8) , intent(out) :: eccf  ! Earth-sun distance factor (ie. (1/r)**2)

    real(rk8) :: lambm  ! Lambda m, mean long of perihelion (rad)
    real(rk8) :: lmm    ! Intermediate argument involving lambm
    real(rk8) :: lamb   ! Lambda, the earths long of perihelion
    real(rk8) :: invrho ! Inverse normalized sun/earth distance
    real(rk8) :: sinl   ! Sine of lmm

    ! Compute eccentricity factor and solar declination using
    ! day value where a round day (such as 213.0) refers to 0z at
    ! Greenwich longitude.
    !
    ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    ! 35:2362-2367.
    !
    ! To get the earths true longitude (position in orbit; lambda in Berger
    ! 1978) which is necessary to find the eccentricity factor and declination,
    ! must first calculate the mean longitude (lambda m in Berger 1978) at
    ! the present day.  This is done by adding to lambm0 (the mean longitude
    ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    ! an increment (delta lambda m in Berger 1978) that is the number of
    ! days past or before (a negative increment) the vernal equinox divided by
    ! the days in a model year times the 2*mathpi radians in a complete orbit.

    lambm = lambm0 + (calday - vernal_equinox)*d_two*mathpi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(d_two*sinl + eccen*(1.25_rk8*sin(d_two*lmm)  &
          + eccen*((13.0_rk8/12.0_rk8)*sin(3._rk8*lmm) - 0.25_rk8*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually
    ! inverse rho will be used), and thus the eccentricity factor (eccf),
    ! can be calculated from formulas given in Berger 1978.

    invrho = (d_one + eccen*cos(lamb - mvelpp)) / (d_one - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta  = asin(sin(obliqr)*sin(lamb))
    eccf   = invrho*invrho
  end subroutine orb_decl

end module mod_sunorbit
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
