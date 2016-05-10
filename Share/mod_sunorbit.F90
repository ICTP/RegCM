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

  implicit none

  private

  real(rkx) , parameter :: ORB_ECCEN_MIN  =   0.0_rkx ! min value for eccen
  real(rkx) , parameter :: ORB_ECCEN_MAX  =   0.1_rkx ! max value for eccen
  real(rkx) , parameter :: ORB_OBLIQ_MIN  = -90.0_rkx ! min value for obliq
  real(rkx) , parameter :: ORB_OBLIQ_MAX  = +90.0_rkx ! max value for obliq
  real(rkx) , parameter :: ORB_MVELP_MIN  =   0.0_rkx ! min value for mvelp
  real(rkx) , parameter :: ORB_MVELP_MAX  = 360.0_rkx ! max value for mvelp
  real(rkx) , parameter :: ORB_UNDEF_REAL = 1.e36_rkx
  integer(ik4) , parameter :: ORB_UNDEF_INT  = 2000000000

  public :: orb_cosz , orb_decl , orb_params

  contains
  !
  ! FUNCTION to return the cosine of the solar zenith angle.
  ! Assumes 365.0 days/year.
  !
  real(rkx) function orb_cosz(jday,lat,lon,declin)
    implicit none
    real(rkx) , intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real(rkx) , intent(in) :: lat    ! Centered latitude (radians)
    real(rkx) , intent(in) :: lon    ! Centered longitude (radians)
    real(rkx) , intent(in) :: declin ! Solar declination (radians)
    orb_cosz = sin(lat)*sin(declin) - &
               cos(lat)*cos(declin)*cos(jday*2.0_rkx*mathpi + lon)
  end function orb_cosz
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
    real(rkx) , intent(inout) :: eccen   ! orbital eccentricity
    real(rkx) , intent(inout) :: obliq   ! obliquity in degrees
    real(rkx) , intent(inout) :: mvelp   ! moving vernal equinox long
    real(rkx) , intent(out) :: obliqr    ! Earths obliquity in rad
    real(rkx) , intent(out) :: lambm0    ! Mean long of perihelion at
                                         ! vernal equinox (radians)
    real(rkx) , intent(out)   :: mvelpp  ! moving vernal equinox long
                                         ! of perihelion plus pi (rad)
    ! # of elements in series wrt obliquity
    integer(ik4) , parameter :: poblen = 47
    ! # of elements in series wrt eccentricity
    integer(ik4) , parameter :: pecclen = 19
    ! # of elements in series wrt vernal equinox
    integer(ik4) , parameter :: pmvelen = 78
    ! arc sec to deg conversion
    real(rkx) , parameter :: psecdeg = 1.0_rkx/3600.0_rkx

    real(rkx) :: yb4_1950AD         ! number of years before 1950 AD

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    ! amplitudes for obliquity cos series
    real(rkx) , parameter :: obamp(poblen) =  &
          (/   -2462.2214466_rkx, -857.3232075_rkx, -629.3231835_rkx,   &
                -414.2804924_rkx, -311.7632587_rkx,  308.9408604_rkx,   &
                -162.5533601_rkx, -116.1077911_rkx,  101.1189923_rkx,   &
                 -67.6856209_rkx,   24.9079067_rkx,   22.5811241_rkx,   &
                 -21.1648355_rkx,  -15.6549876_rkx,   15.3936813_rkx,   &
                  14.6660938_rkx,  -11.7273029_rkx,   10.2742696_rkx,   &
                   6.4914588_rkx,    5.8539148_rkx,   -5.4872205_rkx,   &
                  -5.4290191_rkx,    5.1609570_rkx,    5.0786314_rkx,   &
                  -4.0735782_rkx,    3.7227167_rkx,    3.3971932_rkx,   &
                  -2.8347004_rkx,   -2.6550721_rkx,   -2.5717867_rkx,   &
                  -2.4712188_rkx,    2.4625410_rkx,    2.2464112_rkx,   &
                  -2.0755511_rkx,   -1.9713669_rkx,   -1.8813061_rkx,   &
                  -1.8468785_rkx,    1.8186742_rkx,    1.7601888_rkx,   &
                  -1.5428851_rkx,    1.4738838_rkx,   -1.4593669_rkx,   &
                   1.4192259_rkx,   -1.1818980_rkx,    1.1756474_rkx,   &
                  -1.1316126_rkx,    1.0896928_rkx/)

    ! rates for obliquity cosine series
    real(rkx) , parameter :: obrate(poblen) = &
            (/  31.609974_rkx, 32.620504_rkx, 24.172203_rkx,   &
                31.983787_rkx, 44.828336_rkx, 30.973257_rkx,   &
                43.668246_rkx, 32.246691_rkx, 30.599444_rkx,   &
                42.681324_rkx, 43.836462_rkx, 47.439436_rkx,   &
                63.219948_rkx, 64.230478_rkx,  1.010530_rkx,   &
                 7.437771_rkx, 55.782177_rkx,  0.373813_rkx,   &
                13.218362_rkx, 62.583231_rkx, 63.593761_rkx,   &
                76.438310_rkx, 45.815258_rkx,  8.448301_rkx,   &
                56.792707_rkx, 49.747842_rkx, 12.058272_rkx,   &
                75.278220_rkx, 65.241008_rkx, 64.604291_rkx,   &
                 1.647247_rkx,  7.811584_rkx, 12.207832_rkx,   &
                63.856665_rkx, 56.155990_rkx, 77.448840_rkx,   &
                 6.801054_rkx, 62.209418_rkx, 20.656133_rkx,   &
                48.344406_rkx, 55.145460_rkx, 69.000539_rkx,   &
                11.071350_rkx, 74.291298_rkx, 11.047742_rkx,   &
                 0.636717_rkx, 12.844549_rkx/)

    ! phases for obliquity cosine series
    real(rkx) , parameter :: obphas(poblen) = &
          (/    251.9025_rkx, 280.8325_rkx, 128.3057_rkx,   &
                292.7252_rkx,  15.3747_rkx, 263.7951_rkx,   &
                308.4258_rkx, 240.0099_rkx, 222.9725_rkx,   &
                268.7809_rkx, 316.7998_rkx, 319.6024_rkx,   &
                143.8050_rkx, 172.7351_rkx,  28.9300_rkx,   &
                123.5968_rkx,  20.2082_rkx,  40.8226_rkx,   &
                123.4722_rkx, 155.6977_rkx, 184.6277_rkx,   &
                267.2772_rkx,  55.0196_rkx, 152.5268_rkx,   &
                 49.1382_rkx, 204.6609_rkx,  56.5233_rkx,   &
                200.3284_rkx, 201.6651_rkx, 213.5577_rkx,   &
                 17.0374_rkx, 164.4194_rkx,  94.5422_rkx,   &
                131.9124_rkx,  61.0309_rkx, 296.2073_rkx,   &
                135.4894_rkx, 114.8750_rkx, 247.0691_rkx,   &
                256.6114_rkx,  32.1008_rkx, 143.6804_rkx,   &
                 16.8784_rkx, 160.6835_rkx,  27.5932_rkx,   &
                348.1074_rkx,  82.6496_rkx/)

    ! Cosine/sine series data for computation of eccentricity and fixed vernal
    ! equinox longitude of perihelion (fvelp): amplitude,
    ! rate (arc seconds/year), phase (degrees).

    ! ampl for eccen/fvelp cos/sin series
    real(rkx) , parameter :: ecamp (pecclen) = &
          (/   0.01860798_rkx,  0.01627522_rkx, -0.01300660_rkx,   &
               0.00988829_rkx, -0.00336700_rkx,  0.00333077_rkx,   &
              -0.00235400_rkx,  0.00140015_rkx,  0.00100700_rkx,   &
               0.00085700_rkx,  0.00064990_rkx,  0.00059900_rkx,   &
               0.00037800_rkx, -0.00033700_rkx,  0.00027600_rkx,   &
               0.00018200_rkx, -0.00017400_rkx, -0.00012400_rkx,   &
               0.00001250_rkx/)

    ! rates for eccen/fvelp cos/sin series
    real(rkx) , parameter :: ecrate(pecclen) = &
          (/    4.2072050_rkx,  7.3460910_rkx, 17.8572630_rkx,  &
               17.2205460_rkx, 16.8467330_rkx,  5.1990790_rkx,  &
               18.2310760_rkx, 26.2167580_rkx,  6.3591690_rkx,  &
               16.2100160_rkx,  3.0651810_rkx, 16.5838290_rkx,  &
               18.4939800_rkx,  6.1909530_rkx, 18.8677930_rkx,  &
               17.4255670_rkx,  6.1860010_rkx, 18.4174410_rkx,  &
                0.6678630_rkx/)

    ! phases for eccen/fvelp cos/sin series
    real(rkx) , parameter :: ecphas(pecclen) = &
          (/    28.620089_rkx, 193.788772_rkx, 308.307024_rkx,  &
               320.199637_rkx, 279.376984_rkx,  87.195000_rkx,  &
               349.129677_rkx, 128.443387_rkx, 154.143880_rkx,  &
               291.269597_rkx, 114.860583_rkx, 332.092251_rkx,  &
               296.414411_rkx, 145.769910_rkx, 337.237063_rkx,  &
               152.092288_rkx, 126.839891_rkx, 210.667199_rkx,  &
                72.108838_rkx/)

    ! Sine series data for computation of moving vernal equinox longitude of
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    ! amplitudes for mvelp sine series
    real(rkx) , parameter :: mvamp (pmvelen) = &
          (/   7391.0225890_rkx, 2555.1526947_rkx, 2022.7629188_rkx,  &
              -1973.6517951_rkx, 1240.2321818_rkx,  953.8679112_rkx,  &
               -931.7537108_rkx,  872.3795383_rkx,  606.3544732_rkx,  &
               -496.0274038_rkx,  456.9608039_rkx,  346.9462320_rkx,  &
               -305.8412902_rkx,  249.6173246_rkx, -199.1027200_rkx,  &
                191.0560889_rkx, -175.2936572_rkx,  165.9068833_rkx,  &
                161.1285917_rkx,  139.7878093_rkx, -133.5228399_rkx,  &
                117.0673811_rkx,  104.6907281_rkx,   95.3227476_rkx,  &
                 86.7824524_rkx,   86.0857729_rkx,   70.5893698_rkx,  &
                -69.9719343_rkx,  -62.5817473_rkx,   61.5450059_rkx,  &
                -57.9364011_rkx,   57.1899832_rkx,  -57.0236109_rkx,  &
                -54.2119253_rkx,   53.2834147_rkx,   52.1223575_rkx,  &
                -49.0059908_rkx,  -48.3118757_rkx,  -45.4191685_rkx,  &
                -42.2357920_rkx,  -34.7971099_rkx,   34.4623613_rkx,  &
                -33.8356643_rkx,   33.6689362_rkx,  -31.2521586_rkx,  &
                -30.8798701_rkx,   28.4640769_rkx,  -27.1960802_rkx,  &
                 27.0860736_rkx,  -26.3437456_rkx,   24.7253740_rkx,  &
                 24.6732126_rkx,   24.4272733_rkx,   24.0127327_rkx,  &
                 21.7150294_rkx,  -21.5375347_rkx,   18.1148363_rkx,  &
                -16.9603104_rkx,  -16.1765215_rkx,   15.5567653_rkx,  &
                 15.4846529_rkx,   15.2150632_rkx,   14.5047426_rkx,  &
                -14.3873316_rkx,   13.1351419_rkx,   12.8776311_rkx,  &
                 11.9867234_rkx,   11.9385578_rkx,   11.7030822_rkx,  &
                 11.6018181_rkx,  -11.2617293_rkx,  -10.4664199_rkx,  &
                 10.4333970_rkx,  -10.2377466_rkx,   10.1934446_rkx,  &
                -10.1280191_rkx,   10.0289441_rkx,  -10.0034259_rkx/)

    ! rates for mvelp sine series
    real(rkx) , parameter :: mvrate(pmvelen) = &
          (/    31.609974_rkx, 32.620504_rkx, 24.172203_rkx,   &
                 0.636717_rkx, 31.983787_rkx,  3.138886_rkx,   &
                30.973257_rkx, 44.828336_rkx,  0.991874_rkx,   &
                 0.373813_rkx, 43.668246_rkx, 32.246691_rkx,   &
                30.599444_rkx,  2.147012_rkx, 10.511172_rkx,   &
                42.681324_rkx, 13.650058_rkx,  0.986922_rkx,   &
                 9.874455_rkx, 13.013341_rkx,  0.262904_rkx,   &
                 0.004952_rkx,  1.142024_rkx, 63.219948_rkx,   &
                 0.205021_rkx,  2.151964_rkx, 64.230478_rkx,   &
                43.836462_rkx, 47.439436_rkx,  1.384343_rkx,   &
                 7.437771_rkx, 18.829299_rkx,  9.500642_rkx,   &
                 0.431696_rkx,  1.160090_rkx, 55.782177_rkx,   &
                12.639528_rkx,  1.155138_rkx,  0.168216_rkx,   &
                 1.647247_rkx, 10.884985_rkx,  5.610937_rkx,   &
                12.658184_rkx,  1.010530_rkx,  1.983748_rkx,   &
                14.023871_rkx,  0.560178_rkx,  1.273434_rkx,   &
                12.021467_rkx, 62.583231_rkx, 63.593761_rkx,   &
                76.438310_rkx,  4.280910_rkx, 13.218362_rkx,   &
                17.818769_rkx,  8.359495_rkx, 56.792707_rkx,   &
                 8.448301_rkx,  1.978796_rkx,  8.863925_rkx,   &
                 0.186365_rkx,  8.996212_rkx,  6.771027_rkx,   &
                45.815258_rkx, 12.002811_rkx, 75.278220_rkx,   &
                65.241008_rkx, 18.870667_rkx, 22.009553_rkx,   &
                64.604291_rkx, 11.498094_rkx,  0.578834_rkx,   &
                 9.237738_rkx, 49.747842_rkx,  2.147012_rkx,   &
                 1.196895_rkx,  2.133898_rkx,  0.173168_rkx/)

    ! phases for mvelp sine series
    real(rkx) , parameter :: mvphas(pmvelen) = &
          (/    251.9025_rkx, 280.8325_rkx, 128.3057_rkx,   &
                348.1074_rkx, 292.7252_rkx, 165.1686_rkx,   &
                263.7951_rkx,  15.3747_rkx,  58.5749_rkx,   &
                 40.8226_rkx, 308.4258_rkx, 240.0099_rkx,   &
                222.9725_rkx, 106.5937_rkx, 114.5182_rkx,   &
                268.7809_rkx, 279.6869_rkx,  39.6448_rkx,   &
                126.4108_rkx, 291.5795_rkx, 307.2848_rkx,   &
                 18.9300_rkx, 273.7596_rkx, 143.8050_rkx,   &
                191.8927_rkx, 125.5237_rkx, 172.7351_rkx,   &
                316.7998_rkx, 319.6024_rkx,  69.7526_rkx,   &
                123.5968_rkx, 217.6432_rkx,  85.5882_rkx,   &
                156.2147_rkx,  66.9489_rkx,  20.2082_rkx,   &
                250.7568_rkx,  48.0188_rkx,   8.3739_rkx,   &
                 17.0374_rkx, 155.3409_rkx,  94.1709_rkx,   &
                221.1120_rkx,  28.9300_rkx, 117.1498_rkx,   &
                320.5095_rkx, 262.3602_rkx, 336.2148_rkx,   &
                233.0046_rkx, 155.6977_rkx, 184.6277_rkx,   &
                267.2772_rkx,  78.9281_rkx, 123.4722_rkx,   &
                188.7132_rkx, 180.1364_rkx,  49.1382_rkx,   &
                152.5268_rkx,  98.2198_rkx,  97.4808_rkx,   &
                221.5376_rkx, 168.2438_rkx, 161.1199_rkx,   &
                 55.0196_rkx, 262.6495_rkx, 200.3284_rkx,   &
                201.6651_rkx, 294.6547_rkx,  99.8233_rkx,   &
                213.5577_rkx, 154.1631_rkx, 232.7153_rkx,   &
                138.3034_rkx, 204.6609_rkx, 106.5938_rkx,   &
                250.4676_rkx, 332.3345_rkx,  27.3039_rkx/)

    integer(ik4) :: i       ! Index for series summations
    real(rkx) :: obsum   ! Obliquity series summation
    real(rkx) :: cossum  ! Cos series summation for eccentricity/fvelp
    real(rkx) :: sinsum  ! Sin series summation for eccentricity/fvelp
    real(rkx) :: fvelp   ! Fixed vernal equinox long of perihelion
    real(rkx) :: mvsum   ! mvelp series summation
    real(rkx) :: beta    ! Intermediate argument for lambm0
    real(rkx) :: years   ! Years to time of interest ( pos <=> future)
    real(rkx) :: eccen2  ! eccentricity squared
    real(rkx) :: eccen3  ! eccentricity cubed

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

      yb4_1950AD = 1950.0_rkx - real(iyear_AD,rkx)
      if ( abs(yb4_1950AD) > 1000000.0_rkx )then
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
      obliq = 23.320556_rkx + obsum

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

      if (abs(cossum) <= 1.0e-8_rkx) then
        if (sinsum == d_zero) then
          fvelp = d_zero
        else if (sinsum < d_zero) then
          fvelp = 1.5_rkx*mathpi
        else if (sinsum > d_zero) then
          fvelp = d_half*mathpi
        end if
      else if (cossum < d_zero) then
        fvelp = atan(sinsum/cossum) + mathpi
      else if (cossum > d_zero) then
        if (sinsum < d_zero) then
          fvelp = atan(sinsum/cossum) + 2.0_rkx*mathpi
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
      mvelp = fvelp/degrad + 50.439273_rkx*psecdeg*years + 3.392506_rkx + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp < d_zero)
        mvelp = mvelp + 360.0_rkx
      end do
      do while (mvelp >= 360.0_rkx)
        mvelp = mvelp - 360.0_rkx
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

    mvelpp = (mvelp + 180._rkx)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = sqrt(d_one - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = d_two*((d_half*eccen+.125_rkx*eccen3) * &
            (d_one+beta)*sin(mvelpp) - &
              .250_rkx*eccen2*(d_half+beta)*sin(d_two*mvelpp) + &
              .125_rkx*eccen3*(d_one/3._rkx + beta)*sin(3._rkx*mvelpp))

  end subroutine orb_params
  !
  ! Compute earth/orbit parameters using formula suggested by
  ! Duane Thresher.
  !
  subroutine orb_decl(calday,eccen,mvelpp,lambm0,obliqr,delta,eccf)
    implicit none
    real(rkx) , intent(in) :: calday ! Calendar day, including fraction
    real(rkx) , intent(in) :: eccen  ! Eccentricity
    real(rkx) , intent(in) :: obliqr ! Earths obliquity in radians
    real(rkx) , intent(in) :: lambm0 ! Mean long of perihelion at the
                                     ! vernal equinox (radians)
    real(rkx) , intent(in) :: mvelpp ! moving vernal equinox longitude
                                     ! of perihelion plus pi (radians)
    real(rkx) , intent(out) :: delta ! Solar declination angle in rad
    real(rkx) , intent(out) :: eccf  ! Earth-sun distance factor (ie. (1/r)**2)

    real(rkx) , parameter :: dayspy = 365.0_rkx  ! days per year
    real(rkx) , parameter :: ve     = 80.5_rkx   ! Calday of vernal equinox
                                              ! assumes Jan 1 = calday 1

    real(rkx) :: lambm  ! Lambda m, mean long of perihelion (rad)
    real(rkx) :: lmm    ! Intermediate argument involving lambm
    real(rkx) :: lamb   ! Lambda, the earths long of perihelion
    real(rkx) :: invrho ! Inverse normalized sun/earth distance
    real(rkx) :: sinl   ! Sine of lmm

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

    lambm = lambm0 + (calday - ve)*d_two*mathpi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(d_two*sinl + eccen*(1.25_rkx*sin(d_two*lmm)  &
          + eccen*((13.0_rkx/12.0_rkx)*sin(3._rkx*lmm) - 0.25_rkx*sinl)))

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
