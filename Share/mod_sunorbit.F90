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
  use mod_constants , only : mathpi , degrad

  implicit none

  private

  real(rk8) , parameter :: ORB_ECCEN_MIN  =   0.0D0 ! min value for eccen
  real(rk8) , parameter :: ORB_ECCEN_MAX  =   0.1D0 ! max value for eccen
  real(rk8) , parameter :: ORB_OBLIQ_MIN  = -90.0D0 ! min value for obliq
  real(rk8) , parameter :: ORB_OBLIQ_MAX  = +90.0D0 ! max value for obliq
  real(rk8) , parameter :: ORB_MVELP_MIN  =   0.0D0 ! min value for mvelp
  real(rk8) , parameter :: ORB_MVELP_MAX  = 360.0D0 ! max value for mvelp
  real(rk8) , parameter :: ORB_UNDEF_REAL = 1.D36
  integer(ik4) , parameter :: ORB_UNDEF_INT  = 2000000000

  public :: orb_cosz , orb_decl , orb_params

  contains
  !
  ! FUNCTION to return the cosine of the solar zenith angle.
  ! Assumes 365.0 days/year.
  !
  real(rk8) function orb_cosz(jday,lat,lon,declin)
    implicit none
    real(rk8) , intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real(rk8) , intent(in) :: lat    ! Centered latitude (radians)
    real(rk8) , intent(in) :: lon    ! Centered longitude (radians)
    real(rk8) , intent(in) :: declin ! Solar declination (radians)
    orb_cosz = sin(lat)*sin(declin) - &
               cos(lat)*cos(declin)*cos(jday*2.0D0*mathpi + lon)
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
    real(rk8) , parameter :: psecdeg = 1.0D0/3600.0D0

    real(rk8) :: yb4_1950AD         ! number of years before 1950 AD

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    ! amplitudes for obliquity cos series
    real(rk8) , parameter :: obamp(poblen) =  &
          (/   -2462.2214466D0, -857.3232075D0, -629.3231835D0,   &
                -414.2804924D0, -311.7632587D0,  308.9408604D0,   &
                -162.5533601D0, -116.1077911D0,  101.1189923D0,   &
                 -67.6856209D0,   24.9079067D0,   22.5811241D0,   &
                 -21.1648355D0,  -15.6549876D0,   15.3936813D0,   &
                  14.6660938D0,  -11.7273029D0,   10.2742696D0,   &
                   6.4914588D0,    5.8539148D0,   -5.4872205D0,   &
                  -5.4290191D0,    5.1609570D0,    5.0786314D0,   &
                  -4.0735782D0,    3.7227167D0,    3.3971932D0,   &
                  -2.8347004D0,   -2.6550721D0,   -2.5717867D0,   &
                  -2.4712188D0,    2.4625410D0,    2.2464112D0,   &
                  -2.0755511D0,   -1.9713669D0,   -1.8813061D0,   &
                  -1.8468785D0,    1.8186742D0,    1.7601888D0,   &
                  -1.5428851D0,    1.4738838D0,   -1.4593669D0,   &
                   1.4192259D0,   -1.1818980D0,    1.1756474D0,   &
                  -1.1316126D0,    1.0896928D0/)

    ! rates for obliquity cosine series
    real(rk8) , parameter :: obrate(poblen) = &
            (/  31.609974D0, 32.620504D0, 24.172203D0,   &
                31.983787D0, 44.828336D0, 30.973257D0,   &
                43.668246D0, 32.246691D0, 30.599444D0,   &
                42.681324D0, 43.836462D0, 47.439436D0,   &
                63.219948D0, 64.230478D0,  1.010530D0,   &
                 7.437771D0, 55.782177D0,  0.373813D0,   &
                13.218362D0, 62.583231D0, 63.593761D0,   &
                76.438310D0, 45.815258D0,  8.448301D0,   &
                56.792707D0, 49.747842D0, 12.058272D0,   &
                75.278220D0, 65.241008D0, 64.604291D0,   &
                 1.647247D0,  7.811584D0, 12.207832D0,   &
                63.856665D0, 56.155990D0, 77.448840D0,   &
                 6.801054D0, 62.209418D0, 20.656133D0,   &
                48.344406D0, 55.145460D0, 69.000539D0,   &
                11.071350D0, 74.291298D0, 11.047742D0,   &
                 0.636717D0, 12.844549D0/)

    ! phases for obliquity cosine series
    real(rk8) , parameter :: obphas(poblen) = &
          (/    251.9025D0, 280.8325D0, 128.3057D0,   &
                292.7252D0,  15.3747D0, 263.7951D0,   &
                308.4258D0, 240.0099D0, 222.9725D0,   &
                268.7809D0, 316.7998D0, 319.6024D0,   &
                143.8050D0, 172.7351D0,  28.9300D0,   &
                123.5968D0,  20.2082D0,  40.8226D0,   &
                123.4722D0, 155.6977D0, 184.6277D0,   &
                267.2772D0,  55.0196D0, 152.5268D0,   &
                 49.1382D0, 204.6609D0,  56.5233D0,   &
                200.3284D0, 201.6651D0, 213.5577D0,   &
                 17.0374D0, 164.4194D0,  94.5422D0,   &
                131.9124D0,  61.0309D0, 296.2073D0,   &
                135.4894D0, 114.8750D0, 247.0691D0,   &
                256.6114D0,  32.1008D0, 143.6804D0,   &
                 16.8784D0, 160.6835D0,  27.5932D0,   &
                348.1074D0,  82.6496D0/)

    ! Cosine/sine series data for computation of eccentricity and fixed vernal
    ! equinox longitude of perihelion (fvelp): amplitude,
    ! rate (arc seconds/year), phase (degrees).

    ! ampl for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecamp (pecclen) = &
          (/   0.01860798D0,  0.01627522D0, -0.01300660D0,   &
               0.00988829D0, -0.00336700D0,  0.00333077D0,   &
              -0.00235400D0,  0.00140015D0,  0.00100700D0,   &
               0.00085700D0,  0.00064990D0,  0.00059900D0,   &
               0.00037800D0, -0.00033700D0,  0.00027600D0,   &
               0.00018200D0, -0.00017400D0, -0.00012400D0,   &
               0.00001250D0/)

    ! rates for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecrate(pecclen) = &
          (/    4.2072050D0,  7.3460910D0, 17.8572630D0,  &
               17.2205460D0, 16.8467330D0,  5.1990790D0,  &
               18.2310760D0, 26.2167580D0,  6.3591690D0,  &
               16.2100160D0,  3.0651810D0, 16.5838290D0,  &
               18.4939800D0,  6.1909530D0, 18.8677930D0,  &
               17.4255670D0,  6.1860010D0, 18.4174410D0,  &
                0.6678630D0/)

    ! phases for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecphas(pecclen) = &
          (/    28.620089D0, 193.788772D0, 308.307024D0,  &
               320.199637D0, 279.376984D0,  87.195000D0,  &
               349.129677D0, 128.443387D0, 154.143880D0,  &
               291.269597D0, 114.860583D0, 332.092251D0,  &
               296.414411D0, 145.769910D0, 337.237063D0,  &
               152.092288D0, 126.839891D0, 210.667199D0,  &
                72.108838D0/)

    ! Sine series data for computation of moving vernal equinox longitude of
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    ! amplitudes for mvelp sine series
    real(rk8) , parameter :: mvamp (pmvelen) = &
          (/   7391.0225890D0, 2555.1526947D0, 2022.7629188D0,  &
              -1973.6517951D0, 1240.2321818D0,  953.8679112D0,  &
               -931.7537108D0,  872.3795383D0,  606.3544732D0,  &
               -496.0274038D0,  456.9608039D0,  346.9462320D0,  &
               -305.8412902D0,  249.6173246D0, -199.1027200D0,  &
                191.0560889D0, -175.2936572D0,  165.9068833D0,  &
                161.1285917D0,  139.7878093D0, -133.5228399D0,  &
                117.0673811D0,  104.6907281D0,   95.3227476D0,  &
                 86.7824524D0,   86.0857729D0,   70.5893698D0,  &
                -69.9719343D0,  -62.5817473D0,   61.5450059D0,  &
                -57.9364011D0,   57.1899832D0,  -57.0236109D0,  &
                -54.2119253D0,   53.2834147D0,   52.1223575D0,  &
                -49.0059908D0,  -48.3118757D0,  -45.4191685D0,  &
                -42.2357920D0,  -34.7971099D0,   34.4623613D0,  &
                -33.8356643D0,   33.6689362D0,  -31.2521586D0,  &
                -30.8798701D0,   28.4640769D0,  -27.1960802D0,  &
                 27.0860736D0,  -26.3437456D0,   24.7253740D0,  &
                 24.6732126D0,   24.4272733D0,   24.0127327D0,  &
                 21.7150294D0,  -21.5375347D0,   18.1148363D0,  &
                -16.9603104D0,  -16.1765215D0,   15.5567653D0,  &
                 15.4846529D0,   15.2150632D0,   14.5047426D0,  &
                -14.3873316D0,   13.1351419D0,   12.8776311D0,  &
                 11.9867234D0,   11.9385578D0,   11.7030822D0,  &
                 11.6018181D0,  -11.2617293D0,  -10.4664199D0,  &
                 10.4333970D0,  -10.2377466D0,   10.1934446D0,  &
                -10.1280191D0,   10.0289441D0,  -10.0034259D0/)

    ! rates for mvelp sine series
    real(rk8) , parameter :: mvrate(pmvelen) = &
          (/    31.609974D0, 32.620504D0, 24.172203D0,   &
                 0.636717D0, 31.983787D0,  3.138886D0,   &
                30.973257D0, 44.828336D0,  0.991874D0,   &
                 0.373813D0, 43.668246D0, 32.246691D0,   &
                30.599444D0,  2.147012D0, 10.511172D0,   &
                42.681324D0, 13.650058D0,  0.986922D0,   &
                 9.874455D0, 13.013341D0,  0.262904D0,   &
                 0.004952D0,  1.142024D0, 63.219948D0,   &
                 0.205021D0,  2.151964D0, 64.230478D0,   &
                43.836462D0, 47.439436D0,  1.384343D0,   &
                 7.437771D0, 18.829299D0,  9.500642D0,   &
                 0.431696D0,  1.160090D0, 55.782177D0,   &
                12.639528D0,  1.155138D0,  0.168216D0,   &
                 1.647247D0, 10.884985D0,  5.610937D0,   &
                12.658184D0,  1.010530D0,  1.983748D0,   &
                14.023871D0,  0.560178D0,  1.273434D0,   &
                12.021467D0, 62.583231D0, 63.593761D0,   &
                76.438310D0,  4.280910D0, 13.218362D0,   &
                17.818769D0,  8.359495D0, 56.792707D0,   &
                8.448301D0,  1.978796D0,  8.863925D0,   &
                 0.186365D0,  8.996212D0,  6.771027D0,   &
                45.815258D0, 12.002811D0, 75.278220D0,   &
                65.241008D0, 18.870667D0, 22.009553D0,   &
                64.604291D0, 11.498094D0,  0.578834D0,   &
                 9.237738D0, 49.747842D0,  2.147012D0,   &
                 1.196895D0,  2.133898D0,  0.173168D0/)

    ! phases for mvelp sine series
    real(rk8) , parameter :: mvphas(pmvelen) = &
          (/    251.9025D0, 280.8325D0, 128.3057D0,   &
                348.1074D0, 292.7252D0, 165.1686D0,   &
                263.7951D0,  15.3747D0,  58.5749D0,   &
                 40.8226D0, 308.4258D0, 240.0099D0,   &
                222.9725D0, 106.5937D0, 114.5182D0,   &
                268.7809D0, 279.6869D0,  39.6448D0,   &
                126.4108D0, 291.5795D0, 307.2848D0,   &
                 18.9300D0, 273.7596D0, 143.8050D0,   &
                191.8927D0, 125.5237D0, 172.7351D0,   &
                316.7998D0, 319.6024D0,  69.7526D0,   &
                123.5968D0, 217.6432D0,  85.5882D0,   &
                156.2147D0,  66.9489D0,  20.2082D0,   &
                250.7568D0,  48.0188D0,   8.3739D0,   &
                 17.0374D0, 155.3409D0,  94.1709D0,   &
                221.1120D0,  28.9300D0, 117.1498D0,   &
                320.5095D0, 262.3602D0, 336.2148D0,   &
                233.0046D0, 155.6977D0, 184.6277D0,   &
                267.2772D0,  78.9281D0, 123.4722D0,   &
                188.7132D0, 180.1364D0,  49.1382D0,   &
                152.5268D0,  98.2198D0,  97.4808D0,   &
                221.5376D0, 168.2438D0, 161.1199D0,   &
                 55.0196D0, 262.6495D0, 200.3284D0,   &
                201.6651D0, 294.6547D0,  99.8233D0,   &
                213.5577D0, 154.1631D0, 232.7153D0,   &
                138.3034D0, 204.6609D0, 106.5938D0,   &
                250.4676D0, 332.3345D0,  27.3039D0/)

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

      yb4_1950AD = 1950.0D0 - real(iyear_AD,rk8)
      if ( abs(yb4_1950AD) .gt. 1000000.0D0 )then
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

      obsum = 0.0D0
      do i = 1 , poblen
        obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
                obphas(i))*degrad)
      end do
      obliq = 23.320556D0 + obsum

      ! Summation of cosine and sine series for computation of eccentricity
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation
      ! of moving vernal equinox longitude of perihelion.  Convert the rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.

      cossum = 0.0D0
      do i = 1 , pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      sinsum = 0.0D0
      do i = 1 , pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      ! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen

      ! A series of cases for fvelp, which is in radians.

      if (abs(cossum) .le. 1.0D-8) then
        if (sinsum .eq. 0.0D0) then
          fvelp = 0.0D0
        else if (sinsum .lt. 0.0D0) then
          fvelp = 1.5D0*mathpi
        else if (sinsum .gt. 0.0D0) then
          fvelp = .5D0*mathpi
        end if
      else if (cossum .lt. 0.0D0) then
        fvelp = atan(sinsum/cossum) + mathpi
      else if (cossum .gt. 0.0D0) then
        if (sinsum .lt. 0.0D0) then
          fvelp = atan(sinsum/cossum) + 2.0D0*mathpi
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

      mvsum = 0.0D0
      do i = 1 , pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
                mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273D0*psecdeg*years + 3.392506D0 + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp .lt. 0.0D0)
        mvelp = mvelp + 360.0D0
      end do
      do while (mvelp .ge. 360.0D0)
        mvelp = mvelp - 360.0D0
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

    mvelpp = (mvelp + 180.D0)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = sqrt(1.D0 - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = 2.D0*((.5D0*eccen + .125D0*eccen3)*(1.D0 + beta)*sin(mvelpp)  &
           - .250D0*eccen2*(.5D0    + beta)*sin(2.D0*mvelpp)               &
           + .125D0*eccen3*(1.D0/3.D0 + beta)*sin(3.D0*mvelpp))

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

    real(rk8) , parameter :: dayspy = 365.0D0  ! days per year
    real(rk8) , parameter :: ve     = 80.5D0   ! Calday of vernal equinox
                                              ! assumes Jan 1 = calday 1

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

    lambm = lambm0 + (calday - ve)*2.D0*mathpi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(2.D0*sinl + eccen*(1.25D0*sin(2.D0*lmm)  &
          + eccen*((13.0D0/12.0D0)*sin(3.D0*lmm) - 0.25D0*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually
    ! inverse rho will be used), and thus the eccentricity factor (eccf),
    ! can be calculated from formulas given in Berger 1978.

    invrho = (1.D0 + eccen*cos(lamb - mvelpp)) / (1.D0 - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta  = asin(sin(obliqr)*sin(lamb))
    eccf   = invrho*invrho
  end subroutine orb_decl

end module mod_sunorbit
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
