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

module mod_rad_tracer

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam

  implicit none

  private

  public :: trcmix , trcpth , trcab , trcabn , trcems , trcplk

  real(rkx) , public :: cfc110 , cfc120 , ch40 , co2mmr , n2o0

  contains
  !
  !-----------------------------------------------------------------------
  !
  ! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
  ! CFC12
  !          Code: J.T.Kiehl November 21, 1994
  !
  !-----------------------------------------------------------------------
  !
  !------------------------------input------------------------------------
  !
  ! pmid   - model pressures
  !
  !------------------------------output-----------------------------------
  !
  ! n2o    - nitrous oxide mass mixing ratio
  ! ch4    - methane mass mixing ratio
  ! cfc11  - cfc11 mass mixing ratio
  ! cfc12  - cfc12 mass mixing ratio
  !
  !-----------------------------------------------------------------------
  !
  subroutine trcmix(n1,n2,dlat,xptrop,pmid,n2o,ch4,cfc11,cfc12)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , pointer , dimension(:) , intent(in) :: dlat , xptrop
    real(rkx) , pointer , dimension(:,:) , intent(in) :: pmid
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: cfc11
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: cfc12
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ch4
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: n2o
    !
    ! dlat   - absolute value of latitude in degrees
    ! xn2o   - pressure scale height for n2o
    ! xch4   - pressure scale height for ch4
    ! xcfc11 - pressure scale height for cfc11
    ! xcfc12 - pressure scale height for cfc12
    ! ptrop  - pressure level of tropopause
    ! pratio - pressure divided by ptrop
    !
    real(rkx) :: pratio , xcfc11 , xcfc12 , xch4 , xn2o , alat
    integer(ik4) :: j , k

    xcfc11 = d_zero
    xcfc12 = d_zero
    xch4 = d_zero
    xn2o = d_zero
    do j = n1 , n2
      alat = dlat(j) ! This is absolute value of latitude in degrees
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
      !  set stratospheric scale height factor for gases
      do k = 1 , kz
        if ( pmid(j,k) >= xptrop(j) ) then
          ch4(j,k) = ch40
          n2o(j,k) = n2o0
          cfc11(j,k) = cfc110
          cfc12(j,k) = cfc120
        else
          pratio = pmid(j,k)/xptrop(j)
          ch4(j,k) = ch40*(pratio**xch4)
          n2o(j,k) = n2o0*(pratio**xn2o)
          cfc11(j,k) = cfc110*(pratio**xcfc11)
          cfc12(j,k) = cfc120*(pratio**xcfc12)
        end if
      end do
    end do
  end subroutine trcmix
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
  subroutine trcpth(n1,n2,tnm,pnm,cfc11,cfc12,n2o,ch4,qnm,  &
                    ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212, &
                    uco213,uco221,uco222,uco223,bn2o0,bn2o1,bch4, &
                    uptype)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , pointer , dimension(:,:) :: bch4 , bn2o0 , bn2o1 , pnm ,  &
                   ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,  &
                   uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
    real(rkx) , pointer , dimension(:,:) :: cfc11 , cfc12 , ch4 , &
                                            n2o , qnm , tnm
    intent (in) cfc11 , cfc12 , ch4 , n2o , pnm , qnm , tnm
    intent (inout) bch4 , bn2o0 , bn2o1 , ucfc11 , ucfc12 , uch4 ,    &
                   uco211 , uco212 , uco213 , uco221 , uco222 ,       &
                   uco223 , un2o0 , un2o1 , uptype
    !
    !   co2fac - co2 factor
    !   alpha1 - stimulated emission term
    !   alpha2 - stimulated emission term
    !   rt     - reciprocal of local temperature
    !   rsqrt  - reciprocal of sqrt of temp
    !   pbar   - mean pressure
    !   dpnm   - difference in pressure
    !
    real(rkx) :: diff , alpha1 , alpha2 , dpnm , pbar , rsqrt , rt , co2fac
    integer(ik4) :: j , k
    data diff/1.66_rkx/           ! diffusivity factor

    !-----------------------------------------------------------------------
    !   Calculate path lengths for the trace gases
    !-----------------------------------------------------------------------
    do j = n1 , n2
      ucfc11(j,1) = 1.8_rkx*cfc11(j,1)*pnm(j,1)*regravgts
      ucfc12(j,1) = 1.8_rkx*cfc12(j,1)*pnm(j,1)*regravgts
      un2o0(j,1) = diff*1.02346e5_rkx*n2o(j,1)*pnm(j,1)*regravgts/sqrt(tnm(j,1))
      un2o1(j,1) = diff*2.01909_rkx*un2o0(j,1)*exp(-847.36_rkx/tnm(j,1))
      uch4(j,1) = diff*8.60957e4_rkx*ch4(j,1)*pnm(j,1)* &
                  regravgts/sqrt(tnm(j,1))
      co2fac = diff*co2mmr*pnm(j,1)*regravgts
      alpha1 = (d_one-exp(-1540.0_rkx/tnm(j,1)))**3/sqrt(tnm(j,1))
      alpha2 = (d_one-exp(-1360.0_rkx/tnm(j,1)))**3/sqrt(tnm(j,1))
      uco211(j,1) = 3.42217e3_rkx*co2fac*alpha1*exp(-1849.7_rkx/tnm(j,1))
      uco212(j,1) = 6.02454e3_rkx*co2fac*alpha1*exp(-2782.1_rkx/tnm(j,1))
      uco213(j,1) = 5.53143e3_rkx*co2fac*alpha1*exp(-3723.2_rkx/tnm(j,1))
      uco221(j,1) = 3.88984e3_rkx*co2fac*alpha2*exp(-1997.6_rkx/tnm(j,1))
      uco222(j,1) = 3.67108e3_rkx*co2fac*alpha2*exp(-3843.8_rkx/tnm(j,1))
      uco223(j,1) = 6.50642e3_rkx*co2fac*alpha2*exp(-2989.7_rkx/tnm(j,1))
      bn2o0(j,1) = diff*19.399_rkx*pnm(j,1)**2*n2o(j,1) * &
                   1.02346e5_rkx*regravgts/(sslp*tnm(j,1))
      bn2o1(j,1) = bn2o0(j,1)*exp(-847.36_rkx/tnm(j,1))*2.06646e5_rkx
      bch4(j,1) = diff*2.94449_rkx*ch4(j,1)*pnm(j,1)**2*regravgts * &
                  8.60957e4_rkx/(sslp*tnm(j,1))
      uptype(j,1) = diff*qnm(j,1)*pnm(j,1)**2*exp(1800.0_rkx* &
                   (d_one/tnm(j,1)-d_one/296.0_rkx))*regravgts/sslp
    end do
    do k = 1 , kz
      do j = n1 , n2
        rt = d_one/tnm(j,k)
        rsqrt = sqrt(rt)
        pbar = ((pnm(j,k+1)+pnm(j,k))*d_half)/sslp
        dpnm = (pnm(j,k+1)-pnm(j,k))*regravgts
        alpha1 = diff*rsqrt*(d_one-exp(-1540.0_rkx/tnm(j,k)))**3
        alpha2 = diff*rsqrt*(d_one-exp(-1360.0_rkx/tnm(j,k)))**3
        ucfc11(j,k+1) = ucfc11(j,k) + 1.8_rkx*cfc11(j,k)*dpnm
        ucfc12(j,k+1) = ucfc12(j,k) + 1.8_rkx*cfc12(j,k)*dpnm
        un2o0(j,k+1) = un2o0(j,k) + diff*1.02346e5_rkx*n2o(j,k)*rsqrt*dpnm
        un2o1(j,k+1) = un2o1(j,k) + diff*2.06646e5_rkx*n2o(j,k) * &
                       rsqrt*exp(-847.36_rkx/tnm(j,k))*dpnm
        uch4(j,k+1) = uch4(j,k) + diff*8.60957e4_rkx*ch4(j,k)*rsqrt*dpnm
        uco211(j,k+1) = uco211(j,k) + 1.15_rkx*3.42217e3_rkx*alpha1 * &
                        co2mmr*exp(-1849.7_rkx/tnm(j,k))*dpnm
        uco212(j,k+1) = uco212(j,k) + 1.15_rkx*6.02454e3_rkx*alpha1 * &
                        co2mmr*exp(-2782.1_rkx/tnm(j,k))*dpnm
        uco213(j,k+1) = uco213(j,k) + 1.15_rkx*5.53143e3_rkx*alpha1 * &
                        co2mmr*exp(-3723.2_rkx/tnm(j,k))*dpnm
        uco221(j,k+1) = uco221(j,k) + 1.15_rkx*3.88984e3_rkx*alpha2 * &
                        co2mmr*exp(-1997.6_rkx/tnm(j,k))*dpnm
        uco222(j,k+1) = uco222(j,k) + 1.15_rkx*3.67108e3_rkx*alpha2 * &
                        co2mmr*exp(-3843.8_rkx/tnm(j,k))*dpnm
        uco223(j,k+1) = uco223(j,k) + 1.15_rkx*6.50642e3_rkx*alpha2 * &
                        co2mmr*exp(-2989.7_rkx/tnm(j,k))*dpnm
        bn2o0(j,k+1) = bn2o0(j,k) + diff*19.399_rkx*pbar*rt * &
                       1.02346e5_rkx*n2o(j,k)*dpnm
        bn2o1(j,k+1) = bn2o1(j,k) + diff*19.399_rkx*pbar*rt * &
                       2.06646e5_rkx*exp(-847.36_rkx/tnm(j,k))*n2o(j,k)*dpnm
        bch4(j,k+1) = bch4(j,k) + diff*2.94449_rkx*rt*pbar * &
                      8.60957e4_rkx*ch4(j,k)*dpnm
        uptype(j,k+1) = uptype(j,k) + diff*qnm(j,k)*exp(1800.0_rkx*(d_one / &
                        tnm(j,k)-d_one/296.0_rkx))*pbar*dpnm
      end do
    end do
  end subroutine trcpth
  !
  !----------------------------------------------------------------------
  ! Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
  ! CFC12.
  !
  !             Coded by J.T. Kiehl November 21, 1994
  !-----------------------------------------------------------------------
  !
  !------------------------------Arguments--------------------------------
  ! to3co2  - pressure weighted temperature
  ! pnm     - interface pressures
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
  ! dw      - h2o path length
  ! pnew    - pressure
  ! s2c     - continuum path length
  ! uptype  - p-type h2o path length
  ! dplh2o  - p squared h2o path length
  ! abplnk1 - Planck factor
  ! tco2    - co2 transmission factor
  ! th2o    - h2o transmission factor
  ! to3     - o3 transmission factor
  !
  !------------------------------Output Arguments-------------------------
  !
  ! abstrc  - total trace gas absorptivity
  !
  !-----------------------------------------------------------------------
  !
  subroutine trcab(n1,n2,k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4, &
                   uco211,uco212,uco213,uco221,uco222,uco223,bn2o0,  &
                   bn2o1,bch4,to3co2,pnm,dw,pnew,s2c,uptype,dplh2o,  &
                   abplnk1,tco2,th2o,to3,abstrc)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2 , k1 , k2
    real(rkx) , pointer , dimension(:,:,:) :: abplnk1
    real(rkx) , pointer , dimension(:) :: abstrc , dplh2o , dw , pnew , tco2 ,&
                                        th2o , to3 , to3co2
    real(rkx) , pointer , dimension(:,:) :: bch4 , bn2o0 , bn2o1 , pnm ,  &
             s2c , ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,  &
             uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
    intent (in) abplnk1 , bch4 , bn2o0 , bn2o1 , dplh2o , dw ,      &
                pnew , pnm , s2c , tco2 , th2o , to3 , to3co2 ,     &
                ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 , &
                uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
    intent (inout) abstrc
    !
    !-----------------------------------------------------------------------
    !
    ! sqti    - square root of mean temp
    ! du1     - cfc11 path length
    ! du2     - cfc12 path length
    ! acfc1   - cfc11 absorptivity 798 cm-1
    ! acfc2   - cfc11 absorptivity 846 cm-1
    ! acfc3   - cfc11 absorptivity 933 cm-1
    ! acfc4   - cfc11 absorptivity 1085 cm-1
    ! acfc5   - cfc12 absorptivity 889 cm-1
    ! acfc6   - cfc12 absorptivity 923 cm-1
    ! acfc7   - cfc12 absorptivity 1102 cm-1
    ! acfc8   - cfc12 absorptivity 1161 cm-1
    ! du01    - n2o path length
    ! dbeta01 - n2o pressure factor
    ! dbeta11 -        "
    ! an2o1   - absorptivity of 1285 cm-1 n2o band
    ! du02    - n2o path length
    ! dbeta02 - n2o pressure factor
    ! an2o2   - absorptivity of 589 cm-1 n2o band
    ! du03    - n2o path length
    ! dbeta03 - n2o pressure factor
    ! an2o3   - absorptivity of 1168 cm-1 n2o band
    ! duch4   - ch4 path length
    ! dbetac  - ch4 pressure factor
    ! ach4    - absorptivity of 1306 cm-1 ch4 band
    ! du11    - co2 path length
    ! du12    -       "
    ! du13    -       "
    ! dbetc1  - co2 pressure factor
    ! dbetc2  - co2 pressure factor
    ! aco21   - absorptivity of 1064 cm-1 band
    ! du21    - co2 path length
    ! du22    -       "
    ! du23    -       "
    ! aco22   - absorptivity of 961 cm-1 band
    ! tt      - temp. factor for h2o overlap factor
    ! psi1    -                 "
    ! phi1    -                 "
    ! p1      - h2o overlap factor
    ! w1      -        "
    ! ds2c    - continuum path length
    ! duptyp  - p-type path length
    ! tw      - h2o transmission factor
    ! g1      -         "
    ! g2      -         "
    ! g3      -         "
    ! g4      -         "
    ! ab      - h2o temp. factor
    ! bb      -         "
    ! abp     -         "
    ! bbp     -         "
    ! tcfc3   - transmission for cfc11 band
    ! tcfc4   - transmission for cfc11 band
    ! tcfc6   - transmission for cfc12 band
    ! tcfc7   - transmission for cfc12 band
    ! tcfc8   - transmission for cfc12 band
    ! tlw     - h2o transmission
    ! tch4    - ch4 transmission
    !
    !-----------------------------------------------------------------------
    !
    real(rkx) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4 , tw
    real(rkx) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 ,&
               acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , &
               dbeta01 , dbeta02 , dbeta03 , dbeta11 , dbetac ,       &
               dbetc1 , dbetc2 , du01 , du02 , du03 , du1 , du11 ,    &
               du12 , du13 , du2 , du21 , du22 , du23 , duch4 , p1 ,  &
               phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 ,  &
               tch4 , tlw , w1 , sqti , ds2c , duptyp , tt
    integer(ik4) :: j , l

    data g1 / 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
              0.0304380_rkx , 0.0540398_rkx , 0.0321962_rkx /
    data g2 / 14.48320_rkx , 4.302420_rkx ,  5.23523_rkx , &
               3.25342_rkx , 0.698935_rkx , 16.55990_rkx /
    data g3 / 26.18980_rkx , 18.44760_rkx , 15.36330_rkx , &
              12.19270_rkx ,  9.14992_rkx ,  8.07092_rkx /
    data g4 / 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
              0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx /
    data ab / 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
              2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx /
    data bb /-1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
             -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx /
    data abp / 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
               2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx/
    data bbp /-1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
              -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx /

    do j = n1 , n2
      sqti = sqrt(to3co2(j))
      ! h2o transmission
      tt = abs(to3co2(j)-250.0_rkx)
      ds2c = abs(s2c(j,k1)-s2c(j,k2))
      duptyp = abs(uptype(j,k1)-uptype(j,k2))
      do l = 1 , 6
        psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
        phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
        p1 = pnew(j)*(psi1/phi1)/sslp
        w1 = dw(j)*phi1
        tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1)) - &
                     d_one)-g3(l)*ds2c-g4(l)*duptyp)
      end do
      du1 = abs(ucfc11(j,k1)-ucfc11(j,k2))
      du2 = abs(ucfc12(j,k1)-ucfc12(j,k2))
      ! cfc transmissions
      tcfc3 = exp(-175.005_rkx*du1)
      tcfc4 = exp(-1202.18_rkx*du1)
      tcfc6 = exp(-5786.73_rkx*du2)
      tcfc7 = exp(-2873.51_rkx*du2)
      tcfc8 = exp(-2085.59_rkx*du2)
      ! Absorptivity for CFC11 bands
      acfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*du1))*tw(1)*abplnk1(7,j,k2)
      acfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*du1))*tw(2)*abplnk1(8,j,k2)
      acfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*abplnk1(9,j,k2)
      acfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*abplnk1(10,j,k2)
      ! Absorptivity for CFC12 bands
      acfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*du2))*tw(3)*abplnk1(11,j,k2)
      acfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*abplnk1(12,j,k2)
      acfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*abplnk1(13,j,k2)
      acfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*abplnk1(14,j,k2)
      ! Emissivity for CH4 band 1306 cm-1
      tlw = exp(-d_one*sqrt(dplh2o(j)))
      duch4 = abs(uch4(j,k1)-uch4(j,k2))
      dbetac = abs(bch4(j,k1)-bch4(j,k2))/duch4
      ach4 = 6.00444_rkx*sqti*log(d_one+func(duch4,dbetac)) * &
             tlw*abplnk1(3,j,k2)
      tch4 = d_one/(d_one+0.02_rkx*func(duch4,dbetac))
      ! Absorptivity for N2O bands
      du01 = abs(un2o0(j,k1)-un2o0(j,k2))
      du11 = abs(un2o1(j,k1)-un2o1(j,k2))
      dbeta01 = abs(bn2o0(j,k1)-bn2o0(j,k2))/du01
      dbeta11 = abs(bn2o1(j,k1)-bn2o1(j,k2))/du11
      ! 1285 cm-1 band
      an2o1 = 2.35558_rkx*sqti * &
             log(d_one+func(du01,dbeta01)+func(du11,dbeta11)) * &
             tlw*tch4*abplnk1(4,j,k2)
      du02 = 0.100090_rkx*du01
      du12 = 0.0992746_rkx*du11
      dbeta02 = 0.964282_rkx*dbeta01
      ! 589 cm-1 band
      an2o2 = 2.65581_rkx*sqti * &
              log(d_one+func(du02,dbeta02) + &
              func(du12,dbeta02))*th2o(j)*tco2(j)*abplnk1(5,j,k2)
      du03 = 0.0333767_rkx*du01
      dbeta03 = 0.982143_rkx*dbeta01
      ! 1168 cm-1 band
      an2o3 = 2.54034_rkx*sqti*log(d_one+func(du03,dbeta03)) * &
              tw(6)*tcfc8*abplnk1(6,j,k2)
      ! Emissivity for 1064 cm-1 band of CO2
      du11 = abs(uco211(j,k1)-uco211(j,k2))
      du12 = abs(uco212(j,k1)-uco212(j,k2))
      du13 = abs(uco213(j,k1)-uco213(j,k2))
      dbetc1 = 2.97558_rkx*abs(pnm(j,k1)+pnm(j,k2))/(d_two*sslp*sqti)
      dbetc2 = d_two*dbetc1
      aco21 = 3.7571_rkx*sqti * &
              log(d_one+func(du11,dbetc1)+func(du12,dbetc2) + &
              func(du13,dbetc2))*to3(j)*tw(5)*tcfc4*tcfc7*abplnk1(2,j,k2)
      ! Emissivity for 961 cm-1 band
      du21 = abs(uco221(j,k1)-uco221(j,k2))
      du22 = abs(uco222(j,k1)-uco222(j,k2))
      du23 = abs(uco223(j,k1)-uco223(j,k2))
      aco22 = 3.8443_rkx*sqti * &
              log(d_one+func(du21,dbetc1)+func(du22,dbetc1) + &
              func(du23,dbetc2))*tw(4)*tcfc3*tcfc6*abplnk1(1,j,k2)
      ! total trace gas absorptivity
      abstrc(j) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 +  &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +   &
                  aco21 + aco22
    end do
  end subroutine trcab
  !
  !----------------------------------------------------------------------
  ! Calculate nearest layer absorptivity due to CH4, N2O, CFC11 and CFC12
  !
  !         Coded by J.T. Kiehl November 21, 1994
  !-----------------------------------------------------------------------
  !
  !------------------------------Arguments--------------------------------
  !
  ! tbar   - pressure weighted temperature
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
  ! bplnk  - weighted Planck function for absorptivity
  ! winpl  - fractional path length
  ! pinpl  - pressure factor for subdivided layer
  ! tco2   - co2 transmission
  ! th2o   - h2o transmission
  ! to3    - o3 transmission
  ! dw     - h2o path length
  ! pnew   - pressure factor
  ! s2c    - h2o continuum factor
  ! uptype - p-type path length
  ! up2    - p squared path length
  ! uinpl  - Nearest layer subdivision factor
  !
  !------------------------------Output Arguments-------------------------
  !
  ! abstrc - total trace gas absorptivity
  !
  !-----------------------------------------------------------------------
  !
  subroutine trcabn(n1,n2,k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,  &
                    uco211,uco212,uco213,uco221,uco222,uco223,tbar,    &
                    bplnk,winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,up2, &
                    pnew,abstrc,uinpl)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2 , k2 , kn
    real(rkx) , pointer , dimension(:) :: abstrc , dw , pnew , tco2 , th2o ,  &
                                        to3 , up2
    real(rkx) , pointer , dimension(:,:,:) :: bplnk
    real(rkx) , pointer , dimension(:,:) :: pinpl , tbar , uinpl , winpl
    real(rkx) , pointer , dimension(:,:) :: s2c , ucfc11 , ucfc12 , uch4 , &
             uco211 , uco212 , uco213 , uco221 , uco222 , uco223 ,       &
             un2o0 , un2o1 , uptype
    intent (in) bplnk , dw , pinpl , pnew , s2c , tbar , tco2 , th2o ,    &
                to3 , ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 , &
                uco221 , uco222 , uco223 , uinpl , un2o0 , un2o1 , up2 ,  &
                uptype , winpl
    intent (inout) abstrc
    !
    ! sqti    - square root of mean temp
    ! rsqti   - reciprocal of sqti
    ! du1     - cfc11 path length
    ! du2     - cfc12 path length
    ! acfc1   - absorptivity of cfc11 798 cm-1 band
    ! acfc2   - absorptivity of cfc11 846 cm-1 band
    ! acfc3   - absorptivity of cfc11 933 cm-1 band
    ! acfc4   - absorptivity of cfc11 1085 cm-1 band
    ! acfc5   - absorptivity of cfc11 889 cm-1 band
    ! acfc6   - absorptivity of cfc11 923 cm-1 band
    ! acfc7   - absorptivity of cfc11 1102 cm-1 band
    ! acfc8   - absorptivity of cfc11 1161 cm-1 band
    ! du01    - n2o path length
    ! dbeta01 - n2o pressure factors
    ! dbeta11 -        "
    ! an2o1   - absorptivity of the 1285 cm-1 n2o band
    ! du02    - n2o path length
    ! dbeta02 - n2o pressure factor
    ! an2o2   - absorptivity of the 589 cm-1 n2o band
    ! du03    - n2o path length
    ! dbeta03 - n2o pressure factor
    ! an2o3   - absorptivity of the 1168 cm-1 n2o band
    ! duch4   - ch4 path length
    ! dbetac  - ch4 pressure factor
    ! ach4    - absorptivity of the 1306 cm-1 ch4 band
    ! du11    - co2 path length
    ! du12    -       "
    ! du13    -       "
    ! dbetc1 -  co2 pressure factor
    ! dbetc2 -  co2 pressure factor
    ! aco21  -  absorptivity of the 1064 cm-1 co2 band
    ! du21   -  co2 path length
    ! du22   -        "
    ! du23   -        "
    ! aco22  -  absorptivity of the 961 cm-1 co2 band
    ! tt     -  temp. factor for h2o overlap
    ! psi1   -           "
    ! phi1   -           "
    ! p1     -  factor for h2o overlap
    ! w1     -           "
    ! ds2c   -  continuum path length
    ! duptyp -  p-type path length
    ! tw     -  h2o transmission overlap
    ! g1     -  h2o overlap factor
    ! g2     -          "
    ! g3     -          "
    ! g4     -          "
    ! ab     -  h2o temp. factor
    ! bb     -          "
    ! abp    -          "
    ! bbp    -          "
    ! tcfc3  -  transmission of cfc11 band
    ! tcfc4  -  transmission of cfc11 band
    ! tcfc6  -  transmission of cfc12 band
    ! tcfc7  -          "
    ! tcfc8  -          "
    ! tlw    -  h2o transmission
    ! tch4   -  ch4 transmission
    !
    real(rkx) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4 , tw
    real(rkx) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 ,&
               acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , &
               dbeta01 , dbeta02 , dbeta03 , dbeta11 , dbetac ,       &
               dbetc1 , dbetc2 , du01 , du02 , du03 , du1 , du11 ,    &
               du12 , du13 , du2 , du21 , du22 , du23 , duch4 , p1 ,  &
               phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 ,  &
               tch4 , tlw , w1 , ds2c , duptyp , rsqti , sqti , tt
    integer(ik4) :: j , l

    data g1 / 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
              0.0304380_rkx , 0.0540398_rkx ,  0.0321962_rkx /
    data g2 / 14.483200_rkx ,  4.302420_rkx ,  5.2352300_rkx , &
               3.253420_rkx ,  0.698935_rkx , 16.5599000_rkx/
    data g3 / 26.18980_rkx , 18.44760_rkx , 15.36330_rkx , &
              12.19270_rkx ,  9.14992_rkx ,  8.07092_rkx /
    data g4 / 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
              0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx /
    data ab / 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
              2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx /
    data bb /-1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
             -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx /
    data abp / 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
               2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx /
    data bbp /-1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
              -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx /

    do j = n1 , n2
      sqti = sqrt(tbar(j,kn))
      rsqti = d_one/sqti
      ! h2o transmission
      tt = abs(tbar(j,kn)-250.0_rkx)
      ds2c = abs(s2c(j,k2+1)-s2c(j,k2))*uinpl(j,kn)
      duptyp = abs(uptype(j,k2+1)-uptype(j,k2))*uinpl(j,kn)
      do l = 1 , 6
        psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
        phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
        p1 = pnew(j)*(psi1/phi1)/sslp
        w1 = dw(j)*winpl(j,kn)*phi1
        tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1))-d_one)-g3(l) * &
                     ds2c-g4(l)*duptyp)
      end do
      du1 = abs(ucfc11(j,k2+1)-ucfc11(j,k2))*winpl(j,kn)
      du2 = abs(ucfc12(j,k2+1)-ucfc12(j,k2))*winpl(j,kn)
      ! cfc transmissions
      tcfc3 = exp(-175.005_rkx*du1)
      tcfc4 = exp(-1202.18_rkx*du1)
      tcfc6 = exp(-5786.73_rkx*du2)
      tcfc7 = exp(-2873.51_rkx*du2)
      tcfc8 = exp(-2085.59_rkx*du2)
      ! Absorptivity for CFC11 bands
      acfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*du1))*tw(1)*bplnk(7,j,kn)
      acfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*du1))*tw(2)*bplnk(8,j,kn)
      acfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*bplnk(9,j,kn)
      acfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*bplnk(10,j,kn)
      ! Absorptivity for CFC12 bands
      acfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*du2))*tw(3)*bplnk(11,j,kn)
      acfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*bplnk(12,j,kn)
      acfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*bplnk(13,j,kn)
      acfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*bplnk(14,j,kn)
      ! Emissivity for CH4 band 1306 cm-1
      tlw = exp(-d_one*sqrt(up2(j)))
      duch4 = abs(uch4(j,k2+1)-uch4(j,k2))*winpl(j,kn)
      dbetac = 2.94449_rkx*pinpl(j,kn)*rsqti/sslp
      ach4 = 6.00444_rkx*sqti*log(d_one+func(duch4,dbetac))*tlw*bplnk(3,j,kn)
      tch4 = d_one/(d_one+0.02_rkx*func(duch4,dbetac))
      ! Absorptivity for N2O bands
      du01 = abs(un2o0(j,k2+1)-un2o0(j,k2))*winpl(j,kn)
      du11 = abs(un2o1(j,k2+1)-un2o1(j,k2))*winpl(j,kn)
      dbeta01 = 19.399_rkx*pinpl(j,kn)*rsqti/sslp
      dbeta11 = dbeta01
      ! 1285 cm-1 band
      an2o1 = 2.35558_rkx*sqti * &
              log(d_one+func(du01,dbeta01)+func(du11,dbeta11)) * &
              tlw*tch4*bplnk(4,j,kn)
      du02 = 0.100090_rkx*du01
      du12 = 0.0992746_rkx*du11
      dbeta02 = 0.964282_rkx*dbeta01
      ! 589 cm-1 band
      an2o2 = 2.65581_rkx*sqti * &
              log(d_one+func(du02,dbeta02)+func(du12,dbeta02)) * &
              tco2(j)*th2o(j)*bplnk(5,j,kn)
      du03 = 0.0333767_rkx*du01
      dbeta03 = 0.982143_rkx*dbeta01
      ! 1168 cm-1 band
      an2o3 = 2.54034_rkx*sqti*log(d_one+func(du03,dbeta03))*tw(6) * &
              tcfc8*bplnk(6,j,kn)
      ! Emissivity for 1064 cm-1 band of CO2
      du11 = abs(uco211(j,k2+1)-uco211(j,k2))*winpl(j,kn)
      du12 = abs(uco212(j,k2+1)-uco212(j,k2))*winpl(j,kn)
      du13 = abs(uco213(j,k2+1)-uco213(j,k2))*winpl(j,kn)
      dbetc1 = 2.97558_rkx*pinpl(j,kn)*rsqti/sslp
      dbetc2 = d_two*dbetc1
      aco21 = 3.7571_rkx*sqti * &
              log(d_one+func(du11,dbetc1)+func(du12,dbetc2) + &
              func(du13,dbetc2))*to3(j)*tw(5)*tcfc4*tcfc7*bplnk(2,j,kn)
      ! Emissivity for 961 cm-1 band of co2
      du21 = abs(uco221(j,k2+1)-uco221(j,k2))*winpl(j,kn)
      du22 = abs(uco222(j,k2+1)-uco222(j,k2))*winpl(j,kn)
      du23 = abs(uco223(j,k2+1)-uco223(j,k2))*winpl(j,kn)
      aco22 = 3.8443_rkx*sqti * &
              log(d_one+func(du21,dbetc1)+func(du22,dbetc1) + &
              func(du23,dbetc2))*tw(4)*tcfc3*tcfc6*bplnk(1,j,kn)
      ! total trace gas absorptivity
      abstrc(j) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +  &
                  aco21 + aco22
    end do
  end subroutine trcabn
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
    real(rkx) , pointer , dimension(:,:,:) :: abplnk1 , abplnk2
    real(rkx) , pointer , dimension(:,:) :: emplnk
    real(rkx) , pointer , dimension(:,:) :: tint , tlayr
    real(rkx) , pointer , dimension(:) :: tplnke
    intent (in) tint , tlayr , tplnke
    intent (inout) abplnk1 , abplnk2 , emplnk
    !
    ! wvl   - wavelength index
    ! f1    - Planck function factor
    ! f2    -       "
    ! f3    -       "
    !
    real(rkx) , dimension(14) :: f1 , f2 , f3
    integer(ik4) :: j , k , wvl

    data f1 / 5.85713e8_rkx , 7.94950e8_rkx , 1.47009e9_rkx , 1.40031e9_rkx , &
              1.34853e8_rkx , 1.05158e9_rkx , 3.35370e8_rkx , 3.99601e8_rkx , &
              5.35994e8_rkx , 8.42955e8_rkx , 4.63682e8_rkx , 5.18944e8_rkx , &
              8.83202e8_rkx , 1.03279e9_rkx/
    data f2 / 2.02493e11_rkx , 3.04286e11_rkx , 6.90698e11_rkx , &
              6.47333e11_rkx , 2.85744e10_rkx , 4.41862e11_rkx , &
              9.62780e10_rkx , 1.21618e11_rkx , 1.79905e11_rkx , &
              3.29029e11_rkx , 1.48294e11_rkx , 1.72315e11_rkx , &
              3.50140e11_rkx , 4.31364e11_rkx/
    data f3 / 1383.0_rkx , 1531.0_rkx , 1879.0_rkx , 1849.0_rkx ,  848.0_rkx , &
              1681.0_rkx , 1148.0_rkx , 1217.0_rkx , 1343.0_rkx , 1561.0_rkx , &
              1279.0_rkx , 1328.0_rkx , 1586.0_rkx , 1671.0_rkx/
    !
    ! Calculate emissivity Planck factor
    !
    do wvl = 1 , 14
      do j = n1 , n2
        emplnk(wvl,j) = f1(wvl)/(tplnke(j)**4 * &
                      (exp(f3(wvl)/tplnke(j))-d_one))
      end do
    end do
    !
    ! Calculate absorptivity Planck factor for tint and tlayr temperatures
    !
    do wvl = 1 , 14
      do k = 1 , kzp1
        do j = n1 , n2
          ! non-nearlest layer function
          abplnk1(wvl,j,k) = (f2(wvl)*exp(f3(wvl)/tint(j,k)))        &
                           & /(tint(j,k)**5*                      &
                           & (exp(f3(wvl)/tint(j,k))-d_one)**2)
          ! nearest layer function
          abplnk2(wvl,j,k) = (f2(wvl)*exp(f3(wvl)/tlayr(j,k)))       &
                           & /(tlayr(j,k)**5*                     &
                           & (exp(f3(wvl)/tlayr(j,k))-d_one)**2)
        end do
      end do
    end do
  end subroutine trcplk
  !
  !----------------------------------------------------------------------
  !  Calculate emissivity for CH4, N2O, CFC11 and CFC12 bands.
  !            Coded by J.T. Kiehl November 21, 1994
  !-----------------------------------------------------------------------
  !
  !------------------------------Arguments--------------------------------
  !
  ! co2t   - pressure weighted temperature
  ! pnm    - interface pressure
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
  ! uptype - continuum path length
  ! bn2o0  - pressure factor for n2o
  ! bn2o1  - pressure factor for n2o
  ! bch4   - pressure factor for ch4
  ! emplnk - emissivity Planck factor
  ! th2o   - water vapor overlap factor
  ! tco2   - co2 overlap factor
  ! to3    - o3 overlap factor
  ! s2c    - h2o continuum path length
  ! w      - h2o path length
  ! up2    - pressure squared h2o path length
  ! k      - level index
  !
  !------------------------------Output Arguments-------------------------
  !
  ! emstrc - total trace gas emissivity
  !
  !-----------------------------------------------------------------------
  !
  subroutine trcems(n1,n2,k,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,  &
                    bn2o0,bn2o1,uch4,bch4,uco211,uco212,uco213,uco221, &
                    uco222,uco223,uptype,w,s2c,up2,emplnk,th2o,tco2,   &
                    to3,emstrc)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2 , k
    real(rkx) , pointer , dimension(:,:) :: bch4 , bn2o0 , bn2o1 , co2t , &
             emstrc , pnm , s2c , ucfc11 , ucfc12 , uch4 , uco211 ,     &
             uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,       &
             un2o1 , uptype , w
    real(rkx) , pointer , dimension(:,:) :: emplnk
    real(rkx) , pointer , dimension(:) :: tco2 , th2o , to3 , up2
    intent (in) bch4 , bn2o0 , bn2o1 , co2t , emplnk , pnm , s2c ,    &
                tco2 , th2o , to3 , ucfc11 , ucfc12 , uch4 , uco211 , &
                uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,  &
                un2o1 , up2 , uptype , w
    intent (inout) emstrc
    !
    ! sqti   - square root of mean temp
    ! ecfc1  - emissivity of cfc11 798 cm-1 band
    ! ecfc2  -     "      "    "   846 cm-1 band
    ! ecfc3  -     "      "    "   933 cm-1 band
    ! ecfc4  -     "      "    "   1085 cm-1 band
    ! ecfc5  -     "      "  cfc12 889 cm-1 band
    ! ecfc6  -     "      "    "   923 cm-1 band
    ! ecfc7  -     "      "    "   1102 cm-1 band
    ! ecfc8  -     "      "    "   1161 cm-1 band
    ! u01    - n2o path length
    ! u11    - n2o path length
    ! beta01 - n2o pressure factor
    ! beta11 - n2o pressure factor
    ! en2o1  - emissivity of the 1285 cm-1 N2O band
    ! u02    - n2o path length
    ! u12    - n2o path length
    ! beta02 - n2o pressure factor
    ! en2o2  - emissivity of the 589 cm-1 N2O band
    ! u03    - n2o path length
    ! beta03 - n2o pressure factor
    ! en2o3  - emissivity of the 1168 cm-1 N2O band
    ! betac  - ch4 pressure factor
    ! ech4   - emissivity of 1306 cm-1 CH4 band
    ! betac1 - co2 pressure factor
    ! betac2 - co2 pressure factor
    ! eco21  - emissivity of 1064 cm-1 CO2 band
    ! eco22  - emissivity of 961 cm-1 CO2 band
    ! tt     - temp. factor for h2o overlap factor
    ! psi1   - narrow band h2o temp. factor
    ! phi1   -            "
    ! p1     - h2o line overlap factor
    ! w1     -          "
    ! tw     - h2o transmission overlap
    ! g1     - h2o overlap factor
    ! g2     -          "
    ! g3     -          "
    ! g4     -          "
    ! ab     -          "
    ! bb     -          "
    ! abp    -          "
    ! bbp    -          "
    ! tcfc3  - transmission for cfc11 band
    ! tcfc4  -         "
    ! tcfc6  - transmission for cfc12 band
    ! tcfc7  -          "
    ! tcfc8  -          "
    ! tlw    - h2o overlap factor
    ! tch4   - ch4 overlap factor
    !
    real(rkx) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4 , tw
    real(rkx) :: beta01 , beta02 , beta03 , beta11 , betac , sqti , tt , &
               betac1 , betac2 , ecfc1 , ecfc2 , ecfc3 , ecfc4 ,       &
               ecfc5 , ecfc6 , ecfc7 , ecfc8 , ech4 , eco21 , eco22 ,  &
               en2o1 , en2o2 , en2o3 , p1 , phi1 , psi1 , tcfc3 ,      &
               tcfc4 , tcfc6 , tcfc7 , tcfc8 , tch4 , tlw , u01 ,      &
               u02 , u03 , u11 , u12 , w1
    integer(ik4) :: j , l

    data g1 / 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
              0.0304380_rkx , 0.0540398_rkx , 0.0321962_rkx /
    data g2 / 14.48320_rkx ,  4.302420_rkx ,  5.23523_rkx , &
               3.25342_rkx ,  0.698935_rkx , 16.55990_rkx /
    data g3 / 26.1898_rkx , 18.44760_rkx , 15.36330_rkx , &
              12.1927_rkx ,  9.14992_rkx ,  8.07092_rkx /
    data g4 / 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
              0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx /
    data ab / 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
              2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx /
    data bb / -1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
              -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx /
    data abp / 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
               2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx /
    data bbp / -1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
               -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx /

    do j = n1 , n2
      sqti = sqrt(co2t(j,k))
      ! Transmission for h2o
      tt = abs(co2t(j,k)-250.0_rkx)
      ! transmission due to cfc bands
      do l = 1 , 6
        psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
        phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
        p1 = pnm(j,k)*(psi1/phi1)/sslp
        w1 = w(j,k)*phi1
        tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1)) - &
                      d_one)-g3(l)*s2c(j,k)-g4(l)*uptype(j,k))
      end do
      tcfc3 = exp(-175.005_rkx*ucfc11(j,k))
      tcfc4 = exp(-1202.18_rkx*ucfc11(j,k))
      tcfc6 = exp(-5786.73_rkx*ucfc12(j,k))
      tcfc7 = exp(-2873.51_rkx*ucfc12(j,k))
      tcfc8 = exp(-2085.59_rkx*ucfc12(j,k))
      ! Emissivity for CFC11 bands
      ecfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*ucfc11(j,k)))*tw(1)*emplnk(7,j)
      ecfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*ucfc11(j,k)))*tw(2)*emplnk(8,j)
      ecfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*emplnk(9,j)
      ecfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*emplnk(10,j)
      ! Emissivity for CFC12 bands
      ecfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*ucfc12(j,k)))*tw(3)*emplnk(11,j)
      ecfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*emplnk(12,j)
      ecfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*emplnk(13,j)
      ecfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*emplnk(14,j)
      ! Emissivity for CH4 band 1306 cm-1
      tlw = exp(-d_one*sqrt(up2(j)))
      betac = bch4(j,k)/uch4(j,k)
      ech4 = 6.00444_rkx*sqti*log(d_one+func(uch4(j,k),betac))*tlw*emplnk(3,j)
      tch4 = d_one/(d_one+0.02_rkx*func(uch4(j,k),betac))
      ! Emissivity for N2O bands
      u01 = un2o0(j,k)
      u11 = un2o1(j,k)
      beta01 = bn2o0(j,k)/un2o0(j,k)
      beta11 = bn2o1(j,k)/un2o1(j,k)
      ! 1285 cm-1 band
      en2o1 = 2.35558_rkx*sqti * &
             log(d_one+func(u01,beta01)+func(u11,beta11))*tlw*tch4*emplnk(4,j)
      u02 = 0.100090_rkx*u01
      u12 = 0.0992746_rkx*u11
      beta02 = 0.964282_rkx*beta01
      ! 589 cm-1 band
      en2o2 = 2.65581_rkx*sqti * &
              log(d_one+func(u02,beta02)+func(u12,beta02))*tco2(j) * &
              th2o(j)*emplnk(5,j)
      u03 = 0.0333767_rkx*u01
      beta03 = 0.982143_rkx*beta01
      ! 1168 cm-1 band
      en2o3 = 2.54034_rkx*sqti*log(d_one+func(u03,beta03))*tw(6) * &
             tcfc8*emplnk(6,j)
      ! Emissivity for 1064 cm-1 band of CO2
      betac1 = 2.97558_rkx*pnm(j,k)/(sslp*sqti)
      betac2 = d_two*betac1
      eco21 = 3.7571_rkx*sqti * &
              log(d_one+func(uco211(j,k),betac1) +  &
                        func(uco212(j,k),betac2) +  &
                        func(uco213(j,k),betac2)) * &
              to3(j)*tw(5)*tcfc4*tcfc7*emplnk(2,j)
      ! Emissivity for 961 cm-1 band
      eco22 = 3.8443_rkx*sqti * &
              log(d_one+func(uco221(j,k),betac1) +  &
                        func(uco222(j,k),betac1) +  &
                        func(uco223(j,k),betac2)) * &
              tw(4)*tcfc3*tcfc6*emplnk(1,j)
      ! total trace gas emissivity
      emstrc(j,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 + ecfc6 +  &
                    ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 +   &
                    eco21 + eco22
    end do
  end subroutine trcems

  pure real(rkx) function func(u,b)
    implicit none
    real(rkx) , intent(in) :: u , b
    func = u/sqrt(d_four+u*(d_one+d_one/b))
  end function func

end module mod_rad_tracer
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
