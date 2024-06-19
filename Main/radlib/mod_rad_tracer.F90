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
  use mod_memutil
  use mod_dynparam

  implicit none

  private

  public :: trcmix , trcab , trcabn , trcems

  contains

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
  pure subroutine trcmix(dlat,xptrop,pmid,n2o0,ch40,cfc110,cfc120, &
                         n2o,ch4,cfc11,cfc12)
!$acc routine seq
    implicit none
    !
    ! dlat   - latitude in degrees
    ! xn2o   - pressure scale height for n2o
    ! xch4   - pressure scale height for ch4
    ! xcfc11 - pressure scale height for cfc11
    ! xcfc12 - pressure scale height for cfc12
    ! ptrop  - pressure level of tropopause
    ! pratio - pressure divided by ptrop
    !
    real(rkx) , intent(in) :: pmid , dlat , xptrop
    real(rkx) , intent(in) :: n2o0 , ch40 , cfc110 , cfc120
    real(rkx) , intent(out) :: cfc11 , cfc12 , ch4 , n2o
#ifndef RCEMIP
    real(rkx) :: alat
#endif
    real(rkx) :: pratio , xcfc11 , xcfc12 , xch4 , xn2o

#ifdef RCEMIP
    xn2o = 0.3478_rkx
    xch4 = 0.2353_rkx
    xcfc11 = 0.7273_rkx
    xcfc12 = 0.4000_rkx
#else
    alat = abs(dlat) ! This is absolute value of latitude in degrees
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
    !  set stratospheric scale height factor for gases
    if ( pmid >= xptrop ) then
      ch4 = ch40
      n2o = n2o0
      cfc11 = cfc110
      cfc12 = cfc120
    else
      pratio = pmid/xptrop
      ch4 = ch40*(pratio**xch4)
      n2o = n2o0*(pratio**xn2o)
      cfc11 = cfc110*(pratio**xcfc11)
      cfc12 = cfc120*(pratio**xcfc12)
    end if
  end subroutine trcmix
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
  pure real(rkx) function trcab(dpint,ds2c,duptyp,du1,du2,duch4,dbetac, &
      du01,du11,dbeta01,dbeta11,duco11,duco12,duco13,duco21,duco22,     &
      duco23,dw,pnew,to3co2,dplh2o,tco2,th2o,to3,abplnk1) result(abstrc)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: dpint , ds2c , duptyp , du1 , du2
    real(rkx) , intent(in) :: duch4 , dbetac , du01 , du11
    real(rkx) , intent(in) :: dbeta01 , dbeta11
    real(rkx) , intent(in) :: duco11 , duco12 , duco13
    real(rkx) , intent(in) :: duco21 , duco22 , duco23
    real(rkx) , intent(in) :: to3co2 , dw , pnew , dplh2o
    real(rkx) , intent(in) :: tco2 , th2o , to3
    real(rkx) , dimension(14) , intent(in) :: abplnk1
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
    real(rkx) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 , &
      acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , dbeta02 ,   &
      dbeta03 , dbetc1 , dbetc2 , du02 , du12 , du03 , p1 , phi1 , psi1 ,&
      tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 , tch4 , tlw , w1 , sqti , tt
    real(rkx) , dimension(6) :: tw
    integer(ik4) :: l

    real(rkx) , dimension(6) , parameter :: g1 = &
        [ 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
          0.0304380_rkx , 0.0540398_rkx , 0.0321962_rkx ]
    real(rkx) , dimension(6) , parameter :: g2 = &
        [ 14.48320_rkx , 4.302420_rkx ,  5.23523_rkx , &
           3.25342_rkx , 0.698935_rkx , 16.55990_rkx ]
    real(rkx) , dimension(6) , parameter :: g3 = &
        [ 26.18980_rkx , 18.44760_rkx , 15.36330_rkx , &
          12.19270_rkx ,  9.14992_rkx ,  8.07092_rkx ]
    real(rkx) , dimension(6) , parameter :: g4 = &
        [ 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
          0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx ]
    real(rkx) , dimension(6) , parameter :: ab = &
        [ 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
          2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bb = &
        [ -1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
          -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx ]
    real(rkx) , dimension(6) , parameter :: abp = &
        [ 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
          2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bbp = &
        [ -1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
          -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx ]

    sqti = sqrt(to3co2)
    ! h2o transmission
    tt = abs(to3co2-250.0_rkx)
    do l = 1 , 6
      psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
      phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
      p1 = pnew*(psi1/phi1)/sslp
      w1 = dw*phi1
      tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1)) - &
                   d_one)-g3(l)*ds2c-g4(l)*duptyp)
    end do
    ! cfc transmissions
    tcfc3 = exp(-175.005_rkx*du1)
    tcfc4 = exp(-1202.18_rkx*du1)
    tcfc6 = exp(-5786.73_rkx*du2)
    tcfc7 = exp(-2873.51_rkx*du2)
    tcfc8 = exp(-2085.59_rkx*du2)
    ! Absorptivity for CFC11 bands
    acfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*du1))*tw(1)*abplnk1(7)
    acfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*du1))*tw(2)*abplnk1(8)
    acfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*abplnk1(9)
    acfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*abplnk1(10)
    ! Absorptivity for CFC12 bands
    acfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*du2))*tw(3)*abplnk1(11)
    acfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*abplnk1(12)
    acfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*abplnk1(13)
    acfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*abplnk1(14)
    ! Emissivity for CH4 band 1306 cm-1
    tlw = exp(-d_one*sqrt(dplh2o))
    ach4 = 6.00444_rkx*sqti*log(d_one+func(duch4,dbetac))*tlw*abplnk1(3)
    tch4 = d_one/(d_one+0.02_rkx*func(duch4,dbetac))
    ! Absorptivity for N2O bands
    ! 1285 cm-1 band
    an2o1 = 2.35558_rkx*sqti * log(d_one+func(du01,dbeta01) + &
            func(du11,dbeta11)) * tlw*tch4*abplnk1(4)
    du02 = 0.100090_rkx*du01
    du12 = 0.0992746_rkx*du11
    dbeta02 = 0.964282_rkx*dbeta01
    ! 589 cm-1 band
    an2o2 = 2.65581_rkx*sqti * log(d_one+func(du02,dbeta02) + &
            func(du12,dbeta02))*th2o*tco2*abplnk1(5)
    du03 = 0.0333767_rkx*du01
    dbeta03 = 0.982143_rkx*dbeta01
    ! 1168 cm-1 band
    an2o3 = 2.54034_rkx*sqti*log(d_one+func(du03,dbeta03)) * &
            tw(6)*tcfc8*abplnk1(6)
    ! Emissivity for 1064 cm-1 band of CO2
    dbetc1 = 2.97558_rkx*dpint/(d_two*sslp*sqti)
    dbetc2 = d_two*dbetc1
    aco21 = 3.7571_rkx*sqti * &
            log(d_one+func(duco11,dbetc1)+func(duco12,dbetc2) + &
            func(duco13,dbetc2))*to3*tw(5)*tcfc4*tcfc7*abplnk1(2)
    ! Emissivity for 961 cm-1 band
    aco22 = 3.8443_rkx*sqti * &
            log(d_one+func(duco21,dbetc1)+func(duco22,dbetc1) + &
            func(duco23,dbetc2))*tw(4)*tcfc3*tcfc6*abplnk1(1)
    ! total trace gas absorptivity
    abstrc = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 +  &
             acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +   &
             aco21 + aco22
  end function trcab
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
  pure real(rkx) function trcabn(tbar,dw,pnew,tco2,th2o,to3,up2,     &
      pinpl,winpl,ds2c,duptyp,du1,du2,duch4,du01,du11,duco11,duco12, &
      duco13,duco21,duco22,duco23,bplnk) result(abstrc)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: tbar , dw , pnew , tco2 , th2o , to3 , up2
    real(rkx) , intent(in) :: winpl , pinpl , ds2c , duptyp , du1 , du2
    real(rkx) , intent(in) :: duch4 , du01 , du11 , duco11 , duco12
    real(rkx) , intent(in) :: duco13 , duco21 , duco22 , duco23
    real(rkx) , dimension(14) , intent(in) :: bplnk
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
    real(rkx) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 ,&
      acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , dbeta01 ,  &
      dbeta02 , dbeta03 , dbeta11 , dbetac , dbetc1 , dbetc2 , du02 ,   &
      du03 , p1 , phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 , &
      tch4 , tlw , w1 , rsqti , sqti , tt , du12
    real(rkx) , dimension(6) :: tw
    integer(ik4) :: l

    real(rkx) , dimension(6) , parameter :: g1 = &
     [ 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
       0.0304380_rkx , 0.0540398_rkx ,  0.0321962_rkx ]
    real(rkx) , dimension(6) , parameter :: g2 = &
     [ 14.483200_rkx ,  4.302420_rkx ,  5.2352300_rkx , &
        3.253420_rkx ,  0.698935_rkx , 16.5599000_rkx ]
    real(rkx) , dimension(6) , parameter :: g3 = &
     [ 26.18980_rkx , 18.44760_rkx , 15.36330_rkx , &
       12.19270_rkx ,  9.14992_rkx ,  8.07092_rkx ]
    real(rkx) , dimension(6) , parameter :: g4 = &
     [ 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
       0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx ]
    real(rkx) , dimension(6) , parameter :: ab = &
     [ 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
       2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bb = &
     [ -1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
       -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx ]
    real(rkx) , dimension(6) , parameter :: abp = &
     [ 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
       2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bbp = &
     [ -1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
       -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx ]

    sqti = sqrt(tbar)
    rsqti = d_one/sqti
    ! h2o transmission
    tt = abs(tbar-250.0_rkx)
    do l = 1 , 6
      psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
      phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
      p1 = pnew*(psi1/phi1)/sslp
      w1 = dw*winpl*phi1
      tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1))-d_one)-g3(l) * &
                   ds2c-g4(l)*duptyp)
    end do
    ! cfc transmissions
    tcfc3 = exp(-175.005_rkx*du1)
    tcfc4 = exp(-1202.18_rkx*du1)
    tcfc6 = exp(-5786.73_rkx*du2)
    tcfc7 = exp(-2873.51_rkx*du2)
    tcfc8 = exp(-2085.59_rkx*du2)
    ! Absorptivity for CFC11 bands
    acfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*du1))*tw(1)*bplnk(7)
    acfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*du1))*tw(2)*bplnk(8)
    acfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*bplnk(9)
    acfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*bplnk(10)
    ! Absorptivity for CFC12 bands
    acfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*du2))*tw(3)*bplnk(11)
    acfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*bplnk(12)
    acfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*bplnk(13)
    acfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*bplnk(14)
    ! Emissivity for CH4 band 1306 cm-1
    tlw = exp(-d_one*sqrt(up2))
    dbetac = 2.94449_rkx*pinpl*rsqti/sslp
    ach4 = 6.00444_rkx*sqti*log(d_one+func(duch4,dbetac))*tlw*bplnk(3)
    tch4 = d_one/(d_one+0.02_rkx*func(duch4,dbetac))
    ! Absorptivity for N2O bands
    dbeta01 = 19.399_rkx*pinpl*rsqti/sslp
    dbeta11 = dbeta01
    ! 1285 cm-1 band
    an2o1 = 2.35558_rkx*sqti * &
            log(d_one+func(du01,dbeta01)+func(du11,dbeta11)) * &
            tlw*tch4*bplnk(4)
    du02 = 0.100090_rkx*du01
    du12 = 0.0992746_rkx*du11
    dbeta02 = 0.964282_rkx*dbeta01
    ! 589 cm-1 band
    an2o2 = 2.65581_rkx*sqti * &
            log(d_one+func(du02,dbeta02)+func(du12,dbeta02)) * &
            tco2*th2o*bplnk(5)
    du03 = 0.0333767_rkx*du01
    dbeta03 = 0.982143_rkx*dbeta01
    ! 1168 cm-1 band
    an2o3 = 2.54034_rkx*sqti*log(d_one+func(du03,dbeta03))*tw(6) * &
            tcfc8*bplnk(6)
    ! Emissivity for 1064 cm-1 band of CO2
    dbetc1 = 2.97558_rkx*pinpl*rsqti/sslp
    dbetc2 = d_two*dbetc1
    aco21 = 3.7571_rkx*sqti * &
            log(d_one+func(duco11,dbetc1)+func(duco12,dbetc2) + &
            func(duco13,dbetc2))*to3*tw(5)*tcfc4*tcfc7*bplnk(2)
    ! Emissivity for 961 cm-1 band of co2
    aco22 = 3.8443_rkx*sqti * &
            log(d_one+func(duco21,dbetc1)+func(duco22,dbetc1) + &
            func(duco23,dbetc2))*tw(4)*tcfc3*tcfc6*bplnk(1)
    ! total trace gas absorptivity
    abstrc = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
             acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +  &
             aco21 + aco22
  end function trcabn
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
  pure real(rkx) function trcems(co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,  &
     bn2o0,bn2o1,uch4,bch4,uco211,uco212,uco213,uco221,uco222,uco223, &
     uptype,w,s2c,up2,emplnk,th2o,tco2,to3) result(emstrc)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: bn2o0 , bn2o1
    real(rkx) , intent(in) :: un2o0 , un2o1
    real(rkx) , intent(in) :: bch4 , uch4 , co2t
    real(rkx) , intent(in) :: pnm , s2c
    real(rkx) , intent(in) :: ucfc11 , ucfc12
    real(rkx) , intent(in) :: uco211 , uco212
    real(rkx) , intent(in) :: uco213 , uco221
    real(rkx) , intent(in) :: uco222 , uco223
    real(rkx) , intent(in) :: tco2 , th2o , to3 , up2
    real(rkx) , intent(in) :: uptype , w
    real(rkx) , dimension(14) , intent(in) :: emplnk
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
    real(rkx) :: beta01 , beta02 , beta03 , beta11 , betac , sqti , tt , &
                 betac1 , betac2 , ecfc1 , ecfc2 , ecfc3 , ecfc4 ,       &
                 ecfc5 , ecfc6 , ecfc7 , ecfc8 , ech4 , eco21 , eco22 ,  &
                 en2o1 , en2o2 , en2o3 , p1 , phi1 , psi1 , tcfc3 ,      &
                 tcfc4 , tcfc6 , tcfc7 , tcfc8 , tch4 , tlw , u01 ,      &
                 u02 , u03 , u11 , u12 , w1
    real(rkx) , dimension(6) :: tw
    integer(ik4) :: l

    real(rkx) , dimension(6) , parameter :: g1 = &
      [ 0.0468556_rkx , 0.0397454_rkx , 0.0407664_rkx , &
        0.0304380_rkx , 0.0540398_rkx , 0.0321962_rkx ]
    real(rkx) , dimension(6) , parameter :: g2 = &
      [ 14.48320_rkx ,  4.302420_rkx ,  5.23523_rkx , &
         3.25342_rkx ,  0.698935_rkx , 16.55990_rkx  ]
    real(rkx) , dimension(6) , parameter :: g3 = &
      [ 26.1898_rkx , 18.44760_rkx , 15.36330_rkx , &
        12.1927_rkx ,  9.14992_rkx ,  8.07092_rkx ]
    real(rkx) , dimension(6) , parameter :: g4 = &
      [ 0.0261782_rkx , 0.0369516_rkx , 0.0307266_rkx , &
        0.0243854_rkx , 0.0182932_rkx , 0.0161418_rkx ]
    real(rkx) , dimension(6) , parameter :: ab = &
      [ 3.0857e-2_rkx , 2.3524e-2_rkx , 1.7310e-2_rkx , &
        2.6661e-2_rkx , 2.8074e-2_rkx , 2.2915e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bb = &
      [ -1.3512e-4_rkx ,-6.8320e-5_rkx ,-3.2609e-5_rkx , &
        -1.0228e-5_rkx ,-9.5743e-5_rkx ,-1.0304e-4_rkx ]
    real(rkx) , dimension(6) , parameter :: abp = &
      [ 2.9129e-2_rkx , 2.4101e-2_rkx , 1.9821e-2_rkx , &
        2.6904e-2_rkx , 2.9458e-2_rkx , 1.9892e-2_rkx ]
    real(rkx) , dimension(6) , parameter :: bbp = &
      [ -1.3139e-4_rkx ,-5.5688e-5_rkx ,-4.6380e-5_rkx , &
        -8.0362e-5_rkx ,-1.0115e-4_rkx ,-8.8061e-5_rkx ]

    sqti = sqrt(co2t)
    ! Transmission for h2o
    tt = abs(co2t-250.0_rkx)
    ! transmission due to cfc bands
    do l = 1 , 6
      psi1 = exp(abp(l)*tt+bbp(l)*tt*tt)
      phi1 = exp(ab(l)*tt+bb(l)*tt*tt)
      p1 = pnm*(psi1/phi1)/sslp
      w1 = w*phi1
      tw(l) = exp(-g1(l)*p1*(sqrt(d_one+g2(l)*(w1/p1)) - &
              d_one)-g3(l)*s2c-g4(l)*uptype)
    end do
    tcfc3 = exp(-175.005_rkx*ucfc11)
    tcfc4 = exp(-1202.18_rkx*ucfc11)
    tcfc6 = exp(-5786.73_rkx*ucfc12)
    tcfc7 = exp(-2873.51_rkx*ucfc12)
    tcfc8 = exp(-2085.59_rkx*ucfc12)
    ! Emissivity for CFC11 bands
    ecfc1 = 50.0_rkx*(d_one-exp(-54.09_rkx*ucfc11))*tw(1)*emplnk(7)
    ecfc2 = 60.0_rkx*(d_one-exp(-5130.03_rkx*ucfc11))*tw(2)*emplnk(8)
    ecfc3 = 60.0_rkx*(d_one-tcfc3)*tw(4)*tcfc6*emplnk(9)
    ecfc4 = 100.0_rkx*(d_one-tcfc4)*tw(5)*emplnk(10)
    ! Emissivity for CFC12 bands
    ecfc5 = 45.0_rkx*(d_one-exp(-1272.35_rkx*ucfc12))*tw(3)*emplnk(11)
    ecfc6 = 50.0_rkx*(d_one-tcfc6)*tw(4)*emplnk(12)
    ecfc7 = 80.0_rkx*(d_one-tcfc7)*tw(5)*tcfc4*emplnk(13)
    ecfc8 = 70.0_rkx*(d_one-tcfc8)*tw(6)*emplnk(14)
    ! Emissivity for CH4 band 1306 cm-1
    tlw = exp(-d_one*sqrt(up2))
    betac = bch4/uch4
    ech4 = 6.00444_rkx*sqti*log(d_one+func(uch4,betac))*tlw*emplnk(3)
    tch4 = d_one/(d_one+0.02_rkx*func(uch4,betac))
    ! Emissivity for N2O bands
    u01 = un2o0
    u11 = un2o1
    beta01 = bn2o0/un2o0
    beta11 = bn2o1/un2o1
    ! 1285 cm-1 band
    en2o1 = 2.35558_rkx*sqti * &
           log(d_one+func(u01,beta01)+func(u11,beta11))*tlw*tch4*emplnk(4)
    u02 = 0.100090_rkx*u01
    u12 = 0.0992746_rkx*u11
    beta02 = 0.964282_rkx*beta01
    ! 589 cm-1 band
    en2o2 = 2.65581_rkx*sqti * &
            log(d_one+func(u02,beta02)+func(u12,beta02))*tco2 * &
            th2o*emplnk(5)
    u03 = 0.0333767_rkx*u01
    beta03 = 0.982143_rkx*beta01
    ! 1168 cm-1 band
    en2o3 = 2.54034_rkx*sqti*log(d_one+func(u03,beta03))*tw(6)*tcfc8*emplnk(6)
    ! Emissivity for 1064 cm-1 band of CO2
    betac1 = 2.97558_rkx*pnm/(sslp*sqti)
    betac2 = d_two*betac1
    eco21 = 3.7571_rkx*sqti * &
            log(d_one+func(uco211,betac1) + func(uco212,betac2) + &
                func(uco213,betac2))*to3*tw(5)*tcfc4*tcfc7*emplnk(2)
    ! Emissivity for 961 cm-1 band
    eco22 = 3.8443_rkx*sqti * &
            log(d_one+func(uco221,betac1) + func(uco222,betac1) +  &
                func(uco223,betac2))*tw(4)*tcfc3*tcfc6*emplnk(1)
    ! total trace gas emissivity
    emstrc = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 + ecfc6 +  &
             ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 +   &
             eco21 + eco22
  end function trcems

  pure real(rkx) function func(u,b)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: u , b
    func = u/sqrt(d_four+u*(d_one+d_one/b))
  end function func

end module mod_rad_tracer
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
