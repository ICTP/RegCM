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

      subroutine trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213,
     $                         uco221, uco222, uco223, bn2o0,
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         dplh2o, abplnk1,tco2,   th2o,
     $                         to3,    abstrc)
c----------------------------------------------------------------------
c Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
c CFC12.
c
c             Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcab.F,v 1.2 1995/02/17 21:28:52 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
C
C------------------------------Commons----------------------------------
#     include "include/crdcon.cb"
C
C------------------------------Arguments--------------------------------
      integer k1,k2
      real(kind=8)  to3co2(plond),       ! pressure weighted temperature
     $     pnm(plond,plevp),    ! interface pressures
     $     ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp)    ! pressure factor for ch4
      real(kind=8)  dw(plond),           ! h2o path length
     $     pnew(plond),         ! pressure
     $     s2c(plond,plevp),    ! continuum path length
     $     uptype(plond,plevp), ! p-type h2o path length
     $     dplh2o(plond)        ! p squared h2o path length
      real(kind=8)  abplnk1(14,plond,plevp) ! Planck factor
      real(kind=8)  tco2(plond),         ! co2 transmission factor
     $     th2o(plond),         ! h2o transmission factor
     $     to3(plond)           ! o3 transmission factor
c
c  Output Arguments
c
      real(kind=8)  abstrc(plond)  ! total trace gas absorptivity
c
c  Local Variables
c
      real(kind=8)  sqti(plond),         ! square root of mean temp
     $     du1,                 ! cfc11 path length
     $     du2,                 ! cfc12 path length
     $     acfc1,               ! cfc11 absorptivity 798 cm-1
     $     acfc2,               ! cfc11 absorptivity 846 cm-1
     $     acfc3,               ! cfc11 absorptivity 933 cm-1
     $     acfc4,               ! cfc11 absorptivity 1085 cm-1
     $     acfc5,               ! cfc12 absorptivity 889 cm-1
     $     acfc6,               ! cfc12 absorptivity 923 cm-1
     $     acfc7,               ! cfc12 absorptivity 1102 cm-1
     $     acfc8,               ! cfc12 absorptivity 1161 cm-1
     $     du01,                ! n2o path length
     $     dbeta01,             ! n2o pressure factor
     $     dbeta11,             !         "
     $     an2o1,               ! absorptivity of 1285 cm-1 n2o band
     $     du02,                ! n2o path length
     $     dbeta02,             ! n2o pressure factor
     $     an2o2,               ! absorptivity of 589 cm-1 n2o band
     $     du03                 ! n2o path length
      real(kind=8)  dbeta03,             ! n2o pressure factor
     $     an2o3,               ! absorptivity of 1168 cm-1 n2o band
     $     duch4,               ! ch4 path length
     $     dbetac,              ! ch4 pressure factor
     $     ach4,                ! absorptivity of 1306 cm-1 ch4 band
     $     du11,                ! co2 path length
     $     du12,                !       "
     $     du13,                !       "
     $     dbetc1,              ! co2 pressure factor
     $     dbetc2,              ! co2 pressure factor
     $     aco21,               ! absorptivity of 1064 cm-1 band
     $     du21,                ! co2 path length
     $     du22,                !       "
     $     du23,                !       "
     $     aco22                ! absorptivity of 961 cm-1 band
      real(kind=8)  tt(plond),           ! temp. factor for h2o overlap factor
     $     psi1,                !                 "
     $     phi1,                !                 "
     $     p1,                  ! h2o overlap factor
     $     w1,                  !        "
     $     ds2c(plond),         ! continuum path length
     $     duptyp(plond),       ! p-type path length
     $     tw(plond,6),         ! h2o transmission factor
     $     g1(6),               !         "
     $     g2(6),               !         "
     $     g3(6),               !         "
     $     g4(6),               !         "
     $     ab(6),               ! h2o temp. factor
     $     bb(6),               !         "
     $     abp(6),              !         "
     $     bbp(6)               !         "
      real(kind=8)  tcfc3,               ! transmission for cfc11 band
     $     tcfc4,               ! transmission for cfc11 band
     $     tcfc6,               ! transmission for cfc12 band
     $     tcfc7,               ! transmission for cfc12 band
     $     tcfc8,               ! transmission for cfc12 band
     $     tlw,                 ! h2o transmission
     $     tch4                 ! ch4 transmission
c
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
c------------------------------------------------------------------
      real(kind=8)  func, u, b
      func(u,b) = u/dsqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = dsqrt(to3co2(i))
c h2o transmission
         tt(i) = dabs(to3co2(i) - 250.0)
         ds2c(i) = dabs(s2c(i,k1) - s2c(i,k2))
         duptyp(i) = dabs(uptype(i,k1) - uptype(i,k2))
      end do
c
      do l = 1,6
            do i = 1,plon
               psi1 = dexp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
               phi1 = dexp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
               p1 = pnew(i) * (psi1/phi1) / sslp
               w1 = dw(i) * phi1
               tw(i,l) = dexp(- g1(l)*p1*(dsqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                     - g3(l)*ds2c(i)-g4(l)*duptyp(i))
            end do
      end do
c
      do i = 1,plon
            du1 = dabs(ucfc11(i,k1) - ucfc11(i,k2))
            du2 = dabs(ucfc12(i,k1) - ucfc12(i,k2))
c cfc transmissions
            tcfc3 = dexp(-175.005*du1)
            tcfc4 = dexp(-1202.18*du1)
            tcfc6 = dexp(-5786.73*du2)
            tcfc7 = dexp(-2873.51*du2)
            tcfc8 = dexp(-2085.59*du2)
c  Absorptivity for CFC11 bands
            acfc1 = 50.0*(1.0-dexp(-54.09*du1))*tw(i,1)*abplnk1(7,i,k2)
            acfc2 = 60.0*(1.0-dexp(-5130.03*du1))*tw(i,2)
     $                                            *abplnk1(8,i,k2)
            acfc3 = 60.0*(1.0-tcfc3) * tw(i,4)*tcfc6*abplnk1(9,i,k2)
            acfc4 = 100.0*(1.0-tcfc4)* tw(i,5) * abplnk1(10,i,k2)
c  Absorptivity for CFC12 bands
            acfc5 = 45.0*(1.0 - dexp(-1272.35*du2))*tw(i,3)*
     $                                             abplnk1(11,i,k2)
            acfc6 = 50.0*(1.0 - tcfc6)* tw(i,4)*abplnk1(12,i,k2)
            acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5)*tcfc4*abplnk1(13,i,k2)
            acfc8 = 70.0*(1.0 - tcfc8)* tw(i,6)*abplnk1(14,i,k2)
c  Emissivity for CH4 band 1306 cm-1
            tlw = dexp(-1.0*dsqrt(dplh2o(i)))
            duch4 = dabs(uch4(i,k1) - uch4(i,k2))
            dbetac = dabs(bch4(i,k1) - bch4(i,k2))/duch4
            ach4 = 6.00444*sqti(i)*dlog(1.0 + func(duch4,dbetac)) *
     $           tlw * abplnk1(3,i,k2)
            tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
c  Absorptivity for N2O bands
            du01 = dabs(un2o0(i,k1) - un2o0(i,k2))
            du11 = dabs(un2o1(i,k1) - un2o1(i,k2))
            dbeta01 = dabs(bn2o0(i,k1) - bn2o0(i,k2))/du01
            dbeta11 = dabs(bn2o1(i,k1) - bn2o1(i,k2))/du11
c     1285 cm-1 band
            an2o1 = 2.35558*sqti(i)*dlog(1.0 + func(du01,dbeta01)
     $            +  func(du11,dbeta11)) * tlw* tch4*abplnk1(4,i,k2)
            du02 = 0.100090*du01
            du12 = 0.0992746*du11
            dbeta02 = 0.964282*dbeta01
c     589 cm-1 band
            an2o2 = 2.65581*sqti(i)*dlog(1.0 + func(du02,dbeta02)
     $            +  func(du12,dbeta02)) * th2o(i) * tco2(i) *
     $               abplnk1(5,i,k2)
            du03 = 0.0333767*du01
            dbeta03 = 0.982143*dbeta01
c     1168 cm-1 band
            an2o3 = 2.54034*sqti(i)*dlog(1.0 + func(du03,dbeta03)) *
     $           tw(i,6) * tcfc8 * abplnk1(6,i,k2)
c  Emissivity for 1064 cm-1 band of CO2
            du11 = dabs(uco211(i,k1) - uco211(i,k2))
            du12 = dabs(uco212(i,k1) - uco212(i,k2))
            du13 = dabs(uco213(i,k1) - uco213(i,k2))
            dbetc1 = 2.97558*dabs(pnm(i,k1) + pnm(i,k2))/
     $                                           (2.0*sslp*sqti(i))
            dbetc2 = 2.0 * dbetc1
            aco21 = 3.7571*sqti(i)*dlog(1.0 + func(du11,dbetc1)
     $         + func(du12,dbetc2) + func(du13,dbetc2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * abplnk1(2,i,k2)
c  Emissivity for 961 cm-1 band
            du21 = dabs(uco221(i,k1) - uco221(i,k2))
            du22 = dabs(uco222(i,k1) - uco222(i,k2))
            du23 = dabs(uco223(i,k1) - uco223(i,k2))
            aco22 = 3.8443*sqti(i)*dlog(1.0 + func(du21,dbetc1)
     $         + func(du22,dbetc1) + func(du23,dbetc2))
     $         * tw(i,4) * tcfc3 * tcfc6 * abplnk1(1,i,k2)
c total trace gas absorptivity
            abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6
     $            +  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4
     $            +  aco21 + aco22
      end do
      return
      end
