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
 
      subroutine trcabn(k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,    &
                      & uco212,uco213,uco221,uco222,uco223,tbar,bplnk,  &
                      & winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,up2,    &
                      & pnew,abstrc,uinpl)
!
      use mod_dynparam
      use mod_constants , only : sslp
      implicit none
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
! Dummy arguments
!
      integer :: k2 , kn
      real(8) , dimension(iym1) :: abstrc , dw , pnew , tco2 , th2o ,  &
                                  & to3 , up2
      real(8) , dimension(14,iym1,4) :: bplnk
      real(8) , dimension(iym1,4) :: pinpl , tbar , uinpl , winpl
      real(8) , dimension(iym1,kzp1) :: s2c , ucfc11 , ucfc12 , uch4 ,&
           & uco211 , uco212 , uco213 , uco221 , uco222 , uco223 ,      &
           & un2o0 , un2o1 , uptype
      intent (in) bplnk , dw , k2 , kn , pinpl , pnew , s2c , tbar ,    &
                & tco2 , th2o , to3 , ucfc11 , ucfc12 , uch4 , uco211 , &
                & uco212 , uco213 , uco221 , uco222 , uco223 , uinpl ,  &
                & un2o0 , un2o1 , up2 , uptype , winpl
      intent (out) abstrc
!
! Local variables
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
      real(8) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4
      real(8) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 ,&
               & acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , &
               & b , dbeta01 , dbeta02 , dbeta03 , dbeta11 , dbetac ,   &
               & dbetc1 , dbetc2 , du01 , du02 , du03 , du1 , du11 ,    &
               & du12 , du13 , du2 , du21 , du22 , du23 , duch4 , p1 ,  &
               & phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 ,  &
               & tch4 , tlw , u , w1
      real(8) , dimension(iym1) :: ds2c , duptyp , rsqti , sqti , tt
      real(8) :: func
      integer :: i , l
      real(8) , dimension(iym1,6) :: tw
      data g1/0.0468556 , 0.0397454 , 0.0407664 , 0.0304380 ,           &
         & 0.0540398 , 0.0321962/
      data g2/14.4832 , 4.30242 , 5.23523 , 3.25342 , 0.698935 ,        &
         & 16.5599/
      data g3/26.1898 , 18.4476 , 15.3633 , 12.1927 , 9.14992 , 8.07092/
      data g4/0.0261782 , 0.0369516 , 0.0307266 , 0.0243854 ,           &
         & 0.0182932 , 0.0161418/
      data ab/3.0857E-2 , 2.3524E-2 , 1.7310E-2 , 2.6661E-2 ,           &
         & 2.8074E-2 , 2.2915E-2/
      data bb/ - 1.3512E-4 , -6.8320E-5 , -3.2609E-5 , -1.0228E-5 ,     &
         & -9.5743E-5 , -1.0304E-4/
      data abp/2.9129E-2 , 2.4101E-2 , 1.9821E-2 , 2.6904E-2 ,          &
         & 2.9458E-2 , 1.9892E-2/
      data bbp/ - 1.3139E-4 , -5.5688E-5 , -4.6380E-5 , -8.0362E-5 ,    &
         & -1.0115E-4 , -8.8061E-5/
!------------------------------------------------------------------
      func(u,b) = u/dsqrt(4.0+u*(1.0+1.0/b))
!
      do i = 1 , iym1
        sqti(i) = dsqrt(tbar(i,kn))
        rsqti(i) = 1./sqti(i)
!       h2o transmission
        tt(i) = dabs(tbar(i,kn)-250.0)
        ds2c(i) = dabs(s2c(i,k2+1)-s2c(i,k2))*uinpl(i,kn)
        duptyp(i) = dabs(uptype(i,k2+1)-uptype(i,k2))*uinpl(i,kn)
      end do
!
      do l = 1 , 6
        do i = 1 , iym1
          psi1 = dexp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
          phi1 = dexp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
          p1 = pnew(i)*(psi1/phi1)/sslp
          w1 = dw(i)*winpl(i,kn)*phi1
          tw(i,l) = dexp(-g1(l)*p1*(dsqrt(1.0+g2(l)*(w1/p1))-1.0)-g3(l) &
                  & *ds2c(i)-g4(l)*duptyp(i))
        end do
      end do
!
      do i = 1 , iym1
!
        du1 = dabs(ucfc11(i,k2+1)-ucfc11(i,k2))*winpl(i,kn)
        du2 = dabs(ucfc12(i,k2+1)-ucfc12(i,k2))*winpl(i,kn)
!       cfc transmissions
        tcfc3 = dexp(-175.005*du1)
        tcfc4 = dexp(-1202.18*du1)
        tcfc6 = dexp(-5786.73*du2)
        tcfc7 = dexp(-2873.51*du2)
        tcfc8 = dexp(-2085.59*du2)
!       Absorptivity for CFC11 bands
        acfc1 = 50.0*(1.0-dexp(-54.09*du1))*tw(i,1)*bplnk(7,i,kn)
        acfc2 = 60.0*(1.0-dexp(-5130.03*du1))*tw(i,2)*bplnk(8,i,kn)
        acfc3 = 60.0*(1.0-tcfc3)*tw(i,4)*tcfc6*bplnk(9,i,kn)
        acfc4 = 100.0*(1.0-tcfc4)*tw(i,5)*bplnk(10,i,kn)
!       Absorptivity for CFC12 bands
        acfc5 = 45.0*(1.0-dexp(-1272.35*du2))*tw(i,3)*bplnk(11,i,kn)
        acfc6 = 50.0*(1.0-tcfc6)*tw(i,4)*bplnk(12,i,kn)
        acfc7 = 80.0*(1.0-tcfc7)*tw(i,5)*tcfc4*bplnk(13,i,kn)
        acfc8 = 70.0*(1.0-tcfc8)*tw(i,6)*bplnk(14,i,kn)
!       Emissivity for CH4 band 1306 cm-1
        tlw = dexp(-1.0*dsqrt(up2(i)))
        duch4 = dabs(uch4(i,k2+1)-uch4(i,k2))*winpl(i,kn)
        dbetac = 2.94449*pinpl(i,kn)*rsqti(i)/sslp
        ach4 = 6.00444*sqti(i)*dlog(1.0+func(duch4,dbetac))             &
             & *tlw*bplnk(3,i,kn)
        tch4 = 1.0/(1.0+0.02*func(duch4,dbetac))
!       Absorptivity for N2O bands
        du01 = dabs(un2o0(i,k2+1)-un2o0(i,k2))*winpl(i,kn)
        du11 = dabs(un2o1(i,k2+1)-un2o1(i,k2))*winpl(i,kn)
        dbeta01 = 19.399*pinpl(i,kn)*rsqti(i)/sslp
        dbeta11 = dbeta01
!       1285 cm-1 band
        an2o1 = 2.35558*sqti(i)                                         &
              & *dlog(1.0+func(du01,dbeta01)+func(du11,dbeta11))        &
              & *tlw*tch4*bplnk(4,i,kn)
        du02 = 0.100090*du01
        du12 = 0.0992746*du11
        dbeta02 = 0.964282*dbeta01
!       589 cm-1 band
        an2o2 = 2.65581*sqti(i)                                         &
              & *dlog(1.0+func(du02,dbeta02)+func(du12,dbeta02))*tco2(i)&
              & *th2o(i)*bplnk(5,i,kn)
        du03 = 0.0333767*du01
        dbeta03 = 0.982143*dbeta01
!       1168 cm-1 band
        an2o3 = 2.54034*sqti(i)*dlog(1.0+func(du03,dbeta03))*tw(i,6)    &
              & *tcfc8*bplnk(6,i,kn)
!       Emissivity for 1064 cm-1 band of CO2
        du11 = dabs(uco211(i,k2+1)-uco211(i,k2))*winpl(i,kn)
        du12 = dabs(uco212(i,k2+1)-uco212(i,k2))*winpl(i,kn)
        du13 = dabs(uco213(i,k2+1)-uco213(i,k2))*winpl(i,kn)
        dbetc1 = 2.97558*pinpl(i,kn)*rsqti(i)/sslp
        dbetc2 = 2.0*dbetc1
        aco21 = 3.7571*sqti(i)                                          &
              & *dlog(1.0+func(du11,dbetc1)+func(du12,dbetc2)           &
              & +func(du13,dbetc2))*to3(i)*tw(i,5)                      &
              & *tcfc4*tcfc7*bplnk(2,i,kn)
!       Emissivity for 961 cm-1 band of co2
        du21 = dabs(uco221(i,k2+1)-uco221(i,k2))*winpl(i,kn)
        du22 = dabs(uco222(i,k2+1)-uco222(i,k2))*winpl(i,kn)
        du23 = dabs(uco223(i,k2+1)-uco223(i,k2))*winpl(i,kn)
        aco22 = 3.8443*sqti(i)                                          &
              & *dlog(1.0+func(du21,dbetc1)+func(du22,dbetc1)           &
              & +func(du23,dbetc2))*tw(i,4)*tcfc3*tcfc6*bplnk(1,i,kn)
!       total trace gas absorptivity
        abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 +     &
                  & acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +      &
                  & aco21 + aco22
      end do
      end subroutine trcabn
