!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine trcems(k,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,bn2o0,     &
                      & bn2o1,uch4,bch4,uco211,uco212,uco213,uco221,    &
                      & uco222,uco223,uptype,w,s2c,up2,emplnk,th2o,tco2,&
                      & to3,emstrc)
!
      use mod_regcm_param
      use mod_parrad
      use mod_crdcon
      implicit none
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
! Dummy arguments
!
      integer :: k
      real(8) , dimension(plond,plevp) :: bch4 , bn2o0 , bn2o1 , co2t , &
           & emstrc , pnm , s2c , ucfc11 , ucfc12 , uch4 , uco211 ,     &
           & uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,       &
           & un2o1 , uptype , w
      real(8) , dimension(14,plond) :: emplnk
      real(8) , dimension(plond) :: tco2 , th2o , to3 , up2
      intent (in) bch4 , bn2o0 , bn2o1 , co2t , emplnk , k , pnm , s2c ,&
                & tco2 , th2o , to3 , ucfc11 , ucfc12 , uch4 , uco211 , &
                & uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,  &
                & un2o1 , up2 , uptype , w
      intent (out) emstrc
!
! Local variables
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
      real(8) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4
      real(8) :: b , beta01 , beta02 , beta03 , beta11 , betac ,        &
               & betac1 , betac2 , ecfc1 , ecfc2 , ecfc3 , ecfc4 ,      &
               & ecfc5 , ecfc6 , ecfc7 , ecfc8 , ech4 , eco21 , eco22 , &
               & en2o1 , en2o2 , en2o3 , p1 , phi1 , psi1 , tcfc3 ,     &
               & tcfc4 , tcfc6 , tcfc7 , tcfc8 , tch4 , tlw , u , u01 , &
               & u02 , u03 , u11 , u12 , w1
      real(8) :: func
      integer :: i , l
      real(8) , dimension(plond) :: sqti , tt
      real(8) , dimension(plond,6) :: tw
!
      data g1/0.0468556D0 , 0.0397454D0 , 0.0407664D0 , 0.0304380D0 ,   &
         & 0.0540398D0 , 0.0321962D0/
      data g2/14.4832D0 , 4.30242D0 , 5.23523D0 , 3.25342D0 ,           &
         & 0.698935D0 , 16.5599D0/
      data g3/26.1898D0 , 18.4476D0 , 15.3633D0 , 12.1927D0 ,           &
         & 9.14992D0 , 8.07092D0/
      data g4/0.0261782D0 , 0.0369516D0 , 0.0307266D0 , 0.0243854D0 ,   &
         & 0.0182932D0 , 0.0161418D0/
      data ab/3.0857D-2 , 2.3524D-2 , 1.7310D-2 , 2.6661D-2 ,           &
         & 2.8074D-2 , 2.2915D-2/
      data bb/ - 1.3512D-4 , -6.8320D-5 , -3.2609D-5 , -1.0228D-5 ,     &
         & -9.5743D-5 , -1.0304D-4/
      data abp/2.9129D-2 , 2.4101D-2 , 1.9821D-2 , 2.6904D-2 ,          &
         & 2.9458D-2 , 1.9892D-2/
      data bbp/ - 1.3139D-4 , -5.5688D-5 , -4.6380D-5 , -8.0362D-5 ,    &
         & -1.0115D-4 , -8.8061D-5/
      func(u,b) = u/dsqrt(4.0+u*(1.0+1.0/b))
!
      do i = 1 , plon
        sqti(i) = dsqrt(co2t(i,k))
!       Transmission for h2o
        tt(i) = dabs(co2t(i,k)-250.0)
      end do
!
      do l = 1 , 6
        do i = 1 , plon
          psi1 = dexp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
          phi1 = dexp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
          p1 = pnm(i,k)*(psi1/phi1)/sslp
          w1 = w(i,k)*phi1
          tw(i,l) = dexp(-g1(l)*p1*(dsqrt(1.0+g2(l)*(w1/p1))-1.0)-g3(l) &
                  & *s2c(i,k)-g4(l)*uptype(i,k))
        end do
      end do
!
      do i = 1 , plon
!       transmission due to cfc bands
        tcfc3 = dexp(-175.005*ucfc11(i,k))
        tcfc4 = dexp(-1202.18*ucfc11(i,k))
        tcfc6 = dexp(-5786.73*ucfc12(i,k))
        tcfc7 = dexp(-2873.51*ucfc12(i,k))
        tcfc8 = dexp(-2085.59*ucfc12(i,k))
!       Emissivity for CFC11 bands
        ecfc1 = 50.0*(1.0-dexp(-54.09*ucfc11(i,k)))*tw(i,1)*emplnk(7,i)
        ecfc2 = 60.0*(1.0-dexp(-5130.03*ucfc11(i,k)))*tw(i,2)           &
              & *emplnk(8,i)
        ecfc3 = 60.0*(1.0-tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
        ecfc4 = 100.0*(1.0-tcfc4)*tw(i,5)*emplnk(10,i)
!       Emissivity for CFC12 bands
        ecfc5 = 45.0*(1.0-dexp(-1272.35*ucfc12(i,k)))*tw(i,3)           &
              & *emplnk(11,i)
        ecfc6 = 50.0*(1.0-tcfc6)*tw(i,4)*emplnk(12,i)
        ecfc7 = 80.0*(1.0-tcfc7)*tw(i,5)*tcfc4*emplnk(13,i)
        ecfc8 = 70.0*(1.0-tcfc8)*tw(i,6)*emplnk(14,i)
!       Emissivity for CH4 band 1306 cm-1
        tlw = dexp(-1.0*dsqrt(up2(i)))
        betac = bch4(i,k)/uch4(i,k)
        ech4 = 6.00444*sqti(i)*dlog(1.0+func(uch4(i,k),betac))          &
             & *tlw*emplnk(3,i)
        tch4 = 1.0/(1.0+0.02*func(uch4(i,k),betac))
!       Emissivity for N2O bands
        u01 = un2o0(i,k)
        u11 = un2o1(i,k)
        beta01 = bn2o0(i,k)/un2o0(i,k)
        beta11 = bn2o1(i,k)/un2o1(i,k)
!       1285 cm-1 band
        en2o1 = 2.35558*sqti(i)                                         &
              & *dlog(1.0+func(u01,beta01)+func(u11,beta11))            &
              & *tlw*tch4*emplnk(4,i)
        u02 = 0.100090*u01
        u12 = 0.0992746*u11
        beta02 = 0.964282*beta01
!       589 cm-1 band
        en2o2 = 2.65581*sqti(i)                                         &
              & *dlog(1.0+func(u02,beta02)+func(u12,beta02))*tco2(i)    &
              & *th2o(i)*emplnk(5,i)
        u03 = 0.0333767*u01
        beta03 = 0.982143*beta01
!       1168 cm-1 band
        en2o3 = 2.54034*sqti(i)*dlog(1.0+func(u03,beta03))*tw(i,6)      &
              & *tcfc8*emplnk(6,i)
!       Emissivity for 1064 cm-1 band of CO2
        betac1 = 2.97558*pnm(i,k)/(sslp*sqti(i))
        betac2 = 2.0*betac1
        eco21 = 3.7571*sqti(i)                                          &
              & *dlog(1.0+func(uco211(i,k),betac1)+func(uco212(i,k),    &
              & betac2)+func(uco213(i,k),betac2))*to3(i)*tw(i,5)        &
              & *tcfc4*tcfc7*emplnk(2,i)
!       Emissivity for 961 cm-1 band
        eco22 = 3.8443*sqti(i)                                          &
              & *dlog(1.0+func(uco221(i,k),betac1)+func(uco222(i,k),    &
              & betac1)+func(uco223(i,k),betac2))*tw(i,4)               &
              & *tcfc3*tcfc6*emplnk(1,i)
!       total trace gas emissivity
        emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 + ecfc6 +   &
                    & ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 +    &
                    & eco21 + eco22
      end do
      end subroutine trcems
