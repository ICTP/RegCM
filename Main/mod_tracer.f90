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

      module mod_tracer

      use mod_constants
      use mod_dynparam

      private
!
      public :: trcmix , trcpth , trcab , trcabn , trcems , trcplk
      public :: cfc110 , cfc120 , ch40 , co2mmr , n2o0
!
      real(8) :: cfc110 , cfc120 , ch40 , co2mmr , n2o0
!
      contains

      subroutine trcmix(pmid,alat,coslat,n2o,ch4,cfc11,cfc12)
!
      implicit none
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
! alat   - current latitude in radians
! coslat - cosine of latitude
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
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & pmid
      real(8) , dimension(iym1) :: alat , coslat
      intent (in) alat , coslat , pmid
      intent (out) cfc11 , cfc12 , ch4 , n2o
!
! i      - longitude loop index
! k      - level index
! dlat   - latitude in degrees
! xn2o   - pressure scale height for n2o
! xch4   - pressure scale height for ch4
! xcfc11 - pressure scale height for cfc11
! xcfc12 - pressure scale height for cfc12
! ptrop  - pressure level of tropopause
! pratio - pressure divided by ptrop
!
      real(8) :: dlat , pratio , xcfc11 , xcfc12 , xch4 , xn2o
      integer :: i , k
      real(8) , dimension(iym1) :: ptrop
!
      xcfc11 = d_zero
      xcfc12 = d_zero
      xch4 = d_zero
      xn2o = d_zero
      do i = 1 , iym1
!       set stratospheric scale height factor for gases
        dlat = dabs(57.2958D0*alat(i))
        if ( dlat <= 45.0D0 ) then
          xn2o = 0.3478D0 + 0.00116D0*dlat
          xch4 = 0.2353D0
          xcfc11 = 0.7273D0 + 0.00606D0*dlat
          xcfc12 = 0.4000D0 + 0.00222D0*dlat
        else
          xn2o = 0.4000D0 + 0.013333D0*(dlat-45.0D0)
          xch4 = 0.2353D0 + 0.0225489D0*(dlat-45.0D0)
          xcfc11 = 1.00D0 + 0.013333D0*(dlat-45.0D0)
          xcfc12 = 0.50D0 + 0.024444D0*(dlat-45.0D0)
        end if
!
!       pressure of tropopause
        ptrop(i) = 250.0D2 - 150.0D2*coslat(i)**d_two
!
      end do
!
      do k = 1 , kz
        do i = 1 , iym1
          if ( pmid(i,k) >= ptrop(i) ) then
            ch4(i,k) = ch40
            n2o(i,k) = n2o0
            cfc11(i,k) = cfc110
            cfc12(i,k) = cfc120
          else
            pratio = pmid(i,k)/ptrop(i)
            ch4(i,k) = ch40*(pratio)**xch4
            n2o(i,k) = n2o0*(pratio)**xn2o
            cfc11(i,k) = cfc110*(pratio)**xcfc11
            cfc12(i,k) = cfc120*(pratio)**xcfc12
          end if
        end do
      end do

      end subroutine trcmix
!
      subroutine trcpth(tnm,pnm,cfc11,cfc12,n2o,ch4,qnm,ucfc11,ucfc12,  &
                      & un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,   &
                      & uco222,uco223,bn2o0,bn2o1,bch4,uptype)
!
      implicit none
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
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , pnm ,  &
           & ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,        &
           & uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & qnm , tnm
      intent (in) cfc11 , cfc12 , ch4 , n2o , pnm , qnm , tnm
      intent (inout) bch4 , bn2o0 , bn2o1 , ucfc11 , ucfc12 , uch4 ,    &
                   & uco211 , uco212 , uco213 , uco221 , uco222 ,       &
                   & uco223 , un2o0 , un2o1 , uptype
!
! i      - Longitude index
! k      - Level index
! co2fac - co2 factor
! alpha1 - stimulated emission term
! alpha2 - stimulated emission term
! rt     - reciprocal of local temperature
! rsqrt  - reciprocal of sqrt of temp
! pbar   - mean pressure
! dpnm   - difference in pressure
!
      real(8) , dimension(iym1) :: alpha1 , alpha2 , dpnm , pbar ,     &
                                  & rsqrt , rt
      real(8) , dimension(iym1,1) :: co2fac
      real(8) :: diff
      integer :: i , k
      data diff/1.66D0/           ! diffusivity factor

!-----------------------------------------------------------------------
!     Calculate path lengths for the trace gases
!-----------------------------------------------------------------------
      do i = 1 , iym1
        ucfc11(i,1) = 1.8D0*cfc11(i,1)*pnm(i,1)*rga
        ucfc12(i,1) = 1.8D0*cfc12(i,1)*pnm(i,1)*rga
        un2o0(i,1) = diff*1.02346D5*n2o(i,1)*pnm(i,1)         &
                   & *rga/dsqrt(tnm(i,1))
        un2o1(i,1) = diff*2.01909D0*un2o0(i,1)*dexp(-847.36D0/tnm(i,1))
        uch4(i,1) = diff*8.60957D4*ch4(i,1)*pnm(i,1)*rga/dsqrt(tnm(i,1))
        co2fac(i,1) = diff*co2mmr*pnm(i,1)*rga
        alpha1(i) = (d_one-dexp(-1540.0D0/ &
                     tnm(i,1)))**d_three/dsqrt(tnm(i,1))
        alpha2(i) = (d_one-dexp(-1360.0D0/ &
                     tnm(i,1)))**d_three/dsqrt(tnm(i,1))
        uco211(i,1) = 3.42217D3*co2fac(i,1)*alpha1(i)         &
                    & *dexp(-1849.7D0/tnm(i,1))
        uco212(i,1) = 6.02454D3*co2fac(i,1)*alpha1(i)         &
                    & *dexp(-2782.1D0/tnm(i,1))
        uco213(i,1) = 5.53143D3*co2fac(i,1)*alpha1(i)         &
                    & *dexp(-3723.2D0/tnm(i,1))
        uco221(i,1) = 3.88984D3*co2fac(i,1)*alpha2(i)         &
                    & *dexp(-1997.6D0/tnm(i,1))
        uco222(i,1) = 3.67108D3*co2fac(i,1)*alpha2(i)         &
                    & *dexp(-3843.8D0/tnm(i,1))
        uco223(i,1) = 6.50642D3*co2fac(i,1)*alpha2(i)         &
                    & *dexp(-2989.7D0/tnm(i,1))
        bn2o0(i,1) = diff*19.399D0*pnm(i,1)**d_two*n2o(i,1)   &
                   & *1.02346D5*rga/(sslp*tnm(i,1))
        bn2o1(i,1) = bn2o0(i,1)*dexp(-847.36D0/tnm(i,1))*2.06646D5
        bch4(i,1) = diff*2.94449D0*ch4(i,1)*pnm(i,1)          &
                  & **d_two*rga*8.60957D4/(sslp*tnm(i,1))
        uptype(i,1) = diff*qnm(i,1)*pnm(i,1)                  &
                      **d_two*dexp(1800.0D0*                  &
                     (d_one/tnm(i,1)-d_one/296.0D0))*rga/sslp
      end do
      do k = 1 , kz
        do i = 1 , iym1
          rt(i) = d_one/tnm(i,k)
          rsqrt(i) = dsqrt(rt(i))
          pbar(i) = d_half*(pnm(i,k+1)+pnm(i,k))/sslp
          dpnm(i) = (pnm(i,k+1)-pnm(i,k))*rga
          alpha1(i) = diff*rsqrt(i)* &
                    (d_one-dexp(-1540.0D0/tnm(i,k)))**d_three
          alpha2(i) = diff*rsqrt(i)* &
                    (d_one-dexp(-1360.0D0/tnm(i,k)))**d_three
          ucfc11(i,k+1) = ucfc11(i,k) + 1.8D0*cfc11(i,k)*dpnm(i)
          ucfc12(i,k+1) = ucfc12(i,k) + 1.8D0*cfc12(i,k)*dpnm(i)
          un2o0(i,k+1) = un2o0(i,k) + diff*1.02346D5*n2o(i,k)*rsqrt(i) &
                         *dpnm(i)
          un2o1(i,k+1) = un2o1(i,k) + diff*2.06646D5*n2o(i,k)*rsqrt(i) &
                         *dexp(-847.36D0/tnm(i,k))*dpnm(i)
          uch4(i,k+1) = uch4(i,k) + diff*8.60957D4*ch4(i,k)*rsqrt(i)   &
                        *dpnm(i)
          uco211(i,k+1) = uco211(i,k) + 1.15D0*3.42217D3*alpha1(i)     &
                          *co2mmr*dexp(-1849.7D0/tnm(i,k))*dpnm(i)
          uco212(i,k+1) = uco212(i,k) + 1.15D0*6.02454D3*alpha1(i)     &
                          *co2mmr*dexp(-2782.1D0/tnm(i,k))*dpnm(i)
          uco213(i,k+1) = uco213(i,k) + 1.15D0*5.53143D3*alpha1(i)     &
                          *co2mmr*dexp(-3723.2D0/tnm(i,k))*dpnm(i)
          uco221(i,k+1) = uco221(i,k) + 1.15D0*3.88984D3*alpha2(i)     &
                          *co2mmr*dexp(-1997.6D0/tnm(i,k))*dpnm(i)
          uco222(i,k+1) = uco222(i,k) + 1.15D0*3.67108D3*alpha2(i)     &
                          *co2mmr*dexp(-3843.8D0/tnm(i,k))*dpnm(i)
          uco223(i,k+1) = uco223(i,k) + 1.15D0*6.50642D3*alpha2(i)     &
                          *co2mmr*dexp(-2989.7D0/tnm(i,k))*dpnm(i)
          bn2o0(i,k+1) = bn2o0(i,k) + diff*19.399D0*pbar(i)*rt(i)      &
                         *1.02346D5*n2o(i,k)*dpnm(i)
          bn2o1(i,k+1) = bn2o1(i,k) + diff*19.399D0*pbar(i)*rt(i)      &
                         *2.06646D5*dexp(-847.36D0/tnm(i,k))*n2o(i,k)  &
                         *dpnm(i)
          bch4(i,k+1) = bch4(i,k) + diff*2.94449D0*rt(i)*pbar(i)       &
                        *8.60957D4*ch4(i,k)*dpnm(i)
          uptype(i,k+1) = uptype(i,k) + diff*qnm(i,k)                  &
                          *dexp(1800.0D0*(d_one/                       &
                          tnm(i,k)-d_one/296.0D0))*pbar(i)*dpnm(i)
        end do
      end do

      end subroutine trcpth
!
      subroutine trcab(k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,     &
                     & uco212,uco213,uco221,uco222,uco223,bn2o0,bn2o1,  &
                     & bch4,to3co2,pnm,dw,pnew,s2c,uptype,dplh2o,       &
                     & abplnk1,tco2,th2o,to3,abstrc)
!
      implicit none
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
      integer :: k1 , k2
      real(8) , dimension(14,iym1,kzp1) :: abplnk1
      real(8) , dimension(iym1) :: abstrc , dplh2o , dw , pnew , tco2 ,&
                                  & th2o , to3 , to3co2
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , pnm ,  &
           & s2c , ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,  &
           & uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
      intent (in) abplnk1 , bch4 , bn2o0 , bn2o1 , dplh2o , dw , k1 ,   &
                & k2 , pnew , pnm , s2c , tco2 , th2o , to3 , to3co2 ,  &
                & ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,   &
                & uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
      intent (out) abstrc
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
      real(8) , dimension(6) :: ab , abp , bb , bbp , g1 , g2 , g3 , g4
      real(8) :: acfc1 , acfc2 , acfc3 , acfc4 , acfc5 , acfc6 , acfc7 ,&
               & acfc8 , ach4 , aco21 , aco22 , an2o1 , an2o2 , an2o3 , &
               & dbeta01 , dbeta02 , dbeta03 , dbeta11 , dbetac ,       &
               & dbetc1 , dbetc2 , du01 , du02 , du03 , du1 , du11 ,    &
               & du12 , du13 , du2 , du21 , du22 , du23 , duch4 , p1 ,  &
               & phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 ,  &
               & tch4 , tlw , w1
      real(8) , dimension(iym1) :: ds2c , duptyp , sqti , tt
      integer :: i , l
      real(8) , dimension(iym1,6) :: tw
!
      data g1/0.0468556D0 , 0.0397454D0 , 0.0407664D0 , 0.0304380D0 , &
              0.0540398D0 , 0.0321962D0/
      data g2/14.4832D0 , 4.30242D0 , 5.23523D0 , 3.25342D0 , &
               0.698935D0 , 16.5599D0/
      data g3/26.1898D0 , 18.4476D0 , 15.3633D0 , 12.1927D0 , &
               9.14992D0 , 8.07092D0/
      data g4/0.0261782D0 , 0.0369516D0 , 0.0307266D0 , 0.0243854D0 , &
              0.0182932D0 , 0.0161418D0/
      data ab/3.0857D-2 , 2.3524D-2 , 1.7310D-2 , 2.6661D-2 , &
              2.8074D-2 , 2.2915D-2/
      data bb/ -1.3512D-4 , -6.8320D-5 , -3.2609D-5 , -1.0228D-5 ,  &
               -9.5743D-5 , -1.0304D-4/
      data abp/2.9129D-2 , 2.4101D-2 , 1.9821D-2 , 2.6904D-2 ,  &
               2.9458D-2 , 1.9892D-2/
      data bbp/ -1.3139D-4 , -5.5688D-5 , -4.6380D-5 , -8.0362D-5 , &
                -1.0115D-4 , -8.8061D-5/
!------------------------------------------------------------------
      do i = 1 , iym1
        sqti(i) = dsqrt(to3co2(i))
!       h2o transmission
        tt(i) = dabs(to3co2(i)-250.0D0)
        ds2c(i) = dabs(s2c(i,k1)-s2c(i,k2))
        duptyp(i) = dabs(uptype(i,k1)-uptype(i,k2))
      end do
!
      do l = 1 , 6
        do i = 1 , iym1
          psi1 = dexp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
          phi1 = dexp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
          p1 = pnew(i)*(psi1/phi1)/sslp
          w1 = dw(i)*phi1
          tw(i,l) = dexp(-g1(l)*p1*(dsqrt(d_one+g2(l)*(w1/p1))-&
                    d_one)-g3(l)*ds2c(i)-g4(l)*duptyp(i))
        end do
      end do
!
      do i = 1 , iym1
        du1 = dabs(ucfc11(i,k1)-ucfc11(i,k2))
        du2 = dabs(ucfc12(i,k1)-ucfc12(i,k2))
!       cfc transmissions
        tcfc3 = dexp(-175.005D0*du1)
        tcfc4 = dexp(-1202.18D0*du1)
        tcfc6 = dexp(-5786.73D0*du2)
        tcfc7 = dexp(-2873.51D0*du2)
        tcfc8 = dexp(-2085.59D0*du2)
!       Absorptivity for CFC11 bands
        acfc1 = 50.0D0*(d_one-dexp(-54.09D0*du1))*tw(i,1)*abplnk1(7,i,k2)
        acfc2 = 60.0D0*(d_one-dexp(-5130.03D0*du1))*tw(i,2)*abplnk1(8,i,k2)
        acfc3 = 60.0D0*(d_one-tcfc3)*tw(i,4)*tcfc6*abplnk1(9,i,k2)
        acfc4 = 100.0D0*(d_one-tcfc4)*tw(i,5)*abplnk1(10,i,k2)
!       Absorptivity for CFC12 bands
        acfc5 = 45.0D0*(d_one-dexp(-1272.35D0*du2))*tw(i,3)*abplnk1(11,i,k2)
        acfc6 = 50.0D0*(d_one-tcfc6)*tw(i,4)*abplnk1(12,i,k2)
        acfc7 = 80.0D0*(d_one-tcfc7)*tw(i,5)*tcfc4*abplnk1(13,i,k2)
        acfc8 = 70.0D0*(d_one-tcfc8)*tw(i,6)*abplnk1(14,i,k2)
!       Emissivity for CH4 band 1306 cm-1
        tlw = dexp(-d_one*dsqrt(dplh2o(i)))
        duch4 = dabs(uch4(i,k1)-uch4(i,k2))
        dbetac = dabs(bch4(i,k1)-bch4(i,k2))/duch4
        ach4 = 6.00444D0*sqti(i)*dlog(d_one+func(duch4,dbetac)) &
             & *tlw*abplnk1(3,i,k2)
        tch4 = d_one/(d_one+0.02D0*func(duch4,dbetac))
!       Absorptivity for N2O bands
        du01 = dabs(un2o0(i,k1)-un2o0(i,k2))
        du11 = dabs(un2o1(i,k1)-un2o1(i,k2))
        dbeta01 = dabs(bn2o0(i,k1)-bn2o0(i,k2))/du01
        dbeta11 = dabs(bn2o1(i,k1)-bn2o1(i,k2))/du11
!       1285 cm-1 band
        an2o1 = 2.35558D0*sqti(i)                                  &
              & *dlog(d_one+func(du01,dbeta01)+func(du11,dbeta11)) &
              & *tlw*tch4*abplnk1(4,i,k2)
        du02 = 0.100090D0*du01
        du12 = 0.0992746D0*du11
        dbeta02 = 0.964282D0*dbeta01
!       589 cm-1 band
        an2o2 = 2.65581D0*sqti(i)    &
              & *dlog(d_one+func(du02,dbeta02)+ &
              & func(du12,dbeta02))*th2o(i)*tco2(i)*abplnk1(5,i,k2)
        du03 = 0.0333767D0*du01
        dbeta03 = 0.982143D0*dbeta01
!       1168 cm-1 band
        an2o3 = 2.54034D0*sqti(i)*dlog(d_one+func(du03,dbeta03))* &
                tw(i,6)*tcfc8*abplnk1(6,i,k2)
!       Emissivity for 1064 cm-1 band of CO2
        du11 = dabs(uco211(i,k1)-uco211(i,k2))
        du12 = dabs(uco212(i,k1)-uco212(i,k2))
        du13 = dabs(uco213(i,k1)-uco213(i,k2))
        dbetc1 = 2.97558D0*dabs(pnm(i,k1)+pnm(i,k2))/(d_two*sslp*sqti(i))
        dbetc2 = d_two*dbetc1
        aco21 = 3.7571D0*sqti(i)                                  &
              & *dlog(d_one+func(du11,dbetc1)+func(du12,dbetc2)   &
              & +func(du13,dbetc2))*to3(i)*tw(i,5)                &
              & *tcfc4*tcfc7*abplnk1(2,i,k2)
!       Emissivity for 961 cm-1 band
        du21 = dabs(uco221(i,k1)-uco221(i,k2))
        du22 = dabs(uco222(i,k1)-uco222(i,k2))
        du23 = dabs(uco223(i,k1)-uco223(i,k2))
        aco22 = 3.8443D0*sqti(i)                                   &
              & *dlog(d_one+func(du21,dbetc1)+func(du22,dbetc1)    &
              & +func(du23,dbetc2))*tw(i,4)*tcfc3*tcfc6*abplnk1(1,i,k2)
!       total trace gas absorptivity
        abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 +  &
                  & acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +   &
                  & aco21 + aco22
      end do
      end subroutine trcab
!
      subroutine trcabn(k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,    &
                      & uco212,uco213,uco221,uco222,uco223,tbar,bplnk,  &
                      & winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,up2,    &
                      & pnew,abstrc,uinpl)
!
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
               & dbeta01 , dbeta02 , dbeta03 , dbeta11 , dbetac ,       &
               & dbetc1 , dbetc2 , du01 , du02 , du03 , du1 , du11 ,    &
               & du12 , du13 , du2 , du21 , du22 , du23 , duch4 , p1 ,  &
               & phi1 , psi1 , tcfc3 , tcfc4 , tcfc6 , tcfc7 , tcfc8 ,  &
               & tch4 , tlw , w1
      real(8) , dimension(iym1) :: ds2c , duptyp , rsqti , sqti , tt
      integer :: i , l
      real(8) , dimension(iym1,6) :: tw
      data g1/0.0468556D0 , 0.0397454D0 , 0.0407664D0 , 0.0304380D0 , &
              0.0540398D0 , 0.0321962D0/
      data g2/14.4832D0 , 4.30242D0 , 5.23523D0 , 3.25342D0 , &
               0.698935D0 , 16.5599D0/
      data g3/26.1898D0 , 18.4476D0 , 15.3633D0 , 12.1927D0 , &
               9.14992D0 , 8.07092D0/
      data g4/0.0261782D0 , 0.0369516D0 , 0.0307266D0 , 0.0243854D0 , &
              0.0182932D0 , 0.0161418D0/
      data ab/3.0857D-2 , 2.3524D-2 , 1.7310D-2 , 2.6661D-2 , &
              2.8074D-2 , 2.2915D-2/
      data bb/ -1.3512D-4 , -6.8320D-5 , -3.2609D-5 , -1.0228D-5 , &
               -9.5743D-5 , -1.0304D-4/
      data abp/2.9129D-2 , 2.4101D-2 , 1.9821D-2 , 2.6904D-2 , &
               2.9458D-2 , 1.9892D-2/
      data bbp/ -1.3139D-4 , -5.5688D-5 , -4.6380D-5 , -8.0362D-5 , &
                -1.0115D-4 , -8.8061D-5/
!------------------------------------------------------------------
!
      do i = 1 , iym1
        sqti(i) = dsqrt(tbar(i,kn))
        rsqti(i) = d_one/sqti(i)
!       h2o transmission
        tt(i) = dabs(tbar(i,kn)-250.0D0)
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
          tw(i,l) = dexp(-g1(l)*p1*(dsqrt(d_one+g2(l)* &
                    (w1/p1))-d_one)-g3(l) &
                    *ds2c(i)-g4(l)*duptyp(i))
        end do
      end do
!
      do i = 1 , iym1
!
        du1 = dabs(ucfc11(i,k2+1)-ucfc11(i,k2))*winpl(i,kn)
        du2 = dabs(ucfc12(i,k2+1)-ucfc12(i,k2))*winpl(i,kn)
!       cfc transmissions
        tcfc3 = dexp(-175.005D0*du1)
        tcfc4 = dexp(-1202.18D0*du1)
        tcfc6 = dexp(-5786.73D0*du2)
        tcfc7 = dexp(-2873.51D0*du2)
        tcfc8 = dexp(-2085.59D0*du2)
!       Absorptivity for CFC11 bands
        acfc1 = 50.0D0*(d_one-dexp(-54.09D0*du1))*tw(i,1)*bplnk(7,i,kn)
        acfc2 = 60.0D0*(d_one-dexp(-5130.03D0*du1))*tw(i,2)*bplnk(8,i,kn)
        acfc3 = 60.0D0*(d_one-tcfc3)*tw(i,4)*tcfc6*bplnk(9,i,kn)
        acfc4 = 100.0D0*(d_one-tcfc4)*tw(i,5)*bplnk(10,i,kn)
!       Absorptivity for CFC12 bands
        acfc5 = 45.0D0*(d_one-dexp(-1272.35D0*du2))*tw(i,3)*bplnk(11,i,kn)
        acfc6 = 50.0D0*(d_one-tcfc6)*tw(i,4)*bplnk(12,i,kn)
        acfc7 = 80.0D0*(d_one-tcfc7)*tw(i,5)*tcfc4*bplnk(13,i,kn)
        acfc8 = 70.0D0*(d_one-tcfc8)*tw(i,6)*bplnk(14,i,kn)
!       Emissivity for CH4 band 1306 cm-1
        tlw = dexp(-d_one*dsqrt(up2(i)))
        duch4 = dabs(uch4(i,k2+1)-uch4(i,k2))*winpl(i,kn)
        dbetac = 2.94449D0*pinpl(i,kn)*rsqti(i)/sslp
        ach4 = 6.00444D0*sqti(i)*dlog(d_one+func(duch4,dbetac))  &
             & *tlw*bplnk(3,i,kn)
        tch4 = d_one/(d_one+0.02D0*func(duch4,dbetac))
!       Absorptivity for N2O bands
        du01 = dabs(un2o0(i,k2+1)-un2o0(i,k2))*winpl(i,kn)
        du11 = dabs(un2o1(i,k2+1)-un2o1(i,k2))*winpl(i,kn)
        dbeta01 = 19.399D0*pinpl(i,kn)*rsqti(i)/sslp
        dbeta11 = dbeta01
!       1285 cm-1 band
        an2o1 = 2.35558D0*sqti(i)                                  &
              & *dlog(d_one+func(du01,dbeta01)+func(du11,dbeta11)) &
              & *tlw*tch4*bplnk(4,i,kn)
        du02 = 0.100090D0*du01
        du12 = 0.0992746D0*du11
        dbeta02 = 0.964282D0*dbeta01
!       589 cm-1 band
        an2o2 = 2.65581D0*sqti(i)                                  &
              & *dlog(d_one+func(du02,dbeta02)+func(du12,dbeta02))*&
              & tco2(i)*th2o(i)*bplnk(5,i,kn)
        du03 = 0.0333767D0*du01
        dbeta03 = 0.982143D0*dbeta01
!       1168 cm-1 band
        an2o3 = 2.54034D0*sqti(i)*dlog(d_one+func(du03,dbeta03))*tw(i,6)    &
              & *tcfc8*bplnk(6,i,kn)
!       Emissivity for 1064 cm-1 band of CO2
        du11 = dabs(uco211(i,k2+1)-uco211(i,k2))*winpl(i,kn)
        du12 = dabs(uco212(i,k2+1)-uco212(i,k2))*winpl(i,kn)
        du13 = dabs(uco213(i,k2+1)-uco213(i,k2))*winpl(i,kn)
        dbetc1 = 2.97558D0*pinpl(i,kn)*rsqti(i)/sslp
        dbetc2 = d_two*dbetc1
        aco21 = 3.7571D0*sqti(i)                                   &
              & *dlog(d_one+func(du11,dbetc1)+func(du12,dbetc2)    &
              & +func(du13,dbetc2))*to3(i)*tw(i,5)                 &
              & *tcfc4*tcfc7*bplnk(2,i,kn)
!       Emissivity for 961 cm-1 band of co2
        du21 = dabs(uco221(i,k2+1)-uco221(i,k2))*winpl(i,kn)
        du22 = dabs(uco222(i,k2+1)-uco222(i,k2))*winpl(i,kn)
        du23 = dabs(uco223(i,k2+1)-uco223(i,k2))*winpl(i,kn)
        aco22 = 3.8443D0*sqti(i)                                   &
              & *dlog(d_one+func(du21,dbetc1)+func(du22,dbetc1)    &
              & +func(du23,dbetc2))*tw(i,4)*tcfc3*tcfc6*bplnk(1,i,kn)
!       total trace gas absorptivity
        abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  & acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 +  &
                  & aco21 + aco22
      end do
      end subroutine trcabn
!
      subroutine trcplk(tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
      implicit none
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
      real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
      real(8) , dimension(14,iym1) :: emplnk
      real(8) , dimension(iym1,kzp1) :: tint , tlayr
      real(8) , dimension(iym1) :: tplnke
      intent (in) tint , tlayr , tplnke
      intent (out) abplnk1 , abplnk2 , emplnk
!
! wvl   - wavelength index
! f1    - Planck function factor
! f2    -       "
! f3    -       "
!
      real(8) , dimension(14) :: f1 , f2 , f3
      integer :: i , k , wvl
!
      data f1/5.85713D8 , 7.94950D8 , 1.47009D9 , 1.40031D9 , &
              1.34853D8 , 1.05158D9 , 3.35370D8 , 3.99601D8 , &
              5.35994D8 , 8.42955D8 , 4.63682D8 , 5.18944D8 , &
              8.83202D8 , 1.03279D9/
      data f2/2.02493D11 , 3.04286D11 , 6.90698D11 , 6.47333D11 , &
              2.85744D10 , 4.41862D11 , 9.62780D10 , 1.21618D11 , &
              1.79905D11 , 3.29029D11 , 1.48294D11 , 1.72315D11 , &
              3.50140D11 , 4.31364D11/
      data f3/1383.0D0 , 1531.0D0 , 1879.0D0 , 1849.0D0 ,  848.0D0 , &
              1681.0D0 , 1148.0D0 , 1217.0D0 , 1343.0D0 , 1561.0D0 , &
              1279.0D0 , 1328.0D0 , 1586.0D0 , 1671.0D0/
!
!     Calculate emissivity Planck factor
!
      do wvl = 1 , 14
        do i = 1 , iym1
          emplnk(wvl,i) = f1(wvl)                       &
                        & /(tplnke(i)**d_four*           &
                        (dexp(f3(wvl)/tplnke(i))-d_one))
        end do
      end do
!
!     Calculate absorptivity Planck factor for tint and tlayr temperatures
!
      do wvl = 1 , 14
        do k = 1 , kzp1
          do i = 1 , iym1
!           non-nearlest layer function
            abplnk1(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tint(i,k)))        &
                             & /(tint(i,k)**5.0D0*                      &
                             & (dexp(f3(wvl)/tint(i,k))-d_one)**d_two)
!           nearest layer function
            abplnk2(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tlayr(i,k)))       &
                             & /(tlayr(i,k)**5.0D0*                     &
                             & (dexp(f3(wvl)/tlayr(i,k))-d_one)**d_two)
          end do
        end do
      end do

      end subroutine trcplk
!
      subroutine trcems(k,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,bn2o0,     &
                      & bn2o1,uch4,bch4,uco211,uco212,uco213,uco221,    &
                      & uco222,uco223,uptype,w,s2c,up2,emplnk,th2o,tco2,&
                      & to3,emstrc)
!
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
      integer :: k
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2t , &
           & emstrc , pnm , s2c , ucfc11 , ucfc12 , uch4 , uco211 ,     &
           & uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,       &
           & un2o1 , uptype , w
      real(8) , dimension(14,iym1) :: emplnk
      real(8) , dimension(iym1) :: tco2 , th2o , to3 , up2
      intent (in) bch4 , bn2o0 , bn2o1 , co2t , emplnk , k , pnm , s2c ,&
                & tco2 , th2o , to3 , ucfc11 , ucfc12 , uch4 , uco211 , &
                & uco212 , uco213 , uco221 , uco222 , uco223 , un2o0 ,  &
                & un2o1 , up2 , uptype , w
      intent (out) emstrc
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
      real(8) :: beta01 , beta02 , beta03 , beta11 , betac ,        &
               & betac1 , betac2 , ecfc1 , ecfc2 , ecfc3 , ecfc4 ,      &
               & ecfc5 , ecfc6 , ecfc7 , ecfc8 , ech4 , eco21 , eco22 , &
               & en2o1 , en2o2 , en2o3 , p1 , phi1 , psi1 , tcfc3 ,     &
               & tcfc4 , tcfc6 , tcfc7 , tcfc8 , tch4 , tlw , u01 , &
               & u02 , u03 , u11 , u12 , w1
      integer :: i , l
      real(8) , dimension(iym1) :: sqti , tt
      real(8) , dimension(iym1,6) :: tw
!
      data g1/0.0468556D0 , 0.0397454D0 , 0.0407664D0 , 0.0304380D0 , &
              0.0540398D0 , 0.0321962D0/
      data g2/14.4832D0 ,    4.30242D0 ,  5.23523D0 , 3.25342D0 ,     &
               0.698935D0 , 16.5599D0/
      data g3/26.1898D0 , 18.4476D0 , 15.3633D0 , 12.1927D0 ,         &
               9.14992D0 , 8.07092D0/
      data g4/0.0261782D0 , 0.0369516D0 , 0.0307266D0 , 0.0243854D0 , &
              0.0182932D0 , 0.0161418D0/
      data ab/3.0857D-2 , 2.3524D-2 , 1.7310D-2 , 2.6661D-2 ,         &
              2.8074D-2 , 2.2915D-2/
      data bb/ -1.3512D-4 , -6.8320D-5 , -3.2609D-5 , -1.0228D-5 ,    &
               -9.5743D-5 , -1.0304D-4/
      data abp/2.9129D-2 , 2.4101D-2 , 1.9821D-2 , 2.6904D-2 ,       &
               2.9458D-2 , 1.9892D-2/
      data bbp/ -1.3139D-4 , -5.5688D-5 , -4.6380D-5 , -8.0362D-5 ,    &
                -1.0115D-4 , -8.8061D-5/
!
      do i = 1 , iym1
        sqti(i) = dsqrt(co2t(i,k))
!       Transmission for h2o
        tt(i) = dabs(co2t(i,k)-250.0D0)
      end do
!
      do l = 1 , 6
        do i = 1 , iym1
          psi1 = dexp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
          phi1 = dexp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
          p1 = pnm(i,k)*(psi1/phi1)/sslp
          w1 = w(i,k)*phi1
          tw(i,l) = dexp(-g1(l)*p1*(dsqrt(d_one+g2(l)*(w1/p1))-&
                        d_one)-g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
        end do
      end do
!
      do i = 1 , iym1
!       transmission due to cfc bands
        tcfc3 = dexp(-175.005D0*ucfc11(i,k))
        tcfc4 = dexp(-1202.18D0*ucfc11(i,k))
        tcfc6 = dexp(-5786.73D0*ucfc12(i,k))
        tcfc7 = dexp(-2873.51D0*ucfc12(i,k))
        tcfc8 = dexp(-2085.59D0*ucfc12(i,k))
!       Emissivity for CFC11 bands
        ecfc1 = 50.0D0*(d_one-dexp(-54.09D0*ucfc11(i,k)))*tw(i,1)* &
                        emplnk(7,i)
        ecfc2 = 60.0D0*(d_one-dexp(-5130.03D0*ucfc11(i,k)))*tw(i,2)  &
              & *emplnk(8,i)
        ecfc3 = 60.0D0*(d_one-tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
        ecfc4 = 100.0D0*(d_one-tcfc4)*tw(i,5)*emplnk(10,i)
!       Emissivity for CFC12 bands
        ecfc5 = 45.0D0*(d_one-dexp(-1272.35D0*ucfc12(i,k)))*tw(i,3)  &
              & *emplnk(11,i)
        ecfc6 = 50.0D0*(d_one-tcfc6)*tw(i,4)*emplnk(12,i)
        ecfc7 = 80.0D0*(d_one-tcfc7)*tw(i,5)*tcfc4*emplnk(13,i)
        ecfc8 = 70.0D0*(d_one-tcfc8)*tw(i,6)*emplnk(14,i)
!       Emissivity for CH4 band 1306 cm-1
        tlw = dexp(-d_one*dsqrt(up2(i)))
        betac = bch4(i,k)/uch4(i,k)
        ech4 = 6.00444D0*sqti(i)*dlog(d_one+func(uch4(i,k),betac))   &
             & *tlw*emplnk(3,i)
        tch4 = d_one/(d_one+0.02D0*func(uch4(i,k),betac))
!       Emissivity for N2O bands
        u01 = un2o0(i,k)
        u11 = un2o1(i,k)
        beta01 = bn2o0(i,k)/un2o0(i,k)
        beta11 = bn2o1(i,k)/un2o1(i,k)
!       1285 cm-1 band
        en2o1 = 2.35558D0*sqti(i)                                      &
              & *dlog(d_one+func(u01,beta01)+func(u11,beta11))         &
              & *tlw*tch4*emplnk(4,i)
        u02 = 0.100090D0*u01
        u12 = 0.0992746D0*u11
        beta02 = 0.964282D0*beta01
!       589 cm-1 band
        en2o2 = 2.65581D0*sqti(i)                                      &
              & *dlog(d_one+func(u02,beta02)+func(u12,beta02))*tco2(i) &
              & *th2o(i)*emplnk(5,i)
        u03 = 0.0333767D0*u01
        beta03 = 0.982143D0*beta01
!       1168 cm-1 band
        en2o3 = 2.54034D0*sqti(i)*dlog(d_one+func(u03,beta03))*tw(i,6) &
              & *tcfc8*emplnk(6,i)
!       Emissivity for 1064 cm-1 band of CO2
        betac1 = 2.97558D0*pnm(i,k)/(sslp*sqti(i))
        betac2 = d_two*betac1
        eco21 = 3.7571D0*sqti(i)                                       &
              & *dlog(d_one+func(uco211(i,k),betac1)+func(uco212(i,k), &
              & betac2)+func(uco213(i,k),betac2))*to3(i)*tw(i,5)       &
              & *tcfc4*tcfc7*emplnk(2,i)
!       Emissivity for 961 cm-1 band
        eco22 = 3.8443D0*sqti(i)                                       &
              & *dlog(d_one+func(uco221(i,k),betac1)+func(uco222(i,k), &
              & betac1)+func(uco223(i,k),betac2))*tw(i,4)              &
              & *tcfc3*tcfc6*emplnk(1,i)
!       total trace gas emissivity
        emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 + ecfc6 +  &
                    & ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 +   &
                    & eco21 + eco22
      end do
      end subroutine trcems
!
      function func(u,b)
        implicit none
        real(8) :: func
        real(8) , intent(in) :: u , b
        func = u/dsqrt(d_four+u*(d_one+d_one/b))
      end function func
!
      end module mod_tracer
