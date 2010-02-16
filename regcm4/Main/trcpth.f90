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
 
      subroutine trcpth(tnm,pnm,cfc11,cfc12,n2o,ch4,qnm,ucfc11,ucfc12,  &
                      & un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,   &
                      & uco222,uco223,bn2o0,bn2o1,bch4,uptype)
!
      use mod_regcm_param
      use mod_crdcon
      use mod_tracer
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
! Dummy arguments
!
      real(8) , dimension(ixm1,kxp1) :: bch4 , bn2o0 , bn2o1 , pnm ,  &
           & ucfc11 , ucfc12 , uch4 , uco211 , uco212 , uco213 ,        &
           & uco221 , uco222 , uco223 , un2o0 , un2o1 , uptype
      real(8) , dimension(ixm1,kx) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & qnm , tnm
      intent (in) cfc11 , cfc12 , ch4 , n2o , pnm , qnm , tnm
      intent (inout) bch4 , bn2o0 , bn2o1 , ucfc11 , ucfc12 , uch4 ,    &
                   & uco211 , uco212 , uco213 , uco221 , uco222 ,       &
                   & uco223 , un2o0 , un2o1 , uptype
!
! Local variables
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
      real(8) , dimension(ixm1) :: alpha1 , alpha2 , dpnm , pbar ,     &
                                  & rsqrt , rt
      real(8) , dimension(ixm1,1) :: co2fac
      real(8) :: diff
      integer :: i , k
      data diff/1.66/           ! diffusivity factor

!-----------------------------------------------------------------------
!     Calculate path lengths for the trace gases
!-----------------------------------------------------------------------
      do i = 1 , ixm1
        ucfc11(i,1) = 1.8*cfc11(i,1)*pnm(i,1)*rga
        ucfc12(i,1) = 1.8*cfc12(i,1)*pnm(i,1)*rga
        un2o0(i,1) = diff*1.02346E5*n2o(i,1)*pnm(i,1)                   &
                   & *rga/dsqrt(tnm(i,1))
        un2o1(i,1) = diff*2.01909*un2o0(i,1)*dexp(-847.36/tnm(i,1))
        uch4(i,1) = diff*8.60957E4*ch4(i,1)*pnm(i,1)*rga/dsqrt(tnm(i,1))
        co2fac(i,1) = diff*co2mmr*pnm(i,1)*rga
        alpha1(i) = (1.0-dexp(-1540.0/tnm(i,1)))**3.0/dsqrt(tnm(i,1))
        alpha2(i) = (1.0-dexp(-1360.0/tnm(i,1)))**3.0/dsqrt(tnm(i,1))
        uco211(i,1) = 3.42217E3*co2fac(i,1)*alpha1(i)                   &
                    & *dexp(-1849.7/tnm(i,1))
        uco212(i,1) = 6.02454E3*co2fac(i,1)*alpha1(i)                   &
                    & *dexp(-2782.1/tnm(i,1))
        uco213(i,1) = 5.53143E3*co2fac(i,1)*alpha1(i)                   &
                    & *dexp(-3723.2/tnm(i,1))
        uco221(i,1) = 3.88984E3*co2fac(i,1)*alpha2(i)                   &
                    & *dexp(-1997.6/tnm(i,1))
        uco222(i,1) = 3.67108E3*co2fac(i,1)*alpha2(i)                   &
                    & *dexp(-3843.8/tnm(i,1))
        uco223(i,1) = 6.50642E3*co2fac(i,1)*alpha2(i)                   &
                    & *dexp(-2989.7/tnm(i,1))
        bn2o0(i,1) = diff*19.399*pnm(i,1)**2.0*n2o(i,1)                 &
                   & *1.02346E5*rga/(sslp*tnm(i,1))
        bn2o1(i,1) = bn2o0(i,1)*dexp(-847.36/tnm(i,1))*2.06646E5
        bch4(i,1) = diff*2.94449*ch4(i,1)*pnm(i,1)                      &
                  & **2.0*rga*8.60957E4/(sslp*tnm(i,1))
        uptype(i,1) = diff*qnm(i,1)*pnm(i,1)                            &
                    & **2.0*dexp(1800.0*(1.0/tnm(i,1)-1.0/296.0))       &
                    & *rga/sslp
      end do
      do k = 1 , kx
        do i = 1 , ixm1
          rt(i) = 1./tnm(i,k)
          rsqrt(i) = dsqrt(rt(i))
          pbar(i) = 0.5*(pnm(i,k+1)+pnm(i,k))/sslp
          dpnm(i) = (pnm(i,k+1)-pnm(i,k))*rga
          alpha1(i) = diff*rsqrt(i)*(1.0-dexp(-1540.0/tnm(i,k)))**3.0
          alpha2(i) = diff*rsqrt(i)*(1.0-dexp(-1360.0/tnm(i,k)))**3.0
          ucfc11(i,k+1) = ucfc11(i,k) + 1.8*cfc11(i,k)*dpnm(i)
          ucfc12(i,k+1) = ucfc12(i,k) + 1.8*cfc12(i,k)*dpnm(i)
          un2o0(i,k+1) = un2o0(i,k) + diff*1.02346E5*n2o(i,k)*rsqrt(i)  &
                       & *dpnm(i)
          un2o1(i,k+1) = un2o1(i,k) + diff*2.06646E5*n2o(i,k)*rsqrt(i)  &
                       & *dexp(-847.36/tnm(i,k))*dpnm(i)
          uch4(i,k+1) = uch4(i,k) + diff*8.60957E4*ch4(i,k)*rsqrt(i)    &
                      & *dpnm(i)
          uco211(i,k+1) = uco211(i,k) + 1.15*3.42217E3*alpha1(i)        &
                        & *co2mmr*dexp(-1849.7/tnm(i,k))*dpnm(i)
          uco212(i,k+1) = uco212(i,k) + 1.15*6.02454E3*alpha1(i)        &
                        & *co2mmr*dexp(-2782.1/tnm(i,k))*dpnm(i)
          uco213(i,k+1) = uco213(i,k) + 1.15*5.53143E3*alpha1(i)        &
                        & *co2mmr*dexp(-3723.2/tnm(i,k))*dpnm(i)
          uco221(i,k+1) = uco221(i,k) + 1.15*3.88984E3*alpha2(i)        &
                        & *co2mmr*dexp(-1997.6/tnm(i,k))*dpnm(i)
          uco222(i,k+1) = uco222(i,k) + 1.15*3.67108E3*alpha2(i)        &
                        & *co2mmr*dexp(-3843.8/tnm(i,k))*dpnm(i)
          uco223(i,k+1) = uco223(i,k) + 1.15*6.50642E3*alpha2(i)        &
                        & *co2mmr*dexp(-2989.7/tnm(i,k))*dpnm(i)
          bn2o0(i,k+1) = bn2o0(i,k) + diff*19.399*pbar(i)*rt(i)         &
                       & *1.02346E5*n2o(i,k)*dpnm(i)
          bn2o1(i,k+1) = bn2o1(i,k) + diff*19.399*pbar(i)*rt(i)         &
                       & *2.06646E5*dexp(-847.36/tnm(i,k))*n2o(i,k)     &
                       & *dpnm(i)
          bch4(i,k+1) = bch4(i,k) + diff*2.94449*rt(i)*pbar(i)          &
                      & *8.60957E4*ch4(i,k)*dpnm(i)
          uptype(i,k+1) = uptype(i,k) + diff*qnm(i,k)                   &
                        & *dexp(1800.0*(1.0/tnm(i,k)-1.0/296.0))*pbar(i)&
                        & *dpnm(i)
        end do
      end do

      end subroutine trcpth
