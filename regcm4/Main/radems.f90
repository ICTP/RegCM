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
 
      subroutine radems(s2c,s2t,w,tplnke,plh2o,pnm,plco2,tint,tint4,    &
                      & tlayr,tlayr4,plol,plos,ucfc11,ucfc12,un2o0,     &
                      & un2o1,uch4,uco211,uco212,uco213,uco221,uco222,  &
                      & uco223,uptype,bn2o0,bn2o1,bch4,co2em,co2eml,    &
                      & co2t,h2otr,abplnk1,abplnk2,jslc)
!
!-----------------------------------------------------------------------
!
! Compute emissivity for H2O, CO2, O3
!
! H2O  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
!
!
! CO2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019
!
! O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R. Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! Computes individual emissivities, accounting for band overlap, and
! sums to obtain the total.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_crdcae
      use mod_constants , only : rsslp , rga , dpfco2 , dpfo3
      use mod_radbuf
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! s2c     - H2o continuum path length
! s2t     - Tmp and prs wghted h2o path length
! w       - H2o path length
! tplnke  - Layer planck temperature
! plh2o   - H2o prs wghted path length
! pnm     - Model interface pressure
! plco2   - Prs wghted path of co2
! tint    - Model interface temperatures
! tint4   - Tint to the 4th power
! tlayr   - K-1 model layer temperature
! tlayr4  - Tlayr to the 4th power
! plol    - Pressure wghtd ozone path
! plos    - Ozone path
!
!     Trace gas variables
!
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
! uptype  - p-type continuum path length
!
!     Output arguments
!
! co2em   - Layer co2 normalzd plnck funct drvtv
! co2eml  - Intrfc co2 normalzd plnck func drvtv
! co2t    - Tmp and prs weighted path length
! h2otr   - H2o transmission over o3 band
! emplnk  - emissivity Planck factor
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
! emstrc  - total trace gas emissivity
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(14,ixm1,kxp1) :: abplnk1 , abplnk2
      real(8) , dimension(ixm1,kxp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , h2otr , plco2 , plh2o , plol , plos , pnm , s2c ,   &
           & s2t , tint , tint4 , tlayr , tlayr4 , ucfc11 , ucfc12 ,    &
           & uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,        &
           & uco223 , un2o0 , un2o1 , uptype , w
      real(8) , dimension(ixm1,kx) :: co2eml
      real(8) , dimension(ixm1) :: tplnke
      intent (in) jslc , plco2 , plh2o , plol , plos , s2t , tint4 ,    &
                & tlayr4
      intent (out) co2em , co2eml
      intent (inout) co2t , h2otr
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! iband   - H2o band index
!
!     Local variables for H2O:
!
! h2oems  - H2o emissivity
! tpathe  - Used to compute h2o emissivity
! a       - Eq(2) in table A3a of R&D
! corfac  - Correction factors in table A3b rotation band absorptivity
! dtp     - Path temperature minus 300 K used in
! dtx     - Planck temperature minus 250 K
! dty     - Path temperature minus 250 K
! dtz     - Planck temperature minus 300 K
! emis    - Total emissivity (h2o+co2+o3)
! rsum    - Eq(1) in table A2 of R&D
! term1   - Equation(5) in table A3a of R&D(1986)
! term2   - Delta a(Te) in table A3a of R&D(1986)
! term3   - B(T) function for rotation and vibration-rotation band emissivity
! term4   - Equation(6) in table A3a of R&D(1986)
! term5   - Delta a(Tp) in table A3a of R&D(1986)
! term6   - B(T) function for window region
! term7   - Kl_inf(i) in eq(8) of table A3a of R&D
! term8   - Delta kl_inf(i) in eq(8)
! term9   - B(T) function for 500-800 cm-1 region
! tr1     - Equation(6) in table A2 for 650-800
! tr2     - Equation(6) in table A2 for 500-650
! tr3     - Equation(4) in table A2 for 650-800
! tr4     - Equation(4),table A2 of R&D for 500-650 
! tr7     - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! tr8     - Equation (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! uc      - Y + 0.002U in eq(8) of table A2 of R&D
! pnew    - Effective pressure for h2o linewidth
! trline  - Transmission due to H2O lines in window
! k21     - Exponential coefficient used to calc rot band transmissivity
!           in the 650-800 cm-1 region (tr1)
! k22     - Exponential coefficient used to calc rot band transmissivity
!           in the 500-650 cm-1 region (tr2)
! u       - Pressure weighted H2O path length
! uc1     - H2o continuum pathlength 500-800 cm-1
! r80257  - Conversion factor for h2o pathlength
! a11     - A1 in table A3b for rotation band emiss
! a31     - A3 in table A3b for rotation band emiss
! a21     - First part in numerator of A2 table A3b
! a22     - Second part in numertor of A2 table A3b
! a23     - Denominator of A2 table A3b (rot band)
! t1t4    - Eq(3) in table A3a of R&D
! t2t5    - Eq(4) in table A3a of R&D
! fwk     - Equation(33) in R&D far wing correction
! a41     - Numerator in A2 in Vib-rot (table A3b)
! a51     - Denominator in A2 in Vib-rot(table A3b)
! a61     - A3 factor for Vib-rot band in table A3b
! phi     - Eq(11) in table A3a of R&D
! psi     - Eq(12) in table A3a of R&D
! ubar    - H2o scaled path comment eq(10) table A2
! g1      - Part of eq(10) table A2
! pbar    - H2o scaled pres comment eq(10) table A2
! g3      - Part of eq(10) table A2
! g2      - Part of arguement in eq(10) in table A2
! g4      - Arguement in exp() in eq(10) table A2
! cf812   - Eq(11) in table A2 of R&D
! troco2  - H2o overlap factor for co2 absorption
!
!     Local variables for CO2:
!
! co2ems  - Co2 emissivity
! co2plk  - Used to compute co2 emissivity
! xsum    - Used to calculate path temperature
! t1i     - Co2 hot band temperature factor
! sqti    - Sqrt of temperature
! pi      - Pressure used in co2 mean line width
! et      - Co2 hot band factor
! et2     - Co2 hot band factor
! et4     - Co2 hot band factor
! omet    - Co2 stimulated emission term
! ex      - Part of co2 planck function
! f1co2   - Co2 weak band factor
! f2co2   - Co2 weak band factor
! f3co2   - Co2 weak band factor
! t1co2   - Overlap factor weak bands strong band
! sqwp    - Sqrt of co2 pathlength
! f1sqwp  - Main co2 band factor
! oneme   - Co2 stimulated emission term
! alphat  - Part of the co2 stimulated emiss term
! wco2    - Consts used to define co2 pathlength
! posqt   - Effective pressure for co2 line width
! rbeta7  - Inverse of co2 hot band line width par
! rbeta8  - Inverse of co2 hot band line width par
! rbeta9  - Inverse of co2 hot band line width par
! rbeta13 - Inverse of co2 hot band line width par
! tpath   - Path temp used in co2 band model
! tmp1    - Co2 band factor
! tmp2    - Co2 band factor
! tmp3    - Co2 band factor
! tlayr5  - Temperature factor in co2 Planck func
! rsqti   - Reciprocal of sqrt of temperature
! exm1sq  - Part of co2 Planck function
! u7      - Absorber amount for various co2 band systems
! u8      - Absorber amount for various co2 band systems
! u9      - Absorber amount for various co2 band systems
! u13     - Absorber amount for various co2 band systems
! r250    - Inverse 250K
! r300    - Inverse 300K
! rsslp   - Inverse standard sea-level pressure
!
!     Local variables for O3:
!
! o3ems   - Ozone emissivity
! dbvtt   - Tmp drvtv of planck fctn for tplnke
! te      - Temperature factor
! u1      - Path length factor
! u2      - Path length factor
! phat    - Effecitive path length pressure
! tlocal  - Local planck function temperature
! tcrfac  - Scaled temperature factor
! beta    - Absorption funct factor voigt effect
! realnu  - Absorption function factor
! o3bndi  - Band absorption factor
!
!     Transmission terms for various spectral intervals:
!
! trem4   - H2o   800 - 1000 cm-1
! trem6   - H2o  1000 - 1200 cm-1
! absbnd  - Proportional to co2 band absorptance
! tco2    - co2 overlap factor
! th2o    - h2o overlap factor
! to3     - o3 overlap factor
!
!
! Local variables
!
      real(8) , dimension(ixm1) :: a , co2plk , corfac , dbvtt , dtp ,  &
                                  & dtx , dty , dtz , k21 , k22 , pnew ,&
                                  & rsum , tco2 , th2o , to3 , tpathe , &
                                  & tr1 , tr2 , tr3 , tr4 , tr7 , tr8 , &
                                  & trem4 , trem6 , u , uc , uc1 , xsum
      real(8) :: a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,        &
               & absbnd , alphat , beta , cf812 , et , et2 , et4 , ex , &
               & exm1sq , f1co2 , f1sqwp , f2co2 , f3co2 , fwk , g1 ,   &
               & g2 , g3 , g4 , o3bndi , omet , oneme , pbar , phat ,   &
               & phi , pi , posqt , psi , r250 , r300 , r80257 ,        &
               & rbeta13 , rbeta7 , rbeta8 , rbeta9 , realnu , rsqti ,  &
               & sqti , sqwp , t , t1co2 , t1i , t1t4 , t2t5 ,          &
               & tcrfac , te , tlayr5 , tlocal , tmp1 , tmp2 , tmp3 ,   &
               & tpath , u1 , u13 , u2 , u7 , u8 , u9 , ubar , ux , vx ,&
               & wco2
      real(8) , dimension(ixm1,kxp1) :: co2ems , emstrc , h2oems ,    &
           & o3ems , troco2
      real(8) :: dbvt , fo3
      real(8) , dimension(ixm1,4) :: emis , term1 , term2 , term3 ,    &
                                    & term4 , term5
      real(8) , dimension(14,ixm1) :: emplnk
      integer :: i , iband , k , k1
      real(8) , dimension(ixm1,2) :: term6 , term7 , term8 , term9 ,   &
                                    & trline
!
!---------------------------Statement functions-------------------------
!
!     Statement functions
!     Derivative of planck function at 9.6 micro-meter wavelength, and
!     an absorption function factor:
!
!
      dbvt(t) = (-2.8911366682D-4+(2.3771251896D-6+1.1305188929D-10*t)  &
              & *t)/(1.D0+(-6.1364820707D-3+1.5550319767D-5*t)*t)
!
      fo3(ux,vx) = ux/dsqrt(4.D0+ux*(1.D0+vx))
!
!-----------------------------------------------------------------------
!
!     Initialize
!
      r80257 = 1.D0/8.0257D-04
!
      r250 = 1.D0/250.D0
      r300 = 1.D0/300.D0
!
!     Planck function for co2
!
      do i = 1 , ixm1
        ex = dexp(960.D0/tplnke(i))
        co2plk(i) = 5.D8/((tplnke(i)**4)*(ex-1.))
        co2t(i,1) = tplnke(i)
        xsum(i) = co2t(i,1)*pnm(i,1)
      end do
      k = 1
      do k1 = kxp1 , 2 , -1
        k = k + 1
        do i = 1 , ixm1
          xsum(i) = xsum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
          ex = dexp(960./tlayr(i,k1))
          tlayr5 = tlayr(i,k1)*tlayr4(i,k1)
          co2eml(i,k1-1) = 1.2E11*ex/(tlayr5*(ex-1.)**2)
          co2t(i,k) = xsum(i)/pnm(i,k)
        end do
      end do
!     bndfct = 2.d0*22.18/(dsqrt(196.d0)*300.)
!
!     Initialize planck function derivative for O3
!
      do i = 1 , ixm1
        dbvtt(i) = dbvt(tplnke(i))
      end do
!
!     Calculate trace gas Planck functions
!
      call trcplk(tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
!     Interface loop
!
      do k1 = 1 , kxp1
!
!       H2O emissivity
!
!       emis(i,1)     0 -  800 cm-1   rotation band
!       emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!       emis(i,3)   800 - 1200 cm-1   window
!       emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
!       For the p type continuum
!
        do i = 1 , ixm1
          uc(i) = s2c(i,k1) + 2.E-3*plh2o(i,k1)
          u(i) = plh2o(i,k1)
          pnew(i) = u(i)/w(i,k1)
!
!         Apply scaling factor for 500-800 continuum
!
          uc1(i) = (s2c(i,k1)+1.7E-3*plh2o(i,k1))*(1.+2.*s2c(i,k1))     &
                 & /(1.+15.*s2c(i,k1))
          tpathe(i) = s2t(i,k1)/plh2o(i,k1)
        end do
        do i = 1 , ixm1
          dtx(i) = tplnke(i) - 250.
          dty(i) = tpathe(i) - 250.
          dtz(i) = dtx(i) - 50.
          dtp(i) = dty(i) - 50.
        end do
        do iband = 1 , 3 , 2
          do i = 1 , ixm1
            term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i)     &
                           & *(1.+c1(iband)*dtx(i))
            term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i)     &
                           & *(1.+c2(iband)*dtx(i)*(1.+c3(iband)*dtx(i))&
                           & )
            term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i)     &
                           & *(1.+c4(iband)*dtx(i)*(1.+c5(iband)*dtx(i))&
                           & )
            term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i)     &
                           & *(1.+c6(iband)*dty(i))
            term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i)     &
                           & *(1.+c7(iband)*dty(i))
          end do
        end do
!
!       emis(i,1)     0 -  800 cm-1   rotation band
!
        do i = 1 , ixm1
          a11 = .37 - 3.33E-5*dtz(i) + 3.33E-6*dtz(i)*dtz(i)
          a31 = 1.07 - 1.00E-3*dtp(i) + 1.475E-5*dtp(i)*dtp(i)
          a21 = 1.3870 + 3.80E-3*dtz(i) - 7.8E-6*dtz(i)*dtz(i)
          a22 = 1.0 - 1.21E-3*dtp(i) - 5.33E-6*dtp(i)*dtp(i)
          a23 = 0.9 + 2.62*dsqrt(u(i))
          corfac(i) = a31*(a11+((a21*a22)/a23))
          t1t4 = term1(i,1)*term4(i,1)
          t2t5 = term2(i,1)*term5(i,1)
          a(i) = t1t4 + t2t5/(1.+t2t5*dsqrt(u(i))*corfac(i))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
          emis(i,1) = (1.-rsum(i))*term3(i,1)
!         trem1(i)  = rsum(i)
!
!         emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!
          a41 = 1.75 - 3.96E-3*dtz(i)
          a51 = 1.00 + 1.3*dsqrt(u(i))
          a61 = 1.00 + 1.25E-3*dtp(i) + 6.25E-5*dtp(i)*dtp(i)
          corfac(i) = .3*(1.+(a41)/(a51))*a61
          t1t4 = term1(i,3)*term4(i,3)
          t2t5 = term2(i,3)*term5(i,3)
          a(i) = t1t4 + t2t5/(1.+t2t5*dsqrt(u(i))*corfac(i))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
          emis(i,2) = (1.-rsum(i))*term3(i,3)
!         trem7(i) = rsum(i)
        end do
!
!       Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
        do k = 1 , 2
          do i = 1 , ixm1
            phi = a1(k)*(dty(i)+15.) + a2(k)*(dty(i)+15.)**2
            psi = b1(k)*(dty(i)+15.) + b2(k)*(dty(i)+15.)**2
            phi = dexp(phi)
            psi = dexp(psi)
            ubar = w(i,k1)*phi
            ubar = (ubar*1.66)*r80257
            pbar = pnew(i)*(psi/phi)
            cf812 = cfa1 + ((1.-cfa1)/(1.+ubar*pbar*10.))
            g1 = (realk(k)*pbar)/(2.*st(k))
            g2 = 1. + (ubar*4.0*st(k)*cf812)/pbar
            g3 = dsqrt(g2) - 1.
            g4 = g1*g3
            trline(i,k) = dexp(-g4)
          end do
        end do
        do i = 1 , ixm1
          term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
          term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
          term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
          term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
        end do
!
!       emis(i,3)   800 - 1200 cm-1   window
!
        do i = 1 , ixm1
          term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i)                   &
                     & *(1.+c8*dtx(i)*(1.+c10*dtx(i)                    &
                     & *(1.+c12*dtx(i)*(1.+c14*dtx(i)))))
          trem4(i) = dexp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))        &
                   & *trline(i,2)
          trem6(i) = dexp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))        &
                   & *trline(i,1)
          emis(i,3) = term6(i,1)*(1.-.5*trem4(i)-.5*trem6(i))
!
!         emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
          k21(i) = term7(i,1) + term8(i,1)                              &
                 & /(1.+(c30+c31*(dty(i)-10.)*(dty(i)-10.))*dsqrt(u(i)))
          k22(i) = term7(i,2) + term8(i,2)                              &
                 & /(1.+(c28+c29*(dty(i)-10.))*dsqrt(u(i)))
          term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i)                   &
                     & *(1.+c18*dtx(i)*(1.+c20*dtx(i)                   &
                     & *(1.+c22*dtx(i)*(1.+c24*dtx(i)))))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          tr1(i) = dexp(-(k21(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
          tr2(i) = dexp(-(k22(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
          tr3(i) = dexp(-((coefh(1,1)+coefh(2,1)*dtx(i))*uc1(i)))
          tr4(i) = dexp(-((coefh(1,2)+coefh(2,2)*dtx(i))*uc1(i)))
          tr7(i) = tr1(i)*tr3(i)
          tr8(i) = tr2(i)*tr4(i)
          emis(i,4) = term9(i,1)*.5*(tr1(i)-tr7(i)+tr2(i)-tr8(i))
          h2oems(i,k1) = emis(i,1) + emis(i,2) + emis(i,3) + emis(i,4)
          troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
          th2o(i) = tr8(i)
!         trem2(i)     = troco2(i,k1)
        end do
!
!       CO2 emissivity for 15 micron band system
!
        do i = 1 , ixm1
          t1i = dexp(-480./co2t(i,k1))
          sqti = dsqrt(co2t(i,k1))
          rsqti = 1./sqti
          et = t1i
          et2 = et*et
          et4 = et2*et2
          omet = 1. - 1.5*et2
          f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
          sqwp = dsqrt(plco2(i,k1))
          f1sqwp = f1co2*sqwp
          t1co2 = 1./(1.+245.18*omet*sqwp*rsqti)
          oneme = 1. - et2
          alphat = oneme**3*rsqti
          wco2 = 2.5221*co2vmr*pnm(i,k1)*rga
          u7 = 4.9411E4*alphat*et2*wco2
          u8 = 3.9744E4*alphat*et4*wco2
          u9 = 1.0447E5*alphat*et4*et2*wco2
          u13 = 2.8388E3*alphat*et4*wco2
!
          tpath = co2t(i,k1)
          tlocal = tplnke(i)
          tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
          pi = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
          posqt = pi/(2.*sqti)
          rbeta7 = 1./(5.3288*posqt)
          rbeta8 = 1./(10.6576*posqt)
          rbeta9 = rbeta7
          rbeta13 = rbeta9
          f2co2 = (u7/dsqrt(4.+u7*(1.+rbeta7)))                         &
                & + (u8/dsqrt(4.+u8*(1.+rbeta8)))                       &
                & + (u9/dsqrt(4.+u9*(1.+rbeta9)))
          f3co2 = u13/dsqrt(4.+u13*(1.+rbeta13))
          tmp1 = dlog(1.+f1sqwp)
          tmp2 = dlog(1.+f2co2)
          tmp3 = dlog(1.+f3co2)
          absbnd = (tmp1+2.*t1co2*tmp2+2.*tmp3)*sqti
          tco2(i) = 1.0/(1.0+10.0*(u7/dsqrt(4.+u7*(1.+rbeta7))))
          co2ems(i,k1) = troco2(i,k1)*absbnd*co2plk(i)
          ex = dexp(960./tint(i,k1))
          exm1sq = (ex-1.)**2
          co2em(i,k1) = 1.2E11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
!         trem3(i) = 1. - bndfct*absbnd
        end do
!
!       O3 emissivity
!
        do i = 1 , ixm1
          h2otr(i,k1) = dexp(-12.*s2c(i,k1))
          te = (co2t(i,k1)/293.)**.7
          u1 = 18.29*plos(i,k1)/te
          u2 = .5649*plos(i,k1)/te
          phat = plos(i,k1)/plol(i,k1)
          tlocal = tplnke(i)
          tcrfac = dsqrt(tlocal*r250)*te
          beta = (1./.3205D0)*((1./phat)+(dpfo3*tcrfac))
          realnu = (1./beta)*te
          o3bndi = 74.*te*(tplnke(i)/375.)                              &
                 & *dlog(1.+fo3(u1,realnu)+fo3(u2,realnu))
          o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
          to3(i) = 1.0/(1.+0.1*fo3(u1,realnu)+0.1*fo3(u2,realnu))
!         trem5(i)    = 1.-(o3bndi/(1060-980.))
        end do
!
!       Calculate trace gas emissivities
!
        call trcems(k1,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,bn2o0,bn2o1,  &
                  & uch4,bch4,uco211,uco212,uco213,uco221,uco222,uco223,&
                  & uptype,w,s2c,u,emplnk,th2o,tco2,to3,emstrc)
!
!       Total emissivity:
!
        do i = 1 , ixm1
          emstot(i,k1,jslc) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1) &
                            & + emstrc(i,k1)
        end do
      end do            ! End of interface loop
!
      end subroutine radems
