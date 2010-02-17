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
 
      subroutine radabs(pbr,pnm,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,    &
                      & plco2,plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,&
                      & ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212,   &
                      & uco213,uco221,uco222,uco223,uptype,bn2o0,bn2o1, &
                      & bch4,abplnk1,abplnk2,jslc)

!-----------------------------------------------------------------------
!
! Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
!
! h2o  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
! co2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019.
!            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
!            are also included.
!
! o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R.Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of CFC11
!
! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of CFC12
!
!
! Computes individual absorptivities for non-adjacent layers, accounting
! for band overlap, and sums to obtain the total; then, computes the
! nearest layer contribution.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Briegleb, May 1996
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_crdcae
      use mod_constants , only : rga , sslp , dpfco2 , dpfo3
      use mod_radbuf
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pbr      - Prssr at mid-levels (dynes/cm2)
! pnm      - Prssr at interfaces (dynes/cm2)
! co2em    - Co2 emissivity function
! co2eml   - Co2 emissivity function
! tplnka   - Planck fnctn level temperature
! s2c      - H2o continuum path length
! s2t      - H2o tmp and prs wghted path
! w        - H2o prs wghted path
! h2otr    - H2o trnsmssn fnct for o3 overlap
! plco2    - Co2 prs wghted path length
! plh2o    - H2o prs wfhted path length
! co2t     - Tmp and prs wghted path length
! tint     - Interface temperatures
! tlayr    - K-1 level temperatures
! plol     - Ozone prs wghted path length
! plos     - Ozone path length
! pmln     - Ln(pmidm1)
! piln     - Ln(pintm1)
!
!     Trace gas variables
!
! ucfc11   - CFC11 path length
! ucfc12   - CFC12 path length
! un2o0    - N2O path length
! un2o1    - N2O path length (hot band)
! uch4     - CH4 path length
! uco211   - CO2 9.4 micron band path length
! uco212   - CO2 9.4 micron band path length
! uco213   - CO2 9.4 micron band path length
! uco221   - CO2 10.4 micron band path length
! uco222   - CO2 10.4 micron band path length
! uco223   - CO2 10.4 micron band path length
! uptype   - continuum path length
! bn2o0    - pressure factor for n2o
! bn2o1    - pressure factor for n2o
! bch4     - pressure factor for ch4
! abplnk1   - non-nearest layer Plack factor
! abplnk2   - nearest layer factor
! abstrc    - total trace gas absorptivity
! bplnk     - Planck functions for sub-divided layers
!
!     Output arguments (radbuf)
!
! Nearest layer absorptivities
! Non-adjacent layer absorptivites
! Total emissivity
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(14,ixm1,kxp1) :: abplnk1 , abplnk2
      real(8) , dimension(ixm1,kxp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , h2otr , piln , plco2 , plh2o , plol , plos , pnm ,  &
           & s2c , s2t , tint , tlayr , tplnka , ucfc11 , ucfc12 ,      &
           & uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,        &
           & uco223 , un2o0 , un2o1 , uptype , w
      real(8) , dimension(ixm1,kx) :: co2eml , pbr , pmln
      intent (in) abplnk2 , co2em , co2eml , co2t , h2otr , jslc , pbr ,&
                & piln , plco2 , plh2o , plol , plos , pmln , s2t ,     &
                & tint , tlayr , tplnka , w
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! k1       - Level index
! k2       - Level index
! kn       - Nearest level index
! iband    - Band  index
! pnew     - Effective pressure for H2O vapor linewidth
! trline   - Transmission due to H2O lines in window
! u        - Pressure weighted H2O path length
! tbar     - Mean layer temperature
! emm      - Mean co2 emissivity
! o3emm    - Mean o3 emissivity
! o3bndi   - Ozone band parameter
! temh2o   - Mean layer temperature equivalent to tbar
! k21      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 650-800 cm-1 region (tr1)
! k22      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 500-650 cm-1 region (tr2)
! uc1      - H2o continuum pathlength in 500-800 cm-1
! to3h2o   - H2o trnsmsn for overlap with o3
! pi       - For co2 absorptivity computation
! sqti     - Used to store sqrt of mean temperature
! et       - Co2 hot band factor
! et2      - Co2 hot band factor squared
! et4      - Co2 hot band factor to fourth power
! omet     - Co2 stimulated emission term
! f1co2    - Co2 central band factor
! f2co2    - Co2 weak band factor
! f3co2    - Co2 weak band factor
! t1co2    - Overlap factr weak bands on strong band
! sqwp     - Sqrt of co2 pathlength
! f1sqwp   - Main co2 band factor
! oneme    - Co2 stimulated emission term
! alphat   - Part of the co2 stimulated emission term
! wco2     - Constants used to define co2 pathlength
! posqt    - Effective pressure for co2 line width
! u7       - Co2 hot band path length
! u8       - Co2 hot band path length
! u9       - Co2 hot band path length
! u13      - Co2 hot band path length
! rbeta7   - Inverse of co2 hot band line width par
! rbeta8   - Inverse of co2 hot band line width par
! rbeta9   - Inverse of co2 hot band line width par
! rbeta13  - Inverse of co2 hot band line width par
! tpatha   - For absorptivity computation
! a        - Eq(2) in table A3a of R&D
! abso     - Absorptivity for various gases/bands
! dtp      - Path temp minus 300 K used in h2o rotation band absorptivity
! dtx      - Planck temperature minus 250 K
! dty      - Path temperature minus 250 K
! dtz      - Planck temperature minus 300 K
! term1    - Equation(5) in table A3a of R&D(1986)
! term2    - Delta a(Te) in table A3a of R&D(1986)
! term3    - DB/dT function for rotation and
! term4    - Equation(6) in table A3a of R&D(1986)
! term5    - Delta a(Tp) in table A3a of R&D(1986)
! term6    - DB/dT function for window region
! term7    - Kl_inf(i) in eq(8) of table A3a of R&D
! term8    - Delta kl_inf(i) in eq(8)
! term9    - DB/dT function for 500-800 cm-1 region
! tr1      - Eqn(6) in table A2 of R&D for 650-800
! tr10     - Eqn (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! tr2      - Eqn(6) in table A2 of R&D for 500-650
! tr5      - Eqn(4) in table A2 of R&D for 650-800
! tr6      - Eqn(4) in table A2 of R&D for 500-650
! tr9      - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! uc       - Y + 0.002U in eq(8) of table A2 of R&D
! sqrtu    - Sqrt of pressure weighted h20 pathlength
! fwk      - Equation(33) in R&D far wing correction
! fwku     - GU term in eqs(1) and (6) in table A2
! r2st     - 1/(2*beta) in eq(10) in table A2
! dtyp15   - DeltaTp in eqs(11) & (12) in table A3a
! dtyp15sq - (DeltaTp)^2 in eqs(11) & (12) table A3a
! to3co2   - P weighted temp in ozone band model
! dpnm     - Pressure difference between two levels
! pnmsq    - Pressure squared
! dw       - Amount of h2o between two levels
! uinpl    - Nearest layer subdivision factor
! winpl    - Nearest layer subdivision factor
! zinpl    - Nearest layer subdivision factor
! pinpl    - Nearest layer subdivision factor
! dplh2o   - Difference in press weighted h2o amount
! r80257   - Conversion factor for h2o pathlength
! r293     - 1/293
! r250     - 1/250
! r3205    - Line width factor for o3 (see R&Di)
! r300     - 1/300
! rsslp    - Reciprocal of sea level pressure
! r2sslp   - 1/2 of rsslp
! ds2c     - Y in eq(7) in table A2 of R&D
! a11      - A1 in table A3b for rotation band absorptivity
! a31      - A3 in table A3b for rotation band absorptivity
! a21      - First part in numerator of A2 in table A3b
! a22      - Second part in numerator of A2 in table A3b
! a23      - Denominator of A2 in table A3b (rotation band)
! t1t4     - Eq(3) in table A3a of R&D
! t2t5     - Eq(4) in table A3a of R&D
! rsum     - Eq(1) in table A2 of R&D
! a41      - Numerator in A2 in Vib-rot abstivity(table A3b)
! a51      - Denominator in A2 in Vib-rot (table A3b)
! a61      - A3 factor for Vib-rot band in table A3b
! phi      - Eq(11) in table A3a of R&D
! psi      - Eq(12) in table A3a of R&D
! cf812    - Eq(11) in table A2 of R&D
! ubar     - H2o scaled path see comment for eq(10) table A2
! pbar     - H2o scaled pres see comment for eq(10) table A2
! g4       - Arguement in exp() in eq(10) table A2
! dplos    - Ozone pathlength eq(A2) in R&Di
! dplol    - Presure weighted ozone pathlength
! beta     - Local interface temperature (includes Voigt line correction factor)
! rphat    - Effective pressure for ozone beta
! tcrfac   - Ozone temperature factor table 1 R&Di
! tmp1     - Ozone band factor see eq(A1) in R&Di
! u1       - Effective ozone pathlength eq(A2) in R&Di
! realnu   - 1/beta factor in ozone band model eq(A1)
! tmp2     - Ozone band factor see eq(A1) in R&Di
! u2       - Effective ozone pathlength eq(A2) in R&Di
! rsqti    - Reciprocal of sqrt of path temperature
! tpath    - Path temperature used in co2 band model
! tmp3     - Weak band factor see K&B
! rdpnmsq  - Reciprocal of difference in press^2
! rdpnm    - Reciprocal of difference in press
! p1       - Mean pressure factor
! p2       - Mean pressure factor
! dtym10   - T - 260 used in eq(9) and (10) table A3a
! dplco2   - Co2 pathlength
! corfac   - Correction factors in table A3b
! g2       - Part of arguement in eq(10) in table A2
! te       - A_0 T factor in ozone model table 1 of R&Di
! denom    - Denominator in eq(8) of table A3a of R&D
! trab2    - Transmission terms for H2o  500 -  800 cm-1
! trab4    - Transmission terms for H2o  800 - 1000 cm-1
! trab6    - Transmission terms for H2o 1000 - 1200 cm-1
! absbnd   - Proportional to co2 band absorptance
! dbvtit   - Intrfc drvtv plnck fnctn for o3
! dbvtly   - Level drvtv plnck fnctn for o3
! dbvt     - Planck fnctn tmp derivative for o3
!
!
! Local variables
!
      real(8) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,    &
               & absbnd , alphat , beta , cf812 , corfac , denom ,      &
               & dplco2 , dplol , dplos , ds2c , dtym10 , et , et2 ,    &
               & et4 , f1co2 , g2 , g4 , k21 , k22 , o3bndi , omet ,    &
               & oneme , p1 , p2 , pbar , phi , pi , posqt , psi ,      &
               & r250 , r293 , r2sslp , r300 , r3205 , r80257 ,         &
               & rbeta13 , rbeta8 , rbeta9 , rdpnm , rdpnmsq , realnu , &
               & rphat , rsqti , rsslp , rsum , sqwp , t , t1t4 , t2t5 ,&
               & tcrfac , te , tlocal , tmp1 , tmp2 , tmp3 , tpath ,    &
               & tr1 , tr2 , tr5 , tr6 , u1 , u13 , u2 , u8 , u9 ,      &
               & ubar , wco2
      real(8) , dimension(ixm1,6) :: abso
      real(8) , dimension(ixm1) :: abstrc , dplh2o , dpnm , dtp , dtx ,&
                                  & dty , dtyp15 , dtyp15sq , dtz , dw ,&
                                  & f1sqwp , f2co2 , f3co2 , fwk ,      &
                                  & fwku , pnew , rbeta7 , sqrtu ,      &
                                  & sqti , t1co2 , tco2 , th2o , to3 ,  &
                                  & to3co2 , to3h2o , tpatha , tr10 ,   &
                                  & tr9 , trab2 , trab4 , trab6 , u ,   &
                                  & u7 , uc , uc1
      real(8) , dimension(14,ixm1,4) :: bplnk
      real(8) :: dbvt
      real(8) , dimension(ixm1,kxp1) :: dbvtit , pnmsq , term6 , term9
      real(8) , dimension(ixm1,kx) :: dbvtly
      real(8) , dimension(ixm1,4) :: emm , o3emm , pinpl , tbar ,      &
                                    & temh2o , term1 , term2 , term3 ,  &
                                    & term4 , term5 , uinpl , winpl ,   &
                                    & zinpl
      integer :: i , iband , k , k1 , k2 , kn , wvl
      real(8) , dimension(2) :: r2st
      real(8) , dimension(ixm1,2) :: term7 , term8 , trline
!
!--------------------------Statement function---------------------------
!
      dbvt(t) = (-2.8911366682E-4+(2.3771251896E-6+1.1305188929E-10*t)  &
              & *t)/(1.0+(-6.1364820707E-3+1.5550319767E-5*t)*t)
!
!-----------------------------------------------------------------------
!
!     Initialize
!
      do k = 1 , kx
        do i = 1 , ixm1
          dbvtly(i,k) = dbvt(tlayr(i,k+1))
          dbvtit(i,k) = dbvt(tint(i,k))
        end do
      end do
      do i = 1 , ixm1
        dbvtit(i,kxp1) = dbvt(tint(i,kx + 1))
      end do
!
      r80257 = 1./8.0257D-04
      r293 = 1./293.D0
      r250 = 1./250.D0
      r3205 = 1./.3205D0
      r300 = 1./300.D0
      rsslp = 1./sslp
      r2sslp = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
!     bndfct  = 2.0*22.18d0/(dsqrt(196.d0)*300.)
!
!     Non-adjacent layer absorptivity:
!
!     abso(i,1)     0 -  800 cm-1   h2o rotation band
!     abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!     abso(i,3)   800 - 1200 cm-1   h2o window
!     abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!     abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!     abso(i,6)   co2 15  micrometer band system
!
      do k = 1 , kxp1
        do i = 1 , ixm1
          pnmsq(i,k) = pnm(i,k)**2
          dtx(i) = tplnka(i,k) - 250.
          term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)                   &
                     & *(1.+c9*dtx(i)*(1.+c11*dtx(i)                    &
                     & *(1.+c13*dtx(i)*(1.+c15*dtx(i)))))
          term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)                   &
                     & *(1.+c19*dtx(i)*(1.+c21*dtx(i)                   &
                     & *(1.+c23*dtx(i)*(1.+c25*dtx(i)))))
        end do
      end do
!
!     Non-nearest layer level loops
!
      do k1 = kxp1 , 1 , -1
        do k2 = kxp1 , 1 , -1
          if ( k1.ne.k2 ) then
            do i = 1 , ixm1
              dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
              u(i) = dabs(dplh2o(i))
              sqrtu(i) = dsqrt(u(i))
              ds2c = dabs(s2c(i,k1)-s2c(i,k2))
              dw(i) = dabs(w(i,k1)-w(i,k2))
              uc1(i) = (ds2c+1.7E-3*u(i))*(1.+2.*ds2c)/(1.+15.*ds2c)
              uc(i) = ds2c + 2.E-3*u(i)
            end do
            do i = 1 , ixm1
              pnew(i) = u(i)/dw(i)
              tpatha(i) = (s2t(i,k1)-s2t(i,k2))/dplh2o(i)
              dtx(i) = tplnka(i,k2) - 250.
              dty(i) = tpatha(i) - 250.
              dtyp15(i) = dty(i) + 15.
              dtyp15sq(i) = dtyp15(i)**2
              dtz(i) = dtx(i) - 50.
              dtp(i) = dty(i) - 50.
            end do
            do iband = 2 , 4 , 2
              do i = 1 , ixm1
                term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i) &
                               & *(1.+c1(iband)*dtx(i))
                term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i) &
                               & *(1.+c2(iband)*dtx(i)                  &
                               & *(1.+c3(iband)*dtx(i)))
                term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i) &
                               & *(1.+c4(iband)*dtx(i)                  &
                               & *(1.+c5(iband)*dtx(i)))
                term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i) &
                               & *(1.+c6(iband)*dty(i))
                term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i) &
                               & *(1.+c7(iband)*dty(i))
              end do
            end do
!
!           abso(i,1)     0 -  800 cm-1   h2o rotation band
!
            do i = 1 , ixm1
              a11 = 0.44 + 3.380E-4*dtz(i) - 1.520E-6*dtz(i)*dtz(i)
              a31 = 1.05 - 6.000E-3*dtp(i) + 3.000E-6*dtp(i)*dtp(i)
              a21 = 1.00 + 1.717E-3*dtz(i) - 1.133E-5*dtz(i)*dtz(i)
              a22 = 1.00 + 4.443E-3*dtp(i) + 2.750E-5*dtp(i)*dtp(i)
              a23 = 1.00 + 3.600*sqrtu(i)
              corfac = a31*(a11+((2.*a21*a22)/a23))
              t1t4 = term1(i,2)*term4(i,2)
              t2t5 = term2(i,2)*term5(i,2)
              a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
              fwk(i) = fwcoef + fwc1/(1.+fwc2*u(i))
              fwku(i) = fwk(i)*u(i)
              rsum = dexp(-a*(sqrtu(i)+fwku(i)))
              abso(i,1) = (1.-rsum)*term3(i,2)
!             trab1(i)  = rsum
            end do
!
!           abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
            do i = 1 , ixm1
              a41 = 1.75 - 3.960E-03*dtz(i)
              a51 = 1.00 + 1.3*sqrtu(i)
              a61 = 1.00 + 1.250E-03*dtp(i) + 6.250E-05*dtp(i)*dtp(i)
              corfac = .29*(1.+a41/a51)*a61
              t1t4 = term1(i,4)*term4(i,4)
              t2t5 = term2(i,4)*term5(i,4)
              a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
              rsum = dexp(-a*(sqrtu(i)+fwku(i)))
              abso(i,2) = (1.-rsum)*term3(i,4)
!             trab7(i)  = rsum
            end do
!
!           Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
            do k = 1 , 2
              do i = 1 , ixm1
                phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
                psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
                ubar = dw(i)*phi*1.66*r80257
                pbar = pnew(i)*(psi/phi)
                cf812 = cfa1 + (1.-cfa1)/(1.+ubar*pbar*10.)
                g2 = 1. + ubar*4.0*st(k)*cf812/pbar
                g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-1.)
                trline(i,k) = dexp(-g4)
              end do
            end do
            do i = 1 , ixm1
              term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)               &
                         & *(1.+c16*dty(i))
              term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)               &
                         & *(1.+c17*dty(i))
              term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)               &
                         & *(1.+c26*dty(i))
              term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)               &
                         & *(1.+c27*dty(i))
            end do
!
!           abso(i,3)   800 - 1200 cm-1   h2o window
!           abso(i,4)   500 -  800 cm-1   h2o rotation band overlap
!           with co2
            do i = 1 , ixm1
              k21 = term7(i,1) + term8(i,1)                             &
                  & /(1.+(c30+c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
              k22 = term7(i,2) + term8(i,2)                             &
                  & /(1.+(c28+c29*(dty(i)-10.))*sqrtu(i))
              tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
              tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
              tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
              tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
              tr9(i) = tr1*tr5
              tr10(i) = tr2*tr6
              th2o(i) = tr10(i)
              trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
              trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
              trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
              abso(i,3) = term6(i,k2)                                   &
                        & *(1.-.5*trab4(i)*trline(i,2)-.5*trab6(i)      &
                        & *trline(i,1))
              abso(i,4) = term9(i,k2)*.5*(tr1-tr9(i)+tr2-tr10(i))
            end do
            if ( k2.lt.k1 ) then
              do i = 1 , ixm1
                to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
              end do
            else
              do i = 1 , ixm1
                to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
              end do
            end if
!
!           abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
            do i = 1 , ixm1
              dpnm(i) = pnm(i,k1) - pnm(i,k2)
              to3co2(i) = (pnm(i,k1)*co2t(i,k1)-pnm(i,k2)*co2t(i,k2))   &
                        & /dpnm(i)
              te = (to3co2(i)*r293)**.7
              dplos = plos(i,k1) - plos(i,k2)
              dplol = plol(i,k1) - plol(i,k2)
              u1 = 18.29*dabs(dplos)/te
              u2 = .5649*dabs(dplos)/te
              rphat = dplol/dplos
              tlocal = tint(i,k2)
              tcrfac = dsqrt(tlocal*r250)*te
              beta = r3205*(rphat+dpfo3*tcrfac)
              realnu = te/beta
              tmp1 = u1/dsqrt(4.+u1*(1.+realnu))
              tmp2 = u2/dsqrt(4.+u2*(1.+realnu))
              o3bndi = 74.*te*dlog(1.+tmp1+tmp2)
              abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
              to3(i) = 1.0/(1.+0.1*tmp1+0.1*tmp2)
!             trab5(i)  = 1.-(o3bndi/(1060-980.))
            end do
!
!           abso(i,6)      co2 15  micrometer band system
!
            do i = 1 , ixm1
              sqwp = dsqrt(dabs(plco2(i,k1)-plco2(i,k2)))
              et = dexp(-480./to3co2(i))
              sqti(i) = dsqrt(to3co2(i))
              rsqti = 1./sqti(i)
              et2 = et*et
              et4 = et2*et2
              omet = 1. - 1.5*et2
              f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
              f1sqwp(i) = f1co2*sqwp
              t1co2(i) = 1./(1.+(245.18*omet*sqwp*rsqti))
              oneme = 1. - et2
              alphat = oneme**3*rsqti
              pi = dabs(dpnm(i))
              wco2 = 2.5221*co2vmr*pi*rga
              u7(i) = 4.9411E4*alphat*et2*wco2
              u8 = 3.9744E4*alphat*et4*wco2
              u9 = 1.0447E5*alphat*et4*et2*wco2
              u13 = 2.8388E3*alphat*et4*wco2
              tpath = to3co2(i)
              tlocal = tint(i,k2)
              tcrfac = dsqrt(tlocal*r250*tpath*r300)
              posqt = ((pnm(i,k2)+pnm(i,k1))*r2sslp+dpfco2*tcrfac)*rsqti
              rbeta7(i) = 1./(5.3228*posqt)
              rbeta8 = 1./(10.6576*posqt)
              rbeta9 = rbeta7(i)
              rbeta13 = rbeta9
              f2co2(i) = (u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i))))         &
                       & + (u8/dsqrt(4.+u8*(1.+rbeta8)))                &
                       & + (u9/dsqrt(4.+u9*(1.+rbeta9)))
              f3co2(i) = u13/dsqrt(4.+u13*(1.+rbeta13))
            end do
            if ( k2.ge.k1 ) then
              do i = 1 , ixm1
                sqti(i) = dsqrt(tlayr(i,k2))
              end do
            end if
!
            do i = 1 , ixm1
              tmp1 = dlog(1.+f1sqwp(i))
              tmp2 = dlog(1.+f2co2(i))
              tmp3 = dlog(1.+f3co2(i))
              absbnd = (tmp1+2.*t1co2(i)*tmp2+2.*tmp3)*sqti(i)
              abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
              tco2(i) = 1./(1.0+10.0*(u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)&
                      & ))))
!             trab3(i)  = 1. - bndfct*absbnd
            end do
!
!           Calculate absorptivity due to trace gases
!
            call trcab(k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,     &
                     & uco212,uco213,uco221,uco222,uco223,bn2o0,bn2o1,  &
                     & bch4,to3co2,pnm,dw,pnew,s2c,uptype,u,abplnk1,    &
                     & tco2,th2o,to3,abstrc)
!
!           Sum total absorptivity
!
            do i = 1 , ixm1
              abstot(i,k1,k2,jslc) = abso(i,1) + abso(i,2) + abso(i,3)  &
                                   & + abso(i,4) + abso(i,5) + abso(i,6)&
                                   & + abstrc(i)
            end do
          end if
        end do
      end do          ! End of non-nearest layer level loops
!
!     Non-adjacent layer absorptivity:
!
!     abso(i,1)     0 -  800 cm-1   h2o rotation band
!     abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!     abso(i,3)   800 - 1200 cm-1   h2o window
!     abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!     abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!     abso(i,6)   co2 15  micrometer band system
!
!     Nearest layer level loop
!
      do k2 = kx , 1 , -1
        do i = 1 , ixm1
          tbar(i,1) = 0.5*(tint(i,k2+1)+tlayr(i,k2+1))
          emm(i,1) = 0.5*(co2em(i,k2+1)+co2eml(i,k2))
          tbar(i,2) = 0.5*(tlayr(i,k2+1)+tint(i,k2))
          emm(i,2) = 0.5*(co2em(i,k2)+co2eml(i,k2))
          tbar(i,3) = 0.5*(tbar(i,2)+tbar(i,1))
          emm(i,3) = emm(i,1)
          tbar(i,4) = tbar(i,3)
          emm(i,4) = emm(i,2)
          o3emm(i,1) = 0.5*(dbvtit(i,k2+1)+dbvtly(i,k2))
          o3emm(i,2) = 0.5*(dbvtit(i,k2)+dbvtly(i,k2))
          o3emm(i,3) = o3emm(i,1)
          o3emm(i,4) = o3emm(i,2)
          temh2o(i,1) = tbar(i,1)
          temh2o(i,2) = tbar(i,2)
          temh2o(i,3) = tbar(i,1)
          temh2o(i,4) = tbar(i,2)
          dpnm(i) = pnm(i,k2+1) - pnm(i,k2)
        end do
!----------------------------------------------------------
!       Weighted Planck functions for trace gases
!
        do wvl = 1 , 14
          do i = 1 , ixm1
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1)+abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2)+abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
          end do
        end do
!---------------------------------------------------------
        do i = 1 , ixm1
          rdpnmsq = 1./(pnmsq(i,k2+1)-pnmsq(i,k2))
          rdpnm = 1./dpnm(i)
          p1 = .5*(pbr(i,k2)+pnm(i,k2+1))
          p2 = .5*(pbr(i,k2)+pnm(i,k2))
          uinpl(i,1) = (pnmsq(i,k2+1)-p1**2)*rdpnmsq
          uinpl(i,2) = -(pnmsq(i,k2)-p2**2)*rdpnmsq
          uinpl(i,3) = -(pnmsq(i,k2)-p1**2)*rdpnmsq
          uinpl(i,4) = (pnmsq(i,k2+1)-p2**2)*rdpnmsq
          winpl(i,1) = (.5*(pnm(i,k2+1)-pbr(i,k2)))*rdpnm
          winpl(i,2) = (.5*(-pnm(i,k2)+pbr(i,k2)))*rdpnm
          winpl(i,3) = (.5*(pnm(i,k2+1)+pbr(i,k2))-pnm(i,k2))*rdpnm
          winpl(i,4) = (.5*(-pnm(i,k2)-pbr(i,k2))+pnm(i,k2+1))*rdpnm
          tmp1 = 1./(piln(i,k2+1)-piln(i,k2))
          tmp2 = piln(i,k2+1) - pmln(i,k2)
          tmp3 = piln(i,k2) - pmln(i,k2)
          zinpl(i,1) = (.5*tmp2)*tmp1
          zinpl(i,2) = (-.5*tmp3)*tmp1
          zinpl(i,3) = (.5*tmp2-tmp3)*tmp1
          zinpl(i,4) = (tmp2-.5*tmp3)*tmp1
          pinpl(i,1) = 0.5*(p1+pnm(i,k2+1))
          pinpl(i,2) = 0.5*(p2+pnm(i,k2))
          pinpl(i,3) = 0.5*(p1+pnm(i,k2))
          pinpl(i,4) = 0.5*(p2+pnm(i,k2+1))
        end do
        do kn = 1 , 4
          do i = 1 , ixm1
            u(i) = uinpl(i,kn)*dabs(plh2o(i,k2)-plh2o(i,k2+1))
            sqrtu(i) = dsqrt(u(i))
            dw(i) = dabs(w(i,k2)-w(i,k2+1))
            pnew(i) = u(i)/(winpl(i,kn)*dw(i))
            ds2c = dabs(s2c(i,k2)-s2c(i,k2+1))
            uc1(i) = uinpl(i,kn)*ds2c
            uc1(i) = (uc1(i)+1.7E-3*u(i))*(1.+2.*uc1(i))/(1.+15.*uc1(i))
            uc(i) = uinpl(i,kn)*ds2c + 2.E-3*u(i)
          end do
          do i = 1 , ixm1
            dtx(i) = temh2o(i,kn) - 250.
            dty(i) = tbar(i,kn) - 250.
            dtyp15(i) = dty(i) + 15.
            dtyp15sq(i) = dtyp15(i)**2
            dtz(i) = dtx(i) - 50.
            dtp(i) = dty(i) - 50.
          end do
          do iband = 2 , 4 , 2
            do i = 1 , ixm1
              term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i)   &
                             & *(1.+c1(iband)*dtx(i))
              term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i)   &
                             & *(1.+c2(iband)*dtx(i)                    &
                             & *(1.+c3(iband)*dtx(i)))
              term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i)   &
                             & *(1.+c4(iband)*dtx(i)                    &
                             & *(1.+c5(iband)*dtx(i)))
              term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i)   &
                             & *(1.+c6(iband)*dty(i))
              term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i)   &
                             & *(1.+c7(iband)*dty(i))
            end do
          end do
!
!         abso(i,1)     0 -  800 cm-1   h2o rotation band
!
          do i = 1 , ixm1
            a11 = 0.44 + 3.380E-4*dtz(i) - 1.520E-6*dtz(i)*dtz(i)
            a31 = 1.05 - 6.000E-3*dtp(i) + 3.000E-6*dtp(i)*dtp(i)
            a21 = 1.00 + 1.717E-3*dtz(i) - 1.133E-5*dtz(i)*dtz(i)
            a22 = 1.00 + 4.443E-3*dtp(i) + 2.750E-5*dtp(i)*dtp(i)
            a23 = 1.00 + 3.600*sqrtu(i)
            corfac = a31*(a11+((2.*a21*a22)/a23))
            t1t4 = term1(i,2)*term4(i,2)
            t2t5 = term2(i,2)*term5(i,2)
            a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
            fwk(i) = fwcoef + fwc1/(1.+fwc2*u(i))
            fwku(i) = fwk(i)*u(i)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,1) = (1.-rsum)*term3(i,2)
!           trab1(i) = rsum
          end do
!
!         abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
          do i = 1 , ixm1
            a41 = 1.75 - 3.960E-03*dtz(i)
            a51 = 1.00 + 1.3*sqrtu(i)
            a61 = 1.00 + 1.250E-03*dtp(i) + 6.250E-05*dtp(i)*dtp(i)
            corfac = .29*(1.+a41/a51)*a61
            t1t4 = term1(i,4)*term4(i,4)
            t2t5 = term2(i,4)*term5(i,4)
            a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,2) = (1.-rsum)*term3(i,4)
!           trab7(i) = rsum
          end do
!
!         Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
          do k = 1 , 2
            do i = 1 , ixm1
              phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
              psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
              ubar = dw(i)*phi*winpl(i,kn)*1.66*r80257
              pbar = pnew(i)*(psi/phi)
              cf812 = cfa1 + (1.-cfa1)/(1.+ubar*pbar*10.)
              g2 = 1. + ubar*4.0*st(k)*cf812/pbar
              g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-1.)
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
!         abso(i,3)   800 - 1200 cm-1   h2o window
!         abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with
!         co2
          do i = 1 , ixm1
            dtym10 = dty(i) - 10.
            denom = 1. + (c30+c31*dtym10*dtym10)*sqrtu(i)
            k21 = term7(i,1) + term8(i,1)/denom
            denom = 1. + (c28+c29*dtym10)*sqrtu(i)
            k22 = term7(i,2) + term8(i,2)/denom
            term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)                 &
                       & *(1.+c19*dtx(i)*(1.+c21*dtx(i)                 &
                       & *(1.+c23*dtx(i)*(1.+c25*dtx(i)))))
            tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
            tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
            tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
            tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i) = tr1*tr5
            tr10(i) = tr2*tr6
            trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
            th2o(i) = tr10(i)
            trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
            trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
            term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)                 &
                       & *(1.+c9*dtx(i)*(1.+c11*dtx(i)                  &
                       & *(1.+c13*dtx(i)*(1.+c15*dtx(i)))))
            abso(i,3) = term6(i,2)                                      &
                      & *(1.-.5*trab4(i)*trline(i,2)-.5*trab6(i)        &
                      & *trline(i,1))
            abso(i,4) = term9(i,2)*.5*(tr1-tr9(i)+tr2-tr10(i))
          end do
!
!         abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
!
          do i = 1 , ixm1
            te = (tbar(i,kn)*r293)**.7
            dplos = dabs(plos(i,k2+1)-plos(i,k2))
            u1 = zinpl(i,kn)*18.29*dplos/te
            u2 = zinpl(i,kn)*.5649*dplos/te
            tlocal = tbar(i,kn)
            tcrfac = dsqrt(tlocal*r250)*te
            beta = r3205*(pinpl(i,kn)*rsslp+dpfo3*tcrfac)
            realnu = te/beta
            tmp1 = u1/dsqrt(4.+u1*(1.+realnu))
            tmp2 = u2/dsqrt(4.+u2*(1.+realnu))
            o3bndi = 74.*te*dlog(1.+tmp1+tmp2)
            abso(i,5) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
            to3(i) = 1.0/(1.+0.1*tmp1+0.1*tmp2)
!           trab5(i) = 1.-(o3bndi/(1060-980.))
          end do
!
!         abso(i,6)   co2 15  micrometer band system
!
          do i = 1 , ixm1
            dplco2 = plco2(i,k2+1) - plco2(i,k2)
            sqwp = dsqrt(uinpl(i,kn)*dplco2)
            et = dexp(-480./tbar(i,kn))
            sqti(i) = dsqrt(tbar(i,kn))
            rsqti = 1./sqti(i)
            et2 = et*et
            et4 = et2*et2
            omet = (1.-1.5*et2)
            f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
            f1sqwp(i) = f1co2*sqwp
            t1co2(i) = 1./(1.+(245.18*omet*sqwp*rsqti))
            oneme = 1. - et2
            alphat = oneme**3*rsqti
            pi = dabs(dpnm(i))*winpl(i,kn)
            wco2 = 2.5221*co2vmr*pi*rga
            u7(i) = 4.9411E4*alphat*et2*wco2
            u8 = 3.9744E4*alphat*et4*wco2
            u9 = 1.0447E5*alphat*et4*et2*wco2
            u13 = 2.8388E3*alphat*et4*wco2
            tpath = tbar(i,kn)
            tlocal = tbar(i,kn)
            tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
            posqt = (pinpl(i,kn)*rsslp+dpfco2*tcrfac)*rsqti
            rbeta7(i) = 1./(5.3228*posqt)
            rbeta8 = 1./(10.6576*posqt)
            rbeta9 = rbeta7(i)
            rbeta13 = rbeta9
            f2co2(i) = u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)))             &
                     & + u8/dsqrt(4.+u8*(1.+rbeta8))                    &
                     & + u9/dsqrt(4.+u9*(1.+rbeta9))
            f3co2(i) = u13/dsqrt(4.+u13*(1.+rbeta13))
            tmp1 = dlog(1.+f1sqwp(i))
            tmp2 = dlog(1.+f2co2(i))
            tmp3 = dlog(1.+f3co2(i))
            absbnd = (tmp1+2.*t1co2(i)*tmp2+2.*tmp3)*sqti(i)
            abso(i,6) = trab2(i)*emm(i,kn)*absbnd
            tco2(i) = 1.0/(1.0+10.0*u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)))&
                    & )
!           trab3(i) = 1. - bndfct*absbnd
          end do
!
!         Calculate trace gas absorptivity for nearest layer
!
          call trcabn(k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,      &
                    & uco212,uco213,uco221,uco222,uco223,tbar,bplnk,    &
                    & winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,u,pnew,   &
                    & abstrc,uinpl)
!
!         Total next layer absorptivity:
!
          do i = 1 , ixm1
            absnxt(i,k2,kn,jslc) = abso(i,1) + abso(i,2) + abso(i,3)    &
                                 & + abso(i,4) + abso(i,5) + abso(i,6)  &
                                 & + abstrc(i)
          end do
        end do
      end do                    !  end of nearest layer level loop
!
      end subroutine radabs
