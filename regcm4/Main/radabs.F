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

      subroutine radabs(pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  ,
     $                  s2c    ,s2t     ,w        ,h2otr   ,plco2   ,
     $                  plh2o  ,co2t    ,tint     ,tlayr   ,plol    ,
     $                  plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  ,
     $                  un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  ,
     $                  uco213 ,uco221  ,uco222   ,uco223  ,uptype  ,
     $                  bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 ,
     $                  jslc  )
C-----------------------------------------------------------------------
C
C Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
C
C h2o  ....  Uses nonisothermal emissivity for water vapor from
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
C
C co2  ....  Uses absorptance parameterization of the 15 micro-meter
C            (500 - 800 cm-1) band system of Carbon Dioxide, from
C            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
C            of the Absorptance Due to the 15 micro-meter Band System
C            of Carbon Dioxide Jouranl of Geophysical Research,
C            vol. 96., D5, pp 9013-9019.
C            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
C            are also included.
C
C o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
C            band system of ozone, from Ramanathan, V. and R.Dickinson,
C            1979: The Role of stratospheric ozone in the zonal and
C            seasonal radiative energy balance of the earth-troposphere
C            system. Journal of the Atmospheric Sciences, Vol. 36,
C            pp 1084-1104
C
C ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
C
C n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
C            bands of nitrous oxide
C
C cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
C            micron bands of CFC11
C
C cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
C            micron bands of CFC12
C
C
C Computes individual absorptivities for non-adjacent layers, accounting
C for band overlap, and sums to obtain the total; then, computes the
C nearest layer contribution.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C
C-----------------------------------------------------------------------
c
c $Id: radabs.F,v 1.1.16.1 1998/08/17 21:41:08 zender Exp $
c
      implicit none
      integer jslc
c
C------------------------------Parameters-------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
C
C------------------------------Commons----------------------------------
#     include "include/crdcae.cb"
#     include "include/crdcon.cb"
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(kind=8)  pbr(plond,plev),           ! Prssr at mid-levels (dynes/cm2)
     $     pnm(plond,plevp),          ! Prssr at interfaces (dynes/cm2)
     $     co2em(plond,plevp),        ! Co2 emissivity function
     $     co2eml(plond,plev),        ! Co2 emissivity function
     $     tplnka(plond,plevp),       ! Planck fnctn level temperature
     $     s2c(plond,plevp),          ! H2o continuum path length
     $     s2t(plond,plevp),          ! H2o tmp and prs wghted path
     $     w(plond,plevp),            ! H2o prs wghted path
     $     h2otr(plond,plevp),        ! H2o trnsmssn fnct for o3 overlap
     $     plco2(plond,plevp),        ! Co2 prs wghted path length
     $     plh2o(plond,plevp),        ! H2o prs wfhted path length
     $     co2t(plond,plevp),         ! Tmp and prs wghted path length
     $     tint(plond,plevp),         ! Interface temperatures
     $     tlayr(plond,plevp),        ! K-1 level temperatures
     $     plol(plond,plevp),         ! Ozone prs wghted path length
     $     plos(plond,plevp)          ! Ozone path length
      real(kind=8)  pmln(plond,plev),          ! Ln(pmidm1)
     $     piln(plond,plevp)          ! Ln(pintm1)
c
c   Trace gas variables
c
      real(kind=8)  ucfc11(plond,plevp), ! CFC11 path length
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
     $     uptype(plond,plevp), ! continuum path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp)    ! pressure factor for ch4
      real(kind=8)  abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
      real(kind=8)  abstrc(plond)            ! total trace gas absorptivity
      real(kind=8)  bplnk(14,plond,4) ! Planck functions for sub-divided layers
C
C Output arguments
C
      real(kind=8)  absnxt,              ! Nearest layer absorptivities
     $     abstot,              ! Non-adjacent layer absorptivites
     $     emstot               ! Total emissivity
#ifdef MPP1
      common /radbuf/ absnxt(plond,plev,4,jxp)
     &      , abstot(plond,plevp,plevp,jxp),emstot(plond,plevp,jxp)
#else
      common /radbuf/ absnxt(plond,plev,4,jlx)
     &      , abstot(plond,plevp,plevp,jlx),emstot(plond,plevp,jlx)
#endif
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        k1,           ! Level index
     $        k2,           ! Level index
     $        kn,           ! Nearest level index
     $        iband         ! Band  index
C
      real(kind=8)  pnew(plond),     ! Effective pressure for H2O vapor linewidth
     $     trline(plond,2), ! Transmission due to H2O lines in window
     $     u(plond),        ! Pressure weighted H2O path length
     $     tbar(plond,4),   ! Mean layer temperature
     $     emm(plond,4),    ! Mean co2 emissivity
     $     o3emm(plond,4),  ! Mean o3 emissivity
     $     o3bndi,          ! Ozone band parameter
     $     temh2o(plond,4), ! Mean layer temperature equivalent to tbar
     $     k21,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 650-800
C                           !  cm-1 region (tr1)
     $     k22,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 500-650
C                           !  cm-1 region (tr2)
     $     uc1(plond)       ! H2o continuum pathlength in 500-800 cm-1
      real(kind=8)  to3h2o(plond),   ! H2o trnsmsn for overlap with o3
     $     pi,              ! For co2 absorptivity computation
     $     sqti(plond),     ! Used to store sqrt of mean temperature
     $     et,              ! Co2 hot band factor
     $     et2,             ! Co2 hot band factor squared
     $     et4,             ! Co2 hot band factor to fourth power
     $     omet,            ! Co2 stimulated emission term
     $     f1co2,           ! Co2 central band factor
     $     f2co2(plond),    ! Co2 weak band factor
     $     f3co2(plond),    ! Co2 weak band factor
     $     t1co2(plond),    ! Overlap factr weak bands on strong band
     $     sqwp,            ! Sqrt of co2 pathlength
     $     f1sqwp(plond)    ! Main co2 band factor
      real(kind=8)  oneme,           ! Co2 stimulated emission term
     $     alphat,          ! Part of the co2 stimulated emission term
     $     wco2,            ! Constants used to define co2 pathlength
     $     posqt,           ! Effective pressure for co2 line width
     $     u7(plond),       ! Co2 hot band path length
     $     u8,              ! Co2 hot band path length
     $     u9,              ! Co2 hot band path length
     $     u13,             ! Co2 hot band path length
     $     rbeta7(plond),   ! Inverse of co2 hot band line width par
     $     rbeta8,          ! Inverse of co2 hot band line width par
     $     rbeta9,          ! Inverse of co2 hot band line width par
     $     rbeta13          ! Inverse of co2 hot band line width par
      real(kind=8)  tpatha(plond),     ! For absorptivity computation
     $     a,                 ! Eq(2) in table A3a of R&D
     $     abso(plond,6),     ! Absorptivity for various gases/bands
     $     dtp(plond),        ! Path temp minus 300 K used in h2o
C                             !  rotation band absorptivity
     $     dtx(plond),        ! Planck temperature minus 250 K
     $     dty(plond),        ! Path temperature minus 250 K
     $     dtz(plond),        ! Planck temperature minus 300 K
     $     term1(plond,4),    ! Equation(5) in table A3a of R&D(1986)
     $     term2(plond,4)     ! Delta a(Te) in table A3a of R&D(1986)
      real(kind=8)  term3(plond,4),    ! DB/dT function for rotation and
C                             !  vibration-rotation band absorptivity
     $     term4(plond,4),    ! Equation(6) in table A3a of R&D(1986)
     $     term5(plond,4),    ! Delta a(Tp) in table A3a of R&D(1986)
     $     term6(plond,plevp),! DB/dT function for window region
     $     term7(plond,2),    ! Kl_inf(i) in eq(8) of table A3a of R&D
     $     term8(plond,2),    ! Delta kl_inf(i) in eq(8)
     $     term9(plond,plevp),! DB/dT function for 500-800 cm-1 region
     $     tr1,               ! Eqn(6) in table A2 of R&D for 650-800
     $     tr10(plond),       ! Eqn (6) times eq(4) in table A2
C                             !  of R&D for 500-650 cm-1 region
     $     tr2                ! Eqn(6) in table A2 of R&D for 500-650
      real(kind=8)  tr5,               ! Eqn(4) in table A2 of R&D for 650-800
     $     tr6,               ! Eqn(4) in table A2 of R&D for 500-650
     $     tr9(plond),        ! Equation (6) times eq(4) in table A2
C                             !  of R&D for 650-800 cm-1 region
     $     uc(plond)          ! Y + 0.002U in eq(8) of table A2 of R&D
      real(kind=8)  sqrtu(plond),      ! Sqrt of pressure weighted h20 pathlength
     $     fwk(plond),        ! Equation(33) in R&D far wing correction
     $     fwku(plond),       ! GU term in eqs(1) and (6) in table A2
     $     r2st(2),           ! 1/(2*beta) in eq(10) in table A2
     $     dtyp15(plond),     ! DeltaTp in eqs(11) & (12) in table A3a
     $     dtyp15sq(plond),   ! (DeltaTp)^2 in eqs(11) & (12) table A3a
     $     to3co2(plond),     ! P weighted temp in ozone band model
     $     dpnm(plond),       ! Pressure difference between two levels
     $     pnmsq(plond,plevp),! Pressure squared
     $     dw(plond),         ! Amount of h2o between two levels
     $     uinpl(plond,4),    ! Nearest layer subdivision factor
     $     winpl(plond,4),    ! Nearest layer subdivision factor
     $     zinpl(plond,4),    ! Nearest layer subdivision factor
     $     pinpl(plond,4),    ! Nearest layer subdivision factor
     $     dplh2o(plond)      ! Difference in press weighted h2o amount
      real(kind=8)  r80257,            ! Conversion factor for h2o pathlength
     $     r293,              ! 1/293
     $     r250,              ! 1/250
     $     r3205,             ! Line width factor for o3 (see R&Di)
     $     r300,              ! 1/300
     $     rsslp,             ! Reciprocal of sea level pressure
     $     r2sslp             ! 1/2 of rsslp
      real(kind=8)   ds2c,     ! Y in eq(7) in table A2 of R&D
     $      a11,      ! A1 in table A3b for rotation band absorptivity
     $      a31,      ! A3 in table A3b for rotation band absorptivity
     $      a21,      ! First part in numerator of A2 in table A3b
     $      a22,      ! Second part in numerator of A2 in table A3b
     $      a23,      ! Denominator of A2 in table A3b (rotation band)
     $      t1t4,     ! Eq(3) in table A3a of R&D
     $      t2t5,     ! Eq(4) in table A3a of R&D
     $      rsum,     ! Eq(1) in table A2 of R&D
     $      a41,      ! Numerator in A2 in Vib-rot abstivity(table A3b)
     $      a51,      ! Denominator in A2 in Vib-rot (table A3b)
     $      a61       ! A3 factor for Vib-rot band in table A3b
      real(kind=8)   phi,      ! Eq(11) in table A3a of R&D
     $      psi,      ! Eq(12) in table A3a of R&D
     $      cf812,    ! Eq(11) in table A2 of R&D
     $      ubar,     ! H2o scaled path see comment for eq(10) table A2
     $      pbar,     ! H2o scaled pres see comment for eq(10) table A2
     $      g4        ! Arguement in exp() in eq(10) table A2
C
      real(kind=8)   dplos,    ! Ozone pathlength eq(A2) in R&Di
     $      dplol,    ! Presure weighted ozone pathlength
     $      tlocal,   ! Local interface temperature
     $      beta,     ! Ozone mean line parameter eq(A3) in R&Di
C                       (includes Voigt line correction factor)
     $      rphat,    ! Effective pressure for ozone beta
     $      tcrfac,   ! Ozone temperature factor table 1 R&Di
     $      tmp1,     ! Ozone band factor see eq(A1) in R&Di
     $      u1,       ! Effective ozone pathlength eq(A2) in R&Di
     $      realnu,   ! 1/beta factor in ozone band model eq(A1)
     $      tmp2,     ! Ozone band factor see eq(A1) in R&Di
     $      u2,       ! Effective ozone pathlength eq(A2) in R&Di
     $      rsqti     ! Reciprocal of sqrt of path temperature
C
      real(kind=8)   tpath,    ! Path temperature used in co2 band model
     $      tmp3,     ! Weak band factor see K&B
     $      rdpnmsq,  ! Reciprocal of difference in press^2
     $      rdpnm,    ! Reciprocal of difference in press
     $      p1,       ! Mean pressure factor
     $      p2,       ! Mean pressure factor
     $      dtym10,   ! T - 260 used in eq(9) and (10) table A3a
     $      dplco2,   ! Co2 pathlength
     $      corfac,   ! Correction factors in table A3b
     $      g2,       ! Part of arguement in eq(10) in table A2
     $      te,       ! A_0 T factor in ozone model table 1 of R&Di
     $      denom     ! Denominator in eq(8) of table A3a of R&D
c
c
      real(kind=8)  th2o(plond),
     $     tco2(plond),
     $     to3(plond)
      integer wvl
C
C
C Transmission terms for various spectral intervals:
C
C     real(kind=8)  trab1(plond)   ! H2o     0 -  800 cm-1
      real(kind=8)  trab2(plond)   ! H2o   500 -  800 cm-1
C     real(kind=8)  trab3(plond)   ! Co2   band system
      real(kind=8)  trab4(plond)   ! H2o   800 - 1000 cm-1
C     real(kind=8)  trab5(plond)   ! 9.6 micrometer band
      real(kind=8)  trab6(plond)   ! H2o  1000 - 1200 cm-1
C     real(kind=8)  trab7(plond)   ! H2o  1200 - 2200 cm-1
C
      real(kind=8)  bndfct, ! Band absorptance parameter for co2
     $     absbnd  ! Proportional to co2 band absorptance
C
      real(kind=8)  dbvtit(plond,plevp),      ! Intrfc drvtv plnck fnctn for o3
     $     dbvtly(plond,plev)        ! Level drvtv plnck fnctn for o3
C
C--------------------------Statement function---------------------------
C
      real(kind=8)  dbvt,t     ! Planck fnctn tmp derivative for o3
C
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
C
C-----------------------------------------------------------------------
C
C Initialize
C
      do k=1,plev
         do i=1,plon
            dbvtly(i,k) = dbvt(tlayr(i,k+1))
            dbvtit(i,k) = dbvt(tint(i,k))
         end do
      end do
      do i=1,plon
         dbvtit(i,plevp) = dbvt(tint(i,plevp))
      end do
C
      r80257  = 1./8.0257d-04
      r293    = 1./293.d0
      r250    = 1./250.d0
      r3205   = 1./.3205d0
      r300    = 1./300.d0
      rsslp   = 1./sslp
      r2sslp  = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
      bndfct  = 2.0*22.18d0/(dsqrt(196.d0)*300.)
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
      do k=1,plevp
         do i=1,plon
            pnmsq(i,k) = pnm(i,k)**2
            dtx(i) = tplnka(i,k) - 250.
            term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                   (1. +  c9*dtx(i)*(1. + c11*dtx(i)*
     $                   (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
            term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                   (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                   (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
         end do
      end do
C
C Non-nearest layer level loops
C
      do 200 k1=plevp,1,-1
         do 100 k2=plevp,1,-1
            if(k1.eq.k2) go to 100
            do i=1,plon
               dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
               u(i)      = dabs(dplh2o(i))
               sqrtu(i)  = dsqrt(u(i))
               ds2c      = dabs(s2c(i,k1) - s2c(i,k2))
               dw(i)     = dabs(w(i,k1) - w(i,k2))
               uc1(i)    = (ds2c + 1.7e-3*u(i))*(1. +  2.*ds2c)/
     $                                          (1. + 15.*ds2c)
               uc(i)     = ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               pnew(i)   = u(i)/dw(i)
               tpatha(i) = (s2t(i,k1) - s2t(i,k2))/dplh2o(i)
               dtx(i)      = tplnka(i,k2) - 250.
               dty(i)      = tpatha(i)    - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac  = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4    = term1(i,2)*term4(i,2)
               t2t5    = term2(i,2)*term5(i,2)
               a       = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)  = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i) = fwk(i)*u(i)
               rsum    = dexp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i)  = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = dexp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i)  = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = dexp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = dexp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(dsqrt(g2) - 1.)
                  trline(i,k) = dexp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               k21    = term7(i,1) + term8(i,1)/
     $             (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
               k22    = term7(i,2) + term8(i,2)/
     $             (1. + (c28 + c29*(dty(i)-10.))*sqrtu(i))
               tr1    = dexp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2    = dexp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5    = dexp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6    = dexp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)   = tr1*tr5
               tr10(i)  = tr2*tr6
               th2o(i) = tr10(i)
               trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
               trab4(i) = dexp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i) = dexp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               abso(i,3) = term6(i,k2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4) = term9(i,k2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
            end do
            if(k2.lt.k1) then
               do i=1,plon
                  to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
               end do
            else
               do i=1,plon
                  to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
               end do
            end if
C
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C
            do i=1,plon
               dpnm(i)  = pnm(i,k1) - pnm(i,k2)
               to3co2(i)=(pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/
     $                   dpnm(i)
               te       = (to3co2(i)*r293)**.7
               dplos    = plos(i,k1) - plos(i,k2)
               dplol    = plol(i,k1) - plol(i,k2)
               u1       = 18.29*dabs(dplos)/te
               u2       = .5649*dabs(dplos)/te
               rphat    = dplol/dplos
               tlocal   = tint(i,k2)
               tcrfac   = dsqrt(tlocal*r250)*te
               beta     = r3205*(rphat + dpfo3*tcrfac)
               realnu   = te/beta
               tmp1     = u1/dsqrt(4. + u1*(1. + realnu))
               tmp2     = u2/dsqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*dlog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
               to3(i)   = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i)  = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)      co2 15  micrometer band system
C
            do i=1,plon
               sqwp      = dsqrt(dabs(plco2(i,k1) - plco2(i,k2)))
               et        = dexp(-480./to3co2(i))
               sqti(i)   = dsqrt(to3co2(i))
               rsqti     = 1./sqti(i)
               et2       = et*et
               et4       = et2*et2
               omet      = 1. - 1.5*et2
               f1co2     = 899.70*omet*
     $                   (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i) = f1co2*sqwp
               t1co2(i)  = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme     = 1. - et2
               alphat    = oneme**3*rsqti
               pi        = dabs(dpnm(i))
               wco2      =  2.5221*co2vmr*pi*rga
               u7(i)     =  4.9411e4*alphat*et2*wco2
               u8        =  3.9744e4*alphat*et4*wco2
               u9        =  1.0447e5*alphat*et4*et2*wco2
               u13       = 2.8388e3*alphat*et4*wco2
               tpath     = to3co2(i)
               tlocal    = tint(i,k2)
               tcrfac    = dsqrt(tlocal*r250*tpath*r300)
               posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp +
     $                     dpfco2*tcrfac)*rsqti
               rbeta7(i) = 1./(5.3228*posqt)
               rbeta8    = 1./(10.6576*posqt)
               rbeta9    = rbeta7(i)
               rbeta13   = rbeta9
               f2co2(i)  = (u7(i)/dsqrt(4. + u7(i)*(1. + rbeta7(i)))) +
     $                     (u8/dsqrt(4. + u8*(1. + rbeta8))) +
     $                     (u9/dsqrt(4. + u9*(1. + rbeta9)))
               f3co2(i)  = u13/dsqrt(4. + u13*(1. + rbeta13))
            end do
            if (k2.ge.k1) then
               do i=1,plon
                  sqti(i) = dsqrt(tlayr(i,k2))
               end do
            end if
C
            do i=1,plon
               tmp1      = dlog(1. + f1sqwp(i))
               tmp2      = dlog(1. + f2co2(i))
               tmp3      = dlog(1. + f3co2(i))
               absbnd    = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
               tco2(i)=1./(1.0+10.0*(u7(i)
     &                         /dsqrt(4.+u7(i)*(1.+rbeta7(i)))))
C               trab3(i)  = 1. - bndfct*absbnd
            end do
c
c     Calculate absorptivity due to trace gases
c
            call trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213,
     $                         uco221, uco222, uco223, bn2o0,
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         u,      abplnk1,tco2,   th2o,
     $                         to3,    abstrc)
C
C Sum total absorptivity
C
            do i=1,plon
               abstot(i,k1,k2,jslc) =
     &                           abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  100    continue
  200 continue        ! End of non-nearest layer level loops
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
C Nearest layer level loop
C
      do 500 k2=plev,1,-1
         do i=1,plon
            tbar(i,1)   = 0.5*(tint(i,k2+1) + tlayr(i,k2+1))
            emm(i,1)    = 0.5*(co2em(i,k2+1) + co2eml(i,k2))
            tbar(i,2)   = 0.5*(tlayr(i,k2+1) + tint(i,k2))
            emm(i,2)    = 0.5*(co2em(i,k2) + co2eml(i,k2))
            tbar(i,3)   = 0.5*(tbar(i,2) + tbar(i,1))
            emm(i,3)    = emm(i,1)
            tbar(i,4)   = tbar(i,3)
            emm(i,4)    = emm(i,2)
            o3emm(i,1)  = 0.5*(dbvtit(i,k2+1) + dbvtly(i,k2))
            o3emm(i,2)  = 0.5*(dbvtit(i,k2) + dbvtly(i,k2))
            o3emm(i,3)  = o3emm(i,1)
            o3emm(i,4)  = o3emm(i,2)
            temh2o(i,1) = tbar(i,1)
            temh2o(i,2) = tbar(i,2)
            temh2o(i,3) = tbar(i,1)
            temh2o(i,4) = tbar(i,2)
            dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
         end do
c----------------------------------------------------------
c  Weighted Planck functions for trace gases
c
      do wvl = 1,14
         do i = 1,plon
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1) +
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2) +
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
         end do
      end do
c---------------------------------------------------------
         do i=1,plon
            rdpnmsq    = 1./(pnmsq(i,k2+1) - pnmsq(i,k2))
            rdpnm      = 1./dpnm(i)
            p1         = .5*(pbr(i,k2) + pnm(i,k2+1))
            p2         = .5*(pbr(i,k2) + pnm(i,k2  ))
            uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
            uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
            uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
            uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
            winpl(i,1) = (.5*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
            winpl(i,2) = (.5*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
            winpl(i,3) = (.5*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*
     $                   rdpnm
            winpl(i,4) = (.5*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*
     $                   rdpnm
            tmp1       = 1./(piln(i,k2+1) - piln(i,k2))
            tmp2       = piln(i,k2+1) - pmln(i,k2)
            tmp3       = piln(i,k2  ) - pmln(i,k2)
            zinpl(i,1) = (.5*tmp2          )*tmp1
            zinpl(i,2) = (        - .5*tmp3)*tmp1
            zinpl(i,3) = (.5*tmp2 -    tmp3)*tmp1
            zinpl(i,4) = (   tmp2 - .5*tmp3)*tmp1
            pinpl(i,1) = 0.5*(p1 + pnm(i,k2+1))
            pinpl(i,2) = 0.5*(p2 + pnm(i,k2  ))
            pinpl(i,3) = 0.5*(p1 + pnm(i,k2  ))
            pinpl(i,4) = 0.5*(p2 + pnm(i,k2+1))
         end do
         do 400 kn=1,4
            do i=1,plon
               u(i)     = uinpl(i,kn)*dabs(plh2o(i,k2) - plh2o(i,k2+1))
               sqrtu(i) = dsqrt(u(i))
               dw(i)    = dabs(w(i,k2) - w(i,k2+1))
               pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
               ds2c     = dabs(s2c(i,k2) - s2c(i,k2+1))
               uc1(i)   = uinpl(i,kn)*ds2c
               uc1(i)   = (uc1(i) + 1.7e-3*u(i))*(1. +  2.*uc1(i))/
     $                                           (1. + 15.*uc1(i))
               uc(i)    = uinpl(i,kn)*ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               dtx(i)      = temh2o(i,kn) - 250.
               dty(i)      = tbar(i,kn) - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac    = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4      = term1(i,2)*term4(i,2)
               t2t5      = term2(i,2)*term5(i,2)
               a         = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)    = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i)   = fwk(i)*u(i)
               rsum      = dexp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i) = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = dexp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i) = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = dexp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = dexp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*winpl(i,kn)*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(dsqrt(g2) - 1.)
                  trline(i,k) = dexp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               dtym10     = dty(i) - 10.
               denom      = 1. + (c30 + c31*dtym10*dtym10)*sqrtu(i)
               k21        = term7(i,1) + term8(i,1)/denom
               denom      = 1. + (c28 + c29*dtym10       )*sqrtu(i)
               k22        = term7(i,2) + term8(i,2)/denom
               term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                     (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                      (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
               tr1     = dexp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2     = dexp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5     = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
               tr6     = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)  = tr1*tr5
               tr10(i) = tr2*tr6
               trab2(i)= 0.65*tr9(i) + 0.35*tr10(i)
               th2o(i) = tr10(i)
               trab4(i)= dexp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i)= dexp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                     (1. + c9*dtx(i)*(1. + c11*dtx(i)*
     $                     (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
               abso(i,3)  = term6(i,2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4)  = term9(i,2)*.5*(tr1-tr9(i)+tr2-tr10(i))
            end do
C
C abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
C
            do i=1,plon
               te        = (tbar(i,kn)*r293)**.7
               dplos     = dabs(plos(i,k2+1) - plos(i,k2))
               u1        = zinpl(i,kn)*18.29*dplos/te
               u2        = zinpl(i,kn)*.5649*dplos/te
               tlocal    = tbar(i,kn)
               tcrfac    = dsqrt(tlocal*r250)*te
               beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
               realnu    = te/beta
               tmp1      = u1/dsqrt(4. + u1*(1. + realnu))
               tmp2      = u2/dsqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*dlog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*o3emm(i,kn)*
     $                     (h2otr(i,k2+1)/h2otr(i,k2))
               to3(i)    = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i) = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)   co2 15  micrometer band system
C
            do 300 i=1,plon
               dplco2   = plco2(i,k2+1) - plco2(i,k2)
               sqwp     = dsqrt(uinpl(i,kn)*dplco2)
               et       = dexp(-480./tbar(i,kn))
               sqti(i)  = dsqrt(tbar(i,kn))
               rsqti    = 1./sqti(i)
               et2      = et*et
               et4      = et2*et2
               omet     = (1. - 1.5*et2)
               f1co2    = 899.70*omet*
     $                    (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i)= f1co2*sqwp
               t1co2(i) = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme    = 1. - et2
               alphat   = oneme**3*rsqti
               pi       = dabs(dpnm(i))*winpl(i,kn)
               wco2     = 2.5221*co2vmr*pi*rga
               u7(i)    = 4.9411e4*alphat*et2*wco2
               u8       = 3.9744e4*alphat*et4*wco2
               u9       = 1.0447e5*alphat*et4*et2*wco2
               u13      = 2.8388e3*alphat*et4*wco2
               tpath    = tbar(i,kn)
               tlocal   = tbar(i,kn)
               tcrfac   = dsqrt((tlocal*r250)*(tpath*r300))
               posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
               rbeta7(i)= 1./(5.3228*posqt)
               rbeta8   = 1./(10.6576*posqt)
               rbeta9   = rbeta7(i)
               rbeta13  = rbeta9
               f2co2(i) = u7(i)/dsqrt(4. + u7(i)*(1. + rbeta7(i))) +
     $                    u8/dsqrt(4. + u8*(1. + rbeta8)) +
     $                    u9/dsqrt(4. + u9*(1. + rbeta9))
               f3co2(i) = u13/dsqrt(4. + u13*(1. + rbeta13))
               tmp1     = dlog(1. + f1sqwp(i))
               tmp2     = dlog(1. + f2co2(i))
               tmp3     = dlog(1. + f3co2(i))
               absbnd   = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6)= trab2(i)*emm(i,kn)*absbnd
               tco2(i)=1.0/(1.0+ 10.0*u7(i)
     &                          /dsqrt(4.+u7(i)*(1.+rbeta7(i))))
C               trab3(i) = 1. - bndfct*absbnd
  300       continue
c
c   Calculate trace gas absorptivity for nearest layer
c
            call trcabn(k2, kn, ucfc11, ucfc12, un2o0,  un2o1,
     $                          uch4,   uco211, uco212, uco213,
     $                          uco221, uco222, uco223,
     $                          tbar,   bplnk,
     $                          winpl,  pinpl,  tco2,   th2o,
     $                          to3,    uptype, dw,     s2c,
     $                          u,      pnew,   abstrc     ,uinpl)
C
C Total next layer absorptivity:
C
            do i=1,plon
               absnxt(i,k2,kn,jslc) =
     &                           abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  400    continue
  500 continue                  !  end of nearest layer level loop
C
      return
      end
