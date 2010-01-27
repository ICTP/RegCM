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
C
      subroutine zengocn(u,ts,t,q,hgt,zi,ps,qs,ak,grav,rgas,cpa,u10,tau,
     &                   alh,ash,dth,dqh,ustar,zo)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C
C     U. of Arizona (UA)   Bulk Aerodynamic Algorithm
C     12/22/97
C     --------
C     Add the `zo/L' term to the equations, because in a global modeling
C     environment (particularly over land), the absolute value of L could
C     be smaller than zo (either zom or zoh)   (11/04/98)
C     Add restriction: -100 <= z/L <= 2        (02/23/00)
C
C     Reference: Zeng et al. 1998, Intercomparison of bulk aerodynamic
C     algorithms for the computation of sea surface fluxes
C     using the TOGA COARE and TAO data. J. Climate,
C     11, 2628-2644.
C
C     For additional information, contact
C     Prof. Xubin Zeng
C     Department of Atmospheric Science
C     PAS Building, #81
C     The University of Arizona
C     Tucson, AZ 85721
C     USA
C     Tel:520-621-4782
C     Email:xubin@gogo.atmo.arizona.edu
C
C     input:
C       u   = sqrt(u_x^2 + u_y^2): wind speed in m/s at hu (m) height
C       ts: surface temperature in (deg C)
C       t:  air temperature in (deg C) at ht (m) height
C       q: air specific humidity in (kg/kg) at hq (m) height
C     output:
C       u10: wind speed at 10 meter (m/s)
C       tau: wind stress (N/m2)
C       alh: latent heat flux (W/m2)
C       ash: sensible heat flux (W/m2)
C       dth: air-surface potential temperature difference
C       dqh: air-surface specific humidity difference
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      implicit none
C
C     Input
C
      real(kind=8) u , ts , t , q , hgt , zi , ps , ak , grav , rgas , 
     &             cpa
C
C     Output
C
      real(kind=8) u10 , tau , alh , ash , dth , dqh , ustar , zo
C
C     Local variables
C
      real(kind=8) zbeta , tok , pr , z10 , hu , ht , hq , th , thv , 
     &             dthv , rho , qsat , qs , xlv , visa , wc , um , rb , 
     &             zeta , obu , zot , zoq , zetam , zetat , psi , 
     &             tstar , qstar , thvstar
      integer i
C
C***********************************************************************
C
      zbeta = 1.   ! -  (in computing W_*)
      tok = 273.16 ! Celsius to Kelvin
      pr = 0.71    ! =nu/thermal diffusivity (the Prandtl number)
      z10 = 10.    ! m  (reference height)
C
      hu = hgt
      ht = hgt
      hq = hgt
C
      th = (t+tok)*(1000./ps)**(rgas/cpa)
                                        !potential T
      dth = t + 0.0098*ht - ts
      qs = qsat(ts,ps)*0.98
      qs = .62197*qs/(ps-0.378*qs)
                                 ! in kg/kg
      dqh = q - qs
      thv = th*(1.+0.61*q)
                          ! virtual potential T
      dthv = dth*(1.+0.61*q) + 0.61*th*dqh
      rho = ps*100./(rgas*(ts+tok)*(1.+0.61*qs))
                                               ! density
      xlv = (2.501-0.00237*ts)*1.E+6
                                   !J/kg
C
C     Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep.
C     89-11
C
      visa = 1.326E-5*(1+6.542E-3*t+8.301E-6*t*t-4.84E-9*t*t*t)
C
C     initial values of u* and convective velocity
C
      ustar = 0.06
      wc = 0.5
      if ( dthv>=0. ) then
        um = max(u,0.1D0)
      else
        um = sqrt(u*u+wc*wc)
      end if
C
C     loop to obtain initial and good ustar and zo
C
      do i = 1 , 5
        zo = 0.013*ustar*ustar/grav + 0.11*visa/ustar
        ustar = ak*um/dlog(hu/zo)
      end do
C
      rb = grav*hu*dthv/(thv*um*um)
      if ( rb>=0. ) then         ! neutral or stable
        zeta = rb*dlog(hu/zo)/(1.-5.*dmin1(rb,0.19D0))
        zeta = dmin1(2.D0,dmax1(zeta,1.D-6))
      else                      !unstable
        zeta = rb*dlog(hu/zo)
        zeta = dmax1(-100.D0,dmin1(zeta,-1.D-6))
      end if
      obu = hu/zeta
C
C     main iterations (2-10 iterations would be fine)
C
      do i = 1 , 10
        call rough(zo,zot,zoq,ustar,visa,grav)
C
C       wind
C
        zeta = hu/obu
        zetam = 1.574
        if ( zeta<-zetam ) then  ! zeta < -1
          ustar = ak*um/(dlog(-zetam*obu/zo)-psi(1,-zetam)+psi(1,zo/obu)
     &            +1.14*((-zeta)**0.333-(zetam)**0.333))
        else if ( zeta<0. ) then  ! -1 <= zeta < 0
          ustar = ak*um/(dlog(hu/zo)-psi(1,zeta)+psi(1,zo/obu))
        else if ( zeta<=1. ) then !  0 <= zeta <= 1
          ustar = ak*um/(dlog(hu/zo)+5.*zeta-5.*zo/obu)
        else                   !  1 < zeta, phi=5+zeta
          ustar = ak*um/(dlog(obu/zo)+5.-5.*zo/obu+
     &            (5.*dlog(zeta)+zeta-1.))
        end if
C
C       temperature
C
        zeta = ht/obu
        zetat = 0.465
        if ( zeta<-zetat ) then  ! zeta < -1
          tstar = ak*dth/(dlog(-zetat*obu/zot)-psi(2,-zetat)
     &            +psi(2,zot/obu)
     &            +0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
        else if ( zeta<0. ) then  ! -1 <= zeta < 0
          tstar = ak*dth/(dlog(ht/zot)-psi(2,zeta)+psi(2,zot/obu))
        else if ( zeta<=1. ) then !  0 <= ztea <= 1
          tstar = ak*dth/(dlog(ht/zot)+5.*zeta-5.*zot/obu)
        else                   !  1 < zeta, phi=5+zeta
          tstar = ak*dth/(dlog(obu/zot)+5.-5.*zot/obu+
     &            (5.*dlog(zeta)+zeta-1.))
        end if
C
C       humidity
C
        zeta = hq/obu
        zetat = 0.465
        if ( zeta<-zetat ) then  ! zeta < -1
          qstar = ak*dqh/(dlog(-zetat*obu/zoq)-psi(2,-zetat)
     &            +psi(2,zoq/obu)
     &            +0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
        else if ( zeta<0. ) then  ! -1 <= zeta < 0
          qstar = ak*dqh/(dlog(hq/zoq)-psi(2,zeta)+psi(2,zoq/obu))
        else if ( zeta<=1. ) then !  0 <= ztea <= 1
          qstar = ak*dqh/(dlog(hq/zoq)+5.*zeta-5.*zoq/obu)
        else                   !  1 < zeta, phi=5+zeta
          qstar = ak*dqh/(dlog(obu/zoq)+5.-5.*zoq/obu+
     &            (5.*dlog(zeta)+zeta-1.))
        end if
        thvstar = tstar*(1.+0.61*q) + 0.61*th*qstar
C
        zeta = ak*grav*thvstar*hu/(ustar**2*thv)
        if ( zeta>=0 ) then     !neutral or stable
          um = max(u,0.1D0)
          zeta = dmin1(2.D0,max(zeta,1.D-6))
        else                   !unstable
          wc = zbeta*(-grav*ustar*thvstar*zi/thv)**0.333
          um = sqrt(u*u+wc*wc)
          zeta = dmax1(-100.D0,min(zeta,-1.D-6))
        end if
        obu = hu/zeta
      end do
C
C--------------------------------------------------------------
C
      tau = rho*ustar*ustar*u/um
      alh = -rho*xlv*qstar*ustar
      ash = -rho*cpa*tstar*ustar
C
C     x and y components of tau:
C     taux=rho*ustar*ustar*u_x/um
C     tauy=rho*ustar*ustar*u_y/um
C     10-meter wind (without w_* part)
C
      zeta = z10/obu
      if ( zeta<0. ) then
        u10 = u + (ustar/ak)*(dlog(z10/hu)-(psi(1,zeta)-psi(1,hu/obu)))
      else
        u10 = u + (ustar/ak)*(dlog(z10/hu)+5.*zeta-5.*hu/obu)
      end if
      end subroutine zengocn
C
C stability function for rb < 0
C
      function psi(k,zeta)
      implicit none
C
      integer k
      real(kind=8) zeta
C
      real(kind=8) psi
C
      real(kind=8) chik
C
      chik = (1.-16*zeta)**0.25
      if ( k==1 ) then
        psi = 2.*dlog((1.+chik)*0.5) + dlog((1.+chik*chik)*0.5)
     &        - 2.*atan(chik) + 2.*atan(1.)
      else
        psi = 2.*dlog((1.+chik*chik)*0.5)
      end if
      end function psi
C
C Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532
C p in mb, t in C, and qsat in mb
C
      function qsat(t,p)
      implicit none
C
      real(kind=8) qsat
      real(kind=8) t , p
C
      qsat = (1.0007+3.46E-6*p)*6.1121*exp(17.502*t/(240.97+t))
C
      end function qsat
C
C  our formulation for zo,zot,zoq
C
      subroutine rough(zo,zot,zoq,ustar,visa,g)
C
      implicit none
C
      real(kind=8) zo , zot , zoq , ustar , visa , g
      real(kind=8) re , xq , xt
CIm
C     zo=0.013*ustar*ustar/g+0.11*visa/ustar
C     zo=0.013*ustar*ustar/g
      zo = 0.0065*ustar*ustar/g
CIm_
      re = ustar*zo/visa
      xq = 2.67*re**0.25 - 2.57
      xt = xq
      zoq = zo/exp(xq)
      zot = zo/exp(xt)
      end subroutine rough
