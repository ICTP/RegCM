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
!
      subroutine zengocn(u,ts,t,q,hgt,zi,ps,qs,u10,tau,                 &
                       & alh,ash,dth,dqh,ustar,zo)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
!     U. of Arizona (UA)   Bulk Aerodynamic Algorithm
!     12/22/97
!     --------
!     Add the `zo/L' term to the equations, because in a global modeling
!     environment (particularly over land), the absolute value of L could
!     be smaller than zo (either zom or zoh)   (11/04/98)
!     Add restriction: -100 <= z/L <= 2        (02/23/00)
!
!     Reference: Zeng et al. 1998, Intercomparison of bulk aerodynamic
!     algorithms for the computation of sea surface fluxes
!     using the TOGA COARE and TAO data. J. Climate,
!     11, 2628-2644.
!
!     For additional information, contact
!     Prof. Xubin Zeng
!     Department of Atmospheric Science
!     PAS Building, #81
!     The University of Arizona
!     Tucson, AZ 85721
!     USA
!     Tel:520-621-4782
!     Email:xubin@gogo.atmo.arizona.edu
!
!     input:
!       u   = sqrt(u_x^2 + u_y^2): wind speed in m/s at hu (m) height
!       ts: surface temperature in (deg C)
!       t:  air temperature in (deg C) at ht (m) height
!       q: air specific humidity in (kg/kg) at hq (m) height
!     output:
!       u10: wind speed at 10 meter (m/s)
!       tau: wind stress (N/m2)
!       alh: latent heat flux (W/m2)
!       ash: sensible heat flux (W/m2)
!       dth: air-surface potential temperature difference
!       dqh: air-surface specific humidity difference
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      use mod_constants , only : gti , rgti , rgas , rovcp , vonkar ,   &
                                 cpd , rwwd , tmelt
      implicit none
!
! Dummy arguments
!
      real(kind=8) , intent (in) :: hgt , q , t , u , zi , ts , ps
      real(kind=8) , intent (out) :: alh , ash , tau , u10
      real(kind=8) , intent (inout) :: dqh , dth , qs , ustar , zo
!
! Local variables
!
      real(kind=8) :: dthv , hq , ht , hu , obu , pr , qstar , rb ,     &
               & rho , th , thv , thvstar , tstar , um , visa , zot ,   &
               & wc , xlv , z10 , zbeta , zeta , zetam , zetat , zoq
      real(kind=8), external :: psi, qsat
      integer :: i
!
!***********************************************************************
!
      zbeta = 1.   ! -  (in computing W_*)
      pr = 0.71    ! =nu/thermal diffusivity (the Prandtl number)
      z10 = 10.    ! m  (reference height)
!
      hu = hgt
      ht = hgt
      hq = hgt
!
      th = (t+tmelt)*(1000./ps)**rovcp
      ! potential T
      dth = t + 0.0098*ht - ts
      qs = qsat(ts,ps)*0.98
      qs = rwwd*qs/(ps-0.378*qs)
      ! in kg/kg
      dqh = q - qs
      thv = th*(1.+0.61*q)
      ! virtual potential T
      dthv = dth*(1.+0.61*q) + 0.61*th*dqh
      rho = ps*100./(rgas*(ts+tmelt)*(1.+0.61*qs))
      ! density
      xlv = (2.501-0.00237*ts)*1.E+6
      ! J/kg
!
!     Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep.
!     89-11
!
      visa = 1.326E-5*(1+6.542E-3*t+8.301E-6*t*t-4.84E-9*t*t*t)
!
!     initial values of u* and convective velocity
!
      ustar = 0.06
      wc = 0.5
      if ( dthv.ge.0. ) then
        um = max(u,0.1D0)
      else
        um = sqrt(u*u+wc*wc)
      end if
!
!     loop to obtain initial and good ustar and zo
!
      do i = 1 , 5
        zo = 0.013*ustar*ustar*rgti + 0.11*visa/ustar
        ustar = vonkar*um/dlog(hu/zo)
      end do
!
      rb = gti*hu*dthv/(thv*um*um)
      if ( rb.ge.0. ) then       ! neutral or stable
        zeta = rb*dlog(hu/zo)/(1.-5.*dmin1(rb,0.19D0))
        zeta = dmin1(2.D0,dmax1(zeta,1.D-6))
      else                      !unstable
        zeta = rb*dlog(hu/zo)
        zeta = dmax1(-100.D0,dmin1(zeta,-1.D-6))
      end if
      obu = hu/zeta
!
!     main iterations (2-10 iterations would be fine)
!
      do i = 1 , 10
        call rough(zo,zot,zoq,ustar,visa,gti)
!
!       wind
!
        zeta = hu/obu
        zetam = 1.574
        if ( zeta.lt.-zetam ) then
                                 ! zeta < -1
          ustar = vonkar*um/(dlog(-zetam*obu/zo)-psi(1,-zetam)+         &
                & psi(1,zo/obu)+1.14*((-zeta)**0.333-(zetam)**0.333))
        else if ( zeta.lt.0. ) then
                                  ! -1 <= zeta < 0
          ustar = vonkar*um/(dlog(hu/zo)-psi(1,zeta)+psi(1,zo/obu))
        else if ( zeta.le.1. ) then
                                  !  0 <= zeta <= 1
          ustar = vonkar*um/(dlog(hu/zo)+5.*zeta-5.*zo/obu)
        else                   !  1 < zeta, phi=5+zeta
          ustar = vonkar*um/(dlog(obu/zo)+5.-5.*zo/obu+                 &
                & (5.*dlog(zeta)+zeta-1.))
        end if
!
!       temperature
!
        zeta = ht/obu
        zetat = 0.465
        if ( zeta.lt.-zetat ) then
                                 ! zeta < -1
          tstar = vonkar*dth/(dlog(-zetat*obu/zot)-psi(2,-zetat)        &
                & +psi(2,zot/obu)                                       &
                & +0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
        else if ( zeta.lt.0. ) then
                                  ! -1 <= zeta < 0
          tstar = vonkar*dth/(dlog(ht/zot)-psi(2,zeta)+psi(2,zot/obu))
        else if ( zeta.le.1. ) then
                                  !  0 <= ztea <= 1
          tstar = vonkar*dth/(dlog(ht/zot)+5.*zeta-5.*zot/obu)
        else                   !  1 < zeta, phi=5+zeta
          tstar = vonkar*dth/(dlog(obu/zot)+5.-5.*zot/obu+              &
                & (5.*dlog(zeta)+zeta-1.))
        end if
!
!       humidity
!
        zeta = hq/obu
        zetat = 0.465
        if ( zeta.lt.-zetat ) then
                                 ! zeta < -1
          qstar = vonkar*dqh/(dlog(-zetat*obu/zoq)-psi(2,-zetat)        &
                & +psi(2,zoq/obu)                                       &
                & +0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
        else if ( zeta.lt.0. ) then
                                  ! -1 <= zeta < 0
          qstar = vonkar*dqh/(dlog(hq/zoq)-psi(2,zeta)+psi(2,zoq/obu))
        else if ( zeta.le.1. ) then
                                  !  0 <= ztea <= 1
          qstar = vonkar*dqh/(dlog(hq/zoq)+5.*zeta-5.*zoq/obu)
        else                   !  1 < zeta, phi=5+zeta
          qstar = vonkar*dqh/(dlog(obu/zoq)+5.-5.*zoq/obu+              &
                & (5.*dlog(zeta)+zeta-1.))
        end if
        thvstar = tstar*(1.+0.61*q) + 0.61*th*qstar
!
        zeta = vonkar*gti*thvstar*hu/(ustar**2*thv)
        if ( zeta.ge.0 ) then   !neutral or stable
          um = max(u,0.1D0)
          zeta = dmin1(2.D0,max(zeta,1.D-6))
        else                   !unstable
          wc = zbeta*(-gti*ustar*thvstar*zi/thv)**0.333
          um = sqrt(u*u+wc*wc)
          zeta = dmax1(-100.D0,min(zeta,-1.D-6))
        end if
        obu = hu/zeta
      end do
!
!--------------------------------------------------------------
!
      tau = rho*ustar*ustar*u/um
      alh = -rho*xlv*qstar*ustar
      ash = -rho*cpd*tstar*ustar
!
!     x and y components of tau:
!     taux=rho*ustar*ustar*u_x/um
!     tauy=rho*ustar*ustar*u_y/um
!     10-meter wind (without w_* part)
!
      zeta = z10/obu
      if ( zeta.lt.0. ) then
        u10 = u + (ustar/vonkar)*(dlog(z10/hu)-(psi(1,zeta)-            &
              & psi(1,hu/obu)))
      else
        u10 = u + (ustar/vonkar)*(dlog(z10/hu)+5.*zeta-5.*hu/obu)
      end if
      end subroutine zengocn
!
! stability function for rb < 0
!
      function psi(k,zeta)
      implicit none
!
! Dummy arguments
!
      integer , intent(in) :: k
      real(kind=8) , intent(in) :: zeta
      real(kind=8) :: psi
!
! Local variables
!
      real(kind=8) :: chik
!
      chik = (1.-16*zeta)**0.25
      if ( k.eq.1 ) then
        psi = 2.*dlog((1.+chik)*0.5) + dlog((1.+chik*chik)*0.5)         &
            & - 2.*atan(chik) + 2.*atan(1.)
      else
        psi = 2.*dlog((1.+chik*chik)*0.5)
      end if
      end function psi

!
! Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532
! p in mb, t in C, and qsat in mb
!
      function qsat(t,p)
      implicit none
!
! Dummy arguments
!
      real(kind=8) , intent (in) :: p , t
      real(kind=8) :: qsat
!
      qsat = (1.0007+3.46E-6*p)*6.1121*exp(17.502*t/(240.97+t))
!
      end function qsat

!
!  our formulation for zo,zot,zoq
!
      subroutine rough(zo,zot,zoq,ustar,visa,g)
!
      implicit none
!
! Dummy arguments
!
      real(kind=8) , intent (in) :: g , ustar , visa
      real(kind=8) , intent (out) :: zoq , zot
      real(kind=8) , intent (inout) :: zo
!
! Local variables
!
      real(kind=8) :: re , xq , xt
!
!Im
!     zo=0.013*ustar*ustar/g+0.11*visa/ustar
!     zo=0.013*ustar*ustar/g
      zo = 0.0065*ustar*ustar/g
!Im_
      re = ustar*zo/visa
      xq = 2.67*re**0.25 - 2.57
      xt = xq
      zoq = zo/exp(xq)
      zot = zo/exp(xt)
      end subroutine rough
