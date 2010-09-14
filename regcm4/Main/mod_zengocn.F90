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

      module mod_zengocn
!
! Ocean flux model
! Implement Zeng and Beljaars, GRL , 2005, ZB2005
!
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_pbldim
      use mod_slice
      use mod_bats
      use mod_date
      use mod_bdycod
#ifdef CLM
      use clm_varsur , only : landmask
#endif
!
      private
!
      public :: zengocndrv
!
      contains
!
      subroutine zengocndrv(j , ng , istart , iend , k)
!
      implicit none
!
      integer , intent (in) :: j , ng , istart , iend , k
!
      real(kind=8) :: dqh , dth , facttq , lh , psurf , q995 , qs , sh ,&
               & t995 , tau , tsurf , ustar , uv10 , uv995 , z995 , zi ,&
               & zo
      integer :: i , n
#ifdef CLM
      integer :: jj
#endif
!     Implement Zeng and Beljaars, GRL , 2005, ZB2005
!     Account for SST diurnal evoluation warm layer/ skin temperature
!     scheme
!     real(8) :: lwds , lwus
      real(8) :: rs , rd , td , tdelta , delta
      real(8) :: q , ustarw , fd , l , phidl , aa , bb , cc , lamb
      real(8) :: dtstend , dts , fs , tskin , dtsst
      real(8) , parameter :: a1 = 0.28D+00
      real(8) , parameter :: a2 = 0.27D+00
      real(8) , parameter :: a3 = 0.45D+00
      real(8) , parameter :: b1 = 71.5D+00
      real(8) , parameter :: b2 = 2.8D+00
      real(8) , parameter :: b3 = 0.07D+00
      real(8) , parameter :: alphaw = 0.207D-06
      real(8) , parameter :: nuw = 1.004D-06
      real(8) , parameter :: kw = 0.60
      real(8) , parameter :: nu = 0.3
      real(8) , parameter :: d = 3 ! reference depth for bulk SST
!
#ifdef CLM
      jj = (jxp*myid) + j
#endif
      do i = istart , iend
        do n = 1 , ng
#ifdef CLM
          if ( ocld2d(n,i,j).lt.0.5 .or. landmask(jj,i).eq.3 ) then
#else
          if ( ocld2d(n,i,j).lt.0.5 ) then
#endif
            uv995 = sqrt(ubx3d(i,k,j)**2+vbx3d(i,k,j)**2)
            tsurf = tgb(i,j) - tzero
            t995 = tb3d(i,k,j) - tzero
            q995 = qvb3d(i,k,j)/(1.+qvb3d(i,k,j))
            z995 = za(i,k,j)
            zi = zpbl(i,j)
            psurf = (psb(i,j)+r8pt)*10.
            call zengocn(uv995,tsurf,t995,q995,z995,zi,psurf,qs,        &
                       & uv10,tau,lh,sh,dth,dqh,ustar,zo)
            if (idcsst == 1) then
!             time step considered for the integration of prognostic skin
!             temperature , equal to BATS time step
              dtsst = dtbat
!             handle the first call of the scheme
              if ( .not.firstcall(i,j) ) then
                deltas(i,j) = 0.001
                tdeltas(i,j) = tgb(i,j) - 0.001
                firstcall(i,j) = .true.
                td = tdeltas(i,j)
              end if
!             Init local variables
              delta = deltas(i,j)
              tdelta = tdeltas(i,j)
!             td is now the 3m bulk SST from the forcing variable
              td = ts1(i,j)
!
!             deep impact of aod on sst
!             if ( sum(aerext(i,:,j)).le.1 ) then
!               td = ts1(i,j) - sum(aerext(i,:,j))*0.8
!             else if ( sum(aerext(i,:,j)).gt.1 ) then
!               td = ts1(i,j)- 1.*0.8
!             end if
!
!             rs is the net surface sw flux (sw energy absorbed)
              rs = fsw2d(i,j)
!             rd is sw flux at 3m
              rd = rs*(a1*dexp(-d*b1) + a2*dexp(-d*b2) + a3*dexp(-d*b3))
!             ustar water (with air density ==1)
              ustarw = 0.5*ustar*(rhox2d(i,j)/rhoh2o)**0.5
!             lwds =  flwd2d(i,j)
!             lwus =  emsw*sigm*(tsurf+273.16)**4
!             q is the skin cooling term inckude net lw flux from
!             the radiative scheme
!             q = -(lh+sh+(lwus-lwds))
              q = -(lh+sh+flw2d(i,j))
!             fraction of solar radiation abosrbed in the sublayer
              fs = 0.065+11.*delta-(6.6e-5/delta)*(1-dexp(-delta/8.e-4))
!             dts= temperature difference between bulk level and skin level
!                determined from previous time step (via tdelta and td)
              dts = tdelta-td
!             m.o lenght calculation
              if ( dts.gt.0 ) then
                fd = (nu*gti*alphaw/(5*d))**0.5*                       &
                    &  rhoh2o*cpw0*ustarw**2*dts**0.5
              else
                fd = gti*alphaw*(q+rs-rd)
              end if
              l = rhoh2o*cpw0*ustarw**3/(vonkar*fd)
!             calulation of phidl (stability function)
              if ( (d/l).ge.0 ) then
                phidl = 1+5.*(d/l)
              else
                phidl = (1-16.*(d/l))**(-0.5)
              end if
!             prognostic evolution of dts
!             we can split the tendencies ddts/dt = a - b * dts
!             with a and b are ultimately function of dts through q
              aa = (q + rs - rd) / (d * cpw0 * rhoh2o * nu/(nu+1))
              bb = (nu+1) * vonkar * ustarw / (d*phidl)
!             exponential solution
              dtstend = aa - dts*(1-dexp(-bb*dtsst))/dtsst
!             update dts
              dts = dts + dtstend * dtsst
!             update tdelta
              tdelta = dts + td
!             update delta thickness  and cool skin tempearture
              aa = -16.*gti*alphaw*rhoh2o*cpw0*nuw**3./                 &
                  &     (ustarw**4. *kw**2.)
              bb =  aa *(q+rs*fs)
              if ( bb.gt.0 ) then
!               case of cool skin layer correction
                cc= bb**(3./4.)
                lamb=6.*( (1.+(aa*(q+rs*fs))**0.75)**(-0.333))
                delta = lamb*nuw/ustarw
                tskin= delta/(rhoh2o*cpw0*kw)*(q+rs*fs) + tdelta
              else
!               no cool skin layer in this case, tskin = warm layer
!               temperature
                tskin=tdelta
              end if
!             save the temperature difference and skin layer thickness
!             for next time step
              deltas(i,j) = delta
              tdeltas(i,j) = tdelta
              dtskin(i,j) = tskin-td
!             now feedback tskin in surface variable
              tgb(i,j) = tskin
            end if
            tg1d(n,i) = tgb(i,j)
            tgb1d(n,i) = tgb(i,j)
            sent1d(n,i) = sh
            evpr1d(n,i) = lh/wlhv
!           Back out Drag Coefficient
            drag1d(n,i) = ustar**2*rhox2d(i,j)/uv995
            facttq = dlog(z995/2.)/dlog(z995/zo)
            u10m1d(n,i) = ubx3d(i,k,j)*uv10/uv995
            v10m1d(n,i) = vbx3d(i,k,j)*uv10/uv995
            t2m_1d(n,i) = t995 + tzero - dth*facttq
!
            if ( mod(ntime+nint(dtmin*60.),kbats).eq.0 .or.             &
               & (jyear.eq.jyearr .and. ktau.eq.ktaur) ) then
              facttq = dlog(z995/2.)/dlog(z995/zo)
              q2m_1d(n,i) = q995 - dqh*facttq
              tgb2d(n,i,j) = tgb(i,j)
            end if
          end if
        end do
      end do
!
      end subroutine zengocndrv
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
      subroutine zengocn(u,ts,t,q,hgt,zi,ps,qs,u10,tau,                 &
                       & alh,ash,dth,dqh,ustar,zo)
!
      implicit none
!
      real(kind=8) , intent (in) :: hgt , q , t , u , zi , ts , ps
      real(kind=8) , intent (out) :: alh , ash , tau , u10
      real(kind=8) , intent (inout) :: dqh , dth , qs , ustar , zo
!
      real(kind=8) :: dthv , hq , ht , hu , obu , pr , qstar , rb ,     &
               & rho , th , thv , thvstar , tstar , um , visa , zot ,   &
               & wc , xlv , z10 , zbeta , zeta , zetam , zetat , zoq
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
      th = (t+tzero)*(1000./ps)**rovcp
      ! potential T
      dth = t + 0.0098*ht - ts
      qs = qsat(ts,ps)*0.98
      qs = ep2*qs/(ps-0.378*qs)
      ! in kg/kg
      dqh = q - qs
      thv = th*(1.+0.61*q)
      ! virtual potential T
      dthv = dth*(1.+0.61*q) + 0.61*th*dqh
      rho = ps*100./(rgas*(ts+tzero)*(1.+0.61*qs))
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
        call ocnrough(zo,zot,zoq,ustar,visa,gti)
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
      integer , intent(in) :: k
      real(kind=8) , intent(in) :: zeta
      real(kind=8) :: psi
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
      real(kind=8) , intent (in) :: p , t
      real(kind=8) :: qsat
!
      qsat = (1.0007+3.46E-6*p)*6.1121*exp(17.502*t/(240.97+t))
!
      end function qsat

!
!  our formulation for zo,zot,zoq
!
      subroutine ocnrough(zo,zot,zoq,ustar,visa,g)
!
      implicit none
!
      real(kind=8) , intent (in) :: g , ustar , visa
      real(kind=8) , intent (out) :: zoq , zot
      real(kind=8) , intent (inout) :: zo
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
      end subroutine ocnrough
!
      end module mod_zengocn
