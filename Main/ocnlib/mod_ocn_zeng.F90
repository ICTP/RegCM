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

module mod_ocn_zeng
  !
  ! Ocean flux model
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_ocn_internal
  use mod_runparams , only : iocnrough , iocnzoq , ktau
  use mod_runparams , only : iocncpl, ntcpl
  use mod_runparams , only : iwavcpl, zomax, ustarmax

  implicit none

  private

  public :: zengocndrv

  ! Module Constants

  real(rk8) , parameter :: minz = 1.0D-6
  real(rk8) , parameter :: a1 = 0.28D+00
  real(rk8) , parameter :: a2 = 0.27D+00
  real(rk8) , parameter :: a3 = 0.45D+00
  real(rk8) , parameter :: b1 = 71.5D+00
  real(rk8) , parameter :: b2 = 2.8D+00
  real(rk8) , parameter :: b3 = 0.07D+00
  real(rk8) , parameter :: alphaw = 0.207D-06
  real(rk8) , parameter :: nuw = 1.004D-06
  real(rk8) , parameter :: kw = 0.60D0
  real(rk8) , parameter :: nu = 0.3D0
  real(rk8) , parameter :: d = 3.0D0 ! reference depth for bulk SST

  ! nu / thermal diffusivity
  real(rk8) , parameter :: pr = 0.71D0   ! Prandtl number

  real(rk8) , parameter :: z10 = d_10    ! m  (reference height)
  real(rk8) , parameter :: zbeta = d_one ! -  (in computing W_*)

  real(rk8) , parameter :: zetat = 0.465D0
  real(rk8) , parameter :: zetam = 1.574D0

  real(rk8) , parameter :: threedays = 86400.0D0*3.0D0  ! 3 days

  real(rk8) , parameter :: missing_r8 = 1.0D20
  real(rk8) , parameter :: tol = missing_r8/2.0D0
  logical :: flag1 , flag2

  contains
  !
  ! Implement Zeng and Beljaars, GRL , 2005, ZB2005
  ! Account for SST diurnal evoluation warm layer/ skin temperature scheme
  !
  subroutine zengocndrv
    implicit none
    real(rk8) :: dqh , dth , facttq , lh , q995 , qs , sh , zo , &
                 t995 , tau , tsurf , ustar , uv10 , uv995 , z995 , zi
    real(rk8) :: dthv , hq , zh , hu , obu , qstar , rb , xdens , &
                 th , thv , thvstar , tstar , um , visa , zot ,   &
                 xlv , wc , zeta , zoq , wt1 , wt2
    integer(ik4) :: i , nconv
!   real(rk8) :: lwds , lwus
    real(rk8) :: rs , rd , td , tdelta , delta
    real(rk8) :: q , ustarw , fd , l , phidl , aa , bb , lamb
    real(rk8) :: dtstend , dts , fs , tskin_new
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'zengocndrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    wt1 = (threedays-dtocn)/threedays
    wt2 = dtocn/threedays

    do i = iocnbeg , iocnend
      if ( mask(i) /= 1 ) cycle

      ! Get temperature from SST
      if ( ldcsst ) then
        tgrd(i) = tskin(i)
      else
        tgrd(i) = tgb(i)
      end if
      tgbrd(i) = tgb(i)

      uv995 = dsqrt(usw(i)**2+vsw(i)**2)
      tsurf = tgrd(i) - tzero
      t995 = sts(i) - tzero
      q995 = qv(i)
      z995 = ht(i)
      zi = max(d_two,hpbl(i))
      hu = z995
      zh = z995
      hq = z995
      ! potential T
      th = sts(i)*(1.0D5/sfps(i))**rovcp
      dth = tatm(i) - tgrd(i)
      qs = pfqsat(tgrd(i),sfps(i))*0.98D0
      ! in kg/kg
      dqh = q995 - qs
      thv = th*(d_one+0.61D0*q995)
      ! virtual potential T
      dthv = dth*(d_one+0.61D0*q995) + 0.61D0*th*dqh
      xdens = sfps(i)/(rgas*tgrd(i)*(d_one+0.61D0*qs))
      ! density
      xlv = (2.501D0-0.00237D0*tsurf)*1.0D+6
      ! J/kg
      !
      ! Kinematic viscosity of dry air (m2/s)
      !   Andreas (1989) CRREL Rep. 89-11
      !
      visa = 1.326D-5*(d_one+6.542D-3*t995 + 8.301D-6*t995*t995 - &
                       4.84D-9*t995*t995*t995)
      !
      ! initial values of u* and convective velocity
      !
      ustar = 0.06D0
      wc = d_half
      if ( dthv >= d_zero ) then
        um = max(uv995,0.1D0)
      else
        um = dsqrt(uv995*uv995+wc*wc)
      end if
      !
      ! zo comes from wave model
      ! flag1 is used as mask for zo
      !
      flag1 = .true.
      if ( iwavcpl == 1 ) then
        if ( zoo(i) < tol .and. ktau+1 > ntcpl ) then
          zo = zoo(i)
          if ( zo > zomax ) zo = zomax
          flag1 = .false.
        end if
      end if
      !
      ! ustr comes from wave model
      ! flag2 is used as mask for ustr
      !
      flag2 = .true.
      if ( iwavcpl == 1 ) then
        if ( ustr(i) < tol .and. ktau+1 > ntcpl ) then
          ustar = ustr(i)
          if ( ustar > ustarmax ) ustar = ustarmax
          flag2 = .false.
        end if
      end if
      !
      ! loop to obtain initial and good ustar and zo
      !
      do nconv = 1 , 5
        call ocnrough(zo,zot,zoq,ustar,um10(i),wc,visa)
        if ( flag2 .or. ktau+1 <= ntcpl ) then
          ustar = vonkar*um/dlog(hu/zo)
        end if
      end do
      rb = egrav*hu*dthv/(thv*um*um)
      if ( rb >= d_zero ) then       ! neutral or stable
        zeta = rb*log(hu/zo)/(d_one-d_five*min(rb,0.19D0))
        zeta = min(d_two,max(zeta,minz))
      else                           ! unstable
        zeta = rb*log(hu/zo)
        zeta = max(-d_100,min(zeta,-minz))
      end if
      obu = hu/zeta
      wc = ustar * (max(-zi*vonkar/obu,d_zero))**onet
      !
      ! main iterations (2-10 iterations would be fine)
      !
      do nconv = 1 , 10
        call ocnrough(zo,zot,zoq,ustar,um10(i),wc,visa)
        !
        ! wind
        !
        if ( flag2 .or. ktau+1 <= ntcpl ) then
          zeta = hu/obu
          if ( zeta < -zetam ) then      ! zeta < -1
            ustar = vonkar*um/(log(-zetam*obu/zo)-psi(1,-zetam)+ &
                    psi(1,zo/obu)+1.14D0*((-zeta)**onet-(zetam)**onet))
          else if ( zeta < d_zero ) then ! -1 <= zeta < 0
            ustar = vonkar*um/(log(hu/zo) - psi(1,zeta)+psi(1,zo/obu))
          else if ( zeta <= d_one ) then !  0 <= zeta <= 1
            ustar = vonkar*um/(log(hu/zo) + d_five*zeta-d_five*zo/obu)
          else                           !  1 < zeta, phi=5+zeta
            ustar = vonkar*um/(log(obu/zo)+d_five-d_five*zo/obu+  &
                    (d_five*log(zeta)+zeta-d_one))
          end if
        end if
        !
        ! temperature
        !
        zeta = zh/obu
        if ( zeta < -zetat ) then      ! zeta < -1
          tstar = vonkar*dth/ &
                  (log(-zetat*obu/zot)-psi(2,-zetat)+psi(2,zot/obu)+ &
                    0.8D0*((zetat)**(-onet)-(-zeta)**(-onet)))
        else if ( zeta < d_zero ) then ! -1 <= zeta < 0
          tstar = vonkar*dth/(log(zh/zot) - psi(2,zeta)+psi(2,zot/obu))
        else if ( zeta <= d_one ) then !  0 <= ztea <= 1
          tstar = vonkar*dth/(log(zh/zot) + d_five*zeta-d_five*zot/obu)
        else                           !  1 < zeta, phi=5+zeta
          tstar = vonkar*dth/(log(obu/zot) + d_five-d_five*zot/obu+ &
                  (d_five*log(zeta)+zeta-d_one))
        end if
        !
        ! humidity
        !
        zeta = hq/obu
        if ( zeta < -zetat ) then      ! zeta < -1
          qstar = vonkar*dqh/ &
                 (log(-zetat*obu/zoq)-psi(2,-zetat)+psi(2,zoq/obu)+ &
                       0.8D0*((zetat)**(-onet)-(-zeta)**(-onet)))
        else if ( zeta < d_zero ) then ! -1 <= zeta < 0
          qstar = vonkar*dqh/(log(hq/zoq) - psi(2,zeta)+psi(2,zoq/obu))
        else if ( zeta <= d_one ) then !  0 <= ztea <= 1
          qstar = vonkar*dqh/(log(hq/zoq) + d_five*zeta-d_five*zoq/obu)
        else                           !  1 < zeta, phi=5+zeta
          qstar = vonkar*dqh/(log(obu/zoq) + d_five-d_five*zoq/obu+ &
                  (d_five*log(zeta)+zeta-d_one))
        end if
        thvstar = tstar*(d_one+0.61D0*q995) + 0.61D0*th*qstar
        zeta = vonkar*egrav*thvstar*hu/(ustar**2*thv)
        if ( zeta >= d_zero ) then   !neutral or stable
          um = max(uv995,0.1D0)
          zeta = min(d_two,max(zeta,minz))
        else                   !unstable
          wc = zbeta*(-egrav*ustar*thvstar*zi/thv)**onet
          um = dsqrt(uv995*uv995+wc*wc)
          zeta = max(-d_100,min(zeta,-minz))
        end if
        obu = hu/zeta
      end do
      tau = xdens*ustar*ustar*uv995/um
      lh = -xdens*xlv*qstar*ustar
      sh = -xdens*cpd*tstar*ustar
      !
      ! x and y components of tau:
      ! lms%taux=xdens*ustar*ustar*u_x/um
      ! lms%tauy=xdens*ustar*ustar*u_y/um
      ! 10-meter wind (without w_* part)
      !
      zeta = z10/obu
      if ( zeta < d_zero ) then
        uv10 = uv995 + (ustar/vonkar)*(log(z10/hu)- &
                        (psi(1,zeta)-psi(1,hu/obu)))
      else
        uv10 = uv995 + (ustar/vonkar)* &
                       (log(z10/hu)+d_five*zeta-d_five*hu/obu)
      end if
      if ( ldcsst ) then
        ! time step considered for the integration of prognostic skin
        ! temperature , equal to BATS time step
        ! Init local variables
        delta = deltas(i)
        tdelta = tdeltas(i)
        ! td is now the 3m bulk SST from the forcing variable
        td = sst(i)
        !
        ! deep impact of aod on sst
        ! if ( sum(aerext(i)) <= 1 ) then
        !   td = sst(i) - sum(aerext(i))*0.8D0
        ! else if ( sum(aerext(i)) > 1 ) then
        !   td = sst(i)- d_one*0.8D0
        ! end if
        !
        ! rs is the net surface sw flux (sw energy absorbed)
        rs = rswf(i)
        ! rd is sw flux at 3m
        rd = rs*(a1*exp(-d*b1) + a2*exp(-d*b2) + a3*exp(-d*b3))
        ! ustar water (with air density == 1)
        ustarw = d_half*ustar*dsqrt(rhox(i)/rhoh2o)
        ! lwds =  dwrlwf(i)
        ! lwus =  emsw*sigm*(tsurf+273.16)**4
        ! q is the skin cooling term inckude net lw flux from
        ! the radiative scheme
        ! q = -(lh+sh+(lwus-lwds))
        q = -(lh+sh+rlwf(i))
        ! fraction of solar radiation abosrbed in the sublayer
        fs = 0.065D0+11.0D0*delta-(6.6D-5/delta) * &
             (d_one-exp(-delta/8.0D-4))
        ! dts= temperature difference between bulk level and skin level
        ! determined from previous time step (via tdelta and td)
        dts = tdelta-td
        ! m.o lenght calculation
        if ( dts > d_zero ) then
          fd = dsqrt(nu*egrav*alphaw/(d_five*d))*    &
                     rhoh2o*cpw0*ustarw**2*dsqrt(dts)
        else
          fd = egrav*alphaw*(q+rs-rd)
        end if
        l = rhoh2o*cpw0*ustarw**3/(vonkar*fd)
        ! calulation of phidl (stability function)
        if ( (d/l) >= d_zero ) then
          phidl = d_one+d_five*(d/l)
        else
          phidl = (d_one-16.0D0*(d/l))**(-d_half)
        end if
        ! prognostic evolution of dts
        ! we can split the tendencies ddts/dt = a - b * dts
        ! with a and b are ultimately function of dts through q
        aa = (q + rs - rd) / (d * cpw0 * rhoh2o * nu/(nu+d_one))
        bb = (nu+d_one) * vonkar * ustarw / (d*phidl)
        ! exponential solution
        dtstend = aa - dts*(d_one-exp(-bb*dtsst))/dtsst
        ! update dts
        dts = dts + dtstend * dtsst
        ! update tdelta
        tdelta = dts + td
        ! update delta thickness and cool skin tempearture
        aa = -16.0D0*egrav*alphaw*rhoh2o*cpw0*nuw**3/(ustarw**4 * kw**2)
        bb =  aa *(q+rs*fs)
        if ( bb > d_zero ) then
          ! case of cool skin layer correction
          lamb = 6.0D0*((d_one+(aa*(q+rs*fs))**0.75D0)**(-onet))
          delta = lamb*nuw/ustarw
          tskin_new = delta/(rhoh2o*cpw0*kw)*(q+rs*fs) + tdelta
        else
          ! no cool skin layer in this case, tskin_new = warm layer temp
          tskin_new = tdelta
        end if
        ! save the temperature difference and skin layer thickness
        ! for next time step
        deltas(i)   = delta
        tdeltas(i)  = tdelta
        tskin(i)   = tskin_new
        ! now feedback tskin in surface variables
        tgb(i) = tskin_new
        tgrd(i) = tskin_new
        tgbrd(i) = sst(i)
      end if ! dcsst

      sent(i)  = sh
      evpr(i)  = lh*rwlhv
      ! Back out Drag Coefficient
      facttq = log(z995*d_half)/log(z995/zo)
      drag(i) = ustar**2*rhox(i)/uv995
      ustr(i) = ustar
      zoo(i)  = zo
      rhoa(i) = xdens
      u10m(i) = usw(i)*uv10/uv995
      v10m(i) = vsw(i)*uv10/uv995
      taux(i) = tau*(usw(i)/uv995)
      tauy(i) = tau*(vsw(i)/uv995)
      t2m(i)  = t995 + tzero - dth*facttq
      q2m(i)  = q995 - dqh*facttq
      ! We need specific humidity in output
      q2m(i) = q2m(i)/(d_one+q2m(i))
      um10(i) = um10(i) * wt1 + sqrt(u10m(i)**2+v10m(i)**2) * wt2
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains
    !
    ! stability function for rb < 0
    !
    pure real(rk8) function psi(k,zeta)
      implicit none
      integer(ik4) , intent(in) :: k
      real(rk8) , intent(in) :: zeta
      real(rk8) :: chik
      chik = (d_one-16.0D0*zeta)**d_rfour
      if ( k == 1 ) then
        psi = d_two*log((d_one+chik)*d_half) +       &
                    log((d_one+chik*chik)*d_half) -  &
              d_two*atan(chik) + d_two*atan(d_one)
      else
        psi = d_two*log((d_one+chik*chik)*d_half)
      end if
    end function psi
    !
    ! our formulation for zo,zot,zoq
    !
    subroutine ocnrough(zo,zot,zoq,ustar,um10,wc,visa)
      implicit none
      real(rk8) , intent (in) :: ustar , um10 , wc , visa
      real(rk8) , intent (out) :: zo , zoq , zot
      real(rk8) :: cp , charnockog , re , xtq , rt , rq , alph
      ! if surface roughness not provided by wave model
      if ( flag1 .or. ktau+1 <= ntcpl ) then
        ! Wave age. The wind here is the mean last N days wind
        cp = 1.2D0*um10
        ! Smith et al. (1992), Carlsson et al. (2009)
        ! Charnock parameter as power function of the wave age
        ! We consider here dominant wind sea waves
        ! Swell dominated sea would require a wave model...
        charnockog = regrav*0.063D0*(cp/ustar)**(-0.4D0)
        if ( iocnrough == 1 ) then
          zo = 0.0065D0*regrav*ustar*ustar
        else if ( iocnrough == 2 ) then
          zo = 0.013D0*regrav*ustar*ustar + 0.11D0*visa/ustar
        else if ( iocnrough == 3 ) then
          zo = 0.017D0*regrav*ustar*ustar
        else if ( iocnrough == 4 ) then
          ! C.H. Huang, 2012
          ! Modification of the Charnock Wind Stress Formula
          ! to Include the Effects of Free Convection and Swell
          ! Advanced Methods for Practical Applications in Fluid Mechanics
          zo = charnockog*(ustar*ustar*ustar+0.11D0*wc*wc*wc)**twot
        else if ( iocnrough == 5 ) then
          if ( um10 < 10.0D0 ) then
            alph = 0.011D0
          else if ( um10 > 10.0D0 .and. um10 < 18.0D0 ) then
            alph = 0.011D0 + 0.000875*(um10-10.0D0)
          else if ( um10 > 18.0D0 .and. um10 < 25.0D0 ) then
            alph = 0.018D0
          else
            alph = max(2.D-3,0.018D0/(d_one+0.050D0*(ustar-0.02D0)**2 - &
                                            0.018D0*(ustar-0.02D0)**1.6D0))
          end if
          zo = alph*regrav*ustar*ustar + 0.11D0*visa/ustar
        else
          zo = charnockog*ustar*ustar
        end if
      end if
      re = (ustar*zo)/visa
      if ( iocnzoq == 2 ) then
        zoq = min(4.0D-4, 2.0D-4*re**(-3.3D0))
        zot = zoq
      else if ( iocnzoq == 3 ) then
        if ( re <= 0.11D0 ) then
          rt = 0.177D0
          rq = 0.292D0
        else if ( re <= 0.8D0 ) then
          rt = 1.376D0*re**0.929D0
          rq = 1.808D0*re**0.826D0
        else if ( re <= 3.0D0 ) then
          rt = 1.026D0*re**(-0.599D0)
          rq = 1.393D0*re**(-0.528D0)
        else if ( re <= 10.0D0 ) then
          rt = 1.625D0*re**(-1.018D0)
          rq = 1.956D0*re**(-0.870D0)
        else if ( re <= 30.0D0 ) then
          rt = 4.661D0*re**(-1.475D0)
          rq = 4.994D0*re**(-1.297D0)
        else if (re <= 100.0D0) then
          rt = 34.904D0*re**(-2.067D0)
          rq = 30.709D0*re**(-1.845D0)
        else if ( re <= 300.0D0 ) then
          rt = 1667.19D0*re**(-2.907D0)
          rq = 1448.68D0*re**(-2.682D0)
        else if ( re <= 1000.0D0 ) then
          rt = 5.88d5*re**(-3.935D0)
          rq = 2.98d5*re**(-3.616D0)
        end if
        zot = rt*visa/ustar
        zoq = rq*visa/ustar
      else
        xtq = 2.67D0*(re**d_rfour) - 2.57D0
        zoq = zo/exp(xtq)
        zot = zoq
      end if
     end subroutine ocnrough

  end subroutine zengocndrv

end module mod_ocn_zeng

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
