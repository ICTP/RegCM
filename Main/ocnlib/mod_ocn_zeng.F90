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
  use mod_constants
  use mod_service
  use mod_ocn_internal
  use mod_runparams , only : iocnrough , iocnzoq , syncro_cpl
  use mod_runparams , only : iocncpl , iwavcpl
  use mod_runparams , only : zomax , ustarmax

  implicit none

  private

  public :: zengocndrv

  ! Module Constants

  real(rkx) , parameter :: minz = 1.0e-6_rkx
  real(rkx) , parameter :: a1 = 0.28_rkx
  real(rkx) , parameter :: a2 = 0.27_rkx
  real(rkx) , parameter :: a3 = 0.45_rkx
  real(rkx) , parameter :: b1 = 71.5_rkx
  real(rkx) , parameter :: b2 = 2.8_rkx
  real(rkx) , parameter :: b3 = 0.07_rkx
  real(rkx) , parameter :: alphaw = 0.207e-06_rkx
  real(rkx) , parameter :: nuw = 1.004e-06_rkx
  real(rkx) , parameter :: kw = 0.60_rkx
  real(rkx) , parameter :: nu = 0.3_rkx
  real(rkx) , parameter :: d = 3.0_rkx ! reference depth for bulk SST

  ! nu / thermal diffusivity
  real(rkx) , parameter :: pr = 0.71_rkx   ! Prandtl number

  real(rkx) , parameter :: z10 = d_10    ! m  (reference height)
  real(rkx) , parameter :: zbeta = d_one ! -  (in computing W_*)

  real(rkx) , parameter :: zetat = 0.465_rkx
  real(rkx) , parameter :: zetam = 1.574_rkx
  real(rkx) , parameter :: minw = 0.1_rkx

  real(rkx) , parameter :: missing_r8 = 1.0e20_rkx
  real(rkx) , parameter :: tol = missing_r8/2.0_rkx
  logical :: flag1 , flag2

  contains
  !
  ! Implement Zeng and Beljaars, GRL , 2005, ZB2005
  ! Account for SST diurnal evoluation warm layer/ skin temperature scheme
  !
  subroutine zengocndrv
    implicit none
    real(rkx) :: dqh , dth , facttq , lh , qs , sh , zo , &
                 tau , tsurf , ustar , uv10 , zi
    real(rkx) :: t995 , q995 , uv995 , z995
    real(rkx) :: dthv , hq , zh , hu , obu , qstar , xdens ,    &
                 th , thv , thvstar , tstar , um , visa , zot , &
                 wc , zeta , zoq , wt1 , wt2 , tha , nobu , rlv
    integer(ik4) :: i , nconv
!   real(rkx) :: lwds , lwus
    real(rkx) :: rs , rd , td , tdelta , delta
    real(rkx) :: q , ustarw , fd , l , phidl , aa , bb , lamb
    real(rkx) :: dtstend , dts , fs , tskin_new
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

      uv995 = max(sqrt(usw(i)**2+vsw(i)**2),minw)
      tsurf = tgrd(i) - tzero
      t995 = tatm(i) - tzero
      q995 = qv(i)
      z995 = ht(i)
      rlv = wlh(tgrd(i))
      zi = max(z995,hpbl(i))
      hu = z995
      zh = z995
      hq = z995
      ! potential T
      th = tgrd(i)*(p00/sfps(i))**rovcp
      tha = tatm(i)*(p00/patm(i))**rovcp
      dth = tha - th
      qs = pfwsat(tgrd(i),sfps(i))*0.98_rkx
      ! in kg/kg
      dqh = q995 - qs
      ! virtual potential T
      thv = th*(d_one+ep1*q995)
      dthv = dth*(d_one+ep1*q995) + ep1*th*dqh
      ! density
      xdens = sfps(i)/(rgas*tatm(i)*(d_one+ep1*q995))
      ! J/kg
      ! Kinematic viscosity of dry air (m2/s)
      !   Andreas (1989) CRREL Rep. 89-11
      !
      visa = 1.326e-5_rkx*(d_one + 6.542e-3_rkx * t995 + &
                                   8.301e-6_rkx * t995*t995 - &
                                   4.840e-9_rkx * t995*t995*t995)
      !
      ! initial values of u* and convective velocity
      !
      ustar = 0.06_rkx
      wc = d_half
      if ( dthv >= d_zero ) then
        um = uv995
      else
        um = sqrt(uv995*uv995+wc*wc)
      end if
      !
      ! zo comes from wave model
      ! flag1 is used as mask for zo
      !
      flag1 = .true.
      if ( iwavcpl == 1 ) then
        if ( syncro_cpl%lcount > 1 ) then
          if ( zoo(i) < tol .and. syncro_cpl%act( ) ) then
            zo = zoo(i)
            if ( zo > zomax ) zo = zomax
            flag1 = .false.
          end if
        end if
      end if
      !
      ! ustr comes from wave model
      ! flag2 is used as mask for ustr
      !
      flag2 = .true.
      if ( iwavcpl == 1 ) then
        if ( syncro_cpl%lcount > 1 ) then
          if ( ustr(i) < tol .and. syncro_cpl%act( ) ) then
            ustar = ustr(i)
            if ( ustar > ustarmax ) ustar = ustarmax
            flag2 = .false.
          end if
        end if
      end if
      !
      ! loop to obtain initial and good ustar and zo
      !
      do nconv = 1 , 2
        call ocnrough(zo,zot,zoq,ustar,um10(i),wc,visa)
        if ( flag2 ) then
          ustar = vonkar*um/log(hu/zo)
        end if
      end do
      br(i) = egrav*hu*dthv/(thv*um*um)
      if ( br(i) >= d_zero ) then       ! neutral or stable
        zeta = br(i)*log(hu/zo)/(d_one-d_five*min(br(i),0.19_rkx))
        zeta = min(d_two,max(zeta,minz))
      else                              ! unstable
        zeta = br(i)*log(hu/zo)
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
        if ( flag2 ) then
          zeta = hu/obu
          if ( zeta < -zetam ) then      ! zeta < -1
            ram1(i) = (log(-zetam*obu/zo)-psi(1,-zetam)+ &
                    psi(1,zo/obu)+1.14_rkx*((-zeta)**onet-(zetam)**onet))
          else if ( zeta < d_zero ) then ! -1 <= zeta < 0
            ram1(i) = (log(hu/zo) - psi(1,zeta)+psi(1,zo/obu))
          else if ( zeta <= d_one ) then !  0 <= zeta <= 1
            ram1(i) = (log(hu/zo) + d_five*zeta-d_five*zo/obu)
          else                           !  1 < zeta, phi=5+zeta
            ram1(i) = (log(obu/zo)+d_five-d_five*zo/obu+  &
                        (d_five*log(zeta)+zeta-d_one))
          end if
          ustar = vonkar*um/ram1(i)
        else
          ram1(i) = vonkar*um/ustar
        end if
        !
        ! temperature
        !
        zeta = zh/obu
        if ( zeta < -zetat ) then      ! zeta < -1
          rah1(i) = (log(-zetat*obu/zot)-psi(2,-zetat)+psi(2,zot/obu)+ &
                    0.8_rkx*((zetat)**(-onet)-(-zeta)**(-onet)))
        else if ( zeta < d_zero ) then ! -1 <= zeta < 0
          rah1(i) = (log(zh/zot) - psi(2,zeta)+psi(2,zot/obu))
        else if ( zeta <= d_one ) then !  0 <= ztea <= 1
          rah1(i) = (log(zh/zot) + d_five*zeta-d_five*zot/obu)
        else                           !  1 < zeta, phi=5+zeta
          rah1(i) = (log(obu/zot) + d_five-d_five*zot/obu+ &
                  (d_five*log(zeta)+zeta-d_one))
        end if
        tstar = vonkar*dth/rah1(i)
        !
        ! humidity
        !
        zeta = hq/obu
        if ( zeta < -zetat ) then      ! zeta < -1
          qstar = vonkar*dqh/ &
                 (log(-zetat*obu/zoq)-psi(2,-zetat)+psi(2,zoq/obu)+ &
                       0.8_rkx*((zetat)**(-onet)-(-zeta)**(-onet)))
        else if ( zeta < d_zero ) then ! -1 <= zeta < 0
          qstar = vonkar*dqh/(log(hq/zoq) - psi(2,zeta)+psi(2,zoq/obu))
        else if ( zeta <= d_one ) then !  0 <= ztea <= 1
          qstar = vonkar*dqh/(log(hq/zoq) + d_five*zeta-d_five*zoq/obu)
        else                           !  1 < zeta, phi=5+zeta
          qstar = vonkar*dqh/(log(obu/zoq) + d_five-d_five*zoq/obu+ &
                  (d_five*log(zeta)+zeta-d_one))
        end if
        thvstar = tstar*(d_one+ep1*q995) + ep1*th*qstar
        zeta = vonkar*egrav*thvstar*hu/(ustar**2*thv)
        if ( zeta >= d_zero ) then     !  neutral or stable
          um = uv995
          zeta = min(d_two,max(zeta,minz))
        else                           ! unstable
          wc = zbeta*(-egrav*ustar*thvstar*zi/thv)**onet
          um = sqrt(uv995*uv995+wc*wc)
          zeta = max(-d_100,min(zeta,-minz))
        end if
        nobu = hu/zeta
        if ( abs(nobu-obu) < 0.1_rkx ) then
          obu = nobu
          exit
        end if
        obu = nobu
      end do
      tau = xdens*ustar*ustar*uv995/um
      lh = -xdens*rlv*qstar*ustar
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
        sst(i) = tgb(i)
        delta = deltas(i)
        tdelta = tdeltas(i)
        ! td is now the 3m bulk SST from the forcing variable
        td = sst(i)
        !
        ! deep impact of aod on sst
        ! if ( sum(aerext(i)) <= 1 ) then
        !   td = sst(i) - sum(aerext(i))*0.8_rkx
        ! else if ( sum(aerext(i)) > 1 ) then
        !   td = sst(i)- d_one*0.8_rkx
        ! end if
        !
        ! rs is the net surface sw flux (sw energy absorbed)
        rs = rswf(i)
        ! rd is sw flux at 3m
        rd = rs*(a1*exp(-d*b1) + a2*exp(-d*b2) + a3*exp(-d*b3))
        ! ustar water (with air density == 1)
        ustarw = d_half*ustar*sqrt(rhox(i)/rhoh2o)
        ! lwds =  dwrlwf(i)
        ! lwus =  emsw*sigm*(tsurf+273.16)**4
        ! q is the skin cooling term inckude net lw flux from
        ! the radiative scheme
        ! q = -(lh+sh+(lwus-lwds))
        q = -(lh+sh+rlwf(i))
        ! fraction of solar radiation abosrbed in the sublayer
        fs = 0.065_rkx+11.0_rkx*delta-(6.6e-5_rkx/delta) * &
             (d_one-exp(-delta/8.0e-4_rkx))
        ! dts= temperature difference between bulk level and skin level
        ! determined from previous time step (via tdelta and td)
        dts = tdelta-td
        ! m.o lenght calculation
        if ( dts > d_zero ) then
          fd = sqrt(nu*egrav*alphaw/(d_five*d))*    &
                     rhoh2o*cpw0*ustarw**2*sqrt(dts)
        else
          fd = egrav*alphaw*(q+rs-rd)
        end if
        if ( fd > d_zero ) then
          l = rhoh2o*cpw0*ustarw**3/(vonkar*fd)
          ! calulation of phidl (stability function)
          if ( (d/l) >= d_zero ) then
            phidl = d_one+d_five*(d/l)
          else
            phidl = (d_one-16.0_rkx*(d/l))**(-d_half)
          end if
        else
           phidl = d_one
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
        aa = -16.0_rkx*egrav*alphaw*rhoh2o*cpw0*nuw**3/(ustarw**4 * kw**2)
        bb =  aa *(q+rs*fs)
        if ( bb > d_zero ) then
          ! case of cool skin layer correction
          lamb = 6.0_rkx*((d_one+(aa*(q+rs*fs))**0.75_rkx)**(-onet))
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

      sent(i) = sh
      evpr(i) = lh/rlv
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
      t2m(i)  = t995 + tzero - (dth * (sfps(i)/p00)**rovcp) * facttq
      q2m(i)  = q995 - dqh*facttq
      ! We need specific humidity in output
      q2m(i) = q2m(i)/(d_one+q2m(i))
      um10(i) = um10(i) * wt1 + sqrt(u10m(i)**2+v10m(i)**2) * wt2
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains

#include <pfesat.inc>
#include <pfwsat.inc>
#include <wlh.inc>
    !
    ! stability function for rb < 0
    !
    pure real(rkx) function psi(k,zeta)
      implicit none
      integer(ik4) , intent(in) :: k
      real(rkx) , intent(in) :: zeta
      real(rkx) :: chik
      chik = (d_one-16.0_rkx*zeta)**d_rfour
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
      real(rkx) , intent (in) :: ustar , um10 , wc , visa
      real(rkx) , intent (out) :: zo , zoq , zot
      real(rkx) :: cp , charnockog , re , xtq , rt , rq , alph
      ! if surface roughness not provided by wave model
      if ( flag1 ) then
        ! Wave age. The wind here is the mean last N days wind
        cp = 1.2_rkx*um10
        ! Smith et al. (1992), Carlsson et al. (2009)
        ! Charnock parameter as power function of the wave age
        ! We consider here dominant wind sea waves
        ! Swell dominated sea would require a wave model...
        charnockog = regrav*0.063_rkx*(cp/ustar)**(-0.4_rkx)
        if ( iocnrough == 1 ) then
          zo = 0.0065_rkx*regrav*ustar*ustar
        else if ( iocnrough == 2 ) then
          zo = 0.013_rkx*regrav*ustar*ustar + 0.11_rkx*visa/ustar
        else if ( iocnrough == 3 ) then
          zo = 0.017_rkx*regrav*ustar*ustar
        else if ( iocnrough == 4 ) then
          ! C.H. Huang, 2012
          ! Modification of the Charnock Wind Stress Formula
          ! to Include the Effects of Free Convection and Swell
          ! Advanced Methods for Practical Applications in Fluid Mechanics
          zo = charnockog*(ustar*ustar*ustar+0.11_rkx*wc*wc*wc)**twot
        else if ( iocnrough == 5 ) then
          if ( um10 < 10.0_rkx ) then
            alph = 0.011_rkx
          else if ( um10 > 10.0_rkx .and. um10 < 18.0_rkx ) then
            alph = 0.011_rkx + 0.000875*(um10-10.0_rkx)
          else if ( um10 > 18.0_rkx .and. um10 < 25.0_rkx ) then
            alph = 0.018_rkx
          else
            alph = max(2.e-3_rkx,0.018_rkx / &
                 (d_one+0.050_rkx*(ustar-0.02_rkx)**2 - &
                  0.018_rkx*(ustar-0.02_rkx)**1.6_rkx))
          end if
          zo = alph*regrav*ustar*ustar + 0.11_rkx*visa/ustar
        else
          zo = charnockog*ustar*ustar
        end if
      end if
      zo = max(zo,1.0e-8_rkx)
      re = (ustar*zo)/visa
      if ( iocnzoq == 2 ) then
        zoq = min(4.0e-4_rkx, 2.0e-4_rkx*re**(-3.3_rkx))
        zot = zoq
      else if ( iocnzoq == 3 ) then
        if ( re <= 0.11_rkx ) then
          rt = 0.177_rkx
          rq = 0.292_rkx
        else if ( re <= 0.8_rkx ) then
          rt = 1.376_rkx*re**0.929_rkx
          rq = 1.808_rkx*re**0.826_rkx
        else if ( re <= 3.0_rkx ) then
          rt = 1.026_rkx*re**(-0.599_rkx)
          rq = 1.393_rkx*re**(-0.528_rkx)
        else if ( re <= 10.0_rkx ) then
          rt = 1.625_rkx*re**(-1.018_rkx)
          rq = 1.956_rkx*re**(-0.870_rkx)
        else if ( re <= 30.0_rkx ) then
          rt = 4.661_rkx*re**(-1.475_rkx)
          rq = 4.994_rkx*re**(-1.297_rkx)
        else if (re <= 100.0_rkx) then
          rt = 34.904_rkx*re**(-2.067_rkx)
          rq = 30.709_rkx*re**(-1.845_rkx)
        else if ( re <= 300.0_rkx ) then
          rt = 1667.19_rkx*re**(-2.907_rkx)
          rq = 1448.68_rkx*re**(-2.682_rkx)
        else if ( re <= 1000.0_rkx ) then
          rt = 5.88e5_rkx*re**(-3.935_rkx)
          rq = 2.98e5_rkx*re**(-3.616_rkx)
        else
          rt = 1e-10_rkx
          rq = 1e-10_rkx
        end if
        zot = rt*visa/ustar
        zoq = rq*visa/ustar
      else
        xtq = 2.67_rkx*(re**d_rfour) - 2.57_rkx
        zoq = zo/exp(xtq)
        zot = zoq
      end if
     end subroutine ocnrough

  end subroutine zengocndrv

end module mod_ocn_zeng

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
