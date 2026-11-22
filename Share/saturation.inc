!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !
  ! Saturation of water vapour: the single formulation used by the model
  ! physics. Include this file (contains section of a module) instead of
  ! writing a formula in a scheme.
  !
  ! Saturation vapour pressure over a flat surface of liquid water from
  ! the Clausius-Clapeyron equation with a latent heat linear in
  ! temperature, L(T) = L0 - (c_l - c_pv) (T - T0), integrated from the
  ! triple point (Ambaum, 2020, QJRMS, https://doi.org/10.1002/qj.3899):
  !
  !   ln es = ln es0 + a ln(T0/T) + L0/(Rv T0) - L(T)/(Rv T)
  !         = c0 - a ln T - c1/T
  !
  ! with a = (c_l - c_pv)/Rv, c1 = (L0 + (c_l - c_pv) T0)/Rv and
  ! c0 = ln es0 + a ln T0 + L0/(Rv T0) + a (one exp and one log per call).
  !
  ! Over ice the formula of Murphy and Koop (2005, QJRMS 131, 1539,
  ! Eq. 7, https://doi.org/10.1256/qj.04.94), same cost:
  !
  !   ln es = 9.550426 - 5723.265/T + 3.53068 ln T - 0.00728332 T
  !
  ! Mixed phase: between rtice (-23 C) and rtwat (0 C) the vapour pressure
  ! is the blend es = w es_liq + (1-w) es_ice, with the liquid fraction
  ! w = ((T - rtice)/(rtwat - rtice))**2 used by the cloud microphysics
  ! and convection; liquid above 0 C, ice below -23 C.
  !
  ! The derivatives are exact for the vapour pressures at fixed phase
  ! fraction: d ln es/dT = L(T)/(Rv T**2) over liquid, the derivative of
  ! the Murphy-Koop formula over ice, weighted with the phase
  ! contributions to es for the blend.
  !
  ! The functions with the 'l' suffix are over liquid water at every
  ! temperature (sea water, FAO evapotranspiration). The functions with
  ! the 's' suffix switch from liquid to ice at 0 C: surfaces, whose phase
  ! follows their temperature, and the relative humidity diagnostics
  ! (CMIP convention: over liquid water above 0 C, over ice below).
  !
  ! The vapour pressure used for the humidities is capped at 0.15 p.
  !
  ! esatliq(t), esatice(t), esatmix(t)       vapour pressure (Pa)
  ! dlnesliq(t), dlnesice(t), dlnesmix(t)    d ln es / dT (1/K)
  ! satliqfrac(t)                            liquid fraction w
  ! pfesat(t,p)                              mixed phase es (Pa), capped
  ! pfesdt(t)                                mixed phase d es/dT (Pa/K)
  ! pfwsat(t,p), pfqsat(t,p)                 mixed phase saturation mixing
  !                                          ratio, specific humidity
  ! pfmxsat(t,p,ws,leff)                     mixed phase ws and latent heat
  ! pfesatl(t,p), pfwsatl(t,p), pfqsatl(t,p) the same over liquid water
  ! pfesdtl(t)                               d es/dT over liquid water
  ! pfesats(t,p), pfwsats(t,p), pfqsats(t,p) liquid above 0 C, ice below
  ! pfesdts(t)                               d es/dT, liquid/ice at 0 C
  ! pfqsdt(t,p), pfqsdtl(t,p), pfqsdts(t,p)  d qs/dT (1/K), specific
  !                                          humidity: mixed, liquid,
  !                                          liquid/ice at 0 C
  ! tdewliq(e)                               dew point (K) over liquid
  !                                          water for the vapour pressure
  !                                          e (Pa), inverse of esatliq
  !
!DIR$ ATTRIBUTES FORCEINLINE :: esatliq
  pure elemental real(rkx) function esatliq(t) result(es)
!$acc routine seq
    use mod_constants, only : rwat
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rk8), parameter :: es0 = 611.655_rk8 ! Pa, triple point
    real(rk8), parameter :: t0 = 273.16_rk8
    real(rk8), parameter :: l0 = 2.501e6_rk8  ! J kg-1
    real(rk8), parameter :: dc = 2320.0_rk8   ! c_l - c_pv, J kg-1 K-1
    real(rk8), parameter :: rv = real(rwat,rk8)
    real(rk8), parameter :: a = dc/rv
    real(rk8), parameter :: c1 = (l0+dc*t0)/rv
    real(rk8), parameter :: c0 = log(es0) + a*log(t0) + l0/(rv*t0) + a
    real(rk8) :: tt
    tt = real(t,rk8)
    es = real(exp(c0 - a*log(tt) - c1/tt),rkx)
  end function esatliq

!DIR$ ATTRIBUTES FORCEINLINE :: esatice
  pure elemental real(rkx) function esatice(t) result(es)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    ! Murphy and Koop (2005), Eq. 7, valid above 110 K
    real(rk8), parameter :: c0 = 9.550426_rk8
    real(rk8), parameter :: c1 = 5723.265_rk8
    real(rk8), parameter :: a = 3.53068_rk8
    real(rk8), parameter :: b = 0.00728332_rk8
    real(rk8) :: tt
    tt = real(t,rk8)
    es = real(exp(c0 - c1/tt + a*log(tt) - b*tt),rkx)
  end function esatice

!DIR$ ATTRIBUTES FORCEINLINE :: dlnesliq
  pure elemental real(rkx) function dlnesliq(t) result(d)
!$acc routine seq
    use mod_constants, only : rwat
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    d = (2.501e6_rkx-2320.0_rkx*(t-273.16_rkx))/(rwat*t*t)
  end function dlnesliq

!DIR$ ATTRIBUTES FORCEINLINE :: dlnesice
  pure elemental real(rkx) function dlnesice(t) result(d)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    d = 5723.265_rkx/(t*t) + 3.53068_rkx/t - 0.00728332_rkx
  end function dlnesice

!DIR$ ATTRIBUTES FORCEINLINE :: satliqfrac
  pure elemental real(rkx) function satliqfrac(t) result(w)
!$acc routine seq
    use mod_constants, only : rtice, rtwat, rtwat_rtice_r
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    w = ((max(rtice,min(rtwat,t))-rtice)*rtwat_rtice_r)**2
  end function satliqfrac

!DIR$ ATTRIBUTES FORCEINLINE :: esatmix
  pure elemental real(rkx) function esatmix(t) result(es)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx) :: w
    w = satliqfrac(t)
    es = w*esatliq(t) + (1.0_rkx-w)*esatice(t)
  end function esatmix

!DIR$ ATTRIBUTES FORCEINLINE :: dlnesmix
  pure elemental real(rkx) function dlnesmix(t) result(d)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx) :: w, el, ei
    w = satliqfrac(t)
    el = w*esatliq(t)
    ei = (1.0_rkx-w)*esatice(t)
    d = (el*dlnesliq(t) + ei*dlnesice(t))/(el+ei)
  end function dlnesmix

!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    es = min(esatmix(t),0.15_rkx*p)
  end function pfesat

!DIR$ ATTRIBUTES FORCEINLINE :: pfesdt
  pure elemental real(rkx) function pfesdt(t) result(esdt)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx) :: w, el, ei
    w = satliqfrac(t)
    el = w*esatliq(t)
    ei = (1.0_rkx-w)*esatice(t)
    esdt = el*dlnesliq(t) + ei*dlnesice(t)
  end function pfesdt

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure elemental real(rkx) function pfwsat(t,p) result(ws) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesat(t,p)
    ws = ep2*es/(p-es)
  end function pfwsat

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsat
  pure elemental real(rkx) function pfqsat(t,p) result(qs) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesat(t,p)
    qs = ep2*es/(p-(1.0_rkx-ep2)*es)
  end function pfqsat

!DIR$ ATTRIBUTES FORCEINLINE :: pfmxsat
  pure elemental subroutine pfmxsat(t,p,ws,leff)
!$acc routine seq
    use mod_constants, only : wlhv, wlhs
    implicit none
    real(rkx), intent(in)  :: t    ! Temperature (K)
    real(rkx), intent(in)  :: p    ! Pressure (Pa)
    real(rkx), intent(out) :: ws   ! Saturation mixing ratio (kg/kg)
    real(rkx), intent(out) :: leff ! Effective latent heat (J/kg)
    real(rkx) :: w
    w = satliqfrac(t)
    ws = pfwsat(t,p)
    leff = w*wlhv + (1.0_rkx-w)*wlhs
  end subroutine pfmxsat

!DIR$ ATTRIBUTES FORCEINLINE :: pfesatl
  pure elemental real(rkx) function pfesatl(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    es = min(esatliq(t),0.15_rkx*p)
  end function pfesatl

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsatl
  pure elemental real(rkx) function pfwsatl(t,p) result(ws) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesatl(t,p)
    ws = ep2*es/(p-es)
  end function pfwsatl

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsatl
  pure elemental real(rkx) function pfqsatl(t,p) result(qs) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesatl(t,p)
    qs = ep2*es/(p-(1.0_rkx-ep2)*es)
  end function pfqsatl

!DIR$ ATTRIBUTES FORCEINLINE :: pfesdtl
  pure elemental real(rkx) function pfesdtl(t) result(esdt)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    esdt = esatliq(t)*dlnesliq(t)
  end function pfesdtl

!DIR$ ATTRIBUTES FORCEINLINE :: pfesats
  pure elemental real(rkx) function pfesats(t,p) result(es)
!$acc routine seq
    use mod_constants, only : tzero
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    if ( t >= tzero ) then
      es = min(esatliq(t),0.15_rkx*p)
    else
      es = min(esatice(t),0.15_rkx*p)
    end if
  end function pfesats

!DIR$ ATTRIBUTES FORCEINLINE :: pfesdts
  pure elemental real(rkx) function pfesdts(t) result(esdt)
!$acc routine seq
    use mod_constants, only : tzero
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    if ( t >= tzero ) then
      esdt = esatliq(t)*dlnesliq(t)
    else
      esdt = esatice(t)*dlnesice(t)
    end if
  end function pfesdts

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsats
  pure elemental real(rkx) function pfwsats(t,p) result(ws) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesats(t,p)
    ws = ep2*es/(p-es)
  end function pfwsats

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsats
  pure elemental real(rkx) function pfqsats(t,p) result(qs) ! In kg/kg
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesats(t,p)
    qs = ep2*es/(p-(1.0_rkx-ep2)*es)
  end function pfqsats

!DIR$ ATTRIBUTES FORCEINLINE :: tdewliq
  pure elemental real(rkx) function tdewliq(e) result(td)
!$acc routine seq
    use mod_constants, only : rwat
    implicit none
    real(rkx), intent(in) :: e     ! Vapour pressure (Pa)
    real(rk8), parameter :: es0 = 611.655_rk8
    real(rk8), parameter :: t0 = 273.16_rk8
    real(rk8), parameter :: l0 = 2.501e6_rk8
    real(rk8), parameter :: dc = 2320.0_rk8
    real(rk8), parameter :: rv = real(rwat,rk8)
    real(rk8), parameter :: a = dc/rv
    real(rk8), parameter :: c1 = (l0+dc*t0)/rv
    real(rk8), parameter :: c0 = log(es0) + a*log(t0) + l0/(rv*t0) + a
    real(rk8) :: le, tt, f, df
    integer(ik4) :: it
    le = log(max(real(e,rk8),1.0e-10_rk8))
    ! first guess: inverse of the Bolton (1980) formula
    tt = (4780.8_rk8 - 32.19_rk8*(le-log(613.3_rk8))) / &
         (17.502_rk8 - (le-log(613.3_rk8)))
    ! Newton iterations on ln es(T) - ln e = c0 - a ln T - c1/T - ln e
    do it = 1, 3
      f = c0 - a*log(tt) - c1/tt - le
      df = -a/tt + c1/(tt*tt)
      tt = tt - f/df
    end do
    td = real(tt,rkx)
  end function tdewliq

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsdt
  pure elemental real(rkx) function pfqsdt(t,p) result(qsdt)
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesat(t,p)
    qsdt = ep2*p*pfesdt(t)/(p-(1.0_rkx-ep2)*es)**2
  end function pfqsdt

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsdtl
  pure elemental real(rkx) function pfqsdtl(t,p) result(qsdt)
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesatl(t,p)
    qsdt = ep2*p*pfesdtl(t)/(p-(1.0_rkx-ep2)*es)**2
  end function pfqsdtl

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsdts
  pure elemental real(rkx) function pfqsdts(t,p) result(qsdt)
!$acc routine seq
    use mod_constants, only : ep2
    implicit none
    real(rkx), intent(in) :: t     ! Temperature (K)
    real(rkx), intent(in) :: p     ! Pressure (Pa)
    real(rkx) :: es
    es = pfesats(t,p)
    qsdt = ep2*p*pfesdts(t)/(p-(1.0_rkx-ep2)*es)**2
  end function pfqsdts

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
