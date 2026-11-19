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

! Unit test of Share/saturation.inc: vapour pressure over liquid water and
! ice against Murphy and Koop (2005, QJRMS 131, 1539), the mixed phase
! blend, the derivatives against finite differences, the humidities, the
! liquid/ice switch at 0 C, the dew point inverse and the pressure cap.

program test_saturation
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rk8, rkx
  use mod_constants, only : ep2, rk8, rkx, rtice, rtwat, rtwat_rtice_r, rwat, &
    tzero, wlhs, wlhv
  use share_test_utils, only : check, test_end
  implicit none
  real(rkx) :: t, p, err, w, d, fd
  real(rkx), parameter :: h = 1.0e-3_rkx
  integer(ik4) :: i

  ! Triple point
  call check(abs(esatliq(273.16_rkx)/611.655_rkx-1.0_rkx) < 1.0e-12_rkx, &
    'esatliq at the triple point')
  call check(abs(esatice(273.16_rkx)/611.655_rkx-1.0_rkx) < 1.0e-12_rkx, &
    'esatice at the triple point')

  ! Accuracy against Murphy and Koop (2005). With a constant difference
  ! of heat capacities the error grows away from the triple point: 0.23 %
  ! over liquid from -20 to 40 C (0.53 % at -40 C), 0.63 % over ice down
  ! to 220 K, 3.2 % at 190 K.
  err = 0.0_rkx
  do i = 0, 60
    t = 253.15_rkx + real(i,rkx)
    err = max(err,abs(esatliq(t)/mkliq(t)-1.0_rkx))
  end do
  call check(err < 2.5e-3_rkx, 'esatliq vs Murphy-Koop, -20 to 40 C')
  err = 0.0_rkx
  do i = 0, 80
    t = 233.15_rkx + real(i,rkx)
    err = max(err,abs(esatliq(t)/mkliq(t)-1.0_rkx))
  end do
  call check(err < 6.0e-3_rkx, 'esatliq vs Murphy-Koop, -40 to 40 C')
  err = 0.0_rkx
  do i = 0, 53
    t = 220.0_rkx + real(i,rkx)
    err = max(err,abs(esatice(t)/mkice(t)-1.0_rkx))
  end do
  call check(err < 7.0e-3_rkx, 'esatice vs Murphy-Koop, 220 K to 0 C')
  err = 0.0_rkx
  do i = 0, 83
    t = 190.0_rkx + real(i,rkx)
    err = max(err,abs(esatice(t)/mkice(t)-1.0_rkx))
  end do
  call check(err < 3.5e-2_rkx, 'esatice vs Murphy-Koop, 190 K to 0 C')

  ! Mixed phase: liquid above 0 C, ice below -23 C, monotonic blend
  call check(abs(esatmix(280.0_rkx)-esatliq(280.0_rkx)) < 1.0e-9_rkx, &
    'esatmix is liquid above 0 C')
  call check(abs(esatmix(240.0_rkx)-esatice(240.0_rkx)) < 1.0e-9_rkx, &
    'esatmix is ice below -23 C')
  err = 0.0_rkx
  do i = 1, 230
    t = rtice + 0.1_rkx*real(i,rkx)
    if ( esatmix(t) <= esatmix(t-0.1_rkx) ) err = 1.0_rkx
    if ( esatmix(t) < esatice(t)-1.0e-9_rkx .or. &
         esatmix(t) > esatliq(t)+1.0e-9_rkx ) err = 1.0_rkx
  end do
  call check(err < 0.5_rkx, 'esatmix monotonic between ice and liquid')

  ! Derivatives against centered finite differences
  err = 0.0_rkx
  do i = 0, 130
    t = 190.0_rkx + real(i,rkx)
    fd = (log(esatliq(t+h))-log(esatliq(t-h)))/(2.0_rkx*h)
    err = max(err,abs(dlnesliq(t)/fd-1.0_rkx))
    fd = (log(esatice(t+h))-log(esatice(t-h)))/(2.0_rkx*h)
    err = max(err,abs(dlnesice(t)/fd-1.0_rkx))
  end do
  call check(err < 1.0e-6_rkx, 'dlnesliq, dlnesice vs finite differences')
  ! mixed: derivative of the blend at fixed liquid fraction
  err = 0.0_rkx
  do i = 0, 60
    t = 240.0_rkx + real(i,rkx)
    w = satliqfrac(t)
    fd = (w*(esatliq(t+h)-esatliq(t-h)) + &
          (1.0_rkx-w)*(esatice(t+h)-esatice(t-h)))/(2.0_rkx*h)
    err = max(err,abs(pfesdt(t)/fd-1.0_rkx))
    err = max(err,abs(dlnesmix(t)*esatmix(t)/fd-1.0_rkx))
  end do
  call check(err < 1.0e-6_rkx, 'pfesdt, dlnesmix vs finite differences')
  err = 0.0_rkx
  p = 85000.0_rkx
  do i = 0, 40
    t = 274.0_rkx + real(i,rkx)
    fd = (pfqsat(t+h,p)-pfqsat(t-h,p))/(2.0_rkx*h)
    err = max(err,abs(pfqsdt(t,p)/fd-1.0_rkx))
    fd = (pfqsatl(t+h-40.0_rkx,p)-pfqsatl(t-h-40.0_rkx,p))/(2.0_rkx*h)
    err = max(err,abs(pfqsdtl(t-40.0_rkx,p)/fd-1.0_rkx))
    fd = (pfqsats(t+h-60.0_rkx,p)-pfqsats(t-h-60.0_rkx,p))/(2.0_rkx*h)
    err = max(err,abs(pfqsdts(t-60.0_rkx,p)/fd-1.0_rkx))
  end do
  call check(err < 1.0e-5_rkx, 'pfqsdt, pfqsdtl, pfqsdts vs finite differences')

  ! Humidities: q = w/(1+w), consistent with the vapour pressure
  err = 0.0_rkx
  do i = 0, 100
    t = 220.0_rkx + real(i,rkx)
    p = 30000.0_rkx + 700.0_rkx*real(i,rkx)
    w = pfwsat(t,p)
    err = max(err,abs(pfqsat(t,p)/(w/(1.0_rkx+w))-1.0_rkx))
    err = max(err,abs(w/(ep2*pfesat(t,p)/(p-pfesat(t,p)))-1.0_rkx))
    w = pfwsatl(t,p)
    err = max(err,abs(pfqsatl(t,p)/(w/(1.0_rkx+w))-1.0_rkx))
  end do
  call check(err < 1.0e-12_rkx, 'pfwsat, pfqsat, pfwsatl, pfqsatl')

  ! Liquid/ice switch at 0 C
  call check(abs(pfesats(tzero,1.0e5_rkx)-esatliq(tzero)) < 1.0e-9_rkx .and. &
    abs(pfesats(tzero-1.0e-4_rkx,1.0e5_rkx) - &
        esatice(tzero-1.0e-4_rkx)) < 1.0e-9_rkx, &
    'pfesats switches from liquid to ice at 0 C')

  ! Mixed phase latent heat of pfmxsat
  call pfmxsat(260.0_rkx,60000.0_rkx,w,d)
  call check(abs(w-pfwsat(260.0_rkx,60000.0_rkx)) < 1.0e-15_rkx .and. &
    d > wlhv .and. d < wlhs, 'pfmxsat mixing ratio and latent heat')

  ! Dew point inverse
  err = 0.0_rkx
  do i = 0, 120
    t = 200.0_rkx + real(i,rkx)
    err = max(err,abs(tdewliq(esatliq(t))-t))
  end do
  call check(err < 1.0e-6_rkx, 'tdewliq inverse of esatliq, 200 to 320 K')

  ! Pressure cap
  call check(abs(pfesat(330.0_rkx,10000.0_rkx)-1500.0_rkx) < 1.0e-9_rkx, &
    'vapour pressure capped at 0.15 p')

  call test_end('test_saturation')

  contains

  ! Murphy and Koop (2005), Eq. 10 (liquid) and Eq. 7 (ice), Pa
  pure real(rkx) function mkliq(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rk8) :: tt
    tt = real(t,rk8)
    mkliq = real(exp(54.842763_rk8-6763.22_rk8/tt-4.210_rk8*log(tt) + &
      0.000367_rk8*tt+tanh(0.0415_rk8*(tt-218.8_rk8)) * &
      (53.878_rk8-1331.22_rk8/tt-9.44523_rk8*log(tt)+0.014025_rk8*tt)),rkx)
  end function mkliq

  pure real(rkx) function mkice(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rk8) :: tt
    tt = real(t,rk8)
    mkice = real(exp(9.550426_rk8-5723.265_rk8/tt+3.53068_rk8*log(tt) - &
      0.00728332_rk8*tt),rkx)
  end function mkice

#include <saturation.inc>

end program test_saturation
