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

! Unit test of Share/mod_heatindex.F90: the heat index against the NWS
! heat index chart, the cold limit where it equals the air temperature,
! monotonicity in temperature and humidity, and continuity across the
! regime boundaries.

program test_heatindex
  use mod_intkinds
  use mod_realkinds
  use mod_heatindex
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nchart = 5
  ! NWS chart: air temperature (F), relative humidity, heat index (F)
  real(rkx), dimension(nchart), parameter :: tchart = &
    [ 80.0_rkx, 90.0_rkx, 100.0_rkx, 104.0_rkx, 86.0_rkx ]
  real(rkx), dimension(nchart), parameter :: rchart = &
    [ 0.40_rkx, 0.50_rkx, 0.40_rkx, 0.40_rkx, 0.90_rkx ]
  real(rkx), dimension(nchart), parameter :: hchart = &
    [ 80.0_rkx, 95.0_rkx, 109.0_rkx, 119.0_rkx, 112.0_rkx ]
  real(rkx) :: t, rh, h, hp, dmax, emax
  logical :: mono
  integer(ik4) :: i, j

  call check(heatindex(0.0_rkx,0.5_rkx) == 0.0_rkx, &
    'zero for zero temperature')

  emax = 0.0_rkx
  do i = 1, nchart
    t = (tchart(i)-32.0_rkx)/1.8_rkx + 273.15_rkx
    h = (heatindex(t,rchart(i))-273.15_rkx)*1.8_rkx + 32.0_rkx
    emax = max(emax,abs(h-hchart(i)))
  end do
  call check(emax < 2.0_rkx,'NWS chart within 2 F')

  ! Cold air: the body regulates its skin blood flow alone.
  emax = 0.0_rkx
  do j = 0, 10
    rh = 0.1_rkx*j
    do i = 0, 15
      t = 220.0_rkx + real(i,rkx)
      emax = max(emax,abs(heatindex(t,rh)-t))
    end do
  end do
  call check(emax < 0.05_rkx,'equal to the air temperature below 235 K')

  ! Monotonic and continuous in temperature.
  mono = .true.
  dmax = 0.0_rkx
  do j = 0, 10
    rh = 0.1_rkx*j
    hp = heatindex(230.0_rkx,rh)
    do i = 1, 2000
      t = 230.0_rkx + 0.05_rkx*real(i,rkx)
      h = heatindex(t,rh)
      mono = mono .and. h > hp
      dmax = max(dmax,h-hp)
      hp = h
    end do
  end do
  call check(mono,'increasing with temperature')
  call check(dmax < 1.0_rkx,'continuous in temperature')

  ! Non decreasing and continuous in humidity.
  mono = .true.
  dmax = 0.0_rkx
  do i = 0, 20
    t = 230.0_rkx + 5.0_rkx*real(i,rkx)
    hp = heatindex(t,0.0_rkx)
    do j = 1, 200
      rh = 0.005_rkx*real(j,rkx)
      h = heatindex(t,rh)
      mono = mono .and. h > hp - 1.0e-6_rkx
      dmax = max(dmax,h-hp)
      hp = h
    end do
  end do
  call check(mono,'non decreasing with humidity')
  call check(dmax < 1.0_rkx,'continuous in humidity')

  call check(heatindex(305.0_rkx,1.0_rkx) > &
             heatindex(305.0_rkx,0.0_rkx) + 20.0_rkx, &
    'humid air feels hotter than dry air at 305 K')

  call test_end('test_heatindex')
end program test_heatindex
