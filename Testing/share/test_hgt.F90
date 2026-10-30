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

! Unit test of Share/mod_hgt.F90 on a standard atmosphere (constant lapse
! rate lrate from 288.15 K and 101325 Pa at sea level) above a surface at
! 1500 m: geopotential height on pressure levels, height of model levels,
! pressure of height levels, sea level pressure and hydrostatic heights on
! sigma levels.

program test_hgt
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_hgt
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nk = 8, np = 7, nf = 80, ns = 18
  real(rk8), parameter :: t0 = 288.15_rk8, p0 = 101325.0_rk8
  real(rkx), dimension(nk) :: zl, tl, pl
  real(rkx), dimension(1,1,nk) :: h3, t3, p3
  real(rkx), dimension(1,1,np) :: hp
  real(rkx), dimension(np) :: pp
  real(rkx), dimension(1,1,nf) :: hf, tf, pf
  real(rkx), dimension(1,1) :: ps, ht
  real(rk4), dimension(1,1,nk) :: t4
  real(rk4), dimension(1,1) :: ps4, ht4, slp4, pl4
  real(rkx), dimension(ns) :: sig, zex
  real(rkx), dimension(1,1,ns) :: hs, ts
  real(rkx) :: ptop, p
  real(rk8), dimension(4), parameter :: zsfc = &
    [500.0_rk8, 1500.0_rk8, 3000.0_rk8, 5000.0_rk8]
  real(rk8), dimension(4), parameter :: slptol = &
    [1.0_rk8, 1.0_rk8, 10.0_rk8, 50.0_rk8]
  real(rk8) :: err
  integer(ik4) :: k, n
  character(len=80) :: msg

  ! Model levels from the top, 30 m to 10 km above the surface.
  ht = 1500.0_rkx
  ps = pz(1500.0_rk8)
  zl = [10000.0_rkx, 6000.0_rkx, 3000.0_rkx, 1500.0_rkx, &
        700.0_rkx, 300.0_rkx, 100.0_rkx, 30.0_rkx] + 1500.0_rkx
  do k = 1, nk
    tl(k) = tz(real(zl(k),rk8))
    pl(k) = pz(real(zl(k),rk8))
  end do
  ! Pressure levels inside the model column.
  pp = [80000.0_rkx, 75000.0_rkx, 70000.0_rkx, 50000.0_rkx, &
        40000.0_rkx, 30000.0_rkx, 25000.0_rkx]

  h3(1,1,:) = zl
  t3(1,1,:) = tl
  p3(1,1,:) = pl
  call height(hp,h3,t3,ps,p3,ht,1,1,nk,pp,np)
  err = 0.0_rk8
  do n = 1, np
    err = max(err,abs(hp(1,1,n)-zp(real(pp(n),rk8))))
  end do
  call check(err < 5.0_rk8,'height on pressure levels, top down (5 m)')
  h3(1,1,:) = zl(nk:1:-1)
  t3(1,1,:) = tl(nk:1:-1)
  p3(1,1,:) = pl(nk:1:-1)
  call height(hp,h3,t3,ps,p3,ht,1,1,nk,pp,np)
  err = 0.0_rk8
  do n = 1, np
    err = max(err,abs(hp(1,1,n)-zp(real(pp(n),rk8))))
  end do
  call check(err < 5.0_rk8,'height on pressure levels, bottom up (5 m)')

  h3 = 0.0_rkx
  t3(1,1,:) = tl
  p3(1,1,:) = pl
  call htsig(t3,h3,p3,ps,ht,1,1,nk)
  call check(maxval(abs(h3(1,1,:)-zl)) < 10.0_rkx, &
    'height of model levels (10 m)')
  call check(maxval(abs(h3(1,1,4:nk)-zl(4:nk))) < 0.05_rkx, &
    'height of the lowest model levels (5 cm)')

  ! Height levels every 100 m from the bottom.
  do k = 1, nf
    hf(1,1,k) = 1500.0_rkx + 100.0_rkx*real(k,rkx)
    tf(1,1,k) = tz(real(hf(1,1,k),rk8))
  end do
  ! The first layer uses the temperature of its upper level only: a
  ! relative error of about 1.4e-5 for 100 m.
  call psig(tf,hf,pf,ps,ht,1,1,nf)
  err = 0.0_rk8
  do k = 1, nf
    err = max(err,abs(pf(1,1,k)/pz(real(hf(1,1,k),rk8))-1.0_rk8))
  end do
  call check(err < 2.0e-5_rk8,'psig, pressure of height levels')
  call psig1(tf,hf,pf,ps,ht,1,1,nf)
  err = 0.0_rk8
  do k = 1, nf
    err = max(err,abs(pf(1,1,k)/pz(real(hf(1,1,k),rk8))-1.0_rk8))
  end do
  call check(err < 2.0e-3_rk8,'psig1, pressure of height levels')

  ! Sea level pressure from the lowest level, 30 m above the surface.
  do n = 1, 4
    ht4 = real(zsfc(n),rk4)
    ps4 = real(pz(zsfc(n)),rk4)
    pl4 = real(pz(zsfc(n)+30.0_rk8),rk4)
    t4(1,1,nk) = real(tz(zsfc(n)+30.0_rk8),rk4)
    call mslp(t4,pl4,ps4,ht4,slp4,1,1,nk)
    write(msg,'(a,i0,a)') 'sea level pressure, surface at ', &
      nint(zsfc(n)), ' m'
    call check(abs(real(slp4(1,1),rk8)-p0) < slptol(n),trim(msg))
  end do

  ! Hydrostatic heights of sigma levels, with p* = ps - ptop.
  ptop = 5000.0_rkx
  sig = [0.05_rkx, 0.15_rkx, 0.25_rkx, 0.35_rkx, 0.45_rkx, 0.55_rkx, &
         0.65_rkx, 0.73_rkx, 0.79_rkx, 0.84_rkx, 0.88_rkx, 0.915_rkx, &
         0.945_rkx, 0.965_rkx, 0.98_rkx, 0.99_rkx, 0.995_rkx, 0.9985_rkx]
  ht = 500.0_rkx
  ps = pz(500.0_rk8)
  do k = 1, ns
    p = sig(k)*(ps(1,1)-ptop) + ptop
    zex(k) = real(zp(real(p,rk8)),rkx)
    ts(1,1,k) = tz(real(zex(k),rk8))
  end do
  ps = ps - ptop
  call hydrost(hs,ts,ht,ps,ptop,sig,1,1,ns)
  call check(maxval(abs(hs(1,1,:)-zex)) < 5.0_rkx, &
    'hydrostatic height of sigma levels (5 m)')

  call test_end('test_hgt')

  contains

  real(rkx) function tz(z)
    implicit none
    real(rk8), intent(in) :: z
    tz = real(t0 - real(lrate,rk8)*z,rkx)
  end function tz

  real(rkx) function pz(z)
    implicit none
    real(rk8), intent(in) :: z
    pz = real(p0*(real(tz(z),rk8)/t0)**(egrav/(rgas*lrate)),rkx)
  end function pz

  real(rk8) function zp(p)
    implicit none
    real(rk8), intent(in) :: p
    zp = (t0/lrate)*(1.0_rk8-(p/p0)**(rgas*lrate/egrav))
  end function zp

end program test_hgt
