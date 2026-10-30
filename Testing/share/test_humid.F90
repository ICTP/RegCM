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

! Unit test of Share/mod_humid.F90: relative humidity to mixing ratio on
! pressure, sigma and height levels against the same saturation mixing
! ratio at the exact pressure of a standard atmosphere, the inverse
! conversion, and specific humidity to mixing ratio and back.

program test_humid
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_humid
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nk = 8
  real(rk8), parameter :: t0 = 288.15_rk8, p0 = 101325.0_rk8
  real(rkx), dimension(1,1,nk) :: z, t, q, w
  real(rkx), dimension(nk) :: pex, qex, rh, sig
  real(rkx), dimension(1,1) :: ps
  real(rk8), dimension(1,1,nk) :: t8, q8, p8
  real(rkx) :: ptop
  real(rk8) :: ht
  integer(ik4) :: k

  ! Heights above a surface at 500 m, from the top, and relative humidity
  ! from 10 to 100 %.
  ht = 500.0_rk8
  z(1,1,:) = [12000.0_rkx, 8000.0_rkx, 5000.0_rkx, 3000.0_rkx, &
              1500.0_rkx, 500.0_rkx, 150.0_rkx, 30.0_rkx]
  rh = [0.1_rkx, 0.2_rkx, 0.4_rkx, 0.5_rkx, 0.7_rkx, 0.8_rkx, &
        0.9_rkx, 1.0_rkx]
  ps = real(pz(ht)/1000.0_rk8,rkx)
  do k = 1, nk
    t(1,1,k) = real(t0-lrate*(ht+real(z(1,1,k),rk8)),rkx)
    pex(k) = real(pz(ht+real(z(1,1,k),rk8)),rkx)
  end do
  qex = rh
  call rh2mxr(t(1,1,:),qex,pex,nk)

  q(1,1,:) = rh
  call rh2mxr_z(t,q,z,ps,1,1,nk)
  call check(maxval(abs(q(1,1,:)/qex-1.0_rkx)) < 2.0e-3_rkx, &
    'rh2mxr_z, height levels')

  ! Sigma levels with ptop = 5 kPa at the same pressures.
  ptop = 5.0_rkx
  sig = (pex/1000.0_rkx-ptop)/(ps(1,1)-ptop)
  q(1,1,:) = rh
  call rh2mxr(t,q,ps-ptop,ptop,sig,1,1,nk)
  call check(maxval(abs(q(1,1,:)/qex-1.0_rkx)) < 1.0e-6_rkx, &
    'rh2mxr, sigma levels')

  q(1,1,:) = rh
  call rh2mxr(t,q,pex,1,1,nk)
  call check(maxval(abs(q(1,1,:)/qex-1.0_rkx)) < 1.0e-12_rkx, &
    'rh2mxr, pressure levels')

  ! Back to relative humidity.
  t8 = real(t,rk8)
  q8(1,1,:) = real(qex,rk8)
  p8(1,1,:) = real(pex,rk8)
  call mxr2rh(t8,q8,p8,1,1,nk)
  call check(maxval(abs(q8(1,1,:)-real(rh,rk8))) < 1.0e-6_rk8, &
    'mxr2rh is the inverse of rh2mxr')

  ! Specific humidity to mixing ratio and back.
  w(1,1,:) = qex
  call mxr2sph(w,1,1,nk)
  call check(all(w(1,1,:) < qex),'specific humidity below mixing ratio')
  call sph2mxr(w,1,1,nk)
  call check(maxval(abs(w(1,1,:)-qex)) < 1.0e-12_rkx, &
    'sph2mxr is the inverse of mxr2sph')

  call test_end('test_humid')

  contains

  real(rk8) function pz(zz)
    implicit none
    real(rk8), intent(in) :: zz
    pz = p0*((t0-lrate*zz)/t0)**(egrav/(rgas*lrate))
  end function pz

end program test_humid
