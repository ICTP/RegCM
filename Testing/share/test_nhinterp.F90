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

! Unit test of Share/mod_nhinterp.F90 on a column over 3000 m of terrain:
! the reference surface pressure of nhbase is the pressure of the
! reference atmosphere at the terrain height, the reference heights of the
! levels are the hydrostatic heights of the reference temperature, and an
! atmosphere equal to the reference one is interpolated by nhinterp to the
! same values (a field linear in height) with no pressure perturbation
! (nhpp).

program test_nhinterp
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_nhinterp
  use share_test_utils
  implicit none
  integer(ik4), parameter :: kx = 12
  real(rkx), parameter :: pbase = 101325.0_rkx, tlp = 47.70_rkx
  real(rkx), parameter :: ts0 = 288.15_rkx, ptop = 20.0_rkx
  real(rkx), parameter :: ter = 3000.0_rkx
  real(rkx), pointer, contiguous, dimension(:) :: sigf => null()
  real(rkx), pointer, contiguous, dimension(:) :: sigh => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: ht => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: ps0 => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: ps => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pr0 => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: t0 => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rho0 => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: z0 => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: tv => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: f => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null()
  real(rkx) :: zsfc, err, zz, p, pl, pu, tl, tu
  real(rkx), dimension(kx) :: ztrue
  integer(ik4) :: k, n

  call memory_init
  call getmem(sigf,1,kx+1,'sigf')
  call getmem(sigh,1,kx,'sigh')
  call getmem(ht,1,1,1,1,'ht')
  call getmem(ps0,1,1,1,1,'ps0')
  call getmem(ps,1,1,1,1,'ps')
  call getmem(pr0,1,1,1,1,1,kx,'pr0')
  call getmem(t0,1,1,1,1,1,kx,'t0')
  call getmem(rho0,1,1,1,1,1,kx,'rho0')
  call getmem(z0,1,1,1,1,1,kx,'z0')
  call getmem(tv,1,1,1,1,1,kx,'tv')
  call getmem(f,1,1,1,1,1,kx,'f')
  call getmem(pp,1,1,1,1,1,kx,'pp')
  sigf = [0.0_rkx, 0.1_rkx, 0.2_rkx, 0.3_rkx, 0.4_rkx, 0.5_rkx, &
          0.6_rkx, 0.7_rkx, 0.8_rkx, 0.88_rkx, 0.94_rkx, 0.98_rkx, &
          1.0_rkx]
  do k = 1, kx
    sigh(k) = 0.5_rkx*(sigf(k)+sigf(k+1))
  end do
  ht = ter

  ! ptop in cb, the base pressure in Pa
  call nhsetup(ptop/1000.0_rkx,pbase,tlp,ts0)
  call nhbase(1,1,1,1,kx,sigh,ht,ps0,pr0,t0,rho0,z0)
  zsfc = zref(ps0(1,1)+ptop)
  call check(abs(zsfc-ter) < 0.01_rkx, &
    'reference surface pressure at the terrain height')
  err = 0.0_rkx
  do k = 1, kx
    err = max(err,abs(t0(1,1,k) - &
              max(ts0+tlp*log(pr0(1,1,k)/pbase),tiso)))
    err = max(err,abs(rho0(1,1,k)-pr0(1,1,k)/(rgas*t0(1,1,k))))
  end do
  call check(err < 1.0e-6_rkx,'reference temperature and density')

  ! Heights above the surface by numerical integration of the reference
  ! temperature (without the isothermal limit, as the formula).
  err = 0.0_rkx
  do k = 1, kx
    zz = 0.0_rkx
    do n = 1, 2000
      pl = (ps0(1,1)+ptop) + (pr0(1,1,k)-(ps0(1,1)+ptop)) * &
           real(n-1,rkx)/2000.0_rkx
      pu = (ps0(1,1)+ptop) + (pr0(1,1,k)-(ps0(1,1)+ptop)) * &
           real(n,rkx)/2000.0_rkx
      tl = ts0+tlp*log(pl/pbase)
      tu = ts0+tlp*log(pu/pbase)
      zz = zz + rgas/egrav*0.5_rkx*(tl+tu)*log(pl/pu)
    end do
    ztrue(k) = zz
    err = max(err,abs(z0(1,1,k)-zz))
  end do
  call check(err < 1.0_rkx,'reference heights over 3000 m of terrain')

  ! An atmosphere equal to the reference one: tv on the layers, ps in cb.
  ps(1,1) = ps0(1,1)/1000.0_rkx
  do k = 1, kx
    p = sigh(k)*ps0(1,1) + ptop
    tv(1,1,k) = ts0 + tlp*log(p/pbase)
  end do
  call nhbase(1,1,1,1,kx,sigh,ht,ps0,pr0,t0,rho0,z0)
  call column_heights
  call nhinterp(1,1,1,1,kx,sigh,sigf,f,tv,ps,ps0,1)
  err = 0.0_rkx
  do k = 2, kx-1
    err = max(err,abs(f(1,1,k)-(2.0_rkx*ztrue(k)+5.0_rkx)))
  end do
  call check(err < 2.0_rkx, &
    'nhinterp of a field linear in height (1 m)')

  do k = 1, kx
    p = sigh(k)*ps0(1,1) + ptop
    t0(1,1,k) = ts0 + tlp*log(p/pbase)
  end do
  call nhpp(1,1,1,1,kx,sigf,t0,pr0,t0,t0,ps,ps0,pp)
  call check(maxval(abs(pp)) < 1.0e-3_rkx, &
    'no pressure perturbation for the reference atmosphere')

  call memory_destroy
  call test_end('test_nhinterp')

  contains

  ! Height of the reference atmosphere at pressure p (Pa) above p0
  real(rkx) function zref(p)
    implicit none
    real(rkx), intent(in) :: p
    real(rkx) :: x
    x = log(p/pbase)
    zref = -(0.5_rkx*rgas/egrav*tlp*x*x + rgas/egrav*ts0*x)
  end function zref

  ! The input field: linear in the height of the layers of the column
  subroutine column_heights
    implicit none
    real(rkx) :: zq1, zq2
    integer(ik4) :: kk
    zq2 = 0.0_rkx
    do kk = kx, 1, -1
      zq1 = zq2 - rgas/egrav*tv(1,1,kk) * &
        log((sigf(kk)*ps(1,1)+ptop/1000.0_rkx) / &
            (sigf(kk+1)*ps(1,1)+ptop/1000.0_rkx))
      f(1,1,kk) = 2.0_rkx*0.5_rkx*(zq1+zq2) + 5.0_rkx
      zq2 = zq1
    end do
  end subroutine column_heights

end program test_nhinterp
