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

! Unit test of Share/mod_ensemble.F90: randify adds to each value a noise
! within a fraction of it, of zero mean, leaving zeros unchanged;
! random_pick splits a total in random positive parts that sum to it, of
! equal mean, and gives zeros for a zero total.

program test_ensemble
  use mod_intkinds
  use mod_realkinds
  use mod_ensemble
  use share_test_utils
  implicit none
  integer(ik4), parameter :: ni = 40, nj = 30, nk = 5, nv = 24
  integer(ik4), parameter :: ntry = 20000
  real(rkx), parameter :: frac = 0.1_rkx
  real(rkx), dimension(ni,nj) :: f2, g2
  real(rkx), dimension(ni,nj,nk) :: f3, g3
  real(rk8), pointer, contiguous, dimension(:) :: p8 => null()
  real(rk4), pointer, contiguous, dimension(:) :: p4 => null()
  real(rk8), dimension(nv) :: acc
  integer(ik4) :: i, j, k, n
  logical :: lok

  call fixed_seed
  do j = 1, nj
    do i = 1, ni
      f2(i,j) = real(i-ni/2,rkx) + 0.5_rkx
      do k = 1, nk
        f3(i,j,k) = 280.0_rkx + real(i+j+k,rkx)
      end do
    end do
  end do
  f2(1:3,1) = 0.0_rkx

  g2 = f2
  call randify(g2,frac,ni,nj)
  call check(all(abs(g2-f2) <= frac*abs(f2)+1.0e-12_rkx), &
    'randify 2-D: noise within the fraction, for both signs')
  call check(all(g2(1:3,1) == 0.0_rkx),'randify: zeros unchanged')
  call check(any(g2 /= f2),'randify 2-D changes the field')

  g3 = f3
  call randify(g3,frac,ni,nj,nk)
  call check(all(abs(g3-f3) <= frac*f3+1.0e-12_rkx), &
    'randify 3-D: noise within the fraction')
  call check(abs(sum((g3-f3)/f3))/real(ni*nj*nk,rkx) < 0.01_rkx*frac, &
    'randify 3-D: noise of zero mean')

  allocate(p8(nv),p4(nv))
  lok = .true.
  acc = 0.0_rk8
  do n = 1, ntry
    call random_pick(24.0_rk8,p8,nv)
    if ( any(p8 < 0.0_rk8) .or. abs(sum(p8)-24.0_rk8) > 1.0e-10_rk8 ) &
      lok = .false.
    acc = acc + p8
  end do
  call check(lok,'random_pick: positive parts summing to the total')
  call check(maxval(abs(acc/real(ntry,rk8)-1.0_rk8)) < 0.05_rk8, &
    'random_pick: equal mean of the parts')
  call random_pick(0.0_rk8,p8,nv)
  call check(all(p8 == 0.0_rk8),'random_pick: zero total')
  call random_pick(3.0_rk4,p4,nv)
  call check(abs(sum(p4)-3.0_rk4) < 1.0e-5_rk4 .and. all(p4 >= 0.0_rk4), &
    'random_pick: single precision')
  deallocate(p8,p4)

  call test_end('test_ensemble')
end program test_ensemble
