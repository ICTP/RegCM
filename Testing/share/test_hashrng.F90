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

! Unit test of Share/mod_hashrng.F90: reference values of the unsigned
! 32 bit hash (computed independently with exact integer arithmetic),
! the range of the uniform numbers, their mean, variance and histogram,
! the correlation between neighbouring points and successive steps, and
! independence from the order in which the points are visited.

program test_hashrng
  use mod_intkinds
  use mod_realkinds
  use mod_hashrng
  use share_test_utils
  implicit none
  integer(ik4), parameter :: seed = 19590209
  integer(ik4), parameter :: nj = 64, ni = 64, nk = 16, nbin = 16
  integer(ik4), dimension(nbin) :: hist
  real(rkx), dimension(nj,ni,nk) :: r, rb
  real(rk8) :: mean, var, chi2, cj, ct, rmin, rmax
  integer(ik4) :: i, j, k, n
  integer(ik8) :: step, key, key1

  call check(hash_uint32(hash_key(seed,0_ik8),1,2,2) == 2983617404_ik8, &
    'reference value 1')
  call check(hash_uint32(hash_key(seed,1_ik8),1,2,2) == 3330568734_ik8, &
    'reference value 2')
  call check(hash_uint32(hash_key(seed,123456_ik8),18,117,127) == &
    1961980565_ik8,'reference value 3')
  call check(hash_uint32(hash_key(1,8589934597_ik8),3,4,5) == 1923184515_ik8, &
    'reference value with a counter above 2**32')

  step = 1000_ik8
  key = hash_key(seed,step)
  key1 = hash_key(seed,step+1_ik8)
  do k = 1, nk
    do i = 1, ni
      do j = 1, nj
        r(j,i,k) = hash_uniform(key,k,i,j)
      end do
    end do
  end do
  ! Visit the points in the opposite order, as another decomposition
  ! or a GPU would
  do concurrent ( j = nj:1:-1, i = ni:1:-1, k = nk:1:-1 )
    rb(j,i,k) = hash_uniform(key,k,i,j)
  end do
  call check(all(r == rb),'independent of the visiting order')

  rmin = minval(r)
  rmax = maxval(r)
  call check(rmin >= 0.0_rkx .and. rmax < 1.0_rkx,'values in [0,1)')

  n = nj*ni*nk
  mean = sum(real(r,rk8))/n
  var = sum((real(r,rk8)-mean)**2)/n
  call check(abs(mean-0.5_rk8) < 0.01_rk8,'mean 1/2')
  call check(abs(var-1.0_rk8/12.0_rk8) < 0.003_rk8,'variance 1/12')

  hist = 0
  do k = 1, nk
    do i = 1, ni
      do j = 1, nj
        hist(int(r(j,i,k)*nbin)+1) = hist(int(r(j,i,k)*nbin)+1) + 1
      end do
    end do
  end do
  chi2 = sum((real(hist,rk8)-real(n,rk8)/nbin)**2)/(real(n,rk8)/nbin)
  ! 15 degrees of freedom: the 0.999 quantile is 37.7
  call check(chi2 < 37.7_rk8,'uniform histogram')

  ! Neighbouring points along j and the next time step
  do concurrent ( j = 1:nj, i = 1:ni, k = 1:nk )
    rb(j,i,k) = hash_uniform(key1,k,i,j)
  end do
  cj = sum((real(r(1:nj-1,:,:),rk8)-mean) * &
           (real(r(2:nj,:,:),rk8)-mean))/(var*(nj-1)*ni*nk)
  ct = sum((real(r,rk8)-mean)*(real(rb,rk8)-mean))/(var*n)
  call check(abs(cj) < 0.03_rk8,'no correlation between neighbours')
  call check(abs(ct) < 0.03_rk8,'no correlation between time steps')

  call check(hash_uniform(hash_key(seed+1,step),1,1,1) /= &
             hash_uniform(key,1,1,1),'seeds give different streams')

  call test_end('test_hashrng')
end program test_hashrng
