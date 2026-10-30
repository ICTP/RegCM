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

! Unit test of Share/mod_wavelet.F90: one level Haar transform with soft
! thresholding of the detail coefficients. With a zero threshold the field
! is reproduced; a field constant on 2x2 blocks keeps its value for any
! threshold; a checkerboard of amplitude a has diagonal detail coefficients
! 2a and is removed by a threshold of 2a; odd sizes are handled.

program test_wavelet
  use mod_intkinds
  use mod_realkinds
  use mod_wavelet
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nx = 20, ny = 16
  real(rkx), dimension(nx,ny) :: f, g, base
  real(rkx), dimension(nx-1,ny-3) :: fo, go
  real(rkx), parameter :: amp = 0.3_rkx
  integer(ik4) :: i, j

  do j = 1, ny
    do i = 1, nx
      f(i,j) = sin(0.3_rkx*real(i,rkx)) + cos(0.2_rkx*real(j,rkx)) + &
               0.01_rkx*real(mod(i*j,7),rkx)
      base(i,j) = real((i+1)/2 + 3*((j+1)/2),rkx)
    end do
  end do

  call wavelet_denoise(nx,ny,f,g,0.0_rkx)
  call check(maxval(abs(g-f)) < 1.0e-12_rkx,'zero threshold reproduces')

  call wavelet_denoise(nx,ny,base,g,100.0_rkx)
  call check(maxval(abs(g-base)) < 1.0e-12_rkx, &
    'a field constant on 2x2 blocks is kept')

  do j = 1, ny
    do i = 1, nx
      f(i,j) = base(i,j) + amp*real(1-2*mod(i+j,2),rkx)
    end do
  end do
  call wavelet_denoise(nx,ny,f,g,2.0_rkx*amp)
  call check(maxval(abs(g-base)) < 1.0e-12_rkx,'checkerboard removed')
  call wavelet_denoise(nx,ny,f,g,amp)
  call check(maxval(abs(abs(g-base)-0.5_rkx*amp)) < 1.0e-12_rkx, &
    'half the threshold: checkerboard halved')

  do j = 1, ny-3
    do i = 1, nx-1
      fo(i,j) = sin(0.3_rkx*real(i,rkx))*cos(0.25_rkx*real(j,rkx))
    end do
  end do
  call wavelet_denoise(nx-1,ny-3,fo,go,0.0_rkx)
  call check(maxval(abs(go-fo)) < 1.0e-12_rkx,'odd sizes, zero threshold')
  go = 0.0_rkx
  call wavelet_denoise(nx-1,ny-3,fo,go,0.05_rkx)
  call check(maxval(abs(go-fo)) < 0.2_rkx .and. &
             abs(sum(go)-sum(fo)) < 0.2_rkx,'odd sizes, smooth field kept')

  call test_end('test_wavelet')
end program test_wavelet
