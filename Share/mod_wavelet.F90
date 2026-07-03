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

module mod_wavelet

  use mod_realkinds
  use mod_intkinds
  use mod_constants

  implicit none

  private

  public :: wavelet_denoise

  contains

  subroutine wavelet_denoise(nx, ny, image_in, image_out, threshold)
    implicit none
    integer(ik4), intent(in) :: nx, ny
    real(rkx), dimension(nx, ny), intent(in)  :: image_in
    real(rkx), dimension(nx, ny), intent(out) :: image_out
    real(rkx), intent(in) :: threshold

    ! Internal dimensions (guaranteed to be even)
    integer(ik4) :: nx_e, ny_e, hx, hy
    integer(ik4) :: i, j

    ! Dynamic internal arrays
    real(rkx), allocatable, dimension(:,:) :: img_e, coeff, temp, out_e
    real(rkx) :: val

    nx_e = nx + mod(nx, 2)
    ny_e = ny + mod(ny, 2)

    hx = nx_e / 2
    hy = ny_e / 2

    allocate(img_e(nx_e, ny_e))
    allocate(coeff(nx_e, ny_e))
    allocate(temp(nx_e, ny_e))
    allocate(out_e(nx_e, ny_e))

    img_e(1:nx, 1:ny) = image_in

    ! apply Edge-Reflective Padding if dimensions are odd
    ! Replicate the last column into the extra padded column
    ! Handle the isolated bottom-right corner pixel
    if (mod(nx, 2) /= 0) then
      img_e(nx_e, 1:ny) = image_in(nx, 1:ny)
    end if
    if (mod(ny, 2) /= 0) then
      img_e(1:nx, ny_e) = image_in(1:nx, ny)
    end if
    if (mod(nx, 2) /= 0 .and. mod(ny, 2) /= 0) then
      img_e(nx_e, ny_e) = image_in(nx, ny)
    end if

    ! ================================================
    ! FORWARD 2D HAAR WAVELET TRANSFORM (On Even Grid)
    ! ================================================
    do j = 1, ny_e
      do i = 1, hx
        temp(i, j)    = (img_e(2*i-1, j) + img_e(2*i, j)) * rsqrt2
        temp(i+hx, j) = (img_e(2*i-1, j) - img_e(2*i, j)) * rsqrt2
      end do
    end do
    do j = 1, hy
      do i = 1, nx_e
        coeff(i, j)    = (temp(i, 2*j-1) + temp(i, 2*j)) * rsqrt2
        coeff(i, j+hy) = (temp(i, 2*j-1) - temp(i, 2*j)) * rsqrt2
      end do
    end do

    ! =================
    ! SOFT THRESHOLDING
    ! =================
    do j = 1, ny_e
      do i = 1, nx_e
        ! Skip Low-Frequency Appoximation
        if ( i <= hx .and. j <= hy ) cycle
        val = coeff(i, j)
        if ( abs(val) <= threshold ) then
          coeff(i, j) = 0.0_rkx
        else
          coeff(i, j) = sign(1.0_rkx, val) * (abs(val) - threshold)
        end if
      end do
    end do

    ! =================================
    ! INVERSE 2D HAAR WAVELET TRANSFORM
    ! =================================
    do j = 1, hy
      do i = 1, nx_e
        temp(i, 2*j-1) = (coeff(i, j) + coeff(i, j+hy)) * rsqrt2
        temp(i, 2*j)   = (coeff(i, j) - coeff(i, j+hy)) * rsqrt2
      end do
    end do
    do j = 1, ny_e
      do i = 1, hx
        out_e(2*i-1, j) = (temp(i, j) + temp(i+hx, j)) * rsqrt2
        out_e(2*i, j)   = (temp(i, j) - temp(i+hx, j)) * rsqrt2
      end do
    end do

    ! ====
    ! CROP
    ! ====
    image_out = out_e(1:nx, 1:ny)
    deallocate(img_e, coeff, temp, out_e)

  end subroutine wavelet_denoise

end module mod_wavelet

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
