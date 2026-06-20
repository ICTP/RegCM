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

module mod_lowpass

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_memutil

  implicit none

  integer(ik4) :: nlon, nlat, km, lm
  real(rkx), pointer, dimension(:,:), contiguous :: bvx, bvy

  private

  public :: lowpass_init, lowpass_filter

  contains

  subroutine lowpass_init(lat,lon)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: lat, lon
    real(rkx), dimension(:), allocatable :: px, py
    real(rkx) :: dx, dy
    integer(ik4) :: i, j, k, l

    nlon = size(lon,1)
    nlat = size(lat,2)
    km = nlon/40
    lm = nlat/40
    dx = mathpi/real(nlon-1,rkx)
    dy = mathpi/real(nlat-1,rkx)
    allocate(px(2*km), py(2*lm))
    call getmem(bvx,1,nlon,1,2*km,'lowpass::bvx')
    call getmem(bvy,1,nlat,1,2*lm,'lowpass::bvy')
    do k = 1, 2*km
      px(k) = exp(-(real(k,rkx)/real(km,rkx))**2)
    end do
    do l = 1, 2*lm
      py(l) = exp(-(real(l,rkx)/real(lm,rkx))**2)
    end do
    do k = 1, 2*km
      do i = 1, nlon
        bvx(i,k) = sqrt(2.0_rkx/real(nlon-1,rkx)*px(k))*sin(k*(i-2)*dx)
      end do
    end do
    do l = 1, 2*lm
      do j = 1, nlat
        bvy(j,l) = sqrt(2.0_rkx/real(nlat-1,rkx)*py(l))*sin(l*(j-2)*dy)
      end do
    end do
    deallocate(px,py)
  end subroutine lowpass_init

  subroutine lowpass_filter(f)
    implicit none
    real(rkx), dimension(:,:), intent(inout) :: f
    real(rkx), dimension(:,:), allocatable :: x, sx, sy
    integer(ik4) :: i, j, k, l

    allocate(sx(nlat,2*km), sy(nlon,2*lm))
    allocate(x(nlon,nlat))
    x(:,:) = f(:,:)
    do k = 1, 2*km
      do j = 1, nlat
        sx(j,k) = 0.0_rkx
        do i = 2, nlon-1
          sx(j,k) = sx(j,k) + x(i,j)*bvx(i,k)
        end do
      end do
    end do
    x(:,:) = 0.0_rkx
    do k = 1, 2*km
      do j = 1, nlat
        do i = 1, nlon
          x(i,j) = x(i,j) + sx(j,k)*bvx(i,k)
        end do
      end do
    end do
    do l = 1, 2*lm
      do i = 1, nlon  ! y-transform
        sy(i,l) = 0.
        do j = 2, nlat-1
          sy(i,l) = sy(i,l) + x(i,j)*bvy(j,l)
        end do
      end do
    end do
    x(:,:) = 0.0_rkx
    do l = 1, 2*lm
      do j = 1, nlat
        do i = 1, nlon
          x(i,j) = x(i,j) + sy(i,l)*bvy(j,l)
        end do
      end do
    end do
    f(:,:) = f(:,:) - 0.250_rkx * (f(:,:)-x(:,:))
    deallocate(x,sx,sy)
  end subroutine lowpass_filter

end module mod_lowpass

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
