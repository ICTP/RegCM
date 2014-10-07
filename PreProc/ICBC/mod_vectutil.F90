!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_vectutil

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  contains

  subroutine p1p2(pd,px,ni,nj)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rk8) , intent(in) , dimension(ni,nj) :: px
    real(rk8) , intent(out) , dimension(ni,nj) :: pd
!
    integer(ik4) :: i , j , ni1 , nj1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    ni1 = ni - 1
    nj1 = nj - 1
!
    do j = 2 , nj1
      do i = 2 , ni1
        pd(i,j) = d_rfour*(px(i,j)+px(i-1,j)+px(i,j-1)+px(i-1,j-1))
      end do
    end do
!
    do i = 2 , ni1
      pd(i,1) = d_half*(px(i,1)+px(i-1,1))
      pd(i,nj) = d_half*(px(i,nj1)+px(i-1,nj1))
    end do
!
    do j = 2 , nj1
      pd(1,j) = d_half*(px(1,j)+px(1,j-1))
      pd(ni,j) = d_half*(px(ni1,j)+px(ni1,j-1))
    end do
!
    pd(1,1) = px(1,1)
    pd(1,nj) = px(1,nj1)
    pd(ni,1) = px(ni1,1)
    pd(ni,nj) = px(ni1,nj1)
!
  end subroutine p1p2
!
  subroutine p1p2_band(pd,px,ni,nj)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rk8) , intent(in) , dimension(ni,nj) :: px
    real(rk8) , intent(out) , dimension(ni,nj) :: pd
!
    integer(ik4) :: i , j , nj1 , im1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    nj1 = nj - 1
!
    do j = 2 , nj1
      do i = 1 , ni
        im1 = i-1
        if (im1 == 0) im1 = ni
        pd(i,j) = d_rfour*(px(i,j)+px(im1,j)+px(i,j-1)+px(im1,j-1))
      end do
    end do
!
    do i = 1 , ni
      im1 = i-1
      if (im1 == 0) im1 = ni
      pd(i,1) = d_half*(px(i,1)+px(im1,1))
      pd(i,nj) = d_half*(px(i,nj1)+px(im1,nj1))
    end do
!
  end subroutine p1p2_band
!
!-----------------------------------------------------------------------
!
  subroutine top2btm(x,nlon1,nlat1,nlev1)
    implicit none
    integer(ik4) , intent(in) :: nlat1 , nlev1 , nlon1
    real(rk8) , intent(inout) , dimension(nlon1,nlat1,nlev1) :: x
!
    integer(ik4) :: i , j , k , kr
    real(rk8) , dimension(nlev1) :: work
!
    do j = 1 , nlat1
      do i = 1 , nlon1
        do k = 1 , nlev1
          work(k) = x(i,j,k)
        end do
        do k = 1 , nlev1
          kr = nlev1 - k + 1
          x(i,j,k) = work(kr)
        end do
      end do
    end do
  end subroutine top2btm
!
end module mod_vectutil
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
