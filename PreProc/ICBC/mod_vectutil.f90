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

  use m_realkinds
  use m_stdio

  contains

  subroutine mxmn3d(var,cvar,jx,iy,np)
  implicit none
!
  character(2) :: cvar
  integer :: iy , jx , np
  real(sp) , dimension(jx,iy,np) :: var
  intent (in) cvar , iy , jx , np , var
!
  integer :: i , j , k
  real(sp) :: smax , smin
!
  do k = 1 , np
    smax = -1.E8
    smin = 1.E8
    do j = 1 , iy
      do i = 1 , jx
        if ( smax<var(i,j,k) ) smax = var(i,j,k)
        if ( smin>var(i,j,k) ) smin = var(i,j,k)
      end do
    end do
    write (stdout,*) cvar , k , smax , smin
  end do
  end subroutine mxmn3d
!
!-----------------------------------------------------------------------
!
  subroutine p1p2(pd,px,ni,nj)
  implicit none
!
  integer :: ni , nj
  real(sp) , dimension(ni,nj) :: pd , px
  intent (in) ni , nj , px
  intent (out) pd
!
  integer :: i , j , ni1 , nj1
!
!     THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
!     ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
!     SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE
!     I'S.
!
  ni1 = ni - 1
  nj1 = nj - 1
!
  do j = 2 , nj1
    do i = 2 , ni1
      pd(i,j) = 0.25*(px(i,j)+px(i-1,j)+px(i,j-1)+px(i-1,j-1))
    end do
  end do
!
  do i = 2 , ni1
    pd(i,1) = 0.5*(px(i,1)+px(i-1,1))
    pd(i,nj) = 0.5*(px(i,nj1)+px(i-1,nj1))
  end do
!
  do j = 2 , nj1
    pd(1,j) = 0.5*(px(1,j)+px(1,j-1))
    pd(ni,j) = 0.5*(px(ni1,j)+px(ni1,j-1))
  end do
!
  pd(1,1) = px(1,1)
  pd(1,nj) = px(1,nj1)
  pd(ni,1) = px(ni1,1)
  pd(ni,nj) = px(ni1,nj1)
!
  end subroutine p1p2
!
!
  subroutine p1p2_band(pd,px,ni,nj)
  implicit none
!
  integer :: ni , nj
  real(sp) , dimension(ni,nj) :: pd , px
  intent (in) ni , nj , px
  intent (out) pd
!
  integer :: i , j , nj1 , im1
!
!     THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
!     ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
!     SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE
!     I'S.
!
  nj1 = nj - 1
!
  do j = 2 , nj1
    do i = 1 , ni
      im1=i-1
      if(im1 == 0) im1=ni
      pd(i,j) = 0.25*(px(i,j)+px(im1,j)+px(i,j-1)+px(im1,j-1))
    end do
  end do
!
  do i = 1 , ni
    im1=i-1
    if(im1 == 0) im1=ni
    pd(i,1) = 0.5*(px(i,1)+px(im1,1))
    pd(i,nj) = 0.5*(px(i,nj1)+px(im1,nj1))
  end do
!
  end subroutine p1p2_band
!
!-----------------------------------------------------------------------
!
  subroutine top2btm(x,nlon1,nlat1,nlev1)
  implicit none
!
  integer :: nlat1 , nlev1 , nlon1
  real(sp) , dimension(nlon1,nlat1,nlev1) :: x
  intent (in) nlat1 , nlev1 , nlon1
  intent (inout) x
!
  integer :: i , j , k , kr
  real(sp) , dimension(nlev1) :: work
!
  do i = 1 , nlon1
    do j = 1 , nlat1
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
