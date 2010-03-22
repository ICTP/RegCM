!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine param(nx,ny,kz,xlat,xlon,varmin,varmax,xlat1d,xlon1d,  &
                     & xlonmin,xlonmax,xlatmin,xlatmax,iadim,ndim)
 
      implicit none
!
! Dummy arguments
!
      real(4) :: xlatmax , xlatmin , xlonmax , xlonmin
      integer :: kz , ndim , nx , ny
      integer , dimension(ndim) :: iadim
      real(4) , dimension(ndim) :: varmax , varmin
      real(4) , dimension(ny,nx) :: xlat , xlon
      real(4) , dimension(ny) :: xlat1d
      real(4) , dimension(nx) :: xlon1d
      intent (in) kz , ndim , nx , ny , xlat , xlon
      intent (out) iadim , varmax , varmin , xlat1d , xlatmax , xlatmin ,&
                 & xlon1d , xlonmax , xlonmin
!
! Local variables
!
      integer :: i , j
!
      varmin(1) = minval(xlon)
      varmin(2) = minval(xlat)
      varmin(3) = 1 !1050.
      varmax(1) = maxval(xlon)
      varmax(2) = maxval(xlat)
      varmax(3) = kz !1050.
      iadim(1) = nx
      iadim(2) = ny
      iadim(3) = kz
      do i = 1 , nx
        xlon1d(i) = xlon(ny/2,i)
      end do
      do j = 1 , ny
        xlat1d(j) = xlat(j,nx/2)
      end do
      xlonmin = minval(xlon)
      xlonmax = maxval(xlon)
      xlatmin = minval(xlat)
      xlatmax = maxval(xlat)
 
      end subroutine param
