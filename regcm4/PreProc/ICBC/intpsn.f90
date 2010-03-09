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

      subroutine intpsn(psrcm,zrcm,pa,za,tlayer,pt,ni,nj)
      use mod_constants , only : govr
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj
      real :: pt
      real , dimension(ni,nj) :: pa , psrcm , tlayer , za , zrcm
      intent (in) ni , nj , pa , pt , tlayer , za , zrcm
      intent (out) psrcm
!
! Local variables
!
      real :: tb
      integer :: i , j
!
!     EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
!     USE TLAYER CALCULATED IN INTGTB.
!     PSRCM = SURFACE PRESSURE - PTOP
!
      do i = 1 , ni
        do j = 1 , nj
          tb = tlayer(i,j)
          psrcm(i,j) = pa(i,j)*exp(-govr*(zrcm(i,j)-za(i,j))/tb) - pt
        end do
      end do
 
!     PRINT *, 'ZRCM, ZA, PA, PT =', ZRCM(5,5), ZA(5,5), PA(5,5), PT
!     PRINT *, 'TLAYER(5,5), PSRCM(5,5) = ', TLAYER(5,5), PSRCM(5,5)
 
      end subroutine intpsn
