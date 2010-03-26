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

      subroutine calcmse2d(f2d,zs,nx,ny,nfld,nta,nqa,nmse)
 
      use mod_constants , only : cpd , gti , wlhv

      implicit none
!
! Dummy arguments
!
      integer :: nfld , nmse , nqa , nta , nx , ny
      real(4) , dimension(nx,ny,nfld) :: f2d
      real(4) , dimension(nx,ny) :: zs
      intent (in) nfld , nmse , nqa , nta , nx , ny , zs
      intent (inout) f2d
!
! Local variables
!
      integer :: i , j
      real(4) :: qa , ta , za
!
      do j = 1 , ny
        do i = 1 , nx
          ta = f2d(i,j,nta)
          qa = f2d(i,j,nqa)
          za = zs(i,j) + 2.
          f2d(i,j,nmse) = gti*za + cpd*ta + wlhv*qa
        end do
      end do
 
      end subroutine calcmse2d
