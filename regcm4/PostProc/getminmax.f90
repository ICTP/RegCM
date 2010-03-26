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

      subroutine getminmax(f,n1,n2,n3,vmin,vmax,vmisdat)
 
      implicit none
!
! Dummy arguments
!
      integer :: n1 , n2 , n3
      real(4) :: vmax , vmin , vmisdat
      real(4) , dimension(n1,n2,n3) :: f
      intent (in) f , n1 , n2 , n3 , vmisdat
      intent (inout) vmax , vmin
!
! Local variables
!
      integer :: i , j , k
      real(4) :: misdat
!
      if ( vmisdat>0.0 ) then
        misdat = -1.0*vmisdat
      else
        misdat = vmisdat
      end if
      vmin = 1.E30
      vmax = -1.E30
      do k = 1 , n3
        do j = 1 , n2
          do i = 1 , n1
            if ( f(i,j,k)>misdat ) then
              if ( f(i,j,k)>vmax ) vmax = f(i,j,k)
              if ( f(i,j,k)<vmin ) vmin = f(i,j,k)
            end if
          end do
        end do
      end do
 
      end subroutine getminmax
