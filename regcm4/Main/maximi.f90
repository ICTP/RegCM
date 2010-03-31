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
 
      subroutine maximi(array,iy,kz,ks,ke,imax,istart,iend)

      implicit none
!
! Dummy arguments
!
      integer :: iend , istart , iy , ke , ks , kz
      real(8) , dimension(iy,kz) :: array
      integer , dimension(iy) :: imax
      intent (in) array , iend , istart , iy , ke , ks , kz
      intent (out) imax
!
! Local variables
!
      integer :: i , k
      real(8) :: x , xar
!
      do i = istart , iend
        imax(i) = ks
        x = array(i,ks)
        do k = ks , ke
          xar = array(i,k)
          if ( xar.ge.x ) then
            x = xar
            imax(i) = k
          end if
        end do
      end do
!
      end subroutine maximi
