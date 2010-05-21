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
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      program tograds1
      implicit none
!
! Local variables
!
      real(4) , dimension(44,44) :: a
      integer :: i , j , k
      integer(2) , dimension(44,44) :: ia
      real(8) :: offset , xscale
!
      open (10,file='EH_RF1961JAN',form='unformatted',recl=44*44*2+16,  &
          & access='direct')
      open (20,file='JAN61.dat',form='unformatted',recl=44*44*4,        &
           &access='direct')
      do k = 1 , 105
        read (10,rec=k) offset , xscale , ia
        do j = 1 , 44
          do i = 1 , 44
            a(i,j) = ia(i,j)*xscale + offset
          end do
        end do
        write (20,rec=k) a
      end do
      end program tograds1
