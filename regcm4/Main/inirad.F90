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
 
      subroutine inirad
 
      use mod_dynparam
      use mod_o3blk , only : o3data
      use mod_rad , only : o3prof , heatrt
      use mod_date , only : jyear , jyear0 , ktau
      implicit none
!
! Local variables
!
      integer :: i , j , k
!
!     compute ozone mixing ratio distribution
!
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendl
            do i = 1 , iym1
              heatrt(i,k,j) = 0.
              o3prof(i,k,j) = 0.
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , iym1
              heatrt(i,k,j) = 0.
              o3prof(i,k,j) = 0.
            end do
          end do
#endif
        end do
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iym1
            o3prof(i,kzp1,j) = 0.
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , iym1
            o3prof(i,kzp1,j) = 0.
          end do
        end do
#endif
        call o3data
#ifdef MPP1
        if ( myid.eq.0 ) then
          print * , 'ozone profiles'
          do k = 1 , kzp1
            write (6,99001) o3prof(3,k,jxp)
          end do
        end if
#else
        print * , 'ozone profiles'
        do k = 1 , kzp1
          write (6,99001) o3prof(3,k,3)
        end do
#endif
      end if
99001 format (1x,7E12.4)
 
      end subroutine inirad
