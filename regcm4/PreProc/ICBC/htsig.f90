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

      subroutine htsig(t,h,p3d,ps,ht,im,jm,km)
      use mod_constants , only : rgti, rgas
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km
      real , dimension(im,jm,km) :: h , p3d , t
      real , dimension(im,jm) :: ht , ps
      intent (in) ht , im , jm , km , p3d , ps , t
      intent (inout) h
!
! Local variables
!
      real :: tbar
      integer :: i , j , k
!
      do j = 1 , jm
        do i = 1 , im
          if ( ps(i,j)>-9995.0 ) then
            h(i,j,km) = ht(i,j) + rgas*rgti*t(i,j,km)                   &
                      & *log(ps(i,j)/p3d(i,j,km))
          else
            h(i,j,km) = -9999.0
          end if
        end do
      end do
      do k = km - 1 , 1 , -1
        do j = 1 , jm
          do i = 1 , im
            if ( h(i,j,k+1)>-9995.0 ) then
              tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
              h(i,j,k) = h(i,j,k+1)                                     &
                       & + rgas*rgti*tbar*log(p3d(i,j,k+1)/p3d(i,j,k))
            else
              h(i,j,k) = -9999.0
            end if
          end do
        end do
      end do
      end subroutine htsig
