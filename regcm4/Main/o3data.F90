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
 
!----------------------------------------------------------------------
      subroutine o3data
!
      use mod_dynparam
      use mod_param3 , only : ptop , sigma
      use mod_main
      use mod_rad
      use mod_o3blk
      implicit none
!
! Local variables
!
      integer :: i , j , jj , k , kj , klevp1
      real(8) :: pb1 , pb2 , pt1 , pt2
!
      do k = 1 , 31
        ppann(k) = ppsum(k)
      end do
      o3ann(1) = 0.5*(o3sum(1)+o3win(1))
!
      do k = 2 , 31
        o3ann(k) = o3win(k-1) + (o3win(k)-o3win(k-1))                   &
                 & /(ppwin(k)-ppwin(k-1))*(ppsum(k)-ppwin(k-1))
      end do
      do k = 2 , 31
        o3ann(k) = 0.5*(o3ann(k)+o3sum(k))
      end do
      do k = 1 , 31
        o3wrk(k) = o3ann(k)
        ppwrk(k) = ppann(k)
      end do
!
!     calculate half pressure levels for model and data levels
!
      klevp1 = kzp1
!
#ifdef MPP1
      do j = 1 , jendx
#else
      do j = 1 , jxm1
#endif
        do i = 1 , iym1
          do k = klevp1 , 1 , -1
            kj = klevp1 - k + 1
            prlevh(kj) = (sigma(k)*psb(i,j)+ptop)*10.
          end do
          ppwrkh(1) = 1100.
          do k = 2 , 31
            ppwrkh(k) = (ppwrk(k)+ppwrk(k-1))/2.
          end do
          ppwrkh(32) = 0.
          do k = 1 , kz
            o3prof(i,k,j) = 0.
            do jj = 1 , 31
              if ( (-(prlevh(k)-ppwrkh(jj))).ge.0. ) then
                pb1 = 0.
              else
                pb1 = prlevh(k) - ppwrkh(jj)
              end if
              if ( (-(prlevh(k)-ppwrkh(jj+1))).ge.0. ) then
                pb2 = 0.
              else
                pb2 = prlevh(k) - ppwrkh(jj+1)
              end if
              if ( (-(prlevh(k+1)-ppwrkh(jj))).ge.0. ) then
                pt1 = 0.
              else
                pt1 = prlevh(k+1) - ppwrkh(jj)
              end if
              if ( (-(prlevh(k+1)-ppwrkh(jj+1))).ge.0. ) then
                pt2 = 0.
              else
                pt2 = prlevh(k+1) - ppwrkh(jj+1)
              end if
              o3prof(i,k,j) = o3prof(i,k,j) + (pb2-pb1-pt2+pt1)         &
                            & *o3wrk(jj)
            end do
            o3prof(i,k,j) = o3prof(i,k,j)/(prlevh(k)-prlevh(k+1))
          end do
        end do
      end do
!
      end subroutine o3data
