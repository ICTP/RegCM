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
!    MDRCHANTABILITY or FITNDSS FOR A PARTICULAR PURPOSD.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_o3blk
      use mod_dynparam
      implicit none
      real(8) , dimension(31) :: o3ann , o3sum , o3win , o3wrk , ppann ,&
                               & ppsum , ppwin , ppwrk
      real(8) , dimension(32) :: ppwrkh
      real(8) , allocatable , dimension(:) :: prlevh
      private :: o3ann , o3sum , o3win , o3wrk , ppann , ppsum , ppwin ,&
              &  ppwrk , ppwrkh , prlevh
!
      data o3sum/5.297D-8 , 5.852D-8 , 6.579D-8 , 7.505D-8 , 8.577D-8 , &
         & 9.895D-8 , 1.175D-7 , 1.399D-7 , 1.677D-7 , 2.003D-7 ,       &
         & 2.571D-7 , 3.325D-7 , 4.438D-7 , 6.255D-7 , 8.168D-7 ,       &
         & 1.036D-6 , 1.366D-6 , 1.855D-6 , 2.514D-6 , 3.240D-6 ,       &
         & 4.033D-6 , 4.854D-6 , 5.517D-6 , 6.089D-6 , 6.689D-6 ,       &
         & 1.106D-5 , 1.462D-5 , 1.321D-5 , 9.856D-6 , 5.960D-6 ,       &
         & 5.960D-6/
      data ppsum/955.890 , 850.532 , 754.599 , 667.742 , 589.841 ,      &
         & 519.421 , 455.480 , 398.085 , 347.171 , 301.735 , 261.310 ,  &
         & 225.360 , 193.419 , 165.490 , 141.032 , 120.125 , 102.689 ,  &
         & 87.829 , 75.123 , 64.306 , 55.086 , 47.209 , 40.535 ,        &
         & 34.795 , 29.865 , 19.122 , 9.277 , 4.660 , 2.421 , 1.294 ,   &
         & 0.647/
!
      data o3win/4.629D-8 , 4.686D-8 , 5.017D-8 , 5.613D-8 , 6.871D-8 , &
         & 8.751D-8 , 1.138D-7 , 1.516D-7 , 2.161D-7 , 3.264D-7 ,       &
         & 4.968D-7 , 7.338D-7 , 1.017D-6 , 1.308D-6 , 1.625D-6 ,       &
         & 2.011D-6 , 2.516D-6 , 3.130D-6 , 3.840D-6 , 4.703D-6 ,       &
         & 5.486D-6 , 6.289D-6 , 6.993D-6 , 7.494D-6 , 8.197D-6 ,       &
         & 9.632D-6 , 1.113D-5 , 1.146D-5 , 9.389D-6 , 6.135D-6 ,       &
         & 6.135D-6/
      data ppwin/955.747 , 841.783 , 740.199 , 649.538 , 568.404 ,      &
         & 495.815 , 431.069 , 373.464 , 322.354 , 277.190 , 237.635 ,  &
         & 203.433 , 174.070 , 148.949 , 127.408 , 108.915 , 93.114 ,   &
         & 79.551 , 67.940 , 58.072 , 49.593 , 42.318 , 36.138 ,        &
         & 30.907 , 26.362 , 16.423 , 7.583 , 3.620 , 1.807 , 0.938 ,   &
         & 0.469/

      contains

      subroutine allocate_mod_o3blk
        implicit none
        allocate(prlevh(kzp2))
      end subroutine allocate_mod_o3blk
!
!----------------------------------------------------------------------
!
      subroutine o3data
!
      use mod_dynparam
      use mod_runparams , only : r8pt , sigma
      use mod_main
      use mod_rad
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
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
        do i = 1 , iym1
          do k = klevp1 , 1 , -1
            kj = klevp1 - k + 1
            prlevh(kj) = (sigma(k)*psb(i,j)+r8pt)*10.
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
      end module mod_o3blk
