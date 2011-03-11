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

      use mod_constants
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
               & 9.895D-8 , 1.175D-7 , 1.399D-7 , 1.677D-7 , 2.003D-7 , &
               & 2.571D-7 , 3.325D-7 , 4.438D-7 , 6.255D-7 , 8.168D-7 , &
               & 1.036D-6 , 1.366D-6 , 1.855D-6 , 2.514D-6 , 3.240D-6 , &
               & 4.033D-6 , 4.854D-6 , 5.517D-6 , 6.089D-6 , 6.689D-6 , &
               & 1.106D-5 , 1.462D-5 , 1.321D-5 , 9.856D-6 , 5.960D-6 , &
               & 5.960D-6/
      data ppsum      /955.890D0 , 850.532D0 , 754.599D0 , 667.742D0 , &
           589.841D0 , 519.421D0 , 455.480D0 , 398.085D0 , 347.171D0 , &
           301.735D0 , 261.310D0 , 225.360D0 , 193.419D0 , 165.490D0 , &
           141.032D0 , 120.125D0 , 102.689D0 ,  87.829D0 ,  75.123D0 , &
            64.306D0 ,  55.086D0 ,  47.209D0 ,  40.535D0 ,  34.795D0 , &
            29.865D0 ,  19.122D0 ,   9.277D0 ,   4.660D0 ,   2.421D0 , &
             1.294D0 ,   0.647D0/
!
      data o3win/4.629D-8 , 4.686D-8 , 5.017D-8 , 5.613D-8 , 6.871D-8 , &
               & 8.751D-8 , 1.138D-7 , 1.516D-7 , 2.161D-7 , 3.264D-7 , &
               & 4.968D-7 , 7.338D-7 , 1.017D-6 , 1.308D-6 , 1.625D-6 , &
               & 2.011D-6 , 2.516D-6 , 3.130D-6 , 3.840D-6 , 4.703D-6 , &
               & 5.486D-6 , 6.289D-6 , 6.993D-6 , 7.494D-6 , 8.197D-6 , &
               & 9.632D-6 , 1.113D-5 , 1.146D-5 , 9.389D-6 , 6.135D-6 , &
               & 6.135D-6/
      data ppwin      /955.747D0 , 841.783D0 , 740.199D0 , 649.538D0 , &
           568.404D0 , 495.815D0 , 431.069D0 , 373.464D0 , 322.354D0 , &
           277.190D0 , 237.635D0 , 203.433D0 , 174.070D0 , 148.949D0 , &
           127.408D0 , 108.915D0 ,  93.114D0 ,  79.551D0 ,  67.940D0 , &
            58.072D0 ,  49.593D0 ,  42.318D0 ,  36.138D0 ,  30.907D0 , &
            26.362D0 ,  16.423D0 ,   7.583D0 ,   3.620D0 ,   1.807D0 , &
             0.938D0 ,   0.469D0/

      contains

      subroutine allocate_mod_o3blk
        implicit none
        allocate(prlevh(kzp2))
        prlevh = d_zero
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
      integer :: i , j , jj , k , kj , klevp1
      real(8) :: pb1 , pb2 , pt1 , pt2
!
      do k = 1 , 31
        ppann(k) = ppsum(k)
      end do
      o3ann(1) = 0.5D0*(o3sum(1)+o3win(1))
!
      do k = 2 , 31
        o3ann(k) = o3win(k-1) + (o3win(k)-o3win(k-1))                   &
                 & /(ppwin(k)-ppwin(k-1))*(ppsum(k)-ppwin(k-1))
      end do
      do k = 2 , 31
        o3ann(k) = 0.5D0*(o3ann(k)+o3sum(k))
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
            prlevh(kj) = (sigma(k)*sps2%ps(i,j)+r8pt)*d_10
          end do
          ppwrkh(1) = 1100.0D0
          do k = 2 , 31
            ppwrkh(k) = (ppwrk(k)+ppwrk(k-1))/d_two
          end do
          ppwrkh(32) = d_zero
          do k = 1 , kz
            o3prof(i,k,j) = d_zero
            do jj = 1 , 31
              if ( (-(prlevh(k)-ppwrkh(jj))) >= d_zero ) then
                pb1 = d_zero
              else
                pb1 = prlevh(k) - ppwrkh(jj)
              end if
              if ( (-(prlevh(k)-ppwrkh(jj+1))) >= d_zero ) then
                pb2 = d_zero
              else
                pb2 = prlevh(k) - ppwrkh(jj+1)
              end if
              if ( (-(prlevh(k+1)-ppwrkh(jj))) >= d_zero ) then
                pt1 = d_zero
              else
                pt1 = prlevh(k+1) - ppwrkh(jj)
              end if
              if ( (-(prlevh(k+1)-ppwrkh(jj+1))) >= d_zero ) then
                pt2 = d_zero
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
!
      end module mod_o3blk
