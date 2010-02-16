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
 
      subroutine initlk(veg2d,ix1,jx1)
 
      use mod_regcm_param
      use mod_lake
      implicit none
!
! Dummy arguments
!
      integer :: ix1 , jx1
      real(8) , dimension(ix1,jx1) :: veg2d
      intent (in) ix1 , jx1
      intent (out) veg2d
!
! Local variables
!
      integer :: depth , freeze , idata , ilake , j , jlake , lakeset , &
               & n
      real(8) :: eta , hi , hice , hsnow
      real(8) , dimension(400) :: t
!
 
!     ******  lkpts = number of lake points desired

      data hi , hice , hsnow , eta/.01 , 0. , 0. , .5/
 
!     ******  unit number containing lake pt locations, depths
      idata = 40
 
      lcount = 0
      iin = 41
      iout = 42
      numpts = lkpts
 
!     read in vegetation type desired for lake points
!     (18 = mixed woodland (no lake model);  14 = inland water (lake
!     model)
      read (idata,99001) lakeset
      print * , '*** lake points set to bats surface type ' , lakeset
 
!     initialize data for lake model
      do n = 1 , lkpts
        read (idata,99002) ilake , jlake , depth
        freeze = 1
 
!       ******     Initially set eta to mean values
        if ( depth.lt.50 ) then
          eta = .7
        else if ( depth.gt.100 ) then
          eta = .3
        else
          eta = .5
        end if
 
!       ******     Initially set lake points isothermal at 6.0 C (June)
        do j = 1 , depth
          t(j) = 6.0
        end do
 
        veg2d(ilake,jlake) = dble(lakeset)
        write (iin) ilake , jlake , depth , freeze , hi , hice , hsnow ,&
                  & eta , (t(j),j=1,depth)
        print * , ilake , jlake , depth , freeze
 
      end do
 
      rewind (iin)
99001 format (i2)
99002 format (3I4)
 
      end subroutine initlk
