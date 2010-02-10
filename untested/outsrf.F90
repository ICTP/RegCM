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
 
      subroutine outsrf
 
! ******    write bats fields to unit iutbat
 
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_date
      use mod_iunits
      use mod_bats
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Local variables
!
      integer :: i , j , n
!
!     ****** check if at desired output time for bats variables
      write (*,*) 'BATS variables written at ' , idatex , xtime
      call fillbat
      if ( iotyp.eq.2 ) write (iutbat) idatex
      do n = 1 , numbat
        if ( iotyp.eq.1 ) then
          nrcbat = nrcbat + 1
#ifdef MPP1
          write (iutbat,rec=nrcbat)                                     &
                & ((fbat_io(j,i,n),j=1,jxm2),i=1,ixm2)
#else
          write (iutbat,rec=nrcbat) ((fbat(j,i,n),j=1,jxm2),i=1,ixm2)
#endif
        else if ( iotyp.eq.2 ) then
#ifdef MPP1
          write (iutbat) ((fbat_io(j,i,n),j=1,jxm2),i=1,ixm2)
#else
          write (iutbat) ((fbat(j,i,n),j=1,jxm2),i=1,ixm2)
#endif
        else
        end if
      end do
 
      end subroutine outsrf
