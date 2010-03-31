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
 
      subroutine finddate(npos,idate)
      use mod_regcm_param
      use mod_message
      implicit none
!
! Dummy arguments
!
      integer :: idate , npos
      intent (in) idate
      intent (inout) npos
!
! Local variables
!
      integer :: i
!
      i = 0
      do
        i = i + 1
        if ( mdatez(i).eq.idate ) then
          npos = i
          exit
        end if
        if ( i.gt.289276 ) then
 
          write (aline,*) 'error in finddate : ' , npos , idate
          call say
          call fatal(__FILE__,__LINE__,'Date not found')
          exit
        end if
      end do
      end subroutine finddate
