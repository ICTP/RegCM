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
 
      function lenchr(chrstg)

!----------------------------------------------------------------------
!
!l         return position of right-most non-blank, non-null character
!l         in chrstg.
!
!----------------------------------------------------------------------
      implicit none
!
! Dummy arguments
!
      character(*) :: chrstg
      integer :: lenchr
      intent (in) chrstg
!
      lenchr = len(chrstg)
      if ( chrstg(lenchr:lenchr).eq.' ' .or. chrstg(lenchr:lenchr)      &
         & .eq.char(0) ) then
        do
          lenchr = lenchr - 1
          if ( chrstg(lenchr:lenchr).eq.' ' .or. chrstg(lenchr:lenchr)  &
             & .eq.char(0) ) then
            if ( lenchr.gt.0 ) cycle
          end if
          exit
        end do
      end if
      end function lenchr
