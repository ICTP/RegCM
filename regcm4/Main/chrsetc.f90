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
 
      subroutine chrsetc(iichar,ilen,isetch)
      implicit none
!
! Dummy arguments
!
      character(*) :: iichar
      integer :: ilen
      character(1) :: isetch
      intent (in) ilen , isetch
      intent (out) iichar
!
! Local variables
!
      integer :: i
!
!**** this routine sets to the first ilen characters in iichar
!**** to the character variable specified by isetch
!
!**** on input-
!**** ilen    - number of characters to be set in iichar (integer)
!**** isetch - character to which each character in iichar will
!**** be set equal to if input is character. (char)
!
!**** on output-
!**** iichar  - ilen charcters in iichar is set to the character
!**** specified by either isetch.  (character)
!
!
      do i = 1 , ilen
        iichar(i:i) = isetch
      end do
 
      end subroutine chrsetc
