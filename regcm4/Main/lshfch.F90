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
 
      subroutine lshfch(ixchar,lchar,ioff,iflg)
!
!****     this routine shifts the characters in ixchar ioff characters
!****     to the left.  the right side off ixchar is either blank
!****     or cnull filled depending on iflg.  for example, if ixchar
!****     is to start with the 3rd character instead of the 1st ioff
!****     would be set to 2 shifting ixchar to the left.
!
!****     on input-
!****         ixchar  - character input data. (character)
!****         lchar  - number of characters in ixchar. (integer)
!****         ioff   - number of characters ioff is to be offset to
!****                  left. (integer)
!****         iflg   - flag specifying type of fill to be used for
!****                  the ioff rightmost characters in ixchar.  if
!****                      0 - blank fill is used
!****                      1 - cnull fill is used
!
!****     on output-
!****         ixchar  - ixchar is retured shifted ioff characters to the
!****                  left with the ioff rightmost characters filled
!****                  in accordance with iflg.
!
      implicit none
!
! Dummy arguments
!
      integer :: iflg , ioff , lchar
      character(*) :: ixchar
      intent (in) iflg , ioff , lchar
      intent (inout) ixchar
!
! Local variables
!
      integer :: cnull , i , iblnk , ifill , ipc
      character(1) :: temp
!
      data iblnk/32/
      data cnull/0/
!
      ipc = 1
      do i = ioff + 1 , lchar
        temp = ixchar(i:i)
        ixchar(ipc:ipc) = temp
        ipc = ipc + 1
      end do
!
      if ( iflg.eq.0 ) then
        ifill = iblnk
      else if ( iflg.eq.1 ) then
        ifill = cnull
      else
        call fatal(__FILE__,__LINE__,'iflg out of range')
      end if
!
      do i = ipc , lchar
        ixchar(i:i) = char(ifill)
      end do
 
      end subroutine lshfch
