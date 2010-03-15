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

      subroutine fexist(filnam)
      implicit none
!
! Dummy arguments
!
      character(256) :: filnam
      intent (inout) filnam
!
! Local variables
!
      logical :: there
      character(1) :: yesno

 100  continue
      inquire (file=filnam,exist=there)
      if ( there ) then
 150    continue
        print * , ' '
        print * , ' '
        print * , '**************************************************'
        print * , 'FILE ALREADY EXISTS:  ' , trim(filnam)
        print * , 'Do you want to overwrite the existing file? [y/n/q]'
        read (*,*) yesno
        if ( yesno=='y' ) then
          return
        else if ( yesno=='n' ) then
          print * , 'ENTER NEW FILE NAME'
          read (*,*) filnam
          goto 100
        else if ( yesno=='q' ) then
          stop 999
        else
          go to 150
        end if
      end if
 
      end subroutine fexist
