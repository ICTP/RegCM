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

      subroutine fexistnew(filnam,there)
      implicit none
!
! Dummy arguments
!
      character(256) :: filnam
      logical :: there
      intent (inout) filnam , there
!
! Local variables
!
      character(1) :: yesno
!
 100  continue
      inquire (file=filnam,exist=there)
      if ( .not.there ) then
 150    continue
        print * , 'FILE CAN NOT BE OPENED BECAUSE IT DOES NOT EXISTS: ' &
            & , filnam
        print * , 'DO YOU WANT TO CONTINUE? (y/n)'
        read (*,*) yesno
        if ( yesno=='y' ) then
          print * , 'ENTER NEW FILE NAME'
          read (*,'(a)') filnam
          go to 100
        else if ( yesno=='n' ) then
          return
        else
          print * , 'I DO NOT UNDERSTAND YOUR RESPONSE!!!'
          go to 150
        end if
      end if
      print * , 'OPEN NEW FILE:' , filnam
 
      end subroutine fexistnew
