      subroutine fexist(filnam)
      implicit none
!
! Dummy arguments
!
      character(*) :: filnam
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
        print * , 'FILE ALREADY EXISTS:  ' , filnam
        print * , 'Do you want to overwrite the existing file? [y/n/q]'
        read (*,*) yesno
        if ( yesno=='y' ) then
          return
        else if ( yesno=='n' ) then
          print * , 'ENTER NEW FILE NAME'
          read (*,*) filnam
          go to 100
        else if ( yesno=='q' ) then
          stop 999
        else
          go to 150
        end if
      end if
 
      end subroutine fexist
