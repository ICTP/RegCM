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
