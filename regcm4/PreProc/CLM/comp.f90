      subroutine comp(fields,bvoc)
      implicit none
!
! Dummy arguments
!
      integer :: fields
      logical :: bvoc
      intent (in) bvoc
      intent (out) fields
!
! Local variables
!
      integer :: numcompounds
 
      if ( bvoc ) then
 50     continue
        print * , ' '
        print * , ' '
        print * , '********************************************'
        print * , 'Creating biogenic emissions files'
        print * , 'ENTER NUMBER OF SPECIES'
        read (*,*) numcompounds
        fields = 14 + numcompounds
        if ( fields>50 ) then
          stop 999
        else if ( fields>50 ) then
          go to 50
        else
        end if
      else
        fields = 14
      end if
 
      print * , 'producing ' , fields , ' files'
      end subroutine comp
