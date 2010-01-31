      module culookup
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: jptlucu1 = 50000 , jptlucu2 = 370000
!
! COMMON /LOOKUP/
!
      logical :: lookupoverflow
      real(8) , dimension(jptlucu1:jptlucu2) :: tlucua , tlucuaw ,      &
           & tlucub , tlucuc
      end module culookup
