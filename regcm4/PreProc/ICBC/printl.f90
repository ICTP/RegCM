      subroutine printl(a,n1,n2)
      implicit none
!
! Dummy arguments
!
      integer :: n1 , n2
      real , dimension(n1,n2) :: a
      intent (in) a , n1 , n2
!
! Local variables
!
      integer :: inc1 , inc2 , k , l
!
      inc2 = -1
      inc1 = 1
 
      do l = n2 , 1 , inc2
        print 99001 , l , (a(k,l),k=1,17,inc1)
      end do
99001 format (1x,i3,2x,17F7.2)
 
      end subroutine printl
