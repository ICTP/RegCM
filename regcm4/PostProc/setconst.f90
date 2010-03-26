      subroutine setconst(f,const,n1,n2,n3,n4,n5,i1,ni1,i2,ni2)
 
      implicit none
!
! Dummy arguments
!
      real(4) :: const
      integer :: i1 , i2 , n1 , n2 , n3 , n4 , n5 , ni1 , ni2
      real(4) , dimension(n1,n2,n3,n4,n5) :: f
      intent (in) const , i1 , i2 , n1 , n2 , n3 , n4 , n5 , ni1 , ni2
      intent (out) f
!
! Local variables
!
      integer :: i , j , k , l , m
!
      do m = 1 , n5
        do l = 1 , n4
          do k = 1 , n3
            do j = i2 , ni2
              do i = i1 , ni1
                f(i,j,k,l,m) = const
              end do
            end do
          end do
        end do
      end do
 
      end subroutine setconst
