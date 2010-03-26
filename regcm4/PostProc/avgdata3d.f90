      subroutine avgdata3d(favg,f,n1,n2,n3,n4,n5,ihr,vmisdat)
 
      implicit none
!
! Dummy arguments
!
      integer :: ihr , n1 , n2 , n3 , n4 , n5
      real(4) :: vmisdat
      real(4) , dimension(n1,n2,n3,n4) :: f
      real(4) , dimension(n1,n2,n3,n4,n5) :: favg
      intent (in) f , ihr , n1 , n2 , n3 , n4 , n5 , vmisdat
      intent (inout) favg
!
! Local variables
!
      integer :: i , j , k , l
!
      do l = 1 , n4
        do k = 1 , n3
          do j = 1 , n2
            do i = 1 , n1
              if ( f(i,j,k,l)>vmisdat ) then
                favg(i,j,k,l,ihr) = favg(i,j,k,l,ihr) + f(i,j,k,l)
              else
                favg(i,j,k,l,ihr) = vmisdat
              end if
            end do
          end do
        end do
      end do
 
      end subroutine avgdata3d
