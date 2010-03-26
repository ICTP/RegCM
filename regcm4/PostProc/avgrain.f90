      subroutine avgrain(favg,f,f1,nx,ny,nhr,ihr)
 
      implicit none
!
! Dummy arguments
!
      integer :: ihr , nhr , nx , ny
      real(4) , dimension(nx,ny) :: f , f1
      real(4) , dimension(nx,ny,nhr) :: favg
      intent (in) f , ihr , nhr , nx , ny
      intent (inout) f1 , favg
!
! Local variables
!
      integer :: i , j
!
      do j = 1 , ny
        do i = 1 , nx
          favg(i,j,ihr) = favg(i,j,ihr) + (f(i,j)-f1(i,j))
          f1(i,j) = f(i,j)
        end do
      end do
 
      end subroutine avgrain
