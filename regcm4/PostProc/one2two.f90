      subroutine one2two(i,nx,ny,fin,fout)
 
      implicit none
!
! Dummy arguments
!
      integer :: i , nx , ny
      real(4) , dimension(ny) :: fin
      real(4) , dimension(nx,ny) :: fout
      intent (in) i , fin , nx , ny
      intent (out) fout
!
! Local variables
!
      integer :: j
!
      do j = 1 , ny
        fout(i,j) = fin(j)
      end do
 
      end subroutine one2two
