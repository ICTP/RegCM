      subroutine two2three(i,nx,ny,kx,fin,fout)
 
      implicit none
!
! Dummy arguments
!
      integer :: i , nx , ny , kx
      real(4) , dimension(ny,kx) :: fin
      real(4) , dimension(nx,ny,kx) :: fout
      intent (in) i , fin , nx , ny , kx
      intent (out) fout
!
! Local variables
!
      integer :: j , k
!
      do k = 1 , kx
        do j = 1 , ny
          fout(i,j,k) = fin(j,k)
        end do
      end do
 
      end subroutine two2three
