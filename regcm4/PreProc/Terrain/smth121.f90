      subroutine smth121(htgrid,iy,jx,hscr1)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx
      real(4) , dimension(iy,jx) :: hscr1 , htgrid
      intent (in) iy , jx
      intent (inout) hscr1 , htgrid
!
! Local variables
!
      integer :: i , j
!
!     PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
!     2DX WAVES FROM THE FIELDS htgrid
!
      do j = 1 , jx
        do i = 1 , iy
          hscr1(i,j) = htgrid(i,j)
        end do
      end do
      do i = 1 , iy
        do j = 2 , jx - 1
          if ( (htgrid(i,j)<=-.1) .or. (htgrid(i,j)>0.) ) hscr1(i,j)    &
             & = .25*(2.*htgrid(i,j)+htgrid(i,j+1)+htgrid(i,j-1))
        end do
      end do
      do j = 1 , jx
        do i = 2 , iy - 1
          if ( (hscr1(i,j)<=-.1) .or. (hscr1(i,j)>0.) ) htgrid(i,j)     &
             & = .25*(2.*hscr1(i,j)+hscr1(i+1,j)+hscr1(i-1,j))
        end do
      end do
      end subroutine smth121
