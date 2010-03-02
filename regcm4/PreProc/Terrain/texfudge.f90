      subroutine texfudge(fudge,ch,texout,htgrid,iy,jx,char_tex)
      implicit none
!
! Dummy arguments
!
      character(10) :: char_tex
      logical :: fudge
      integer :: iy , jx
      character(1) , dimension(iy,jx) :: ch
      real(4) , dimension(iy,jx) :: htgrid , texout
      intent (in) char_tex , fudge , iy , jx
      intent (out) htgrid
      intent (inout) ch , texout
!
! Local variables
!
      integer :: i , j
!
      if ( fudge ) then
        open (13,file=char_tex,form='formatted')
        do i = iy , 1 , -1
          read (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
        do j = 1 , jx
          do i = 1 , iy
            if ( ch(i,j)==' ' ) then
              texout(i,j) = 14.
            else if ( ch(i,j)=='1' ) then
              texout(i,j) = 1.
            else if ( ch(i,j)=='2' ) then
              texout(i,j) = 2.
            else if ( ch(i,j)=='3' ) then
              texout(i,j) = 3.
            else if ( ch(i,j)=='4' ) then
              texout(i,j) = 4.
            else if ( ch(i,j)=='5' ) then
              texout(i,j) = 5.
            else if ( ch(i,j)=='6' ) then
              texout(i,j) = 6.
            else if ( ch(i,j)=='7' ) then
              texout(i,j) = 7.
            else if ( ch(i,j)=='8' ) then
              texout(i,j) = 8.
            else if ( ch(i,j)=='9' ) then
              texout(i,j) = 9.
            else if ( ch(i,j)=='A' ) then
              texout(i,j) = 10.
            else if ( ch(i,j)=='B' ) then
              texout(i,j) = 11.
            else if ( ch(i,j)=='C' ) then
              texout(i,j) = 12.
            else if ( ch(i,j)=='D' ) then
              texout(i,j) = 13.
            else if ( ch(i,j)=='E' ) then
              texout(i,j) = 14.
            else if ( ch(i,j)=='F' ) then
              texout(i,j) = 15.
            else if ( ch(i,j)=='G' ) then
              texout(i,j) = 16.
            else if ( ch(i,j)=='H' ) then
              texout(i,j) = 17.
            else if ( nint(texout(i,j))==0 ) then
!             ch(i,j) = 'X'
              ch(i,j) = ' '
            else
              write (*,*) 'TEXTURE TYPE exceed the limit'
              stop
            end if
            if ( nint(texout(i,j))==14 ) htgrid(i,j) = 0.0
          end do
        end do
      else
        do j = 1 , jx
          do i = 1 , iy
            if ( nint(texout(i,j))==14 ) then
              ch(i,j) = ' '
            else if ( nint(texout(i,j))==1 ) then
              ch(i,j) = '1'
            else if ( nint(texout(i,j))==2 ) then
              ch(i,j) = '2'
            else if ( nint(texout(i,j))==3 ) then
              ch(i,j) = '3'
            else if ( nint(texout(i,j))==4 ) then
              ch(i,j) = '4'
            else if ( nint(texout(i,j))==5 ) then
              ch(i,j) = '5'
            else if ( nint(texout(i,j))==6 ) then
              ch(i,j) = '6'
            else if ( nint(texout(i,j))==7 ) then
              ch(i,j) = '7'
            else if ( nint(texout(i,j))==8 ) then
              ch(i,j) = '8'
            else if ( nint(texout(i,j))==9 ) then
              ch(i,j) = '9'
            else if ( nint(texout(i,j))==10 ) then
              ch(i,j) = 'A'
            else if ( nint(texout(i,j))==11 ) then
              ch(i,j) = 'B'
            else if ( nint(texout(i,j))==12 ) then
              ch(i,j) = 'C'
            else if ( nint(texout(i,j))==13 ) then
              ch(i,j) = 'D'
            else if ( nint(texout(i,j))==15 ) then
              ch(i,j) = 'F'
            else if ( nint(texout(i,j))==16 ) then
              ch(i,j) = 'G'
            else if ( nint(texout(i,j))==17 ) then
              ch(i,j) = 'H'
            else
              write (*,*) 'TEXTURE TYPE' , nint(texout(i,j)) ,          &
                         &'exceed the limit'
              stop
            end if
          end do
        end do
        open (13,file=char_tex,form='formatted')
        do i = iy , 1 , -1
          write (13,99001) (ch(i,j),j=1,jx)
        end do
        close (13)
      end if
99001 format (132A1)
      end subroutine texfudge
