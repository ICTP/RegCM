!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine lndfudge(fudge,ch,lndout,htgrid,ix,jx,lsmtyp,char_lnd)
      implicit none
!
! Dummy arguments
!
      character(10) :: char_lnd
      logical :: fudge
      integer :: ix , jx
      character(4) :: lsmtyp
      character(1) , dimension(ix,jx) :: ch
      real(4) , dimension(ix,jx) :: htgrid , lndout
      intent (in) char_lnd , fudge , ix , jx , lsmtyp
      intent (inout) ch , htgrid , lndout
!
! Local variables
!
      integer :: i , j
!
      if ( lsmtyp=='BATS' ) then
        if ( fudge ) then
          open (13,file=char_lnd,form='formatted')
          do i = ix , 1 , -1
            read (13,99001) (ch(i,j),j=1,jx)
          end do
          close (13)
          do j = 1 , jx
            do i = 1 , ix
              if ( ch(i,j)==' ' ) then
                lndout(i,j) = 15.
              else if ( ch(i,j)=='1' ) then
                lndout(i,j) = 1.
              else if ( ch(i,j)=='2' ) then
                lndout(i,j) = 2.
              else if ( ch(i,j)=='3' ) then
                lndout(i,j) = 3.
              else if ( ch(i,j)=='4' ) then
                lndout(i,j) = 4.
              else if ( ch(i,j)=='5' ) then
                lndout(i,j) = 5.
              else if ( ch(i,j)=='6' ) then
                lndout(i,j) = 6.
              else if ( ch(i,j)=='7' ) then
                lndout(i,j) = 7.
              else if ( ch(i,j)=='8' ) then
                lndout(i,j) = 8.
              else if ( ch(i,j)=='9' ) then
                lndout(i,j) = 9.
              else if ( ch(i,j)=='A' ) then
                lndout(i,j) = 10.
              else if ( ch(i,j)=='B' ) then
                lndout(i,j) = 11.
              else if ( ch(i,j)=='C' ) then
                lndout(i,j) = 12.
              else if ( ch(i,j)=='D' ) then
                lndout(i,j) = 13.
              else if ( ch(i,j)=='E' ) then
                lndout(i,j) = 14.
              else if ( ch(i,j)=='F' ) then
                lndout(i,j) = 15.
              else if ( ch(i,j)=='G' ) then
                lndout(i,j) = 16.
              else if ( ch(i,j)=='H' ) then
                lndout(i,j) = 17.
              else if ( ch(i,j)=='I' ) then
                lndout(i,j) = 18.
              else if ( ch(i,j)=='J' ) then
                lndout(i,j) = 19.
              else if ( ch(i,j)=='K' ) then
                lndout(i,j) = 20.
              else if ( nint(lndout(i,j))==0 ) then
!               ch(i,j) = 'X'
                ch(i,j) = ' '
              else
                write (*,*) 'LANDUSE MASK exceed the limit'
                stop
              end if
!_fix         if(nint(lndout(i,j)).eq.15) htgrid(i,j) = 0.0
              if ( htgrid(i,j)<0.1 .and. nint(lndout(i,j))==15 )        &
                 & htgrid(i,j) = 0.0
            end do
          end do
        else
          do j = 1 , jx
            do i = 1 , ix
              if ( nint(lndout(i,j))==15 ) then
                ch(i,j) = ' '
              else if ( nint(lndout(i,j))==1 ) then
                ch(i,j) = '1'
              else if ( nint(lndout(i,j))==2 ) then
                ch(i,j) = '2'
              else if ( nint(lndout(i,j))==3 ) then
                ch(i,j) = '3'
              else if ( nint(lndout(i,j))==4 ) then
                ch(i,j) = '4'
              else if ( nint(lndout(i,j))==5 ) then
                ch(i,j) = '5'
              else if ( nint(lndout(i,j))==6 ) then
                ch(i,j) = '6'
              else if ( nint(lndout(i,j))==7 ) then
                ch(i,j) = '7'
              else if ( nint(lndout(i,j))==8 ) then
                ch(i,j) = '8'
              else if ( nint(lndout(i,j))==9 ) then
                ch(i,j) = '9'
              else if ( nint(lndout(i,j))==10 ) then
                ch(i,j) = 'A'
              else if ( nint(lndout(i,j))==11 ) then
                ch(i,j) = 'B'
              else if ( nint(lndout(i,j))==12 ) then
                ch(i,j) = 'C'
              else if ( nint(lndout(i,j))==13 ) then
                ch(i,j) = 'D'
              else if ( nint(lndout(i,j))==14 ) then
                ch(i,j) = 'E'
              else if ( nint(lndout(i,j))==16 ) then
                ch(i,j) = 'G'
              else if ( nint(lndout(i,j))==17 ) then
                ch(i,j) = 'H'
              else if ( nint(lndout(i,j))==18 ) then
                ch(i,j) = 'I'
              else if ( nint(lndout(i,j))==19 ) then
                ch(i,j) = 'J'
              else if ( nint(lndout(i,j))==20 ) then
                ch(i,j) = 'K'
              else
                write (*,*) 'LANDUSE MASK' , nint(lndout(i,j)) ,        &
                           &'exceed the limit'
                stop
              end if
            end do
          end do
          open (13,file=char_lnd,form='formatted')
          do i = ix , 1 , -1
            write (13,99001) (ch(i,j),j=1,jx)
          end do
          close (13)
        end if
      else if ( lsmtyp=='USGS' ) then
        if ( fudge ) then
          open (13,file=char_lnd,form='formatted')
          do i = ix , 1 , -1
            read (13,99001) (ch(i,j),j=1,jx)
          end do
          close (13)
          do j = 1 , jx
            do i = 1 , ix
              if ( ch(i,j)==' ' ) then
                lndout(i,j) = 25.
              else if ( ch(i,j)=='1' ) then
                lndout(i,j) = 1.
              else if ( ch(i,j)=='2' ) then
                lndout(i,j) = 2.
              else if ( ch(i,j)=='3' ) then
                lndout(i,j) = 3.
              else if ( ch(i,j)=='4' ) then
                lndout(i,j) = 4.
              else if ( ch(i,j)=='5' ) then
                lndout(i,j) = 5.
              else if ( ch(i,j)=='6' ) then
                lndout(i,j) = 6.
              else if ( ch(i,j)=='7' ) then
                lndout(i,j) = 7.
              else if ( ch(i,j)=='8' ) then
                lndout(i,j) = 8.
              else if ( ch(i,j)=='9' ) then
                lndout(i,j) = 9.
              else if ( ch(i,j)=='A' ) then
                lndout(i,j) = 10.
              else if ( ch(i,j)=='B' ) then
                lndout(i,j) = 11.
              else if ( ch(i,j)=='C' ) then
                lndout(i,j) = 12.
              else if ( ch(i,j)=='D' ) then
                lndout(i,j) = 13.
              else if ( ch(i,j)=='E' ) then
                lndout(i,j) = 14.
              else if ( ch(i,j)=='F' ) then
                lndout(i,j) = 15.
              else if ( ch(i,j)=='G' ) then
                lndout(i,j) = 16.
              else if ( ch(i,j)=='H' ) then
                lndout(i,j) = 17.
              else if ( ch(i,j)=='I' ) then
                lndout(i,j) = 18.
              else if ( ch(i,j)=='J' ) then
                lndout(i,j) = 19.
              else if ( ch(i,j)=='K' ) then
                lndout(i,j) = 20.
              else if ( ch(i,j)=='L' ) then
                lndout(i,j) = 21.
              else if ( ch(i,j)=='M' ) then
                lndout(i,j) = 22.
              else if ( ch(i,j)=='N' ) then
                lndout(i,j) = 23.
              else if ( ch(i,j)=='O' ) then
                lndout(i,j) = 24.
              else
                write (*,*) 'LANDUSE MASK exceed the limit'
                stop
              end if
              if ( htgrid(i,j)<0.1 .and.                                &
                 & (nint(lndout(i,j))==25 .or. nint(lndout(i,j))==0) )  &
                 & htgrid(i,j) = 0.0
            end do
          end do
        else
          do j = 1 , jx
            do i = 1 , ix
              if ( nint(lndout(i,j))==25 .or. nint(lndout(i,j))==0 )    &
                 & then
                ch(i,j) = ' '
              else if ( nint(lndout(i,j))==1 ) then
                ch(i,j) = '1'
              else if ( nint(lndout(i,j))==2 ) then
                ch(i,j) = '2'
              else if ( nint(lndout(i,j))==3 ) then
                ch(i,j) = '3'
              else if ( nint(lndout(i,j))==4 ) then
                ch(i,j) = '4'
              else if ( nint(lndout(i,j))==5 ) then
                ch(i,j) = '5'
              else if ( nint(lndout(i,j))==6 ) then
                ch(i,j) = '6'
              else if ( nint(lndout(i,j))==7 ) then
                ch(i,j) = '7'
              else if ( nint(lndout(i,j))==8 ) then
                ch(i,j) = '8'
              else if ( nint(lndout(i,j))==9 ) then
                ch(i,j) = '9'
              else if ( nint(lndout(i,j))==10 ) then
                ch(i,j) = 'A'
              else if ( nint(lndout(i,j))==11 ) then
                ch(i,j) = 'B'
              else if ( nint(lndout(i,j))==12 ) then
                ch(i,j) = 'C'
              else if ( nint(lndout(i,j))==13 ) then
                ch(i,j) = 'D'
              else if ( nint(lndout(i,j))==14 ) then
                ch(i,j) = 'E'
              else if ( nint(lndout(i,j))==15 ) then
                ch(i,j) = 'F'
              else if ( nint(lndout(i,j))==16 ) then
                ch(i,j) = 'G'
              else if ( nint(lndout(i,j))==17 ) then
                ch(i,j) = 'H'
              else if ( nint(lndout(i,j))==18 ) then
                ch(i,j) = 'I'
              else if ( nint(lndout(i,j))==19 ) then
                ch(i,j) = 'J'
              else if ( nint(lndout(i,j))==20 ) then
                ch(i,j) = 'K'
              else if ( nint(lndout(i,j))==21 ) then
                ch(i,j) = 'L'
              else if ( nint(lndout(i,j))==22 ) then
                ch(i,j) = 'M'
              else if ( nint(lndout(i,j))==23 ) then
                ch(i,j) = 'N'
              else if ( nint(lndout(i,j))==24 ) then
                ch(i,j) = 'O'
              else
                write (*,*) 'LANDUSE MASK' , nint(lndout(i,j)) ,        &
                           &'exceed the limit'
                stop
              end if
            end do
          end do
          open (13,file=char_lnd,form='formatted')
          do i = ix , 1 , -1
            write (13,99001) (ch(i,j),j=1,jx)
          end do
          close (13)
        end if
      else
        print * , 'LANDUSE LEGEND DOES NOT EXIST'
        stop 'subroutine LNDFUDGE'
      end if
99001 format (132A1)
      end subroutine lndfudge
