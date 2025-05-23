!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_fudge

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message

  private

  public :: lndfudge, texfudge, lakfudge

  character(len=*), parameter :: f99001 = '(132A1)'

  contains

  subroutine lndfudge(fudge,lndout,jx,iy,char_lnd)
    implicit none
    character(len=*) :: char_lnd
    logical :: fudge, there
    integer(ik4) :: iy, jx
    real(rkx), dimension(jx,iy) :: lndout
    intent (in) char_lnd, fudge, iy, jx
    intent (inout) lndout
    integer(ik4) :: iunit
    integer(ik4) :: i, j
    character(len=1), dimension(jx,iy) :: ch

    if ( fudge ) then
      inquire (file=char_lnd,exist=there)
      if ( .not.there ) then
        write(stderr,*) 'Fudging requested for landuse but '// &
                 ' missing input ascii file ',trim(char_lnd)
        write(stderr,*)  'ERROR OPENING ', char_lnd,  &
            ' FILE:  FILE DOES NOT EXIST'
        call die('lndfudge')
      endif
      open (newunit=iunit,file=char_lnd,form='formatted')
      do i = iy, 1, -1
        read (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
      do i = 1, iy
        do j = 1, jx
          if ( ch(j,i) == ' ' ) then
            lndout(j,i) = 15.
          else if ( ch(j,i) == '1' ) then
            lndout(j,i) = 1.
          else if ( ch(j,i) == '2' ) then
            lndout(j,i) = 2.
          else if ( ch(j,i) == '3' ) then
            lndout(j,i) = 3.
          else if ( ch(j,i) == '4' ) then
            lndout(j,i) = 4.
          else if ( ch(j,i) == '5' ) then
            lndout(j,i) = 5.
          else if ( ch(j,i) == '6' ) then
            lndout(j,i) = 6.
          else if ( ch(j,i) == '7' ) then
            lndout(j,i) = 7.
          else if ( ch(j,i) == '8' ) then
            lndout(j,i) = 8.
          else if ( ch(j,i) == '9' ) then
            lndout(j,i) = 9.
          else if ( ch(j,i) == 'A' ) then
            lndout(j,i) = 10.
          else if ( ch(j,i) == 'B' ) then
            lndout(j,i) = 11.
          else if ( ch(j,i) == 'C' ) then
            lndout(j,i) = 12.
          else if ( ch(j,i) == 'D' ) then
            lndout(j,i) = 13.
          else if ( ch(j,i) == 'E' ) then
            lndout(j,i) = 14.
          else if ( ch(j,i) == 'F' ) then
            lndout(j,i) = 15.
          else if ( ch(j,i) == 'G' ) then
            lndout(j,i) = 16.
          else if ( ch(j,i) == 'H' ) then
            lndout(j,i) = 17.
          else if ( ch(j,i) == 'I' ) then
            lndout(j,i) = 18.
          else if ( ch(j,i) == 'J' ) then
            lndout(j,i) = 19.
          else if ( ch(j,i) == 'K' ) then
            lndout(j,i) = 20.
          else if ( ch(j,i) == 'L' ) then
            lndout(j,i) = 21.
          else if ( ch(j,i) == 'M' ) then
            lndout(j,i) = 22.
          else if ( nint(lndout(j,i)) == 0 ) then
            ch(j,i) = ' '
          else
            write (stderr,*) 'LANDUSE MASK exceed the limit'
            call die('lndfudge')
          end if
        end do
      end do
    else
      do i = 1, iy
        do j = 1, jx
          if ( nint(lndout(j,i)) == 15 .or. nint(lndout(j,i)) == 0 ) then
            ch(j,i) = ' '
          else if ( nint(lndout(j,i)) == 1 ) then
            ch(j,i) = '1'
          else if ( nint(lndout(j,i)) == 2 ) then
            ch(j,i) = '2'
          else if ( nint(lndout(j,i)) == 3 ) then
            ch(j,i) = '3'
          else if ( nint(lndout(j,i)) == 4 ) then
            ch(j,i) = '4'
          else if ( nint(lndout(j,i)) == 5 ) then
            ch(j,i) = '5'
          else if ( nint(lndout(j,i)) == 6 ) then
            ch(j,i) = '6'
          else if ( nint(lndout(j,i)) == 7 ) then
            ch(j,i) = '7'
          else if ( nint(lndout(j,i)) == 8 ) then
            ch(j,i) = '8'
          else if ( nint(lndout(j,i)) == 9 ) then
            ch(j,i) = '9'
          else if ( nint(lndout(j,i)) == 10 ) then
            ch(j,i) = 'A'
          else if ( nint(lndout(j,i)) == 11 ) then
            ch(j,i) = 'B'
          else if ( nint(lndout(j,i)) == 12 ) then
            ch(j,i) = 'C'
          else if ( nint(lndout(j,i)) == 13 ) then
            ch(j,i) = 'D'
          else if ( nint(lndout(j,i)) == 14 ) then
            ch(j,i) = 'E'
          else if ( nint(lndout(j,i)) == 16 ) then
            ch(j,i) = 'G'
          else if ( nint(lndout(j,i)) == 17 ) then
            ch(j,i) = 'H'
          else if ( nint(lndout(j,i)) == 18 ) then
            ch(j,i) = 'I'
          else if ( nint(lndout(j,i)) == 19 ) then
            ch(j,i) = 'J'
          else if ( nint(lndout(j,i)) == 20 ) then
            ch(j,i) = 'K'
          else if ( nint(lndout(j,i)) == 21 ) then
            ch(j,i) = 'L'
          else if ( nint(lndout(j,i)) == 22 ) then
            ch(j,i) = 'M'
          else
            write (stderr,*) 'LANDUSE MASK', nint(lndout(j,i)),        &
                        'exceed the limit'
            call die('lndfudge')
          end if
        end do
      end do
      open (newunit=iunit,file=char_lnd,form='formatted',err=100)
      do i = iy, 1, -1
        write (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
    end if
    return
100 write(stderr, *) 'Cannot create file ',trim(char_lnd)
    write (stderr, *)  'Is the directory "dirter" present?'
    write (stderr, *)  'Have you got write privileges on it?'
    call die('lndfudge','Path or permission problem',1)
  end subroutine lndfudge

  subroutine texfudge(fudge,texout,lnduse,jx,iy,char_tex)
    implicit none
    character(len=*) :: char_tex
    logical :: fudge, there
    integer(ik4) :: iy, jx
    real(rkx), dimension(jx,iy) :: texout, lnduse
    intent (in) char_tex, fudge, iy, jx
    intent (inout) texout, lnduse
    integer(ik4) :: i, j
    integer(ik4) :: iunit
    real(rkx) :: oval
    character(len=1), dimension(jx,iy) :: ch

    if ( fudge ) then
      inquire (file=char_tex,exist=there)
      if ( .not.there ) then
        write(stderr,*) 'Fudging requested for texture but '// &
                 'missing input ascii file ',trim(char_tex)
        write(stderr,*)  'ERROR OPENING ', char_tex,   &
                 ' FILE:  FILE DOES NOT EXIST'
        call die('texfudge')
      endif
      open (newunit=iunit,file=char_tex,form='formatted')
      do i = iy, 1, -1
        read (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
      do i = 1, iy
        do j = 1, jx
          oval = texout(j,i)
          if ( ch(j,i)==' ' ) then
            texout(j,i) = 14.
          else if ( ch(j,i) == '1' ) then
            texout(j,i) = 1.
          else if ( ch(j,i) == '2' ) then
            texout(j,i) = 2.
          else if ( ch(j,i) == '3' ) then
            texout(j,i) = 3.
          else if ( ch(j,i) == '4' ) then
            texout(j,i) = 4.
          else if ( ch(j,i) == '5' ) then
            texout(j,i) = 5.
          else if ( ch(j,i) == '6' ) then
            texout(j,i) = 6.
          else if ( ch(j,i) == '7' ) then
            texout(j,i) = 7.
          else if ( ch(j,i) == '8' ) then
            texout(j,i) = 8.
          else if ( ch(j,i) == '9' ) then
            texout(j,i) = 9.
          else if ( ch(j,i) == 'A' ) then
            texout(j,i) = 10.
          else if ( ch(j,i) == 'B' ) then
            texout(j,i) = 11.
          else if ( ch(j,i) == 'C' ) then
            texout(j,i) = 12.
          else if ( ch(j,i) == 'D' ) then
            texout(j,i) = 13.
          else if ( ch(j,i) == 'E' ) then
            texout(j,i) = 14.
          else if ( ch(j,i) == 'F' ) then
            texout(j,i) = 15.
          else if ( ch(j,i) == 'G' ) then
            texout(j,i) = 16.
          else if ( ch(j,i) == 'H' ) then
            texout(j,i) = 17.
          else if ( nint(texout(j,i)) == 0 ) then
            ch(j,i) = ' '
          else
            write (stderr,*) 'TEXTURE TYPE exceed the limit'
            call die('texfudge')
          end if
          if ( nint(oval) /= nint(texout(j,i)) ) then
            if ( nint(texout(j,i)) == 14 ) then
              if ( lnduse(j,i) < 13.5 .or. lnduse(j,i) > 15.5 ) then
                lnduse(j,i) = 14.0
              end if
            else
              if ( lnduse(j,i) > 13.5 .and. lnduse(j,i) < 15.5 ) then
                lnduse(j,i) = 8.0
              end if
            end if
          end if
        end do
      end do
    else
      do i = 1, iy
        do j = 1, jx
          if ( nint(texout(j,i)) == 14 ) then
            ch(j,i) = ' '
          else if ( nint(texout(j,i)) == 1 ) then
            ch(j,i) = '1'
          else if ( nint(texout(j,i)) == 2 ) then
            ch(j,i) = '2'
          else if ( nint(texout(j,i)) == 3 ) then
            ch(j,i) = '3'
          else if ( nint(texout(j,i)) == 4 ) then
            ch(j,i) = '4'
          else if ( nint(texout(j,i)) == 5 ) then
            ch(j,i) = '5'
          else if ( nint(texout(j,i)) == 6 ) then
            ch(j,i) = '6'
          else if ( nint(texout(j,i)) == 7 ) then
            ch(j,i) = '7'
          else if ( nint(texout(j,i)) == 8 ) then
            ch(j,i) = '8'
          else if ( nint(texout(j,i)) == 9 ) then
            ch(j,i) = '9'
          else if ( nint(texout(j,i)) == 10 ) then
            ch(j,i) = 'A'
          else if ( nint(texout(j,i)) == 11 ) then
            ch(j,i) = 'B'
          else if ( nint(texout(j,i)) == 12 ) then
            ch(j,i) = 'C'
          else if ( nint(texout(j,i)) == 13 ) then
            ch(j,i) = 'D'
          else if ( nint(texout(j,i)) == 15 ) then
            ch(j,i) = 'F'
          else if ( nint(texout(j,i)) == 16 ) then
            ch(j,i) = 'G'
          else if ( nint(texout(j,i)) == 17 ) then
            ch(j,i) = 'H'
          else
            write (stderr,*) 'TEXTURE TYPE', nint(texout(j,i)),          &
                        'exceed the limit'
            call die('texfudge')
          end if
        end do
      end do
      open (newunit=iunit,file=char_tex,form='formatted',err=100)
      do i = iy, 1, -1
        write (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
    end if
    return
100 write(stderr, *) 'Cannot create file ',trim(char_tex)
    write (stderr, *)  'Is the directory "dirter" present?'
    write (stderr, *)  'Have you got write privileges on it?'
    call die('lndfudge','Path or permission problem',1)
  end subroutine texfudge

  subroutine lakfudge(fudge,dpth,lnd,jx,iy,char_lak)
    implicit none
    character(len=*), intent(in) :: char_lak
    logical, intent(in) :: fudge
    integer(ik4), intent(in) :: iy, jx
    real(rkx), dimension(jx,iy), intent(inout) :: dpth, lnd
    logical :: there
    integer(ik4) :: i, j
    integer(ik4) :: iunit
    character(len=1), dimension(jx,iy) :: ch

    if ( fudge ) then
      inquire (file=char_lak,exist=there)
      if ( .not.there ) then
        write(stderr,*) 'Fudging requested for lake but '// &
                 'missing input ascii file ',trim(char_lak)
        write(stderr,*)  'ERROR OPENING ', char_lak,  &
            ' FILE:  FILE DOES NOT EXIST'
        call die('lakfudge')
      endif
      open (newunit=iunit,file=char_lak,form='formatted')
      do i = iy, 1, -1
        read (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
      do i = 1, iy
        do j = 1, jx
          if (lnd(j,i) > 13.5 .and. lnd(j,i) < 14.5) then
            if ( ch(j,i)/='L' ) then
              dpth(j,i) = 0.0
              lnd(j,i) = 15.0
             end if
           end if
        end do
      end do
    else
      do i = 1, iy
        do j = 1, jx
          if ( dpth(j,i) > 0.0 ) then
            if (lnd(j,i) > 13.5 .and. lnd(j,i) < 14.5) then
              ch(j,i) = 'L'
            else
              ch(j,i) = '*'
            end if
          else
            ch(j,i) = '-'
          end if
        end do
      end do
      open (newunit=iunit,file=char_lak,form='formatted',err=100)
      do i = iy, 1, -1
        write (iunit,f99001) (ch(j,i),j=1,jx)
      end do
      close (iunit)
    end if
    return
100 write(stderr, *) 'Cannot create file ',trim(char_lak)
    write (stderr, *)  'Is the directory "dirter" present?'
    write (stderr, *)  'Have you got write privileges on it?'
    call die('lndfudge','Path or permission problem',1)
  end subroutine lakfudge

end module mod_fudge
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
