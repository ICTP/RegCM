!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_message

  use mod_stdio
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  character(len=512) :: aline
  character(len=8) :: cline

  character , parameter :: esc = achar(27)
  integer(ik4) , parameter :: maxplot = 10
  character , parameter :: plot(0:maxplot) = (/ &
      ' ', '.', ',', ':', '=', '+', 'o', 'x', 'X', '#', '@' /)

  integer(ik4) :: rows = 24
  integer(ik4) :: cols = 48
  integer(ik4) :: mrows , mcols

  public :: viz_init , viz_clear , viz_plot , viz_done

  public :: setup_mesg , die , aline , say , note , cry , fatal , checkalloc

  integer(ik4) :: iprank = 0

  interface die
    module procedure die0
    module procedure die1
    module procedure die2
    module procedure die4
  end interface die

  contains

  subroutine setup_mesg(ipid)
    implicit none
    integer(ik4) :: ipid
    iprank = ipid
  end subroutine setup_mesg

  subroutine say
    implicit none
    if ( iprank == 0 ) write (stdout,*) trim(aline)
  end subroutine say
 
  subroutine note
    implicit none
    write (aline,*) '------------------ NOTICE -----------------'
    write (stderr,*) ' Processor ' , iprank , ' : ' , trim(aline)
    write (aline,*) '-------------------------------------------'
  end subroutine note
 
  subroutine cry
    implicit none
    if ( iprank == 0 ) then
      write (aline,*) '------------- IMPORTANT NOTICE ------------'
      write (stderr,*) trim(aline)
      write (aline,*) '-------------------------------------------'
    end if
  end subroutine cry
 
  subroutine fatal(filename,line,str)
    implicit none
    character(*) , intent(in) :: filename , str
    integer(ik4) , intent(in) :: line
    write (cline,'(i8)') line
    write (stderr,*) '-------------- FATAL CALLED ---------------'
    if ( line > 0 ) then
      write (aline,*) 'Fatal in file: '//filename//' at line: '//trim(cline)
      write (stderr,*) trim(aline)
    end if
    write (stderr,*) str
    write (stderr,*) '-------------------------------------------'
    call die(filename,trim(cline),1)
  end subroutine fatal

  subroutine checkalloc(ival,filename,line,arg)
    implicit none
    integer(ik4) , intent(in) :: ival , line
    character(*) , intent(in) :: filename , arg
    if ( ival /= 0 ) then
      write (cline,'(i8)') line
      write (stderr,*) 'Memory error in allocating ', arg
      call die(filename,trim(cline),ival)
    end if
  end subroutine checkalloc

  subroutine die0(msg)
    implicit none
    character (len=*) , intent(in) :: msg
    external :: myabort
    if ( iprank == 0 ) write (stderr,*) msg
    call myabort
  end subroutine die0

  subroutine die1(msg,msg1)
    implicit none
    character (len=*) , intent(in) :: msg , msg1
    external :: myabort
    if ( iprank == 0 ) write (stderr,*) msg , ' : ', msg1
    call myabort
  end subroutine die1

  subroutine die2(msg,msg1,ier1)
    implicit none
    character (len=*) , intent(in) :: msg , msg1
    integer(ik4) , intent(in) :: ier1
    external :: myabort
    if ( iprank == 0 ) write (stderr,*) msg , ' : ', msg1 , ': ', ier1
    call myabort
  end subroutine die2

  subroutine die4(msg,msg1,ier1,msg2,ier2)
    implicit none
    character (len=*) , intent(in) :: msg , msg1 , msg2
    integer(ik4) , intent(in) :: ier1 , ier2
    external :: myabort
    if ( iprank == 0 ) write (stderr,*) msg , ' : ', msg1 , &
                           ': ', ier1 , ' : ', msg2 , ': ', ier2
    call myabort
  end subroutine die4

! small module to do text mode graphics using ANSI terminal escape sequences
! Copyright (c) 2013 Axel Kohlmeyer <akohlmey@gmail.com>

  ! initialize package and set size of plot area
  subroutine viz_init(x,y)
    integer(ik4) , intent(in) , optional :: x , y
    if ( present(x) ) cols = x
    if ( present(y) ) rows = y
    mcols = cols
    mrows = rows
    ! hide the cursor to reduce flicker
    write(stdout,fmt='(a1,a)',advance='no') esc,'[?25l'
  end subroutine viz_init

  ! restore settings to something sane.
  subroutine viz_done
    ! place cursor in last line
    call viz_pos(1,mrows)
    ! set forground color to principal color
    write(stdout,fmt='(a1,a,a1)',advance='no') esc,'[30m'
    ! re-enable the cursor
    write(stdout,fmt='(a1,a)', advance='no') esc,'[?25h'
    ! and call for a reset
    write(stdout,fmt='(a1,a)') esc,'[0m'
  end subroutine viz_done

  ! clear the screen
  subroutine viz_clear
    write(stdout,fmt='(a1,a)',advance='no') esc,'[2J'
  end subroutine viz_clear

  subroutine add_code(code,c,n)
    character(len=1) , intent(inout) , dimension(:) :: code
    character(len=1) , intent(in) :: c
    integer(ik4) , intent(inout) :: n
    n = n+1
    code(n) = c
  end subroutine add_code

  ! position cursor at a given location within screen.
  ! top left corner is (1,1)
  subroutine viz_pos(x,y)
    integer(ik4) , intent(in) :: x , y
    integer :: i , n
    character(len=1) :: code(7)
    n = 0

    i = y
    if (i < 1) i = 1
    if (i > mrows) i = mrows
    call add_code(code,'[',n)
    if (i > 9) call add_code(code,achar(48+i/10),n)
    call add_code(code,achar(48+mod(i,10)),n)
    call add_code(code,';',n)

    i = x
    if (i < 1) i = 1
    if (i > mcols) i = mcols
    if (i > 9) call add_code(code,achar(48+i/10),n)
    call add_code(code,achar(48+mod(i,10)),n)
    call add_code(code,'H',n)
    write(stdout,fmt='(a1,a)',advance='no') esc,code(1:n)
  end subroutine viz_pos

  subroutine viz_plot(val,rmax)
    real(rk8) , intent(in) , dimension(:,:) :: val
    real(rk8) , intent(in) , optional :: rmax
    integer(ik4) :: nx , ny , i , j , k , l , m , n , dx , dy
    real(rk8) :: vmax , scalef , tmp

    nx = size(val,1)
    ny = size(val,2)

    mcols = min(cols,nx)
    mrows = min(rows,ny)

    ! set blocksize for averaging
    dx = nx / mcols
    dy = ny / mrows

    ! set or determine scaling factor for data points
    vmax = 1.0D-30
    if ( present(rmax) ) then
      vmax = abs(rmax)
    else
      ! find absolute maximum value for scaling
      do j = 1 , mrows
        do i = 1 , mcols
          ! average over cells
          tmp = 0.0D0
          n = 0
          do k = (j-1)*dy+1 , j*dy
            do l = (i-1)*dx+1 , i*dx
              tmp = tmp + val(l,k)
              n = n + 1
            end do
          end do
          tmp = abs(tmp)/dble(n)
          if ( vmax < tmp ) vmax = tmp
        end do
      end do
    end if
    scalef = dble(maxplot)/vmax

    ! now plot
    do j = mrows , 1 , -1
      call viz_pos(1,mrows-j)
      do i = 1 , mcols
        ! average over cells
        tmp = 0.0D0
        n = 0
        do k = (j-1)*dy+1 , j*dy
          do l = (i-1)*dx+1 , i*dx
            tmp = tmp + val(l,k)
            n = n + 1
          end do
        end do
        ! convert absolute value into character
        m = int(scalef*abs(tmp)/dble(n)+0.5D0)
        m = max(0,min(maxplot,m))
        if ( tmp < 0.0D0 ) then
          if ( m > 5 ) then
            write(stdout,fmt='(a1,a,a1)',advance='no') esc,'[36m',plot(m)
          else
            write(stdout,fmt='(a1,a,a1)',advance='no') esc,'[34m',plot(m)
          end if
        else
          if ( m > 5 ) then
            write(stdout,fmt='(a1,a,a1)',advance='no') esc,'[31m',plot(m)
          else
            write(stdout,fmt='(a1,a,a1)',advance='no') esc,'[30m',plot(m)
          end if
        end if
      end do
      write(stdout,fmt='(a1,a)',advance='no') esc,'[30m|'
    end do
  end subroutine viz_plot

end module mod_message
