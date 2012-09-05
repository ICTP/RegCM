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

  private
!
  character(512) :: aline
  character(8) :: cline

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
!
    write (cline,'(i8)') line
    write (aline,*) '-------------- FATAL CALLED ---------------'
    call say
    if ( line > 0 ) then
      write (aline,*) 'Fatal in file: '//filename//' at line: '//trim(cline)
      call say
    end if
    write (aline,*) str
    call say
    write (aline,*) '-------------------------------------------'
    call say
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

end module mod_message
