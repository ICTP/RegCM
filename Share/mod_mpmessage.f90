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

module mod_mpmessage

  use mpi
  use mod_stdio
  use mod_dynparam , only : mycomm , myid
  use netcdf

  private
!
  character(512) :: aline
  character(8) :: cline

  public :: die , aline , say , note , cry , fatal , checkalloc , checkncerr

  interface die
    module procedure die0
    module procedure die1
    module procedure die2
    module procedure die4
  end interface die

  contains

  subroutine say
    implicit none
    if ( myid == 0 ) write (stdout,*) trim(aline)
  end subroutine say
 
  subroutine note
    implicit none
    write (aline,*) '------------------ NOTICE -----------------'
    write (stderr,*) ' Processor ' , myid , ' : ' , trim(aline)
    write (aline,*) '-------------------------------------------'
  end subroutine note
 
  subroutine cry
    implicit none
    if ( myid == 0 ) then
      write (aline,*) '------------- IMPORTANT NOTICE ------------'
      write (stderr,*) trim(aline)
      write (aline,*) '-------------------------------------------'
    end if
  end subroutine cry
 
  subroutine fatal(filename,line,str)
    use mpi
    implicit none
!
    character(*) , intent(in) :: filename , str
    integer , intent(in) :: line
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
    integer , intent(in) :: ival , line
    character(*) , intent(in) :: filename , arg
    if ( ival /= 0 ) then
      write (cline,'(i8)') line
      write (stderr,*) 'Memory error in allocating ', arg
      call die(filename,trim(cline),ival)
    end if
  end subroutine checkalloc

  subroutine checkncerr(ival,filename,line,arg)
    implicit none
    integer , intent(in) :: ival , line
    character(*) , intent(in) :: filename , arg
    if ( ival /= nf90_noerr ) then
      write (cline,'(i8)') line
      write (stderr,*) nf90_strerror(ival)
      call die(filename,trim(cline)//':'//arg,ival)
    end if
  end subroutine checkncerr

  subroutine die0(msg)
    implicit none
    character (len=*) , intent(in) :: msg
    integer :: ierr
    if ( myid == 0 ) write (stderr,*) msg
    call mpi_abort(mycomm,1,ierr)
  end subroutine die0

  subroutine die1(msg,msg1)
    implicit none
    character (len=*) , intent(in) :: msg , msg1
    integer :: ierr
    if ( myid == 0 ) write (stderr,*) msg , ' : ', msg1
    call mpi_abort(mycomm,1,ierr)
  end subroutine die1

  subroutine die2(msg,msg1,ier1)
    implicit none
    character (len=*) , intent(in) :: msg , msg1
    integer , intent(in) :: ier1
    integer :: ierr
    if ( myid == 0 ) write (stderr,*) msg , ' : ', msg1 , ': ', ier1
    call mpi_abort(mycomm,1,ierr)
  end subroutine die2

  subroutine die4(msg,msg1,ier1,msg2,ier2)
    implicit none
    character (len=*) , intent(in) :: msg , msg1 , msg2
    integer , intent(in) :: ier1 , ier2
    integer :: ierr
    if ( myid == 0 ) write (stderr,*) msg , ' : ', msg1 , &
                           ': ', ier1 , ' : ', msg2 , ': ', ier2
    call mpi_abort(mycomm,1,ierr)
  end subroutine die4

end module mod_mpmessage
