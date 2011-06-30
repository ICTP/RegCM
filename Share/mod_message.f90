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

  use m_die
  use m_stdio

  private
!
  character(512) :: aline
  character(8) :: cline

  public :: aline , say , note , cry , fatal , checkalloc

  contains

  subroutine say(myid)
    implicit none
    integer , intent(in) :: myid
    if ( myid == 0 ) write (stdout,*) trim(aline)
  end subroutine say
 
  subroutine note(myid)
    implicit none
    integer , intent(in) :: myid
    write (stderr,*) ' Processor ' , myid , ' : ' , trim(aline)
  end subroutine note
 
  subroutine cry(myid)
    implicit none
    integer , intent(in) :: myid
    if ( myid == 0 ) write (stderr,*) trim(aline)
  end subroutine cry
 
  subroutine fatal(filename,line,str)
    use mpi
    implicit none
!
    character(*) , intent(in) :: filename , str
    integer , intent(in) :: line
!
!
    write (cline,'(i8)') line
    write (aline,*) '-------------- FATAL CALLED ---------------'
    call cry(0)
    if ( line > 0 ) then
      write (aline,*) 'Fatal in file: '//filename//' at line: '//trim(cline)
      call cry(0)
    end if
    write (aline,*) str
    call cry(0)
    write (aline,*) '-------------------------------------------'
    call cry(0)
    call die(filename,trim(cline),1)
  end subroutine fatal

  subroutine checkalloc(ival,filename,line,arg)
    implicit none
    integer , intent(in) :: ival , line
    character(*) , intent(in) :: filename , arg
    write (cline,'(i8)') line
    if ( ival /= 0 ) then
      write (stderr,*) 'Memory error in allocating ', arg
      call die(filename,trim(cline),ival)
    end if
  end subroutine checkalloc

end module mod_message
