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

      use mod_dynparam , only : myid

      implicit none
!
      character(512) :: aline

      contains

      subroutine say
      implicit none
      if ( myid == 0 ) print * , trim(aline)
      end subroutine say
 
      subroutine warning
      implicit none
      print * , ' Processor ' , myid , ' : ' , trim(aline)
      end subroutine warning
 
      subroutine fatal(filename,line,str)
      use mpi
      implicit none
!
      character(*) :: filename , str
      integer :: line
      intent (in) filename , line , str
!
      character(8) :: cline
      integer :: ierr
!
      write (cline,'(i6)') line
      write (aline,*) '-------------- FATAL CALLED ---------------'
      call say
      if ( line > 0 ) then
        write (aline,*) 'Fatal in file: '//filename//' at line: '//     &
                      & trim(cline)
        call say
      end if
      write (aline,*) str
      call say
      write (aline,*) '-------------------------------------------'
      call say
      call mpi_abort(mpi_comm_world,1,ierr)
      end subroutine fatal

      end module mod_message
