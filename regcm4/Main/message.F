C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      logical function active( level )
        implicit none
#       include "parameter.inc"
#       include "include/message.cb"
        integer , intent(in) :: level
        active = ( level .le. debug_level )
        return
      end function active

      subroutine say
        implicit none
#       include "parameter.inc"
#       include "include/message.cb"
#ifdef MPP1
        if ( myid .eq. 0 ) print*, trim(aline)
#else
        print*, trim(aline)
#endif
      end subroutine say

      subroutine warning
        implicit none
#       include "parameter.inc"
#       include "include/message.cb"
#ifdef MPP1
        print*, ' Processor ', myid,' : ',trim(aline)
#else
        print*, trim(aline)
#endif
      end subroutine warning

      subroutine fatal( filename, line, str )
        implicit none
#       include "parameter.inc"
#       include "include/message.cb"
#ifdef MPP1
        include 'mpif.h'
        INTEGER ierr
#endif
        character*(*) filename
        character*(*) str
        integer , intent (in) :: line
        character(len=8) :: cline

        write(cline,'(i6)') line
        write(aline, *) '-------------- FATAL CALLED ---------------'
        call say
        if ( line > 0 ) then
          write(aline, *) 'Fatal in file: '//filename//
     &                    ' at line: '//TRIM(cline)
          call say
        endif
        write(aline, *) str
        CALL say
        write(aline, *) '-------------------------------------------'
        CALL say
#ifdef MPP1
        CALL mpi_abort(MPI_COMM_WORLD,1,ierr)
#else
        CALL abort
#endif
      end subroutine fatal

      subroutine vprntv( a, n, nam )
        implicit none
#       include "parameter.inc"

        integer, intent(in) :: n
        character(len=8), intent(in) :: nam
        real(kind=8), dimension(n), intent(in) :: a
#ifdef MPP1
        if ( myid .eq. 0 ) print 1, nam, a
#else
        print 1, nam, a
#endif
1       format('0',a8,1x,1p,11g11.3,1x,/,9x,1p,11g11.3)
      return
      end subroutine vprntv

      subroutine vprntm( a, n1, n2, nam )
        implicit none
#       include "parameter.inc"

        integer, intent(in) :: n1, n2
        character(len=8), intent(in) :: nam
        real(kind=8), dimension(n1,n2), intent(in) :: a
        integer :: k, l
#ifdef MPP1
        if ( myid .eq. 0 ) then
          print 10, nam
          do k = 1, n1
            print 11, k, (a(k,l),l=1,n2)
          end do
        end if
#else
        print 10, nam
        do k = 1, n1
          print 11, k, (a(k,l),l=1,n2)
        end do
#endif
   10   format(1h1,a8,/)
   11   format(1x,i3,5x,1p,11g11.3,1x,/,9x,1p,11g11.3)
      return
      end subroutine vprntm
