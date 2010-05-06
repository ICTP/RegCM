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

#ifdef MPP1
      use mod_dynparam , only : myid
#endif
      use mod_dynparam , only : debug_level

      implicit none
!
      character(512) :: aline

      contains

      function active(level)
      implicit none
!
! Dummy arguments
!
      integer :: level
      logical :: active
      intent (in) level
!
      active = (level.le.debug_level)
      end function active
 
      subroutine say
      implicit none
#ifdef MPP1
      if ( myid.eq.0 ) print * , trim(aline)
#else
      print * , trim(aline)
#endif
      end subroutine say
 
      subroutine warning
      implicit none
#ifdef MPP1
      print * , ' Processor ' , myid , ' : ' , trim(aline)
#else
      print * , trim(aline)
#endif
      end subroutine warning
 
      subroutine fatal(filename,line,str)
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
      character(*) :: filename , str
      integer :: line
      intent (in) filename , line , str
!
! Local variables
!
      character(8) :: cline
#ifdef MPP1
      integer :: ierr
#endif
!
      write (cline,'(i6)') line
      write (aline,*) '-------------- FATAL CALLED ---------------'
      call say
      if ( line.gt.0 ) then
        write (aline,*) 'Fatal in file: '//filename//' at line: '//     &
                      & trim(cline)
        call say
      end if
      write (aline,*) str
      call say
      write (aline,*) '-------------------------------------------'
      call say
#ifdef MPP1
      call mpi_abort(mpi_comm_world,1,ierr)
#else
      call abort
#endif
      end subroutine fatal
 
      subroutine vprntv(a,n,nam)
      implicit none
!
! Dummy arguments
!
      integer :: n
      character(8) :: nam
      real(8) , dimension(n) :: a
      intent (in) a , n , nam
 
#ifdef MPP1
      if ( myid.eq.0 ) print 99001 , nam , a
#else
      print 99001 , nam , a
#endif
99001 format ('0',a8,1x,1p,11G11.3,1x,/,9x,1p,11G11.3)
      end subroutine vprntv
 
      subroutine vprntm(a,n1,n2,nam)
      implicit none
!
! Dummy arguments
!
      integer :: n1 , n2
      character(8) :: nam
      real(8) , dimension(n1,n2) :: a
      intent (in) a , n1 , n2 , nam
!
! Local variables
!
      integer :: k , l
 
#ifdef MPP1
      if ( myid.eq.0 ) then
        print 99001 , nam
        do k = 1 , n1
          print 99002 , k , (a(k,l),l=1,n2)
        end do
      end if
#else
      print 99001 , nam
      do k = 1 , n1
        print 99002 , k , (a(k,l),l=1,n2)
      end do
#endif
99001 format ('1',a8,/)
99002 format (1x,i3,5x,1p,11G11.3,1x,/,9x,1p,11G11.3)
      end subroutine vprntm

      end module mod_message
