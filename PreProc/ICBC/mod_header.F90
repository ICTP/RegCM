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

      module mod_header
      use m_stdio

      contains

      subroutine header(myname)
      implicit none 
      character (len=*) , intent(in) :: myname
      integer :: ihost , idir
      integer :: hostnm
      integer :: getcwd
      character (len=24) :: cdata='?'
      character (len=32) :: hostname='?' 
      character (len=32) :: user='?' 
      character (len=256) :: directory='?'
      integer , parameter :: nrite=stdout
  

      write (nrite,99002)  myname 
      write (nrite,99001)  SVN_REV, __DATE__ , __TIME__   

#ifdef IBM
        hostname='ibm platform '
        user= 'Unknown'
        call fdate_(cdata)
#else
        ihost = hostnm(hostname)
        call getlog(user)
        call fdate(cdata)
#endif 

        idir = getcwd(directory)

        write (nrite,*) ": this run start at    : ",cdata
        write (nrite,*) ": it is submitted by   : ",trim(user)
        write (nrite,*) ": it is running on     : ",trim(hostname)
        write (nrite,*) ": in directory         : ",trim(directory)
        write (nrite,*) "                      " 
      return 
99002 format(/,1x,' This is ',A,' part of the RegCM version 4')
99001 format(2x,' SVN Revision: ',A,' compiled at: data : ',A,          &
        &    '  time: ',A,/)
      end subroutine header

      subroutine finaltime(myid)
        implicit none
        integer , intent (in) :: myid
        character (len=24) :: cdata='?'

        if ( myid.eq. 0 ) then
#ifdef IBM
          call fdate_(cdata)
#else
          call fdate(cdata)
#endif 
          write ( stdout,* ) 'Input ready at : ', cdata
        end if

        return
      end subroutine finaltime

      end module mod_header
