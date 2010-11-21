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

      subroutine header(myid)
      implicit none 
!
! Dummy arguments
!
      integer , intent(in) :: myid
!
! Local variables:
!
      integer :: hostnm
      integer :: ihost, idir
      integer :: getcwd
      character (len=32) :: cdata='?'
      character (len=32) :: hostname='?' 
      character (len=32) :: user='?' 
      character (len=128) :: directory='?'
!
      if (myid.eq.1)  then 
        write (stdout,"(/,2x, &
          'This is Terrain part of RegCM package version 4 ')")
        write (stdout,100)  SVN_REV, __DATE__ , __TIME__   
100     format(2x,' SVN Revision: ',a,' compiled at: data : ',          &
           &   a,'  time: ',a,/)

#ifdef IBM
        hostname='ibm platform '
        user= 'Unknown'
        call fdate_(cdata)
#else
        Ihost = hostnm(hostname)
        call getlog(user)
        call fdate(cdata)
#endif 

        Idir=getcwd(directory)

        write(stdout,*) ": this run start at  : ",trim(cdata)
        write(stdout,*) ": it is submitted by : ",trim(user)
        write(stdout,*) ": it is running on   : ",trim(hostname)
        write(stdout,*) ": in directory       : ",trim(directory)
        write(stdout,*) "                     " 
      end if 
      end subroutine header

      end module mod_header
