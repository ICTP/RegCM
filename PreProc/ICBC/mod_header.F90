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

  use mod_intkinds
  use mod_stdio

  contains

  subroutine header(myname)
    implicit none 
    character (len=*) , intent(in) :: myname
    integer(ik4) :: ihost , idir
    integer(ik4) :: hostnm
    integer(ik4) :: getcwd
    character (len=32) :: cdata = '?'
    character (len=5) :: czone = '?'
    character (len=32) :: hostname = '?' 
    character (len=32) :: user = '?' 
    character (len=256) :: directory = '?'
    integer(ik4) , parameter :: nrite = stdout
    integer(ik4) , dimension(8) :: tval

    write (nrite,99002)  myname 
    write (nrite,99001)  SVN_REV, __DATE__ , __TIME__   

#ifdef IBM
    hostname = 'ibm platform '
    user = 'Unknown'
#else
    ihost = hostnm(hostname)
    call getlog(user)
#endif 
    call date_and_time(zone=czone,values=tval)
    idir = getcwd(directory)

    write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
          tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
    write (nrite,*) ": this run start at    : ",trim(cdata)
    write (nrite,*) ": it is submitted by   : ",trim(user)
    write (nrite,*) ": it is running on     : ",trim(hostname)
    write (nrite,*) ": in directory         : ",trim(directory)
    write (nrite,*) "                      " 
99002 format(/,1x,' This is ',A,' part of the RegCM version 4')
99001 format(2x,' SVN Revision: ',A,' compiled at: data : ',A,'  time: ',A,/)
  end subroutine header

  subroutine finaltime(myid)
    implicit none
    integer(ik4) , intent (in) :: myid
    character (len=24) :: cdata = '?'

    if ( myid ==  0 ) then
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
