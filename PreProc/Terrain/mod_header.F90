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

  subroutine header(myid)
  implicit none
!
  integer(ik4) , intent(in) :: myid
!
  integer(ik4) :: hostnm
  integer(ik4) :: ihost, idir
  integer(ik4) :: getcwd
  integer(ik4) , dimension(8) :: tval
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  character (len=32) :: hostname='?'
  character (len=32) :: user='?'
  character (len=128) :: directory='?'
  character (len=*) , parameter :: f99001 = &
    '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  if (myid.eq.1)  then
    write (stdout, "(/,2x,'This is Terrain part of RegCM package version 4 ')")
    write (stdout,f99001)  GIT_VER, __DATE__ , __TIME__

#ifdef IBM
    hostname='ibm platform '
    user= 'Unknown'
#else
    ihost = hostnm(hostname)
    call getlog(user)
#endif
    call date_and_time(zone=czone,values=tval)
    idir = getcwd(directory)

    write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
       tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
    write(stdout,*) ": this run start at  : ",trim(cdata)
    write(stdout,*) ": it is submitted by : ",trim(user)
    write(stdout,*) ": it is running on   : ",trim(hostname)
    write(stdout,*) ": in directory       : ",trim(directory)
    write(stdout,*) "                     "
  end if
  end subroutine header

end module mod_header
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
