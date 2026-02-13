!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_header

  use mod_intkinds
  use mod_stdio
  use mod_posix, only : hostname

  implicit none

  public

  contains

  subroutine header(myname)
    implicit none
    character (len=*), intent(in) :: myname
    character (len=32) :: cdata
    character (len=5) :: czone
    character (len=32) :: hostnm
    character (len=32) :: user
    character (len=256) :: directory
    integer(ik4), dimension(8) :: tval
    character(len=*), parameter :: f99001 = &
           '(/,1x," This is ",A," part of the RegCM version 5")'
    character(len=*), parameter :: f99002 = &
           '(2x," SVN Revision: ",A," compiled at: data : ",A,"  time: ",A,/)'

    cdata = '?'
    czone = '?'
    hostnm = '?'
    user = '?'
    directory = '?'

    write (stdout,f99001)  myname
#ifdef NAGFOR
    write (stdout,f99001)  GIT_VER, '1900-01-01', '00:00:00'
#else
    write (stdout,f99002)  GIT_VER, __DATE__, __TIME__
#endif

#ifdef IBM
    hostnm = 'ibm platform '
#else
    call hostname(hostnm)
#endif
    call date_and_time(zone=czone,values=tval)
    call get_environment_variable('PWD',directory)
    call get_environment_variable('USER',user)

    write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
          tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
    write (stdout,*) ": this run start at    : ",trim(cdata)
    write (stdout,*) ": it is submitted by   : ",trim(user)
    write (stdout,*) ": it is running on     : ",trim(hostnm)
    write (stdout,*) ": in directory         : ",trim(directory)
    write (stdout,*) "                      "
  end subroutine header

  subroutine finaltime(myid)
    implicit none
    integer(ik4), intent (in) :: myid
    character (len=24) :: cdata
#ifdef __INTEL_COMPILER
    external :: fdate
#endif

    cdata = '?'
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
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
