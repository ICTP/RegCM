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
  use mod_realkinds
  use mod_constants
  use mod_mppparam
  use mod_date
  use mod_dynparam , only : nproc
  use mod_stdio

  implicit none

  private

  public :: whoami , header , checktime , finaltime

  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  integer(ik4) , dimension(8) :: tval
  real(rk4) :: start_time , last_time

  contains

  subroutine whoami(myid)
    implicit none
    integer(ik4) , intent(in) :: myid
    character(len=*) , parameter :: f99001 = &
        '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

    if ( myid == iocpu ) then
      call cpu_time(start_time)
      last_time = start_time
      write (stdout,"(/,2x,'This is RegCM trunk')")
      write (stdout,f99001)  GIT_VER, __DATE__ , __TIME__
    end if

  end subroutine whoami

  subroutine header(myid,nproc)
    implicit none
    integer(ik4) , intent(in) :: myid , nproc
    integer(ik4) :: ihost , idir
    integer(ik4) :: hostnm
    integer(ik4) :: getcwd
    character (len=32) :: hostname='?'
    character (len=32) :: user='?'
    character (len=256) :: directory='?'

    if ( myid == iocpu )  then
#ifdef IBM
      hostname='ibm platform '
      user= 'Unknown'
#else
      ihost = hostnm(hostname)
      call getlog(user)
#endif
      idir = getcwd(directory)
      call date_and_time(zone=czone,values=tval)
      write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
            tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
      write (stdout,*) ": this run start at  : ",trim(cdata)
      write (stdout,*) ": it is submitted by : ",trim(user)
      write (stdout,*) ": it is running on   : ",trim(hostname)
      write (stdout,*) ": it is using        : ",nproc, &
                       '  processors'
      write (stdout,*) ": in directory       : ",trim(directory)
      write (stdout,*) "                      "
    end if
  end subroutine header

  subroutine checktime(myid,ctime)
    implicit none
    integer(ik4) , intent(in) :: myid
    character(len=*) , intent(in) :: ctime
    integer(ik4) :: iunit
    real(rk4) :: check_time
    if ( myid == iocpu ) then
      call cpu_time(check_time)
      write (stdout,*) 'Elapsed seconds of run for this month : ', &
                (check_time-last_time)
      open(newunit=iunit,file=trim(ctime)//'.txt',form='formatted', &
           status='replace',action='write')
      write(iunit,*) 'Elapsed seconds of run for this month : ', &
            (check_time-last_time)
      close(iunit)
      last_time = check_time
    end if
  end subroutine checktime

  subroutine finaltime(myid)
    implicit none
    integer(ik4) , intent (in) :: myid
    real(rkx) :: finish_time

    if ( myid == iocpu ) then
      call cpu_time(finish_time)
      call date_and_time(zone=czone,values=tval)
      write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
            tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
      write (stdout,*) ': this run stops at  : ', trim(cdata)
      write (stdout,*) ': Run has been completed using ', nproc, ' processors.'
      write (stdout,*) ': Total elapsed seconds of run : ', &
                (finish_time - start_time)
    end if
  end subroutine finaltime

end module mod_header
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
