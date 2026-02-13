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
  use mod_realkinds
  use mod_constants
  use mod_mppparam
  use mod_date
  use mod_dynparam, only : nproc
  use mod_posix, only : hostname
  use mod_stdio

  implicit none

  private

  public :: whoami, header, checktime, finaltime

  character (len=32) :: cdata
  character (len=5) :: czone
  integer(ik4), dimension(8) :: tval
  real(rk4) :: start_time, last_time

  contains

  subroutine whoami(myid)
    implicit none
    integer(ik4), intent(in) :: myid
    character(len=*), parameter :: f99001 = &
        '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

    if ( myid == iocpu ) then
      call cpu_time(start_time)
      last_time = start_time
      write (stdout,"(/,2x,'This is RegCM development version')")
      write (stdout,f99001)  GIT_VER, __DATE__, __TIME__
    end if

  end subroutine whoami

  subroutine header(myid,nproc)
    implicit none
    integer(ik4), intent(in) :: myid, nproc
    character (len=32) :: hostnm
    character (len=32) :: user
    character (len=256) :: directory

    cdata = '?'
    czone = '?'
    hostnm = '?'
    user = '?'
    directory = '?'

    if ( myid == iocpu )  then
#ifdef IBM
      hostnm='ibm platform '
#else
      call hostname(hostnm)
#endif
      call date_and_time(zone=czone,values=tval)
      call get_environment_variable('PWD',directory)
      call get_environment_variable('USER',user)

      write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
            tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
      write (stdout,*) ": this run start at  : ",trim(cdata)
      write (stdout,*) ": it is submitted by : ",trim(user)
      write (stdout,*) ": it is running on   : ",trim(hostnm)
      write (stdout,*) ": it is using        : ",nproc, &
                       '  processors'
      write (stdout,*) ": in directory       : ",trim(directory)
      write (stdout,*) "                      "
    end if
  end subroutine header

  subroutine checktime(myid,ctime,period)
    implicit none
    integer(ik4), intent(in) :: myid
    character(len=*), intent(in) :: ctime, period
    integer(ik4) :: iunit
    real(rk4) :: check_time
    if ( myid == iocpu ) then
      call cpu_time(check_time)
      write (stdout,*) 'Elapsed seconds of run for this ',trim(period), &
            ' : ', (check_time-last_time)
      open(newunit=iunit,file=trim(ctime)//'.txt',form='formatted', &
           status='replace',action='write')
      write(iunit,*) 'Elapsed seconds of run for this ',trim(period), &
            ' : ', (check_time-last_time)
      close(iunit)
      last_time = check_time
    end if
  end subroutine checktime

  subroutine finaltime(myid)
    implicit none
    integer(ik4), intent (in) :: myid
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
