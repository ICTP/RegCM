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

  private

  public :: whoami , header , finaltime

  integer(ik4) , parameter :: nrite=6
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  integer(ik4) :: clock_count , clock_rate , clock_max
  integer(ik4) , dimension(8) :: tval
  real(rk8) :: start_time

  contains

  subroutine whoami(myid)
    implicit none 
    integer(ik4) , intent(in) :: myid

    if ( myid == iocpu ) then 
      call system_clock(clock_count,clock_rate,clock_max)
      start_time = dble(clock_count)
      write (nrite,"(/,2x,'This is RegCM trunk')")
      write (nrite,99001)  SVN_REV, __DATE__ , __TIME__   
    end if

99001   format(2x,' SVN Revision: ',a,' compiled at: data : ',a,  &
    &    '  time: ',a,/)
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
      write (nrite,*) ": this run start at  : ",trim(cdata)
      write (nrite,*) ": it is submitted by : ",trim(user)
      write (nrite,*) ": it is running on   : ",trim(hostname)
      write (nrite,*) ": it is using        : ",nproc, &
                       '  processors' 
      write (nrite,*) ": in directory       : ",trim(directory)
      write (nrite,*) "                      " 
    end if 
  end subroutine header

  subroutine finaltime(myid)
    implicit none
    integer(ik4) , intent (in) :: myid
    real(rk8) :: finish_time

    if ( myid == iocpu ) then
      call system_clock(clock_count,clock_rate,clock_max)
      finish_time = dble(clock_count)
      call date_and_time(zone=czone,values=tval)
      write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
            tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
      write (nrite,*) ': this run stops at  : ', trim(cdata)
      write (nrite,*) ': Total elapsed seconds of run : ', &
                (finish_time - start_time)/dble(clock_rate)
    end if
  end subroutine finaltime

end module mod_header
