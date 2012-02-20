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

  use mod_constants
  use mod_realkinds
  use mod_date

  private

  public :: whoami , header , finaltime

  integer , parameter :: nrite=6
  character (len=24) :: cdata='?'
  integer , dimension(8) :: timearr
  real(dp) :: start_time

  contains

  subroutine whoami(myid)
    implicit none 
    integer , intent(in) :: myid

    if ( myid == 0 )  then 
      call date_and_time(values=timearr)
      start_time = julianday(timearr(1), timearr(2), timearr(3)) + &
                   timearr(5)*secph+ timearr(6)*secpm + &
                   dble(timearr(7)) + d_r1000*timearr(8)
      write (nrite,"(/,2x,'This is RegCM trunk')")
      write (nrite,99001)  SVN_REV, __DATE__ , __TIME__   
    end if

99001   format(2x,' SVN Revision: ',a,' compiled at: data : ',a,  &
    &    '  time: ',a,/)
  end subroutine whoami

  subroutine header(myid,nproc)
    implicit none 
    integer , intent(in) :: myid , nproc
    integer :: ihost , idir
    integer :: hostnm
    integer :: getcwd
    character (len=32) :: hostname='?' 
    character (len=32) :: user='?' 
    character (len=256) :: directory='?'
  
    if ( myid == 0 )  then 
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
      write (nrite,*) ": this run start at  : ",cdata
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
    integer , intent (in) :: myid
    real(dp) :: finish_time

    if ( myid ==  0 ) then
#ifdef IBM
      call fdate_(cdata)
#else
      call fdate(cdata)
#endif 
      call date_and_time(values=timearr)
      finish_time = julianday(timearr(1), timearr(2), timearr(3)) + &
                    timearr(5)*secph+ timearr(6)*secpm + &
                    dble(timearr(7)) + d_r1000*timearr(8)
      write (nrite,*) ': this run stops at  : ', cdata
      write (nrite,*) ': Total elapsed seconds of run : ', &
                      finish_time - start_time
    end if
  end subroutine finaltime

end module mod_header
