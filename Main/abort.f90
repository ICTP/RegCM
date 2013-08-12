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

  subroutine myabort
    use mod_dynparam , only : mycomm
    use mod_intkinds
    implicit none
    include 'mpif.h'
    integer(ik4) :: ierr , myid
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    call date_and_time(date,time,zone)
    call mpi_comm_rank(mycomm, myid, ierr)
    write (6,*) 'Abort called by computing node ', myid, 'at ', &
            date(1:4),'-',date(5:6),'-',date(7:8),' ', &
            time(1:2),':',time(3:4),':',time(5:10),' ',&
            zone
    call mpi_abort(mycomm,1,ierr)
  end subroutine myabort

