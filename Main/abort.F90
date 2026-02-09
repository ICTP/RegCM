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

  subroutine myabort
    use mod_dynparam, only : mycomm
    use mod_intkinds
    use mpi
    implicit none (type, external)
    integer(ik4) :: ierr, myid
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    call date_and_time(date,time,zone)
    call mpi_comm_rank(mycomm, myid, ierr)
    write (0,*) 'Abort called by computing node ', myid, 'at ', &
            date(1:4),'-',date(5:6),'-',date(7:8),' ', &
            time(1:2),':',time(3:4),':',time(5:10),' ',&
            zone
    write(0,*) 'Execution terminated because of runtime error'
    call mpi_abort(mycomm,1,ierr)
  end subroutine myabort

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
