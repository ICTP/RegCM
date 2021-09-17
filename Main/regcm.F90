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

program regcm
!
!**********************************************************************
!
!     Used module declarations
!
!**********************************************************************
!
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_constants
  use mod_dynparam
  use mod_regcm_interface
  use mod_runparams
#ifndef MPI_SERIAL
  use mpi
#endif

  implicit none

  real(rk8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer(ik4) :: ierr , iprov
#ifdef MPI_SERIAL
  include 'mpif.h'
  integer(ik4) , parameter :: mpi_thread_single = 0
#endif
!
!**********************************************************************
!
! Model Initialization
!
!**********************************************************************
!
  call mpi_init_thread(mpi_thread_single,iprov,ierr)
  if ( ierr /= mpi_success ) then
    write(stderr,*) 'Cannot initilize MPI'
    stop
  end if
  call RCM_initialize()
!
!**********************************************************************
!
! Model Run
!
!**********************************************************************
!
  timestr = d_zero
  tdif = idate2 - idate1
  timeend = tohours(tdif) * secph

  call RCM_run(timestr, timeend)
!
!**********************************************************************
!
! Model Finalize
!
!**********************************************************************
!
  call RCM_finalize()
  call mpi_finalize(ierr)
!
end program regcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
