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
#ifdef OASIS
  use mod_oasis_interface
#endif

  implicit none (type, external)

  real(rk8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer(ik4) :: ierr, iprov
#ifdef OASIS
  integer :: localCommunicator
#endif
#ifdef MPI_SERIAL
  include 'mpif.h'
  integer(ik4), parameter :: mpi_thread_single = 0
#endif
!
!**********************************************************************
!
! Model Initialization
!
!**********************************************************************
!
#ifndef OASIS
  call mpi_init_thread(mpi_thread_funneled,iprov,ierr)
  if ( ierr /= mpi_success ) then
    write(stderr,*) 'Cannot initilize MPI'
    stop
  end if
  call RCM_initialize()
#else
  !
  ! OASIS Initialization
  !
  call oasisxregcm_init(localCommunicator)
  call RCM_initialize(localCommunicator)
#endif
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
#ifndef OASIS
  call mpi_finalize(ierr)
#else
  !
  ! OASIS Finalization
  !
  call oasisxregcm_finalize
#endif
!
end program regcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
