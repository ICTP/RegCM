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
  use mod_dynparam
  use mod_regcm_interface
  use mod_runparams

  implicit none
  include 'mpif.h'

  real(rk8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer(ik4) :: ierr
!
!**********************************************************************
!
! Model Initialization
!
!**********************************************************************
!
  call mpi_init(ierr)
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
