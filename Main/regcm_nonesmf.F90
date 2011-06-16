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
  use mod_interface
  use mod_runparams
  use mod_interface
  use mpi
!
  implicit none
!
  real(8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer :: ierr
  logical :: first
!
!**********************************************************************
!
! Model Initialization
!
!**********************************************************************
!
  call mpi_init(ierr)
!
  call RCM_initialize()
!
!**********************************************************************
!
! Model Run 
!
!**********************************************************************
!
  first = .true.
  timestr = d_zero
  tdif = idate2 - idate1
  timeend = tdif%hours() * secph

  call RCM_run(timestr, timeend, first)
!
!**********************************************************************
!
! Model Finalize 
!
!**********************************************************************
!
  call RCM_finalize()
!
end program regcm
