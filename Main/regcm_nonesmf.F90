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
  use mod_memutil
  use mod_runparams
  use mod_date
  use mod_message
  use mod_ncio
  use mod_output
  use mod_split
  use mod_bdycod
  use mod_che_semdde
  use mod_init
  use mod_header
  use mod_param
  use mod_tendency
  use mod_tstep
  use mod_service
  use mod_interface
#ifdef CHEMTEST
  use mod_chem
#endif
  use mod_mppio
  use mpi
#ifdef CLM
  use perf_mod
  use spmdMod, only: mpicom
#endif
!
  real(8) :: timestr, timeend
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
  first = .TRUE.
  timestr = 0.0D0
  timeend = idatediff(idate2,idate1)*3600.0D0

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
