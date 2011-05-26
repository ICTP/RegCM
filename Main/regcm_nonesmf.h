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
!     Model Initialization
!
!**********************************************************************
!
      call mpi_init(ierr)
!
      call RCM_initialize()
!
!**********************************************************************
!
!     Model Run 
!
!**********************************************************************
!
      first = .TRUE.
      timestr = 0.0
      timeend = idatediff(idate2,idate1)*3600.0
      call RCM_run(timestr, timeend, first)
!
!**********************************************************************
!
!     Model Finalize 
!
!**********************************************************************
!
      call RCM_finalize()
!
      end program regcm
