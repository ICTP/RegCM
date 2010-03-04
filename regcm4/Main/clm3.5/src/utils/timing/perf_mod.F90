module perf_mod

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for controlling the performance
!          timer logic.
! 
! Author:  P. Worley, January 2007
!
! $Id$
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_sys_mod,       only: shr_sys_abort
   use shr_kind_mod,      only: shr_kind_cl, shr_kind_r8

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public t_initf
   public t_profile_onf
   public t_barrier_onf
   public t_single_filef
   public t_stampf
   public t_startf
   public t_stopf
   public t_enablef
   public t_disablef
   public t_adj_detailf
   public t_barrierf
   public t_prf
   public t_finalizef

!-----------------------------------------------------------------------
! Private interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   private perf_defaultopts
   private perf_setopts

!-----------------------------------------------------------------------
!- include statements --------------------------------------------------
!-----------------------------------------------------------------------
#ifdef HAVE_PAPI
#include <f77papi.h>
#endif
#include <mpif.h>  
#include "gptl.inc"

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
   logical, parameter :: def_perf_single_file = .true.          ! default
   logical, private   :: perf_single_file = def_perf_single_file
                         ! flag indicating whether the performance timer
                         ! output should be written to a single file 
                         ! (per component communicator) or to a 
                         ! separate file for each process

   logical, parameter :: def_timing_initialized = .false.       ! default
   logical, private   :: timing_initialized = def_timing_initialized
                         ! flag indicating whether timing library has
                         ! been initialized

   logical, parameter :: def_timing_disable = .false.           ! default
   logical, private   :: timing_disable = def_timing_disable
                         ! flag indicating whether timers are disabled

   logical, parameter :: def_timing_barrier = .false.           ! default
   logical, private   :: timing_barrier = def_timing_barrier
                         ! flag indicating whether the mpi_barrier in
                         ! t_barrierf should be called

   integer, parameter :: def_timer_depth_limit = 99999          ! default
   integer, private   :: timer_depth_limit = def_timer_depth_limit
                         ! integer indicating maximum number of levels of
                         ! timer nesting 

   integer, parameter :: def_timing_detail_limit = 0            ! default
   integer, private   :: timing_detail_limit = def_timer_depth_limit
                         ! integer indicating maximum detail level to
                         ! profile

   integer, parameter :: init_timing_disable_depth = 0          ! init
   integer, private   :: timing_disable_depth = init_timing_disable_depth
                         ! integer indicating depth of t_disablef calls

   integer, parameter :: init_timing_detail = 0                 ! init
   integer, private   :: cur_timing_detail = init_timing_detail
                         ! current timing detail level

#ifdef UNICOSMP
   integer, parameter :: def_perf_timer = GPTLrtc               ! default
#else
   integer, parameter :: def_perf_timer = GPTLmpiwtime          ! default
#endif
   integer, private   :: perf_timer = def_perf_timer            ! default
                         ! integer indicating which timer to use
                         ! (as defined in gptl.inc)

!=======================================================================
contains
!=======================================================================

!
!========================================================================
!
   subroutine perf_defaultopts(timing_disable_out, &
                               perf_timer_out, &
                               timer_depth_limit_out, &
                               timing_detail_limit_out, &
                               timing_barrier_out, &
                               perf_single_file_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
   ! timers disable/enable option
   logical, intent(out), optional :: timing_disable_out
   ! performance timer option
   integer, intent(out), optional :: perf_timer_out
   ! timer depth limit option
   integer, intent(out), optional :: timer_depth_limit_out
   ! timer detail limit option
   integer, intent(out), optional :: timing_detail_limit_out
   ! timing barrier enable/disable option
   logical, intent(out), optional :: timing_barrier_out
   ! timing single / multple output file options
   logical, intent(out), optional :: perf_single_file_out
!-----------------------------------------------------------------------
   if ( present(timing_disable_out) ) then
      timing_disable_out = def_timing_disable
   endif
   if ( present(perf_timer_out) ) then
      perf_timer_out = def_perf_timer
   endif
   if ( present(timer_depth_limit_out) ) then
      timer_depth_limit_out = def_timer_depth_limit
   endif
   if ( present(timing_detail_limit_out) ) then
      timing_detail_limit_out = def_timing_detail_limit
   endif
   if ( present(timing_barrier_out) ) then
      timing_barrier_out = def_timing_barrier
   endif
   if ( present(perf_single_file_out) ) then
      perf_single_file_out = def_perf_single_file
   endif
!
   return
   end subroutine perf_defaultopts
!
!========================================================================
!
   subroutine perf_setopts(mastertask, &
                           LogPrint, &
                           timing_disable_in, &
                           perf_timer_in, &
                           timer_depth_limit_in, &
                           timing_detail_limit_in, &
                           timing_barrier_in, &
                           perf_single_file_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments----------------------------
!
   ! master process?
   logical, intent(in) :: mastertask
   ! Print out to log file?
   logical, intent(IN) :: LogPrint        
   ! timers disable/enable option
   logical, intent(in), optional :: timing_disable_in
   ! performance timer option
   integer, intent(in), optional :: perf_timer_in
   ! timer depth limit option
   integer, intent(in), optional :: timer_depth_limit_in
   ! timer detail limit option
   integer, intent(in), optional :: timing_detail_limit_in
   ! timing barrier enable/disable option
   logical, intent(in), optional :: timing_barrier_in
   ! timing single / multple output file options
   logical, intent(in), optional :: perf_single_file_in
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! error return
!-----------------------------------------------------------------------
   if ( .not. timing_initialized ) then

      if ( present(timing_disable_in) ) then
         timing_disable = timing_disable_in
         if (timing_disable) then
            ierr = GPTLdisable()
         else 
            ierr = GPTLenable()
         endif
      endif
      if ( present(perf_timer_in) ) then
         if ((perf_timer_in .eq. GPTLgettimeofday) .or. &
             (perf_timer_in .eq. GPTLnanotime) .or. &
             (perf_timer_in .eq. GPTLrtc) .or. &
             (perf_timer_in .eq. GPTLmpiwtime) .or. &
             (perf_timer_in .eq. GPTLclockgettime) .or. &
             (perf_timer_in .eq. GPTLpapitime)) then
            perf_timer = perf_timer_in
!pw            call GPTLsetutr(perf_timer)
         else
            if (mastertask) then
               write(6,*) 'PERF_SETOPTS: illegal timer requested=',&
                          perf_timer_in, '. Request ignored.'
            endif
         endif
      endif
      if ( present(timer_depth_limit_in) ) then
         timer_depth_limit = timer_depth_limit_in
      endif
      if ( present(timing_detail_limit_in) ) then
         timing_detail_limit = timing_detail_limit_in
      endif
      if ( present(timing_barrier_in) ) then
         timing_barrier = timing_barrier_in
      endif
      if ( present(perf_single_file_in) ) then
         perf_single_file = perf_single_file_in
      endif
!
      if (mastertask .and. LogPrint) then
         write(6,*) 'PERF_SETOPTS:  Using profile_disable=', timing_disable, &             
               ' profile_timer=', perf_timer, &
               ' profile_depth_limit=', timer_depth_limit, &    
               ' profile_detail_limit=', timing_detail_limit, &
               ' profile_barrier=', timing_barrier, &
               ' profile_single_file=', perf_single_file 
      endif                                                                               
!
#ifdef DEBUG
   else
      write(6,*) 'PERF_SETOPTS: timing library already initialized. Request ignored.'
#endif
   endif
!
   return
   end subroutine perf_setopts

!
!========================================================================
!
   logical function t_profile_onf()
!----------------------------------------------------------------------- 
! Purpose: Return flag indicating whether profiling is currently active.
!          Part of workaround to implement FVbarrierclock before
!          communicators exposed in Pilgrim. Does not check level of
!          event nesting.
! Author: P. Worley 
!-----------------------------------------------------------------------

   if ((.not. timing_initialized) .or. &
       (timing_disable_depth > 0) .or. &
       (cur_timing_detail > timing_detail_limit)) then
      t_profile_onf = .false.
   else
      t_profile_onf = .true.
   endif

   end function t_profile_onf
!
!========================================================================
!
   logical function t_barrier_onf()
!----------------------------------------------------------------------- 
! Purpose: Return timing_barrier. Part of workaround to implement 
!          FVbarrierclock before communicators exposed in Pilgrim. 
! Author: P. Worley 
!-----------------------------------------------------------------------

   t_barrier_onf = timing_barrier

   end function t_barrier_onf
!
!========================================================================
!
   logical function t_single_filef()
!----------------------------------------------------------------------- 
! Purpose: Return perf_single_file. Used to control output of other
!          performance data, only spmdstats currently.
! Author: P. Worley 
!-----------------------------------------------------------------------

   t_single_filef = perf_single_file

   end function t_single_filef
!
!========================================================================
!
   subroutine t_stampf(wall, usr, sys)
!----------------------------------------------------------------------- 
! Purpose: Record wallclock, user, and system times (seconds).
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Output arguments-----------------------------
!
   real(shr_kind_r8), intent(out) :: wall ! wallclock time
   real(shr_kind_r8), intent(out) :: usr  ! user time
   real(shr_kind_r8), intent(out) :: sys  ! system time
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if ((.not. timing_initialized) .or. &
       (timing_disable_depth > 0)) then
      wall = 0.0
      usr = 0.0
      sys = 0.0
   else
      ierr = GPTLstamp(wall, usr, sys)
   endif

   return
   end subroutine t_stampf
!
!========================================================================
!
   subroutine t_startf(event)
!----------------------------------------------------------------------- 
! Purpose: Start an event timer
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   character(len=*), intent(in) :: event  ! performance timer event name
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      ierr = GPTLstart(event)

   endif

   return
   end subroutine t_startf
!
!========================================================================
!
   subroutine t_stopf(event)
!----------------------------------------------------------------------- 
! Purpose: Stop an event timer
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   character(len=*), intent(in) :: event  ! performance timer event name
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      ierr = GPTLstop(event)

   endif

   return
   end subroutine t_stopf
!
!========================================================================
!
   subroutine t_enablef()
!----------------------------------------------------------------------- 
! Purpose: Enable t_startf, t_stopf, t_stampf, and t_barrierf. Ignored
!          in threaded regions.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!---------------------------Externals-----------------------------------
!
#if ( defined _OPENMP )
   logical omp_in_parallel
   external omp_in_parallel
#endif
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

#if ( defined _OPENMP )
   if (omp_in_parallel()) return
#endif

   if (timing_disable_depth > 0) then
      if (timing_disable_depth .eq. 1) then
         ierr = GPTLenable()
      endif
      timing_disable_depth = timing_disable_depth - 1
   endif

   return
   end subroutine t_enablef
!
!========================================================================
!
   subroutine t_disablef()
!----------------------------------------------------------------------- 
! Purpose: Disable t_startf, t_stopf, t_stampf, and t_barrierf. Ignored
!          in threaded regions.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!---------------------------Externals-----------------------------------
!
#if ( defined _OPENMP )
   logical omp_in_parallel
   external omp_in_parallel
#endif
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

#if ( defined _OPENMP )
   if (omp_in_parallel()) return
#endif

   if (timing_disable_depth .eq. 0) then
      ierr = GPTLdisable()
   endif
   timing_disable_depth = timing_disable_depth + 1

   return
   end subroutine t_disablef
!
!========================================================================
!
   subroutine t_adj_detailf(detail_adjustment)
!----------------------------------------------------------------------- 
! Purpose: Modify current detail level. Ignored in threaded regions.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   integer, intent(in) :: detail_adjustment ! user defined increase or
                                            ! decrease in detail level
!
!---------------------------Externals-----------------------------------
!
#if ( defined _OPENMP )
   logical omp_in_parallel
   external omp_in_parallel
#endif
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

#if ( defined _OPENMP )
   if (omp_in_parallel()) return
#endif

   cur_timing_detail = cur_timing_detail + detail_adjustment

   return
   end subroutine t_adj_detailf
!
!========================================================================
!
   subroutine t_barrierf(event, mpicom)
!----------------------------------------------------------------------- 
! Purpose: Call (and time) mpi_barrier. Ignored inside OpenMP
!          threaded regions. Note that barrier executed even if
!          event not recorded because of level of timer event nesting.
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Use statements------------------------------
   use shr_mpi_mod,       only: shr_mpi_barrier
!---------------------------Input arguments-----------------------------
   ! mpi communicator id
   integer, intent(in), optional :: mpicom
   ! performance timer event name
   character(len=*), intent(in), optional :: event
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!---------------------------Externals-----------------------------------
!
#if ( defined _OPENMP )
   logical omp_in_parallel
   external omp_in_parallel
#endif
!
!-----------------------------------------------------------------------
!
#if ( defined _OPENMP )
   if (omp_in_parallel()) return
#endif
   if ((timing_initialized) .and. &
       (timing_disable_depth .eq. 0) .and. &
       (cur_timing_detail .le. timing_detail_limit)) then

      if (timing_barrier) then

         if ( present (event) ) then
            ierr = GPTLstart(event)
         endif

         if ( present (mpicom) ) then
            call shr_mpi_barrier(mpicom, 'T_BARRIERF: bad mpi communicator')
         else
            call shr_mpi_barrier(MPI_COMM_WORLD, 'T_BARRIERF: bad mpi communicator')
         endif

         if ( present (event) ) then
            ierr = GPTLstop(event)
         endif

      endif

   endif

   return
   end subroutine t_barrierf
!
!========================================================================
!
   subroutine t_prf(filename, mpicom)
!----------------------------------------------------------------------- 
! Purpose: Write out performance timer data
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Use statements------------------------------
   use shr_file_mod,      only: shr_file_getUnit, shr_file_freeUnit
!
!---------------------------Input arguments-----------------------------
!
   ! performance timer output file name
   character(len=*), intent(in), optional :: filename
   ! mpi communicator id
   integer, intent(in), optional :: mpicom
!
!---------------------------Local workspace-----------------------------
!
   integer  i                     ! loop index
   integer  me                    ! communicator local process id
   integer  npes                  ! local communicator group size
   integer  gme                   ! global process id
   integer  ierr                  ! MPI error return
   integer  signal                ! send/recv variable for single
                                  ! output file logic
   integer  str_length            ! string length
   integer  unitn                 ! file unit number
   integer status (MPI_STATUS_SIZE)    ! Status of message
   character(len=5) cme                ! string representation of process id
   character(len=SHR_KIND_CL+12) fname ! timing output filename
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

!$OMP MASTER
   call mpi_comm_rank(MPI_COMM_WORLD, gme, ierr)
   if ( present(mpicom) ) then
      call mpi_comm_size(mpicom, npes, ierr)
!pw      if (ierr .eq. MPI_ERR_COMM) then
!pw         call shr_sys_abort('T_PRF: bad mpi communicator')
!pw      endif
      call mpi_comm_rank(mpicom, me, ierr)
   else
      call mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
      me = gme
   endif

   do i=1,SHR_KIND_CL+12
     fname(i:i) = " "
   enddo

   unitn = shr_file_getUnit()

   ! If a single timing output file, take turns writing to it.
   if (perf_single_file) then

      if ( present(filename) ) then
         str_length = min(SHR_KIND_CL,len_trim(filename))
         fname(1:str_length) = filename(1:str_length)
      else
         fname(1:10) = "timing_all"
      endif

      signal = 0
      if (me .eq. 0) then

         open( unitn, file=trim(fname), status='UNKNOWN' )
         write( unitn, 101) me, gme
 101     format(/,"************ PROCESS ",I5," (",I5,") ************",/)
         close( unitn )

         ierr = GPTLpr(0, trim(fname))

      else

         call mpi_recv (signal, 1, mpi_integer, me-1, me-1, mpicom, status, ierr)
         if (ierr /= mpi_success) then
            write(6,*) 'T_PRF: mpi_recv failed ierr=',ierr
            call shr_sys_abort()
         end if

         open( unitn, file=trim(fname), status='OLD', position='APPEND' )
         write( unitn, 101) me, gme
         close( unitn )

         ierr = GPTLpr(0, trim(fname))

      endif

      if (me+1 < npes) &
         call mpi_send (signal, 1, mpi_integer, me+1, me, mpicom, ierr)

   else

      write(cme,'(i5.5)') me
      if ( present(filename) ) then
         str_length = min(SHR_KIND_CL-6,len_trim(filename))
         fname(1:str_length) = filename(1:str_length)
      else
         str_length = 6
         fname(1:10) = "timing"
      endif
      fname(str_length+1:str_length+1) = '.'
      fname(str_length+2:str_length+6) = cme

      open( unitn, file=trim(fname), status='UNKNOWN' )
      write( unitn, 101) me, gme
      close( unitn )

      ierr = GPTLpr(0, trim(fname))

   endif

   call shr_file_freeUnit( unitn )
!$OMP END MASTER

   return
   end subroutine t_prf
!
!========================================================================
!
   subroutine t_initf(NLFilename, LogPrint, mpicom, MasterTask)
!----------------------------------------------------------------------- 
! Purpose:  Set default values of runtime timing options 
!           before namelist prof_inparm is read,
!           read namelist (and broadcast, if SPMD),
!           then initialize timing library.
! Author:   P. Worley (based on shr_inputinfo_mod and runtime_opts)
!-----------------------------------------------------------------------
!---------------------------Use statements- ----------------------------
   use shr_file_mod,      only: shr_file_getUnit, shr_file_freeUnit
   use shr_mpi_mod,       only: shr_mpi_bcast
!-----------------------------------------------------------------------
!---------------------------Include statements--------------------------
!-----------------------------------------------------------------------
#ifdef HAVE_PAPI
#include <f77papi.h>
#endif
!
!---------------------------Input arguments-----------------------------
!
   character(len=*),   intent(IN) :: NLFilename      ! Name-list filename
   logical, optional,  intent(IN) :: LogPrint        ! If print out to log file
   integer, optional,  intent(IN) :: mpicom          ! MPI communicator
   logical, optional,  intent(IN) :: MasterTask      ! If MPI master task
!
!---------------------------Local workspace-----------------------------
!
   character(len=*), parameter    :: subname = '(T_INITF) '
   logical                        :: MasterTask2     ! If MPI master task
   logical                        :: LogPrint2       ! If print to stdout
!
!---------------------------Namelist -----------------------------------
!
   logical profile_disable
   logical profile_barrier
   logical profile_single_file
   integer profile_depth_limit
   integer profile_detail_limit
   integer profile_timer
   namelist /prof_inparm/ profile_disable, profile_barrier, &
                          profile_single_file, profile_depth_limit, &
                          profile_detail_limit, profile_timer
!
!---------------------------Local workspace-----------------------------
!
   integer  me                    ! communicator local process id
   integer  ierr                  ! error return
   integer  unitn                 ! file unit number
!-----------------------------------------------------------------------
    if ( timing_initialized ) then
#ifdef DEBUG
       write(6,*) 'T_INITF: timing library already initialized. Request ignored.'
#endif
       return
    endif

!$OMP MASTER
    if ( present(MasterTask) .and. present(mpicom) )then
       call mpi_comm_rank(mpicom, me, ierr)
!pw      if (ierr .eq. MPI_ERR_COMM) then
!pw         call shr_sys_abort('T_INITF: bad mpi communicator')
!pw      endif
       if (me .eq. 0) then
          MasterTask2 = .true.
       else
          MasterTask2 = .false.
       endif
    else
       MasterTask2 = .true.
    end if

    if ( present(LogPrint) ) then
       LogPrint2 = LogPrint
    else
       LogPrint2 = .true.
    endif

    ! Set defaults, then override with user-specified input
    call perf_defaultopts(timing_disable_out=profile_disable, &
                          perf_timer_out=profile_timer, &
                          timer_depth_limit_out=profile_depth_limit, &
                          timing_detail_limit_out=profile_detail_limit, &
                          timing_barrier_out=profile_barrier, &
                          perf_single_file_out=profile_single_file )
    if ( MasterTask2 ) then
       unitn = shr_file_getUnit()
       write(6,*) 'Read in prof_inparm namelist from: '//trim(NLFilename)
       ierr = 1
       open( unitn, file=trim(NLFilename), status='old', iostat=ierr )
       if (ierr .eq. 0) then
          ierr = 1
          do while( ierr > 0 )
             read(unitn, nml=prof_inparm, iostat=ierr)
          end do
          close(unitn)
       endif
       call shr_file_freeUnit( unitn )
    endif
    ! This logic assumes that there will be only one MasterTask
    ! per communicator, and that this MasterTask is process 0.
    if ( present(MasterTask) .and. present(mpicom) )then
       call shr_mpi_bcast( profile_disable,    MPICom )
       call shr_mpi_bcast( profile_barrier,    MPICom )
       call shr_mpi_bcast( profile_single_file,  MPICom )
       call shr_mpi_bcast( profile_depth_limit, MPICom )
       call shr_mpi_bcast( profile_detail_limit, MPICom )
       call shr_mpi_bcast( profile_timer,        MPICom )
    end if
    call perf_setopts    (MasterTask2, LogPrint2, &
                          timing_disable_in=profile_disable, &
                          perf_timer_in=profile_timer, &
                          timer_depth_limit_in=profile_depth_limit, &
                          timing_detail_limit_in=profile_detail_limit, &
                          timing_barrier_in=profile_barrier, &
                          perf_single_file_in=profile_single_file )
!$OMP END MASTER
!$OMP BARRIER

   if (timing_disable) return

!$OMP MASTER
   !
   ! Set options and initialize timing library.  
   ! 
   ! For logical settings, 2nd arg 0 
   ! to gptlsetoption means disable, non-zero means enable
   !
   ! Turn off CPU timing (expensive)
   !
   if (gptlsetoption (gptlcpu, 0) < 0) call shr_sys_abort (subname//':: gptlsetoption')
   !
   ! Set max timer depth
   !
   if (gptlsetoption (gptldepthlimit, timer_depth_limit) < 0) &
     call shr_sys_abort (subname//':: gptlsetoption')
   !
   ! Next 2 calls only work if PAPI is enabled.  These examples enable counting
   ! of total cycles and floating point ops, respectively
   !
#ifdef HAVE_PAPI
!pw should make these runtime options in the future
   if (gptlsetoption (PAPI_TOT_CYC, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
   if (gptlsetoption (PAPI_FP_OPS, 1) < 0) call shr_sys_abort (subname//':: gptlsetoption')
#endif
   !
   ! Initialize the timing lib.  This call must occur after all gptlsetoption
   ! calls and before all other timing lib calls.
   !
   if (gptlinitialize () < 0) call shr_sys_abort (subname//':: gptlinitialize')
   timing_initialized = .true.
!$OMP END MASTER
!$OMP BARRIER

   return
   end subroutine t_initf
!
!========================================================================
!
   subroutine t_finalizef()
!----------------------------------------------------------------------- 
! Purpose: shut down timing library
! Author: P. Worley 
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer  ierr                  ! GPTL error return
!
!-----------------------------------------------------------------------
!
   if (.not. timing_initialized) return

!$OMP MASTER
   ierr = GPTLfinalize()
   timing_initialized = .false.
!$OMP END MASTER
!$OMP BARRIER

   return
   end subroutine t_finalizef

end module perf_mod
