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
!

module mod_interface

!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
  use mod_memutil
  use mod_runparams
  use mod_message
  use mod_ncio
  use mod_output
  use mod_split
  use mod_bdycod
  use mod_che_semdde
  use mod_init
  use mod_header
  use mod_params
  use mod_tendency
  use mod_tstep
  use mod_service
#ifdef CHEMTEST
  use mod_chem
#endif
  use mod_mppio
  use mpi
#ifdef CLM
  use perf_mod
  use spmdMod, only: mpicom
#endif
  implicit none
!
  private
  public :: RCM_initialize
  public :: RCM_run
  public :: RCM_finalize

  real(8) :: dtinc

  contains
 
  subroutine RCM_initialize(mpiCommunicator)
!
!**********************************************************************
!
!     Imported variable declarations 
!
!**********************************************************************
!
    integer, intent(in), optional :: mpiCommunicator
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
    integer :: MyComm, ncpu, ierr
    character(256) :: namelistfile, prgname
! 
!**********************************************************************
!
!**********************************************************************
!
!     MPI Initialization
!
!**********************************************************************
!
    if (present(mpiCommunicator)) then
      MyComm = mpiCommunicator
    else
      MyComm = MPI_COMM_WORLD
    end if
    call mpi_comm_rank(MyComm, myid, ierr)
    call mpi_comm_size(MyComm, ncpu, ierr)
!
    call whoami(myid)
!
#ifdef DEBUG 
    call activate_debug()
#endif
!
!**********************************************************************
!
!     Read input global namelist
!
!**********************************************************************
!
    if ( myid == 0 ) then
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr/=0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
    end if
!
    call memory_init
!
    if ( myid == 0 ) then
#ifdef BAND
      call init_mod_ncio(.true.)
#else
      call init_mod_ncio(.false.)
#endif
    end if

!
!**********************************************************************
!
!     MPI Initialization for model
!
!**********************************************************************
!
    call broadcast_params

    call set_nproc(ncpu)

    if ( ncpu /= nproc ) then
      write (aline,*) 'The number of CPU is not well set'
      call say
      write (aline,*) 'NCPU = ' , ncpu , '    nproc =' , nproc
      call say
      call fatal(__FILE__,__LINE__,'CPU Count mismatch')
    end if
!      print * , "process" , myid , "of" , nproc
    call mpi_barrier(mpi_comm_world,ierr)
!     starttime= MPI_WTIME()
    if ( myid > 0 ) then
      iwest = myid - 1
    else
#ifdef BAND
      iwest = nproc-1
#else
      iwest = mpi_proc_null
#endif
    end if
    if ( myid < nproc-1 ) then
      ieast = myid + 1
    else
#ifdef BAND
      ieast = 0
#else
      ieast = mpi_proc_null
#endif
    end if
    if ( jxp < 3 ) then
      write (aline,*) 'The number of jxp must be greater than 2'
      call say
      write (aline,*) 'jxp = ' , jxp , '   jx = ' , jx
      call say
      call fatal(__FILE__,__LINE__,'Domain too small')
    end if
    if ( jxp*nproc /= jx ) then
      write (aline,*) 'jx should be divided by nproc'
      call say
      write (aline,*) 'jx = ' , jx , '   nproc = ' , nproc
      call say
      call fatal(__FILE__,__LINE__,                                   &
               & 'Domain dimension not multiple of' //                &
               & ' processor number')
    end if
    jbegin = 1
    jendl = jxp
    jendx = jxp
    jendm = jxp
#ifndef BAND
    if ( myid == 0 ) jbegin = 2
    if ( myid == nproc-1 ) then
      jendx = jxp - 1
      jendm = jxp - 2
    end if
#endif
!
!**********************************************************************
!
!     RegCM V4 printout header
!
!**********************************************************************
!
    call header(myid,nproc)
!
!**********************************************************************
!
!     Parameter Setup
!
!**********************************************************************
!
    call param
    dtinc = dt
!
!**********************************************************************
!
!     Read initial data
!
!**********************************************************************
!
! this below enable debugging
#ifdef DEBUG 
    call start_debug()
#endif 
!
    call init
#ifdef CHEMTEST
    if ( ichem == 1 ) then
      call init_chem
    end if
#endif
!
!**********************************************************************
!
!     Read Boundary conditions
!
!**********************************************************************
!
    call bdyin
#ifdef CHEMTEST
    if ( ichem == 1 ) then
      call bdyin_chem
    end if
#endif
!
    call spinit
! 
    if ( ichem == 1 ) call chsrfem
!
!**********************************************************************
!
!     Write initial state to output
!
!**********************************************************************
!
    call output
!
!**********************************************************************
!
!     Clean up and logging
!
!**********************************************************************
!
    call free_mpp_initspace
    call time_print(6,'inizialization phase')
    call time_reset()
!
    return

  end subroutine RCM_initialize
!
!=======================================================================
!                                                                      !
!     This routine runs RegCM model from specified starting (TimeStr)  !
!     to ending (TimeEnd) time-steps.                                  !
!                                                                      !
!=======================================================================
!
  subroutine RCM_run(timestr, timeend, first)
    implicit none
!
!**********************************************************************
!
!     Imported variable declarations 
!
!**********************************************************************
!
    real(8), intent(in) :: timestr   ! starting time-step
    real(8), intent(in) :: timeend   ! ending   time-step
    logical, intent(in) :: first
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
    real(8) :: extime
    integer :: iexec
!
!**********************************************************************
!
!     Initialization of parameters 
!
!**********************************************************************
!
    if ( first ) then
      extime = d_zero
      iexec  = 1
    end if
!
!**********************************************************************
!
!     Model run 
!
!**********************************************************************
!     
    do while ( extime >= timestr .and. extime < timeend)
!
!     Read in boundary conditions if needed
!
      if ( idatex == bdydate1 .and. idatex < idate2 ) then
        call bdyin
#ifdef CHEMTEST
        if ( ichem == 1 ) call bdyin_chem
#endif
      end if
!
!     Refined start
!
      if ( .not.ifrest ) then
        if ( rfstrt ) then
          if ( (ktau == 0) .or. dtinc /= deltmx ) then
            call tstep(extime,dtinc)
            write (aline, 99001) extime , dtinc , dt , dt2 ,          &
                               & dtmin , ktau , idatex%year
            call say
          end if
        end if
      end if
!
!     Compute tendencies
!
      call tend(iexec)
!
!     Split modes
!
      call splitf
!
!     Write output for this timestep if requested
!
      if (ifrest) done_restart = .true.
      call output
!
!     Increment time
!
      extime = extime + dtinc
      if (debug_level > 3) then
        if (myid == 0) then
          write(6,'(a,a,f12.2)') 'Simulation time: ', idatex%tostring(), extime
        end if
      end if
!
    end do

!this below close down debug 
#ifdef DEBUG
    call stop_debug()
#endif 
    call time_print(6,'evolution phase')
!
!**********************************************************************
!
!     Formats
!
!**********************************************************************
!
99001 format (6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,       &
        & ' dt = ',f7.1,' dt2 = ',f7.1,' dtmin = ',f6.1,' ktau = ', &
        & i7,' in year ',i4)

  end subroutine RCM_run

  subroutine RCM_finalize
    implicit none
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
    integer :: ierr
    type(rcm_time_interval) :: tdif
!
!**********************************************************************
!
!     Simulation completed
!
!**********************************************************************
!
    call release_mod_ncio
!
    write (aline, 99002) xtime , ktau , idatex%year
    call say
!
!**********************************************************************
!
!     Set length of next run (auto-restart option)
!
!**********************************************************************
!
    tdif = idate2 - idate1
    idate1 = idate2
    idate2 = idate1 + tdif
    write (aline, *) ' *** new max DATE will be ' , idate2%tostring()
    call say
!
!**********************************************************************
!
!     Finalize the components 
!
!**********************************************************************
!
    if ( myid == 0 ) then
      call for_next
    end if
    call mpi_finalize(ierr)
#ifdef CLM
    call t_prf('timing_all',mpicom)
    call t_finalizef()
#endif
!
    call memory_destroy
    call finaltime(myid)
!
    if ( myid == 0 ) then
      print *, 'RegCM V4 simulation successfully reached end'
    end if
!
!**********************************************************************
!
!     Formats
!
!**********************************************************************
!
99002 format (                                                          &
     & ' ***** restart file for next run is written at time     = ',&
     & f10.2,' minutes, ktau = ',i7,' in year ',i4)

    contains
!
!**********************************************************************
!
!    Subroutine to write file restparam.nl with an hint 
!    for restarting the model
!
!**********************************************************************
!
    subroutine for_next
      implicit none
! 
!**********************************************************************
!
!     Open file and write data
!
!**********************************************************************
!
      open(99, file='restparam.nl')
      write (99,99001) '&restartparam'
      if ( idate1 < globidate2 ) then
        write (99,99001) 'ifrest  = .true. '
      else
        write (99,99001) 'ifrest  = .false.'
      end if
      write (99,99002) 'idate0  = ' , idate0%toidate()
      write (99,99002) 'idate1  = ' , idate1%toidate()
      write (99,99002) 'idate2  = ' , globidate2%toidate()
      write (99,99002) '/'
! 
!**********************************************************************
!
!     Close file
!
!**********************************************************************
!
      close (99)
!
!**********************************************************************
!
!     Formats
!
!**********************************************************************
!
99001 format (1x,a)
99002 format (1x,a,i10,',')
!
    end subroutine for_next
!
  end subroutine RCM_finalize
!
end module mod_interface
