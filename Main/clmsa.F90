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
program clmsa

#ifdef CLM45
  use mod_realkinds
  use mod_intkinds
  use mod_date
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_service
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_header
  use mod_clm_params
  use mod_service
  use mod_date
  use mod_dynparam
  use mod_atm_stub
  use mod_clm_regcm
  use mod_ncio
#ifndef MPI_SERIAL
  use mpi
#endif
  implicit none

  real(rk8) :: extime
  real(rk8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer(ik4) :: ierr , iprov
#ifdef MPI_SERIAL
  include 'mpif.h'
  integer(ik4) , parameter :: mpi_thread_single = 0
#endif
  data extime /0.0_rk8/

  call mpi_init_thread(mpi_thread_single,iprov,ierr)
  if ( ierr /= mpi_success ) then
    write(stderr,*) 'Cannot initilize MPI'
    stop
  end if

  call CLM_initialize()

  timestr = d_zero
  tdif = idate2 - idate1
  timeend = tohours(tdif) * secph

  call CLM_run(timestr, timeend)

  call CLM_finalize()
  call mpi_finalize(ierr)

  contains

  subroutine CLM_initialize(mpiCommunicator)
    implicit none
    integer, intent(in), optional :: mpiCommunicator
    integer(ik4) :: ierr
    if (present(mpiCommunicator)) then
      call mpi_comm_dup(mpiCommunicator,mycomm,ierr)
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'Cannot get communicator!')
      end if
    else
      call mpi_comm_dup(MPI_COMM_WORLD,mycomm,ierr)
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'Cannot get communicator!')
      end if
    end if
    call mpi_comm_rank(mycomm, myid, ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'mpi_comm_rank Failure!')
    end if
    call mpi_comm_size(mycomm, nproc, ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'mpi_comm_size Failure!')
    end if
#ifndef MPI_SERIAL
#ifdef DEBUG
    call mpi_comm_set_errhandler(mycomm, mpi_errors_return, ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'mpi_comm_set_errhandler Failure!')
    end if
#endif
#endif

    call whoami(myid)
    call setup_mesg(myid)

#ifdef DEBUG
    call activate_debug()
#endif
    !
    ! Read input global namelist
    !
    if ( myid == iocpu ) then
      call get_command_argument(0,value=prgname)
      call get_command_argument(1,value=namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr /= 0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
    end if

    call broadcast_params

    call memory_init

    call header(myid,nproc)
    call set_nproc
    call setup_model_indexes

#ifdef DEBUG
    call start_debug()
#endif
    !
    ! Parameter Setup
    !
    call param
    !
    ! INIT
    !
    call init_bdy
    call atmval
    if ( sfbcread == 1 ) then
      call read_clmbc(pptc,solar,flwd,totc)
      swdif = solar * 0.75_rkx * totc
      swdir = solar - swdif
      lwdif = flwd * 0.75_rkx * totc
      lwdir = flwd - lwdif
    end if
    call initsaclm45(lm)
    !
    ! Clean up and logging
    !
#ifdef DEBUG
    call time_print(6,'inizialization phase')
    call time_reset()
#endif
  end subroutine CLM_initialize
  !
  !=======================================================================
  !                                                                      !
  !     This routine runs CLM45 model from specified starting (TimeStr)  !
  !     to ending (TimeEnd) time-steps.                                  !
  !                                                                      !
  !=======================================================================
  !
  subroutine CLM_run(timestr, timeend)
    implicit none
    real(rk8) , intent(in) :: timestr   ! starting time-step
    real(rk8) , intent(in) :: timeend   ! ending   time-step

    do while ( extime >= timestr .and. extime < timeend )

      call atmval
      call runsaclm45(lm)
      call rcmtimer%advance( )

      if ( .not. rcmtimer%reached_endtime ) then
        if ( alarm_in_bdy%act( ) ) then
          !
          ! Read in new boundary conditions
          !
          call bdyin
        end if
        if ( sfbcread == 1 ) then
          if ( alarm_hour%act( ) ) then
            call read_clmbc(pptc,solar,flwd,totc)
            swdif = solar * 0.75_rkx * totc
            swdir = solar - swdif
            lwdif = flwd * 0.75_rkx * totc
            lwdir = flwd - lwdif
          end if
        end if
      end if
      !
      ! Increment execution time and boundary time
      !
      extime = extime + real(dtsrf,rk8)
      if ( debug_level > 3 ) then
        if ( myid == italk ) then
          write(6,'(a,a,f12.2)') 'Simulation time: ', rcmtimer%str( ), extime
        end if
      end if
    end do

#ifdef DEBUG
    call time_print(6,'evolution phase')
    call stop_debug()
#endif

  end subroutine CLM_run

  subroutine CLM_finalize
    implicit none

    if ( myid == italk ) then
      write(stdout,*) 'Final time ', trim(rcmtimer%str( )) , ' reached.'
    end if

    call close_icbc
    if ( sfbcread == 1 ) then
      call close_clmbc
    end if

    call rcmtimer%dismiss( )
    call memory_destroy
    call finaltime(myid)

    if ( myid == italk ) then
      write(stdout,*) 'CLM45 simulation successfully reached end'
    end if
  end subroutine CLM_finalize

#else
  write(0,*) 'This programs is enabled only if CLM45 is compiled in.'
#endif

end program clmsa
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
