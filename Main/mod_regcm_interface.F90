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
module mod_regcm_interface

  use mod_memutil
  use mod_service
  use mod_che_interface
  use mod_lm_interface
  use mod_atm_interface
  use mod_pbl_interface
  use mod_rad_interface
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_ncio
  use mod_ncout
  use mod_output
  use mod_split
  use mod_bdycod
  use mod_init
  use mod_header
  use mod_params
  use mod_tendency
  use mod_service
  use mod_moloch
#ifdef CPL
  use mod_update, only: rcm_get, rcm_put
#endif
  use mpi
  implicit none

  private
  public :: RCM_initialize
  public :: RCM_run
  public :: RCM_finalize
  public :: atm_model

  real(rk8) :: extime

  type atm_model
    character(len=5) :: model_name = 'RegCM'
    character(len=31) :: model_longname = 'The ICTP Regional Climate Model'
  end type atm_model

  data extime /0.0_rk8/
  contains

  subroutine RCM_initialize(mpiCommunicator)
    implicit none
    integer, intent(in), optional :: mpiCommunicator
    integer(ik4) :: ierr
    !
    ! MPI Initialization
    !
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
    ! Read IC and BC data.
    !
    if ( .not. ifrest ) then
      if ( irceideal /= 1 ) then
        call init_bdy
      end if
    end if
    !
    ! Initialize data (from IC or restart)
    !
    call init
    !
    ! Initialize split explicit scheme ( hydrostatic )
    !
    if ( idynamic == 1 ) call spinit
    !
    ! Setup the output files
    !
    call init_output_streams(do_parallel_netcdf_out)
    call output
    !
    ! Setup valid BC's
    !
    if ( irceideal /= 1 ) call bdyval
    !
    ! Clean up and logging
    !
#ifdef DEBUG
    call time_print(6,'inizialization phase')
    call time_reset()
#endif
  end subroutine RCM_initialize
  !
  !=======================================================================
  !                                                                      !
  !     This routine runs RegCM model from specified starting (TimeStr)  !
  !     to ending (TimeEnd) time-steps.                                  !
  !                                                                      !
  !=======================================================================
  !
  subroutine RCM_run(timestr, timeend)
    implicit none
    real(rk8) , intent(in) :: timestr   ! starting time-step
    real(rk8) , intent(in) :: timeend   ! ending   time-step

    do while ( extime >= timestr .and. extime < timeend )
      !
      ! Retrieve information from the driver
      !
#ifdef CPL
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        if ( rcmtimer%integrating( ) ) then
          call rcm_get(myid)
        end if
      end if
#endif
      !
      ! Compute tendencies
      !
      if ( idynamic == 3 ) then
        call moloch
      else
        call tend
      end if
      !
      ! Write output for this timestep if requested
      !
      call output
      !
      ! Boundary code
      !
      if ( .not. rcmtimer%reached_endtime ) then
        if ( alarm_in_bdy%act( ) ) then
          !
          ! Read in new boundary conditions
          !
          if ( irceideal /= 1 ) call bdyin
        end if
        !
        ! fill up the boundary values for xxb and xxa variables:
        !
        if ( irceideal /= 1 ) call bdyval
      end if
      !
      ! Send information to the driver
      !
#ifdef CPL
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call rcm_put(myid)
      end if
#endif
      !
      ! Increment execution time and boundary time
      !
      extime = extime + real(dtsec,rk8)
      if ( debug_level > 3 ) then
        if ( myid == italk ) then
          write(6,'(a,a,f12.2)') 'Simulation time: ', rcmtimer%str( )
        end if
      end if

    end do

#ifdef DEBUG
    call time_print(6,'evolution phase')
    call stop_debug()
#endif

  end subroutine RCM_run

  subroutine RCM_finalize
    implicit none

    if ( myid == italk ) then
      write(stdout,*) 'Final time ', trim(rcmtimer%str( )) , ' reached.'
    end if

    call close_icbc
    if ( ichem == 1 ) call close_chbc( )
    call dispose_output_streams
    call checktime(myid,trim(dirout)//pthsep//trim(prestr)//trim(domname)// &
                       '.'//tochar10(lastout))

#ifdef CLM
    call t_prf('timing_all',mpicom)
    call t_finalizef()
#endif

    if ( iclimao3 == 1 ) then
      call closeo3
    end if

    if ( iclimaaer == 1 ) then
      call closeaerosol
    end if

    call rcmtimer%dismiss( )
    call memory_destroy
    call finaltime(myid)

    if ( myid == italk ) then
      write(stdout,*) 'RegCM V5 simulation successfully reached end'
    end if
  end subroutine RCM_finalize

end module mod_regcm_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
