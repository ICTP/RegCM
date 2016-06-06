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
  use mod_cloud_s1
#ifdef CPL
  use mod_update, only: rcm_get, rcm_put
#endif
  use mpi
  implicit none

  private
  public :: RCM_initialize
  public :: RCM_run
  public :: RCM_finalize

  real(rkx) :: extime

  data extime /d_zero/
  contains

  subroutine RCM_initialize(mpiCommunicator)
    implicit none
    integer, intent(in), optional :: mpiCommunicator
    integer(ik4) :: ierr
    !
    ! MPI Initialization
    !
    if (present(mpiCommunicator)) then
      mycomm = mpiCommunicator
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
    call mpi_comm_set_errhandler(mycomm, mpi_errors_return, ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'mpi_comm_set_errhandler Failure!')
    end if
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
    if ( .not. ifrest ) call init_bdy
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
    call bdyval
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
    real(rkx) , intent(in) :: timestr   ! starting time-step
    real(rkx) , intent(in) :: timeend   ! ending   time-step
    character(len=32) :: appdat

#ifdef DEBUG
    ! if ( ipptls == 2 ) call grid_nc_create('qqxp',cross,zqxn,qqxp)
    ! call grid_nc_create('qxatm',cross,atm1%qx,nc_4d)
#endif
    do while ( extime >= timestr .and. extime < timeend )
#ifdef DEBUG
      ! call grid_nc_write(nc_4d)
#endif
      !
      ! Retrieve information from the driver
      !
#ifdef CPL
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        if (ktau > 0) then
          call rcm_get(myid)
        end if
      end if
#endif
      !
      ! Compute tendencies
      !
      call tend
      !
      ! Boundary code
      !
      if ( ktau /= mtau ) then
        if ( nbdytime == kbdy ) then
          !
          ! Read in new boundary conditions
          !
          call bdyin
        else
          xbctime = xbctime + dtsec
        end if
        !
        ! fill up the boundary values for xxb and xxa variables:
        !
        call bdyval
      end if
      !
      ! Write output for this timestep if requested
      !
      call output
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
      extime = extime + dtsec
      if ( debug_level > 3 ) then
        if ( myid == italk ) then
          appdat = tochar(idatex)
          write(6,'(a,a,f12.2)') 'Simulation time: ', appdat, extime
        end if
      end if
    end do

#ifdef DEBUG
    call stop_debug()
    ! if ( ipptls == 2 ) call grid_nc_destroy(qqxp)
    ! call grid_nc_destroy(nc_4d)
    call time_print(6,'evolution phase')
#endif

  end subroutine RCM_run

  subroutine RCM_finalize
    implicit none
    character(len=32) :: appdat

    if ( myid == italk ) then
      appdat = tochar(idate2)
      write(stdout,*) 'Restart file for next run is written at time ',appdat
    end if

    call close_icbc
    if ( ichem == 1 ) call close_chbc
    call dispose_output_streams

#ifdef CLM
    call t_prf('timing_all',mpicom)
    call t_finalizef()
#endif

    call memory_destroy
    call finaltime(myid)

    if ( myid == italk ) then
      write(stdout,*) 'RegCM V4 simulation successfully reached end'
    end if
  end subroutine RCM_finalize

end module mod_regcm_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
