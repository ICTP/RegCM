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
!
#ifdef CLM45
program clmsa
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
  use mod_ipcc_scenario, only : update_scenario
  use mod_ncio
  use mpi_f08
  implicit none

  real(rk8) :: extime
  real(rk8) :: timestr, timeend
  type(rcm_time_interval) :: tdif
  integer(ik4) :: ierr, iprov
  ! ERA5 net surface shortwave, from which the downward one is recovered
  real(rkx), pointer, contiguous, dimension(:,:) :: ssnet => null( )
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
    type(mpi_comm) :: fcom
    if (present(mpiCommunicator)) then
      fcom%mpi_val = mpiCommunicator
      call mpi_comm_dup(fcom,mycomm,ierr)
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
#ifdef DEBUG
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
    else
      call init_random_number_generator( )
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
      call getmem(ssnet,jci1,jci2,ici1,ici2,'clmsa:ssnet')
      call read_clmbc(pptc,ssnet,flwd,totc)
    end if
    call initsaclm45(lm)
    if ( sfbcread == 1 ) then
      call first_albedo
      call split_solar
    end if
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
    real(rk8), intent(in) :: timestr   ! starting time-step
    real(rk8), intent(in) :: timeend   ! ending   time-step

    do while ( extime >= timestr .and. extime < timeend )

      call update_scenario(rcmtimer%year,rcmtimer%month)
      call atmval
      call runsaclm45(lm)
      call albedoclm45(lm,lms)
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
            call read_clmbc(pptc,ssnet,flwd,totc)
            call split_solar
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
      write(stdout,*) 'Final time ', trim(rcmtimer%str( )), ' reached.'
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
  !
  ! Albedos of the initial state, with the defaults of the coupled model
  ! where CLM has none (sun below the horizon)
  !
  subroutine first_albedo
    implicit none
    where ( lm%ldmsk1 == 1 )
      lms%swdiralb = 0.16_rkx
      lms%swdifalb = 0.16_rkx
      lms%lwdiralb = 0.32_rkx
      lms%lwdifalb = 0.32_rkx
    end where
    call albedoclm45(lm,lms)
  end subroutine first_albedo
  !
  ! The forcing gives the net surface shortwave. The downward one is the
  ! net divided by the absorbed fraction, from the CLM albedos of the
  ! last step (land points only). The downward shortwave is split in
  ! visible (sw) and near-infrared (lw) bands, half each, and each band
  ! in direct and diffuse with the empirical direct fractions of the CESM
  ! data atmosphere (CDEPS). The fractions depend on the downward flux,
  ! so the absorbed fraction is iterated.
  !
  subroutine split_solar
    implicit none
    integer(ik4) :: i, j, n, nl, iter
    integer(ik4), parameter :: niter = 4
    ! Highest albedo allowed in the division
    real(rkx), parameter :: maxalb = 0.95_rkx
    real(rkx) :: s, rvis, rnir, avisd, avisi, anird, aniri, absf
    do i = ici1, ici2
      do j = jci1, jci2
        ! Interpolation may leave round-off negatives
        ssnet(j,i) = max(ssnet(j,i),0.0_rkx)
        avisd = 0.0_rkx
        avisi = 0.0_rkx
        anird = 0.0_rkx
        aniri = 0.0_rkx
        nl = 0
        do n = 1, nnsg
          if ( lm%ldmsk1(n,j,i) == 1 ) then
            nl = nl + 1
            avisd = avisd + lms%swdiralb(n,j,i)
            avisi = avisi + lms%swdifalb(n,j,i)
            anird = anird + lms%lwdiralb(n,j,i)
            aniri = aniri + lms%lwdifalb(n,j,i)
          end if
        end do
        if ( nl > 0 ) then
          avisd = min(avisd/real(nl,rkx),maxalb)
          avisi = min(avisi/real(nl,rkx),maxalb)
          anird = min(anird/real(nl,rkx),maxalb)
          aniri = min(aniri/real(nl,rkx),maxalb)
        end if
        solar(j,i) = ssnet(j,i)
        do iter = 1, niter
          s = solar(j,i) * 0.5_rkx
          rvis = vis_direct(s)
          rnir = nir_direct(s)
          absf = 0.5_rkx * (rvis*(1.0_rkx-avisd) +           &
                            (1.0_rkx-rvis)*(1.0_rkx-avisi) + &
                            rnir*(1.0_rkx-anird) +           &
                            (1.0_rkx-rnir)*(1.0_rkx-aniri))
          solar(j,i) = ssnet(j,i) / absf
        end do
        s = solar(j,i) * 0.5_rkx
        rvis = vis_direct(s)
        rnir = nir_direct(s)
        swdir(j,i) = rvis * s
        swdif(j,i) = (1.0_rkx - rvis) * s
        lwdir(j,i) = rnir * s
        lwdif(j,i) = (1.0_rkx - rnir) * s
      end do
    end do
  end subroutine split_solar

  pure real(rkx) function vis_direct(s)
    implicit none
    real(rkx), intent(in) :: s
    vis_direct = min(0.99_rkx,max(0.17639_rkx + 0.00380_rkx*s - &
                 9.0039e-06_rkx*s**2 + 8.1351e-09_rkx*s**3,0.01_rkx))
  end function vis_direct

  pure real(rkx) function nir_direct(s)
    implicit none
    real(rkx), intent(in) :: s
    nir_direct = min(0.99_rkx,max(0.29548_rkx + 0.00504_rkx*s - &
                 1.4957e-05_rkx*s**2 + 1.4881e-08_rkx*s**3,0.01_rkx))
  end function nir_direct

end program clmsa

#else
program clmsa
  implicit none
  write(0,*) 'This programs is enabled only if CLM45 is compiled in.'
end program clmsa
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
