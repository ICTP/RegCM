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
  use mod_sun, only : zenitm
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
    logical :: lsaved
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
    ! atmval advances the boundary interpolation time. This call only
    ! provides the fields for the initialization: rewind, or the first
    ! interval of the run is interpolated one step ahead (and a restarted
    ! run does not reproduce the continuous one).
    xbctime = d_zero
    if ( sfbcread == 1 ) then
      call getmem(ssnet,jci1,jci2,ici1,ici2,'clmsa:ssnet')
      call read_clmbc(pptc,ssnet,flwd,totc)
    end if
    call initsaclm45(lm)
    if ( sfbcread == 1 ) then
      ! On restart, the albedos kept by albedoclm45 through the night are
      ! taken from the save file written along with the CLM restart.
      lsaved = .false.
      if ( ifrest ) lsaved = read_sa_albedo(tochar10(rcmtimer%idate))
      if ( .not. lsaved ) call first_albedo
      call split_solar
    end if
    !
    ! The restart flag forces the solar orbit computation at every call of
    ! zenitm. After the initialization it must follow the same (daily)
    ! update as in a continuous run, or a restart does not reproduce it.
    !
    doing_restart = .false.
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
    logical :: lrstwr
    character(len=16) :: rstdate

    do while ( extime >= timestr .and. extime < timeend )

      call update_scenario(rcmtimer%year,rcmtimer%month)
      !
      ! Solar declination (used by CLM for the daylength and the
      ! photosynthetic capacity) and zenith angle, updated as in the
      ! coupled model. Without this the declination stays the one of the
      ! start date of the run.
      !
      call zenitm(mddom%xlat,mddom%xlon,coszrs)
      call atmval
      call runsaclm45(lm,lrstwr,rstdate)
      call albedoclm45(lm,lms)
      if ( sfbcread == 1 .and. lrstwr ) then
        call write_sa_albedo(rstdate)
      end if
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
  ! Save file of the standalone driver. albedoclm45 keeps through the night
  ! the last albedos with the sun above the horizon, and split_solar uses
  ! them when the forcing has sunlight and CLM has not (at sunrise). This
  ! state is not CLM state (the coupled model keeps it in its SAV file),
  ! so it is saved here, at the same time and date of the CLM restart.
  !
  character(len=256) function sa_savefile(cdate)
    implicit none
    character(len=*), intent(in) :: cdate
    sa_savefile = trim(dirout)//trim(domname)//'.clmsa.sav.'// &
                  trim(cdate)//'.nc'
  end function sa_savefile

  subroutine sa_check(istat,what,fname)
    use netcdf
    implicit none
    integer(ik4), intent(in) :: istat
    character(len=*), intent(in) :: what, fname
    if ( istat /= nf90_noerr ) then
      write(stderr,*) trim(what), ' ', trim(fname), ' : ', &
                      trim(nf90_strerror(istat))
      call fatal(__FILE__,__LINE__,'CLMSA SAVE FILE')
    end if
  end subroutine sa_check

  subroutine sa_albedo_pointer(n,p)
    implicit none
    integer(ik4), intent(in) :: n
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(out) :: p
    select case (n)
      case (1)
        p => lms%swdiralb
      case (2)
        p => lms%swdifalb
      case (3)
        p => lms%lwdiralb
      case default
        p => lms%lwdifalb
    end select
  end subroutine sa_albedo_pointer

  subroutine write_sa_albedo(cdate)
    use netcdf
    implicit none
    character(len=*), intent(in) :: cdate
    character(len=8), dimension(4), parameter :: vname = &
      [character(len=8) :: 'swdiralb', 'swdifalb', 'lwdiralb', 'lwdifalb']
    real(rkx), pointer, contiguous, dimension(:,:,:) :: gbuf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: p => null( )
    character(len=256) :: fname
    integer(ik4) :: ncid, istat, n
    integer(ik4), dimension(3) :: idims
    integer(ik4), dimension(4) :: ivar
    fname = sa_savefile(cdate)
    if ( myid == iocpu ) then
      call getmem(gbuf,1,nnsg,jcross1,jcross2,icross1,icross2,'clmsa:gbuf')
      istat = nf90_create(fname,ior(nf90_clobber,nf90_netcdf4),ncid)
      call sa_check(istat,'Cannot create',fname)
      istat = nf90_def_dim(ncid,'nnsg',nnsg,idims(1))
      call sa_check(istat,'Cannot define nnsg in',fname)
      istat = nf90_def_dim(ncid,'jx',jcross2-jcross1+1,idims(2))
      call sa_check(istat,'Cannot define jx in',fname)
      istat = nf90_def_dim(ncid,'iy',icross2-icross1+1,idims(3))
      call sa_check(istat,'Cannot define iy in',fname)
      do n = 1, 4
        istat = nf90_def_var(ncid,vname(n),nf90_double,idims,ivar(n))
        call sa_check(istat,'Cannot define '//vname(n)//' in',fname)
      end do
      istat = nf90_put_att(ncid,nf90_global,'title', &
        'clmsa albedos kept through the night, for restart')
      istat = nf90_put_att(ncid,nf90_global,'date',trim(cdate))
      istat = nf90_enddef(ncid)
      call sa_check(istat,'Cannot end definition of',fname)
    end if
    do n = 1, 4
      call sa_albedo_pointer(n,p)
      call subgrid_collect(p,gbuf,jci1,jci2,ici1,ici2)
      if ( myid == iocpu ) then
        istat = nf90_put_var(ncid,ivar(n),gbuf)
        call sa_check(istat,'Cannot write '//vname(n)//' in',fname)
      end if
    end do
    if ( myid == iocpu ) then
      istat = nf90_close(ncid)
      call sa_check(istat,'Cannot close',fname)
      call relmem(gbuf)
      if ( myid == italk ) then
        write(stdout,*) 'Write clmsa save file ', trim(fname)
      end if
    end if
  end subroutine write_sa_albedo
  !
  ! Returns .false. (and the caller falls back on the default albedos) if
  ! the save file is not there, as for restart files of older versions.
  !
  logical function read_sa_albedo(cdate) result(lread)
    use netcdf
    implicit none
    character(len=*), intent(in) :: cdate
    character(len=8), dimension(4), parameter :: vname = &
      [character(len=8) :: 'swdiralb', 'swdifalb', 'lwdiralb', 'lwdifalb']
    real(rkx), pointer, contiguous, dimension(:,:,:) :: gbuf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: p => null( )
    character(len=256) :: fname
    integer(ik4) :: ncid, istat, n, ivar, idim, ilen
    integer(ik4), dimension(3) :: nexp
    character(len=4), dimension(3), parameter :: dname = &
      [character(len=4) :: 'nnsg', 'jx', 'iy']
    fname = sa_savefile(cdate)
    if ( myid == iocpu ) then
      inquire(file=fname,exist=lread)
    end if
    call bcast(lread)
    if ( .not. lread ) then
      if ( myid == italk ) then
        write(stdout,*) 'No clmsa save file ', trim(fname), &
                        ': default albedos until the first sunrise.'
      end if
      return
    end if
    if ( myid == iocpu ) then
      call getmem(gbuf,1,nnsg,jcross1,jcross2,icross1,icross2,'clmsa:gbuf')
      istat = nf90_open(fname,nf90_nowrite,ncid)
      call sa_check(istat,'Cannot open',fname)
      nexp = [nnsg, jcross2-jcross1+1, icross2-icross1+1]
      do n = 1, 3
        istat = nf90_inq_dimid(ncid,dname(n),idim)
        call sa_check(istat,'Missing dimension '//dname(n)//' in',fname)
        istat = nf90_inquire_dimension(ncid,idim,len=ilen)
        call sa_check(istat,'Cannot read dimension '//dname(n)//' in',fname)
        if ( ilen /= nexp(n) ) then
          write(stderr,*) 'Dimension ', dname(n), ' is ', ilen, &
                          ', expected ', nexp(n), ' in ', trim(fname)
          call fatal(__FILE__,__LINE__,'CLMSA SAVE FILE')
        end if
      end do
    end if
    do n = 1, 4
      if ( myid == iocpu ) then
        istat = nf90_inq_varid(ncid,vname(n),ivar)
        call sa_check(istat,'Missing '//vname(n)//' in',fname)
        istat = nf90_get_var(ncid,ivar,gbuf)
        call sa_check(istat,'Cannot read '//vname(n)//' in',fname)
      end if
      call sa_albedo_pointer(n,p)
      call subgrid_distribute(gbuf,p,jci1,jci2,ici1,ici2)
    end do
    if ( myid == iocpu ) then
      istat = nf90_close(ncid)
      call sa_check(istat,'Cannot close',fname)
      call relmem(gbuf)
    end if
    where ( lm%ldmsk1 == 1 )
      lms%swalb = max(lms%swdiralb,lms%swdifalb)
      lms%lwalb = max(lms%lwdiralb,lms%lwdifalb)
    end where
    if ( myid == italk ) then
      write(stdout,*) 'Read clmsa save file ', trim(fname)
    end if
  end function read_sa_albedo
  !
  ! The forcing gives the net surface shortwave (ERA5 ssr). If the ERA5
  ! downward shortwave (ssrd) becomes available for download, clmbc should
  ! write it and this recovery is no longer needed. The downward one is the
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
