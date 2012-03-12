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
  use mod_mppgrid
  use mod_service
  use mod_che_interface
  use mod_atm_interface
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_ncio
  use mod_output
  use mod_split
  use mod_bdycod
  use mod_init
  use mod_header
  use mod_params
  use mod_tendency
  use mod_tstep
  use mod_service
  use mod_mppio
  use mpi
  use mod_sun
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

  real(dp) :: dtinc
  real(dp) :: extime

  data extime /d_zero/
  contains
 
  subroutine RCM_initialize(mpiCommunicator)
    implicit none
    integer, intent(in), optional :: mpiCommunicator
!
    integer :: ncpu, ierr
    character(256) :: namelistfile, prgname
! 
!**********************************************************************
!
!   MPI Initialization
!
!**********************************************************************
!
    if (present(mpiCommunicator)) then
      mycomm = mpiCommunicator
    else
      mycomm = MPI_COMM_WORLD
    end if
    call mpi_comm_rank(mycomm, myid, ierr)
    call mpi_comm_size(mycomm, ncpu, ierr)
!
    call whoami(myid)
!
#ifdef DEBUG 
    call activate_debug()
#endif
!
!**********************************************************************
!
!   Read input global namelist
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
    call broadcast_params

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
!   MPI Initialization for model
!
!**********************************************************************
!
    call set_nproc(ncpu)

    if ( ncpu /= nproc ) then
      write (aline,*) 'The number of CPU is not well set'
      call say
      write (aline,*) 'NCPU = ' , ncpu , '    nproc =' , nproc
      call say
      call fatal(__FILE__,__LINE__,'CPU Count mismatch')
    end if
!      print * , "process" , myid , "of" , nproc
    call mpi_barrier(mycomm,ierr)
!     starttime= MPI_WTIME()
    if ( nproc == 1 ) then
#ifdef BAND
      iwest = myid
      ieast = myid
#else
      ieast = mpi_proc_null
      iwest = mpi_proc_null
#endif
    else
      if ( myid == 0 ) then
#ifdef BAND
        iwest = nproc-1
#else
        iwest = mpi_proc_null
#endif
        ieast = myid+1
      else if ( myid == nproc-1 ) then
#ifdef BAND
        ieast = 0
#else
        ieast = mpi_proc_null
#endif
        iwest = myid-1
      else
        ieast = myid+1
        iwest = myid-1
      end if
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
!   RegCM V4 printout header
!
!**********************************************************************
!
    call header(myid,nproc)
!
#ifdef DEBUG 
    call start_debug()
#endif 
!
!**********************************************************************
!
!   Parameter Setup
!
!**********************************************************************
!
    call param
    dtinc = dt
!
!**********************************************************************
!
!   Read initial data and boundary conditions
!
!**********************************************************************
!
    !
    ! Calculate solar declination angle at startup
    !
    if (myid == 0) then
      write (6,*) 'Calculate solar declination angle at ',toint10(idatex)
    end if
    call solar1
    call init_bdy
!
!**********************************************************************
!
!   Initialize data (from IC or restart)
!
!**********************************************************************
!
    call init
!
    if ( ichem == 1 ) then
      call start_chem(ice1,ice2,jce1,jce2,ifrest,bdydate1,bdydate2)
    end if
!
    if ( .not. ifrest ) then
      if ( ichem == 1 ) then
        call chem_bdyin(150D00, bdydate1, bdydate2)
      end if
    end if
!
!**********************************************************************
!
!   Initialize split explicit scheme
!
!**********************************************************************
!
    call spinit
!
!**********************************************************************
!
!   Initialize emission dataset
!
!**********************************************************************
!
    if ( ichem == 1 ) call chem_emission(xmonth)
!
!**********************************************************************
!
!   Setup the output files
!
!**********************************************************************
!
    call output
!
!**********************************************************************
!
!   Set the boundary conditions for this timestep
!
!**********************************************************************
!
    call bdyval(xbctime)
    if ( ichem == 1 ) then
      call chem_bdyval(xbctime,nbdytime,dtbdys,ktau,ifrest)
    end if
!
!**********************************************************************
!
!   Calculate Zenital Angle
!
!**********************************************************************
!
    call zenitm(coszrs,jci1,jci2,ici1,ici2)
!
!**********************************************************************
!
!   Clean up and logging
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
  subroutine RCM_run(timestr, timeend)
    implicit none
    real(dp) , intent(in) :: timestr   ! starting time-step
    real(dp) , intent(in) :: timeend   ! ending   time-step
    character(len=32) :: appdat
!
#ifdef DEBUG
    call deco1d_nc_create('psa',cross,sfs%psa,psa)
    call deco1d_nc_create('psb',cross,sfs%psb,psb)
    call deco1d_nc_create('ua',dot,atm1%u,uax)
    call deco1d_nc_create('va',dot,atm1%v,vax)
    call deco1d_nc_create('ta',cross,atm1%t,tax)
    call deco1d_nc_create('qa',cross,atm1%qv,qax)

    call deco1d_nc_create('tten',cross,aten%t,taten)
    call deco1d_nc_create('heatrt',cross,heatrt,aheat)

    call deco1d_nc_write(psa)
    call deco1d_nc_write(psb)
    call deco1d_nc_write(uax)
    call deco1d_nc_write(vax)
    call deco1d_nc_write(tax)
    call deco1d_nc_write(qax)
#endif

    do while ( extime >= timestr .and. extime < timeend)
      !
      ! Refined start
      !
      if ( .not. ifrest ) then
        if ( rfstrt ) then
          if ( (ktau == 0) .or. dtinc /= deltmx ) then
            call tstep(extime,dtinc)
            write (aline, 99001) extime , dtinc , dt , dt2 ,          &
                                 dtsec , ktau , xyear
            call say
          end if
        end if
      end if
      !
      ! Compute tendencies
      !
      call tend
      !
      ! Split modes
      !
      call splitf
      !
      ! Boundary code (do not execute at the end of run)
      !
      if ( ktau /= mtau ) then
        if ( nbdytime == 0 ) then
          !
          ! recalculate solar declination angle if reading bdy
          !
          if (myid == 0) then
            write (6,*) 'Calculate solar declination angle at ',toint10(idatex)
          end if
          call solar1
          !
          ! Read in new boundary conditions
          !
          call bdyin
          if ( ichem == 1 ) call chem_bdyin(150D00, bdydate1, bdydate2)
        end if
        !
        ! fill up the boundary values for xxb and xxa variables:
        !
        call bdyval(xbctime)
        if ( ichem == 1 ) then
          call chem_bdyval(xbctime,nbdytime,dtbdys,ktau,ifrest)
        end if
      end if
      !
      ! Write output for this timestep if requested
      !
      call output
      !
      ! Increment execution time
      !
      extime = extime + dtinc
      if (debug_level > 3) then
        if (myid == 0) then
          appdat = tochar(idatex)
          write(6,'(a,a,f12.2)') 'Simulation time: ', appdat, extime
        end if
      end if
    end do

#ifdef DEBUG
    call stop_debug()
    call deco1d_nc_write(psa)
    call deco1d_nc_write(psb)
    call deco1d_nc_write(uax)
    call deco1d_nc_write(vax)
    call deco1d_nc_write(tax)
    call deco1d_nc_write(qax)
    call deco1d_nc_destroy(psa)
    call deco1d_nc_destroy(psb)
    call deco1d_nc_destroy(uax)
    call deco1d_nc_destroy(vax)
    call deco1d_nc_destroy(tax)
    call deco1d_nc_destroy(qax)
    call deco1d_nc_destroy(taten)
#endif
    call time_print(6,'evolution phase')
!
99001 format (6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,       &
        & ' dt = ',f7.1,' dt2 = ',f7.1,' dtsec = ',f6.1,' ktau = ', &
        & i7,' in year ',i4)

  end subroutine RCM_run

  subroutine RCM_finalize
    implicit none
    character(len=32) :: appdat
!
    call release_mod_ncio
!
    appdat = tochar(idate2)
    write (aline, 99002) appdat
    call say
!
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
99002 format ('Restart file for next run is written at time = ',a)
!
  end subroutine RCM_finalize
!
end module mod_regcm_interface
