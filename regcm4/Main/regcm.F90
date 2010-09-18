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
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! This is the Regional Climatic Model RegCM from ICTP Trieste
!
!  References :
!
!     1. model:
!
!     Anthes, R. A., and T. T. Warner, 1978: Development of
!     hydrodynamic models suitable for air pollution and
!     other mesometeorological studies. Mon. Wea. Rev.,
!     106, 1045-1078.
!
!     2. cumulus parameterization :
!
!     Anthes, R. A., 1977: A cumulus parameterization scheme
!     utilizing a one-dimensional cloud model. Mon. Wea.
!     Rev., 105, 270-286.
!
!     Kuo, Y.-H., 1983: A diagnostic case study of the effects
!     of deep extratropical convection on the large-scale
!     temperature and moisture structure. PH.D. thesis,
!     Department of Meteorology, the Pennsylvania State
!     University, 222 pp.
!
!     Grell,
!
!     3. explicit moisture :
!
!     Hsie, E.-Y., 1983: Frontogenesis in a moist atmosphere.
!     PH.D. thesis, Department of Meteorology, the
!     Pennsylvania State University, 251 pp.
!
!     4. pbl parameterization :
!
!     Holtslag, De Bruijn and Pan - MWR - 8/90
!
!     5. radiation parameterization :
!
!     CCM2 radiation column model, bruce briegleb, jan '92
!
!     CCM3 radiation column model, NCAR/TN-422+PPR, Description
!     of the NCAR CCM,      J. T. Kiehl, J. Hack et al.,
!     introduced by   Filippo Giorgi, N. Keiichi, Yun Qian
!
!     CCM3.6.6 code introduced by Xunqiang Bi, 2000
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      program regcm
!
      use mod_dynparam
      use mod_date
      use mod_runparams
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
      use service_mod
#ifdef MPP1
      use mpi
#ifdef CLM
      use perf_mod
      use spmdMod, only: mpicom
#endif
#endif
      implicit none
!
      real(8) :: dtinc , extime
      integer :: iexec , iexecn
      integer :: nhours
      integer :: ierr
#ifdef MPP1
      integer :: ncpu
#endif
      character(256) :: namelistfile, prgname
! 
!**********************************************************************
!
#ifdef MPP1
!**********************************************************************
!
!     MPI Initialization
!
!**********************************************************************
!
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,ncpu,ierr)

#endif
!

      call activate_debug()
!**********************************************************************
!
!     Read input global namelist
!
!**********************************************************************
!
#ifdef MPP1
      if ( myid.eq.0 ) then
#endif

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
!
#ifdef BAND
      call init_mod_ncio(.true.)
#else
      call init_mod_ncio(.false.)
#endif
!
#ifdef MPP1
      end if
#endif

#ifdef MPP1
!
!**********************************************************************
!
!     MPI Initialization for model
!
!**********************************************************************
!
      call broadcast_params

      call set_nproc(ncpu)

      if ( ncpu.ne.nproc ) then
        write (aline,*) 'The number of CPU is not well set'
        call say
        write (aline,*) 'NCPU = ' , ncpu , '    nproc =' , nproc
        call say
        call fatal(__FILE__,__LINE__,'CPU Count mismatch')
      end if
!      print * , "process" , myid , "of" , nproc
      call mpi_barrier(mpi_comm_world,ierr)
!     starttime= MPI_WTIME()
      if ( myid.gt.0 ) then
        iwest = myid - 1
      else
#ifdef BAND
        iwest = nproc-1
#else
        iwest = mpi_proc_null
#endif
      end if
      if ( myid.lt.nproc-1 ) then
        ieast = myid + 1
      else
#ifdef BAND
        ieast = 0
#else
        ieast = mpi_proc_null
#endif
      end if
      if ( jxp.lt.2 ) then
        write (aline,*) 'The number of jxp must be greater than 1'
        call say
        write (aline,*) 'jxp = ' , jxp , '   jx = ' , jx
        call say
        call fatal(__FILE__,__LINE__,'Domain too small')
      end if
      if ( jxp*nproc.ne.jx ) then
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
      if ( myid.eq.0 ) jbegin = 2
      if ( myid.eq.nproc-1 ) then
        jendx = jxp - 1
        jendm = jxp - 2
      end if
#endif
#else
      myid = 0
      nproc= 1 
#endif

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
      extime = 0.
      dtinc  = 0.
      iexec  = 1
      iexecn = 1
      call param
!
!**********************************************************************
!
!     Read initial data and setup boundary condition
!
!**********************************************************************
!
      call init
! 
      call bdyin
! 
      call spinit(sigma,kzp1)
! 
      if ( ichem.eq.1 ) call chsrfem
! 
!**********************************************************************
!
!     Write initial state to output
!
!**********************************************************************
!
      call output
      call time_print(6,'inizialization phase')
!
!**********************************************************************
!
!     Time Loop : begin Forecast
!
!**********************************************************************
!
      do while ( nnnnnn.lt.nnnend )
!
!       Read in boundary conditions if needed
!
        if ( nnnnnn.gt.nnbase ) call bdyin
!
!       Refined start
!
        if ( .not.ifrest ) then
          if ( rfstrt ) then
            if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or.                 &
               & dtinc.ne.deltmx ) then
              call tstep(extime,dtinc,deltmx)
              write (aline, 99001) extime , dtinc , dt , dt2 ,          &
                                   & dtmin , ktau , jyear
              call say
            end if
          end if
        end if
!
!       Compute tendencies
!
        call tend(iexec)
!
!       Split modes
!
        call splitf
!
!       Write output for this timestep
!
        call output
!
!       Increment time
!
        extime = extime + dtinc

      end do
      call time_print(6,'evolution phase')
!
!**********************************************************************
!
!     Simulation completed
!
!**********************************************************************
!
      call release_mod_ncio
!
      write (aline, 99002) xtime , ktau , jyear
      call say
!


!**********************************************************************
!
!     Set length of next run (auto-restart option)
!
!**********************************************************************
!
      nhours = idatediff(idate2,idate1)
      idate1 = idate2
      call addhours(idate2,nhours)
      write (aline, *) ' *** new max DATE will be ' , idate2
      call say
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        call for_next
      end if
      call mpi_finalize(ierr)
#else
      call for_next
#endif
#ifdef CLM
      call t_prf('timing_all',mpicom)
      call t_finalizef()
#endif
!
      call finaltime(myid)
!
      if ( myid.eq.0 ) then
        print *, 'RegCM V4 simulation successfully reached end'
      end if
!
! Formats
!
99001 format (6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,       &
            & ' dt = ',f7.1,' dt2 = ',f7.1,' dtmin = ',f6.1,' ktau = ', &
            & i7,' in year ',i4)
99002 format (                                                          &
         & ' ***** restart file for next run is written at time     = ',&
         & f10.2,' minutes, ktau = ',i7,' in year ',i4)

      contains
!
! Subroutine to write file restparam.nl with an hint for restarting the
!    model
!
      subroutine for_next
! 
      open(99, file='restparam.nl')
      write (99,99001) '&restartparam'
      if ( idate1.lt.globidate2 ) then
        write (99,99001) 'ifrest  = .true. '
      else
        write (99,99001) 'ifrest  = .false.'
      end if
      write (99,99002) 'idate0  = ' , idate0
      write (99,99002) 'idate1  = ' , idate1
      write (99,99002) 'idate2  = ' , globidate2
      write (99,99002) '/'
      close (99)
!
99001 format (1x,a)
99002 format (1x,a,i10,',')
!
      end subroutine for_next
!
      end program regcm
