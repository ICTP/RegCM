!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      program regcm
!
      use mod_regcm_param
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3 , only : ptop , deltmx , sigma
      use mod_message
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Local variables
!
      real(8) :: dtinc , extime
      integer :: iexec , iexecn
#ifdef MPP1
      integer :: ierr , ncpu
#else
      integer :: myid
#endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!     ---before runing this model, the following data sets must be
!     prepared:
!     1. initial data -- output from "ccm to mm4" routine,
!     2. lateral boundary data
!
!     ---for each forecast, the parameters needed to be changed are
!     in include file "parame", and namelist "mm4.in".
!     information on data disposal set in inlcude file "dispose.f"
!
!
!     parameters :
!
!     ix, jx, kx : dimensions of arrays in y, x, z directions
!     for large domain.
!
!     nx (=7) : seven "j" (north-south) slices are needed for
!     tau+1 variables.
!
!     nspgx (=6) : number of cross-point slices on the boundary
!     affected by sponge boundary conditions.
!
!     nspgd (=6) : number of dot-point slices on the boundary
!     affected by sponge boundary conditions.
!
!     common blocks :
!
!     /main/    : stores all the 3d prognostic variables in two
!     timesteps and all the 2d variavles and constants
!
!     /bdycod/  : stores boundary values and tendencies of p*u,
!     p*v, p*t, p*qv and p*, and outmost 2 slices of
!     u and v for large domain.
!
!
!     /cvaria/  : stores the prognostic variables at tau+1,
!     decoupled variables, diagnostic variables and
!     working spaces needed in the model.
!
!     /param1/  : stores the parameters specified in subroutine
!     "param" for large domain only.
!
!
!     /param2/  : stores the logicals and constants specified in
!     subroutine "param".
!
!     /param3/  : stores the indexes and constants specified in
!     subroutine "param".
!
!     /iunits/  : stores the input/output unit numbers specified in
!     subroutine "param".
!
!     /pmoist/  : stores the parameters and constants related to
!     moisture calculations.
!
!     /pbldim/  : stores the parameters and constants related to
!     the boundary layer
!
!     *** comments ***
!     1. all the variables stored in main common
!     block must be saved for restart.
!     2. the variables stored in common blocks:
!     main   ,
!     /bdycod/, /param1/
!     are equivalent to large arrays for simplicity
!     to transfer the variables through arguments.
!
!     references :
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!**********************************************************************

#ifdef MPP1
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,ncpu,ierr)
      if ( ncpu.ne.nproc ) then
        write (aline,*) 'The number of CPU is not well set'
        call say
        write (aline,*) 'NCPU = ' , ncpu , '    nproc =' , nproc
        call say
        call fatal(__FILE__,__LINE__,'CPU Count mismatch')
      end if
      print * , "process" , myid , "of" , nproc
      call mpi_barrier(mpi_comm_world,ierr)
!     starttime= MPI_WTIME()
      if ( myid.gt.0 ) then
        iwest = myid - 1
      else
        iwest = mpi_proc_null
      end if
      if ( myid.lt.nproc-1 ) then
        ieast = myid + 1
      else
        ieast = mpi_proc_null
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
                 & 'Domain dimension not multiple of'//                 &
                  &' processor number')
      end if
      jbegin = 1
      jendl = jxp
      jendx = jxp
      jendm = jxp
      if ( myid.eq.0 ) jbegin = 2
      if ( myid.eq.nproc-1 ) then
        jendx = jxp - 1
        jendm = jxp - 2
      end if
#else
      myid = 1
#endif
!!
      call header(myid)
!
!-----set up parameters:
!
      extime = 0.
      dtinc = 0.
      iexec = 1
      iexecn = 1
!
      call param
!
!-----read in initial data:
!
      call init
 
      call bdyin
!
      call spinit(ptop,sigma,kxp1)
!
!chem2
      if ( ichem.eq.1 ) call chsrfem
!chem2_
 
      call output(iexec)
      do
!
!-----begin forecast:
!
!
!-----read in boundary conditions:
!
        if ( nnnnnn.gt.nnbase ) call bdyin
!
!.....refined start:
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
        call tend(iexec)
!
        call splitf
!
!-----output:
!
        call output(iexec)
!
        extime = extime + dtinc
!       print*,nnnnnn,nnnend
        if ( nnnnnn.ge.nnnend ) then
          write (aline, 99002) xtime , ktau , jyear
          call say
!
!-----set length of next run (auto-restart option)
!
!         xchar = 't'
          idate1 = idate2
          idate2 = mdatez(nnnend+nslice)
          write (aline, *) ' *** new max DATE will be ' , idate2
          call say

#ifdef MPP1
          if ( myid.eq.0 ) then
            call for_next
          end if
          call mpi_finalize(ierr)
#else
          call for_next
#endif
!         endtime = MPI_WTIME()
!         print *,"The Program took  ",endtime-starttime," secondes"
          stop 99999
        end if
      end do

99001 format (6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,       &
             &' dt = ',f7.1,' dt2 = ',f7.1,' dtmin = ',f6.1,' ktau = ', &
             &i7,' in year ',i4)
99002 format (                                                          &
          &' ***** restart file for next run is written at time     = ',&
          &f10.2,' minutes, ktau = ',i7,' in year ',i4)

      end program regcm
