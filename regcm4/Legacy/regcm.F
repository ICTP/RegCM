C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      program regcm
      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c                                                                     c
c     ---before runing this model, the following data sets must be    c
c        prepared:                                                    c
c        1. initial data -- output from "ccm to mm4" routine,         c
c        2. lateral boundary data                                     c
c                                                                     c
c     ---for each forecast, the parameters needed to be changed are   c
c        in include file "parame", and namelist "mm4.in".             c
c        information on data disposal set in inlcude file "dispose.f" c
c                                                                     c
c                                                                     c
c     parameters :                                                    c
c                                                                     c
c        ix, jx, kx : dimensions of arrays in y, x, z directions      c
c                     for large domain.                               c
c                                                                     c
c        nx (=7) : seven "j" (north-south) slices are needed for      c
c                  tau+1 variables.                                   c
c                                                                     c
c        nspgx (=6) : number of cross-point slices on the boundary    c
c                     affected by sponge boundary conditions.         c
c                                                                     c
c        nspgd (=6) : number of dot-point slices on the boundary      c
c                     affected by sponge boundary conditions.         c
c                                                                     c
c     common blocks :                                                 c
c                                                                     c
c       /main/    : stores all the 3d prognostic variables in two     c
c                   timesteps and all the 2d variavles and constants  c
c                                                                     c
c       /bdycod/  : stores boundary values and tendencies of p*u,     c
c                   p*v, p*t, p*qv and p*, and outmost 2 slices of    c
c                   u and v for large domain.                         c
c                                                                     c
c                                                                     c
c       /cvaria/  : stores the prognostic variables at tau+1,         c
c                   decoupled variables, diagnostic variables and     c
c                   working spaces needed in the model.               c
c                                                                     c
c       /param1/  : stores the parameters specified in subroutine     c
c                   "param" for large domain only.                    c
c                                                                     c
c                                                                     c
c       /param2/  : stores the logicals and constants specified in    c
c                   subroutine "param".                               c
c                                                                     c
c       /param3/  : stores the indexes and constants specified in     c
c                   subroutine "param".                               c
c                                                                     c
c       /iunits/  : stores the input/output unit numbers specified in c
c                   subroutine "param".                               c
c                                                                     c
c       /pmoist/  : stores the parameters and constants related to    c
c                   moisture calculations.                            c
c                                                                     c
c       /pbldim/  : stores the parameters and constants related to    c
c                   the boundary layer                                c
c                                                                     c
c                   *** comments ***                                  c
c                   1. all the variables stored in main common        c
c                      block must be saved for restart.               c
c                   2. the variables stored in common blocks:         c
c                      main   ,                                       c
c                     /bdycod/, /param1/                              c
c                      are equivalent to large arrays for simplicity  c
c                      to transfer the variables through arguments.   c
c                                                                     c
c     references :                                                    c
c                                                                     c
c        1. model:                                                    c
c                                                                     c
c           Anthes, R. A., and T. T. Warner, 1978: Development of     c
c                hydrodynamic models suitable for air pollution and   c
c                other mesometeorological studies. Mon. Wea. Rev.,    c
c                106, 1045-1078.                                      c
c                                                                     c
c        2. cumulus parameterization :                                c
c                                                                     c
c           Anthes, R. A., 1977: A cumulus parameterization scheme    c
c                utilizing a one-dimensional cloud model. Mon. Wea.   c
c                Rev., 105, 270-286.                                  c
c                                                                     c
c           Kuo, Y.-H., 1983: A diagnostic case study of the effects  c
c                of deep extratropical convection on the large-scale  c
c                temperature and moisture structure. PH.D. thesis,    c
c                Department of Meteorology, the Pennsylvania State    c
c                University, 222 pp.                                  c
c                                                                     c
c           Grell,                                                    c
c                                                                     c
c        3. explicit moisture :                                       c
c                                                                     c
c           Hsie, E.-Y., 1983: Frontogenesis in a moist atmosphere.   c
c                PH.D. thesis, Department of Meteorology, the         c
c                Pennsylvania State University, 251 pp.               c
c                                                                     c
c        4. pbl parameterization :                                    c
c                                                                     c
c           Holtslag, De Bruijn and Pan - MWR - 8/90                  c
c                                                                     c
c        5. radiation parameterization :                              c
c                                                                     c
c           CCM2 radiation column model, bruce briegleb, jan '92      c
c                                                                     c
c           CCM3 radiation column model, NCAR/TN-422+PPR, Description c
c             of the NCAR CCM,      J. T. Kiehl, J. Hack et al.,      c
c             introduced by   Filippo Giorgi, N. Keiichi, Yun Qian    c
c                                                                     c
c           CCM3.6.6 code introduced by Xunqiang Bi, 2000             c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#     include "parameter.inc"

#     include "include/date.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/message.cb"

#ifdef MPP1
c     real(kind=8)  starttime,endtime
      include 'mpif.h'
      integer NCPU
      integer ierr
#endif
c
      INTEGER         MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
c
      real(kind=8)  extime,dtinc
      integer iexec,iexecn
c
c**********************************************************************
#ifdef MPP1
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,NCPU,ierr)
      if (NCPU.ne.nproc) then
        write(aline,*) 'The number of CPU is not well set'
        call say
        write(aline,*) 'NCPU = ', NCPU, '    nproc =',nproc
        call say
        call fatal(__FILE__,__LINE__,'CPU Count mismatch')
      end if
      print *,"process",myid,"of",nproc
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     starttime= MPI_WTIME()
      if (myid .gt. 0) then
        iwest=myid-1
      else
        iwest=MPI_PROC_NULL
      end if
      if (myid .lt. nproc-1) then
        ieast=myid+1
      else
        ieast=MPI_PROC_NULL
      end if
      if (jxp.lt.2) then
        write(aline,*) 'The number of jxp must be greater than 1'
        call say
        write(aline,*) 'jxp = ', jxp, '   mjx = ',mjx
        call say
        call fatal(__FILE__,__LINE__,'Domain too small')
      end if
      if (jxp*nproc .ne. mjx) then
        write(aline,*) 'mjx should be divided by nproc'
        call say
        write(aline,*) 'mjx = ',mjx, '   nproc = ',nproc
        call say
        call fatal(__FILE__,__LINE__,'Domain dimension not multiple of'
     &             //' processor number')
      end if
      jbegin = 1
      jendl = jxp
      jendx = jxp
      jendm = jxp
      if (myid .eq. 0) jbegin = 2
      if (myid .eq. nproc-1) then
        jendx = jxp-1
        jendm = jxp-2
      end if
#else  
      integer :: myid 
      myid = 1
#endif
c! 
      call header(myid)
c
c-----set up parameters:
c
      extime = 0.
      dtinc = 0.
      iexec = 1
      iexecn = 1
c
      call param
c
c-----read in initial data:
c
      call init

      call bdyin
c
      call spinit( ptop, sigma, kxp1)
c
Cchem2  
      if (ichem.eq.1) then 
       call chsrfem 
      end if
Cchem2_

      call output(iexec)
c
c-----begin forecast:
c
10    continue
c
c-----read in boundary conditions:
c
      if (NNNNNN .gt. NNBASE) then
       call bdyin
      end if
c
c.....refined start:
c
      if (.not. ifrest) then
       if (rfstrt) then
        if ((jyear.eq.jyear0.and.ktau.eq.0) .or. dtinc.ne.deltmx) then
         call tstep(extime,dtinc,deltmx)
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         print 30, extime,dtinc,dt,dt2,dtmin,ktau,jyear
30       format(6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,
     1          ' dt = ',f7.1,' dt2 = ',f7.1,' dtmin = ',f6.1,
     2          ' ktau = ',i7,' in year ',i4)
#ifdef MPP1
         end if
#endif
        end if
       end if
      end if
c
      call tend(iexec)
c
      call splitf( gnuhf )
c
c-----output:
c
      call output(iexec)
c
      extime=extime+dtinc
c     print*,nnnnnn,nnnend
      if (NNNNNN.lt.NNNEND) goto 10
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      print 60, xtime,ktau,jyear
60    format(//' ***** restart file for next run is written at time
     1= ',f10.2,' minutes, ktau = ',i7,' in year ',i4)
#ifdef MPP1
      end if
#endif
c
c-----set length of next run (auto-restart option)
c
c     xchar = 't'
      IDATE1= IDATE2
      IDATE2= MDATEZ(NNNEND+nslice)
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      print*,' *** new max DATE will be ',IDATE2
      call for_next

#ifdef MPP1
      end if
c     endtime = MPI_WTIME()
c     print *,"The Program took  ",endtime-starttime," secondes"
      call MPI_FINALIZE(ierr)
#endif
      stop 99999
      end program regcm
