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

#ifdef MPP1

#ifdef CLM

      module mod_mtrxclm

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_date
      use mod_clm
      use mod_bats
      use mod_mppio
      use mod_main
      use mod_pbldim
      use mod_slice
      use mod_bats
      use mod_vecbats
      use mod_drag
      use mod_zengocn

      private

      public :: mtrxclm
      public :: initclm
      public :: interfclm
      public :: albedoclm

      contains

      subroutine mtrxclm
!
!=======================================================================
!  based on: clm version 3.5
!=======================================================================
! Matrix CLM was written by Jason Bell for coupling the clm3.0 land
! surface model to RegCM3.
!! 2008.5.8    A Tawfik Revised to work with RegCM3 and clm3.5
!
! mtrxclm must be called outside of the jslc loop. CLM is passed full
! arrays of variables and returns likewise. This is done to avoid
! altering the functionality of CLM as much as possible.
!=======================================================================
!=======================================================================
!  si version  - water fluxes are generally calculated in kg/m**2/s.
!  1000 kg/m**2/s = 1 m/s  - converted to energy units for display:
!                            1 kg/m**2/s = 2.5 x 10.e6  watts/m**2.
!  note also  1 kg/m**2/s = 1 mm/m**2/s so fluxes can be thought of
!                            in mm/s.
!=======================================================================
! include clm modules
      use atmdrvMod       , only : rcmdrv
      use clm_comp        , only : clm_run1, clm_run2
!
      implicit none

!---------------------------------------------------------------------
 
      call interfclm(1)
      call rcmdrv()
      call clm_run1(r2cdoalb,r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
      call clm_run2(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
      call interfclm(2)
      end subroutine mtrxclm
!
      subroutine initclm(instep)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Routine for initializing clm3 in regcm3
!
! written by Jason Bell 6/2004 for clm2. Updated by J.Bell 4/2005.
! Modified by: Ahmed Tawfik 8/2009
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ==========================================
! Possible atmospheric input fields to clm:
! ==========================================
! Name     Description                              Required/Optional
! -----------------------------------------------------------------------------
! TBOT     temperature (K)                          Required
! WIND     wind:sqrt(u**2+v**2) (m/s)               Required
! QBOT     specific humidity (kg/kg)                Required
! Tdew     dewpoint temperature (K)                 Alternative to Q
! RH       relative humidity (percent)              Alternative to Q
! ZBOT     reference height (m)                     optional
! PSRF     surface pressure (Pa)                    optional
! FSDS     total incident solar radiation (W/m**2)  Required
! FSDSdir  direct incident solar radiation (W/m**2) optional (replaces FSDS)
! FSDSdif  diffuse incident solar rad (W/m**2)      optional (replaces FSDS)
!clm2 start
! FSDSsdir  direct incident solar radiation (W/m**2) optional (replaces FSDSdir)
! FSDSsdif  diffuse incident solar rad (W/m**2)      optional (replaces FSDSdif)
! FSDSldir  direct incident solar radiation (W/m**2) optional (replaces FSDSdir)
! FSDSldif  diffuse incident solar rad (W/m**2)      optional (replaces FSDSdif)
!clm2 end
! FLDS     incident longwave radiation (W/m**2)     optional
! PRECTmms total precipitation (mm H2O / sec)       Required
! PRECCmms convective precipitation (mm H2O / sec)  optional (replaces PRECT)
! PRECLmms large-scale precipitation (mm H2O / sec) optional (replaces PRECT)
 
!
      use initializeMod
      use shr_orb_mod
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use clm_varpar,    only : lsmlon , lsmlat
      use clm_varsur,    only : landmask , landfrac , satbrt_clm
      use clm_varsur,    only : r2cimask , init_tgb , r2coutfrq
      use clm_varsur,    only : clm2bats_veg , ht_rcm
      use clm_varsur,    only : clm_fracveg
      use clm_varsur,    only : slmo
      use atmdrvMod
      use spmdMod,       only : masterproc, iam, spmd_init
      use program_offMod
      use clm_comp 
      use clmtype
      use perf_mod
      use mpi
! 
      implicit none
!
! Dummy arguments
!
      integer :: instep
      intent (out) instep
!
! Local variables
!
      integer :: ci , cj , i , ii , j , je , jj , js , n , ierr
      real(8) , dimension(jxp,iy) :: r2cflwd , r2cpsb , r2cqb ,         &
                & r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
                & r2csolsd , r2ctb , r2cuxb , r2cvxb , r2cxlat ,        &
                & r2cxlatd , r2cxlon , r2cxlond , r2czga
      real(r8) , dimension(jxp*iy) :: work_in
      real(r8) , dimension(jx*iy) :: work_out
!
!----------------------------------------------------------------------
!     About the dimension ordering:
!     regcm: ix=lat,jx=lon, arrays are lat by lon
!     clm: i=lon, j=lat, arrays are lon by lat
!----------------------------------------------------------------------
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Initialize run control variables for clm
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( ifrest ) then         !abt added all these indices
        r2cnsrest = 1
      else
        r2cnsrest = 0
      end if
!     land surface timestep
      r2cdtime = abatm
!     start date and time
      r2cstart_ymd = idate1/100
      r2cstart_tod = mod(idate1,100)
!     stop date and time
      r2cstop_ymd = idate2/100
      r2cstop_tod = mod(idate2,100)
!     calendar type (GREGORIAN not available in regcm)
      r2cclndr = 'NO_LEAP'
!     don't write to NCAR Mass Store
      r2cmss_irt = 0
!     clm output frequency
      r2coutfrq = clmfrq
!     radiation calculation frequency
!     regcm: radfrq is in minutes
!     clm: irad is (+) iterations or (-) hours
!     clm: hours gets converted to seconds then divided by dtime
      r2cirad = (radfrq*60/r2cdtime)
!     write output
      if ( ifbat ) then
        r2cwrtdia = .true.
      else
        r2cwrtdia = .false.
      end if
!     Set grid spacing resolution
      r2cdx = dx
!     Set gridcell area
      r2carea = (dx/1000)*(dx/1000)
!     Set landmask method
      r2cimask = imask
!     Set elevation and BATS landuse type (abt added)
      if ( .not.allocated(ht_rcm) ) allocate(ht_rcm(iy,jx))
      if ( .not.allocated(satbrt_clm) ) allocate(satbrt_clm(iy,jx))
      if ( .not.allocated(init_tgb) ) allocate(init_tgb(iy,jx))
      if ( .not.allocated(clm2bats_veg) ) allocate(clm2bats_veg(jx,iy))
      if ( .not.allocated(clm_fracveg) ) allocate(clm_fracveg(iy,jx))
      if ( myid==0 ) then
        do j = 1 , jx
          do i = 1 , iy
            ht_rcm(i,j) = ht_io(i,j)
            satbrt_clm(i,j) = satbrt_io(i,j)
            init_tgb(i,j) = ts0_io(i,j)
            clm_fracveg(i,j) = 0
          end do
        end do
      end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     End of clm run control variable initialization
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Initialize clm input variables:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Assign regcm values to the passed variables
!     Flip lat and lon for clm
!     regcm writes uneven # of j values to arrays. for now fix by
!     copying neighboring values
 
!     Some vars have not been init'd yet
      if ( .not. ifrest ) then
!       Rainfall
        pptc(:,:) = 0.0
        pptnc(:,:) = 0.0
!       Radiation
        sols2d(:,:) = 0.0
        soll2d(:,:) = 0.0
        solsd2d(:,:) = 0.0
        solld2d(:,:) = 0.0
        flwd2d(:,:) = 0.0
!       Albedo
!       Set initial albedos to clm dry soil values for mid-colored soils
        aldirs2d(:,:) = 0.16
        aldirl2d(:,:) = 0.32
        aldifs2d(:,:) = 0.16
        aldifl2d(:,:) = 0.32
      end if
 
      do j = 1 , jxp
 
!       clm3 currently works on all iy,jx instead of 2-ilx and 2-jlx so
!       copy neighboring values for now

        do i = 1 , iy
 
!         10/05 treat all variables identically. Copy i=2 to i=1 and
!         i=ilx to i=iy. For myid = 0, copy j=2 to j=1. For myid =
!         nproc-1, copy j=jendx to j=jxp.

#ifdef BAND
          cj = j
#else
          if ( myid==0 .and. j==1 ) then
            cj = 2
          else if ( myid==(nproc-1) .and. j==jxp ) then
            cj = jxp - 1
          else
            cj = j
          end if
#endif
          if ( i==1 ) then
            ci = 2
          else if ( i==iy ) then
            ci = iy - 1
          else
            ci = i
          end if
 
!         xlat,xlon in degrees
          r2cxlatd(j,i) = mddom%xlat(ci,cj)
          r2cxlond(j,i) = mddom%xlong(ci,cj)
!         xlat,xlon in radians
          r2cxlat(j,i) = mddom%xlat(ci,cj)*degrad
          r2cxlon(j,i) = mddom%xlong(ci,cj)*degrad
 
          if ( .not.ifrest ) then
!           T(K) at bottom layer
            r2ctb(j,i) = tb3d(ci,kz,cj)
!           Specific Humidity
            r2cqb(j,i) = qvb3d(ci,kz,cj)/(1.+qvb3d(ci,kz,cj))
!           Reference Height (m)
!           abt               r2czga(j,i) = za3d(ci,kz,cj)
            r2czga(j,i) = za(ci,kz,cj)
!           Surface winds
            r2cuxb(j,i) = ubx3d(ci,kz,cj)
            r2cvxb(j,i) = vbx3d(ci,kz,cj)
!           Surface Pressure in Pa from hPa
            r2cpsb(j,i) = (sps2%ps(ci,cj)+ptop)*1000.
!
            r2crnc(j,i) = pptc(ci,cj)
            r2crnnc(j,i) = pptnc(ci,cj)
            r2csols(j,i) = sols2d(ci,cj)
            r2csoll(j,i) = soll2d(ci,cj)
            r2csolsd(j,i) = solsd2d(ci,cj)
            r2csolld(j,i) = solld2d(ci,cj)
            r2cflwd(j,i) = flwd2d(ci,cj)
 
          end if
 
        end do
      end do
 
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c    1. Copy 2d (jxp,iy) arrays to 1d work_in (jx*iy) array.
!c    2. Gather jxp values of each nproc work_in array and fill
!c    work_out(jx*iy) array.
!c    3. Copy 1d work_out array to 2d (jx,iy) array for passing
!c    to clm.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( .not.ifrest ) then
!       TGB
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2ctb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2ctb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       QB
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2cqb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2cqb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       ZGA
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2czga(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2czga_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       PSB
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2cpsb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2cpsb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       UXB
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2cuxb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2cuxb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       VXB
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2cvxb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2cvxb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       RNC
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2crnc(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2crnc_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       RNNC
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2crnnc(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2crnnc_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLS
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2csols(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2csols_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLL
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2csoll(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2csoll_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLSD
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2csolsd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2csolsd_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLLD
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2csolld(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2csolld_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       LONGWAVE RAD.
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            work_in(ii) = r2cflwd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,&
                         & jxp*iy,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , iy
            r2cflwd_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
 
      end if  ! if not restart
 
 
!     XLAT in radians
      ii = 1
      do j = 1 , jxp
        do i = 1 , iy
          work_in(ii) = r2cxlat(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,  &
                       & jxp*iy,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , iy
          r2cxlat_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLON in radians
      ii = 1
      do j = 1 , jxp
        do i = 1 , iy
          work_in(ii) = r2cxlon(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,  &
                       & jxp*iy,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , iy
          r2cxlon_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLAT in degrees
      ii = 1
      do j = 1 , jxp
        do i = 1 , iy
          work_in(ii) = r2cxlatd(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,  &
                       & jxp*iy,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , iy
          r2cxlatd_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLON in degrees
      ii = 1
      do j = 1 , jxp
        do i = 1 , iy
          work_in(ii) = r2cxlond(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*iy,mpi_double_precision,work_out,  &
                       & jxp*iy,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , iy
          r2cxlond_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
 
!     Set grid edges
      r2cedgen = r2cxlatd_all(1,1)
      r2cedges = r2cxlatd_all(1,1)
      r2cedgee = r2cxlond_all(1,1)
      r2cedgew = r2cxlond_all(1,1)
 
!     clm vars are (lon x lat), regcm (lat x lon)
      do j = 1 , jx
        do i = 1 , iy
          r2cedgen = dmax1(r2cxlatd_all(j,i),r2cedgen)
          r2cedges = dmin1(r2cxlatd_all(j,i),r2cedges)
          r2cedgee = dmax1(r2cxlond_all(j,i),r2cedgee)
          r2cedgew = dmin1(r2cxlond_all(j,i),r2cedgew)
        end do
      end do
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     End Initialization of atm variables:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Get orbital parameters for use in initialize_
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      r2cnstep = ktau
 
      lsmlon = jx   !abt changed clm_varpar_init also
      lsmlat = iy
 
!!!!! Program_off initializes MPI,timing,and surface variables
 
      call program_off(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
      call clm_init0()
      call clm_init1()
      call clm_init2()
 
! -----------------------------------------------------------------
! Initialize "external" atmospheric forcing for CLM
! -----------------------------------------------------------------
 
! Read atmospheric forcing dataset one time to obtain the longitudes
! and latitudes of the atmospheric dataset, as well as the edges. When
! coupled to atm model, these are input variables. If no
! atmospheric data files are provided, model uses dummy atmospheric
! forcing and sets atmospheric grid to land grid.
 
      if ( masterproc ) write (6,*)                                     &
                                  &'Attempting to make atmospheric grid'
      call rcmdrv_init()
      if ( masterproc ) write (6,*) 'Successfully make atmospheric grid'
 
!jlb  12/05: orbital params calc'd in solar1_clm called from init.f and
!     tend.f Initializes grid and surface variables (e.g.veg type, soil
!     text....)
      call mpi_bcast(landmask,size(landmask),mpi_integer,0,             &
                   & mpi_comm_world,ierr)
 
      call mpi_bcast(landfrac,size(landfrac),mpi_real,0,mpi_comm_world, &
                   & ierr)
 
!     Initialize radiation and atmosphere variables
      if ( .not.ifrest ) then
        instep = ktau
        call rcmdrv()
      end if !end ifrest test
 
!     Initialize ocld2d now that clm has determined the land sea mask
!     Initialize accumulation variables at zero
 
      js = (jxp*myid) + 1
      je = jxp*(myid+1)
      jj = 0
      do j = js , je
        jj = jj + 1
        do i = 1 , iy - 1
 
          if ( .not.ifrest ) then
            do n = 1 , nnsg
              ocld2d(n,i,jj) = dble(landmask(j,i))
              tgb2d(n,i,jj) = sts2%tg(i,jj)
              taf2d(n,i,jj) = sts2%tg(i,jj)
              tlef2d(n,i,jj) = sts2%tg(i,jj)
              dew2d(n,i,jj) = 0.
              sag2d(n,i,jj) = 0.
              scv2d(n,i,jj) = dmax1(snowc(n,i,jj),0.D0)
              sice2d(n,i,jj) = 0.
              fsw2d(i,jj) = 0.
              flw2d(i,jj) = 0.
              sabv2d(i,jj) = 0.
              sol2d(i,jj) = 0.
              gwet2d(n,i,jj) = 0.5
              fswa2d(i,jj) = 0.
              flwa2d(i,jj) = 0.
              sena2d(n,i,jj) = 0.
              evpa2d(n,i,jj) = 0.
              prca2d(i,jj) = 0.
              prnca2d(i,jj) = 0.
              rnos2d(n,i,jj) = 0.
              rno2d(n,i,jj) = 0.
              svga2d(i,jj) = 0.
              sina2d(i,jj) = 0.
              ircp2d(n,i,jj) = 0.
            end do
          end if !end ifrest test
 
          if ( landfrac(j,i)/=1. .and. landfrac(j,i)/=0. )              &
             &    landmask(j,i) = 3.

          ! Set some clm land surface/vegetation variables to the ones
          ! used in RegCM.  Make sure all are consistent  

          mddom%satbrt(i,jj) = clm2bats_veg(j,i)
          if ( clm2bats_veg(j,i).eq.0 ) mddom%satbrt(i,jj) = 15
          do n = 1 , nnsg
            satbrt1(n,i,jj) = clm2bats_veg(j,i)
            if ( clm2bats_veg(j,i).eq.0 ) satbrt1(n,i,jj) = 15
          end do
          if ( mddom%satbrt(i,jj).gt.13.9 .and. &
               mddom%satbrt(i,jj).lt.15.1 ) then
            veg2d(i,jj)  = 0
            do n = 1 , nnsg
              veg2d1(n,i,jj)  = 0
            end do
          else
            veg2d(i,jj) = mddom%satbrt(i,jj)
            do n = 1 , nnsg
              veg2d1(n,i,jj)  = mddom%satbrt(i,jj)
            end do
          end if
!          svegfrac2d(i,jj) = clm_fracveg(i,j) 
          do n = 1 , nnsg
            if ( veg2d(i,jj).eq.0 .and. ocld2d(n,i,jj).eq.1 ) then
              veg2d(i,jj)     =  2
              veg2d1(n,i,jj)  =  2
              satbrt1(n,i,jj) =  2
              mddom%satbrt(i,jj)    =  2
            end if
          end do
          if ( landfrac(j,i).ne.1 .and. landfrac(j,i).ne.0 ) then
            landmask(j,i) = 3
          end if
        end do
      end do
 
!     deallocate some variables used in CLM initialization only
      if ( allocated(ht_rcm) ) deallocate(ht_rcm)
      if ( allocated(satbrt_clm) ) deallocate(satbrt_clm)
      if ( allocated(init_tgb) ) deallocate(init_tgb)
      if ( allocated(clm2bats_veg) ) deallocate(clm2bats_veg)
      if ( allocated(clm_fracveg) ) deallocate(clm_fracveg)
 
      end subroutine initclm
!
      subroutine albedoclm(j,iemiss)
 
      use clm_varsur,     only : landfrac
      implicit none
!
! Dummy arguments
!
      integer :: iemiss , j
      intent (in) iemiss , j
!
! Local variables
!
      real(8) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               & albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               & conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               & dralbs , fsol1 , fsol2 , sfac , sical0 , sical1 , sl , &
               & sl2 , sli , snal0 , snal1 , tdiff , tdiffs , wet , x
      real(8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 & aldifs_s , aldirl_s , aldirs_s
      real(8) :: fseas
      integer :: kolour , n , i, jj
!
!     Albedo calculates fragmented albedos (direct and diffuse) in
!     wavelength regions split at 0.7um.
!
!     CM hands albedos to radiation package which computes
!     fsw1d(i) = net solar absorbed over full grid square
!     sabveg(i) = vegetation absorbed (full solar spectrum)
!     solis(i) = shortwave  solar incident
!
!     Here these are calculated at the end of albedo - they use only
!     direct albedos for now
!
!     in both versions :  lftemp uses sabveg
!     tgrund uses sabveg & fsw1d(i) to get
!     ground absorbed solar
!     photosynthesis uses solis - see subrouts
!     stomat and co2 (carbon)
!
!     For sea, sea-ice veg albedos are not set
!     these albedos are not treated as arrays here
!
!     (depuv/10.0)= the ratio of upper soil layer to total
!     root depth; used to compute "wet" for soil albedo
!
!     =================================================================
!     1. set initial parameters
!     =================================================================
      fseas(x) = dmax1(0.D0,1.D0-0.0016D0*dmax1(298.D0-x,0.D0)**2)
!
!     1.1 constants
!
!     Solar flux partitioned at wavelength of 0.7micr
      fsol1 = 0.5D0
      fsol2 = 0.5D0
!     Short and long wave albedo for new snow
      snal0 = 0.95D0
      snal1 = 0.65D0
!     Short and long wave albedo for sea ice
      sical0 = 0.6D0
      sical1 = 0.4D0
 
!
!     In depth, wt is frac of grid square covered by snow;
!     depends on average snow depth, vegetation, etc.
!
      call depth
 
!     1.2 set pointers
!     ***************************************************
!     *    set n"x"k params here  in ccm but not needed *
!     ***************************************************
!
!     1.3  set default vegetation and albedo
!     do loop 50 in ccm not used here )
 
      do i = 2 , iym1
        czen(i) = dmax1(coszrs(i),0.D0)
        czeta = czen(i)
        do n = 1 , nnsg
          albgs = 0.0D0
          albgl = 0.0D0
          albgsd = 0.0D0
          albgld = 0.0D0
          albs = 0.0D0
          albl = 0.0D0
          albsd = 0.0D0
          albld = 0.0D0
 
          albvs_s(n) = 0.0D0
          albvl_s(n) = 0.0D0
          ldoc1d(n,i) = ocld2d(n,i,j)
 
!================================================================
!         2.   get albedo over land
!================================================================
!         can't use pointer "nalbk" here because not set - use nldock
!         instead tgb1d(i) used instead of tbelow
!
          if (iseaice == 1) then
            if ( ldoc1d(n,i).gt.1.5 ) then
              tdiffs=ts1d(n,i)-tzero
              tdiff=dmax1(tdiffs,0.d0)
              tdiffs=dmin1(tdiff,20.d0)
              albgl=sical1-1.1e-2*tdiffs
              albgs=sical0-2.45e-2*tdiffs
              albg=fsol1*albgs+fsol2*albgl
              albgsd=albgs
              albgld=albgl
            end if
          else if ( ldoc1d(n,i).gt.0.1D0 .and. sice1d(n,i).eq.0.D0 ) then
            sfac = 1.D0 - fseas(tgb1d(n,i))
!           **********  ccm tests here on land mask for veg and soils
!c          data *** reduces albedo at low temps !!!!!should respond to
!c          moisture too the following card inactivated (commented out)
!c          (pat, 27 oct 86)
!           veg1d(i)=vegc(lveg(i))-seasf(lveg(i))*sfac
            albs = albvgs(lveg(n,i))
            albl = albvgl(lveg(n,i))
 
!----------------------------------------------------------------------
            if ( (lveg(n,i).lt.12) .or. (lveg(n,i).gt.15) ) then
 
!             2.1  bare soil albedos
!             (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,i))
              wet = ssw1d(n,i)/depuv(lveg(n,i))
              alwet = dmax1((11.D0-40.D0*wet),0.D0)*0.01D0
              alwet = dmin1(alwet,solour(kolour))
              albg = solour(kolour) + alwet
!             if((lveg(n,i).eq.8)) albg=0.40      !Laura, cambiato il
!             DESERTO
              albgs = albg
              albgl = 2.D0*albg
!             **********            higher nir albedos
!             **********              set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl
 
!             Dec. 15   albzn=0.85+1./(1.+10.*czen(i))
!             Dec. 12, 2008
              albzn = 1.0D0
!             Dec. 15, 2008
 
!             **********            leafless hardwood canopy: no or
!             inverse zen dep
              if ( lveg(n,i).eq.5 .and. sfac.lt.0.1 ) albzn = 1.
!             **********            multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
 
!             **********            albedo over vegetation after zenith
!             angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
 
            else if ( lveg(n,i).eq.12 ) then
 
!             2.2   permanent ice sheet
              albgs = 0.8D0
              albgsd = 0.8D0
              albgl = 0.55D0
              albgld = 0.55D0
            else
 
!             2.3  inland water, swamps, rice paddies etc.
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgsd = albg
              albgl = albg
              albgld = albg
            end if
 
          else if ( sice1d(n,i).gt.0.D0 ) then
!====================================================================
!           3.  get albedo over sea ice
!====================================================================
!           **********          albedo depends on wave-length and ts.
!           the ts **********          dependence accounts for melt
!           water puddles.
            tdiffs = ts1d(n,i) - tzero
            tdiff = dmax1(tdiffs,0.D0)
            tdiffs = dmin1(tdiff,20.D0)
            albgl = sical1 - 1.1E-2*tdiffs
            albgs = sical0 - 2.45E-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          else
          end if
! ===================================================================
!         4.  correct for snow cover
! ===================================================================
          if ( scv1d(n,i).gt.0.0D0 ) then
!           **********            snow albedo depends on  snow-age,
!           zenith angle, **********            and thickness of snow
 
!           **********            zenith angle set in zenitm
!           **********            snow albedoes for visible and ir
!           solar rad **********            visible albedo depends on
!           snow age **********            age gives reduction of
!           visible rad snow albedo **********              due to age
            cons = 0.2D0
            conn = 0.5D0
            age = (1.D0-1.D0/(1.D0+sag1d(n,i)))
!           **********            sl helps control albedo zenith
!           dependence
            sl = 2.0D0
            sli = 1.D0/sl
            sl2 = 2.D0*sl
!           **********            snal0= new snow albedo for vis rad,
!           sol zen le 6 **********            snal1= new snow albedo
!           for long-wave rad
            dfalbs = snal0*(1.D0-cons*age)
!           **********            czf corrects albedo of new snow for
!           solar zenith
            cf1 = ((1.D0+sli)/(1.D0+sl2*czen(i))-sli)
            cff = dmax1(cf1,0.D0)
            czf = 0.4D0*cff*(1.D0-dfalbs)
            dralbs = dfalbs + czf
            dfalbl = snal1*(1.D0-conn*age)
            czf = 0.4D0*cff*(1.D0-dfalbl)
            dralbl = dfalbl + czf
 
            if ( veg1d(n,i).gt.0.001D0 ) then
!             **********            effective albedo over vegetation
!             with snow
              albl = (1.D0-wt(n,i))*albl + dralbl*wt(n,i)
              albld = (1.D0-wt(n,i))*albld + dfalbl*wt(n,i)
              albs = (1.D0-wt(n,i))*albs + dralbs*wt(n,i)
              albsd = (1.D0-wt(n,i))*albsd + dfalbs*wt(n,i)
            end if
 
!----------------------------------------------------------------------
!           4.1  compute albedo for snow on bare ground
!----------------------------------------------------------------------
            albgs = (1.D0-scvk(n,i))*albgs + dralbs*scvk(n,i)
            albgl = (1.D0-scvk(n,i))*albgl + dralbl*scvk(n,i)
            albgsd = (1.D0-scvk(n,i))*albgsd + dfalbs*scvk(n,i)
            albgld = (1.D0-scvk(n,i))*albgld + dfalbl*scvk(n,i)
          end if
 
!=====================================================================
!         5.  albedo over open ocean
!=====================================================================
          if ( ldoc1d(n,i).eq.0.D0 ) then
!           *********   ocean albedo depends on zenith angle
            if ( czeta.ge.0.0D0 ) then
!             **********   albedo independent of wavelength
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            end if
          end if
 
!
!         ***************not part of albedo in the ccm ****************
!
          aldirs_s(n) = (1.D0-veg1d(n,i))*albgs + veg1d(n,i)*albs
          aldirl_s(n) = (1.D0-veg1d(n,i))*albgl + veg1d(n,i)*albl
          aldifs_s(n) = (1.D0-veg1d(n,i))*albgsd + veg1d(n,i)*albsd
          aldifl_s(n) = (1.D0-veg1d(n,i))*albgld + veg1d(n,i)*albld
        end do
        albvs(i) = albvs_s(1)
        albvl(i) = albvl_s(1)
        aldirs(i) = aldirs_s(1)
        aldirl(i) = aldirl_s(1)
        aldifs(i) = aldifs_s(1)
        aldifl(i) = aldifl_s(1)
        if ( iemiss.eq.1 ) emiss1d(i) = emiss2d(1,i,j)
        do n = 2 , nnsg
          albvs(i) = albvs(i) + albvs_s(n)
          albvl(i) = albvl(i) + albvl_s(n)
          aldirs(i) = aldirs(i) + aldirs_s(n)
          aldirl(i) = aldirl(i) + aldirl_s(n)
          aldifs(i) = aldifs(i) + aldifs_s(n)
          aldifl(i) = aldifl(i) + aldifl_s(n)
          if ( iemiss.eq.1 ) emiss1d(i) = emiss1d(i) + emiss2d(n,i,j)
        end do
        albvs(i) = albvs(i)/dble(nnsg)
        albvl(i) = albvl(i)/dble(nnsg)
        aldirs(i) = aldirs(i)/dble(nnsg)
        aldirl(i) = aldirl(i)/dble(nnsg)
        aldifs(i) = aldifs(i)/dble(nnsg)
        aldifl(i) = aldifl(i)/dble(nnsg)
        if ( iemiss.eq.1 ) emiss1d(i) = emiss1d(i)/dble(nnsg)
 
!       ******   fsw1d(i),sabveg(i),solis(i) computed in colmod (radcsw)

!       ****** Section Below added for albedo to be corrected by CLM
!       ****** calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
!       ****** to be equal to 1 which can cause a FPE.  To avoid this
!       ****** use albedo calculated with BATS method when albedo=1
#ifdef MPP1
        jj = j+(jxp*myid)
#else
        jj = j
#endif
        if(ocld2d(1,i,j)>=1 .and. aldirs2d(i,j)/=1.) then
           aldirs(i) = aldirs2d(i,j)*landfrac(jj,i) +                   &
                       aldirs(i)*(1-landfrac(jj,i))
           aldirl(i) = aldirl2d(i,j)*landfrac(jj,i) +                   &
                       aldirl(i)*(1-landfrac(jj,i))
           aldifs(i) = aldifs2d(i,j)*landfrac(jj,i) +                   &
                       aldifs(i)*(1-landfrac(jj,i))
           aldifl(i) = aldifl2d(i,j)*landfrac(jj,i) +                   &
                       aldifl(i)*(1-landfrac(jj,i))
           albvs(i)  = aldirs2d(i,j)*landfrac(jj,i) +                   &
                       albvs(i) *(1-landfrac(jj,i))
           albvl(i)  = aldirl2d(i,j)*landfrac(jj,i) +                   &
                       albvl(i) *(1-landfrac(jj,i)) 
        end if

      end do   ! end of i loop
  
      end subroutine albedoclm
!
      subroutine interfclm(ivers)

!=======================================================================
!l  built for clm version 3.0
!=======================================================================
! ivers = 1 : regcm -> clm
! ivers = 2 : clm -> regcm
!
      use clm_varsur,    only : landmask, landfrac
      use clmtype
      use clm_varsur,    only : c2r_allout,omap_i,omap_j
      use mpi
      implicit none
!
! Dummy arguments
!
      integer :: ivers
      intent (in) ivers
!
! Local variables
!
      real(8) :: mmpd , wpm2
      integer :: ci , cj , counter , i , ii , iii , j , jj , kk ,       &
                & locid , n , nn1 , nnn , nout
      real(8) , dimension(jxp,iy) :: r2cflwd , r2cpsb , r2cqb ,    &
                & r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
                & r2csolsd , r2ctb , r2cuxb , r2cvxb , r2czga
      real(4) :: real_4
      real(8) , dimension(jxp*iy*13) :: workin
      real(8) , dimension(jx*iy*13) :: workout
      integer :: ierr
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     About the dimension ordering:
!     regcm: ix=lat,jx=lon, arrays are lat by lon
!     clm: i=lon, j=lat, arrays are lon by lat
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      if ( ivers==1 ) then
 
        locid = 1
!       clm3 currently works on all iy,jx instead of 2-ilx and 2-jlx so
!       copy neighboring values for now
        do j = 1 , jxp
          do i = 1 , iy
 
!jlb        10/05 treat all variables identically. Copy i=2 to i=1 and
!           i=ilx to i=iy. For myid = 0, copy j=2 to j=1. For myid =
!           nproc-1, copy j=jendx to j=jxp.
#ifdef BAND
            cj = j
#else
            if ( myid==0 .and. j==1 ) then
              cj = 2
            else if ( myid==(nproc-1) .and. j==jxp ) then
              cj = jxp - 1
            else
              cj = j
            end if
#endif
            if ( i==1 ) then
              ci = 2
            else if ( i==iy ) then
              ci = iy - 1
            else
              ci = i
            end if
 
!           T(K) at bottom layer
            r2ctb(j,i) = tb3d(ci,kz,cj)
!           Specific Humidity ?
            r2cqb(j,i) = qvb3d(ci,kz,cj)/(1+qvb3d(ci,kz,cj))
!           Reference Height (m)
            r2czga(j,i) = za(ci,kz,cj)
!           Surface winds
            r2cuxb(j,i) = ubx3d(ci,kz,cj)
            r2cvxb(j,i) = vbx3d(ci,kz,cj)
!           Surface Pressure in Pa from cbar
            r2cpsb(j,i) = (sps2%ps(ci,cj)+ptop)*1000.
!           Rainfall
            r2crnc(j,i) = pptc(ci,cj)
            r2crnnc(j,i) = pptnc(ci,cj)
!           Incident Solar Radiation
            r2csols(j,i) = sols2d(ci,cj)
            r2csoll(j,i) = soll2d(ci,cj)
            r2csolsd(j,i) = solsd2d(ci,cj)
            r2csolld(j,i) = solld2d(ci,cj)
            r2cflwd(j,i) = flwd2d(ci,cj)
 
          end do
        end do
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      cc
!c      1. Copy 2d (jxp,iy) arrays to 1d work_in (jx*iy) array.      cc
!c      2. Gather jxp values of each nproc work_in array and fill     cc
!c      work_out(jx*iy) array.                                   cc
!c      3. Copy 1d work_out array to 2d (jx,iy) array for passing    cc
!c      to clm.                                                   cc
!c      abt updated below 1/09                                        cc
!c      UPDATE:  copy all r2c vars to one large array; this allows    cc
!c      for one MPI_ALLGATHER call instead of several        cc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
        ii = 1
        do j = 1 , jxp
          do i = 1 , iy
            workin(ii) = r2ctb(j,i)
            workin(ii+(jxp*iy)) = r2cqb(j,i)
            workin(ii+(2*jxp*iy)) = r2czga(j,i)
            workin(ii+(3*jxp*iy)) = r2cpsb(j,i)
            workin(ii+(4*jxp*iy)) = r2cuxb(j,i)
            workin(ii+(5*jxp*iy)) = r2cvxb(j,i)
            workin(ii+(6*jxp*iy)) = r2crnc(j,i)
            workin(ii+(7*jxp*iy)) = r2crnnc(j,i)
            workin(ii+(8*jxp*iy)) = r2csols(j,i)
            workin(ii+(9*jxp*iy)) = r2csoll(j,i)
            workin(ii+(10*jxp*iy)) = r2csolsd(j,i)
            workin(ii+(11*jxp*iy)) = r2csolld(j,i)
            workin(ii+(12*jxp*iy)) = r2cflwd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(workin,13*jxp*iy,mpi_double_precision,       &
                         & workout,13*jxp*iy,mpi_double_precision,      &
                         & mpi_comm_world,ierr)
 
        ii = 1
        kk = 1
        counter = 1
        do j = 1 , jx
          do i = 1 , iy
            r2ctb_all(j,i) = workout(ii)
            r2cqb_all(j,i) = workout(ii+(jxp*iy))
            r2czga_all(j,i) = workout(ii+(2*jxp*iy))
            r2cpsb_all(j,i) = workout(ii+(3*jxp*iy))
            r2cuxb_all(j,i) = workout(ii+(4*jxp*iy))
            r2cvxb_all(j,i) = workout(ii+(5*jxp*iy))
            r2crnc_all(j,i) = workout(ii+(6*jxp*iy))
            r2crnnc_all(j,i) = workout(ii+(7*jxp*iy))
            r2csols_all(j,i) = workout(ii+(8*jxp*iy))
            r2csoll_all(j,i) = workout(ii+(9*jxp*iy))
            r2csolsd_all(j,i) = workout(ii+(10*jxp*iy))
            r2csolld_all(j,i) = workout(ii+(11*jxp*iy))
            r2cflwd_all(j,i) = workout(ii+(12*jxp*iy))
            ii = ii + 1
            counter = counter + 1
          end do
          if ( counter>jxp*iy ) then
            kk = kk + 1
            counter = 1
            ii = jxp*(kk-1)*13*iy + 1
          end if
        end do

      else if ( ivers==2 ) then ! end of ivers = 1

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c 1. Copy the parts of 2d (jx,iy) clm arrays that contain
!c     data to 1d work_in (jx*iy) array.
!c 2. Gather jxp values of each nproc work_in array and fill
!c     work_out (jx*iy) array.
!c 3. Copy 1d work_out array to 2d (jx,iy) array for passing
!c     to nproc 2d (jxp,iy) arrays.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!JLB    3/06: NEED TO CREATE MASK/FILTER FOR EACH CPU. CLM DOES NOT
!       ASSIGN GRID CELLS IN QUITE THE ROUND ROBIN FASHION INDICATED.
!       TGB
!       3/06 -> New routine utilizing c2rprocmap in case of non-round
!       robin cpu assignment in clm.
!abt    updated below
!       1/09 -> update use MPI_ALLGATHER call located in clm_atmlnd.F90
!       module This module saves each land surface variable regcm needs
!       into a large array (c2r_all).  gathers that array (c2r_allout)
!       and the corresponding grid point on the i/j grid to omap_i/j
!       Finally, below c2r_allout vars are placed into the proper place
!       on the grid for each c2r variable
!       NOTE: CLM refers to i = lon while REGCM refers to j = lon
!
        iii = 0
        jj = 1
        if ( aertyp/='AER00D0' ) then
          nout = 22
        else
          nout = 20
        end if
        do nn1 = 1 , nproc
          do ii = 1 , c2rngc(nn1)
            kk = c2rngc(nn1)
            j = omap_i(jj)
            i = omap_j(jj)
 
            c2rtgb(j,i) = c2r_allout(ii+iii)
            c2rsnowc(j,i) = c2r_allout(ii+kk+iii)
            c2rsenht(j,i) = c2r_allout(ii+(2*kk)+iii)
            c2rlatht(j,i) = c2r_allout(ii+(3*kk)+iii)
            c2ruvdrag(j,i) = c2r_allout(ii+(4*kk)+iii)
            c2ralbdirs(j,i) = c2r_allout(ii+(5*kk)+iii)
            c2ralbdirl(j,i) = c2r_allout(ii+(6*kk)+iii)
            c2ralbdifs(j,i) = c2r_allout(ii+(7*kk)+iii)
            c2ralbdifl(j,i) = c2r_allout(ii+(8*kk)+iii)
            c2rtgbb(j,i) = c2r_allout(ii+(9*kk)+iii)
            c2r2mt(j,i) = c2r_allout(ii+(10*kk)+iii)
            c2r2mq(j,i) = c2r_allout(ii+(11*kk)+iii)
            c2ru10(j,i) = c2r_allout(ii+(12*kk)+iii)
            c2rtlef(j,i) = c2r_allout(ii+(13*kk)+iii)
            c2rsm10cm(j,i) = c2r_allout(ii+(14*kk)+iii)
            c2rsm1m(j,i) = c2r_allout(ii+(15*kk)+iii)
            c2rsmtot(j,i) = c2r_allout(ii+(16*kk)+iii)
            c2rinfl(j,i) = c2r_allout(ii+(17*kk)+iii)
            c2rro_sur(j,i) = c2r_allout(ii+(18*kk)+iii)
            c2rro_sub(j,i) = c2r_allout(ii+(19*kk)+iii)
            if ( aertyp/='AER00D0' ) then
              c2rfracsno(j,i) = c2r_allout(ii+(20*kk)+iii)
              c2rfvegnosno(j,i) = c2r_allout(ii+(21*kk)+iii)
            end if
            jj = jj + 1
          end do
          iii = iii + c2rngc(nn1)*nout
        end do
 
        deallocate(c2r_allout)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      cc
!c      Fill nproc 2d (jxp,iy) arrays from full 2d (jx,iy) clm data. cc
!c      cc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
        if ( jyear==jyear0 .and. ktau<=1 ) then
          mmpd = 86400./dtbat
          wpm2 = 1./dtbat
        else if ( jyear==jyear0 .and. dble(ktau*dtmin)<=batfrq*60.+     &
                & 0.01 ) then
          mmpd = 24./(batfrq-dtmin/60.)
          wpm2 = 1./((batfrq-dtmin/60.)*3600.)
        else
          mmpd = 24./batfrq
          wpm2 = 1./(batfrq*3600.)
        end if
 
 
        jj = 0
        do j = jbegin , jendx
          jj = (jxp*myid) + j
 
          call interf(1 , j , kz , 2 , iym1 , nnsg)
          if ( iocnflx==2 ) call zengocndrv(j, nnsg , 2 , iym1 , kz)
 
          do i = 2 , iym1
            ci = i
            sfsta%uvdrag(i,j) = 0.0
            sfsta%hfx(i,j) = 0.0
            sfsta%qfx(i,j) = 0.0
            sts2%tg(i,j) = 0.0
            sts1%tg(i,j) = 0.0
            sfsta%tgbb(i,j) = 0.0
!chem2
            ssw2da(i,j) = 0.0
            sdeltk2d(i,j) = 0.0
            sdelqk2d(i,j) = 0.0
            sfracv2d(i,j) = 0.0
            sfracb2d(i,j) = 0.0
            sfracs2d(i,j) = 0.0
!chem2_
            if ( landmask(jj,ci)==1 ) then
              sts2%tg(i,j) = c2rtgb(jj,ci)
              sts1%tg(i,j) = c2rtgb(jj,ci)
              sfsta%hfx(i,j) = c2rsenht(jj,ci)
              sfsta%qfx(i,j) = c2rlatht(jj,ci)
              sfsta%uvdrag(i,j) = c2ruvdrag(jj,ci)
              sfsta%tgbb(i,j) = c2rtgbb(jj,ci)
 
              if ( r2cdoalb ) coszrs2d(i,j) = c2rcosz(jj,ci)
              if ( i<=iy-1 ) then
                aldirs2d(i,j) = c2ralbdirs(jj,ci)
                aldirl2d(i,j) = c2ralbdirl(jj,ci)
                aldifs2d(i,j) = c2ralbdifs(jj,ci)
                aldifl2d(i,j) = c2ralbdifl(jj,ci)
              end if
 
              do n = 1 , nnsg
                snowc(n,i,j) = c2rsnowc(jj,ci)
                tg2d(n,i,j) = c2rtgb(jj,ci)
                tgb2d(n,i,j) = c2rtgb(jj,ci)
                !supposed to be lower soil layer temp not tgrnd
                taf2d(n,i,j) = c2r2mt(jj,ci)
                tlef2d(n,i,j) = c2rtlef(jj,ci)
                swt2d(n,i,j) = c2rsmtot(jj,ci)
                srw2d(n,i,j) = c2rsm1m(jj,ci)
                ssw2d(n,i,j) = c2rsm10cm(jj,ci)
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)    !snow age
                scv2d(n,i,j) = c2rsnowc(jj,ci)
                sice2d(n,i,j) = sice1d(n,i)  ! sea ice
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*sfsta%qfx(i,j)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sfsta%hfx(i,j)
                rnos2d(n,i,j) = c2rro_sur(jj,ci)*dtbat
                rno2d(n,i,j) = (c2rro_sub(jj,ci)+c2rro_sur(jj,ci))*dtbat
 
                ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
                sfracb2d(i,j) = sfracb2d(i,j)                           &
                              & + 1 - (c2rfvegnosno(jj,ci)+             &
                              & c2rfracsno(jj,ci))
                sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
              end do
 
              !abt added for 2m humidity when landmask = 1 or 3
              q2d(i,j) = c2r2mq(jj,ci)
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
            else if ( landmask(jj,ci)==0 ) then !ocean
 
              do n = 1 , nnsg
                sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
                sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
                sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
                sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
                sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)
!chem2
                ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
                sfracb2d(i,j) = sfracb2d(i,j) + (1.-sigf(n,i))          &
                              & *(1.-scvk(n,i))
                sfracs2d(i,j) = sfracs2d(i,j) + sigf(n,i)*wt(n,i)       &
                              & + (1.-sigf(n,i))*scvk(n,i)
!               svegfrac2d(i,j)= svegfrac2d(i,j)+veg1d(n,i)
!chem2_
 
                if ( iocnflx==1 .or.                                    &
                   & (iocnflx==2 .and. ocld2d(n,i,j)>=0.5) ) then
                  sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                     &
                            & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)&
                            & *tlef1d(n,i)**4)**0.25
                else
                  sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
                end if
                ssw1d(n,i) = -1.E34
                rsw1d(n,i) = -1.E34
                tsw1d(n,i) = -1.E34
                rno1d(n,i) = -1.E34
                rnos1d(n,i) = -1.E34
                scv1d(n,i) = -1.E34
              end do
 
              sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)/float(nnsg)
              sfsta%hfx(i,j) = sfsta%hfx(i,j)/float(nnsg)
              sfsta%qfx(i,j) = sfsta%qfx(i,j)/float(nnsg)
              sts2%tg(i,j) = sts2%tg(i,j)/float(nnsg)
              sts1%tg(i,j) = sts1%tg(i,j)/float(nnsg)
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)/float(nnsg)
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
              do n = 1 , nnsg
                snowc(n,i,j) = scv1d(n,i)
!fix
                tg2d(n,i,j) = tg1d(n,i)
!fix_
                tgb2d(n,i,j) = tgb1d(n,i)
!               taf2d(n,i,j)=taf1d(n,i)
                taf2d(n,i,j) = t2m_1d(n,i)
                !note taf2d is not temp in canopy but 2m temp
                tlef2d(n,i,j) = tlef1d(n,i)
                swt2d(n,i,j) = tsw1d(n,i)
                srw2d(n,i,j) = rsw1d(n,i)
                ssw2d(n,i,j) = ssw1d(n,i)
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)
                scv2d(n,i,j) = scv1d(n,i)
                sice2d(n,i,j) = sice1d(n,i)
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
                if ( rnos2d(n,i,j)>-1.E10 .and. rnos1d(n,i)>-1.E10 )    &
                   & then
                  rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)/tau1*dtbat
                else
                  rnos2d(n,i,j) = -1.E34
                end if
                if ( rno2d(n,i,j)>-1.E10 .and. rnos1d(n,i)>-1.E10 .and. &
                   & rno1d(n,i)>-1.E10 ) then
                  rno2d(n,i,j) = rno2d(n,i,j) + (rno1d(n,i)-rnos1d(n,i))&
                               & /tau1*dtbat
                else
                  rno2d(n,i,j) = -1.E34
                end if
              end do
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
 
            else if ( landmask(jj,ci)==3 ) then
            !gridcell with some % land and ocean
 
              do n = 1 , nnsg
                sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
                sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
                sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
                sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
                sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)
!chem2
                ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
                sfracb2d(i,j) = sfracb2d(i,j)                           &
                              & + 1 - (c2rfvegnosno(jj,ci)+             &
                              & c2rfracsno(jj,ci))
                sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
 
                ssw2da(i,j) = ssw2da(i,j)*landfrac(jj,ci)               &
                            & + (1-landfrac(jj,ci))*ssw1d(n,i)
                sdeltk2d(i,j) = sdeltk2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*sigf(n,i)
                sfracb2d(i,j) = sfracb2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*(1.-sigf(n,i))    &
                              & *(1.-scvk(n,i))
                sfracs2d(i,j) = sfracs2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))                   &
                              & *(sigf(n,i)*wt(n,i)+(1.-sigf(n,i))      &
                              & *scvk(n,i))
!chem2_
 
                if ( iocnflx==1 .or.                                    &
                   & (iocnflx==2 .and. ocld2d(n,i,j)>=0.5) ) then
                  sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                     &
                            & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)&
                            & *tlef1d(n,i)**4)**0.25
                else
                  sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
                end if
!               ssw1d(n,i)=-1.e34
!               rsw1d(n,i)=-1.e34
!               tsw1d(n,i)=-1.e34
!               rno1d(n,i)=-1.e34
!               rnos1d(n,i)=-1.e34
!               scv1d(n,i)=-1.e34
              end do
 
              sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)*(1-landfrac(jj,ci)) &
                          & + c2ruvdrag(jj,ci)*landfrac(jj,ci)
              sfsta%hfx(i,j) = sfsta%hfx(i,j)*(1-landfrac(jj,ci)) +     &
                               c2rsenht(jj,ci)*landfrac(jj,ci)
              sfsta%qfx(i,j) = sfsta%qfx(i,j)*(1-landfrac(jj,ci)) +     &
                               c2rlatht(jj,ci)*landfrac(jj,ci)
              sts2%tg(i,j) = sts2%tg(i,j)*(1-landfrac(jj,ci)) +         &
                             c2rtgb(jj,ci)*landfrac(jj,ci)
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)*(1-landfrac(jj,ci)) +   &
                                c2rtgbb(jj,ci)*landfrac(jj,ci)
              sts1%tg(i,j) = sts1%tg(i,j)*(1-landfrac(jj,ci)) +         &
                             c2rtgb(jj,ci)*landfrac(jj,ci)
 
              sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)/float(nnsg)
              sfsta%hfx(i,j) = sfsta%hfx(i,j)/float(nnsg)
              sfsta%qfx(i,j) = sfsta%qfx(i,j)/float(nnsg)
              sts2%tg(i,j) = sts2%tg(i,j)/float(nnsg)
              sts1%tg(i,j) = sts1%tg(i,j)/float(nnsg)
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)/float(nnsg)
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
              do n = 1 , nnsg
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)
                scv2d(n,i,j) = scv1d(n,i)
                sice2d(n,i,j) = sice1d(n,i)
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
!abt            added below for the landfraction method
                snowc(n,i,j) = c2rsnowc(jj,ci)*landfrac(jj,ci)          &
                             & + scv1d(n,i)*(1-landfrac(jj,ci))
                tg2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci) + tg1d(n,i) &
                            & *(1-landfrac(jj,ci))
                tgb2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci)            &
                             & + tgb1d(n,i)*(1-landfrac(jj,ci))
                taf2d(n,i,j) = c2r2mt(jj,ci)*landfrac(jj,ci)            &
                             & + t2m_1d(n,i)*(1-landfrac(jj,ci))
                !note taf2d is 2m temp not temp in foilage
                tlef2d(n,i,j) = c2rtlef(jj,ci)*landfrac(jj,ci)          &
                              & + tlef1d(n,i)*(1-landfrac(jj,ci))
                swt2d(n,i,j) = c2rsmtot(jj,ci)*landfrac(jj,ci)          &
                             & + tsw1d(n,i)*(1-landfrac(jj,ci))
                srw2d(n,i,j) = c2rsm1m(jj,ci)*landfrac(jj,ci)           &
                             & + rsw1d(n,i)*(1-landfrac(jj,ci))
                ssw2d(n,i,j) = c2rsm10cm(jj,ci)*landfrac(jj,ci)         &
                             & + ssw1d(n,i)*(1-landfrac(jj,ci))
                q2d(i,j) = c2r2mq(jj,ci)*landfrac(jj,ci) + q2m_1d(n,i)  &
                         & *(1-landfrac(jj,ci))
 
 
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*sfsta%qfx(i,j)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sfsta%hfx(i,j)
                rnos2d(n,i,j) = c2rro_sur(jj,ci)*dtbat
                rno2d(n,i,j) = c2rro_sub(jj,ci)*dtbat + c2rro_sur(jj,ci)&
                             & *dtbat
!abt            above
              end do
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
 
            else
               !landmask
            end if
          end do !i loop
 
!!!!!!!!!!!!!!!! addition from new RegCM !!!!!!!!!!!!!!!!!!!
 
          do i = 2 , iym1
            ci = i
 
            u10m_o(j,i-1) = 0.0
            v10m_o(j,i-1) = 0.0
            tg_o(j,i-1) = 0.0
            t2m_o(j,i-1) = 0.0
 
            do n = 1 , nnsg
              if ( ocld2d(n,i,j)>0.5 ) then
                u10m_s(n,j,i-1) = ubx3d(i,kz,j)
                v10m_s(n,j,i-1) = vbx3d(i,kz,j)
                tg_s(n,j,i-1) = tg2d(n,i,j)
                t2m_s(n,j,i-1) = taf2d(n,i,j)
!               abt            u10m_o(j,i-1)= u10m_o(j,i-1)+ u10m1d(n,i)
!               abt            v10m_o(j,i-1)= v10m_o(j,i-1)+ v10m1d(n,i)
                u10m_o(j,i-1) = u10m_o(j,i-1) + ubx3d(i,kz,j)
                v10m_o(j,i-1) = v10m_o(j,i-1) + vbx3d(i,kz,j)
                t2m_o(j,i-1) = t2m_o(j,i-1) + taf2d(n,i,j)
                tg_o(j,i-1) = tg_o(j,i-1) + tg2d(n,i,j)
              else if ( ocld2d(n,i,j)<0.5 ) then
                tg_s(n,j,i-1) = tg1d(n,i)
                u10m_s(n,j,i-1) = u10m1d(n,i)
                v10m_s(n,j,i-1) = v10m1d(n,i)
                t2m_s(n,j,i-1) = t2m_1d(n,i)
 
                u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
                v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
                t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
                tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
              else
              end if
            end do
 
            u10m_o(j,i-1) = u10m_o(j,i-1)/float(nnsg)
            v10m_o(j,i-1) = v10m_o(j,i-1)/float(nnsg)
            t2m_o(j,i-1) = t2m_o(j,i-1)/float(nnsg)
            tg_o(j,i-1) = tg_o(j,i-1)/float(nnsg)
 
            tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
            tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
            t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
            t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
            w10x_o(j,i-1) = amax1(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**2+  &
                          & v10m_o(j,i-1)**2))
            real_4 = (sps2%ps(i,j)+ptop)*10.
            psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
 
          end do
             !i loop
 
          if ( mod(ntime+nint(dtmin*60.),kbats)==0 .or.                 &
             & (jyear==jyearr .and. ktau==ktaur) ) then
            if ( jyear==jyear0 .and. ktau<=1 ) then
              mmpd = 86400./dtbat
              wpm2 = 1./dtbat
            else if ( jyear==jyear0 .and. dble(ktau*dtmin)<=batfrq*60.+ &
                    & 0.01 ) then
              mmpd = 24./(batfrq-dtmin/60.)
              wpm2 = 1./((batfrq-dtmin/60.)*3600.)
            else
              mmpd = 24./batfrq
              wpm2 = 1./(batfrq*3600.)
            end if
 
            do i = 2 , iym1
              ci = i
 
              drag_o(j,i-1) = 0.0
              q2m_o(j,i-1) = 0.0
              evpa_o(j,i-1) = 0.0
              sena_o(j,i-1) = 0.0
              do n = 1 , nnsg
                if ( ocld2d(n,i,j)>=0.5 ) then
                  q2m_s(n,j,i-1) = q2d(i,j)
                  drag_s(n,j,i-1) = sfsta%uvdrag(i,j)
                  evpa_s(n,j,i-1) = evpa2d(n,ci,j)*mmpd
                  sena_s(n,j,i-1) = sena2d(n,ci,j)*wpm2
                  tpr_s(n,j,i-1) = (prnca2d(ci,j)+prca2d(ci,j))*mmpd
                  prcv_s(n,j,i-1) = prca2d(ci,j)*mmpd
                  ps_s(n,j,i-1) = p1d(n,i)*0.01
 
                  q2m_o(j,i-1) = q2m_o(j,i-1) + q2d(i,j)
                  drag_o(j,i-1) = drag_o(j,i-1) + sfsta%uvdrag(i,j)
                  evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,ci,j)
                  sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,ci,j)
                else if ( ocld2d(n,i,j)<=0.5 ) then
                  q2m_s(n,j,i-1) = q2m_1d(n,i)
                  drag_s(n,j,i-1) = drag1d(n,i)
                  evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
                  sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
                  tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
                  prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
                  ps_s(n,j,i-1) = p1d(n,i)*0.01
 
                  q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
                  drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
                  evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
                  sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
                else
                end if
              end do
              tpr_o(j,i-1) = (prnca2d(ci,j)+prca2d(ci,j))*mmpd
              q2m_o(j,i-1) = q2m_o(j,i-1)/float(nnsg)
              drag_o(j,i-1) = drag_o(j,i-1)/float(nnsg)
              evpa_o(j,i-1) = evpa_o(j,i-1)/float(nnsg)*mmpd
              sena_o(j,i-1) = sena_o(j,i-1)/float(nnsg)*wpm2
              flwa_o(j,i-1) = flwa2d(ci,j)*wpm2
              fswa_o(j,i-1) = fswa2d(ci,j)*wpm2
              flwd_o(j,i-1) = flwda2d(ci,j)*wpm2
              sina_o(j,i-1) = sina2d(ci,j)*wpm2
              prcv_o(j,i-1) = prca2d(ci,j)*mmpd
              ps_o(j,i-1) = (sps2%ps(i,j)+ptop)*10.
              zpbl_o(j,i-1) = sfsta%zpbl(i,j)
 
              tlef_o(j,i-1) = 0.0
              ssw_o(j,i-1) = 0.0
              rsw_o(j,i-1) = 0.0
              rnos_o(j,i-1) = 0.0
              scv_o(j,i-1) = 0.0
              nnn = 0
              do n = 1 , nnsg
!abt           if(ocld2d(n,ci,j).ge.0.5) then
                if ( ocld2d(n,ci,j)>=0.5 .and. landmask(jj,ci)/=3 ) then
                  tlef_o(j,i-1) = tlef_o(j,i-1) + c2rtlef(jj,ci)
                  ssw_o(j,i-1) = ssw_o(j,i-1) + c2rsm10cm(jj,ci)
                  rsw_o(j,i-1) = rsw_o(j,i-1) + c2rsm1m(jj,ci)
                  rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,ci,j)
                  scv_o(j,i-1) = scv_o(j,i-1) + c2rsnowc(jj,ci)
                  tlef_s(n,j,i-1) = c2rtlef(jj,ci)
                  ssw_s(n,j,i-1) = c2rsm10cm(jj,ci)
                  rsw_s(n,j,i-1) = c2rsm1m(jj,ci)
                  rnos_s(n,j,i-1) = rnos2d(n,ci,j)*mmpd
                  scv_s(n,j,i-1) = c2rsnowc(jj,ci)
                  nnn = nnn + 1
                else
                  tlef_s(n,j,i-1) = -1.E34
                  ssw_s(n,j,i-1) = -1.E34
                  rsw_s(n,j,i-1) = -1.E34
                  rnos_s(n,j,i-1) = -1.E34
                  scv_s(n,j,i-1) = -1.E34
                end if
              end do
              if ( nnn>=max0(nnsg/2,1) ) then
                tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
                ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
                rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
                rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
                scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
              else
                tlef_o(j,i-1) = -1.E34
                ssw_o(j,i-1) = -1.E34
                rsw_o(j,i-1) = -1.E34
                rnos_o(j,i-1) = -1.E34
                scv_o(j,i-1) = -1.E34
              end if
!             ******    reset accumulation arrays to zero
              do n = 1 , nnsg
                evpa2d(n,ci,j) = 0.
                rnos2d(n,ci,j) = 0.
                sena2d(n,ci,j) = 0.
              end do
              prnca2d(ci,j) = 0.
              prca2d(ci,j) = 0.
              flwa2d(ci,j) = 0.
              flwda2d(ci,j) = 0.
              fswa2d(ci,j) = 0.
              svga2d(ci,j) = 0.
              sina2d(ci,j) = 0.
 
            end do  ! end of i loop
          end if    ! end if jyear eq jyearr
        end do      ! end of j loop
      else          ! end if ivers = 2
      end if
 
      end subroutine interfclm
!
      end module mod_mtrxclm

#endif

#endif
