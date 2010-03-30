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

#ifdef MPP1
#ifdef CLM
      subroutine initclm(nstep)

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
 
      use mod_regcm_param
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
!
      use mpi
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_clm
      use mod_bats
      use mod_mppio
      use mod_constants
      use mod_main
      use mod_pbldim
      use mod_slice
! 
      implicit none
!
! Dummy arguments
!
      integer :: nstep
      intent (out) nstep
!
! Local variables
!
      integer :: ci , cj , i , ii , j , je , jj , js , n , ierr
      real(8) , dimension(jxp,ix) :: r2cflwd , r2cpsb , r2cqb ,         &
                & r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
                & r2csolsd , r2ctb , r2cuxb , r2cvxb , r2cxlat ,        &
                & r2cxlatd , r2cxlon , r2cxlond , r2czga
      real(r8) , dimension(jxp*ix) :: work_in
      real(r8) , dimension(jx*ix) :: work_out
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
      if ( .not.allocated(ht_rcm) ) allocate(ht_rcm(ix,jx))
      if ( .not.allocated(satbrt_clm) ) allocate(satbrt_clm(ix,jx))
      if ( .not.allocated(init_tgb) ) allocate(init_tgb(ix,jx))
      if ( .not.allocated(clm2bats_veg) ) allocate(clm2bats_veg(jx,ix))
      if ( .not.allocated(clm_fracveg) ) allocate(clm_fracveg(ix,jx))
      if ( myid==0 ) then
        do j = 1 , jx
          do i = 1 , ix
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
 
!       clm3 currently works on all ix,jx instead of 2-ilx and 2-jlx so
!       copy neighboring values for now

        do i = 1 , ix
 
!         10/05 treat all variables identically. Copy i=2 to i=1 and
!         i=ilx to i=ix. For myid = 0, copy j=2 to j=1. For myid =
!         nproc-1, copy j=jendx to j=jxp.

          if ( myid==0 .and. j==1 ) then
            cj = 2
          else if ( myid==(nproc-1) .and. j==jxp ) then
            cj = jxp - 1
          else
            cj = j
          end if
          if ( i==1 ) then
            ci = 2
          else if ( i==ix ) then
            ci = ix - 1
          else
            ci = i
          end if
 
!         xlat,xlon in degrees
          r2cxlatd(j,i) = xlat(ci,cj)
          r2cxlond(j,i) = xlong(ci,cj)
!         xlat,xlon in radians
          r2cxlat(j,i) = xlat(ci,cj)*degrad
          r2cxlon(j,i) = xlong(ci,cj)*degrad
 
          if ( .not.ifrest ) then
!           T(K) at bottom layer
            r2ctb(j,i) = tb3d(ci,kx,cj)
!           Specific Humidity
            r2cqb(j,i) = qvb3d(ci,kx,cj)/(1.+qvb3d(ci,kx,cj))
!           Reference Height (m)
!           abt               r2czga(j,i) = za3d(ci,kx,cj)
            r2czga(j,i) = za(ci,kx,cj)
!           Surface winds
            r2cuxb(j,i) = ubx3d(ci,kx,cj)
            r2cvxb(j,i) = vbx3d(ci,kx,cj)
!           Surface Pressure in Pa from hPa
            r2cpsb(j,i) = (psb(ci,cj)+ptop)*1000.
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
!c    1. Copy 2d (jxp,ix) arrays to 1d work_in (jx*ix) array.
!c    2. Gather jxp values of each nproc work_in array and fill
!c    work_out(jx*ix) array.
!c    3. Copy 1d work_out array to 2d (jx,ix) array for passing
!c    to clm.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( .not.ifrest ) then
!       TGB
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2ctb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2ctb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       QB
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2cqb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2cqb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       ZGA
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2czga(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2czga_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       PSB
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2cpsb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2cpsb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       UXB
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2cuxb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2cuxb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       VXB
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2cvxb(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2cvxb_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       RNC
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2crnc(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2crnc_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       RNNC
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2crnnc(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2crnnc_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLS
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2csols(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2csols_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLL
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2csoll(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2csoll_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLSD
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2csolsd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2csolsd_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       SOLLD
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2csolld(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2csolld_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
!       LONGWAVE RAD.
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            work_in(ii) = r2cflwd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,&
                         & jxp*ix,mpi_double_precision,mpi_comm_world,  &
                         & ierr)
        ii = 1
        do j = 1 , jx
          do i = 1 , ix
            r2cflwd_all(j,i) = work_out(ii)
            ii = ii + 1
          end do
        end do
 
      end if  ! if not restart
 
 
!     XLAT in radians
      ii = 1
      do j = 1 , jxp
        do i = 1 , ix
          work_in(ii) = r2cxlat(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,  &
                       & jxp*ix,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , ix
          r2cxlat_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLON in radians
      ii = 1
      do j = 1 , jxp
        do i = 1 , ix
          work_in(ii) = r2cxlon(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,  &
                       & jxp*ix,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , ix
          r2cxlon_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLAT in degrees
      ii = 1
      do j = 1 , jxp
        do i = 1 , ix
          work_in(ii) = r2cxlatd(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,  &
                       & jxp*ix,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , ix
          r2cxlatd_all(j,i) = work_out(ii)
          ii = ii + 1
        end do
      end do
!     XLON in degrees
      ii = 1
      do j = 1 , jxp
        do i = 1 , ix
          work_in(ii) = r2cxlond(j,i)
          ii = ii + 1
        end do
      end do
      call mpi_allgather(work_in,jxp*ix,mpi_double_precision,work_out,  &
                       & jxp*ix,mpi_double_precision,mpi_comm_world,    &
                       & ierr)
      ii = 1
      do j = 1 , jx
        do i = 1 , ix
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
        do i = 1 , ix
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
      lsmlat = ix
 
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
        nstep = ktau
        call rcmdrv()
      end if !end ifrest test
 
!     Initialize ocld2d now that clm has determined the land sea mask
!     Initialize accumulation variables at zero
 
      js = (jxp*myid) + 1
      je = jxp*(myid+1)
      jj = 0
      do j = js , je
        jj = jj + 1
        do i = 1 , ix - 1
 
          if ( .not.ifrest ) then
            do n = 1 , nnsg
              ocld2d(n,i,jj) = dble(landmask(j,i))
              tgb2d(n,i,jj) = tgb(i,jj)
              taf2d(n,i,jj) = tgb(i,jj)
              tlef2d(n,i,jj) = tgb(i,jj)
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
 
          satbrt(i,jj) = clm2bats_veg(j,i)
          if ( clm2bats_veg(j,i)==0 ) satbrt(i,jj) = 15
          if ( satbrt(i,jj).gt.13.9 .and. satbrt(i,jj).lt.15.1 ) then
            veg2d(i,jj) = 0
          else
            veg2d(i,jj) = satbrt(i,jj)
          end if
          svegfrac2d(i,jj) = clm_fracveg(i,j)

          if ( landfrac(j,i)/=1. .and. landfrac(j,i)/=0. )              &
             &    landmask(j,i) = 3.

          ! Set some clm land surface/vegetation variables to the ones
          ! used in RegCM.  Make sure all are consistent  

          satbrt(i,jj) = clm2bats_veg(j,i)
          if ( clm2bats_veg(j,i).eq.0 ) satbrt(i,jj) = 15
          do n = 1 , nnsg
            satbrt1(n,i,jj) = clm2bats_veg(j,i)
            if ( clm2bats_veg(j,i).eq.0 ) satbrt1(n,i,jj) = 15
          end do
          if ( satbrt(i,jj).gt.13.9 .and. satbrt(i,jj).lt.15.1 ) then
            veg2d(i,jj)  = 0
            do n = 1 , nnsg
              veg2d1(n,i,jj)  = 0
            end do
          else
            veg2d(i,jj) = satbrt(i,jj)
            do n = 1 , nnsg
              veg2d1(n,i,jj)  = satbrt(i,jj)
            end do
          end if
          svegfrac2d(i,jj) = clm_fracveg(i,j) 
          do n = 1 , nnsg
            if ( veg2d(i,jj).eq.0 .and. ocld2d(n,i,jj).eq.1 ) then
              veg2d(i,jj)     =  2
              veg2d1(n,i,jj)  =  2
              satbrt1(n,i,jj) =  2
              satbrt(i,jj)    =  2
            end if
          end do
          if ( landfrac(j,i).ne.1 .and. landfrac(j,i).ne.0 ) then
            landmask(j,i) = 3
          endif
        end do
      end do
 
!     deallocate some variables used in CLM initialization only
      if ( allocated(ht_rcm) ) deallocate(ht_rcm)
      if ( allocated(satbrt_clm) ) deallocate(satbrt_clm)
      if ( allocated(init_tgb) ) deallocate(init_tgb)
      if ( allocated(clm2bats_veg) ) deallocate(clm2bats_veg)
      if ( allocated(clm_fracveg) ) deallocate(clm_fracveg)
 
      end subroutine initclm

#endif
#endif
