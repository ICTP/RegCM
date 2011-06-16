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

#ifdef CLM

module mod_mtrxclm

  use mod_runparams
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
!
  subroutine mtrxclm
!
    use atmdrvMod , only : rcmdrv
    use clm_comp , only : clm_run1 , clm_run2
!
    implicit none

    call interfclm(1)
    call rcmdrv()
    call clm_run1(r2cdoalb,r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call clm_run2(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call interfclm(2)
!
  end subroutine mtrxclm
!
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine initclm
!
  use initializeMod
  use shr_orb_mod
  use shr_kind_mod,  only : r8 => shr_kind_r8
  use clm_varpar,    only : lsmlon , lsmlat
  use clm_varsur,    only : landmask , landfrac
  use clm_varsur,    only : r2cimask , init_tgb , r2coutfrq
  use clm_varsur,    only : clm2bats_veg , ht_rcm
  use clm_varsur,    only : clm_fracveg
  use clm_varsur,    only : slmo
  use atmdrvMod
  use program_offMod
  use clm_comp 
  use clmtype
  use perf_mod
  use mpi
! 
  implicit none
!
  integer :: ci , cj , i , ii , j , jj , n , ierr
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
  r2cdtime = idint(dtsrf)
!     start date and time
  r2cstart_ymd = idate1%year*10000+idate1%month*100+idate1%day
  r2cstart_tod = idate1%second_of_day
!     stop date and time
  r2cstop_ymd = idate2%year*10000+idate2%month*100+idate2%day
  r2cstop_tod = idate2%second_of_day
!     calendar type (GREGORIAN not available in regcm)
  if ( ical == noleap ) then
    r2cclndr = 'NO_LEAP'
  else if ( ical == gregorian ) then
    r2cclndr = 'GREGORIAN'
  else
    call fatal(__FILE__,__LINE__,'CLM supports only gregorian and noleap')
  end if
!     don't write to NCAR Mass Store
  r2cmss_irt = 0
!     clm output frequency
  r2coutfrq = idint(clmfrq)
!     radiation calculation frequency
!     regcm: dtrad is in minutes
!     clm: irad is (+) iterations or (-) hours
!     clm: hours gets converted to seconds then divided by dtime
  r2cirad = idint(dtrad*minph/r2cdtime)
!     write output
  if ( ifsrf ) then
    r2cwrtdia = .true.
  else
    r2cwrtdia = .false.
  end if
!     Set grid spacing resolution
  r2cdx = dx
!     Set gridcell area
  r2carea = (dx*d_r1000)*(dx*d_r1000)
!     Set landmask method
  r2cimask = imask
!     Set elevation and BATS landuse type (abt added)
  if ( .not.allocated(ht_rcm) ) allocate(ht_rcm(iy,jx))
  if ( .not.allocated(init_tgb) ) allocate(init_tgb(iy,jx))
  if ( .not.allocated(clm2bats_veg) ) allocate(clm2bats_veg(jx,iy))
  if ( .not.allocated(clm_fracveg) ) allocate(clm_fracveg(iy,jx))
  if ( myid==0 ) then
    do j = 1 , jx
      do i = 1 , iy
        ht_rcm(i,j)      = mddom_io%ht(i,j)
        init_tgb(i,j)  = ts0_io(i,j)
        clm_fracveg(i,j) = d_zero
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
    pptc(:,:)  = d_zero
    pptnc(:,:) = d_zero
!       Radiation
    sols2d(:,:)  = d_zero
    soll2d(:,:)  = d_zero
    solsd2d(:,:) = d_zero
    solld2d(:,:) = d_zero
    flwd2d(:,:)  = d_zero
!       Albedo
!       Set initial albedos to clm dry soil values for mid-colored soils
    aldirs2d(:,:) = 0.16D0
    aldifs2d(:,:) = 0.16D0
    aldirl2d(:,:) = 0.32D0
    aldifl2d(:,:) = 0.32D0
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
      r2cxlond(j,i) = mddom%xlon(ci,cj)
!         xlat,xlon in radians
      r2cxlat(j,i) = mddom%xlat(ci,cj)*degrad
      r2cxlon(j,i) = mddom%xlon(ci,cj)*degrad
 
      if ( .not.ifrest ) then
!           T(K) at bottom layer
        r2ctb(j,i) = tb3d(ci,kz,cj)
!           Specific Humidity
        r2cqb(j,i) = qvb3d(ci,kz,cj)/(d_one+qvb3d(ci,kz,cj))
!           Reference Height (m)
        r2czga(j,i) = za(ci,kz,cj)
!           Surface winds
        r2cuxb(j,i) = ubx3d(ci,kz,cj)
        r2cvxb(j,i) = vbx3d(ci,kz,cj)
!           Surface Pressure in Pa from hPa
        r2cpsb(j,i)   = (sps2%ps(ci,cj)+r8pt)*d_1000
        r2crnc(j,i)   = pptc(ci,cj)
        r2crnnc(j,i)  = pptnc(ci,cj)
        r2csols(j,i)  = sols2d(ci,cj)
        r2csoll(j,i)  = soll2d(ci,cj)
        r2csolsd(j,i) = solsd2d(ci,cj)
        r2csolld(j,i) = solld2d(ci,cj)
        r2cflwd(j,i)  = flwd2d(ci,cj)
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
 
  if ( myid == 0 ) &
     write (6,*) 'Attempting to make atmospheric grid'
  call rcmdrv_init()
  if ( myid == 0 ) &
    write (6,*) 'Successfully make atmospheric grid'
 
!     Initialize radiation and atmosphere variables

  if ( .not. ifrest ) then
    call rcmdrv()
  end if !end ifrest test
 
!     Correct landmask
  do j = 1 , jx
    do i = 1 , iy
      if ( dabs(landfrac(j,i)-d_one) >= 0.1D0 .and. &
           dabs(landfrac(j,i)) >= 0.1D0 ) then
        landmask(j,i) = 3
      end if
    end do
  end do

!     Initialize ocld2d now that clm has determined the land sea mask
!     Initialize accumulation variables at zero
 
  if ( .not.ifrest ) then
    do j = 1 , jxp
      jj = myid*jxp + j
      do i = 1 , iym1
        do n = 1 , nnsg
          ocld2d(n,i,j) = landmask(jj,i)
          tgb2d(n,i,j) = sts2%tg(i,j)
          taf2d(n,i,j) = sts2%tg(i,j)
          tlef2d(n,i,j) = sts2%tg(i,j)
          dew2d(n,i,j) = d_zero
          sag2d(n,i,j) = d_zero
          scv2d(n,i,j) = dmax1(scv2d(n,i,j),d_zero)
          sice2d(n,i,j) = d_zero
          gwet2d(n,i,j) = d_half
          sena2d(n,i,j) = d_zero
          evpa2d(n,i,j) = d_zero
          rnos2d(n,i,j) = d_zero
          rno2d(n,i,j) = d_zero
          ircp2d(n,i,j) = d_zero
        end do
        fsw2d(i,j)   = d_zero
        flw2d(i,j)   = d_zero
        sabv2d(i,j)  = d_zero
        sol2d(i,j)   = d_zero
        fswa2d(i,j)  = d_zero
        flwa2d(i,j)  = d_zero
        prca2d(i,j)  = d_zero
        prnca2d(i,j) = d_zero
        svga2d(i,j)  = d_zero
        sina2d(i,j)  = d_zero
 
        ! Set some clm land surface/vegetation variables to the ones
        ! used in RegCM.  Make sure all are consistent  

        mddom%lndcat(i,j) = clm2bats_veg(jj,i)
        if ( clm2bats_veg(jj,i) < 0.1D0 ) mddom%lndcat(i,j) = 15.0D0
        do n = 1 , nnsg
          lndcat1(n,i,j) = clm2bats_veg(jj,i)
          if ( clm2bats_veg(jj,i) < 0.1D0 ) lndcat1(n,i,j) = 15.0D0
        end do

        veg2d(i,j) = idnint(mddom%lndcat(i,j))
        do n = 1 , nnsg
          veg2d1(n,i,j)  = idnint(lndcat1(n,i,j))
        end do

        if ( ( veg2d(i,j) == 14 .or. veg2d(i,j) == 15 ) .and. &
               ldmsk(i,j) /= 0 ) then
          veg2d(i,j)        =  2
          mddom%lndcat(i,j) =  d_two
        end if
        do n = 1 , nnsg
          if ( ( veg2d1(n,i,j) == 14 .or. veg2d1(n,i,j) == 15 ) .and. &
               ocld2d(n,i,j) /= 0 ) then
            veg2d1(n,i,j)     =  2
            lndcat1(n,i,j)    =  d_two
          end if
        end do
      end do
    end do
    ! Save CLM modified landuse for restart
    lndcat2d(:,:) = mddom%lndcat(:,:)
  end if !end ifrest test

!     deallocate some variables used in CLM initialization only
  if ( allocated(ht_rcm) )       deallocate(ht_rcm)
  if ( allocated(init_tgb) )     deallocate(init_tgb)
  if ( allocated(clm2bats_veg) ) deallocate(clm2bats_veg)
  if ( allocated(clm_fracveg) )  deallocate(clm_fracveg)
 
  end subroutine initclm
!
  subroutine albedoclm(j,iemiss)
 
  use clm_varsur , only : landfrac
  implicit none
!
  integer , intent(in) :: iemiss , j
!
  integer :: i , jj
!
  call albedov(j,iemiss)
! 
!     ****** Section Below added for albedo to be corrected by CLM
!     ****** calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
!     ****** to be equal to 1 which can cause a FPE.  To avoid this
!     ****** use albedo calculated with BATS method when albedo=1
!
  do i = 2 , iym1
    jj = j+(jxp*myid)
    if (ocld2d(1,i,j) /= 0 .and. &
        (d_one-aldirs2d(i,j)) > 1.0D-10) then
      aldirs(i) = aldirs2d(i,j)*landfrac(jj,i) + &
                  aldirs(i)*(d_one-landfrac(jj,i))
      aldirl(i) = aldirl2d(i,j)*landfrac(jj,i) + &
                  aldirl(i)*(d_one-landfrac(jj,i))
      aldifs(i) = aldifs2d(i,j)*landfrac(jj,i) + &
                  aldifs(i)*(d_one-landfrac(jj,i))
      aldifl(i) = aldifl2d(i,j)*landfrac(jj,i) + &
                  aldifl(i)*(d_one-landfrac(jj,i))
      albvs(i)  = aldirs2d(i,j)*landfrac(jj,i) + &
                  albvs(i) *(d_one-landfrac(jj,i))
      albvl(i)  = aldirl2d(i,j)*landfrac(jj,i) + &
                  albvl(i) *(d_one-landfrac(jj,i)) 
    end if
    aldirs_o(j,i-1) = real(aldirs(i))
    aldifs_o(j,i-1) = real(aldifs(i))
  end do
 
  end subroutine albedoclm
!
  subroutine interfclm(ivers)

!=======================================================================
!l  built for clm version 3.0
!=======================================================================
! ivers = 1 : regcm -> clm
! ivers = 2 : clm -> regcm
!
  use clmtype
  use clm_varsur , only : landmask, landfrac
  use clm_varsur , only : c2r_allout,omap_i,omap_j
  use mpi
  implicit none
!
  integer , intent(in) :: ivers
!
  real(8) :: mmpd , wpm2
  integer :: ci , cj , counter , i , ii , iii , j , jj , kk ,       &
            & n , nn1 , nnn , nout
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
        r2cqb(j,i) = qvb3d(ci,kz,cj)/(d_one+qvb3d(ci,kz,cj))
!           Reference Height (m)
        r2czga(j,i) = za(ci,kz,cj)
!           Surface winds
        r2cuxb(j,i) = ubx3d(ci,kz,cj)
        r2cvxb(j,i) = vbx3d(ci,kz,cj)
!           Surface Pressure in Pa from cbar
        r2cpsb(j,i) = (sps2%ps(ci,cj)+r8pt)*d_1000
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
 
    if ( ktau==0 ) then
      mmpd = secpd/dtbat
      wpm2 = d_one/dtbat
    else if ( dble(ktau*dtmin) <= srffrq*minph+0.01D0 ) then
      mmpd = houpd/(srffrq-dtmin/minph)
      wpm2 = d_one/((srffrq-dtmin/minph)*secph)
    else
      mmpd = houpd/srffrq
      wpm2 = d_one/(srffrq*secph)
    end if
 
    do j = jbegin , jendx
      jj = (jxp*myid) + j
 
      call interf(1,j,2,iym1)

      if ( iocnflx==2 ) then
        call zengocndrv(j,2,iym1)
      end if
 
      do i = 2 , iym1
        ci = i
        sfsta%uvdrag(i,j) = d_zero
        sfsta%hfx(i,j) = d_zero
        sfsta%qfx(i,j) = d_zero
        sts2%tg(i,j) = d_zero
        sts1%tg(i,j) = d_zero
        sfsta%tgbb(i,j) = d_zero

        if ( ichem == 1) then
          ssw2da(i,j) = d_zero
          sdeltk2d(i,j) = d_zero
          sdelqk2d(i,j) = d_zero
          sfracv2d(i,j) = d_zero
          sfracb2d(i,j) = d_zero
          sfracs2d(i,j) = d_zero
        end if

        if ( landmask(jj,ci)==1 ) then
          sts2%tg(i,j) = c2rtgb(jj,ci)
          sts1%tg(i,j) = c2rtgb(jj,ci)
          sfsta%hfx(i,j) = c2rsenht(jj,ci)
          sfsta%qfx(i,j) = c2rlatht(jj,ci)
          sfsta%uvdrag(i,j) = c2ruvdrag(jj,ci)
          sfsta%tgbb(i,j) = c2rtgbb(jj,ci)
 
          if ( i<=iym1 ) then
            aldirs2d(i,j) = c2ralbdirs(jj,ci)
            aldirl2d(i,j) = c2ralbdirl(jj,ci)
            aldifs2d(i,j) = c2ralbdifs(jj,ci)
            aldifl2d(i,j) = c2ralbdifl(jj,ci)
          end if
 
          do n = 1 , nnsg
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
 
            if ( ichem == 1 ) then
              ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
              sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
              sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
              sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
              sfracb2d(i,j) = sfracb2d(i,j) + d_one -               &
                           & (c2rfvegnosno(jj,ci)+c2rfracsno(jj,ci))
              sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
            end if
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
          pptnc(i,j) = d_zero
          pptc(i,j) = d_zero
        else if ( landmask(jj,ci)==0 ) then !ocean
 
          do n = 1 , nnsg
            sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
            sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
            sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
            sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
            sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)

            if ( ichem == 1 ) then
              ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
              sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
              sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
              sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
              sfracb2d(i,j) = sfracb2d(i,j) + (d_one-sigf(n,i))    &
                            & *(d_one-scvk(n,i))
              sfracs2d(i,j) = sfracs2d(i,j) + sigf(n,i)*wt(n,i) &
                            & + (d_one-sigf(n,i))*scvk(n,i)
            end if
 
            if ( iocnflx==1 .or.                                    &
               & (iocnflx==2 .and. ocld2d(n,i,j) /= 0) ) then
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                     &
                        & + ((d_one-veg1d(n,i))*tg1d(n,i)**d_four+   &
                        & veg1d(n,i)*tlef1d(n,i)**d_four)**d_rfour
            else
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
            end if
            ssw1d(n,i)  = dmissval
            rsw1d(n,i)  = dmissval
            tsw1d(n,i)  = dmissval
            rno1d(n,i)  = dmissval
            rnos1d(n,i) = dmissval
            scv1d(n,i)  = dmissval
          end do
 
          do n = 1 , nnsg
            tg2d(n,i,j) = tg1d(n,i)
            tgb2d(n,i,j) = tgb1d(n,i)
            taf2d(n,i,j) = t2m1d(n,i)
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
            if ( dabs(rnos1d(n,i)) > 1.0D-10 ) then
              rnos2d(n,i,j) = rnos2d(n,i,j) + &
                              rnos1d(n,i)/secpd*dtbat
            end if
            if ( dabs(rnos1d(n,i)) > 1.0D-10 .and. &
               & dabs(rno1d(n,i))  > 1.0D-10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + &
                      (rno1d(n,i)-rnos1d(n,i))/secpd*dtbat
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
          pptnc(i,j) = d_zero
          pptc(i,j) = d_zero
 
        else if ( landmask(jj,ci)==3 ) then
        !gridcell with some % land and ocean
 
          do n = 1 , nnsg
            sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
            sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
            sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
            sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
            sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)

            if ( ichem == 1 ) then
              ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
              sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
              sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
              sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
              sfracb2d(i,j) = sfracb2d(i,j)                         &
                            & + d_one - (c2rfvegnosno(jj,ci)+       &
                            & c2rfracsno(jj,ci))
              sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
              ssw2da(i,j) = ssw2da(i,j)*landfrac(jj,ci)             &
                          & + (d_one-landfrac(jj,ci))*ssw1d(n,i)
              sdeltk2d(i,j) = sdeltk2d(i,j)*landfrac(jj,ci)         &
                            & + (d_one-landfrac(jj,ci))*delt1d(n,i)
              sdelqk2d(i,j) = sdelqk2d(i,j)*landfrac(jj,ci)         &
                            & + (d_one-landfrac(jj,ci))*delq1d(n,i)
              sfracv2d(i,j) = sfracv2d(i,j)*landfrac(jj,ci)         &
                            & + (d_one-landfrac(jj,ci))*sigf(n,i)
              sfracb2d(i,j) = sfracb2d(i,j)*landfrac(jj,ci)         &
                            & + (d_one-landfrac(jj,ci))*            &
                            & (d_one-sigf(n,i))*(d_one-scvk(n,i))
              sfracs2d(i,j) = sfracs2d(i,j)*landfrac(jj,ci)         &
                            & + (d_one-landfrac(jj,ci))             &
                            & *(sigf(n,i)*wt(n,i)+(d_one-sigf(n,i)) &
                            & *scvk(n,i))
            end if
 
            if ( iocnflx==1 .or.                                    &
               & (iocnflx==2 .and. ocld2d(n,i,j) /= 0) ) then
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                     &
                        & + ((d_one-veg1d(n,i))*tg1d(n,i)**d_four+   &
                        &   veg1d(n,i)*tlef1d(n,i)**d_four)**d_rfour
            else
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
            end if
          end do
 
          sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)* &
                        (d_one-landfrac(jj,ci))+ &
                        c2ruvdrag(jj,ci)*landfrac(jj,ci)
          sfsta%hfx(i,j) = sfsta%hfx(i,j)*(d_one-landfrac(jj,ci)) + &
                           c2rsenht(jj,ci)*landfrac(jj,ci)
          sfsta%qfx(i,j) = sfsta%qfx(i,j)*(d_one-landfrac(jj,ci)) + &
                           c2rlatht(jj,ci)*landfrac(jj,ci)
          sts2%tg(i,j) = sts2%tg(i,j)*(d_one-landfrac(jj,ci)) +     &
                         c2rtgb(jj,ci)*landfrac(jj,ci)
          sfsta%tgbb(i,j) = sfsta%tgbb(i,j)* &
                            (d_one-landfrac(jj,ci)) +   &
                            c2rtgbb(jj,ci)*landfrac(jj,ci)
          sts1%tg(i,j) = sts1%tg(i,j)*(d_one-landfrac(jj,ci)) +     &
                         c2rtgb(jj,ci)*landfrac(jj,ci)
 
          do n = 1 , nnsg
            dew2d(n,i,j) = ldew1d(n,i)
            sag2d(n,i,j) = sag1d(n,i)
            scv2d(n,i,j) = scv1d(n,i)
            sice2d(n,i,j) = sice1d(n,i)
            gwet2d(n,i,j) = gwet1d(n,i)
            ircp2d(n,i,j) = ircp1d(n,i)
!abt            added below for the landfraction method
            scv2d(n,i,j) = c2rsnowc(jj,ci)*landfrac(jj,ci)          &
                         & + scv1d(n,i)*(d_one-landfrac(jj,ci))
            tg2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci) + tg1d(n,i) &
                        & *(d_one-landfrac(jj,ci))
            tgb2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci)            &
                         & + tgb1d(n,i)*(d_one-landfrac(jj,ci))
            taf2d(n,i,j) = c2r2mt(jj,ci)*landfrac(jj,ci)            &
                         & + t2m1d(n,i)*(d_one-landfrac(jj,ci))
            !note taf2d is 2m temp not temp in foilage
            tlef2d(n,i,j) = c2rtlef(jj,ci)*landfrac(jj,ci)          &
                          & + tlef1d(n,i)*(d_one-landfrac(jj,ci))
            swt2d(n,i,j) = c2rsmtot(jj,ci)*landfrac(jj,ci)          &
                         & + tsw1d(n,i)*(d_one-landfrac(jj,ci))
            srw2d(n,i,j) = c2rsm1m(jj,ci)*landfrac(jj,ci)           &
                         & + rsw1d(n,i)*(d_one-landfrac(jj,ci))
            ssw2d(n,i,j) = c2rsm10cm(jj,ci)*landfrac(jj,ci)         &
                         & + ssw1d(n,i)*(d_one-landfrac(jj,ci))
            q2d(i,j) = c2r2mq(jj,ci)*landfrac(jj,ci) + q2m1d(n,i)  &
                     & *(d_one-landfrac(jj,ci))
 
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*sfsta%qfx(i,j)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sfsta%hfx(i,j)
            rnos2d(n,i,j) = c2rro_sur(jj,ci)*dtbat
            rno2d(n,i,j) = c2rro_sub(jj,ci)*dtbat + c2rro_sur(jj,ci)&
                         & *dtbat
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
          pptnc(i,j) = d_zero
          pptc(i,j) = d_zero
 
        end if
      end do !i loop
 
!         Fill output arrays if needed
 
      if ( mod(ntime+idnint(dtmin*minph),nsrffrq) == 0 .or. ktau == 0 ) then
 
        do i = 2 , iym1
          ci = i
 
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
 
          do n = 1 , nnsg
            if ( ocld2d(n,i,j) /= 0 ) then
              u10m_s(n,j,i-1) = real(ubx3d(i,kz,j))
              v10m_s(n,j,i-1) = real(vbx3d(i,kz,j))
              tg_s(n,j,i-1) = real(tg2d(n,i,j))
              t2m_s(n,j,i-1) = real(taf2d(n,i,j))
              u10m_o(j,i-1) = u10m_o(j,i-1) + real(ubx3d(i,kz,j))
              v10m_o(j,i-1) = v10m_o(j,i-1) + real(vbx3d(i,kz,j))
              t2m_o(j,i-1) = t2m_o(j,i-1) + real(taf2d(n,i,j))
              tg_o(j,i-1) = tg_o(j,i-1) + real(tg2d(n,i,j))
            else if ( ocld2d(n,i,j) == 0 ) then
              tg_s(n,j,i-1) = real(tg1d(n,i))
              u10m_s(n,j,i-1) = real(u10m1d(n,i))
              v10m_s(n,j,i-1) = real(v10m1d(n,i))
              t2m_s(n,j,i-1) = real(t2m1d(n,i))
              u10m_o(j,i-1) = u10m_o(j,i-1) + real(u10m1d(n,i))
              v10m_o(j,i-1) = v10m_o(j,i-1) + real(v10m1d(n,i))
              t2m_o(j,i-1) = t2m_o(j,i-1) + real(t2m1d(n,i))
              tg_o(j,i-1) = tg_o(j,i-1) + real(tg1d(n,i))
            end if
          end do
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1), &
                   sqrt(u10m_o(j,i-1)**2.0+v10m_o(j,i-1)**2.0))
          real_4 = real((sps2%ps(i,j)+r8pt)*d_10)
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
 
          drag_o(j,i-1) = 0.0
          q2m_o(j,i-1) = 0.0
          evpa_o(j,i-1) = 0.0
          sena_o(j,i-1) = 0.0
          do n = 1 , nnsg
            if ( ocld2d(n,i,j) /= 0 ) then
              q2m_s(n,j,i-1) = real(q2d(i,j))
              drag_s(n,j,i-1) = real(sfsta%uvdrag(i,j))
              evpa_s(n,j,i-1) = real(evpa2d(n,ci,j)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,ci,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(ci,j)+prca2d(ci,j))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(ci,j)*mmpd)
              ps_s(n,j,i-1) = real(p1d(n,i)*0.01D0)
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + real(q2d(i,j))
              drag_o(j,i-1) = drag_o(j,i-1) + real(sfsta%uvdrag(i,j))
              evpa_o(j,i-1) = evpa_o(j,i-1) + real(evpa2d(n,ci,j))
              sena_o(j,i-1) = sena_o(j,i-1) + real(sena2d(n,ci,j))
            else if ( ocld2d(n,i,j) == 0 ) then
              q2m_s(n,j,i-1) = real(q2m1d(n,i))
              drag_s(n,j,i-1) = real(drag1d(n,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,i,j)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(i,j)*mmpd)
              ps_s(n,j,i-1) = real(p1d(n,i)*0.01D0)
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + real(q2m1d(n,i))
              drag_o(j,i-1) = drag_o(j,i-1) + real(drag1d(n,i))
              evpa_o(j,i-1) = evpa_o(j,i-1) + real(evpa2d(n,i,j))
              sena_o(j,i-1) = sena_o(j,i-1) + real(sena2d(n,i,j))
            end if
          end do
          tpr_o(j,i-1) = real((prnca2d(ci,j)+prca2d(ci,j))*mmpd)
          evpa_o(j,i-1) = evpa_o(j,i-1)*real(mmpd)
          sena_o(j,i-1) = sena_o(j,i-1)*real(wpm2)
          flwa_o(j,i-1) = real(flwa2d(ci,j)*wpm2)
          fswa_o(j,i-1) = real(fswa2d(ci,j)*wpm2)
          flwd_o(j,i-1) = real(flwda2d(ci,j)*wpm2)
          sina_o(j,i-1) = real(sina2d(ci,j)*wpm2)
          prcv_o(j,i-1) = real(prca2d(ci,j)*mmpd)
          ps_o(j,i-1) = real((sps2%ps(i,j)+r8pt)*d_10)
          zpbl_o(j,i-1) = real(sfsta%zpbl(i,j))
 
          tlef_o(j,i-1) = 0.0
          ssw_o(j,i-1) = 0.0
          rsw_o(j,i-1) = 0.0
          rnos_o(j,i-1) = 0.0
          scv_o(j,i-1) = 0.0
          nnn = 0
          do n = 1 , nnsg
            if ( ocld2d(n,ci,j) /= 0 .and. landmask(jj,ci)/=3 ) then
              tlef_o(j,i-1) = tlef_o(j,i-1) + real(c2rtlef(jj,ci))
              ssw_o(j,i-1) = ssw_o(j,i-1) + real(c2rsm10cm(jj,ci))
              rsw_o(j,i-1) = rsw_o(j,i-1) + real(c2rsm1m(jj,ci))
              rnos_o(j,i-1) = rnos_o(j,i-1) + real(rnos2d(n,ci,j))
              scv_o(j,i-1) = scv_o(j,i-1) + real(c2rsnowc(jj,ci))
              tlef_s(n,j,i-1) = real(c2rtlef(jj,ci))
              ssw_s(n,j,i-1) = real(c2rsm10cm(jj,ci))
              rsw_s(n,j,i-1) = real(c2rsm1m(jj,ci))
              rnos_s(n,j,i-1) = real(rnos2d(n,ci,j)*mmpd)
              scv_s(n,j,i-1) = real(c2rsnowc(jj,ci))
              nnn = nnn + 1
            else
              tlef_s(n,j,i-1) = smissval
              ssw_s(n,j,i-1) = smissval
              rsw_s(n,j,i-1) = smissval
              rnos_s(n,j,i-1) = smissval
              scv_s(n,j,i-1) = smissval
            end if
          end do
          if ( nnn>=max0(nnsg/2,1) ) then
            tlef_o(j,i-1) = tlef_o(j,i-1)/real(nnn)
            ssw_o(j,i-1) = ssw_o(j,i-1)/real(nnn)
            rsw_o(j,i-1) = rsw_o(j,i-1)/real(nnn)
            rnos_o(j,i-1) = rnos_o(j,i-1)/real(nnn)*real(mmpd)
            scv_o(j,i-1) = scv_o(j,i-1)/real(nnn)
          else
            tlef_o(j,i-1) = smissval
            ssw_o(j,i-1) = smissval
            rsw_o(j,i-1) = smissval
            rnos_o(j,i-1) = smissval
            scv_o(j,i-1) = smissval
          end if
!               ******    reset accumulation arrays to zero
          do n = 1 , nnsg
            evpa2d(n,ci,j) = d_zero
            rnos2d(n,ci,j) = d_zero
            sena2d(n,ci,j) = d_zero
          end do
          prnca2d(ci,j) = d_zero
          prca2d(ci,j) = d_zero
          flwa2d(ci,j) = d_zero
          flwda2d(ci,j) = d_zero
          fswa2d(ci,j) = d_zero
          svga2d(ci,j) = d_zero
          sina2d(ci,j) = d_zero
 
        end do  ! end of i loop
      end if
    end do      ! end of j loop
  end if        ! end if ivers = 2
 
  end subroutine interfclm
!
end module mod_mtrxclm

#endif
