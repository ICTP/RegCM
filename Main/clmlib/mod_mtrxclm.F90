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

  use mod_dynparam
  use mod_mpmessage
  use mod_clm
  use mod_bats_common
  use mod_bats_mtrxbats
  use mod_bats_drag
  use mod_bats_zengocn

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
  subroutine mtrxclm(ktau)
!
    use atmdrvMod , only : rcmdrv
    use clm_comp , only : clm_run1 , clm_run2
!
    implicit none
    integer(8) , intent(in) :: ktau

    call interfclm(1,ktau)
    call rcmdrv()
    call clm_run1(r2cdoalb,r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call clm_run2(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call interfclm(1,ktau)
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
subroutine initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf)
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
  logical , intent(in) :: ifrest
  type(rcm_time_and_date) , intent(in) :: idate1 , idate2
  real(8) , intent(in) :: dtrad , dtsrf , dx
!
  integer :: ci , cj , i , ii , j , jj , n , ierr
  real(8) , dimension(jxp,iy) :: r2cflwd , r2cpsb , r2cqb ,         &
              r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
              r2csolsd , r2ctb , r2cuxb , r2cvxb , r2cxlat ,        &
              r2cxlatd , r2cxlon , r2cxlond , r2czga
  real(r8) , dimension(jxp*iy) :: work_in
  real(r8) , dimension(jx*iy) :: work_out
  integer :: year , month , day , hour
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
  call split_idate(idate1,year,month,day,hour)
  r2cstart_ymd = year*10000+month*100+day
  r2cstart_tod = idate1%second_of_day
!     stop date and time
  call split_idate(idate2,year,month,day,hour)
  r2cstop_ymd = year*10000+month*100+day
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
        ht_rcm(i,j)    = htf(i,j)
        init_tgb(i,j)  = tsf(i,j)
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
      r2cxlatd(j,i) = xlat(ci,cj)
      r2cxlond(j,i) = xlon(ci,cj)
!         xlat,xlon in radians
      r2cxlat(j,i) = xlat(ci,cj)*degrad
      r2cxlon(j,i) = xlon(ci,cj)*degrad
 
      if ( .not.ifrest ) then
!           T(K) at bottom layer
        r2ctb(j,i) = tatm(cj,ci,kz)
!           Specific Humidity
        r2cqb(j,i) = qvatm(cj,ci,kz)/(d_one+qvatm(cj,ci,kz))
!           Reference Height (m)
        r2czga(j,i) = hgt(ci,kz,cj)
!           Surface winds
        r2cuxb(j,i) = uatm(cj,ci,kz)
        r2cvxb(j,i) = vatm(cj,ci,kz)
!           Surface Pressure in Pa from hPa
        r2cpsb(j,i)   = (sfps(ci,cj)+ptop)*d_1000
        r2crnc(j,i)   = pptc(cj,ci)
        r2crnnc(j,i)  = pptnc(cj,ci)
        r2csols(j,i)  = sols2d(ci,cj)
        r2csoll(j,i)  = soll2d(ci,cj)
        r2csolsd(j,i) = solsd2d(ci,cj)
        r2csolld(j,i) = solld2d(ci,cj)
        r2cflwd(j,i)  = flwd2d(cj,ci)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                       jxp*iy,mpi_double_precision,mycomm,  &
                       ierr)
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
                     jxp*iy,mpi_double_precision,mycomm,    &
                     ierr)
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
                     jxp*iy,mpi_double_precision,mycomm,    &
                     ierr)
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
                     jxp*iy,mpi_double_precision,mycomm,    &
                     ierr)
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
                     jxp*iy,mpi_double_precision,mycomm,    &
                     ierr)
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
          ocld2d(n,j,i) = landmask(jj,i)
          tgbrd(n,j,i) = tground2(j,i)
          taf(n,j,i) = tground2(j,i)
          tlef(n,j,i) = tground2(j,i)
          dew2d(n,j,i) = d_zero
          snag(n,j,i) = d_zero
          sncv(n,j,i) = dmax1(sncv(n,j,i),d_zero)
          sfice(n,j,i) = d_zero
          gwet(n,j,i) = d_half
          sena2d(n,i,j) = d_zero
          evpa2d(n,j,i) = d_zero
          srfrno(n,j,i) = d_zero
          runoff(n,j,i) = d_zero
          ircp(n,j,i) = d_zero
        end do
        fsw2d(j,i)   = d_zero
        flw2d(j,i)   = d_zero
        sabveg(j,i)  = d_zero
        fswa(j,i)  = d_zero
        flwa(j,i)  = d_zero
        prca2d(j,i)  = d_zero
        prnca2d(j,i) = d_zero
        svga(j,i)  = d_zero
        sina(j,i)  = d_zero
 
        ! Set some clm land surface/vegetation variables to the ones
        ! used in RegCM.  Make sure all are consistent  

        lndcat(i,j) = clm2bats_veg(jj,i)
        if ( clm2bats_veg(jj,i) < 0.1D0 ) lndcat(i,j) = 15.0D0
        do n = 1 , nnsg
          lndcat1(n,j,i) = clm2bats_veg(jj,i)
          if ( clm2bats_veg(jj,i) < 0.1D0 ) lndcat1(n,j,i) = 15.0D0
        end do

        veg2d(j,i) = idnint(lndcat(i,j))
        do n = 1 , nnsg
          veg2d1(n,j,i)  = idnint(lndcat1(n,j,i))
        end do

        if ( ( veg2d(j,i) == 14 .or. veg2d(j,i) == 15 ) .and. &
               ldmsk(j,i) /= 0 ) then
          veg2d(j,i)        =  2
          lndcat(i,j) =  d_two
        end if
        do n = 1 , nnsg
          if ( ( veg2d1(n,j,i) == 14 .or. veg2d1(n,j,i) == 15 ) .and. &
               ocld2d(n,j,i) /= 0 ) then
            veg2d1(n,j,i)     =  2
            lndcat1(n,j,i)    =  d_two
          end if
        end do
      end do
    end do
    ! Save CLM modified landuse for restart
    lndcat2d(:,:) = lndcat(:,:)
  end if !end ifrest test

!     deallocate some variables used in CLM initialization only
  if ( allocated(ht_rcm) )       deallocate(ht_rcm)
  if ( allocated(init_tgb) )     deallocate(init_tgb)
  if ( allocated(clm2bats_veg) ) deallocate(clm2bats_veg)
  if ( allocated(clm_fracveg) )  deallocate(clm_fracveg)
 
  end subroutine initclm
!
  subroutine albedoclm(imon,jstart,jend,istart,iend)
 
  use clm_varsur , only : landfrac
  implicit none
  integer , intent(in) :: jstart , jend , istart , iend
  integer , intent(in) :: imon
  integer :: i , j , jj
!
  call albedov(imon,jstart,jend,istart,iend)
! 
!     ****** Section Below added for albedo to be corrected by CLM
!     ****** calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
!     ****** to be equal to 1 which can cause a FPE.  To avoid this
!     ****** use albedo calculated with BATS method when albedo=1
!
  do i = istart , iend
    do j = jstart , jend
      jj = j+(jxp*myid)
      if (ocld2d(1,i,j) /= 0 .and. &
          (d_one-aldirs2d(i,j)) > 1.0D-10) then
        aldirs(j,i) = aldirs2d(i,j)*landfrac(jj,i) + &
                    aldirs(j,i)*(d_one-landfrac(jj,i))
        aldirl(j,i) = aldirl2d(i,j)*landfrac(jj,i) + &
                    aldirl(j,i)*(d_one-landfrac(jj,i))
        aldifs(j,i) = aldifs2d(i,j)*landfrac(jj,i) + &
                    aldifs(j,i)*(d_one-landfrac(jj,i))
        aldifl(j,i) = aldifl2d(i,j)*landfrac(jj,i) + &
                    aldifl(j,i)*(d_one-landfrac(jj,i))
        albvs(j,i)  = aldirs2d(i,j)*landfrac(jj,i) + &
                    albvs(j,i) *(d_one-landfrac(jj,i))
        albvl(j,i)  = aldirl2d(i,j)*landfrac(jj,i) + &
                    albvl(j,i) *(d_one-landfrac(jj,i)) 
      end if
      aldirs_o(j,i-1) = real(aldirs(j,i))
      aldifs_o(j,i-1) = real(aldifs(j,i))
    end do
  end do
 
  end subroutine albedoclm
!
  subroutine interfclm(ivers,ktau)

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
  integer(8) , intent(in) :: ktau
!
  real(8) :: mmpd , wpm2
  integer :: ci , cj , counter , i , ii , iii , j , jj , kk ,       &
              n , nn1 , nnn , nout
  real(8) , dimension(jxp,iy) :: r2cflwd , r2cpsb , r2cqb ,    &
              r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
              r2csolsd , r2ctb , r2cuxb , r2cvxb , r2czga
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
        r2ctb(j,i) = tatm(cj,ci,kz)
!           Specific Humidity ?
        r2cqb(j,i) = qvatm(cj,ci,kz)/(d_one+qvatm(cj,ci,kz))
!           Reference Height (m)
        r2czga(j,i) = hgt(ci,kz,cj)
!           Surface winds
        r2cuxb(j,i) = uatm(cj,ci,kz)
        r2cvxb(j,i) = vatm(cj,ci,kz)
!           Surface Pressure in Pa from cbar
        r2cpsb(j,i) = (sfps(ci,cj)+ptop)*d_1000
!           Rainfall
        r2crnc(j,i) = pptc(cj,ci)
        r2crnnc(j,i) = pptnc(cj,ci)
!           Incident Solar Radiation
        r2csols(j,i) = sols2d(ci,cj)
        r2csoll(j,i) = soll2d(ci,cj)
        r2csolsd(j,i) = solsd2d(ci,cj)
        r2csolld(j,i) = solld2d(ci,cj)
        r2cflwd(j,i) = flwd2d(cj,ci)
 
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
                       workout,13*jxp*iy,mpi_double_precision,      &
                       mycomm,ierr)
 
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
    else if ( ktau+1 == kbats ) then
      mmpd = houpd/(srffrq-xdtsec/secph)
      wpm2 = d_one/((srffrq-xdtsec/secph)*secph)
    else
      mmpd = houpd/srffrq
      wpm2 = d_one/(srffrq*secph)
    end if
 
    call interf(1,jbegin,jendx,2,iym1,ktau)

    if ( iocnflx==2 ) then
      call zengocndrv(jbegin,jendx,2,iym1,ktau)
    end if
 
    do j = jbegin , jendx
      jj = (jxp*myid) + j
      do i = 2 , iym1
        ci = i
        uvdrag(j,i) = d_zero
        hfx(j,i) = d_zero
        qfx(j,i) = d_zero
        tground2(j,i) = d_zero
        tground1(j,i) = d_zero
        tgbb(j,i) = d_zero

        if ( lchem ) then
          ssw2da(j,i) = d_zero
          sdeltk2d(j,i) = d_zero
          sdelqk2d(j,i) = d_zero
          sfracv2d(j,i) = d_zero
          sfracb2d(j,i) = d_zero
          sfracs2d(j,i) = d_zero
        end if

        if ( landmask(jj,ci)==1 ) then
          tground2(j,i) = c2rtgb(jj,ci)
          tground1(j,i) = c2rtgb(jj,ci)
          hfx(j,i) = c2rsenht(jj,ci)
          qfx(j,i) = c2rlatht(jj,ci)
          uvdrag(j,i) = c2ruvdrag(jj,ci)
          tgbb(j,i) = c2rtgbb(jj,ci)
 
          if ( i<=iym1 ) then
            aldirs2d(i,j) = c2ralbdirs(jj,ci)
            aldirl2d(i,j) = c2ralbdirl(jj,ci)
            aldifs2d(i,j) = c2ralbdifs(jj,ci)
            aldifl2d(i,j) = c2ralbdifl(jj,ci)
          end if
 
          do n = 1 , nnsg
            tgrd(n,j,i) = c2rtgb(jj,ci)
            tgbrd(n,j,i) = c2rtgb(jj,ci)
            !supposed to be lower soil layer temp not tgrnd
            taf(n,j,i) = c2r2mt(jj,ci)
            tlef(n,j,i) = c2rtlef(jj,ci)
            tsw(n,j,i) = c2rsmtot(jj,ci)
            rsw(n,j,i) = c2rsm1m(jj,ci)
            ssw(n,j,i) = c2rsm10cm(jj,ci)
            dew2d(n,j,i) = ldew(n,j,i)
            sncv(n,j,i) = c2rsnowc(jj,ci)
            evpa2d(n,j,i) = evpa2d(n,j,i) + dtbat*qfx(j,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*hfx(j,i)
            srfrno(n,j,i) = c2rro_sur(jj,ci)*dtbat
            runoff(n,j,i) = (c2rro_sub(jj,ci)+c2rro_sur(jj,ci))*dtbat
 
            if ( lchem ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw2d(n,i,j)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jj,ci)
              sfracb2d(j,i) = sfracb2d(j,i) + d_one -               &
                             (c2rfvegnosno(jj,ci)+c2rfracsno(jj,ci))
              sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jj,ci)
            end if
          end do
 
          !abt added for 2m humidity when landmask = 1 or 3
          q2d(i,j) = c2r2mq(jj,ci)
!
!             quantities stored on 2d surface array for bats use only
!
          prca2d(j,i) = prca2d(j,i) + dtbat*pptc(j,i)
          prnca2d(j,i) = prnca2d(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd2d(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
        else if ( landmask(jj,ci)==0 ) then !ocean
 
          do n = 1 , nnsg
            uvdrag(j,i) = uvdrag(j,i) + drag(n,j,i)
            hfx(j,i) = hfx(j,i) + sent(n,j,i)
            qfx(j,i) = qfx(j,i) + evpr(n,j,i)
            tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
            tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

            if ( lchem  ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + sigf(n,j,i)
              sfracb2d(j,i) = sfracb2d(j,i) + (d_one-sigf(n,j,i))    &
                              *(d_one-scvk(n,j,i))
              sfracs2d(j,i) = sfracs2d(j,i) + sigf(n,j,i)*wt(n,j,i) &
                              + (d_one-sigf(n,j,i))*scvk(n,j,i)
            end if
 
            if ( iocnflx==1 .or.                                    &
                 (iocnflx==2 .and. ocld2d(n,j,i) /= 0) ) then
              tgbb(j,i) = tgbb(j,i)                     &
                          + ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four+   &
                          lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
            else
              tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
            end if
            ssw(n,j,i)  = dmissval
            rsw(n,j,i)  = dmissval
            tsw(n,j,i)  = dmissval
            trnof(n,j,i)  = dmissval
            srnof(n,j,i) = dmissval
            sncv(n,j,i)  = dmissval
          end do
 
          do n = 1 , nnsg
            taf(n,j,i) = t2m(n,j,i)
            dew2d(n,j,i) = ldew(n,j,i)
            sncv(n,j,i) = sncv(n,j,i)
            evpa2d(n,j,i) = evpa2d(n,j,i) + dtbat*evpr(n,j,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent(n,j,i)
            if ( dabs(trnof(n,j,i)) > 1.0D-10 ) then
              rnos2d(n,i,j) = rnos2d(n,i,j) + &
                              trnof(n,j,i)/secpd*dtbat
            end if
            if ( dabs(srnof(n,j,i)) > 1.0D-10 .and. &
                 dabs(trnof(n,j,i))  > 1.0D-10 ) then
              runoff(n,j,i) = runoff(n,j,i) + &
                      (trnof(n,j,i)-srnof(n,j,i))/secpd*dtbat
            end if
          end do
!
!             quantities stored on 2d surface array for bats use only
!
          prca2d(j,i) = prca2d(j,i) + dtbat*pptc(j,i)
          prnca2d(j,i) = prnca2d(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd2d(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
 
        else if ( landmask(jj,ci)==3 ) then
        !gridcell with some % land and ocean
 
          do n = 1 , nnsg
            uvdrag(j,i) = uvdrag(j,i) + drag(n,j,i)
            hfx(j,i) = hfx(j,i) + sent(n,j,i)
            qfx(j,i) = qfx(j,i) + evpr(n,j,i)
            tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
            tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

            if ( lchem ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw2d(n,i,j)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jj,ci)
              sfracb2d(j,i) = sfracb2d(j,i)                         &
                              + d_one - (c2rfvegnosno(jj,ci)+       &
                              c2rfracsno(jj,ci))
              sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jj,ci)
              ssw2da(j,i) = ssw2da(j,i)*landfrac(jj,ci)             &
                            + (d_one-landfrac(jj,ci))*ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i)*landfrac(jj,ci)         &
                              + (d_one-landfrac(jj,ci))*delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i)*landfrac(jj,ci)         &
                              + (d_one-landfrac(jj,ci))*delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i)*landfrac(jj,ci)         &
                              + (d_one-landfrac(jj,ci))*sigf(n,j,i)
              sfracb2d(j,i) = sfracb2d(j,i)*landfrac(jj,ci)         &
                              + (d_one-landfrac(jj,ci))*            &
                              (d_one-sigf(n,j,i))*(d_one-scvk(n,j,i))
              sfracs2d(j,i) = sfracs2d(j,i)*landfrac(jj,ci)         &
                              + (d_one-landfrac(jj,ci))             &
                              *(sigf(n,j,i)*wt(n,j,i)+(d_one-sigf(n,j,i)) &
                              *scvk(n,j,i))
            end if
 
            if ( iocnflx==1 .or.                                    &
                 (iocnflx==2 .and. ocld2d(n,j,i) /= 0) ) then
              tgbb(j,i) = tgbb(j,i) + &
                        ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four+ &
                            lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
            else
              tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
            end if
          end do
 
          uvdrag(j,i) = uvdrag(j,i)* &
                        (d_one-landfrac(jj,ci))+ &
                        c2ruvdrag(jj,ci)*landfrac(jj,ci)
          hfx(j,i) = hfx(j,i)*(d_one-landfrac(jj,ci)) + &
                           c2rsenht(jj,ci)*landfrac(jj,ci)
          qfx(j,i) = qfx(j,i)*(d_one-landfrac(jj,ci)) + &
                           c2rlatht(jj,ci)*landfrac(jj,ci)
          tground2(j,i) = tground2(j,i)*(d_one-landfrac(jj,ci)) +     &
                         c2rtgb(jj,ci)*landfrac(jj,ci)
          tgbb(j,i) = tgbb(j,i)* (d_one-landfrac(jj,ci)) +   &
                            c2rtgbb(jj,ci)*landfrac(jj,ci)
          tground1(j,i) = tground1(j,i)*(d_one-landfrac(jj,ci)) +     &
                         c2rtgb(jj,ci)*landfrac(jj,ci)
 
          do n = 1 , nnsg
            dew2d(n,j,i) = ldew(n,j,i)
!abt            added below for the landfraction method
            sncv(n,j,i) = c2rsnowc(jj,ci)*landfrac(jj,ci)          &
                           + sncv(n,j,i)*(d_one-landfrac(jj,ci))
            tgrd(n,j,i) = c2rtgb(jj,ci)*landfrac(jj,ci) + tgrd(n,j,i) &
                          *(d_one-landfrac(jj,ci))
            tgbrd(n,j,i) = c2rtgb(jj,ci)*landfrac(jj,ci)            &
                           + tgbrd(n,j,i)*(d_one-landfrac(jj,ci))
            taf(n,j,i) = c2r2mt(jj,ci)*landfrac(jj,ci)            &
                           + t2m(n,j,i)*(d_one-landfrac(jj,ci))
            !note taf is 2m temp not temp in foilage
            tlef(n,j,i) = c2rtlef(jj,ci)*landfrac(jj,ci)          &
                            + tlef(n,j,i)*(d_one-landfrac(jj,ci))
            tsw(n,j,i) = c2rsmtot(jj,ci)*landfrac(jj,ci)          &
                           + tsw(n,j,i)*(d_one-landfrac(jj,ci))
            rsw(n,j,i) = c2rsm1m(jj,ci)*landfrac(jj,ci)           &
                           + rsw(n,j,i)*(d_one-landfrac(jj,ci))
            ssw(n,j,i) = c2rsm10cm(jj,ci)*landfrac(jj,ci)         &
                           + ssw(n,j,i)*(d_one-landfrac(jj,ci))
            q2d(i,j) = c2r2mq(jj,ci)*landfrac(jj,ci) + q2m(n,j,i)  &
                       *(d_one-landfrac(jj,ci))
 
            evpa2d(n,j,i) = evpa2d(n,j,i) + dtbat*qfx(j,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*hfx(j,i)
            srfrno(n,j,i) = c2rro_sur(jj,ci)*dtbat
            runoff(n,j,i) = c2rro_sub(jj,ci)*dtbat + c2rro_sur(jj,ci)*dtbat
          end do
!
!             quantities stored on 2d surface array for bats use only
!
          prca2d(j,i) = prca2d(j,i) + dtbat*pptc(j,i)
          prnca2d(j,i) = prnca2d(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd2d(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
 
        end if
      end do !i loop
 
!         Fill output arrays if needed
 
      if ( mod(ktau+1,kbats) == 0 .or. ktau == 0 ) then
 
        do i = 2 , iym1
          ci = i
 
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
 
          do n = 1 , nnsg
            if ( ocld2d(n,j,i) /= 0 ) then
              u10m_s(n,j,i-1) = real(uatm(j,i,kz))
              v10m_s(n,j,i-1) = real(vatm(j,i,kz))
              tg_s(n,j,i-1) = real(tgrd(n,j,i))
              t2m_s(n,j,i-1) = real(taf(n,j,i))
              u10m_o(j,i-1) = u10m_o(j,i-1) + real(uatm(j,i,kz))
              v10m_o(j,i-1) = v10m_o(j,i-1) + real(vatm(j,i,kz))
              t2m_o(j,i-1) = t2m_o(j,i-1) + real(taf(n,j,i))
              tg_o(j,i-1) = tg_o(j,i-1) + real(tgrd(n,j,i))
            else if ( ocld2d(n,j,i) == 0 ) then
              tg_s(n,j,i-1) = real(tgrd(n,j,i))
              u10m_s(n,j,i-1) = real(u10m(n,j,i))
              v10m_s(n,j,i-1) = real(v10m(n,j,i))
              t2m_s(n,j,i-1) = real(t2m(n,j,i))
              u10m_o(j,i-1) = u10m_o(j,i-1) + real(u10m(n,j,i))
              v10m_o(j,i-1) = v10m_o(j,i-1) + real(v10m(n,j,i))
              t2m_o(j,i-1) = t2m_o(j,i-1) + real(t2m(n,j,i))
              tg_o(j,i-1) = tg_o(j,i-1) + real(tgrd(n,j,i))
            end if
          end do
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1), &
                   sqrt(u10m_o(j,i-1)**2.0+v10m_o(j,i-1)**2.0))
          real_4 = real((sfps(i,j)+ptop)*d_10)
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
 
          drag_o(j,i-1) = 0.0
          q2m_o(j,i-1) = 0.0
          evpa_o(j,i-1) = 0.0
          sena_o(j,i-1) = 0.0
          do n = 1 , nnsg
            if ( ocld2d(n,j,i) /= 0 ) then
              q2m_s(n,j,i-1) = real(q2d(i,j))
              drag_s(n,j,i-1) = real(uvdrag(j,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,j,ci)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,ci,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(j,ci)+prca2d(j,ci))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(j,ci)*mmpd)
              ps_s(n,j,i-1) = real(sfcp(n,j,i)*0.01D0)
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + real(q2d(i,j))
              drag_o(j,i-1) = drag_o(j,i-1) + real(uvdrag(j,i))
              evpa_o(j,i-1) = evpa_o(j,i-1) + real(evpa2d(n,j,ci))
              sena_o(j,i-1) = sena_o(j,i-1) + real(sena2d(n,ci,j))
            else if ( ocld2d(n,j,i) == 0 ) then
              q2m_s(n,j,i-1) = real(q2m(n,j,i))
              drag_s(n,j,i-1) = real(drag(n,j,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,j,i)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(j,i)+prca2d(j,i))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(j,i)*mmpd)
              ps_s(n,j,i-1) = real(sfcp(n,j,i)*0.01D0)
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + real(q2m(n,j,i))
              drag_o(j,i-1) = drag_o(j,i-1) + real(drag(n,j,i))
              evpa_o(j,i-1) = evpa_o(j,i-1) + real(evpa2d(n,j,i))
              sena_o(j,i-1) = sena_o(j,i-1) + real(sena2d(n,i,j))
            end if
          end do
          tpr_o(j,i-1) = real((prnca2d(j,ci)+prca2d(j,ci))*mmpd)
          evpa_o(j,i-1) = evpa_o(j,i-1)*real(mmpd)
          sena_o(j,i-1) = sena_o(j,i-1)*real(wpm2)
          flwa_o(j,i-1) = real(flwa(j,ci)*wpm2)
          fswa_o(j,i-1) = real(fswa(j,ci)*wpm2)
          flwd_o(j,i-1) = real(flwda(j,ci)*wpm2)
          sina_o(j,i-1) = real(sina(j,ci)*wpm2)
          prcv_o(j,i-1) = real(prca2d(j,ci)*mmpd)
          ps_o(j,i-1) = real((sfps(i,j)+ptop)*d_10)
          zpbl_o(j,i-1) = real(zpbl(i,j))
 
          tlef_o(j,i-1) = 0.0
          ssw_o(j,i-1) = 0.0
          rsw_o(j,i-1) = 0.0
          rnos_o(j,i-1) = 0.0
          scv_o(j,i-1) = 0.0
          nnn = 0
          do n = 1 , nnsg
            if ( ocld2d(n,j,ci) /= 0 .and. landmask(jj,ci)/=3 ) then
              tlef_o(j,i-1) = tlef_o(j,i-1) + real(c2rtlef(jj,ci))
              ssw_o(j,i-1) = ssw_o(j,i-1) + real(c2rsm10cm(jj,ci))
              rsw_o(j,i-1) = rsw_o(j,i-1) + real(c2rsm1m(jj,ci))
              rnos_o(j,i-1) = rnos_o(j,i-1) + real(srfrno(n,j,ci))
              scv_o(j,i-1) = scv_o(j,i-1) + real(c2rsnowc(jj,ci))
              tlef_s(n,j,i-1) = real(c2rtlef(jj,ci))
              ssw_s(n,j,i-1) = real(c2rsm10cm(jj,ci))
              rsw_s(n,j,i-1) = real(c2rsm1m(jj,ci))
              rnos_s(n,j,i-1) = real(srfrno(n,j,ci)*mmpd)
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
            evpa2d(n,j,ci) = d_zero
            srfrno(n,j,ci) = d_zero
            sena2d(n,ci,j) = d_zero
          end do
          prnca2d(j,ci) = d_zero
          prca2d(j,ci) = d_zero
          flwa(j,ci) = d_zero
          flwda(j,ci) = d_zero
          fswa(j,ci) = d_zero
          svga(j,ci) = d_zero
          sina(j,ci) = d_zero
 
        end do  ! end of i loop
      end if
    end do      ! end of j loop
  end if        ! end if ivers = 2
 
  end subroutine interfclm
!
end module mod_mtrxclm

#endif
