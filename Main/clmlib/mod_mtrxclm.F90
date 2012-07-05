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
  use mod_runparams , only : idate0
  use mod_mppparam , only : iqv
  use mod_mpmessage
  use mod_service
  use mod_mppparam
  use mod_date
  use mod_clm
  use mod_bats_common
  use mod_bats_mtrxbats
  use mod_bats_drag
  use mod_bats_zengocn

  use clm_time_manager , only : get_curr_calday
  use shr_orb_mod , only : shr_orb_cosz , shr_orb_decl , &
                           shr_orb_params
  private

  public :: mtrxclm
  public :: initclm
  public :: interfclm
  public :: solar_clm
  public :: zenit_clm
  public :: albedoclm

  interface fill_frame
    module procedure fill_frame2d , fill_frame3d , fill_frame4d
  end interface fill_frame

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
    call interfclm(2,ktau)
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
  use clm_varpar,    only : lsmlon , lsmlat
  use clm_varsur,    only : landmask , landfrac , satbrt_clm
  use clm_varsur,    only : r2cimask , init_tgb , r2coutfrq
  use clm_varsur,    only : clm2bats_veg , ht_rcm
  use clm_varsur,    only : clm_fracveg
  use atmdrvMod
  use program_offMod
  use clm_comp 
  use clmtype
  use perf_mod
! 
  implicit none
!
  logical , intent(in) :: ifrest
  type(rcm_time_and_date) , intent(in) :: idate1 , idate2
  real(8) , intent(in) :: dtrad , dtsrf , dx
!
  integer :: i , j , jj , n
  integer :: year , month , day , hour
  !
  ! Initialize run control variables for clm
  !
  if ( ifrest ) then
    r2cnsrest = 1
  else
    r2cnsrest = 0
  end if
  ! land surface timestep
  r2cdtime = idint(dtsrf)
  ! start date and time
  call split_idate(idate1,year,month,day,hour)
  r2cstart_ymd = year*10000+month*100+day
  r2cstart_tod = idate1%second_of_day
  ! stop date and time
  call split_idate(idate2,year,month,day,hour)
  r2cstop_ymd = year*10000+month*100+day
  r2cstop_tod = idate2%second_of_day
  r2cclndr = calstr(idate0%calendar)
  r2cmss_irt = 0
  ! clm output frequency
  r2coutfrq = idint(clmfrq)
  ! radiation calculation frequency
  ! regcm: dtrad is in minutes
  ! clm: irad is (+) iterations or (-) hours
  ! clm: hours gets converted to seconds then divided by dtime
  r2cirad = idint(dtrad*minph/r2cdtime)
  ! write output
  if ( ifsrf ) then
    r2cwrtdia = .true.
  else
    r2cwrtdia = .false.
  end if
  ! Set grid spacing resolution
  r2cdx = dx
  ! Set gridcell area
  r2carea = (dx*d_r1000)*(dx*d_r1000)
  ! Set landmask method
  r2cimask = imask
  ! Set elevation and BATS landuse type (abt added)
  if ( .not.allocated(ht_rcm) ) allocate(ht_rcm(jx,iy))
  if ( .not.allocated(init_tgb) ) allocate(init_tgb(jx,iy))
  if ( .not.allocated(satbrt_clm) ) allocate(satbrt_clm(jx,iy))
  if ( .not.allocated(clm_fracveg) ) allocate(clm_fracveg(jx,iy))
  if ( .not.allocated(clm2bats_veg) ) allocate(clm2bats_veg(jx,iy))
  if ( myid==0 ) then
    do j = 1 , jx
      do i = 1 , iy
        ht_rcm(j,i)    = htf(j,i)
        satbrt_clm(j,i) = lndcatf(j,i)
        init_tgb(j,i)  = tsf(j,i)
        clm_fracveg(j,i) = d_zero
      end do
    end do
  end if
  !
  ! End of clm run control variable initialization
  !
  ! Assign regcm values to the passed variables
  ! regcm writes uneven # of j values to arrays. for now fix by
  ! copying neighboring values
  !
  if ( .not. ifrest ) then
    ! Rainfall
    pptc(:,:)  = d_zero
    pptnc(:,:) = d_zero
    ! Radiation
    sols2d(:,:)  = d_zero
    soll2d(:,:)  = d_zero
    solsd2d(:,:) = d_zero
    solld2d(:,:) = d_zero
    flwd(:,:)  = d_zero
    ! Albedo
    ! Set initial albedos to clm dry soil values for mid-colored soils
    aldirs(:,:) = 0.16D0
    aldifs(:,:) = 0.16D0
    aldirl(:,:) = 0.32D0
    aldifl(:,:) = 0.32D0
  end if
 
  call fill_frame(xlat,r2cxlatd)
  call fill_frame(xlon,r2cxlond)
  r2cxlat = r2cxlatd*degrad
  r2cxlon = r2cxlond*degrad
  if ( .not.ifrest ) then
    call fill_frame(tatm,r2ctb)
    call fill_frame(qxatm,r2cqb,iqv)
    r2cqb = r2cqb/(d_one+r2cqb)
    call fill_frame(hgt,r2czga)
    call fill_frame(uatm,r2cuxb)
    call fill_frame(vatm,r2cvxb)
    call fill_frame(sfps,r2cpsb)
    r2cpsb = (r2cpsb+ptop)*d_1000
    call fill_frame(pptc,r2crnc)
    call fill_frame(pptnc,r2crnnc)
    call fill_frame(sols2d,r2csols)
    call fill_frame(soll2d,r2csoll)
    call fill_frame(solsd2d,r2csolsd)
    call fill_frame(solld2d,r2csolld)
    call fill_frame(flwd,r2cflwd)
  end if

  !
  !    Gather jxp values of each nproc work_in array and fill
  !    work_out(jx*iy) array.
  !    3. Copy 1d work_out array to 2d (jx,iy) array for passing
  !    to clm.
  !
  if ( .not.ifrest ) then
    call grid_fill(r2ctb,r2ctb_all)
    call grid_fill(r2cqb,r2cqb_all)
    call grid_fill(r2czga,r2czga_all)
    call grid_fill(r2cpsb,r2cpsb_all)
    call grid_fill(r2cuxb,r2cuxb_all)
    call grid_fill(r2cvxb,r2cvxb_all)
    call grid_fill(r2crnc,r2crnc_all)
    call grid_fill(r2crnnc,r2crnnc_all)
    call grid_fill(r2csols,r2csols_all)
    call grid_fill(r2csoll,r2csoll_all)
    call grid_fill(r2csolsd,r2csolsd_all)
    call grid_fill(r2csolld,r2csolld_all)
    call grid_fill(r2cflwd,r2cflwd_all)
  end if  ! if not restart

  call grid_fill(r2cxlat,r2cxlat_all)
  call grid_fill(r2cxlon,r2cxlon_all)
  call grid_fill(r2cxlatd,r2cxlatd_all)
  call grid_fill(r2cxlond,r2cxlond_all)
  !
  ! Set grid edges
  !
  r2cedgen = r2cxlatd_all(1,1)
  r2cedges = r2cxlatd_all(1,1)
  r2cedgee = r2cxlond_all(1,1)
  r2cedgew = r2cxlond_all(1,1)
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
 
  ! Initialize radiation and atmosphere variables
  if ( .not. ifrest ) then
    call rcmdrv()
  end if !end ifrest test
 
  ! Correct landmask
  do j = 1 , jx
    do i = 1 , iy
      if ( dabs(landfrac(j,i)-d_one) >= 0.1D0 .and. &
           dabs(landfrac(j,i)) >= 0.1D0 ) then
        landmask(j,i) = 3
      end if
    end do
  end do

  ! Initialize ldmsk1 now that clm has determined the land sea mask
  ! Initialize accumulation variables at zero
 
  if ( .not.ifrest ) then
    do n = 1 , nnsg
      do i = ici1 , ici2
        do j = jci1 , jci2
          jj = myid*jxp + j
          ldmsk1(n,j,i) = landmask(jj,i)
          tgbrd(n,j,i) = tground2(j,i)
          taf(n,j,i) = tground2(j,i)
          tlef(n,j,i) = tground2(j,i)
          ldew(n,j,i) = d_zero
          snag(n,j,i) = d_zero
          sncv(n,j,i) = dmax1(sncv(n,j,i),d_zero)
          sfice(n,j,i) = d_zero
          gwet(n,j,i) = d_half
          sena(n,j,i) = d_zero
          evpa(n,j,i) = d_zero
          srfrna(n,j,i) = d_zero
          runoff(n,j,i) = d_zero
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        jj = myid*jxp + j
        fsw(j,i)   = d_zero
        flw(j,i)   = d_zero
        sabveg(j,i)  = d_zero
        fswa(j,i)  = d_zero
        flwa(j,i)  = d_zero
        prca(j,i)  = d_zero
        prnca(j,i) = d_zero
        svga(j,i)  = d_zero
        sina(j,i)  = d_zero
 
        ! Set some clm land surface/vegetation variables to the ones
        ! used in RegCM.  Make sure all are consistent  

        lndcat(j,i) = clm2bats_veg(jj,i)
        if ( clm2bats_veg(jj,i) < 0.1D0 ) lndcat(j,i) = 15.0D0
        do n = 1 , nnsg
          lndcat1(n,j,i) = clm2bats_veg(jj,i)
          if ( clm2bats_veg(jj,i) < 0.1D0 ) lndcat1(n,j,i) = 15.0D0
        end do

        iveg(j,i) = idnint(lndcat(j,i))
        do n = 1 , nnsg
          iveg1(n,j,i) = idnint(lndcat1(n,j,i))
        end do

        if ( ( iveg(j,i) == 14 .or. iveg(j,i) == 15 ) .and. &
               ldmsk(j,i) /= 0 ) then
          iveg(j,i)   =  2
          lndcat(j,i) =  d_two
        end if
        do n = 1 , nnsg
          if ( ( iveg1(n,j,i) == 14 .or. iveg1(n,j,i) == 15 ) .and. &
               ldmsk1(n,j,i) /= 0 ) then
            iveg1(n,j,i)   =  2
            lndcat1(n,j,i) =  d_two
          end if
        end do
      end do
    end do
    ! Save CLM modified landuse for restart
    lndcat2d(jci1:jci2,ici1:ici2) = lndcat(jci1:jci2,ici1:ici2)
  end if !end ifrest test

  ! deallocate some variables used in CLM initialization only
  if ( allocated(ht_rcm) )       deallocate(ht_rcm)
  if ( allocated(init_tgb) )     deallocate(init_tgb)
  if ( allocated(clm2bats_veg) ) deallocate(clm2bats_veg)
  if ( allocated(clm_fracveg) )  deallocate(clm_fracveg)
 
  end subroutine initclm
!
  subroutine albedoclm(imon)
 
  use clm_varsur , only : landfrac
  implicit none
  integer , intent(in) :: imon
  integer :: i , j , jj
!
  call albedobats(imon)
! 
!     ****** Section Below added for albedo to be corrected by CLM
!     ****** calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
!     ****** to be equal to 1 which can cause a FPE.  To avoid this
!     ****** use albedo calculated with BATS method when albedo=1
!
  do i = ici1 , ici2
    do j = jci1 , jci2
      jj = j+(jxp*myid)
      if (ldmsk1(1,j,i) /= 0 .and. &
          (d_one-aldirs(j,i)) > 1.0D-10) then
        aldirs(j,i) = aldirs(j,i)*landfrac(jj,i) + &
                    aldirs(j,i)*(d_one-landfrac(jj,i))
        aldirl(j,i) = aldirl(j,i)*landfrac(jj,i) + &
                    aldirl(j,i)*(d_one-landfrac(jj,i))
        aldifs(j,i) = aldifs(j,i)*landfrac(jj,i) + &
                    aldifs(j,i)*(d_one-landfrac(jj,i))
        aldifl(j,i) = aldifl(j,i)*landfrac(jj,i) + &
                    aldifl(j,i)*(d_one-landfrac(jj,i))
        albvs(j,i)  = aldirs(j,i)*landfrac(jj,i) + &
                    albvs(j,i) *(d_one-landfrac(jj,i))
        albvl(j,i)  = aldirl(j,i)*landfrac(jj,i) + &
                    albvl(j,i) *(d_one-landfrac(jj,i)) 
      end if
      aldirs_o(j,i) = real(aldirs(j,i))
      aldifs_o(j,i) = real(aldifs(j,i))
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
  implicit none
!
  integer , intent(in) :: ivers
  integer(8) , intent(in) :: ktau
!
  real(8) :: mmpd , wpm2
  integer :: i , ii , iii , j , jj , kk , n , nn1 , nnn , nout
  real(4) :: real_4
!
  if ( ivers == 1 ) then
 
    call fill_frame(tatm,r2ctb)
    call fill_frame(qxatm,r2cqb,iqv)
    r2cqb = r2cqb/(d_one+r2cqb)
    call fill_frame(hgt,r2czga)
    call fill_frame(uatm,r2cuxb)
    call fill_frame(vatm,r2cvxb)
    call fill_frame(sfps,r2cpsb)
    r2cpsb = (r2cpsb+ptop)*d_1000
    call fill_frame(pptc,r2crnc)
    call fill_frame(pptnc,r2crnnc)
    call fill_frame(sols2d,r2csols)
    call fill_frame(soll2d,r2csoll)
    call fill_frame(solsd2d,r2csolsd)
    call fill_frame(solld2d,r2csolld)
    call fill_frame(flwd,r2cflwd)
 
    call grid_fill(r2ctb,r2ctb_all)
    call grid_fill(r2cqb,r2cqb_all)
    call grid_fill(r2czga,r2czga_all)
    call grid_fill(r2cpsb,r2cpsb_all)
    call grid_fill(r2cuxb,r2cuxb_all)
    call grid_fill(r2cvxb,r2cvxb_all)
    call grid_fill(r2crnc,r2crnc_all)
    call grid_fill(r2crnnc,r2crnnc_all)
    call grid_fill(r2csols,r2csols_all)
    call grid_fill(r2csoll,r2csoll_all)
    call grid_fill(r2csolsd,r2csolsd_all)
    call grid_fill(r2csolld,r2csolld_all)
    call grid_fill(r2cflwd,r2cflwd_all)

  else if ( ivers == 2 ) then ! end of ivers = 1
    iii = 0
    jj = 1
    if ( ichem == 1 ) then
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
        if ( ichem == 1 ) then
          c2rfracsno(j,i) = c2r_allout(ii+(20*kk)+iii)
          c2rfvegnosno(j,i) = c2r_allout(ii+(21*kk)+iii)
        end if
        jj = jj + 1
      end do
      iii = iii + c2rngc(nn1)*nout
    end do
 
    if ( ktau == 0 .and. debug_level == 2 ) then
      mmpd = secpd/dtbat
      wpm2 = d_one/dtbat
    else if ( ktau+1 == kbats .and. debug_level == 2 ) then
      mmpd = houpd/(srffrq-xdtsec/secph)
      wpm2 = d_one/((srffrq-xdtsec/secph)*secph)
    else
      mmpd = houpd/srffrq
      wpm2 = d_one/(srffrq*secph)
    end if
 
    call interf(1,ktau)

    if ( iocnflx == 2 ) then
      call zengocndrv(ktau)
    end if
 
    do i = ici1 , ici2
      do j = jci1 , jci2

        jj = (jxp*myid) + j

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

        if ( landmask(jj,i)==1 ) then
          tground2(j,i) = c2rtgb(jj,i)
          tground1(j,i) = c2rtgb(jj,i)
          hfx(j,i) = c2rsenht(jj,i)
          qfx(j,i) = c2rlatht(jj,i)
          uvdrag(j,i) = c2ruvdrag(jj,i)
          tgbb(j,i) = c2rtgbb(jj,i)
 
          if ( i<=iym1 ) then
            aldirs(j,i) = c2ralbdirs(jj,i)
            aldirl(j,i) = c2ralbdirl(jj,i)
            aldifs(j,i) = c2ralbdifs(jj,i)
            aldifl(j,i) = c2ralbdifl(jj,i)
          end if
 
          do n = 1 , nnsg
            tgrd(n,j,i) = c2rtgb(jj,i)
            tgbrd(n,j,i) = c2rtgb(jj,i)
            !supposed to be lower soil layer temp not tgrnd
            taf(n,j,i) = c2r2mt(jj,i)
            tlef(n,j,i) = c2rtlef(jj,i)
            tsw(n,j,i) = c2rsmtot(jj,i)
            rsw(n,j,i) = c2rsm1m(jj,i)
            ssw(n,j,i) = c2rsm10cm(jj,i)
            ldew(n,j,i) = ldew(n,j,i)
            sncv(n,j,i) = c2rsnowc(jj,i)
            evpa(n,j,i) = evpa(n,j,i) + dtbat*qfx(j,i)
            sena(n,j,i) = sena(n,j,i) + dtbat*hfx(j,i)
            srfrna(n,j,i) = c2rro_sur(jj,i)*dtbat
            runoff(n,j,i) = (c2rro_sub(jj,i)+c2rro_sur(jj,i))*dtbat
 
            if ( lchem ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jj,i)
              sfracb2d(j,i) = sfracb2d(j,i) + d_one -               &
                             (c2rfvegnosno(jj,i)+c2rfracsno(jj,i))
              sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jj,i)
            end if
          end do
 
          !abt added for 2m humidity when landmask = 1 or 3
          q2d(j,i) = c2r2mq(jj,i)
          prca(j,i) = prca(j,i) + dtbat*pptc(j,i)
          prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
        else if ( landmask(jj,i) == 0 ) then !ocean
 
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
                 (iocnflx==2 .and. ldmsk1(n,j,i) /= 0) ) then
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
            sncv(n,j,i) = sncv(n,j,i)
            evpa(n,j,i) = evpa(n,j,i) + dtbat*evpr(n,j,i)
            sena(n,j,i) = sena(n,j,i) + dtbat*sent(n,j,i)
            if ( dabs(trnof(n,j,i)) > 1.0D-10 ) then
              srfrna(n,j,i) = srfrna(n,j,i) + &
                              trnof(n,j,i)/secpd*dtbat
            end if
            if ( dabs(srnof(n,j,i)) > 1.0D-10 .and. &
                 dabs(trnof(n,j,i))  > 1.0D-10 ) then
              runoff(n,j,i) = runoff(n,j,i) + &
                      (trnof(n,j,i)-srnof(n,j,i))/secpd*dtbat
            end if
          end do
          prca(j,i) = prca(j,i) + dtbat*pptc(j,i)
          prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
 
        else if ( landmask(jj,i) == 3 ) then
          !gridcell with some % land and ocean
          do n = 1 , nnsg
            uvdrag(j,i) = uvdrag(j,i) + drag(n,j,i)
            hfx(j,i) = hfx(j,i) + sent(n,j,i)
            qfx(j,i) = qfx(j,i) + evpr(n,j,i)
            tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
            tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

            if ( lchem ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jj,i)
              sfracb2d(j,i) = sfracb2d(j,i)                         &
                              + d_one - (c2rfvegnosno(jj,i)+       &
                              c2rfracsno(jj,i))
              sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jj,i)
              ssw2da(j,i) = ssw2da(j,i)*landfrac(jj,i)             &
                            + (d_one-landfrac(jj,i))*ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i)*landfrac(jj,i)         &
                              + (d_one-landfrac(jj,i))*delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i)*landfrac(jj,i)         &
                              + (d_one-landfrac(jj,i))*delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i)*landfrac(jj,i)         &
                              + (d_one-landfrac(jj,i))*sigf(n,j,i)
              sfracb2d(j,i) = sfracb2d(j,i)*landfrac(jj,i)         &
                              + (d_one-landfrac(jj,i))*            &
                              (d_one-sigf(n,j,i))*(d_one-scvk(n,j,i))
              sfracs2d(j,i) = sfracs2d(j,i)*landfrac(jj,i)         &
                              + (d_one-landfrac(jj,i))             &
                              *(sigf(n,j,i)*wt(n,j,i)+(d_one-sigf(n,j,i)) &
                              *scvk(n,j,i))
            end if
 
            if ( iocnflx==1 .or.                                    &
                 (iocnflx==2 .and. ldmsk1(n,j,i) /= 0) ) then
              tgbb(j,i) = tgbb(j,i) + &
                        ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four+ &
                            lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
            else
              tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
            end if
          end do
 
          uvdrag(j,i) = uvdrag(j,i)* &
                        (d_one-landfrac(jj,i))+ &
                        c2ruvdrag(jj,i)*landfrac(jj,i)
          hfx(j,i) = hfx(j,i)*(d_one-landfrac(jj,i)) + &
                           c2rsenht(jj,i)*landfrac(jj,i)
          qfx(j,i) = qfx(j,i)*(d_one-landfrac(jj,i)) + &
                           c2rlatht(jj,i)*landfrac(jj,i)
          tground2(j,i) = tground2(j,i)*(d_one-landfrac(jj,i)) +     &
                         c2rtgb(jj,i)*landfrac(jj,i)
          tgbb(j,i) = tgbb(j,i)* (d_one-landfrac(jj,i)) +   &
                            c2rtgbb(jj,i)*landfrac(jj,i)
          tground1(j,i) = tground1(j,i)*(d_one-landfrac(jj,i)) +     &
                         c2rtgb(jj,i)*landfrac(jj,i)
 
          do n = 1 , nnsg
            !abt added below for the landfraction method
            sncv(n,j,i) = c2rsnowc(jj,i)*landfrac(jj,i)          &
                           + sncv(n,j,i)*(d_one-landfrac(jj,i))
            tgrd(n,j,i) = c2rtgb(jj,i)*landfrac(jj,i) + tgrd(n,j,i) &
                          *(d_one-landfrac(jj,i))
            tgbrd(n,j,i) = c2rtgb(jj,i)*landfrac(jj,i)            &
                           + tgbrd(n,j,i)*(d_one-landfrac(jj,i))
            taf(n,j,i) = c2r2mt(jj,i)*landfrac(jj,i)            &
                           + t2m(n,j,i)*(d_one-landfrac(jj,i))
            !note taf is 2m temp not temp in foilage
            tlef(n,j,i) = c2rtlef(jj,i)*landfrac(jj,i)          &
                            + tlef(n,j,i)*(d_one-landfrac(jj,i))
            tsw(n,j,i) = c2rsmtot(jj,i)*landfrac(jj,i)          &
                           + tsw(n,j,i)*(d_one-landfrac(jj,i))
            rsw(n,j,i) = c2rsm1m(jj,i)*landfrac(jj,i)           &
                           + rsw(n,j,i)*(d_one-landfrac(jj,i))
            ssw(n,j,i) = c2rsm10cm(jj,i)*landfrac(jj,i)         &
                           + ssw(n,j,i)*(d_one-landfrac(jj,i))
            q2d(j,i) = c2r2mq(jj,i)*landfrac(jj,i) + q2m(n,j,i)  &
                       *(d_one-landfrac(jj,i))
 
            evpa(n,j,i) = evpa(n,j,i) + dtbat*qfx(j,i)
            sena(n,j,i) = sena(n,j,i) + dtbat*hfx(j,i)
            srfrna(n,j,i) = c2rro_sur(jj,i)*dtbat
            runoff(n,j,i) = c2rro_sub(jj,i)*dtbat + c2rro_sur(jj,i)*dtbat
          end do
          prca(j,i) = prca(j,i) + dtbat*pptc(j,i)
          prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
        end if
      end do
    end do
 
      ! Fill output arrays if needed
 
    if ( mod(ktau+1,kbats) == 0 .or. ktau == 0 ) then
 
      do i = ici1 , ici2
        do j = jci1 , jci2

          jj = (jxp*myid) + j

          u10m_o(j,i) = 0.0
          v10m_o(j,i) = 0.0
          tg_o(j,i) = 0.0
          t2m_o(j,i) = 0.0
 
          do n = 1 , nnsg
            if ( ldmsk1(n,j,i) /= 0 ) then
              u10m_s(n,j,i) = real(uatm(j,i,kz))
              v10m_s(n,j,i) = real(vatm(j,i,kz))
              tg_s(n,j,i) = real(tgrd(n,j,i))
              t2m_s(n,j,i) = real(taf(n,j,i))
              u10m_o(j,i) = u10m_o(j,i) + real(uatm(j,i,kz))
              v10m_o(j,i) = v10m_o(j,i) + real(vatm(j,i,kz))
              t2m_o(j,i) = t2m_o(j,i) + real(taf(n,j,i))
              tg_o(j,i) = tg_o(j,i) + real(tgrd(n,j,i))
            else if ( ldmsk1(n,j,i) == 0 ) then
              tg_s(n,j,i) = real(tgrd(n,j,i))
              u10m_s(n,j,i) = real(u10m(n,j,i))
              v10m_s(n,j,i) = real(v10m(n,j,i))
              t2m_s(n,j,i) = real(t2m(n,j,i))
              u10m_o(j,i) = u10m_o(j,i) + real(u10m(n,j,i))
              v10m_o(j,i) = v10m_o(j,i) + real(v10m(n,j,i))
              t2m_o(j,i) = t2m_o(j,i) + real(t2m(n,j,i))
              tg_o(j,i) = tg_o(j,i) + real(tgrd(n,j,i))
            end if
          end do
 
          tgmx_o(j,i) = amax1(tgmx_o(j,i),tg_o(j,i))
          tgmn_o(j,i) = amin1(tgmn_o(j,i),tg_o(j,i))
          t2mx_o(j,i) = amax1(t2mx_o(j,i),t2m_o(j,i))
          t2mn_o(j,i) = amin1(t2mn_o(j,i),t2m_o(j,i))
          w10x_o(j,i) = amax1(w10x_o(j,i), &
                   sqrt(u10m_o(j,i)**2.0+v10m_o(j,i)**2.0))
          real_4 = real((pptnc(j,i)+pptc(j,i)))
          pcpx_o(j,i) = amax1(pcpx_o(j,i),real_4)
          pcpa_o(j,i) = pcpa_o(j,i) + real_4/fdaysrf
          tavg_o(j,i) = tavg_o(j,i)+t2m_o(j,i)/fdaysrf
          real_4 = real((sfps(j,i)+ptop)*d_10)
          psmn_o(j,i) = amin1(psmn_o(j,i),real_4)
          pptnc(j,i) = d_zero
          pptc(j,i) = d_zero
          if ( fsw(j,i) > 120.0D0 ) then
            sund_o(j,i) = sund_o(j,i) + real(dtbat)
            sunt_o(j,i) = sunt_o(j,i) + real(dtbat)
          end if
 
          drag_o(j,i) = 0.0
          q2m_o(j,i) = 0.0
          evpa_o(j,i) = 0.0
          sena_o(j,i) = 0.0
          do n = 1 , nnsg
            if ( ldmsk1(n,j,i) /= 0 ) then
              q2m_s(n,j,i) = real(q2d(j,i))
              drag_s(n,j,i) = real(uvdrag(j,i))
              evpa_s(n,j,i) = real(evpa(n,j,i)*mmpd)
              sena_s(n,j,i) = real(sena(n,j,i)*wpm2)
              tpr_s(n,j,i) = real((prnca(j,i)+prca(j,i))*mmpd)
              prcv_s(n,j,i) = real(prca(j,i)*mmpd)
              ps_s(n,j,i) = real(sfcp(n,j,i)*0.01D0)
 
              q2m_o(j,i) = q2m_o(j,i) + real(q2d(j,i))
              drag_o(j,i) = drag_o(j,i) + real(uvdrag(j,i))
              evpa_o(j,i) = evpa_o(j,i) + real(evpa(n,j,i))
              sena_o(j,i) = sena_o(j,i) + real(sena(n,j,i))
            else if ( ldmsk1(n,j,i) == 0 ) then
              q2m_s(n,j,i) = real(q2m(n,j,i))
              drag_s(n,j,i) = real(drag(n,j,i))
              evpa_s(n,j,i) = real(evpa(n,j,i)*mmpd)
              sena_s(n,j,i) = real(sena(n,j,i)*wpm2)
              tpr_s(n,j,i) = real((prnca(j,i)+prca(j,i))*mmpd)
              prcv_s(n,j,i) = real(prca(j,i)*mmpd)
              ps_s(n,j,i) = real(sfcp(n,j,i)*0.01D0)
 
              q2m_o(j,i) = q2m_o(j,i) + real(q2m(n,j,i))
              drag_o(j,i) = drag_o(j,i) + real(drag(n,j,i))
              evpa_o(j,i) = evpa_o(j,i) + real(evpa(n,j,i))
              sena_o(j,i) = sena_o(j,i) + real(sena(n,j,i))
            end if
          end do
          tpr_o(j,i) = real((prnca(j,i)+prca(j,i))*mmpd)
          evpa_o(j,i) = evpa_o(j,i)*real(mmpd)
          sena_o(j,i) = sena_o(j,i)*real(wpm2)
          flwa_o(j,i) = real(flwa(j,i)*wpm2)
          fswa_o(j,i) = real(fswa(j,i)*wpm2)
          flwd_o(j,i) = real(flwda(j,i)*wpm2)
          sina_o(j,i) = real(sina(j,i)*wpm2)
          prcv_o(j,i) = real(prca(j,i)*mmpd)
          ps_o(j,i) = real((sfps(j,i)+ptop)*d_10)
          zpbl_o(j,i) = real(hpbl(j,i))
 
          tlef_o(j,i) = 0.0
          ssw_o(j,i) = 0.0
          rsw_o(j,i) = 0.0
          rnos_o(j,i) = 0.0
          scv_o(j,i) = 0.0
          nnn = 0
          do n = 1 , nnsg
            if ( ldmsk1(n,j,i) /= 0 .and. landmask(jj,i)/=3 ) then
              tlef_o(j,i) = tlef_o(j,i) + real(c2rtlef(jj,i))
              ssw_o(j,i) = ssw_o(j,i) + real(c2rsm10cm(jj,i))
              rsw_o(j,i) = rsw_o(j,i) + real(c2rsm1m(jj,i))
              ! Correct unit of measure of runoff coming from CLM
              rnos_o(j,i) = rnos_o(j,i) + real(srfrna(n,j,i)*d_r1000)
              scv_o(j,i) = scv_o(j,i) + real(c2rsnowc(jj,i))
              tlef_s(n,j,i) = real(c2rtlef(jj,i))
              ssw_s(n,j,i) = real(c2rsm10cm(jj,i))
              rsw_s(n,j,i) = real(c2rsm1m(jj,i))
              rnos_s(n,j,i) = real(srfrna(n,j,i)*mmpd)
              scv_s(n,j,i) = real(c2rsnowc(jj,i))
              nnn = nnn + 1
            else
              tlef_s(n,j,i) = smissval
              ssw_s(n,j,i) = smissval
              rsw_s(n,j,i) = smissval
              rnos_s(n,j,i) = smissval
              scv_s(n,j,i) = smissval
            end if
          end do
          if ( nnn>=max0(nnsg/2,1) ) then
            tlef_o(j,i) = tlef_o(j,i)/real(nnn)
            ssw_o(j,i) = ssw_o(j,i)/real(nnn)
            rsw_o(j,i) = rsw_o(j,i)/real(nnn)
            rnos_o(j,i) = rnos_o(j,i)/real(nnn)*real(mmpd)
            scv_o(j,i) = scv_o(j,i)/real(nnn)
          else
            tlef_o(j,i) = smissval
            ssw_o(j,i) = smissval
            rsw_o(j,i) = smissval
            rnos_o(j,i) = smissval
            scv_o(j,i) = smissval
          end if
          ! reset accumulation arrays to zero
          do n = 1 , nnsg
            evpa(n,j,i) = d_zero
            srfrna(n,j,i) = d_zero
            sena(n,j,i) = d_zero
          end do
          prnca(j,i) = d_zero
          prca(j,i) = d_zero
          flwa(j,i) = d_zero
          flwda(j,i) = d_zero
          fswa(j,i) = d_zero
          svga(j,i) = d_zero
          sina(j,i) = d_zero
 
        end do
      end do
    end if
  end if  ! end if ivers = 2
 
  end subroutine interfclm
!
  subroutine fill_frame2d(a,b)
    implicit none
    real(dp) , pointer , intent(in) , dimension(:,:) :: a
    real(dp) , pointer , intent(out) , dimension(:,:) :: b
    b(jci1:jci2,ici1:ici2) = a(jci1:jci2,ici1:ici2)
    if ( ma%has_bdyleft ) then
      b(jce1,ici1:ici2) = a(jci1,ici1:ici2)
    end if
    if ( ma%has_bdyright ) then
      b(jce2,ici1:ici2) = a(jci2,ici1:ici2)
      b(jde2,ici1:ici2) = a(jci2,ici1:ici2)
    end if
    if ( ma%has_bdybottom ) then
      b(jci1:jci2,ice1) = a(jci1:jci2,ici1)
    end if
    if ( ma%has_bdytop ) then
      b(jci1:jci2,ice2) = a(jci1:jci2,ici2)
      b(jci1:jci2,ide2) = a(jci1:jci2,ici2)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdybottom ) then
      b(jce1,ice1) = a(jci1,ici1)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdytop ) then
      b(jce1,ice2) = a(jci1,ici2)
      b(jce1,ide2) = a(jci1,ici2)
    end if
    if ( ma%has_bdyright .and. ma%has_bdybottom ) then
      b(jce2,ice1) = a(jci2,ici1)
      b(jde2,ice1) = a(jci2,ici1)
      b(jce2,ide1) = a(jci2,ici1)
      b(jde2,ide1) = a(jci2,ici1)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      b(jce2,ice2) = a(jci2,ici2)
      b(jde2,ice2) = a(jci2,ici2)
      b(jce2,ide2) = a(jci2,ici2)
      b(jde2,ide2) = a(jci2,ici2)
    end if
  end subroutine fill_frame2d

  subroutine fill_frame3d(a,b)
    implicit none
    real(dp) , pointer , intent(in) , dimension(:,:,:) :: a
    real(dp) , pointer , intent(out) , dimension(:,:) :: b
    b(jci1:jci2,ici1:ici2) = a(jci1:jci2,ici1:ici2,kz)
    if ( ma%has_bdyleft ) then
      b(jce1,ici1:ici2) = a(jci1,ici1:ici2,kz)
    end if
    if ( ma%has_bdyright ) then
      b(jce2,ici1:ici2) = a(jci2,ici1:ici2,kz)
      b(jde2,ici1:ici2) = a(jci2,ici1:ici2,kz)
    end if
    if ( ma%has_bdybottom ) then
      b(jci1:jci2,ice1) = a(jci1:jci2,ici1,kz)
    end if
    if ( ma%has_bdytop ) then
      b(jci1:jci2,ice2) = a(jci1:jci2,ici2,kz)
      b(jci1:jci2,ide2) = a(jci1:jci2,ici2,kz)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdybottom ) then
      b(jce1,ice1) = a(jci1,ici1,kz)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdytop ) then
      b(jce1,ice2) = a(jci1,ici2,kz)
      b(jce1,ide2) = a(jci1,ici2,kz)
    end if
    if ( ma%has_bdyright .and. ma%has_bdybottom ) then
      b(jce2,ice1) = a(jci2,ici1,kz)
      b(jce2,ide1) = a(jci2,ici1,kz)
      b(jde2,ice1) = a(jci2,ici1,kz)
      b(jde2,ide1) = a(jci2,ici1,kz)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      b(jce2,ice2) = a(jci2,ici2,kz)
      b(jce2,ide2) = a(jci2,ici2,kz)
      b(jde2,ice2) = a(jci2,ici2,kz)
      b(jde2,ide2) = a(jci2,ici2,kz)
    end if
  end subroutine fill_frame3d

  subroutine fill_frame4d(a,b,l)
    implicit none
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: a
    integer , intent(in) :: l
    real(dp) , pointer , intent(out) , dimension(:,:) :: b
    b(jci1:jci2,ici1:ici2) = a(jci1:jci2,ici1:ici2,kz,l)
    if ( ma%has_bdyleft ) then
      b(jce1,ici1:ici2) = a(jci1,ici1:ici2,kz,l)
    end if
    if ( ma%has_bdyright ) then
      b(jce2,ici1:ici2) = a(jci2,ici1:ici2,kz,l)
      b(jde2,ici1:ici2) = a(jci2,ici1:ici2,kz,l)
    end if
    if ( ma%has_bdybottom ) then
      b(jci1:jci2,ice1) = a(jci1:jci2,ici1,kz,l)
    end if
    if ( ma%has_bdytop ) then
      b(jci1:jci2,ice2) = a(jci1:jci2,ici2,kz,l)
      b(jci1:jci2,ide2) = a(jci1:jci2,ici2,kz,l)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdybottom ) then
      b(jce1,ice1) = a(jci1,ici1,kz,l)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdytop ) then
      b(jce1,ice2) = a(jci1,ici2,kz,l)
      b(jce1,ide2) = a(jci1,ici2,kz,l)
    end if
    if ( ma%has_bdyright .and. ma%has_bdybottom ) then
      b(jce2,ice1) = a(jci2,ici1,kz,l)
      b(jce2,ide1) = a(jci2,ici1,kz,l)
      b(jde2,ice1) = a(jci2,ici1,kz,l)
      b(jde2,ide1) = a(jci2,ici1,kz,l)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      b(jce2,ice2) = a(jci2,ici2,kz,l)
      b(jce2,ide2) = a(jci2,ici2,kz,l)
      b(jde2,ice2) = a(jci2,ici2,kz,l)
      b(jde2,ide2) = a(jci2,ici2,kz,l)
    end if
  end subroutine fill_frame4d

  subroutine solar_clm(idatex,calday,declin,xyear)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idatex
    integer , intent(in)  :: xyear
    real(dp) , intent(out) :: calday , declin
    real(dp) :: decdeg
    real(dp) :: mvelp , obliq
    integer :: iyear_ad
    logical :: log_print
    character (len=64) :: subroutine_name='solar_clm'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    calday = yeardayfrac(idatex)
!
    log_print = .false.
    iyear_ad = xyear
!   Get eccen,obliq,mvelp,obliqr,lambm0,mvelpp
    call shr_orb_params(iyear_ad,r2ceccen,obliq,mvelp,r2cobliqr,      &
                        r2clambm0,r2cmvelpp,log_print)
!
!   Get declin,eccf
    call shr_orb_decl(calday,r2ceccen,r2cmvelpp,r2clambm0,r2cobliqr,  &
                      declin,r2ceccf)

!   convert declin to degrees
    decdeg = declin/degrad
    write (aline, 99001) calday, decdeg
    call say

    call time_end(subroutine_name,idindx)

99001 format (11x,'*** Day ',f12.4,' solar declination angle = ',f12.8,&
        &   ' degrees.')
!
  end subroutine solar_clm

  subroutine zenit_clm(coszrs)
    implicit none
    real(dp) , pointer , intent(out), dimension(:,:) :: coszrs
!
    integer :: i , j
    real(dp) :: cldy , declinp1 , xxlon
    real(dp) :: xxlat
    integer :: idindx=0
    character (len=64) :: subroutine_name='zenitm_clm'
!
    call time_begin(subroutine_name,idindx)
    cldy = get_curr_calday()
    call shr_orb_decl(cldy,r2ceccen,r2cmvelpp,r2clambm0, &
                      r2cobliqr,declinp1,r2ceccf)
    do i = ici1 , ici2
      do j = jci1 , jci2
        xxlat = xlat(j,i)*degrad
        xxlon = xlon(j,i)*degrad
        coszrs(j,i) = shr_orb_cosz(cldy,xxlat,xxlon,declinp1)
        coszrs(j,i) = dmax1(0.0D0,coszrs(j,i))
        coszrs(j,i) = dmin1(1.0D0,coszrs(j,i))
      end do
    end do
    call time_end(subroutine_name,idindx)
  end subroutine zenit_clm

end module mod_mtrxclm

#endif
