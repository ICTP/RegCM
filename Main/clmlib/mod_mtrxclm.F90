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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams , only : idate0 , iqv , solcon , clmfrq , iocnflx , &
      imask , ilawrence_albedo , ichem , rcmtimer , syncro_srf , &
      dx , dtsrf , dtrad , igaschem , iaerosol , chtrname , idate1 ,    &
      idate2 , eccen , obliqr , lambm0 , mvelpp
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_mppparam
  use mod_date
  use mod_clm
  use mod_regcm_types
  use mod_stdio
  use mod_sunorbit

  use clm_time_manager , only : get_curr_calday
  use shr_orb_mod , only : shr_orb_cosz , shr_orb_decl , &
                           shr_orb_params
  use clm_varsur,    only : clm_fracveg

  implicit none

  private

  public :: mtrxclm
  public :: initclm
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
  subroutine mtrxclm(lm,lms)
    use atmdrvMod , only : rcmdrv
    use clm_comp , only : clm_run1 , clm_run2
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mtrxclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    call interfclm(lm,lms,1)
    call rcmdrv()
    call clm_run1(r2cdoalb,eccen,obliqr,lambm0,mvelpp)
    call clm_run2(eccen,obliqr,lambm0,mvelpp)
    call interfclm(lm,lms,2)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
  subroutine initclm(lm,lms)
    use initializeMod
    use shr_orb_mod
    use clm_varpar,    only : lsmlon , lsmlat
    use clm_varsur,    only : landmask , landfrac , satbrt_clm
    use clm_varsur,    only : r2cimask , init_tgb , init_snow , r2coutfrq
    use clm_varsur,    only : clm2bats_veg , ht_rcm
    use clm_varsur,    only : r2cilawrence_albedo
    use clm_varsur,    only : cgaschem, caerosol, numdays
    use atmdrvMod
    use program_offMod
    use clm_comp
    use clmtype
    use perf_mod

#if (defined VOC)
    use clm_varpar ,   only : nvoc
#endif
    use clm_drydep,    only : seq_drydep_init
!
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
!
    integer(ik4) :: i , j , n
    integer(ik4) :: year , month , day , hour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Initialize run control variables for clm
    !
    numdays = dayspy
    r2comm = mycomm
    if ( rcmtimer%integrating( ) ) then
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
    r2cirad = idint(dtrad/r2cdtime)
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
    ! Lawrence modifications
    r2cilawrence_albedo = ilawrence_albedo

    !chemistry fields
    cgaschem = igaschem
    caerosol = iaerosol

    ! Set elevation and BATS landuse type (abt added)
    allocate(ht_rcm(jx,iy))
    allocate(init_tgb(jx,iy))
    allocate(init_snow(jx,iy))
    allocate(satbrt_clm(jx,iy))
    allocate(clm_fracveg(jx,iy))
    allocate(clm2bats_veg(jx,iy))
    clm_fracveg(:,:) = d_zero

    if ( igaschem == 1 ) then
#if (defined VOC)
      lm%voc_em0(:,:) = d_zero
      lm%voc_em1(:,:) = d_zero
      lm%voc_em2(:,:) = d_zero
#endif
      lm%dep_vels(:,:,:) = d_zero
    end if

    call grid_fill(lm%ht,ht_rcm)
    call grid_fill(lm%lndcat,satbrt_clm)
    call grid_fill(lm%snowam,init_snow)

    !
    ! End of clm run control variable initialization
    !
    ! Assign regcm values to the passed variables
    ! regcm writes uneven # of j values to arrays. for now fix by
    ! copying neighboring values
    !
    call fill_frame(lm%xlat,r2cxlatd)
    call fill_frame(lm%xlon,r2cxlond)

    r2cxlat = r2cxlatd*degrad
    r2cxlon = r2cxlond*degrad

    if ( rcmtimer%start( ) ) then
      !
      !    Gather values of each nproc work_in array and fill
      !      work_out(jx*iy) array.
      !    3. Copy 1d work_out array to 2d (jx,iy) array for passing
      !    to clm.
      !
      call fill_frame(lm%tatm,r2ctb)
      call fill_frame(lm%qvatm,r2cqb)
      r2cqb = r2cqb/(d_one+r2cqb)
      call fill_frame(lm%hgt,r2czga)
      call fill_frame(lm%uatm,r2cuxb)
      call fill_frame(lm%vatm,r2cvxb)
      call fill_frame(lm%sfps,r2cpsb)
      call fill_frame(lm%cprate,r2crnc)
      call fill_frame(lm%ncprate,r2crnnc)
      r2crnc = r2crnc * syncro_srf%rw
      r2crnnc = r2crnnc * syncro_srf%rw
      call fill_frame(lm%swdir,r2csols)
      call fill_frame(lm%lwdir,r2csoll)
      call fill_frame(lm%swdif,r2csolsd)
      call fill_frame(lm%lwdif,r2csolld)
      call fill_frame(lm%dwrlwf,r2cflwd)

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
    do i = 1 , iy
      do j = 1 , jx
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
    if ( igaschem == 1 ) call seq_drydep_init(ntr, chtrname)
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

    if ( myid == iocpu ) write (stdout,*) 'Attempting to make atmospheric grid'
    call rcmdrv_init()
    if ( myid == iocpu ) write (stdout,*) 'Successfully  make atmospheric grid'

    ! Initialize radiation and atmosphere variables
    if ( rcmtimer%start( ) ) then
      call rcmdrv()
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ichem == 1 ) then
            lm%svegfrac2d(j,i) = clm_fracveg(j,i)
          end if
          if ( rcmtimer%start( ) ) then
            ! Set initial albedos to clm dry soil values for mid-colored soils
            lms%swdiralb(n,j,i) = 0.16_rkx
            lms%swdifalb(n,j,i) = 0.16_rkx
            lms%lwdiralb(n,j,i) = 0.32_rkx
            lms%lwdifalb(n,j,i) = 0.32_rkx
            lms%swalb(n,j,i) = (lms%swdiralb(n,j,i)+lms%swdifalb(n,j,i)) * &
                clm_fracveg(j,i)
            lms%lwalb(n,j,i) = (lms%lwdiralb(n,j,i)+lms%lwdifalb(n,j,i)) * &
                clm_fracveg(j,i)
          end if
        end do
      end do
    end do

    ! Correct landmask
    do i = 1 , iy
      do j = 1 , jx
        if ( dabs(landfrac(j,i)-d_one) >= 0.1_rkx .and. &
             dabs(landfrac(j,i))       >= 0.1_rkx ) then
          landmask(j,i) = 3
        end if
      end do
    end do

    ! Initialize ldmsk1 now that clm has determined the land sea mask
    ! Initialize accumulation variables at zero

    if ( rcmtimer%start( ) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            if ( landmask(j,i) == 3 ) then
              lm%ldmsk1(n,j,i) = 5
            else
              lm%ldmsk1(n,j,i) = landmask(j,i)
            end if
            lms%tgbrd(n,j,i) = lm%tg(j,i)
            lms%taf(n,j,i)   = lm%tg(j,i)
            lms%tlef(n,j,i)  = lm%tg(j,i)
            lms%snag(n,j,i)  = d_zero
            lms%sncv(n,j,i)  = dmax1(lms%sncv(n,j,i),d_zero)
            lms%sfice(n,j,i) = d_zero
            lms%ldew(n,j,i)  = d_zero
            lms%gwet(n,j,i)  = d_half
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          !
          ! Set some clm land surface/vegetation variables to the ones
          ! used in RegCM.  Make sure all are consistent
          !
          lm%lndcat(j,i) = clm2bats_veg(j,i)
          if ( clm2bats_veg(j,i) < 0.1_rkx ) lm%lndcat(j,i) = 15.0_rkx
          do n = 1 , nnsg
            lm%lndcat1(n,j,i) = clm2bats_veg(j,i)
            if ( clm2bats_veg(j,i) < 0.1_rkx ) lm%lndcat1(n,j,i) = 15.0_rkx
          end do
        end do
      end do
    end if

    ! deallocate some variables used in CLM initialization only
    deallocate(ht_rcm)
    deallocate(init_tgb)
    deallocate(init_snow)
    deallocate(satbrt_clm)
    deallocate(clm2bats_veg)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initclm
!
  subroutine albedoclm(lm,lms)
    use clm_varsur , only : landfrac
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'albedoclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    !
    ! Section Below added for albedo to be corrected by CLM
    ! calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
    ! to be equal to 1 which can cause a FPE.  To avoid this
    ! use albedo calculated with BATS method when albedo=1
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          !
          ! Use over land CLM calculated albedo (good when < 1)
          !
          if ( (d_one-c2ralbdirs(j,i)) > dlowval ) then
            lms%swdiralb(n,j,i) = c2ralbdirs(j,i)
          end if
          if ( (d_one-c2ralbdirl(j,i)) > dlowval ) then
            lms%lwdiralb(n,j,i) = c2ralbdirl(j,i)
          end if
          if ( (d_one-c2ralbdifs(j,i)) > dlowval ) then
            lms%swdifalb(n,j,i) = c2ralbdifs(j,i)
          end if
          if ( (d_one-c2ralbdifl(j,i)) > dlowval ) then
            lms%lwdifalb(n,j,i) = c2ralbdifl(j,i)
          end if
          lms%swalb(n,j,i) = (lms%swdiralb(n,j,i)+lms%swdifalb(n,j,i)) * &
                  clm_fracveg(j,i)
          lms%lwalb(n,j,i) = (lms%lwdiralb(n,j,i)+lms%lwdifalb(n,j,i)) * &
                  clm_fracveg(j,i)
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine albedoclm
!
  subroutine interfclm(lm,lms,ivers)
    use clmtype
    use clm_varsur , only : landmask , landfrac
    use clm_varsur , only : c2r_allout , omap_i , omap_j
    use clm_varsur,  only : cgaschem, caerosol

#if (defined VOC)
    use clm_varpar,  only : nvoc
#endif
    use clm_drydep,  only : c2r_depout, n_drydep

    implicit none
    !
    ! ivers = 1 : regcm -> clm
    ! ivers = 2 : clm -> regcm
    !
    integer(ik4) , intent(in) :: ivers
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms

    integer(ik4) :: i , j , ic , jc , ib , kk , n
    integer(ik4) :: idep , icpu , nout
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interfclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    if ( ivers == 1 ) then

      call fill_frame(lm%tatm,r2ctb)
      call fill_frame(lm%qvatm,r2cqb)
      r2cqb = r2cqb/(d_one+r2cqb)
      call fill_frame(lm%hgt,r2czga)
      call fill_frame(lm%uatm,r2cuxb)
      call fill_frame(lm%vatm,r2cvxb)
      call fill_frame(lm%sfps,r2cpsb)
      call fill_frame(lm%cprate,r2crnc)
      call fill_frame(lm%ncprate,r2crnnc)
      r2crnc = r2crnc * syncro_srf%rw
      r2crnnc = r2crnnc * syncro_srf%rw
      call fill_frame(lm%swdir,r2csols)
      call fill_frame(lm%lwdir,r2csoll)
      call fill_frame(lm%swdif,r2csolsd)
      call fill_frame(lm%lwdif,r2csolld)
      call fill_frame(lm%dwrlwf,r2cflwd)

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

      ic   = 0
      jc   = 1
      idep = 0
      nout = 20
      if ( ichem == 1 ) then  !Aerosol and/or Chem schemes on
#if (defined VOC)
        if ( cgaschem == 1 ) nout = nout + 3
#endif
        if ( caerosol == 1 ) nout = nout + 2
      end if
      do icpu = 1 , nproc
        kk = c2rngc(icpu)
        do ib = 1 , kk
          j = omap_i(jc)
          i = omap_j(jc)
          c2rtgb(j,i)     = c2r_allout(ib+( 0*kk)+ic)
          c2rsnowc(j,i)   = c2r_allout(ib+( 1*kk)+ic)
          c2rsenht(j,i)   = c2r_allout(ib+( 2*kk)+ic)
          c2rlatht(j,i)   = c2r_allout(ib+( 3*kk)+ic)
          c2ruvdrag(j,i)  = c2r_allout(ib+( 4*kk)+ic)
          c2ralbdirs(j,i) = c2r_allout(ib+( 5*kk)+ic)
          c2ralbdirl(j,i) = c2r_allout(ib+( 6*kk)+ic)
          c2ralbdifs(j,i) = c2r_allout(ib+( 7*kk)+ic)
          c2ralbdifl(j,i) = c2r_allout(ib+( 8*kk)+ic)
          c2rtgbb(j,i)    = c2r_allout(ib+( 9*kk)+ic)
          c2r2mt(j,i)     = c2r_allout(ib+(10*kk)+ic)
          c2r2mq(j,i)     = c2r_allout(ib+(11*kk)+ic)
          c2ru10(j,i)     = c2r_allout(ib+(12*kk)+ic)
          c2rtlef(j,i)    = c2r_allout(ib+(13*kk)+ic)
          c2rsm10cm(j,i)  = c2r_allout(ib+(14*kk)+ic)
          c2rsm1m(j,i)    = c2r_allout(ib+(15*kk)+ic)
          c2rsmtot(j,i)   = c2r_allout(ib+(16*kk)+ic)
          c2rinfl(j,i)    = c2r_allout(ib+(17*kk)+ic)
          c2rro_sur(j,i)  = c2r_allout(ib+(18*kk)+ic)
          c2rro_sub(j,i)  = c2r_allout(ib+(19*kk)+ic)
          if ( ichem == 1 ) then
            if ( cgaschem == 1 .and. caerosol == 1 ) then
              c2rfracsno(j,i)   = c2r_allout(ib+(20*kk)+ic)
              c2rfvegnosno(j,i) = c2r_allout(ib+(21*kk)+ic)
              !**** Dry deposition velocities from CLM4
#if (defined VOC)
              lm%voc_em0(j,i)    = c2r_allout(ib+(22*kk)+ic)
              lm%voc_em1(j,i)    = c2r_allout(ib+(23*kk)+ic)
              lm%voc_em2(j,i)    = c2r_allout(ib+(24*kk)+ic)
#endif
              do n = 1 , ntr
                lm%dep_vels(j,i,n) = c2r_depout(ib+(kk*(nout-1))+idep)
              end do
            else if ( cgaschem == 1 .and. caerosol /= 1 ) then
              !**** Dry deposition velocities from CLM
#if (defined VOC)
              lm%voc_em0(j,i)  = c2r_allout(ib+(20*kk)+ic)
              lm%voc_em1(j,i)   = c2r_allout(ib+(21*kk)+ic)
              lm%voc_em2(j,i)   = c2r_allout(ib+(22*kk)+ic)
#endif
              do n = 1 , ntr
                lm%dep_vels(j,i,n) = c2r_depout(ib+(kk*(nout-1))+idep)
              end do
            else if ( cgaschem /= 1 .and. caerosol == 1 ) then
              c2rfracsno(j,i)   = c2r_allout(ib+(20*kk)+ic)
              c2rfvegnosno(j,i) = c2r_allout(ib+(21*kk)+ic)
            end if
          end if
          jc = jc + 1
        end do
        ic = ic + c2rngc(icpu)*nout
        if ( ichem == 1 .and. cgaschem == 1 ) then
          idep = idep + c2rngc(icpu)*n_drydep
        end if
      end do

      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( landmask(j,i) == 1 .or. landmask(j,i) == 3 ) then
            do n = 1 , nnsg
              lms%tgbb(n,j,i)   = c2rtgbb(j,i)
              lms%drag(n,j,i)   = c2ruvdrag(j,i)
              lms%prcp(n,j,i)   = r2crnc(j,i) + r2crnnc(j,i)
              lms%tgrd(n,j,i)   = c2rtgb(j,i)
              lms%tgbrd(n,j,i)  = c2rtgb(j,i)
              lms%evpr(n,j,i)   = c2rlatht(j,i)
              lms%sent(n,j,i)   = c2rsenht(j,i)
              lms%taf(n,j,i)    = c2r2mt(j,i)
              lms%t2m(n,j,i)    = c2r2mt(j,i)
              lms%u10m(n,j,i)   = lm%uatm(j,i)/dlog(lm%hgt(j,i)*d_r10)
              lms%v10m(n,j,i)   = lm%vatm(j,i)/dlog(lm%hgt(j,i)*d_r10)
              lms%sfcp(n,j,i)   = lm%sfps(j,i)
              lms%tlef(n,j,i)   = c2rtlef(j,i)
              lms%tsw(n,j,i)    = c2rsmtot(j,i)
              lms%rsw(n,j,i)    = c2rsm1m(j,i)
              lms%ssw(n,j,i)    = c2rsm10cm(j,i)
              lms%sncv(n,j,i)   = c2rsnowc(j,i)
              lms%srnof(n,j,i)  = c2rro_sur(j,i)
              lms%trnof(n,j,i)  = (c2rro_sub(j,i)+c2rro_sur(j,i))
              lms%q2m(n,j,i)    = c2r2mq(j,i)
              lms%deltat(n,j,i) = lms%tgbrd(n,j,i)-lm%tatm(j,i)
              lms%deltaq(n,j,i) = (lm%qvatm(j,i)/(d_one+lm%qvatm(j,i))) - &
                      c2r2mq(j,i)
            end do
          end if
        end do
      end do

    end if  ! end if ivers = 2
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine interfclm
!
  subroutine fill_frame2d(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
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
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
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
      b(jde2,ice1) = a(jci2,ici1,kz)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      b(jce2,ice2) = a(jci2,ici2,kz)
      b(jde2,ice2) = a(jci2,ici2,kz)
      b(jce2,ide2) = a(jci2,ici2,kz)
      b(jde2,ide2) = a(jci2,ici2,kz)
    end if
  end subroutine fill_frame3d

  subroutine fill_frame4d(a,b,l)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
    integer(ik4) , intent(in) :: l
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
      b(jde2,ice1) = a(jci2,ici1,kz,l)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      b(jce2,ice2) = a(jci2,ici2,kz,l)
      b(jde2,ice2) = a(jci2,ici2,kz,l)
      b(jce2,ide2) = a(jci2,ici2,kz,l)
      b(jde2,ide2) = a(jci2,ici2,kz,l)
    end if
  end subroutine fill_frame4d

end module mod_mtrxclm

#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
