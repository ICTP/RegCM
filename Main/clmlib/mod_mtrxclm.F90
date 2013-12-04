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
  use mod_runparams , only : idate0 , iqv , solcon , clmfrq , &
                             imask , ilawrence_albedo , ichem , ksrf , xmonth
  use mod_mpmessage
  use mod_service
  use mod_mppparam
  use mod_date
  use mod_clm
  use mod_bats_common
  use mod_bats_mtrxbats
  use mod_bats_drag
  use mod_bats_zengocn
  use mod_outvars

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
    integer(ik8) , intent(in) :: ktau
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mtrxclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    call interfclm(1,ktau)
    call rcmdrv()
    call clm_run1(r2cdoalb,r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call clm_run2(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
    call interfclm(2,ktau)
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
  subroutine initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf,igases,iaeros,ctracers)
!
    use initializeMod
    use shr_orb_mod
    use clm_varpar,    only : lsmlon , lsmlat
    use clm_varsur,    only : landmask , landfrac , satbrt_clm
    use clm_varsur,    only : r2cimask , init_tgb , init_snow , r2coutfrq
    use clm_varsur,    only : clm2bats_veg , ht_rcm
    use clm_varsur,    only : clm_fracveg , r2cilawrence_albedo
    use clm_varsur,    only : cgaschem, caerosol
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
!
    logical , intent(in)                   :: ifrest
    type(rcm_time_and_date) , intent(in)   :: idate1 , idate2
    real(rk8) , intent(in)                 :: dtrad , dtsrf , dx
    integer(ik4) , intent(in), optional         :: igases
    integer(ik4) , intent(in), optional         :: iaeros
    character(len=6), intent(in), optional :: ctracers(*) 
!
    integer(ik4) :: i , j , ig , jg , n
    integer(ik4) :: year , month , day , hour
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Initialize run control variables for clm
    !
    r2comm = mycomm
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
    ! Lawrence modifications
    r2cilawrence_albedo = ilawrence_albedo

    !chemistry fields
    cgaschem = igases
    caerosol = iaeros

    ! Set elevation and BATS landuse type (abt added)
    allocate(ht_rcm(jx,iy))
    allocate(init_tgb(jx,iy))
    allocate(init_snow(jx,iy))
    allocate(satbrt_clm(jx,iy))
    allocate(clm_fracveg(jx,iy))
    allocate(clm2bats_veg(jx,iy))
    clm_fracveg(:,:) = d_zero

#if (defined VOC)
    voc_em(:,:) = d_zero
    voc_em1(:,:) = d_zero
    voc_em2(:,:) = d_zero
#endif
    if ( igases == 1 ) then
      dep_vels(:,:,:) = d_zero
    end if

    call grid_fill(ht,ht_rcm)
    call grid_fill(lndcat,satbrt_clm)
    call grid_fill(tground1,init_tgb)
    call grid_fill(snowam,init_snow)

    !
    ! End of clm run control variable initialization
    !
    ! Assign regcm values to the passed variables
    ! regcm writes uneven # of j values to arrays. for now fix by
    ! copying neighboring values
    !
    if ( .not. ifrest ) then
      ! Rainfall
      cprate(:,:)  = d_zero
      ncprate(:,:) = d_zero
      ! Radiation
      sols2d(:,:)  = d_zero
      soll2d(:,:)  = d_zero
      solsd2d(:,:) = d_zero
      solld2d(:,:) = d_zero
      dwrlwf(:,:)  = d_zero
      ! Albedo
      ! Set initial albedos to clm dry soil values for mid-colored soils
      swdiralb(:,:) = 0.16D0
      swdifalb(:,:) = 0.16D0
      lwdiralb(:,:) = 0.32D0
      lwdifalb(:,:) = 0.32D0
    end if

    call fill_frame(xlat,r2cxlatd)
    call fill_frame(xlon,r2cxlond)

    r2cxlat = r2cxlatd*degrad
    r2cxlon = r2cxlond*degrad

    if ( .not.ifrest ) then
      !
      !    Gather values of each nproc work_in array and fill
      !      work_out(jx*iy) array.
      !    3. Copy 1d work_out array to 2d (jx,iy) array for passing
      !    to clm.
      !
      call fill_frame(tatm,r2ctb)
      call fill_frame(qvatm,r2cqb)
      r2cqb = r2cqb/(d_one+r2cqb)
      call fill_frame(hgt,r2czga)
      call fill_frame(uatm,r2cuxb)
      call fill_frame(vatm,r2cvxb)
      call fill_frame(sfps,r2cpsb)
      r2cpsb = (r2cpsb+ptop)*d_1000
      call fill_frame(cprate,r2crnc)
      call fill_frame(ncprate,r2crnnc)
      call fill_frame(sols2d,r2csols)
      call fill_frame(soll2d,r2csoll)
      call fill_frame(solsd2d,r2csolsd)
      call fill_frame(solld2d,r2csolld)
      call fill_frame(dwrlwf,r2cflwd)

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
    if ( igases == 1 ) call seq_drydep_init(ntr, ctracers)
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
    if ( .not. ifrest ) then
      call rcmdrv()
    end if !end ifrest test

    ! Correct landmask
    do i = 1 , iy
      do j = 1 , jx
        if ( dabs(landfrac(j,i)-d_one) >= 0.1D0 .and. &
             dabs(landfrac(j,i))       >= 0.1D0 ) then
          landmask(j,i) = 3
        end if
      end do
    end do

    ! Initialize ldmsk1 now that clm has determined the land sea mask
    ! Initialize accumulation variables at zero

    if ( .not. ifrest ) then
      do n = 1 , nnsg
        do i = ici1 , ici2
          ig = global_cross_istart+i-1
          do j = jci1 , jci2
            jg = global_cross_jstart+j-1
            ldmsk1(n,j,i) = landmask(jg,ig)
            tgbrd(n,j,i)  = tground2(j,i)
            taf(n,j,i)    = tground2(j,i)
            tlef(n,j,i)   = tground2(j,i)
            ldew(n,j,i)   = d_zero
            snag(n,j,i)   = d_zero
            sncv(n,j,i)   = dmax1(sncv(n,j,i),d_zero)
            sfice(n,j,i)  = d_zero
            gwet(n,j,i)   = d_half
          end do
        end do
      end do
      do i = ici1 , ici2
        ig = global_cross_istart+i-1
        do j = jci1 , jci2
          jg = global_cross_jstart+j-1
          rswf(j,i)    = d_zero
          rlwf(j,i)    = d_zero
          vegswab(j,i) = d_zero
          !
          ! Set some clm land surface/vegetation variables to the ones
          ! used in RegCM.  Make sure all are consistent
          !
          lndcat(j,i) = clm2bats_veg(jg,ig)
          if ( clm2bats_veg(jg,ig) < 0.1D0 ) lndcat(j,i) = 15.0D0
          do n = 1 , nnsg
            lndcat1(n,j,i) = clm2bats_veg(jg,ig)
            if ( clm2bats_veg(jg,ig) < 0.1D0 ) lndcat1(n,j,i) = 15.0D0
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
    deallocate(ht_rcm)
    deallocate(init_tgb)
    deallocate(init_snow)
    deallocate(satbrt_clm)
    deallocate(clm2bats_veg)
    deallocate(clm_fracveg)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initclm
!
  subroutine albedoclm
    use clm_varsur , only : landfrac
    implicit none
    integer(ik4) :: i , j , ig , jg
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'albedoclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    !
    ! Calculate DEFAULT albedo to be the BATS one
    !
    call albedobats
    !
    ! Section Below added for albedo to be corrected by CLM
    ! calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
    ! to be equal to 1 which can cause a FPE.  To avoid this
    ! use albedo calculated with BATS method when albedo=1
    !
    do i = ici1 , ici2
      ig = global_cross_istart+i-1
      do j = jci1 , jci2
        jg = global_cross_jstart+j-1
        if ( ( ldmsk(j,i) == 0 .and. landfrac(jg,ig) > d_zero ) ) then
          !
          ! Correct "Mixed points" from CLM for their water fraction
          ! in the albedo (good when < 1)
          !
          if ( (d_one-c2ralbdirs(jg,ig)) > dlowval ) then
            swdiralb(j,i) = swdiralb(j,i)*(d_one-landfrac(jg,ig)) + &
                          c2ralbdirs(j,i)*landfrac(jg,ig)
          end if
          if ( (d_one-c2ralbdirl(jg,ig)) > dlowval ) then
            lwdiralb(j,i) = lwdiralb(j,i)*(d_one-landfrac(jg,ig)) + &
                          c2ralbdirl(j,i)*landfrac(jg,ig)
          end if
          if ( (d_one-c2ralbdifs(jg,ig)) > dlowval ) then 
            swdifalb(j,i) = swdifalb(j,i)*(d_one-landfrac(jg,ig)) + &
                          c2ralbdifs(jg,ig)*landfrac(jg,ig)
          end if
          if ( (d_one-c2ralbdifl(jg,ig)) > dlowval ) then
            lwdifalb(j,i) = lwdifalb(j,i)*(d_one-landfrac(jg,ig)) + &
                          c2ralbdifl(jg,ig)
          end if
        else if (ldmsk(j,i) /= 0 ) then
          !
          ! Use over land CLM calculated albedo (good when < 1)
          !
          if ( (d_one-c2ralbdirs(jg,ig)) > dlowval ) then
            swdiralb(j,i) = c2ralbdirs(jg,ig)
          end if
          if ( (d_one-c2ralbdirl(jg,ig)) > dlowval ) then
            lwdiralb(j,i) = c2ralbdirl(jg,ig)
          end if
          if ( (d_one-c2ralbdifs(jg,ig)) > dlowval ) then 
            swdifalb(j,i) = c2ralbdifs(jg,ig)
          end if
          if ( (d_one-c2ralbdifl(jg,ig)) > dlowval ) then
            lwdifalb(j,i) = c2ralbdifl(jg,ig)
          end if
        end if
        swalb(j,i) = swdiralb(j,i)+swdifalb(j,i)
        lwalb(j,i) = lwdiralb(j,i)+lwdifalb(j,i)
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine albedoclm
!
  subroutine interfclm(ivers,ktau)
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
    integer(ik8) , intent(in) :: ktau
!
    integer(ik4) :: i , j , ic , jc , ib , jg , ig , kk , n
    integer(ik4) :: idep , icpu , nout
    real(rk8) :: fact
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interfclm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    if ( ivers == 1 ) then

      call fill_frame(tatm,r2ctb)
      call fill_frame(qvatm,r2cqb)
      r2cqb = r2cqb/(d_one+r2cqb)
      call fill_frame(hgt,r2czga)
      call fill_frame(uatm,r2cuxb)
      call fill_frame(vatm,r2cvxb)
      call fill_frame(sfps,r2cpsb)
      r2cpsb = (r2cpsb+ptop)*d_1000
      call fill_frame(cprate,r2crnc)
      call fill_frame(ncprate,r2crnnc)
      call fill_frame(sols2d,r2csols)
      call fill_frame(soll2d,r2csoll)
      call fill_frame(solsd2d,r2csolsd)
      call fill_frame(solld2d,r2csolld)
      call fill_frame(dwrlwf,r2cflwd)

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
        do ib = 1 , c2rngc(icpu)
          kk = c2rngc(icpu)
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
              !**** Dry deposition velocities from CLM4
              do n = 1 , ntr
                dep_vels(j,i,n) = c2r_depout(ib+(kk*(n-1))+idep)
              end do
              c2rfracsno(j,i)   = c2r_allout(ib+(20*kk)+ic)
              c2rfvegnosno(j,i) = c2r_allout(ib+(21*kk)+ic)
#if (defined VOC)
              voc_em(j,i)       = c2r_allout(ib+(22*kk)+ic)
              voc_em1(j,i)      = c2r_allout(ib+(23*kk)+ic)
              voc_em2(j,i)      = c2r_allout(ib+(24*kk)+ic)
#endif
            else if ( cgaschem == 1 .and. caerosol /= 1 ) then
              !**** Dry deposition velocities from CLM4
              do n = 1 , ntr
                dep_vels(j,i,n) = c2r_depout(ib+(kk*(n-1))+idep)
              end do
#if (defined VOC)
              voc_em(j,i)    = c2r_allout(ib+(20*kk)+ic)
              voc_em1(j,i)   = c2r_allout(ib+(21*kk)+ic)
              voc_em2(j,i)   = c2r_allout(ib+(22*kk)+ic)
#endif
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

      call interf(1,ktau)

      if ( iocnflx == 2 ) then
        call zengocndrv
      else if ( iocnflx == 1 ) then
        call dragc
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( ldmsk(j,i) == 0 ) then
              tgrd(:,j,i) = tground2(j,i)
              drag(:,j,i) = cdrx(:,j,i)*vspda(:,j,i)*rhs(:,j,i)
              tlef(:,j,i) = sts(:,j,i)
              delq(:,j,i) =  qs(:,j,i) - qgrd(:,j,i)
              delt(:,j,i) = sts(:,j,i) - tgrd(:,j,i)
              evpr(:,j,i) = -drag(:,j,i)*delq(:,j,i)
              sent(:,j,i) = -drag(:,j,i)*cpd*delt(:,j,i)
              fact = z10fra(1,j,i)/zlgocn(1,j,i)
              u10m(:,j,i) = usw(j,i)*(d_one-fact)
              v10m(:,j,i) = vsw(j,i)*(d_one-fact)
              fact = z2fra(1,j,i)/zlgocn(1,j,i)
              t2m(:,j,i) = sts(:,j,i) - delt(:,j,i)*fact
              q2m(:,j,i) = qs(:,j,i) - delq(:,j,i)*fact
            end if
          end do
        end do
      end if

      do i = ici1 , ici2
        ig = global_cross_istart+i-1
        do j = jci1 , jci2
          jg = global_cross_jstart+j-1
          uvdrag(j,i)   = d_zero
          hfx(j,i)      = d_zero
          qfx(j,i)      = d_zero
          tground2(j,i) = d_zero
          tground1(j,i) = d_zero
          tgbb(j,i)     = d_zero

          if ( ichem == 1 ) then
            ssw2da(j,i) = d_zero
            sdeltk2d(j,i) = d_zero
            sdelqk2d(j,i) = d_zero
            sfracv2d(j,i) = d_zero
            sfracb2d(j,i) = d_zero
            sfracs2d(j,i) = d_zero
          end if

          if ( landmask(jg,ig) == 1 ) then
            tground2(j,i) = c2rtgb(jg,ig)
            tground1(j,i) = c2rtgb(jg,ig)
            hfx(j,i)      = c2rsenht(jg,ig)
            qfx(j,i)      = c2rlatht(jg,ig)
            uvdrag(j,i)   = c2ruvdrag(jg,ig)
            tgbb(j,i)     = c2rtgbb(jg,ig)

            do n = 1 , nnsg
              tgrd(n,j,i)   = c2rtgb(jg,ig)
              tgbrd(n,j,i)  = c2rtgb(jg,ig)
              ! supposed to be lower soil layer temp not tgrnd
              taf(n,j,i)    = c2r2mt(jg,ig)
              t2m(n,j,i)    = c2r2mt(jg,ig)
              u10m(n,j,i)   = uatm(j,i)
              v10m(n,j,i)   = vatm(j,i)
              tlef(n,j,i)   = c2rtlef(jg,ig)
              tsw(n,j,i)    = c2rsmtot(jg,ig)
              rsw(n,j,i)    = c2rsm1m(jg,ig)
              ssw(n,j,i)    = c2rsm10cm(jg,ig)
              sncv(n,j,i)   = c2rsnowc(jg,ig)
              srnof(n,j,i)  = c2rro_sur(jg,ig)*dtbat
              trnof(n,j,i)  = (c2rro_sub(jg,ig)+c2rro_sur(jg,ig))*dtbat
              q2m(n,j,i)    = c2r2mq(jg,ig)

              if ( ichem == 1 ) then
                ssw2da(j,i)   = ssw2da(j,i) + ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jg,ig)
                sfracb2d(j,i) = sfracb2d(j,i) + d_one -               &
                               (c2rfvegnosno(jg,ig)+c2rfracsno(jg,ig))
                sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jg,ig)
              end if
            end do
          else if ( landmask(jg,ig) == 0 ) then !ocean
            do n = 1 , nnsg
              uvdrag(j,i)   = uvdrag(j,i) + drag(n,j,i)
              hfx(j,i)      = hfx(j,i) + sent(n,j,i)
              qfx(j,i)      = qfx(j,i) + evpr(n,j,i)
              tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
              tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

              if ( ichem == 1  ) then
                ssw2da(j,i)   = ssw2da(j,i) + ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i) + sigf(n,j,i)
                sfracb2d(j,i) = sfracb2d(j,i) + (d_one-sigf(n,j,i))    &
                                *(d_one-scvk(n,j,i))
                sfracs2d(j,i) = sfracs2d(j,i) + sigf(n,j,i)*wt(n,j,i) &
                                + (d_one-sigf(n,j,i))*scvk(n,j,i)
              end if

              if ( ldmsk1(n,j,i) /= 0 ) then
                tgbb(j,i) = tgbb(j,i) + &
                     ((d_one-lncl(n,j,i))*tgrd(n,j,i)**4 +   &
                     lncl(n,j,i)*tlef(n,j,i)**4)**d_rfour
              else
                tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
              end if
              ssw(n,j,i)   = dmissval
              rsw(n,j,i)   = dmissval
              tsw(n,j,i)   = dmissval
              trnof(n,j,i) = dmissval
              srnof(n,j,i) = dmissval
              sncv(n,j,i)  = dmissval
            end do

            do n = 1 , nnsg
              taf(n,j,i)  = t2m(n,j,i)
              sncv(n,j,i) = sncv(n,j,i)
            end do
          else if ( landmask(jg,ig) == 3 ) then
            !gridcell with some % land and ocean
            do n = 1 , nnsg
              uvdrag(j,i)   = uvdrag(j,i) + drag(n,j,i)
              hfx(j,i)      = hfx(j,i) + sent(n,j,i)
              qfx(j,i)      = qfx(j,i) + evpr(n,j,i)
              tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
              tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

              if ( ichem == 1 ) then
                ssw2da(j,i)   = ssw2da(j,i) + ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jg,ig)
                sfracb2d(j,i) = sfracb2d(j,i) + d_one - &
                               (c2rfvegnosno(jg,ig) + c2rfracsno(jg,ig))
                sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jg,ig)
                ssw2da(j,i) = ssw2da(j,i)*landfrac(jg,ig)             &
                              + (d_one-landfrac(jg,ig))*ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i)*landfrac(jg,ig)         &
                                + (d_one-landfrac(jg,ig))*delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i)*landfrac(jg,ig)         &
                                + (d_one-landfrac(jg,ig))*delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i)*landfrac(jg,ig)         &
                                + (d_one-landfrac(jg,ig))*sigf(n,j,i)
                sfracb2d(j,i) = sfracb2d(j,i)*landfrac(jg,ig)         &
                                + (d_one-landfrac(jg,ig))*            &
                                (d_one-sigf(n,j,i))*(d_one-scvk(n,j,i))
                sfracs2d(j,i) = sfracs2d(j,i)*landfrac(jg,ig)         &
                                + (d_one-landfrac(jg,ig))             &
                                *(sigf(n,j,i)*wt(n,j,i)+(d_one-sigf(n,j,i)) &
                                *scvk(n,j,i))
              end if

              if ( ldmsk1(n,j,i) /= 0 ) then
                tgbb(j,i) = tgbb(j,i) + &
                          ((d_one-lncl(n,j,i))*tgrd(n,j,i)**4+ &
                              lncl(n,j,i)*tlef(n,j,i)**4)**d_rfour
              else
                tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
              end if
            end do

            uvdrag(j,i)   = uvdrag(j,i)*(d_one-landfrac(jg,ig)) + &
                            c2ruvdrag(jg,ig)*landfrac(jg,ig)
            hfx(j,i)      = hfx(j,i)*(d_one-landfrac(jg,ig)) + &
                            c2rsenht(jg,ig)*landfrac(jg,ig)
            qfx(j,i)      = qfx(j,i)*(d_one-landfrac(jg,ig)) + &
                            c2rlatht(jg,ig)*landfrac(jg,ig)
            tground2(j,i) = tground2(j,i)*(d_one-landfrac(jg,ig)) + &
                            c2rtgb(jg,ig)*landfrac(jg,ig)
            tgbb(j,i)     = tgbb(j,i)* (d_one-landfrac(jg,ig)) +   &
                            c2rtgbb(jg,ig)*landfrac(jg,ig)
            tground1(j,i) = tground1(j,i)*(d_one-landfrac(jg,ig)) + &
                            c2rtgb(jg,ig)*landfrac(jg,ig)

            do n = 1 , nnsg
              !abt added below for the landfraction method
              sncv(n,j,i)   = c2rsnowc(jg,ig)*landfrac(jg,ig)          &
                             + sncv(n,j,i)*(d_one-landfrac(jg,ig))
              tgrd(n,j,i)   = c2rtgb(jg,ig)*landfrac(jg,ig) + tgrd(n,j,i) &
                             *(d_one-landfrac(jg,ig))
              tgbrd(n,j,i)  = c2rtgb(jg,ig)*landfrac(jg,ig)            &
                             + tgbrd(n,j,i)*(d_one-landfrac(jg,ig))
              taf(n,j,i)    = c2r2mt(jg,ig)*landfrac(jg,ig)            &
                             + t2m(n,j,i)*(d_one-landfrac(jg,ig))
              !note taf is 2m temp not temp in foilage
              tlef(n,j,i)   = c2rtlef(jg,ig)*landfrac(jg,ig)          &
                             + tlef(n,j,i)*(d_one-landfrac(jg,ig))
              tsw(n,j,i)    = c2rsmtot(jg,ig)*landfrac(jg,ig)          &
                             + tsw(n,j,i)*(d_one-landfrac(jg,ig))
              rsw(n,j,i)    = c2rsm1m(jg,ig)*landfrac(jg,ig)           &
                             + rsw(n,j,i)*(d_one-landfrac(jg,ig))
              ssw(n,j,i)    = c2rsm10cm(jg,ig)*landfrac(jg,ig)         &
                             + ssw(n,j,i)*(d_one-landfrac(jg,ig))
              q2m(n,j,i)    = c2r2mq(jg,ig)*landfrac(jg,ig) + q2m(n,j,i)  &
                             *(d_one-landfrac(jg,ig))
              srnof(n,j,i)  = c2rro_sur(jg,ig)*dtbat
              trnof(n,j,i) = c2rro_sub(jg,ig)*dtbat + c2rro_sur(jg,ig)*dtbat
            end do
          end if
        end do
      end do

      ! Those are needed for output purposes

      ! Fill accumulators

      if ( ktau > 0 ) then
        if ( ifatm ) then
          if ( associated(atm_tgb_out) ) &
            atm_tgb_out = atm_tgb_out + sum(tgbrd,1)*rdnnsg
          if ( associated(atm_tsw_out) ) &
            atm_tsw_out = atm_tsw_out + sum(tsw,1)*rdnnsg
        end if
        if ( ifsrf ) then
          if ( associated(srf_evp_out) ) &
            srf_evp_out = srf_evp_out + sum(evpr,1)*rdnnsg
          if ( associated(srf_tpr_out) ) &
            srf_tpr_out = srf_tpr_out + totpr
          if ( associated(srf_prcv_out) ) &
            srf_prcv_out = srf_prcv_out + cprate
          ! Reset accumulation from precip and cumulus
          ncprate = d_zero
          cprate  = d_zero
          if ( associated(srf_zpbl_out) ) &
            srf_zpbl_out = srf_zpbl_out + hpbl
          if ( associated(srf_scv_out) ) &
            srf_scv_out = srf_scv_out + sum(sncv,1)*rdnnsg
          if ( associated(srf_sund_out) ) then
            where( rswf > 120.0D0 )
              srf_sund_out = srf_sund_out + dtbat
            end where
          end if
          if ( associated(srf_runoff_out) ) then
            srf_runoff_out(:,:,1) = srf_runoff_out(:,:,1) + sum(srnof,1)*rdnnsg
            srf_runoff_out(:,:,2) = srf_runoff_out(:,:,2) + sum(trnof,1)*rdnnsg
          end if
          if ( associated(srf_sena_out) ) &
            srf_sena_out = srf_sena_out + sum(sent,1)*rdnnsg
          if ( associated(srf_flw_out) ) &
            srf_flw_out = srf_flw_out + rlwf
          if ( associated(srf_fsw_out) ) &
            srf_fsw_out = srf_fsw_out + rswf
          if ( associated(srf_fld_out) ) &
            srf_fld_out = srf_fld_out + dwrlwf
          if ( associated(srf_sina_out) ) &
            srf_sina_out = srf_sina_out + sinc
        end if
        if ( ifsub ) then
          call reorder_add_subgrid(sfcp,sub_ps_out)
          if ( associated(sub_evp_out) ) &
            call reorder_add_subgrid(evpr,sub_evp_out)
          if ( associated(sub_scv_out) ) &
            call reorder_add_subgrid(sncv,sub_scv_out,mask=ldmsk1)
          if ( associated(sub_sena_out) ) &
            call reorder_add_subgrid(sent,sub_sena_out)
          if ( associated(sub_runoff_out) ) then
            call reorder_add_subgrid(srnof,sub_runoff_out,1,ldmsk1)
            call reorder_add_subgrid(trnof,sub_runoff_out,2,ldmsk1)
          end if
        end if
        if ( ifsts ) then
          if ( associated(sts_tgmax_out) ) &
            sts_tgmax_out = max(sts_tgmax_out,sum(tgrd,1)*rdnnsg)
          if ( associated(sts_tgmin_out) ) &
            sts_tgmin_out = min(sts_tgmin_out,sum(tgrd,1)*rdnnsg)
          if ( associated(sts_t2max_out) ) &
            sts_t2max_out(:,:,1) = max(sts_t2max_out(:,:,1),sum(t2m,1)*rdnnsg)
          if ( associated(sts_t2min_out) ) &
            sts_t2min_out(:,:,1) = min(sts_t2min_out(:,:,1),sum(t2m,1)*rdnnsg)
          if ( associated(sts_t2min_out) ) &
            sts_t2avg_out(:,:,1) = sts_t2avg_out(:,:,1) + sum(t2m,1)*rdnnsg
          if ( associated(sts_w10max_out) ) &
            sts_w10max_out(:,:,1) = max(sts_w10max_out(:,:,1), &
              sqrt(sum((u10m**2+v10m**2),1)*rdnnsg))
          if ( associated(sts_pcpmax_out) ) &
            sts_pcpmax_out = max(sts_pcpmax_out,totpr)
          if ( associated(sts_pcpavg_out) ) &
            sts_pcpavg_out = sts_pcpavg_out + totpr
          if ( associated(sts_psmin_out) ) &
            sts_psmin_out = min(sts_psmin_out, &
              (sfps(jci1:jci2,ici1:ici2)+ptop)*d_10)
          if ( associated(sts_sund_out) ) then
            where( rswf > 120.0D0 )
              sts_sund_out = sts_sund_out + dtbat
            end where
          end if
        end if
      end if

      ! Those are for the output, but collected only at POINT in time

      if ( mod(ktau+1,ksrf) == 0 ) then

        call albedoclm

        if ( ifsrf ) then
          if ( associated(srf_uvdrag_out) ) &
            srf_uvdrag_out = uvdrag
          if ( associated(srf_tg_out) ) &
            srf_tg_out = tground1(jci1:jci2,ici1:ici2)
          if ( associated(srf_tlef_out) ) then
            where ( sum(ldmsk1,1) > nnsg/2 )
              srf_tlef_out = sum(tlef,1)*rdnnsg
            elsewhere
              srf_tlef_out = dmissval
            end where
          end if
          if ( associated(srf_aldirs_out) ) &
            srf_aldirs_out = swdiralb
          if ( associated(srf_aldifs_out) ) &
            srf_aldifs_out = swdifalb
          if ( associated(srf_seaice_out) ) &
            srf_seaice_out = sum(sfice,1)*rdnnsg*d_r1000
          if ( associated(srf_t2m_out) ) &
            srf_t2m_out(:,:,1) = sum(t2m,1)*rdnnsg
          if ( associated(srf_q2m_out) ) &
            srf_q2m_out(:,:,1) = sum(q2m,1)*rdnnsg
          if ( associated(srf_u10m_out) ) &
            srf_u10m_out(:,:,1) = sum(u10m,1)*rdnnsg
          if ( associated(srf_v10m_out) ) &
            srf_v10m_out(:,:,1) = sum(v10m,1)*rdnnsg
          if ( associated(srf_smw_out) ) then
            srf_smw_out(:,:,1) = sum(ssw,1)*rdnnsg
            srf_smw_out(:,:,2) = sum(rsw,1)*rdnnsg
          end if
        end if

        if ( ifsub ) then
          if ( associated(sub_uvdrag_out) ) &
            call reorder_subgrid(drag,sub_uvdrag_out)
          if ( associated(sub_tg_out) ) &
            call reorder_subgrid(tgrd,sub_tg_out)
          if ( associated(sub_tlef_out) ) &
            call reorder_subgrid(tlef,sub_tlef_out,mask=ldmsk1)
          if ( associated(sub_u10m_out) ) &
            call reorder_subgrid(u10m,sub_u10m_out)
          if ( associated(sub_v10m_out) ) &
            call reorder_subgrid(v10m,sub_v10m_out)
          if ( associated(sub_t2m_out) ) &
            call reorder_subgrid(t2m,sub_t2m_out)
          if ( associated(sub_q2m_out) ) &
            call reorder_subgrid(q2m,sub_q2m_out)
          if ( associated(sub_smw_out) ) then
            call reorder_subgrid(ssw,sub_smw_out,1,ldmsk1)
            call reorder_subgrid(rsw,sub_smw_out,2,ldmsk1)
          end if
        end if

      end if ! IF output time

    end if  ! end if ivers = 2
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine interfclm
!
  subroutine fill_frame2d(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
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
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
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
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
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

  subroutine solar_clm(idatex,calday,declin,xyear)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idatex
    integer(ik4) , intent(in)  :: xyear
    real(rk8) , intent(out) :: calday , declin
    real(rk8) :: decdeg
    real(rk8) :: mvelp , obliq
    integer(ik4) :: iyear_ad
    logical :: log_print
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'solar_clm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    calday = yeardayfrac(idatex)
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
    if ( myid == italk ) then
      write (stdout,'(a,f12.2,a,f12.8,a)') ' JDay ', calday , &
        ' solar declination angle = ', decdeg , ' degrees'
      write(stdout, '(18x,a,f12.4,a)') ' solar TSI irradiance    = ' , &
        solcon, ' W/m^2'
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
!
  end subroutine solar_clm

  subroutine zenit_clm(coszrs)
    implicit none
    real(rk8) , pointer , intent(out), dimension(:,:) :: coszrs
!
    integer(ik4) :: i , j
    real(rk8) :: cldy , declinp1 , xxlon
    real(rk8) :: xxlat
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'zenitm_clm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine zenit_clm

end module mod_mtrxclm

#endif
