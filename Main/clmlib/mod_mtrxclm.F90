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
  use mod_runparams , only : idate0 , iqv
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
  subroutine initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf,igases,iaeros,ctracers)
!
    use initializeMod
    use shr_orb_mod
    use clm_varpar,    only : lsmlon , lsmlat
    use clm_varsur,    only : landmask , landfrac , satbrt_clm
    use clm_varsur,    only : r2cimask , init_tgb , r2coutfrq
    use clm_varsur,    only : clm2bats_veg , ht_rcm
    use clm_varsur,    only : clm_fracveg
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
    real(8) , intent(in)                   :: dtrad , dtsrf , dx
    integer , intent(in), optional         :: igases
    integer , intent(in), optional         :: iaeros
    character(len=6), intent(in), optional :: ctracers(*) 
!
    integer :: i , j , ig , jg , n , mpierr
    integer :: year , month , day , hour
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

    !chemistry fields
    cgaschem = igases
    caerosol = iaeros

    ! Set elevation and BATS landuse type (abt added)
    if ( .not.allocated(ht_rcm) )       allocate(ht_rcm(jx,iy))
    if ( .not.allocated(init_tgb) )     allocate(init_tgb(jx,iy))
    if ( .not.allocated(satbrt_clm) )   allocate(satbrt_clm(jx,iy))
    if ( .not.allocated(clm_fracveg) )  allocate(clm_fracveg(jx,iy))
    if ( .not.allocated(clm2bats_veg) ) allocate(clm2bats_veg(jx,iy))
    clm_fracveg(:,:) = d_zero

#if (defined VOC)
    voc_em(:,:) = d_zero
#endif
    if ( igases == 1 ) then
      dep_vels(:,:,:) = d_zero
    end if

    if ( myid == iocpu ) then
      ! Broadcast of those in CLM code.
      do i = 1 , iy
        do j = 1 , jx
          ht_rcm(j,i)      = htf(j,i)
          satbrt_clm(j,i)  = lndcatf(j,i)
          init_tgb(j,i)    = tsf(j,i)
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
      flwd(:,:)    = d_zero
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
      !
      !    Gather values of each nproc work_in array and fill
      !      work_out(jx*iy) array.
      !    3. Copy 1d work_out array to 2d (jx,iy) array for passing
      !    to clm.
      !
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

    if ( myid == iocpu ) write (6,*) 'Attempting to make atmospheric grid'
    call rcmdrv_init()
    if ( myid == iocpu ) write (6,*) 'Successfully  make atmospheric grid'

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
          ig = global_istart+i-1
          do j = jci1 , jci2
            jg = global_jstart+j-1
            ldmsk1(n,j,i) = landmask(jg,ig)
            tgbrd(n,j,i)  = tground2(j,i)
            taf(n,j,i)    = tground2(j,i)
            tlef(n,j,i)   = tground2(j,i)
            ldew(n,j,i)   = d_zero
            snag(n,j,i)   = d_zero
            sncv(n,j,i)   = dmax1(sncv(n,j,i),d_zero)
            sfice(n,j,i)  = d_zero
            gwet(n,j,i)   = d_half
            sena(n,j,i)   = d_zero
            evpa(n,j,i)   = d_zero
            srfrna(n,j,i) = d_zero
            runoff(n,j,i) = d_zero
          end do
        end do
      end do
      do i = ici1 , ici2
        ig = global_istart+i-1
        do j = jci1 , jci2
          jg = global_jstart+j-1
          fsw(j,i)    = d_zero
          flw(j,i)    = d_zero
          sabveg(j,i) = d_zero
          fswa(j,i)   = d_zero
          flwa(j,i)   = d_zero
          prca(j,i)   = d_zero
          prnca(j,i)  = d_zero
          svga(j,i)   = d_zero
          sina(j,i)   = d_zero
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
    integer :: i , j , ig , jg
!
    call albedobats(imon)
!
!     ****** Section Below added for albedo to be corrected by CLM
!     ****** calculated albedo.  NOTE: for cosz<=0 CLM assigns albedo
!     ****** to be equal to 1 which can cause a FPE.  To avoid this
!     ****** use albedo calculated with BATS method when albedo=1
!
    do i = ici1 , ici2
      ig = global_istart+i-1
      do j = jci1 , jci2
        jg = global_jstart+j-1
        if (ldmsk1(1,j,i) /= 0 .and. (d_one-aldirs(j,i)) > 1.0D-10 ) then
          aldirs(j,i) = aldirs(j,i)*landfrac(jg,ig) + &
                        aldirs(j,i)*(d_one-landfrac(jg,ig))
          aldirl(j,i) = aldirl(j,i)*landfrac(jg,ig) + &
                        aldirl(j,i)*(d_one-landfrac(jg,ig))
          aldifs(j,i) = aldifs(j,i)*landfrac(jg,ig) + &
                        aldifs(j,i)*(d_one-landfrac(jg,ig))
          aldifl(j,i) = aldifl(j,i)*landfrac(jg,ig) + &
                        aldifl(j,i)*(d_one-landfrac(jg,ig))
          albvs(j,i)  = aldirs(j,i)*landfrac(jg,ig) + &
                        albvs(j,i) *(d_one-landfrac(jg,ig))
          albvl(j,i)  = aldirl(j,i)*landfrac(jg,ig) + &
                        albvl(j,i) *(d_one-landfrac(jg,ig))
        end if
        fbat(j,i,aldirs_o) = real(aldirs(j,i))
        fbat(j,i,aldifs_o) = real(aldifs(j,i))
      end do
    end do
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
    integer , intent(in) :: ivers
    integer(8) , intent(in) :: ktau
!
    real(8) :: mmpd , wpm2
    integer :: i , j , ic , jc , ib , jg , ig , kk , n , icpu , nnn , nout
    integer :: idep, iddep
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

      ic   = 0
      jc   = 1
      idep = 0
      nout = 20
      if ( ichem == 1 ) then  !Aerosol and/or Chem schemes on
#if (defined VOC)
        if ( cgaschem == 1 ) nout = nout + 1
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
              do iddep = 1 , ntr
                dep_vels(j,i,iddep) = c2r_depout(ib+(kk*(iddep-1))+idep)
              end do
              c2rfracsno(j,i)   = c2r_allout(ib+(20*kk)+ic)
              c2rfvegnosno(j,i) = c2r_allout(ib+(21*kk)+ic)
#if (defined VOC)
              voc_em(j,i)       = c2r_allout(ib+(22*kk)+ic)
#endif
            else if ( cgaschem == 1 .and. caerosol /= 1 ) then
              !**** Dry deposition velocities from CLM4
              do iddep = 1 , ntr
                dep_vels(j,i,iddep) = c2r_depout(ib+(kk*(iddep-1))+idep)
              end do
#if (defined VOC)
              voc_em(j,i)    = c2r_allout(ib+(20*kk)+ic)
#endif
            else if ( cgaschem /= 1 .and. caerosol == 1 ) then
              c2rfracsno(j,i)   = c2r_allout(ib+(20*kk)+ic)
              c2rfvegnosno(j,i) = c2r_allout(ib+(21*kk)+ic)
            end if
          end if
          jc = jc + 1
        end do
        ic = ic + c2rngc(icpu)*nout
        if ( cgaschem == 1 ) idep = idep + c2rngc(icpu)*n_drydep
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
        ig = global_istart+i-1
        do j = jci1 , jci2
          jg = global_jstart+j-1
          uvdrag(j,i)   = d_zero
          hfx(j,i)      = d_zero
          qfx(j,i)      = d_zero
          tground2(j,i) = d_zero
          tground1(j,i) = d_zero
          tgbb(j,i)     = d_zero

          if ( lchem ) then
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

            aldirs(j,i)   = c2ralbdirs(jg,ig)
            aldirl(j,i)   = c2ralbdirl(jg,ig)
            aldifs(j,i)   = c2ralbdifs(jg,ig)
            aldifl(j,i)   = c2ralbdifl(jg,ig)

            do n = 1 , nnsg
              tgrd(n,j,i)   = c2rtgb(jg,ig)
              tgbrd(n,j,i)  = c2rtgb(jg,ig)
              ! supposed to be lower soil layer temp not tgrnd
              taf(n,j,i)    = c2r2mt(jg,ig)
              tlef(n,j,i)   = c2rtlef(jg,ig)
              tsw(n,j,i)    = c2rsmtot(jg,ig)
              rsw(n,j,i)    = c2rsm1m(jg,ig)
              ssw(n,j,i)    = c2rsm10cm(jg,ig)
              ldew(n,j,i)   = ldew(n,j,i)
              sncv(n,j,i)   = c2rsnowc(jg,ig)
              evpa(n,j,i)   = evpa(n,j,i) + dtbat*qfx(j,i)
              sena(n,j,i)   = sena(n,j,i) + dtbat*hfx(j,i)
              srfrna(n,j,i) = c2rro_sur(jg,ig)*dtbat
              runoff(n,j,i) = (c2rro_sub(jg,ig)+c2rro_sur(jg,ig))*dtbat

              if ( lchem ) then
                ssw2da(j,i)   = ssw2da(j,i) + ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i) + c2rfvegnosno(jg,ig)
                sfracb2d(j,i) = sfracb2d(j,i) + d_one -               &
                               (c2rfvegnosno(jg,ig)+c2rfracsno(jg,ig))
                sfracs2d(j,i) = sfracs2d(j,i) + c2rfracsno(jg,ig)
              end if
            end do
            !abt added for 2m humidity when landmask = 1 or 3
            q2d(j,i)   = c2r2mq(jg,ig)
            prca(j,i)  = prca(j,i) + dtbat*pptc(j,i)
            prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
            flwa(j,i)  = flwa(j,i) + dtbat*flw(j,i)
            flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
            fswa(j,i)  = fswa(j,i) + dtbat*fsw(j,i)
            svga(j,i)  = svga(j,i) + dtbat*sabveg(j,i)
            sina(j,i)  = sina(j,i) + dtbat*sinc(j,i)
          else if ( landmask(jg,ig) == 0 ) then !ocean
            do n = 1 , nnsg
              uvdrag(j,i)   = uvdrag(j,i) + drag(n,j,i)
              hfx(j,i)      = hfx(j,i) + sent(n,j,i)
              qfx(j,i)      = qfx(j,i) + evpr(n,j,i)
              tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
              tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

              if ( lchem  ) then
                ssw2da(j,i)   = ssw2da(j,i) + ssw(n,j,i)
                sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
                sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
                sfracv2d(j,i) = sfracv2d(j,i) + sigf(n,j,i)
                sfracb2d(j,i) = sfracb2d(j,i) + (d_one-sigf(n,j,i))    &
                                *(d_one-scvk(n,j,i))
                sfracs2d(j,i) = sfracs2d(j,i) + sigf(n,j,i)*wt(n,j,i) &
                                + (d_one-sigf(n,j,i))*scvk(n,j,i)
              end if

              if ( iocnflx == 1 .or.                                    &
                  (iocnflx == 2 .and. ldmsk1(n,j,i) /= 0) ) then
                tgbb(j,i) = tgbb(j,i) + &
                     ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four +   &
                     lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
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
              evpa(n,j,i) = evpa(n,j,i) + dtbat*evpr(n,j,i)
              sena(n,j,i) = sena(n,j,i) + dtbat*sent(n,j,i)
              if ( dabs(trnof(n,j,i)) > 1.0D-10 ) then
                srfrna(n,j,i) = srfrna(n,j,i) + trnof(n,j,i)/secpd*dtbat
              end if
              if ( dabs(srnof(n,j,i)) > 1.0D-10 .and. &
                   dabs(trnof(n,j,i)) > 1.0D-10 ) then
                runoff(n,j,i) = runoff(n,j,i) + &
                        (trnof(n,j,i)-srnof(n,j,i))/secpd*dtbat
              end if
            end do
            prca(j,i)  = prca(j,i) + dtbat*pptc(j,i)
            prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
            flwa(j,i)  = flwa(j,i) + dtbat*flw(j,i)
            flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
            fswa(j,i)  = fswa(j,i) + dtbat*fsw(j,i)
            svga(j,i)  = svga(j,i) + dtbat*sabveg(j,i)
            sina(j,i)  = sina(j,i) + dtbat*sinc(j,i)
          else if ( landmask(jg,ig) == 3 ) then
            !gridcell with some % land and ocean
            do n = 1 , nnsg
              uvdrag(j,i)   = uvdrag(j,i) + drag(n,j,i)
              hfx(j,i)      = hfx(j,i) + sent(n,j,i)
              qfx(j,i)      = qfx(j,i) + evpr(n,j,i)
              tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
              tground1(j,i) = tground1(j,i) + tgrd(n,j,i)

              if ( lchem ) then
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

              if ( iocnflx == 1 .or.                                    &
                  (iocnflx == 2 .and. ldmsk1(n,j,i) /= 0) ) then
                tgbb(j,i) = tgbb(j,i) + &
                          ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four+ &
                              lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
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
              q2d(j,i)      = c2r2mq(jg,ig)*landfrac(jg,ig) + q2m(n,j,i)  &
                             *(d_one-landfrac(jg,ig))
              evpa(n,j,i)   = evpa(n,j,i) + dtbat*qfx(j,i)
              sena(n,j,i)   = sena(n,j,i) + dtbat*hfx(j,i)
              srfrna(n,j,i) = c2rro_sur(jg,ig)*dtbat
              runoff(n,j,i) = c2rro_sub(jg,ig)*dtbat + c2rro_sur(jg,ig)*dtbat
            end do
            prca(j,i)  = prca(j,i) + dtbat*pptc(j,i)
            prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
            flwa(j,i)  = flwa(j,i) + dtbat*flw(j,i)
            flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
            fswa(j,i)  = fswa(j,i) + dtbat*fsw(j,i)
            svga(j,i)  = svga(j,i) + dtbat*sabveg(j,i)
            sina(j,i)  = sina(j,i) + dtbat*sinc(j,i)
          end if
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          fbat(j,i,u10m_o) = 0.0
          fbat(j,i,v10m_o) = 0.0
          fbat(j,i,tg_o)   = 0.0
          fbat(j,i,t2m_o)  = 0.0
          fbat(j,i,q2m_o)  = 0.0
          do n = 1 , nnsg
            if ( ldmsk1(n,j,i) /= 0 ) then
              fsub(n,j,i,q2m_s)  = real(q2d(j,i))
              fsub(n,j,i,u10m_s) = real(uatm(j,i,kz))
              fsub(n,j,i,v10m_s) = real(vatm(j,i,kz))
              fsub(n,j,i,tg_s)   = real(tgrd(n,j,i))
              fsub(n,j,i,t2m_s)  = real(taf(n,j,i))
              fbat(j,i,q2m_o)    = fbat(j,i,q2m_o) + real(q2d(j,i))
              fbat(j,i,u10m_o)   = fbat(j,i,u10m_o) + real(uatm(j,i,kz))
              fbat(j,i,v10m_o)   = fbat(j,i,v10m_o) + real(vatm(j,i,kz))
              fbat(j,i,t2m_o)    = fbat(j,i,t2m_o) + real(taf(n,j,i))
              fbat(j,i,tg_o)     = fbat(j,i,tg_o) + real(tgrd(n,j,i))
            else if ( ldmsk1(n,j,i) == 0 ) then
              fsub(n,j,i,q2m_s)  = real(q2m(n,j,i))
              fsub(n,j,i,tg_s)   = real(tgrd(n,j,i))
              fsub(n,j,i,u10m_s) = real(u10m(n,j,i))
              fsub(n,j,i,v10m_s) = real(v10m(n,j,i))
              fsub(n,j,i,t2m_s)  = real(t2m(n,j,i))
              fbat(j,i,q2m_o)    = fbat(j,i,q2m_o) + real(q2m(n,j,i))
              fbat(j,i,u10m_o)   = fbat(j,i,u10m_o) + real(u10m(n,j,i))
              fbat(j,i,v10m_o)   = fbat(j,i,v10m_o) + real(v10m(n,j,i))
              fbat(j,i,t2m_o)    = fbat(j,i,t2m_o) + real(t2m(n,j,i))
              fbat(j,i,tg_o)     = fbat(j,i,tg_o) + real(tgrd(n,j,i))
            end if
          end do
          fbat(j,i,tgmx_o) = amax1(fbat(j,i,tgmx_o),fbat(j,i,tg_o))
          fbat(j,i,tgmn_o) = amin1(fbat(j,i,tgmn_o),fbat(j,i,tg_o))
          fbat(j,i,t2mx_o) = amax1(fbat(j,i,t2mx_o),fbat(j,i,t2m_o))
          fbat(j,i,t2mn_o) = amin1(fbat(j,i,t2mn_o),fbat(j,i,t2m_o))
          fbat(j,i,w10x_o) = amax1(fbat(j,i,w10x_o), &
                        sqrt(fbat(j,i,u10m_o)**2.0+fbat(j,i,v10m_o)**2.0))
          real_4      = real((pptnc(j,i)+pptc(j,i)))
          fbat(j,i,pcpx_o) = amax1(fbat(j,i,pcpx_o),real_4)
          fbat(j,i,pcpa_o) = fbat(j,i,pcpa_o) + real_4/fdaysrf
          fbat(j,i,tavg_o) = fbat(j,i,tavg_o)+fbat(j,i,t2m_o)/fdaysrf
          real_4      = real((sfps(j,i)+ptop)*d_10)
          fbat(j,i,psmn_o) = amin1(fbat(j,i,psmn_o),real_4)
          if ( fsw(j,i) > 120.0D0 ) then
            fbat(j,i,sund_o) = fbat(j,i,sund_o) + real(dtbat)
            fbat(j,i,sunt_o) = fbat(j,i,sunt_o) + real(dtbat)
          end if
          pptnc(j,i) = d_zero
          pptc(j,i)  = d_zero
        end do
      end do

      ! Fill output arrays if needed

      if ( mod(ktau+1,kbats) == 0 .or. ktau == 0 ) then

        do i = ici1 , ici2
          ig = global_istart+i-1
          do j = jci1 , jci2
            jg = global_jstart+j-1
            fbat(j,i,drag_o) = 0.0
            fbat(j,i,evpa_o) = 0.0
            fbat(j,i,sena_o) = 0.0
            do n = 1 , nnsg
              if ( ldmsk1(n,j,i) /= 0 ) then
                fsub(n,j,i,drag_s) = real(uvdrag(j,i))
                fsub(n,j,i,evpa_s) = real(evpa(n,j,i)*mmpd)
                fsub(n,j,i,sena_s) = real(sena(n,j,i)*wpm2)
                fsub(n,j,i,tpr_s)  = real((prnca(j,i)+prca(j,i))*mmpd)
                fsub(n,j,i,prcv_s) = real(prca(j,i)*mmpd)
                fsub(n,j,i,ps_s)   = real(sfcp(n,j,i)*0.01D0)
                fbat(j,i,drag_o)   = fbat(j,i,drag_o) + real(uvdrag(j,i))
                fbat(j,i,evpa_o)   = fbat(j,i,evpa_o) + real(evpa(n,j,i))
                fbat(j,i,sena_o)   = fbat(j,i,sena_o) + real(sena(n,j,i))
              else if ( ldmsk1(n,j,i) == 0 ) then
                fsub(n,j,i,drag_s) = real(drag(n,j,i))
                fsub(n,j,i,evpa_s) = real(evpa(n,j,i)*mmpd)
                fsub(n,j,i,sena_s) = real(sena(n,j,i)*wpm2)
                fsub(n,j,i,tpr_s)  = real((prnca(j,i)+prca(j,i))*mmpd)
                fsub(n,j,i,prcv_s) = real(prca(j,i)*mmpd)
                fsub(n,j,i,ps_s)   = real(sfcp(n,j,i)*0.01D0)
                fbat(j,i,drag_o)   = fbat(j,i,drag_o) + real(drag(n,j,i))
                fbat(j,i,evpa_o)   = fbat(j,i,evpa_o) + real(evpa(n,j,i))
                fbat(j,i,sena_o)   = fbat(j,i,sena_o) + real(sena(n,j,i))
              end if
            end do
            fbat(j,i,tpr_o)  = real((prnca(j,i)+prca(j,i))*mmpd)
            fbat(j,i,evpa_o) = fbat(j,i,evpa_o)*real(mmpd)
            fbat(j,i,sena_o) = fbat(j,i,sena_o)*real(wpm2)
            fbat(j,i,flwa_o) = real(flwa(j,i)*wpm2)
            fbat(j,i,fswa_o) = real(fswa(j,i)*wpm2)
            fbat(j,i,flwd_o) = real(flwda(j,i)*wpm2)
            fbat(j,i,sina_o) = real(sina(j,i)*wpm2)
            fbat(j,i,prcv_o) = real(prca(j,i)*mmpd)
            fbat(j,i,ps_o) = real((sfps(j,i)+ptop)*d_10)
            fbat(j,i,zpbl_o) = real(hpbl(j,i))
            fbat(j,i,tlef_o) = 0.0
            fbat(j,i,ssw_o)  = 0.0
            fbat(j,i,rsw_o)  = 0.0
            fbat(j,i,rnos_o) = 0.0
            fbat(j,i,scv_o)  = 0.0
            nnn = 0
            do n = 1 , nnsg
              if ( ldmsk1(n,j,i) /= 0 .and. landmask(jg,ig) /= 3 ) then
                fbat(j,i,tlef_o)   = fbat(j,i,tlef_o) + real(c2rtlef(jg,ig))
                fbat(j,i,ssw_o)    = fbat(j,i,ssw_o) + real(c2rsm10cm(jg,ig))
                fbat(j,i,rsw_o)    = fbat(j,i,rsw_o) + real(c2rsm1m(jg,ig))
                ! Correct unit of measure of runoff coming from CLM
                fbat(j,i,rnos_o)   = fbat(j,i,rnos_o) + &
                                     real(srfrna(n,j,i)*d_r1000)
                fbat(j,i,scv_o)    = fbat(j,i,scv_o) + real(c2rsnowc(jg,ig))
                fsub(n,j,i,tlef_s) = real(c2rtlef(jg,ig))
                fsub(n,j,i,ssw_s)  = real(c2rsm10cm(jg,ig))
                fsub(n,j,i,rsw_s)  = real(c2rsm1m(jg,ig))
                fsub(n,j,i,rnos_s) = real(srfrna(n,j,i)*mmpd)
                fsub(n,j,i,scv_s)  = real(c2rsnowc(jg,ig))
                nnn = nnn + 1
              else
                fsub(n,j,i,tlef_s) = smissval
                fsub(n,j,i,ssw_s)  = smissval
                fsub(n,j,i,rsw_s)  = smissval
                fsub(n,j,i,rnos_s) = smissval
                fsub(n,j,i,scv_s)  = smissval
              end if
            end do
            if ( nnn >= max0(nnsg/2,1) ) then
              fbat(j,i,tlef_o) = fbat(j,i,tlef_o)/real(nnn)
              fbat(j,i,ssw_o)  = fbat(j,i,ssw_o)/real(nnn)
              fbat(j,i,rsw_o)  = fbat(j,i,rsw_o)/real(nnn)
              fbat(j,i,rnos_o) = fbat(j,i,rnos_o)/real(nnn)*real(mmpd)
              fbat(j,i,scv_o)  = fbat(j,i,scv_o)/real(nnn)
            else
              fbat(j,i,tlef_o) = smissval
              fbat(j,i,ssw_o)  = smissval
              fbat(j,i,rsw_o)  = smissval
              fbat(j,i,rnos_o) = smissval
              fbat(j,i,scv_o)  = smissval
            end if
            ! reset accumulation arrays to zero
            do n = 1 , nnsg
              evpa(n,j,i)   = d_zero
              srfrna(n,j,i) = d_zero
              sena(n,j,i)   = d_zero
            end do
            prnca(j,i) = d_zero
            prca(j,i)  = d_zero
            flwa(j,i)  = d_zero
            flwda(j,i) = d_zero
            fswa(j,i)  = d_zero
            svga(j,i)  = d_zero
            sina(j,i)  = d_zero
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
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: a
    real(dp) , pointer , intent(out) , dimension(:,:) :: b
    integer , intent(in) :: l
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
