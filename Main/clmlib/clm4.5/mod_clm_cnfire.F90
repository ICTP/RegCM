module mod_clm_cnfire
#ifdef CN
  !
  ! module for fire dynamics
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)"
  ! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the
  ! 20th Century transient simulations at f19_g16 with
  ! (newfire05_clm45sci15_clm4_0_58)
  ! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
  ! climatological lightning data.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_mpmessage
  use mod_runparams , only : dtsrf , idatex , idate1 , ktau , dtsec
  use mod_dynparam
  use mod_stdio
  use mod_date
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_nchelper
  use mod_clm_subgridave , only : p2c
  use mod_clm_varpar , only : nlevdecomp , ndecomp_pools
  use mod_clm_varpar , only : maxpatch_pft , max_pft_per_col
  use mod_clm_varcon , only : dzsoi_decomp , rpi , tfrz , secspday
  use mod_clm_domain , only : ldomain
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varctl , only : fsurdat , inst_name
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_pftvarcon , only : fsr_pft , fd_pft , noveg
  use mod_clm_pftvarcon , only : nc4_grass , nc3crop , ndllf_evr_tmp_tree
  use mod_clm_pftvarcon , only : nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree
  use mod_clm_pftvarcon , only : nbrdlf_evr_shrub
  use mod_clm_pftvarcon , only : cc_leaf , cc_lstem , cc_dstem , cc_other
  use mod_clm_pftvarcon , only : fm_leaf , fm_lstem , fm_dstem , fm_other
  use mod_clm_pftvarcon , only : fm_root , fm_lroot , fm_droot
  use mod_clm_pftvarcon , only : nc3crop , lf_flab , lf_fcel , lf_flig
  use mod_clm_pftvarcon , only : fr_flab , fr_fcel , fr_flig
  use mod_clm_varctl , only : inst_name
  use mod_clm_histfile , only : hist_addfld1d

  implicit none

  private

  save

  public :: CNFireInit    ! Initialization of CNFire
  public :: CNFireInterp  ! Interpolate fire data
  public :: CNFireArea    ! Calculate fire area
  public :: CNFireFluxes  ! Calculate fire fluxes

  ! position datasets for dynamic human population density
  private :: hdm_init
  ! interpolates between two years of human pop. density file data
  private :: hdm_interp
  ! position datasets for Lightning
  private :: lnfm_init
  ! interpolates between two years of Lightning file data
  private :: lnfm_interp

  real(rkx) , pointer , dimension(:) :: forc_lnfm  ! Lightning frequency
  real(rkx) , pointer , dimension(:) :: forc_hdm   ! Human population density
  real(rkx) , pointer , dimension(:) :: hdm_p1 , hdm_p2
  real(rkx) , pointer , dimension(:) :: lnfm_p1 , lnfm_p2
  real(rkx) , parameter :: secsphr = 3600._rkx  ! Seconds in an hour
  real(rkx) , parameter :: borealat = 40._rkx   ! Latitude for boreal peat fires

  ! Human population density input data stream
  type(clm_filetype) :: sdat_hdm
  integer(ik4) :: ipoprec
  ! Lightning input data stream
  type(clm_filetype) :: sdat_lnfm
  integer(ik4) :: ilnfmrec

  contains
  !
  ! Initialize CN Fire module
  !
  subroutine CNFireInit( begg, endg )
    implicit none
    integer(ik4) , intent(in) :: begg , endg   ! gridcell index bounds
    call hdm_init(   begg, endg )
    call lnfm_init(  begg, endg )
    call CNFireInterp( )
  end subroutine CNFireInit
  !
  ! Interpolate CN Fire datasets
  !
  subroutine CNFireInterp()
    implicit none
    call hdm_interp()
    call lnfm_interp()
  end subroutine CNFireInterp
  !
  ! Computes column-level burned area in each timestep
  !
  subroutine CNFireArea (num_soilc, filter_soilc, num_soilp, filter_soilp)
    implicit none
    ! number of soil columns in filter
    integer(ik4) , intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4) , dimension(:) , intent(in) :: filter_soilc
    ! number of soil pfts in filter
    integer(ik4) , intent(in) :: num_soilp
    ! filter for soil pfts
    integer(ik4) , dimension(:) , intent(in) :: filter_soilp
    ! 10-day running mean of tot. precipitation
    real(rkx) , pointer , dimension(:) :: prec10
    ! 60-day running mean of tot. precipitation
    real(rkx) , pointer , dimension(:) :: prec60
    ! decrease of pft weight (0-1) on the col. for timestep
    real(rkx) , pointer , dimension(:) :: lfpftd
    ! pft weight on the column
    real(rkx) , pointer , dimension(:) :: wtcol
    ! vegetation type for this pft
    integer(ik4) , pointer , dimension(:) :: ivt
    ! (gC/m2) dead coarse root C
    real(rkx) , pointer , dimension(:) :: deadcrootc
    ! (gC/m2) dead coarse root C storage
    real(rkx) , pointer , dimension(:) :: deadcrootc_storage
    ! (gC/m2) dead coarse root C transfer
    real(rkx) , pointer , dimension(:) :: deadcrootc_xfer
    ! (gC/m2) fine root C
    real(rkx) , pointer , dimension(:) :: frootc
    ! (gC/m2) fine root C storage
    real(rkx) , pointer , dimension(:) :: frootc_storage
    ! (gC/m2) fine root C transfer
    real(rkx) , pointer , dimension(:) :: frootc_xfer
    ! (gC/m2) live coarse root C
    real(rkx) , pointer , dimension(:) :: livecrootc
    ! (gC/m2) live coarse root C storage
    real(rkx) , pointer , dimension(:) :: livecrootc_storage
    ! (gC/m2) live coarse root C transfer
    real(rkx) , pointer , dimension(:) :: livecrootc_xfer
    ! (gC/m2) total vegetation carbon, excluding cpool
    real(rkx) , pointer , dimension(:) :: totvegc
    ! root zone soil wetness
    real(rkx) , pointer , dimension(:) :: btran2
    ! pft's column index
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! (gC/m2) leaf C
    real(rkx) , pointer , dimension(:) :: leafc
    ! (gC/m2) leaf C storage
    real(rkx) , pointer , dimension(:) :: leafc_storage
    ! (gC/m2) leaf C transfer
    real(rkx) , pointer , dimension(:) :: leafc_xfer
    ! (gC/m2) live stem C
    real(rkx) , pointer , dimension(:) :: livestemc
    ! (gC/m2) live stem C storage
    real(rkx) , pointer , dimension(:) :: livestemc_storage
    ! (gC/m2) live stem C transfer
    real(rkx) , pointer , dimension(:) :: livestemc_xfer
    ! (gC/m2) dead stem C
    real(rkx) , pointer , dimension(:) :: deadstemc
    ! (gC/m2) dead stem C storage
    real(rkx) , pointer , dimension(:) :: deadstemc_storage
    ! (gC/m2) dead stem C transfer
    real(rkx) , pointer , dimension(:) :: deadstemc_xfer
    ! burn date for crop
    integer(ik4) , pointer , dimension(:) :: burndate

    ! column-level
    ! fractional area with water table at surface
    real(rkx) , pointer , dimension(:) :: fsat
    ! conversion area frac. of BET+BDT that haven't burned before
    real(rkx) , pointer , dimension(:) :: lfc
    ! column's weight relative to corresponding gridcell
    real(rkx) , pointer , dimension(:) :: cwtgcell
    ! annual decreased fraction coverage of BET+BDT on gridcell
    real(rkx) , pointer , dimension(:) :: dtrotr_col
    ! pft weight of BET on the gridcell (0-1)
    real(rkx) , pointer , dimension(:) :: trotr1_col
    ! pft weight of BDT on the gridcell (0-1)
    real(rkx) , pointer , dimension(:) :: trotr2_col
    ! 10-day running mean of tot. precipitation
    real(rkx) , pointer , dimension(:) :: prec10_col
    ! 60-day running mean of tot. precipitation
    real(rkx) , pointer , dimension(:) :: prec60_col
    ! number of pfts on the column
    integer(ik4) , pointer , dimension(:) :: npfts
    ! pft index array
    integer(ik4) , pointer , dimension(:) :: pfti
    ! gridcell of corresponding column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! soil water as frac. of whc for top 0.05 m
    real(rkx) , pointer , dimension(:) :: wf
    ! soil water as frac. of whc for top 0.17 m
    real(rkx) , pointer , dimension(:) :: wf2
    ! soil T for top 0.17 m
    real(rkx) , pointer , dimension(:) :: tsoi17
    ! gdp data
    real(rkx) , pointer , dimension(:) :: gdp_lf
    ! peatland fraction data
    real(rkx) , pointer , dimension(:) :: peatf_lf
    ! proscribed crop fire time
    integer(ik4), pointer , dimension(:) :: abm_lf
    ! (gC/m2) total lit C (column-level mean)
    real(rkx) , pointer , dimension(:) :: totlitc
    ! fire spread rate at column level
    real(rkx) , pointer , dimension(:) :: fsr_col
    ! fire duration at column level
    real(rkx) , pointer , dimension(:) :: fd_col
    ! root carbon
    real(rkx) , pointer , dimension(:) :: rootc_col
    ! burned area fraction for cropland
    real(rkx) , pointer , dimension(:) :: baf_crop
    ! burned area fraction for peatland
    real(rkx) , pointer , dimension(:) :: baf_peatf
    ! total burned area out of conversion
    real(rkx) , pointer , dimension(:) :: fbac
    ! burned area out of conversion region due to land use fire
    real(rkx) , pointer , dimension(:) :: fbac1
    ! cropland fraction in veg column
    real(rkx) , pointer , dimension(:) :: cropf_col
    ! transpiration wetness factor (0 to 1)
    real(rkx) , pointer , dimension(:) :: btran_col
    ! fractional coverage of non-crop PFTs
    real(rkx) , pointer , dimension(:) :: wtlf
    ! fractional coverage of non-crop and non-bare-soil PFTs
    real(rkx) , pointer , dimension(:) :: lfwt
    ! totvegc at column level
    real(rkx) , pointer , dimension(:) :: totvegc_col
    ! leaf carbon at column level
    real(rkx) , pointer , dimension(:) :: leafc_col
    ! gdp limitation factor for nfire
    real(rkx) , pointer , dimension(:) :: lgdp_col
    ! gdp limitation factor for baf per fire
    real(rkx) , pointer , dimension(:) :: lgdp1_col
    ! pop limitation factor for baf per fire
    real(rkx) , pointer , dimension(:) :: lpop_col
    ! fuel avalability factor for Reg.C
    real(rkx) , pointer , dimension(:) :: fuelc
    ! fuel avalability factor for Reg.A
    real(rkx) , pointer , dimension(:) :: fuelc_crop
    ! (gC/m3)  vert.-resolved decomposing c pools
    real(rkx) , pointer , dimension(:,:,:) :: decomp_cpools_vr

    ! grid-level
    ! latitude (degrees)
    real(rkx) , pointer , dimension(:) :: latdeg
    ! rain
    real(rkx) , pointer , dimension(:) :: forc_rain
    ! snow
    real(rkx) , pointer , dimension(:) :: forc_snow
    ! relative humidity
    real(rkx) , pointer , dimension(:) :: forc_rh
    ! atmospheric temperature (Kelvin)
    real(rkx) , pointer , dimension(:) :: forc_t
    !atmospheric wind speed (m/s)
    real(rkx) , pointer , dimension(:) :: forc_wind
    ! column-level
    ! fire counts (count/km2/timestep), valid only in Reg. C
    real(rkx) , pointer , dimension(:) :: nfire
    ! fractional area burned in this timestep
    real(rkx) , pointer , dimension(:) :: farea_burned
    ! TRUE => pool is a cwd pool
    logical , pointer , dimension(:) :: is_cwd
    ! lower threshold of fuel mass (gC/m2) for ignition
    real(rkx) , parameter :: lfuel =  110._rkx
    ! upper threshold of fuel mass(gC/m2) for ignition
    real(rkx) , parameter :: ufuel = 1050._rkx
    ! g(W) when W=0 m/s
    real(rkx) , parameter :: g0 = 0.05_rkx
    ! a1 parameter for cropland fire in Li et. al. 2013 (was different in paper)
    real(rkx) , parameter :: cropfire_a1 = 0.153_rkx
    ! c parameter for peatland fire in Li et. al. 2013
    ! boreal peat fires (was different in paper)
    real(rkx) , parameter :: boreal_peatfire_c = 2.1d-5
    ! non-boreal peat fires (was different in paper)
    real(rkx) , parameter :: non_boreal_peatfire_c = 0.0005d00

    integer(ik4) :: g , l , c , p , pi , j , fc , fp
    integer(ik4) :: kyr , kmo , kda ! index variables
    real(rkx) :: dt        ! time step variable (s)
    real(rkx) :: m         ! top-layer soil moisture (proportion)
    real(rkx) ::cli       !
    real(rkx) , parameter ::cli_scale = 1370.0_rkx
    real(rkx) ::cri       !
    real(rkx) :: fb        ! availability of fuel
    real(rkx) :: fhd       ! impact of hd on agricultural fire
    real(rkx) :: fgdp      ! impact of gdp on agricultural fire
    real(rkx) :: fire_m    ! combustability of fuel on fire occurrence
    real(rkx) :: spread_m  ! combustability of fuel on fire spread
    real(rkx) :: Lb_lf     ! length-to-breadth ratio added by Lifang
    integer(ik4) :: i_cwd     ! cwd pool
    real(rkx) :: lh       !
    real(rkx) :: fs       !
    real(rkx) :: ig       !
    real(rkx) :: hdmlf    ! human density

    wtcol              => clm3%g%l%c%p%wtcol
    ivt                => clm3%g%l%c%p%itype
    prec60             => clm3%g%l%c%p%pps%prec60
    prec10             => clm3%g%l%c%p%pps%prec10
    deadcrootc         => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer    => clm3%g%l%c%p%pcs%deadcrootc_xfer
    frootc             => clm3%g%l%c%p%pcs%frootc
    frootc_storage     => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer        => clm3%g%l%c%p%pcs%frootc_xfer
    livecrootc         => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer    => clm3%g%l%c%p%pcs%livecrootc_xfer
    totvegc            => clm3%g%l%c%p%pcs%totvegc
    btran2             => clm3%g%l%c%p%pps%btran2
    pcolumn            => clm3%g%l%c%p%column
    leafc              => clm3%g%l%c%p%pcs%leafc
    leafc_storage      => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer         => clm3%g%l%c%p%pcs%leafc_xfer
    livestemc          => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage  => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer     => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc          => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage  => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer     => clm3%g%l%c%p%pcs%deadstemc_xfer
    lfpftd             => clm3%g%l%c%p%pps%lfpftd
    burndate           => clm3%g%l%c%p%pps%burndate

    ! assign local pointers to derived type members (column-level)
    cwtgcell         => clm3%g%l%c%wtgcell
    npfts            => clm3%g%l%c%npfts
    pfti             => clm3%g%l%c%pfti
    wf               => clm3%g%l%c%cps%wf
    wf2              => clm3%g%l%c%cps%wf2
    tsoi17           => clm3%g%l%c%ces%tsoi17
    farea_burned     => clm3%g%l%c%cps%farea_burned
    baf_crop         => clm3%g%l%c%cps%baf_crop
    baf_peatf        => clm3%g%l%c%cps%baf_peatf
    fbac             => clm3%g%l%c%cps%fbac
    fbac1            => clm3%g%l%c%cps%fbac1
    cropf_col        => clm3%g%l%c%cps%cropf_col
    gdp_lf           => clm3%g%l%c%cps%gdp_lf
    peatf_lf         => clm3%g%l%c%cps%peatf_lf
    abm_lf           => clm3%g%l%c%cps%abm_lf
    nfire            => clm3%g%l%c%cps%nfire
    totlitc          => clm3%g%l%c%ccs%totlitc
    fsr_col          => clm3%g%l%c%cps%fsr_col
    fd_col           => clm3%g%l%c%cps%fd_col
    rootc_col        => clm3%g%l%c%ccs%rootc_col
    totvegc_col      => clm3%g%l%c%ccs%totvegc_col
    leafc_col        => clm3%g%l%c%ccs%leafc_col
    lgdp_col         => clm3%g%l%c%cps%lgdp_col
    lgdp1_col        => clm3%g%l%c%cps%lgdp1_col
    lpop_col         => clm3%g%l%c%cps%lpop_col
    fuelc            => clm3%g%l%c%ccs%fuelc
    fuelc_crop       => clm3%g%l%c%ccs%fuelc_crop
    btran_col        => clm3%g%l%c%cps%btran_col
    wtlf             => clm3%g%l%c%cps%wtlf
    lfwt             => clm3%g%l%c%cps%lfwt
    cgridcell        => clm3%g%l%c%gridcell
    trotr1_col       => clm3%g%l%c%cps%trotr1_col
    trotr2_col       => clm3%g%l%c%cps%trotr2_col
    dtrotr_col       => clm3%g%l%c%cps%dtrotr_col
    prec60_col       => clm3%g%l%c%cps%prec60_col
    prec10_col       => clm3%g%l%c%cps%prec10_col
    lfc              => clm3%g%l%c%cps%lfc
    fsat             => clm3%g%l%c%cws%fsat
    is_cwd           => decomp_cascade_con%is_cwd
    decomp_cpools_vr => clm3%g%l%c%ccs%decomp_cpools_vr

    !assign local pointers to derived type members (grid-level)
    forc_rh    => clm_a2l%forc_rh
    forc_wind  => clm_a2l%forc_wind
    forc_t     => clm_a2l%forc_t
    forc_rain  => clm_a2l%forc_rain
    forc_snow  => clm_a2l%forc_snow
    latdeg     => clm3%g%latdeg

    !pft to column average
    call p2c(num_soilc, filter_soilc, prec10,  prec10_col)
    call p2c(num_soilc, filter_soilc, prec60,  prec60_col)
    call p2c(num_soilc, filter_soilc,totvegc, totvegc_col)
    call p2c(num_soilc, filter_soilc,leafc, leafc_col)
    ! Get model step size
    dt = dtsrf
    !
    ! On first time-step, just set area burned to zero and exit
    !
    if ( ktau == 0 ) then
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        farea_burned(c) = 0._rkx
        baf_crop(c)     = 0._rkx
        baf_peatf(c)    = 0._rkx
        fbac(c)         = 0._rkx
        fbac1(c)        = 0._rkx
      end do
      return
    end if
    !
    ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil
    ! vegetation (lfwt) in vegetated column
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      cropf_col(c) = 0._rkx
      lfwt(c)      = 0._rkx
    end do
    do pi = 1 , max_pft_per_col
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        if ( pi <=  npfts(c) ) then
          p = pfti(c) + pi - 1
          ! For crop veg types
          if ( ivt(p) > nc4_grass ) then
            cropf_col(c) = cropf_col(c) + wtcol(p)
          end if
          ! For natural vegetation (non-crop)
          if ( ivt(p) >= ndllf_evr_tmp_tree .and. ivt(p) <= nc4_grass ) then
            lfwt(c) = lfwt(c) + wtcol(p)
          end if
        end if
      end do
    end do
    !
    ! Calculate crop fuel
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      fuelc_crop(c)=0._rkx
    end do
    do pi = 1 , max_pft_per_col
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        if ( pi <=  npfts(c) ) then
          p = pfti(c) + pi - 1
          ! For crop PFT's
          if ( ivt(p) > nc4_grass .and. wtcol(p) > 0._rkx .and. &
                  leafc_col(c) > 0._rkx ) then
            fuelc_crop(c) = fuelc_crop(c) + (leafc(p) + leafc_storage(p) + &
                      leafc_xfer(p))*wtcol(p)/cropf_col(c)     + &
                      totlitc(c)*leafc(p)/leafc_col(c)*wtcol(p)/cropf_col(c)
          end if
        end if
      end do
    end do
    !
    ! Calculate noncrop column variables
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      fsr_col(c)   = 0._rkx
      fd_col(c)    = 0._rkx
      rootc_col(c) = 0._rkx
      lgdp_col(c)  = 0._rkx
      lgdp1_col(c) = 0._rkx
      lpop_col(c)  = 0._rkx
      btran_col(c) = 0._rkx
      wtlf(c)      = 0._rkx
#ifdef DYNPFT
      trotr1_col(c) = 0._rkx
      trotr2_col(c) = 0._rkx
      dtrotr_col(c) = 0._rkx
#endif
    end do
    do pi = 1 , max_pft_per_col
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        g = cgridcell(c)
        if ( pi <= npfts(c) ) then
          p = pfti(c) + pi - 1
          ! For non-crop -- natural vegetation and bare-soil
          if ( ivt(p) < nc3crop .and. cropf_col(c) < 1.0_rkx ) then
            if ( .not. is_nan(btran2(p)) .and. &
                    btran2(p) <= 1._rkx ) then
              btran_col(c) = btran_col(c)+btran2(p)*wtcol(p)
              wtlf(c)      = wtlf(c)+wtcol(p)
            end if
#ifdef DYNPFT
            if ( ivt(p) == nbrdlf_evr_trp_tree .and. wtcol(p) > 0._rkx ) then
              trotr1_col(c) = trotr1_col(c) + wtcol(p)*cwtgcell(c)
            end if
            if ( ivt(p) == nbrdlf_dcd_trp_tree .and. wtcol(p) > 0._rkx ) then
              trotr2_col(c) = trotr2_col(c) + wtcol(p)*cwtgcell(c)
            end if
            if ( ivt(p) == nbrdlf_evr_trp_tree .or. &
                 ivt(p) == nbrdlf_dcd_trp_tree ) then
              if ( lfpftd(p) > 0._rkx ) then
                dtrotr_col(c) = dtrotr_col(c)+lfpftd(p)*cwtgcell(c)
              end if
            end if
#endif
            rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                           frootc_xfer(p) + deadcrootc(p) +                &
                           deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                           livecrootc(p)+livecrootc_storage(p) +           &
                           livecrootc_xfer(p))*wtcol(p)

            fsr_col(c) = fsr_col(c) + &
                    fsr_pft(ivt(p))*wtcol(p)/(1.0_rkx-cropf_col(c))

            if ( lfwt(c) /= 0.0_rkx ) then
              hdmlf = forc_hdm(g)

              ! all these constants are in Li et al. BG (2012a,b;2013)

              if ( hdmlf > 0.1_rkx ) then
                ! For NOT bare-soil
                if ( ivt(p) /= noveg ) then
                  ! For shrub and grass (crop already excluded above)
                  if ( ivt(p) >= nbrdlf_evr_shrub ) then !for shurb and grass
                    lgdp_col(c)  = lgdp_col(c) + (0.1_rkx + 0.9_rkx*    &
                                      exp(-1._rkx*rpi* &
                                      (gdp_lf(c)/8._rkx)**0.5_rkx))*wtcol(p) &
                                      /(1.0_rkx - cropf_col(c))
                    lgdp1_col(c) = lgdp1_col(c) + (0.2_rkx + 0.8_rkx*   &
                                     exp(-1._rkx*rpi* &
                                     (gdp_lf(c)/7._rkx)))*wtcol(p)/lfwt(c)
                    lpop_col(c)  = lpop_col(c) + (0.2_rkx + 0.8_rkx*    &
                                     exp(-1._rkx*rpi* &
                                   (hdmlf/450._rkx)**0.5_rkx))*wtcol(p)/lfwt(c)
                  else   ! for trees
                    if ( gdp_lf(c) > 20._rkx ) then
                      lgdp_col(c) = lgdp_col(c) + &
                              0.39_rkx*wtcol(p)/(1.0_rkx - cropf_col(c))
                    else
                      lgdp_col(c) = lgdp_col(c)+wtcol(p)/(1.0_rkx - cropf_col(c))
                    end if
                    if ( gdp_lf(c) > 20._rkx ) then
                      lgdp1_col(c) = lgdp1_col(c)+0.62_rkx*wtcol(p)/lfwt(c)
                    else
                      if ( gdp_lf(c) > 8._rkx ) then
                        lgdp1_col(c)=lgdp1_col(c)+0.83_rkx*wtcol(p)/lfwt(c)
                      else
                        lgdp1_col(c)=lgdp1_col(c)+wtcol(p)/lfwt(c)
                      end if
                    end if
                    lpop_col(c) = lpop_col(c) + (0.4_rkx + 0.6_rkx*    &
                                        exp(-1._rkx*rpi* &
                                        (hdmlf/125._rkx)))*wtcol(p)/lfwt(c)
                  end if
                end if
              else
                lgdp_col(c)  = lgdp_col(c)+wtcol(p)/(1.0_rkx - cropf_col(c))
                lgdp1_col(c) = lgdp1_col(c)+wtcol(p)/lfwt(c)
                lpop_col(c)  = lpop_col(c)+wtcol(p)/lfwt(c)
              end if
            end if
            fd_col(c) = fd_col(c) + &
                    fd_pft(ivt(p))*wtcol(p)*secsphr/(1.0_rkx-cropf_col(c))
          end if
        end if
      end do
    end do

#ifdef DYNPFT
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      if ( dtrotr_col(c) > 0._rkx ) then
        if ( date_is(idatex,1,1) .and. time_is(idatex,0,dtsrf) ) then
          lfc(c) = 0._rkx
        end if
        if ( date_is(idatex,1,1) .and. &
             time_is(idatex,int(dtsrf+dt/2.0_rkx)) ) then
          lfc(c) = dtrotr_col(c)*dayspy*secspday/dt
        end if
      else
        lfc(c) = 0._rkx
      end if
    end do
#endif
    !
    ! calculate burned area fraction in cropland
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      baf_crop(c) = 0._rkx
    end do

    do fp = 1 , num_soilp
      p = filter_soilp(fp)
      if ( date_is(idatex,1,1) .and. time_is(idatex,0,dtsrf) ) then
        burndate(p) = 10000 ! init. value; actual range [0 365]
      end if
    end do

    call split_idate(idatex,kyr,kmo,kda)
    do pi = 1 , max_pft_per_col
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = cgridcell(c)
        hdmlf = forc_hdm(g)
        if ( pi <=  npfts(c) ) then
          p = pfti(c) + pi - 1
          ! For crop
          if ( forc_t(g) >= tfrz .and. &
               ivt(p) > nc4_grass .and.  &
               kmo == abm_lf(c) .and. &
               forc_rain(g)+forc_snow(g) == 0._rkx  .and. &
               burndate(p) >= 999 .and. &
               wtcol(p) > 0._rkx ) then ! catch  crop burn time
            ! calculate human density impact on ag. fire
            fhd = 0.04_rkx+0.96_rkx*exp(-1._rkx*rpi*(hdmlf/350._rkx)**0.5_rkx)
            ! calculate impact of GDP on ag. fire
            fgdp = 0.01_rkx+0.99_rkx*exp(-1._rkx*rpi*(gdp_lf(c)/10._rkx))
            ! calculate burned area
            fb   = max(0.0_rkx,min(1.0_rkx,(fuelc_crop(c)-lfuel)/(ufuel-lfuel)))
            ! crop fire only for generic crop types at this time
            ! managed crops are treated as grasses if crop model is turned on
            ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND
            ! CURRENTLY DOES NOT!
            !  As such results are only valid for a time-step of a half-hour.
            baf_crop(c) = baf_crop(c) + cropfire_a1*fb*fhd*fgdp*wtcol(p)
            if ( fb*fhd*fgdp*wtcol(p) > 0._rkx ) then
              burndate(p) = kda
            end if
          end if
        end if
      end do
    end do
    !
    ! calculate peatland fire
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)
      ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND
      ! CURRENTLY DOES NOT!
      ! As such results are only valid for a time-step of a half-hour.
      if ( latdeg(g)<borealat ) then
        baf_peatf(c) = non_boreal_peatfire_c*max(0._rkx, &
                       min(1._rkx,(4.0_rkx-prec60_col(c)*secspday)/ &
                       4.0_rkx))**2*peatf_lf(c)*(1._rkx-fsat(c))
      else
        baf_peatf(c) = boreal_peatfire_c*exp(-rpi*(max(wf2(c),0._rkx)/0.3_rkx))* &
          max(0._rkx,min(1._rkx,(tsoi17(c)-tfrz)/10._rkx))*peatf_lf(c)* &
          (1._rkx-fsat(c))
      end if
    end do
    !
    ! calculate other fires
    !

    ! Set the number of timesteps for e-folding.
    ! When the simulation has run fewer than this number of steps,
    ! re-scale the e-folding time to get a stable early estimate.

    ! find which pool is the cwd pool
    i_cwd = 0
    do l = 1 , ndecomp_pools
      if ( is_cwd(l) ) then
        i_cwd = l
      end if
    end do

    !
    ! begin column loop to calculate fractional area affected by fire
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)
      hdmlf = forc_hdm(g)
      if ( cropf_col(c) < 1.0 ) then
        fuelc(c) = totlitc(c)+totvegc_col(c) - &
                   rootc_col(c)-fuelc_crop(c)*cropf_col(c)
        do j = 1 , nlevdecomp
          fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j)
        end do
        fuelc(c) = fuelc(c)/(1._rkx-cropf_col(c))
        fb       = max(0.0_rkx,min(1.0_rkx,(fuelc(c)-lfuel)/(ufuel-lfuel)))
        m        = max(0._rkx,wf(c))
        fire_m   = exp(-rpi *(m/0.69_rkx)**2)*(1.0_rkx - max(0._rkx, &
                   min(1._rkx,(forc_rh(g)-30._rkx)/(70._rkx-30._rkx))))*  &
                   min(1._rkx,exp(rpi*(forc_t(g)-tfrz)/10._rkx))
        lh       = 0.0035_rkx*6.8_rkx*hdmlf**(0.43_rkx)/30._rkx/24._rkx
        fs       = 1._rkx-(0.01_rkx+0.98_rkx*exp(-0.025_rkx*hdmlf))
        ig       = (lh+forc_lnfm(g)*0.25_rkx)*(1._rkx-fs)*(1._rkx-cropf_col(c))
        nfire(c) = ig/secsphr*dt*fb*fire_m*lgdp_col(c) !fire counts/km2/timestep
        Lb_lf    = 1._rkx+10.0_rkx*(1._rkx-exp(-0.06_rkx*forc_wind(g)))
        if ( wtlf(c) > 0.0_rkx )then
          spread_m = (1.0_rkx - max(0._rkx,min(1._rkx,(btran_col(c)/wtlf(c)-0.3_rkx)/ &
                     (0.7_rkx-0.3_rkx))))*(1.0-max(0._rkx, &
                     min(1._rkx,(forc_rh(g)-30._rkx)/(70._rkx-30._rkx))))
        else
          spread_m = 0.0_rkx
        end if
        farea_burned(c) = min(1._rkx,(g0*spread_m*fsr_col(c)* &
             fd_col(c)/1000._rkx)**2*lgdp1_col(c)* &
             lpop_col(c)*nfire(c)*rpi*Lb_lf+ &
             baf_crop(c)+baf_peatf(c))  ! fraction (0-1) per timestep
        !
        ! if landuse change data is used, calculate deforestation fires and
        ! add it in the total of burned area fraction
        !
#ifdef DYNPFT
        if ( trotr1_col(c)+trotr2_col(c) > 0.6_rkx ) then
          if ( (date_is(idatex,1,1) .and. time_is(idatex,0,dtsrf)) .or. &
                dtrotr_col(c) <=0._rkx ) then
            fbac1(c)        = 0._rkx
            farea_burned(c) = baf_crop(c)+baf_peatf(c)
          else
            cri = (4.0_rkx*trotr1_col(c)+1.8_rkx*trotr2_col(c)) / &
                  (trotr1_col(c)+trotr2_col(c))
            cli = (max(0._rkx,min(1._rkx,(cri-prec60_col(c) * &
                    secspday)/cri))**0.5)* &
                  (max(0._rkx,min(1._rkx,(cri-prec10_col(c) * &
                    secspday)/cri))**0.5)* &
                   max(0.0005_rkx,min(1._rkx,19._rkx*dtrotr_col(c) * &
                    dayspy*secspday/dt-0.001_rkx))* &
                   max(0._rkx,min(1._rkx,(0.25_rkx-(forc_rain(g) + &
                   forc_snow(g))*secsphr)/0.25_rkx))
            ! NOTE: THIS SHOULD TAKE INTO ACCOUNT THE TIME-STEP AND
            ! CURRENTLY DOES NOT!
            !  As such results are only valid for a time-step of a half-hour.
            farea_burned(c) = cli/cli_scale +baf_crop(c)+baf_peatf(c)
            ! burned area out of conversion region due to land use fire
            fbac1(c) = max(0._rkx,cli/cli_scale - 2.0_rkx*lfc(c))
          end if
          ! total burned area out of conversion
          fbac(c) = fbac1(c)+baf_crop(c)+baf_peatf(c)
        else
          fbac(c) = farea_burned(c)
        end if
#else
        farea_burned(c) = min(1._rkx,baf_crop(c)+baf_peatf(c))
#endif
      end if

#if (defined NOFIRE)
      ! zero out the fire area if NOFIRE flag is on
      farea_burned(c) = 0._rkx
      baf_crop(c)     = 0._rkx
      baf_peatf(c)    = 0._rkx
      fbac(c)         = 0._rkx
      fbac1(c)        = 0._rkx
      ! with NOFIRE, tree carbon is still removed in landuse change
      ! regions by the landuse code
#endif
    end do  ! end of column loop
  end subroutine CNFireArea
  !
  ! Fire effects routine for coupled carbon-nitrogen code (CN).
  ! Relies primarily on estimate of fractional area burned in this
  ! timestep, from CNFireArea().
  !
  ! Total fire carbon emissions (g C/m2 land area/yr)
  !   = avg(COL_FIRE_CLOSS)*seconds_per_year +
  !     avg(SOMC_FIRE)*seconds_per_year +
  !     avg(LF_CONV_CFLUX)*seconds_per_year *
  !         min(1.0,avg(LFC2)/dt*seconds_per_year)*0.8
  ! where dt is the time step size (sec),avg means the temporal average
  ! in a year seconds_per_year is the number of seconds in a year.
  !
  subroutine CNFireFluxes (num_soilc, filter_soilc, num_soilp, filter_soilp)
    implicit none
    ! number of soil columns in filter
    integer(ik4) , intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4) , intent(in) , dimension(:) :: filter_soilc
    ! number of soil pfts in filter
    integer(ik4) , intent(in) :: num_soilp
    ! filter for soil pfts
    integer(ik4) , intent(in) , dimension(:) :: filter_soilp
#if (defined CNDV)
    ! number of individuals (#/m2)
    real(rkx) , pointer , dimension(:) :: nind
#endif
    ! woody lifeform (1=woody, 0=not woody)
    real(rkx) , pointer , dimension(:) :: woody
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    ! pft vegetation type
    integer(ik4) , pointer :: ivt(:)
    ! pft weight relative to column
    real(rkx) , pointer :: wtcol(:)
    ! latitude (degrees)
    real(rkx) , pointer :: latdeg(:)
    ! gridcell of corresponding column
    integer(ik4) , pointer :: cgridcell(:)
    ! number of pfts for each column
    integer(ik4) , pointer :: npfts(:)
    ! beginning pft index for each column
    integer(ik4) , pointer :: pfti(:)
    ! pft's column index
    integer(ik4) , pointer :: pcolumn(:)
    ! timestep fractional area burned (proportion)
    real(rkx) , pointer :: farea_burned(:)
    ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
    real(rkx) , pointer :: fire_mortality_c_to_cwdc(:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rkx) , pointer :: decomp_cpools_vr(:,:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rkx) , pointer :: decomp_npools_vr(:,:,:)
    ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
    real(rkx) , pointer :: fire_mortality_n_to_cwdn(:,:)
    ! conversion area frac. of BET+BDT that haven't burned before
    real(rkx) , pointer :: lfc(:)
    ! conversion area frac. of BET+BDT that burned this timestep
    real(rkx) , pointer :: lfc2(:)
    ! burned area out of conversion region due to land use fire
    real(rkx) , pointer :: fbac1(:)
    ! baf for cropland
    real(rkx) , pointer :: baf_crop(:)
    ! baf for peatlabd
    real(rkx) , pointer :: baf_peatf(:)
    ! (gC/m2) ann max leaf C
    real(rkx) , pointer :: leafcmax(:)
    ! total burned area out of conversion
    real(rkx) , pointer :: fbac(:)

    ! annual decreased fraction coverage of BET+BDT (0-1) on the gridcell
    real(rkx) , pointer :: dtrotr_col(:)
    ! pft weight of BET on the gridcell (0-1)
    real(rkx) , pointer :: trotr1_col(:)
    ! pft weight of BDT on the gridcell (0-1)
    real(rkx) , pointer :: trotr2_col(:)

    ! (gC/m2) total soil organic matter carbon
    real(rkx) , pointer :: totsomc(:)
    ! (gC/m2/s)fire carbon emissions due to peat burning
    real(rkx) , pointer :: somc_fire(:)

    ! (gC/m2) leaf C
    real(rkx) , pointer :: leafc(:)
    ! (gC/m2) leaf C storage
    real(rkx) , pointer :: leafc_storage(:)
    ! (gC/m2) leaf C transfer
    real(rkx) , pointer :: leafc_xfer(:)
    ! (gC/m2) live stem C
    real(rkx) , pointer :: livestemc(:)
    ! (gC/m2) live stem C storage
    real(rkx) , pointer :: livestemc_storage(:)
    ! (gC/m2) live stem C transfer
    real(rkx) , pointer :: livestemc_xfer(:)

    ! (gC/m2) dead stem C
    real(rkx) , pointer :: deadstemc(:)
    ! (gC/m2) dead stem C storage
    real(rkx) , pointer :: deadstemc_storage(:)
    ! (gC/m2) dead stem C transfer
    real(rkx) , pointer :: deadstemc_xfer(:)
    ! (gC/m2) fine root C
    real(rkx) , pointer :: frootc(:)
    ! (gC/m2) fine root C storage
    real(rkx) , pointer :: frootc_storage(:)
    ! (gC/m2) fine root C transfer
    real(rkx) , pointer :: frootc_xfer(:)
    ! (gC/m2) dead coarse root C
    real(rkx) , pointer :: deadcrootc(:)
    ! (gC/m2) dead coarse root C storage
    real(rkx) , pointer :: deadcrootc_storage(:)
    ! (gC/m2) dead coarse root C transfer
    real(rkx) , pointer :: deadcrootc_xfer(:)
    ! (gC/m2) live coarse root C
    real(rkx) , pointer :: livecrootc(:)
    ! (gC/m2) live coarse root C storage
    real(rkx) , pointer :: livecrootc_storage(:)
    ! (gC/m2) live coarse root C transfer
    real(rkx) , pointer :: livecrootc_xfer(:)
    ! (gC/m2) growth respiration storage
    real(rkx) , pointer :: gresp_storage(:)
    ! (gC/m2) growth respiration transfer
    real(rkx) , pointer :: gresp_xfer(:)

    ! (gN/m2) leaf N
    real(rkx) , pointer :: leafn(:)
    ! (gN/m2) leaf N storage
    real(rkx) , pointer :: leafn_storage(:)
    ! (gN/m2) leaf N transfer
    real(rkx) , pointer :: leafn_xfer(:)
    ! (gN/m2) live stem N
    real(rkx) , pointer :: livestemn(:)
    ! (gN/m2) live stem N storage
    real(rkx) , pointer :: livestemn_storage(:)
    ! (gN/m2) live stem N transfer
    real(rkx) , pointer :: livestemn_xfer(:)
    ! (gN/m2) dead stem N
    real(rkx) , pointer :: deadstemn(:)
    ! (gN/m2) dead stem N storage
    real(rkx) , pointer :: deadstemn_storage(:)
    ! (gN/m2) dead stem N transfer
    real(rkx) , pointer :: deadstemn_xfer(:)
    ! (gN/m2) fine root N
    real(rkx) , pointer :: frootn(:)
    ! (gN/m2) fine root N storage
    real(rkx) , pointer :: frootn_storage(:)
    ! (gN/m2) fine root N transfer
    real(rkx) , pointer :: frootn_xfer(:)
    ! (gN/m2) live coarse root N
    real(rkx) , pointer :: livecrootn(:)
    ! (gN/m2) live coarse root N storage
    real(rkx) , pointer :: livecrootn_storage(:)
    ! (gN/m2) live coarse root N transfer
    real(rkx) , pointer :: livecrootn_xfer(:)
    ! (gN/m2) dead coarse root N
    real(rkx) , pointer :: deadcrootn(:)
    ! (gN/m2) dead coarse root N storage
    real(rkx) , pointer :: deadcrootn_storage(:)
    ! (gN/m2) dead coarse root N transfer
    real(rkx) , pointer :: deadcrootn_xfer(:)
    ! (gN/m2) plant pool of retranslocated N
    real(rkx) , pointer :: retransn(:)

    ! (gC/m2/s) fire C emissions from leafc
    real(rkx) , pointer :: m_leafc_to_fire(:)
    ! (gC/m2/s) fire C emissions from leafc_storage
    real(rkx) , pointer :: m_leafc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from leafc_xfer
    real(rkx) , pointer :: m_leafc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from livestemc
    real(rkx) , pointer :: m_livestemc_to_fire(:)
    ! (gC/m2/s) fire C emissions from livestemc_storage
    real(rkx) , pointer :: m_livestemc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from livestemc_xfer
    real(rkx) , pointer :: m_livestemc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rkx) , pointer :: m_deadstemc_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadstemc_storage
    real(rkx) , pointer :: m_deadstemc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rkx) , pointer :: m_deadstemc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from frootc
    real(rkx) , pointer :: m_frootc_to_fire(:)
    ! (gC/m2/s) fire C emissions from frootc_storage
    real(rkx) , pointer :: m_frootc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from frootc_xfer
    real(rkx) , pointer :: m_frootc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from livecrootc
    real(rkx) , pointer :: m_livecrootc_to_fire(:)
    ! (gC/m2/s) fire C emissions from livecrootc_storage
    real(rkx) , pointer :: m_livecrootc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from livecrootc_xfer
    real(rkx) , pointer :: m_livecrootc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadcrootc
    real(rkx) , pointer :: m_deadcrootc_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadcrootc_storage
    real(rkx) , pointer :: m_deadcrootc_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from deadcrootc_xfer
    real(rkx) , pointer :: m_deadcrootc_xfer_to_fire(:)
    ! (gC/m2/s) fire C emissions from gresp_storage
    real(rkx) , pointer :: m_gresp_storage_to_fire(:)
    ! (gC/m2/s) fire C emissions from gresp_xfer
    real(rkx) , pointer :: m_gresp_xfer_to_fire(:)
    ! (gC/m3/s) vertically-resolved decomposing C fire loss
    real(rkx) , pointer :: m_decomp_cpools_to_fire_vr(:,:,:)

    ! (gN/m2/s) fire N emissions from leafn
    real(rkx) , pointer :: m_leafn_to_fire(:)
    ! (gN/m2/s) fire N emissions from leafn_storage
    real(rkx) , pointer :: m_leafn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from leafn_xfer
    real(rkx) , pointer :: m_leafn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from livestemn
    real(rkx) , pointer :: m_livestemn_to_fire(:)
    ! (gN/m2/s) fire N emissions from livestemn_storage
    real(rkx) , pointer :: m_livestemn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from livestemn_xfer
    real(rkx) , pointer :: m_livestemn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadstemn
    real(rkx) , pointer :: m_deadstemn_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadstemn_storage
    real(rkx) , pointer :: m_deadstemn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadstemn_xfer
    real(rkx) , pointer :: m_deadstemn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from frootn
    real(rkx) , pointer :: m_frootn_to_fire(:)
    ! (gN/m2/s) fire N emissions from frootn_storage
    real(rkx) , pointer :: m_frootn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from frootn_xfer
    real(rkx) , pointer :: m_frootn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
    real(rkx) , pointer :: m_livecrootn_to_fire(:)
    ! (gN/m2/s) fire N emissions from livecrootn_storage
    real(rkx) , pointer :: m_livecrootn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from livecrootn_xfer
    real(rkx) , pointer :: m_livecrootn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadcrootn
    real(rkx) , pointer :: m_deadcrootn_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadcrootn_storage
    real(rkx) , pointer :: m_deadcrootn_storage_to_fire(:)
    ! (gN/m2/s) fire N emissions from deadcrootn_xfer
    real(rkx) , pointer :: m_deadcrootn_xfer_to_fire(:)
    ! (gN/m2/s) fire N emissions from retransn
    real(rkx) , pointer :: m_retransn_to_fire(:)
    ! vertically-resolved decomposing N fire loss (gN/m3/s)
    real(rkx) , pointer :: m_decomp_npools_to_fire_vr(:,:,:)

    ! (gC/m2/s) C transfers from various C pools to litter and
    ! cwd pools due to fire mortality
    real(rkx) , pointer :: m_leafc_to_litter_fire(:)
    real(rkx) , pointer :: m_leafc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_leafc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemc_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemc_to_deadstemc_fire(:)
    real(rkx) , pointer :: m_deadstemc_to_litter_fire(:)
    real(rkx) , pointer :: m_deadstemc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_deadstemc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_frootc_to_litter_fire(:)
    real(rkx) , pointer :: m_frootc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_frootc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootc_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootc_to_deadcrootc_fire(:)
    real(rkx) , pointer :: m_deadcrootc_to_litter_fire(:)
    real(rkx) , pointer :: m_deadcrootc_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_deadcrootc_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_gresp_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_gresp_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_c_to_litr_met_fire(:,:)
    real(rkx) , pointer :: m_c_to_litr_cel_fire(:,:)
    real(rkx) , pointer :: m_c_to_litr_lig_fire(:,:)

    ! (gN/m2/s) N transfers from various C pools to litter and
    ! cwd pools due to fire mortality
    real(rkx) , pointer :: m_leafn_to_litter_fire(:)
    real(rkx) , pointer :: m_leafn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_leafn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemn_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livestemn_to_deadstemn_fire(:)
    real(rkx) , pointer :: m_deadstemn_to_litter_fire(:)
    real(rkx) , pointer :: m_deadstemn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_deadstemn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_frootn_to_litter_fire(:)
    real(rkx) , pointer :: m_frootn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_frootn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootn_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_livecrootn_to_deadcrootn_fire(:)
    real(rkx) , pointer :: m_deadcrootn_to_litter_fire(:)
    real(rkx) , pointer :: m_deadcrootn_storage_to_litter_fire(:)
    real(rkx) , pointer :: m_deadcrootn_xfer_to_litter_fire(:)
    real(rkx) , pointer :: m_retransn_to_litter_fire(:)
    real(rkx) , pointer :: m_n_to_litr_met_fire(:,:)
    real(rkx) , pointer :: m_n_to_litr_cel_fire(:,:)
    real(rkx) , pointer :: m_n_to_litr_lig_fire(:,:)

    ! TRUE => pool is a cwd pool
    logical , pointer  :: is_cwd(:)
    ! TRUE => pool is a litter pool
    logical , pointer  :: is_litter(:)
    ! (1/m) profile of fine roots
    real(rkx) , pointer :: froot_prof(:,:)
    ! (1/m) profile of coarse roots
    real(rkx) , pointer :: croot_prof(:,:)
    ! (1/m) profile of stems
    real(rkx) , pointer :: stem_prof(:,:)
    ! (1/m) profile of leaves
    real(rkx) , pointer :: leaf_prof(:,:)
    integer(ik4) :: g , c , p , j , l , pi
    integer(ik4) :: fp , fc           ! filter indices
    real(rkx) :: f                    ! rate for fire effects (1/s)
    real(rkx) :: dt                   ! time step variable (s)

    ! assign local pointers

#if (defined CNDV)
    nind                     => clm3%g%l%c%p%pdgvs%nind
#endif
    pcolumn                   => clm3%g%l%c%p%column
    cgridcell                 => clm3%g%l%c%gridcell
    farea_burned              => clm3%g%l%c%cps%farea_burned
    woody                     => pftcon%woody
    fire_mortality_c_to_cwdc  => clm3%g%l%c%ccf%fire_mortality_c_to_cwdc
    fire_mortality_n_to_cwdn  => clm3%g%l%c%cnf%fire_mortality_n_to_cwdn

    lfc                       => clm3%g%l%c%cps%lfc
    lfc2                      => clm3%g%l%c%cps%lfc2
    fbac1                     => clm3%g%l%c%cps%fbac1
    fbac                      => clm3%g%l%c%cps%fbac
    baf_crop                  => clm3%g%l%c%cps%baf_crop
    baf_peatf                 => clm3%g%l%c%cps%baf_peatf
    leafcmax                  => clm3%g%l%c%p%pcs%leafcmax
    latdeg                    => clm3%g%latdeg
    wtcol                     => clm3%g%l%c%p%wtcol
    pfti                      => clm3%g%l%c%pfti

    ivt                       => clm3%g%l%c%p%itype
    npfts                     => clm3%g%l%c%npfts

    trotr1_col                => clm3%g%l%c%cps%trotr1_col
    trotr2_col                => clm3%g%l%c%cps%trotr2_col
    dtrotr_col                => clm3%g%l%c%cps%dtrotr_col

    somc_fire                 => clm3%g%l%c%ccf%somc_fire
    totsomc                   => clm3%g%l%c%ccs%totsomc
    decomp_cpools_vr          => clm3%g%l%c%ccs%decomp_cpools_vr
    decomp_npools_vr          => clm3%g%l%c%cns%decomp_npools_vr

    leafc                     => clm3%g%l%c%p%pcs%leafc
    leafc_storage             => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                => clm3%g%l%c%p%pcs%leafc_xfer
    livestemc                 => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage         => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer            => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc                 => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage         => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer            => clm3%g%l%c%p%pcs%deadstemc_xfer
    frootc                    => clm3%g%l%c%p%pcs%frootc
    frootc_storage            => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer               => clm3%g%l%c%p%pcs%frootc_xfer
    livecrootc                => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage        => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer           => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc                => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage        => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer           => clm3%g%l%c%p%pcs%deadcrootc_xfer
    gresp_storage             => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                => clm3%g%l%c%p%pcs%gresp_xfer

    leafn                 => clm3%g%l%c%p%pns%leafn
    leafn_storage         => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer            => clm3%g%l%c%p%pns%leafn_xfer
    livestemn             => clm3%g%l%c%p%pns%livestemn
    livestemn_storage     => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer        => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn             => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage     => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer        => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                => clm3%g%l%c%p%pns%frootn
    frootn_storage        => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer           => clm3%g%l%c%p%pns%frootn_xfer
    livecrootn            => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage    => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer       => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn            => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage    => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer       => clm3%g%l%c%p%pns%deadcrootn_xfer
    retransn              => clm3%g%l%c%p%pns%retransn
    pactive               => clm3%g%l%c%p%active

    m_leafc_to_fire         => clm3%g%l%c%p%pcf%m_leafc_to_fire
    m_leafc_storage_to_fire => clm3%g%l%c%p%pcf%m_leafc_storage_to_fire
    m_leafc_xfer_to_fire    => clm3%g%l%c%p%pcf%m_leafc_xfer_to_fire
    m_livestemc_to_fire     => clm3%g%l%c%p%pcf%m_livestemc_to_fire
    m_livestemc_storage_to_fire => clm3%g%l%c%p%pcf%m_livestemc_storage_to_fire
    m_livestemc_xfer_to_fire => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_fire
    m_deadstemc_to_fire  => clm3%g%l%c%p%pcf%m_deadstemc_to_fire
    m_deadstemc_storage_to_fire => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_fire
    m_deadstemc_xfer_to_fire => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_fire
    m_frootc_to_fire => clm3%g%l%c%p%pcf%m_frootc_to_fire
    m_frootc_storage_to_fire => clm3%g%l%c%p%pcf%m_frootc_storage_to_fire
    m_frootc_xfer_to_fire  => clm3%g%l%c%p%pcf%m_frootc_xfer_to_fire
    m_livecrootc_to_fire   => clm3%g%l%c%p%pcf%m_livecrootc_to_fire
    m_livecrootc_storage_to_fire => &
            clm3%g%l%c%p%pcf%m_livecrootc_storage_to_fire
    m_livecrootc_xfer_to_fire  => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_fire
    m_deadcrootc_to_fire       => clm3%g%l%c%p%pcf%m_deadcrootc_to_fire
    m_deadcrootc_storage_to_fire => &
            clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_fire
    m_deadcrootc_xfer_to_fire  => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_fire
    m_gresp_storage_to_fire    => clm3%g%l%c%p%pcf%m_gresp_storage_to_fire
    m_gresp_xfer_to_fire       => clm3%g%l%c%p%pcf%m_gresp_xfer_to_fire

    m_leafn_to_fire            => clm3%g%l%c%p%pnf%m_leafn_to_fire
    m_leafn_storage_to_fire    => clm3%g%l%c%p%pnf%m_leafn_storage_to_fire
    m_leafn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire
    m_livestemn_to_fire        => clm3%g%l%c%p%pnf%m_livestemn_to_fire
    m_livestemn_storage_to_fire => clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire
    m_livestemn_xfer_to_fire    => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire
    m_deadstemn_to_fire         => clm3%g%l%c%p%pnf%m_deadstemn_to_fire
    m_deadstemn_storage_to_fire => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire
    m_deadstemn_xfer_to_fire    => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire
    m_frootn_to_fire            => clm3%g%l%c%p%pnf%m_frootn_to_fire
    m_frootn_storage_to_fire    => clm3%g%l%c%p%pnf%m_frootn_storage_to_fire
    m_frootn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire
    m_livecrootn_to_fire        => clm3%g%l%c%p%pnf%m_livecrootn_to_fire
    m_livecrootn_storage_to_fire=> clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire
    m_livecrootn_xfer_to_fire   => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire
    m_deadcrootn_to_fire        => clm3%g%l%c%p%pnf%m_deadcrootn_to_fire
    m_deadcrootn_storage_to_fire=> clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_xfer_to_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire
    m_retransn_to_fire          => clm3%g%l%c%p%pnf%m_retransn_to_fire

    m_leafc_to_litter_fire    => clm3%g%l%c%p%pcf%m_leafc_to_litter_fire
    m_leafc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_leafc_storage_to_litter_fire
    m_leafc_xfer_to_litter_fire    => &
            clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter_fire
    m_livestemc_to_litter_fire  => clm3%g%l%c%p%pcf%m_livestemc_to_litter_fire
    m_livestemc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter_fire
    m_livestemc_xfer_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter_fire
    m_livestemc_to_deadstemc_fire => &
            clm3%g%l%c%p%pcf%m_livestemc_to_deadstemc_fire
    m_deadstemc_to_litter_fire => clm3%g%l%c%p%pcf%m_deadstemc_to_litter_fire
    m_deadstemc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter_fire
    m_deadstemc_xfer_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter_fire
    m_frootc_to_litter_fire => clm3%g%l%c%p%pcf%m_frootc_to_litter_fire
    m_frootc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_frootc_storage_to_litter_fire
    m_frootc_xfer_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter_fire
    m_livecrootc_to_litter_fire => clm3%g%l%c%p%pcf%m_livecrootc_to_litter_fire
    m_livecrootc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter_fire
    m_livecrootc_xfer_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter_fire
    m_livecrootc_to_deadcrootc_fire => &
            clm3%g%l%c%p%pcf%m_livecrootc_to_deadcrootc_fire
    m_deadcrootc_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_deadcrootc_to_litter_fire
    m_deadcrootc_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter_fire
    m_deadcrootc_xfer_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter_fire
    m_gresp_storage_to_litter_fire => &
            clm3%g%l%c%p%pcf%m_gresp_storage_to_litter_fire
    m_gresp_xfer_to_litter_fire => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter_fire
    m_decomp_cpools_to_fire_vr => clm3%g%l%c%ccf%m_decomp_cpools_to_fire_vr
    m_c_to_litr_met_fire       => clm3%g%l%c%ccf%m_c_to_litr_met_fire
    m_c_to_litr_cel_fire       => clm3%g%l%c%ccf%m_c_to_litr_cel_fire
    m_c_to_litr_lig_fire       => clm3%g%l%c%ccf%m_c_to_litr_lig_fire

    m_leafn_to_litter_fire     => clm3%g%l%c%p%pnf%m_leafn_to_litter_fire
    m_leafn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_leafn_storage_to_litter_fire
    m_leafn_xfer_to_litter_fire => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter_fire
    m_livestemn_to_litter_fire  => clm3%g%l%c%p%pnf%m_livestemn_to_litter_fire
    m_livestemn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter_fire
    m_livestemn_xfer_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter_fire
    m_livestemn_to_deadstemn_fire => &
            clm3%g%l%c%p%pnf%m_livestemn_to_deadstemn_fire
    m_deadstemn_to_litter_fire => clm3%g%l%c%p%pnf%m_deadstemn_to_litter_fire
    m_deadstemn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter_fire
    m_deadstemn_xfer_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter_fire
    m_frootn_to_litter_fire => clm3%g%l%c%p%pnf%m_frootn_to_litter_fire
    m_frootn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_frootn_storage_to_litter_fire
    m_frootn_xfer_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter_fire
    m_livecrootn_to_litter_fire => clm3%g%l%c%p%pnf%m_livecrootn_to_litter_fire
    m_livecrootn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter_fire
    m_livecrootn_xfer_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter_fire
    m_livecrootn_to_deadcrootn_fire => &
            clm3%g%l%c%p%pnf%m_livecrootn_to_deadcrootn_fire
    m_deadcrootn_to_litter_fire => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter_fire
    m_deadcrootn_storage_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter_fire
    m_deadcrootn_xfer_to_litter_fire => &
            clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter_fire
    m_retransn_to_litter_fire => clm3%g%l%c%p%pnf%m_retransn_to_litter_fire
    m_decomp_npools_to_fire_vr => clm3%g%l%c%cnf%m_decomp_npools_to_fire_vr
    m_n_to_litr_met_fire   => clm3%g%l%c%cnf%m_n_to_litr_met_fire
    m_n_to_litr_cel_fire   => clm3%g%l%c%cnf%m_n_to_litr_cel_fire
    m_n_to_litr_lig_fire   => clm3%g%l%c%cnf%m_n_to_litr_lig_fire

    is_cwd                 => decomp_cascade_con%is_cwd
    is_litter              => decomp_cascade_con%is_litter
    croot_prof             => clm3%g%l%c%p%pps%croot_prof
    stem_prof              => clm3%g%l%c%p%pps%stem_prof
    froot_prof             => clm3%g%l%c%p%pps%froot_prof
    leaf_prof              => clm3%g%l%c%p%pps%leaf_prof

    ! Get model step size
    ! calculate burned area fraction per sec
    dt = dtsrf
    !
    ! pft loop
    !
    do fp = 1 , num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! get the column-level fractional area burned for this timestep
      ! and convert to a rate per second
      ! For non-crop (bare-soil and natural vegetation)
      if ( ivt(p) < nc3crop ) then
#ifdef DYNPFT
        f = (fbac(c)-baf_crop(c))/ dt
#else
        f = (farea_burned(c))/ dt
#endif
      else
        f = baf_crop(c) / dt
      end if

      ! apply this rate to the pft state variables to get flux rates
      ! biomass burning
      ! carbon fluxes
      m_leafc_to_fire(p) = leafc(p) * f * cc_leaf(ivt(p))
      m_leafc_storage_to_fire(p) = leafc_storage(p) * f * cc_other(ivt(p))
      m_leafc_xfer_to_fire(p) = leafc_xfer(p) * f * cc_other(ivt(p))
      m_livestemc_to_fire(p) = livestemc(p) * f * cc_lstem(ivt(p))
      m_livestemc_storage_to_fire(p) = livestemc_storage(p) * &
              f * cc_other(ivt(p))
      m_livestemc_xfer_to_fire(p) = livestemc_xfer(p) * f * cc_other(ivt(p))
      m_deadstemc_to_fire(p) = deadstemc(p) * f * cc_dstem(ivt(p))
      m_deadstemc_storage_to_fire(p) = deadstemc_storage(p) * &
              f * cc_other(ivt(p))
      m_deadstemc_xfer_to_fire(p) = deadstemc_xfer(p) * f * cc_other(ivt(p))
      m_frootc_to_fire(p) =  frootc(p) * f * 0._rkx
      m_frootc_storage_to_fire(p) = frootc_storage(p) * f * cc_other(ivt(p))
      m_frootc_xfer_to_fire(p) = frootc_xfer(p) * f * cc_other(ivt(p))
      m_livecrootc_to_fire(p) = livecrootc(p) * f * 0._rkx
      m_livecrootc_storage_to_fire(p) = livecrootc_storage(p) * &
              f * cc_other(ivt(p))
      m_livecrootc_xfer_to_fire(p) = livecrootc_xfer(p) * f * cc_other(ivt(p))
      m_deadcrootc_to_fire(p) =  deadcrootc(p) * f * 0._rkx
      m_deadcrootc_storage_to_fire(p) = deadcrootc_storage(p) * &
              f *  cc_other(ivt(p))
      m_deadcrootc_xfer_to_fire(p) = deadcrootc_xfer(p) * f * cc_other(ivt(p))
      m_gresp_storage_to_fire(p) = gresp_storage(p) * f * cc_other(ivt(p))
      m_gresp_xfer_to_fire(p) = gresp_xfer(p) * f * cc_other(ivt(p))

      ! nitrogen fluxes
      m_leafn_to_fire(p) = leafn(p) * f * cc_leaf(ivt(p))
      m_leafn_storage_to_fire(p) = leafn_storage(p) * f * cc_other(ivt(p))
      m_leafn_xfer_to_fire(p) = leafn_xfer(p) * f * cc_other(ivt(p))
      m_livestemn_to_fire(p) = livestemn(p) * f * cc_lstem(ivt(p))
      m_livestemn_storage_to_fire(p) = livestemn_storage(p) * &
              f * cc_other(ivt(p))
      m_livestemn_xfer_to_fire(p) = livestemn_xfer(p) * f * cc_other(ivt(p))
      m_deadstemn_to_fire(p) = deadstemn(p) * f * cc_dstem(ivt(p))
      m_deadstemn_storage_to_fire(p) = deadstemn_storage(p) * &
              f * cc_other(ivt(p))
      m_deadstemn_xfer_to_fire(p) = deadstemn_xfer(p) * f * cc_other(ivt(p))
      m_frootn_to_fire(p) = frootn(p) * f * 0._rkx
      m_frootn_storage_to_fire(p) = frootn_storage(p) * f * cc_other(ivt(p))
      m_frootn_xfer_to_fire(p) = frootn_xfer(p) * f * cc_other(ivt(p))
      m_livecrootn_to_fire(p) = livecrootn(p) * f * 0._rkx
      m_livecrootn_storage_to_fire(p) = livecrootn_storage(p) * &
              f * cc_other(ivt(p))
      m_livecrootn_xfer_to_fire(p) = livecrootn_xfer(p) * f * cc_other(ivt(p))
      m_deadcrootn_to_fire(p) = deadcrootn(p) * f * 0._rkx
      m_deadcrootn_xfer_to_fire(p) = deadcrootn_xfer(p) * f * cc_other(ivt(p))
      m_deadcrootn_storage_to_fire(p) = deadcrootn_storage(p) * &
              f * cc_other(ivt(p))
      m_retransn_to_fire(p) = retransn(p) * f * cc_other(ivt(p))

      ! mortality due to fire
      ! carbon bool
      m_leafc_to_litter_fire(p) = leafc(p) * f * &
                 (1._rkx - cc_leaf(ivt(p))) * fm_leaf(ivt(p))
      m_leafc_storage_to_litter_fire(p) = leafc_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_leafc_xfer_to_litter_fire(p) = leafc_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemc_to_litter_fire(p) = livestemc(p) * f * &
                 (1._rkx - cc_lstem(ivt(p))) * fm_droot(ivt(p))
      m_livestemc_storage_to_litter_fire(p) = livestemc_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemc_xfer_to_litter_fire(p) = livestemc_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemc_to_deadstemc_fire(p) = livestemc(p) * f * &
                 (1._rkx - cc_lstem(ivt(p))) * (fm_lstem(ivt(p))-fm_droot(ivt(p)))
      m_deadstemc_to_litter_fire(p) = deadstemc(p) * f * &
                 (1._rkx - cc_dstem(ivt(p))) * fm_droot(ivt(p))
      m_deadstemc_storage_to_litter_fire(p) = deadstemc_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_deadstemc_xfer_to_litter_fire(p) = deadstemc_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_frootc_to_litter_fire(p) = frootc(p) * f * fm_root(ivt(p))
      m_frootc_storage_to_litter_fire(p) = frootc_storage(p) * f * &
                       fm_other(ivt(p))
      m_frootc_xfer_to_litter_fire(p) = frootc_xfer(p) * f * fm_other(ivt(p))
      m_livecrootc_to_litter_fire(p) = livecrootc(p) * f * fm_droot(ivt(p))
      m_livecrootc_storage_to_litter_fire(p) = livecrootc_storage(p) * f * &
                       fm_other(ivt(p))
      m_livecrootc_xfer_to_litter_fire(p) = livecrootc_xfer(p) * f * &
                       fm_other(ivt(p))
      m_livecrootc_to_deadcrootc_fire(p) = livecrootc(p) * f * &
                       (fm_lroot(ivt(p))-fm_droot(ivt(p)))
      m_deadcrootc_to_litter_fire(p) = deadcrootc(p) * f * fm_droot(ivt(p))
      m_deadcrootc_storage_to_litter_fire(p) = deadcrootc_storage(p) * f * &
                       fm_other(ivt(p))
      m_deadcrootc_xfer_to_litter_fire(p) = deadcrootc_xfer(p) * f * &
                       fm_other(ivt(p))
      m_gresp_storage_to_litter_fire(p) = gresp_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_gresp_xfer_to_litter_fire(p) = gresp_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))

      ! nitrogen bool
      m_leafn_to_litter_fire(p) = leafn(p) * f * &
                 (1._rkx - cc_leaf(ivt(p))) * fm_leaf(ivt(p))
      m_leafn_storage_to_litter_fire(p) = leafn_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_leafn_xfer_to_litter_fire(p) = leafn_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemn_to_litter_fire(p) = livestemn(p) * f * &
                 (1._rkx - cc_lstem(ivt(p))) * fm_droot(ivt(p))
      m_livestemn_storage_to_litter_fire(p) = livestemn_storage(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemn_xfer_to_litter_fire(p) = livestemn_xfer(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))
      m_livestemn_to_deadstemn_fire(p) = livestemn(p) * f * &
                 (1._rkx - cc_lstem(ivt(p))) * (fm_lstem(ivt(p))-fm_droot(ivt(p)))
      m_frootn_to_litter_fire(p) = frootn(p) * f * fm_root(ivt(p))
      m_frootn_storage_to_litter_fire(p) = frootn_storage(p) * f * &
                      fm_other(ivt(p))
      m_frootn_xfer_to_litter_fire(p) = frootn_xfer(p) * f * fm_other(ivt(p))
      m_livecrootn_to_litter_fire(p) = livecrootn(p) * f * fm_droot(ivt(p))
      m_livecrootn_storage_to_litter_fire(p) = livecrootn_storage(p) * f * &
                      fm_other(ivt(p))
      m_livecrootn_xfer_to_litter_fire(p) = livecrootn_xfer(p) * f * &
                      fm_other(ivt(p))
      m_livecrootn_to_deadcrootn_fire(p) = livecrootn(p) * f * &
                      (fm_lroot(ivt(p))-fm_droot(ivt(p)))
      m_deadcrootn_to_litter_fire(p) = deadcrootn(p) * f * fm_droot(ivt(p))
      m_deadcrootn_storage_to_litter_fire(p) = deadcrootn_storage(p) * f * &
                      fm_other(ivt(p))
      m_deadcrootn_xfer_to_litter_fire(p) = deadcrootn_xfer(p) * f * &
                      fm_other(ivt(p))
      m_retransn_to_litter_fire(p) = retransn(p) * f * &
                 (1._rkx - cc_other(ivt(p))) * fm_other(ivt(p))

#if (defined CNDV)
      if ( woody(ivt(p)) == 1._rkx ) then
        if ( livestemc(p)+deadstemc(p) > 0._rkx ) then
          nind(p) = nind(p)*(1._rkx-1._rkx*fm_droot(ivt(p))*f)
        else
          nind(p) = 0._rkx
        end if
      end if
      leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
      if ( ivt(p) == noveg ) leafcmax(p) = 0._rkx
#endif
    end do  ! end of pfts loop
    !
    ! fire-affected carbon to litter and cwd
    !
    do j = 1 , nlevdecomp
      do pi = 1 , max_pft_per_col
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if ( pactive(p) ) then
              fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                 m_deadstemc_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
              fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                 m_deadcrootc_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
              fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                 m_deadstemn_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
              fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                 m_deadcrootn_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)

              fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                 m_livestemc_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
              fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                 m_livecrootc_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
              fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                 m_livestemn_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
              fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                 m_livecrootn_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)

              m_c_to_litr_met_fire(c,j)=m_c_to_litr_met_fire(c,j) + &
                ((m_leafc_to_litter_fire(p)*lf_flab(ivt(p)) &
                +m_leafc_storage_to_litter_fire(p) + &
                m_leafc_xfer_to_litter_fire(p) + &
                m_gresp_storage_to_litter_fire(p) &
                +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                (m_frootc_to_litter_fire(p)*fr_flab(ivt(p)) &
                +m_frootc_storage_to_litter_fire(p) + &
                m_frootc_xfer_to_litter_fire(p))*froot_prof(p,j) &
                +(m_livestemc_storage_to_litter_fire(p) + &
                m_livestemc_xfer_to_litter_fire(p) &
                +m_deadstemc_storage_to_litter_fire(p) + &
                m_deadstemc_xfer_to_litter_fire(p))* stem_prof(p,j)&
                +(m_livecrootc_storage_to_litter_fire(p) + &
                m_livecrootc_xfer_to_litter_fire(p) &
                +m_deadcrootc_storage_to_litter_fire(p) + &
                m_deadcrootc_xfer_to_litter_fire(p))*croot_prof(p,j))* wtcol(p)
                m_c_to_litr_cel_fire(c,j)=m_c_to_litr_cel_fire(c,j) + &
                (m_leafc_to_litter_fire(p)*lf_fcel(ivt(p))*leaf_prof(p,j) + &
                m_frootc_to_litter_fire(p)*&
                fr_fcel(ivt(p))*froot_prof(p,j))* wtcol(p)
                m_c_to_litr_lig_fire(c,j)=m_c_to_litr_lig_fire(c,j) + &
                (m_leafc_to_litter_fire(p)*lf_flig(ivt(p))*leaf_prof(p,j) + &
                m_frootc_to_litter_fire(p)* &
                fr_flig(ivt(p))*froot_prof(p,j))* wtcol(p)

              m_n_to_litr_met_fire(c,j)=m_n_to_litr_met_fire(c,j) + &
                ((m_leafn_to_litter_fire(p)*lf_flab(ivt(p)) &
                +m_leafn_storage_to_litter_fire(p) + &
                m_leafn_xfer_to_litter_fire(p)+m_retransn_to_litter_fire(p)) &
                *leaf_prof(p,j) +(m_frootn_to_litter_fire(p)*fr_flab(ivt(p)) &
                +m_frootn_storage_to_litter_fire(p) + &
                m_frootn_xfer_to_litter_fire(p))*froot_prof(p,j) &
                +(m_livestemn_storage_to_litter_fire(p) + &
                m_livestemn_xfer_to_litter_fire(p) &
                +m_deadstemn_storage_to_litter_fire(p) + &
                m_deadstemn_xfer_to_litter_fire(p))* stem_prof(p,j)&
                +(m_livecrootn_storage_to_litter_fire(p) + &
                m_livecrootn_xfer_to_litter_fire(p) &
                +m_deadcrootn_storage_to_litter_fire(p) + &
                m_deadcrootn_xfer_to_litter_fire(p))* croot_prof(p,j))* wtcol(p)
                m_n_to_litr_cel_fire(c,j)=m_n_to_litr_cel_fire(c,j) + &
                (m_leafn_to_litter_fire(p)*lf_fcel(ivt(p))*leaf_prof(p,j) + &
                m_frootn_to_litter_fire(p)*&
                fr_fcel(ivt(p))*froot_prof(p,j))* wtcol(p)
                m_n_to_litr_lig_fire(c,j)=m_n_to_litr_lig_fire(c,j) + &
                (m_leafn_to_litter_fire(p)*lf_flig(ivt(p))*leaf_prof(p,j) + &
                m_frootn_to_litter_fire(p)*fr_flig(ivt(p))*&
                froot_prof(p,j))* wtcol(p)
            end if
          end if
        end do
      end do
    end do
    !
    ! vertically-resolved decomposing C/N fire loss
    ! column loop
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      f = farea_burned(c) / dt

      ! apply this rate to the column state variables to get flux rates

      do j = 1 , nlevdecomp
        ! carbon fluxes
        do l = 1 , ndecomp_pools
          if ( is_litter(l) ) then
            m_decomp_cpools_to_fire_vr(c,j,l) = &
                    decomp_cpools_vr(c,j,l) * f * 0.4_rkx
          end if
          if ( is_cwd(l) ) then
            m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                                               (f-baf_crop(c)/dt) * 0.2_rkx
          end if
        end do

        ! nitrogen fluxes
        do l = 1 , ndecomp_pools
          if ( is_litter(l) ) then
            m_decomp_npools_to_fire_vr(c,j,l) = &
                    decomp_npools_vr(c,j,l) * f * 0.4_rkx
          end if
          if ( is_cwd(l) ) then
            m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                                           (f-baf_crop(c)/ dt) * 0.2_rkx
          end if
        end do
      end do
    end do  ! end of column loop
    !
    ! carbon loss due to deforestation fires
    !
#ifdef DYNPFT
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      lfc2(c)=0._rkx
      if ( date_is(idatex,1,1) .and. time_is(idatex,0,dtsrf) ) then
        if ( trotr1_col(c)+trotr2_col(c) > 0.6_rkx .and. &
             dtrotr_col(c) > 0._rkx .and. &
             lfc(c) > 0._rkx .and. fbac1(c) == 0._rkx ) then
          lfc2(c) = max(0._rkx,min(lfc(c),(farea_burned(c)-baf_crop(c) - &
                    baf_peatf(c))/2.0))/(dtrotr_col(c)*dayspy*secspday/dt)
          lfc(c)  = lfc(c)-max(0._rkx,min(lfc(c),&
                  (farea_burned(c)-baf_crop(c) - baf_peatf(c))/2.0_rkx))
        end if
      end if
    end do
#endif
    !
    ! Carbon loss due to peat fires
    !
    ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
    ! soil carbon b/c clm4 soil carbon was very low in peatland areas
    ! Fang Li has not checked clm45 soil carbon in peatland areas
    !
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)
      if ( latdeg(g) < borealat) then
        somc_fire(c)= totsomc(c)*baf_peatf(c)/dt*6.0_rkx/33.9_rkx
      else
        somc_fire(c)= baf_peatf(c)/dt*2.2e3_rkx
      end if
    end do

    ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
    ! They will be added here in proportion to the carbon emission
    ! Emission factors differ for various fire types
  end subroutine CNFireFluxes
  !
  ! Initialize data stream information for population density.
  !
  subroutine hdm_init( begg, endg )
    implicit none
    integer(ik4) , intent(in) :: begg , endg ! gridcell index bounds
    integer(ik4) :: yr , mon , day , ih

    allocate( forc_hdm(begg:endg) )
    allocate( hdm_p1(begg:endg) )
    allocate( hdm_p2(begg:endg) )

    call split_idate(idatex,yr,mon,day,ih)
    ipoprec = yr - 1850 + 1
    call clm_openfile(fsurdat,sdat_hdm)
    call clm_readvar(sdat_hdm,'HDM',hdm_p1,gcomm_gridcell,ipoprec)
    call clm_readvar(sdat_hdm,'HDM',hdm_p2,gcomm_gridcell,ipoprec+1)

    ! Add history fields
    call hist_addfld1d (fname='HDM', units='counts/km^2',      &
          avgflag='A', long_name='human population density',   &
          ptr_lnd=forc_hdm, default='inactive')

  end subroutine hdm_init
  !
  ! Interpolate data stream information for population density.
  !
  subroutine hdm_interp( )
    implicit none
    integer(ik4) :: yr , mon , day , ih , ip
    real(rkx) :: w1 , w2

    call split_idate(idatex,yr,mon,day,ih)
    ip = yr - 1850 + 1
    if ( ip /= ipoprec ) then
      ipoprec = ipoprec+1
      call clm_readvar(sdat_hdm,'HDM',hdm_p1,gcomm_gridcell,ipoprec)
      call clm_readvar(sdat_hdm,'HDM',hdm_p2,gcomm_gridcell,ipoprec+1)
    end if
    w1 = d_one - (yeardayfrac(idatex)/dayspy)
    w2 = d_one - w1
    forc_hdm(:) = hdm_p1(:)*w1 + hdm_p2(:)*w2
  end subroutine hdm_interp
  !
  ! Initialize data stream information for Lightning.
  !
  subroutine lnfm_init( begg, endg )
    implicit none
    integer(ik4) , intent(in) :: begg , endg   ! gridcell index bounds
    integer(ik4) :: yr , mon , day , ih

    allocate( forc_lnfm(begg:endg) )
    allocate( lnfm_p1(begg:endg) )
    allocate( lnfm_p2(begg:endg) )

    call split_idate(idatex,yr,mon,day,ih)
    ilnfmrec = int(yeardayfrac(idate1)/dayspy*365.0_rkx*8.0_rkx+1.0_rkx)+ih/3

    call clm_openfile(fsurdat,sdat_lnfm)
    call clm_readvar(sdat_hdm,'LNFM',lnfm_p1,gcomm_gridcell,ilnfmrec)
    if ( ilnfmrec+1 > 365*8 ) then
      call clm_readvar(sdat_hdm,'LNFM',lnfm_p2,gcomm_gridcell,1)
    else
      call clm_readvar(sdat_hdm,'LNFM',lnfm_p2,gcomm_gridcell,ilnfmrec+1)
    end if

    ! Add history fields
    call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
         avgflag='A', long_name='Lightning frequency',        &
         ptr_lnd=forc_lnfm, default='inactive')

  end subroutine lnfm_init
  !
  ! Interpolate data stream information for Lightning.
  !
  subroutine lnfm_interp( )
    implicit none
    integer(ik4) :: ip
    real(rkx) :: w1 , w2
    integer(ik4) :: yr , mon , day , ih

    call split_idate(idatex,yr,mon,day,ih)
    ip = int(yeardayfrac(idatex)/dayspy*365.0_rkx*8.0_rkx+1.0_rkx)+ih/3
    if ( ip /= ilnfmrec ) then
      ilnfmrec = ilnfmrec+1
      if ( ilnfmrec > 365*8 ) ilnfmrec = 1
      call clm_readvar(sdat_hdm,'LNFM',lnfm_p1,gcomm_gridcell,ilnfmrec)
      if ( ilnfmrec+1 > 365*8 ) then
        call clm_readvar(sdat_hdm,'LNFM',lnfm_p2,gcomm_gridcell,1)
      else
        call clm_readvar(sdat_hdm,'LNFM',lnfm_p2,gcomm_gridcell,ilnfmrec+1)
      end if
    end if
    w1 = d_one - (yeardayfrac(idatex)-dble(ip/8))
    w2 = d_one - w1
    forc_lnfm(:) = lnfm_p1(:)*w1 + lnfm_p2(:)*w2
  end subroutine lnfm_interp
#endif

end module mod_clm_cnfire
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
