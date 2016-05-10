module mod_clm_cninitimevar
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save

  public :: CNiniTimeVar

  contains
  !
  ! Initializes time varying variables used only in
  ! coupled carbon-nitrogen mode (CN):
  !
  subroutine CNiniTimeVar()
#ifdef CN
    use mod_clm_type
    use mod_clm_atmlnd  , only: clm_a2l
    use mod_clm_varcon  , only: istsoil, zsoi
    use mod_clm_varpar  , only: nlevgrnd, nlevsoi, nlevdecomp, ndecomp_pools, nlevdecomp_full
    use mod_clm_varcon  , only: istcrop, c13ratio, c14ratio
    use mod_clm_varctl  , only: use_c13, use_c14
    use mod_clm_pftvarcon   , only: noveg
    use mod_clm_pftvarcon   , only: npcropmin
    use mod_clm_decomp   , only: get_proc_bounds
    use mod_clm_surfrd   , only: crop_prog
!
    implicit none
!
    real(rkx), pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
    real(rkx), pointer :: woody(:)     ! binary flag for woody lifeform (1=woody, 0=not woody)
    real(rkx), pointer :: leafcn(:)    ! leaf C:N (gC/gN)
    real(rkx), pointer :: deadwdcn(:)  ! dead wood (xylem and heartwood) C:N (gC/gN)
    integer(ik4) , pointer :: ivt(:)       ! pft vegetation type
    logical , pointer :: lakpoi(:)    ! true => landunit is a lake point
    integer(ik4) , pointer :: plandunit(:) ! landunit index associated with each pft
    integer(ik4) , pointer :: clandunit(:) ! landunit index associated with each column
    integer(ik4) , pointer :: itypelun(:)  ! landunit type
!
! local pointers to implicit out arguments
!
    real(rkx), pointer :: forc_hgt_u_pft(:)    !observational height of wind at pft-level [m]
    real(rkx), pointer :: annsum_counter(:) ! seconds since last annual accumulator turnover
    real(rkx), pointer :: cannsum_npp(:)    ! annual sum of NPP, averaged from pft-level (gC/m2/yr)
    real(rkx), pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
    real(rkx), pointer :: sminn(:)              ! (gN/m2) soil mineral N
    real(rkx), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
    real(rkx), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
    real(rkx), pointer :: decomp_cpools(:,:)         ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_cpools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(rkx), pointer :: decomp_npools(:,:)         ! (gC/m2)  decomposing (litter, cwd, soil) N pools
    real(rkx), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rkx), pointer :: decomp_npools_1m(:,:)           ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
    real(rkx), pointer :: sminn_vr(:,:)              ! (gN/m3) soil mineral N
    real(rkx), pointer :: col_ctrunc(:)              ! (gC/m2) column-level sink for C truncation (diagnostic)
    real(rkx), pointer :: col_ctrunc_vr(:,:)         ! (gC/m3) column-level sink for C truncation (prognostic)
    real(rkx), pointer :: col_ntrunc_vr(:,:)         ! (gN/m3) column-level sink for N truncation
    real(rkx), pointer :: nfixation_prof(:,:)     ! (1/m) profile for N fixation additions
    real(rkx), pointer :: ndep_prof(:,:)     ! (1/m) profile for N fixation additions

    real(rkx), pointer :: fpi_vr(:,:)
    real(rkx), pointer :: alt(:)
    real(rkx), pointer :: altmax(:)
    real(rkx), pointer :: altmax_lastyear(:)
    integer(ik4), pointer :: alt_indx(:)
    integer(ik4), pointer :: altmax_indx(:)
    integer(ik4), pointer :: altmax_lastyear_indx(:)
    real(rkx), pointer :: som_adv_coef(:,:)
    real(rkx), pointer :: som_diffus_coef(:,:)
#ifdef NITRIF_DENITRIF
    real(rkx), pointer :: smin_nh4_vr(:,:)           ! (gN/m3) soil mineral NH4 pool
    real(rkx), pointer :: smin_no3_vr(:,:)           ! (gN/m3) soil mineral NO3 pool
    real(rkx), pointer :: smin_nh4(:)           ! (gN/m2) soil mineral NH4 pool
    real(rkx), pointer :: smin_no3(:)           ! (gN/m2) soil mineral NO3 pool
#endif

    real(rkx), pointer :: leafc(:)              ! (gC/m2) leaf C
    real(rkx), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
    real(rkx), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
    real(rkx), pointer :: grainc(:)             ! (gC/m2) grain C
    real(rkx), pointer :: grainc_storage(:)     ! (gC/m2) grain C storage
    real(rkx), pointer :: grainc_xfer(:)        ! (gC/m2) grain C transfer
    real(rkx), pointer :: frootc(:)             ! (gC/m2) fine root C
    real(rkx), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
    real(rkx), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
    real(rkx), pointer :: livestemc(:)          ! (gC/m2) live stem C
    real(rkx), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
    real(rkx), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
    real(rkx), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
    real(rkx), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
    real(rkx), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rkx), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
    real(rkx), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
    real(rkx), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
    real(rkx), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
    real(rkx), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
    real(rkx), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
    real(rkx), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
    real(rkx), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
    real(rkx), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
    real(rkx), pointer :: leafn(:)              ! (gN/m2) leaf N
    real(rkx), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
    real(rkx), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
    real(rkx), pointer :: grainn(:)             ! (gN/m2) grain N
    real(rkx), pointer :: grainn_storage(:)     ! (gN/m2) grain N storage
    real(rkx), pointer :: grainn_xfer(:)        ! (gN/m2) grain N transfer
    real(rkx), pointer :: frootn(:)             ! (gN/m2) fine root N
    real(rkx), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
    real(rkx), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
    real(rkx), pointer :: livestemn(:)          ! (gN/m2) live stem N
    real(rkx), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rkx), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rkx), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
    real(rkx), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
    real(rkx), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
    real(rkx), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
    real(rkx), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
    real(rkx), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
    real(rkx), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rkx), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
    real(rkx), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
    real(rkx), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
    real(rkx), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
    real(rkx), pointer :: psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: c13_psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: c13_psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: c14_psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: c14_psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: laisun(:)             ! sunlit projected leaf area index
    real(rkx), pointer :: laisha(:)             ! shaded projected leaf area index
    real(rkx), pointer :: dormant_flag(:)       ! dormancy flag
    real(rkx), pointer :: days_active(:)        ! number of days since last dormancy
    real(rkx), pointer :: onset_flag(:)         ! onset flag
    real(rkx), pointer :: onset_counter(:)      ! onset days counter
    real(rkx), pointer :: onset_gddflag(:)      ! onset flag for growing degree day sum
    real(rkx), pointer :: onset_fdd(:)          ! onset freezing degree days counter
    real(rkx), pointer :: onset_gdd(:)          ! onset growing degree days
    real(rkx), pointer :: onset_swi(:)          ! onset soil water index
    real(rkx), pointer :: offset_flag(:)        ! offset flag
    real(rkx), pointer :: offset_counter(:)     ! offset days counter
    real(rkx), pointer :: offset_fdd(:)         ! offset freezing degree days counter
    real(rkx), pointer :: offset_swi(:)         ! offset soil water index
    real(rkx), pointer :: lgsf(:)               ! long growing season factor [0-1]
    real(rkx), pointer :: bglfr(:)              ! background litterfall rate (1/s)
    real(rkx), pointer :: bgtr(:)               ! background transfer rate (1/s)
    real(rkx), pointer :: dayl(:)               ! daylength (seconds)
    real(rkx), pointer :: prev_dayl(:)          ! daylength from previous timestep (seconds)
    real(rkx), pointer :: annavg_t2m(:)         ! annual average 2m air temperature (K)
    real(rkx), pointer :: tempavg_t2m(:)        ! temporary average 2m air temperature (K)
    real(rkx), pointer :: gpp(:)                ! GPP flux before downregulation (gC/m2/s)
    real(rkx), pointer :: availc(:)             ! C flux available for allocation (gC/m2/s)
    real(rkx), pointer :: xsmrpool_recover(:)   ! C flux assigned to recovery of negative cpool (gC/m2/s)
    real(rkx), pointer :: xsmrpool_c13ratio(:)  ! C flux assigned to recovery of negative cpool (gC/m2/s)
    real(rkx), pointer :: alloc_pnow(:)         ! fraction of current allocation to display as new growth (DIM)
    real(rkx), pointer :: c_allometry(:)        ! C allocation index (DIM)
    real(rkx), pointer :: n_allometry(:)        ! N allocation index (DIM)
    real(rkx), pointer :: plant_ndemand(:)      ! N flux required to support initial GPP (gN/m2/s)
    real(rkx), pointer :: tempsum_potential_gpp(:) ! temporary annual sum of plant_ndemand
    real(rkx), pointer :: annsum_potential_gpp(:)  ! annual sum of plant_ndemand
    real(rkx), pointer :: tempmax_retransn(:)   ! temporary max of retranslocated N pool (gN/m2)
    real(rkx), pointer :: annmax_retransn(:)    ! annual max of retranslocated N pool (gN/m2)
    real(rkx), pointer :: avail_retransn(:)     ! N flux available from retranslocation pool (gN/m2/s)
    real(rkx), pointer :: plant_nalloc(:)       ! total allocated N flux (gN/m2/s)
    real(rkx), pointer :: plant_calloc(:)       ! total allocated C flux (gC/m2/s)
    real(rkx), pointer :: excess_cflux(:)       ! C flux not allocated due to downregulation (gC/m2/s)
    real(rkx), pointer :: downreg(:)            ! fractional reduction in GPP due to N limitation (DIM)
    real(rkx), pointer :: tempsum_npp(:)        ! temporary annual sum of NPP
    real(rkx), pointer :: annsum_npp(:)         ! annual sum of NPP
#if (defined CNDV)
    real(rkx), pointer :: tempsum_litfall(:)    ! temporary annual sum of litfall
    real(rkx), pointer :: annsum_litfall(:)     ! annual sum of litfall
#endif
    real(rkx), pointer :: rc13_canair(:)        !C13O2/C12O2 in canopy air
    real(rkx), pointer :: rc13_psnsun(:)        !C13O2/C12O2 in sunlit canopy psn flux
    real(rkx), pointer :: rc13_psnsha(:)        !C13O2/C12O2 in shaded canopy psn flux
    real(rkx), pointer :: alphapsnsun(:)        !sunlit 13c fractionation ([])
    real(rkx), pointer :: alphapsnsha(:)        !shaded 13c fractionation ([])
    real(rkx), pointer :: qflx_drain(:)         ! sub-surface runoff (mm H2O /s)
    real(rkx), pointer :: qflx_surf(:)          ! surface runoff (mm H2O /s)
    real(rkx), pointer :: qflx_irrig(:)         !irrigation flux (mm H2O/s)

    ! fire-related variables changed by F. Li and S. Levis
    real(rkx), pointer :: wf(:)                 ! soil moisture in top 0.05 m
    real(rkx), pointer :: wf2(:)
    real(rkx), pointer :: nfire(:)              ! fire counts/km2/timestep
    real(rkx), pointer :: baf_crop(:)          ! burned area fraction in crop
    real(rkx), pointer :: baf_peatf(:)         ! burned area fraction in peatland
    real(rkx), pointer :: fbac(:)
    real(rkx), pointer :: fbac1(:)
    real(rkx), pointer :: farea_burned(:)       ! timestep fractional area burned (proportion)

    real(rkx), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
    real(rkx), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
    real(rkx), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
    real(rkx), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
    real(rkx), pointer :: totlitc_1m(:)         ! (gC/m2) total litter carbon to 1 meter
    real(rkx), pointer :: totsomc_1m(:)         ! (gC/m2) total soil organic matter carbon to 1 meter

    real(rkx), pointer :: woodc(:)              ! (gC/m2) pft-level wood C
    real(rkx), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
    real(rkx), pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg
    real(rkx), pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
    real(rkx), pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
    real(rkx), pointer :: totlitn_1m(:)         ! (gN/m2) total litter nitrogen to 1 meter
    real(rkx), pointer :: totsomn_1m(:)         ! (gN/m2) total soil organic matter nitrogen to 1 meter
    real(rkx), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
    real(rkx), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
    real(rkx), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
    real(rkx), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
    real(rkx), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
    real(rkx), pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
    real(rkx), pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
    real(rkx), pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
    real(rkx), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
    real(rkx), pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
    real(rkx), pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
    real(rkx), pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
    !!! C13
    real(rkx), pointer :: cwdc13(:)               ! (gC/m2) coarse woody debris C
    real(rkx), pointer :: decomp_c13pools(:,:)              ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_c13pools_vr(:,:,:)         ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_c13pools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(rkx), pointer :: c13_col_ctrunc_vr(:,:)            ! (gC/m3) C truncation term
    real(rkx), pointer :: leafc13(:)              ! (gC/m2) leaf C
    real(rkx), pointer :: leafc13_storage(:)      ! (gC/m2) leaf C storage
    real(rkx), pointer :: leafc13_xfer(:)         ! (gC/m2) leaf C transfer
    real(rkx), pointer :: frootc13(:)             ! (gC/m2) fine root C
    real(rkx), pointer :: frootc13_storage(:)     ! (gC/m2) fine root C storage
    real(rkx), pointer :: frootc13_xfer(:)        ! (gC/m2) fine root C transfer
    real(rkx), pointer :: livestemc13(:)          ! (gC/m2) live stem C
    real(rkx), pointer :: livestemc13_storage(:)  ! (gC/m2) live stem C storage
    real(rkx), pointer :: livestemc13_xfer(:)     ! (gC/m2) live stem C transfer
    real(rkx), pointer :: deadstemc13(:)          ! (gC/m2) dead stem C
    real(rkx), pointer :: deadstemc13_storage(:)  ! (gC/m2) dead stem C storage
    real(rkx), pointer :: deadstemc13_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rkx), pointer :: livecrootc13(:)         ! (gC/m2) live coarse root C
    real(rkx), pointer :: livecrootc13_storage(:) ! (gC/m2) live coarse root C storage
    real(rkx), pointer :: livecrootc13_xfer(:)    ! (gC/m2) live coarse root C transfer
    real(rkx), pointer :: deadcrootc13(:)         ! (gC/m2) dead coarse root C
    real(rkx), pointer :: deadcrootc13_storage(:) ! (gC/m2) dead coarse root C storage
    real(rkx), pointer :: deadcrootc13_xfer(:)    ! (gC/m2) dead coarse root C transfer
    real(rkx), pointer :: c13_gresp_storage(:)    ! (gC/m2) growth respiration storage
    real(rkx), pointer :: c13_gresp_xfer(:)       ! (gC/m2) growth respiration transfer
    real(rkx), pointer :: c13pool(:)              ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: c13xsmrpool(:)          ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: c13_pft_ctrunc(:)       ! (gC/m2) C truncation term
    real(rkx), pointer :: totvegc13(:)            ! (gC/m2) total vegetation carbon, excluding cpool

    !!! C14
    real(rkx), pointer :: cwdc14(:)               ! (gC/m2) coarse woody debris C
    real(rkx), pointer :: decomp_c14pools(:,:)              ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_c14pools_vr(:,:,:)         ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_c14pools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(rkx), pointer :: c14_col_ctrunc_vr(:,:)            ! (gC/m3) C truncation term
    real(rkx), pointer :: leafc14(:)              ! (gC/m2) leaf C
    real(rkx), pointer :: leafc14_storage(:)      ! (gC/m2) leaf C storage
    real(rkx), pointer :: leafc14_xfer(:)         ! (gC/m2) leaf C transfer
    real(rkx), pointer :: frootc14(:)             ! (gC/m2) fine root C
    real(rkx), pointer :: frootc14_storage(:)     ! (gC/m2) fine root C storage
    real(rkx), pointer :: frootc14_xfer(:)        ! (gC/m2) fine root C transfer
    real(rkx), pointer :: livestemc14(:)          ! (gC/m2) live stem C
    real(rkx), pointer :: livestemc14_storage(:)  ! (gC/m2) live stem C storage
    real(rkx), pointer :: livestemc14_xfer(:)     ! (gC/m2) live stem C transfer
    real(rkx), pointer :: deadstemc14(:)          ! (gC/m2) dead stem C
    real(rkx), pointer :: deadstemc14_storage(:)  ! (gC/m2) dead stem C storage
    real(rkx), pointer :: deadstemc14_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rkx), pointer :: livecrootc14(:)         ! (gC/m2) live coarse root C
    real(rkx), pointer :: livecrootc14_storage(:) ! (gC/m2) live coarse root C storage
    real(rkx), pointer :: livecrootc14_xfer(:)    ! (gC/m2) live coarse root C transfer
    real(rkx), pointer :: deadcrootc14(:)         ! (gC/m2) dead coarse root C
    real(rkx), pointer :: deadcrootc14_storage(:) ! (gC/m2) dead coarse root C storage
    real(rkx), pointer :: deadcrootc14_xfer(:)    ! (gC/m2) dead coarse root C transfer
    real(rkx), pointer :: c14_gresp_storage(:)    ! (gC/m2) growth respiration storage
    real(rkx), pointer :: c14_gresp_xfer(:)       ! (gC/m2) growth respiration transfer
    real(rkx), pointer :: c14pool(:)              ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: c14xsmrpool(:)          ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: c14_pft_ctrunc(:)       ! (gC/m2) C truncation term
    real(rkx), pointer :: totvegc14(:)            ! (gC/m2) total vegetation carbon, excluding cpool
    real(rkx), pointer :: rc14_atm(:)             !C14O2/C12O2 in atmosphere

    ! dynamic landuse variables
    real(rkx), pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
    real(rkx), pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
    real(rkx), pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
    real(rkx), pointer :: totprodc(:)           ! (gC/m2) total wood product C

    !!! C13
    real(rkx), pointer :: seedc13(:)              ! (gC/m2) column-level pool for seeding new PFTs
    real(rkx), pointer :: prod10c13(:)          ! (gC/m2) wood product C13 pool, 10-year lifespan
    real(rkx), pointer :: prod100c13(:)         ! (gC/m2) wood product C13 pool, 100-year lifespan
    real(rkx), pointer :: totprodc13(:)         ! (gC/m2) total wood product C13
    !!! C14
    real(rkx), pointer :: seedc14(:)              ! (gC/m2) column-level pool for seeding new PFTs
    real(rkx), pointer :: prod10c14(:)          ! (gC/m2) wood product C14 pool, 10-year lifespan
    real(rkx), pointer :: prod100c14(:)         ! (gC/m2) wood product C14 pool, 100-year lifespan
    real(rkx), pointer :: totprodc14(:)         ! (gC/m2) total wood product C14

    real(rkx), pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
    real(rkx), pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
    real(rkx), pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
    real(rkx), pointer :: totprodn(:)           ! (gN/m2) total wood product N
    real(rkx), pointer :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
    real(rkx), pointer :: initial_stock(:)                  ! initial concentration for seeding at spinup

    ! crop
    real(rkx), pointer :: fert_counter(:)
    real(rkx), pointer :: fert(:)
    real(rkx), pointer :: soyfixn(:)
    real(rkx), pointer :: grain_flag(:)
!
    integer(ik4) :: l,c,p,j,k   ! indices
    integer(ik4) :: begp, endp  ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc  ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl  ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg  ! per-proc gridcell ending gridcell indices

    ! assign local pointers at the gridcell level

    ! assign local pointers at the landunit level
    lakpoi                         => clm3%g%l%lakpoi
    itypelun                       => clm3%g%l%itype

    ! assign local pointers at the column level
    clandunit                      => clm3%g%l%c%landunit
    annsum_counter                 => clm3%g%l%c%cps%annsum_counter
    cannsum_npp                    => clm3%g%l%c%cps%cannsum_npp
    cannavg_t2m                    => clm3%g%l%c%cps%cannavg_t2m

     !fire related variables changed by F. Li and S. Levis
    wf                             => clm3%g%l%c%cps%wf
    wf2                            => clm3%g%l%c%cps%wf2
    nfire                          => clm3%g%l%c%cps%nfire
    baf_crop                       => clm3%g%l%c%cps%baf_crop
    baf_peatf                      => clm3%g%l%c%cps%baf_peatf
    fbac                           => clm3%g%l%c%cps%fbac
    fbac1                          => clm3%g%l%c%cps%fbac1
    farea_burned                   => clm3%g%l%c%cps%farea_burned

    qflx_drain                     => clm3%g%l%c%cwf%qflx_drain
    qflx_surf                      => clm3%g%l%c%cwf%qflx_surf
    decomp_cpools                  => clm3%g%l%c%ccs%decomp_cpools
    decomp_cpools_1m               => clm3%g%l%c%ccs%decomp_cpools_1m
    decomp_cpools_vr               => clm3%g%l%c%ccs%decomp_cpools_vr
    decomp_npools                  => clm3%g%l%c%cns%decomp_npools
    decomp_npools_vr               => clm3%g%l%c%cns%decomp_npools_vr
    decomp_npools_1m               => clm3%g%l%c%cns%decomp_npools_1m
    nfixation_prof                 => clm3%g%l%c%cps%nfixation_prof
    ndep_prof                      => clm3%g%l%c%cps%ndep_prof
    qflx_irrig                     => clm3%g%l%c%cwf%qflx_irrig

    ! dynamic landuse variables
    seedc                          => clm3%g%l%c%ccs%seedc
    prod10c                        => clm3%g%l%c%ccs%prod10c
    prod100c                       => clm3%g%l%c%ccs%prod100c
    totprodc                       => clm3%g%l%c%ccs%totprodc
    seedn                          => clm3%g%l%c%cns%seedn
    prod10n                        => clm3%g%l%c%cns%prod10n
    prod100n                       => clm3%g%l%c%cns%prod100n
    totprodn                       => clm3%g%l%c%cns%totprodn
    sminn                          => clm3%g%l%c%cns%sminn
    col_ctrunc                     => clm3%g%l%c%ccs%col_ctrunc
    sminn_vr                       => clm3%g%l%c%cns%sminn_vr
    col_ctrunc_vr                  => clm3%g%l%c%ccs%col_ctrunc_vr
    col_ntrunc_vr                  => clm3%g%l%c%cns%col_ntrunc_vr

    fpi_vr                         => clm3%g%l%c%cps%fpi_vr
    alt                            => clm3%g%l%c%cps%alt
    altmax                         => clm3%g%l%c%cps%altmax
    altmax_lastyear                => clm3%g%l%c%cps%altmax_lastyear
    som_adv_coef                   => clm3%g%l%c%cps%som_adv_coef
    som_diffus_coef                => clm3%g%l%c%cps%som_diffus_coef
    alt_indx                       => clm3%g%l%c%cps%alt_indx
    altmax_indx                    => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx           => clm3%g%l%c%cps%altmax_lastyear_indx
#ifdef NITRIF_DENITRIF
    smin_nh4_vr                    => clm3%g%l%c%cns%smin_nh4_vr
    smin_no3_vr                    => clm3%g%l%c%cns%smin_no3_vr
    smin_nh4                       => clm3%g%l%c%cns%smin_nh4
    smin_no3                       => clm3%g%l%c%cns%smin_no3
#endif

    totcolc                        => clm3%g%l%c%ccs%totcolc
    cwdc                           => clm3%g%l%c%ccs%cwdc
    totecosysc                     => clm3%g%l%c%ccs%totecosysc
    totlitc                        => clm3%g%l%c%ccs%totlitc
    totsomc                        => clm3%g%l%c%ccs%totsomc
    totlitc_1m                     => clm3%g%l%c%ccs%totlitc_1m
    totsomc_1m                     => clm3%g%l%c%ccs%totsomc_1m

    totcoln                        => clm3%g%l%c%cns%totcoln
    cwdn                           => clm3%g%l%c%cns%cwdn
    totecosysn                     => clm3%g%l%c%cns%totecosysn
    totlitn                        => clm3%g%l%c%cns%totlitn
    totsomn                        => clm3%g%l%c%cns%totsomn
    totlitn_1m                     => clm3%g%l%c%cns%totlitn_1m
    totsomn_1m                     => clm3%g%l%c%cns%totsomn_1m
    if ( use_c13 ) then
      seedc13                      => clm3%g%l%c%cc13s%seedc
      prod10c13                    => clm3%g%l%c%cc13s%prod10c
      prod100c13                   => clm3%g%l%c%cc13s%prod100c
      totprodc13                   => clm3%g%l%c%cc13s%totprodc
      cwdc13                       => clm3%g%l%c%cc13s%cwdc
      decomp_c13pools              => clm3%g%l%c%cc13s%decomp_cpools
      decomp_c13pools_vr           => clm3%g%l%c%cc13s%decomp_cpools_vr
      c13_col_ctrunc_vr            => clm3%g%l%c%cc13s%col_ctrunc_vr
      decomp_c13pools_1m           => clm3%g%l%c%cc13s%decomp_cpools_1m
      c13_psnsun                   => clm3%g%l%c%p%pc13f%psnsun
      c13_psnsha                   => clm3%g%l%c%p%pc13f%psnsha
      xsmrpool_c13ratio            => clm3%g%l%c%p%pepv%xsmrpool_c13ratio
      alphapsnsun                  => clm3%g%l%c%p%pps%alphapsnsun
      alphapsnsha                  => clm3%g%l%c%p%pps%alphapsnsha
      leafc13                      => clm3%g%l%c%p%pc13s%leafc
      leafc13_storage              => clm3%g%l%c%p%pc13s%leafc_storage
      leafc13_xfer                 => clm3%g%l%c%p%pc13s%leafc_xfer
      frootc13                     => clm3%g%l%c%p%pc13s%frootc
      frootc13_storage             => clm3%g%l%c%p%pc13s%frootc_storage
      frootc13_xfer                => clm3%g%l%c%p%pc13s%frootc_xfer
      livestemc13                  => clm3%g%l%c%p%pc13s%livestemc
      livestemc13_storage          => clm3%g%l%c%p%pc13s%livestemc_storage
      livestemc13_xfer             => clm3%g%l%c%p%pc13s%livestemc_xfer
      deadstemc13                  => clm3%g%l%c%p%pc13s%deadstemc
      deadstemc13_storage          => clm3%g%l%c%p%pc13s%deadstemc_storage
      deadstemc13_xfer             => clm3%g%l%c%p%pc13s%deadstemc_xfer
      livecrootc13                 => clm3%g%l%c%p%pc13s%livecrootc
      livecrootc13_storage         => clm3%g%l%c%p%pc13s%livecrootc_storage
      livecrootc13_xfer            => clm3%g%l%c%p%pc13s%livecrootc_xfer
      deadcrootc13                 => clm3%g%l%c%p%pc13s%deadcrootc
      deadcrootc13_storage         => clm3%g%l%c%p%pc13s%deadcrootc_storage
      deadcrootc13_xfer            => clm3%g%l%c%p%pc13s%deadcrootc_xfer
      c13_gresp_storage            => clm3%g%l%c%p%pc13s%gresp_storage
      c13_gresp_xfer               => clm3%g%l%c%p%pc13s%gresp_xfer
      c13pool                      => clm3%g%l%c%p%pc13s%cpool
      c13xsmrpool                  => clm3%g%l%c%p%pc13s%xsmrpool
      c13_pft_ctrunc               => clm3%g%l%c%p%pc13s%pft_ctrunc
      totvegc13                    => clm3%g%l%c%p%pc13s%totvegc
      rc13_canair                  => clm3%g%l%c%p%pepv%rc13_canair
      rc13_psnsun                  => clm3%g%l%c%p%pepv%rc13_psnsun
      rc13_psnsha                  => clm3%g%l%c%p%pepv%rc13_psnsha
    end if
    if ( use_c14 ) then
      seedc14                      => clm3%g%l%c%cc14s%seedc
      prod10c14                    => clm3%g%l%c%cc14s%prod10c
      prod100c14                   => clm3%g%l%c%cc14s%prod100c
      totprodc14                   => clm3%g%l%c%cc14s%totprodc
      cwdc14                       => clm3%g%l%c%cc14s%cwdc
      decomp_c14pools              => clm3%g%l%c%cc14s%decomp_cpools
      decomp_c14pools_vr           => clm3%g%l%c%cc14s%decomp_cpools_vr
      c14_col_ctrunc_vr            => clm3%g%l%c%cc14s%col_ctrunc_vr
      decomp_c14pools_1m           => clm3%g%l%c%cc14s%decomp_cpools_1m
      c14_psnsun                   => clm3%g%l%c%p%pc14f%psnsun
      c14_psnsha                   => clm3%g%l%c%p%pc14f%psnsha
      leafc14                      => clm3%g%l%c%p%pc14s%leafc
      leafc14_storage              => clm3%g%l%c%p%pc14s%leafc_storage
      leafc14_xfer                 => clm3%g%l%c%p%pc14s%leafc_xfer
      frootc14                     => clm3%g%l%c%p%pc14s%frootc
      frootc14_storage             => clm3%g%l%c%p%pc14s%frootc_storage
      frootc14_xfer                => clm3%g%l%c%p%pc14s%frootc_xfer
      livestemc14                  => clm3%g%l%c%p%pc14s%livestemc
      livestemc14_storage          => clm3%g%l%c%p%pc14s%livestemc_storage
      livestemc14_xfer             => clm3%g%l%c%p%pc14s%livestemc_xfer
      deadstemc14                  => clm3%g%l%c%p%pc14s%deadstemc
      deadstemc14_storage          => clm3%g%l%c%p%pc14s%deadstemc_storage
      deadstemc14_xfer             => clm3%g%l%c%p%pc14s%deadstemc_xfer
      livecrootc14                 => clm3%g%l%c%p%pc14s%livecrootc
      livecrootc14_storage         => clm3%g%l%c%p%pc14s%livecrootc_storage
      livecrootc14_xfer            => clm3%g%l%c%p%pc14s%livecrootc_xfer
      deadcrootc14                 => clm3%g%l%c%p%pc14s%deadcrootc
      deadcrootc14_storage         => clm3%g%l%c%p%pc14s%deadcrootc_storage
      deadcrootc14_xfer            => clm3%g%l%c%p%pc14s%deadcrootc_xfer
      c14_gresp_storage            => clm3%g%l%c%p%pc14s%gresp_storage
      c14_gresp_xfer               => clm3%g%l%c%p%pc14s%gresp_xfer
      c14pool                      => clm3%g%l%c%p%pc14s%cpool
      c14xsmrpool                  => clm3%g%l%c%p%pc14s%xsmrpool
      c14_pft_ctrunc               => clm3%g%l%c%p%pc14s%pft_ctrunc
      totvegc14                    => clm3%g%l%c%p%pc14s%totvegc
      rc14_atm                     => clm3%g%l%c%p%pepv%rc14_atm
    end if
    ! crop
    soyfixn                        => clm3%g%l%c%p%pnf%soyfixn
    fert                           => clm3%g%l%c%p%pnf%fert
    fert_counter                   => clm3%g%l%c%p%pepv%fert_counter
    grain_flag                     => clm3%g%l%c%p%pepv%grain_flag
    ! assign local pointers at the pft level
    ivt                            => clm3%g%l%c%p%itype
    plandunit                      => clm3%g%l%c%p%landunit
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    grainc                         => clm3%g%l%c%p%pcs%grainc
    grainc_storage                 => clm3%g%l%c%p%pcs%grainc_storage
    grainc_xfer                    => clm3%g%l%c%p%pcs%grainc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    cpool                          => clm3%g%l%c%p%pcs%cpool
    xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
    forc_hgt_u_pft                 => clm3%g%l%c%p%pps%forc_hgt_u_pft
    woodc                          => clm3%g%l%c%p%pcs%woodc
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    grainn                         => clm3%g%l%c%p%pns%grainn
    grainn_storage                 => clm3%g%l%c%p%pns%grainn_storage
    grainn_xfer                    => clm3%g%l%c%p%pns%grainn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn
    npool                          => clm3%g%l%c%p%pns%npool
    psnsun                         => clm3%g%l%c%p%pcf%psnsun
    psnsha                         => clm3%g%l%c%p%pcf%psnsha
    laisun                         => clm3%g%l%c%p%pps%laisun
    laisha                         => clm3%g%l%c%p%pps%laisha
    dormant_flag                   => clm3%g%l%c%p%pepv%dormant_flag
    days_active                    => clm3%g%l%c%p%pepv%days_active
    onset_flag                     => clm3%g%l%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%l%c%p%pepv%onset_counter
    onset_gddflag                  => clm3%g%l%c%p%pepv%onset_gddflag
    onset_fdd                      => clm3%g%l%c%p%pepv%onset_fdd
    onset_gdd                      => clm3%g%l%c%p%pepv%onset_gdd
    onset_swi                      => clm3%g%l%c%p%pepv%onset_swi
    offset_flag                    => clm3%g%l%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%l%c%p%pepv%offset_counter
    offset_fdd                     => clm3%g%l%c%p%pepv%offset_fdd
    offset_swi                     => clm3%g%l%c%p%pepv%offset_swi
    lgsf                           => clm3%g%l%c%p%pepv%lgsf
    bglfr                          => clm3%g%l%c%p%pepv%bglfr
    bgtr                           => clm3%g%l%c%p%pepv%bgtr
    dayl                           => clm3%g%l%c%p%pepv%dayl
    prev_dayl                      => clm3%g%l%c%p%pepv%prev_dayl
    annavg_t2m                     => clm3%g%l%c%p%pepv%annavg_t2m
    tempavg_t2m                    => clm3%g%l%c%p%pepv%tempavg_t2m
    gpp                            => clm3%g%l%c%p%pepv%gpp
    availc                         => clm3%g%l%c%p%pepv%availc
    xsmrpool_recover               => clm3%g%l%c%p%pepv%xsmrpool_recover
    alloc_pnow                     => clm3%g%l%c%p%pepv%alloc_pnow
    c_allometry                    => clm3%g%l%c%p%pepv%c_allometry
    n_allometry                    => clm3%g%l%c%p%pepv%n_allometry
    plant_ndemand                  => clm3%g%l%c%p%pepv%plant_ndemand
    tempsum_potential_gpp          => clm3%g%l%c%p%pepv%tempsum_potential_gpp
    annsum_potential_gpp           => clm3%g%l%c%p%pepv%annsum_potential_gpp
    tempmax_retransn               => clm3%g%l%c%p%pepv%tempmax_retransn
    annmax_retransn                => clm3%g%l%c%p%pepv%annmax_retransn
    avail_retransn                 => clm3%g%l%c%p%pepv%avail_retransn
    plant_nalloc                   => clm3%g%l%c%p%pepv%plant_nalloc
    plant_calloc                   => clm3%g%l%c%p%pepv%plant_calloc
    excess_cflux                   => clm3%g%l%c%p%pepv%excess_cflux
    downreg                        => clm3%g%l%c%p%pepv%downreg
    tempsum_npp                    => clm3%g%l%c%p%pepv%tempsum_npp
    annsum_npp                     => clm3%g%l%c%p%pepv%annsum_npp
#if (defined CNDV)
    tempsum_litfall                => clm3%g%l%c%p%pepv%tempsum_litfall
    annsum_litfall                 => clm3%g%l%c%p%pepv%annsum_litfall
#endif
    dispvegc                       => clm3%g%l%c%p%pcs%dispvegc
    pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
    storvegc                       => clm3%g%l%c%p%pcs%storvegc
    totpftc                        => clm3%g%l%c%p%pcs%totpftc
    totvegc                        => clm3%g%l%c%p%pcs%totvegc
    prev_frootc_to_litter          => clm3%g%l%c%p%pepv%prev_frootc_to_litter
    prev_leafc_to_litter           => clm3%g%l%c%p%pepv%prev_leafc_to_litter
    dispvegn                       => clm3%g%l%c%p%pns%dispvegn
    pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
    storvegn                       => clm3%g%l%c%p%pns%storvegn
    totpftn                        => clm3%g%l%c%p%pns%totpftn
    totvegn                        => clm3%g%l%c%p%pns%totvegn

    ! assign local pointers for ecophysiological constants
    evergreen                      => pftcon%evergreen
    woody                          => pftcon%woody
    leafcn                         => pftcon%leafcn
    deadwdcn                       => pftcon%deadwdcn
    ! decomposoition parameters
    initial_cn_ratio    => decomp_cascade_con%initial_cn_ratio
    initial_stock       => decomp_cascade_con%initial_stock

    ! Determine subgrid bounds on this processor
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to CNVegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information
    ! to set this properly (e.g., pft-level displacement height and roughness
    ! length). So leave at 30m.
    do p = begp, endp
      forc_hgt_u_pft(p) = 30._rkx
    end do

    ! initialize column-level variables
    do c = begc, endc
      l = clandunit(c)
      if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then

        ! column physical state variables
        annsum_counter(c) = 0._rkx
        cannsum_npp(c)    = 0._rkx
        cannavg_t2m(c)    = 280._rkx
        ! fire related variables changed by F. Li and S. Levis
        ! it needs to be non zero so the first time step has no fires
        wf(c) = 1.0_rkx
        wf2(c) = 1.0_rkx
        nfire(c) = 0._rkx
        baf_crop(c) = 0._rkx
        baf_peatf(c) = 0._rkx
        fbac(c) = 0._rkx
        fbac1(c) = 0._rkx
        farea_burned(c) = 0._rkx

        ! needed for CNNLeaching
        qflx_drain(c) = 0._rkx
        qflx_surf(c) = 0._rkx
        qflx_irrig(c) = 0._rkx

        ! column carbon state variable initialization
        do j = 1, nlevdecomp
          do k = 1, ndecomp_pools
            if (zsoi(j) < 0.3 ) then  !! only initialize upper soil column
              decomp_cpools_vr(c,j,k) = initial_stock(k)
            else
              decomp_cpools_vr(c,j,k) = 0._rkx
            end if
          end do
          col_ctrunc_vr(c,j) = 0._rkx
        end do
        if ( nlevdecomp > 1 ) then
          do j = nlevdecomp+1, nlevdecomp_full
            do k = 1, ndecomp_pools
              decomp_cpools_vr(c,j,k) = 0._rkx
            end do
            col_ctrunc_vr(c,j) = 0._rkx
          end do
        end if
        do k = 1, ndecomp_pools
          decomp_cpools(c,k) = initial_stock(k)
          decomp_cpools_1m(c,k) = initial_stock(k)
        end do

        do j = 1, nlevdecomp_full
          ! initialize fpi_vr so that levels below nlevsoi are not nans
          fpi_vr(c,j) = 0._rkx
          som_adv_coef(c,j) = 0._rkx
          som_diffus_coef(c,j) = 0._rkx
          ! here initialize the profiles for converting to
          ! vertically resolved carbon pools
          nfixation_prof(c,j) = 0._rkx
          ndep_prof(c,j) = 0._rkx
        end do

        ! and define alt variables to be zero
        alt(c) = 0._rkx
        altmax(c) = 0._rkx
        altmax_lastyear(c) = 0._rkx
        alt_indx(c) = 0
        altmax_indx(c) = 0
        altmax_lastyear_indx = 0

        cwdc(c)       = 0._rkx
        col_ctrunc(c) = 0._rkx
        totlitc(c)    = 0._rkx
        totsomc(c)    = 0._rkx
        totlitc_1m(c) = 0._rkx
        totsomc_1m(c) = 0._rkx
        totecosysc(c) = 0._rkx
        totcolc(c)    = 0._rkx

        if ( use_c13 ) then
          do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
              decomp_c13pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c13ratio
            end do
            c13_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c13ratio
          end do
          if ( nlevdecomp > 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
              do k = 1, ndecomp_pools
                decomp_c13pools_vr(c,j,k) = 0._rkx
              end do
              c13_col_ctrunc_vr(c,j) = 0._rkx
            end do
          end if
          cwdc13(c) = cwdc(c) * c13ratio
          do k = 1, ndecomp_pools
            decomp_c13pools(c,k) = decomp_cpools(c,k) * c13ratio
            decomp_c13pools_1m(c,k) = decomp_cpools_1m(c,k) * c13ratio
          end do
        end if

        if ( use_c14 ) then
          do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
              decomp_c14pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c14ratio
            end do
            c14_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c14ratio
          end do
          if ( nlevdecomp > 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
              do k = 1, ndecomp_pools
                decomp_c14pools_vr(c,j,k) = 0._rkx
              end do
              c14_col_ctrunc_vr(c,j) = 0._rkx
            end do
          end if
          cwdc14(c) = cwdc(c) * c14ratio
          do k = 1, ndecomp_pools
            decomp_c14pools(c,k) = decomp_cpools(c,k) * c14ratio
            decomp_c14pools_1m(c,k) = decomp_cpools_1m(c,k) * c14ratio
          end do
        end if

        ! column nitrogen state variables
        sminn(c) = 0._rkx
        do j = 1, nlevdecomp
          do k = 1, ndecomp_pools
            decomp_npools_vr(c,j,k) = &
                    decomp_cpools_vr(c,j,k) / initial_cn_ratio(k)
          end do
          sminn_vr(c,j) = 0._rkx
          col_ntrunc_vr(c,j) = 0._rkx
        end do
        if ( nlevdecomp > 1 ) then
          do j = nlevdecomp+1, nlevdecomp_full
            do k = 1, ndecomp_pools
              decomp_npools_vr(c,j,k) = 0._rkx
            end do
            sminn_vr(c,j) = 0._rkx
            col_ntrunc_vr(c,j) = 0._rkx
          end do
        end if
        do k = 1, ndecomp_pools
          decomp_npools(c,k) = decomp_cpools(c,k) / initial_cn_ratio(k)
          decomp_npools_1m(c,k) = decomp_cpools_1m(c,k) / initial_cn_ratio(k)
        end do

#ifdef NITRIF_DENITRIF
        do j = 1, nlevdecomp_full
          smin_nh4_vr(c,j) = 0._rkx
          smin_no3_vr(c,j) = 0._rkx
        end do
        smin_nh4(c) = 0._rkx
        smin_no3(c) = 0._rkx
#endif
        totlitn(c)    = 0._rkx
        totsomn(c)    = 0._rkx
        totlitn_1m(c) = 0._rkx
        totsomn_1m(c) = 0._rkx
        totecosysn(c) = 0._rkx
        totcoln(c)    = 0._rkx
        cwdn(c)       = 0._rkx

        ! dynamic landcover state variables
        seedc(c)  = 0._rkx
        prod10c(c)    = 0._rkx
        prod100c(c)   = 0._rkx
        totprodc(c)   = 0._rkx

        if ( use_c13 ) then
          seedc13(c)    = 0._rkx
          prod10c13(c)  = 0._rkx
          prod100c13(c) = 0._rkx
          totprodc13(c) = 0._rkx
        end if

        if ( use_c14 ) then
          seedc14(c)    = 0._rkx
          prod10c14(c)  = 0._rkx
          prod100c14(c) = 0._rkx
          totprodc14(c) = 0._rkx
        end if

        seedn(c)      = 0._rkx
        prod10n(c)    = 0._rkx
        prod100n(c)   = 0._rkx
        totprodn(c)   = 0._rkx

        ! also initialize dynamic landcover fluxes so that they have
        ! real values on first timestep, prior to calling pftdyn_cnbal
        clm3%g%l%c%ccf%dwt_seedc_to_leaf(c) = 0._rkx
        clm3%g%l%c%ccf%dwt_seedc_to_deadstem(c) = 0._rkx
        clm3%g%l%c%ccf%dwt_conv_cflux(c) = 0._rkx
        clm3%g%l%c%ccf%lf_conv_cflux(c) = 0._rkx
        clm3%g%l%c%ccf%dwt_prod10c_gain(c) = 0._rkx
        clm3%g%l%c%ccf%prod10c_loss(c) = 0._rkx
        clm3%g%l%c%ccf%dwt_prod100c_gain(c) = 0._rkx
        clm3%g%l%c%ccf%prod100c_loss(c) = 0._rkx
        do j = 1, nlevdecomp_full
          clm3%g%l%c%ccf%dwt_frootc_to_litr_met_c(c,j) = 0._rkx
          clm3%g%l%c%ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._rkx
          clm3%g%l%c%ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._rkx
          clm3%g%l%c%ccf%dwt_livecrootc_to_cwdc(c,j) = 0._rkx
          clm3%g%l%c%ccf%dwt_deadcrootc_to_cwdc(c,j) = 0._rkx
        end do
        clm3%g%l%c%ccf%dwt_closs(c) = 0._rkx

        if ( use_c13 ) then
          clm3%g%l%c%cc13f%dwt_seedc_to_leaf(c) = 0._rkx
          clm3%g%l%c%cc13f%dwt_seedc_to_deadstem(c) = 0._rkx
          clm3%g%l%c%cc13f%dwt_conv_cflux(c) = 0._rkx
          clm3%g%l%c%cc13f%dwt_prod10c_gain(c) = 0._rkx
          clm3%g%l%c%cc13f%prod10c_loss(c) = 0._rkx
          clm3%g%l%c%cc13f%dwt_prod100c_gain(c) = 0._rkx
          clm3%g%l%c%cc13f%prod100c_loss(c) = 0._rkx
          do j = 1, nlevdecomp_full
            clm3%g%l%c%cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._rkx
            clm3%g%l%c%cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._rkx
            clm3%g%l%c%cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._rkx
            clm3%g%l%c%cc13f%dwt_livecrootc_to_cwdc(c,j) = 0._rkx
            clm3%g%l%c%cc13f%dwt_deadcrootc_to_cwdc(c,j) = 0._rkx
          end do
          clm3%g%l%c%cc13f%dwt_closs(c) = 0._rkx
        end if

        if ( use_c14 ) then
          clm3%g%l%c%cc14f%dwt_seedc_to_leaf(c) = 0._rkx
          clm3%g%l%c%cc14f%dwt_seedc_to_deadstem(c) = 0._rkx
          clm3%g%l%c%cc14f%dwt_conv_cflux(c) = 0._rkx
          clm3%g%l%c%cc14f%dwt_prod10c_gain(c) = 0._rkx
          clm3%g%l%c%cc14f%prod10c_loss(c) = 0._rkx
          clm3%g%l%c%cc14f%dwt_prod100c_gain(c) = 0._rkx
          clm3%g%l%c%cc14f%prod100c_loss(c) = 0._rkx
          do j = 1 , nlevdecomp_full
            clm3%g%l%c%cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._rkx
            clm3%g%l%c%cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._rkx
            clm3%g%l%c%cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._rkx
            clm3%g%l%c%cc14f%dwt_livecrootc_to_cwdc(c,j) = 0._rkx
            clm3%g%l%c%cc14f%dwt_deadcrootc_to_cwdc(c,j) = 0._rkx
          end do
          clm3%g%l%c%cc14f%dwt_closs(c) = 0._rkx
        end if

        clm3%g%l%c%cnf%dwt_seedn_to_leaf(c) = 0._rkx
        clm3%g%l%c%cnf%dwt_seedn_to_deadstem(c) = 0._rkx
        clm3%g%l%c%cnf%dwt_conv_nflux(c) = 0._rkx
        clm3%g%l%c%cnf%dwt_prod10n_gain(c) = 0._rkx
        clm3%g%l%c%cnf%prod10n_loss(c) = 0._rkx
        clm3%g%l%c%cnf%dwt_prod100n_gain(c) = 0._rkx
        clm3%g%l%c%cnf%prod100n_loss(c) = 0._rkx
        do j = 1 , nlevdecomp_full
          clm3%g%l%c%cnf%dwt_frootn_to_litr_met_n(c,j) = 0._rkx
          clm3%g%l%c%cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._rkx
          clm3%g%l%c%cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._rkx
          clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn(c,j) = 0._rkx
          clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn(c,j) = 0._rkx
        end do
        clm3%g%l%c%cnf%dwt_nloss(c) = 0._rkx
      end if
    end do

    ! initialize pft-level variables
    do p = begp , endp
      l = plandunit(p)
      if ( itypelun(l) == istsoil .or. itypelun(l) == istcrop ) then
        ! carbon state variables
        if (ivt(p) == noveg) then
          leafc(p) = 0._rkx
          leafc_storage(p) = 0._rkx
        else
          if (evergreen(ivt(p)) == 1._rkx) then
            leafc(p) = 1._rkx
            leafc_storage(p) = 0._rkx
          else if (ivt(p) >= npcropmin) then ! prognostic crop types
            leafc(p) = 0._rkx
            leafc_storage(p) = 0._rkx
          else
            leafc(p) = 0._rkx
            leafc_storage(p) = 1._rkx
          end if
        end if
        leafc_xfer(p) = 0._rkx
        if ( crop_prog )then
          grainc(p) = 0._rkx
          grainc_storage(p) = 0._rkx
          grainc_xfer(p) = 0._rkx
          fert(p) = 0._rkx
          soyfixn(p) = 0._rkx
        end if
        frootc(p) = 0._rkx
        frootc_storage(p) = 0._rkx
        frootc_xfer(p) = 0._rkx
        livestemc(p) = 0._rkx
        livestemc_storage(p) = 0._rkx
        livestemc_xfer(p) = 0._rkx

        ! tree types need to be initialized with some stem mass so that
        ! roughness length is not zero in canopy flux calculation

        if (woody(ivt(p)) == 1._rkx) then
          deadstemc(p) = 0.1_rkx
        else
          deadstemc(p) = 0._rkx
        end if

        deadstemc_storage(p) = 0._rkx
        deadstemc_xfer(p) = 0._rkx
        livecrootc(p) = 0._rkx
        livecrootc_storage(p) = 0._rkx
        livecrootc_xfer(p) = 0._rkx
        deadcrootc(p) = 0._rkx
        deadcrootc_storage(p) = 0._rkx
        deadcrootc_xfer(p) = 0._rkx
        gresp_storage(p) = 0._rkx
        gresp_xfer(p) = 0._rkx
        cpool(p) = 0._rkx
        xsmrpool(p) = 0._rkx
        pft_ctrunc(p) = 0._rkx
        dispvegc(p) = 0._rkx
        storvegc(p) = 0._rkx
        totpftc(p)  = 0._rkx
        ! calculate totvegc explicitly so that it is available for the isotope
        ! code on the first time step.
        totvegc(p)  = leafc(p) + leafc_storage(p) + &
                leafc_xfer(p) + frootc(p) +  &
                frootc_storage(p) + frootc_xfer(p) + &
                livestemc(p) + livestemc_storage(p) +  &
                livestemc_xfer(p) + deadstemc(p) + &
                deadstemc_storage(p) + deadstemc_xfer(p) +  &
                livecrootc(p) + livecrootc_storage(p) + &
                livecrootc_xfer(p) + deadcrootc(p) +  &
                deadcrootc_storage(p) + deadcrootc_xfer(p) + &
                gresp_storage(p) +  gresp_xfer(p) + cpool(p)

        if ( crop_prog )then
          totvegc(p) = totvegc(p) + grainc(p) + &
                  grainc_storage(p) + grainc_xfer(p)
        end if

        woodc(p)    = 0._rkx

        if ( use_c13 ) then
          leafc13(p)               = leafc(p)               * c13ratio
          leafc13_storage(p)       = leafc_storage(p)       * c13ratio
          leafc13_xfer(p)          = leafc_xfer(p)          * c13ratio
          frootc13(p)              = frootc(p)              * c13ratio
          frootc13_storage(p)      = frootc_storage(p)      * c13ratio
          frootc13_xfer(p)         = frootc_xfer(p)         * c13ratio
          livestemc13(p)           = livestemc(p)           * c13ratio
          livestemc13_storage(p)   = livestemc_storage(p)   * c13ratio
          livestemc13_xfer(p)      = livestemc_xfer(p)      * c13ratio
          deadstemc13(p)           = deadstemc(p)           * c13ratio
          deadstemc13_storage(p)   = deadstemc_storage(p)   * c13ratio
          deadstemc13_xfer(p)      = deadstemc_xfer(p)      * c13ratio
          livecrootc13(p)          = livecrootc(p)          * c13ratio
          livecrootc13_storage(p)  = livecrootc_storage(p)  * c13ratio
          livecrootc13_xfer(p)     = livecrootc_xfer(p)     * c13ratio
          deadcrootc13(p)          = deadcrootc(p)          * c13ratio
          deadcrootc13_storage(p)  = deadcrootc_storage(p)  * c13ratio
          deadcrootc13_xfer(p)     = deadcrootc_xfer(p)     * c13ratio
          c13_gresp_storage(p)     = gresp_storage(p)       * c13ratio
          c13_gresp_xfer(p)        = gresp_xfer(p)          * c13ratio
          c13pool(p)               = cpool(p)               * c13ratio
          c13xsmrpool(p)           = xsmrpool(p)            * c13ratio
          c13_pft_ctrunc(p)        = pft_ctrunc(p)          * c13ratio

          ! calculate totvegc explicitly so that it is available for the isotope
          ! code on the first time step.
          totvegc13(p)  = leafc13(p) + leafc13_storage(p) + &
                  leafc13_xfer(p) + frootc13(p) +  &
                  frootc13_storage(p) + frootc13_xfer(p) + &
                  livestemc13(p) + livestemc13_storage(p) +  &
                  livestemc13_xfer(p) + deadstemc13(p) + &
                  deadstemc13_storage(p) + deadstemc13_xfer(p) +  &
                  livecrootc13(p) + livecrootc13_storage(p) + &
                  livecrootc13_xfer(p) + deadcrootc13(p) +  &
                  deadcrootc13_storage(p) + deadcrootc13_xfer(p) + &
                  c13_gresp_storage(p) +  &
                  c13_gresp_xfer(p) + c13pool(p)
        end if

        if ( use_c14 ) then
          leafc14(p)               = leafc(p)               * c14ratio
          leafc14_storage(p)       = leafc_storage(p)       * c14ratio
          leafc14_xfer(p)          = leafc_xfer(p)          * c14ratio
          frootc14(p)              = frootc(p)              * c14ratio
          frootc14_storage(p)      = frootc_storage(p)      * c14ratio
          frootc14_xfer(p)         = frootc_xfer(p)         * c14ratio
          livestemc14(p)           = livestemc(p)           * c14ratio
          livestemc14_storage(p)   = livestemc_storage(p)   * c14ratio
          livestemc14_xfer(p)      = livestemc_xfer(p)      * c14ratio
          deadstemc14(p)           = deadstemc(p)           * c14ratio
          deadstemc14_storage(p)   = deadstemc_storage(p)   * c14ratio
          deadstemc14_xfer(p)      = deadstemc_xfer(p)      * c14ratio
          livecrootc14(p)          = livecrootc(p)          * c14ratio
          livecrootc14_storage(p)  = livecrootc_storage(p)  * c14ratio
          livecrootc14_xfer(p)     = livecrootc_xfer(p)     * c14ratio
          deadcrootc14(p)          = deadcrootc(p)          * c14ratio
          deadcrootc14_storage(p)  = deadcrootc_storage(p)  * c14ratio
          deadcrootc14_xfer(p)     = deadcrootc_xfer(p)     * c14ratio
          c14_gresp_storage(p)     = gresp_storage(p)       * c14ratio
          c14_gresp_xfer(p)        = gresp_xfer(p)          * c14ratio
          c14pool(p)               = cpool(p)               * c14ratio
          c14xsmrpool(p)           = xsmrpool(p)            * c14ratio
          c14_pft_ctrunc(p)        = pft_ctrunc(p)          * c14ratio

          ! calculate totvegc explicitly so that it is available for the isotope
          ! code on the first time step.
          totvegc14(p)  = leafc14(p) + leafc14_storage(p) + &
                  leafc14_xfer(p) + frootc14(p) +  &
                  frootc14_storage(p) + frootc14_xfer(p) + &
                  livestemc14(p) + livestemc14_storage(p) +  &
                  livestemc14_xfer(p) + deadstemc14(p) + &
                  deadstemc14_storage(p) + deadstemc14_xfer(p) +  &
                  livecrootc14(p) + livecrootc14_storage(p) + &
                  livecrootc14_xfer(p) + deadcrootc14(p) +  &
                  deadcrootc14_storage(p) + deadcrootc14_xfer(p) + &
                  c14_gresp_storage(p) +  &
                  c14_gresp_xfer(p) + c14pool(p)
          rc14_atm(p) = c14ratio
        end if

        ! nitrogen state variables
        if (ivt(p) == noveg) then
          leafn(p) = 0._rkx
          leafn_storage(p) = 0._rkx
        else
          leafn(p) = leafc(p) / leafcn(ivt(p))
          leafn_storage(p) = leafc_storage(p) / leafcn(ivt(p))
        end if

        leafn_xfer(p) = 0._rkx
        if ( crop_prog )then
          grainn(p) = 0._rkx
          grainn_storage(p) = 0._rkx
          grainn_xfer(p) = 0._rkx
        end if
        frootn(p) = 0._rkx
        frootn_storage(p) = 0._rkx
        frootn_xfer(p) = 0._rkx
        livestemn(p) = 0._rkx
        livestemn_storage(p) = 0._rkx
        livestemn_xfer(p) = 0._rkx

        ! tree types need to be initialized with some stem mass so that
        ! roughness length is not zero in canopy flux calculation

        if (woody(ivt(p)) == 1._rkx) then
          deadstemn(p) = deadstemc(p) / deadwdcn(ivt(p))
        else
          deadstemn(p) = 0._rkx
        end if

        deadstemn_storage(p) = 0._rkx
        deadstemn_xfer(p) = 0._rkx
        livecrootn(p) = 0._rkx
        livecrootn_storage(p) = 0._rkx
        livecrootn_xfer(p) = 0._rkx
        deadcrootn(p) = 0._rkx
        deadcrootn_storage(p) = 0._rkx
        deadcrootn_xfer(p) = 0._rkx
        retransn(p) = 0._rkx
        npool(p) = 0._rkx
        pft_ntrunc(p) = 0._rkx
        dispvegn(p) = 0._rkx
        storvegn(p) = 0._rkx
        totvegn(p)  = 0._rkx
        totpftn(p)  = 0._rkx

        ! initialization for psnsun and psnsha required for
        ! proper arbitrary initialization of allocation routine
        ! in initial ecosysdyn call

        psnsun(p) = 0._rkx
        psnsha(p) = 0._rkx

        if ( use_c13 ) then
          c13_psnsun(p) = 0._rkx
          c13_psnsha(p) = 0._rkx
        end if

        if ( use_c14 ) then
          c14_psnsun(p) = 0._rkx
          c14_psnsha(p) = 0._rkx
        end if

        laisun(p) = 0._rkx
        laisha(p) = 0._rkx

        ! ecophysiological variables
        ! phenology variables
        dormant_flag(p) = 1._rkx
        days_active(p) = 0._rkx
        onset_flag(p) = 0._rkx
        onset_counter(p) = 0._rkx
        onset_gddflag(p) = 0._rkx
        onset_fdd(p) = 0._rkx
        onset_gdd(p) = 0._rkx
        onset_swi(p) = 0.0_rkx
        offset_flag(p) = 0._rkx
        offset_counter(p) = 0._rkx
        offset_fdd(p) = 0._rkx
        offset_swi(p) = 0._rkx
        lgsf(p) = 0._rkx
        bglfr(p) = 0._rkx
        bgtr(p) = 0._rkx
        annavg_t2m(p) = 280._rkx
        tempavg_t2m(p) = 0._rkx
        fert_counter(p) = 0._rkx
        grain_flag(p) = 0._rkx

        ! non-phenology variables
        gpp(p) = 0._rkx
        availc(p) = 0._rkx
        xsmrpool_recover(p) = 0._rkx
        alloc_pnow(p) = 1._rkx
        c_allometry(p) = 0._rkx
        n_allometry(p) = 0._rkx
        plant_ndemand(p) = 0._rkx
        tempsum_potential_gpp(p) = 0._rkx
        annsum_potential_gpp(p) = 0._rkx
        tempmax_retransn(p) = 0._rkx
        annmax_retransn(p) = 0._rkx
        avail_retransn(p) = 0._rkx
        plant_nalloc(p) = 0._rkx
        plant_calloc(p) = 0._rkx
        excess_cflux(p) = 0._rkx
        downreg(p) = 0._rkx
        prev_leafc_to_litter(p) = 0._rkx
        prev_frootc_to_litter(p) = 0._rkx
        tempsum_npp(p) = 0._rkx
        annsum_npp(p) = 0._rkx
#if (defined CNDV)
        tempsum_litfall(p) = 0._rkx
        annsum_litfall(p) = 0._rkx
#endif

        if ( use_c13 ) then
          xsmrpool_c13ratio(p) = c13ratio
          rc13_canair(p) = 0._rkx
          rc13_psnsun(p) = 0._rkx
          rc13_psnsha(p) = 0._rkx
          alphapsnsun(p) = 0._rkx
          alphapsnsha(p) = 0._rkx
        end if
      end if   ! end of if-istsoil block
    end do   ! end of loop over pfts
#endif
  end subroutine CNiniTimeVar

end module mod_clm_cninitimevar
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
