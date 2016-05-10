module mod_clm_cncstateupdate1
#ifdef CN
  !
  ! Module for carbon state variable update, non-mortality fluxes.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_runparams , only : dtsrf
  use mod_clm_varpar , only : ndecomp_cascade_transitions, nlevdecomp

  implicit none

  save

  private

  public:: CStateUpdate1
  public:: CStateUpdate0

  contains
  !
  ! On the radiation time step, update cpool carbon state
  !
  subroutine CStateUpdate0(num_soilp, filter_soilp, isotope)
    use mod_clm_type
    implicit none
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts
    character(len=*), intent(in) :: isotope     ! 'bulk', 'c13' or 'c14'

    real(rkx), pointer :: psnshade_to_cpool(:)
    real(rkx), pointer :: psnsun_to_cpool(:)

    real(rkx), pointer :: cpool(:)  ! (gC/m2) temporary photosynthate C pool

    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    integer(ik4) :: p     ! indices
    integer(ik4) :: fp   ! lake filter indices
    real(rkx):: dt      ! radiation time step (seconds)

    ! select which isotope
    select case (isotope)
      case ('bulk')
        pcisof => clm3%g%l%c%p%pcf
        pcisos => clm3%g%l%c%p%pcs
      case ('c14')
        pcisof => clm3%g%l%c%p%pc14f
        pcisos => clm3%g%l%c%p%pc14s
      case ('c13')
        pcisof => clm3%g%l%c%p%pc13f
        pcisos => clm3%g%l%c%p%pc13s
      case default
        call fatal(__FILE__,__LINE__, &
          'CNCIsoStateUpdate1Mod: iso must be bulk, c13 or c14')
    end select

    ! assign local pointers at the pft level
    cpool              => pcisos%cpool
    psnshade_to_cpool  => pcisof%psnshade_to_cpool
    psnsun_to_cpool    => pcisof%psnsun_to_cpool

    ! set time steps
    dt = dtsrf

    ! pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)
      ! gross photosynthesis fluxes
      cpool(p) = cpool(p) + psnsun_to_cpool(p)*dt
      cpool(p) = cpool(p) + psnshade_to_cpool(p)*dt
    end do
  end subroutine CStateUpdate0
  !
  ! On the radiation time step, update all the prognostic carbon state
  ! variables (except for gap-phase mortality and fire fluxes)
  !
  subroutine CStateUpdate1(num_soilc, filter_soilc, num_soilp, &
                           filter_soilp, isotope)
    use mod_clm_type
    use mod_clm_varpar , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    use mod_clm_pftvarcon , only : npcropmin
    implicit none
    integer(ik4), intent(in) :: num_soilc       ! number of soil columns filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts
    character(len=*), intent(in) :: isotope     ! 'bulk', 'c13' or 'c14'

    ! binary flag for woody lifeform (1=woody, 0=not woody)
    real(rkx), pointer :: woody(:)
    integer(ik4) , pointer :: ivt(:)       ! pft vegetation type
    integer(ik4) , pointer :: harvdate(:)  ! harvest date
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rkx), pointer :: xsmrpool_to_atm(:)
    real(rkx), pointer :: deadcrootc_xfer_to_deadcrootc(:)
    real(rkx), pointer :: deadstemc_xfer_to_deadstemc(:)
    real(rkx), pointer :: frootc_xfer_to_frootc(:)
    real(rkx), pointer :: leafc_xfer_to_leafc(:)
    real(rkx), pointer :: livecrootc_xfer_to_livecrootc(:)
    real(rkx), pointer :: livestemc_xfer_to_livestemc(:)
    real(rkx), pointer :: cpool_to_xsmrpool(:)
    real(rkx), pointer :: cpool_to_deadcrootc(:)
    real(rkx), pointer :: cpool_to_deadcrootc_storage(:)
    real(rkx), pointer :: cpool_to_deadstemc(:)
    real(rkx), pointer :: cpool_to_deadstemc_storage(:)
    real(rkx), pointer :: cpool_to_frootc(:)
    real(rkx), pointer :: cpool_to_frootc_storage(:)
    real(rkx), pointer :: cpool_to_gresp_storage(:)
    real(rkx), pointer :: cpool_to_leafc(:)
    real(rkx), pointer :: cpool_to_leafc_storage(:)
    real(rkx), pointer :: cpool_to_livecrootc(:)
    real(rkx), pointer :: cpool_to_livecrootc_storage(:)
    real(rkx), pointer :: cpool_to_livestemc(:)
    real(rkx), pointer :: cpool_to_livestemc_storage(:)
    real(rkx), pointer :: deadcrootc_storage_to_xfer(:)
    real(rkx), pointer :: deadstemc_storage_to_xfer(:)
    real(rkx), pointer :: frootc_storage_to_xfer(:)
    real(rkx), pointer :: frootc_to_litter(:)
    real(rkx), pointer :: gresp_storage_to_xfer(:)
    real(rkx), pointer :: leafc_storage_to_xfer(:)
    real(rkx), pointer :: leafc_to_litter(:)
    real(rkx), pointer :: livecrootc_storage_to_xfer(:)
    real(rkx), pointer :: livecrootc_to_deadcrootc(:)
    real(rkx), pointer :: livestemc_storage_to_xfer(:)
    real(rkx), pointer :: livestemc_to_deadstemc(:)
    real(rkx), pointer :: livestem_curmr(:)
    real(rkx), pointer :: froot_curmr(:)
    real(rkx), pointer :: leaf_curmr(:)
    real(rkx), pointer :: livecroot_curmr(:)
    real(rkx), pointer :: grain_curmr(:)
    real(rkx), pointer :: livestem_xsmr(:)
    real(rkx), pointer :: froot_xsmr(:)
    real(rkx), pointer :: leaf_xsmr(:)
    real(rkx), pointer :: livecroot_xsmr(:)
    real(rkx), pointer :: grain_xsmr(:)
    real(rkx), pointer :: cpool_deadcroot_gr(:)
    real(rkx), pointer :: cpool_deadcroot_storage_gr(:)
    real(rkx), pointer :: cpool_deadstem_gr(:)
    real(rkx), pointer :: cpool_deadstem_storage_gr(:)
    real(rkx), pointer :: cpool_froot_gr(:)
    real(rkx), pointer :: cpool_froot_storage_gr(:)
    real(rkx), pointer :: cpool_leaf_gr(:)
    real(rkx), pointer :: cpool_leaf_storage_gr(:)
    real(rkx), pointer :: cpool_livecroot_gr(:)
    real(rkx), pointer :: cpool_livecroot_storage_gr(:)
    ! live stem growth respiration (gC/m2/s)
    real(rkx), pointer :: cpool_livestem_gr(:)
    ! live stem growth respiration to storage (gC/m2/s)
    real(rkx), pointer :: cpool_livestem_storage_gr(:)
    ! dead coarse root growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_deadcroot_gr(:)
    ! dead stem growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_deadstem_gr(:)
    ! fine root  growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_froot_gr(:)
    ! leaf growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_leaf_gr(:)
    ! live coarse root growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_livecroot_gr(:)
    ! live stem growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_livestem_gr(:)
    ! allocation to grain C (gC/m2/s)
    real(rkx), pointer :: cpool_to_grainc(:)
    ! allocation to grain C storage (gC/m2/s)
    real(rkx), pointer :: cpool_to_grainc_storage(:)
    ! grain C shift storage to transfer (gC/m2/s)
    real(rkx), pointer :: grainc_storage_to_xfer(:)
    ! live stem C litterfall (gC/m2/s)
    real(rkx), pointer :: livestemc_to_litter(:)
    ! grain C to food (gC/m2/s)
    real(rkx), pointer :: grainc_to_food(:)
    ! grain C growth from storage (gC/m2/s)
    real(rkx), pointer :: grainc_xfer_to_grainc(:)
    ! grain growth respiration (gC/m2/s)
    real(rkx), pointer :: cpool_grain_gr(:)
    ! grain growth respiration to storage (gC/m2/s)
    real(rkx), pointer :: cpool_grain_storage_gr(:)
    ! grain growth respiration from storage (gC/m2/s)
    real(rkx), pointer :: transfer_grain_gr(:)

    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rkx), pointer :: decomp_cpools_vr(:,:,:)
    ! (gC/m3/timestep)  change in decomposing c pools. Used to update
    ! concentrations concurrently with vertical transport equation
    real(rkx), pointer :: decomp_cpools_sourcesink(:,:,:)
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rkx), pointer :: decomp_cascade_hr_vr(:,:,:)
    ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
    real(rkx), pointer :: decomp_cascade_ctransfer_vr(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer :: cascade_receiver_pool(:)
    real(rkx), pointer :: grainc(:)         ! (gC/m2) grain C
    real(rkx), pointer :: grainc_storage(:) ! (gC/m2) grain C storage
    real(rkx), pointer :: grainc_xfer(:)    ! (gC/m2) grain C transfer
    real(rkx), pointer :: cpool(:)      ! (gC/m2) temporary photosynthate C pool
    real(rkx), pointer :: xsmrpool(:)   ! (gC/m2) execss maint resp C pool
    real(rkx), pointer :: deadcrootc(:) ! (gC/m2) dead coarse root C
    ! (gC/m2) dead coarse root C storage
    real(rkx), pointer :: deadcrootc_storage(:)
    ! (gC/m2) dead coarse root C transfer
    real(rkx), pointer :: deadcrootc_xfer(:)
    real(rkx), pointer :: deadstemc(:)         ! (gC/m2) dead stem C
    real(rkx), pointer :: deadstemc_storage(:) ! (gC/m2) dead stem C storage
    real(rkx), pointer :: deadstemc_xfer(:)    ! (gC/m2) dead stem C transfer
    real(rkx), pointer :: frootc(:)            ! (gC/m2) fine root C
    real(rkx), pointer :: frootc_storage(:)    ! (gC/m2) fine root C storage
    real(rkx), pointer :: frootc_xfer(:)       ! (gC/m2) fine root C transfer
    ! (gC/m2) growth respiration storage
    real(rkx), pointer :: gresp_storage(:)
    ! (gC/m2) growth respiration transfer
    real(rkx), pointer :: gresp_xfer(:)
    real(rkx), pointer :: leafc(:)           ! (gC/m2) leaf C
    real(rkx), pointer :: leafc_storage(:)   ! (gC/m2) leaf C storage
    real(rkx), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
    real(rkx), pointer :: livecrootc(:)      ! (gC/m2) live coarse root C
    ! (gC/m2) live coarse root C storage
    real(rkx), pointer :: livecrootc_storage(:)
    ! (gC/m2) live coarse root C transfer
    real(rkx), pointer :: livecrootc_xfer(:)
    real(rkx), pointer :: livestemc(:)          ! (gC/m2) live stem C
    real(rkx), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
    real(rkx), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
    ! C fluxes associated with phenology (litterfall and crop)
    ! to litter metabolic pool (gC/m3/s)
    real(rkx), pointer :: phenology_c_to_litr_met_c(:,:)
    ! C fluxes associated with phenology (litterfall and crop)
    ! to litter cellulose pool (gC/m3/s)
    real(rkx), pointer :: phenology_c_to_litr_cel_c(:,:)
    ! C fluxes associated with phenology (litterfall and crop)
    ! to litter lignin pool (gC/m3/s)
    real(rkx), pointer :: phenology_c_to_litr_lig_c(:,:)

    real(rkx), pointer :: dwt_seedc_to_leaf(:)
    real(rkx), pointer :: dwt_seedc_to_deadstem(:)
    real(rkx), pointer :: dwt_frootc_to_litr_met_c(:,:)
    real(rkx), pointer :: dwt_frootc_to_litr_cel_c(:,:)
    real(rkx), pointer :: dwt_frootc_to_litr_lig_c(:,:)
    ! (gC/m2/s) live coarse root to CWD due to landcover change
    real(rkx), pointer :: dwt_livecrootc_to_cwdc(:,:)
    ! (gC/m2/s) dead coarse root to CWD due to landcover change
    real(rkx), pointer :: dwt_deadcrootc_to_cwdc(:,:)
    real(rkx), pointer :: seedc(:)

    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    type(column_cflux_type), pointer :: ccisof
    type(column_cstate_type), pointer :: ccisos
    integer(ik4) :: c,p,j,k     ! indices
    integer(ik4) :: fp,fc   ! lake filter indices
    real(rkx):: dt      ! radiation time step (seconds)

    ! select which isotope
    select case (isotope)
      case ('bulk')
        pcisof => clm3%g%l%c%p%pcf
        pcisos => clm3%g%l%c%p%pcs
        ccisof => clm3%g%l%c%ccf
        ccisos => clm3%g%l%c%ccs
      case ('c14')
        pcisof => clm3%g%l%c%p%pc14f
        pcisos => clm3%g%l%c%p%pc14s
        ccisof => clm3%g%l%c%cc14f
        ccisos => clm3%g%l%c%cc14s
      case ('c13')
        pcisof => clm3%g%l%c%p%pc13f
        pcisos => clm3%g%l%c%p%pc13s
        ccisof => clm3%g%l%c%cc13f
        ccisos => clm3%g%l%c%cc13s
      case default
        call fatal(__FILE__,__LINE__,&
           'CNCIsoStateUpdate1Mod: iso must be bulk, c13 or c14')
    end select

    ! assign local pointers
    woody => pftcon%woody

    ! assign local pointers at the column level
    decomp_cpools_vr            => ccisos%decomp_cpools_vr
    decomp_cpools_sourcesink    => ccisof%decomp_cpools_sourcesink
    decomp_cascade_hr_vr        => ccisof%decomp_cascade_hr_vr
    decomp_cascade_ctransfer_vr => ccisof%decomp_cascade_ctransfer_vr
    cascade_donor_pool          => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool       => decomp_cascade_con%cascade_receiver_pool
    phenology_c_to_litr_met_c   => ccisof%phenology_c_to_litr_met_c
    phenology_c_to_litr_cel_c   => ccisof%phenology_c_to_litr_cel_c
    phenology_c_to_litr_lig_c   => ccisof%phenology_c_to_litr_lig_c

    ! new pointers for dynamic landcover
    dwt_seedc_to_leaf              => ccisof%dwt_seedc_to_leaf
    dwt_seedc_to_deadstem          => ccisof%dwt_seedc_to_deadstem
    dwt_frootc_to_litr_met_c       => ccisof%dwt_frootc_to_litr_met_c
    dwt_frootc_to_litr_cel_c       => ccisof%dwt_frootc_to_litr_cel_c
    dwt_frootc_to_litr_lig_c       => ccisof%dwt_frootc_to_litr_lig_c
    dwt_livecrootc_to_cwdc         => ccisof%dwt_livecrootc_to_cwdc
    dwt_deadcrootc_to_cwdc         => ccisof%dwt_deadcrootc_to_cwdc
    seedc                          => ccisos%seedc

    ! assign local pointers at the pft level
    ivt                            => clm3%g%l%c%p%itype
    cpool_deadcroot_gr             => pcisof%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => pcisof%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => pcisof%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => pcisof%cpool_deadstem_storage_gr
    cpool_froot_gr                 => pcisof%cpool_froot_gr
    cpool_froot_storage_gr         => pcisof%cpool_froot_storage_gr
    cpool_leaf_gr                  => pcisof%cpool_leaf_gr
    cpool_leaf_storage_gr          => pcisof%cpool_leaf_storage_gr
    cpool_livecroot_gr             => pcisof%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => pcisof%cpool_livecroot_storage_gr
    cpool_livestem_gr              => pcisof%cpool_livestem_gr
    cpool_livestem_storage_gr      => pcisof%cpool_livestem_storage_gr
    cpool_to_xsmrpool              => pcisof%cpool_to_xsmrpool
    cpool_to_deadcrootc            => pcisof%cpool_to_deadcrootc
    cpool_to_deadcrootc_storage    => pcisof%cpool_to_deadcrootc_storage
    cpool_to_deadstemc             => pcisof%cpool_to_deadstemc
    cpool_to_deadstemc_storage     => pcisof%cpool_to_deadstemc_storage
    cpool_to_frootc                => pcisof%cpool_to_frootc
    cpool_to_frootc_storage        => pcisof%cpool_to_frootc_storage
    cpool_to_gresp_storage         => pcisof%cpool_to_gresp_storage
    cpool_to_leafc                 => pcisof%cpool_to_leafc
    cpool_to_leafc_storage         => pcisof%cpool_to_leafc_storage
    cpool_to_livecrootc            => pcisof%cpool_to_livecrootc
    cpool_to_livecrootc_storage    => pcisof%cpool_to_livecrootc_storage
    cpool_to_livestemc             => pcisof%cpool_to_livestemc
    cpool_to_livestemc_storage     => pcisof%cpool_to_livestemc_storage
    deadcrootc_storage_to_xfer     => pcisof%deadcrootc_storage_to_xfer
    deadcrootc_xfer_to_deadcrootc  => pcisof%deadcrootc_xfer_to_deadcrootc
    deadstemc_storage_to_xfer      => pcisof%deadstemc_storage_to_xfer
    deadstemc_xfer_to_deadstemc    => pcisof%deadstemc_xfer_to_deadstemc
    froot_curmr                    => pcisof%froot_curmr
    froot_xsmr                     => pcisof%froot_xsmr
    frootc_storage_to_xfer         => pcisof%frootc_storage_to_xfer
    frootc_to_litter               => pcisof%frootc_to_litter
    frootc_xfer_to_frootc          => pcisof%frootc_xfer_to_frootc
    gresp_storage_to_xfer          => pcisof%gresp_storage_to_xfer
    leaf_curmr                     => pcisof%leaf_curmr
    leaf_xsmr                      => pcisof%leaf_xsmr
    leafc_storage_to_xfer          => pcisof%leafc_storage_to_xfer
    leafc_to_litter                => pcisof%leafc_to_litter
    leafc_xfer_to_leafc            => pcisof%leafc_xfer_to_leafc
    livecroot_curmr                => pcisof%livecroot_curmr
    livecroot_xsmr                 => pcisof%livecroot_xsmr
    livecrootc_storage_to_xfer     => pcisof%livecrootc_storage_to_xfer
    livecrootc_to_deadcrootc       => pcisof%livecrootc_to_deadcrootc
    livecrootc_xfer_to_livecrootc  => pcisof%livecrootc_xfer_to_livecrootc
    livestem_curmr                 => pcisof%livestem_curmr
    livestem_xsmr                  => pcisof%livestem_xsmr
    livestemc_storage_to_xfer      => pcisof%livestemc_storage_to_xfer
    livestemc_to_deadstemc         => pcisof%livestemc_to_deadstemc
    livestemc_xfer_to_livestemc    => pcisof%livestemc_xfer_to_livestemc
    transfer_deadcroot_gr          => pcisof%transfer_deadcroot_gr
    transfer_deadstem_gr           => pcisof%transfer_deadstem_gr
    transfer_froot_gr              => pcisof%transfer_froot_gr
    transfer_leaf_gr               => pcisof%transfer_leaf_gr
    transfer_livecroot_gr          => pcisof%transfer_livecroot_gr
    transfer_livestem_gr           => pcisof%transfer_livestem_gr
    harvdate                       => clm3%g%l%c%p%pps%harvdate
    xsmrpool_to_atm                => pcisof%xsmrpool_to_atm
    cpool_grain_gr                 => pcisof%cpool_grain_gr
    cpool_grain_storage_gr         => pcisof%cpool_grain_storage_gr
    cpool_to_grainc                => pcisof%cpool_to_grainc
    cpool_to_grainc_storage        => pcisof%cpool_to_grainc_storage
    livestemc_to_litter            => pcisof%livestemc_to_litter
    grain_curmr                    => clm3%g%l%c%p%pcf%grain_curmr
    grain_xsmr                     => clm3%g%l%c%p%pcf%grain_xsmr
    grainc_storage_to_xfer         => pcisof%grainc_storage_to_xfer
    grainc_to_food                 => pcisof%grainc_to_food
    grainc_xfer_to_grainc          => pcisof%grainc_xfer_to_grainc
    transfer_grain_gr              => pcisof%transfer_grain_gr
    grainc                         => pcisos%grainc
    grainc_storage                 => pcisos%grainc_storage
    grainc_xfer                    => pcisos%grainc_xfer
    cpool                          => pcisos%cpool
    xsmrpool                       => pcisos%xsmrpool
    deadcrootc                     => pcisos%deadcrootc
    deadcrootc_storage             => pcisos%deadcrootc_storage
    deadcrootc_xfer                => pcisos%deadcrootc_xfer
    deadstemc                      => pcisos%deadstemc
    deadstemc_storage              => pcisos%deadstemc_storage
    deadstemc_xfer                 => pcisos%deadstemc_xfer
    frootc                         => pcisos%frootc
    frootc_storage                 => pcisos%frootc_storage
    frootc_xfer                    => pcisos%frootc_xfer
    gresp_storage                  => pcisos%gresp_storage
    gresp_xfer                     => pcisos%gresp_xfer
    leafc                          => pcisos%leafc
    leafc_storage                  => pcisos%leafc_storage
    leafc_xfer                     => pcisos%leafc_xfer
    livecrootc                     => pcisos%livecrootc
    livecrootc_storage             => pcisos%livecrootc_storage
    livecrootc_xfer                => pcisos%livecrootc_xfer
    livestemc                      => pcisos%livestemc
    livestemc_storage              => pcisos%livestemc_storage
    livestemc_xfer                 => pcisos%livestemc_xfer

    ! set time steps
    dt = dtsrf

    ! column level fluxes

    ! column loop
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      ! seeding fluxes, from dynamic landcover
      seedc(c) = seedc(c) - dwt_seedc_to_leaf(c) * dt
      seedc(c) = seedc(c) - dwt_seedc_to_deadstem(c) * dt
    end do

    ! plant to litter fluxes
    do j = 1 , nlevdecomp
     ! column loop
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       ! phenology and dynamic land cover fluxes
       decomp_cpools_sourcesink(c,j,i_met_lit) = &
         ( phenology_c_to_litr_met_c(c,j) + dwt_frootc_to_litr_met_c(c,j) ) *dt
       decomp_cpools_sourcesink(c,j,i_cel_lit) = &
         ( phenology_c_to_litr_cel_c(c,j) + dwt_frootc_to_litr_cel_c(c,j) ) *dt
       decomp_cpools_sourcesink(c,j,i_lig_lit) = &
         ( phenology_c_to_litr_lig_c(c,j) + dwt_frootc_to_litr_lig_c(c,j) ) *dt
       decomp_cpools_sourcesink(c,j,i_cwd) = &
         ( dwt_livecrootc_to_cwdc(c,j) + dwt_deadcrootc_to_cwdc(c,j) ) *dt
      end do
    end do

    ! litter and SOM HR fluxes
    do k = 1 , ndecomp_cascade_transitions
      do j = 1 , nlevdecomp
        ! column loop
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          decomp_cpools_sourcesink(c,j,cascade_donor_pool(k)) = &
            decomp_cpools_sourcesink(c,j,cascade_donor_pool(k)) - &
            ( decomp_cascade_hr_vr(c,j,k) + &
              decomp_cascade_ctransfer_vr(c,j,k)) *dt
        end do
      end do
    end do
    do k = 1 , ndecomp_cascade_transitions
      if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
        do j = 1 , nlevdecomp
          ! column loop
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            decomp_cpools_sourcesink(c,j,cascade_receiver_pool(k)) = &
              decomp_cpools_sourcesink(c,j,cascade_receiver_pool(k)) + &
              decomp_cascade_ctransfer_vr(c,j,k)*dt
          end do
        end do
      end if
    end do

    ! pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      ! phenology: transfer growth fluxes
      leafc(p) = leafc(p) + leafc_xfer_to_leafc(p)*dt
      leafc_xfer(p) = leafc_xfer(p) - leafc_xfer_to_leafc(p)*dt
      frootc(p) = frootc(p) + frootc_xfer_to_frootc(p)*dt
      frootc_xfer(p) = frootc_xfer(p) - frootc_xfer_to_frootc(p)*dt
      if (woody(ivt(p)) == 1._rkx) then
        livestemc(p) = livestemc(p) + livestemc_xfer_to_livestemc(p)*dt
        livestemc_xfer(p) = livestemc_xfer(p) - &
          livestemc_xfer_to_livestemc(p)*dt
        deadstemc(p) = deadstemc(p) + deadstemc_xfer_to_deadstemc(p)*dt
        deadstemc_xfer(p) = deadstemc_xfer(p) - &
          deadstemc_xfer_to_deadstemc(p)*dt
        livecrootc(p) = livecrootc(p) + livecrootc_xfer_to_livecrootc(p)*dt
        livecrootc_xfer(p) = livecrootc_xfer(p) - &
          livecrootc_xfer_to_livecrootc(p)*dt
        deadcrootc(p) = deadcrootc(p) + deadcrootc_xfer_to_deadcrootc(p)*dt
        deadcrootc_xfer(p) = deadcrootc_xfer(p) - &
          deadcrootc_xfer_to_deadcrootc(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        ! lines here for consistency; the transfer terms are zero
        livestemc(p) = livestemc(p) + livestemc_xfer_to_livestemc(p)*dt
        livestemc_xfer(p) = livestemc_xfer(p) - &
          livestemc_xfer_to_livestemc(p)*dt
        grainc(p) = grainc(p) + grainc_xfer_to_grainc(p)*dt
        grainc_xfer(p) = grainc_xfer(p) - grainc_xfer_to_grainc(p)*dt
      end if

      ! phenology: litterfall fluxes
      leafc(p) = leafc(p) - leafc_to_litter(p)*dt
      frootc(p) = frootc(p) - frootc_to_litter(p)*dt

      ! livewood turnover fluxes
      if ( woody(ivt(p)) == 1._rkx ) then
        livestemc(p)  = livestemc(p)  - livestemc_to_deadstemc(p)*dt
        deadstemc(p)  = deadstemc(p)  + livestemc_to_deadstemc(p)*dt
        livecrootc(p) = livecrootc(p) - livecrootc_to_deadcrootc(p)*dt
        deadcrootc(p) = deadcrootc(p) + livecrootc_to_deadcrootc(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        livestemc(p)  = livestemc(p)  - livestemc_to_litter(p)*dt
        grainc(p)     = grainc(p)     - grainc_to_food(p)*dt
      end if

      ! maintenance respiration fluxes from cpool
      cpool(p) = cpool(p) - cpool_to_xsmrpool(p)*dt
      cpool(p) = cpool(p) - leaf_curmr(p)*dt
      cpool(p) = cpool(p) - froot_curmr(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        cpool(p) = cpool(p) - livestem_curmr(p)*dt
        cpool(p) = cpool(p) - livecroot_curmr(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        cpool(p) = cpool(p) - livestem_curmr(p)*dt
        cpool(p) = cpool(p) - grain_curmr(p)*dt
      end if

      ! maintenance respiration fluxes from xsmrpool
      xsmrpool(p) = xsmrpool(p) + cpool_to_xsmrpool(p)*dt
      xsmrpool(p) = xsmrpool(p) - leaf_xsmr(p)*dt
      xsmrpool(p) = xsmrpool(p) - froot_xsmr(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        xsmrpool(p) = xsmrpool(p) - livestem_xsmr(p)*dt
        xsmrpool(p) = xsmrpool(p) - livecroot_xsmr(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        xsmrpool(p) = xsmrpool(p) - livestem_xsmr(p)*dt
        xsmrpool(p) = xsmrpool(p) - grain_xsmr(p)*dt
        if ( harvdate(p) < 999 ) then ! beginning at harvest, send to atm
          xsmrpool_to_atm(p) = xsmrpool_to_atm(p) + xsmrpool(p)/dt
          xsmrpool(p) = xsmrpool(p) - xsmrpool_to_atm(p)*dt
        end if
      end if

      ! allocation fluxes
      cpool(p) = cpool(p) - cpool_to_leafc(p)*dt
      leafc(p) = leafc(p) + cpool_to_leafc(p)*dt
      cpool(p) = cpool(p) - cpool_to_leafc_storage(p)*dt
      leafc_storage(p) = leafc_storage(p) + cpool_to_leafc_storage(p)*dt
      cpool(p) = cpool(p) - cpool_to_frootc(p)*dt
      frootc(p) = frootc(p) + cpool_to_frootc(p)*dt
      cpool(p) = cpool(p) - cpool_to_frootc_storage(p)*dt
      frootc_storage(p) = frootc_storage(p) + cpool_to_frootc_storage(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        cpool(p) = cpool(p) - cpool_to_livestemc(p)*dt
        livestemc(p) = livestemc(p) + cpool_to_livestemc(p)*dt
        cpool(p) = cpool(p) - cpool_to_livestemc_storage(p)*dt
        livestemc_storage(p) = livestemc_storage(p) + &
          cpool_to_livestemc_storage(p)*dt
        cpool(p) = cpool(p) - cpool_to_deadstemc(p)*dt
        deadstemc(p) = deadstemc(p) + cpool_to_deadstemc(p)*dt
        cpool(p) = cpool(p) - cpool_to_deadstemc_storage(p)*dt
        deadstemc_storage(p) = deadstemc_storage(p) + &
          cpool_to_deadstemc_storage(p)*dt
        cpool(p) = cpool(p) - cpool_to_livecrootc(p)*dt
        livecrootc(p) = livecrootc(p) + cpool_to_livecrootc(p)*dt
        cpool(p) = cpool(p) - cpool_to_livecrootc_storage(p)*dt
        livecrootc_storage(p) = livecrootc_storage(p) + &
          cpool_to_livecrootc_storage(p)*dt
        cpool(p) = cpool(p) - cpool_to_deadcrootc(p)*dt
        deadcrootc(p) = deadcrootc(p) + cpool_to_deadcrootc(p)*dt
        cpool(p) = cpool(p) - cpool_to_deadcrootc_storage(p)*dt
        deadcrootc_storage(p) = deadcrootc_storage(p) + &
          cpool_to_deadcrootc_storage(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        cpool(p) = cpool(p) - cpool_to_livestemc(p)*dt
        livestemc(p) = livestemc(p) + cpool_to_livestemc(p)*dt
        cpool(p) = cpool(p) - cpool_to_livestemc_storage(p)*dt
        livestemc_storage(p) = livestemc_storage(p) + &
          cpool_to_livestemc_storage(p)*dt
        cpool(p) = cpool(p) - cpool_to_grainc(p)*dt
        grainc(p) = grainc(p) + cpool_to_grainc(p)*dt
        cpool(p) = cpool(p) - cpool_to_grainc_storage(p)*dt
        grainc_storage(p) = grainc_storage(p) + cpool_to_grainc_storage(p)*dt
      end if

      ! growth respiration fluxes for current growth
      cpool(p) = cpool(p) - cpool_leaf_gr(p)*dt
      cpool(p) = cpool(p) - cpool_froot_gr(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        cpool(p) = cpool(p) - cpool_livestem_gr(p)*dt
        cpool(p) = cpool(p) - cpool_deadstem_gr(p)*dt
        cpool(p) = cpool(p) - cpool_livecroot_gr(p)*dt
        cpool(p) = cpool(p) - cpool_deadcroot_gr(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        cpool(p) = cpool(p) - cpool_livestem_gr(p)*dt
        cpool(p) = cpool(p) - cpool_grain_gr(p)*dt
      end if

      ! growth respiration for transfer growth
      gresp_xfer(p) = gresp_xfer(p) - transfer_leaf_gr(p)*dt
      gresp_xfer(p) = gresp_xfer(p) - transfer_froot_gr(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        gresp_xfer(p) = gresp_xfer(p) - transfer_livestem_gr(p)*dt
        gresp_xfer(p) = gresp_xfer(p) - transfer_deadstem_gr(p)*dt
        gresp_xfer(p) = gresp_xfer(p) - transfer_livecroot_gr(p)*dt
        gresp_xfer(p) = gresp_xfer(p) - transfer_deadcroot_gr(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        gresp_xfer(p) = gresp_xfer(p) - transfer_livestem_gr(p)*dt
        gresp_xfer(p) = gresp_xfer(p) - transfer_grain_gr(p)*dt
      end if

      ! growth respiration at time of storage
      cpool(p) = cpool(p) - cpool_leaf_storage_gr(p)*dt
      cpool(p) = cpool(p) - cpool_froot_storage_gr(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        cpool(p) = cpool(p) - cpool_livestem_storage_gr(p)*dt
        cpool(p) = cpool(p) - cpool_deadstem_storage_gr(p)*dt
        cpool(p) = cpool(p) - cpool_livecroot_storage_gr(p)*dt
        cpool(p) = cpool(p) - cpool_deadcroot_storage_gr(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        cpool(p) = cpool(p) - cpool_livestem_storage_gr(p)*dt
        cpool(p) = cpool(p) - cpool_grain_storage_gr(p)*dt
      end if

      ! growth respiration stored for release during transfer growth
      cpool(p)         = cpool(p)         - cpool_to_gresp_storage(p)*dt
      gresp_storage(p) = gresp_storage(p) + cpool_to_gresp_storage(p)*dt

      ! move storage pools into transfer pools
      leafc_storage(p) = leafc_storage(p) - leafc_storage_to_xfer(p)*dt
      leafc_xfer(p) = leafc_xfer(p) + leafc_storage_to_xfer(p)*dt
      frootc_storage(p) = frootc_storage(p) - frootc_storage_to_xfer(p)*dt
      frootc_xfer(p) = frootc_xfer(p) + frootc_storage_to_xfer(p)*dt
      if ( woody(ivt(p)) == 1._rkx ) then
        livestemc_storage(p) = livestemc_storage(p) - &
          livestemc_storage_to_xfer(p)*dt
        livestemc_xfer(p) = livestemc_xfer(p) + &
          livestemc_storage_to_xfer(p)*dt
        deadstemc_storage(p) = deadstemc_storage(p) - &
          deadstemc_storage_to_xfer(p)*dt
        deadstemc_xfer(p) = deadstemc_xfer(p) + &
          deadstemc_storage_to_xfer(p)*dt
        livecrootc_storage(p) = livecrootc_storage(p) - &
          livecrootc_storage_to_xfer(p)*dt
        livecrootc_xfer(p) = livecrootc_xfer(p) + &
          livecrootc_storage_to_xfer(p)*dt
        deadcrootc_storage(p) = deadcrootc_storage(p) - &
          deadcrootc_storage_to_xfer(p)*dt
        deadcrootc_xfer(p) = deadcrootc_xfer(p) + &
          deadcrootc_storage_to_xfer(p)*dt
        gresp_storage(p) = gresp_storage(p) - gresp_storage_to_xfer(p)*dt
        gresp_xfer(p) = gresp_xfer(p) + gresp_storage_to_xfer(p)*dt
      end if
      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        ! lines here for consistency; the transfer terms are zero
        livestemc_storage(p) = livestemc_storage(p) - &
          livestemc_storage_to_xfer(p)*dt
        livestemc_xfer(p) = livestemc_xfer(p) + livestemc_storage_to_xfer(p)*dt
        grainc_storage(p) = grainc_storage(p) - grainc_storage_to_xfer(p)*dt
        grainc_xfer(p) = grainc_xfer(p) + grainc_storage_to_xfer(p)*dt
      end if
    end do ! end of pft loop
  end subroutine CStateUpdate1

#endif

end module mod_clm_cncstateupdate1
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
