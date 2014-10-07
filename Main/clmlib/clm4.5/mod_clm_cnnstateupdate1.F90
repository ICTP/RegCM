module mod_clm_cnnstateupdate1
#ifdef CN
  !
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams , only : dtsrf

  implicit none

  save

  private

  public:: NStateUpdate1

  contains
  !
  ! On the radiation time step, update all the prognostic nitrogen state
  ! variables (except for gap-phase mortality and fire fluxes)
  !
  subroutine NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_varpar , only : nlevdecomp, ndecomp_pools, &
            ndecomp_cascade_transitions
    use mod_clm_varpar , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#ifdef NITRIF_DENITRIF
    use mod_clm_varcon , only : nitrif_n2o_loss_frac
#endif
    use mod_clm_pftvarcon , only : npcropmin, nc3crop
    use mod_clm_surfrd , only : crop_prog
    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    integer(ik4) , pointer :: ivt(:) ! pft vegetation type
    ! binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8), pointer :: woody(:)
    real(rk8), pointer :: ndep_to_sminn(:)
    ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
    real(rk8), pointer :: nfix_to_sminn(:)
    real(rk8), pointer :: fert_to_sminn(:)
    real(rk8), pointer :: soyfixn_to_sminn(:)
    real(rk8), pointer :: sminn_to_denit_excess_vr(:,:)
    ! vertically-resolved denitrification along decomp cascade (gN/m3/s)
    real(rk8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)
    real(rk8), pointer :: sminn_to_plant_vr(:,:)
    real(rk8), pointer :: supplement_to_sminn_vr(:,:)
    real(rk8), pointer :: deadcrootn_storage_to_xfer(:)
    real(rk8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
    real(rk8), pointer :: deadstemn_storage_to_xfer(:)
    real(rk8), pointer :: deadstemn_xfer_to_deadstemn(:)
    real(rk8), pointer :: frootn_storage_to_xfer(:)
    real(rk8), pointer :: frootn_to_litter(:)
    real(rk8), pointer :: frootn_xfer_to_frootn(:)
    real(rk8), pointer :: frootn_to_retransn(:)
    real(rk8), pointer :: leafn_storage_to_xfer(:)
    real(rk8), pointer :: leafn_to_litter(:)
    real(rk8), pointer :: leafn_to_retransn(:)
    real(rk8), pointer :: leafn_xfer_to_leafn(:)
    real(rk8), pointer :: livecrootn_storage_to_xfer(:)
    real(rk8), pointer :: livecrootn_to_deadcrootn(:)
    real(rk8), pointer :: livecrootn_to_retransn(:)
    real(rk8), pointer :: livecrootn_xfer_to_livecrootn(:)
    real(rk8), pointer :: livestemn_storage_to_xfer(:)
    real(rk8), pointer :: livestemn_to_deadstemn(:)
    real(rk8), pointer :: livestemn_to_retransn(:)
    real(rk8), pointer :: livestemn_xfer_to_livestemn(:)
    real(rk8), pointer :: npool_to_deadcrootn(:)
    real(rk8), pointer :: npool_to_deadcrootn_storage(:)
    real(rk8), pointer :: npool_to_deadstemn(:)
    real(rk8), pointer :: npool_to_deadstemn_storage(:)
    real(rk8), pointer :: npool_to_frootn(:)
    real(rk8), pointer :: npool_to_frootn_storage(:)
    real(rk8), pointer :: npool_to_leafn(:)
    real(rk8), pointer :: npool_to_leafn_storage(:)
    real(rk8), pointer :: npool_to_livecrootn(:)
    real(rk8), pointer :: npool_to_livecrootn_storage(:)
    ! allocation to live stem N (gN/m2/s)
    real(rk8), pointer :: npool_to_livestemn(:)
    ! allocation to live stem N storage (gN/m2/s)
    real(rk8), pointer :: npool_to_livestemn_storage(:)
    ! deployment of retranslocated N (gN/m2/s)
    real(rk8), pointer :: retransn_to_npool(:)
    ! deployment of soil mineral N uptake (gN/m2/s)
    real(rk8), pointer :: sminn_to_npool(:)
    ! grain N shift storage to transfer (gN/m2/s)
    real(rk8), pointer :: grainn_storage_to_xfer(:)
    ! grain N to food (gN/m2/s)
    real(rk8), pointer :: grainn_to_food(:)
    ! grain N growth from storage (gN/m2/s)
    real(rk8), pointer :: grainn_xfer_to_grainn(:)
    ! livestem N to litter (gN/m2/s)
    real(rk8), pointer :: livestemn_to_litter(:)
    ! allocation to grain N (gN/m2/s)
    real(rk8), pointer :: npool_to_grainn(:)
    ! allocation to grain N storage (gN/m2/s)
    real(rk8), pointer :: npool_to_grainn_storage(:)

    real(rk8), pointer :: sminn_vr(:,:)    ! (gN/m3) soil mineral N
#ifdef NITRIF_DENITRIF
    real(rk8), pointer :: smin_no3_vr(:,:) ! (gN/m3) soil NO3
    real(rk8), pointer :: smin_nh4_vr(:,:) ! (gN/m3) soil NH4
    real(rk8), pointer :: f_nit_vr(:,:)    ! (gN/m3/s) soil nitrification flux
    real(rk8), pointer :: f_denit_vr(:,:)  ! (gN/m3/s) soil denitrification flux
    real(rk8), pointer :: actual_immob_no3_vr(:,:)  ! (gN/m3/s)
    real(rk8), pointer :: actual_immob_nh4_vr(:,:)  ! (gN/m3/s)
    real(rk8), pointer :: smin_no3_to_plant_vr(:,:) ! (gN/m3/s)
    real(rk8), pointer :: smin_nh4_to_plant_vr(:,:) ! (gN/m3/s)
    real(rk8), pointer :: gross_nmin_vr(:,:)        ! (gN/m3/s)
#endif
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer :: decomp_npools_vr(:,:,:)
    ! (gC/m3)  change in decomposing N pools over a timestep.
    ! Used to update concentrations concurrently with vertical transport
    real(rk8), pointer :: decomp_npools_sourcesink(:,:,:)
    ! vert-res transfer of N from donor to receiver pool along
    ! decomp. cascade (gN/m3/s)
    real(rk8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)
    ! vert-res mineral N flux for transition along decomposition
    ! cascade (gN/m3/s)
    real(rk8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer :: cascade_receiver_pool(:)
    ! profile over which N deposition is distributed through column (1/m)
    real(rk8), pointer :: ndep_prof(:,:)
    ! profile over which N fixation is distributed through column (1/m)
    real(rk8), pointer :: nfixation_prof(:,:)
    real(rk8), pointer :: grainn(:)             ! (gN/m2) grain N
    real(rk8), pointer :: grainn_storage(:)     ! (gN/m2) grain N storage
    real(rk8), pointer :: grainn_xfer(:)        ! (gN/m2) grain N transfer
    real(rk8), pointer :: frootn(:)             ! (gN/m2) fine root N
    real(rk8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
    real(rk8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
    real(rk8), pointer :: leafn(:)              ! (gN/m2) leaf N
    real(rk8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
    real(rk8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
    real(rk8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
    real(rk8), pointer :: livecrootn_storage(:) ! (gN/m2) live corse root N strg
    real(rk8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live corse root N trfr
    real(rk8), pointer :: livestemn(:)          ! (gN/m2) live stem N
    real(rk8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rk8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rk8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rk8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead corse root N strg
    real(rk8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead corse root N trfr
    real(rk8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
    real(rk8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
    real(rk8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
    real(rk8), pointer :: retransn(:)           ! (gN/m2) plant pool of
                                                ! retranslocated N
    real(rk8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
    ! N fluxes associated with phenology (litterfall and crop) to
    ! litter metabolic pool (gN/m3/s)
    real(rk8), pointer :: phenology_n_to_litr_met_n(:,:)
    ! N fluxes associated with phenology (litterfall and crop) to
    ! litter cellulose pool (gN/m3/s)
    real(rk8), pointer :: phenology_n_to_litr_cel_n(:,:)
    ! N fluxes associated with phenology (litterfall and crop) to
    ! litter lignin pool (gN/m3/s)
    real(rk8), pointer :: phenology_n_to_litr_lig_n(:,:)

    real(rk8), pointer :: dwt_seedn_to_leaf(:)
    real(rk8), pointer :: dwt_seedn_to_deadstem(:)
    real(rk8), pointer :: dwt_frootn_to_litr_met_n(:,:)
    real(rk8), pointer :: dwt_frootn_to_litr_cel_n(:,:)
    real(rk8), pointer :: dwt_frootn_to_litr_lig_n(:,:)
    real(rk8), pointer :: dwt_livecrootn_to_cwdn(:,:)
    real(rk8), pointer :: dwt_deadcrootn_to_cwdn(:,:)
    real(rk8), pointer :: seedn(:)

    integer(ik4) :: c,p,j,k  ! indices
    integer(ik4) :: fp,fc    ! lake filter indices
    real(rk8):: dt       ! radiation time step (seconds)

    ! assign local pointers
    woody      => pftcon%woody

    ! assign local pointers at the column level
    ndep_to_sminn    => clm3%g%l%c%cnf%ndep_to_sminn
    nfix_to_sminn    => clm3%g%l%c%cnf%nfix_to_sminn
    fert_to_sminn    => clm3%g%l%c%cnf%fert_to_sminn
    soyfixn_to_sminn => clm3%g%l%c%cnf%soyfixn_to_sminn
#ifndef NITRIF_DENITRIF
    sminn_to_denit_excess_vr         => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
    sminn_to_denit_decomp_cascade_vr => &
            clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr
#else
    smin_no3_vr          => clm3%g%l%c%cns%smin_no3_vr
    smin_nh4_vr          => clm3%g%l%c%cns%smin_nh4_vr
    f_nit_vr             => clm3%g%l%c%cnf%f_nit_vr
    f_denit_vr           => clm3%g%l%c%cnf%f_denit_vr
    actual_immob_no3_vr  => clm3%g%l%c%cnf%actual_immob_no3_vr
    actual_immob_nh4_vr  => clm3%g%l%c%cnf%actual_immob_nh4_vr
    smin_no3_to_plant_vr => clm3%g%l%c%cnf%smin_no3_to_plant_vr
    smin_nh4_to_plant_vr => clm3%g%l%c%cnf%smin_nh4_to_plant_vr
    gross_nmin_vr        => clm3%g%l%c%cnf%gross_nmin_vr
#endif
    sminn_to_plant_vr            => clm3%g%l%c%cnf%sminn_to_plant_vr
    decomp_cascade_sminn_flux_vr => &
            clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr
    decomp_cascade_ntransfer_vr => &
            clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr
    cascade_donor_pool        => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool     => decomp_cascade_con%cascade_receiver_pool
    supplement_to_sminn_vr    => clm3%g%l%c%cnf%supplement_to_sminn_vr
    decomp_npools_vr          => clm3%g%l%c%cns%decomp_npools_vr
    decomp_npools_sourcesink  => clm3%g%l%c%cnf%decomp_npools_sourcesink
    sminn_vr                  => clm3%g%l%c%cns%sminn_vr
    ndep_prof                 => clm3%g%l%c%cps%ndep_prof
    nfixation_prof            => clm3%g%l%c%cps%nfixation_prof
    phenology_n_to_litr_met_n => clm3%g%l%c%cnf%phenology_n_to_litr_met_n
    phenology_n_to_litr_cel_n => clm3%g%l%c%cnf%phenology_n_to_litr_cel_n
    phenology_n_to_litr_lig_n => clm3%g%l%c%cnf%phenology_n_to_litr_lig_n

    ! new pointers for dynamic landcover
    dwt_seedn_to_leaf        => clm3%g%l%c%cnf%dwt_seedn_to_leaf
    dwt_seedn_to_deadstem    => clm3%g%l%c%cnf%dwt_seedn_to_deadstem
    dwt_frootn_to_litr_met_n => clm3%g%l%c%cnf%dwt_frootn_to_litr_met_n
    dwt_frootn_to_litr_cel_n => clm3%g%l%c%cnf%dwt_frootn_to_litr_cel_n
    dwt_frootn_to_litr_lig_n => clm3%g%l%c%cnf%dwt_frootn_to_litr_lig_n
    dwt_livecrootn_to_cwdn   => clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn
    dwt_deadcrootn_to_cwdn   => clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn
    seedn => clm3%g%l%c%cns%seedn

    ! assign local pointers at the pft level
    ivt                         => clm3%g%l%c%p%itype
    deadcrootn_storage_to_xfer  => clm3%g%l%c%p%pnf%deadcrootn_storage_to_xfer
    deadcrootn_xfer_to_deadcrootn  => &
            clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn
    deadstemn_storage_to_xfer      => clm3%g%l%c%p%pnf%deadstemn_storage_to_xfer
    deadstemn_xfer_to_deadstemn    => &
            clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn
    frootn_storage_to_xfer     => clm3%g%l%c%p%pnf%frootn_storage_to_xfer
    frootn_to_litter           => clm3%g%l%c%p%pnf%frootn_to_litter
    frootn_to_retransn         => clm3%g%l%c%p%pnf%frootn_to_retransn
    frootn_xfer_to_frootn      => clm3%g%l%c%p%pnf%frootn_xfer_to_frootn
    leafn_storage_to_xfer      => clm3%g%l%c%p%pnf%leafn_storage_to_xfer
    leafn_to_litter            => clm3%g%l%c%p%pnf%leafn_to_litter
    leafn_to_retransn          => clm3%g%l%c%p%pnf%leafn_to_retransn
    leafn_xfer_to_leafn        => clm3%g%l%c%p%pnf%leafn_xfer_to_leafn
    livecrootn_storage_to_xfer => clm3%g%l%c%p%pnf%livecrootn_storage_to_xfer
    livecrootn_to_deadcrootn   => clm3%g%l%c%p%pnf%livecrootn_to_deadcrootn
    livecrootn_to_retransn     => clm3%g%l%c%p%pnf%livecrootn_to_retransn
    livecrootn_xfer_to_livecrootn  => &
            clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn
    livestemn_storage_to_xfer  => clm3%g%l%c%p%pnf%livestemn_storage_to_xfer
    livestemn_to_deadstemn     => clm3%g%l%c%p%pnf%livestemn_to_deadstemn
    livestemn_to_retransn      => clm3%g%l%c%p%pnf%livestemn_to_retransn
    livestemn_xfer_to_livestemn => &
            clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn
    npool_to_deadcrootn         => clm3%g%l%c%p%pnf%npool_to_deadcrootn
    npool_to_deadcrootn_storage => clm3%g%l%c%p%pnf%npool_to_deadcrootn_storage
    npool_to_deadstemn          => clm3%g%l%c%p%pnf%npool_to_deadstemn
    npool_to_deadstemn_storage  => clm3%g%l%c%p%pnf%npool_to_deadstemn_storage
    npool_to_frootn             => clm3%g%l%c%p%pnf%npool_to_frootn
    npool_to_frootn_storage     => clm3%g%l%c%p%pnf%npool_to_frootn_storage
    npool_to_leafn              => clm3%g%l%c%p%pnf%npool_to_leafn
    npool_to_leafn_storage      => clm3%g%l%c%p%pnf%npool_to_leafn_storage
    npool_to_livecrootn         => clm3%g%l%c%p%pnf%npool_to_livecrootn
    npool_to_livecrootn_storage => clm3%g%l%c%p%pnf%npool_to_livecrootn_storage
    npool_to_livestemn          => clm3%g%l%c%p%pnf%npool_to_livestemn
    npool_to_livestemn_storage  => clm3%g%l%c%p%pnf%npool_to_livestemn_storage
    retransn_to_npool           => clm3%g%l%c%p%pnf%retransn_to_npool
    sminn_to_npool              => clm3%g%l%c%p%pnf%sminn_to_npool
    grainn_storage_to_xfer      => clm3%g%l%c%p%pnf%grainn_storage_to_xfer
    grainn_to_food              => clm3%g%l%c%p%pnf%grainn_to_food
    grainn_xfer_to_grainn       => clm3%g%l%c%p%pnf%grainn_xfer_to_grainn
    livestemn_to_litter         => clm3%g%l%c%p%pnf%livestemn_to_litter
    npool_to_grainn             => clm3%g%l%c%p%pnf%npool_to_grainn
    npool_to_grainn_storage     => clm3%g%l%c%p%pnf%npool_to_grainn_storage
    grainn                      => clm3%g%l%c%p%pns%grainn
    grainn_storage              => clm3%g%l%c%p%pns%grainn_storage
    grainn_xfer                 => clm3%g%l%c%p%pns%grainn_xfer
    deadcrootn                  => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage          => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer             => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                   => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage           => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer              => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                      => clm3%g%l%c%p%pns%frootn
    frootn_storage              => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                 => clm3%g%l%c%p%pns%frootn_xfer
    leafn                       => clm3%g%l%c%p%pns%leafn
    leafn_storage               => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                  => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                  => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage          => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer             => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                   => clm3%g%l%c%p%pns%livestemn
    livestemn_storage           => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer              => clm3%g%l%c%p%pns%livestemn_xfer
    npool                       => clm3%g%l%c%p%pns%npool
    retransn                    => clm3%g%l%c%p%pns%retransn

    ! set time steps
    dt = dtsrf

    ! column-level fluxes

    ! column loop
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      ! seeding fluxes, from dynamic landcover
      seedn(c) = seedn(c) - dwt_seedn_to_leaf(c) * dt
      seedn(c) = seedn(c) - dwt_seedn_to_deadstem(c) * dt
    end do

    do j = 1 , nlevdecomp
      ! column loop
      do fc = 1 , num_soilc
        c = filter_soilc(fc)

#ifndef NITRIF_DENITRIF
        ! N deposition and fixation
        sminn_vr(c,j) = sminn_vr(c,j) + ndep_to_sminn(c)*dt * ndep_prof(c,j)
        sminn_vr(c,j) = sminn_vr(c,j) + nfix_to_sminn(c)*dt * &
                nfixation_prof(c,j)
#else
        ! N deposition and fixation (put all into NH4 pool)
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + ndep_to_sminn(c)*dt * &
                ndep_prof(c,j)
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + &
                nfix_to_sminn(c)*dt * nfixation_prof(c,j)
#endif

        ! plant to litter fluxes
        ! phenology and dynamic landcover fluxes
        decomp_npools_sourcesink(c,j,i_met_lit) = &
              ( phenology_n_to_litr_met_n(c,j) + &
                dwt_frootn_to_litr_met_n(c,j) ) *dt
        decomp_npools_sourcesink(c,j,i_cel_lit) = &
              ( phenology_n_to_litr_cel_n(c,j) + &
              dwt_frootn_to_litr_cel_n(c,j) ) *dt
        decomp_npools_sourcesink(c,j,i_lig_lit) = &
              ( phenology_n_to_litr_lig_n(c,j) + &
              dwt_frootn_to_litr_lig_n(c,j) ) *dt
        decomp_npools_sourcesink(c,j,i_cwd) = &
              ( dwt_livecrootn_to_cwdn(c,j) + &
              dwt_deadcrootn_to_cwdn(c,j) )*dt
      end do
    end do

    ! repeating N dep and fixation for crops
    if ( crop_prog )then
      do j = 1 , nlevdecomp
        ! column loop
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
#ifndef NITRIF_DENITRIF
          ! N deposition and fixation
          sminn_vr(c,j) = sminn_vr(c,j) + fert_to_sminn(c)*dt * ndep_prof(c,j)
          sminn_vr(c,j) = sminn_vr(c,j) + &
                  soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
#else
          ! N deposition and fixation (put all into NH4 pool)
          smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + &
                  fert_to_sminn(c)*dt * ndep_prof(c,j)
          smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + &
                  soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
#endif
        end do
      end do
    end if

    ! decomposition fluxes
    do k = 1 , ndecomp_cascade_transitions
      do j = 1 , nlevdecomp
        ! column loop
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = &
                 decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                 decomp_cascade_ntransfer_vr(c,j,k) * dt
        end do
      end do
    end do
    do k = 1 , ndecomp_cascade_transitions
      if ( cascade_receiver_pool(k) .ne. 0 ) then  ! skip terminal transitions
        do j = 1, nlevdecomp
          ! column loop
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) = &
                    decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) + &
                    (decomp_cascade_ntransfer_vr(c,j,k) + &
                     decomp_cascade_sminn_flux_vr(c,j,k)) * dt
          end do
        end do
      else  ! terminal transitions
        do j = 1 , nlevdecomp
          ! column loop
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = &
                   decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                   decomp_cascade_sminn_flux_vr(c,j,k) * dt
          end do
        end do
      end if
    end do

#ifndef NITRIF_DENITRIF
    ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes
    ! and denitrification fluxes
    do k = 1 , ndecomp_cascade_transitions
      if ( cascade_receiver_pool(k) .ne. 0 ) then  ! skip terminal transitions
        do j = 1 , nlevdecomp
          ! column loop
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            sminn_vr(c,j) = sminn_vr(c,j) - &
               (sminn_to_denit_decomp_cascade_vr(c,j,k) + &
                decomp_cascade_sminn_flux_vr(c,j,k))* dt
          end do
        end do
      else
        do j = 1 , nlevdecomp
          ! column loop
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            sminn_vr(c,j) = sminn_vr(c,j) - &
               sminn_to_denit_decomp_cascade_vr(c,j,k)* dt
               sminn_vr(c,j)  = sminn_vr(c,j) + &
               decomp_cascade_sminn_flux_vr(c,j,k)* dt
          end do
        end do
      endif
    end do

    do j = 1 , nlevdecomp
      ! column loop
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        ! "bulk denitrification"
        sminn_vr(c,j) = sminn_vr(c,j) - sminn_to_denit_excess_vr(c,j) * dt

        ! total plant uptake from mineral N
        sminn_vr(c,j) = sminn_vr(c,j) - sminn_to_plant_vr(c,j)*dt

        ! flux that prevents N limitation (when Carbon_only is set)
        sminn_vr(c,j) = sminn_vr(c,j) + supplement_to_sminn_vr(c,j)*dt
      end do
    end do
#else
    !-------------    NITRIF_DENITRIF --------------------

    do j = 1 , nlevdecomp
      ! column loop
      do fc = 1 , num_soilc
        c = filter_soilc(fc)

        ! mineralization fluxes (divert a fraction of this stream to
        ! nitrification flux, add the rest to NH4 pool)
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + gross_nmin_vr(c,j)*dt

        ! immobilization fluxes
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - actual_immob_nh4_vr(c,j)*dt
        smin_no3_vr(c,j) = smin_no3_vr(c,j) - actual_immob_no3_vr(c,j)*dt

        ! plant uptake fluxes
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - smin_nh4_to_plant_vr(c,j)*dt
        smin_no3_vr(c,j) = smin_no3_vr(c,j) - smin_no3_to_plant_vr(c,j)*dt

        ! Account for nitrification fluxes
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - f_nit_vr(c,j) * dt
        smin_no3_vr(c,j) = smin_no3_vr(c,j) + f_nit_vr(c,j) * dt * &
                (1.D0 - nitrif_n2o_loss_frac)
        ! Account for denitrification fluxes
        smin_no3_vr(c,j) = smin_no3_vr(c,j) - f_denit_vr(c,j) * dt

        ! flux that prevents N limitation (when Carbon_only is set;
        !  put all into NH4)
        smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + supplement_to_sminn_vr(c,j)*dt

        ! update diagnostic total
        sminn_vr(c,j) = smin_nh4_vr(c,j) + smin_no3_vr(c,j)
      end do ! end of column loop
    end do
#endif

    ! pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      ! phenology: transfer growth fluxes
      leafn(p)       = leafn(p)       + leafn_xfer_to_leafn(p)*dt
      leafn_xfer(p)  = leafn_xfer(p)  - leafn_xfer_to_leafn(p)*dt
      frootn(p)      = frootn(p)      + frootn_xfer_to_frootn(p)*dt
      frootn_xfer(p) = frootn_xfer(p) - frootn_xfer_to_frootn(p)*dt
      if (woody(ivt(p)) == 1.0D0) then
        livestemn(p)       = livestemn(p)       + &
                livestemn_xfer_to_livestemn(p)*dt
        livestemn_xfer(p)  = livestemn_xfer(p)  - &
                livestemn_xfer_to_livestemn(p)*dt
        deadstemn(p)       = deadstemn(p)       + &
                deadstemn_xfer_to_deadstemn(p)*dt
        deadstemn_xfer(p)  = deadstemn_xfer(p)  - &
                deadstemn_xfer_to_deadstemn(p)*dt
        livecrootn(p)      = livecrootn(p)      + &
                livecrootn_xfer_to_livecrootn(p)*dt
        livecrootn_xfer(p) = livecrootn_xfer(p) - &
                livecrootn_xfer_to_livecrootn(p)*dt
        deadcrootn(p)      = deadcrootn(p)      + &
                deadcrootn_xfer_to_deadcrootn(p)*dt
        deadcrootn_xfer(p) = deadcrootn_xfer(p) - &
                deadcrootn_xfer_to_deadcrootn(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
        ! lines here for consistency; the transfer terms are zero
        livestemn(p)       = livestemn(p)      + &
                 livestemn_xfer_to_livestemn(p)*dt
        livestemn_xfer(p)  = livestemn_xfer(p) - &
                 livestemn_xfer_to_livestemn(p)*dt
        grainn(p)          = grainn(p)         + &
                 grainn_xfer_to_grainn(p)*dt
        grainn_xfer(p)     = grainn_xfer(p)    - &
                 grainn_xfer_to_grainn(p)*dt
      end if

      ! phenology: litterfall and retranslocation fluxes
      leafn(p)    = leafn(p)    - leafn_to_litter(p)*dt
      frootn(p)   = frootn(p)   - frootn_to_litter(p)*dt
      leafn(p)    = leafn(p)    - leafn_to_retransn(p)*dt
      retransn(p) = retransn(p) + leafn_to_retransn(p)*dt

      ! live wood turnover and retranslocation fluxes
      if (woody(ivt(p)) == 1.D0) then
        livestemn(p)  = livestemn(p)  - livestemn_to_deadstemn(p)*dt
        deadstemn(p)  = deadstemn(p)  + livestemn_to_deadstemn(p)*dt
        livestemn(p)  = livestemn(p)  - livestemn_to_retransn(p)*dt
        retransn(p)   = retransn(p)   + livestemn_to_retransn(p)*dt
        livecrootn(p) = livecrootn(p) - livecrootn_to_deadcrootn(p)*dt
        deadcrootn(p) = deadcrootn(p) + livecrootn_to_deadcrootn(p)*dt
        livecrootn(p) = livecrootn(p) - livecrootn_to_retransn(p)*dt
        retransn(p)   = retransn(p)   + livecrootn_to_retransn(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
        frootn(p)     = frootn(p)     - frootn_to_retransn(p)*dt
        retransn(p)   = retransn(p)   + frootn_to_retransn(p)*dt
        livestemn(p)  = livestemn(p)  - livestemn_to_litter(p)*dt
        livestemn(p)  = livestemn(p)  - livestemn_to_retransn(p)*dt
        retransn(p)   = retransn(p)   + livestemn_to_retransn(p)*dt
        grainn(p)     = grainn(p)     - grainn_to_food(p)*dt
      end if

      ! uptake from soil mineral N pool
      npool(p) = npool(p) + sminn_to_npool(p)*dt

      ! deployment from retranslocation pool
      npool(p)    = npool(p)    + retransn_to_npool(p)*dt
      retransn(p) = retransn(p) - retransn_to_npool(p)*dt

      ! allocation fluxes
      npool(p)           = npool(p)          - npool_to_leafn(p)*dt
      leafn(p)           = leafn(p)          + npool_to_leafn(p)*dt
      npool(p)           = npool(p)          - npool_to_leafn_storage(p)*dt
      leafn_storage(p)   = leafn_storage(p)  + npool_to_leafn_storage(p)*dt
      npool(p)           = npool(p)          - npool_to_frootn(p)*dt
      frootn(p)          = frootn(p)         + npool_to_frootn(p)*dt
      npool(p)           = npool(p)          - npool_to_frootn_storage(p)*dt
      frootn_storage(p)  = frootn_storage(p) + npool_to_frootn_storage(p)*dt
      if (woody(ivt(p)) == 1.D0) then
        npool(p)              = npool(p)              - &
                npool_to_livestemn(p)*dt
        livestemn(p)          = livestemn(p)          + &
                npool_to_livestemn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_livestemn_storage(p)*dt
        livestemn_storage(p)  = livestemn_storage(p)  + &
                npool_to_livestemn_storage(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_deadstemn(p)*dt
        deadstemn(p)          = deadstemn(p)          + &
                npool_to_deadstemn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_deadstemn_storage(p)*dt
        deadstemn_storage(p)  = deadstemn_storage(p)  + &
                npool_to_deadstemn_storage(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_livecrootn(p)*dt
        livecrootn(p)         = livecrootn(p)         + &
                npool_to_livecrootn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_livecrootn_storage(p)*dt
        livecrootn_storage(p) = livecrootn_storage(p) + &
                npool_to_livecrootn_storage(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_deadcrootn(p)*dt
        deadcrootn(p)         = deadcrootn(p)         + &
                npool_to_deadcrootn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_deadcrootn_storage(p)*dt
        deadcrootn_storage(p) = deadcrootn_storage(p) + &
                npool_to_deadcrootn_storage(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
        npool(p)              = npool(p)              - &
                npool_to_livestemn(p)*dt
        livestemn(p)          = livestemn(p)          + &
                npool_to_livestemn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_livestemn_storage(p)*dt
        livestemn_storage(p)  = livestemn_storage(p)  + &
                npool_to_livestemn_storage(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_grainn(p)*dt
        grainn(p)             = grainn(p)             + &
                npool_to_grainn(p)*dt
        npool(p)              = npool(p)              - &
                npool_to_grainn_storage(p)*dt
        grainn_storage(p)     = grainn_storage(p)     + &
                npool_to_grainn_storage(p)*dt
      end if

      ! move storage pools into transfer pools
      leafn_storage(p)  = leafn_storage(p)  - leafn_storage_to_xfer(p)*dt
      leafn_xfer(p)     = leafn_xfer(p)     + leafn_storage_to_xfer(p)*dt
      frootn_storage(p) = frootn_storage(p) - frootn_storage_to_xfer(p)*dt
      frootn_xfer(p)    = frootn_xfer(p)    + frootn_storage_to_xfer(p)*dt
      if (woody(ivt(p)) == 1.D0) then
        livestemn_storage(p)  = livestemn_storage(p)  - &
                livestemn_storage_to_xfer(p)*dt
        livestemn_xfer(p)     = livestemn_xfer(p)     + &
                livestemn_storage_to_xfer(p)*dt
        deadstemn_storage(p)  = deadstemn_storage(p)  - &
                deadstemn_storage_to_xfer(p)*dt
        deadstemn_xfer(p)     = deadstemn_xfer(p)     + &
                deadstemn_storage_to_xfer(p)*dt
        livecrootn_storage(p) = livecrootn_storage(p) - &
                livecrootn_storage_to_xfer(p)*dt
        livecrootn_xfer(p)    = livecrootn_xfer(p)    + &
                livecrootn_storage_to_xfer(p)*dt
        deadcrootn_storage(p) = deadcrootn_storage(p) - &
                deadcrootn_storage_to_xfer(p)*dt
        deadcrootn_xfer(p)    = deadcrootn_xfer(p)    + &
                deadcrootn_storage_to_xfer(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
        ! lines here for consistency; the transfer terms are zero
        livestemn_storage(p)  = livestemn_storage(p) - &
                livestemn_storage_to_xfer(p)*dt
        livestemn_xfer(p)     = livestemn_xfer(p)    + &
                livestemn_storage_to_xfer(p)*dt
        grainn_storage(p)     = grainn_storage(p)    - &
                grainn_storage_to_xfer(p)*dt
        grainn_xfer(p)        = grainn_xfer(p)       + &
                grainn_storage_to_xfer(p)*dt
      end if
    end do
  end subroutine NStateUpdate1

#endif

end module mod_clm_cnnstateupdate1
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
