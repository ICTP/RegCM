module mod_clm_cnnstateupdate2
#ifdef CN
  !
  ! Module for nitrogen state variable update, mortality fluxes.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams, only : dtsrf
  use mod_clm_varpar, only : nlevsoi, nlevdecomp

  implicit none
  save
  private

  public:: NStateUpdate2
  public:: NStateUpdate2h

  contains
  !
  ! On the radiation time step, update all the prognostic nitrogen state
  ! variables affected by gap-phase mortality fluxes
  !
  subroutine NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    implicit none
    integer(ik4), intent(in) :: num_soilc   ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_met_n(:,:)
    ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_cel_n(:,:)
    ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_lig_n(:,:)
    ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_cwdn(:,:)
    real(rk8), pointer, contiguous :: m_deadcrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_retransn_to_litter(:)

    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous :: decomp_npools_vr(:,:,:)
    real(rk8), pointer, contiguous :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rk8), pointer, contiguous :: deadcrootn_storage(:) ! (gN/m2) dead corse root N strg
    real(rk8), pointer, contiguous :: deadcrootn_xfer(:)    ! (gN/m2) dead corse root N trfr
    real(rk8), pointer, contiguous :: deadstemn(:)          ! (gN/m2) dead stem N
    real(rk8), pointer, contiguous :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
    real(rk8), pointer, contiguous :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
    real(rk8), pointer, contiguous :: frootn(:)             ! (gN/m2) fine root N
    real(rk8), pointer, contiguous :: frootn_storage(:)     ! (gN/m2) fine root N storage
    real(rk8), pointer, contiguous :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
    real(rk8), pointer, contiguous :: leafn(:)              ! (gN/m2) leaf N
    real(rk8), pointer, contiguous :: leafn_storage(:)      ! (gN/m2) leaf N storage
    real(rk8), pointer, contiguous :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
    real(rk8), pointer, contiguous :: livecrootn(:)         ! (gN/m2) live coarse root N
    real(rk8), pointer, contiguous :: livecrootn_storage(:) ! (gN/m2) live corse root N strg
    real(rk8), pointer, contiguous :: livecrootn_xfer(:)    ! (gN/m2) live corse root N trfr
    real(rk8), pointer, contiguous :: livestemn(:)          ! (gN/m2) live stem N
    real(rk8), pointer, contiguous :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rk8), pointer, contiguous :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rk8), pointer, contiguous :: retransn(:)           ! (gN/m2) plant pool of
                                                ! retranslocated N

    integer(ik4) :: c,p,j ! indices
    integer(ik4) :: fp,fc ! lake filter indices
    real(rk8):: dt        ! radiation time step (seconds)

    ! assign local pointers at the column level
    gap_mortality_n_to_litr_met_n  => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_met_n
    gap_mortality_n_to_litr_cel_n  => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_cel_n
    gap_mortality_n_to_litr_lig_n  => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_lig_n
    gap_mortality_n_to_cwdn        => clm3%g%l%c%cnf%gap_mortality_n_to_cwdn
    decomp_npools_vr               => clm3%g%l%c%cns%decomp_npools_vr
    ! assign local pointers at the pft level
    m_deadcrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
    m_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
    m_deadcrootn_xfer_to_litter    => &
            clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
    m_deadstemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
    m_deadstemn_to_litter          => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
    m_deadstemn_xfer_to_litter     => &
            clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
    m_frootn_storage_to_litter     => &
            clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
    m_frootn_to_litter             => clm3%g%l%c%p%pnf%m_frootn_to_litter
    m_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
    m_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
    m_leafn_to_litter              => clm3%g%l%c%p%pnf%m_leafn_to_litter
    m_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
    m_livecrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
    m_livecrootn_to_litter         => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
    m_livecrootn_xfer_to_litter    => &
            clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
    m_livestemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
    m_livestemn_to_litter          => clm3%g%l%c%p%pnf%m_livestemn_to_litter
    m_livestemn_xfer_to_litter     => &
            clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
    m_retransn_to_litter           => clm3%g%l%c%p%pnf%m_retransn_to_litter
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn

    ! set time steps
    dt = dtsrf

    do j = 1, nlevdecomp
      ! column loop
      do fc = 1, num_soilc
        c = filter_soilc(fc)

        ! column-level nitrogen fluxes from gap-phase mortality
        decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + &
                gap_mortality_n_to_litr_met_n(c,j) * dt
        decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + &
                gap_mortality_n_to_litr_cel_n(c,j) * dt
        decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + &
                gap_mortality_n_to_litr_lig_n(c,j) * dt
        decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + &
                gap_mortality_n_to_cwdn(c,j)  * dt

      end do ! end of column loop
    end do

    ! pft loop
    do fp = 1, num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from gap-phase mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - m_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter(p)  * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - m_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - &
              m_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - &
              m_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - &
              m_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - &
              m_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - &
              m_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - &
              m_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - &
              m_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - &
              m_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - &
              m_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - &
              m_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - &
              m_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - &
              m_deadcrootn_xfer_to_litter(p) * dt

    end do
  end subroutine NStateUpdate2
  !
  ! Update all the prognostic nitrogen state
  ! variables affected by harvest mortality fluxes
  !
  subroutine NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous :: harvest_n_to_litr_met_n(:,:)
    ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous :: harvest_n_to_litr_cel_n(:,:)
    ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous :: harvest_n_to_litr_lig_n(:,:)
    ! N fluxes associated with harvest to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous :: harvest_n_to_cwdn(:,:)
    real(rk8), pointer, contiguous :: hrv_deadcrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_deadcrootn_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_deadcrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_deadstemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_deadstemn_to_prod10n(:)
    real(rk8), pointer, contiguous :: hrv_deadstemn_to_prod100n(:)
    real(rk8), pointer, contiguous :: hrv_deadstemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_frootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_frootn_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_frootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_leafn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_leafn_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_leafn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livecrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livecrootn_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livecrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livestemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livestemn_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_livestemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: hrv_retransn_to_litter(:)

    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous :: decomp_npools_vr(:,:,:)
    real(rk8), pointer, contiguous :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rk8), pointer, contiguous :: deadcrootn_storage(:) ! (gN/m2) dead corse root N strg
    real(rk8), pointer, contiguous :: deadcrootn_xfer(:)    ! (gN/m2) dead corse root N trfr
    real(rk8), pointer, contiguous :: deadstemn(:)          ! (gN/m2) dead stem N
    real(rk8), pointer, contiguous :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
    real(rk8), pointer, contiguous :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
    real(rk8), pointer, contiguous :: frootn(:)             ! (gN/m2) fine root N
    real(rk8), pointer, contiguous :: frootn_storage(:)     ! (gN/m2) fine root N storage
    real(rk8), pointer, contiguous :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
    real(rk8), pointer, contiguous :: leafn(:)              ! (gN/m2) leaf N
    real(rk8), pointer, contiguous :: leafn_storage(:)      ! (gN/m2) leaf N storage
    real(rk8), pointer, contiguous :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
    real(rk8), pointer, contiguous :: livecrootn(:)         ! (gN/m2) live corse root N
    real(rk8), pointer, contiguous :: livecrootn_storage(:) ! (gN/m2) live corse root N strg
    real(rk8), pointer, contiguous :: livecrootn_xfer(:)    ! (gN/m2) live corse root N trfr
    real(rk8), pointer, contiguous :: livestemn(:)          ! (gN/m2) live stem N
    real(rk8), pointer, contiguous :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rk8), pointer, contiguous :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rk8), pointer, contiguous :: retransn(:)           ! (gN/m2) plant pool of
                                                ! retranslocated N

    integer(ik4) :: c,p,j ! indices
    integer(ik4) :: fp,fc ! lake filter indices
    real(rk8):: dt        ! radiation time step (seconds)

    ! assign local pointers at the column level
    harvest_n_to_litr_met_n => clm3%g%l%c%cnf%harvest_n_to_litr_met_n
    harvest_n_to_litr_cel_n => clm3%g%l%c%cnf%harvest_n_to_litr_cel_n
    harvest_n_to_litr_lig_n => clm3%g%l%c%cnf%harvest_n_to_litr_lig_n
    harvest_n_to_cwdn       => clm3%g%l%c%cnf%harvest_n_to_cwdn
    decomp_npools_vr        => clm3%g%l%c%cns%decomp_npools_vr

    ! assign local pointers at the pft level
    hrv_deadcrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
    hrv_deadcrootn_to_litter  => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
    hrv_deadcrootn_xfer_to_litter    => &
            clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter
    hrv_deadstemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
    hrv_deadstemn_to_prod10n  => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
    hrv_deadstemn_to_prod100n => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
    hrv_deadstemn_xfer_to_litter => &
            clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
    hrv_frootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
    hrv_frootn_to_litter       => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
    hrv_frootn_xfer_to_litter  => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
    hrv_leafn_storage_to_litter => &
            clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
    hrv_leafn_to_litter        => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
    hrv_leafn_xfer_to_litter   => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
    hrv_livecrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
    hrv_livecrootn_to_litter   => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
    hrv_livecrootn_xfer_to_litter    => &
            clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
    hrv_livestemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
    hrv_livestemn_to_litter => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
    hrv_livestemn_xfer_to_litter => &
            clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
    hrv_retransn_to_litter     => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
    deadcrootn                 => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage         => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer            => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                  => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage          => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer             => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                     => clm3%g%l%c%p%pns%frootn
    frootn_storage             => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                => clm3%g%l%c%p%pns%frootn_xfer
    leafn                      => clm3%g%l%c%p%pns%leafn
    leafn_storage              => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                 => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                 => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage         => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer            => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                  => clm3%g%l%c%p%pns%livestemn
    livestemn_storage          => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer             => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                   => clm3%g%l%c%p%pns%retransn

    ! set time steps
    dt = dtsrf

    do j = 1, nlevdecomp
      ! column loop
      do fc = 1,num_soilc
        c = filter_soilc(fc)

        ! column-level nitrogen fluxes from harvest mortality
        decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + &
                harvest_n_to_litr_met_n(c,j) * dt
        decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + &
                harvest_n_to_litr_cel_n(c,j) * dt
        decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + &
                harvest_n_to_litr_lig_n(c,j) * dt
        decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + &
                harvest_n_to_cwdn(c,j)  * dt

      end do ! end of column loop
    end do

    ! pft loop
    do fp = 1, num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from harvest mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - hrv_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - hrv_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - hrv_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod10n(p) * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod100n(p)* dt
      livecrootn(p) = livecrootn(p) - hrv_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - hrv_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - hrv_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - &
              hrv_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - &
              hrv_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - &
              hrv_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - &
              hrv_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - &
              hrv_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - &
              hrv_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - &
              hrv_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - &
              hrv_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - &
              hrv_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - &
              hrv_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - &
              hrv_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - &
              hrv_deadcrootn_xfer_to_litter(p) * dt

    end do
  end subroutine NStateUpdate2h

#endif

end module mod_clm_cnnstateupdate2
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
