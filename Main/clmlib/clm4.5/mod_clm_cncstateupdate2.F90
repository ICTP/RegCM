module mod_clm_cncstateupdate2
#ifdef CN
  !
  ! Module for carbon state variable update, mortality fluxes.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams , only : dtsrf
  use mod_mpmessage

  implicit none

  save

  private

  public:: CStateUpdate2
  public:: CStateUpdate2h

  contains
  !
  ! On the radiation time step, update all the prognostic carbon state
  ! variables affected by gap-phase mortality fluxes
  !
  subroutine CStateUpdate2(num_soilc, filter_soilc, &
                           num_soilp, filter_soilp, isotope)
    use mod_clm_type
    use mod_clm_varpar , only : nlevdecomp
    use mod_clm_varpar , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    implicit none
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(:)
    ! number of soil pfts in filter
    integer(ik4), intent(in) :: num_soilp
    ! filter for soil pfts
    integer(ik4), intent(in) :: filter_soilp(:)
    character(len=*), intent(in) :: isotope  ! 'bulk', 'c13' or 'c14'

    ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
    real(rk8), pointer :: gap_mortality_c_to_litr_met_c(:,:)
    ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
    real(rk8), pointer :: gap_mortality_c_to_litr_cel_c(:,:)
    ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
    real(rk8), pointer :: gap_mortality_c_to_litr_lig_c(:,:)
    ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
    real(rk8), pointer :: gap_mortality_c_to_cwdc(:,:)
    real(rk8), pointer :: m_deadcrootc_storage_to_litter(:)
    real(rk8), pointer :: m_deadcrootc_to_litter(:)
    real(rk8), pointer :: m_deadcrootc_xfer_to_litter(:)
    real(rk8), pointer :: m_deadstemc_storage_to_litter(:)
    real(rk8), pointer :: m_deadstemc_to_litter(:)
    real(rk8), pointer :: m_deadstemc_xfer_to_litter(:)
    real(rk8), pointer :: m_frootc_storage_to_litter(:)
    real(rk8), pointer :: m_frootc_to_litter(:)
    real(rk8), pointer :: m_frootc_xfer_to_litter(:)
    real(rk8), pointer :: m_gresp_storage_to_litter(:)
    real(rk8), pointer :: m_gresp_xfer_to_litter(:)
    real(rk8), pointer :: m_leafc_storage_to_litter(:)
    real(rk8), pointer :: m_leafc_to_litter(:)
    real(rk8), pointer :: m_leafc_xfer_to_litter(:)
    real(rk8), pointer :: m_livecrootc_storage_to_litter(:)
    real(rk8), pointer :: m_livecrootc_to_litter(:)
    real(rk8), pointer :: m_livecrootc_xfer_to_litter(:)
    real(rk8), pointer :: m_livestemc_storage_to_litter(:)
    real(rk8), pointer :: m_livestemc_to_litter(:)
    real(rk8), pointer :: m_livestemc_xfer_to_litter(:)

    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer :: decomp_cpools_vr(:,:,:)
    ! (gC/m2) dead coarse root C
    real(rk8), pointer :: deadcrootc(:)
    ! (gC/m2) dead coarse root C storage
    real(rk8), pointer :: deadcrootc_storage(:)
    !(gC/m2) dead coarse root C transfer
    real(rk8), pointer :: deadcrootc_xfer(:)
    real(rk8), pointer :: deadstemc(:)         ! (gC/m2) dead stem C
    real(rk8), pointer :: deadstemc_storage(:) ! (gC/m2) dead stem C storage
    real(rk8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
    real(rk8), pointer :: frootc(:)          ! (gC/m2) fine root C
    real(rk8), pointer :: frootc_storage(:)  ! (gC/m2) fine root C storage
    real(rk8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
    real(rk8), pointer :: gresp_storage(:) ! (gC/m2) growth respiration storage
    real(rk8), pointer :: gresp_xfer(:)    ! (gC/m2) growth respiration transfer
    real(rk8), pointer :: leafc(:)          ! (gC/m2) leaf C
    real(rk8), pointer :: leafc_storage(:)  ! (gC/m2) leaf C storage
    real(rk8), pointer :: leafc_xfer(:)     ! (gC/m2) leaf C transfer
    real(rk8), pointer :: livecrootc(:)     ! (gC/m2) live coarse root C
    ! (gC/m2) live coarse root C storage
    real(rk8), pointer :: livecrootc_storage(:)
    !(gC/m2) live coarse root C transfer
    real(rk8), pointer :: livecrootc_xfer(:)
    real(rk8), pointer :: livestemc(:)       ! (gC/m2) live stem C
    ! (gC/m2) live stem C storage
    real(rk8), pointer :: livestemc_storage(:)
    real(rk8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer

    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    type(column_cflux_type), pointer :: ccisof
    type(column_cstate_type), pointer :: ccisos
    integer(ik4) :: c,p,j      ! indices
    integer(ik4) :: fp,fc    ! lake filter indices
    real(rk8):: dt       ! radiation time step (seconds)

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
        call fatal(__FILE__,__LINE__, &
         'CNCIsoStateUpdate2Mod: iso must be bulk, c13 or c14')
    end select

    ! assign local pointers at the column level
    gap_mortality_c_to_litr_met_c  => ccisof%gap_mortality_c_to_litr_met_c
    gap_mortality_c_to_litr_cel_c  => ccisof%gap_mortality_c_to_litr_cel_c
    gap_mortality_c_to_litr_lig_c  => ccisof%gap_mortality_c_to_litr_lig_c
    gap_mortality_c_to_cwdc        => ccisof%gap_mortality_c_to_cwdc
    decomp_cpools_vr               => ccisos%decomp_cpools_vr

    ! assign local pointers at the pft level
    m_deadcrootc_storage_to_litter => pcisof%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_litter         => pcisof%m_deadcrootc_to_litter
    m_deadcrootc_xfer_to_litter    => pcisof%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_litter  => pcisof%m_deadstemc_storage_to_litter
    m_deadstemc_to_litter          => pcisof%m_deadstemc_to_litter
    m_deadstemc_xfer_to_litter     => pcisof%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_litter     => pcisof%m_frootc_storage_to_litter
    m_frootc_to_litter             => pcisof%m_frootc_to_litter
    m_frootc_xfer_to_litter        => pcisof%m_frootc_xfer_to_litter
    m_gresp_storage_to_litter      => pcisof%m_gresp_storage_to_litter
    m_gresp_xfer_to_litter         => pcisof%m_gresp_xfer_to_litter
    m_leafc_storage_to_litter      => pcisof%m_leafc_storage_to_litter
    m_leafc_to_litter              => pcisof%m_leafc_to_litter
    m_leafc_xfer_to_litter         => pcisof%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_litter => pcisof%m_livecrootc_storage_to_litter
    m_livecrootc_to_litter         => pcisof%m_livecrootc_to_litter
    m_livecrootc_xfer_to_litter    => pcisof%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_litter  => pcisof%m_livestemc_storage_to_litter
    m_livestemc_to_litter          => pcisof%m_livestemc_to_litter
    m_livestemc_xfer_to_litter     => pcisof%m_livestemc_xfer_to_litter
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

    ! column level carbon fluxes from gap-phase mortality
    do j = 1 , nlevdecomp
      ! column loop
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        ! column gap mortality fluxes
        decomp_cpools_vr(c,j,i_met_lit) = decomp_cpools_vr(c,j,i_met_lit) + &
          gap_mortality_c_to_litr_met_c(c,j) * dt
        decomp_cpools_vr(c,j,i_cel_lit) = decomp_cpools_vr(c,j,i_cel_lit) + &
          gap_mortality_c_to_litr_cel_c(c,j) * dt
        decomp_cpools_vr(c,j,i_lig_lit) = decomp_cpools_vr(c,j,i_lig_lit) + &
          gap_mortality_c_to_litr_lig_c(c,j) * dt
        decomp_cpools_vr(c,j,i_cwd) = decomp_cpools_vr(c,j,i_cwd) + &
          gap_mortality_c_to_cwdc(c,j) * dt
      end do
    end do

    ! pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      ! pft-level carbon fluxes from gap-phase mortality
      ! displayed pools
      leafc(p) = leafc(p) - m_leafc_to_litter(p) * dt
      frootc(p) = frootc(p) - m_frootc_to_litter(p) * dt
      livestemc(p) = livestemc(p) - m_livestemc_to_litter(p) * dt
      deadstemc(p) = deadstemc(p) - m_deadstemc_to_litter(p) * dt
      livecrootc(p) = livecrootc(p) - m_livecrootc_to_litter(p) * dt
      deadcrootc(p) = deadcrootc(p) - m_deadcrootc_to_litter(p) * dt

      ! storage pools
      leafc_storage(p) = leafc_storage(p) - m_leafc_storage_to_litter(p) * dt
      frootc_storage(p) = frootc_storage(p) - &
        m_frootc_storage_to_litter(p) * dt
      livestemc_storage(p) = livestemc_storage(p) - &
        m_livestemc_storage_to_litter(p) * dt
      deadstemc_storage(p) = deadstemc_storage(p) - &
        m_deadstemc_storage_to_litter(p) * dt
      livecrootc_storage(p) = livecrootc_storage(p) - &
        m_livecrootc_storage_to_litter(p) * dt
      deadcrootc_storage(p) = deadcrootc_storage(p) - &
        m_deadcrootc_storage_to_litter(p) * dt
      gresp_storage(p) = gresp_storage(p) - m_gresp_storage_to_litter(p) * dt

      ! transfer pools
      leafc_xfer(p) = leafc_xfer(p) - m_leafc_xfer_to_litter(p) * dt
      frootc_xfer(p) = frootc_xfer(p) - m_frootc_xfer_to_litter(p) * dt
      livestemc_xfer(p) = livestemc_xfer(p) - &
        m_livestemc_xfer_to_litter(p) * dt
      deadstemc_xfer(p) = deadstemc_xfer(p) - &
        m_deadstemc_xfer_to_litter(p) * dt
      livecrootc_xfer(p) = livecrootc_xfer(p) - &
        m_livecrootc_xfer_to_litter(p) * dt
      deadcrootc_xfer(p) = deadcrootc_xfer(p) - &
        m_deadcrootc_xfer_to_litter(p) * dt
      gresp_xfer(p) = gresp_xfer(p) - m_gresp_xfer_to_litter(p) * dt
    end do ! end of pft loop
  end subroutine CStateUpdate2
  !
  ! Update all the prognostic carbon state
  ! variables affected by harvest mortality fluxes
  !
  subroutine CStateUpdate2h(num_soilc, filter_soilc, &
                            num_soilp, filter_soilp, isotope)
    use mod_clm_type
    use mod_clm_varpar , only : nlevdecomp
    use mod_clm_varpar , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    implicit none
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts
    character(len=*), intent(in) :: isotope     ! 'bulk', 'c13' or 'c14'

    ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
    real(rk8), pointer :: harvest_c_to_litr_met_c(:,:)
    ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
    real(rk8), pointer :: harvest_c_to_litr_cel_c(:,:)
    ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
    real(rk8), pointer :: harvest_c_to_litr_lig_c(:,:)
    ! C fluxes associated with harvest to CWD pool (gC/m3/s)
    real(rk8), pointer :: harvest_c_to_cwdc(:,:)
    real(rk8), pointer :: hrv_deadcrootc_storage_to_litter(:)
    real(rk8), pointer :: hrv_deadcrootc_to_litter(:)
    real(rk8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_deadstemc_storage_to_litter(:)
    real(rk8), pointer :: hrv_deadstemc_to_prod10c(:)
    real(rk8), pointer :: hrv_deadstemc_to_prod100c(:)
    real(rk8), pointer :: hrv_deadstemc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_frootc_storage_to_litter(:)
    real(rk8), pointer :: hrv_frootc_to_litter(:)
    real(rk8), pointer :: hrv_frootc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_gresp_storage_to_litter(:)
    real(rk8), pointer :: hrv_gresp_xfer_to_litter(:)
    real(rk8), pointer :: hrv_leafc_storage_to_litter(:)
    real(rk8), pointer :: hrv_leafc_to_litter(:)
    real(rk8), pointer :: hrv_leafc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_livecrootc_storage_to_litter(:)
    real(rk8), pointer :: hrv_livecrootc_to_litter(:)
    real(rk8), pointer :: hrv_livecrootc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_livestemc_storage_to_litter(:)
    real(rk8), pointer :: hrv_livestemc_to_litter(:)
    real(rk8), pointer :: hrv_livestemc_xfer_to_litter(:)
    real(rk8), pointer :: hrv_xsmrpool_to_atm(:)

    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer :: decomp_cpools_vr(:,:,:)
    ! (gC/m2) dead coarse root C
    real(rk8), pointer :: deadcrootc(:)
    ! (gC/m2) dead coarse root C storage
    real(rk8), pointer :: deadcrootc_storage(:)
    ! (gC/m2) dead coarse root C transfer
    real(rk8), pointer :: deadcrootc_xfer(:)
    real(rk8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
    real(rk8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
    real(rk8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rk8), pointer :: frootc(:)             ! (gC/m2) fine root C
    real(rk8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
    real(rk8), pointer :: frootc_xfer(:)   ! (gC/m2) fine root C transfer
    real(rk8), pointer :: gresp_storage(:) ! (gC/m2) growth respiration storage
    real(rk8), pointer :: gresp_xfer(:)    ! (gC/m2) growth respiration transfer
    real(rk8), pointer :: leafc(:)         ! (gC/m2) leaf C
    real(rk8), pointer :: leafc_storage(:) ! (gC/m2) leaf C storage
    real(rk8), pointer :: leafc_xfer(:)    ! (gC/m2) leaf C transfer
    real(rk8), pointer :: livecrootc(:)    ! (gC/m2) live coarse root C
    ! (gC/m2) live coarse root C storage
    real(rk8), pointer :: livecrootc_storage(:)
    ! (gC/m2) live coarse root C transfer
    real(rk8), pointer :: livecrootc_xfer(:)
    real(rk8), pointer :: livestemc(:)         ! (gC/m2) live stem C
    real(rk8), pointer :: livestemc_storage(:) ! (gC/m2) live stem C storage
    real(rk8), pointer :: livestemc_xfer(:)    ! (gC/m2) live stem C transfer
    ! (gC/m2) abstract C pool to meet excess MR demand
    real(rk8), pointer :: xsmrpool(:)

    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    type(column_cflux_type), pointer :: ccisof
    type(column_cstate_type), pointer :: ccisos
    integer(ik4) :: c,p,j    ! indices
    integer(ik4) :: fp,fc    ! lake filter indices
    real(rk8):: dt       ! radiation time step (seconds)

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
        call fatal(__FILE__,__LINE__, &
          'CNCIsoStateUpdate2Mod: iso must be bulk, c13 or c14')
    end select

    ! assign local pointers at the column level
    harvest_c_to_litr_met_c          => ccisof%harvest_c_to_litr_met_c
    harvest_c_to_litr_cel_c          => ccisof%harvest_c_to_litr_cel_c
    harvest_c_to_litr_lig_c          => ccisof%harvest_c_to_litr_lig_c
    harvest_c_to_cwdc                => ccisof%harvest_c_to_cwdc
    decomp_cpools_vr                 => ccisos%decomp_cpools_vr

    ! assign local pointers at the pft level
    hrv_deadcrootc_storage_to_litter => pcisof%hrv_deadcrootc_storage_to_litter
    hrv_deadcrootc_to_litter         => pcisof%hrv_deadcrootc_to_litter
    hrv_deadcrootc_xfer_to_litter    => pcisof%hrv_deadcrootc_xfer_to_litter
    hrv_deadstemc_storage_to_litter  => pcisof%hrv_deadstemc_storage_to_litter
    hrv_deadstemc_to_prod10c         => pcisof%hrv_deadstemc_to_prod10c
    hrv_deadstemc_to_prod100c        => pcisof%hrv_deadstemc_to_prod100c
    hrv_deadstemc_xfer_to_litter     => pcisof%hrv_deadstemc_xfer_to_litter
    hrv_frootc_storage_to_litter     => pcisof%hrv_frootc_storage_to_litter
    hrv_frootc_to_litter             => pcisof%hrv_frootc_to_litter
    hrv_frootc_xfer_to_litter        => pcisof%hrv_frootc_xfer_to_litter
    hrv_gresp_storage_to_litter      => pcisof%hrv_gresp_storage_to_litter
    hrv_gresp_xfer_to_litter         => pcisof%hrv_gresp_xfer_to_litter
    hrv_leafc_storage_to_litter      => pcisof%hrv_leafc_storage_to_litter
    hrv_leafc_to_litter              => pcisof%hrv_leafc_to_litter
    hrv_leafc_xfer_to_litter         => pcisof%hrv_leafc_xfer_to_litter
    hrv_livecrootc_storage_to_litter => pcisof%hrv_livecrootc_storage_to_litter
    hrv_livecrootc_to_litter         => pcisof%hrv_livecrootc_to_litter
    hrv_livecrootc_xfer_to_litter    => pcisof%hrv_livecrootc_xfer_to_litter
    hrv_livestemc_storage_to_litter  => pcisof%hrv_livestemc_storage_to_litter
    hrv_livestemc_to_litter          => pcisof%hrv_livestemc_to_litter
    hrv_livestemc_xfer_to_litter     => pcisof%hrv_livestemc_xfer_to_litter
    hrv_xsmrpool_to_atm              => pcisof%hrv_xsmrpool_to_atm
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
    xsmrpool                       => pcisos%xsmrpool

    ! set time steps
    dt = dtsrf

    ! column level carbon fluxes from harvest mortality
    do j = 1 , nlevdecomp
      ! column loop
      do fc = 1 , num_soilc
        c = filter_soilc(fc)

        ! column harvest fluxes
        decomp_cpools_vr(c,j,i_met_lit) = decomp_cpools_vr(c,j,i_met_lit) + &
          harvest_c_to_litr_met_c(c,j) * dt
        decomp_cpools_vr(c,j,i_cel_lit) = decomp_cpools_vr(c,j,i_cel_lit) + &
          harvest_c_to_litr_cel_c(c,j) * dt
        decomp_cpools_vr(c,j,i_lig_lit) = decomp_cpools_vr(c,j,i_lig_lit) + &
          harvest_c_to_litr_lig_c(c,j) * dt
        decomp_cpools_vr(c,j,i_cwd) = decomp_cpools_vr(c,j,i_cwd) + &
          harvest_c_to_cwdc(c,j)  * dt

        ! wood to product pools - states updated in CNWoodProducts()
      end do
    end do

    ! pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      ! pft-level carbon fluxes from harvest mortality
      ! displayed pools
      leafc(p) = leafc(p) - hrv_leafc_to_litter(p) * dt
      frootc(p) = frootc(p) - hrv_frootc_to_litter(p) * dt
      livestemc(p) = livestemc(p) - hrv_livestemc_to_litter(p) * dt
      deadstemc(p) = deadstemc(p) - hrv_deadstemc_to_prod10c(p) * dt
      deadstemc(p) = deadstemc(p) - hrv_deadstemc_to_prod100c(p) * dt
      livecrootc(p) = livecrootc(p) - hrv_livecrootc_to_litter(p) * dt
      deadcrootc(p) = deadcrootc(p) - hrv_deadcrootc_to_litter(p) * dt

      ! xsmrpool
      xsmrpool(p) = xsmrpool(p) - hrv_xsmrpool_to_atm(p) * dt

      ! storage pools
      leafc_storage(p) = leafc_storage(p) - &
        hrv_leafc_storage_to_litter(p) * dt
      frootc_storage(p) = frootc_storage(p) - &
        hrv_frootc_storage_to_litter(p) * dt
      livestemc_storage(p) = livestemc_storage(p) - &
        hrv_livestemc_storage_to_litter(p)  * dt
      deadstemc_storage(p) = deadstemc_storage(p) - &
        hrv_deadstemc_storage_to_litter(p)  * dt
      livecrootc_storage(p) = livecrootc_storage(p) - &
        hrv_livecrootc_storage_to_litter(p) * dt
      deadcrootc_storage(p) = deadcrootc_storage(p) - &
        hrv_deadcrootc_storage_to_litter(p) * dt
      gresp_storage(p) = gresp_storage(p) - &
        hrv_gresp_storage_to_litter(p) * dt

      ! transfer pools
      leafc_xfer(p) = leafc_xfer(p) - hrv_leafc_xfer_to_litter(p) * dt
      frootc_xfer(p) = frootc_xfer(p) - hrv_frootc_xfer_to_litter(p) * dt
      livestemc_xfer(p) = livestemc_xfer(p) - &
        hrv_livestemc_xfer_to_litter(p) * dt
      deadstemc_xfer(p) = deadstemc_xfer(p) - &
        hrv_deadstemc_xfer_to_litter(p) * dt
      livecrootc_xfer(p) = livecrootc_xfer(p) - &
        hrv_livecrootc_xfer_to_litter(p) * dt
      deadcrootc_xfer(p) = deadcrootc_xfer(p) - &
        hrv_deadcrootc_xfer_to_litter(p) * dt
      gresp_xfer(p) = gresp_xfer(p) - hrv_gresp_xfer_to_litter(p) * dt

    end do ! end of pft loop
  end subroutine CStateUpdate2h

#endif

end module mod_clm_cncstateupdate2
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
