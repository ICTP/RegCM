module mod_clm_cngresp
#ifdef CN
  !
  ! Module for growth respiration fluxes,
  ! for coupled carbon-nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds

  implicit none

  save

  private

  public :: CNGResp

  contains
  !
  ! On the radiation time step, update all the prognostic carbon state
  ! variables
  !
  subroutine CNGResp(num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_pftvarcon , only : npcropmin, grperc, grpnow
    implicit none
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    integer(ik4) , pointer :: ivt(:)  ! pft vegetation type
    real(rk8), pointer :: cpool_to_leafc(:)
    real(rk8), pointer :: cpool_to_leafc_storage(:)
    real(rk8), pointer :: cpool_to_frootc(:)
    real(rk8), pointer :: cpool_to_frootc_storage(:)
    real(rk8), pointer :: cpool_to_livestemc(:)
    real(rk8), pointer :: cpool_to_livestemc_storage(:)
    real(rk8), pointer :: cpool_to_deadstemc(:)
    real(rk8), pointer :: cpool_to_deadstemc_storage(:)
    real(rk8), pointer :: cpool_to_livecrootc(:)
    real(rk8), pointer :: cpool_to_livecrootc_storage(:)
    ! allocation to dead coarse root C (gC/m2/s)
    real(rk8), pointer :: cpool_to_deadcrootc(:)
    ! allocation to dead coarse root C storage (gC/m2/s)
    real(rk8), pointer :: cpool_to_deadcrootc_storage(:)
    ! allocation to grain C (gC/m2/s)
    real(rk8), pointer :: cpool_to_grainc(:)
    ! allocation to grain C storage (gC/m2/s)
    real(rk8), pointer :: cpool_to_grainc_storage(:)
    ! grain C growth from storage (gC/m2/s)
    real(rk8), pointer :: grainc_xfer_to_grainc(:)
    ! leaf C growth from storage (gC/m2/s)
    real(rk8), pointer :: leafc_xfer_to_leafc(:)
    ! fine root C growth from storage (gC/m2/s)
    real(rk8), pointer :: frootc_xfer_to_frootc(:)
    ! live stem C growth from storage (gC/m2/s)
    real(rk8), pointer :: livestemc_xfer_to_livestemc(:)
    ! dead stem C growth from storage (gC/m2/s)
    real(rk8), pointer :: deadstemc_xfer_to_deadstemc(:)
    ! live coarse root C growth from storage (gC/m2/s)
    real(rk8), pointer :: livecrootc_xfer_to_livecrootc(:)
    ! dead coarse root C growth from storage (gC/m2/s)
    real(rk8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
    ! binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8), pointer :: woody(:)

    real(rk8), pointer :: cpool_grain_gr(:)
    real(rk8), pointer :: cpool_grain_storage_gr(:)
    real(rk8), pointer :: transfer_grain_gr(:)
    real(rk8), pointer :: cpool_leaf_gr(:)
    real(rk8), pointer :: cpool_leaf_storage_gr(:)
    real(rk8), pointer :: transfer_leaf_gr(:)
    real(rk8), pointer :: cpool_froot_gr(:)
    real(rk8), pointer :: cpool_froot_storage_gr(:)
    real(rk8), pointer :: transfer_froot_gr(:)
    real(rk8), pointer :: cpool_livestem_gr(:)
    real(rk8), pointer :: cpool_livestem_storage_gr(:)
    real(rk8), pointer :: transfer_livestem_gr(:)
    real(rk8), pointer :: cpool_deadstem_gr(:)
    real(rk8), pointer :: cpool_deadstem_storage_gr(:)
    real(rk8), pointer :: transfer_deadstem_gr(:)
    real(rk8), pointer :: cpool_livecroot_gr(:)
    real(rk8), pointer :: cpool_livecroot_storage_gr(:)
    real(rk8), pointer :: transfer_livecroot_gr(:)
    real(rk8), pointer :: cpool_deadcroot_gr(:)
    real(rk8), pointer :: cpool_deadcroot_storage_gr(:)
    real(rk8), pointer :: transfer_deadcroot_gr(:)

    integer(ik4) :: p      ! indices
    integer(ik4) :: fp     ! lake filter pft index

    ! Assign local pointers to derived type arrays (in)
    ivt                         => clm3%g%l%c%p%itype
    cpool_to_leafc              => clm3%g%l%c%p%pcf%cpool_to_leafc
    cpool_to_leafc_storage      => clm3%g%l%c%p%pcf%cpool_to_leafc_storage
    cpool_to_frootc             => clm3%g%l%c%p%pcf%cpool_to_frootc
    cpool_to_frootc_storage     => clm3%g%l%c%p%pcf%cpool_to_frootc_storage
    cpool_to_livestemc          => clm3%g%l%c%p%pcf%cpool_to_livestemc
    cpool_to_livestemc_storage  => clm3%g%l%c%p%pcf%cpool_to_livestemc_storage
    cpool_to_deadstemc          => clm3%g%l%c%p%pcf%cpool_to_deadstemc
    cpool_to_deadstemc_storage  => clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage
    cpool_to_livecrootc         => clm3%g%l%c%p%pcf%cpool_to_livecrootc
    cpool_to_livecrootc_storage => clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage
    cpool_to_deadcrootc         => clm3%g%l%c%p%pcf%cpool_to_deadcrootc
    cpool_to_deadcrootc_storage => clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage
    cpool_to_grainc             => clm3%g%l%c%p%pcf%cpool_to_grainc
    cpool_to_grainc_storage     => clm3%g%l%c%p%pcf%cpool_to_grainc_storage
    grainc_xfer_to_grainc       => clm3%g%l%c%p%pcf%grainc_xfer_to_grainc
    leafc_xfer_to_leafc         => clm3%g%l%c%p%pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc       => clm3%g%l%c%p%pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc => clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc => clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc => &
            clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc => &
            clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc
    woody => pftcon%woody

    ! Assign local pointers to derived type arrays (out)
    cpool_grain_gr                => clm3%g%l%c%p%pcf%cpool_grain_gr
    cpool_grain_storage_gr        => clm3%g%l%c%p%pcf%cpool_grain_storage_gr
    transfer_grain_gr             => clm3%g%l%c%p%pcf%transfer_grain_gr
    cpool_leaf_gr                 => clm3%g%l%c%p%pcf%cpool_leaf_gr
    cpool_leaf_storage_gr         => clm3%g%l%c%p%pcf%cpool_leaf_storage_gr
    transfer_leaf_gr              => clm3%g%l%c%p%pcf%transfer_leaf_gr
    cpool_froot_gr                => clm3%g%l%c%p%pcf%cpool_froot_gr
    cpool_froot_storage_gr        => clm3%g%l%c%p%pcf%cpool_froot_storage_gr
    transfer_froot_gr             => clm3%g%l%c%p%pcf%transfer_froot_gr
    cpool_livestem_gr             => clm3%g%l%c%p%pcf%cpool_livestem_gr
    cpool_livestem_storage_gr     => clm3%g%l%c%p%pcf%cpool_livestem_storage_gr
    transfer_livestem_gr          => clm3%g%l%c%p%pcf%transfer_livestem_gr
    cpool_deadstem_gr             => clm3%g%l%c%p%pcf%cpool_deadstem_gr
    cpool_deadstem_storage_gr     => clm3%g%l%c%p%pcf%cpool_deadstem_storage_gr
    transfer_deadstem_gr          => clm3%g%l%c%p%pcf%transfer_deadstem_gr
    cpool_livecroot_gr            => clm3%g%l%c%p%pcf%cpool_livecroot_gr
    cpool_livecroot_storage_gr    => clm3%g%l%c%p%pcf%cpool_livecroot_storage_gr
    transfer_livecroot_gr         => clm3%g%l%c%p%pcf%transfer_livecroot_gr
    cpool_deadcroot_gr            => clm3%g%l%c%p%pcf%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr    => clm3%g%l%c%p%pcf%cpool_deadcroot_storage_gr
    transfer_deadcroot_gr         => clm3%g%l%c%p%pcf%transfer_deadcroot_gr

    ! Loop through pfts
    ! start pft loop
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      if ( ivt(p) >= npcropmin ) then ! skip 2 generic crops
        cpool_livestem_gr(p)         = cpool_to_livestemc(p) * grperc(ivt(p))
        cpool_livestem_storage_gr(p) = cpool_to_livestemc_storage(p) * &
                                       grperc(ivt(p)) * grpnow(ivt(p))
        transfer_livestem_gr(p)      = livestemc_xfer_to_livestemc(p) * &
                                       grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
        cpool_grain_gr(p)            = cpool_to_grainc(p) * grperc(ivt(p))
        cpool_grain_storage_gr(p)    = cpool_to_grainc_storage(p) * &
                                       grperc(ivt(p)) * grpnow(ivt(p))
        transfer_grain_gr(p)         = grainc_xfer_to_grainc(p) * &
                                       grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
      end if

      ! leaf and fine root growth respiration
      cpool_leaf_gr(p)          = cpool_to_leafc(p) * grperc(ivt(p))
      cpool_leaf_storage_gr(p)  = cpool_to_leafc_storage(p) * grperc(ivt(p)) * &
                                  grpnow(ivt(p))
      transfer_leaf_gr(p)       = leafc_xfer_to_leafc(p) * grperc(ivt(p)) * &
                                  (1.D0 - grpnow(ivt(p)))
      cpool_froot_gr(p)         = cpool_to_frootc(p) * grperc(ivt(p))
      cpool_froot_storage_gr(p) = cpool_to_frootc_storage(p) * &
                                  grperc(ivt(p)) * grpnow(ivt(p))
      transfer_froot_gr(p)      = frootc_xfer_to_frootc(p) * grperc(ivt(p)) * &
                                  (1.D0 - grpnow(ivt(p)))

      if ( woody(ivt(p)) == 1.D0 ) then
        cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc(ivt(p))
        cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * &
                                        grperc(ivt(p)) * grpnow(ivt(p))
        transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * &
                                        grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
        cpool_deadstem_gr(p)          = cpool_to_deadstemc(p) * grperc(ivt(p))
        cpool_deadstem_storage_gr(p)  = cpool_to_deadstemc_storage(p) * &
                                        grperc(ivt(p)) * grpnow(ivt(p))
        transfer_deadstem_gr(p)       = deadstemc_xfer_to_deadstemc(p) * &
                                        grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
        cpool_livecroot_gr(p)         = cpool_to_livecrootc(p) * grperc(ivt(p))
        cpool_livecroot_storage_gr(p) = cpool_to_livecrootc_storage(p) * &
                                        grperc(ivt(p)) * grpnow(ivt(p))
        transfer_livecroot_gr(p)      = livecrootc_xfer_to_livecrootc(p) * &
                                        grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
        cpool_deadcroot_gr(p)         = cpool_to_deadcrootc(p) * grperc(ivt(p))
        cpool_deadcroot_storage_gr(p) = cpool_to_deadcrootc_storage(p) * &
                                        grperc(ivt(p)) * grpnow(ivt(p))
        transfer_deadcroot_gr(p)      = deadcrootc_xfer_to_deadcrootc(p) * &
                                        grperc(ivt(p)) * (1.D0 - grpnow(ivt(p)))
      end if
    end do
  end subroutine CNGResp

#endif

end module mod_clm_cngresp
