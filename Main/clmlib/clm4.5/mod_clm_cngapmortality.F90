module mod_clm_cngapmortality
#ifdef CN
  !
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam, only : dayspy

  implicit none

  save

  private

  public :: CNGapMortality

  contains
  !
  ! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
  !
  subroutine CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_varcon   , only : secspday
    use mod_clm_pftvarcon, only : npcropmin
    implicit none
    integer(ik4), intent(in) :: num_soilc  ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! column filter for soil points
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! pft filter for soil points

    integer(ik4), pointer, contiguous :: ivt(:) ! pft vegetation type
    real(rk8), pointer, contiguous :: woody(:)   ! binary flag for woody lifeform
                                     ! (1=woody, 0=not woody)
    real(rk8), pointer, contiguous :: leafc(:)              ! (gC/m2) leaf C
    real(rk8), pointer, contiguous :: frootc(:)             ! (gC/m2) fine root C
    real(rk8), pointer, contiguous :: livestemc(:)          ! (gC/m2) live stem C
    real(rk8), pointer, contiguous :: deadstemc(:)          ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous :: livecrootc(:)         ! (gC/m2) live coarse root C
    real(rk8), pointer, contiguous :: deadcrootc(:)         ! (gC/m2) dead coarse root C
    real(rk8), pointer, contiguous :: leafc_storage(:)      ! (gC/m2) leaf C storage
    real(rk8), pointer, contiguous :: frootc_storage(:)     ! (gC/m2) fine root C storage
    real(rk8), pointer, contiguous :: livestemc_storage(:)  ! (gC/m2) live stem C storage
    real(rk8), pointer, contiguous :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
    real(rk8), pointer, contiguous :: livecrootc_storage(:) ! (gC/m2) live corse root C stge
    real(rk8), pointer, contiguous :: deadcrootc_storage(:) ! (gC/m2) dead corse root C stge
    real(rk8), pointer, contiguous :: gresp_storage(:)      ! (gC/m2) growth respir. stge
    real(rk8), pointer, contiguous :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
    real(rk8), pointer, contiguous :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
    real(rk8), pointer, contiguous :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
    real(rk8), pointer, contiguous :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rk8), pointer, contiguous :: livecrootc_xfer(:)    ! (gC/m2) live corse root C trfr
    real(rk8), pointer, contiguous :: deadcrootc_xfer(:)    ! (gC/m2) dead corse root C trfr
    real(rk8), pointer, contiguous :: gresp_xfer(:)         ! (gC/m2) growth resp trfr
    real(rk8), pointer, contiguous :: leafn(:)              ! (gN/m2) leaf N
    real(rk8), pointer, contiguous :: frootn(:)             ! (gN/m2) fine root N
    real(rk8), pointer, contiguous :: livestemn(:)          ! (gN/m2) live stem N
    real(rk8), pointer, contiguous :: deadstemn(:)          ! (gN/m2) dead stem N
    real(rk8), pointer, contiguous :: livecrootn(:)         ! (gN/m2) live coarse root N
    real(rk8), pointer, contiguous :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rk8), pointer, contiguous :: retransn(:)           ! (gN/m2) plant pool of
                                                ! retranslocated N
    real(rk8), pointer, contiguous :: leafn_storage(:)      ! (gN/m2) leaf N storage
    real(rk8), pointer, contiguous :: frootn_storage(:)     ! (gN/m2) fine root N storage
    real(rk8), pointer, contiguous :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rk8), pointer, contiguous :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
    real(rk8), pointer, contiguous :: livecrootn_storage(:) ! (gN/m2) live corse root N stge
    real(rk8), pointer, contiguous :: deadcrootn_storage(:) ! (gN/m2) dead corse root N stge
    real(rk8), pointer, contiguous :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
    real(rk8), pointer, contiguous :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
    real(rk8), pointer, contiguous :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rk8), pointer, contiguous :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
    real(rk8), pointer, contiguous :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N trfr
    real(rk8), pointer, contiguous :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N trfr
#if (defined CNDV)
    real(rk8), pointer, contiguous :: greffic(:)
    real(rk8), pointer, contiguous :: heatstress(:)
    ! number of individuals (#/m2) added by F. Li and S. Levis
    real(rk8), pointer, contiguous :: nind(:)
#endif

    real(rk8), pointer, contiguous :: m_leafc_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_gresp_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_gresp_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_retransn_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_xfer_to_litter(:)

    integer(ik4) :: p  ! pft index
    integer(ik4) :: fp ! pft filter index
    real(rk8):: am     ! rate for fractional mortality (1/yr)
    real(rk8):: m      ! rate for fractional mortality (1/s)
#ifdef CNDV
    ! asymptotic max mortality rate (/yr)
    real(rk8):: mort_max
    !coeff of growth efficiency in mortality equation
    real(rk8), parameter :: k_mort = 0.3_rk8
    logical, dimension(num_soilp) :: iswoody
#endif

    ! assign local pointers
    woody                          => pftcon%woody

    ! assign local pointers to pft-level arrays
    ivt                        => clm3%g%l%c%p%itype
    leafc                      => clm3%g%l%c%p%pcs%leafc
    frootc                     => clm3%g%l%c%p%pcs%frootc
    livestemc                  => clm3%g%l%c%p%pcs%livestemc
    deadstemc                  => clm3%g%l%c%p%pcs%deadstemc
    livecrootc                 => clm3%g%l%c%p%pcs%livecrootc
    deadcrootc                 => clm3%g%l%c%p%pcs%deadcrootc
    leafc_storage              => clm3%g%l%c%p%pcs%leafc_storage
    frootc_storage             => clm3%g%l%c%p%pcs%frootc_storage
    livestemc_storage          => clm3%g%l%c%p%pcs%livestemc_storage
    deadstemc_storage          => clm3%g%l%c%p%pcs%deadstemc_storage
    livecrootc_storage         => clm3%g%l%c%p%pcs%livecrootc_storage
    deadcrootc_storage         => clm3%g%l%c%p%pcs%deadcrootc_storage
    gresp_storage              => clm3%g%l%c%p%pcs%gresp_storage
    leafc_xfer                 => clm3%g%l%c%p%pcs%leafc_xfer
    frootc_xfer                => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc_xfer             => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc_xfer             => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc_xfer            => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc_xfer            => clm3%g%l%c%p%pcs%deadcrootc_xfer
    gresp_xfer                 => clm3%g%l%c%p%pcs%gresp_xfer
    leafn                      => clm3%g%l%c%p%pns%leafn
    frootn                     => clm3%g%l%c%p%pns%frootn
    livestemn                  => clm3%g%l%c%p%pns%livestemn
    deadstemn                  => clm3%g%l%c%p%pns%deadstemn
    livecrootn                 => clm3%g%l%c%p%pns%livecrootn
    deadcrootn                 => clm3%g%l%c%p%pns%deadcrootn
    retransn                   => clm3%g%l%c%p%pns%retransn
    leafn_storage              => clm3%g%l%c%p%pns%leafn_storage
    frootn_storage             => clm3%g%l%c%p%pns%frootn_storage
    livestemn_storage          => clm3%g%l%c%p%pns%livestemn_storage
    deadstemn_storage          => clm3%g%l%c%p%pns%deadstemn_storage
    livecrootn_storage         => clm3%g%l%c%p%pns%livecrootn_storage
    deadcrootn_storage         => clm3%g%l%c%p%pns%deadcrootn_storage
    leafn_xfer                 => clm3%g%l%c%p%pns%leafn_xfer
    frootn_xfer                => clm3%g%l%c%p%pns%frootn_xfer
    livestemn_xfer             => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn_xfer             => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn_xfer            => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn_xfer            => clm3%g%l%c%p%pns%deadcrootn_xfer
    m_leafc_to_litter          => clm3%g%l%c%p%pcf%m_leafc_to_litter
    m_frootc_to_litter         => clm3%g%l%c%p%pcf%m_frootc_to_litter
    m_livestemc_to_litter      => clm3%g%l%c%p%pcf%m_livestemc_to_litter
    m_deadstemc_to_litter      => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
    m_livecrootc_to_litter     => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
    m_deadcrootc_to_litter     => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
    m_leafc_storage_to_litter  => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
    m_frootc_storage_to_litter => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
    m_livestemc_storage_to_litter  => &
            clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
    m_deadstemc_storage_to_litter  => &
            clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
    m_livecrootc_storage_to_litter => &
            clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
    m_deadcrootc_storage_to_litter => &
            clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
    m_gresp_storage_to_litter  => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
    m_leafc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
    m_frootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
    m_livestemc_xfer_to_litter => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
    m_deadstemc_xfer_to_litter => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
    m_livecrootc_xfer_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
    m_deadcrootc_xfer_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
    m_gresp_xfer_to_litter      => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
    m_leafn_to_litter           => clm3%g%l%c%p%pnf%m_leafn_to_litter
    m_frootn_to_litter          => clm3%g%l%c%p%pnf%m_frootn_to_litter
    m_livestemn_to_litter       => clm3%g%l%c%p%pnf%m_livestemn_to_litter
    m_deadstemn_to_litter       => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
    m_livecrootn_to_litter      => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
    m_deadcrootn_to_litter      => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
    m_retransn_to_litter        => clm3%g%l%c%p%pnf%m_retransn_to_litter
    m_leafn_storage_to_litter   => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
    m_frootn_storage_to_litter  => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
    m_livestemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
    m_deadstemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
    m_livecrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
    m_deadcrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
    m_leafn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
    m_frootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
    m_livestemn_xfer_to_litter => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
    m_deadstemn_xfer_to_litter => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
    m_livecrootn_xfer_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
    m_deadcrootn_xfer_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
#if (defined CNDV)
    greffic     => clm3%g%l%c%p%pdgvs%greffic
    heatstress  => clm3%g%l%c%p%pdgvs%heatstress
    nind        => clm3%g%l%c%p%pdgvs%nind     ! F. Li and S. Levis
#endif

    ! set the mortality rate based on annual rate
    am = 0.02_rk8

#if (defined CNDV)
    do fp = 1, num_soilp
      p = filter_soilp(fp)
      if ( abs(woody(ivt(p))-1._rk8) < epsilon(1.0) ) then
        iswoody(fp) = .true.
      else
        iswoody(fp) = .false.
      end if
    end do
#endif

    ! pft loop
    do fp = 1, num_soilp
      p = filter_soilp(fp)

#if (defined CNDV)
      ! Stress mortality from lpj's subr Mortality.
      if ( iswoody(fp) ) then
        if ( ivt(p) == 8 ) then
          mort_max = 0.03_rk8 ! BDT boreal
        else
          mort_max = 0.01_rk8 ! original value for all pfts
        end if
        ! heatstress and greffic calculated in Establishment once/yr
        ! Mortality rate inversely related to growth efficiency
        ! (Prentice et al 1993)
        am = mort_max / (1.0_rk8 + k_mort * greffic(p))
        am = min(1.0_rk8, am + heatstress(p))
      else ! lpj didn't set this for grasses; cn does
        ! set the mortality rate based on annual rate
        am = 0.02_rk8
      end if
#endif
      m = am/(dayspy * secspday)

      ! pft-level gap mortality carbon fluxes
      ! displayed pools
      m_leafc_to_litter(p)               = leafc(p)               * m
      m_frootc_to_litter(p)              = frootc(p)              * m
      m_livestemc_to_litter(p)           = livestemc(p)           * m
      m_deadstemc_to_litter(p)           = deadstemc(p)           * m
      m_livecrootc_to_litter(p)          = livecrootc(p)          * m
      m_deadcrootc_to_litter(p)          = deadcrootc(p)          * m

      ! storage pools
      m_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
      m_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
      m_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
      m_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
      m_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
      m_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
      m_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

      ! transfer pools
      m_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
      m_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
      m_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
      m_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
      m_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
      m_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
      m_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

      ! pft-level gap mortality nitrogen fluxes
      ! displayed pools
      m_leafn_to_litter(p)               = leafn(p)               * m
      m_frootn_to_litter(p)              = frootn(p)              * m
      m_livestemn_to_litter(p)           = livestemn(p)           * m
      m_deadstemn_to_litter(p)           = deadstemn(p)           * m
      m_livecrootn_to_litter(p)          = livecrootn(p)          * m
      m_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
      if (ivt(p) < npcropmin) m_retransn_to_litter(p) = retransn(p) * m

      ! storage pools
      m_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
      m_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
      m_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
      m_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
      m_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
      m_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

      ! transfer pools
      m_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
      m_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
      m_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
      m_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
      m_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
      m_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m

#if (defined CNDV)
      ! added by F. Li and S. Levis
      if ( iswoody(fp) ) then
        if ( livestemc(p)+deadstemc(p) > 0.0_rk8 ) then
          nind(p) = nind(p) * (1.0_rk8-m)
        else
          nind(p) = 0.0_rk8
        end if
      end if
#endif
    end do ! end of pft loop

    ! gather all pft-level litterfall fluxes to the column
    ! for litter C and N inputs

    call CNGapPftToColumn(num_soilc, filter_soilc)

  end subroutine CNGapMortality
  !
  ! called in the middle of CNGapMoratlity to gather all pft-level
  ! gap mortality fluxes to the column level and assign them to the
  ! three litter pools
  !
  subroutine CNGapPftToColumn (num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varpar, only : maxpatch_pft, nlevdecomp
    implicit none
    integer(ik4), intent(in) :: num_soilc  ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! soil column filter

    ! true=>do computations on this pft (see reweightMod for details)
    logical, pointer, contiguous :: pactive(:)
    integer(ik4), pointer, contiguous :: ivt(:)      ! pft vegetation type
    real(rk8), pointer, contiguous :: wtcol(:)    ! pft weight relative to column (0-1)
    real(rk8), pointer, contiguous :: lf_flab(:)  ! leaf litter labile fraction
    real(rk8), pointer, contiguous :: lf_fcel(:)  ! leaf litter cellulose fraction
    real(rk8), pointer, contiguous :: lf_flig(:)  ! leaf litter lignin fraction
    real(rk8), pointer, contiguous :: fr_flab(:)  ! fine root litter labile fraction
    real(rk8), pointer, contiguous :: fr_fcel(:)  ! fine root litter cellulose fraction
    real(rk8), pointer, contiguous :: fr_flig(:)  ! fine root litter lignin fraction
    integer(ik4), pointer, contiguous :: npfts(:)    ! number of pfts for each column
    integer(ik4), pointer, contiguous :: pfti(:)     ! beginning pft index for each column
    real(rk8), pointer, contiguous :: m_leafc_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_gresp_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootc_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_gresp_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_to_litter(:)
    real(rk8), pointer, contiguous :: m_retransn_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_storage_to_litter(:)
    real(rk8), pointer, contiguous :: m_leafn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_frootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livestemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadstemn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_livecrootn_xfer_to_litter(:)
    real(rk8), pointer, contiguous :: m_deadcrootn_xfer_to_litter(:)
    ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_c_to_litr_met_c(:,:)
    ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_c_to_litr_cel_c(:,:)
    ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_c_to_litr_lig_c(:,:)
    ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_c_to_cwdc(:,:)
    ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_met_n(:,:)
    ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_cel_n(:,:)
    ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_litr_lig_n(:,:)
    ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous :: gap_mortality_n_to_cwdn(:,:)

    real(rk8), pointer, contiguous :: leaf_prof(:,:)  ! (1/m) profile of leaves
    real(rk8), pointer, contiguous :: froot_prof(:,:) ! (1/m) profile of fine roots
    real(rk8), pointer, contiguous :: croot_prof(:,:) ! (1/m) profile of coarse roots
    real(rk8), pointer, contiguous :: stem_prof(:,:)  ! (1/m) profile of stems

    integer(ik4) :: fc,c,pi,p,j  ! indices

    ! assign local pointers
    lf_flab   => pftcon%lf_flab
    lf_fcel   => pftcon%lf_fcel
    lf_flig   => pftcon%lf_flig
    fr_flab   => pftcon%fr_flab
    fr_fcel   => pftcon%fr_fcel
    fr_flig   => pftcon%fr_flig

    ! assign local pointers to column-level arrays
    npfts                          => clm3%g%l%c%npfts
    pfti                           => clm3%g%l%c%pfti

    ! assign local pointers to pft-level arrays
    pactive                        => clm3%g%l%c%p%active
    ivt                            => clm3%g%l%c%p%itype
    wtcol                          => clm3%g%l%c%p%wtcol
    m_leafc_to_litter              => clm3%g%l%c%p%pcf%m_leafc_to_litter
    m_frootc_to_litter             => clm3%g%l%c%p%pcf%m_frootc_to_litter
    m_livestemc_to_litter          => clm3%g%l%c%p%pcf%m_livestemc_to_litter
    m_deadstemc_to_litter          => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
    m_livecrootc_to_litter         => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
    m_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
    m_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
    m_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
    m_livestemc_storage_to_litter  => &
            clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
    m_deadstemc_storage_to_litter  => &
            clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
    m_livecrootc_storage_to_litter => &
            clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
    m_deadcrootc_storage_to_litter => &
            clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
    m_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
    m_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
    m_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
    m_livestemc_xfer_to_litter => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
    m_deadstemc_xfer_to_litter => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
    m_livecrootc_xfer_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
    m_deadcrootc_xfer_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
    m_gresp_xfer_to_litter      => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
    m_leafn_to_litter           => clm3%g%l%c%p%pnf%m_leafn_to_litter
    m_frootn_to_litter          => clm3%g%l%c%p%pnf%m_frootn_to_litter
    m_livestemn_to_litter       => clm3%g%l%c%p%pnf%m_livestemn_to_litter
    m_deadstemn_to_litter       => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
    m_livecrootn_to_litter      => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
    m_deadcrootn_to_litter      => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
    m_retransn_to_litter        => clm3%g%l%c%p%pnf%m_retransn_to_litter
    m_leafn_storage_to_litter   => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
    m_frootn_storage_to_litter  => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
    m_livestemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
    m_deadstemn_storage_to_litter  => &
            clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
    m_livecrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
    m_deadcrootn_storage_to_litter => &
            clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
    m_leafn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
    m_frootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
    m_livestemn_xfer_to_litter => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
    m_deadstemn_xfer_to_litter => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
    m_livecrootn_xfer_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
    m_deadcrootn_xfer_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
    gap_mortality_c_to_litr_met_c => &
            clm3%g%l%c%ccf%gap_mortality_c_to_litr_met_c
    gap_mortality_c_to_litr_cel_c => &
            clm3%g%l%c%ccf%gap_mortality_c_to_litr_cel_c
    gap_mortality_c_to_litr_lig_c => &
            clm3%g%l%c%ccf%gap_mortality_c_to_litr_lig_c
    gap_mortality_c_to_cwdc        => clm3%g%l%c%ccf%gap_mortality_c_to_cwdc
    gap_mortality_n_to_litr_met_n => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_met_n
    gap_mortality_n_to_litr_cel_n => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_cel_n
    gap_mortality_n_to_litr_lig_n => &
            clm3%g%l%c%cnf%gap_mortality_n_to_litr_lig_n
    gap_mortality_n_to_cwdn => clm3%g%l%c%cnf%gap_mortality_n_to_cwdn
    leaf_prof               => clm3%g%l%c%p%pps%leaf_prof
    froot_prof              => clm3%g%l%c%p%pps%froot_prof
    croot_prof              => clm3%g%l%c%p%pps%croot_prof
    stem_prof               => clm3%g%l%c%p%pps%stem_prof

    do j = 1, nlevdecomp
      do pi = 1, maxpatch_pft
        do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1
            if (pactive(p)) then
              ! leaf gap mortality carbon fluxes
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j) + &
                      m_leafc_to_litter(p) * lf_flab(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)
              gap_mortality_c_to_litr_cel_c(c,j) = &
                      gap_mortality_c_to_litr_cel_c(c,j) + &
                      m_leafc_to_litter(p) * lf_fcel(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)
              gap_mortality_c_to_litr_lig_c(c,j) = &
                      gap_mortality_c_to_litr_lig_c(c,j) + &
                      m_leafc_to_litter(p) * lf_flig(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)

              ! fine root gap mortality carbon fluxes
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j) + &
                      m_frootc_to_litter(p) * fr_flab(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)
              gap_mortality_c_to_litr_cel_c(c,j) = &
                      gap_mortality_c_to_litr_cel_c(c,j) + &
                      m_frootc_to_litter(p) * fr_fcel(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)
              gap_mortality_c_to_litr_lig_c(c,j) = &
                      gap_mortality_c_to_litr_lig_c(c,j) + &
                      m_frootc_to_litter(p) * fr_flig(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)

              ! wood gap mortality carbon fluxes
              gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                       (m_livestemc_to_litter(p) + &
                       m_deadstemc_to_litter(p)) * wtcol(p) * stem_prof(p,j)
              gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                       (m_livecrootc_to_litter(p) + &
                       m_deadcrootc_to_litter(p)) * wtcol(p) * croot_prof(p,j)

              ! storage gap mortality carbon fluxes
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j) + &
                      (m_leafc_storage_to_litter(p) + &
                      m_gresp_storage_to_litter(p)) * wtcol(p) * leaf_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j)     + &
                      m_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j)  + &
                      (m_livestemc_storage_to_litter(p) + &
                      m_deadstemc_storage_to_litter(p)) * &
                      wtcol(p) * stem_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j) + &
                      (m_livecrootc_storage_to_litter(p) + &
                      m_deadcrootc_storage_to_litter(p)) * &
                      wtcol(p) * croot_prof(p,j)

              ! transfer gap mortality carbon fluxes
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j)      + &
                      (m_leafc_xfer_to_litter(p) + &
                      m_gresp_xfer_to_litter(p))     * wtcol(p) * leaf_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j)     + &
                      m_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j)  + &
                      (m_livestemc_xfer_to_litter(p) + &
                      m_deadstemc_xfer_to_litter(p)) * wtcol(p) * stem_prof(p,j)
              gap_mortality_c_to_litr_met_c(c,j) = &
                      gap_mortality_c_to_litr_met_c(c,j) + &
                      (m_livecrootc_xfer_to_litter(p) + &
                      m_deadcrootc_xfer_to_litter(p)) * &
                      wtcol(p) * croot_prof(p,j)

              ! leaf gap mortality nitrogen fluxes
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      m_leafn_to_litter(p) * lf_flab(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)
              gap_mortality_n_to_litr_cel_n(c,j) = &
                      gap_mortality_n_to_litr_cel_n(c,j) + &
                      m_leafn_to_litter(p) * lf_fcel(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)
              gap_mortality_n_to_litr_lig_n(c,j) = &
                      gap_mortality_n_to_litr_lig_n(c,j) + &
                      m_leafn_to_litter(p) * lf_flig(ivt(p)) * &
                      wtcol(p) * leaf_prof(p,j)

              ! fine root litter nitrogen fluxes
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      m_frootn_to_litter(p) * fr_flab(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)
              gap_mortality_n_to_litr_cel_n(c,j) = &
                      gap_mortality_n_to_litr_cel_n(c,j) + &
                      m_frootn_to_litter(p) * fr_fcel(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)
              gap_mortality_n_to_litr_lig_n(c,j) = &
                      gap_mortality_n_to_litr_lig_n(c,j) + &
                      m_frootn_to_litter(p) * fr_flig(ivt(p)) * &
                      wtcol(p) * froot_prof(p,j)

              ! wood gap mortality nitrogen fluxes
              gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j)  + &
                       (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p)) * &
                       wtcol(p) * stem_prof(p,j)
              gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                       (m_livecrootn_to_litter(p) + &
                       m_deadcrootn_to_litter(p)) * &
                       wtcol(p) * croot_prof(p,j)

              ! retranslocated N pool gap mortality fluxes
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

              ! storage gap mortality nitrogen fluxes
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j)      + &
                      m_leafn_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j)     + &
                      m_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      (m_livestemn_storage_to_litter(p) + &
                      m_deadstemn_storage_to_litter(p)) * &
                      wtcol(p) * stem_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      (m_livecrootn_storage_to_litter(p) + &
                      m_deadcrootn_storage_to_litter(p)) * &
                      wtcol(p) * croot_prof(p,j)

              ! transfer gap mortality nitrogen fluxes
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j)      + &
                      m_leafn_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j)     + &
                      m_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j)  = &
                      gap_mortality_n_to_litr_met_n(c,j)  + &
                      (m_livestemn_xfer_to_litter(p) + &
                      m_deadstemn_xfer_to_litter(p)) * wtcol(p) * stem_prof(p,j)
              gap_mortality_n_to_litr_met_n(c,j) = &
                      gap_mortality_n_to_litr_met_n(c,j) + &
                      (m_livecrootn_xfer_to_litter(p) + &
                      m_deadcrootn_xfer_to_litter(p)) * &
                      wtcol(p) * croot_prof(p,j)
            end if
          end if
        end do
      end do
    end do
  end subroutine CNGapPftToColumn

#endif

end module mod_clm_cngapmortality
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
