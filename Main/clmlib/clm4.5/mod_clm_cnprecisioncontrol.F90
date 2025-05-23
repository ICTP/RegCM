module mod_clm_cnprecisioncontrol

#ifdef CN
  !
  ! controls on very low values in critical state variables
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_clm_varpar, only : ndecomp_pools

  implicit none

  save

  private

  public:: CNPrecisionControl

  contains
  !
  ! On the radiation time step, force leaf and deadstem c and n to 0 if
  ! they get too small.
  !
  subroutine CNPrecisionControl(num_soilc,filter_soilc,num_soilp,filter_soilp)
    use mod_clm_type
    use mod_clm_varctl, only : use_c13, use_c14
    use mod_clm_varpar, only : nlevdecomp
    use mod_clm_pftvarcon, only : nc3crop
    use mod_clm_surfrd, only : crop_prog
    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    ! (gC/m3) column-level sink for C truncation
    real(rk8), pointer, contiguous :: col_ctrunc_vr(:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous :: decomp_cpools_vr(:,:,:)
    ! (gN/m3) column-level sink for N truncation
    real(rk8), pointer, contiguous :: col_ntrunc_vr(:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous :: decomp_npools_vr(:,:,:)
    ! (gC/m2) temporary photosynthate C pool
    real(rk8), pointer, contiguous :: cpool(:)
    real(rk8), pointer, contiguous :: deadcrootc(:)         ! (gC/m2) dead coarse root C
    real(rk8), pointer, contiguous :: deadcrootc_storage(:) ! (gC/m2) dead cors root C strg
    real(rk8), pointer, contiguous :: deadcrootc_xfer(:)    ! (gC/m2) dead cors root C trnf
    real(rk8), pointer, contiguous :: deadstemc(:)          ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
    real(rk8), pointer, contiguous :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
    real(rk8), pointer, contiguous :: frootc(:)             ! (gC/m2) fine root C
    real(rk8), pointer, contiguous :: frootc_storage(:)     ! (gC/m2) fine root C storage
    real(rk8), pointer, contiguous :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
    real(rk8), pointer, contiguous :: gresp_storage(:)      ! (gC/m2) growth resp. strg
    real(rk8), pointer, contiguous :: gresp_xfer(:)         ! (gC/m2) growth resp. trnfr
    real(rk8), pointer, contiguous :: leafc(:)              ! (gC/m2) leaf C
    real(rk8), pointer, contiguous :: leafc_storage(:)      ! (gC/m2) leaf C storage
    real(rk8), pointer, contiguous :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
    real(rk8), pointer, contiguous :: livecrootc(:)         ! (gC/m2) live coarse root C
    real(rk8), pointer, contiguous :: livecrootc_storage(:) ! (gC/m2) live cors root C strg
    real(rk8), pointer, contiguous :: livecrootc_xfer(:)    ! (gC/m2) live cors root C trnf
    real(rk8), pointer, contiguous :: livestemc(:)          ! (gC/m2) live stem C
    real(rk8), pointer, contiguous :: livestemc_storage(:)  ! (gC/m2) live stem C storage
    real(rk8), pointer, contiguous :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
    real(rk8), pointer, contiguous :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C
                                                ! truncation
    real(rk8), pointer, contiguous :: xsmrpool(:)     ! (gC/m2) execss maint resp C pool
    real(rk8), pointer, contiguous :: grainc(:)       ! (gC/m2) grain C
    real(rk8), pointer, contiguous :: grainc_storage(:)  ! (gC/m2) grain C storage
    real(rk8), pointer, contiguous :: grainc_xfer(:)     ! (gC/m2) grain C transfer

    !!! C13
    ! (gC/m3) column-level sink for C truncation
    real(rk8), pointer, contiguous :: c13_col_ctrunc_vr(:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous :: decomp_c13pools_vr(:,:,:)
    ! (gC/m2) temporary photosynthate C pool
    real(rk8), pointer, contiguous :: c13_cpool(:)
    ! (gC/m2) dead coarse root C
    real(rk8), pointer, contiguous :: c13_deadcrootc(:)
    ! (gC/m2) dead coarse root C storage
    real(rk8), pointer, contiguous :: c13_deadcrootc_storage(:)
    ! (gC/m2) dead coarse root C transfer
    real(rk8), pointer, contiguous :: c13_deadcrootc_xfer(:)
    real(rk8), pointer, contiguous :: c13_deadstemc(:)          ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous :: c13_deadstemc_storage(:)  ! (gC/m2) dead stem C stor
    real(rk8), pointer, contiguous :: c13_deadstemc_xfer(:)     ! (gC/m2) dead stem C trnsf
    real(rk8), pointer, contiguous :: c13_frootc(:)             ! (gC/m2) fine root C
    real(rk8), pointer, contiguous :: c13_frootc_storage(:)     ! (gC/m2) fine root C stor
    real(rk8), pointer, contiguous :: c13_frootc_xfer(:)        ! (gC/m2) fine root C trnsf
    real(rk8), pointer, contiguous :: c13_gresp_storage(:) ! (gC/m2) growth respiration strg
    real(rk8), pointer, contiguous :: c13_gresp_xfer(:)    ! (gC/m2) growth respiration trnf
    real(rk8), pointer, contiguous :: c13_leafc(:)           ! (gC/m2) leaf C
    real(rk8), pointer, contiguous :: c13_leafc_storage(:)   ! (gC/m2) leaf C storage
    real(rk8), pointer, contiguous :: c13_leafc_xfer(:)      ! (gC/m2) leaf C transfer
    real(rk8), pointer, contiguous :: c13_livecrootc(:)      ! (gC/m2) live coarse root C
    ! (gC/m2) live coarse root C storage
    real(rk8), pointer, contiguous :: c13_livecrootc_storage(:)
    ! (gC/m2) live coarse root C transfer
    real(rk8), pointer, contiguous :: c13_livecrootc_xfer(:)
    real(rk8), pointer, contiguous :: c13_livestemc(:)          ! (gC/m2) live stem C
    real(rk8), pointer, contiguous :: c13_livestemc_storage(:)  ! (gC/m2) live stem C stor
    real(rk8), pointer, contiguous :: c13_livestemc_xfer(:)     ! (gC/m2) live stem C transf
    ! (gC/m2) pft-level sink for C truncation
    real(rk8), pointer, contiguous :: c13_pft_ctrunc(:)

    !!! C14
    ! (gC/m3) column-level sink for C truncation
    real(rk8), pointer, contiguous :: c14_col_ctrunc_vr(:,:)
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous :: decomp_c14pools_vr(:,:,:)
    ! (gC/m2) temporary photosynthate C pool
    real(rk8), pointer, contiguous :: c14_cpool(:)
    ! (gC/m2) dead coarse root C
    real(rk8), pointer, contiguous :: c14_deadcrootc(:)
    ! (gC/m2) dead coarse root C storage
    real(rk8), pointer, contiguous :: c14_deadcrootc_storage(:)
    ! (gC/m2) dead coarse root C transfer
    real(rk8), pointer, contiguous :: c14_deadcrootc_xfer(:)
    real(rk8), pointer, contiguous :: c14_deadstemc(:)         ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous :: c14_deadstemc_storage(:) ! (gC/m2) dead stem C storage
    real(rk8), pointer, contiguous :: c14_deadstemc_xfer(:)    ! (gC/m2) dead stem C transf
    real(rk8), pointer, contiguous :: c14_frootc(:)            ! (gC/m2) fine root C
    real(rk8), pointer, contiguous :: c14_frootc_storage(:)    ! (gC/m2) fine root C storage
    real(rk8), pointer, contiguous :: c14_frootc_xfer(:)       ! (gC/m2) fine root C transf
    real(rk8), pointer, contiguous :: c14_gresp_storage(:) ! (gC/m2) growth respiration strg
    real(rk8), pointer, contiguous :: c14_gresp_xfer(:)    ! (gC/m2) growth respiration trnf
    real(rk8), pointer, contiguous :: c14_leafc(:)             ! (gC/m2) leaf C
    real(rk8), pointer, contiguous :: c14_leafc_storage(:)     ! (gC/m2) leaf C storage
    real(rk8), pointer, contiguous :: c14_leafc_xfer(:)        ! (gC/m2) leaf C transfer
    real(rk8), pointer, contiguous :: c14_livecrootc(:)        ! (gC/m2) live coarse root C
    ! (gC/m2) live coarse root C storage
    real(rk8), pointer, contiguous :: c14_livecrootc_storage(:)
    ! (gC/m2) live coarse root C transfer
    real(rk8), pointer, contiguous :: c14_livecrootc_xfer(:)
    real(rk8), pointer, contiguous :: c14_livestemc(:)         ! (gC/m2) live stem C
    real(rk8), pointer, contiguous :: c14_livestemc_storage(:) ! (gC/m2) live stem C storage
    real(rk8), pointer, contiguous :: c14_livestemc_xfer(:)    ! (gC/m2) live stem C trans
    ! (gC/m2) pft-level sink for C truncation
    real(rk8), pointer, contiguous :: c14_pft_ctrunc(:)

    real(rk8), pointer, contiguous :: deadcrootn(:)         ! (gN/m2) dead coarse root N
    real(rk8), pointer, contiguous :: deadcrootn_storage(:) ! (gN/m2) dead cors root N strg
    real(rk8), pointer, contiguous :: deadcrootn_xfer(:)    ! (gN/m2) dead cors root N trnf
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
    real(rk8), pointer, contiguous :: livecrootn_storage(:) ! (gN/m2) live cors root N strg
    real(rk8), pointer, contiguous :: livecrootn_xfer(:)    ! (gN/m2) live cors root N trnf
    real(rk8), pointer, contiguous :: grainn(:)             ! (gC/m2) grain N
    real(rk8), pointer, contiguous :: grainn_storage(:)     ! (gC/m2) grain N storage
    real(rk8), pointer, contiguous :: grainn_xfer(:)        ! (gC/m2) grain N transfer
    real(rk8), pointer, contiguous :: livestemn(:)          ! (gN/m2) live stem N
    real(rk8), pointer, contiguous :: livestemn_storage(:)  ! (gN/m2) live stem N storage
    real(rk8), pointer, contiguous :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
    real(rk8), pointer, contiguous :: npool(:)              ! (gN/m2) temporary plant N pool
    real(rk8), pointer, contiguous :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for
                                                ! N truncation
    real(rk8), pointer, contiguous :: retransn(:) ! (gN/m2) plant pool of retranslocated N
#ifdef NITRIF_DENITRIF
    real(rk8), pointer, contiguous :: smin_no3_vr(:,:)   ! (gN/m3) soil mineral NO3
    real(rk8), pointer, contiguous :: smin_nh4_vr(:,:)   ! (gN/m3) soil mineral NH4
#endif
    integer(ik4), pointer, contiguous :: ivt(:)   ! pft vegetation type

    integer(ik4) :: c,p,j,k    ! indices
    integer(ik4) :: fp,fc    ! lake filter indices
    real(rk8):: pc,pn    ! truncation terms for pft-level corrections
    real(rk8):: cc,cn    ! truncation terms for column-level corrections
    !!! C13
    real(rk8):: pc13     ! truncation terms for pft-level corrections
    real(rk8):: cc13     ! truncation terms for column-level corrections
    !!! C14
    real(rk8):: pc14     ! truncation terms for pft-level corrections
    real(rk8):: cc14     ! truncation terms for column-level corrections

    real(rk8):: ccrit    ! critical carbon state value for truncation
    real(rk8):: ncrit    ! critical nitrogen state value for truncation

    ! assign local pointers at the column level
    col_ctrunc_vr                  => clm3%g%l%c%ccs%col_ctrunc_vr
    decomp_cpools_vr               => clm3%g%l%c%ccs%decomp_cpools_vr
    col_ntrunc_vr                  => clm3%g%l%c%cns%col_ntrunc_vr
    decomp_npools_vr               => clm3%g%l%c%cns%decomp_npools_vr
    ! assign local pointers at the pft level
    ivt                            => clm3%g%l%c%p%itype
    cpool                          => clm3%g%l%c%p%pcs%cpool
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
    xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
    grainc                         => clm3%g%l%c%p%pcs%grainc
    grainc_storage                 => clm3%g%l%c%p%pcs%grainc_storage
    grainc_xfer                    => clm3%g%l%c%p%pcs%grainc_xfer

    if ( use_c13 ) then
      c13_col_ctrunc_vr        => clm3%g%l%c%cc13s%col_ctrunc_vr
      decomp_c13pools_vr       => clm3%g%l%c%cc13s%decomp_cpools_vr
      c13_cpool                => clm3%g%l%c%p%pc13s%cpool
      c13_deadcrootc           => clm3%g%l%c%p%pc13s%deadcrootc
      c13_deadcrootc_storage   => clm3%g%l%c%p%pc13s%deadcrootc_storage
      c13_deadcrootc_xfer      => clm3%g%l%c%p%pc13s%deadcrootc_xfer
      c13_deadstemc            => clm3%g%l%c%p%pc13s%deadstemc
      c13_deadstemc_storage    => clm3%g%l%c%p%pc13s%deadstemc_storage
      c13_deadstemc_xfer       => clm3%g%l%c%p%pc13s%deadstemc_xfer
      c13_frootc               => clm3%g%l%c%p%pc13s%frootc
      c13_frootc_storage       => clm3%g%l%c%p%pc13s%frootc_storage
      c13_frootc_xfer          => clm3%g%l%c%p%pc13s%frootc_xfer
      c13_gresp_storage        => clm3%g%l%c%p%pc13s%gresp_storage
      c13_gresp_xfer           => clm3%g%l%c%p%pc13s%gresp_xfer
      c13_leafc                => clm3%g%l%c%p%pc13s%leafc
      c13_leafc_storage        => clm3%g%l%c%p%pc13s%leafc_storage
      c13_leafc_xfer           => clm3%g%l%c%p%pc13s%leafc_xfer
      c13_livecrootc           => clm3%g%l%c%p%pc13s%livecrootc
      c13_livecrootc_storage   => clm3%g%l%c%p%pc13s%livecrootc_storage
      c13_livecrootc_xfer      => clm3%g%l%c%p%pc13s%livecrootc_xfer
      c13_livestemc            => clm3%g%l%c%p%pc13s%livestemc
      c13_livestemc_storage    => clm3%g%l%c%p%pc13s%livestemc_storage
      c13_livestemc_xfer       => clm3%g%l%c%p%pc13s%livestemc_xfer
      c13_pft_ctrunc           => clm3%g%l%c%p%pc13s%pft_ctrunc
    end if

    if ( use_c14 ) then
      c14_col_ctrunc_vr        => clm3%g%l%c%cc14s%col_ctrunc_vr
      decomp_c14pools_vr       => clm3%g%l%c%cc14s%decomp_cpools_vr
      c14_cpool                => clm3%g%l%c%p%pc14s%cpool
      c14_deadcrootc           => clm3%g%l%c%p%pc14s%deadcrootc
      c14_deadcrootc_storage   => clm3%g%l%c%p%pc14s%deadcrootc_storage
      c14_deadcrootc_xfer      => clm3%g%l%c%p%pc14s%deadcrootc_xfer
      c14_deadstemc            => clm3%g%l%c%p%pc14s%deadstemc
      c14_deadstemc_storage    => clm3%g%l%c%p%pc14s%deadstemc_storage
      c14_deadstemc_xfer       => clm3%g%l%c%p%pc14s%deadstemc_xfer
      c14_frootc               => clm3%g%l%c%p%pc14s%frootc
      c14_frootc_storage       => clm3%g%l%c%p%pc14s%frootc_storage
      c14_frootc_xfer          => clm3%g%l%c%p%pc14s%frootc_xfer
      c14_gresp_storage        => clm3%g%l%c%p%pc14s%gresp_storage
      c14_gresp_xfer           => clm3%g%l%c%p%pc14s%gresp_xfer
      c14_leafc                => clm3%g%l%c%p%pc14s%leafc
      c14_leafc_storage        => clm3%g%l%c%p%pc14s%leafc_storage
      c14_leafc_xfer           => clm3%g%l%c%p%pc14s%leafc_xfer
      c14_livecrootc           => clm3%g%l%c%p%pc14s%livecrootc
      c14_livecrootc_storage   => clm3%g%l%c%p%pc14s%livecrootc_storage
      c14_livecrootc_xfer      => clm3%g%l%c%p%pc14s%livecrootc_xfer
      c14_livestemc            => clm3%g%l%c%p%pc14s%livestemc
      c14_livestemc_storage    => clm3%g%l%c%p%pc14s%livestemc_storage
      c14_livestemc_xfer       => clm3%g%l%c%p%pc14s%livestemc_xfer
      c14_pft_ctrunc           => clm3%g%l%c%p%pc14s%pft_ctrunc
    end if

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
    grainn                     => clm3%g%l%c%p%pns%grainn
    grainn_storage             => clm3%g%l%c%p%pns%grainn_storage
    grainn_xfer                => clm3%g%l%c%p%pns%grainn_xfer
    livestemn                  => clm3%g%l%c%p%pns%livestemn
    livestemn_storage          => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer             => clm3%g%l%c%p%pns%livestemn_xfer
    npool                      => clm3%g%l%c%p%pns%npool
    pft_ntrunc                 => clm3%g%l%c%p%pns%pft_ntrunc
    retransn                   => clm3%g%l%c%p%pns%retransn

#ifdef NITRIF_DENITRIF
   smin_nh4_vr => clm3%g%l%c%cns%smin_nh4_vr
   smin_no3_vr => clm3%g%l%c%cns%smin_no3_vr
#endif

    ! set the critical carbon state value for truncation (gC/m2)
    ccrit = 1.e-8_rk8
    ! set the critical nitrogen state value for truncation (gN/m2)
    ncrit = 1.e-8_rk8

    ! pft loop
    do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! initialize the pft-level C and N truncation terms
      pc = 0._rk8

      if ( use_c13 ) then
        pc13 = 0._rk8
      end if
      if ( use_c14 ) then
        pc14 = 0._rk8
      end if

      pn = 0._rk8

      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate C, C13, and N components

      ! leaf C and N
      if (abs(leafc(p)) < ccrit) then
        pc = pc + leafc(p)
        leafc(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_leafc(p)
          c13_leafc(p) = 0._rk8
        end if

        if ( use_c14 ) then
          pc14 = pc14 + c14_leafc(p)
          c14_leafc(p) = 0._rk8
        end if
        pn = pn + leafn(p)
        leafn(p) = 0._rk8
      end if

      ! leaf storage C and N
      if (abs(leafc_storage(p)) < ccrit) then
        pc = pc + leafc_storage(p)
        leafc_storage(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_leafc_storage(p)
          c13_leafc_storage(p) = 0._rk8
        end if

        if ( use_c14 ) then
          pc14 = pc14 + c14_leafc_storage(p)
          c14_leafc_storage(p) = 0._rk8
        end if

        pn = pn + leafn_storage(p)
        leafn_storage(p) = 0._rk8
      end if

      ! leaf transfer C and N
      if (abs(leafc_xfer(p)) < ccrit) then
        pc = pc + leafc_xfer(p)
        leafc_xfer(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_leafc_xfer(p)
          c13_leafc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_leafc_xfer(p)
          c14_leafc_xfer(p) = 0._rk8
        end if

        pn = pn + leafn_xfer(p)
        leafn_xfer(p) = 0._rk8
      end if

      ! froot C and N
      if (abs(frootc(p)) < ccrit) then
        pc = pc + frootc(p)
        frootc(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_frootc(p)
          c13_frootc(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_frootc(p)
          c14_frootc(p) = 0._rk8
        end if

        pn = pn + frootn(p)
        frootn(p) = 0._rk8
      end if

      ! froot storage C and N
      if (abs(frootc_storage(p)) < ccrit) then
        pc = pc + frootc_storage(p)
        frootc_storage(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_frootc_storage(p)
          c13_frootc_storage(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_frootc_storage(p)
          c14_frootc_storage(p) = 0._rk8
        end if

        pn = pn + frootn_storage(p)
        frootn_storage(p) = 0._rk8
      end if

      ! froot transfer C and N
      if (abs(frootc_xfer(p)) < ccrit) then
        pc = pc + frootc_xfer(p)
        frootc_xfer(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_frootc_xfer(p)
          c13_frootc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_frootc_xfer(p)
          c14_frootc_xfer(p) = 0._rk8
        end if

        pn = pn + frootn_xfer(p)
        frootn_xfer(p) = 0._rk8
      end if

      if ( crop_prog .and. ivt(p) >= nc3crop )then
        ! grain C and N
        if (abs(grainc(p)) < ccrit) then
          pc = pc + grainc(p)
          grainc(p) = 0._rk8
          pn = pn + grainn(p)
          grainn(p) = 0._rk8
        end if

        ! grain storage C and N
        if (abs(grainc_storage(p)) < ccrit) then
          pc = pc + grainc_storage(p)
          grainc_storage(p) = 0._rk8
          pn = pn + grainn_storage(p)
          grainn_storage(p) = 0._rk8
        end if

        ! grain transfer C and N
        if (abs(grainc_xfer(p)) < ccrit) then
          pc = pc + grainc_xfer(p)
          grainc_xfer(p) = 0._rk8
          pn = pn + grainn_xfer(p)
          grainn_xfer(p) = 0._rk8
        end if
      end if

      ! livestem C and N
      if (abs(livestemc(p)) < ccrit) then
        pc = pc + livestemc(p)
        livestemc(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_livestemc(p)
          c13_livestemc(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_livestemc(p)
          c14_livestemc(p) = 0._rk8
        end if

        pn = pn + livestemn(p)
        livestemn(p) = 0._rk8
      end if

      ! livestem storage C and N
      if (abs(livestemc_storage(p)) < ccrit) then
        pc = pc + livestemc_storage(p)
        livestemc_storage(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_livestemc_storage(p)
          c13_livestemc_storage(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_livestemc_storage(p)
          c14_livestemc_storage(p) = 0._rk8
        end if
        pn = pn + livestemn_storage(p)
        livestemn_storage(p) = 0._rk8
      end if

      ! livestem transfer C and N
      if (abs(livestemc_xfer(p)) < ccrit) then
        pc = pc + livestemc_xfer(p)
        livestemc_xfer(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_livestemc_xfer(p)
          c13_livestemc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_livestemc_xfer(p)
          c14_livestemc_xfer(p) = 0._rk8
        end if
        pn = pn + livestemn_xfer(p)
        livestemn_xfer(p) = 0._rk8
      end if

      ! deadstem C and N
      if (abs(deadstemc(p)) < ccrit) then
        pc = pc + deadstemc(p)
        deadstemc(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadstemc(p)
          c13_deadstemc(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadstemc(p)
          c14_deadstemc(p) = 0._rk8
        end if
        pn = pn + deadstemn(p)
        deadstemn(p) = 0._rk8
      end if

      ! deadstem storage C and N
      if (abs(deadstemc_storage(p)) < ccrit) then
        pc = pc + deadstemc_storage(p)
        deadstemc_storage(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadstemc_storage(p)
          c13_deadstemc_storage(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadstemc_storage(p)
          c14_deadstemc_storage(p) = 0._rk8
        end if
        pn = pn + deadstemn_storage(p)
        deadstemn_storage(p) = 0._rk8
      end if

      ! deadstem transfer C and N
      if (abs(deadstemc_xfer(p)) < ccrit) then
        pc = pc + deadstemc_xfer(p)
        deadstemc_xfer(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadstemc_xfer(p)
          c13_deadstemc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadstemc_xfer(p)
          c14_deadstemc_xfer(p) = 0._rk8
        end if
        pn = pn + deadstemn_xfer(p)
        deadstemn_xfer(p) = 0._rk8
      end if

      ! livecroot C and N
      if (abs(livecrootc(p)) < ccrit) then
        pc = pc + livecrootc(p)
        livecrootc(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_livecrootc(p)
          c13_livecrootc(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_livecrootc(p)
          c14_livecrootc(p) = 0._rk8
        end if
        pn = pn + livecrootn(p)
        livecrootn(p) = 0._rk8
      end if

      ! livecroot storage C and N
      if (abs(livecrootc_storage(p)) < ccrit) then
        pc = pc + livecrootc_storage(p)
        livecrootc_storage(p) = 0._rk8

        if ( use_c13 ) then
          pc13 = pc13 + c13_livecrootc_storage(p)
          c13_livecrootc_storage(p) = 0._rk8
        end if

        if ( use_c14 ) then
          pc14 = pc14 + c14_livecrootc_storage(p)
          c14_livecrootc_storage(p) = 0._rk8
        end if

        pn = pn + livecrootn_storage(p)
        livecrootn_storage(p) = 0._rk8
      end if

      ! livecroot transfer C and N
      if (abs(livecrootc_xfer(p)) < ccrit) then
        pc = pc + livecrootc_xfer(p)
        livecrootc_xfer(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_livecrootc_xfer(p)
          c13_livecrootc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_livecrootc_xfer(p)
          c14_livecrootc_xfer(p) = 0._rk8
        end if
        pn = pn + livecrootn_xfer(p)
        livecrootn_xfer(p) = 0._rk8
      end if

      ! deadcroot C and N
      if (abs(deadcrootc(p)) < ccrit) then
        pc = pc + deadcrootc(p)
        deadcrootc(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadcrootc(p)
          c13_deadcrootc(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadcrootc(p)
          c14_deadcrootc(p) = 0._rk8
        end if
        pn = pn + deadcrootn(p)
        deadcrootn(p) = 0._rk8
      end if

      ! deadcroot storage C and N
      if (abs(deadcrootc_storage(p)) < ccrit) then
        pc = pc + deadcrootc_storage(p)
        deadcrootc_storage(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadcrootc_storage(p)
          c13_deadcrootc_storage(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadcrootc_storage(p)
          c14_deadcrootc_storage(p) = 0._rk8
        end if
        pn = pn + deadcrootn_storage(p)
        deadcrootn_storage(p) = 0._rk8
      end if

      ! deadcroot transfer C and N
      if (abs(deadcrootc_xfer(p)) < ccrit) then
        pc = pc + deadcrootc_xfer(p)
        deadcrootc_xfer(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_deadcrootc_xfer(p)
          c13_deadcrootc_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_deadcrootc_xfer(p)
          c14_deadcrootc_xfer(p) = 0._rk8
        end if
        pn = pn + deadcrootn_xfer(p)
        deadcrootn_xfer(p) = 0._rk8
      end if

      ! gresp_storage (C only)
      if (abs(gresp_storage(p)) < ccrit) then
        pc = pc + gresp_storage(p)
        gresp_storage(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_gresp_storage(p)
          c13_gresp_storage(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_gresp_storage(p)
          c14_gresp_storage(p) = 0._rk8
        end if
      end if

      ! gresp_xfer (C only)
      if (abs(gresp_xfer(p)) < ccrit) then
        pc = pc + gresp_xfer(p)
        gresp_xfer(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_gresp_xfer(p)
          c13_gresp_xfer(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_gresp_xfer(p)
          c14_gresp_xfer(p) = 0._rk8
        end if
      end if

      ! cpool (C only)
      if (abs(cpool(p)) < ccrit) then
        pc = pc + cpool(p)
        cpool(p) = 0._rk8
        if ( use_c13 ) then
          pc13 = pc13 + c13_cpool(p)
          c13_cpool(p) = 0._rk8
        end if
        if ( use_c14 ) then
          pc14 = pc14 + c14_cpool(p)
          c14_cpool(p) = 0._rk8
        end if
      end if

      if ( crop_prog .and. ivt(p) >= nc3crop )then
        ! xsmrpool (C only)
        if (abs(xsmrpool(p)) < ccrit) then
          pc = pc + xsmrpool(p)
          xsmrpool(p) = 0._rk8
        end if
      end if

      ! retransn (N only)
      if (abs(retransn(p)) < ncrit) then
        pn = pn + retransn(p)
        retransn(p) = 0._rk8
      end if

      ! npool (N only)
      if (abs(npool(p)) < ncrit) then
        pn = pn + npool(p)
        npool(p) = 0._rk8
      end if

      pft_ctrunc(p) = pft_ctrunc(p) + pc
      if ( use_c13 ) then
        c13_pft_ctrunc(p) = c13_pft_ctrunc(p) + pc13
      end if

      if ( use_c14 ) then
        c14_pft_ctrunc(p) = c14_pft_ctrunc(p) + pc14
      end if
      pft_ntrunc(p) = pft_ntrunc(p) + pn
    end do ! end of pft loop

    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      do j = 1, nlevdecomp
        ! initialize the column-level C and N truncation terms
        cc = 0._rk8
        if ( use_c13 ) then
          cc13 = 0._rk8
        end if
        if ( use_c14 ) then
          cc14 = 0._rk8
        end if
        cn = 0._rk8

        ! do tests on state variables for precision control
        ! for linked C-N state variables, perform precision test on
        ! the C component, but truncate both C and N components

        ! all decomposing pools C and N
        do k = 1, ndecomp_pools

          if (abs(decomp_cpools_vr(c,j,k)) < ccrit) then
            cc = cc + decomp_cpools_vr(c,j,k)
            decomp_cpools_vr(c,j,k) = 0._rk8
            if ( use_c13 ) then
              cc13 = cc13 + decomp_c13pools_vr(c,j,k)
              decomp_c13pools_vr(c,j,k) = 0._rk8
            end if
            if ( use_c14 ) then
              cc14 = cc14 + decomp_c14pools_vr(c,j,k)
              decomp_c14pools_vr(c,j,k) = 0._rk8
            end if
            cn = cn + decomp_npools_vr(c,j,k)
            decomp_npools_vr(c,j,k) = 0._rk8
          end if
        end do

        ! not doing precision control on soil mineral N, since it will
        ! be getting the N truncation flux anyway.

        col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) + cc
        if ( use_c13 ) then
          c13_col_ctrunc_vr(c,j) = c13_col_ctrunc_vr(c,j) + cc13
        end if
        if ( use_c14 ) then
          c14_col_ctrunc_vr(c,j) = c14_col_ctrunc_vr(c,j) + cc14
        end if
        col_ntrunc_vr(c,j) = col_ntrunc_vr(c,j) + cn
      end do
    end do   ! end of column loop

#ifdef NITRIF_DENITRIF
    ! remove small negative perturbations for stability purposes,
    ! if any should arise.

    do fc = 1, num_soilc
      c = filter_soilc(fc)
      do j = 1, nlevdecomp
        if (abs(smin_no3_vr(c,j)) < ncrit/1e4_rk8) then
          if ( smin_no3_vr(c,j) < 0._rk8 ) then
            write(stderr, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
            write(stderr, *) 'smin_no3_vr(c,j), c, j: ', smin_no3_vr(c,j), c, j
            smin_no3_vr(c,j) = 0._rk8
          end if
        end if
        if (abs(smin_nh4_vr(c,j)) < ncrit/1e4_rk8) then
          if ( smin_nh4_vr(c,j) < 0._rk8 ) then
            write(stderr, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
            write(stderr, *) 'smin_nh4_vr(c,j), c, j: ', smin_nh4_vr(c,j), c, j
            smin_nh4_vr(c,j) = 0._rk8
          end if
        end if
      end do
    end do
#endif
  end subroutine CNPrecisionControl

#endif

end module mod_clm_cnprecisioncontrol
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
