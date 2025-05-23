module mod_clm_cninispecial
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save

  public :: CNiniSpecial

  contains
  !
  ! One-time initialization of CN variables for special landunits
  !
  subroutine CNiniSpecial ()
#ifdef CN
    use mod_clm_pftvarcon, only : noveg
    use mod_clm_decomp, only : get_proc_bounds
    use mod_clm_varcon, only : spval
    use mod_clm_varpar, only : nlevdecomp_full
    use mod_clm_varctl, only : use_c13, use_c14
    use mod_clm_type
    use mod_clm_cnsetvalue
    use mod_clm_surfrd, only : crop_prog
    implicit none
    ! landunit index of corresponding column
    integer(ik4), pointer, contiguous :: clandunit(:)
    ! landunit index of corresponding pft
    integer(ik4), pointer, contiguous :: plandunit(:)
    ! BOOL: true=>landunit is wetland,ice,lake, or urban
    logical, pointer, contiguous :: ifspecial(:)

    integer(ik4) :: fc,fp,l,c,p,j ! indices
    integer(ik4) :: begp, endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg    ! per-proc gridcell ending gridcell indices
    integer(ik4) :: num_specialc  ! number of good values in specialc filter
    integer(ik4) :: num_specialp  ! number of good values in specialp filter
    integer(ik4), allocatable :: specialc(:) ! special landunit filter - columns
    integer(ik4), allocatable :: specialp(:) ! special landunit filter - pfts

    ! assign local pointers at the landunit level
    ifspecial => clm3%g%l%ifspecial

    ! assign local pointers at the column level
    clandunit => clm3%g%l%c%landunit

    ! assign local pointers at the pft level
    plandunit => clm3%g%l%c%p%landunit

    ! Determine subgrid bounds on this processor
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! allocate special landunit filters
    allocate(specialc(endc-begc+1))
    allocate(specialp(endp-begp+1))

    ! fill special landunit filters
    num_specialc = 0
    do c = begc, endc
      l = clandunit(c)
      if (ifspecial(l)) then
        num_specialc = num_specialc + 1
        specialc(num_specialc) = c
      end if
    end do

    num_specialp = 0
    do p = begp, endp
      l = plandunit(p)
      if (ifspecial(l)) then
        num_specialp = num_specialp + 1
        specialp(num_specialp) = p
      end if
    end do

    ! initialize column-level fields
    call CNSetCps(num_specialc, specialc, spval, clm3%g%l%c%cps)
    call CNSetCcs(num_specialc, specialc, 0._rk8, clm3%g%l%c%ccs)
    call CNSetCns(num_specialc, specialc, 0._rk8, clm3%g%l%c%cns)
    call CNSetCcf(num_specialc, specialc, 0._rk8, clm3%g%l%c%ccf)
    call CNSetCnf(num_specialc, specialc, 0._rk8, clm3%g%l%c%cnf)

    if ( use_c13 ) then
      call CNSetCcs(num_specialc, specialc, 0._rk8, clm3%g%l%c%cc13s)
      call CNSetCcf(num_specialc, specialc, 0._rk8, clm3%g%l%c%cc13f)
    end if

    if ( use_c14 ) then
      call CNSetCcs(num_specialc, specialc, 0._rk8, clm3%g%l%c%cc14s)
      call CNSetCcf(num_specialc, specialc, 0._rk8, clm3%g%l%c%cc14f)
    end if

    ! initialize column-average pft fields
    call CNSetPps(num_specialc, specialc, spval, clm3%g%l%c%cps%pps_a)
    call CNSetPcs(num_specialc, specialc, 0._rk8, clm3%g%l%c%ccs%pcs_a)
    call CNSetPns(num_specialc, specialc, 0._rk8, clm3%g%l%c%cns%pns_a)
    call CNSetPcf(num_specialc, specialc, 0._rk8, clm3%g%l%c%ccf%pcf_a)
    call CNSetPnf(num_specialc, specialc, 0._rk8, clm3%g%l%c%cnf%pnf_a)

    ! initialize pft-level fields
    call CNSetPepv(num_specialp, specialp, spval, clm3%g%l%c%p%pepv)
    call CNSetPps(num_specialp, specialp, spval, clm3%g%l%c%p%pps)
    call CNSetPcs(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pcs)
    call CNSetPns(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pns)
    call CNSetPcf(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pcf)
    call CNSetPnf(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pnf)

    if ( use_c13 ) then
      call CNSetPcs(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pc13s)
      call CNSetPcf(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pc13f)
    end if

    if ( use_c14 ) then
      call CNSetPcs(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pc14s)
      call CNSetPcf(num_specialp, specialp, 0._rk8, clm3%g%l%c%p%pc14f)
    end if

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for SurfaceAlbedo and biogeophysics
    ! also set pcf%psnsun and pcf%psnsha to 0 (not included in CNSetPcf())

    do fp = 1,num_specialp
      p = specialp(fp)
      clm3%g%l%c%p%pps%tlai(p) = 0._rk8
      clm3%g%l%c%p%pps%tsai(p) = 0._rk8
      clm3%g%l%c%p%pps%elai(p) = 0._rk8
      clm3%g%l%c%p%pps%esai(p) = 0._rk8
      clm3%g%l%c%p%pps%htop(p) = 0._rk8
      clm3%g%l%c%p%pps%hbot(p) = 0._rk8
      clm3%g%l%c%p%pps%fwet(p) = 0._rk8
      clm3%g%l%c%p%pps%fdry(p) = 0._rk8
      clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = 0
      clm3%g%l%c%p%pps%frac_veg_nosno(p) = 0
      clm3%g%l%c%p%pcf%psnsun(p) = 0._rk8
      clm3%g%l%c%p%pcf%psnsha(p) = 0._rk8
      if (crop_prog) clm3%g%l%c%p%pnf%fert(p)   = 0._rk8

      if ( use_c13 ) then
        clm3%g%l%c%p%pc13f%psnsun(p) = 0._rk8
        clm3%g%l%c%p%pc13f%psnsha(p) = 0._rk8
      end if

      if ( use_c14 ) then
        clm3%g%l%c%p%pc14f%psnsun(p) = 0._rk8
        clm3%g%l%c%p%pc14f%psnsha(p) = 0._rk8
      end if
    end do

    do fc = 1,num_specialc
      c = specialc(fc)
      clm3%g%l%c%cwf%qflx_irrig(c) = 0._rk8
      clm3%g%l%c%ccf%pcf_a%psnsun(c) = 0._rk8
      clm3%g%l%c%ccf%pcf_a%psnsha(c) = 0._rk8

      if ( use_c13 ) then
        clm3%g%l%c%cc13f%pcf_a%psnsun(c) = 0._rk8
        clm3%g%l%c%cc13f%pcf_a%psnsha(c) = 0._rk8
      end if

      if ( use_c14 ) then
        clm3%g%l%c%cc14f%pcf_a%psnsun(c) = 0._rk8
        clm3%g%l%c%cc14f%pcf_a%psnsha(c) = 0._rk8
      end if

      ! adding dynpft code
      clm3%g%l%c%ccs%seedc(c) = 0._rk8
      clm3%g%l%c%ccs%prod10c(c) = 0._rk8
      clm3%g%l%c%ccs%prod100c(c) = 0._rk8
      clm3%g%l%c%ccs%totprodc(c) = 0._rk8

      if ( use_c13 ) then
        clm3%g%l%c%cc13s%seedc(c) = 0._rk8
        clm3%g%l%c%cc13s%prod10c(c) = 0._rk8
        clm3%g%l%c%cc13s%prod100c(c) = 0._rk8
        clm3%g%l%c%cc13s%totprodc(c) = 0._rk8
      end if

      if ( use_c14 ) then
        clm3%g%l%c%cc14s%seedc(c) = 0._rk8
        clm3%g%l%c%cc14s%prod10c(c) = 0._rk8
        clm3%g%l%c%cc14s%prod100c(c) = 0._rk8
        clm3%g%l%c%cc14s%totprodc(c) = 0._rk8
      end if

      clm3%g%l%c%cns%seedn(c) = 0._rk8
      clm3%g%l%c%cns%prod10n(c) = 0._rk8
      clm3%g%l%c%cns%prod100n(c) = 0._rk8
      clm3%g%l%c%cns%totprodn(c) = 0._rk8
      clm3%g%l%c%ccf%dwt_seedc_to_leaf(c) = 0._rk8
      clm3%g%l%c%ccf%dwt_seedc_to_deadstem(c) = 0._rk8
      clm3%g%l%c%ccf%dwt_conv_cflux(c) = 0._rk8
      clm3%g%l%c%ccf%lf_conv_cflux(c) = 0._rk8    ! F. Li and S. Levis
      clm3%g%l%c%ccf%dwt_prod10c_gain(c) = 0._rk8
      clm3%g%l%c%ccf%prod10c_loss(c) = 0._rk8
      clm3%g%l%c%ccf%dwt_prod100c_gain(c) = 0._rk8
      clm3%g%l%c%ccf%prod100c_loss(c) = 0._rk8
      clm3%g%l%c%ccf%dwt_closs(c) = 0._rk8
      clm3%g%l%c%ccf%landuseflux(c) = 0._rk8
      clm3%g%l%c%ccf%landuptake(c) = 0._rk8
      do j = 1, nlevdecomp_full
        clm3%g%l%c%ccf%dwt_frootc_to_litr_met_c(c,j) = 0._rk8
        clm3%g%l%c%ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._rk8
        clm3%g%l%c%ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._rk8
        clm3%g%l%c%ccf%dwt_livecrootc_to_cwdc(c,j) = 0._rk8
        clm3%g%l%c%ccf%dwt_deadcrootc_to_cwdc(c,j) = 0._rk8
      end do
      if ( use_c13 ) then
        clm3%g%l%c%cc13f%dwt_seedc_to_leaf(c) = 0._rk8
        clm3%g%l%c%cc13f%dwt_seedc_to_deadstem(c) = 0._rk8
        clm3%g%l%c%cc13f%dwt_conv_cflux(c) = 0._rk8
        clm3%g%l%c%cc13f%dwt_prod10c_gain(c) = 0._rk8
        clm3%g%l%c%cc13f%prod10c_loss(c) = 0._rk8
        clm3%g%l%c%cc13f%dwt_prod100c_gain(c) = 0._rk8
        clm3%g%l%c%cc13f%prod100c_loss(c) = 0._rk8
        do j = 1, nlevdecomp_full
          clm3%g%l%c%cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._rk8
          clm3%g%l%c%cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._rk8
          clm3%g%l%c%cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._rk8
          clm3%g%l%c%cc13f%dwt_livecrootc_to_cwdc(c,j) = 0._rk8
          clm3%g%l%c%cc13f%dwt_deadcrootc_to_cwdc(c,j) = 0._rk8
        end do
        clm3%g%l%c%cc13f%dwt_closs(c) = 0._rk8
      end if

      if ( use_c14 ) then
        clm3%g%l%c%cc14f%dwt_seedc_to_leaf(c) = 0._rk8
        clm3%g%l%c%cc14f%dwt_seedc_to_deadstem(c) = 0._rk8
        clm3%g%l%c%cc14f%dwt_conv_cflux(c) = 0._rk8
        clm3%g%l%c%cc14f%dwt_prod10c_gain(c) = 0._rk8
        clm3%g%l%c%cc14f%prod10c_loss(c) = 0._rk8
        clm3%g%l%c%cc14f%dwt_prod100c_gain(c) = 0._rk8
        clm3%g%l%c%cc14f%prod100c_loss(c) = 0._rk8
        do j = 1, nlevdecomp_full
          clm3%g%l%c%cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._rk8
          clm3%g%l%c%cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._rk8
          clm3%g%l%c%cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._rk8
          clm3%g%l%c%cc14f%dwt_livecrootc_to_cwdc(c,j) = 0._rk8
          clm3%g%l%c%cc14f%dwt_deadcrootc_to_cwdc(c,j) = 0._rk8
        end do
        clm3%g%l%c%cc14f%dwt_closs(c) = 0._rk8
      end if

      clm3%g%l%c%cnf%dwt_seedn_to_leaf(c) = 0._rk8
      clm3%g%l%c%cnf%dwt_seedn_to_deadstem(c) = 0._rk8
      clm3%g%l%c%cnf%dwt_conv_nflux(c) = 0._rk8
      clm3%g%l%c%cnf%dwt_prod10n_gain(c) = 0._rk8
      clm3%g%l%c%cnf%prod10n_loss(c) = 0._rk8
      clm3%g%l%c%cnf%dwt_prod100n_gain(c) = 0._rk8
      clm3%g%l%c%cnf%prod100n_loss(c) = 0._rk8
      do j = 1, nlevdecomp_full
        clm3%g%l%c%cnf%dwt_frootn_to_litr_met_n(c,j) = 0._rk8
        clm3%g%l%c%cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._rk8
        clm3%g%l%c%cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._rk8
        clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn(c,j) = 0._rk8
        clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn(c,j) = 0._rk8
      end do
      clm3%g%l%c%cnf%dwt_nloss(c) = 0._rk8
    end do

    ! deallocate special landunit filters
    deallocate(specialc)
    deallocate(specialp)

#endif
  end subroutine CNiniSpecial

end module mod_clm_cninispecial
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
