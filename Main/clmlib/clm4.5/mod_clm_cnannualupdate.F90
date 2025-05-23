module mod_clm_cnannualupdate
#ifdef CN
  !
  ! Module for updating annual summation variables
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams, only : dtsrf
  use mod_clm_time_manager, only : is_end_curr_year

  implicit none

  save

  private

  public:: CNAnnualUpdate

  contains
  !
  ! On the radiation time step, update annual summation variables
  !
  subroutine CNAnnualUpdate(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                            num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_varcon , only: secspday
    use mod_clm_subgridave, only: p2c
    implicit none
    integer(ik4), intent(in) :: lbc, ubc   ! column bounds
    integer(ik4), intent(in) :: lbp, ubp   ! pft bounds
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)
    ! number of soil pfts in filter
    integer(ik4), intent(in) :: num_soilp
    ! filter for soil pfts
    integer(ik4), intent(in) :: filter_soilp(ubp-lbp+1)

    integer(ik4), pointer, contiguous :: pcolumn(:) ! index into column level quantities

    ! seconds since last annual accumulator turnover
    real(rk8), pointer, contiguous :: annsum_counter(:)
    ! temporary annual sum of potential GPP
    real(rk8), pointer, contiguous :: tempsum_potential_gpp(:)
    ! annual sum of potential GPP
    real(rk8), pointer, contiguous :: annsum_potential_gpp(:)
    ! temporary annual max of retranslocated N pool (gN/m2)
    real(rk8), pointer, contiguous :: tempmax_retransn(:)
    ! annual max of retranslocated N pool (gN/m2)
    real(rk8), pointer, contiguous :: annmax_retransn(:)
    ! temporary average 2m air temperature (K)
    real(rk8), pointer, contiguous :: tempavg_t2m(:)
    ! annual average 2m air temperature (K)
    real(rk8), pointer, contiguous :: annavg_t2m(:)
    ! temporary sum NPP (gC/m2/yr)
    real(rk8), pointer, contiguous :: tempsum_npp(:)
    ! annual sum NPP (gC/m2/yr)
    real(rk8), pointer, contiguous :: annsum_npp(:)
    ! column annual sum NPP (gC/m2/yr)
    real(rk8), pointer, contiguous :: cannsum_npp(:)
    !annual average of 2m air temperature, averaged from pft-level (K)
    real(rk8), pointer, contiguous :: cannavg_t2m(:)
#if (defined CNDV)
    ! temporary sum litfall (gC/m2/yr)
    real(rk8), pointer, contiguous :: tempsum_litfall(:)
    ! annual sum litfall (gC/m2/yr)
    real(rk8), pointer, contiguous :: annsum_litfall(:)
#endif

    integer(ik4) :: c,p     ! indices
    integer(ik4) :: fp,fc   ! lake filter indices
    logical :: eoy

    ! assign local pointers to derived type arrays
    annsum_counter        => clm3%g%l%c%cps%annsum_counter
    tempsum_potential_gpp => clm3%g%l%c%p%pepv%tempsum_potential_gpp
    annsum_potential_gpp  => clm3%g%l%c%p%pepv%annsum_potential_gpp
    tempmax_retransn      => clm3%g%l%c%p%pepv%tempmax_retransn
    annmax_retransn       => clm3%g%l%c%p%pepv%annmax_retransn
    tempavg_t2m           => clm3%g%l%c%p%pepv%tempavg_t2m
    annavg_t2m            => clm3%g%l%c%p%pepv%annavg_t2m
    tempsum_npp           => clm3%g%l%c%p%pepv%tempsum_npp
    annsum_npp            => clm3%g%l%c%p%pepv%annsum_npp
    cannsum_npp           => clm3%g%l%c%cps%cannsum_npp
    cannavg_t2m           => clm3%g%l%c%cps%cannavg_t2m
#if (defined CNDV)
    tempsum_litfall       => clm3%g%l%c%p%pepv%tempsum_litfall
    annsum_litfall        => clm3%g%l%c%p%pepv%annsum_litfall
#endif
    pcolumn               => clm3%g%l%c%p%column

    eoy = is_end_curr_year( )

    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      annsum_counter(c) = annsum_counter(c) + dtsrf
    end do

    if ( num_soilc > 0 ) then

      if ( eoy ) then
        ! pft loop
        do fp = 1, num_soilp
          p = filter_soilp(fp)
          ! update annual plant ndemand accumulator
          annsum_potential_gpp(p)  = tempsum_potential_gpp(p)
          tempsum_potential_gpp(p) = 0._rk8

          ! update annual total N retranslocation accumulator
          annmax_retransn(p)  = tempmax_retransn(p)
          tempmax_retransn(p) = 0._rk8

          ! update annual average 2m air temperature accumulator
          annavg_t2m(p)  = tempavg_t2m(p)
          tempavg_t2m(p) = 0._rk8

          ! update annual NPP accumulator, convert to annual total
          annsum_npp(p) = tempsum_npp(p) * dtsrf
          tempsum_npp(p) = 0._rk8

#if (defined CNDV)
          ! update annual litfall accumulator, convert to annual total
          annsum_litfall(p) = tempsum_litfall(p) * dtsrf
          tempsum_litfall(p) = 0._rk8
#endif
        end do

        ! use p2c routine to get selected column-average
        ! pft-level fluxes and states
        call p2c(num_soilc, filter_soilc, annsum_npp, cannsum_npp)
        call p2c(num_soilc, filter_soilc, annavg_t2m, cannavg_t2m)
      end if
    end if

    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if ( eoy ) annsum_counter(c) = 0._rk8
    end do
  end subroutine CNAnnualUpdate

#endif

end module mod_clm_cnannualupdate

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
