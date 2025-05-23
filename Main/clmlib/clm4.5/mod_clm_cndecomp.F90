module mod_clm_cndecomp
#ifdef CN
  !
  ! Module holding routines used in litter and soil decomposition model
  ! for coupled carbon-nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_varcon, only : dzsoi_decomp
#ifndef CENTURY_DECOMP
  use mod_clm_cndecompcascadebgc, only : decomp_rate_constants
#else
  use mod_clm_cndecompcascadecentury, only : decomp_rate_constants
#endif
#ifdef NITRIF_DENITRIF
  use mod_clm_cnnitrifdenitrif, only : nitrif_denitrif
#endif
  use mod_clm_cnverticalprofile, only : decomp_vertprofiles

  implicit none

  save

  private

  public :: CNDecompAlloc

  contains

  subroutine CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                           num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_cnallocation, only : CNAllocation
    use mod_clm_varpar, only : nlevdecomp, &
      ndecomp_cascade_transitions, ndecomp_pools
    use mod_clm_subgridave, only : p2c
    implicit none
    integer(ik4), intent(in) :: lbp, ubp    ! pft-index bounds
    integer(ik4), intent(in) :: lbc, ubc    ! column-index bounds
    integer(ik4), intent(in) :: num_soilc   ! number of soil columns in filter
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)
    integer(ik4), intent(in) :: num_soilp   ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts

    ! all c pools involved in decomposition
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous :: decomp_cpools_vr(:,:,:)
    ! all n pools involved in decomposition
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous :: decomp_npools_vr(:,:,:)

    ! index into landunit level quantities
    integer(ik4), pointer, contiguous :: clandunit(:)
    ! landunit type
    integer(ik4), pointer, contiguous :: itypelun(:)
    ! pft level
    ! fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer, contiguous :: rootfr(:,:)

    ! fraction of potential immobilization (no units)
    real(rk8), pointer, contiguous :: fpi_vr(:,:)
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer, contiguous :: decomp_cascade_hr_vr(:,:,:)
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer, contiguous :: decomp_cascade_ctransfer_vr(:,:,:)
    ! vert-res transfer of N from donor to receiver pool along
    ! decomp. cascade (gN/m3/s)
    real(rk8), pointer, contiguous :: decomp_cascade_ntransfer_vr(:,:,:)
    ! vert-res mineral N flux for transition along
    ! decomposition cascade (gN/m3/s)
    real(rk8), pointer, contiguous :: decomp_cascade_sminn_flux_vr(:,:,:)

    real(rk8), pointer, contiguous :: potential_immob_vr(:,:)
#ifndef NITRIF_DENITRIF
    real(rk8), pointer, contiguous :: sminn_to_denit_decomp_cascade_vr(:,:,:)
    real(rk8), pointer, contiguous :: sminn_to_denit_excess_vr(:,:)
#endif
    real(rk8), pointer, contiguous :: gross_nmin_vr(:,:)
    real(rk8), pointer, contiguous :: net_nmin_vr(:,:)
    ! gross rate of N mineralization (gN/m2/s)
    real(rk8), pointer, contiguous :: gross_nmin(:)
    ! net rate of N mineralization (gN/m2/s)
    real(rk8), pointer, contiguous :: net_nmin(:)
    ! For methane code
#ifdef LCH4
    ! fraction of potential SOM + LITTER heterotrophic respiration
    real(rk8), pointer, contiguous :: fphr(:,:)
    ! fraction by which decomposition is limited by moisture availability
    real(rk8), pointer, contiguous :: w_scalar(:,:)
#endif
    ! rate constant for decomposition (1./sec)
    real(rk8), pointer, contiguous :: decomp_k(:,:,:)
    ! respired fraction in decomposition step (frac)
    real(rk8), pointer, contiguous :: rf_decomp_cascade(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer, contiguous :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer, contiguous :: cascade_receiver_pool(:)
    ! what fraction of C leaving a given pool passes through
    ! a given transition (frac)
    real(rk8), pointer, contiguous :: pathfrac_decomp_cascade(:,:,:)
    ! TRUE => pool has fixed C:N ratio
    logical,  pointer, contiguous :: floating_cn_ratio_decomp_pools(:)

    integer(ik4) :: c,j,k,l    !indices
    integer(ik4) :: fc           !lake filter column index
    !potential C loss from one pool to another
    real(rk8):: p_decomp_cpool_loss(lbc:ubc,1:nlevdecomp, &
                                            1:ndecomp_cascade_transitions)
    !potential mineral N flux, from one pool to another
    real(rk8):: pmnf_decomp_cascade(lbc:ubc,1:nlevdecomp, &
                                            1:ndecomp_cascade_transitions)
    !potential N immobilization
    real(rk8):: immob(lbc:ubc,1:nlevdecomp)
    real(rk8):: ratio        !temporary variable
    real(rk8):: dnp          !denitrification proportion
    real(rk8):: cn_decomp_pools(lbc:ubc,1:nlevdecomp,1:ndecomp_pools)
    ! c:n ratio for initialization of pools
    real(rk8), pointer, contiguous :: initial_cn_ratio(:)
    integer(ik4), parameter :: i_atm = 0
    ! maximum annual depth of thaw
    integer(ik4), pointer, contiguous :: altmax_indx(:)
    ! prior year maximum annual depth of thaw
    integer(ik4), pointer, contiguous :: altmax_lastyear_indx(:)

    ! For methane code
#ifndef NITRIF_DENITRIF
    real(rk8):: phr_vr(lbc:ubc,1:nlevdecomp) !potential HR (gC/m3/s)
#else
    real(rk8), pointer, contiguous :: phr_vr(:,:)        !potential HR (gC/m3/s)
#endif
#ifdef LCH4
    real(rk8):: hrsum(lbc:ubc,1:nlevdecomp)  !sum of HR (gC/m2/s)
#endif

    decomp_cpools_vr             => clm3%g%l%c%ccs%decomp_cpools_vr
    decomp_cascade_hr_vr         => clm3%g%l%c%ccf%decomp_cascade_hr_vr
    decomp_cascade_ctransfer_vr  => clm3%g%l%c%ccf%decomp_cascade_ctransfer_vr
    decomp_npools_vr             => clm3%g%l%c%cns%decomp_npools_vr
    decomp_cascade_ntransfer_vr  => clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr
    decomp_cascade_sminn_flux_vr => clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr
    fpi_vr                => clm3%g%l%c%cps%fpi_vr
    potential_immob_vr    => clm3%g%l%c%cnf%potential_immob_vr

    decomp_k                   => clm3%g%l%c%ccf%decomp_k
    rf_decomp_cascade          => clm3%g%l%c%cps%rf_decomp_cascade
    cascade_donor_pool         => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool      => decomp_cascade_con%cascade_receiver_pool
    pathfrac_decomp_cascade    => clm3%g%l%c%cps%pathfrac_decomp_cascade
    floating_cn_ratio_decomp_pools => &
      decomp_cascade_con%floating_cn_ratio_decomp_pools
    initial_cn_ratio       => decomp_cascade_con%initial_cn_ratio
    altmax_indx            => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx   => clm3%g%l%c%cps%altmax_lastyear_indx

#ifndef NITRIF_DENITRIF
    sminn_to_denit_decomp_cascade_vr => &
      clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr
    sminn_to_denit_excess_vr => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
#else
    phr_vr           => clm3%g%l%c%ccf%phr_vr
#endif
    gross_nmin_vr    => clm3%g%l%c%cnf%gross_nmin_vr
    net_nmin_vr      => clm3%g%l%c%cnf%net_nmin_vr
    gross_nmin       => clm3%g%l%c%cnf%gross_nmin
    net_nmin         => clm3%g%l%c%cnf%net_nmin
    ! For methane code
#ifdef LCH4
    fphr             => clm3%g%l%c%cch4%fphr
    w_scalar         => clm3%g%l%c%ccf%w_scalar
#endif

    rootfr           => clm3%g%l%c%p%pps%rootfr
    clandunit        => clm3%g%l%c%landunit
    itypelun         => clm3%g%l%itype

    call decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)

    ! set initial values for potential C and N fluxes
    p_decomp_cpool_loss(:,:,:) = 0._rk8
    pmnf_decomp_cascade(:,:,:) = 0._rk8

    ! column loop to calculate potential decomp rates and total immobilization
    ! demand.

    ! calculate c:n ratios of applicable pools
    do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
        do j = 1, nlevdecomp
          do fc = 1, num_soilc
            c = filter_soilc(fc)
            if ( decomp_npools_vr(c,j,l) > 0._rk8 ) then
              cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / &
                decomp_npools_vr(c,j,l)
            end if
          end do
        end do
      else
        do j = 1,nlevdecomp
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
          end do
        end do
      end if
    end do


    ! calculate the non-nitrogen-limited fluxes
    ! these fluxes include the  "/ dt" term to put them on a
    ! per second basis, since the rate constants have been
    ! calculated on a per timestep basis.

    do k = 1, ndecomp_cascade_transitions
      do j = 1, nlevdecomp
        do fc = 1, num_soilc
          c = filter_soilc(fc)

          if ( decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 .and. &
               decomp_k(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
           p_decomp_cpool_loss(c,j,k) = &
             decomp_cpools_vr(c,j,cascade_donor_pool(k)) * &
             decomp_k(c,j,cascade_donor_pool(k)) * &
             pathfrac_decomp_cascade(c,j,k)
           if ( .not. &
             floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then
             !! not transition of cwd to litter

             if (cascade_receiver_pool(k) /= i_atm ) then
               ! not 100% respiration
               ratio = 0._rk8
               if ( decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
                 ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k)) / &
                   cn_decomp_pools(c,j,cascade_donor_pool(k))
               endif
               pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * &
                 (1.0_rk8 - rf_decomp_cascade(c,j,k) - ratio) / &
                  cn_decomp_pools(c,j,cascade_receiver_pool(k)) )
             else   ! 100% respiration
               pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / &
                 cn_decomp_pools(c,j,cascade_donor_pool(k))
             end if
           else   ! CWD -> litter
             pmnf_decomp_cascade(c,j,k) = 0._rk8
           end if
         end if
       end do
     end do
   end do

   ! Sum up all the potential immobilization fluxes (positive pmnf flux)
   ! and all the mineralization fluxes (negative pmnf flux)
   do j = 1, nlevdecomp
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       immob(c,j) = 0._rk8
     end do
   end do
   do k = 1, ndecomp_cascade_transitions
     do j = 1, nlevdecomp
       do fc = 1, num_soilc
         c = filter_soilc(fc)
         if (pmnf_decomp_cascade(c,j,k) > 0._rk8) then
           immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
         else
           gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
         end if
       end do
     end do
   end do

   do j = 1, nlevdecomp
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       potential_immob_vr(c,j) = immob(c,j)
     end do
   end do

   ! Add up potential hr for methane calculations
   do j = 1, nlevdecomp
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       phr_vr(c,j) = 0._rk8
     end do
   end do
   do k = 1, ndecomp_cascade_transitions
     do j = 1, nlevdecomp
       do fc = 1, num_soilc
         c = filter_soilc(fc)
         phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * &
           p_decomp_cpool_loss(c,j,k)
       end do
     end do
   end do

   call decomp_vertprofiles(lbp, ubp, lbc,ubc,num_soilc, &
     filter_soilc,num_soilp,filter_soilp)

#ifdef NITRIF_DENITRIF
   ! calculate nitrification and denitrification rates
   call nitrif_denitrif(lbc, ubc, num_soilc, filter_soilc)
#endif

   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.

   call CNAllocation(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp, &
        filter_soilp)

   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N

   dnp = 0.01_rk8

   ! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
     if ( floating_cn_ratio_decomp_pools(l) ) then
       do j = 1, nlevdecomp
         do fc = 1, num_soilc
           c = filter_soilc(fc)
           if ( decomp_npools_vr(c,j,l) > 0._rk8 ) then
             cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / &
               decomp_npools_vr(c,j,l)
           end if
         end do
       end do
     else
       do j = 1, nlevdecomp
         do fc = 1, num_soilc
           c = filter_soilc(fc)
           cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
         end do
       end do
     end if
   end do

   ! upon return from CNAllocation, the fraction of potential immobilization
   ! has been set (cps%fpi_vr). now finish the decomp calculations.
   ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
   ! Also calculate denitrification losses as a simple proportion
   ! of mineralization flux.

   do k = 1, ndecomp_cascade_transitions
     do j = 1, nlevdecomp
       do fc = 1, num_soilc
         c = filter_soilc(fc)

         if ( decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
           if ( pmnf_decomp_cascade(c,j,k) > 0._rk8 ) then
             p_decomp_cpool_loss(c,j,k) = &
               p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
             pmnf_decomp_cascade(c,j,k) = &
               pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
#ifndef NITRIF_DENITRIF
             sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._rk8
           else
             sminn_to_denit_decomp_cascade_vr(c,j,k) = -dnp * &
               pmnf_decomp_cascade(c,j,k)
#endif
           end if
           decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * &
             p_decomp_cpool_loss(c,j,k)
           decomp_cascade_ctransfer_vr(c,j,k) = (1._rk8 - &
             rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
           if ( decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 .and. &
                cascade_receiver_pool(k) /= i_atm ) then
             decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / &
               cn_decomp_pools(c,j,cascade_donor_pool(k))
           else
             decomp_cascade_ntransfer_vr(c,j,k) = 0._rk8
           endif
           if ( cascade_receiver_pool(k) /= 0 ) then
             decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
           else  ! keep sign convention negative for terminal pools
             decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
           endif
           net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
         else
           decomp_cascade_ntransfer_vr(c,j,k) = 0._rk8
#ifndef NITRIF_DENITRIF
           sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._rk8
#endif
           decomp_cascade_sminn_flux_vr(c,j,k) = 0._rk8
         end if
       end do
     end do
   end do

#ifdef LCH4
   ! Calculate total fraction of potential HR, for methane code
   do j = 1,nlevdecomp
     do fc = 1,num_soilc
       c = filter_soilc(fc)
       hrsum(c,j) = 0._rk8
     end do
   end do
   do k = 1, ndecomp_cascade_transitions
     do j = 1, nlevdecomp
       do fc = 1, num_soilc
         c = filter_soilc(fc)
         hrsum(c,j) = hrsum(c,j) + &
           rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
       end do
     end do
   end do

   ! Nitrogen limitation / (low)-moisture limitation
   do j = 1, nlevdecomp
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (phr_vr(c,j) > 0._rk8) then
         fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
         ! Prevent overflow errors for 0 respiration
         fphr(c,j) = max(fphr(c,j), 0.01_rk8)
       else
         fphr(c,j) = 1._rk8
       end if
     end do
   end do
#endif

   ! vertically integrate net and gross mineralization fluxes
   ! for diagnostic output
   do j = 1, nlevdecomp
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
       gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)
     end do
   end do
 end subroutine CNDecompAlloc

#endif

end module mod_clm_cndecomp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
