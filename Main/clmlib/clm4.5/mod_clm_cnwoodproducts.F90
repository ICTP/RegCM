module mod_clm_cnwoodproducts
#ifdef CN
  !
  ! Calculate loss fluxes from wood products pools, and update product
  ! pool state variables
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams, only : dtsrf
  use mod_clm_decomp, only : get_proc_bounds
  use mod_clm_varcon, only : istsoil

  implicit none

  save

  private

  public:: CNWoodProducts

  contains
  !
  ! Update all loss fluxes from wood product pools, and update product
  ! pool state variables for both loss and gain terms.
  ! Gain terms are calculated in pftdyn_cnbal() for gains associated
  ! with changes in landcover, and in CNHarvest(), for gains associated
  ! with wood harvest.
  !
  subroutine CNWoodProducts(num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varctl, only : use_c13, use_c14
    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns

    integer(ik4) :: fc  ! lake filter indices
    integer(ik4) :: c   ! indices
    real(rk8):: dt      ! time step (seconds)
    type(column_type), pointer :: cptr  ! pointer to column derived subtype
    real(rk8) :: kprod10   ! decay constant for 10-year product pool
    real(rk8) :: kprod100  ! decay constant for 100-year product pool

    cptr => clm3%g%l%c

    ! calculate column-level losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial
    ! state over 10 and 100 years, respectively, using a discrete-time
    ! fractional decay algorithm.
    kprod10 = 7.2e-9_rk8
    kprod100 = 7.2e-10_rk8

    do fc = 1, num_soilc
      c = filter_soilc(fc)

      ! calculate fluxes (1/sec)
      cptr%ccf%prod10c_loss(c)    = cptr%ccs%prod10c(c)    * kprod10
      cptr%ccf%prod100c_loss(c)   = cptr%ccs%prod100c(c)   * kprod100

      if ( use_c13 ) then
        cptr%cc13f%prod10c_loss(c)  = cptr%cc13s%prod10c(c)  * kprod10
        cptr%cc13f%prod100c_loss(c) = cptr%cc13s%prod100c(c) * kprod100
      end if

      if ( use_c14 ) then
        cptr%cc14f%prod10c_loss(c)  = cptr%cc14s%prod10c(c)  * kprod10
        cptr%cc14f%prod100c_loss(c) = cptr%cc14s%prod100c(c) * kprod100
      end if

      cptr%cnf%prod10n_loss(c)    = cptr%cns%prod10n(c)    * kprod10
      cptr%cnf%prod100n_loss(c)   = cptr%cns%prod100n(c)   * kprod100
    end do

    ! set time steps
    dt = dtsrf

    ! update wood product state variables
    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)

      ! column-level fluxes

      ! fluxes into wood product pools, from landcover change
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    + &
              cptr%ccf%dwt_prod10c_gain(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   + &
              cptr%ccf%dwt_prod100c_gain(c)*dt

      if ( use_c13 ) then
        cptr%cc13s%prod10c(c)  = cptr%cc13s%prod10c(c)  + &
                cptr%cc13f%dwt_prod10c_gain(c)*dt
        cptr%cc13s%prod100c(c) = cptr%cc13s%prod100c(c) + &
                cptr%cc13f%dwt_prod100c_gain(c)*dt
      end if

      if ( use_c14 ) then
        cptr%cc14s%prod10c(c)  = cptr%cc14s%prod10c(c)  + &
                cptr%cc14f%dwt_prod10c_gain(c)*dt
        cptr%cc14s%prod100c(c) = cptr%cc14s%prod100c(c) + &
                cptr%cc14f%dwt_prod100c_gain(c)*dt
      end if

      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    + &
              cptr%cnf%dwt_prod10n_gain(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   + &
              cptr%cnf%dwt_prod100n_gain(c)*dt

      ! fluxes into wood product pools, from harvest
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    + &
              cptr%ccf%hrv_deadstemc_to_prod10c(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   + &
              cptr%ccf%hrv_deadstemc_to_prod100c(c)*dt

      if ( use_c13 ) then
        cptr%cc13s%prod10c(c)  = cptr%cc13s%prod10c(c)  + &
                cptr%cc13f%hrv_deadstemc_to_prod10c(c)*dt
        cptr%cc13s%prod100c(c) = cptr%cc13s%prod100c(c) + &
                cptr%cc13f%hrv_deadstemc_to_prod100c(c)*dt
      end if

      if ( use_c14 ) then
        cptr%cc14s%prod10c(c)  = cptr%cc14s%prod10c(c)  + &
                cptr%cc14f%hrv_deadstemc_to_prod10c(c)*dt
        cptr%cc14s%prod100c(c) = cptr%cc14s%prod100c(c) + &
                cptr%cc14f%hrv_deadstemc_to_prod100c(c)*dt
      end if

      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    + &
              cptr%cnf%hrv_deadstemn_to_prod10n(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   + &
              cptr%cnf%hrv_deadstemn_to_prod100n(c)*dt

      ! fluxes out of wood product pools, from decomposition
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    - &
              cptr%ccf%prod10c_loss(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   - &
              cptr%ccf%prod100c_loss(c)*dt

      if ( use_c13 ) then
        cptr%cc13s%prod10c(c)  = cptr%cc13s%prod10c(c)  - &
                cptr%cc13f%prod10c_loss(c)*dt
        cptr%cc13s%prod100c(c) = cptr%cc13s%prod100c(c) - &
                cptr%cc13f%prod100c_loss(c)*dt
      end if

      if ( use_c14 ) then
        cptr%cc14s%prod10c(c)  = cptr%cc14s%prod10c(c)  - &
                cptr%cc14f%prod10c_loss(c)*dt
        cptr%cc14s%prod100c(c) = cptr%cc14s%prod100c(c) - &
                cptr%cc14f%prod100c_loss(c)*dt
      end if

      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    - &
              cptr%cnf%prod10n_loss(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   - &
              cptr%cnf%prod100n_loss(c)*dt

    end do ! end of column loop

  end subroutine CNWoodProducts

#endif

end module mod_clm_cnwoodproducts
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
