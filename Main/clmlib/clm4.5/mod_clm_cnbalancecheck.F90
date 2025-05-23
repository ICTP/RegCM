module mod_clm_cnbalancecheck
#ifdef CN

  !
  ! Module for carbon mass balance checking.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams, only : dtsrf
  use mod_mpmessage
  use mod_stdio

  implicit none

  save

  private

  public :: BeginCBalance
  public :: BeginNBalance
  public :: CBalanceCheck
  public :: NBalanceCheck

  contains
  !
  ! On the radiation time step, calculate the beginning carbon balance
  ! for mass conservation checks.
  !
  subroutine BeginCBalance(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    implicit none
    integer(ik4), intent(in) :: lbc, ubc  ! column bounds
    ! number of soil columns filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)
    ! (gC/m2) total column carbon, incl veg and cpool
    real(rk8), pointer, contiguous :: totcolc(:)

    ! carbon mass, beginning of time step (gC/m**2)
    real(rk8), pointer, contiguous :: col_begcb(:)
    integer(ik4) :: c     ! indices
    integer(ik4) :: fc   ! lake filter indices

    ! assign local pointers at the column level
    col_begcb => clm3%g%l%c%ccbal%begcb
    totcolc   => clm3%g%l%c%ccs%totcolc

    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      ! calculate beginning column-level carbon balance,
      ! for mass conservation check
      col_begcb(c) = totcolc(c)
    end do ! end of columns loop

  end subroutine BeginCBalance
  !
  ! On the radiation time step, calculate the beginning nitrogen
  ! balance for mass conservation checks.
  !
  subroutine BeginNBalance(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    implicit none
    integer(ik4), intent(in) :: lbc, ubc  ! column bounds
    ! number of soil columns filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)

    ! (gN/m2) total column nitrogen, incl veg
    real(rk8), pointer, contiguous :: totcoln(:)

    ! nitrogen mass, beginning of time step (gN/m**2)
    real(rk8), pointer, contiguous :: col_begnb(:)

    integer(ik4) :: c   ! indices
    integer(ik4) :: fc  ! lake filter indices

    ! assign local pointers at the column level
    col_begnb => clm3%g%l%c%cnbal%begnb
    totcoln   => clm3%g%l%c%cns%totcoln

    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      ! calculate beginning column-level nitrogen balance,
      ! for mass conservation check
      col_begnb(c) = totcoln(c)
    end do ! end of columns loop
  end subroutine BeginNBalance
  !
  ! On the radiation time step, perform carbon mass conservation
  ! check for column and pft
  !
  subroutine CBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    implicit none
    integer(ik4), intent(in) :: lbc, ubc  ! column bounds
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)

    ! (gC/m2) total column carbon, incl veg and cpool
    real(rk8), pointer, contiguous :: totcolc(:)
    ! (gC/m2/s) gross primary production
    real(rk8), pointer, contiguous :: gpp(:)
    ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
    real(rk8), pointer, contiguous :: er(:)
    ! (gC/m2/s) total column-level fire C loss
    real(rk8), pointer, contiguous :: col_fire_closs(:)
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous :: col_hrv_xsmrpool_to_atm(:)
    ! (gC/m2/s) total carbon loss from product pools and conversion
    real(rk8), pointer, contiguous :: dwt_closs(:)
    ! (gC/m2/s) total wood product carbon loss
    real(rk8), pointer, contiguous :: product_closs(:)
    ! total SOM C loss from vertical transport (gC/m^2/s)
    real(rk8), pointer, contiguous :: som_c_leached(:)
    ! grid pointer
    integer(ik4), pointer, contiguous :: gridcell(:)

    ! (gC/m2/s) total column-level carbon inputs (for balance check)
    real(rk8), pointer, contiguous :: col_cinputs(:)
    ! (gC/m2/s) total column-level carbon outputs (for balance check)
    real(rk8), pointer, contiguous :: col_coutputs(:)
    ! carbon mass, beginning of time step (gC/m**2)
    real(rk8), pointer, contiguous :: col_begcb(:)
    ! carbon mass, end of time step (gC/m**2)
    real(rk8), pointer, contiguous :: col_endcb(:)
    ! carbon balance error for the timestep (gC/m**2)
    real(rk8), pointer, contiguous :: col_errcb(:)

    integer(ik4) :: c,err_index,g  ! indices
    integer(ik4) :: fc   ! lake filter indices
    logical :: err_found ! error flag
    real(rk8):: dt       ! radiation time step (seconds)

    ! assign local pointers to column-level arrays
    totcolc                  => clm3%g%l%c%ccs%totcolc
    gpp                      => clm3%g%l%c%ccf%pcf_a%gpp
    er                       => clm3%g%l%c%ccf%er
    col_fire_closs           => clm3%g%l%c%ccf%col_fire_closs
    col_hrv_xsmrpool_to_atm  => clm3%g%l%c%ccf%pcf_a%hrv_xsmrpool_to_atm
    dwt_closs                => clm3%g%l%c%ccf%dwt_closs
    product_closs            => clm3%g%l%c%ccf%product_closs
    gridcell                 => clm3%g%l%c%gridcell

    col_cinputs              => clm3%g%l%c%ccf%col_cinputs
    col_coutputs             => clm3%g%l%c%ccf%col_coutputs
    col_begcb                => clm3%g%l%c%ccbal%begcb
    col_endcb                => clm3%g%l%c%ccbal%endcb
    col_errcb                => clm3%g%l%c%ccbal%errcb
    som_c_leached            => clm3%g%l%c%ccf%som_c_leached

    ! set time steps
    dt = dtsrf

    err_found = .false.
    ! column loop
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      ! calculate the total column-level carbon storage,
      ! for mass conservation check
      col_endcb(c) = totcolc(c)
      ! calculate total column-level inputs
      col_cinputs(c) = gpp(c)
      ! calculate total column-level outputs
      ! er = ar + hr, col_fire_closs includes pft-level fire losses
      col_coutputs(c) = er(c) + col_fire_closs(c) + dwt_closs(c) + &
        product_closs(c) + col_hrv_xsmrpool_to_atm(c)
      ! subtract leaching flux
      col_coutputs(c) = col_coutputs(c) - som_c_leached(c)
      ! calculate the total column-level carbon balance error for this time step
      col_errcb(c) = (col_cinputs(c) - col_coutputs(c))*dt - &
           (col_endcb(c) - col_begcb(c))
      ! check for significant errors

      if ( abs(col_errcb(c)) > 1.0e-3_rk8 .and. 1==2) then
        err_found = .true.
        err_index = c
      end if
    end do ! end of columns loop
    if (err_found) then
      c = err_index
      g = gridcell(c)
      write(stderr,*) 'Latdeg, Londeg =', clm3%g%latdeg(g), clm3%g%londeg(g)
      write(stderr,*) 'Column = ', c
      write(stderr,*) 'column cbalance error = ', col_errcb(c)
      write(stderr,*) 'delta store = ',col_endcb(c)-col_begcb(c)
      write(stderr,*) 'begcb       = ',col_begcb(c)
      write(stderr,*) 'endcb       = ',col_endcb(c)
      write(stderr,*) 'inputs      = ',col_cinputs(c)*dt
      write(stderr,*) 'outputs     = ',col_coutputs(c)*dt
      write(stderr,*) 'ER          = ',er(c)*dt
      write(stderr,*) 'Fire Loss   = ',col_fire_closs(c)*dt
      write(stderr,*) 'DWT Loss    = ',dwt_closs(c)*dt
      write(stderr,*) 'Product CL  = ',product_closs(c)*dt
      write(stderr,*) 'To Atm      = ',col_hrv_xsmrpool_to_atm(c)*dt
      write(stderr,*) 'Leached     = ',som_c_leached(c)*dt
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
  end subroutine CBalanceCheck
  !
  ! On the radiation time step, perform nitrogen mass conservation check
  ! for column and pft
  !
  subroutine NBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_surfrd, only : crop_prog
    implicit none
    integer(ik4), intent(in) :: lbc, ubc ! column bounds
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)

    ! (gN/m2) total column nitrogen, incl veg
    real(rk8), pointer, contiguous :: totcoln(:)
    ! atmospheric N deposition to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous :: ndep_to_sminn(:)
    ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous :: nfix_to_sminn(:)
    real(rk8), pointer, contiguous :: fert_to_sminn(:)
    real(rk8), pointer, contiguous :: soyfixn_to_sminn(:)
    ! supplemental N supply (gN/m2/s)
    real(rk8), pointer, contiguous :: supplement_to_sminn(:)
    ! total rate of denitrification (gN/m2/s)
    real(rk8), pointer, contiguous :: denit(:)
#ifndef NITRIF_DENITRIF
    ! soil mineral N pool loss to leaching (gN/m2/s)
    real(rk8), pointer, contiguous :: sminn_leached(:)
#else
    ! soil mineral NO3 pool loss to leaching (gN/m2/s)
    real(rk8), pointer, contiguous :: smin_no3_leached(:)
    ! soil mineral NO3 pool loss to runoff (gN/m2/s)
    real(rk8), pointer, contiguous :: smin_no3_runoff(:)
    ! flux of N2o from nitrification [gN/m^2/s]
    real(rk8), pointer, contiguous :: f_n2o_nit(:)
#endif
    ! total column-level fire N loss (gN/m2/s)
    real(rk8), pointer, contiguous :: col_fire_nloss(:)
    ! (gN/m2/s) total nitrogen loss from product pools and conversion
    real(rk8), pointer, contiguous :: dwt_nloss(:)
    ! (gN/m2/s) total wood product nitrogen loss
    real(rk8), pointer, contiguous :: product_nloss(:)
    ! total SOM N loss from vertical transport
    real(rk8), pointer, contiguous :: som_n_leached(:)

    ! column-level N inputs (gN/m2/s)
    real(rk8), pointer, contiguous :: col_ninputs(:)
    ! column-level N outputs (gN/m2/s)
    real(rk8), pointer, contiguous :: col_noutputs(:)
    ! nitrogen mass, beginning of time step (gN/m**2)
    real(rk8), pointer, contiguous :: col_begnb(:)
    ! nitrogen mass, end of time step (gN/m**2)
    real(rk8), pointer, contiguous :: col_endnb(:)
    ! nitrogen balance error for the timestep (gN/m**2)
    real(rk8), pointer, contiguous :: col_errnb(:)
    ! grid pointer
    integer(ik4), pointer, contiguous :: gridcell(:)

    integer(ik4) :: c,err_index,g    ! indices
    integer(ik4) :: fc   ! lake filter indices
    logical :: err_found ! error flag
    real(rk8):: dt       ! radiation time step (seconds)

    ! assign local pointers to column-level arrays

    totcoln               => clm3%g%l%c%cns%totcoln
    ndep_to_sminn         => clm3%g%l%c%cnf%ndep_to_sminn
    nfix_to_sminn         => clm3%g%l%c%cnf%nfix_to_sminn
    fert_to_sminn         => clm3%g%l%c%cnf%fert_to_sminn
    soyfixn_to_sminn      => clm3%g%l%c%cnf%soyfixn_to_sminn
    supplement_to_sminn   => clm3%g%l%c%cnf%supplement_to_sminn
    denit                 => clm3%g%l%c%cnf%denit
#ifndef NITRIF_DENITRIF
    sminn_leached         => clm3%g%l%c%cnf%sminn_leached
#else
    smin_no3_leached      => clm3%g%l%c%cnf%smin_no3_leached
    smin_no3_runoff       => clm3%g%l%c%cnf%smin_no3_runoff
    f_n2o_nit             => clm3%g%l%c%cnf%f_n2o_nit
#endif
    col_fire_nloss        => clm3%g%l%c%cnf%col_fire_nloss
    dwt_nloss             => clm3%g%l%c%cnf%dwt_nloss
    product_nloss         => clm3%g%l%c%cnf%product_nloss
    som_n_leached         => clm3%g%l%c%cnf%som_n_leached
    gridcell              => clm3%g%l%c%gridcell

    col_ninputs           => clm3%g%l%c%cnf%col_ninputs
    col_noutputs          => clm3%g%l%c%cnf%col_noutputs
    col_begnb             => clm3%g%l%c%cnbal%begnb
    col_endnb             => clm3%g%l%c%cnbal%endnb
    col_errnb             => clm3%g%l%c%cnbal%errnb

    ! set time steps
    dt = dtsrf

    err_found = .false.
    ! column loop
    do fc = 1, num_soilc
      c=filter_soilc(fc)

      ! calculate the total column-level nitrogen storage, for
      ! mass conservation check
      col_endnb(c) = totcoln(c)

      ! calculate total column-level inputs

      col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + &
        supplement_to_sminn(c)
      if (crop_prog) col_ninputs(c) = col_ninputs(c) + &
                        fert_to_sminn(c) + soyfixn_to_sminn(c)

      ! calculate total column-level outputs

      col_noutputs(c) = denit(c) + col_fire_nloss(c) + &
        dwt_nloss(c) + product_nloss(c)

#ifndef NITRIF_DENITRIF
      col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
#else
      col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

      col_noutputs(c) = col_noutputs(c) + &
        smin_no3_leached(c) + smin_no3_runoff(c)
#endif

      col_noutputs(c) = col_noutputs(c) - som_n_leached(c)

      ! calculate the total column-level nitrogen balance error
      ! for this time step

      col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
            (col_endnb(c) - col_begnb(c))

      if ( abs(col_errnb(c)) > 1e-4_rk8 .and. 1==2) then
        err_found = .true.
        err_index = c
      end if
    end do ! end of columns loop
    if (err_found) then
      c = err_index
      g = gridcell(c)
      write(stderr,*) 'column nbalance error = ', col_errnb(c), c
      write(stderr,*) 'Latdeg,Londeg=', clm3%g%latdeg(g), clm3%g%londeg(g)
      write(stderr,*) 'begnb       = ',col_begnb(c)
      write(stderr,*) 'endnb       = ',col_endnb(c)
      write(stderr,*) 'delta store = ',col_endnb(c)-col_begnb(c)
      write(stderr,*) 'input mass  = ',col_ninputs(c)*dt
      write(stderr,*) 'output mass = ',col_noutputs(c)*dt
      write(stderr,*) 'net flux    = ',(col_ninputs(c)-col_noutputs(c))*dt
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
  end subroutine NBalanceCheck

#endif

end module mod_clm_cnbalancecheck
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
