module mod_clm_driverinit
  !
  ! Initialization of clm driver variables needed from previous timestep
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_clm_varpar , only : nlevsno
  use mod_clm_subgridave , only : p2c
  use mod_clm_varcon , only : h2osno_max , rair , cpair , grav
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_domain , only : ldomain
  use mod_clm_qsat , only : Qsat

  implicit none

  private

  save

  public :: clm_driverInit

  contains
  !
  ! Initialization of clm driver variables needed from previous timestep
  !
  subroutine clm_driverInit(lbc, ubc, lbp, ubp, &
             num_nolakec, filter_nolakec)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc  ! column-index bounds
    integer(ik4) , intent(in) :: lbp , ubp  ! pft-index bounds
    ! number of column non-lake points in column filter
    integer(ik4) , intent(in) :: num_nolakec 
    ! column filter for non-lake points
    integer(ik4) , intent(in) :: filter_nolakec(ubc-lbc+1)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: snl(:) ! number of snow layers
    real(rk8) , pointer :: h2osno(:) ! snow water (mm H2O)
    ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer :: frac_veg_nosno_alb(:)
    ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (pft-level)
    integer(ik4) , pointer :: frac_veg_nosno(:)
    real(rk8) , pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(rk8) , pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    logical , pointer :: do_capsnow(:)     ! true => do snow capping
    ! snow water (mm H2O) at previous time step
    real(rk8) , pointer :: h2osno_old(:)
    ! fraction of ice relative to the tot water
    real(rk8) , pointer :: frac_iceold(:,:)
    integer(ik4) :: g , l , c , p , f , j  ! indices
    ! heat flux from beneath soil/ice column (W/m**2)
    real(rk8) , pointer :: eflx_bot(:)
    ! atmospheric temperature (Kelvin)
    real(rk8) , pointer :: forc_t(:)
    ! atmospheric potential temperature (Kelvin)
    real(rk8) , pointer :: forc_th(:)
    ! atmospheric specific humidity (kg/kg)
    real(rk8) , pointer :: forc_q(:)
    real(rk8) , pointer :: forc_pbot(:)  ! atmospheric pressure (Pa)
    real(rk8) , pointer :: forc_rho(:)   ! atmospheric density (kg/m**3)
    integer(ik4) , pointer :: cgridcell(:) ! column's gridcell
    integer(ik4) , pointer :: clandunit(:) ! column's landunit
    integer(ik4) , pointer :: plandunit(:) ! pft's landunit
    integer(ik4) , pointer :: ityplun(:)   ! landunit type

    ! Assign local pointers to derived type members (landunit-level)

    ityplun            => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    snl                => clm3%g%l%c%cps%snl
    h2osno             => clm3%g%l%c%cws%h2osno
    h2osno_old         => clm3%g%l%c%cws%h2osno_old
    do_capsnow         => clm3%g%l%c%cps%do_capsnow
    frac_iceold        => clm3%g%l%c%cps%frac_iceold
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno
    eflx_bot           => clm3%g%l%c%cef%eflx_bot
    forc_t             => clm3%g%l%c%ces%forc_t
    forc_th            => clm3%g%l%c%ces%forc_th
    forc_q             => clm3%g%l%c%cws%forc_q
    forc_pbot          => clm3%g%l%c%cps%forc_pbot
    forc_rho           => clm3%g%l%c%cps%forc_rho
    clandunit          => clm3%g%l%c%landunit
    cgridcell          => clm3%g%l%c%gridcell

    ! Assign local pointers to derived type members (pft-level)

    pactive            => clm3%g%l%c%p%active
    plandunit          => clm3%g%l%c%p%landunit

    do c = lbc , ubc

      l = clandunit(c)
      g = cgridcell(c)

      ! Initialize column forcing

      forc_t(c)    = clm_a2l%forc_t(g)
      forc_th(c)   = clm_a2l%forc_th(g)
      forc_q(c)    = clm_a2l%forc_q(g)
      forc_pbot(c) = clm_a2l%forc_pbot(g)
      forc_rho(c)  = clm_a2l%forc_rho(g)

      ! Save snow mass at previous time step
      h2osno_old(c) = h2osno(c)

      ! Decide whether to cap snow
      if (h2osno(c) > h2osno_max) then
        do_capsnow(c) = .true.
      else
        do_capsnow(c) = .false.
      end if
      eflx_bot(c) = 0.D0
    end do

    ! Initialize fraction of vegetation not covered by snow (pft-level)

    do p = lbp , ubp
      if (pactive(p)) then
        frac_veg_nosno(p) = frac_veg_nosno_alb(p)
      else
        frac_veg_nosno(p) = 0
      end if
    end do

    ! Initialize set of previous time-step variables
    ! Ice fraction of snow at previous time step

    do j = -nlevsno+1 , 0
      do f = 1 , num_nolakec
        c = filter_nolakec(f)
        if (j >= snl(c) + 1) then
          frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
        end if
      end do
    end do

  end subroutine clm_driverInit

end module mod_clm_driverinit
