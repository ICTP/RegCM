module mod_clm_cnndynamics
#ifdef CN
  !
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam , only : dayspy
  use mod_runparams , only : dtsrf
  use mod_clm_varcon , only : dzsoi_decomp, zisoi
  implicit none

  private

  save

  public :: CNNDeposition
  public :: CNNFixation
  public :: CNNLeaching
  public :: CNNFert
  public :: CNSoyfix

#ifndef NITRIF_DENITRIF
  ! (days) time over which to exponentially relax the npp flux for N
  ! fixation term (if <= 0. or >= 365; use old annual method)
  real(rkx), public :: nfix_timeconst = 0._rkx
#else
  ! (days) time over which to exponentially relax the npp flux for N
  ! fixation term (if <= 0. or >= 365; use old annual method)
  real(rkx), public :: nfix_timeconst = 10._rkx
#endif

  contains
  !
  ! On the radiation time step, update the nitrogen deposition rate
  ! from atmospheric forcing. For now it is assumed that all the atmospheric
  ! N deposition goes to the soil mineral N pool.
  ! This could be updated later to divide the inputs between mineral N absorbed
  ! directly into the canopy and mineral N entering the soil pool.
  !
  subroutine CNNDeposition( lbc, ubc )
    use mod_clm_type
    use mod_clm_atmlnd , only : clm_a2l
    implicit none
    integer(ik4), intent(in) :: lbc, ubc        ! column bounds
    real(rkx) , pointer :: forc_ndep(:)  ! nitrogen deposition rate (gN/m2/s)
    integer(ik4) , pointer :: gridcell(:) ! index into gridcell level quantities
    real(rkx), pointer :: ndep_to_sminn(:)
    integer(ik4) :: g , c  ! indices

    ! Assign local pointers to derived type arrays (in)
    forc_ndep     => clm_a2l%forc_ndep
    gridcell      => clm3%g%l%c%gridcell

    ! Assign local pointers to derived type arrays (out)
    ndep_to_sminn => clm3%g%l%c%cnf%ndep_to_sminn

    ! Loop through columns
    do c = lbc , ubc
      g = gridcell(c)
      ndep_to_sminn(c) = forc_ndep(g)
    end do
  end subroutine CNNDeposition
  !
  ! On the radiation time step, update the nitrogen fixation rate
  ! as a function of annual total NPP. This rate gets updated once per year.
  ! All N fixation goes to the soil mineral N pool.
  !
  subroutine CNNFixation(num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varcon      , only: secspday, spval
    implicit none
    integer(ik4), intent(in) :: num_soilc  ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    real(rkx), pointer :: cannsum_npp(:) ! nitrogen deposition rate (gN/m2/s)
    ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
    real(rkx), pointer :: nfix_to_sminn(:)
    ! (gC/m2/s) lagged net primary production
    real(rkx), pointer :: col_lag_npp(:)
    integer(ik4)  :: c,fc     ! indices
    real(rkx) :: t            ! temporary

    ! Assign local pointers to derived type arrays (in)
    cannsum_npp   => clm3%g%l%c%cps%cannsum_npp

    ! Assign local pointers to derived type arrays (out)
    nfix_to_sminn => clm3%g%l%c%cnf%nfix_to_sminn

    if ( nfix_timeconst > 0._rkx .and. nfix_timeconst < 500._rkx ) then
      col_lag_npp => clm3%g%l%c%cps%col_lag_npp
    end if

    if ( nfix_timeconst > 0._rkx .and. nfix_timeconst < 500._rkx ) then
      ! use exponential relaxation with time constant nfix_timeconst
      ! for NPP - NFIX relation
      ! Loop through columns
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        if (col_lag_npp(c) /= spval) then
          ! need to put npp in units of gC/m^2/year here first
          t = (1.8_rkx * (1._rkx - exp(-0.003_rkx * col_lag_npp(c) * &
                  (secspday * dayspy))))/(secspday * dayspy)
          nfix_to_sminn(c) = max(0._rkx,t)
        else
          nfix_to_sminn(c) = 0._rkx
        end if
      end do
    else
      ! use annual-mean values for NPP-NFIX relation
      do fc = 1 , num_soilc
        c = filter_soilc(fc)

        ! the value 0.001666 is set to give 100 TgN/yr when global
        ! NPP = 60 PgC/yr.  (Cleveland et al., 1999)
        ! Convert from gN/m2/yr -> gN/m2/s
        !t = cannsum_npp(c) * 0.001666_rkx / (secspday * dayspy)
        t = (1.8_rkx * (1._rkx - exp(-0.003_rkx * cannsum_npp(c)))) / &
                (secspday * dayspy)
        nfix_to_sminn(c) = max(0._rkx,t)
        ! PET 2/14/05: commenting out the dependence on NPP, and
        ! forcing Nfix to global constant = 0.4 gN/m2/yr
        !nfix_to_sminn(c) = 0.4 / (secspday*dayspy)

      end do
    end if
  end subroutine CNNFixation
  !
  ! On the radiation time step, update the nitrogen leaching rate
  ! as a function of soluble mineral N and total soil water outflow.
  !
  subroutine CNNLeaching(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varpar      , only : nlevdecomp, nlevsoi
    implicit none
    integer(ik4), intent(in) :: lbc, ubc  ! column bounds
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    ! liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rkx), pointer :: h2osoi_liq(:,:)
    real(rkx), pointer :: qflx_drain(:) ! sub-surface runoff (mm H2O /s)
    !!! awaiting_new_frozen_hydrolgy
    !!! sub-surface runoff from perched wt (mm H2O /s)
    !!! real(rkx), pointer :: qflx_drain_perched(:)
    real(rkx), pointer :: qflx_surf(:)   ! surface runoff (mm H2O /s)
    real(rkx), pointer :: sminn_vr(:,:)  ! (gN/m3) soil mineral N

#ifndef NITRIF_DENITRIF
    ! rate of mineral N leaching (gN/m3/s)
    real(rkx), pointer :: sminn_leached_vr(:,:)
#else
    ! rate of mineral NO3 leaching (gN/m3/s)
    real(rkx), pointer :: smin_no3_leached_vr(:,:)
    ! rate of mineral NO3 loss with runoff (gN/m3/s)
    real(rkx), pointer :: smin_no3_runoff_vr(:,:)
    real(rkx), pointer :: smin_no3_vr(:,:)
#endif
    real(rkx), pointer :: dz(:,:)  !layer thickness (m)

    integer(ik4)  :: j,c,fc  ! indices
    real(rkx) :: dt                 ! radiation time step (seconds)
    real(rkx) :: tot_water(lbc:ubc) ! total column liquid water (kg water/m2)
    ! liquid water to shallow surface depth (kg water/m2)
    real(rkx) :: surface_water(lbc:ubc)
#ifndef NITRIF_DENITRIF
    real(rkx) :: sf         ! soluble fraction of mineral N (unitless)
#else
    real(rkx) :: sf_no3     ! soluble fraction of NO3 (unitless)
#endif
    real(rkx) :: disn_conc  ! dissolved mineral N concentration

    ! (m) depth over which runoff mixes with soil water for N loss to runoff
    real(rkx), parameter :: depth_runoff_Nloss = 0.05
    real(rkx) :: drain_tot(lbc:ubc)  ! total drainage flux (mm H2O /s)

    ! Assign local pointers to derived type arrays (in)
    h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq
    qflx_drain  => clm3%g%l%c%cwf%qflx_drain
    !!! awaiting_new_frozen_hydrolgy qflx_drain_perched
    !!!!                      => clm3%g%l%c%cwf%qflx_drain_perched
    qflx_surf => clm3%g%l%c%cwf%qflx_surf
    sminn_vr => clm3%g%l%c%cns%sminn_vr
    ! Assign local pointers to derived type arrays (out)
#ifndef NITRIF_DENITRIF
    sminn_leached_vr => clm3%g%l%c%cnf%sminn_leached_vr
#else
    smin_no3_leached_vr => clm3%g%l%c%cnf%smin_no3_leached_vr
    smin_no3_runoff_vr  => clm3%g%l%c%cnf%smin_no3_runoff_vr
    smin_no3_vr         => clm3%g%l%c%cns%smin_no3_vr
#endif
    dz => clm3%g%l%c%cps%dz

    ! set time steps
    dt = dtsrf

#ifndef NITRIF_DENITRIF
    ! Assume that 10% of the soil mineral N is in a soluble form
    sf = 0.1_rkx
#else
    ! Assume that 100% of the soil NO3 is in a soluble form
    sf_no3 = 1.0_rkx
#endif

    ! calculate the total soil water
    tot_water(lbc:ubc) = 0._rkx
    do j = 1 , nlevsoi
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
      end do
    end do

    ! for runoff calculation; calculate total water to a given depth
    surface_water(lbc:ubc) = 0._rkx
    do j = 1,nlevsoi
      if ( zisoi(j) <= depth_runoff_Nloss)  then
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
        end do
      else if ( zisoi(j-1) < depth_runoff_Nloss)  then
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * &
                  ( (depth_runoff_Nloss - zisoi(j-1)) / dz(c,j))
        end do
      end if
    end do

    !!! awaiting_new_frozen_hydrolgy ! Loop through columns
    !!! awaiting_new_frozen_hydrolgy do fc = 1,num_soilc
    !!! awaiting_new_frozen_hydrolgy c = filter_soilc(fc)
    !!! awaiting_new_frozen_hydrolgy drain_tot(c) =
    !!!                  qflx_drain(c) + qflx_drain_perched(c)
    !!! awaiting_new_frozen_hydrolgy end do
    ! Loop through columns
    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      drain_tot(c) = qflx_drain(c)
    end do

#ifndef NITRIF_DENITRIF
    do j = 1 , nlevdecomp
      ! Loop through columns
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
#ifndef VERTSOILC
        ! calculate the dissolved mineral N concentration (gN/kg water)
        ! assumes that 10% of mineral nitrogen is soluble
        disn_conc = 0._rkx
        if (tot_water(c) > 0._rkx) then
          disn_conc = (sf * sminn_vr(c,j) ) / tot_water(c)
        end if

        ! calculate the N leaching flux as a function of the dissolved
        ! concentration and the sub-surface drainage flux
        sminn_leached_vr(c,j) = disn_conc * drain_tot(c)
#else
        ! calculate the dissolved mineral N concentration (gN/kg water)
        ! assumes that 10% of mineral nitrogen is soluble
        disn_conc = 0._rkx
        if (h2osoi_liq(c,j) > 0._rkx) then
          disn_conc = (sf * sminn_vr(c,j) * dz(c,j) )/(h2osoi_liq(c,j) )
        end if

        ! calculate the N leaching flux as a function of the dissolved
        ! concentration and the sub-surface drainage flux
        sminn_leached_vr(c,j) = disn_conc * drain_tot(c) * &
                h2osoi_liq(c,j) / ( tot_water(c) * dz(c,j) )

#endif
        ! limit the flux based on current sminn state
        ! only let at most the assumed soluble fraction
        ! of sminn be leached on any given timestep
        sminn_leached_vr(c,j) = min(sminn_leached_vr(c,j), &
                (sf * sminn_vr(c,j))/dt)

        ! limit the flux to a positive value
        sminn_leached_vr(c,j) = max(sminn_leached_vr(c,j), 0._rkx)
      end do
    end do
#else
    ! NITRIF_NITRIF
    do j = 1 , nlevdecomp
      ! Loop through columns
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
#ifndef VERTSOILC
        ! calculate the dissolved mineral N concentration (gN/kg water)
        ! assumes that 10% of mineral nitrogen is soluble
        disn_conc = 0._rkx
        if (tot_water(c) > 0._rkx) then
          disn_conc = (sf_no3 * smin_no3_vr(c,j) )/tot_water(c)
        end if

        ! calculate the N leaching flux as a function of the dissolved
        ! concentration and the sub-surface drainage flux
        smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c)
#else
        ! calculate the dissolved mineral N concentration (gN/kg water)
        ! assumes that 10% of mineral nitrogen is soluble
        disn_conc = 0._rkx
        if (h2osoi_liq(c,j) > 0._rkx) then
          disn_conc = (sf_no3 * smin_no3_vr(c,j) * dz(c,j) )/(h2osoi_liq(c,j) )
        end if
        !
        ! calculate the N leaching flux as a function of the dissolved
        ! concentration and the sub-surface drainage flux
        smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c) * &
                h2osoi_liq(c,j) / ( tot_water(c) * dz(c,j) )
        !
        ! ensure that leaching rate isn't larger than soil N pool
        smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), &
                smin_no3_vr(c,j) / dt )
        !
        ! limit the leaching flux to a positive value
        smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._rkx)
        !
        !
        ! calculate the N loss from surface runoff, assuming a shallow
        ! mixing of surface waters into soil and removal based on runoff
        if ( zisoi(j) <= depth_runoff_Nloss )  then
          smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                  h2osoi_liq(c,j) / ( surface_water(c) * dz(c,j) )
        else if ( zisoi(j-1) < depth_runoff_Nloss )  then
          smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
             h2osoi_liq(c,j) * ((depth_runoff_Nloss - zisoi(j-1)) / &
             dz(c,j)) / ( surface_water(c) * (depth_runoff_Nloss-zisoi(j-1) ))
        else
          smin_no3_runoff_vr(c,j) = 0._rkx
        end if
        !
        ! ensure that runoff rate isn't larger than soil N pool
        smin_no3_runoff_vr(c,j) = min(smin_no3_runoff_vr(c,j), &
                smin_no3_vr(c,j) / dt - smin_no3_leached_vr(c,j))
        !
        ! limit the flux to a positive value
        smin_no3_runoff_vr(c,j) = max(smin_no3_runoff_vr(c,j), 0._rkx)
#endif
        ! limit the flux based on current smin_no3 state
        ! only let at most the assumed soluble fraction
        ! of smin_no3 be leached on any given timestep
        smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), &
                (sf_no3 * smin_no3_vr(c,j))/dt)

        ! limit the flux to a positive value
        smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._rkx)

      end do
    end do
#endif
  end subroutine CNNLeaching
  !
  ! On the radiation time step, update the nitrogen fertilizer for crops
  ! All fertilizer goes into the soil mineral N pool.
  !
  subroutine CNNFert(num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_subgridave , only : p2c
    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns

    real(rkx), pointer :: fert(:)  ! nitrogen fertilizer rate (gN/m2/s)

    real(rkx), pointer :: fert_to_sminn(:)

    ! Assign local pointers to derived type arrays (in)
    fert          => clm3%g%l%c%p%pnf%fert

    ! Assign local pointers to derived type arrays (out)
    fert_to_sminn => clm3%g%l%c%cnf%fert_to_sminn
    call p2c(num_soilc,filter_soilc,fert,fert_to_sminn)
  end subroutine CNNFert
  !
  ! This routine handles the fixation of nitrogen for soybeans based on
  ! the EPICPHASE model M. Cabelguenne et al.,
  ! Agricultural systems 60: 175-196, 1999
  ! N-fixation is based on soil moisture, plant growth phase, and
  ! availibility of nitrogen in the soil root zone.
  !
  subroutine CNSoyfix (num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_pftvarcon , only : nsoybean
    use mod_clm_subgridave , only : p2c

    implicit none
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts

    integer(ik4) , pointer :: ivt(:)         ! pft vegetation type
    integer(ik4) , pointer :: pcolumn(:)     ! pft's column index
    real(rkx), pointer :: fpg(:)   ! fraction of potential gpp (no units)
    real(rkx), pointer :: wf(:)    ! soil water as frac. of whc for top 0.5 m
    ! N flux required to support initial GPP (gN/m2/s)
    real(rkx), pointer :: plant_ndemand(:)
    real(rkx), pointer :: sminn(:)        ! (kgN/m2) soil mineral N
    real(rkx), pointer :: hui(:)          ! =gdd since planting (gddplant)
    real(rkx), pointer :: gddmaturity(:)  ! gdd needed to harvest
    logical , pointer :: croplive(:)  ! true if planted and not harvested

    real(rkx), pointer :: soyfixn(:)  ! nitrogen fixed to each soybean crop
    real(rkx), pointer :: soyfixn_to_sminn(:)

    integer(ik4) :: fp,p,c
    ! soil water factor, nitrogen factor, growth stage factor
    real(rkx):: fxw,fxn,fxg,fxr
    ! difference between nitrogen supply and demand
    real(rkx):: soy_ndemand
    real(rkx):: GDDfrac
    real(rkx):: sminnthreshold1, sminnthreshold2
    real(rkx):: GDDfracthreshold1, GDDfracthreshold2
    real(rkx):: GDDfracthreshold3, GDDfracthreshold4

    ! Assign local pointers to derived type arrays (in)
    ivt           => clm3%g%l%c%p%itype
    pcolumn       => clm3%g%l%c%p%column
    fpg           => clm3%g%l%c%cps%fpg
    wf            => clm3%g%l%c%cps%wf
    plant_ndemand => clm3%g%l%c%p%pepv%plant_ndemand
    sminn         => clm3%g%l%c%cns%sminn
    hui           => clm3%g%l%c%p%pps%gddplant
    gddmaturity   => clm3%g%l%c%p%pps%gddmaturity
    croplive      => clm3%g%l%c%p%pps%croplive

    ! Assign local pointers to derived type arrays (out)
    soyfixn           => clm3%g%l%c%p%pnf%soyfixn
    soyfixn_to_sminn  => clm3%g%l%c%cnf%soyfixn_to_sminn

    sminnthreshold1 = 30._rkx
    sminnthreshold2 = 10._rkx
    GDDfracthreshold1 = 0.15_rkx
    GDDfracthreshold2 = 0.30_rkx
    GDDfracthreshold3 = 0.55_rkx
    GDDfracthreshold4 = 0.75_rkx

    do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! if soybean currently growing then calculate fixation

      if ( ivt(p) == nsoybean .and. croplive(p) ) then

        ! difference between supply and demand

        if ( fpg(c) < 1._rkx ) then
          soy_ndemand = 0._rkx
          soy_ndemand = plant_ndemand(p) - plant_ndemand(p)*fpg(c)

          ! fixation depends on nitrogen, soil water, and growth stage

          ! soil water factor

          fxw = 0._rkx
          fxw = wf(c)/0.85_rkx

          ! soil nitrogen factor (Beth says: CHECK UNITS)

          if ( sminn(c) > sminnthreshold1 ) then
            fxn = 0._rkx
          else if ( sminn(c) > sminnthreshold2 .and. &
                    sminn(c) <= sminnthreshold1 ) then
            fxn = 1.5_rkx - .005_rkx * (sminn(c) * 10._rkx)
          else if (sminn(c) <= sminnthreshold2) then
            fxn = 1._rkx
          end if

          ! growth stage factor
          ! slevis: to replace GDDfrac, assume...
          ! Beth's crit_offset_gdd_def is similar to my gddmaturity
          ! Beth's ac_gdd (base 5C) similar to my hui=gddplant (base 10
          ! for soy)
          ! Ranges below are not firm. Are they lit. based or tuning based?

          GDDfrac = hui(p) / gddmaturity(p)

          if ( GDDfrac > GDDfracthreshold1 .and. &
               GDDfrac <= GDDfracthreshold2 ) then
            fxg = 6.67_rkx * GDDfrac - 1._rkx
          else if ( GDDfrac > GDDfracthreshold2 .and. &
                    GDDfrac <= GDDfracthreshold3 ) then
            fxg = 1._rkx
          else if ( GDDfrac > GDDfracthreshold3 .and. &
                    GDDfrac <= GDDfracthreshold4 ) then
            fxg = 3.75_rkx - 5._rkx * GDDfrac
          end if

          ! calculate the nitrogen fixed by the soybean
          fxr = min(1._rkx, fxw, fxn) * fxg
          fxr = max(0._rkx, fxr)
          soyfixn(p) =  fxr * soy_ndemand
          soyfixn(p) = min(soyfixn(p), soy_ndemand)
        else ! if nitrogen demand met, no fixation
          soyfixn(p) = 0._rkx
        end if
      else ! if not live soybean, no fixation
        soyfixn(p) = 0._rkx
      end if
    end do

    call p2c(num_soilc,filter_soilc,soyfixn,soyfixn_to_sminn)

  end subroutine CNSoyfix

#endif

end module mod_clm_cnndynamics
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
