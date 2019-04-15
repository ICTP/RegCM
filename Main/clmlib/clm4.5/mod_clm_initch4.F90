module mod_clm_initch4
#ifdef LCH4
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_mppparam
  use mod_clm_type
  use mod_clm_varpar , only : nlevsoi , nlevgrnd
  use mod_clm_varpar , only : ngases , nlevdecomp
  use mod_clm_varcon , only : istsoil , istdlak , spval , istcrop
  use mod_clm_decomp , only : get_proc_bounds , get_proc_global
  use mod_clm_ch4varcon , only : allowlakeprod
  !
  ! Contains time constant (and flux / diagnostic vars) and time-varying
  ! (state vars) initialization code for CH4 scheme.
  !
  implicit none

  private

  save

  public :: initch4

  private :: initTimeConst_ch4 ! Set constant parameters.
  private :: makearbinit_ch4   ! Set time-variable parameters for spin up.

  contains
  !
  ! Calls initTimeConst_ch4.
  ! Calls makearbinit_ch4 with logical arbinit. If arbinit == .true. OR
  ! if initial conditions file does not contain methane and oxygen
  ! concentrations, then initializes time varying values. This
  ! allows back-compatibility with initial condition files that have
  ! not been spun up with the new lake code.
  ! In future versions, this could be phased out.
  !
  subroutine initch4( arbinit )
    implicit none
    logical , intent(in) :: arbinit ! Whether mkarbinit has been called.
    call initTimeConst_ch4()
    ! Attn EK
    ! For now
    call makearbinit_ch4(arbinit)
    ! For future versions always using initial condition files spun up
    ! with the new ch4 code: if (arbinit) call makearbinit_ch4(arbinit)
  end subroutine initch4
  !
  ! If arbinit == .true., or if methane & oxygen concentrations (or lake
  ! soil org matter) have not been initialized, then sets time varying values.
  ! Initializes the following time varying variables:
  !
  ! conc_ch4_sat, conc_ch4_unsat, conc_o2_sat, conc_o2_unsat, lake_soilc,
  ! o2stress, finunduated
  !
  subroutine makearbinit_ch4( arbinit )
    implicit none
    logical , intent(in) :: arbinit ! Whether mkarbinit has been called.
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: cellorg
    ! finundated from previous timestep
    real(rk8) , pointer , dimension(:) :: fsat_bef
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_sat
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_unsat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_sat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_unsat
    ! total soil organic matter found in level (g C / m^3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: lake_soilc
    ! time-lagged surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_surf_lag
    ! time-lagged fractional inundated area
    real(rk8) , pointer , dimension(:) :: finundated_lag
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    ! & methanotrophs (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2stress_unsat
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    ! & methanotrophs (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2stress_sat
    ! inundated gridcell fractional area (excluding dedicated wetland columns)
    real(rk8) , pointer , dimension(:) :: finundated
    ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
    real(rk8) , pointer , dimension(:,:) :: layer_sat_lag
    integer(ik4) :: j , l , c , p      ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg  ! per-proc gridcell ending gridcell indices

    if ( myid == italk ) then
      write (stdout,*) &
              'Setting initial data to non-spun up values for CH4 Mod'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit          => clm3%g%l%c%landunit
    conc_ch4_sat       => clm3%g%l%c%cch4%conc_ch4_sat
    conc_ch4_unsat     => clm3%g%l%c%cch4%conc_ch4_unsat
    conc_o2_sat        => clm3%g%l%c%cch4%conc_o2_sat
    conc_o2_unsat      => clm3%g%l%c%cch4%conc_o2_unsat
    lake_soilc         => clm3%g%l%c%cch4%lake_soilc
    cellorg            => clm3%g%l%c%cps%cellorg
    qflx_surf_lag      => clm3%g%l%c%cch4%qflx_surf_lag
    finundated_lag     => clm3%g%l%c%cch4%finundated_lag
    o2stress_sat       => clm3%g%l%c%cch4%o2stress_sat
    o2stress_unsat     => clm3%g%l%c%cch4%o2stress_unsat
    finundated         => clm3%g%l%c%cws%finundated
    fsat_bef           => clm3%g%l%c%cch4%fsat_bef
    layer_sat_lag      => clm3%g%l%c%cch4%layer_sat_lag


    ! Assign local pointers to derived subtypes components (pft-level)

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    do c = begc , endc
      l = clandunit(c)
      if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
        do j = 1 , nlevsoi
          if ( conc_ch4_sat(c,j) == spval .or. arbinit )   &
            conc_ch4_sat(c,j)   = 0._rk8
          if ( conc_ch4_unsat(c,j) == spval .or. arbinit ) &
            conc_ch4_unsat(c,j) = 0._rk8
          if ( conc_o2_sat(c,j) == spval .or. arbinit )    &
            conc_o2_sat(c,j)    = 0._rk8
          if ( conc_o2_unsat(c,j) == spval .or. arbinit )  &
            conc_o2_unsat(c,j)  = 0._rk8
          if ( o2stress_sat(c,j) == spval .or. arbinit )   &
            o2stress_sat(c,j) = 1._rk8
          if ( o2stress_unsat(c,j) == spval .or. arbinit ) &
            o2stress_unsat(c,j) = 1._rk8
          if ( layer_sat_lag(c,j) == spval .or. arbinit )  &
            layer_sat_lag(c,j) = 1._rk8
        end do
        if ( qflx_surf_lag(c) == spval .or. arbinit ) qflx_surf_lag(c) = 0._rk8
        if ( finundated_lag(c) == spval .or. arbinit ) finundated_lag(c) = 0._rk8
        ! finundated will be used to calculate soil decomposition if
        ! anoxia is used
        if ( fsat_bef(c) == spval .or. arbinit ) then
          finundated(c) = 0._rk8
        else
          finundated(c) = fsat_bef(c)
        end if
      else if ( ltype(l) == istdlak ) then
        do j = 1 , nlevsoi
          if ( conc_ch4_sat(c,j) == spval .or. arbinit ) &
            conc_ch4_sat(c,j)   = 0._rk8
          if ( conc_o2_sat(c,j) == spval .or. arbinit )  &
            conc_o2_sat(c,j)    = 0._rk8
          if ( lake_soilc(c,j) == spval .or. arbinit )   &
            lake_soilc(c,j) = 580._rk8*cellorg(c,j)
          ! Need to convert from kg/m^3 organic matter to g C / m^3
          ! (org matter is defined to be 58% C)
        end do
      end if

      ! Set values for all columns equal to zero below nlevsoi
      do j = nlevsoi+1 , nlevgrnd
        conc_ch4_sat(c,j) = 0._rk8
        conc_ch4_unsat(c,j) = 0._rk8
        conc_o2_sat(c,j) = 0._rk8
        conc_o2_unsat(c,j) = 0._rk8
        lake_soilc(c,j) = 0._rk8
        o2stress_sat(c,j) = 1._rk8
        o2stress_unsat(c,j) = 1._rk8
        layer_sat_lag(c,j) = 1._rk8
      end do
    end do
  end subroutine makearbinit_ch4
  !
  ! Initialize variables for ch4 code that will not be input from
  ! restart/inic file. Also set values for inactive CH4 columns to
  ! spval so that they will not be averaged in history file.
  subroutine initTimeConst_ch4
    implicit none
    ! landunit index of column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! gridcell index of column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! landunit type index
    integer(ik4) , pointer , dimension(:) :: ltype
    ! volumetric soil water at saturation (porosity)
    real(rk8) , pointer , dimension(:,:) :: watsat
    ! column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: cellorg
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_sat
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_unsat
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_lake
    ! Total column CH4 aerenchyma (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_aere_sat
    ! Total column CH4 aerenchyma (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_aere_unsat
    ! CH4 ebullition to atmosphere (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_sat
    ! CH4 ebullition to atmosphere (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_unsat
    ! CH4 ebullition to atmosphere (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_lake
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_sat
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_unsat
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_lake
    ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_sat
    ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_unsat
    ! production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_lake
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_ebul_total_sat
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_ebul_total_unsat
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_ebul_depth_sat
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_ebul_depth_unsat
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_aere_depth_sat
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_aere_depth_unsat
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_aere_depth_sat
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_aere_depth_unsat
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_tran_depth_sat
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_tran_depth_unsat
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_oxid_depth_sat
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_oxid_depth_unsat
    ! O2 consumption during decomposition in each soil layer(nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: o2_decomp_depth_sat
    ! O2 consumption during decomposition in each soil layer(nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: o2_decomp_depth_unsat
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_aere_depth_sat
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_aere_depth_unsat
    ! CO2 production during decomposition in each soil layer(nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: co2_decomp_depth_sat
    ! CO2 production during decomposition in each soil layer(nlevsoi) (mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: co2_decomp_depth_unsat
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_oxid_depth_sat
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_oxid_depth_unsat
    ! CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
    real(rk8) , pointer , dimension(:) :: ch4_dfsat_flux
    ! depth of water table for unsaturated fraction (m)
    real(rk8) , pointer , dimension(:) :: zwt_ch4_unsat
    ! tracer conductance for boundary layer [m/s]
    real(rk8) , pointer , dimension(:) :: grnd_ch4_cond
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_lake
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_lake
    ! inundated gridcell fractional area
    real(rk8) , pointer , dimension(:) :: finundated
    ! fraction of potential HR
    real(rk8) , pointer , dimension(:,:) :: fphr
    ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
    real(rk8) , pointer , dimension(:) :: sif
    ! column-averaged root fraction
    real(rk8) , pointer , dimension(:,:) :: rootfr
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    ! & methanotrophs (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2stress_unsat
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    ! & methanotrophs (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2stress_sat
    ! Ratio of methane available to the total per-timestep methane
    ! sinks (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4stress_unsat
    ! Ratio of methane available to the total per-timestep methane
    ! sinks (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4stress_sat
    ! total methane in soil column (g C / m^2)
    real(rk8) , pointer , dimension(:) :: totcolch4
    ! To avoid rare pathologies with allowlakeprod switching between restarts
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_sat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_sat
    ! total soil organic matter found in level (g C / m^3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: lake_soilc
    integer(ik4)  :: g , l , c , p , j        ! indices
    integer(ik4)  :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4)  :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4)  :: begl , endl ! per-proc beginning and ending ldunit indices
    integer(ik4)  :: begg , endg ! per-proc gridcell ending gridcell indices

    watsat     => clm3%g%l%c%cps%watsat

    ch4_surf_diff_sat => clm3%g%l%c%cch4%ch4_surf_diff_sat
    ch4_surf_diff_unsat => clm3%g%l%c%cch4%ch4_surf_diff_unsat
    ch4_surf_diff_lake => clm3%g%l%c%cch4%ch4_surf_diff_lake
    ch4_surf_aere_sat => clm3%g%l%c%cch4%ch4_surf_aere_sat
    ch4_surf_aere_unsat => clm3%g%l%c%cch4%ch4_surf_aere_unsat
    ch4_surf_ebul_sat => clm3%g%l%c%cch4%ch4_surf_ebul_sat
    ch4_surf_ebul_unsat => clm3%g%l%c%cch4%ch4_surf_ebul_unsat
    ch4_surf_ebul_lake => clm3%g%l%c%cch4%ch4_surf_ebul_lake
    ch4_oxid_depth_sat => clm3%g%l%c%cch4%ch4_oxid_depth_sat
    ch4_oxid_depth_unsat => clm3%g%l%c%cch4%ch4_oxid_depth_unsat
    ch4_oxid_depth_lake => clm3%g%l%c%cch4%ch4_oxid_depth_lake
    ch4_prod_depth_sat => clm3%g%l%c%cch4%ch4_prod_depth_sat
    ch4_prod_depth_unsat => clm3%g%l%c%cch4%ch4_prod_depth_unsat
    ch4_prod_depth_lake => clm3%g%l%c%cch4%ch4_prod_depth_lake
    ch4_ebul_total_sat => clm3%g%l%c%cch4%ch4_ebul_total_sat
    ch4_ebul_total_unsat => clm3%g%l%c%cch4%ch4_ebul_total_unsat
    ch4_ebul_depth_sat => clm3%g%l%c%cch4%ch4_ebul_depth_sat
    ch4_ebul_depth_unsat => clm3%g%l%c%cch4%ch4_ebul_depth_unsat
    ch4_aere_depth_sat => clm3%g%l%c%cch4%ch4_aere_depth_sat
    ch4_aere_depth_unsat => clm3%g%l%c%cch4%ch4_aere_depth_unsat
    ch4_tran_depth_sat => clm3%g%l%c%cch4%ch4_tran_depth_sat
    ch4_tran_depth_unsat => clm3%g%l%c%cch4%ch4_tran_depth_unsat
    co2_aere_depth_sat => clm3%g%l%c%cch4%co2_aere_depth_sat
    co2_aere_depth_unsat => clm3%g%l%c%cch4%co2_aere_depth_unsat
    o2_oxid_depth_sat => clm3%g%l%c%cch4%o2_oxid_depth_sat
    o2_oxid_depth_unsat => clm3%g%l%c%cch4%o2_oxid_depth_unsat
    o2_decomp_depth_sat => clm3%g%l%c%cch4%o2_decomp_depth_sat
    o2_decomp_depth_unsat => clm3%g%l%c%cch4%o2_decomp_depth_unsat
    o2_aere_depth_sat => clm3%g%l%c%cch4%o2_aere_depth_sat
    o2_aere_depth_unsat => clm3%g%l%c%cch4%o2_aere_depth_unsat
    co2_decomp_depth_sat => clm3%g%l%c%cch4%co2_decomp_depth_sat
    co2_decomp_depth_unsat => clm3%g%l%c%cch4%co2_decomp_depth_unsat
    co2_oxid_depth_sat => clm3%g%l%c%cch4%co2_oxid_depth_sat
    co2_oxid_depth_unsat => clm3%g%l%c%cch4%co2_oxid_depth_unsat
    ch4_dfsat_flux => clm3%g%l%c%cch4%ch4_dfsat_flux
    zwt_ch4_unsat => clm3%g%l%c%cch4%zwt_ch4_unsat
    grnd_ch4_cond => clm3%g%l%c%cps%pps_a%grnd_ch4_cond
    conc_ch4_lake => clm3%g%l%c%cch4%conc_ch4_lake
    conc_o2_lake => clm3%g%l%c%cch4%conc_o2_lake
    finundated    => clm3%g%l%c%cws%finundated
    fphr       => clm3%g%l%c%cch4%fphr
    sif           => clm3%g%l%c%cch4%sif
    rootfr        => clm3%g%l%c%cps%pps_a%rootfr
    o2stress_unsat   => clm3%g%l%c%cch4%o2stress_unsat
    o2stress_sat     => clm3%g%l%c%cch4%o2stress_sat
    ch4stress_unsat  => clm3%g%l%c%cch4%ch4stress_unsat
    ch4stress_sat    => clm3%g%l%c%cch4%ch4stress_sat
    totcolch4           => clm3%g%l%c%cch4%totcolch4
    conc_ch4_sat => clm3%g%l%c%cch4%conc_ch4_sat
    conc_o2_sat => clm3%g%l%c%cch4%conc_o2_sat
    lake_soilc  => clm3%g%l%c%cch4%lake_soilc
    cellorg            => clm3%g%l%c%cps%cellorg

    if ( myid == italk ) then
      write (stdout,*)'Attempting to initialize non-state variables for CH4 Mod'
    end if

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype           => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit       => clm3%g%l%c%landunit
    cgridcell       => clm3%g%l%c%gridcell

    do c = begc , endc
      l = clandunit(c)

      ! First set levels from nlevsoi+1 to nlevgrnd = 0

      ch4_prod_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_prod_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_prod_depth_lake(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_oxid_depth_lake(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      o2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_tran_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_tran_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_aere_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      co2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_ebul_depth_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4_ebul_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      conc_ch4_lake(c,nlevsoi+1:nlevgrnd) = 0._rk8
      conc_o2_lake(c,nlevsoi+1:nlevgrnd) = 0._rk8
      rootfr(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4stress_unsat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      ch4stress_sat(c,nlevsoi+1:nlevgrnd) = 0._rk8
      fphr(c,nlevdecomp+1:nlevgrnd) = 0._rk8

      if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
        conc_ch4_lake(c,:) = spval
        conc_o2_lake(c,:) = spval
        ch4_surf_diff_lake(c) = spval
        ch4_surf_ebul_lake(c) = spval
        ch4_prod_depth_lake(c,:) = spval
        ch4_oxid_depth_lake(c,:) = spval
      else if ( ltype(l) == istdlak .and. allowlakeprod ) then
        ch4_prod_depth_unsat(c,:) = spval
        ch4_oxid_depth_unsat(c,:) = spval
        o2_oxid_depth_unsat(c,:) = spval
        o2_decomp_depth_unsat(c,:) = spval
        o2_aere_depth_unsat(c,:) = spval
        co2_decomp_depth_unsat(c,:) = spval
        co2_oxid_depth_unsat(c,:) = spval
        ch4_aere_depth_unsat(c,:) = spval
        ch4_tran_depth_unsat(c,:) = spval
        co2_aere_depth_unsat(c,:) = spval
        ch4_surf_aere_unsat(c) = spval
        ch4_ebul_depth_unsat(c,:) = spval
        ch4_ebul_total_unsat(c) = spval
        ch4_surf_ebul_unsat(c) = spval
        ch4_surf_diff_unsat(c) = spval
        ch4_dfsat_flux(c) = spval
        zwt_ch4_unsat(c) = spval
        finundated(c) = spval
        fphr(c,:) = spval
        sif(c)       = spval
        rootfr(c,:) = spval
        o2stress_unsat(c,:) = spval
        ch4stress_unsat(c,:) = spval
      else  ! Inactive CH4 columns
        ch4_prod_depth_sat(c,:) = spval
        ch4_prod_depth_unsat(c,:) = spval
        ch4_prod_depth_lake(c,:) = spval
        ch4_oxid_depth_sat(c,:) = spval
        ch4_oxid_depth_unsat(c,:) = spval
        ch4_oxid_depth_lake(c,:) = spval
        o2_oxid_depth_sat(c,:) = spval
        o2_oxid_depth_unsat(c,:) = spval
        o2_decomp_depth_sat(c,:) = spval
        o2_decomp_depth_unsat(c,:) = spval
        o2_aere_depth_sat(c,:) = spval
        o2_aere_depth_unsat(c,:) = spval
        co2_decomp_depth_sat(c,:) = spval
        co2_decomp_depth_unsat(c,:) = spval
        co2_oxid_depth_sat(c,:) = spval
        co2_oxid_depth_unsat(c,:) = spval
        ch4_aere_depth_sat(c,:) = spval
        ch4_aere_depth_unsat(c,:) = spval
        ch4_tran_depth_sat(c,:) = spval
        ch4_tran_depth_unsat(c,:) = spval
        co2_aere_depth_sat(c,:) = spval
        co2_aere_depth_unsat(c,:) = spval
        ch4_surf_aere_sat(c) = spval
        ch4_surf_aere_unsat(c) = spval
        ch4_ebul_depth_sat(c,:) = spval
        ch4_ebul_depth_unsat(c,:) = spval
        ch4_ebul_total_sat(c) = spval
        ch4_ebul_total_unsat(c) = spval
        ch4_surf_ebul_sat(c) = spval
        ch4_surf_ebul_unsat(c) = spval
        ch4_surf_ebul_lake(c) = spval
        ch4_surf_diff_sat(c) = spval
        ch4_surf_diff_unsat(c) = spval
        ch4_surf_diff_lake(c) = spval
        ch4_dfsat_flux(c) = spval
        zwt_ch4_unsat(c) = spval
        grnd_ch4_cond(c) = spval
        conc_ch4_lake(c,:) = spval
        conc_o2_lake(c,:) = spval
        finundated(c) = spval
        fphr(c,:) = spval
        sif(c)  = spval
        rootfr(c,:) = spval
        o2stress_unsat(c,:) = spval
        o2stress_sat(c,:) = spval
        ch4stress_unsat(c,:) = spval
        ch4stress_sat(c,:) = spval
        ! Set to zero for inactive columns so that this can be used
        ! as an appropriate area-weighted gridcell average soil methane
        ! content.
        totcolch4(c) = 0._rk8
      end if

      if (ltype(l) == istdlak .and. .not. allowlakeprod) then
        ! To avoid rare pathologies with allowlakeprod switching
        ! between restarts
        conc_ch4_sat(c,:) = 0._rk8
        conc_o2_sat(c,:) = 0._rk8
        do j = 1 , nlevsoi
          lake_soilc(c,j) = 580._rk8*cellorg(c,j)
        end do
      end if
    end do

    if ( myid == italk ) then
      write (stdout,*)'Successfully initialized non-state variables for CH4 Mod'
    end if
  end subroutine initTimeConst_ch4

#endif

end module mod_clm_initch4
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
