module mod_clm_cndecompcascadecentury
#ifdef CN

#ifdef CENTURY_DECOMP
  !
  ! Module that sets the coeffiecients used in the decomposition cascade
  ! submodel.  This uses the CENTURY parameters
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams , only : dtsrf
  use mod_dynparam , only : dayspy
  use mod_mpmessage
  use mod_clm_varpar , only : nlevsoi, nlevgrnd, nlevdecomp
  use mod_clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
  use mod_clm_varpar , only : nsompools
  use mod_clm_varpar , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use mod_clm_varctl , only : spinup_state
  use mod_clm_varcon , only : tfrz , zsoi , rpi
#ifdef LCH4
  use mod_clm_varctl , only : anoxia
  use mod_clm_ch4varcon , only : mino2lim
#endif

  implicit none

  save

  private

  public:: init_decompcascade, decomp_rate_constants

#if (defined VERTSOILC)
  ! (meters) e-folding depth for reduction in decomposition
  ! [set to large number for depth-independance]
  real(rk8), public :: decomp_depth_efolding = 0.50D0
#endif

  ! do we normalize the century decomp. rates so that they match the
  ! CLM Q10 at a given tep?
  logical, public :: normalize_q10_to_century_tfunc = .true.
  ! reference temperature for normalizaion (degrees C)
  real(rk8), public :: normalization_tref = 15.0D0
  logical, public :: use_century_tfunc = .false.
#ifdef LCH4
  ! true ==> weight anoxia by inundated fraction
  logical,  public :: anoxia_wtsat = .false.
#endif
  ! separate q10 for frozen soil respiration rates.
  !  default to same as above zero rates
  real(rk8), public :: froz_q10 = 1.50D0
  ! used here and in ch4Mod
  integer(ik4), public :: nlev_soildecomp_standard

  !! parameters for AD spinup
  ! multipliers for soil decomp during accelerated spinup
  real(rk8), public, parameter :: spinup_vector(nsompools) = &
           (/ 1.00D0, 15.00D0, 675.00D0 /)

  contains
  !
  ! initialize rate constants and decomposition pathways following the
  ! decomposition cascade of the CENTURY model.
  ! written by C. Koven
  !
  subroutine init_decompcascade(begc, endc)
    use mod_clm_type
    implicit none
    ! column level
    ! per-proc beginning and ending column indices
    integer(ik4)  :: begc, endc

    !-- properties of each pathway along decomposition cascade
    ! name of transition
    character(len=8), pointer :: cascade_step_name(:)
    ! respired fraction in decomposition step (frac)
    real(rk8), pointer :: rf_decomp_cascade(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer :: cascade_receiver_pool(:)
    ! what fraction of C leaving a given pool passes through a given
    ! transition (frac)
    real(rk8), pointer :: pathfrac_decomp_cascade(:,:,:)
    !-- properties of each decomposing pool
    ! TRUE => pool has fixed C:N ratio
    logical,  pointer :: floating_cn_ratio_decomp_pools(:)
    ! name of pool for restart files
    character(len=8), pointer :: decomp_pool_name_restart(:)
    ! name of pool for history files
    character(len=8), pointer :: decomp_pool_name_history(:)
    ! name of pool for netcdf long names
    character(len=20), pointer :: decomp_pool_name_long(:)
    ! name of pool for netcdf short names
    character(len=8), pointer :: decomp_pool_name_short(:)
    logical, pointer :: is_litter(:)  ! TRUE => pool is a litter pool
    logical, pointer :: is_soil(:)    ! TRUE => pool is a soil pool
    logical, pointer :: is_cwd(:)     ! TRUE => pool is a cwd pool
    ! c:n ratio for initialization of pools
    real(rk8), pointer :: initial_cn_ratio(:)
    ! initial concentration for seeding at spinup
    real(rk8), pointer :: initial_stock(:)
    logical, pointer :: is_metabolic(:) ! TRUE => pool is metabolic material
    logical, pointer :: is_cellulose(:) ! TRUE => pool is cellulose
    logical, pointer :: is_lignin(:)    ! TRUE => pool is lignin
    real(rk8), pointer :: cellclay(:,:)  ! column 3D clay
    real(rk8), pointer :: cellsand(:,:)  ! column 3D sand
    ! factor for AD spinup associated with each pool
    real(rk8), pointer :: spinup_factor(:)
    real(rk8) :: rf_l1s1
    real(rk8) :: rf_l2s1
    real(rk8) :: rf_l3s2
    real(rk8) :: rf_s1s2(begc:endc,1:nlevdecomp)
    real(rk8) :: rf_s1s3(begc:endc,1:nlevdecomp)
    real(rk8) :: rf_s2s1
    real(rk8) :: rf_s2s3
    real(rk8) :: rf_s3s1
    real(rk8) :: rf_cwdl2
    real(rk8) :: rf_cwdl3
    real(rk8):: cwd_fcel
    real(rk8):: cwd_flig
    real(rk8) :: cn_s1
    real(rk8) :: cn_s2
    real(rk8) :: cn_s3
    real(rk8) :: cn_s4
    real(rk8) :: f_s1s2(begc:endc,1:nlevdecomp)
    real(rk8) :: f_s1s3(begc:endc,1:nlevdecomp)
    real(rk8) :: f_s2s1
    real(rk8) :: f_s2s3

    integer(ik4) :: i_litr1
    integer(ik4) :: i_litr2
    integer(ik4) :: i_litr3
    integer(ik4) :: i_soil1
    integer(ik4) :: i_soil2
    integer(ik4) :: i_soil3
    integer(ik4) :: i_l1s1
    integer(ik4) :: i_l2s1
    integer(ik4) :: i_l3s2
    integer(ik4) :: i_s1s2
    integer(ik4) :: i_s1s3
    integer(ik4) :: i_s2s1
    integer(ik4) :: i_s2s3
    integer(ik4) :: i_s3s1
    integer(ik4) :: i_cwdl2
    integer(ik4) :: i_cwdl3

    integer(ik4) :: c, j     ! indices
    real(rk8) :: t       ! temporary variable

    cascade_step_name       => decomp_cascade_con%cascade_step_name
    rf_decomp_cascade       => clm3%g%l%c%cps%rf_decomp_cascade
    cascade_donor_pool      => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool   => decomp_cascade_con%cascade_receiver_pool
    pathfrac_decomp_cascade => clm3%g%l%c%cps%pathfrac_decomp_cascade
    floating_cn_ratio_decomp_pools => &
            decomp_cascade_con%floating_cn_ratio_decomp_pools
    decomp_pool_name_restart => decomp_cascade_con%decomp_pool_name_restart
    decomp_pool_name_history => decomp_cascade_con%decomp_pool_name_history
    decomp_pool_name_long   => decomp_cascade_con%decomp_pool_name_long
    decomp_pool_name_short  => decomp_cascade_con%decomp_pool_name_short
    is_litter           => decomp_cascade_con%is_litter
    is_soil             => decomp_cascade_con%is_soil
    is_cwd              => decomp_cascade_con%is_cwd
    initial_cn_ratio    => decomp_cascade_con%initial_cn_ratio
    initial_stock       => decomp_cascade_con%initial_stock
    is_metabolic        => decomp_cascade_con%is_metabolic
    is_cellulose        => decomp_cascade_con%is_cellulose
    is_lignin           => decomp_cascade_con%is_lignin
    cellclay            => clm3%g%l%c%cps%cellclay
    cellsand            => clm3%g%l%c%cps%cellsand
    spinup_factor       => decomp_cascade_con%spinup_factor

    !------- time-constant coefficients ---------- !
    ! set soil organic matter compartment C:N ratios
    cn_s1 = 8.00D0
    cn_s2 = 11.00D0
    cn_s3 = 11.00D0

    ! set respiration fractions for fluxes between compartments
    rf_l1s1 = 0.550D0
    rf_l2s1 = 0.50D0
    rf_l3s2 = 0.50D0
    rf_s2s1 = 0.550D0
    rf_s2s3 = 0.550D0
    rf_s3s1 = 0.550D0
    rf_cwdl2 = 0.0D0
    rf_cwdl3 = 0.0D0

    ! set the cellulose and lignin fractions for coarse woody debris
    cwd_fcel = 0.760D0
    cwd_flig = 0.240D0

    ! set path fractions
    f_s2s1 = 0.420D0/(0.450D0)
    f_s2s3 = 0.030D0/(0.450D0)

    ! some of these are dependent on the soil texture properties
    do c = begc , endc
      do j = 1 , nlevdecomp
        t = 0.850D0 - 0.680D0 * 0.010D0 * (100.0D0 - cellsand(c,j))
        f_s1s2(c,j) = 1.0D0 - .0040D0 / (1.0D0 - t)
        f_s1s3(c,j) = .0040D0 / (1.0D0 - t)
        rf_s1s2(c,j) = t
        rf_s1s3(c,j) = t
      end do
    end do

    !-------------------  list of pools and their attributes  ------------

    i_litr1 = i_met_lit
    floating_cn_ratio_decomp_pools(i_litr1) = .true.
    decomp_pool_name_restart(i_litr1) = 'litr1'
    decomp_pool_name_history(i_litr1) = 'LITR1'
    decomp_pool_name_long(i_litr1) = 'litter 1'
    decomp_pool_name_short(i_litr1) = 'L1'
    is_litter(i_litr1) = .true.
    is_soil(i_litr1) = .false.
    is_cwd(i_litr1) = .false.
    initial_cn_ratio(i_litr1) = 90.0D0
    initial_stock(i_litr1) = 0.0D0
    is_metabolic(i_litr1) = .true.
    is_cellulose(i_litr1) = .false.
    is_lignin(i_litr1) = .false.

    i_litr2 = i_cel_lit
    floating_cn_ratio_decomp_pools(i_litr2) = .true.
    decomp_pool_name_restart(i_litr2) = 'litr2'
    decomp_pool_name_history(i_litr2) = 'LITR2'
    decomp_pool_name_long(i_litr2) = 'litter 2'
    decomp_pool_name_short(i_litr2) = 'L2'
    is_litter(i_litr2) = .true.
    is_soil(i_litr2) = .false.
    is_cwd(i_litr2) = .false.
    initial_cn_ratio(i_litr2) = 90.0D0
    initial_stock(i_litr2) = 0.0D0
    is_metabolic(i_litr2) = .false.
    is_cellulose(i_litr2) = .true.
    is_lignin(i_litr2) = .false.

    i_litr3 = i_lig_lit
    floating_cn_ratio_decomp_pools(i_litr3) = .true.
    decomp_pool_name_restart(i_litr3) = 'litr3'
    decomp_pool_name_history(i_litr3) = 'LITR3'
    decomp_pool_name_long(i_litr3) = 'litter 3'
    decomp_pool_name_short(i_litr3) = 'L3'
    is_litter(i_litr3) = .true.
    is_soil(i_litr3) = .false.
    is_cwd(i_litr3) = .false.
    initial_cn_ratio(i_litr3) = 90.0D0
    initial_stock(i_litr3) = 0.0D0
    is_metabolic(i_litr3) = .false.
    is_cellulose(i_litr3) = .false.
    is_lignin(i_litr3) = .true.

    ! CWD
    floating_cn_ratio_decomp_pools(i_cwd) = .true.
    decomp_pool_name_restart(i_cwd) = 'cwd'
    decomp_pool_name_history(i_cwd) = 'CWD'
    decomp_pool_name_long(i_cwd) = 'coarse woody debris'
    decomp_pool_name_short(i_cwd) = 'CWD'
    is_litter(i_cwd) = .false.
    is_soil(i_cwd) = .false.
    is_cwd(i_cwd) = .true.
    initial_cn_ratio(i_cwd) = 90.0D0
    initial_stock(i_cwd) = 0.0D0
    is_metabolic(i_cwd) = .false.
    is_cellulose(i_cwd) = .false.
    is_lignin(i_cwd) = .false.

    i_soil1 = 5
    floating_cn_ratio_decomp_pools(i_soil1) = .false.
    decomp_pool_name_restart(i_soil1) = 'soil1'
    decomp_pool_name_history(i_soil1) = 'SOIL1'
    decomp_pool_name_long(i_soil1) = 'soil 1'
    decomp_pool_name_short(i_soil1) = 'S1'
    is_litter(i_soil1) = .false.
    is_soil(i_soil1) = .true.
    is_cwd(i_soil1) = .false.
    initial_cn_ratio(i_soil1) = cn_s1
    initial_stock(i_soil1) = 20.0D0
    is_metabolic(i_soil1) = .false.
    is_cellulose(i_soil1) = .false.
    is_lignin(i_soil1) = .false.

    i_soil2 = 6
    floating_cn_ratio_decomp_pools(i_soil2) = .false.
    decomp_pool_name_restart(i_soil2) = 'soil2'
    decomp_pool_name_history(i_soil2) = 'SOIL2'
    decomp_pool_name_long(i_soil2) = 'soil 2'
    decomp_pool_name_short(i_soil2) = 'S2'
    is_litter(i_soil2) = .false.
    is_soil(i_soil2) = .true.
    is_cwd(i_soil2) = .false.
    initial_cn_ratio(i_soil2) = cn_s2
    initial_stock(i_soil2) = 20.0D0
    is_metabolic(i_soil2) = .false.
    is_cellulose(i_soil2) = .false.
    is_lignin(i_soil2) = .false.

    i_soil3 = 7
    floating_cn_ratio_decomp_pools(i_soil3) = .false.
    decomp_pool_name_restart(i_soil3) = 'soil3'
    decomp_pool_name_history(i_soil3) = 'SOIL3'
    decomp_pool_name_long(i_soil3) = 'soil 3'
    decomp_pool_name_short(i_soil3) = 'S3'
    is_litter(i_soil3) = .false.
    is_soil(i_soil3) = .true.
    is_cwd(i_soil3) = .false.
    initial_cn_ratio(i_soil3) = cn_s3
    initial_stock(i_soil3) = 20.0D0
    is_metabolic(i_soil3) = .false.
    is_cellulose(i_soil3) = .false.
    is_lignin(i_soil3) = .false.

    spinup_factor(i_litr1) = 1.0D0
    spinup_factor(i_litr2) = 1.0D0
    spinup_factor(i_litr3) = 1.0D0
    spinup_factor(i_cwd) = 1.0D0
    spinup_factor(i_soil1) = spinup_vector(1)
    spinup_factor(i_soil2) = spinup_vector(2)
    spinup_factor(i_soil3) = spinup_vector(3)

    ! list of transitions and their time-independent coefficients
    i_l1s1 = 1
    cascade_step_name(i_l1s1) = 'L1S1'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = rf_l1s1
    cascade_donor_pool(i_l1s1) = i_litr1
    cascade_receiver_pool(i_l1s1) = i_soil1
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = 1.00D0

    i_l2s1 = 2
    cascade_step_name(i_l2s1) = 'L2S1'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s1) = rf_l2s1
    cascade_donor_pool(i_l2s1) = i_litr2
    cascade_receiver_pool(i_l2s1) = i_soil1
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s1)= 1.00D0

    i_l3s2 = 3
    cascade_step_name(i_l3s2) = 'L3S2'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s2) = rf_l3s2
    cascade_donor_pool(i_l3s2) = i_litr3
    cascade_receiver_pool(i_l3s2) = i_soil2
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s2) = 1.00D0

    i_s1s2 = 4
    cascade_step_name(i_s1s2) = 'S1S2'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = &
            rf_s1s2(begc:endc,1:nlevdecomp)
    cascade_donor_pool(i_s1s2) = i_soil1
    cascade_receiver_pool(i_s1s2) = i_soil2
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = &
            f_s1s2(begc:endc,1:nlevdecomp)

    i_s1s3 = 5
    cascade_step_name(i_s1s3) = 'S1S3'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s3) = &
            rf_s1s3(begc:endc,1:nlevdecomp)
    cascade_donor_pool(i_s1s3) = i_soil1
    cascade_receiver_pool(i_s1s3) = i_soil3
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s3) = &
            f_s1s3(begc:endc,1:nlevdecomp)

    i_s2s1 = 6
    cascade_step_name(i_s2s1) = 'S2S1'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s1) = rf_s2s1
    cascade_donor_pool(i_s2s1) = i_soil2
    cascade_receiver_pool(i_s2s1) = i_soil1
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s1) = f_s2s1

    i_s2s3 = 7
    cascade_step_name(i_s2s3) = 'S2S3'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = rf_s2s3
    cascade_donor_pool(i_s2s3) = i_soil2
    cascade_receiver_pool(i_s2s3) = i_soil3
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = f_s2s3

    i_s3s1 = 8
    cascade_step_name(i_s3s1) = 'S3S1'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s1) = rf_s3s1
    cascade_donor_pool(i_s3s1) = i_soil3
    cascade_receiver_pool(i_s3s1) = i_soil1
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s1) = 1.00D0

    i_cwdl2 = 9
    cascade_step_name(i_cwdl2) = 'CWDL2'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
    cascade_donor_pool(i_cwdl2) = i_cwd
    cascade_receiver_pool(i_cwdl2) = i_litr2
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = cwd_fcel

    i_cwdl3 = 10
    cascade_step_name(i_cwdl3) = 'CWDL3'
    rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
    cascade_donor_pool(i_cwdl3) = i_cwd
    cascade_receiver_pool(i_cwdl3) = i_litr3
    pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl3) = cwd_flig

  end subroutine init_decompcascade
  !
  ! calculate rate constants and decomposition pathways for the CENTURY
  ! decomposition cascade model
  ! written by C. Koven based on original CLM4 decomposition cascade
  !
  subroutine decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varcon, only: secspday
    implicit none
    integer(ik4), intent(in) :: lbc, ubc        ! column bounds
    ! number of soil columns in filter
    integer(ik4), intent(in) :: num_soilc
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    ! column level
    ! rate constant for decomposition (1./sec)
    real(rk8), pointer :: decomp_k(:,:,:)
    real(rk8), pointer :: t_scalar(:,:)  ! soil temperature scalar for decomp
    real(rk8), pointer :: w_scalar(:,:)  ! soil water scalar for decomp
    ! fraction by which decomposition is limited by anoxia
    real(rk8), pointer :: o_scalar(:,:)

    real(rk8), pointer :: dz(:,:)        ! soil layer thickness (m)
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer :: t_soisno(:,:)
    real(rk8), pointer :: sucsat(:,:)    ! minimum soil suction (mm)
    ! soil water potential in each soil layer (MPa)
    real(rk8), pointer :: soilpsi(:,:)
#ifdef LCH4
    ! Ratio of oxygen available to that demanded by roots, aerobes, &
    ! methanotrophs (nlevsoi)
    real(rk8), pointer :: o2stress_unsat(:,:)
    ! Ratio of oxygen available to that demanded by roots, aerobes, &
    ! methanotrophs (nlevsoi)
    real(rk8), pointer :: o2stress_sat(:,:)
    real(rk8), pointer :: finundated(:)  ! fractional inundated area
#endif
    integer(ik4), pointer :: alt_indx(:) ! current depth of thaw

    real(rk8):: frw(lbc:ubc)     ! rooting fraction weight
    real(rk8), pointer:: fr(:,:) ! column-level rooting fraction by soil depth
    real(rk8):: minpsi, maxpsi   ! limits for soil water scalar for decomp
    real(rk8):: psi              ! temporary soilpsi for water scalar
    ! real(rk8):: w_scalar(lbc:ubc,1:nlevdecomp)  !soil water scalar for decomp
    real(rk8):: rate_scalar      ! combined rate scalar for decomp
    ! decomposition rate constant litter 1 (1/sec)
    real(rk8):: k_l1
    ! decomposition rate constant litter 2 and litter 3 (1/sec)
    real(rk8):: k_l2_l3
    ! decomposition rate constant SOM 1 (1/sec)
    real(rk8):: k_s1
    ! decomposition rate constant SOM 2 (1/sec)
    real(rk8):: k_s2
    ! decomposition rate constant SOM 3 (1/sec)
    real(rk8):: k_s3
    real(rk8):: k_frag     ! fragmentation rate constant CWD (1/sec)
    real(rk8):: tau_l1     ! turnover time of  litter 1 (yr)
    real(rk8):: tau_l2_l3  ! turnover time of  litter 2 and litter 3 (yr)
    real(rk8):: tau_l3     ! turnover time of  litter 3 (yr)
    real(rk8):: tau_s1     ! turnover time of  SOM 1 (yr)
    real(rk8):: tau_s2     ! turnover time of  SOM 2 (yr)
    real(rk8):: tau_s3     ! turnover time of  SOM 3 (yr)
    real(rk8):: tau_cwd    ! corrected fragmentation rate constant CWD
    real(rk8):: cwd_fcel   ! cellulose fraction of coarse woody debris
    real(rk8):: cwd_flig   ! lignin fraction of coarse woody debris
    real(rk8):: cwdc_loss  ! fragmentation rate for CWD carbon (gC/m2/s)
    real(rk8):: cwdn_loss  ! fragmentation rate for CWD nitrogen (gN/m2/s)

    integer(ik4) :: i_litr1
    integer(ik4) :: i_litr2
    integer(ik4) :: i_litr3
    integer(ik4) :: i_soil1
    integer(ik4) :: i_soil2
    integer(ik4) :: i_soil3
    integer(ik4) :: c, fc, j, k, l
    real(rk8) :: q10 = 1.50D0
    real(rk8) :: catanf    ! hyperbolic temperature function from CENTURY
    real(rk8) :: catanf_30 ! reference rate at 30C
    real(rk8) :: t1        ! temperature argument

    ! factor by which to offset the decomposition rates frm century to
    ! a q10 formulation
    real(rk8) :: normalization_factor

#if (defined VERTSOILC)
    real(rk8) :: depth_scalar(lbc:ubc,1:nlevdecomp)
#endif

    !----- CENTURY T response function
    catanf(t1) = 11.750D0 +(29.70D0 / rpi) * &
            atan( rpi * 0.0310D0  * ( t1 - 15.40D0 ))

    ! Assign local pointers to derived type arrays
    t_soisno    => clm3%g%l%c%ces%t_soisno
    sucsat      => clm3%g%l%c%cps%sucsat
    soilpsi     => clm3%g%l%c%cps%soilpsi
    dz          => clm3%g%l%c%cps%dz
    t_scalar    => clm3%g%l%c%ccf%t_scalar
    w_scalar    => clm3%g%l%c%ccf%w_scalar
    o_scalar    => clm3%g%l%c%ccf%o_scalar
    decomp_k    => clm3%g%l%c%ccf%decomp_k
#ifdef LCH4
    o2stress_sat    => clm3%g%l%c%cch4%o2stress_sat
    o2stress_unsat  => clm3%g%l%c%cch4%o2stress_unsat
    finundated      => clm3%g%l%c%cws%finundated
#endif
    alt_indx        => clm3%g%l%c%cps%alt_indx

    if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
      call fatal(__FILE__,__LINE__, &
        'ERROR: cannot have both use_century_tfunc and &
        &normalize_q10_to_century_tfunc set as true' )
    endif

    ! tau (yrs) at reference temperature
    ! the aboveground parameters from century
    ! tau_l1 = 1./14.8
    ! tau_l2 = 1./3.9
    ! tau_s1 = 1./6.0
    ! tau_l2 = 1./0.2
    ! tau_l3 = 1./.0045

    ! the belowground parameters from century
    tau_l1 = 1.0D0/18.5D0
    tau_l2_l3 = 1.0D0/4.9D0
    tau_s1 = 1.0D0/7.3D0
    tau_s2 = 1.0D0/0.2D0
    tau_s3 = 1.0D0/.0045D0

    ! century leaves wood decomposition rates open, within range of 0-0.5 yr^-1
    tau_cwd  = 1.0D0/0.3D0

    ! translate to per-second time constant
    k_l1 = 1.0D0 / (secspday * dayspy * tau_l1)
    k_l2_l3 = 1.0D0 / (secspday * dayspy * tau_l2_l3)
    k_s1 = 1.0D0 / (secspday * dayspy * tau_s1)
    k_s2 = 1.0D0 / (secspday * dayspy * tau_s2)
    k_s3 = 1.0D0 / (secspday * dayspy * tau_s3)
    k_frag = 1.0D0 / (secspday * dayspy * tau_cwd)

    ! calc ref rate
    catanf_30 = catanf(30.0D0)
    ! The following code implements the acceleration part of the AD
    ! spinup algorithm

    if ( spinup_state == 1 ) then
      k_s1 = k_s1 * spinup_vector(1)
      k_s2 = k_s2 * spinup_vector(2)
      k_s3 = k_s3 * spinup_vector(3)
    end if

    i_litr1 = 1
    i_litr2 = 2
    i_litr3 = 3
    i_soil1 = 5
    i_soil2 = 6
    i_soil3 = 7

    !--- time dependent coefficients-----!
    if ( nlevdecomp == 1 ) then
      ! calculate function to weight the temperature and water potential
      ! scalars for decomposition control.

      ! the following normalizes values in fr so that they
      ! sum to 1.0 across top nlevdecomp levels on a column
      frw(lbc:ubc) = 0.0D0
      nlev_soildecomp_standard=5
      allocate(fr(lbc:ubc,nlev_soildecomp_standard))
      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          frw(c) = frw(c) + dz(c,j)
        end do
      end do
      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if (frw(c) /= 0.0D0) then
            fr(c,j) = dz(c,j) / frw(c)
          else
            fr(c,j) = 0.0D0
          end if
        end do
      end do

      if ( .not. use_century_tfunc ) then
        ! calculate rate constant scalar for soil temperature
        ! assuming that the base rate constants are assigned for non-moisture
        ! limiting conditions at 25 C.

        do j = 1 , nlev_soildecomp_standard
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            if (j==1) t_scalar(c,:) = 0.0D0
            !! use separate (possibly equal) t funcs above and below
            !! freezing point
            !! t_scalar(c,1) = t_scalar(c,1) + &
            !!   (q10**((t_soisno(c,j)-(tfrz+25.0D0))/10.0D0))*fr(c,j)
            if (t_soisno(c,j) >= tfrz) then
              t_scalar(c,1) = t_scalar(c,1) + &
                      (q10**((t_soisno(c,j)-(tfrz+25.0D0))/10.0D0))*fr(c,j)
            else
              t_scalar(c,1) = t_scalar(c,1) + (q10**(-25.0D0/10.0D0))* &
                      (froz_q10**((t_soisno(c,j)-tfrz)/10.0D0))*fr(c,j)
            end if
          end do
        end do
      else
        ! original century uses an arctangent function to calculate the
        ! temperature dependence of decomposition
        do j = 1 , nlev_soildecomp_standard
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            if (j==1) t_scalar(c,:) = 0.0D0
            t_scalar(c,1) = t_scalar(c,1) + &
                    max(catanf(t_soisno(c,j)-tfrz)/catanf_30*fr(c,j),0.010D0)
          end do
        end do
      end if

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987.
      ! Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983.
      ! Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

      minpsi = -10.00D0;

      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if (j==1) w_scalar(c,:) = 0.0D0
          maxpsi = sucsat(c,j) * (-9.8D-6)
          psi = min(soilpsi(c,j),maxpsi)
          ! decomp only if soilpsi is higher than minpsi
          if (psi > minpsi) then
            w_scalar(c,1) = w_scalar(c,1) + &
                    (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
          end if
        end do
      end do

#ifdef LCH4
      if (anoxia_wtsat) then ! Adjust for saturated fraction if unfrozen
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( alt_indx(c) >= nlev_soildecomp_standard .and. &
               t_soisno(c,1) > tfrz) then
            w_scalar(c,1) = w_scalar(c,1) * &
                    (1.0D0 - finundated(c)) + finundated(c)
          end if
        end do
      end if
#endif

#ifdef LCH4
      ! Calculate ANOXIA
      if (anoxia) then
        ! Check for anoxia w/o LCH4 now done in controlMod.

        do j = 1 , nlev_soildecomp_standard
          do fc = 1 , num_soilc
            c = filter_soilc(fc)

            if (j==1) o_scalar(c,:) = 0.0D0

            if (.not. anoxia_wtsat) then
              o_scalar(c,1) = o_scalar(c,1) + &
                      fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
            else
              o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * &
                (max(o2stress_unsat(c,j), mino2lim)*(1.0D0 - finundated(c)) + &
                 max(o2stress_sat(c,j), mino2lim)*finundated(c) )
            end if
          end do
        end do
      else
        o_scalar(lbc:ubc,1:nlevdecomp) = 1.0D0
      end if
#else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1.0D0
#endif

      deallocate(fr)

    else

      if ( .not. use_century_tfunc ) then
        ! calculate rate constant scalar for soil temperature
        ! assuming that the base rate constants are assigned for non-moisture
        ! limiting conditions at 25 C.
        ! Peter Thornton: 3/13/09
        ! Replaced the Lloyd and Taylor function with a Q10 formula,
        ! with Q10 = 1.5 as part of the modifications made to improve
        ! the seasonal cycle of atmospheric CO2 concentration in global
        ! simulations. This does not impact the base rates at 25 C,
        ! which are calibrated from microcosm studies.

        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            !! use separate (possibly equal) t funcs above and below
            !! freezing point
            !! t_scalar(c,j)= (q10**((t_soisno(c,j)-(tfrz+25.0D0))/10.0D0))
            if (t_soisno(c,j) >= tfrz) then
              t_scalar(c,j)= (q10**((t_soisno(c,j)-(tfrz+25.0D0))/10.0D0))
            else
              t_scalar(c,j)= (q10**(-25.0D0/10.0D0)) * &
                      (froz_q10**((t_soisno(c,j)-tfrz)/10.0D0))
            end if
          end do
        end do
      else
        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            t_scalar(c,j)= max(catanf(t_soisno(c,j)-tfrz)/catanf_30, 0.010D0)
          end do
        end do
      endif

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987.
      ! Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983.
      ! Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

      minpsi = -10.00D0;
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          maxpsi = sucsat(c,j) * (-9.8D-6)
          psi = min(soilpsi(c,j),maxpsi)
          ! decomp only if soilpsi is higher than minpsi
          if (psi > minpsi) then
            w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
          else
            w_scalar(c,j) = 0.0D0
          end if
#ifdef LCH4
          if (anoxia_wtsat .and. t_soisno(c,j) > tfrz) then
            ! wet area will have w_scalar of 1 if unfrozen
            w_scalar(c,j) = w_scalar(c,j) * &
                    (1.0D0 - finundated(c)) + finundated(c)
          end if
#endif
        end do
      end do
#ifdef LCH4
      ! Calculate ANOXIA
      ! Check for anoxia w/o LCH4 now done in controlMod.
      if (anoxia) then
        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc(fc)

            if (.not. anoxia_wtsat) then
              o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
            else
              o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim) * &
                      (1.0D0 - finundated(c)) + &
                       max(o2stress_sat(c,j), mino2lim) * finundated(c)
            end if
          end do
        end do
      else
        o_scalar(lbc:ubc,1:nlevdecomp) = 1.0D0
      end if
#else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1.0D0
#endif
    end if

    if ( normalize_q10_to_century_tfunc ) then
      ! scale all decomposition rates by a constant to compensate for
      ! offset between original CENTURY temp func and Q10
      normalization_factor = (catanf(normalization_tref)/catanf_30) / &
              (q10**((normalization_tref-25.0D0)/10.0D0))
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          t_scalar(c,j) = t_scalar(c,j) * normalization_factor
        end do
      end do
    endif

#if (defined VERTSOILC)
    ! add a term to reduce decomposition rate at depth
    ! for now used a fixed e-folding depth
    do j = 1 , nlevdecomp
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding)
      end do
    end do
#endif

#if (defined VERTSOILC)
    do j = 1 , nlevdecomp
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * &
                w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
      end do
    end do
#else
    do j = 1 , nlevdecomp
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
        decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * &
                w_scalar(c,j) * o_scalar(c,j)
      end do
    end do
#endif
  end subroutine decomp_rate_constants

#endif

#endif

end module mod_clm_cndecompcascadecentury
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
