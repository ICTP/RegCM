module mod_clm_accflds
  !
  ! This module contains subroutines that initialize, update and extract
  ! the user-specified fields over user-defined intervals. Each interval
  ! and accumulation type is unique to each field processed.
  ! Subroutine [initAccumFlds] defines the fields to be processed
  ! and the type of accumulation. Subroutine [updateAccumFlds] does
  ! the actual accumulation for a given field. Fields are accumulated
  ! by calls to subroutine [update_accum_field]. To accumulate a field,
  ! it must first be defined in subroutine [initAccumFlds] and then
  ! accumulated by calls to [updateAccumFlds].
  ! Four types of accumulations are possible:
  !   o average over time interval
  !   o running mean over time interval
  !   o running accumulation over time interval
  ! Time average fields are only valid at the end of the averaging interval.
  ! Running means are valid once the length of the simulation exceeds the
  ! averaging interval. Accumulated fields are continuously accumulated.
  ! The trigger value "-99999." resets the accumulation to zero.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_date
  use mod_mpmessage
  use mod_runparams
  use mod_clm_type
  use mod_clm_decomp , only : get_proc_bounds , get_proc_global
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_accumul , only : init_accum_field , print_accum_fields
  use mod_clm_varcon , only : secspday , tfrz , spval
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_time_manager , only : is_end_curr_day
  use mod_clm_accumul , only : update_accum_field , extract_accum_field
  use mod_clm_pftvarcon , only : nwcereal , nwcerealirrig , mxtmp , baset
  use mod_clm_time_manager , only : get_start_date
  use mod_clm_pftvarcon , only : ndllf_dcd_brl_tree
  use mod_clm_varctl , only : nsrest , nsrStartup , nextdate

  implicit none

  private

  save

  ! Initialization accumulator fields
  public :: initAccFlds
  ! Initialize clmtype variables obtained from accum fields
  public :: initAccClmtype
  ! Update accumulator fields
  public :: updateAccFlds

  contains
  !
  ! Initializes accumulator and sets up array of accumulated fields
  !
  subroutine initAccFlds()
    implicit none
    integer(ik4) , parameter :: not_used = bigint

    ! Hourly average of 2m temperature.

    call init_accum_field(fname='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600.0_rk8/dtsrf), &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! Hourly average of Urban 2m temperature.

    call init_accum_field(fname='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600.0_rk8/dtsrf), &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! Hourly average of Rural 2m temperature.

    call init_accum_field(fname='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600.0_rk8/dtsrf), &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 24hr average of vegetation temperature (heald, 04/06)
    call init_accum_field (fname='T_VEG24', units='K', &
         desc='24hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 240hr average of vegetation temperature (heald, 04/06)
    call init_accum_field (fname='T_VEG240', units='K', &
         desc='240hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 24hr average of direct solar radiation (heald, 04/06)
    call init_accum_field (fname='FSD24', units='W/m2', &
         desc='24hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 240hr average of direct solar radiation (heald, 04/06)
    call init_accum_field (fname='FSD240', units='W/m2', &
         desc='240hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 24hr average of diffuse solar radiation (heald, 04/06)
    call init_accum_field (fname='FSI24', units='W/m2', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 240hr average of diffuse solar radiation (heald, 04/06)
    call init_accum_field (fname='FSI240', units='W/m2', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 24hr average of fraction of canopy that is sunlit (heald, 04/06)
    call init_accum_field (fname='FSUN24', units='fraction', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! 240hr average of fraction of canopy that is sunlit (heald, 04/06)
    call init_accum_field (fname='FSUN240', units='fraction', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! Average of LAI from previous and current timestep (heald, 04/06)
    call init_accum_field (fname='LAIP', units='m2/m2', &
         desc='leaf area index average over timestep', &
         accum_type='runmean', accum_period=1, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! The following is a running mean.
    ! The accumulation period is set to -10 for a 10-day running mean.

    call init_accum_field (fname='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=tfrz+20.0_rk8)

#if (defined CNDV)
    ! 30-day average of 2m temperature.

    call init_accum_field (fname='TDA', units='K', &
         desc='30-day average of 2-m temperature', &
         accum_type='timeavg', accum_period=-30, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! The following are running means.
    ! The accumulation period is set to -365 for a 365-day running mean.

    call init_accum_field (fname='PREC365', units='MM H2O/S', &
         desc='365-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-365, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    ! The following are accumulated fields.
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field).
    ! Hence, [accper] is not valid.

    call init_accum_field (fname='AGDDTW', units='K', &
         desc='growing degree-days base twmax', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    call init_accum_field (fname='AGDD', units='K', &
         desc='growing degree-days base 5C', &
         accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)
#endif

      call init_accum_field (fname='PREC60', units='MM H2O/S', &
         desc='60-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-60, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

         call init_accum_field (fname='PREC10', units='MM H2O/S', &
         desc='10-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.0_rk8)

    if ( crop_prog )then
       ! 10-day average of min 2m temperature.

       call init_accum_field (fname='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=tfrz)

       ! 5-day average of min 2m temperature.

       call init_accum_field (fname='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=tfrz)

       ! All GDD summations are relative to the planting date
       ! (Kucharik & Brye 2003)

       call init_accum_field (fname='GDD0', units='K', &
            desc='growing degree-days base 0C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0.0_rk8)

       call init_accum_field (fname='GDD8', units='K', &
            desc='growing degree-days base 8C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0.0_rk8)

       call init_accum_field (fname='GDD10', units='K', &
            desc='growing degree-days base 10C from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0.0_rk8)

       call init_accum_field (fname='GDDPLANT', units='K', &
            desc='growing degree-days from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0.0_rk8)

       call init_accum_field (fname='GDDTSOI', units='K', &
            desc='growing degree-days from planting (top two soil layers)', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0.0_rk8)
    end if

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds
  !
  ! Update and/or extract accumulated fields
  !
  subroutine updateAccFlds()
    implicit none
    integer(ik4) , pointer :: itype(:)      ! pft vegetation
    ! index into gridcell level quantities
    integer(ik4) , pointer :: pgridcell(:)
    integer(ik4) , pointer :: plandunit(:)
    real(rk8), pointer :: forc_t(:)     ! atmospheric temperature (Kelvin)
    real(rk8), pointer :: forc_rain(:)  ! rain rate [mm/s]
    real(rk8), pointer :: forc_snow(:)  ! snow rate [mm/s]
    ! 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: t_ref2m(:)
    ! Urban 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: t_ref2m_u(:)
    ! Rural 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: t_ref2m_r(:)
    logical , pointer :: urbpoi(:)     ! true => landunit is an urban point
    logical , pointer :: ifspecial(:)  ! true => landunit is not vegetated
    ! landunit index associated with each pft
    real(rk8), pointer :: vf(:)         ! vernalization factor
    real(rk8), pointer :: t_soisno(:,:) ! soil temperature (K)
    real(rk8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8), pointer :: watsat(:,:)
    real(rk8), pointer :: dz(:,:)       ! layer thickness depth (m)
    real(rk8), pointer :: latdeg(:)     ! latitude (radians)
    logical , pointer :: croplive(:)   ! Flag, true if planted, not harvested
    integer(ik4) , pointer :: pcolumn(:)    ! index into column level quantities
    ! heald (04/06): variables to be accumulated for VOC emissions
    real(rk8), pointer :: t_veg(:)  ! pft vegetation temperature (Kelvin)
    ! direct beam radiation (visible only)
    real(rk8), pointer :: forc_solad(:,:)
    ! diffuse radiation     (visible only)
    real(rk8), pointer :: forc_solai(:,:)
    real(rk8), pointer :: fsun(:) ! sunlit fraction of canopy
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer :: elai(:)
    ! heald (04/06): accumulated variables for VOC emissions
    ! 24hr average vegetation temperature (K)
    real(rk8), pointer :: t_veg24(:)
    ! 240hr average vegetation temperature (Kelvin)
    real(rk8), pointer :: t_veg240(:)
    real(rk8), pointer :: fsd24(:)  ! 24hr average of direct beam radiation
    real(rk8), pointer :: fsd240(:) ! 240hr average of direct beam radiation
    real(rk8), pointer :: fsi24(:)  ! 24hr average of diffuse beam radiation
    real(rk8), pointer :: fsi240(:) ! 240hr average of diffuse beam radiation
    real(rk8), pointer :: fsun24(:) ! 24hr average of sunlit fraction of canopy
    real(rk8), pointer :: fsun240(:)! 240hr average of sunlit fraction of canopy
    real(rk8), pointer :: elai_p(:) ! leaf area index average over timestep

    ! daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min(:)
    ! daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max(:)
    ! instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst(:)
    ! instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst(:)
    ! Urban daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min_u(:)
    ! Rural daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min_r(:)
    ! Urban daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max_u(:)
    ! Rural daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max_r(:)
    ! Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst_u(:)
    ! Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst_r(:)
    ! Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst_u(:)
    ! Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst_r(:)
    ! 10-day running mean of the 2 m temperature (K)
    real(rk8), pointer :: t10(:)
#if (defined CNDV)
    ! 30-day average temperature (Kelvin)
    real(rk8), pointer :: t_mo(:)
    ! annual min of t_mo (Kelvin)
    real(rk8), pointer :: t_mo_min(:)
    ! 365-day running mean of tot. precipitation
    real(rk8), pointer :: prec365(:)
    ! accumulated growing degree days above twmax
    real(rk8), pointer :: agddtw(:)
    ! accumulated growing degree days above 5
    real(rk8), pointer :: agdd(:)
    ! upper limit of temperature of the warmest month
    real(rk8), pointer :: twmax(:)
#endif
    ! 10-day running mean of tot. precipitation
    real(rk8), pointer :: prec10(:)
    ! 60-day running mean of tot. precipitation
    real(rk8), pointer :: prec60(:)
    ! growing degree-days base 0C'
    real(rk8), pointer :: gdd0(:)
    ! growing degree-days base 8C from planting
    real(rk8), pointer :: gdd8(:)
    ! growing degree-days base 10C from planting
    real(rk8), pointer :: gdd10(:)
    ! growing degree-days from planting
    real(rk8), pointer :: gddplant(:)
    ! growing degree-days from planting (top two soil layers)
    real(rk8), pointer :: gddtsoi(:)
    ! 10-day running mean of min 2-m temperature
    real(rk8), pointer :: a10tmin(:)
    ! 5-day running mean of min 2-m temperature
    real(rk8), pointer :: a5tmin(:)
    integer(ik4) :: g , l , c , p ! indices
    integer(ik4) :: itypveg  ! vegetation type
    integer(ik4) :: year     ! year (0, ...) for nstep
    integer(ik4) :: month    ! month (1, ..., 12) for nstep
    integer(ik4) :: day      ! day of month (1, ..., 31) for nstep
    integer(ik4) :: secs     ! seconds into current date for nstep
    logical :: end_cd        ! temporary for is_end_curr_day() value
    integer(ik4) :: ier    ! error status
    integer(ik4) :: begp, endp !  per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc !  per-proc beginning and ending column indices
    integer(ik4) :: begl, endl !  per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg !  per-proc gridcell ending gridcell indices
    real(rk8), pointer :: rbufslp(:)      ! temporary single level - pft level
    integer(ik8) :: kkincr

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_t     => clm_a2l%forc_t
    forc_rain  => clm_a2l%forc_rain
    forc_snow  => clm_a2l%forc_snow
    forc_solad => clm_a2l%forc_solad     ! (heald 04/06)
    forc_solai => clm_a2l%forc_solai     ! (heald 04/06)

    ! Assign local pointers to derived subtypes components (landunit-level)
    ifspecial  => clm3%g%l%ifspecial
    urbpoi     => clm3%g%l%urbpoi

    ! Assign local pointers to derived subtypes components (pft-level)

    itype            => clm3%g%l%c%p%itype
    pgridcell        => clm3%g%l%c%p%gridcell
    t_ref2m          => clm3%g%l%c%p%pes%t_ref2m
    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_ref2m_u        => clm3%g%l%c%p%pes%t_ref2m_u
    t_ref2m_r        => clm3%g%l%c%p%pes%t_ref2m_r
    t_ref2m_max_u    => clm3%g%l%c%p%pes%t_ref2m_max_u
    t_ref2m_max_r    => clm3%g%l%c%p%pes%t_ref2m_max_r
    t_ref2m_min_u    => clm3%g%l%c%p%pes%t_ref2m_min_u
    t_ref2m_min_r    => clm3%g%l%c%p%pes%t_ref2m_min_r
    t_ref2m_max_inst_u => clm3%g%l%c%p%pes%t_ref2m_max_inst_u
    t_ref2m_max_inst_r => clm3%g%l%c%p%pes%t_ref2m_max_inst_r
    t_ref2m_min_inst_u => clm3%g%l%c%p%pes%t_ref2m_min_inst_u
    t_ref2m_min_inst_r => clm3%g%l%c%p%pes%t_ref2m_min_inst_r
    plandunit        => clm3%g%l%c%p%landunit
    t10              => clm3%g%l%c%p%pes%t10
    a10tmin          => clm3%g%l%c%p%pes%a10tmin
    a5tmin           => clm3%g%l%c%p%pes%a5tmin
#if (defined CNDV)
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    t_mo_min         => clm3%g%l%c%p%pdgvs%t_mo_min
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd
    twmax            => dgv_pftcon%twmax
#endif
    prec60           => clm3%g%l%c%p%pps%prec60
    prec10           => clm3%g%l%c%p%pps%prec10
    gdd0             => clm3%g%l%c%p%pps%gdd0
    gdd8             => clm3%g%l%c%p%pps%gdd8
    gdd10            => clm3%g%l%c%p%pps%gdd10
    gddplant         => clm3%g%l%c%p%pps%gddplant
    gddtsoi          => clm3%g%l%c%p%pps%gddtsoi
    vf               => clm3%g%l%c%p%pps%vf
    t_soisno         => clm3%g%l%c%ces%t_soisno
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    watsat           => clm3%g%l%c%cps%watsat
    dz               => clm3%g%l%c%cps%dz
    latdeg           => clm3%g%latdeg
    croplive         => clm3%g%l%c%p%pps%croplive
    pcolumn          => clm3%g%l%c%p%column
    t_veg24          => clm3%g%l%c%p%pvs%t_veg24           ! (heald 04/06)
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240          ! (heald 04/06)
    fsd24            => clm3%g%l%c%p%pvs%fsd24             ! (heald 04/06)
    fsd240           => clm3%g%l%c%p%pvs%fsd240            ! (heald 04/06)
    fsi24            => clm3%g%l%c%p%pvs%fsi24             ! (heald 04/06)
    fsi240           => clm3%g%l%c%p%pvs%fsi240            ! (heald 04/06)
    fsun24           => clm3%g%l%c%p%pvs%fsun24            ! (heald 04/06)
    fsun240          => clm3%g%l%c%p%pvs%fsun240           ! (heald 04/06)
    elai_p           => clm3%g%l%c%p%pvs%elai_p            ! (heald 04/06)
    t_veg            => clm3%g%l%c%p%pes%t_veg             ! (heald 04/06)
    fsun             => clm3%g%l%c%p%pps%fsun              ! (heald 04/06)
    elai             => clm3%g%l%c%p%pps%elai              ! (heald 04/06)

    ! Determine calendar information

    call curr_date (nextdate, year, month, day, secs)

    ! Don't do any accumulation if nstep is zero
    ! (only applies to coupled or cam mode)

    if ( ktau < ntsrf ) return

    kkincr = ktau / ntsrf

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if ( ier/=0 ) then
      write(stderr,*)'update_accum_hist allocation error for rbuf1dp'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', t_ref2m, kkincr)
    call extract_accum_field ('TREFAV', rbufslp, kkincr)
    end_cd = is_end_curr_day()
    do p = begp , endp
      if (rbufslp(p) /= spval) then
        t_ref2m_max_inst(p) = max(rbufslp(p), t_ref2m_max_inst(p))
        t_ref2m_min_inst(p) = min(rbufslp(p), t_ref2m_min_inst(p))
      end if
      if (end_cd) then
        t_ref2m_max(p) = t_ref2m_max_inst(p)
        t_ref2m_min(p) = t_ref2m_min_inst(p)
        t_ref2m_max_inst(p) = -spval
        t_ref2m_min_inst(p) =  spval
      else if (secs <= int(dtsrf)) then
        t_ref2m_max(p) = spval
        t_ref2m_min(p) = spval
      end if
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', t_ref2m_u, kkincr)
    call extract_accum_field ('TREFAV_U', rbufslp, kkincr)
    do p = begp , endp
      l = plandunit(p)
      if (rbufslp(p) /= spval) then
        t_ref2m_max_inst_u(p) = max(rbufslp(p), t_ref2m_max_inst_u(p))
        t_ref2m_min_inst_u(p) = min(rbufslp(p), t_ref2m_min_inst_u(p))
      end if
      if (end_cd) then
        if (urbpoi(l)) then
          t_ref2m_max_u(p) = t_ref2m_max_inst_u(p)
          t_ref2m_min_u(p) = t_ref2m_min_inst_u(p)
          t_ref2m_max_inst_u(p) = -spval
          t_ref2m_min_inst_u(p) =  spval
        end if
      else if (secs <= int(dtsrf)) then
        t_ref2m_max_u(p) = spval
        t_ref2m_min_u(p) = spval
      end if
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', t_ref2m_r, kkincr)
    call extract_accum_field ('TREFAV_R', rbufslp, kkincr)
    do p = begp , endp
      l = plandunit(p)
      if (rbufslp(p) /= spval) then
        t_ref2m_max_inst_r(p) = max(rbufslp(p), t_ref2m_max_inst_r(p))
        t_ref2m_min_inst_r(p) = min(rbufslp(p), t_ref2m_min_inst_r(p))
      end if
      if (end_cd) then
        if (.not.(ifspecial(l))) then
          t_ref2m_max_r(p) = t_ref2m_max_inst_r(p)
          t_ref2m_min_r(p) = t_ref2m_min_inst_r(p)
          t_ref2m_max_inst_r(p) = -spval
          t_ref2m_min_inst_r(p) =  spval
        end if
      else if (secs <= int(dtsrf)) then
        t_ref2m_max_r(p) = spval
        t_ref2m_min_r(p) = spval
      end if
    end do

    ! Accumulate and extract T_VEG24 & T_VEG240 (heald 04/06)
    do p = begp , endp
      rbufslp(p) = t_veg(p)
    end do
    call update_accum_field  ('T_VEG24', rbufslp, kkincr)
    call extract_accum_field ('T_VEG24', t_veg24, kkincr)
    call update_accum_field  ('T_VEG240', rbufslp, kkincr)
    call extract_accum_field ('T_VEG240', t_veg240, kkincr)

    ! Accumulate and extract forc_solad24 & forc_solad240 (heald 04/06)
    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_solad(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp, kkincr)
    call extract_accum_field ('FSD240', fsd240, kkincr)
    call update_accum_field  ('FSD24', rbufslp, kkincr)
    call extract_accum_field ('FSD24', fsd24, kkincr)

    ! Accumulate and extract forc_solai24 & forc_solai240 (heald 04/06)
    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_solai(g,1)
    end do
    call update_accum_field  ('FSI24', rbufslp, kkincr)
    call extract_accum_field ('FSI24', fsi24, kkincr)
    call update_accum_field  ('FSI240', rbufslp, kkincr)
    call extract_accum_field ('FSI240', fsi240, kkincr)

    ! Accumulate and extract fsun24 & fsun240 (heald 04/06)
    do p = begp , endp
      rbufslp(p) = fsun(p)
    end do
    call update_accum_field  ('FSUN24', rbufslp, kkincr)
    call extract_accum_field ('FSUN24', fsun24, kkincr)
    call update_accum_field  ('FSUN240', rbufslp, kkincr)
    call extract_accum_field ('FSUN240', fsun240, kkincr)

    ! Accumulate and extract elai_p (heald 04/06)
    do p = begp , endp
      rbufslp(p) = elai(p)
    end do
    call update_accum_field  ('LAIP', rbufslp, kkincr)
    call extract_accum_field ('LAIP', elai_p, kkincr)

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    call update_accum_field  ('T10', t_ref2m, kkincr)
    call extract_accum_field ('T10', t10, kkincr)

#if (defined CNDV)
    ! Accumulate and extract TDA
    ! (accumulates TBOT as 30-day average)
    ! Also determine t_mo_min

    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_t(g)
    end do
    call update_accum_field  ('TDA', rbufslp, kkincr)
    call extract_accum_field ('TDA', rbufslp, kkincr)
    do p = begp , endp
      t_mo(p) = rbufslp(p)
      t_mo_min(p) = min(t_mo_min(p), rbufslp(p))
    end do

    ! Accumulate and extract PREC365
    ! (accumulates total precipitation as 365-day running mean)

    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC365', rbufslp, kkincr)
    call extract_accum_field ('PREC365', prec365, kkincr)

    ! Accumulate growing degree days based on 10-day running mean temperature.
    ! The trigger to reset the accumulated values to zero is -99999.

    ! Accumulate and extract AGDDTW (gdd base twmax, which is 23 deg C
    ! for boreal woody pfts)
    if ( date_is(nextdate,1,1) .and. time_is(nextdate,0) ) then
      do p = begp , endp
        rbufslp(p) = -99999.0_rk8
      end do
    else
      do p = begp , endp
        rbufslp(p) = max(0.0_rk8, (t10(p) - tfrz - twmax(ndllf_dcd_brl_tree)) &
                     * dtsrf/secspday)
      end do
    end if

    call update_accum_field  ('AGDDTW', rbufslp, kkincr)
    call extract_accum_field ('AGDDTW', agddtw, kkincr)

    ! Accumulate and extract AGDD

    do p = begp , endp
      rbufslp(p) = max(0.0_rk8, (t_ref2m(p) - (tfrz + 5.0_rk8)) &
            * dtsrf/secspday)
    end do
    call update_accum_field  ('AGDD', rbufslp, kkincr)
    call extract_accum_field ('AGDD', agdd, kkincr)
#endif

    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC60', rbufslp, kkincr)
    call extract_accum_field ('PREC60', prec60, kkincr)

    ! Accumulate and extract PREC10
    ! (accumulates total precipitation as 10-day running mean)
    do p = begp , endp
      g = pgridcell(p)
      rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC10', rbufslp, kkincr)
    call extract_accum_field ('PREC10', prec10, kkincr)

    if ( crop_prog ) then
      ! Accumulate and extract TDM10

      do p = begp , endp
        rbufslp(p) = min(t_ref2m_min(p),t_ref2m_min_inst(p)) !slevis: ok choice?
        if (rbufslp(p) > 1.e30_rk8) rbufslp(p) = tfrz !and were 'min'&
      end do                                       !'min_inst' not initialized?
      call update_accum_field  ('TDM10', rbufslp, kkincr)
      call extract_accum_field ('TDM10', a10tmin, kkincr)

      ! Accumulate and extract TDM5

      do p = begp , endp
        rbufslp(p) = min(t_ref2m_min(p),t_ref2m_min_inst(p)) !slevis: ok choice?
        if (rbufslp(p) > 1.e30_rk8) rbufslp(p) = tfrz !and were 'min'&
      end do    !'min_inst' not initialized?
      call update_accum_field  ('TDM5', rbufslp, kkincr)
      call extract_accum_field ('TDM5', a5tmin, kkincr)

      ! Accumulate and extract GDD0

      if ( date_is(nextdate,1,1) .and. time_is(nextdate,0) ) then
        do p = begp , endp
          rbufslp(p) = -99999.0_rk8
        end do
      else
        do p = begp , endp
          g = pgridcell(p)
          if ( ( month > 3 .and. month < 10 .and. latdeg(g) >= 0.0_rk8) .or. &
               ( (month > 9 .or.  month < 4) .and. latdeg(g) < 0.0_rk8) ) then
            rbufslp(p) = max(0.0_rk8, min(26.0_rk8, t_ref2m(p)-tfrz)) &
                             * dtsrf/secspday
          else
            ! keeps gdd unchanged at other times (eg, through Dec in NH)
            rbufslp(p) = 0.0_rk8
          end if
        end do
      end if
      call update_accum_field  ('GDD0', rbufslp, kkincr)
      call extract_accum_field ('GDD0', gdd0, kkincr)

      ! Accumulate and extract GDD8

      if ( date_is(nextdate,1,1) .and. time_is(nextdate,0) ) then
        do p = begp , endp
          rbufslp(p) = -99999.0_rk8
        end do
      else
        do p = begp , endp
          g = pgridcell(p)
          if ( ( month > 3 .and. month < 10 .and. latdeg(g) >= 0.0_rk8) .or. &
               ( (month > 9 .or.  month < 4) .and. latdeg(g) < 0.0_rk8) ) then
            rbufslp(p) = max(0.0_rk8, min(30.0_rk8, &
                         t_ref2m(p)-(tfrz + 8.0_rk8))) * dtsrf/secspday
          else
            ! keeps gdd unchanged at other times (eg, through Dec in NH)
            rbufslp(p) = 0.0_rk8
          end if
        end do
      end if
      call update_accum_field  ('GDD8', rbufslp, kkincr)
      call extract_accum_field ('GDD8', gdd8, kkincr)

      ! Accumulate and extract GDD10

      if ( date_is(nextdate,1,1) .and. time_is(nextdate,0) ) then
        do p = begp , endp
          rbufslp(p) = -99999.0_rk8
        end do
      else
        do p = begp , endp
          g = pgridcell(p)
          if ( ( month > 3 .and. month < 10 .and. latdeg(g) >= 0.0_rk8) .or. &
               ( (month > 9 .or.  month < 4) .and. latdeg(g) < 0.0_rk8) ) then
            rbufslp(p) = max(0.0_rk8, min(30.0_rk8, &
                         t_ref2m(p)-(tfrz + 10.0_rk8))) * dtsrf/secspday
          else
            ! keeps gdd unchanged at other times (eg, through Dec in NH)
            rbufslp(p) = 0.0_rk8
          end if
        end do
      end if
      call update_accum_field  ('GDD10', rbufslp, kkincr)
      call extract_accum_field ('GDD10', gdd10, kkincr)

      ! Accumulate and extract GDDPLANT

      do p = begp , endp
        if (croplive(p)) then ! relative to planting date
          itypveg = itype(p)
          rbufslp(p) = max(0.0_rk8, min(mxtmp(itypveg), &
                           t_ref2m(p)-(tfrz + baset(itypveg)))) &
                          * dtsrf/secspday
          if ( itypveg == nwcereal .or. itypveg == nwcerealirrig) then
            rbufslp(p) = rbufslp(p)*vf(p)
          end if
        else
          rbufslp(p) = -99999.0_rk8
        end if
      end do
      call update_accum_field  ('GDDPLANT', rbufslp, kkincr)
      call extract_accum_field ('GDDPLANT', gddplant, kkincr)

      ! Accumulate and extract GDDTSOI
      ! In agroibis this variable is calculated
      ! to 0.05 m, so here we use the top two soil layers

      do p = begp , endp
        if (croplive(p)) then ! relative to planting date
          itypveg = itype(p)
          c = pcolumn(p)
          rbufslp(p) = max(0.0_rk8, min(mxtmp(itypveg), &
             ((t_soisno(c,1)*dz(c,1)+t_soisno(c,2)*dz(c,2))/ &
             (dz(c,1)+dz(c,2))) - (tfrz + baset(itypveg)))) * dtsrf/secspday
          if (itypveg == nwcereal .or. itypveg == nwcerealirrig) then
            rbufslp(p) = rbufslp(p)*vf(p)
          end if
        else
          rbufslp(p) = -99999.0_rk8
        end if
      end do
      call update_accum_field  ('GDDTSOI', rbufslp, kkincr)
      call extract_accum_field ('GDDTSOI', gddtsoi, kkincr)

    end if

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds
  !
  ! Initialize clmtype variables that are associated with
  ! time accumulated fields. This routine is called in an initial run
  ! at nstep=0 for cam and csm mode.
  ! This routine is also always called for a restart run and
  ! therefore must be called after the restart file is read in
  ! and the accumulated fields are obtained.
  !
  subroutine initAccClmtype
    implicit none
    ! daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min(:)
    ! daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max(:)
    ! instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst(:)
    ! instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst(:)
    ! Urban daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min_u(:)
    ! Rural daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_min_r(:)
    ! Urban daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max_u(:)
    ! Rural daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer :: t_ref2m_max_r(:)
    ! Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst_u(:)
    ! Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_min_inst_r(:)
    ! Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst_u(:)
    ! Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer :: t_ref2m_max_inst_r(:)
    ! 10-day running mean of the 2 m temperature (K)
    real(rk8), pointer :: t10(:)
#ifdef CNDV
    ! 30-day average temperature (Kelvin)
    real(rk8), pointer :: t_mo(:)
    ! 365-day running mean of tot. precipitation
    real(rk8), pointer :: prec365(:)
    ! accumulated growing degree days above twmax
    real(rk8), pointer :: agddtw(:)
    ! accumulated growing degree days above 5
    real(rk8), pointer :: agdd(:)
#endif
    ! 60-day running mean of tot. precipitation
    real(rk8), pointer :: prec60(:)
    ! 10-day running mean of tot. precipitation
    real(rk8), pointer :: prec10(:)
    ! growing degree-days base 0C'
    real(rk8), pointer :: gdd0(:)
    ! growing degree-days base 8C from planting
    real(rk8), pointer :: gdd8(:)
    ! growing degree-days base 10C from planting
    real(rk8), pointer :: gdd10(:)
    ! growing degree-days from planting
    real(rk8), pointer :: gddplant(:)
    ! growing degree-days from planting (top two soil layers)
    real(rk8), pointer :: gddtsoi(:)
    ! 10-day running mean of min 2-m temperature
    real(rk8), pointer :: a10tmin(:)
    ! 5-day running mean of min 2-m temperature
    real(rk8), pointer :: a5tmin(:)
    ! heald (04/06): accumulated variables for VOC emissions
    real(rk8), pointer :: t_veg24(:) ! 24hr average vegetation temperature (K)
    ! 240hr average vegetation temperature (Kelvin)
    real(rk8), pointer :: t_veg240(:)
    real(rk8), pointer :: fsd24(:)   ! 24hr average of direct beam radiation
    real(rk8), pointer :: fsd240(:)  ! 240hr average of direct beam radiation
    real(rk8), pointer :: fsi24(:)   ! 24hr average of diffuse beam radiation
    real(rk8), pointer :: fsi240(:)  ! 240hr average of diffuse beam radiation
    real(rk8), pointer :: fsun24(:)  ! 24hr average of sunlit fraction of canopy
    ! 240hr average of sunlit fraction of canopy
    real(rk8), pointer :: fsun240(:)
    real(rk8), pointer :: elai_p(:)  ! leaf area index average over timestep
    integer(ik4) :: p          ! indices
    integer(ik4) :: ier        ! error status
    integer(ik4) :: begp, endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(rk8), pointer :: rbufslp(:)  ! temporary
    character(len=32) :: subname = 'initAccClmtype'  ! subroutine name
    integer(ik8) :: kkincr

    ! Assign local pointers to derived subtypes components (pft-level)

    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_ref2m_max_inst_u => clm3%g%l%c%p%pes%t_ref2m_max_inst_u
    t_ref2m_max_inst_r => clm3%g%l%c%p%pes%t_ref2m_max_inst_r
    t_ref2m_min_inst_u => clm3%g%l%c%p%pes%t_ref2m_min_inst_u
    t_ref2m_min_inst_r => clm3%g%l%c%p%pes%t_ref2m_min_inst_r
    t_ref2m_max_u      => clm3%g%l%c%p%pes%t_ref2m_max_u
    t_ref2m_max_r      => clm3%g%l%c%p%pes%t_ref2m_max_r
    t_ref2m_min_u      => clm3%g%l%c%p%pes%t_ref2m_min_u
    t_ref2m_min_r      => clm3%g%l%c%p%pes%t_ref2m_min_r
    t10              => clm3%g%l%c%p%pes%t10
    a10tmin          => clm3%g%l%c%p%pes%a10tmin
    a5tmin           => clm3%g%l%c%p%pes%a5tmin
#if (defined CNDV)
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd
#endif
    prec60           => clm3%g%l%c%p%pps%prec60
    prec10           => clm3%g%l%c%p%pps%prec10
    gdd0             => clm3%g%l%c%p%pps%gdd0
    gdd8             => clm3%g%l%c%p%pps%gdd8
    gdd10            => clm3%g%l%c%p%pps%gdd10
    gddplant         => clm3%g%l%c%p%pps%gddplant
    gddtsoi          => clm3%g%l%c%p%pps%gddtsoi
    ! heald (04/06): accumulated variables for VOC emissions
    t_veg24          => clm3%g%l%c%p%pvs%t_veg24
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240
    fsd24            => clm3%g%l%c%p%pvs%fsd24
    fsd240           => clm3%g%l%c%p%pvs%fsd240
    fsi24            => clm3%g%l%c%p%pvs%fsi24
    fsi240           => clm3%g%l%c%p%pvs%fsi240
    fsun24           => clm3%g%l%c%p%pvs%fsun24
    fsun240          => clm3%g%l%c%p%pvs%fsun240
    elai_p           => clm3%g%l%c%p%pvs%elai_p

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Determine time step

    kkincr = ktau / ntsrf

    ! Initialize 2m ref temperature max and min values

    if ( nsrest == nsrStartup ) then
      ! Why not restart? These vars are not in clmr.
      do p = begp , endp
        t_ref2m_max(p) = spval
        t_ref2m_min(p) = spval
        t_ref2m_max_inst(p) = -spval
        t_ref2m_min_inst(p) =  spval
        t_ref2m_max_u(p) = spval
        t_ref2m_max_r(p) = spval
        t_ref2m_min_u(p) = spval
        t_ref2m_min_r(p) = spval
        t_ref2m_max_inst_u(p) = -spval
        t_ref2m_max_inst_r(p) = -spval
        t_ref2m_min_inst_u(p) =  spval
        t_ref2m_min_inst_r(p) =  spval
      end do
    end if

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
      write(stderr,*) &
        'extract_accum_hist allocation error for rbufslp in '//subname
       call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('T_VEG24', rbufslp, kkincr)
    do p = begp , endp
      t_veg24(p) = rbufslp(p)
    end do

    call extract_accum_field ('T_VEG240', rbufslp, kkincr)
    do p = begp , endp
      t_veg240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD24', rbufslp, kkincr)
    do p = begp , endp
      fsd24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD240', rbufslp, kkincr)
    do p = begp , endp
      fsd240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI24', rbufslp, kkincr)
    do p = begp , endp
      fsi24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI240', rbufslp, kkincr)
    do p = begp , endp
      fsi240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN24', rbufslp, kkincr)
    do p = begp , endp
      fsun24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN240', rbufslp, kkincr)
    do p = begp , endp
      fsun240(p) = rbufslp(p)
    end do

    call extract_accum_field ('LAIP', rbufslp, kkincr)
    do p = begp , endp
      elai_p(p) = rbufslp(p)
    end do

    if ( crop_prog )then

      call extract_accum_field ('GDD0', rbufslp, kkincr)
      do p = begp , endp
        gdd0(p) = rbufslp(p)
      end do

      call extract_accum_field ('GDD8', rbufslp, kkincr)
      do p = begp , endp
        gdd8(p) = rbufslp(p)
      end do

      call extract_accum_field ('GDD10', rbufslp, kkincr)
      do p = begp , endp
        gdd10(p) = rbufslp(p)
      end do

      call extract_accum_field ('GDDPLANT', rbufslp, kkincr)
      do p = begp , endp
        gddplant(p) = rbufslp(p)
      end do

      call extract_accum_field ('GDDTSOI', rbufslp, kkincr)
      do p = begp , endp
        gddtsoi(p) = rbufslp(p)
      end do

      call extract_accum_field ('TDM10', rbufslp, kkincr)
      do p = begp , endp
        a10tmin(p) = rbufslp(p)
      end do

      call extract_accum_field ('TDM5', rbufslp, kkincr)
      do p = begp , endp
        a5tmin(p) = rbufslp(p)
      end do

    end if

    call extract_accum_field ('T10', rbufslp, kkincr)
    do p = begp , endp
      t10(p) = rbufslp(p)
    end do

#if (defined CNDV)

    call extract_accum_field ('TDA', rbufslp, kkincr)
    do p = begp , endp
      t_mo(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC365', rbufslp, kkincr)
    do p = begp , endp
      prec365(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDDTW', rbufslp, kkincr)
    do p = begp , endp
      agddtw(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD', rbufslp, kkincr)
    do p = begp , endp
      agdd(p) = rbufslp(p)
    end do

#endif
    call extract_accum_field ('PREC60', rbufslp, kkincr)
    do p = begp , endp
      prec60(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC10', rbufslp, kkincr)
    do p = begp , endp
      prec10(p) = rbufslp(p)
    end do

    deallocate(rbufslp)

  end subroutine initAccClmtype

end module mod_clm_accflds
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
