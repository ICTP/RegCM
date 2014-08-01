module mod_clm_initslake
  !
  ! Contains time constant and time-varying initialization code for S Lake
  ! scheme.
  ! Note that this is called after the restart file is read to facilitate
  ! using existing spin-up files that lack the new lake fields.
  ! A check is done to see whether mkarbinit was called and/or state
  ! variable values are missing (spval).
  ! Once this version is officialized, this is no longer necessary.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_clm_type
  use mod_clm_varpar , only : nlevsno , nlevlak , nlevgrnd , nlevsoi
  use mod_clm_varcon , only : bdsno , denice , denh2o , spval , tfrz , tkwat
  use mod_clm_varcon , only : istdlak , zlak , dzlak , zsoi , dzsoi , zisoi
  use mod_clm_varcon , only : denh2o , denice , secspday , spval
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_snicar , only : snw_rds_min
  use mod_clm_slakecon , only : lsadz , fcrit , minz0lake , depthcrit
  use mod_clm_slakecon , only : mixfact , pudz , lakepuddling
  use mod_clm_slakecon , only : lake_puddle_thick , deepmixing_depthcrit
  use mod_clm_slakecon , only : deepmixing_mixfact , lake_use_old_fcrit_minz0
  use mod_clm_slakecon , only : lake_melt_icealb , alblakwi
  use mod_clm_decomp , only : get_proc_bounds , get_proc_global
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varctl , only : nsrest

  implicit none

  private

  save

  public :: initSLake ! driver

  private :: initTimeConst        ! Set constant parameters (and h2osoi_vol).
  private :: makearbinit          ! Set time-variable parameters for spin up.
  private :: snow_depth2levLake   ! Reset snow layers over S lakes if necessary.

  contains
  !
  ! Calls initTimeConst.
  ! Calls makearbinit with logical arbinit. If arbinit == .true. OR
  ! if initial conditions file does not contain SLake soil and snow state,
  ! then initializes time varying values. This allows back-compatibility
  ! with initial condition files that have not been spun up with the new
  ! lake code. In future versions, this could be phased out.
  !
  subroutine initSLake( arbinit )
    implicit none
    logical , intent(in) :: arbinit ! Whether mkarbinit has been called.
    call initTimeConst()
    ! Attn EK
    ! For now
    call makearbinit(arbinit)
    ! For future versions always using initial condition files spun up
    ! with the new lake code: if (arbinit) call makearbinit(arbinit)
  end subroutine initSLake
  !
  ! If arbinit == .true., or if soil and snow have not been initialized,
  ! then sets time varying values.
  ! Initializes the following time varying variables:
  ! water      : h2osoi_liq, h2osoi_ice, h2osoi_vol
  ! snow       : snl, dz, z, zi
  ! temperature: t_soisno
  ! h2osno, snow_depth, can water var, t_lake, & t_grnd should have been
  ! initialized already from regular mkarbinit.
  !
  subroutine makearbinit( arbinit )
    implicit none
    logical , intent(in) :: arbinit ! Whether mkarbinit has been called.
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! layer thickness depth (m)
    real(rk8) , pointer , dimension(:,:) :: dz
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: watsat
    ! ice lens (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    ! liquid water (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    ! lake temperature (Kelvin)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: t_lake
    ! ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    ! mass fraction of lake layer that is frozen
    real(rk8) , pointer , dimension(:,:) :: lake_icefrac
    ! top level eddy conductivity (W/mK)
    real(rk8) , pointer , dimension(:) :: savedtke1
    ! friction velocity (m/s)
    real(rk8) , pointer , dimension(:) :: ust_lake
    !roughness length over ground, momentum [m]
    real(rk8) , pointer , dimension(:) :: z0mg
    ! New SNICAR variables
    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(rk8) , pointer , dimension(:,:) :: snw_rds
    ! snow grain size, top (col) [microns]
    real(rk8) , pointer , dimension(:) :: snw_rds_top
    ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(rk8) , pointer , dimension(:) :: sno_liq_top
    integer(ik4) :: j , l , c   ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending ldunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    logical :: localarbinit

    localarbinit = arbinit

    if ( myid == italk ) then
      write (stdout,*) 'Setting initial data to non-spun up values for '// &
              'lake points'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype
    lakpoi     => clm3%g%l%lakpoi

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    dz         => clm3%g%l%c%cps%dz
    watsat     => clm3%g%l%c%cps%watsat
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
    h2osno     => clm3%g%l%c%cws%h2osno
    t_soisno   => clm3%g%l%c%ces%t_soisno
    t_lake     => clm3%g%l%c%ces%t_lake
    t_grnd     => clm3%g%l%c%ces%t_grnd
    lake_icefrac => clm3%g%l%c%cws%lake_icefrac
    savedtke1  => clm3%g%l%c%cps%savedtke1
    ust_lake   => clm3%g%l%c%cps%ust_lake
    z0mg       => clm3%g%l%c%cps%z0mg
    ! New SNICAR variables
    snw_rds          => clm3%g%l%c%cps%snw_rds
    snw_rds_top      => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top      => clm3%g%l%c%cps%sno_liq_top

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, snow_depth and snowage has valid values everywhere
    ! canopy water (pft level)

    ! Note (Attn EK): h2osno, snow_depth, and t_lake are currently set
    ! in mkarbinit. Those could be moved here for modularity when old
    ! lake model is removed.

    ! Set snow layer number, depth and thickness

    ! Reset snow layers over lakes for S lake code.
    call snow_depth2levLake(begc,endc,arbinit)
    ! This will now ONLY run if arbinit, or initial condition file appears
    ! to be coming from old lake model, to preserve bit-for-bit results
    ! during restarts.

    ! Set snow/soil temperature, note:
    ! t_soisno has valid values over non-lakes and S lakes.
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    do c = begc , endc

      l = clandunit(c)

      if ( lakpoi(l) ) then
        if (snl(c) < 0) then
          do j = snl(c)+1 , 0
            if ( arbinit .or. t_soisno(c,j) == spval ) t_soisno(c,j) = 250.D0
          end do
        end if
        do j = 1 , nlevgrnd
          if ( arbinit .or. t_soisno(c,j) == spval ) &
            t_soisno(c,j) = t_lake(c,nlevlak)
        end do
        do j = 1 , nlevlak
          if ( arbinit .or. lake_icefrac(c,j) == spval ) then
            ! initialization of SNICAR vars is fragile;
            ! set localarbinit = .true. if any layer is spval
            if ( lake_icefrac(c,j)==spval ) localarbinit = .true.
            lake_icefrac(c,j) = 0.D0
            ! Always initialize with no ice to prevent excessive
            ! ice sheets from forming when
            ! starting with old lake model that has unrealistically
            ! cold lake temperatures.
            ! Keep lake temperature as is, and the energy deficit
            ! below freezing (which is no smaller than it would have
            ! been with prognostic ice, as the temperature would then
            ! have been higher and more heat would have flowed out of
            ! the lake) will be converted to ice in the first timestep.
          end if
        end do
        if ( arbinit .or. savedtke1(c) == spval ) savedtke1(c) = tkwat
        ! For prognostic roughness length
        if ( arbinit .or. ust_lake(c) == spval ) then
          ust_lake(c) = 0.1D0
          z0mg(c) = 0.0004D0
        end if
      end if
    end do

    do j = 1 , nlevgrnd
      do c = begc , endc
        l = clandunit(c)
        if ( lakpoi(l) ) then
          if ( arbinit .or. h2osoi_vol(c,j) == spval ) then
            if ( j <= nlevsoi ) then ! soil
              h2osoi_vol(c,j) = watsat(c,j)
              ! Reset h2osoi_liq & h2osoi_ice
              h2osoi_liq(c,j) = spval
              h2osoi_ice(c,j) = spval
            else ! bedrock
              h2osoi_vol(c,j) = 0.D0
            end if
          end if
          ! soil layers
          if ( arbinit .or. h2osoi_ice(c,j) == spval .or. &
               h2osoi_liq(c,j) == spval) then
            if ( t_soisno(c,j) <= tfrz ) then
              h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
              h2osoi_liq(c,j) = 0.D0
            else
              h2osoi_ice(c,j) = 0.D0
              h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
            endif
          end if
        end if
      end do
    end do

    ! Set snow
    ! If localarbinit only!

    do j = -nlevsno+1 , 0
      do c = begc , endc
        l = clandunit(c)
        if ( lakpoi(l) ) then
          if ( j > snl(c) ) then
            if ( arbinit .or. h2osoi_ice(c,j) == spval) then
              h2osoi_ice(c,j) = dz(c,j)*bdsno
            end if
            if ( arbinit .or. h2osoi_liq(c,j) == spval ) then
              h2osoi_liq(c,j) = 0.D0
            end if
          end if
        end if
      end do
    end do

    do c = begc , endc
      l = clandunit(c)
      if ( lakpoi(l) .and. snl(c) < 0 .and. localarbinit ) then
        ! SNICAR fields not initialized in regular mkarbinit b/c there
        ! lakes have snl==0
        ! The rest of the fields will be initialized to zero in any case
        snw_rds(c,snl(c)+1:0)        = snw_rds_min
        snw_rds(c,-nlevsno+1:snl(c)) = 0.D0
        snw_rds_top(c)               = snw_rds_min
        sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / &
                (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
      end if
    end do
  end subroutine makearbinit
  !
  ! Create snow layers and interfaces given snow depth for lakes.
  ! Note that cps%zi(0) is set in routine iniTimeConst.
  !
  subroutine snow_depth2levLake(lbc,ubc,arbinit)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! column bounds
    logical , intent(in) :: arbinit
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! snow height (m)
    real(rk8) , pointer , dimension(:) :: snow_depth
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! mass fraction of lake layer that is frozen
    real(rk8) , pointer , dimension(:,:) :: lake_icefrac
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! layer depth  (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: z
    ! layer thickness depth (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: dz
    ! interface depth (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: zi
    integer(ik4) :: c , l , j      !indices

    ! Assign local pointers to derived subtypes components (landunit-level)

    lakpoi => clm3%g%l%lakpoi

    ! Assign local pointers to derived type members (column-level)

    clandunit => clm3%g%l%c%landunit
    snow_depth    => clm3%g%l%c%cps%snow_depth
    snl       => clm3%g%l%c%cps%snl
    zi        => clm3%g%l%c%cps%zi
    dz        => clm3%g%l%c%cps%dz
    z         => clm3%g%l%c%cps%z
    lake_icefrac => clm3%g%l%c%cws%lake_icefrac


    ! Determine snow levels and interfaces for lake points
    ! lsadz (default 0.03m) is added to min & max thicknesses to not
    ! stress lake numerics at 30 min. timestep.

    do c = lbc , ubc
      l = clandunit(c)
      if ( lakpoi(l) .and. (arbinit .or. lake_icefrac(c,1) == spval) ) then
        if ( snow_depth(c) < 0.01D0 + lsadz ) then
          snl(c) = 0
          dz(c,-nlevsno+1:0) = 0.D0
          z (c,-nlevsno+1:0) = 0.D0
          zi(c,-nlevsno+0:0) = 0.D0
        else
          if ( (snow_depth(c) >= 0.01D0 + lsadz) .and. &
               (snow_depth(c) <= 0.03D0 + lsadz) ) then
            snl(c) = -1
            dz(c,0)  = snow_depth(c)
          else if ( (snow_depth(c) > 0.03D0 + 2.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.04D0 + 2.D0*lsadz) ) then
            snl(c) = -2
            dz(c,-1) = snow_depth(c)/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.04D0 + 2.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.07D0 + 2.D0*lsadz) ) then
            snl(c) = -2
            dz(c,-1) = 0.02D0 + lsadz
            dz(c, 0) = snow_depth(c) - dz(c,-1)
          else if ( (snow_depth(c) > 0.07D0 + 3.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.12D0 + 3.D0*lsadz) ) then
            snl(c) = -3
            dz(c,-2) = 0.02D0 + lsadz
            dz(c,-1) = (snow_depth(c) - 0.02D0 - lsadz)/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.12D0 + 3.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.18D0 + 3.D0*lsadz) ) then
            snl(c) = -3
            dz(c,-2) = 0.02D0 + lsadz
            dz(c,-1) = 0.05D0 + lsadz
            dz(c, 0) = snow_depth(c) - dz(c,-2) - dz(c,-1)
          else if ( (snow_depth(c) > 0.18D0 + 4.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.29D0 + 4.D0*lsadz) ) then
            snl(c) = -4
            dz(c,-3) = 0.02D0 + lsadz
            dz(c,-2) = 0.05D0 + lsadz
            dz(c,-1) = (snow_depth(c) - dz(c,-3) - dz(c,-2))/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.29D0 + 4.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.41D0 + 4.D0*lsadz) ) then
            snl(c) = -4
            dz(c,-3) = 0.02D0 + lsadz
            dz(c,-2) = 0.05D0 + lsadz
            dz(c,-1) = 0.11D0 + lsadz
            dz(c, 0) = snow_depth(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
          else if ( (snow_depth(c) > 0.41D0 + 5.D0*lsadz) .and. &
                    (snow_depth(c) <= 0.64D0 + 5.D0*lsadz) ) then
            snl(c) = -5
            dz(c,-4) = 0.02D0 + lsadz
            dz(c,-3) = 0.05D0 + lsadz
            dz(c,-2) = 0.11D0 + lsadz
            dz(c,-1) = (snow_depth(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( snow_depth(c) > 0.64D0 + 5.D0*lsadz ) then
            snl(c) = -5
            dz(c,-4) = 0.02D0 + lsadz
            dz(c,-3) = 0.05D0 + lsadz
            dz(c,-2) = 0.11D0 + lsadz
            dz(c,-1) = 0.23D0 + lsadz
            dz(c, 0)=snow_depth(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
          end if
        end if
      end if
    end do

    ! The following loop is currently not vectorized
    do c = lbc , ubc
      l = clandunit(c)
      if ( lakpoi(l) .and. (arbinit .or. lake_icefrac(c,1) == spval) ) then
        do j = 0 , snl(c)+1 , -1
          z(c,j)    = zi(c,j) - 0.5D0*dz(c,j)
          zi(c,j-1) = zi(c,j) - dz(c,j)
        end do
      end if
    end do
  end subroutine snow_depth2levLake
  !
  ! Initialize time invariant clm variables for S Lake code (and h2osoi_vol).
  !
  subroutine initTimeConst
    implicit none
    ! landunit index of column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! gridcell index of column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! landunit type index
    integer(ik4) , pointer , dimension(:) :: ltype
    real(rk8) , pointer , dimension(:,:) :: cellsand  ! column 3D sand
    real(rk8) , pointer , dimension(:,:) :: cellclay  ! column 3D clay
    real(rk8) , pointer , dimension(:,:) :: cellorg   ! column 3D org
    ! liquid water (kg/m2) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    ! ice lens (kg/m2) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    ! variable lake depth (m)
    real(rk8) , pointer , dimension(:) :: lakedepth
    ! layer depth (m)
    real(rk8) , pointer , dimension(:,:) :: z
    ! interface level below a "z" level (m)
    real(rk8) , pointer , dimension(:,:) :: zi
    ! layer thickness depth (m)
    real(rk8) , pointer , dimension(:,:) :: dz
    ! layer thickness for lake (m)
    real(rk8) , pointer , dimension(:,:) :: dz_lake
    ! layer depth for lake (m)
    real(rk8) , pointer , dimension(:,:) :: z_lake
    ! Clapp and Hornberger "b" (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: bsw
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: watsat
    ! btran parameter for btran=0
    real(rk8) , pointer , dimension(:,:) :: watdry
    ! btran parameter for btran = 1
    real(rk8) , pointer , dimension(:,:) :: watopt
    ! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: hksat
    ! minimum soil suction (mm) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: sucsat
    ! heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: csol
    ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tkmg
    ! thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tkdry
    ! thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tksatu
    ! volumetric soil water at field capacity (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: watfc
    !volumetric soil water [m3/m3]  (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    integer(ik4) :: lev
#ifdef EXTRALAKELAYERS
    integer(ik4) :: i
#endif
    integer(ik4) :: g , l , c  ! indices
    real(rk8) :: bd    ! bulk density of dry soil material [kg/m^3]
    real(rk8) :: tkm   ! mineral conductivity
    real(rk8) :: xksat ! maximum hydraulic conductivity of soil [mm/s]
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numg   ! total number of gridcells across all processors
    integer(ik4) :: numl   ! total number of landunits across all processors
    integer(ik4) :: numc   ! total number of columns across all processors
    integer(ik4) :: nump   ! total number of pfts across all processors
    real(rk8) :: depthratio ! ratio of lake depth to standard deep lake depth
    real(rk8) :: clay    ! temporary
    real(rk8) :: sand    ! temporary
    ! New CLM 4 variables
    real(rk8) :: om_frac            ! organic matter fraction
    real(rk8) :: om_watsat = 0.9D0  ! porosity of organic soil
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(rk8) :: om_hksat  = 0.1D0
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(rk8) :: om_tkm    = 0.25D0
    ! saturated suction for organic matter (Letts, 2000)
    real(rk8) :: om_sucsat = 10.3D0
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(rk8) :: om_csol = 2.5D0
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(rk8) :: om_tkd = 0.05D0
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(rk8) :: om_b = 2.7D0
    ! organic matter (kg/m3) where soil is assumed to act like peat
    real(rk8) :: organic_max = 130.D0
    ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(rk8) :: csol_bedrock = 2.0D6
    ! percolation threshold
    real(rk8) :: pc = 0.5D0
    ! percolation exponent
    real(rk8) :: pcbeta = 0.139D0
    ! Note: These constants should be shared with iniTimeConst.
    real(rk8) :: perc_frac   ! "percolating" fraction of organic soil
    real(rk8) :: perc_norm   ! normalize to 1 when 100% organic soil
    real(rk8) :: uncon_hksat ! series conductivity of mineral/organic soil
    real(rk8) :: uncon_frac  ! fraction of "unconnected" soil

    ! timestep used for nominal lsadz
    real(rk8) , parameter :: dtime_lsadz = 1800.D0
    ! nominal minimum snow layer thickness
    ! used in SnowHydrology, etc. This should be a global constant.
    real(rk8) , parameter :: snodzmin = 0.01D0

    if ( myid == italk ) then
      write (stdout,*) 'Attempting to initialize time invariant '// &
              'variables for lakes'
    end if
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call get_proc_global(numg,numl,numc,nump)

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype           => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit       => clm3%g%l%c%landunit
    cgridcell       => clm3%g%l%c%gridcell
    z               => clm3%g%l%c%cps%z
    dz              => clm3%g%l%c%cps%dz
    zi              => clm3%g%l%c%cps%zi
    bsw             => clm3%g%l%c%cps%bsw
    watsat          => clm3%g%l%c%cps%watsat
    watdry          => clm3%g%l%c%cps%watdry
    watopt          => clm3%g%l%c%cps%watopt
    hksat           => clm3%g%l%c%cps%hksat
    sucsat          => clm3%g%l%c%cps%sucsat
    tkmg            => clm3%g%l%c%cps%tkmg
    tksatu          => clm3%g%l%c%cps%tksatu
    tkdry           => clm3%g%l%c%cps%tkdry
    csol            => clm3%g%l%c%cps%csol
    dz_lake         => clm3%g%l%c%cps%dz_lake
    z_lake          => clm3%g%l%c%cps%z_lake
    cellsand        => clm3%g%l%c%cps%cellsand
    cellclay        => clm3%g%l%c%cps%cellclay
    cellorg         => clm3%g%l%c%cps%cellorg
    lakedepth       => clm3%g%l%c%cps%lakedepth
    watfc           => clm3%g%l%c%cps%watfc
    h2osoi_liq      => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice      => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_vol      => clm3%g%l%c%cws%h2osoi_vol


    ! Set SLakeCon constants according to namelist fields
    if ( lake_use_old_fcrit_minz0 ) then
      ! critical dimensionless fetch for Charnock parameter.
      ! From Vickers & Mahrt 1997 but converted to use u instead of u*
      ! (Form used in Subin et al. 2011)
      fcrit   = 22.D0
      ! (m) Minimum allowed roughness length for unfrozen lakes.
      ! (Used in Subin et al. 2011)
      minz0lake = 1.D-5
    else
      fcrit   = 100.D0  ! Vickers & Mahrt 1997
      ! (m) Minimum allowed roughness length for unfrozen lakes.
      ! Now set low so it is only to avoid floating point exceptions.
      minz0lake = 1.D-10
    end if

    if ( lakepuddling ) then
      ! (m) Minimum total ice thickness required to allow lake puddling.
      ! Default is 0.2m.
      ! This option has not been extensively tested.
      ! This option turns on lake_no_melt_icealb, as the decrease in albedo
      ! will be based on whether there is water over nice, not purely a
      ! function of ice top temperature.
      pudz = lake_puddle_thick
    end if

    ! (m) Depth beneath which to increase mixing.
    ! See discussion in Subin et al. 2011
    depthcrit = deepmixing_depthcrit
    ! Mixing increase factor. ! Defaults are 25 m, increase by 10.
    ! Note some other namelists will be used directly in lake physics
    ! during model integration.
    mixfact = deepmixing_mixfact

    ! Adjust lsadz (added to min. snow layer thickness) if dtsrf /= 1800 s.
    if ( dtsrf /= dtime_lsadz ) then
      lsadz = (lsadz + snodzmin) * (dtsrf/dtime_lsadz)**0.5D0  - snodzmin
      lsadz = max(lsadz, 0.D0)
      if ( myid == italk ) then
        write(stdout,*) 'CLM timestep appears to differ from 1800s. '
        write(stdout,*) 'Adjusting minimum snow layer thickness over lakes '// &
                'to maintain stability.'
        write(stdout,*) 'New minimum thickness is ', lsadz + snodzmin, '.'
      end if
    end if

    ! Set alblakwi
    alblakwi(:) = lake_melt_icealb(:)


#ifndef EXTRALAKELAYERS

    ! Lake layers

    dzlak(1) = 0.1D0
    dzlak(2) = 1.D0
    dzlak(3) = 2.D0
    dzlak(4) = 3.D0
    dzlak(5) = 4.D0
    dzlak(6) = 5.D0
    dzlak(7) = 7.D0
    dzlak(8) = 7.D0
    dzlak(9) = 10.45D0
    dzlak(10)= 10.45D0

    zlak(1) =  0.05D0
    zlak(2) =  0.6D0
    zlak(3) =  2.1D0
    zlak(4) =  4.6D0
    zlak(5) =  8.1D0
    zlak(6) = 12.6D0
    zlak(7) = 18.6D0
    zlak(8) = 25.6D0
    zlak(9) = 34.325D0
    zlak(10)= 44.775D0
#else
    dzlak(1) =0.1D0
    dzlak(2) =0.25D0
    dzlak(3) =0.25D0
    dzlak(4) =0.25D0
    dzlak(5) =0.25D0
    dzlak(6) =0.5D0
    dzlak(7) =0.5D0
    dzlak(8) =0.5D0
    dzlak(9) =0.5D0
    dzlak(10) =0.75D0
    dzlak(11) =0.75D0
    dzlak(12) =0.75D0
    dzlak(13) =0.75D0
    dzlak(14) =2D0
    dzlak(15) =2D0
    dzlak(16) =2.5D0
    dzlak(17) =2.5D0
    dzlak(18) =3.5D0
    dzlak(19) =3.5D0
    dzlak(20) =3.5D0
    dzlak(21) =3.5D0
    dzlak(22) =5.225D0
    dzlak(23) =5.225D0
    dzlak(24) =5.225D0
    dzlak(25) =5.225D0

    zlak(1) = dzlak(1)/2.D0
    do i = 2 , nlevlak
      zlak(i) = zlak(i-1) + (dzlak(i-1)+dzlak(i))/2.D0
    end do
#endif
    ! --------------------------------------------------------------------
    ! Initialize soil and lake levels
    ! Initialize soil thermal and hydraulic properties [color is not needed]
    ! --------------------------------------------------------------------

    ! Column level initialization
    do c = begc , endc

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)

      ! Soil hydraulic and thermal properties
      if ( ltype(l)==istdlak ) then
        do lev = 1 , nlevgrnd
          if ( lev <= nlevsoi )then
            clay    = cellclay(c,lev)
            sand    = cellsand(c,lev)
            om_frac = (cellorg(c,lev)/organic_max)**2.D0
          else
            clay    = cellclay(c,nlevsoi)
            sand    = cellsand(c,nlevsoi)
            om_frac = 0.0D0
          end if
          watsat(c,lev) = 0.489D0 - 0.00126D0*sand
          bsw(c,lev)    = 2.91 + 0.159*clay
          sucsat(c,lev) = 10.D0 * ( 10.D0**(1.88D0-0.0131D0*sand) )
          bd            = (1.D0-watsat(c,lev))*2.7D3
          watsat(c,lev) = (1.D0 - om_frac)*watsat(c,lev) + om_watsat*om_frac
          tkm           = (1.D0-om_frac)*(8.80D0*sand+2.92D0*clay) / &
                  (sand+clay)+om_tkm*om_frac ! W/(m K)
          bsw(c,lev)    = (1.D0-om_frac)*(2.91D0 + 0.159D0*clay) + &
                  om_frac*om_b
          sucsat(c,lev) = (1.D0-om_frac)*sucsat(c,lev) + om_sucsat*om_frac
          xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

          ! perc_frac is zero unless perf_frac greater than percolation
          ! threshold
          if ( om_frac > pc ) then
            perc_norm=(1.D0 - pc)**(-pcbeta)
            perc_frac=perc_norm*(om_frac - pc)**pcbeta
          else
            perc_frac=0.D0
          endif
          ! uncon_frac is fraction of mineral soil plus fraction of
          ! "nonpercolating" organic soil
          uncon_frac = (1.D0-om_frac)+(1.D0-perc_frac)*om_frac
          ! uncon_hksat is series addition of mineral/organic conductivites
          if ( om_frac .lt. 1.D0 ) then
            uncon_hksat=uncon_frac/((1.D0-om_frac)/xksat &
                   +((1.D0-perc_frac)*om_frac)/om_hksat)
          else
            uncon_hksat = 0.D0
          end if
          hksat(c,lev) = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat
          tkmg(c,lev) = tkm ** (1.D0- watsat(c,lev))
          tksatu(c,lev) = tkmg(c,lev)*0.57D0**watsat(c,lev)
          tkdry(c,lev) = ((0.135D0*bd + 64.7D0) / &
                  (2.7D3 - 0.947D0*bd))*(1.D0-om_frac) + om_tkd*om_frac
          csol(c,lev) = ((1.D0-om_frac)*(2.128D0*sand+2.385D0*clay) / &
                  (sand+clay) + om_csol*om_frac)*1.D6  ! J/(m3 K)
          if ( lev .gt. nlevsoi ) then
            csol(c,lev) = csol_bedrock
          end if
          watdry(c,lev) = watsat(c,lev) * &
                  (316230.D0/sucsat(c,lev)) ** (-1.D0/bsw(c,lev))
          watopt(c,lev) = watsat(c,lev) * &
                  (158490.D0/sucsat(c,lev)) ** (-1.D0/bsw(c,lev))
          !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
          ! water content at field capacity, defined as hk = 0.1 mm/day
          ! used eqn (7.70) in CLM3 technote with
          ! k = 0.1 (mm/day) / (# seconds/day)
          watfc(c,lev) = watsat(c,lev) * &
                  (0.1D0 / (hksat(c,lev)*secspday))** &
                  (1.D0/(2.D0*bsw(c,lev)+3.D0))
        end do
      endif

      ! Define levels
      if ( ltype(l) == istdlak ) then
        if ( lakedepth(c) == spval ) then
          lakedepth(c) = zlak(nlevlak) + 0.5D0*dzlak(nlevlak)
          z_lake(c,1:nlevlak) = zlak(1:nlevlak)
          dz_lake(c,1:nlevlak) = dzlak(1:nlevlak)
        else if ( lakedepth(c) > 1.D0 .and. lakedepth(c) < 5000.D0 ) then
          depthratio = lakedepth(c) / (zlak(nlevlak) + 0.5D0*dzlak(nlevlak))
          z_lake(c,1) = zlak(1)
          dz_lake(c,1) = dzlak(1)
          dz_lake(c,2:nlevlak-1) = dzlak(2:nlevlak-1)*depthratio
          dz_lake(c,nlevlak) = dzlak(nlevlak)*depthratio - &
                              (dz_lake(c,1) - dzlak(1)*depthratio)
          do lev = 2 , nlevlak
            z_lake(c,lev) = z_lake(c,lev-1) + &
                    (dz_lake(c,lev-1)+dz_lake(c,lev))/2.D0
          end do
        else if ( lakedepth(c) > 0.D0 .and. lakedepth(c) <= 1.D0 ) then
          dz_lake(c,:) = lakedepth(c) / nlevlak;
          z_lake(c,1) = dz_lake(c,1) / 2.D0;
          do lev = 2 , nlevlak
            z_lake(c,lev) = z_lake(c,lev-1) + &
                    (dz_lake(c,lev-1)+dz_lake(c,lev))/2.D0
          end do
        else
          write(stderr,*) 'Bad lake depth: lakedepth, c: ', lakedepth(c), c
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
        z(c,1:nlevgrnd) = zsoi(1:nlevgrnd)
        zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
        dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
      end if

      ! Initialize h2osoi_vol. Moved here from BiogeophysRest to avoid
      ! problems with dz not being defined for lakes there, and to avoid
      ! adding dz to the restart file.
      ! Note (Attn EK): if this subroutine is ever called before
      ! SLakeRest, then this code should be moved there.

      if ( ltype(l) == istdlak .and. h2osoi_liq(c,1) /= spval ) then
        ! If not spval then this was available on the restart file
        ! and it was not an old restart file.
        ! Otherwise it will be set in makearbinit.
        do lev = 1 , nlevgrnd
          h2osoi_vol(c,lev) = h2osoi_liq(c,lev)/(dz(c,lev)*denh2o) + &
                              h2osoi_ice(c,lev)/(dz(c,lev)*denice)
        end do
      end if
    end do

    if ( myid == italk ) then
      write (stdout,*) 'Successfully initialized time invariant '// &
                'variables for lakes'
    end if

  end subroutine initTimeConst

end module mod_clm_initslake
