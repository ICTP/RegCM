module mod_clm_type

  use mod_intkinds
  use mod_realkinds
  use mod_clm_domain

  !
  ! Define derived type hierarchy. Includes declaration of
  ! the clm derived type and 1d mapping arrays.
  !
  ! --------------------------------------------------------
  ! gridcell types can have values of
  ! --------------------------------------------------------
  !   1 => default
  ! --------------------------------------------------------
  ! landunits types can have values of (see clm_varcon.F90)
  ! --------------------------------------------------------
  !   1  => (istsoil)    soil (vegetated or bare soil landunit)
  !   2  => (istice)     land ice
  !   3  => (istdlak)    deep lake
  !   4  => (istslak) shall lake (not currently implemented;
  !          SLake implementation has variable depth)
  !   5  => (istwet)     wetland
  !   6  => (isturb)     urban
  !   8  => (istcrop)    crop (only for crop configuration)
  ! --------------------------------------------------------
  ! column types can have values of
  ! --------------------------------------------------------
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istice)           land ice
  !   3  => (istdlak)          deep lake
  !   4  => (istslak)          shallow lake
  !   5  => (istwet)           wetland
  !   61 => (icol_roof)        urban roof
  !   62 => (icol_sunwall)     urban sunwall
  !   63 => (icol_shadewall)   urban shadewall
  !   64 => (icol_road_imperv) urban impervious road
  !   65 => (icol_road_perv)   urban pervious road
  ! --------------------------------------------------------
  ! pft types can have values of
  ! --------------------------------------------------------
  !   0  => not vegetated
  !   1  => needleleaf evergreen temperate tree
  !   2  => needleleaf evergreen boreal tree
  !   3  => needleleaf deciduous boreal tree
  !   4  => broadleaf evergreen tropical tree
  !   5  => broadleaf evergreen temperate tree
  !   6  => broadleaf deciduous tropical tree
  !   7  => broadleaf deciduous temperate tree
  !   8  => broadleaf deciduous boreal tree
  !   9  => broadleaf evergreen shrub
  !   10 => broadleaf deciduous temperate shrub
  !   11 => broadleaf deciduous boreal shrub
  !   12 => c3 arctic grass
  !   13 => c3 non-arctic grass
  !   14 => c4 grass
  !   15 => c3_crop
  !   16 => c3_irrigated
  !   17 => corn
  !   18 => irrigated corn
  !   19 => spring temperate cereal
  !   20 => irrigated spring temperate cereal
  !   21 => winter temperate cereal
  !   22 => irrigated winter temperate cereal
  !   23 => soybean
  !   24 => irrigated soybean
  ! --------------------------------------------------------
  !

  implicit none

  private

  save

  !----------------------------------------------------
  ! Begin definition of conservation check structures
  !----------------------------------------------------
  ! energy balance structure
  !----------------------------------------------------
  type, public :: energy_balance_type
    ! Soil/lake energy conservation error (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errsoi => null()
    ! Surface energy conservation error (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errseb => null()
    ! Solar radiation conservation error (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errsol => null()
    ! Longwave radiation conservation error (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errlon => null()
  end type energy_balance_type

  !----------------------------------------------------
  ! water balance structure
  !----------------------------------------------------
  type, public :: water_balance_type
    !water mass begining of the time step
    real(rk8), pointer, contiguous, dimension(:) :: begwb => null()
    !water mass end of the time step
    real(rk8), pointer, contiguous, dimension(:) :: endwb => null()
    !water conservation error (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: errh2o => null()
  end type water_balance_type

  !----------------------------------------------------
  ! carbon balance structure
  !----------------------------------------------------
  type, public :: carbon_balance_type
    !carbon mass, beginning of time step (gC/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: begcb => null()
    !carbon mass, end of time step (gC/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: endcb => null()
    !carbon balance error for the timestep (gC/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errcb => null()
  end type carbon_balance_type

  !----------------------------------------------------
  ! nitrogen balance structure
  !----------------------------------------------------
  type, public :: nitrogen_balance_type
    !nitrogen mass, beginning of time step (gN/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: begnb => null()
    !nitrogen mass, end of time step (gN/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: endnb => null()
    !nitrogen balance error for the timestep (gN/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: errnb => null()
  end type nitrogen_balance_type

  !----------------------------------------------------
  ! End definition of conservation check structures
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Begin definition of structures defined at the pft_type level
  !----------------------------------------------------

  !----------------------------------------------------
  ! pft physical state variables structure
  !----------------------------------------------------
  type, public :: pft_pstate_type
    !60-day running mean of tot. precipitation (mm/s)
    ! added by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: prec60 => null()
    !10-day running mean of tot. precipitation (mm/s)
    ! added by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: prec10 => null()
    !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4), pointer, contiguous, dimension(:) :: frac_veg_nosno => null()
    !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4), pointer, contiguous, dimension(:) :: frac_veg_nosno_alb => null()
    !vegetation emissivity
    real(rk8), pointer, contiguous, dimension(:) :: emv => null()
    !roughness length over vegetation, momentum [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0mv => null()
    !roughness length over vegetation, sensible heat [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0hv => null()
    !roughness length over vegetation, latent heat [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0qv => null()
    !fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: rootfr => null()
    !effective fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: rootr => null()
    !root resistance by layer (0-1)  (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: rresis => null()
    !Maximum allowed dew [mm]
    real(rk8), pointer, contiguous, dimension(:) :: dewmx => null()
    !sunlit stomatal resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: rssun => null()
    !relative humidity of the canopy air vs leaf
    real(rk8), pointer, contiguous, dimension(:) :: rhal => null()
    !vpd of the canopy air vs leaf
    real(rk8), pointer, contiguous, dimension(:) :: vpdal => null()
    !shaded stomatal resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: rssha => null()
    !canopy layer: sunlit leaf stomatal resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:,:) :: rssha_z => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: rssun_z => null()
    !canopy layer: shaded leaf stomatal resistance (s/m)
    !sunlit projected leaf area index
    real(rk8), pointer, contiguous, dimension(:) :: laisun => null()
    !shaded projected leaf area index
    real(rk8), pointer, contiguous, dimension(:) :: laisha => null()
    !sunlit leaf area for canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: laisun_z => null()
    !shaded leaf area for canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: laisha_z => null()
    !transpiration wetness factor (0 to 1)
    real(rk8), pointer, contiguous, dimension(:) :: btran => null()
    ! root zone soil wetness factor (0 to 1) added by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: btran2 => null()
    !sunlit fraction of canopy
    real(rk8), pointer, contiguous, dimension(:) :: fsun => null()
    !one-sided leaf area index, no burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: tlai => null()
    !one-sided stem area index, no burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: tsai => null()
    !one-sided leaf area index with burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: elai => null()
    !one-sided stem area index with burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: esai => null()
    !fraction of canopy that is wet (0 to 1)
    real(rk8), pointer, contiguous, dimension(:) :: fwet => null()
    !fraction of foliage that is green and dry [-] (new)
    real(rk8), pointer, contiguous, dimension(:) :: fdry => null()
    !change in t_veg, last iteration (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: dt_veg => null()
    !canopy top (m)
    real(rk8), pointer, contiguous, dimension(:) :: htop => null()
    !canopy bottom (m)
    real(rk8), pointer, contiguous, dimension(:) :: hbot => null()
    !momentum roughness length (m)
    real(rk8), pointer, contiguous, dimension(:) :: z0m => null()
    !displacement height (m)
    real(rk8), pointer, contiguous, dimension(:) :: displa => null()
    !surface albedo (direct)                              (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albd => null()
    !surface albedo (diffuse)                             (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albi => null()
    !flux absorbed by canopy per unit direct flux         (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabd => null()
    !flux absorbed by sunlit canopy per unit direct flux  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabd_sun => null()
    !flux absorbed by shaded canopy per unit direct flux  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabd_sha => null()
    !flux absorbed by canopy per unit diffuse flux        (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabi => null()
    !flux absorbed by sunlit canopy per unit diffuse flux (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabi_sun => null()
    !flux absorbed by shaded canopy per unit diffuse flux (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: fabi_sha => null()
    !down direct flux below canopy per unit direct flx    (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: ftdd => null()
    !down diffuse flux below canopy per unit direct flx   (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: ftid => null()
    !down diffuse flux below canopy per unit diffuse flx  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: ftii => null()
    ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(rk8), pointer, contiguous, dimension(:) :: vcmaxcintsun => null()
    ! leaf to canopy scaling coefficient, shaded leaf vcmax
    real(rk8), pointer, contiguous, dimension(:) :: vcmaxcintsha => null()
    !number of canopy layers
    integer(ik4), pointer, contiguous, dimension(:) :: ncan => null()
    !number of canopy layers, above snow for radiative transfer
    integer(ik4), pointer, contiguous, dimension(:) :: nrad => null()
    !absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: fabd_sun_z => null()
    !absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: fabd_sha_z => null()
    !absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: fabi_sun_z => null()
    !absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: fabi_sha_z => null()
    !sunlit fraction of canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: fsun_z => null()
    !tlai increment for canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: tlai_z => null()
    !tsai increment for canopy layer
    real(rk8), pointer, contiguous, dimension(:,:) :: tsai_z => null()
    !10-m wind (m/s) (for dust model)
    real(rk8), pointer, contiguous, dimension(:) :: u10 => null()
    !10-m wind (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: u10_clm => null()
    !atmospheric wind speed plus convective velocity (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: va => null()
    !Bulk Richardson number
    real(rk8), pointer, contiguous, dimension(:) :: br1 => null()
    !aerodynamical resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: ram1 => null()
    !Thermal resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: rah1 => null()
    !aerodynamical resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: ram1_lake => null()
    ! crop burn date
    integer(ik4), pointer, contiguous, dimension(:) :: burndate => null()
    !fractional humidity at leaf surface (dimensionless)
    real(rk8), pointer, contiguous, dimension(:) :: rh_leaf => null()
    !fractional humidity of canopy air (dimensionless)
    real(rk8), pointer, contiguous, dimension(:) :: rhaf => null()
    !friction velocity (m/s) (for dust model)
    real(rk8), pointer, contiguous, dimension(:) :: fv => null()
    !wind forcing height (10m+z0m+d) (m)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_u_pft => null()
    !temperature forcing height (10m+z0m+d) (m)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_t_pft => null()
    !specific humidity forcing height (10m+z0m+d) (m)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_q_pft => null()
    ! decrease of pft weight (0-1) on the column for the timestep
    ! added by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: lfpftd => null()
    ! Variables for prognostic crop model
    ! cold hardening index?
    real(rk8), pointer, contiguous, dimension(:) :: hdidx => null()
    ! cumulative vernalization d?ependence?
    real(rk8), pointer, contiguous, dimension(:) :: cumvd => null()
    ! max hgt attained by a crop during yr (m)
    real(rk8), pointer, contiguous, dimension(:) :: htmx => null()
    ! vernalization factor for cereal
    real(rk8), pointer, contiguous, dimension(:) :: vf => null()
    ! growing degree days (gdd) needed to harvest (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gddmaturity => null()
    ! growing degree-days base  0C from planting  (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd0 => null()
    ! growing degree-days base  8C from planting  (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd8 => null()
    ! growing degree-days base 10C from planting  (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd10 => null()
    ! 20-year average of gdd0                     (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd020 => null()
    ! 20-year average of gdd8                     (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd820 => null()
    ! 20-year average of gdd10                    (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gdd1020 => null()
    ! accum gdd past planting date for crop       (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gddplant => null()
    ! growing degree-days from planting (top two soil layers) (ddays)
    real(rk8), pointer, contiguous, dimension(:) :: gddtsoi => null()
    ! heat unit index needed from planting to leaf emergence
    real(rk8), pointer, contiguous, dimension(:) :: huileaf => null()
    ! heat unit index needed to reach vegetative maturity
    real(rk8), pointer, contiguous, dimension(:) :: huigrain => null()
    ! saved leaf allocation coefficient from phase 2
    real(rk8), pointer, contiguous, dimension(:) :: aleafi => null()
    ! saved stem allocation coefficient from phase 2
    real(rk8), pointer, contiguous, dimension(:) :: astemi => null()
    ! leaf allocation coefficient
    real(rk8), pointer, contiguous, dimension(:) :: aleaf => null()
    ! stem allocation coefficient
    real(rk8), pointer, contiguous, dimension(:) :: astem => null()
    ! Flag, true if planted, not harvested
    logical, pointer, contiguous, dimension(:) :: croplive => null()
    ! Flag, true if planted
    logical, pointer, contiguous, dimension(:) :: cropplant => null()
    ! harvest date
    ! cropplant and harvdate could be 2D to facilitate rotation
    integer(ik4), pointer, contiguous, dimension(:) :: harvdate => null()
    ! date of planting
    integer(ik4), pointer, contiguous, dimension(:) :: idop => null()
    ! 1: max allowed lai; 0: not at max
    integer(ik4), pointer, contiguous, dimension(:) :: peaklai => null()
    !deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
    real(rk8), pointer, contiguous, dimension(:) :: vds => null()
    !sunlit 13c fractionation ([])
    real(rk8), pointer, contiguous, dimension(:) :: alphapsnsun => null()
    !shaded 13c fractionation ([])
    real(rk8), pointer, contiguous, dimension(:) :: alphapsnsha => null()
    ! sand fraction
    real(rk8), pointer, contiguous, dimension(:) :: sandfrac => null()
    ! clay fraction
    real(rk8), pointer, contiguous, dimension(:) :: clayfrac => null()
    ! for dry deposition of chemical tracers
    ! difference between lai month one and month two
    real(rk8), pointer, contiguous, dimension(:) :: mlaidiff => null()
    ! aerodynamical resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: rb1 => null()
    ! 12 months of monthly lai from input data set
    real(rk8), pointer, contiguous, dimension(:,:) :: annlai => null()

    ! New variable for methane code
#ifdef LCH4
    !tracer conductance for boundary layer [m/s]
    real(rk8), pointer, contiguous, dimension(:) :: grnd_ch4_cond
    !tracer conductance for canopy [m/s]
    real(rk8), pointer, contiguous, dimension(:) :: canopy_cond
#endif
    ! and vertical profiles for calculating fluxes
    ! (1/m) profile of leaves
    real(rk8), pointer, contiguous, dimension(:,:) :: leaf_prof => null()
    ! (1/m) profile of fine roots
    real(rk8), pointer, contiguous, dimension(:,:) :: froot_prof => null()
    ! (1/m) profile of coarse roots
    real(rk8), pointer, contiguous, dimension(:,:) :: croot_prof => null()
    ! (1/m) profile of stems
    real(rk8), pointer, contiguous, dimension(:,:) :: stem_prof => null()
  end type pft_pstate_type

  type, public :: pft_psynstate_type
    ! true if C3 and false if C4
    logical, pointer, contiguous, dimension(:) :: c3flag => null()
    ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: ac => null()
    ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: aj => null()
    ! product-limited (C3) or CO2-limited (C4) gross photosynthesis
    ! (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: ap => null()
    ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: ag => null()
    ! net leaf photosynthesis (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: an => null()
    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: vcmax_z => null()
    ! CO2 compensation point (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: cp => null()
    ! Michaelis-Menten constant for CO2 (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: kc => null()
    ! Michaelis-Menten constant for O2 (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: ko => null()
    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(rk8), pointer, contiguous, dimension(:) :: qe => null()
    ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: tpu_z => null()
    ! initial slope of CO2 response curve (C4 plants)
    real(rk8), pointer, contiguous, dimension(:,:) :: kp_z => null()
    ! empirical curvature parameter for ac, aj photosynthesis co-limitation
    real(rk8), pointer, contiguous, dimension(:) :: theta_cj => null()
    ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(rk8), pointer, contiguous, dimension(:) :: bbb => null()
    ! Ball-Berry slope of conductance-photosynthesis relationship
    real(rk8), pointer, contiguous, dimension(:) :: mbb => null()
    ! leaf stomatal conductance (umol H2O/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gs_mol => null()
    ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(rk8), pointer, contiguous, dimension(:) :: gb_mol => null()
  end type pft_psynstate_type

  !----------------------------------------------------
  ! pft ecophysiological constants structure
  !----------------------------------------------------
  type, public :: pft_epc_type
    !value for not vegetated
    integer(ik4), pointer, contiguous, dimension(:) :: noveg => null()
    !tree or not?
    integer(ik4), pointer, contiguous, dimension(:) :: tree => null()
    !soil water potential at full stomatal opening (mm)
    real(rk8), pointer, contiguous, dimension(:) :: smpso => null()
    !soil water potential at full stomatal closure (mm)
    real(rk8), pointer, contiguous, dimension(:) :: smpsc => null()
    !foliage nitrogen limitation factor (-)
    real(rk8), pointer, contiguous, dimension(:) :: fnitr => null()
    !foliage nitrogen (%)
    real(rk8), pointer, contiguous, dimension(:) :: foln => null()
    !characteristic leaf dimension (m)
    real(rk8), pointer, contiguous, dimension(:) :: dleaf => null()
    !photosynthetic pathway: 0. = c4, 1. = c3
    real(rk8), pointer, contiguous, dimension(:) :: c3psn => null()
    !leaf/stem orientation index
    real(rk8), pointer, contiguous, dimension(:) :: xl => null()
    !leaf reflectance: 1=vis, 2=nir   (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: rhol => null()
    !stem reflectance: 1=vis, 2=nir   (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: rhos => null()
    !leaf transmittance: 1=vis, 2=nir (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: taul => null()
    !stem transmittance: 1=vis, 2=nir (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: taus => null()
    !ratio of momentum roughness length to canopy top height (-)
    real(rk8), pointer, contiguous, dimension(:) :: z0mr => null()
    !ratio of displacement height to canopy top height (-)
    real(rk8), pointer, contiguous, dimension(:) :: displar => null()
    !CLM rooting distribution parameter [1/m]
    real(rk8), pointer, contiguous, dimension(:) :: roota_par => null()
    !CLM rooting distribution parameter [1/m]
    real(rk8), pointer, contiguous, dimension(:) :: rootb_par => null()
    ! new variables for CN code
    !wood density (gC/m3)
    real(rk8), pointer, contiguous, dimension(:) :: dwood => null()
    !specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous, dimension(:) :: slatop => null()
    !dSLA/dLAI, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous, dimension(:) :: dsladlai => null()
    !leaf C:N (gC/gN)
    real(rk8), pointer, contiguous, dimension(:) :: leafcn => null()
    !fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
    real(rk8), pointer, contiguous, dimension(:) :: flnr => null()
    !binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8), pointer, contiguous, dimension(:) :: woody => null()
    !leaf litter C:N (gC/gN)
    real(rk8), pointer, contiguous, dimension(:) :: lflitcn => null()
    !fine root C:N (gC/gN)
    real(rk8), pointer, contiguous, dimension(:) :: frootcn => null()
    !live wood (phloem and ray parenchyma) C:N (gC/gN)
    real(rk8), pointer, contiguous, dimension(:) :: livewdcn => null()
    !dead wood (xylem and heartwood) C:N (gC/gN)
    real(rk8), pointer, contiguous, dimension(:) :: deadwdcn => null()
    !grain C:N (gC/gN) for prognostic crop model
    real(rk8), pointer, contiguous, dimension(:) :: graincn => null()
    !allocation parameter: new fine root C per new leaf C (gC/gC)
    real(rk8), pointer, contiguous, dimension(:) :: froot_leaf => null()
    !allocation parameter: new stem c per new leaf C (gC/gC)
    real(rk8), pointer, contiguous, dimension(:) :: stem_leaf => null()
    !allocation parameter: new coarse root C per new stem C (gC/gC)
    real(rk8), pointer, contiguous, dimension(:) :: croot_stem => null()
    ! allocation parameter: fraction of new wood that is live
    ! (phloem and ray parenchyma) (no units)
    real(rk8), pointer, contiguous, dimension(:) :: flivewd => null()
    ! allocation parameter: fraction of allocation that goes to currently
    ! displayed growth, remainder to storage
    real(rk8), pointer, contiguous, dimension(:) :: fcur => null()
    !leaf litter labile fraction
    real(rk8), pointer, contiguous, dimension(:) :: lf_flab => null()
    !leaf litter cellulose fraction
    real(rk8), pointer, contiguous, dimension(:) :: lf_fcel => null()
    !leaf litter lignin fraction
    real(rk8), pointer, contiguous, dimension(:) :: lf_flig => null()
    !fine root litter labile fraction
    real(rk8), pointer, contiguous, dimension(:) :: fr_flab => null()
    !fine root litter cellulose fraction
    real(rk8), pointer, contiguous, dimension(:) :: fr_fcel => null()
    !fine root litter lignin fraction
    real(rk8), pointer, contiguous, dimension(:) :: fr_flig => null()
    !leaf longevity (yrs)
    real(rk8), pointer, contiguous, dimension(:) :: leaf_long => null()
    !binary flag for evergreen leaf habit (0 or 1)
    real(rk8), pointer, contiguous, dimension(:) :: evergreen => null()
    !binary flag for stress-deciduous leaf habit (0 or 1)
    real(rk8), pointer, contiguous, dimension(:) :: stress_decid => null()
    !binary flag for seasonal-deciduous leaf habit (0 or 1)
    real(rk8), pointer, contiguous, dimension(:) :: season_decid => null()

    ! fire variables added by F. Li and S. Levis
    ! combustion completeness factors (0 to 1)
    !combustion completeness factor for leaf
    real(rk8), pointer, contiguous, dimension(:) :: cc_leaf => null()
    !combustion completeness factor for live stem
    real(rk8), pointer, contiguous, dimension(:) :: cc_lstem => null()
    !combustion completeness factor for dead stem
    real(rk8), pointer, contiguous, dimension(:) :: cc_dstem => null()
    !combustion completeness factor for other plant tissues
    real(rk8), pointer, contiguous, dimension(:) :: cc_other => null()
    !  mortality factors (0 to 1)
    !fire-related mortality factor for leaf
    real(rk8), pointer, contiguous, dimension(:) :: fm_leaf => null()
    !fire-related mortality factor for live stem
    real(rk8), pointer, contiguous, dimension(:) :: fm_lstem => null()
    !fire-related mortality factor for dead stem
    real(rk8), pointer, contiguous, dimension(:) :: fm_dstem => null()
    !fire-related mortality factor for other plant tissues
    real(rk8), pointer, contiguous, dimension(:) :: fm_other => null()
    !fire-related mortality factor for fine roots
    real(rk8), pointer, contiguous, dimension(:) :: fm_root => null()
    !fire-related mortality factor for live roots
    real(rk8), pointer, contiguous, dimension(:) :: fm_lroot => null()
    !fire-related mortality factor for dead roots
    real(rk8), pointer, contiguous, dimension(:) :: fm_droot => null()

    !CLM rooting distribution parameter for C and N inputs [unitless]
    real(rk8), pointer, contiguous, dimension(:) :: rootprof_beta => null()
    ! new variables for crop code
    ! fertilizer applied
    real(rk8), pointer, contiguous, dimension(:) :: fertnitro => null()
    ! C:N during grain fill; leaf
    real(rk8), pointer, contiguous, dimension(:) :: fleafcn => null()
    ! C:N during grain fill; froot
    real(rk8), pointer, contiguous, dimension(:) :: ffrootcn => null()
    ! C:N during grain fill; stem
    real(rk8), pointer, contiguous, dimension(:) :: fstemcn => null()
  end type pft_epc_type

  type, public :: decomp_cascade_type
    !-- properties of each pathway along decomposition cascade
    ! name of transition
    character(len=8), pointer, contiguous, dimension(:) :: cascade_step_name => null()
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer, contiguous, dimension(:) :: cascade_donor_pool => null()
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer, contiguous, dimension(:) :: cascade_receiver_pool => null()
    !-- properties of each decomposing pool
    ! TRUE => pool has fixed C:N ratio
    logical,  pointer, contiguous, dimension(:) :: floating_cn_ratio_decomp_pools => null()
    ! name of pool for restart files
    character(len=8), pointer, contiguous, dimension(:) :: decomp_pool_name_restart => null()
    ! name of pool for history files
    character(len=8), pointer, contiguous, dimension(:) :: decomp_pool_name_history => null()
    ! name of pool for netcdf long names
    character(len=20), pointer, contiguous, dimension(:) :: decomp_pool_name_long => null()
    ! name of pool for netcdf short names
    character(len=8), pointer, contiguous, dimension(:) :: decomp_pool_name_short => null()
    ! TRUE => pool is a litter pool
    logical, pointer, contiguous, dimension(:) :: is_litter => null()
    ! TRUE => pool is a soil pool
    logical, pointer, contiguous, dimension(:) :: is_soil => null()
    ! TRUE => pool is a cwd pool
    logical, pointer, contiguous, dimension(:) :: is_cwd => null()
    ! c:n ratio for initialization of pools
    real(rk8), pointer, contiguous, dimension(:) :: initial_cn_ratio => null()
    ! initial concentration for seeding at spinup
    real(rk8), pointer, contiguous, dimension(:) :: initial_stock => null()
    ! TRUE => pool is metabolic material
    logical, pointer, contiguous, dimension(:) :: is_metabolic => null()
    ! TRUE => pool is cellulose
    logical, pointer, contiguous, dimension(:) :: is_cellulose => null()
    ! TRUE => pool is lignin
    logical, pointer, contiguous, dimension(:) :: is_lignin => null()
    ! factor by which to scale AD and relevant processes by
    real(rk8), pointer, contiguous, dimension(:) :: spinup_factor => null()
  end type decomp_cascade_type

#if (defined CNDV)
  !----------------------------------------------------
  ! pft DGVM-specific ecophysiological constants structure
  !----------------------------------------------------
  type, public :: pft_dgvepc_type
    !tree maximum crown area [m2]
    real(rk8), pointer, contiguous, dimension(:) :: crownarea_max => null()
    !minimum coldest monthly mean temperature [units?]
    real(rk8), pointer, contiguous, dimension(:) :: tcmin => null()
    !maximum coldest monthly mean temperature [units?]
    real(rk8), pointer, contiguous, dimension(:) :: tcmax => null()
    !minimum growing degree days (at or above 5 C)
    real(rk8), pointer, contiguous, dimension(:) :: gddmin => null()
    !upper limit of temperature of the warmest month [units?]
    real(rk8), pointer, contiguous, dimension(:) :: twmax => null()
    !parameter in allometric equation
    real(rk8), pointer, contiguous, dimension(:) :: reinickerp => null()
    !parameter in allometric
    real(rk8), pointer, contiguous, dimension(:) :: allom1 => null()
    !parameter in allometric
    real(rk8), pointer, contiguous, dimension(:) :: allom2 => null()
    !parameter in allometric
    real(rk8), pointer, contiguous, dimension(:) :: allom3 => null()
  end type pft_dgvepc_type
#endif

  !----------------------------------------------------
  ! pft ecophysiological variables structure
  !----------------------------------------------------
  type, public :: pft_epv_type
    !dormancy flag
    real(rk8), pointer, contiguous, dimension(:) :: dormant_flag => null()
    !number of days since last dormancy
    real(rk8), pointer, contiguous, dimension(:) :: days_active => null()
    !onset flag
    real(rk8), pointer, contiguous, dimension(:) :: onset_flag => null()
    !onset days counter
    real(rk8), pointer, contiguous, dimension(:) :: onset_counter => null()
    !onset flag for growing degree day sum
    real(rk8), pointer, contiguous, dimension(:) :: onset_gddflag => null()
    !onset freezing degree days counter
    real(rk8), pointer, contiguous, dimension(:) :: onset_fdd => null()
    !onset growing degree days
    real(rk8), pointer, contiguous, dimension(:) :: onset_gdd => null()
    !onset soil water index
    real(rk8), pointer, contiguous, dimension(:) :: onset_swi => null()
    !offset flag
    real(rk8), pointer, contiguous, dimension(:) :: offset_flag => null()
    !offset days counter
    real(rk8), pointer, contiguous, dimension(:) :: offset_counter => null()
    !offset freezing degree days counter
    real(rk8), pointer, contiguous, dimension(:) :: offset_fdd => null()
    !offset soil water index
    real(rk8), pointer, contiguous, dimension(:) :: offset_swi => null()
    !>0 fertilize; <=0 not
    real(rk8), pointer, contiguous, dimension(:) :: fert_counter => null()
    !1: grain fill stage; 0: not
    real(rk8), pointer, contiguous, dimension(:) :: grain_flag => null()
    !long growing season factor [0-1]
    real(rk8), pointer, contiguous, dimension(:) :: lgsf => null()
    !background litterfall rate (1/s)
    real(rk8), pointer, contiguous, dimension(:) :: bglfr => null()
    !background transfer growth rate (1/s)
    real(rk8), pointer, contiguous, dimension(:) :: bgtr => null()
    !daylength (seconds)
    real(rk8), pointer, contiguous, dimension(:) :: dayl => null()
    !daylength from previous timestep (seconds)
    real(rk8), pointer, contiguous, dimension(:) :: prev_dayl => null()
    !annual average 2m air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: annavg_t2m => null()
    !temporary average 2m air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: tempavg_t2m => null()
    !GPP flux before downregulation (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: gpp => null()
    !C flux available for allocation (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: availc => null()
    !C flux assigned to recovery of negative cpool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: xsmrpool_recover => null()
    !C13/C(12+13) ratio for xsmrpool (proportion)
    real(rk8), pointer, contiguous, dimension(:) :: xsmrpool_c13ratio => null()
    !fraction of current allocation to display as new growth (DIM)
    real(rk8), pointer, contiguous, dimension(:) :: alloc_pnow => null()
    !C allocation index (DIM)
    real(rk8), pointer, contiguous, dimension(:) :: c_allometry => null()
    !N allocation index (DIM)
    real(rk8), pointer, contiguous, dimension(:) :: n_allometry => null()
    !N flux required to support initial GPP (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: plant_ndemand => null()
    !temporary annual sum of potential GPP
    real(rk8), pointer, contiguous, dimension(:) :: tempsum_potential_gpp => null()
    !annual sum of potential GPP
    real(rk8), pointer, contiguous, dimension(:) :: annsum_potential_gpp => null()
    !temporary annual max of retranslocated N pool (gN/m2)
    real(rk8), pointer, contiguous, dimension(:) :: tempmax_retransn => null()
    !annual max of retranslocated N pool (gN/m2)
    real(rk8), pointer, contiguous, dimension(:) :: annmax_retransn => null()
    !N flux available from retranslocation pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: avail_retransn => null()
    !total allocated N flux (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: plant_nalloc => null()
    !total allocated C flux (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: plant_calloc => null()
    !C flux not allocated due to downregulation (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: excess_cflux => null()
    !fractional reduction in GPP due to N limitation (DIM)
    real(rk8), pointer, contiguous, dimension(:) :: downreg => null()
    !previous timestep leaf C litterfall flux (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: prev_leafc_to_litter => null()
    !previous timestep froot C litterfall flux (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: prev_frootc_to_litter => null()
    !temporary annual sum of NPP (gC/m2/yr)
    real(rk8), pointer, contiguous, dimension(:) :: tempsum_npp => null()
    !annual sum of NPP (gC/m2/yr)
    real(rk8), pointer, contiguous, dimension(:) :: annsum_npp => null()
#if (defined CNDV)
    !temporary annual sum of litfall (gC/m2/yr)
    real(rk8), pointer, contiguous, dimension(:) :: tempsum_litfall
    !annual sum of litfall (gC/m2/yr)
    real(rk8), pointer, contiguous, dimension(:) :: annsum_litfall
#endif
    !C13O2/C12O2 in canopy air
    real(rk8), pointer, contiguous, dimension(:) :: rc13_canair => null()
    !C13O2/C12O2 in sunlit canopy psn flux
    real(rk8), pointer, contiguous, dimension(:) :: rc13_psnsun => null()
    !C13O2/C12O2 in shaded canopy psn flux
    real(rk8), pointer, contiguous, dimension(:) :: rc13_psnsha => null()
    !C14O2/C12O2 in atmosphere
    real(rk8), pointer, contiguous, dimension(:) :: rc14_atm => null()
  end type pft_epv_type

  !----------------------------------------------------
  ! pft energy state variables structure
  !----------------------------------------------------
  type, public :: pft_estate_type
    !2 m height surface air temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m => null()
    ! samy : gridded Q10 maintenance respiration at PFT level
    real(rk8), pointer, contiguous, dimension(:) :: q10m => null()
    !daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min => null()
    !daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max => null()
    !instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min_inst => null()
    !instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max_inst => null()
    !2 m height surface specific humidity (kg/kg)
    real(rk8), pointer, contiguous, dimension(:) :: q_ref2m => null()
    !Urban 2 m height surface air temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_u => null()
    !Rural 2 m height surface air temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_r => null()
    !Urban daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min_u => null()
    !Rural daily minimum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min_r => null()
    !Urban daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max_u => null()
    !Rural daily maximum of average 2 m height surface air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max_r => null()
    !Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min_inst_u => null()
    !Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_min_inst_r => null()
    !Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max_inst_u => null()
    !Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m_max_inst_r => null()
    ! 10-day running mean of min 2-m temperature
    real(rk8), pointer, contiguous, dimension(:) :: a10tmin => null()
    ! 5-day running mean of min 2-m temperature
    real(rk8), pointer, contiguous, dimension(:) :: a5tmin => null()
    !10-day running mean of the 2 m temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t10 => null()
    !2 m height surface relative humidity (%)
    real(rk8), pointer, contiguous, dimension(:) :: rh_ref2m => null()
    !Urban 2 m height surface relative humidity (%)
    real(rk8), pointer, contiguous, dimension(:) :: rh_ref2m_u => null()
    !Rural 2 m height surface relative humidity (%)
    real(rk8), pointer, contiguous, dimension(:) :: rh_ref2m_r => null()
    !vegetation temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_veg => null()
    !intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(rk8), pointer, contiguous, dimension(:) :: thm => null()
  end type pft_estate_type

  !----------------------------------------------------
  ! pft water state variables structure
  !----------------------------------------------------
  type, public :: pft_wstate_type
    !canopy water (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: h2ocan => null()
  end type pft_wstate_type

  !----------------------------------------------------
  ! pft carbon state variables structure
  !----------------------------------------------------
  type, public :: pft_cstate_type
    ! (gC/m2) ann max leaf C
    real(rk8), pointer, contiguous, dimension(:) :: leafcmax => null()
    ! variables for prognostic crop model
    ! (gC/m2) grain C
    real(rk8), pointer, contiguous, dimension(:) :: grainc => null()
    ! (gC/m2) grain C storage
    real(rk8), pointer, contiguous, dimension(:) :: grainc_storage => null()
    ! (gC/m2) grain C transfer
    real(rk8), pointer, contiguous, dimension(:) :: grainc_xfer => null()
    !
    ! (gC/m2) leaf C
    real(rk8), pointer, contiguous, dimension(:) :: leafc => null()
    ! (gC/m2) leaf C storage
    real(rk8), pointer, contiguous, dimension(:) :: leafc_storage => null()
    ! (gC/m2) leaf C transfer
    real(rk8), pointer, contiguous, dimension(:) :: leafc_xfer => null()
    ! (gC/m2) fine root C
    real(rk8), pointer, contiguous, dimension(:) :: frootc => null()
    ! (gC/m2) fine root C storage
    real(rk8), pointer, contiguous, dimension(:) :: frootc_storage => null()
    ! (gC/m2) fine root C transfer
    real(rk8), pointer, contiguous, dimension(:) :: frootc_xfer => null()
    ! (gC/m2) live stem C
    real(rk8), pointer, contiguous, dimension(:) :: livestemc => null()
    ! (gC/m2) live stem C storage
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_storage => null()
    ! (gC/m2) live stem C transfer
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_xfer => null()
    ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous, dimension(:) :: deadstemc => null()
    ! (gC/m2) dead stem C storage
    real(rk8), pointer, contiguous, dimension(:) :: deadstemc_storage => null()
    ! (gC/m2) dead stem C transfer
    real(rk8), pointer, contiguous, dimension(:) :: deadstemc_xfer => null()
    ! (gC/m2) live coarse root C
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc => null()
    ! (gC/m2) live coarse root C storage
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc_storage => null()
    ! (gC/m2) live coarse root C transfer
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc_xfer => null()
    ! (gC/m2) dead coarse root C
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootc => null()
    ! (gC/m2) dead coarse root C storage
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootc_storage => null()
    ! (gC/m2) dead coarse root C transfer
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootc_xfer => null()
    ! (gC/m2) growth respiration storage
    real(rk8), pointer, contiguous, dimension(:) :: gresp_storage => null()
    ! (gC/m2) growth respiration transfer
    real(rk8), pointer, contiguous, dimension(:) :: gresp_xfer => null()
    ! (gC/m2) temporary photosynthate C pool
    real(rk8), pointer, contiguous, dimension(:) :: cpool => null()
    ! (gC/m2) abstract C pool to meet excess MR demand
    real(rk8), pointer, contiguous, dimension(:) :: xsmrpool => null()
    ! (gC/m2) pft-level sink for C truncation
    real(rk8), pointer, contiguous, dimension(:) :: pft_ctrunc => null()
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gC/m2) displayed veg carbon, excluding storage and cpool
    real(rk8), pointer, contiguous, dimension(:) :: dispvegc => null()
    ! (gC/m2) stored vegetation carbon, excluding cpool
    real(rk8), pointer, contiguous, dimension(:) :: storvegc => null()
    ! (gC/m2) total vegetation carbon, excluding cpool
    real(rk8), pointer, contiguous, dimension(:) :: totvegc => null()
    ! (gC/m2) total pft-level carbon, including cpool
    real(rk8), pointer, contiguous, dimension(:) :: totpftc => null()
#if (defined CN)
    ! (gC/m2) wood C
    real(rk8), pointer, contiguous, dimension(:) :: woodc
#endif
  end type pft_cstate_type

  !----------------------------------------------------
  ! pft nitrogen state variables structure
  !----------------------------------------------------
  type, public :: pft_nstate_type
    ! variables for prognostic crop model
    ! (gN/m2) grain N
    real(rk8), pointer, contiguous, dimension(:) :: grainn => null()
    ! (gN/m2) grain N storage
    real(rk8), pointer, contiguous, dimension(:) :: grainn_storage => null()
    ! (gN/m2) grain N transfer
    real(rk8), pointer, contiguous, dimension(:) :: grainn_xfer => null()
    !
    ! (gN/m2) leaf N
    real(rk8), pointer, contiguous, dimension(:) :: leafn => null()
    ! (gN/m2) leaf N storage
    real(rk8), pointer, contiguous, dimension(:) :: leafn_storage => null()
    ! (gN/m2) leaf N transfer
    real(rk8), pointer, contiguous, dimension(:) :: leafn_xfer => null()
    ! (gN/m2) fine root N
    real(rk8), pointer, contiguous, dimension(:) :: frootn => null()
    ! (gN/m2) fine root N storage
    real(rk8), pointer, contiguous, dimension(:) :: frootn_storage => null()
    ! (gN/m2) fine root N transfer
    real(rk8), pointer, contiguous, dimension(:) :: frootn_xfer => null()
    ! (gN/m2) live stem N
    real(rk8), pointer, contiguous, dimension(:) :: livestemn => null()
    ! (gN/m2) live stem N storage
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_storage => null()
    ! (gN/m2) live stem N transfer
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_xfer => null()
    ! (gN/m2) dead stem N
    real(rk8), pointer, contiguous, dimension(:) :: deadstemn => null()
    ! (gN/m2) dead stem N storage
    real(rk8), pointer, contiguous, dimension(:) :: deadstemn_storage => null()
    ! (gN/m2) dead stem N transfer
    real(rk8), pointer, contiguous, dimension(:) :: deadstemn_xfer => null()
    ! (gN/m2) live coarse root N
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn => null()
    ! (gN/m2) live coarse root N storage
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_storage => null()
    ! (gN/m2) live coarse root N transfer
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_xfer => null()
    ! (gN/m2) dead coarse root N
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootn => null()
    ! (gN/m2) dead coarse root N storage
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootn_storage => null()
    ! (gN/m2) dead coarse root N transfer
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootn_xfer => null()
    ! (gN/m2) plant pool of retranslocated N
    real(rk8), pointer, contiguous, dimension(:) :: retransn => null()
    ! (gN/m2) temporary plant N pool
    real(rk8), pointer, contiguous, dimension(:) :: npool => null()
    ! (gN/m2) pft-level sink for N truncation
    real(rk8), pointer, contiguous, dimension(:) :: pft_ntrunc => null()
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gN/m2) displayed veg nitrogen, excluding storage
    real(rk8), pointer, contiguous, dimension(:) :: dispvegn => null()
    ! (gN/m2) stored vegetation nitrogen
    real(rk8), pointer, contiguous, dimension(:) :: storvegn => null()
    ! (gN/m2) total vegetation nitrogen
    real(rk8), pointer, contiguous, dimension(:) :: totvegn => null()
    ! (gN/m2) total pft-level nitrogen
    real(rk8), pointer, contiguous, dimension(:) :: totpftn => null()
  end type pft_nstate_type

  !----------------------------------------------------
  ! pft VOC state variables structure
  !----------------------------------------------------
  type, public :: pft_vstate_type
   ! 24hr average vegetation temperature (K)
   real(rk8), pointer, contiguous, dimension(:) :: t_veg24 => null()
   ! 240hr average vegetation temperature (Kelvin)
   real(rk8), pointer, contiguous, dimension(:) :: t_veg240 => null()
   ! 24hr average of direct beam radiation
   real(rk8), pointer, contiguous, dimension(:) :: fsd24 => null()
   ! 240hr average of direct beam radiation
   real(rk8), pointer, contiguous, dimension(:) :: fsd240 => null()
   ! 24hr average of diffuse beam radiation
   real(rk8), pointer, contiguous, dimension(:) :: fsi24 => null()
   ! 240hr average of diffuse beam radiation
   real(rk8), pointer, contiguous, dimension(:) :: fsi240 => null()
   ! 24hr average of sunlit fraction of canopy
   real(rk8), pointer, contiguous, dimension(:) :: fsun24 => null()
   ! 240hr average of sunlit fraction of canopy
   real(rk8), pointer, contiguous, dimension(:) :: fsun240 => null()
   ! leaf area index average over timestep
   real(rk8), pointer, contiguous, dimension(:) :: elai_p => null()
  end type pft_vstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! pft DGVM state variables structure
  !----------------------------------------------------
  type, public :: pft_dgvstate_type
    !accumulated growing degree days above twmax
    real(rk8), pointer, contiguous, dimension(:) :: agddtw => null()
    !accumulated growing degree days above 5
    real(rk8), pointer, contiguous, dimension(:) :: agdd => null()
    !30-day average temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_mo => null()
    !annual min of t_mo (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_mo_min => null()
    !365-day running mean of tot. precipitation
    real(rk8), pointer, contiguous, dimension(:) :: prec365 => null()
    !whether PFT present in patch
    logical, pointer, contiguous, dimension(:) :: present => null()
    !if .false. then exclude seasonal decid pfts from tropics
    logical, pointer, contiguous, dimension(:) :: pftmayexist => null()
    !number of individuals (#/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: nind => null()
    !individual leaf mass
    real(rk8), pointer, contiguous, dimension(:) :: lm_ind => null()
    !LAI per individual
    real(rk8), pointer, contiguous, dimension(:) :: lai_ind => null()
    !foliar projective cover increment (fraction)
    real(rk8), pointer, contiguous, dimension(:) :: fpcinc => null()
    !foliar projective cover on gridcell (fraction)
    real(rk8), pointer, contiguous, dimension(:) :: fpcgrid => null()
    !last yr's fpcgrid
    real(rk8), pointer, contiguous, dimension(:) :: fpcgridold => null()
    !area that each individual tree takes up (m^2)
    real(rk8), pointer, contiguous, dimension(:) :: crownarea => null()
    real(rk8), pointer, contiguous, dimension(:) :: greffic => null()
    real(rk8), pointer, contiguous, dimension(:) :: heatstress => null()
  end type pft_dgvstate_type
#endif

  !----------------------------------------------------
  ! pft energy flux variables structure
  !----------------------------------------------------
  type, public :: pft_eflux_type
    !solar radiation absorbed by soil (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: sabg_soil => null()
    !solar radiation absorbed by snow (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: sabg_snow => null()
    !fsno weighted sum (needed by balancecheck, because fsno changes midway)
    real(rk8), pointer, contiguous, dimension(:) :: sabg_chk => null()
    !solar radiation absorbed by ground (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: sabg => null()
    !solar radiation absorbed by vegetation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: sabv => null()
    !solar radiation absorbed (total) (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsa => null()
    !urban solar radiation absorbed (total) (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsa_u => null()
    !rural solar radiation absorbed (total) (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsa_r => null()
    !solar radiation reflected (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr => null()
    !absorbed PAR for sunlit leaves in canopy layer (W/m**2)
    real(rk8), pointer, contiguous, dimension(:,:) :: parsun_z => null()
    !absorbed PAR for shaded leaves in canopy layer (W/m**2)
    real(rk8), pointer, contiguous, dimension(:,:) :: parsha_z => null()
    !downward longwave radiation below the canopy [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: dlrad => null()
    !upward longwave radiation above the canopy [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: ulrad => null()
    !total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_tot => null()
    !urban total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_tot_u => null()
    !rural total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_tot_r => null()
    !ground evaporation heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_grnd => null()
    !soil heat flux (W/m**2) [+ = into soil]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_soil_grnd => null()
    !urban soil heat flux (W/m**2) [+ = into soil]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_soil_grnd_u => null()
    !rural soil heat flux (W/m**2) [+ = into soil]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_soil_grnd_r => null()
    !total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_tot => null()
    !urban total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_tot_u => null()
    !rural total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_tot_r => null()
    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_grnd => null()
    !sensible heat flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_snow => null()
    !sensible heat flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_soil => null()
    !sensible heat flux from surface water (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_h2osfc => null()
    !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_veg => null()
    !veg evaporation heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_vege => null()
    !veg transpiration heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_vegt => null()
    !sensible heat flux from domestic heating/cooling sources of waste
    ! heat (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_wasteheat_pft => null()
    !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_heat_from_ac_pft => null()
    !traffic sensible heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_traffic_pft => null()
    !total anthropogenic heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_anthro => null()
    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rk8), pointer, contiguous, dimension(:) :: cgrnd => null()
    !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
    real(rk8), pointer, contiguous, dimension(:) :: cgrndl => null()
    !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
    real(rk8), pointer, contiguous, dimension(:) :: cgrnds => null()
    !net heat flux into ground (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_gnet => null()
    ! New lake field
    !net heat flux into lake / snow surface, excluding light
    ! transmission (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_grnd_lake => null()
    !derivative of net ground heat flux wrt soil temp (W/m**2 K)
    real(rk8), pointer, contiguous, dimension(:) :: dgnetdT => null()
    !emitted infrared (longwave) radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lwrad_out => null()
    !net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lwrad_net => null()
    !urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lwrad_net_u => null()
    !rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lwrad_net_r => null()
    !net radiation (W/m**2) [+ = to sfc]
    real(rk8), pointer, contiguous, dimension(:) :: netrad => null()
    !incident direct beam vis solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_vis_d => null()
    !incident direct beam nir solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_nir_d => null()
    !incident diffuse vis solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_vis_i => null()
    !incident diffuse nir solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_nir_i => null()
    !reflected direct beam vis solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_vis_d => null()
    !reflected direct beam nir solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_nir_d => null()
    !reflected diffuse vis solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_vis_i => null()
    !reflected diffuse nir solar radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_nir_i => null()
    !incident direct beam vis solar radiation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_vis_d_ln => null()
    !incident diffuse beam vis solar radiation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_vis_i_ln => null()
    !absorbed par by vegetation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: parveg_ln => null()
    !incident direct beam nir solar radiation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsds_nir_d_ln => null()
    !reflected direct beam vis solar radiation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_vis_d_ln => null()
    !reflected direct beam nir solar radiation at local noon (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_nir_d_ln => null()
    ! absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
    real(rk8), pointer, contiguous, dimension(:,:) :: sabg_lyr => null()
    ! (rural) shortwave radiation penetrating top soisno layer [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sabg_pen => null()
    ! surface forcing of snow with all aerosols (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_aer => null()
    ! surface forcing of snow with BC (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_bc => null()
    ! surface forcing of snow with OC (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_oc => null()
    ! surface forcing of snow with dust (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_dst => null()
    ! surface forcing of snow with all aerosols, averaged only when
    ! snow is present (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_aer_sno => null()
    ! surface forcing of snow with BC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_bc_sno => null()
    ! surface forcing of snow with OC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_oc_sno => null()
    ! surface forcing of snow with dust, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: sfc_frc_dst_sno => null()
    ! reflected direct beam vis solar radiation from snow (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_sno_vd => null()
    ! reflected direct beam NIR solar radiation from snow (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_sno_nd => null()
    ! reflected diffuse vis solar radiation from snow (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_sno_vi => null()
    ! reflected diffuse NIR solar radiation from snow (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_sno_ni => null()
    ! incident visible, direct radiation on snow (for history files)  [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: fsds_sno_vd => null()
    ! incident near-IR, direct radiation on snow (for history files)  [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: fsds_sno_nd => null()
    ! incident visible, diffuse radiation on snow (for history files) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: fsds_sno_vi => null()
    ! incident near-IR, diffuse radiation on snow (for history files) [W/m2]
    real(rk8), pointer, contiguous, dimension(:) :: fsds_sno_ni => null()
  end type pft_eflux_type

  !----------------------------------------------------
  ! pft momentum flux variables structure
  !----------------------------------------------------
  type, public :: pft_mflux_type
    !wind (shear) stress: e-w (kg/m/s**2)
    real(rk8), pointer, contiguous, dimension(:) ::  taux => null()
    !wind (shear) stress: n-s (kg/m/s**2)
    real(rk8), pointer, contiguous, dimension(:) ::  tauy => null()
  end type pft_mflux_type

  !----------------------------------------------------
  ! pft water flux variables structure
  !----------------------------------------------------
  type, public :: pft_wflux_type
    !interception of precipitation [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_prec_intr => null()
    !water onto ground including canopy runoff [kg/(m2 s)]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_prec_grnd => null()
    !rain on ground after interception (mm H2O/s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_rain_grnd => null()
    !snow on ground after interception (mm H2O/s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snow_grnd => null()
    !excess snowfall due to snow capping (mm H2O /s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snwcp_ice => null()
    !excess rainfall due to snow capping (mm H2O /s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snwcp_liq => null()
    !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_veg => null()
    !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_tran_veg => null()
    !evaporation from leaves and stems
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_can => null()
    !soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_soi => null()
    !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_tot => null()
    !ground surface evaporation rate (mm H2O/s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_grnd => null()
    !ground surface dew formation (mm H2O /s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_dew_grnd => null()
    !sublimation rate from snow pack (mm H2O /s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_sub_snow => null()
    !surface dew added to snow pack (mm H2O /s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_dew_snow => null()
    !snow evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_ev_snow => null()
    !soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_ev_soil => null()
    !h2osfc evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_ev_h2osfc => null()
  end type pft_wflux_type

  !----------------------------------------------------
  ! pft carbon flux variables structure
  !----------------------------------------------------
  type, public :: pft_cflux_type
    !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsun => null()
    !shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsha => null()
    !canopy layer: sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:,:) :: psnsun_z => null()
    !canopy layer: shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:,:) :: psnsha_z => null()
    !intracellular sunlit leaf CO2 (Pa)
    real(rk8), pointer, contiguous, dimension(:,:) :: cisun_z => null()
    !intracellular shaded leaf CO2 (Pa)
    real(rk8), pointer, contiguous, dimension(:,:) :: cisha_z => null()
    !sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:) :: lmrsun => null()
    !shaded leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:) :: lmrsha => null()
    !canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: lmrsun_z => null()
    !canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: lmrsha_z => null()
    !photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer, contiguous, dimension(:) :: fpsn => null()
    !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
    real(rk8), pointer, contiguous, dimension(:) :: fco2 => null()
    !Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsun_wc => null()
    !Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsha_wc => null()
    !Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer, contiguous, dimension(:) :: fpsn_wc => null()
    !RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsun_wj => null()
    !RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsha_wj => null()
    !RuBP-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer, contiguous, dimension(:) :: fpsn_wj => null()
    !product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsun_wp => null()
    !product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsha_wp => null()
    !product-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer, contiguous, dimension(:) :: fpsn_wp => null()
    ! new variables for CN code
    ! gap mortality fluxes
    ! leaf C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_to_litter => null()
    ! leaf C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_storage_to_litter => null()
    ! leaf C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_xfer_to_litter => null()
    ! fine root C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_to_litter => null()
    ! fine root C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_storage_to_litter => null()
    ! fine root C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_xfer_to_litter => null()
    ! live stem C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_to_litter => null()
    ! live stem C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_storage_to_litter => null()
    ! live stem C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_xfer_to_litter => null()
    ! dead stem C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_to_litter => null()
    ! dead stem C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_storage_to_litter => null()
    ! dead stem C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_xfer_to_litter => null()
    ! live coarse root C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_to_litter => null()
    ! live coarse root C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_storage_to_litter => null()
    ! live coarse root C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_xfer_to_litter => null()
    ! dead coarse root C mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_to_litter => null()
    ! dead coarse root C storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_storage_to_litter => null()
    ! dead coarse root C transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_xfer_to_litter => null()
    ! growth respiration storage mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_storage_to_litter => null()
    ! growth respiration transfer mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_xfer_to_litter => null()
    ! harvest mortality fluxes
    ! leaf C harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafc_to_litter => null()
    ! leaf C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafc_storage_to_litter => null()
    ! leaf C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafc_xfer_to_litter => null()
    ! fine root C harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootc_to_litter => null()
    ! fine root C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootc_storage_to_litter => null()
    ! fine root C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootc_xfer_to_litter => null()
    ! live stem C harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemc_to_litter => null()
    ! live stem C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemc_storage_to_litter => null()
    ! live stem C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemc_xfer_to_litter => null()
    ! dead stem C harvest to 10-year product pool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_to_prod10c => null()
    ! dead stem C harvest to 100-year product pool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_to_prod100c => null()
    ! dead stem C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_storage_to_litter => null()
    ! dead stem C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_xfer_to_litter => null()
    ! live coarse root C harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootc_to_litter => null()
    ! live coarse root C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootc_storage_to_litter => null()
    ! live coarse root C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootc_xfer_to_litter => null()
    ! dead coarse root C harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootc_to_litter => null()
    ! dead coarse root C storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootc_storage_to_litter => null()
    ! dead coarse root C transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootc_xfer_to_litter => null()
    ! growth respiration storage harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_gresp_storage_to_litter => null()
    ! growth respiration transfer harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_gresp_xfer_to_litter => null()
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_xsmrpool_to_atm => null()
    ! PFT-level fire C fluxes added by F. Li and S. Levis
    ! (gC/m2/s) fire C emissions from leafc
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_to_fire => null()
    ! (gC/m2/s) fire C emissions from leafc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from leafc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from livestemc
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_to_fire => null()
    ! (gC/m2/s) fire C emissions from livestemc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from livestemc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadstemc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from frootc
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_to_fire => null()
    ! (gC/m2/s) fire C emissions from frootc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from frootc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from livecrootc
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_to_fire => null()
    ! (gC/m2/s) fire C emissions from livecrootc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from livecrootc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadcrootc
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadcrootc_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from deadcrootc_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_xfer_to_fire => null()
    ! (gC/m2/s) fire C emissions from gresp_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_storage_to_fire => null()
    ! (gC/m2/s) fire C emissions from gresp_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_xfer_to_fire => null()
    ! (gC/m2/s) from leafc to litter c due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_to_litter_fire => null()
    ! (gC/m2/s) from leafc_storage to litter C  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_storage_to_litter_fire => null()
    ! (gC/m2/s) from leafc_xfer to litter C  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafc_xfer_to_litter_fire => null()
    ! (gC/m2/s) from livestemc to litter C  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_to_litter_fire => null()
    ! (gC/m2/s) from livestemc_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_storage_to_litter_fire => null()
    !(gC/m2/s) from livestemc_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_xfer_to_litter_fire => null()
    !(gC/m2/s) from livestemc to deadstemc due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemc_to_deadstemc_fire => null()
    !(gC/m2/s) from deadstemc to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_to_litter_fire => null()
    !(gC/m2/s) from deadstemc_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_storage_to_litter_fire => null()
    !(gC/m2/s) from deadstemc_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemc_xfer_to_litter_fire => null()
    !(gC/m2/s) from frootc to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_to_litter_fire => null()
    !(gC/m2/s) from frootc_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_storage_to_litter_fire => null()
    !(gC/m2/s) from frootc_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootc_xfer_to_litter_fire => null()
    !(gC/m2/s) from livecrootc to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_to_litter_fire => null()
    !(gC/m2/s) from livecrootc_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_storage_to_litter_fire => null()
    !(gC/m2/s) from livecrootc_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_xfer_to_litter_fire => null()
    !(gC/m2/s) from livecrootc to deadstemc due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootc_to_deadcrootc_fire => null()
    !(gC/m2/s) from deadcrootc to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_to_litter_fire => null()
    !(gC/m2/s) from deadcrootc_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_storage_to_litter_fire => null()
    !(gC/m2/s) from deadcrootc_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootc_xfer_to_litter_fire => null()
    !(gC/m2/s) from gresp_storage to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_storage_to_litter_fire => null()
    !(gC/m2/s) from gresp_xfer to litter C due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_gresp_xfer_to_litter_fire => null()
    ! phenology fluxes from transfer pools
    ! grain C growth from storage for prognostic crop(gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainc_xfer_to_grainc => null()
    ! leaf C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafc_xfer_to_leafc => null()
    ! fine root C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootc_xfer_to_frootc => null()
    ! live stem C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_xfer_to_livestemc => null()
    ! dead stem C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadstemc_xfer_to_deadstemc => null()
    ! live coarse root C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc_xfer_to_livecrootc => null()
    ! dead coarse root C growth from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootc_xfer_to_deadcrootc => null()
    ! leaf and fine root litterfall
    ! leaf C litterfall (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafc_to_litter => null()
    ! fine root C litterfall (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootc_to_litter => null()
    ! live stem C litterfall (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_to_litter => null()
    ! grain C to food for prognostic crop(gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainc_to_food => null()
    ! maintenance respiration fluxes
    ! leaf maintenance respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leaf_mr => null()
    ! fine root maintenance respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: froot_mr => null()
    ! live stem maintenance respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestem_mr => null()
    ! live coarse root maintenance respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecroot_mr => null()
    ! crop grain or organs maint. respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grain_mr => null()
    ! leaf maintenance respiration from current GPP (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leaf_curmr => null()
    ! fine root maintenance respiration from current GPP (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: froot_curmr => null()
    ! live stem maintenance respiration from current GPP (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestem_curmr => null()
    ! live coarse root maintenance respiration from current GPP (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecroot_curmr => null()
    ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grain_curmr => null()
    ! leaf maintenance respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leaf_xsmr => null()
    ! fine root maintenance respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: froot_xsmr => null()
    ! live stem maintenance respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestem_xsmr => null()
    ! live coarse root maintenance respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecroot_xsmr => null()
    ! crop grain or organs maint. respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grain_xsmr => null()
    ! photosynthesis fluxes
    ! C fixation from sunlit canopy (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: psnsun_to_cpool => null()
    ! C fixation from shaded canopy (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: psnshade_to_cpool => null()
    ! allocation fluxes, from current GPP
    ! allocation to maintenance respiration storage pool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_xsmrpool => null()
    ! allocation to grain C for prognostic crop(gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_grainc => null()
    ! allocation to grain C storage for prognostic crop(gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_grainc_storage => null()
    ! allocation to leaf C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_leafc => null()
    ! allocation to leaf C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_leafc_storage => null()
    ! allocation to fine root C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_frootc => null()
    ! allocation to fine root C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_frootc_storage => null()
    ! allocation to live stem C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_livestemc => null()
    ! allocation to live stem C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_livestemc_storage => null()
    ! allocation to dead stem C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_deadstemc => null()
    ! allocation to dead stem C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_deadstemc_storage => null()
    ! allocation to live coarse root C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_livecrootc => null()
    ! allocation to live coarse root C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_livecrootc_storage => null()
    ! allocation to dead coarse root C (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_deadcrootc => null()
    ! allocation to dead coarse root C storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_deadcrootc_storage => null()
    ! allocation to growth respiration storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_to_gresp_storage => null()
    ! growth respiration fluxes
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: xsmrpool_to_atm => null()
    ! leaf growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_leaf_gr => null()
    ! leaf growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_leaf_storage_gr => null()
    ! leaf growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_leaf_gr => null()
    ! fine root growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_froot_gr => null()
    ! fine root  growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_froot_storage_gr => null()
    ! fine root  growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_froot_gr => null()
    ! live stem growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_livestem_gr => null()
    ! live stem growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_livestem_storage_gr => null()
    ! live stem growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_livestem_gr => null()
    ! dead stem growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_deadstem_gr => null()
    ! dead stem growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_deadstem_storage_gr => null()
    ! dead stem growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_deadstem_gr => null()
    ! live coarse root growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_livecroot_gr => null()
    ! live coarse root growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_livecroot_storage_gr => null()
    ! live coarse root growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_livecroot_gr => null()
    ! dead coarse root growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_deadcroot_gr => null()
    ! dead coarse root growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_deadcroot_storage_gr => null()
    ! dead coarse root growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_deadcroot_gr => null()
    ! growth respiration for prognostic crop model
    ! grain growth respiration (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_grain_gr => null()
    ! grain growth respiration to storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: cpool_grain_storage_gr => null()
    ! grain growth respiration from storage (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: transfer_grain_gr => null()
    ! annual turnover of storage to transfer pools
    ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainc_storage_to_xfer => null()
    ! leaf C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafc_storage_to_xfer => null()
    ! fine root C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootc_storage_to_xfer => null()
    ! live stem C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_storage_to_xfer => null()
    ! dead stem C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadstemc_storage_to_xfer => null()
    ! live coarse root C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc_storage_to_xfer => null()
    ! dead coarse root C shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootc_storage_to_xfer => null()
    ! growth respiration shift storage to transfer (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: gresp_storage_to_xfer => null()
    ! turnover of livewood to deadwood
    ! live stem C turnover (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemc_to_deadstemc => null()
    ! live coarse root C turnover (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootc_to_deadcrootc => null()
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) gross primary production
    real(rk8), pointer, contiguous, dimension(:) :: gpp => null()
    ! (gC/m2/s) maintenance respiration
    real(rk8), pointer, contiguous, dimension(:) :: mr => null()
    ! (gC/m2/s) growth resp for new growth displayed in this timestep
    real(rk8), pointer, contiguous, dimension(:) :: current_gr => null()
    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
    real(rk8), pointer, contiguous, dimension(:) :: transfer_gr => null()
    ! (gC/m2/s) growth resp for growth sent to storage for later display
    real(rk8), pointer, contiguous, dimension(:) :: storage_gr => null()
    ! (gC/m2/s) total growth respiration
    real(rk8), pointer, contiguous, dimension(:) :: gr => null()
    ! (gC/m2/s) autotrophic respiration (MR + GR)
    real(rk8), pointer, contiguous, dimension(:) :: ar => null()
    ! (gC/m2/s) root respiration (fine root MR + total root GR)
    real(rk8), pointer, contiguous, dimension(:) :: rr => null()
    ! (gC/m2/s) net primary production
    real(rk8), pointer, contiguous, dimension(:) :: npp => null()
    ! (gC/m2/s) aboveground NPP
    real(rk8), pointer, contiguous, dimension(:) :: agnpp => null()
    ! (gC/m2/s) belowground NPP
    real(rk8), pointer, contiguous, dimension(:) :: bgnpp => null()
    ! (gC/m2/s) litterfall (leaves and fine roots)
    real(rk8), pointer, contiguous, dimension(:) :: litfall => null()
    ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
    real(rk8), pointer, contiguous, dimension(:) :: vegfire => null()
    ! (gC/m2/s) pft-level wood harvest (to product pools)
    real(rk8), pointer, contiguous, dimension(:) :: wood_harvestc => null()
    ! (gC/m2/s) pft-level carbon inputs (for balance checking)
    real(rk8), pointer, contiguous, dimension(:) :: pft_cinputs => null()
    ! (gC/m2/s) pft-level carbon outputs (for balance checking)
    real(rk8), pointer, contiguous, dimension(:) :: pft_coutputs => null()
#if (defined CN)
    ! CLAMP summary (diagnostic) variables, not involved in mass balance
    ! (gC/m2/s) pft-level fine root C alloc
    real(rk8), pointer, contiguous, dimension(:) :: frootc_alloc
    ! (gC/m2/s) pft-level fine root C loss
    real(rk8), pointer, contiguous, dimension(:) :: frootc_loss
    ! (gC/m2/s) pft-level leaf C alloc
    real(rk8), pointer, contiguous, dimension(:) :: leafc_alloc
    ! (gC/m2/s) pft-level leaf C loss
    real(rk8), pointer, contiguous, dimension(:) :: leafc_loss
    ! (gC/m2/s) pft-level wood C alloc
    real(rk8), pointer, contiguous, dimension(:) :: woodc_alloc
    ! (gC/m2/s) pft-level wood C loss
    real(rk8), pointer, contiguous, dimension(:) :: woodc_loss
#endif
    ! new variables for fire code
    ! (gC/m2/s) total pft-level fire C loss
    real(rk8), pointer, contiguous, dimension(:) :: pft_fire_closs => null()
! For aerenchyma calculations in CH4 code
#if (defined LCH4)
    ! (gC/m2/s) annual average aboveground NPP
    real(rk8), pointer, contiguous, dimension(:) :: annavg_agnpp
    ! (gC/m2/s) annual average belowground NPP
    real(rk8), pointer, contiguous, dimension(:) :: annavg_bgnpp
    ! (gC/m2/s) temp. average aboveground NPP
    real(rk8), pointer, contiguous, dimension(:) :: tempavg_agnpp
    ! (gC/m2/s) temp. average belowground NPP
    real(rk8), pointer, contiguous, dimension(:) :: tempavg_bgnpp
#endif
  end type pft_cflux_type

  !----------------------------------------------------
  ! pft nitrogen flux variables structure
  !----------------------------------------------------
  type, public :: pft_nflux_type
    ! new variables for CN code
    ! gap mortality fluxes
    ! leaf N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_to_litter => null()
    ! fine root N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_to_litter => null()
    ! leaf N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_storage_to_litter => null()
    ! fine root N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_storage_to_litter => null()
    ! live stem N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_storage_to_litter => null()
    ! dead stem N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_storage_to_litter => null()
    ! live coarse root N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_storage_to_litter => null()
    ! dead coarse root N storage mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_storage_to_litter => null()
    ! leaf N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_xfer_to_litter => null()
    ! fine root N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_xfer_to_litter => null()
    ! live stem N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_xfer_to_litter => null()
    ! dead stem N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_xfer_to_litter => null()
    ! live coarse root N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_xfer_to_litter => null()
    ! dead coarse root N transfer mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_xfer_to_litter => null()
    ! live stem N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_to_litter => null()
    ! dead stem N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_to_litter => null()
    ! live coarse root N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_to_litter => null()
    ! dead coarse root N mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_to_litter => null()
    ! retranslocated N pool mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: m_retransn_to_litter => null()
    ! harvest mortality fluxes
    ! leaf N harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafn_to_litter => null()
    ! fine root N harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootn_to_litter => null()
    ! leaf N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafn_storage_to_litter => null()
    ! fine root N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootn_storage_to_litter => null()
    ! live stem N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemn_storage_to_litter => null()
    ! dead stem N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_storage_to_litter => null()
    ! live coarse root N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootn_storage_to_litter => null()
    ! dead coarse root N storage harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootn_storage_to_litter => null()
    ! leaf N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_leafn_xfer_to_litter => null()
    ! fine root N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_frootn_xfer_to_litter => null()
    ! live stem N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemn_xfer_to_litter => null()
    ! dead stem N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_xfer_to_litter => null()
    ! live coarse root N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootn_xfer_to_litter => null()
    ! dead coarse root N transfer harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootn_xfer_to_litter => null()
    ! live stem N harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livestemn_to_litter => null()
    ! dead stem N harvest to 10-year product pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_to_prod10n => null()
    ! dead stem N harvest to 100-year product pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_to_prod100n => null()
    ! live coarse root N harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_livecrootn_to_litter => null()
    ! dead coarse root N harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadcrootn_to_litter => null()
    ! retranslocated N pool harvest mortality (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_retransn_to_litter => null()
    ! PFT-level fire N fluxes added by F. Li and S. Levis
    ! (gN/m2/s) fire N emissions from leafn
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_to_fire => null()
    ! (gN/m2/s) fire N emissions from leafn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from leafn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from livestemn
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_to_fire => null()
    ! (gN/m2/s) fire N emissions from livestemn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from livestemn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadstemn
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadstemn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadstemn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from frootn
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_to_fire => null()
    ! (gN/m2/s) fire N emissions from frootn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from frootn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_to_fire => null()
    ! (gN/m2/s) fire N emissions from livecrootn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from livecrootn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadcrootn
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadcrootn_storage
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_storage_to_fire => null()
    ! (gN/m2/s) fire N emissions from deadcrootn_xfer
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_xfer_to_fire => null()
    ! (gN/m2/s) fire N emissions from retransn
    real(rk8), pointer, contiguous, dimension(:) :: m_retransn_to_fire => null()
    ! (gN/m2/s) from leafn to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_to_litter_fire => null()
    ! (gN/m2/s) from leafn_storage to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_storage_to_litter_fire => null()
    ! (gN/m2/s) from leafn_xfer to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_leafn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from livestemn to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_to_litter_fire => null()
    ! (gN/m2/s) from livestemn_storage to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_storage_to_litter_fire => null()
    ! (gN/m2/s) from livestemn_xfer to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from livestemn to deadstemn N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livestemn_to_deadstemn_fire => null()
    ! (gN/m2/s) from deadstemn to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_to_litter_fire => null()
    ! (gN/m2/s) from deadstemn_storage to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_storage_to_litter_fire => null()
    ! (gN/m2/s) from deadstemn_xfer to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadstemn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from frootn to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_to_litter_fire => null()
    ! (gN/m2/s) from frootn_storage to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_storage_to_litter_fire => null()
    ! (gN/m2/s) from frootn_xfer to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_frootn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from livecrootn to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_to_litter_fire => null()
    ! (gN/m2/s) from livecrootn_storage to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_storage_to_litter_fire => null()
    ! (gN/m2/s) from livecrootn_xfer to litter N  due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_livecrootn_to_deadcrootn_fire => null()
    ! (gN/m2/s) from deadcrootn to deadcrootn due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_to_litter_fire => null()
    ! (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_storage_to_litter_fire => null()
    ! (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_deadcrootn_xfer_to_litter_fire => null()
    ! (gN/m2/s) from retransn to deadcrootn due to fire
    real(rk8), pointer, contiguous, dimension(:) :: m_retransn_to_litter_fire => null()
    ! phenology fluxes from transfer pool
    ! grain N growth from storage for prognostic crop model (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainn_xfer_to_grainn => null()
    ! leaf N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafn_xfer_to_leafn => null()
    ! fine root N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootn_xfer_to_frootn => null()
    ! live stem N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_xfer_to_livestemn => null()
    ! dead stem N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadstemn_xfer_to_deadstemn => null()
    ! live coarse root N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_xfer_to_livecrootn => null()
    ! dead coarse root N growth from storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootn_xfer_to_deadcrootn => null()
    ! litterfall fluxes
    ! livestem N to litter (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_to_litter => null()
    ! grain N to food for prognostic crop (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainn_to_food => null()
    ! leaf N litterfall (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafn_to_litter => null()
    ! leaf N to retranslocated N pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafn_to_retransn => null()
    ! fine root N to retranslocated N pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootn_to_retransn => null()
    ! fine root N litterfall (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootn_to_litter => null()
    ! allocation fluxes
    ! deployment of retranslocated N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: retransn_to_npool => null()
    ! deployment of soil mineral N uptake (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: sminn_to_npool => null()
    ! allocation to grain N for prognostic crop (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_grainn => null()
    ! allocation to grain N storage for prognostic crop (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_grainn_storage => null()
    ! allocation to leaf N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_leafn => null()
    ! allocation to leaf N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_leafn_storage => null()
    ! allocation to fine root N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_frootn => null()
    ! allocation to fine root N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_frootn_storage => null()
    ! allocation to live stem N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_livestemn => null()
    ! allocation to live stem N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_livestemn_storage => null()
    ! allocation to dead stem N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_deadstemn => null()
    ! allocation to dead stem N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_deadstemn_storage => null()
    ! allocation to live coarse root N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_livecrootn => null()
    ! allocation to live coarse root N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_livecrootn_storage => null()
    ! allocation to dead coarse root N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_deadcrootn => null()
    ! allocation to dead coarse root N storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: npool_to_deadcrootn_storage => null()
    ! annual turnover of storage to transfer pools
    ! grain N shift storage to transfer for prognostic crop (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: grainn_storage_to_xfer => null()
    ! leaf N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: leafn_storage_to_xfer => null()
    ! fine root N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: frootn_storage_to_xfer => null()
    ! live stem N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_storage_to_xfer => null()
    ! dead stem N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadstemn_storage_to_xfer => null()
    ! live coarse root N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_storage_to_xfer => null()
    ! dead coarse root N shift storage to transfer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: deadcrootn_storage_to_xfer => null()
    ! applied fertilizer (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: fert => null()
    ! soybean fixed N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: soyfixn => null()
    ! turnover of livewood to deadwood, with retranslocation
    ! live stem N turnover (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_to_deadstemn => null()
    ! live stem N to retranslocated N pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livestemn_to_retransn => null()
    ! live coarse root N turnover (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_to_deadcrootn => null()
    ! live coarse root N to retranslocated N pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: livecrootn_to_retransn => null()
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! total N deployed to growth and storage (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ndeploy => null()
    ! total N inputs to pft-level (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: pft_ninputs => null()
    ! total N outputs from pft-level (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: pft_noutputs => null()
    ! total N losses to wood product pools (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: wood_harvestn => null()
    ! new variables for fire code
    ! total pft-level fire N loss (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: pft_fire_nloss => null()
  end type pft_nflux_type

  !----------------------------------------------------
  ! pft VOC fluxes structure for history output
  !----------------------------------------------------
  type, public :: megan_out_type
    !(n_megan_comps) MEGAN flux [ug C m-2 h-1]
    real(rk8), pointer, contiguous, dimension(:) :: flux_out => null()
  end type megan_out_type

  !----------------------------------------------------
  ! pft VOC flux variables structure
  !----------------------------------------------------
  type, public :: pft_vflux_type
    !total VOC flux into atmosphere [moles/m2/sec]
    real(rk8), pointer, contiguous, dimension(:) :: vocflx_tot => null()
    !(num_mech_comps) MEGAN flux [moles/m2/sec]
    real(rk8), pointer, contiguous, dimension(:,:) :: vocflx => null()
    real(rk8), pointer, contiguous, dimension(:) :: Eopt_out => null()   !Eopt coefficient
    real(rk8), pointer, contiguous, dimension(:) :: topt_out => null()   !topt coefficient
    real(rk8), pointer, contiguous, dimension(:) :: alpha_out => null()  !alpha coefficient
    real(rk8), pointer, contiguous, dimension(:) :: cp_out => null()     !cp coefficient
    real(rk8), pointer, contiguous, dimension(:) :: paru_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: par24u_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: par240u_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: para_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: par24a_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: par240a_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gamma_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaL_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaT_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaP_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaA_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaS_out => null()
    real(rk8), pointer, contiguous, dimension(:) :: gammaC_out => null()
    ! points to output fluxes
    type(megan_out_type), pointer, dimension(:) :: meg => null()
  end type pft_vflux_type

  !----------------------------------------------------
  ! pft dry dep velocity variables structure
  !----------------------------------------------------
  type, public :: pft_depvd_type
    real(rk8), pointer, contiguous, dimension(:,:) :: drydepvel => null()
  end type pft_depvd_type

  !----------------------------------------------------
  ! pft dust flux variables structure
  !----------------------------------------------------
  type, public :: pft_dflux_type
    !Fraction of bare ground emitting dust
    real(rk8), pointer, contiguous, dimension(:) :: lnd_frc_mbl_dst => null()
    !(ndst)  !surface dust emission (kg/m**2/s) [ + = to atm]
    real(rk8), pointer, contiguous, dimension(:,:) :: flx_mss_vrt_dst => null()
    !total dust flux into atmosphere
    real(rk8), pointer, contiguous, dimension(:) :: flx_mss_vrt_dst_tot => null()
    !(ndst) turbulent deposition velocity (m/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: vlc_trb => null()
    !turbulent deposition velocity 1(m/s)
    real(rk8), pointer, contiguous, dimension(:) :: vlc_trb_1 => null()
    !turbulent deposition velocity 2(m/s)
    real(rk8), pointer, contiguous, dimension(:) :: vlc_trb_2 => null()
    !turbulent deposition velocity 3(m/s)
    real(rk8), pointer, contiguous, dimension(:) :: vlc_trb_3 => null()
    !turbulent deposition velocity 4(m/s)
    real(rk8), pointer, contiguous, dimension(:) :: vlc_trb_4 => null()
  end type pft_dflux_type

  !----------------------------------------------------
  ! End definition of structures defined at the pft_type level
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Begin definition of structures defined at the column_type level
  !----------------------------------------------------

  !----------------------------------------------------
  ! column physical state variables structure
  !----------------------------------------------------
  type, public :: column_pstate_type
    !pft-level pstate variables averaged to the column
    type(pft_pstate_type) :: pps_a
    !number of snow layers
    integer(ik4), pointer, contiguous, dimension(:) :: snl => null()
    !soil color class
    integer(ik4), pointer, contiguous, dimension(:) :: isoicol => null()
#ifdef HAMSTER_ALBEDO
    real(rk8), pointer, contiguous, dimension(:,:) :: hamster_alb
#endif
    !F. Li and S. Levis
    ! global real gdp data (k US$/capita)
    real(rk8), pointer, contiguous, dimension(:) :: gdp_lf => null()
    ! global peatland fraction data (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: peatf_lf => null()
    ! global peak month of crop fire emissions
    integer(ik4), pointer, contiguous, dimension(:) :: abm_lf => null()
    !gdp limitation factor for fire occurrence (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: lgdp_col => null()
    !gdp limitation factor for fire spreading (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: lgdp1_col => null()
    !pop limitation factor for fire spreading (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: lpop_col => null()
    !Clapp and Hornberger "b" (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: bsw => null()
    !volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: watsat => null()
    !btran parameter for btran=0
    real(rk8), pointer, contiguous, dimension(:,:) :: watdry => null()
    !btran parameter for btran = 1
    real(rk8), pointer, contiguous, dimension(:,:) :: watopt => null()
    !hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: hksat => null()
    !mineral hksat
    real(rk8), pointer, contiguous, dimension(:,:) :: hksat_min => null()
    !thermal conductivity
    real(rk8), pointer, contiguous, dimension(:,:) :: tk_hist => null()
    !heat capacity
    real(rk8), pointer, contiguous, dimension(:,:) :: cv_hist => null()
    !minimum soil suction (mm) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: sucsat => null()
    !decay factor (m)
    real(rk8), pointer, contiguous, dimension(:) :: hkdepth => null()
    !maximum saturated fraction for a gridcell
    real(rk8), pointer, contiguous, dimension(:) :: wtfact => null()
    !fractional impermeability (-)
    real(rk8), pointer, contiguous, dimension(:,:) :: fracice => null()
    !heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: csol => null()
    !thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: tkmg => null()
    !thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: tkdry => null()
    !thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: tksatu => null()
    !restriction for min of soil potential (mm) (new)
    real(rk8), pointer, contiguous, dimension(:) :: smpmin => null()
    !threshold soil moisture based on clay content
    real(rk8), pointer, contiguous, dimension(:) :: gwc_thr => null()
    ![frc] Mass fraction clay limited to 0.20
    real(rk8), pointer, contiguous, dimension(:) :: mss_frc_cly_vld => null()
    !basin factor
    real(rk8), pointer, contiguous, dimension(:) :: mbl_bsn_fct => null()
    !true => do snow capping
    logical, pointer, contiguous, dimension(:) :: do_capsnow => null()
    !snow height of snow covered area (m)
    real(rk8), pointer, contiguous, dimension(:) :: snow_depth => null()
    ! gridcell averaged snow height (m)
    real(rk8), pointer, contiguous, dimension(:) :: snowdp => null()
    !fraction of ground covered by snow (0 to 1)
    real(rk8), pointer, contiguous, dimension(:) :: frac_sno => null()
    !fraction of ground covered by snow (0 to 1)
    real(rk8), pointer, contiguous, dimension(:) :: frac_sno_eff => null()
    !fractional area with surface water greater than zero
    real(rk8), pointer, contiguous, dimension(:) :: frac_h2osfc => null()
    !temporay fractional area with surface water greater than zero
    real(rk8), pointer, contiguous, dimension(:) :: frac_h2osfc_temp => null()
    !gridcell topographic standard deviation (m)
    real(rk8), pointer, contiguous, dimension(:) :: topo_std => null()
    !gridcell topographic index
    real(rk8), pointer, contiguous, dimension(:) :: topo_ndx => null()
    !gridcell topographic slope
    real(rk8), pointer, contiguous, dimension(:) :: topo_slope => null()
    ! microtopography pdf sigma (m)
    real(rk8), pointer, contiguous, dimension(:) :: micro_sigma => null()
    ! level at which h2osfc "percolates"
    real(rk8), pointer, contiguous, dimension(:) :: h2osfc_thresh => null()
    ! SCA shape parameter
    real(rk8), pointer, contiguous, dimension(:) :: n_melt => null()
    !interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: zi => null()
    !layer thickness (m)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: dz => null()
    !layer depth (m) (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: z => null()
    !fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: frac_iceold => null()
    !flag for melting (=1), freezing (=2), Not=0 (new) (-nlevsno+1:nlevgrnd)
    integer(ik4), pointer, contiguous, dimension(:,:) :: imelt => null()
    !effective porosity = porosity - vol_ice (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: eff_porosity => null()
    !ground emissivity
    real(rk8), pointer, contiguous, dimension(:) :: emg => null()
    !roughness length over ground, momentum [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0mg => null()
    !roughness length over ground, sensible heat [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0hg => null()
    !roughness length over ground, latent heat [m]
    real(rk8), pointer, contiguous, dimension(:) :: z0qg => null()
    !latent heat of vapor of water (or sublimation) [j/kg]
    real(rk8), pointer, contiguous, dimension(:) :: htvp => null()
    !coefficient of convective velocity [-]
    real(rk8), pointer, contiguous, dimension(:) :: beta => null()
    !convective boundary height [m]
    real(rk8), pointer, contiguous, dimension(:) :: zii => null()
    !ground albedo (direct) (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgrd => null()
    !ground albedo (diffuse) (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgri => null()
    !effective fraction of roots in each soil layer (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: rootr_column => null()
    !fraction of roots in each soil layer for urban pervious road
    real(rk8), pointer, contiguous, dimension(:,:) :: rootfr_road_perv => null()
    !effective fraction of roots in each soil layer of urban pervious road
    real(rk8), pointer, contiguous, dimension(:,:) :: rootr_road_perv => null()
    !soil water as frac. of whc for top 0.05 m (0-1)
    ! (only comment changed by F. Li and S. Levis)
    real(rk8), pointer, contiguous, dimension(:) :: wf => null()
    !soil water as frac. of whc for top 0.17 m (0-1)
    ! added by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: wf2 => null()
    !irrigation rate
!   real(rk8), pointer, contiguous, dimension(:) :: xirrig
    !maximum daylength for this column (s)
    real(rk8), pointer, contiguous, dimension(:) :: max_dayl => null()
#if (defined VICHYDRO)
    !b infiltration parameter
    real(rk8), pointer, contiguous, dimension(:) :: b_infil
    !fracton of Dsmax where non-linear baseflow begins
    real(rk8), pointer, contiguous, dimension(:) :: ds
    !max. velocity of baseflow (mm/day)
    real(rk8), pointer, contiguous, dimension(:) :: dsmax
    !fraction of maximum soil moisutre where non-liear base flow occurs
    real(rk8), pointer, contiguous, dimension(:) :: wsvic
    !baseflow exponent (Qb)
    real(rk8), pointer, contiguous, dimension(:) :: c_param
    !pore-size distribution related paramter(Q12)
    real(rk8), pointer, contiguous, dimension(:,:) :: expt
    !Saturated hydrologic conductivity
    real(rk8), pointer, contiguous, dimension(:,:) :: ksat
    !soil moisture dissusion parameter
    real(rk8), pointer, contiguous, dimension(:,:) :: phi_s
    !layer depth of upper layer
    real(rk8), pointer, contiguous, dimension(:,:) :: depth
    !soil porisity (1-bulk_density/soil_density)
    real(rk8), pointer, contiguous, dimension(:,:) :: porosity
    !max layer moist + ice (mm)
    real(rk8), pointer, contiguous, dimension(:,:) :: max_moist
    !fraction of VIC layers in CLM layers
    real(rk8), pointer, contiguous, dimension(:,:,:) :: vic_clm_fract
#endif
    ! new variables for CN code
    ! solar declination angle (radians)
    real(rk8), pointer, contiguous, dimension(:) :: decl => null()
    ! cosine of solar zenith angle
    real(rk8), pointer, contiguous, dimension(:) :: coszen => null()
    ! soil water potential in each soil layer (MPa)
    real(rk8), pointer, contiguous, dimension(:,:) :: soilpsi => null()
    ! bulk density of dry soil material [kg/m^3]
    real(rk8), pointer, contiguous, dimension(:,:) :: bd => null()
    ! fraction of potential immobilization (no units)
    real(rk8), pointer, contiguous, dimension(:,:) :: fpi_vr => null()
    ! fraction of potential immobilization (no units)
    real(rk8), pointer, contiguous, dimension(:) :: fpi => null()
    ! respired fraction in decomposition step (frac)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: rf_decomp_cascade => null()
    ! what fraction of C leaving a given pool passes through a given
    ! transition (frac)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: pathfrac_decomp_cascade => null()
    ! (1/m) profile for N fixation additions
    real(rk8), pointer, contiguous, dimension(:,:) :: nfixation_prof => null()
    ! (1/m) profile for N fixation additions
    real(rk8), pointer, contiguous, dimension(:,:) :: ndep_prof => null()
    ! current depth of thaw
    real(rk8), pointer, contiguous, dimension(:) :: alt => null()
    ! maximum annual depth of thaw
    real(rk8), pointer, contiguous, dimension(:) :: altmax => null()
    ! prior year maximum annual depth of thaw
    real(rk8), pointer, contiguous, dimension(:) :: altmax_lastyear => null()
    ! current depth of thaw
    integer(ik4), pointer, contiguous, dimension(:) :: alt_indx => null()
    ! maximum annual depth of thaw
    integer(ik4), pointer, contiguous, dimension(:) :: altmax_indx => null()
    ! prior year maximum annual depth of thaw
    integer(ik4), pointer, contiguous, dimension(:) :: altmax_lastyear_indx => null()
    ! SOM advective flux (m/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: som_adv_coef => null()
    ! SOM diffusivity due to bio/cryo-turbation (m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: som_diffus_coef => null()
#ifdef NITRIF_DENITRIF
    ! maximumn monthly-mean soil temperature
    real(rk8), pointer, contiguous, dimension(:,:) :: tmean_monthly_max_vr
    ! monthly-mean soil temperature
    real(rk8), pointer, contiguous, dimension(:,:) :: tmean_monthly_vr
#endif
    !fraction of potential gpp (no units)
    real(rk8), pointer, contiguous, dimension(:) :: fpg => null()
    !seconds since last annual accumulator turnover
    real(rk8), pointer, contiguous, dimension(:) :: annsum_counter => null()
    !annual sum of NPP, averaged from pft-level (gC/m2/yr)
    real(rk8), pointer, contiguous, dimension(:) :: cannsum_npp => null()
    ! (gC/m2/s) lagged net primary production
    real(rk8), pointer, contiguous, dimension(:) :: col_lag_npp => null()
    !annual average of 2m air temperature, averaged from pft-level (K)
    real(rk8), pointer, contiguous, dimension(:) :: cannavg_t2m => null()
    !volumetric soil water at field capacity (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: watfc => null()
    ! F. Li and S. Levis
    ! fire counts (count/km2/timestep), valid only in Reg. C
    real(rk8), pointer, contiguous, dimension(:) :: nfire => null()
    ! fire spread rate in pft level (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_pft => null()
    ! fire spread rate at column level (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: fsr_col => null()
    ! fire duration at column level (hr)
    real(rk8), pointer, contiguous, dimension(:) :: fd_col => null()
    ! fire duration in pft level    (hr)
    real(rk8), pointer, contiguous, dimension(:) :: fd_pft => null()
    !60-day running mean of tot. precipitation (mm/s)
    real(rk8), pointer, contiguous, dimension(:) :: prec60_col => null()
    !10-day running mean of tot. precipitation (mm/s)
    real(rk8), pointer, contiguous, dimension(:) :: prec10_col => null()
    ! conversion area fraction of BET and BDT that haven't burned before (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: lfc => null()
    ! conversion area fraction of BET and BDT that burned in this
    ! timestep ((timestep)-1)
    real(rk8), pointer, contiguous, dimension(:) :: lfc2 => null()
    ! annual decreased fraction coverage of BET on the gridcell (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: dtrotr_col => null()
    ! pft weight of BET and BDT on the gridcell(0-1)
    real(rk8), pointer, contiguous, dimension(:) :: trotr1_col => null()
    ! pft weight of BDT on the gridcell (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: trotr2_col => null()
    ! crop fraction in veg column (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: cropf_col => null()
    ! baf for cropland per time step(0-1)
    real(rk8), pointer, contiguous, dimension(:) :: baf_crop => null()
    ! baf for peatland per time step (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: baf_peatf => null()
    ! total burned area out of conversion (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: fbac => null()
    ! burned area out of conversion region due to land use fire (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: fbac1 => null()
    ! btran2 at column level (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: btran_col => null()
    ! fractional coverage of non-crop PFTs (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: wtlf => null()
    ! fractional coverage of non-crop and non-bare-soil PFTs (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: lfwt => null()
    !timestep fractional area burned (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: farea_burned => null()
    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: albsnd_hst => null()
    ! snow albedo, diffuse, for history files (col,bnd) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: albsni_hst => null()
    ! soil albedo: direct (col,bnd) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: albsod => null()
    ! soil albedo: diffuse (col,bnd) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: albsoi => null()
    ! absorbed flux per unit incident direct flux: VIS (col,lyr) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: flx_absdv => null()
    ! absorbed flux per unit incident direct flux: NIR (col,lyr) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: flx_absdn => null()
    ! absorbed flux per unit incident diffuse flux: VIS (col,lyr) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: flx_absiv => null()
    ! absorbed flux per unit incident diffuse flux: NIR (col,lyr) [frc]
    real(rk8), pointer, contiguous, dimension(:,:) :: flx_absin => null()
    ! snow grain radius (col,lyr) [m^-6, microns]
    real(rk8), pointer, contiguous, dimension(:,:) :: snw_rds => null()
    ! snow grain radius, top layer (col) [m^-6, microns]
    real(rk8), pointer, contiguous, dimension(:) :: snw_rds_top => null()
    ! snow liquid water fraction (mass), top layer (col) [fraction]
    real(rk8), pointer, contiguous, dimension(:) :: sno_liq_top => null()
    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_bcpho => null()
    ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_bcphi => null()
    ! total mass of BC in snow (pho+phi) (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_bctot => null()
    ! column-integrated mass of total BC (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_bc_col => null()
    ! top-layer mass of total BC (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_bc_top => null()
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_ocpho => null()
    ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_ocphi => null()
    ! total mass of OC in snow (pho+phi) (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_octot => null()
    ! column-integrated mass of total OC (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_oc_col => null()
    ! top-layer mass of total OC (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_oc_top => null()
    ! mass of dust species 1 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_dst1 => null()
    ! mass of dust species 2 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_dst2 => null()
    ! mass of dust species 3 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_dst3 => null()
    ! mass of dust species 4 in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_dst4 => null()
    ! total mass of dust in snow (col,lyr) [kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_dsttot => null()
    ! column-integrated mass of dust in snow (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_dst_col => null()
    ! top-layer mass of dust in snow (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: mss_dst_top => null()
    ! top-layer mass of snow (col) [kg]
    real(rk8), pointer, contiguous, dimension(:) :: h2osno_top => null()
    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_bcphi => null()
    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_bcpho => null()
    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_ocphi => null()
    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_ocpho => null()
    ! mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_dst1 => null()
    ! mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_dst2 => null()
    ! mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_dst3 => null()
    ! mass concentration of dust species 4 in snow (col,lyr) [kg/kg]
    real(rk8), pointer, contiguous, dimension(:,:) :: mss_cnc_dst4 => null()
    ! pure snow ground direct albedo (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgrd_pur => null()
    ! pure snow ground diffuse albedo (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgri_pur => null()
    ! ground direct albedo without BC  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgrd_bc => null()
    ! ground diffuse albedo without BC (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgri_bc => null()
    ! ground direct albedo without OC  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgrd_oc => null()
    ! ground diffuse albedo without OC (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgri_oc => null()
    ! ground direct albedo without dust  (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgrd_dst => null()
    ! ground diffuse albedo without dust (numrad)
    real(rk8), pointer, contiguous, dimension(:,:) :: albgri_dst => null()
    ! temperature gradient in top layer  [K m-1]
    real(rk8), pointer, contiguous, dimension(:) :: dTdz_top => null()
    ! temperature of top snow layer [K]
    real(rk8), pointer, contiguous, dimension(:) :: snot_top => null()
    ! new variables for S Lake code
    ! surface friction velocity (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: ws => null()
    ! coefficient for calculation of decay of eddy diffusivity with depth
    real(rk8), pointer, contiguous, dimension(:) :: ks => null()
    ! lake layer thickness (m)  (1:nlevlak)
    real(rk8), pointer, contiguous, dimension(:,:) :: dz_lake => null()
    ! layer depth for lake (m)
    real(rk8), pointer, contiguous, dimension(:,:) :: z_lake => null()
    ! top level eddy conductivity from previous timestep (W/mK)
    real(rk8), pointer, contiguous, dimension(:) :: savedtke1 => null()
    ! sand value for gridcell containing column (1:nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: cellsand => null()
    ! clay value for gridcell containing column (1:nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: cellclay => null()
    ! organic matter for gridcell containing column (1:nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: cellorg => null()
    ! variable lake depth (m)
    real(rk8), pointer, contiguous, dimension(:) :: lakedepth => null()
    ! lake extinction coefficient from surface data (1/m)
    real(rk8), pointer, contiguous, dimension(:) :: etal => null()
    ! lake fetch from surface data (m)
    real(rk8), pointer, contiguous, dimension(:) :: lakefetch => null()
    ! friction velocity (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: ust_lake => null()
    ! end new variables for S Lake code
    ! New variables for finundated in methane code
#ifdef LCH4
    ! coefficient for determining finundated (m)
    !real(rk8), pointer, contiguous, dimension(:) :: zwt0
    ! maximum inundated fraction for a gridcell (for methane code)
    !real(rk8), pointer, contiguous, dimension(:) :: f0
    ! coefficient for determining finundated (m)
    !real(rk8), pointer, contiguous, dimension(:) :: p3
    real(rk8), pointer, contiguous, dimension(:) :: k
    real(rk8), pointer, contiguous, dimension(:) :: q
    real(rk8), pointer, contiguous, dimension(:) :: v
    real(rk8), pointer, contiguous, dimension(:) :: maxf
    ! added by Lei Meng for pH effects of methane production
    ! pH values
    real(rk8), pointer, contiguous, dimension(:) :: pH
#endif
    ! End New variables for methane code
    ! current irrigation rate [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: irrig_rate => null()
    ! number of time steps for which we still need to irrigate today
    ! (if 0, ignore irrig_rate)
    integer(ik4), pointer , dimension(:) :: n_irrig_steps_left => null()
    ! surface atm pressure, downscaled to column (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_pbot => null()
    ! surface air density, downscaled to column (kg/m^3)
    real(rk8), pointer, contiguous, dimension(:) :: forc_rho => null()
#ifdef CN
    real(rk8), pointer, contiguous, dimension(:) :: q10
    real(rk8), pointer, contiguous, dimension(:) :: ndep
#endif
  end type column_pstate_type

  !----------------------------------------------------
  ! column energy state variables structure
  !----------------------------------------------------
  type, public :: column_estate_type
    !pft-level energy state variables averaged to the column
    type(pft_estate_type) :: pes_a
    !ground temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_grnd => null()
    !Urban ground temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_grnd_u => null()
    !Rural ground temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_grnd_r => null()
    !change in t_grnd, last iteration (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: dt_grnd => null()
    !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: t_soisno => null()
    !soil temperature in top 10cm of soil (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_soi_10cm => null()
    !soil temperature in top 17cm of soil (Kelvin) by F. Li and S. Levis
    real(rk8), pointer, contiguous, dimension(:) :: tsoi17 => null()
    !lake temperature (Kelvin)  (1:nlevlak)
    real(rk8), pointer, contiguous, dimension(:,:) :: t_lake => null()
    !soil/snow temperature before update (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: tssbef => null()
    !virtual potential temperature (kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: thv => null()
    !soil heat content (MJ/m2)
    real(rk8), pointer, contiguous, dimension(:) :: hc_soi => null()
    !soil plus snow heat content (MJ/m2)
    real(rk8), pointer, contiguous, dimension(:) :: hc_soisno => null()
    !atm temperature, downscaled to column (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: forc_t => null()
    !atm potl temperature, downscaled to column (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: forc_th => null()
    !surface water temperature
    real(rk8), pointer, contiguous, dimension(:) :: t_h2osfc => null()
    !surface water temperature from time-step before
    real(rk8), pointer, contiguous, dimension(:) :: t_h2osfc_bef => null()
  end type column_estate_type

  !----------------------------------------------------
  ! column water state variables structure
  !----------------------------------------------------
  type, public :: column_wstate_type
    !pft-level water state variables averaged to the column
    type(pft_wstate_type) :: pws_a
    !surface water (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: h2osfc => null()
    !ground specific humidity [kg/kg]
    real(rk8), pointer, contiguous, dimension(:) :: qg_snow => null()
    !ground specific humidity [kg/kg]
    real(rk8), pointer, contiguous, dimension(:) :: qg_soil => null()
    !ground specific humidity [kg/kg]
    real(rk8), pointer, contiguous, dimension(:) :: qg_h2osfc => null()
    !initial snow water
    real(rk8), pointer, contiguous, dimension(:,:) :: swe_old => null()
    !snow water (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: h2osno => null()
    !imbalance in snow water (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: errh2osno => null()
    !snow sources (mm H2O/s)
    real(rk8), pointer, contiguous, dimension(:) :: snow_sources => null()
    !snow sinks (mm H2O/s)
    real(rk8), pointer, contiguous, dimension(:) :: snow_sinks => null()
    !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: h2osoi_liq => null()
    !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: h2osoi_ice => null()
    !liquid water + ice lens in top 10cm of soil (kg/m2)
    real(rk8), pointer, contiguous, dimension(:) :: h2osoi_liqice_10cm => null()
    !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: h2osoi_vol => null()
    !snow mass for previous time step (kg/m2) (new)
    real(rk8), pointer, contiguous, dimension(:) :: h2osno_old => null()
    !ground specific humidity [kg/kg]
    real(rk8), pointer, contiguous, dimension(:) :: qg => null()
    !d(qg)/dT
    real(rk8), pointer, contiguous, dimension(:) :: dqgdT => null()
    !average snow ice lens
    real(rk8), pointer, contiguous, dimension(:) :: snowice => null()
    !average snow liquid water
    real(rk8), pointer, contiguous, dimension(:) :: snowliq => null()
    !factor that reduces ground saturated specific humidity (-)
    real(rk8), pointer, contiguous, dimension(:) :: soilalpha => null()
    !factor that reduces ground evaporation L&P1992(-)
    real(rk8), pointer, contiguous, dimension(:) :: soilbeta => null()
    !urban factor that reduces ground saturated specific humidity (-)
    real(rk8), pointer, contiguous, dimension(:) :: soilalpha_u => null()
    !water table depth
    real(rk8), pointer, contiguous, dimension(:) :: zwt => null()
    !frost table depth
    real(rk8), pointer, contiguous, dimension(:) :: frost_table => null()
    !perched water table depth
    real(rk8), pointer, contiguous, dimension(:) :: zwt_perched => null()
    !integrated snowfall (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: int_snow => null()
    !fractional impermeable area
    real(rk8), pointer, contiguous, dimension(:) :: fcov => null()
    !water in the unconfined aquifer (mm)
    real(rk8), pointer, contiguous, dimension(:) :: wa => null()
    !aquifer recharge rate (mm/s)
    real(rk8), pointer, contiguous, dimension(:) :: qcharge => null()
    !soil matric potential (mm)
    real(rk8), pointer, contiguous, dimension(:,:) :: smp_l => null()
    !hydraulic conductivity (mm/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: hk_l => null()
    !fractional area with water table at surface
    real(rk8), pointer, contiguous, dimension(:) :: fsat => null()
    !atm specific humidity, downscaled to column (kg/kg)
    real(rk8), pointer, contiguous, dimension(:) :: forc_q => null()
#if (defined VICHYDRO)
    !soil moisture (kg/m2) for VIC soil layers
    real(rk8), pointer, contiguous, dimension(:,:) :: moist
    !soil ice (kg/m2) for VIC soil layers
    real(rk8), pointer, contiguous, dimension(:,:) :: ice
    !volumetric soil moisture for VIC soil layers
    real(rk8), pointer, contiguous, dimension(:,:) :: moist_vol
    !maximum infiltration rate calculated by VIC
    real(rk8), pointer, contiguous, dimension(:) :: max_infil
    !average saturation in top soil layers in VIC
    real(rk8), pointer, contiguous, dimension(:) :: i_0
#endif
#ifdef LCH4
    !fractional inundated area (excluding dedicated wetland columns)
    real(rk8), pointer, contiguous, dimension(:) :: finundated
#endif
    ! new variables for S Lake code
    ! mass fraction of lake layer that is frozen
    real(rk8), pointer, contiguous, dimension(:,:) :: lake_icefrac => null()
    ! ice thickness (m) (integrated if lakepuddling)
    real(rk8), pointer, contiguous, dimension(:) :: lake_icethick => null()
    ! end new variables for S Lake code
  end type column_wstate_type

  !----------------------------------------------------
  ! column carbon state variables structure
  !----------------------------------------------------
  type, public :: column_cstate_type
    !pft-level carbon state variables averaged to the column
    type(pft_cstate_type) :: pcs_a
    ! NOTE: the soilc variable is used by the original CLM C-cycle code,
    ! and is not used by the CN code
    !soil carbon (kg C /m**2)
    real(rk8), pointer, contiguous, dimension(:) :: soilc => null()
    ! all c pools involved in decomposition
    ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cpools_vr => null()
    ! (gC/m3) vertically-resolved column-level sink for C truncation
    real(rk8), pointer, contiguous, dimension(:,:) :: col_ctrunc_vr => null()
    !fire-related variables added by F. Li and S. Levis
    !root carbon at column level (gC/m2)
    real(rk8), pointer, contiguous, dimension(:) :: rootc_col => null()
    !column-level totvegc (gC/m2)
    real(rk8), pointer, contiguous, dimension(:) :: totvegc_col => null()
    !column-level leafc (gC/m2)
    real(rk8), pointer, contiguous, dimension(:) :: leafc_col => null()
    ! fuel avalability factor for Reg.C (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: fuelc => null()
    ! fuel avalability factor for Reg.A (0-1)
    real(rk8), pointer, contiguous, dimension(:) :: fuelc_crop => null()
    ! pools for dynamic landcover
    ! (gC/m2) column-level pool for seeding new PFTs
    real(rk8), pointer, contiguous, dimension(:) :: seedc => null()
    ! (gC/m2) wood product C pool, 10-year lifespan
    real(rk8), pointer, contiguous, dimension(:) :: prod10c => null()
    ! (gC/m2) wood product C pool, 100-year lifespan
    real(rk8), pointer, contiguous, dimension(:) :: prod100c => null()
    ! (gC/m2) total wood product C
    real(rk8), pointer, contiguous, dimension(:) :: totprodc => null()
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cpools => null()
    ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cpools_1m => null()
    ! (gC/m2) Diagnostic: coarse woody debris C
    real(rk8), pointer, contiguous, dimension(:) :: cwdc => null()
    ! (gC/m2) column-level sink for C truncation
    real(rk8), pointer, contiguous, dimension(:) :: col_ctrunc => null()
    ! (gC/m2) total litter carbon
    real(rk8), pointer, contiguous, dimension(:) :: totlitc => null()
    ! (gC/m2) total soil organic matter carbon
    real(rk8), pointer, contiguous, dimension(:) :: totsomc => null()
    ! (gC/m2) total litter carbon to 1 meter
    real(rk8), pointer, contiguous, dimension(:) :: totlitc_1m => null()
    ! (gC/m2) total soil organic matter carbon to 1 meter
    real(rk8), pointer, contiguous, dimension(:) :: totsomc_1m => null()
    ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
    real(rk8), pointer, contiguous, dimension(:) :: totecosysc => null()
    ! (gC/m2) total column carbon, incl veg and cpool
    real(rk8), pointer, contiguous, dimension(:) :: totcolc => null()
  end type column_cstate_type

#ifdef LCH4
  !----------------------------------------------------
  ! column methane variables structure
  !----------------------------------------------------
  type, public :: column_ch4_type
    ! new variables for CH4 code
    ! column-level methane fluxes
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_prod_depth_sat => null()
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_prod_depth_unsat => null()
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_prod_depth_lake => null()
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_oxid_depth_sat => null()
    !CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_oxid_depth_unsat => null()
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_oxid_depth_lake => null()
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_aere_depth_sat => null()
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_aere_depth_unsat => null()
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_tran_depth_sat => null()
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_tran_depth_unsat => null()
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_ebul_depth_sat => null()
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4_ebul_depth_unsat => null()
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_ebul_total_sat => null()
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_ebul_total_unsat => null()
    ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_aere_sat => null()
    ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_aere_unsat => null()
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_ebul_sat => null()
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_ebul_unsat => null()
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_ebul_lake => null()
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_aere_depth_sat => null()
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_aere_depth_unsat => null()
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_oxid_depth_sat => null()
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_oxid_depth_unsat => null()
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_aere_depth_sat => null()
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_aere_depth_unsat => null()
    !O2 consumption during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_decomp_depth_sat => null()
    !O2 consumption during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2_decomp_depth_unsat => null()
    ! CO2 production during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_decomp_depth_sat => null()
    ! CO2 production during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_decomp_depth_unsat => null()
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_oxid_depth_sat => null()
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: co2_oxid_depth_unsat => null()
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_o2_sat => null()
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_o2_unsat => null()
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_o2_lake => null()
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_ch4_sat => null()
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_ch4_unsat => null()
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: conc_ch4_lake => null()
    ! CH4 surface flux (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_diff_sat => null()
    ! CH4 surface flux (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_diff_unsat => null()
    ! CH4 surface flux (mol/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4_surf_diff_lake => null()
    ! CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
    real(rk8), pointer, contiguous, dimension(:) :: ch4_dfsat_flux => null()
    ! Other variables
    ! depth of water table for unsaturated fraction (m)
    real(rk8), pointer, contiguous, dimension(:) :: zwt_ch4_unsat => null()
    !fsat from previous timestep
    real(rk8), pointer, contiguous, dimension(:) :: fsat_bef => null()
    ! total soil organic matter found in level (g C / m^3) (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: lake_soilc => null()
    !aerodynamic resistance for moisture (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: lake_raw => null()
    ! total methane found in soil column (g C / m^2)
    real(rk8), pointer, contiguous, dimension(:) :: totcolch4 => null()
    ! fraction of potential heterotrophic respiration
    real(rk8), pointer, contiguous, dimension(:,:) :: fphr => null()
    ! seconds since last annual accumulator turnover
    real(rk8), pointer, contiguous, dimension(:) :: annsum_counter => null()
    ! temporary average SOM heterotrophic resp. (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: tempavg_somhr => null()
    ! annual average SOM heterotrophic resp. (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: annavg_somhr => null()
    ! respiration-weighted annual average of finundated
    real(rk8), pointer, contiguous, dimension(:) :: tempavg_finrw => null()
    ! respiration-weighted annual average of finundated
    real(rk8), pointer, contiguous, dimension(:) :: annavg_finrw => null()
    ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
    real(rk8), pointer, contiguous, dimension(:) :: sif => null()
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    !  & methanotrophs (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2stress_unsat => null()
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    ! & methanotrophs (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: o2stress_sat => null()
    ! Ratio of methane available to the total per-timestep methane
    ! sinks (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4stress_unsat => null()
    ! Ratio of methane available to the total per-timestep methane
    ! sinks (nlevsoi)
    real(rk8), pointer, contiguous, dimension(:,:) :: ch4stress_sat => null()
    ! time-lagged surface runoff (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_surf_lag => null()
    ! time-lagged fractional inundated area
    real(rk8), pointer, contiguous, dimension(:) :: finundated_lag => null()
    ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
    real(rk8), pointer, contiguous, dimension(:,:) :: layer_sat_lag => null()
  end type column_ch4_type
#endif


  !----------------------------------------------------
  ! column nitrogen state variables structure
  !----------------------------------------------------
  type, public :: column_nstate_type
    !pft-level nitrogen state variables averaged to the column
    type(pft_nstate_type) :: pns_a
    ! all n pools involved in decomposition
    ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_npools_vr => null()
    ! (gN/m3) vertically-resolved soil mineral N
    real(rk8), pointer, contiguous, dimension(:,:) :: sminn_vr => null()
    ! (gN/m3) vertically-resolved column-level sink for N truncation
    real(rk8), pointer, contiguous, dimension(:,:) :: col_ntrunc_vr => null()
#ifdef NITRIF_DENITRIF
    ! (gN/m3) vertically-resolved soil mineral NO3
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_no3_vr
    ! (gN/m2) soil mineral NO3 pool
    real(rk8), pointer, contiguous, dimension(:) :: smin_no3
    ! (gN/m3) vertically-resolved soil mineral NH4
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_nh4_vr
    ! (gN/m2) soil mineral NH4 pool
    real(rk8), pointer, contiguous, dimension(:) :: smin_nh4
#endif
    ! wood product pools, for dynamic landcover
    ! (gN/m2) column-level pool for seeding new PFTs
    real(rk8), pointer, contiguous, dimension(:) :: seedn => null()
    ! (gN/m2) wood product N pool, 10-year lifespan
    real(rk8), pointer, contiguous, dimension(:) :: prod10n => null()
    ! (gN/m2) wood product N pool, 100-year lifespan
    real(rk8), pointer, contiguous, dimension(:) :: prod100n => null()
    ! (gN/m2) total wood product N
    real(rk8), pointer, contiguous, dimension(:) :: totprodn => null()
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gN/m2)  decomposing (litter, cwd, soil) N pools
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_npools => null()
    ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_npools_1m => null()
    ! (gN/m2) soil mineral N
    real(rk8), pointer, contiguous, dimension(:) :: sminn => null()
    ! (gN/m2) column-level sink for N truncation
    real(rk8), pointer, contiguous, dimension(:) :: col_ntrunc => null()
    ! (gN/m2) Diagnostic: coarse woody debris N
    real(rk8), pointer, contiguous, dimension(:) :: cwdn => null()
    ! (gN/m2) total litter nitrogen
    real(rk8), pointer, contiguous, dimension(:) :: totlitn => null()
    ! (gN/m2) total soil organic matter nitrogen
    real(rk8), pointer, contiguous, dimension(:) :: totsomn => null()
    ! (gN/m2) total litter nitrogen to 1 meter
    real(rk8), pointer, contiguous, dimension(:) :: totlitn_1m => null()
    ! (gN/m2) total soil organic matter nitrogen to 1 meter
    real(rk8), pointer, contiguous, dimension(:) :: totsomn_1m => null()
    ! (gN/m2) total ecosystem nitrogen, incl veg
    real(rk8), pointer, contiguous, dimension(:) :: totecosysn => null()
    ! (gN/m2) total column nitrogen, incl veg
    real(rk8), pointer, contiguous, dimension(:) :: totcoln => null()
  end type column_nstate_type

  !----------------------------------------------------
  ! column VOC state variables structure
  !----------------------------------------------------
  type, public :: column_vstate_type
    !pft-level VOC state variables averaged to the column
    type(pft_vstate_type) :: pvs_a
  end type column_vstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! column DGVM state variables structure
  !----------------------------------------------------
  type, public :: column_dgvstate_type
    type(pft_dgvstate_type) :: pdgvs_a
  end type column_dgvstate_type
#endif

  !----------------------------------------------------
  ! column dust state variables structure
  !----------------------------------------------------
  type, public :: column_dstate_type
    real(rk8), pointer, contiguous, dimension(:) :: dummy_entry
  end type column_dstate_type

  !----------------------------------------------------
  ! column energy flux variables structure
  !----------------------------------------------------
  type, public :: column_eflux_type
    ! pft-level energy flux variables averaged to the column
    type(pft_eflux_type) :: pef_a
    ! snow melt heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_snomelt => null()
    ! urban snow melt heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_snomelt_u => null()
    ! rural snow melt heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_snomelt_r => null()
    ! implicit evaporation for soil temperature equation
    real(rk8), pointer, contiguous, dimension(:) :: eflx_impsoil => null()
    ! ground heat flux between soil layers 1 and 2 (W/m2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_fgr12 => null()
    ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
    real(rk8), pointer, contiguous, dimension(:,:) :: eflx_fgr => null()
    ! Urban variable
    ! heat flux from urban building interior to urban walls, roof (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_building_heat => null()
    ! urban air conditioning flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_urban_ac => null()
    ! urban heating flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_urban_heat => null()
    ! heat flux from beneath the soil or ice column (W/m**2)
    ! positive upward; usually eflx_bot >= 0
    real(rk8), pointer, contiguous, dimension(:) :: eflx_bot => null()
  end type column_eflux_type

  !----------------------------------------------------
  ! column momentum flux variables structure
  !----------------------------------------------------
  type, public :: column_mflux_type
    ! pft-level momentum flux variables averaged to the column
    type(pft_mflux_type) ::  pmf_a
  end type column_mflux_type

  !----------------------------------------------------
  ! column water flux variables structure
  !----------------------------------------------------
  type, public :: column_wflux_type
    ! pft-level water flux variables averaged to the column
    type(pft_wflux_type) :: pwf_a
    ! infiltration (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_infl => null()
    ! surface runoff (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_surf => null()
    ! sub-surface runoff (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_drain => null()
    ! net water input into soil from top (mm/s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_top_soil => null()
    ! conversion of h2osfc to ice
    real(rk8), pointer, contiguous, dimension(:) :: qflx_h2osfc_to_ice => null()
    !surface water runoff
    real(rk8), pointer, contiguous, dimension(:) :: qflx_h2osfc_surf => null()
    !snow falling on surface water
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snow_h2osfc => null()
    ! sub-surface runoff from perched wt (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_drain_perched => null()
    ! flood water flux at column level
    real(rk8), pointer, contiguous, dimension(:) :: qflx_floodc => null()
    ! liquid water + ice from layer above soil to top soil layer or sent
    ! to qflx_qrgwl (mm H2O/s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_sl_top_soil => null()
    ! snow melt (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snomelt => null()
    ! snow melt (net)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snow_melt => null()
    ! qflx_surf at glaciers, wetlands, lakes
    real(rk8), pointer, contiguous, dimension(:) :: qflx_qrgwl => null()
    ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_runoff => null()
    ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_runoff_u => null()
    ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_runoff_r => null()
    ! snow melt [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: qmelt => null()
    ! mass balance correction term for dynamic weights
    real(rk8), pointer, contiguous, dimension(:) :: h2ocan_loss => null()
    ! soil saturation excess [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_rsub_sat => null()
    ! dry (BCPHO+BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_bc_dep_dry => null()
    ! wet (BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_bc_dep_wet => null()
    ! hydrophobic BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_bc_dep_pho => null()
    ! hydrophillic BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_bc_dep_phi => null()
    ! total (dry+wet) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_bc_dep => null()
    ! dry (OCPHO+OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_oc_dep_dry => null()
    ! wet (OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_oc_dep_wet => null()
    ! hydrophobic OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_oc_dep_pho => null()
    ! hydrophillic OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_oc_dep_phi => null()
    ! total (dry+wet) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_oc_dep => null()
    ! dust species 1 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_dry1 => null()
    ! dust species 1 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_wet1 => null()
    ! dust species 2 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_dry2 => null()
    ! dust species 2 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_wet2 => null()
    ! dust species 3 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_dry3 => null()
    ! dust species 3 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_wet3 => null()
    ! dust species 4 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_dry4 => null()
    ! dust species 4 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep_wet4 => null()
    ! total (dry+wet) dust deposition on ground (positive definite) (col) [kg/s]
    real(rk8), pointer, contiguous, dimension(:) :: flx_dst_dep => null()
    ! snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(rk8), pointer, contiguous, dimension(:,:) :: qflx_snofrz_lyr => null()
    ! column-integrated snow freezing rate (positive definite) (col)
    ! [kg m-2 s-1]
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snofrz_col => null()
    !irrigation flux (mm H2O/s)
    real(rk8), pointer, contiguous, dimension(:) :: qflx_irrig => null()
  end type column_wflux_type

  !----------------------------------------------------
  ! column carbon flux variables structure
  !----------------------------------------------------
  type, public :: column_cflux_type
    ! pft-level carbon flux variables averaged to the column
    type(pft_cflux_type) :: pcf_a
    ! phenology: litterfall and crop fluxes
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! metabolic pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_c_to_litr_met_c => null()
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! cellulose pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_c_to_litr_cel_c => null()
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! lignin pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_c_to_litr_lig_c => null()
    ! gap mortality
    ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_c_to_litr_met_c => null()
    ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_c_to_litr_cel_c => null()
    ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_c_to_litr_lig_c => null()
    ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_c_to_cwdc => null()
    ! fire
    ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: fire_mortality_c_to_cwdc => null()
    ! harvest
    ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_c_to_litr_met_c => null()
    ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_c_to_litr_cel_c => null()
    ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_c_to_litr_lig_c => null()
    ! C fluxes associated with harvest to CWD pool (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_c_to_cwdc => null()
    ! new variables for CN code
    ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_to_prod10c => null()
    ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemc_to_prod100c => null()
    ! column-level fire fluxes
    ! vertically-resolved decomposing C fire loss (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: m_decomp_cpools_to_fire_vr => null()
    ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_decomp_cpools_to_fire => null()
    ! C from leaf, froot, xfer and storage C to litter labile C by fire
    ! (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_c_to_litr_met_fire => null()
    ! C from leaf, froot, xfer and storage C to litter cellulose C by fire
    ! (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_c_to_litr_cel_fire => null()
    ! C from leaf, froot, xfer and storage C to litter lignin C by fire
    ! (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_c_to_litr_lig_fire => null()
    ! (gC/m2/s) conversion C flux due to BET and BDT area decreasing
    ! (immediate loss to atm)
    real(rk8), pointer, contiguous, dimension(:) :: lf_conv_cflux => null()
    ! (gC/m2/s) carbon emissions due to peat burning
    real(rk8), pointer, contiguous, dimension(:) :: somc_fire => null()

    ! decomposition fluxes
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cascade_hr_vr => null()
    ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cascade_hr => null()
    ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cascade_ctransfer_vr => null()
    ! vertically-integrated (diagnostic) C transferred along deomposition cascade (gC/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cascade_ctransfer => null()
    ! (gC/m3/timestep)  change in decomposing c pools.
    !  Used to update concentrations concurrently with vertical transport
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cpools_sourcesink => null()
    ! rate constant for decomposition (1./sec)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_k => null()
    ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: hr_vr => null()
    ! fraction by which decomposition is limited by anoxia
    real(rk8), pointer, contiguous, dimension(:,:) :: o_scalar => null()
    ! fraction by which decomposition is limited by moisture availability
    real(rk8), pointer, contiguous, dimension(:,:) :: w_scalar => null()
    ! fraction by which decomposition is limited by temperature
    real(rk8), pointer, contiguous, dimension(:,:) :: t_scalar => null()
    ! total SOM C loss from vertical transport (gC/m^2/s)
    real(rk8), pointer, contiguous, dimension(:) :: som_c_leached => null()
    ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cpools_leached => null()
    ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cpools_transport_tendency => null()
#ifdef NITRIF_DENITRIF
    ! potential hr (not N-limited) (gC/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phr_vr
#endif
    ! dynamic landcover fluxes
#ifdef CN
    ! (gC/m2/s) seed source to PFT-level
    real(rk8), pointer, contiguous, dimension(:) :: dwt_seedc_to_leaf
    ! (gC/m2/s) seed source to PFT-level
    real(rk8), pointer, contiguous, dimension(:) :: dwt_seedc_to_deadstem
    ! (gC/m2/s) conversion C flux (immediate loss to atm)
    real(rk8), pointer, contiguous, dimension(:) :: dwt_conv_cflux
    ! (gC/m2/s) addition to 10-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: dwt_prod10c_gain
    ! (gC/m2/s) addition to 100-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: dwt_prod100c_gain
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootc_to_litr_met_c
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootc_to_litr_cel_c
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootc_to_litr_lig_c
    ! (gC/m3/s) live coarse root to CWD due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_livecrootc_to_cwdc
    ! (gC/m3/s) dead coarse root to CWD due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_deadcrootc_to_cwdc
    ! (gC/m2/s) total carbon loss from product pools and conversion
    real(rk8), pointer, contiguous, dimension(:) :: dwt_closs
    ! (gC/m2/s) dwt_closs+product_closs
    real(rk8), pointer, contiguous, dimension(:) :: landuseflux
    ! (gC/m2/s) nee-landuseflux
    real(rk8), pointer, contiguous, dimension(:) :: landuptake
    ! wood product pool loss fluxes
    ! (gC/m2/s) decomposition loss from 10-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: prod10c_loss
    ! (gC/m2/s) decomposition loss from 100-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: prod100c_loss
    ! (gC/m2/s) total wood product carbon loss
    real(rk8), pointer, contiguous, dimension(:) :: product_closs
#endif
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) litter heterotrophic respiration
    real(rk8), pointer, contiguous, dimension(:) :: lithr => null()
    ! (gC/m2/s) soil organic matter heterotrophic respiration
    real(rk8), pointer, contiguous, dimension(:) :: somhr => null()
    ! (gC/m2/s) total heterotrophic respiration
    real(rk8), pointer, contiguous, dimension(:) :: hr => null()
    ! (gC/m2/s) total soil respiration (HR + root resp)
    real(rk8), pointer, contiguous, dimension(:) :: sr => null()
    ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
    real(rk8), pointer, contiguous, dimension(:) :: er => null()
    ! (gC/m2/s) litter fire losses
    real(rk8), pointer, contiguous, dimension(:) :: litfire => null()
    ! (gC/m2/s) soil organic matter fire losses
    real(rk8), pointer, contiguous, dimension(:) :: somfire => null()
    ! (gC/m2/s) total ecosystem fire losses
    real(rk8), pointer, contiguous, dimension(:) :: totfire => null()
    ! (gC/m2/s) net ecosystem production, excludes fire, landuse,
    ! and harvest flux, positive for sink
    real(rk8), pointer, contiguous, dimension(:) :: nep => null()
    ! (gC/m2/s) net biome production, includes fire, landuse,
    ! and harvest flux, positive for sink
    real(rk8), pointer, contiguous, dimension(:) :: nbp => null()
    ! (gC/m2/s) net ecosystem exchange of carbon, includes fire,
    ! landuse, harvest, and hrv_xsmrpool flux, positive for source
    real(rk8), pointer, contiguous, dimension(:) :: nee => null()
    ! (gC/m2/s) total column-level carbon inputs (for balance check)
    real(rk8), pointer, contiguous, dimension(:) :: col_cinputs => null()
    ! (gC/m2/s) total column-level carbon outputs (for balance check)
    real(rk8), pointer, contiguous, dimension(:) :: col_coutputs => null()

#if (defined CN)
    ! CLAMP summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
    real(rk8), pointer, contiguous, dimension(:) :: cwdc_hr
    ! (gC/m2/s) col-level coarse woody debris C loss
    real(rk8), pointer, contiguous, dimension(:) :: cwdc_loss
    ! (gC/m2/s) col-level litter C loss
    real(rk8), pointer, contiguous, dimension(:) :: litterc_loss
#endif

    ! new variables for fire
    ! (gC/m2/s) total column-level fire C loss
    real(rk8), pointer, contiguous, dimension(:) :: col_fire_closs => null()
  end type column_cflux_type

  !----------------------------------------------------
  ! column nitrogen flux variables structure
  !----------------------------------------------------
  type, public :: column_nflux_type
    !pft-level nitrogen flux variables averaged to the column
    type(pft_nflux_type) :: pnf_a
    ! new variables for CN code
    ! deposition fluxes
    ! atmospheric N deposition to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ndep_to_sminn => null()
    ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: nfix_to_sminn => null()
    ! fertilizer N to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: fert_to_sminn => null()
    ! soybean fixation to soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: soyfixn_to_sminn => null()
    ! phenology: litterfall and crop fluxes
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_n_to_litr_met_n => null()
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_n_to_litr_cel_n => null()
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: phenology_n_to_litr_lig_n => null()
    ! gap mortality
    ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_n_to_litr_met_n => null()
    ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_n_to_litr_cel_n => null()
    ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_n_to_litr_lig_n => null()
    ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gap_mortality_n_to_cwdn => null()
    ! fire
    ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: fire_mortality_n_to_cwdn => null()
    ! harvest
    ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_n_to_litr_met_n => null()
    ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_n_to_litr_cel_n => null()
    ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_n_to_litr_lig_n => null()
    ! N fluxes associated with harvest to CWD pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: harvest_n_to_cwdn => null()
    !
    ! dead stem N harvest mortality to 10-year product pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_to_prod10n => null()
    ! dead stem N harvest mortality to 100-year product pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: hrv_deadstemn_to_prod100n => null()
    ! vertically-resolved decomposing N fire loss (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: m_decomp_npools_to_fire_vr => null()
    ! vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_decomp_npools_to_fire => null()
    ! column-level fire N fluxes added by F. Li and S. Levis
    ! N from leaf, froot, xfer and storage N to litter labile N
    ! by fire (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_n_to_litr_met_fire => null()
    ! N from leaf, froot, xfer and storage N to litter cellulose N
    ! by fire (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_n_to_litr_cel_fire => null()
    ! N from leaf, froot, xfer and storage N to litter lignin N
    ! by fire (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: m_n_to_litr_lig_fire => null()

    ! decomposition fluxes
    ! vert-res transfer of N from donor to receiver pool along decomp.
    ! cascade (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cascade_ntransfer_vr => null()
    ! vert-int (diagnostic) transfer of N from donor to receiver pool
    ! along decomp. cascade (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cascade_ntransfer => null()
    ! vert-res mineral N flux for transition along decomposition
    ! cascade (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_cascade_sminn_flux_vr => null()
    ! vert-int (diagnostic) mineral N flux for transition along
    ! decomposition cascade (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_cascade_sminn_flux => null()
    ! (gN/m3)  change in decomposing n pools (sum of all additions and
    ! subtractions from stateupdate1).  Used to update concentrations
    ! concurrently with vertical transport
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_npools_sourcesink => null()
    ! vertically-resolved immobilization fluxes
    ! vertically-resolved potential N immobilization (gN/m3/s) at each level
    real(rk8), pointer, contiguous, dimension(:,:) :: potential_immob_vr => null()
    ! vert-int (diagnostic) potential N immobilization (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: potential_immob => null()
    ! vertically-resolved actual N immobilization (gN/m3/s) at each level
    real(rk8), pointer, contiguous, dimension(:,:) :: actual_immob_vr => null()
    ! vert-int (diagnostic) actual N immobilization (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: actual_immob => null()
    ! vertically-resolved plant uptake of soil mineral N (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: sminn_to_plant_vr => null()
    ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: sminn_to_plant => null()
    ! vertically-resolved supplemental N supply (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: supplement_to_sminn_vr => null()
    ! vert-int (diagnostic) supplemental N supply (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: supplement_to_sminn => null()
    ! vertically-resolved gross rate of N mineralization (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: gross_nmin_vr => null()
    ! vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: gross_nmin => null()
    ! vertically-resolved net rate of N mineralization (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: net_nmin_vr => null()
    ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: net_nmin => null()
#ifdef NITRIF_DENITRIF
    ! nitrification / denitrification fluxes
    ! (gN/m3/s) soil nitrification flux
    real(rk8), pointer, contiguous, dimension(:,:) :: f_nit_vr
    ! (gN/m3/s) soil denitrification flux
    real(rk8), pointer, contiguous, dimension(:,:) :: f_denit_vr
    ! (gN/m2/s) soil nitrification flux
    real(rk8), pointer, contiguous, dimension(:) :: f_nit
    ! (gN/m2/s) soil denitrification flux
    real(rk8), pointer, contiguous, dimension(:) :: f_denit

    ! (gN/m3/s) potential soil nitrification flux
    real(rk8), pointer, contiguous, dimension(:,:) :: pot_f_nit_vr
    ! (gN/m3/s) potential soil denitrification flux
    real(rk8), pointer, contiguous, dimension(:,:) :: pot_f_denit_vr
    ! (gN/m2/s) potential soil nitrification flux
    real(rk8), pointer, contiguous, dimension(:) :: pot_f_nit
    ! (gN/m2/s) potential soil denitrification flux
    real(rk8), pointer, contiguous, dimension(:) :: pot_f_denit
    ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(rk8), pointer, contiguous, dimension(:,:) :: n2_n2o_ratio_denit_vr
    ! flux of N2o from denitrification [gN/m^3/s]
    real(rk8), pointer, contiguous, dimension(:,:) :: f_n2o_denit_vr
    ! flux of N2o from denitrification [gN/m^2/s]
    real(rk8), pointer, contiguous, dimension(:) :: f_n2o_denit
    ! flux of N2o from nitrification [gN/m^3/s]
    real(rk8), pointer, contiguous, dimension(:,:) :: f_n2o_nit_vr
    ! flux of N2o from nitrification [gN/m^2/s]
    real(rk8), pointer, contiguous, dimension(:) :: f_n2o_nit
    ! Total N2O flux from nitrification and denitrification
    real(rk8), pointer, contiguous, dimension(:,:) :: f_n2o_tot_vr
    real(rk8), pointer, contiguous, dimension(:) :: f_n2o_tot

    ! immobilization / uptake fluxes
    ! vertically-resolved actual immobilization of NO3 (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: actual_immob_no3_vr
    ! vertically-resolved actual immobilization of NH4 (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: actual_immob_nh4_vr
    ! vertically-resolved plant uptake of soil NO3 (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_no3_to_plant_vr
    ! vertically-resolved plant uptake of soil NH4 (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_nh4_to_plant_vr
    ! actual immobilization of NO3 (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: actual_immob_no3
    ! actual immobilization of NH4 (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: actual_immob_nh4
    ! plant uptake of soil NO3 (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: smin_no3_to_plant
    ! plant uptake of soil Nh4 (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: smin_nh4_to_plant
    ! leaching fluxes
    ! vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_no3_leached_vr
    ! soil mineral NO3 pool loss to leaching (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: smin_no3_leached
    ! vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_no3_runoff_vr
    ! soil mineral NO3 pool loss to runoff (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: smin_no3_runoff

    ! NITRIF_DENITRIF diagnostic quantities
    ! (ugN / g soil) soil nitrate concentration
    real(rk8), pointer, contiguous, dimension(:,:) :: smin_no3_massdens_vr
    ! (kg soil / m3) bulk density of soil
    real(rk8), pointer, contiguous, dimension(:,:) :: soil_bulkdensity
    real(rk8), pointer, contiguous, dimension(:,:) :: k_nitr_t_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: k_nitr_ph_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: k_nitr_h2o_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: k_nitr_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: wfps_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: fmax_denit_carbonsubstrate_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: fmax_denit_nitrate_vr
    ! nitrification and denitrification fluxes
    real(rk8), pointer, contiguous, dimension(:,:) :: f_denit_base_vr
    real(rk8), pointer, contiguous, dimension(:,:) :: diffus ! diffusivity (m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: ratio_k1
    real(rk8), pointer, contiguous, dimension(:,:) :: ratio_no3_co2
    real(rk8), pointer, contiguous, dimension(:,:) :: soil_co2_prod
    real(rk8), pointer, contiguous, dimension(:,:) :: fr_WFPS

    real(rk8), pointer, contiguous, dimension(:,:) :: r_psi
    real(rk8), pointer, contiguous, dimension(:,:) :: anaerobic_frac
#else
    ! denitrification fluxes
    ! vertically-resolved denitrification along decomp cascade (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: sminn_to_denit_decomp_cascade_vr
    ! vertically-integrated (diagnostic) denitrification along decomp
    ! cascade (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: sminn_to_denit_decomp_cascade
    ! vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: sminn_to_denit_excess_vr
    ! vertically-integrated (diagnostic) denitrification from excess
    ! mineral N pool (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: sminn_to_denit_excess
    ! leaching fluxes
    ! vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: sminn_leached_vr
    ! soil mineral N pool loss to leaching (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: sminn_leached
#endif
    ! dynamic landcover fluxes
    ! (gN/m2/s) seed source to PFT-level
    real(rk8), pointer, contiguous, dimension(:) :: dwt_seedn_to_leaf => null()
    ! (gN/m2/s) seed source to PFT-level
    real(rk8), pointer, contiguous, dimension(:) :: dwt_seedn_to_deadstem => null()
    ! (gN/m2/s) conversion N flux (immediate loss to atm)
    real(rk8), pointer, contiguous, dimension(:) :: dwt_conv_nflux => null()
    ! (gN/m2/s) addition to 10-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: dwt_prod10n_gain => null()
    ! (gN/m2/s) addition to 100-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: dwt_prod100n_gain => null()
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootn_to_litr_met_n => null()
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootn_to_litr_cel_n => null()
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_frootn_to_litr_lig_n => null()
    ! (gN/m3/s) live coarse root to CWD due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_livecrootn_to_cwdn => null()
    ! (gN/m3/s) dead coarse root to CWD due to landcover change
    real(rk8), pointer, contiguous, dimension(:,:) :: dwt_deadcrootn_to_cwdn => null()
    ! (gN/m2/s) total nitrogen loss from product pools and conversion
    real(rk8), pointer, contiguous, dimension(:) :: dwt_nloss => null()
    ! wood product pool loss fluxes
    ! (gN/m2/s) decomposition loss from 10-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: prod10n_loss => null()
    ! (gN/m2/s) decomposition loss from 100-yr wood product pool
    real(rk8), pointer, contiguous, dimension(:) :: prod100n_loss => null()
    ! (gN/m2/s) total wood product nitrogen loss
    real(rk8), pointer, contiguous, dimension(:) :: product_nloss => null()
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! total rate of denitrification (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: denit => null()
    ! column-level N inputs (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: col_ninputs => null()
    ! column-level N outputs (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: col_noutputs => null()
    ! new variables for fire
    ! total column-level fire N loss (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: col_fire_nloss => null()

    ! total SOM N loss from vertical transport (gN/m^2/s)
    real(rk8), pointer, contiguous, dimension(:) :: som_n_leached => null()
    ! N loss from vertical transport from each decomposing N pool (gN/m^2/s)
    real(rk8), pointer, contiguous, dimension(:,:) :: decomp_npools_leached => null()
    ! N tendency due to vertical transport in decomposing N pools (gN/m^3/s)
    real(rk8), pointer, contiguous, dimension(:,:,:) :: decomp_npools_transport_tendency => null()
  end type column_nflux_type

  !----------------------------------------------------
  ! column dust flux variables structure
  !----------------------------------------------------
  type, public :: column_dflux_type
    !pft-level dust flux variables averaged to the column
    type(pft_dflux_type) :: pdf_a
  end type column_dflux_type

  !----------------------------------------------------
  ! End definition of structures defined at the column_type level
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Begin definition of structures defined at the landunit_type level
  !----------------------------------------------------

  !----------------------------------------------------
  ! landunit physical state variables structure
  ! note - landunit type can be vegetated (includes bare soil), deep lake,
  ! shallow lake, wetland, glacier or urban
  !----------------------------------------------------
  type, public :: landunit_pstate_type
    !column-level physical state variables averaged to landunit
    type(column_pstate_type) :: cps_a
    ! Urban variables
    ! internal building temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_building => null()
    ! maximum internal building temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_building_max => null()
    ! minimum internal building temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: t_building_min => null()
    ! thermal conductivity of urban wall (W/m/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: tk_wall => null()
    ! thermal conductivity of urban roof (W/m/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: tk_roof => null()
    ! thermal conductivity of urban impervious road (W/m/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: tk_improad => null()
    ! heat capacity of urban wall (J/m^3/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: cv_wall => null()
    ! heat capacity of urban roof (J/m^3/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: cv_roof => null()
    ! heat capacity of urban impervious road (J/m^3/K)
    real(rk8), pointer, contiguous, dimension(:,:) :: cv_improad => null()
    ! total thickness of urban wall (m)
    real(rk8), pointer, contiguous, dimension(:) :: thick_wall => null()
    ! total thickness of urban roof (m)
    real(rk8), pointer, contiguous, dimension(:) :: thick_roof => null()
    ! number of impervious road layers (-)
    integer(ik4), pointer, contiguous, dimension(:) :: nlev_improad => null()
    ! view factor of sky for road
    real(rk8), pointer, contiguous, dimension(:) :: vf_sr => null()
    ! view factor of one wall for road
    real(rk8), pointer, contiguous, dimension(:) :: vf_wr => null()
    ! view factor of sky for one wall
    real(rk8), pointer, contiguous, dimension(:) :: vf_sw => null()
    ! view factor of road for one wall
    real(rk8), pointer, contiguous, dimension(:) :: vf_rw => null()
    ! view factor of opposing wall for one wall
    real(rk8), pointer, contiguous, dimension(:) :: vf_ww => null()
    ! urban canopy air temperature (K)
    real(rk8), pointer, contiguous, dimension(:) :: taf => null()
    ! urban canopy air specific humidity (kg/kg)
    real(rk8), pointer, contiguous, dimension(:) :: qaf => null()
    ! direct solar absorbed by roof per unit ground area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_roof_dir => null()
    ! diffuse solar absorbed by roof per unit ground area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_roof_dif => null()
    ! direct  solar absorbed by sunwall per unit wall area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_sunwall_dir => null()
    ! diffuse solar absorbed by sunwall per unit wall area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_sunwall_dif => null()
    ! direct  solar absorbed by shadewall per unit wall area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_shadewall_dir => null()
    ! diffuse solar absorbed by shadewall per unit wall area per unit
    ! incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_shadewall_dif => null()
    ! direct  solar absorbed by impervious road per unit ground area
    ! per unit incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_improad_dir => null()
    ! diffuse solar absorbed by impervious road per unit ground area
    ! per unit incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_improad_dif => null()
    ! direct  solar absorbed by pervious road per unit ground area
    ! per unit incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_perroad_dir => null()
    ! diffuse solar absorbed by pervious road per unit ground area
    ! per unit incident flux
    real(rk8), pointer, contiguous, dimension(:,:) :: sabs_perroad_dif => null()
  end type landunit_pstate_type

  !----------------------------------------------------
  ! landunit energy flux variables structure
  !----------------------------------------------------
  type, public :: landunit_eflux_type
    ! column-level energy flux variables averaged to landunit
    type(column_eflux_type) ::  cef_a
    ! Urban variables
    ! multiplicative traffic factor for sensible heat flux
    ! from urban traffic (-)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_traffic_factor => null()
    ! traffic sensible heat flux (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_traffic => null()
    ! sensible heat flux from domestic heating/cooling sources of waste
    ! heat (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_wasteheat => null()
    ! sensible heat flux to be put back into canyon due to removal
    ! by AC (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_heat_from_ac => null()
  end type landunit_eflux_type

  !----------------------------------------------------
  ! End definition of structures defined at the landunit_type level
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Begin definition of structures defined at the gridcell_type level
  !----------------------------------------------------

  !----------------------------------------------------
  ! gridcell physical state variables structure
  !----------------------------------------------------
  type, public :: gridcell_pstate_type
    !column-level physical state variables averaged to gridcell
    type(column_pstate_type) :: cps_a
  end type gridcell_pstate_type

  !----------------------------------------------------
  ! gridcell energy state variables structure
  !----------------------------------------------------
  type, public :: gridcell_estate_type
    !column-level energy state variables averaged to gridcell
    type(column_estate_type) :: ces_a
    ! initial gridcell total heat content
    real(rk8), pointer, contiguous, dimension(:) :: gc_heat1 => null()
    ! post land cover change total heat content
    real(rk8), pointer, contiguous, dimension(:) :: gc_heat2 => null()
  end type gridcell_estate_type

  !----------------------------------------------------
  ! gridcell water state variables structure
  !----------------------------------------------------
  type, public :: gridcell_wstate_type
    !column-level water state variables averaged to gridcell
    type(column_wstate_type) :: cws_a
    ! initial gridcell total h2o liq content
    real(rk8), pointer, contiguous, dimension(:) :: gc_liq1 => null()
    ! post land cover change total liq content
    real(rk8), pointer, contiguous, dimension(:) :: gc_liq2 => null()
    ! initial gridcell total h2o liq content
    real(rk8), pointer, contiguous, dimension(:) :: gc_ice1 => null()
    ! post land cover change total ice content
    real(rk8), pointer, contiguous, dimension(:) :: gc_ice2 => null()
  end type gridcell_wstate_type

  !----------------------------------------------------
  ! gridcell carbon state variables structure
  !----------------------------------------------------
  type, public :: gridcell_cstate_type
    !column-level carbon state variables averaged to gridcell
    type(column_cstate_type) :: ccs_a
  end type gridcell_cstate_type

#ifdef LCH4
  !----------------------------------------------------
  ! gridcell CH4 flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_ch4_type
    ! Atmospheric conc of CH4, O2, CO2 (mol/m3)
    real(rk8), pointer, contiguous, dimension(:,:) :: c_atm => null()
    ! gridcell CO2 production from CH4 oxidation (g C/m**2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4co2f => null()
    !gridcell average CH4 production (g C/m^2/s)
    real(rk8), pointer, contiguous, dimension(:) :: ch4prodg => null()
    !gridcell average net methane correction to CO2 flux (g C/m^2/s)
    real(rk8), pointer, contiguous, dimension(:) :: nem => null()
  end type gridcell_ch4_type
#endif

  !----------------------------------------------------
  ! gridcell nitrogen state variables structure
  !----------------------------------------------------
  type, public :: gridcell_nstate_type
    !column-level nitrogen state variables averaged to gridcell
    type(column_nstate_type) :: cns_a
  end type gridcell_nstate_type

  !----------------------------------------------------
  ! gridcell VOC state variables structure
  !----------------------------------------------------
  type, public :: gridcell_vstate_type
    !column-level VOC state variables averaged to gridcell
    type(column_vstate_type):: cvs_a
  end type gridcell_vstate_type

  !----------------------------------------------------
  ! gridcell VOC emission factor variables structure (heald)
  !----------------------------------------------------
  type, public :: gridcell_efstate_type
    ! isoprene emission factors
    real(rk8), pointer, contiguous, dimension(:,:) :: efisop => null()
  end type gridcell_efstate_type

  !----------------------------------------------------
  ! gridcell dust state variables structure
  !----------------------------------------------------
  type, public :: gridcell_dstate_type
    !column-level dust state variables averaged to gridcell
    type(column_dstate_type) :: cds_a
  end type gridcell_dstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! gridcell DGVM state variables structure
  !----------------------------------------------------
  type, public :: gridcell_dgvstate_type
    !20-yr running mean of agdd
    real(rk8), pointer, contiguous, dimension(:) :: agdd20 => null()
    !20-yr running mean of tmomin
    real(rk8), pointer, contiguous, dimension(:) :: tmomin20 => null()
    !ann minimum of 10-day running mean (K)
    real(rk8), pointer, contiguous, dimension(:) :: t10min => null()
  end type gridcell_dgvstate_type
#endif

  !----------------------------------------------------
  ! gridcell energy flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_eflux_type
    !column-level energy flux variables averaged to gridcell
    type(column_eflux_type) :: cef_a
    ! total grid-level sensible heat flux
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_totg => null()
    ! dynamic land cover change conversion energy flux
    real(rk8), pointer, contiguous, dimension(:) :: eflx_dynbal => null()
  end type gridcell_eflux_type

  !----------------------------------------------------
  ! gridcell momentum flux variables structure
  !-- -------------------------------------------------
  type, public :: gridcell_mflux_type
    !pft-level momentum flux variables averaged to gridcell
    type(pft_mflux_type) :: pmf_a
  end type gridcell_mflux_type

  !----------------------------------------------------
  ! gridcell water flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_wflux_type
    ! column-level water flux variables averaged to gridcell
    type(column_wflux_type) :: cwf_a
    ! total grid-level liq runoff
    real(rk8), pointer, contiguous, dimension(:) :: qflx_runoffg => null()
    ! total grid-level ice runoff
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snwcp_iceg => null()
    ! liq dynamic land cover change conversion runoff flux
    real(rk8), pointer, contiguous, dimension(:) :: qflx_liq_dynbal => null()
    ! ice dynamic land cover change conversion runoff flux
    real(rk8), pointer, contiguous, dimension(:) :: qflx_ice_dynbal => null()
    ! total grid-level flood water flux
    real(rk8), pointer, contiguous, dimension(:) :: qflx_floodg => null()
  end type gridcell_wflux_type

  !----------------------------------------------------
  ! gridcell carbon flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_cflux_type
    !column-level carbon flux variables averaged to gridcell
    type(column_cflux_type) :: ccf_a
  end type gridcell_cflux_type

  !----------------------------------------------------
  ! gridcell nitrogen flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_nflux_type
    !column-level nitrogen flux variables averaged to gridcell
    type(column_nflux_type) :: cnf_a
  end type gridcell_nflux_type

  !----------------------------------------------------
  ! gridcell dust flux variables structure
  !----------------------------------------------------
  type, public :: gridcell_dflux_type
    ! pft-level dust flux variables averaged to gridcell
    type(pft_dflux_type) :: pdf_a
  end type gridcell_dflux_type

  !----------------------------------------------------
  ! End definition of structures defined at the gridcell_type level
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Begin definition of spatial scaling hierarchy
  !----------------------------------------------------

  !----------------------------------------------------
  ! define the pft structure
  !----------------------------------------------------

  type, public :: pft_type
    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into column level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: column => null()
    !weight (relative to column)
    real(rk8), pointer, contiguous, dimension(:) :: wtcol => null()
    !index into landunit level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: landunit => null()
    !weight (relative to landunit)
    real(rk8), pointer, contiguous, dimension(:) :: wtlunit => null()
    !index into gridcell level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: gridcell => null()
    !weight (relative to gridcell)
    real(rk8), pointer, contiguous, dimension(:) :: wtgcell => null()

    ! topological mapping functionality
    !pft vegetation
    integer(ik4), pointer, contiguous, dimension(:) :: itype => null()
    !m index for laixy(i,j,m),etc.
    integer(ik4), pointer, contiguous, dimension(:) :: mxy => null()
    !true=>do computations on this pft (see reweightMod for details)
    logical, pointer, contiguous, dimension(:) :: active => null()

    ! conservation check structures for the pft level
    type(energy_balance_type)   :: pebal !energy balance structure
    type(water_balance_type)    :: pwbal !water balance structure
    type(carbon_balance_type)   :: pcbal !carbon balance structure
    type(nitrogen_balance_type) :: pnbal !nitrogen balance structure

#if (defined CNDV)
    ! DGVM state variables
    type(pft_dgvstate_type) :: pdgvs     !pft DGVM state variables
#endif

    ! CN ecophysiological variables
    type(pft_epv_type)    :: pepv        !pft ecophysiological variables

    ! state variables defined at the pft level
    type(pft_pstate_type) :: pps         !physical state variables
    type(pft_psynstate_type)::ppsyns     !photosynthesis relevant variables
    type(pft_estate_type) :: pes         !pft energy state
    type(pft_wstate_type) :: pws         !pft water state
    type(pft_cstate_type) :: pcs         !pft carbon state
    type(pft_nstate_type) :: pns         !pft nitrogen state
    type(pft_vstate_type) :: pvs         !pft VOC state

    ! flux variables defined at the pft level
    type(pft_eflux_type)  :: pef         !pft energy flux
    type(pft_mflux_type)  :: pmf         !pft momentum flux
    type(pft_wflux_type)  :: pwf         !pft water flux
    type(pft_cflux_type)  :: pcf         !pft carbon flux
    type(pft_nflux_type)  :: pnf         !pft nitrogen flux
    type(pft_vflux_type)  :: pvf         !pft VOC flux
    type(pft_dflux_type)  :: pdf         !pft dust flux
    type(pft_depvd_type)  :: pdd         !dry dep velocity

    type(pft_cstate_type) :: pc13s       !pft carbon-13 state
    type(pft_cflux_type)  :: pc13f       !pft carbon-13 flux

    type(pft_cstate_type) :: pc14s       !pft carbon-14 state
    type(pft_cflux_type)  :: pc14f       !pft carbon-14 flux
  end type pft_type

  !----------------------------------------------------
  ! define the column structure
  !----------------------------------------------------

  type, public :: column_type
    !plant functional type (pft) data structure
    type(pft_type)   :: p

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into landunit level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: landunit => null()
    !weight (relative to landunit)
    real(rk8), pointer, contiguous, dimension(:) :: wtlunit => null()
    !index into gridcell level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: gridcell => null()
    !weight (relative to gridcell)
    real(rk8), pointer, contiguous, dimension(:) :: wtgcell => null()
    !beginning pft index for each column
    integer(ik4), pointer, contiguous, dimension(:) :: pfti => null()
    !ending pft index for each column
    integer(ik4), pointer, contiguous, dimension(:) :: pftf => null()
    !number of pfts for each column
    integer(ik4), pointer, contiguous, dimension(:) :: npfts => null()

    ! topological mapping functionality
    !column type
    integer(ik4), pointer, contiguous, dimension(:) :: itype => null()
    !true=>do computations on this column (see reweightMod for details)
    logical, pointer, contiguous, dimension(:) :: active => null()

    ! conservation check structures for the column level
    type(energy_balance_type)   :: cebal !energy balance structure
    type(water_balance_type)    :: cwbal !water balance structure
    type(carbon_balance_type)   :: ccbal !carbon balance structure
    type(nitrogen_balance_type) :: cnbal !nitrogen balance structure

    ! state variables defined at the column level
    type(column_pstate_type) :: cps      !column physical state variables
    type(column_estate_type) :: ces      !column energy state
    type(column_wstate_type) :: cws      !column water state
    type(column_cstate_type) :: ccs      !column carbon state
    type(column_nstate_type) :: cns      !column nitrogen state
    type(column_dstate_type) :: cds      !column dust state

   ! flux variables defined at the column level
    type(column_eflux_type) :: cef       !column energy flux
    type(column_mflux_type) :: cmf       !column momentum flux
    type(column_wflux_type) :: cwf       !column water flux
    type(column_cflux_type) :: ccf       !column carbon flux
#ifdef LCH4
    type(column_ch4_type)   :: cch4      !column CH4 variables
#endif
    type(column_nflux_type) :: cnf       !column nitrogen flux
    type(column_dflux_type) :: cdf       !column dust flux

#if (defined CNDV)
    ! dgvm variables defined at the column level
    type (column_dgvstate_type) :: cdgvs !column DGVM structure
#endif

    type(column_cstate_type) :: cc13s    !column carbon-13 state
    type(column_cflux_type)  :: cc13f    !column carbon-13 flux

    type(column_cstate_type) :: cc14s    !column carbon-14 state
    type(column_cflux_type)  :: cc14f    !column carbon-14 flux
  end type column_type

  !----------------------------------------------------
  ! define the geomorphological land unit structure
  !----------------------------------------------------

  type, public :: landunit_type
    ! column data structure (soil/snow/canopy columns)
    type(column_type) :: c

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into gridcell level quantities
    integer(ik4), pointer, contiguous, dimension(:) :: gridcell => null()
    !weight (relative to gridcell)
    real(rk8), pointer, contiguous, dimension(:) :: wtgcell => null()
    !beginning column index per landunit
    integer(ik4), pointer, contiguous, dimension(:) :: coli => null()
    !ending column index for each landunit
    integer(ik4), pointer, contiguous, dimension(:) :: colf => null()
    !number of columns for each landunit
    integer(ik4), pointer, contiguous, dimension(:) :: ncolumns => null()
    !beginning pft index for each landunit
    integer(ik4), pointer, contiguous, dimension(:) :: pfti => null()
    !ending pft index for each landunit
    integer(ik4), pointer, contiguous, dimension(:) :: pftf => null()
    !number of pfts for each landunit
    integer(ik4), pointer, contiguous, dimension(:) :: npfts => null()

    ! Urban canyon related properties
    ! urban landunit canyon height to width ratio (-)
    real(rk8), pointer, contiguous, dimension(:) :: canyon_hwr => null()
    ! urban landunit weight of pervious road column to total road (-)
    real(rk8), pointer, contiguous, dimension(:) :: wtroad_perv => null()
    ! weight of roof with respect to urban landunit (-)
    real(rk8), pointer, contiguous, dimension(:) :: wtlunit_roof => null()

    ! Urban related info MV - this should be moved to land physical state - MV
    ! height of urban roof (m)
    real(rk8), pointer, contiguous, dimension(:) :: ht_roof => null()
    ! height above road at which wind in canyon is to be computed (m)
    real(rk8), pointer, contiguous, dimension(:) :: wind_hgt_canyon => null()
    ! urban landunit momentum roughness length (m)
    real(rk8), pointer, contiguous, dimension(:) :: z_0_town => null()
    ! urban landunit displacement height (m)
    real(rk8), pointer, contiguous, dimension(:) :: z_d_town => null()

    ! topological mapping functionality
    !landunit type
    integer(ik4), pointer, contiguous, dimension(:) :: itype => null()
    !true=>landunit is not vegetated
    logical, pointer, contiguous, dimension(:) :: ifspecial => null()
    !true=>lake point
    logical, pointer, contiguous, dimension(:) :: lakpoi => null()
    !true=>urban point
    logical, pointer, contiguous, dimension(:) :: urbpoi => null()
    !urban density type
    integer(ik4), pointer, contiguous, dimension(:) :: udenstype => null()
    !true=>do computations on this landunit (see reweightMod for details)
    logical, pointer, contiguous, dimension(:) :: active => null()

    ! state variables defined at the land unit level
    ! land unit physical state variables
    type(landunit_pstate_type) :: lps

    ! flux variables defined at the landunit level
    ! average of energy fluxes all columns
    type(landunit_eflux_type) :: lef
  end type landunit_type

  !----------------------------------------------------
  ! define the gridcell structure
  !----------------------------------------------------

  type, public :: gridcell_type
    !geomorphological landunits
    type(landunit_type) :: l

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    ! Beginning landunit index
    integer(ik4), pointer, contiguous, dimension(:) :: luni => null()
    ! Ending landunit index
    integer(ik4), pointer, contiguous, dimension(:) :: lunf => null()
    ! Number of landunit for each gridcell
    integer(ik4), pointer, contiguous, dimension(:) :: nlandunits => null()
    ! Beginning column index
    integer(ik4), pointer, contiguous, dimension(:) :: coli => null()
    ! Ending column index
    integer(ik4), pointer, contiguous, dimension(:) :: colf => null()
    ! Number of columns for each gridcell
    integer(ik4), pointer, contiguous, dimension(:) :: ncolumns => null()
    ! Beginning pft index
    integer(ik4), pointer, contiguous, dimension(:) :: pfti => null()
    ! Ending pft index
    integer(ik4), pointer, contiguous, dimension(:) :: pftf => null()
    ! Number of pfts for each gridcell
    integer(ik4), pointer, contiguous, dimension(:) :: npfts => null()

    ! Topological mapping functionality, local 1d gdc arrays

    ! Global index
    integer(ik4), pointer, contiguous, dimension(:) :: gindex => null()
    ! Total land area, gridcell (km^2)
    real(rk8), pointer, contiguous, dimension(:) :: area => null()
    ! Latitude (radians)
    real(rk8), pointer, contiguous, dimension(:) :: lat => null()
    ! Longitude (radians)
    real(rk8), pointer, contiguous, dimension(:) :: lon => null()
    ! Latitude (degrees)
    real(rk8), pointer, contiguous, dimension(:) :: latdeg => null()
    ! Longitude (degrees)
    real(rk8), pointer, contiguous, dimension(:) :: londeg => null()
    ! "atm" global index
    integer(ik4), pointer, contiguous, dimension(:) :: gindex_a => null()
    ! "atm" latitude (radians) for albedo
    real(rk8), pointer, contiguous, dimension(:) :: lat_a => null()
    ! "atm" longitude (radians) for albedo
    real(rk8), pointer, contiguous, dimension(:) :: lon_a => null()
    ! "atm" latitude (degrees) for albedo
    real(rk8), pointer, contiguous, dimension(:) :: latdeg_a => null()
    ! "atm" longitude (degrees) for albedo
    real(rk8), pointer, contiguous, dimension(:) :: londeg_a => null()

    ! Greenland ice sheet mask
    real(rk8), pointer, contiguous, dimension(:) :: gris_mask => null()
    ! Greenland ice-covered area per gridcell (km^2)
    real(rk8), pointer, contiguous, dimension(:) :: gris_area => null()
    ! Antarctic ice sheet mask
    real(rk8), pointer, contiguous, dimension(:) :: aais_mask => null()
    ! Antarctic ice-covered area per gridcell (km^2)
    real(rk8), pointer, contiguous, dimension(:) :: aais_area => null()
    ! Total water storage (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: tws => null()

#if (defined CNDV)
   ! dgvm variables defined at the gridcell level
    type(gridcell_dgvstate_type) :: gdgvs ! gridcell DGVM structure
#endif

    ! state variables defined at the gridcell level
    type(gridcell_estate_type)  :: ges ! average of energy states all landunits
    type(gridcell_wstate_type)  :: gws ! average of water states all landunits
    type(gridcell_efstate_type) :: gve ! gridcell VOC emission factors

    ! flux variables defined at the gridcell level
    type(gridcell_eflux_type) :: gef   ! average of energy fluxes all landunits
    type(gridcell_wflux_type) :: gwf   ! average of water fluxes all landunits
#ifdef LCH4
    type(gridcell_ch4_type)   :: gch4  ! average of CH4 fluxes all landunits
#endif
  end type gridcell_type

  !----------------------------------------------------
  ! define the top-level (model) structure
  !----------------------------------------------------

  type, public :: model_type
    ! lower level in hierarch
    type(gridcell_type) :: g    ! gridicell data structure
    integer(ik4) :: ngridcells  ! number of gridcells for this process
    real(rk8) :: area           ! total land area for all gridcells (km^2)
  end type model_type

  !----------------------------------------------------
  ! End definition of spatial scaling hierarchy
  !----------------------------------------------------

  !*************************************************************************

  !----------------------------------------------------
  ! Declare single instance of clmtype
  !----------------------------------------------------
  type(model_type), public, target :: clm3

  !----------------------------------------------------
  ! Declare single instance of array of ecophysiological constant types
  !----------------------------------------------------
  type(pft_epc_type), public, target :: pftcon

  !----------------------------------------------------
  ! Declare single instance of array of decomposition cascade constant types
  !----------------------------------------------------
  type(decomp_cascade_type), public, target :: decomp_cascade_con

#if (defined CNDV)
  !----------------------------------------------------
  ! Declare single instance of array of dgvm ecophysiological constant types
  !----------------------------------------------------
  type(pft_dgvepc_type), public, target :: dgv_pftcon
#endif

  ! name of gridcells
  character(len=16), parameter, public :: nameg  = 'gridcell'
  ! name of landunits
  character(len=16), parameter, public :: namel  = 'landunit'
  ! name of columns
  character(len=16), parameter, public :: namec  = 'column'
  ! name of pfts
  character(len=16), parameter, public :: namep  = 'pft'

end module mod_clm_type
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
