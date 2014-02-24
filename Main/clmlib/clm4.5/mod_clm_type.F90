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

  !----------------------------------------------------
  ! Begin definition of conservation check structures
  !----------------------------------------------------
  ! energy balance structure
  !----------------------------------------------------
  type , public :: energy_balance_type
    ! Soil/lake energy conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errsoi
    ! Surface energy conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errseb
    ! Solar radiation conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errsol
    ! Longwave radiation conservation error (W/m**2)
    real(rk8) , pointer , dimension(:) :: errlon
  end type energy_balance_type

  !----------------------------------------------------
  ! water balance structure
  !----------------------------------------------------
  type , public :: water_balance_type
    !water mass begining of the time step
    real(rk8) , pointer , dimension(:) :: begwb
    !water mass end of the time step
    real(rk8) , pointer , dimension(:) :: endwb
    !water conservation error (mm H2O)
    real(rk8) , pointer , dimension(:) :: errh2o
  end type water_balance_type

  !----------------------------------------------------
  ! carbon balance structure
  !----------------------------------------------------
  type , public :: carbon_balance_type
    !carbon mass, beginning of time step (gC/m**2)
    real(rk8) , pointer , dimension(:) :: begcb
    !carbon mass, end of time step (gC/m**2)
    real(rk8) , pointer , dimension(:) :: endcb
    !carbon balance error for the timestep (gC/m**2)
    real(rk8) , pointer , dimension(:) :: errcb
  end type carbon_balance_type

  !----------------------------------------------------
  ! nitrogen balance structure
  !----------------------------------------------------
  type , public :: nitrogen_balance_type
    !nitrogen mass, beginning of time step (gN/m**2)
    real(rk8) , pointer , dimension(:) :: begnb
    !nitrogen mass, end of time step (gN/m**2)
    real(rk8) , pointer , dimension(:) :: endnb
    !nitrogen balance error for the timestep (gN/m**2)
    real(rk8) , pointer , dimension(:) :: errnb
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
  type , public :: pft_pstate_type
    !60-day running mean of tot. precipitation (mm/s)
    ! added by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: prec60
    !10-day running mean of tot. precipitation (mm/s)
    ! added by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: prec10
    !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer , dimension(:) :: frac_veg_nosno
    !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer , dimension(:) :: frac_veg_nosno_alb
    !vegetation emissivity
    real(rk8) , pointer , dimension(:) :: emv
    !roughness length over vegetation, momentum [m]
    real(rk8) , pointer , dimension(:) :: z0mv
    !roughness length over vegetation, sensible heat [m]
    real(rk8) , pointer , dimension(:) :: z0hv
    !roughness length over vegetation, latent heat [m]
    real(rk8) , pointer , dimension(:) :: z0qv
    !fraction of roots in each soil layer  (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: rootfr
    !effective fraction of roots in each soil layer  (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: rootr
    !root resistance by layer (0-1)  (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: rresis
    !Maximum allowed dew [mm]
    real(rk8) , pointer , dimension(:) :: dewmx
    !sunlit stomatal resistance (s/m)
    real(rk8) , pointer , dimension(:) :: rssun
    !relative humidity of the canopy air vs leaf
    real(rk8) , pointer , dimension(:) :: rhal
    !vpd of the canopy air vs leaf
    real(rk8) , pointer , dimension(:) :: vpdal
    !shaded stomatal resistance (s/m)
    real(rk8) , pointer , dimension(:) :: rssha
    !canopy layer: sunlit leaf stomatal resistance (s/m)
    real(rk8) , pointer , dimension(:,:) :: rssha_z
    real(rk8) , pointer , dimension(:,:) :: rssun_z
    !canopy layer: shaded leaf stomatal resistance (s/m)
    !sunlit projected leaf area index
    real(rk8) , pointer , dimension(:) :: laisun
    !shaded projected leaf area index
    real(rk8) , pointer , dimension(:) :: laisha
    !sunlit leaf area for canopy layer
    real(rk8) , pointer , dimension(:,:) :: laisun_z
    !shaded leaf area for canopy layer
    real(rk8) , pointer , dimension(:,:) :: laisha_z
    !transpiration wetness factor (0 to 1)
    real(rk8) , pointer , dimension(:) :: btran
    ! root zone soil wetness factor (0 to 1) added by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: btran2
    !sunlit fraction of canopy
    real(rk8) , pointer , dimension(:) :: fsun
    !one-sided leaf area index, no burying by snow
    real(rk8) , pointer , dimension(:) :: tlai
    !one-sided stem area index, no burying by snow
    real(rk8) , pointer , dimension(:) :: tsai
    !one-sided leaf area index with burying by snow
    real(rk8) , pointer , dimension(:) :: elai
    !one-sided stem area index with burying by snow
    real(rk8) , pointer , dimension(:) :: esai
    !fraction of canopy that is wet (0 to 1)
    real(rk8) , pointer , dimension(:) :: fwet
    !fraction of foliage that is green and dry [-] (new)
    real(rk8) , pointer , dimension(:) :: fdry
    !change in t_veg, last iteration (Kelvin)
    real(rk8) , pointer , dimension(:) :: dt_veg
    !canopy top (m)
    real(rk8) , pointer , dimension(:) :: htop
    !canopy bottom (m)
    real(rk8) , pointer , dimension(:) :: hbot
    !momentum roughness length (m)
    real(rk8) , pointer , dimension(:) :: z0m
    !displacement height (m)
    real(rk8) , pointer , dimension(:) :: displa
    !surface albedo (direct)                              (numrad)
    real(rk8) , pointer , dimension(:,:) :: albd
    !surface albedo (diffuse)                             (numrad)
    real(rk8) , pointer , dimension(:,:) :: albi
    !flux absorbed by canopy per unit direct flux         (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabd
    !flux absorbed by sunlit canopy per unit direct flux  (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabd_sun
    !flux absorbed by shaded canopy per unit direct flux  (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabd_sha
    !flux absorbed by canopy per unit diffuse flux        (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabi
    !flux absorbed by sunlit canopy per unit diffuse flux (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabi_sun
    !flux absorbed by shaded canopy per unit diffuse flux (numrad)
    real(rk8) , pointer , dimension(:,:) :: fabi_sha
    !down direct flux below canopy per unit direct flx    (numrad)
    real(rk8) , pointer , dimension(:,:) :: ftdd
    !down diffuse flux below canopy per unit direct flx   (numrad)
    real(rk8) , pointer , dimension(:,:) :: ftid
    !down diffuse flux below canopy per unit diffuse flx  (numrad)
    real(rk8) , pointer , dimension(:,:) :: ftii
    ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsun
    ! leaf to canopy scaling coefficient, shaded leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsha
    !number of canopy layers
    integer(ik4) , pointer , dimension(:) :: ncan
    !number of canopy layers, above snow for radiative transfer
    integer(ik4) , pointer , dimension(:) :: nrad
    !absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sun_z
    !absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sha_z
    !absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sun_z
    !absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sha_z
    !sunlit fraction of canopy layer
    real(rk8) , pointer , dimension(:,:) :: fsun_z
    !tlai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tlai_z
    !tsai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tsai_z
    !10-m wind (m/s) (for dust model)
    real(rk8) , pointer , dimension(:) :: u10
    !10-m wind (m/s)
    real(rk8) , pointer , dimension(:) :: u10_clm
    !atmospheric wind speed plus convective velocity (m/s)
    real(rk8) , pointer , dimension(:) :: va
    !aerodynamical resistance (s/m)
    real(rk8) , pointer , dimension(:) :: ram1
    !aerodynamical resistance (s/m)
    real(rk8) , pointer , dimension(:) :: ram1_lake
    ! crop burn date
    integer(ik4) , pointer , dimension(:) :: burndate
    !fractional humidity at leaf surface (dimensionless)
    real(rk8) , pointer , dimension(:) :: rh_leaf
    !fractional humidity of canopy air (dimensionless)
    real(rk8) , pointer , dimension(:) :: rhaf
    !friction velocity (m/s) (for dust model)
    real(rk8) , pointer , dimension(:) :: fv
    !wind forcing height (10m+z0m+d) (m)
    real(rk8) , pointer , dimension(:) :: forc_hgt_u_pft
    !temperature forcing height (10m+z0m+d) (m)
    real(rk8) , pointer , dimension(:) :: forc_hgt_t_pft
    !specific humidity forcing height (10m+z0m+d) (m)
    real(rk8) , pointer , dimension(:) :: forc_hgt_q_pft
    ! decrease of pft weight (0-1) on the column for the timestep
    ! added by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: lfpftd
    ! Variables for prognostic crop model
    ! cold hardening index?
    real(rk8) , pointer , dimension(:) :: hdidx
    ! cumulative vernalization d?ependence?
    real(rk8) , pointer , dimension(:) :: cumvd
    ! max hgt attained by a crop during yr (m)
    real(rk8) , pointer , dimension(:) :: htmx
    ! vernalization factor for cereal
    real(rk8) , pointer , dimension(:) :: vf
    ! growing degree days (gdd) needed to harvest (ddays)
    real(rk8) , pointer , dimension(:) :: gddmaturity
    ! growing degree-days base  0C from planting  (ddays)
    real(rk8) , pointer , dimension(:) :: gdd0
    ! growing degree-days base  8C from planting  (ddays)
    real(rk8) , pointer , dimension(:) :: gdd8
    ! growing degree-days base 10C from planting  (ddays)
    real(rk8) , pointer , dimension(:) :: gdd10
    ! 20-year average of gdd0                     (ddays)
    real(rk8) , pointer , dimension(:) :: gdd020
    ! 20-year average of gdd8                     (ddays)
    real(rk8) , pointer , dimension(:) :: gdd820
    ! 20-year average of gdd10                    (ddays)
    real(rk8) , pointer , dimension(:) :: gdd1020
    ! accum gdd past planting date for crop       (ddays)
    real(rk8) , pointer , dimension(:) :: gddplant
    ! growing degree-days from planting (top two soil layers) (ddays)
    real(rk8) , pointer , dimension(:) :: gddtsoi
    ! heat unit index needed from planting to leaf emergence
    real(rk8) , pointer , dimension(:) :: huileaf
    ! heat unit index needed to reach vegetative maturity
    real(rk8) , pointer , dimension(:) :: huigrain
    ! saved leaf allocation coefficient from phase 2
    real(rk8) , pointer , dimension(:) :: aleafi
    ! saved stem allocation coefficient from phase 2
    real(rk8) , pointer , dimension(:) :: astemi
    ! leaf allocation coefficient
    real(rk8) , pointer , dimension(:) :: aleaf
    ! stem allocation coefficient
    real(rk8) , pointer , dimension(:) :: astem
    ! Flag, true if planted, not harvested
    logical , pointer , dimension(:) :: croplive
    ! Flag, true if planted
    logical , pointer , dimension(:) :: cropplant
    ! harvest date
    ! cropplant and harvdate could be 2D to facilitate rotation
    integer(ik4) , pointer , dimension(:) :: harvdate
    ! date of planting
    integer(ik4) , pointer , dimension(:) :: idop
    ! 1: max allowed lai; 0: not at max
    integer(ik4) , pointer , dimension(:) :: peaklai
    !deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
    real(rk8) , pointer , dimension(:) :: vds
    !sunlit 13c fractionation ([])
    real(rk8) , pointer , dimension(:) :: alphapsnsun
    !shaded 13c fractionation ([])
    real(rk8) , pointer , dimension(:) :: alphapsnsha
    ! sand fraction
    real(rk8) , pointer , dimension(:) :: sandfrac
    ! clay fraction
    real(rk8) , pointer , dimension(:) :: clayfrac
    ! for dry deposition of chemical tracers
    ! difference between lai month one and month two
    real(rk8) , pointer , dimension(:) :: mlaidiff
    ! aerodynamical resistance (s/m)
    real(rk8) , pointer , dimension(:) :: rb1
    ! 12 months of monthly lai from input data set
    real(rk8) , pointer , dimension(:,:) :: annlai

    ! New variable for methane code
#ifdef LCH4
    !tracer conductance for boundary layer [m/s]
    real(rk8) , pointer , dimension(:) :: grnd_ch4_cond
    !tracer conductance for canopy [m/s]
    real(rk8) , pointer , dimension(:) :: canopy_cond
#endif
    ! and vertical profiles for calculating fluxes
    ! (1/m) profile of leaves
    real(rk8) , pointer , dimension(:,:) :: leaf_prof
    ! (1/m) profile of fine roots
    real(rk8) , pointer , dimension(:,:) :: froot_prof
    ! (1/m) profile of coarse roots
    real(rk8) , pointer , dimension(:,:) :: croot_prof
    ! (1/m) profile of stems
    real(rk8) , pointer , dimension(:,:) :: stem_prof
  end type pft_pstate_type

  type , public :: pft_psynstate_type
    ! true if C3 and false if C4
    logical , pointer , dimension(:) :: c3flag
    ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: ac
    ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: aj
    ! product-limited (C3) or CO2-limited (C4) gross photosynthesis
    ! (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: ap
    ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: ag
    ! net leaf photosynthesis (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: an
    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: vcmax_z
    ! CO2 compensation point (Pa)
    real(rk8) , pointer , dimension(:) :: cp
    ! Michaelis-Menten constant for CO2 (Pa)
    real(rk8) , pointer , dimension(:) :: kc
    ! Michaelis-Menten constant for O2 (Pa)
    real(rk8) , pointer , dimension(:) :: ko
    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(rk8) , pointer , dimension(:) :: qe
    ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: tpu_z
    ! initial slope of CO2 response curve (C4 plants)
    real(rk8) , pointer , dimension(:,:) :: kp_z
    ! empirical curvature parameter for ac, aj photosynthesis co-limitation
    real(rk8) , pointer , dimension(:) :: theta_cj
    ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(rk8) , pointer , dimension(:) :: bbb
    ! Ball-Berry slope of conductance-photosynthesis relationship
    real(rk8) , pointer , dimension(:) :: mbb
    ! leaf stomatal conductance (umol H2O/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: gs_mol
    ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(rk8) , pointer , dimension(:) :: gb_mol
  end type pft_psynstate_type

  !----------------------------------------------------
  ! pft ecophysiological constants structure
  !----------------------------------------------------
  type , public :: pft_epc_type
    !value for not vegetated
    integer(ik4) , pointer , dimension(:) :: noveg
    !tree or not?
    integer(ik4) , pointer , dimension(:) :: tree
    !soil water potential at full stomatal opening (mm)
    real(rk8) , pointer , dimension(:) :: smpso
    !soil water potential at full stomatal closure (mm)
    real(rk8) , pointer , dimension(:) :: smpsc
    !foliage nitrogen limitation factor (-)
    real(rk8) , pointer , dimension(:) :: fnitr
    !foliage nitrogen (%)
    real(rk8) , pointer , dimension(:) :: foln
    !characteristic leaf dimension (m)
    real(rk8) , pointer , dimension(:) :: dleaf
    !photosynthetic pathway: 0. = c4, 1. = c3
    real(rk8) , pointer , dimension(:) :: c3psn
    !leaf/stem orientation index
    real(rk8) , pointer , dimension(:) :: xl
    !leaf reflectance: 1=vis, 2=nir   (numrad)
    real(rk8) , pointer , dimension(:,:) :: rhol
    !stem reflectance: 1=vis, 2=nir   (numrad)
    real(rk8) , pointer , dimension(:,:) :: rhos
    !leaf transmittance: 1=vis, 2=nir (numrad)
    real(rk8) , pointer , dimension(:,:) :: taul
    !stem transmittance: 1=vis, 2=nir (numrad)
    real(rk8) , pointer , dimension(:,:) :: taus
    !ratio of momentum roughness length to canopy top height (-)
    real(rk8) , pointer , dimension(:) :: z0mr
    !ratio of displacement height to canopy top height (-)
    real(rk8) , pointer , dimension(:) :: displar
    !CLM rooting distribution parameter [1/m]
    real(rk8) , pointer , dimension(:) :: roota_par
    !CLM rooting distribution parameter [1/m]
    real(rk8) , pointer , dimension(:) :: rootb_par
    ! new variables for CN code
    !wood density (gC/m3)
    real(rk8) , pointer , dimension(:) :: dwood
    !specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rk8) , pointer , dimension(:) :: slatop
    !dSLA/dLAI, projected area basis [m^2/gC]
    real(rk8) , pointer , dimension(:) :: dsladlai
    !leaf C:N (gC/gN)
    real(rk8) , pointer , dimension(:) :: leafcn
    !fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
    real(rk8) , pointer , dimension(:) :: flnr
    !binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8) , pointer , dimension(:) :: woody
    !leaf litter C:N (gC/gN)
    real(rk8) , pointer , dimension(:) :: lflitcn
    !fine root C:N (gC/gN)
    real(rk8) , pointer , dimension(:) :: frootcn
    !live wood (phloem and ray parenchyma) C:N (gC/gN)
    real(rk8) , pointer , dimension(:) :: livewdcn
    !dead wood (xylem and heartwood) C:N (gC/gN)
    real(rk8) , pointer , dimension(:) :: deadwdcn
    !grain C:N (gC/gN) for prognostic crop model
    real(rk8) , pointer , dimension(:) :: graincn
    !allocation parameter: new fine root C per new leaf C (gC/gC)
    real(rk8) , pointer , dimension(:) :: froot_leaf
    !allocation parameter: new stem c per new leaf C (gC/gC)
    real(rk8) , pointer , dimension(:) :: stem_leaf
    !allocation parameter: new coarse root C per new stem C (gC/gC)
    real(rk8) , pointer , dimension(:) :: croot_stem
    ! allocation parameter: fraction of new wood that is live
    ! (phloem and ray parenchyma) (no units)
    real(rk8) , pointer , dimension(:) :: flivewd
    ! allocation parameter: fraction of allocation that goes to currently
    ! displayed growth, remainder to storage
    real(rk8) , pointer , dimension(:) :: fcur
    !leaf litter labile fraction
    real(rk8) , pointer , dimension(:) :: lf_flab
    !leaf litter cellulose fraction
    real(rk8) , pointer , dimension(:) :: lf_fcel
    !leaf litter lignin fraction
    real(rk8) , pointer , dimension(:) :: lf_flig
    !fine root litter labile fraction
    real(rk8) , pointer , dimension(:) :: fr_flab
    !fine root litter cellulose fraction
    real(rk8) , pointer , dimension(:) :: fr_fcel
    !fine root litter lignin fraction
    real(rk8) , pointer , dimension(:) :: fr_flig
    !leaf longevity (yrs)
    real(rk8) , pointer , dimension(:) :: leaf_long
    !binary flag for evergreen leaf habit (0 or 1)
    real(rk8) , pointer , dimension(:) :: evergreen
    !binary flag for stress-deciduous leaf habit (0 or 1)
    real(rk8) , pointer , dimension(:) :: stress_decid
    !binary flag for seasonal-deciduous leaf habit (0 or 1)
    real(rk8) , pointer , dimension(:) :: season_decid

    ! fire variables added by F. Li and S. Levis
    ! combustion completeness factors (0 to 1)
    !combustion completeness factor for leaf
    real(rk8) , pointer , dimension(:) :: cc_leaf
    !combustion completeness factor for live stem
    real(rk8) , pointer , dimension(:) :: cc_lstem
    !combustion completeness factor for dead stem
    real(rk8) , pointer , dimension(:) :: cc_dstem
    !combustion completeness factor for other plant tissues
    real(rk8) , pointer , dimension(:) :: cc_other
    !  mortality factors (0 to 1)
    !fire-related mortality factor for leaf
    real(rk8) , pointer , dimension(:) :: fm_leaf
    !fire-related mortality factor for live stem
    real(rk8) , pointer , dimension(:) :: fm_lstem
    !fire-related mortality factor for dead stem
    real(rk8) , pointer , dimension(:) :: fm_dstem
    !fire-related mortality factor for other plant tissues
    real(rk8) , pointer , dimension(:) :: fm_other
    !fire-related mortality factor for fine roots
    real(rk8) , pointer , dimension(:) :: fm_root
    !fire-related mortality factor for live roots
    real(rk8) , pointer , dimension(:) :: fm_lroot
    !fire-related mortality factor for dead roots
    real(rk8) , pointer , dimension(:) :: fm_droot

    !CLM rooting distribution parameter for C and N inputs [unitless]
    real(rk8) , pointer , dimension(:) :: rootprof_beta
    ! new variables for crop code
    ! fertilizer applied
    real(rk8) , pointer , dimension(:) :: fertnitro
    ! C:N during grain fill; leaf
    real(rk8) , pointer , dimension(:) :: fleafcn
    ! C:N during grain fill; froot
    real(rk8) , pointer , dimension(:) :: ffrootcn
    ! C:N during grain fill; stem
    real(rk8) , pointer , dimension(:) :: fstemcn
  end type pft_epc_type

  type , public :: decomp_cascade_type
    !-- properties of each pathway along decomposition cascade
    ! name of transition
    character(len=8) , pointer , dimension(:) :: cascade_step_name
    ! which pool is C taken from for a given decomposition step
    integer(ik4) ,  pointer , dimension(:) :: cascade_donor_pool
    ! which pool is C added to for a given decomposition step
    integer(ik4) ,  pointer , dimension(:) :: cascade_receiver_pool
    !-- properties of each decomposing pool
    ! TRUE => pool has fixed C:N ratio
    logical ,  pointer , dimension(:) :: floating_cn_ratio_decomp_pools
    ! name of pool for restart files
    character(len=8) , pointer , dimension(:) :: decomp_pool_name_restart
    ! name of pool for history files
    character(len=8) , pointer , dimension(:) :: decomp_pool_name_history
    ! name of pool for netcdf long names
    character(len=20) , pointer , dimension(:) :: decomp_pool_name_long
    ! name of pool for netcdf short names
    character(len=8) , pointer , dimension(:) :: decomp_pool_name_short
    ! TRUE => pool is a litter pool
    logical , pointer , dimension(:) :: is_litter
    ! TRUE => pool is a soil pool
    logical , pointer , dimension(:) :: is_soil
    ! TRUE => pool is a cwd pool
    logical , pointer , dimension(:) :: is_cwd
    ! c:n ratio for initialization of pools
    real(rk8) , pointer , dimension(:) :: initial_cn_ratio
    ! initial concentration for seeding at spinup
    real(rk8) , pointer , dimension(:) :: initial_stock
    ! TRUE => pool is metabolic material
    logical , pointer , dimension(:) :: is_metabolic
    ! TRUE => pool is cellulose
    logical , pointer , dimension(:) :: is_cellulose
    ! TRUE => pool is lignin
    logical , pointer , dimension(:) :: is_lignin
    ! factor by which to scale AD and relevant processes by
    real(rk8) , pointer , dimension(:) :: spinup_factor
  end type decomp_cascade_type

#if (defined CNDV)
  !----------------------------------------------------
  ! pft DGVM-specific ecophysiological constants structure
  !----------------------------------------------------
  type , public :: pft_dgvepc_type
    !tree maximum crown area [m2]
    real(rk8) , pointer , dimension(:) :: crownarea_max
    !minimum coldest monthly mean temperature [units?]
    real(rk8) , pointer , dimension(:) :: tcmin
    !maximum coldest monthly mean temperature [units?]
    real(rk8) , pointer , dimension(:) :: tcmax
    !minimum growing degree days (at or above 5 C)
    real(rk8) , pointer , dimension(:) :: gddmin
    !upper limit of temperature of the warmest month [units?]
    real(rk8) , pointer , dimension(:) :: twmax
    !parameter in allometric equation
    real(rk8) , pointer , dimension(:) :: reinickerp
    !parameter in allometric
    real(rk8) , pointer , dimension(:) :: allom1
    !parameter in allometric
    real(rk8) , pointer , dimension(:) :: allom2
    !parameter in allometric
    real(rk8) , pointer , dimension(:) :: allom3
  end type pft_dgvepc_type
#endif

  !----------------------------------------------------
  ! pft ecophysiological variables structure
  !----------------------------------------------------
  type , public :: pft_epv_type
    !dormancy flag
    real(rk8) , pointer , dimension(:) :: dormant_flag
    !number of days since last dormancy
    real(rk8) , pointer , dimension(:) :: days_active
    !onset flag
    real(rk8) , pointer , dimension(:) :: onset_flag
    !onset days counter
    real(rk8) , pointer , dimension(:) :: onset_counter
    !onset flag for growing degree day sum
    real(rk8) , pointer , dimension(:) :: onset_gddflag
    !onset freezing degree days counter
    real(rk8) , pointer , dimension(:) :: onset_fdd
    !onset growing degree days
    real(rk8) , pointer , dimension(:) :: onset_gdd
    !onset soil water index
    real(rk8) , pointer , dimension(:) :: onset_swi
    !offset flag
    real(rk8) , pointer , dimension(:) :: offset_flag
    !offset days counter
    real(rk8) , pointer , dimension(:) :: offset_counter
    !offset freezing degree days counter
    real(rk8) , pointer , dimension(:) :: offset_fdd
    !offset soil water index
    real(rk8) , pointer , dimension(:) :: offset_swi
    !>0 fertilize; <=0 not
    real(rk8) , pointer , dimension(:) :: fert_counter
    !1: grain fill stage; 0: not
    real(rk8) , pointer , dimension(:) :: grain_flag
    !long growing season factor [0-1]
    real(rk8) , pointer , dimension(:) :: lgsf
    !background litterfall rate (1/s)
    real(rk8) , pointer , dimension(:) :: bglfr
    !background transfer growth rate (1/s)
    real(rk8) , pointer , dimension(:) :: bgtr
    !daylength (seconds)
    real(rk8) , pointer , dimension(:) :: dayl
    !daylength from previous timestep (seconds)
    real(rk8) , pointer , dimension(:) :: prev_dayl
    !annual average 2m air temperature (K)
    real(rk8) , pointer , dimension(:) :: annavg_t2m
    !temporary average 2m air temperature (K)
    real(rk8) , pointer , dimension(:) :: tempavg_t2m
    !GPP flux before downregulation (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: gpp
    !C flux available for allocation (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: availc
    !C flux assigned to recovery of negative cpool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: xsmrpool_recover
    !C13/C(12+13) ratio for xsmrpool (proportion)
    real(rk8) , pointer , dimension(:) :: xsmrpool_c13ratio
    !fraction of current allocation to display as new growth (DIM)
    real(rk8) , pointer , dimension(:) :: alloc_pnow
    !C allocation index (DIM)
    real(rk8) , pointer , dimension(:) :: c_allometry
    !N allocation index (DIM)
    real(rk8) , pointer , dimension(:) :: n_allometry
    !N flux required to support initial GPP (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: plant_ndemand
    !temporary annual sum of potential GPP
    real(rk8) , pointer , dimension(:) :: tempsum_potential_gpp
    !annual sum of potential GPP
    real(rk8) , pointer , dimension(:) :: annsum_potential_gpp
    !temporary annual max of retranslocated N pool (gN/m2)
    real(rk8) , pointer , dimension(:) :: tempmax_retransn
    !annual max of retranslocated N pool (gN/m2)
    real(rk8) , pointer , dimension(:) :: annmax_retransn
    !N flux available from retranslocation pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: avail_retransn
    !total allocated N flux (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: plant_nalloc
    !total allocated C flux (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: plant_calloc
    !C flux not allocated due to downregulation (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: excess_cflux
    !fractional reduction in GPP due to N limitation (DIM)
    real(rk8) , pointer , dimension(:) :: downreg
    !previous timestep leaf C litterfall flux (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: prev_leafc_to_litter
    !previous timestep froot C litterfall flux (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: prev_frootc_to_litter
    !temporary annual sum of NPP (gC/m2/yr)
    real(rk8) , pointer , dimension(:) :: tempsum_npp
    !annual sum of NPP (gC/m2/yr)
    real(rk8) , pointer , dimension(:) :: annsum_npp
#if (defined CNDV)
    !temporary annual sum of litfall (gC/m2/yr)
    real(rk8) , pointer , dimension(:) :: tempsum_litfall
    !annual sum of litfall (gC/m2/yr)
    real(rk8) , pointer , dimension(:) :: annsum_litfall
#endif
    !C13O2/C12O2 in canopy air
    real(rk8) , pointer , dimension(:) :: rc13_canair
    !C13O2/C12O2 in sunlit canopy psn flux
    real(rk8) , pointer , dimension(:) :: rc13_psnsun
    !C13O2/C12O2 in shaded canopy psn flux
    real(rk8) , pointer , dimension(:) :: rc13_psnsha
    !C14O2/C12O2 in atmosphere
    real(rk8) , pointer , dimension(:) :: rc14_atm
  end type pft_epv_type

  !----------------------------------------------------
  ! pft energy state variables structure
  !----------------------------------------------------
  type , public :: pft_estate_type
    !2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m
    !daily minimum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min
    !daily maximum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max
    !instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min_inst
    !instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max_inst
    !2 m height surface specific humidity (kg/kg)
    real(rk8) , pointer , dimension(:) :: q_ref2m
    !Urban 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_u
    !Rural 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_r
    !Urban daily minimum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min_u
    !Rural daily minimum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min_r
    !Urban daily maximum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max_u
    !Rural daily maximum of average 2 m height surface air temperature (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max_r
    !Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min_inst_u
    !Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_min_inst_r
    !Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max_inst_u
    !Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(rk8) , pointer , dimension(:) :: t_ref2m_max_inst_r
    ! 10-day running mean of min 2-m temperature
    real(rk8) , pointer , dimension(:) :: a10tmin
    ! 5-day running mean of min 2-m temperature
    real(rk8) , pointer , dimension(:) :: a5tmin
    !10-day running mean of the 2 m temperature (K)
    real(rk8) , pointer , dimension(:) :: t10
    !2 m height surface relative humidity (%)
    real(rk8) , pointer , dimension(:) :: rh_ref2m
    !Urban 2 m height surface relative humidity (%)
    real(rk8) , pointer , dimension(:) :: rh_ref2m_u
    !Rural 2 m height surface relative humidity (%)
    real(rk8) , pointer , dimension(:) :: rh_ref2m_r
    !vegetation temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_veg
    !intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(rk8) , pointer , dimension(:) :: thm
  end type pft_estate_type

  !----------------------------------------------------
  ! pft water state variables structure
  !----------------------------------------------------
  type , public :: pft_wstate_type
    !canopy water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2ocan
  end type pft_wstate_type

  !----------------------------------------------------
  ! pft carbon state variables structure
  !----------------------------------------------------
  type , public :: pft_cstate_type
    ! (gC/m2) ann max leaf C
    real(rk8) , pointer , dimension(:) :: leafcmax
    ! variables for prognostic crop model
    ! (gC/m2) grain C
    real(rk8) , pointer , dimension(:) :: grainc
    ! (gC/m2) grain C storage
    real(rk8) , pointer , dimension(:) :: grainc_storage
    ! (gC/m2) grain C transfer
    real(rk8) , pointer , dimension(:) :: grainc_xfer
    !
    ! (gC/m2) leaf C
    real(rk8) , pointer , dimension(:) :: leafc
    ! (gC/m2) leaf C storage
    real(rk8) , pointer , dimension(:) :: leafc_storage
    ! (gC/m2) leaf C transfer
    real(rk8) , pointer , dimension(:) :: leafc_xfer
    ! (gC/m2) fine root C
    real(rk8) , pointer , dimension(:) :: frootc
    ! (gC/m2) fine root C storage
    real(rk8) , pointer , dimension(:) :: frootc_storage
    ! (gC/m2) fine root C transfer
    real(rk8) , pointer , dimension(:) :: frootc_xfer
    ! (gC/m2) live stem C
    real(rk8) , pointer , dimension(:) :: livestemc
    ! (gC/m2) live stem C storage
    real(rk8) , pointer , dimension(:) :: livestemc_storage
    ! (gC/m2) live stem C transfer
    real(rk8) , pointer , dimension(:) :: livestemc_xfer
    ! (gC/m2) dead stem C
    real(rk8) , pointer , dimension(:) :: deadstemc
    ! (gC/m2) dead stem C storage
    real(rk8) , pointer , dimension(:) :: deadstemc_storage
    ! (gC/m2) dead stem C transfer
    real(rk8) , pointer , dimension(:) :: deadstemc_xfer
    ! (gC/m2) live coarse root C
    real(rk8) , pointer , dimension(:) :: livecrootc
    ! (gC/m2) live coarse root C storage
    real(rk8) , pointer , dimension(:) :: livecrootc_storage
    ! (gC/m2) live coarse root C transfer
    real(rk8) , pointer , dimension(:) :: livecrootc_xfer
    ! (gC/m2) dead coarse root C
    real(rk8) , pointer , dimension(:) :: deadcrootc
    ! (gC/m2) dead coarse root C storage
    real(rk8) , pointer , dimension(:) :: deadcrootc_storage
    ! (gC/m2) dead coarse root C transfer
    real(rk8) , pointer , dimension(:) :: deadcrootc_xfer
    ! (gC/m2) growth respiration storage
    real(rk8) , pointer , dimension(:) :: gresp_storage
    ! (gC/m2) growth respiration transfer
    real(rk8) , pointer , dimension(:) :: gresp_xfer
    ! (gC/m2) temporary photosynthate C pool
    real(rk8) , pointer , dimension(:) :: cpool
    ! (gC/m2) abstract C pool to meet excess MR demand
    real(rk8) , pointer , dimension(:) :: xsmrpool
    ! (gC/m2) pft-level sink for C truncation
    real(rk8) , pointer , dimension(:) :: pft_ctrunc
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gC/m2) displayed veg carbon, excluding storage and cpool
    real(rk8) , pointer , dimension(:) :: dispvegc
    ! (gC/m2) stored vegetation carbon, excluding cpool
    real(rk8) , pointer , dimension(:) :: storvegc
    ! (gC/m2) total vegetation carbon, excluding cpool
    real(rk8) , pointer , dimension(:) :: totvegc
    ! (gC/m2) total pft-level carbon, including cpool
    real(rk8) , pointer , dimension(:) :: totpftc
#if (defined CN)
    ! (gC/m2) wood C
    real(rk8) , pointer , dimension(:) :: woodc
#endif
  end type pft_cstate_type

  !----------------------------------------------------
  ! pft nitrogen state variables structure
  !----------------------------------------------------
  type , public :: pft_nstate_type
    ! variables for prognostic crop model
    ! (gN/m2) grain N
    real(rk8) , pointer , dimension(:) :: grainn
    ! (gN/m2) grain N storage
    real(rk8) , pointer , dimension(:) :: grainn_storage
    ! (gN/m2) grain N transfer
    real(rk8) , pointer , dimension(:) :: grainn_xfer
    !
    ! (gN/m2) leaf N
    real(rk8) , pointer , dimension(:) :: leafn
    ! (gN/m2) leaf N storage
    real(rk8) , pointer , dimension(:) :: leafn_storage
    ! (gN/m2) leaf N transfer
    real(rk8) , pointer , dimension(:) :: leafn_xfer
    ! (gN/m2) fine root N
    real(rk8) , pointer , dimension(:) :: frootn
    ! (gN/m2) fine root N storage
    real(rk8) , pointer , dimension(:) :: frootn_storage
    ! (gN/m2) fine root N transfer
    real(rk8) , pointer , dimension(:) :: frootn_xfer
    ! (gN/m2) live stem N
    real(rk8) , pointer , dimension(:) :: livestemn
    ! (gN/m2) live stem N storage
    real(rk8) , pointer , dimension(:) :: livestemn_storage
    ! (gN/m2) live stem N transfer
    real(rk8) , pointer , dimension(:) :: livestemn_xfer
    ! (gN/m2) dead stem N
    real(rk8) , pointer , dimension(:) :: deadstemn
    ! (gN/m2) dead stem N storage
    real(rk8) , pointer , dimension(:) :: deadstemn_storage
    ! (gN/m2) dead stem N transfer
    real(rk8) , pointer , dimension(:) :: deadstemn_xfer
    ! (gN/m2) live coarse root N
    real(rk8) , pointer , dimension(:) :: livecrootn
    ! (gN/m2) live coarse root N storage
    real(rk8) , pointer , dimension(:) :: livecrootn_storage
    ! (gN/m2) live coarse root N transfer
    real(rk8) , pointer , dimension(:) :: livecrootn_xfer
    ! (gN/m2) dead coarse root N
    real(rk8) , pointer , dimension(:) :: deadcrootn
    ! (gN/m2) dead coarse root N storage
    real(rk8) , pointer , dimension(:) :: deadcrootn_storage
    ! (gN/m2) dead coarse root N transfer
    real(rk8) , pointer , dimension(:) :: deadcrootn_xfer
    ! (gN/m2) plant pool of retranslocated N
    real(rk8) , pointer , dimension(:) :: retransn
    ! (gN/m2) temporary plant N pool
    real(rk8) , pointer , dimension(:) :: npool
    ! (gN/m2) pft-level sink for N truncation
    real(rk8) , pointer , dimension(:) :: pft_ntrunc
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gN/m2) displayed veg nitrogen, excluding storage
    real(rk8) , pointer , dimension(:) :: dispvegn
    ! (gN/m2) stored vegetation nitrogen
    real(rk8) , pointer , dimension(:) :: storvegn
    ! (gN/m2) total vegetation nitrogen
    real(rk8) , pointer , dimension(:) :: totvegn
    ! (gN/m2) total pft-level nitrogen
    real(rk8) , pointer , dimension(:) :: totpftn
  end type pft_nstate_type

  !----------------------------------------------------
  ! pft VOC state variables structure
  !----------------------------------------------------
  type , public :: pft_vstate_type
   ! 24hr average vegetation temperature (K)
   real(rk8) , pointer , dimension(:) :: t_veg24
   ! 240hr average vegetation temperature (Kelvin)
   real(rk8) , pointer , dimension(:) :: t_veg240
   ! 24hr average of direct beam radiation
   real(rk8) , pointer , dimension(:) :: fsd24
   ! 240hr average of direct beam radiation
   real(rk8) , pointer , dimension(:) :: fsd240
   ! 24hr average of diffuse beam radiation
   real(rk8) , pointer , dimension(:) :: fsi24
   ! 240hr average of diffuse beam radiation
   real(rk8) , pointer , dimension(:) :: fsi240
   ! 24hr average of sunlit fraction of canopy
   real(rk8) , pointer , dimension(:) :: fsun24
   ! 240hr average of sunlit fraction of canopy
   real(rk8) , pointer , dimension(:) :: fsun240
   ! leaf area index average over timestep
   real(rk8) , pointer , dimension(:) :: elai_p
  end type pft_vstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! pft DGVM state variables structure
  !----------------------------------------------------
  type , public :: pft_dgvstate_type
    !accumulated growing degree days above twmax
    real(rk8) , pointer , dimension(:) :: agddtw
    !accumulated growing degree days above 5
    real(rk8) , pointer , dimension(:) :: agdd
    !30-day average temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_mo
    !annual min of t_mo (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_mo_min
    !365-day running mean of tot. precipitation
    real(rk8) , pointer , dimension(:) :: prec365
    !whether PFT present in patch
    logical , pointer , dimension(:) :: present
    !if .false. then exclude seasonal decid pfts from tropics
    logical , pointer , dimension(:) :: pftmayexist
    !number of individuals (#/m**2)
    real(rk8) , pointer , dimension(:) :: nind
    !individual leaf mass
    real(rk8) , pointer , dimension(:) :: lm_ind
    !LAI per individual
    real(rk8) , pointer , dimension(:) :: lai_ind
    !foliar projective cover increment (fraction)
    real(rk8) , pointer , dimension(:) :: fpcinc
    !foliar projective cover on gridcell (fraction)
    real(rk8) , pointer , dimension(:) :: fpcgrid
    !last yr's fpcgrid
    real(rk8) , pointer , dimension(:) :: fpcgridold
    !area that each individual tree takes up (m^2)
    real(rk8) , pointer , dimension(:) :: crownarea
    real(rk8) , pointer , dimension(:) :: greffic
    real(rk8) , pointer , dimension(:) :: heatstress
  end type pft_dgvstate_type
#endif

  !----------------------------------------------------
  ! pft energy flux variables structure
  !----------------------------------------------------
  type , public :: pft_eflux_type
    !solar radiation absorbed by soil (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg_soil
    !solar radiation absorbed by snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg_snow
    !fsno weighted sum (needed by balancecheck, because fsno changes midway)
    real(rk8) , pointer , dimension(:) :: sabg_chk
    !solar radiation absorbed by ground (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabg
    !solar radiation absorbed by vegetation (W/m**2)
    real(rk8) , pointer , dimension(:) :: sabv
    !solar radiation absorbed (total) (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsa
    !urban solar radiation absorbed (total) (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsa_u
    !rural solar radiation absorbed (total) (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsa_r
    !solar radiation reflected (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr
    !absorbed PAR for sunlit leaves in canopy layer (W/m**2)
    real(rk8) , pointer , dimension(:,:) :: parsun_z
    !absorbed PAR for shaded leaves in canopy layer (W/m**2)
    real(rk8) , pointer , dimension(:,:) :: parsha_z
    !downward longwave radiation below the canopy [W/m2]
    real(rk8) , pointer , dimension(:) :: dlrad
    !upward longwave radiation above the canopy [W/m2]
    real(rk8) , pointer , dimension(:) :: ulrad
    !total latent heat flux (W/m**2)  [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_tot
    !urban total latent heat flux (W/m**2)  [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_tot_u
    !rural total latent heat flux (W/m**2)  [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_tot_r
    !ground evaporation heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_grnd
    !soil heat flux (W/m**2) [+ = into soil]
    real(rk8) , pointer , dimension(:) :: eflx_soil_grnd
    !urban soil heat flux (W/m**2) [+ = into soil]
    real(rk8) , pointer , dimension(:) :: eflx_soil_grnd_u
    !rural soil heat flux (W/m**2) [+ = into soil]
    real(rk8) , pointer , dimension(:) :: eflx_soil_grnd_r
    !total sensible heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_tot
    !urban total sensible heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_tot_u
    !rural total sensible heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_tot_r
    !sensible heat flux from ground (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_grnd
    !sensible heat flux from snow (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_snow
    !sensible heat flux from soil (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_soil
    !sensible heat flux from surface water (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_h2osfc
    !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_sh_veg
    !veg evaporation heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_vege
    !veg transpiration heat flux (W/m**2) [+ to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lh_vegt
    !sensible heat flux from domestic heating/cooling sources of waste
    ! heat (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_wasteheat_pft
    !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_heat_from_ac_pft
    !traffic sensible heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_traffic_pft
    !total anthropogenic heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_anthro
    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rk8) , pointer , dimension(:) :: cgrnd
    !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
    real(rk8) , pointer , dimension(:) :: cgrndl
    !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
    real(rk8) , pointer , dimension(:) :: cgrnds
    !net heat flux into ground (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_gnet
    ! New lake field
    !net heat flux into lake / snow surface, excluding light
    ! transmission (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_grnd_lake
    !derivative of net ground heat flux wrt soil temp (W/m**2 K)
    real(rk8) , pointer , dimension(:) :: dgnetdT
    !emitted infrared (longwave) radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_out
    !net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_net
    !urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_net_u
    !rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_net_r
    !net radiation (W/m**2) [+ = to sfc]
    real(rk8) , pointer , dimension(:) :: netrad
    !incident direct beam vis solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_vis_d
    !incident direct beam nir solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_nir_d
    !incident diffuse vis solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_vis_i
    !incident diffuse nir solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_nir_i
    !reflected direct beam vis solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_vis_d
    !reflected direct beam nir solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_nir_d
    !reflected diffuse vis solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_vis_i
    !reflected diffuse nir solar radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_nir_i
    !incident direct beam vis solar radiation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_vis_d_ln
    !incident diffuse beam vis solar radiation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_vis_i_ln
    !absorbed par by vegetation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: parveg_ln
    !incident direct beam nir solar radiation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsds_nir_d_ln
    !reflected direct beam vis solar radiation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_vis_d_ln
    !reflected direct beam nir solar radiation at local noon (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_nir_d_ln
    ! absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
    real(rk8) , pointer , dimension(:,:) :: sabg_lyr
    ! (rural) shortwave radiation penetrating top soisno layer [W/m2]
    real(rk8) , pointer , dimension(:) :: sabg_pen
    ! surface forcing of snow with all aerosols (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_aer
    ! surface forcing of snow with BC (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_bc
    ! surface forcing of snow with OC (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_oc
    ! surface forcing of snow with dust (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_dst
    ! surface forcing of snow with all aerosols, averaged only when
    ! snow is present (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_aer_sno
    ! surface forcing of snow with BC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_bc_sno
    ! surface forcing of snow with OC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_oc_sno
    ! surface forcing of snow with dust, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8) , pointer , dimension(:) :: sfc_frc_dst_sno
    ! reflected direct beam vis solar radiation from snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_sno_vd
    ! reflected direct beam NIR solar radiation from snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_sno_nd
    ! reflected diffuse vis solar radiation from snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_sno_vi
    ! reflected diffuse NIR solar radiation from snow (W/m**2)
    real(rk8) , pointer , dimension(:) :: fsr_sno_ni
    ! incident visible, direct radiation on snow (for history files)  [W/m2]
    real(rk8) , pointer , dimension(:) :: fsds_sno_vd
    ! incident near-IR, direct radiation on snow (for history files)  [W/m2]
    real(rk8) , pointer , dimension(:) :: fsds_sno_nd
    ! incident visible, diffuse radiation on snow (for history files) [W/m2]
    real(rk8) , pointer , dimension(:) :: fsds_sno_vi
    ! incident near-IR, diffuse radiation on snow (for history files) [W/m2]
    real(rk8) , pointer , dimension(:) :: fsds_sno_ni
  end type pft_eflux_type

  !----------------------------------------------------
  ! pft momentum flux variables structure
  !----------------------------------------------------
  type , public :: pft_mflux_type
    !wind (shear) stress: e-w (kg/m/s**2)
    real(rk8) , pointer , dimension(:) ::  taux
    !wind (shear) stress: n-s (kg/m/s**2)
    real(rk8) , pointer , dimension(:) ::  tauy
  end type pft_mflux_type

  !----------------------------------------------------
  ! pft water flux variables structure
  !----------------------------------------------------
  type , public :: pft_wflux_type
    !interception of precipitation [mm/s]
    real(rk8) , pointer , dimension(:) :: qflx_prec_intr
    !water onto ground including canopy runoff [kg/(m2 s)]
    real(rk8) , pointer , dimension(:) :: qflx_prec_grnd
    !rain on ground after interception (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_rain_grnd
    !snow on ground after interception (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_snow_grnd
    !excess snowfall due to snow capping (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_ice
    !excess rainfall due to snow capping (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_liq
    !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_evap_veg
    !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_tran_veg
    !evaporation from leaves and stems
    real(rk8) , pointer , dimension(:) :: qflx_evap_can
    !soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_evap_soi
    !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8) , pointer , dimension(:) :: qflx_evap_tot
    !ground surface evaporation rate (mm H2O/s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_evap_grnd
    !ground surface dew formation (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_dew_grnd
    !sublimation rate from snow pack (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_sub_snow
    !surface dew added to snow pack (mm H2O /s) [+]
    real(rk8) , pointer , dimension(:) :: qflx_dew_snow
    !snow evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_ev_snow
    !soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_ev_soil
    !h2osfc evaporation (mm H2O/s) (+ = to atm)
    real(rk8) , pointer , dimension(:) :: qflx_ev_h2osfc
  end type pft_wflux_type

  !----------------------------------------------------
  ! pft carbon flux variables structure
  !----------------------------------------------------
  type , public :: pft_cflux_type
    !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsun
    !shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsha
    !canopy layer: sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:,:) :: psnsun_z
    !canopy layer: shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:,:) :: psnsha_z
    !intracellular sunlit leaf CO2 (Pa)
    real(rk8) , pointer , dimension(:,:) :: cisun_z
    !intracellular shaded leaf CO2 (Pa)
    real(rk8) , pointer , dimension(:,:) :: cisha_z
    !sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:) :: lmrsun
    !shaded leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:) :: lmrsha
    !canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: lmrsun_z
    !canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8) , pointer , dimension(:,:) :: lmrsha_z
    !photosynthesis (umol CO2 /m**2 /s)
    real(rk8) , pointer , dimension(:) :: fpsn
    !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
    real(rk8) , pointer , dimension(:) :: fco2
    !Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsun_wc
    !Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsha_wc
    !Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8) , pointer , dimension(:) :: fpsn_wc
    !RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsun_wj
    !RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsha_wj
    !RuBP-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8) , pointer , dimension(:) :: fpsn_wj
    !product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsun_wp
    !product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8) , pointer , dimension(:) :: psnsha_wp
    !product-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8) , pointer , dimension(:) :: fpsn_wp
    ! new variables for CN code
    ! gap mortality fluxes
    ! leaf C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafc_to_litter
    ! leaf C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafc_storage_to_litter
    ! leaf C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafc_xfer_to_litter
    ! fine root C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootc_to_litter
    ! fine root C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootc_storage_to_litter
    ! fine root C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootc_xfer_to_litter
    ! live stem C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemc_to_litter
    ! live stem C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemc_storage_to_litter
    ! live stem C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemc_xfer_to_litter
    ! dead stem C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemc_to_litter
    ! dead stem C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemc_storage_to_litter
    ! dead stem C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemc_xfer_to_litter
    ! live coarse root C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootc_to_litter
    ! live coarse root C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootc_storage_to_litter
    ! live coarse root C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootc_xfer_to_litter
    ! dead coarse root C mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_to_litter
    ! dead coarse root C storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_storage_to_litter
    ! dead coarse root C transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_xfer_to_litter
    ! growth respiration storage mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_gresp_storage_to_litter
    ! growth respiration transfer mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: m_gresp_xfer_to_litter
    ! harvest mortality fluxes
    ! leaf C harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafc_to_litter
    ! leaf C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafc_storage_to_litter
    ! leaf C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafc_xfer_to_litter
    ! fine root C harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootc_to_litter
    ! fine root C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootc_storage_to_litter
    ! fine root C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootc_xfer_to_litter
    ! live stem C harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemc_to_litter
    ! live stem C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemc_storage_to_litter
    ! live stem C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemc_xfer_to_litter
    ! dead stem C harvest to 10-year product pool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_to_prod10c
    ! dead stem C harvest to 100-year product pool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_to_prod100c
    ! dead stem C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_storage_to_litter
    ! dead stem C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_xfer_to_litter
    ! live coarse root C harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootc_to_litter
    ! live coarse root C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootc_storage_to_litter
    ! live coarse root C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootc_xfer_to_litter
    ! dead coarse root C harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootc_to_litter
    ! dead coarse root C storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootc_storage_to_litter
    ! dead coarse root C transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootc_xfer_to_litter
    ! growth respiration storage harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_gresp_storage_to_litter
    ! growth respiration transfer harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_gresp_xfer_to_litter
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_xsmrpool_to_atm
    ! PFT-level fire C fluxes added by F. Li and S. Levis
    ! (gC/m2/s) fire C emissions from leafc
    real(rk8) , pointer , dimension(:) :: m_leafc_to_fire
    ! (gC/m2/s) fire C emissions from leafc_storage
    real(rk8) , pointer , dimension(:) :: m_leafc_storage_to_fire
    ! (gC/m2/s) fire C emissions from leafc_xfer
    real(rk8) , pointer , dimension(:) :: m_leafc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from livestemc
    real(rk8) , pointer , dimension(:) :: m_livestemc_to_fire
    ! (gC/m2/s) fire C emissions from livestemc_storage
    real(rk8) , pointer , dimension(:) :: m_livestemc_storage_to_fire
    ! (gC/m2/s) fire C emissions from livestemc_xfer
    real(rk8) , pointer , dimension(:) :: m_livestemc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rk8) , pointer , dimension(:) :: m_deadstemc_to_fire
    ! (gC/m2/s) fire C emissions from deadstemc_storage
    real(rk8) , pointer , dimension(:) :: m_deadstemc_storage_to_fire
    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(rk8) , pointer , dimension(:) :: m_deadstemc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from frootc
    real(rk8) , pointer , dimension(:) :: m_frootc_to_fire
    ! (gC/m2/s) fire C emissions from frootc_storage
    real(rk8) , pointer , dimension(:) :: m_frootc_storage_to_fire
    ! (gC/m2/s) fire C emissions from frootc_xfer
    real(rk8) , pointer , dimension(:) :: m_frootc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from livecrootc
    real(rk8) , pointer , dimension(:) :: m_livecrootc_to_fire
    ! (gC/m2/s) fire C emissions from livecrootc_storage
    real(rk8) , pointer , dimension(:) :: m_livecrootc_storage_to_fire
    ! (gC/m2/s) fire C emissions from livecrootc_xfer
    real(rk8) , pointer , dimension(:) :: m_livecrootc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from deadcrootc
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_to_fire
    ! (gC/m2/s) fire C emissions from deadcrootc_storage
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_storage_to_fire
    ! (gC/m2/s) fire C emissions from deadcrootc_xfer
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_xfer_to_fire
    ! (gC/m2/s) fire C emissions from gresp_storage
    real(rk8) , pointer , dimension(:) :: m_gresp_storage_to_fire
    ! (gC/m2/s) fire C emissions from gresp_xfer
    real(rk8) , pointer , dimension(:) :: m_gresp_xfer_to_fire
    ! (gC/m2/s) from leafc to litter c due to fire
    real(rk8) , pointer , dimension(:) :: m_leafc_to_litter_fire
    ! (gC/m2/s) from leafc_storage to litter C  due to fire
    real(rk8) , pointer , dimension(:) :: m_leafc_storage_to_litter_fire
    ! (gC/m2/s) from leafc_xfer to litter C  due to fire
    real(rk8) , pointer , dimension(:) :: m_leafc_xfer_to_litter_fire
    ! (gC/m2/s) from livestemc to litter C  due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemc_to_litter_fire
    ! (gC/m2/s) from livestemc_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemc_storage_to_litter_fire
    !(gC/m2/s) from livestemc_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemc_xfer_to_litter_fire
    !(gC/m2/s) from livestemc to deadstemc due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemc_to_deadstemc_fire
    !(gC/m2/s) from deadstemc to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemc_to_litter_fire
    !(gC/m2/s) from deadstemc_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemc_storage_to_litter_fire
    !(gC/m2/s) from deadstemc_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemc_xfer_to_litter_fire
    !(gC/m2/s) from frootc to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_frootc_to_litter_fire
    !(gC/m2/s) from frootc_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_frootc_storage_to_litter_fire
    !(gC/m2/s) from frootc_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_frootc_xfer_to_litter_fire
    !(gC/m2/s) from livecrootc to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootc_to_litter_fire
    !(gC/m2/s) from livecrootc_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootc_storage_to_litter_fire
    !(gC/m2/s) from livecrootc_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootc_xfer_to_litter_fire
    !(gC/m2/s) from livecrootc to deadstemc due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootc_to_deadcrootc_fire
    !(gC/m2/s) from deadcrootc to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_to_litter_fire
    !(gC/m2/s) from deadcrootc_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_storage_to_litter_fire
    !(gC/m2/s) from deadcrootc_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootc_xfer_to_litter_fire
    !(gC/m2/s) from gresp_storage to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_gresp_storage_to_litter_fire
    !(gC/m2/s) from gresp_xfer to litter C due to fire
    real(rk8) , pointer , dimension(:) :: m_gresp_xfer_to_litter_fire
    ! phenology fluxes from transfer pools
    ! grain C growth from storage for prognostic crop(gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grainc_xfer_to_grainc
    ! leaf C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leafc_xfer_to_leafc
    ! fine root C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: frootc_xfer_to_frootc
    ! live stem C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemc_xfer_to_livestemc
    ! dead stem C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: deadstemc_xfer_to_deadstemc
    ! live coarse root C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootc_xfer_to_livecrootc
    ! dead coarse root C growth from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: deadcrootc_xfer_to_deadcrootc
    ! leaf and fine root litterfall
    ! leaf C litterfall (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leafc_to_litter
    ! fine root C litterfall (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: frootc_to_litter
    ! live stem C litterfall (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemc_to_litter
    ! grain C to food for prognostic crop(gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grainc_to_food
    ! maintenance respiration fluxes
    ! leaf maintenance respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leaf_mr
    ! fine root maintenance respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: froot_mr
    ! live stem maintenance respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestem_mr
    ! live coarse root maintenance respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecroot_mr
    ! crop grain or organs maint. respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grain_mr
    ! leaf maintenance respiration from current GPP (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leaf_curmr
    ! fine root maintenance respiration from current GPP (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: froot_curmr
    ! live stem maintenance respiration from current GPP (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestem_curmr
    ! live coarse root maintenance respiration from current GPP (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecroot_curmr
    ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grain_curmr
    ! leaf maintenance respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leaf_xsmr
    ! fine root maintenance respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: froot_xsmr
    ! live stem maintenance respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestem_xsmr
    ! live coarse root maintenance respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecroot_xsmr
    ! crop grain or organs maint. respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grain_xsmr
    ! photosynthesis fluxes
    ! C fixation from sunlit canopy (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: psnsun_to_cpool
    ! C fixation from shaded canopy (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: psnshade_to_cpool
    ! allocation fluxes, from current GPP
    ! allocation to maintenance respiration storage pool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_xsmrpool
    ! allocation to grain C for prognostic crop(gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_grainc
    ! allocation to grain C storage for prognostic crop(gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_grainc_storage
    ! allocation to leaf C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_leafc
    ! allocation to leaf C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_leafc_storage
    ! allocation to fine root C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_frootc
    ! allocation to fine root C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_frootc_storage
    ! allocation to live stem C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_livestemc
    ! allocation to live stem C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_livestemc_storage
    ! allocation to dead stem C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_deadstemc
    ! allocation to dead stem C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_deadstemc_storage
    ! allocation to live coarse root C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_livecrootc
    ! allocation to live coarse root C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_livecrootc_storage
    ! allocation to dead coarse root C (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_deadcrootc
    ! allocation to dead coarse root C storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_deadcrootc_storage
    ! allocation to growth respiration storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_to_gresp_storage
    ! growth respiration fluxes
    ! excess MR pool harvest mortality (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: xsmrpool_to_atm
    ! leaf growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_leaf_gr
    ! leaf growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_leaf_storage_gr
    ! leaf growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_leaf_gr
    ! fine root growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_froot_gr
    ! fine root  growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_froot_storage_gr
    ! fine root  growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_froot_gr
    ! live stem growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_livestem_gr
    ! live stem growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_livestem_storage_gr
    ! live stem growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_livestem_gr
    ! dead stem growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_deadstem_gr
    ! dead stem growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_deadstem_storage_gr
    ! dead stem growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_deadstem_gr
    ! live coarse root growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_livecroot_gr
    ! live coarse root growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_livecroot_storage_gr
    ! live coarse root growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_livecroot_gr
    ! dead coarse root growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_deadcroot_gr
    ! dead coarse root growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_deadcroot_storage_gr
    ! dead coarse root growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_deadcroot_gr
    ! growth respiration for prognostic crop model
    ! grain growth respiration (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_grain_gr
    ! grain growth respiration to storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: cpool_grain_storage_gr
    ! grain growth respiration from storage (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: transfer_grain_gr
    ! annual turnover of storage to transfer pools
    ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: grainc_storage_to_xfer
    ! leaf C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: leafc_storage_to_xfer
    ! fine root C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: frootc_storage_to_xfer
    ! live stem C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemc_storage_to_xfer
    ! dead stem C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: deadstemc_storage_to_xfer
    ! live coarse root C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootc_storage_to_xfer
    ! dead coarse root C shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: deadcrootc_storage_to_xfer
    ! growth respiration shift storage to transfer (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: gresp_storage_to_xfer
    ! turnover of livewood to deadwood
    ! live stem C turnover (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemc_to_deadstemc
    ! live coarse root C turnover (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootc_to_deadcrootc
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) gross primary production
    real(rk8) , pointer , dimension(:) :: gpp
    ! (gC/m2/s) maintenance respiration
    real(rk8) , pointer , dimension(:) :: mr
    ! (gC/m2/s) growth resp for new growth displayed in this timestep
    real(rk8) , pointer , dimension(:) :: current_gr
    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
    real(rk8) , pointer , dimension(:) :: transfer_gr
    ! (gC/m2/s) growth resp for growth sent to storage for later display
    real(rk8) , pointer , dimension(:) :: storage_gr
    ! (gC/m2/s) total growth respiration
    real(rk8) , pointer , dimension(:) :: gr
    ! (gC/m2/s) autotrophic respiration (MR + GR)
    real(rk8) , pointer , dimension(:) :: ar
    ! (gC/m2/s) root respiration (fine root MR + total root GR)
    real(rk8) , pointer , dimension(:) :: rr
    ! (gC/m2/s) net primary production
    real(rk8) , pointer , dimension(:) :: npp
    ! (gC/m2/s) aboveground NPP
    real(rk8) , pointer , dimension(:) :: agnpp
    ! (gC/m2/s) belowground NPP
    real(rk8) , pointer , dimension(:) :: bgnpp
    ! (gC/m2/s) litterfall (leaves and fine roots)
    real(rk8) , pointer , dimension(:) :: litfall
    ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
    real(rk8) , pointer , dimension(:) :: vegfire
    ! (gC/m2/s) pft-level wood harvest (to product pools)
    real(rk8) , pointer , dimension(:) :: wood_harvestc
    ! (gC/m2/s) pft-level carbon inputs (for balance checking)
    real(rk8) , pointer , dimension(:) :: pft_cinputs
    ! (gC/m2/s) pft-level carbon outputs (for balance checking)
    real(rk8) , pointer , dimension(:) :: pft_coutputs
#if (defined CN)
    ! CLAMP summary (diagnostic) variables, not involved in mass balance
    ! (gC/m2/s) pft-level fine root C alloc
    real(rk8) , pointer , dimension(:) :: frootc_alloc
    ! (gC/m2/s) pft-level fine root C loss
    real(rk8) , pointer , dimension(:) :: frootc_loss
    ! (gC/m2/s) pft-level leaf C alloc
    real(rk8) , pointer , dimension(:) :: leafc_alloc
    ! (gC/m2/s) pft-level leaf C loss
    real(rk8) , pointer , dimension(:) :: leafc_loss
    ! (gC/m2/s) pft-level wood C alloc
    real(rk8) , pointer , dimension(:) :: woodc_alloc
    ! (gC/m2/s) pft-level wood C loss
    real(rk8) , pointer , dimension(:) :: woodc_loss
#endif
    ! new variables for fire code
    ! (gC/m2/s) total pft-level fire C loss
    real(rk8) , pointer , dimension(:) :: pft_fire_closs
! For aerenchyma calculations in CH4 code
#if (defined LCH4)
    ! (gC/m2/s) annual average aboveground NPP
    real(rk8) , pointer , dimension(:) :: annavg_agnpp
    ! (gC/m2/s) annual average belowground NPP
    real(rk8) , pointer , dimension(:) :: annavg_bgnpp
    ! (gC/m2/s) temp. average aboveground NPP
    real(rk8) , pointer , dimension(:) :: tempavg_agnpp
    ! (gC/m2/s) temp. average belowground NPP
    real(rk8) , pointer , dimension(:) :: tempavg_bgnpp
#endif
  end type pft_cflux_type

  !----------------------------------------------------
  ! pft nitrogen flux variables structure
  !----------------------------------------------------
  type , public :: pft_nflux_type
    ! new variables for CN code
    ! gap mortality fluxes
    ! leaf N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafn_to_litter
    ! fine root N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootn_to_litter
    ! leaf N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafn_storage_to_litter
    ! fine root N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootn_storage_to_litter
    ! live stem N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemn_storage_to_litter
    ! dead stem N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemn_storage_to_litter
    ! live coarse root N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootn_storage_to_litter
    ! dead coarse root N storage mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_storage_to_litter
    ! leaf N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_leafn_xfer_to_litter
    ! fine root N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_frootn_xfer_to_litter
    ! live stem N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemn_xfer_to_litter
    ! dead stem N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemn_xfer_to_litter
    ! live coarse root N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootn_xfer_to_litter
    ! dead coarse root N transfer mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_xfer_to_litter
    ! live stem N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livestemn_to_litter
    ! dead stem N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadstemn_to_litter
    ! live coarse root N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_livecrootn_to_litter
    ! dead coarse root N mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_to_litter
    ! retranslocated N pool mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: m_retransn_to_litter
    ! harvest mortality fluxes
    ! leaf N harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafn_to_litter
    ! fine root N harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootn_to_litter
    ! leaf N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafn_storage_to_litter
    ! fine root N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootn_storage_to_litter
    ! live stem N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemn_storage_to_litter
    ! dead stem N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_storage_to_litter
    ! live coarse root N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootn_storage_to_litter
    ! dead coarse root N storage harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootn_storage_to_litter
    ! leaf N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_leafn_xfer_to_litter
    ! fine root N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_frootn_xfer_to_litter
    ! live stem N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemn_xfer_to_litter
    ! dead stem N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_xfer_to_litter
    ! live coarse root N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootn_xfer_to_litter
    ! dead coarse root N transfer harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootn_xfer_to_litter
    ! live stem N harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livestemn_to_litter
    ! dead stem N harvest to 10-year product pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_to_prod10n
    ! dead stem N harvest to 100-year product pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_to_prod100n
    ! live coarse root N harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_livecrootn_to_litter
    ! dead coarse root N harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadcrootn_to_litter
    ! retranslocated N pool harvest mortality (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_retransn_to_litter
    ! PFT-level fire N fluxes added by F. Li and S. Levis
    ! (gN/m2/s) fire N emissions from leafn
    real(rk8) , pointer , dimension(:) :: m_leafn_to_fire
    ! (gN/m2/s) fire N emissions from leafn_storage
    real(rk8) , pointer , dimension(:) :: m_leafn_storage_to_fire
    ! (gN/m2/s) fire N emissions from leafn_xfer
    real(rk8) , pointer , dimension(:) :: m_leafn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from livestemn
    real(rk8) , pointer , dimension(:) :: m_livestemn_to_fire
    ! (gN/m2/s) fire N emissions from livestemn_storage
    real(rk8) , pointer , dimension(:) :: m_livestemn_storage_to_fire
    ! (gN/m2/s) fire N emissions from livestemn_xfer
    real(rk8) , pointer , dimension(:) :: m_livestemn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from deadstemn
    real(rk8) , pointer , dimension(:) :: m_deadstemn_to_fire
    ! (gN/m2/s) fire N emissions from deadstemn_storage
    real(rk8) , pointer , dimension(:) :: m_deadstemn_storage_to_fire
    ! (gN/m2/s) fire N emissions from deadstemn_xfer
    real(rk8) , pointer , dimension(:) :: m_deadstemn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from frootn
    real(rk8) , pointer , dimension(:) :: m_frootn_to_fire
    ! (gN/m2/s) fire N emissions from frootn_storage
    real(rk8) , pointer , dimension(:) :: m_frootn_storage_to_fire
    ! (gN/m2/s) fire N emissions from frootn_xfer
    real(rk8) , pointer , dimension(:) :: m_frootn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
    real(rk8) , pointer , dimension(:) :: m_livecrootn_to_fire
    ! (gN/m2/s) fire N emissions from livecrootn_storage
    real(rk8) , pointer , dimension(:) :: m_livecrootn_storage_to_fire
    ! (gN/m2/s) fire N emissions from livecrootn_xfer
    real(rk8) , pointer , dimension(:) :: m_livecrootn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from deadcrootn
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_to_fire
    ! (gN/m2/s) fire N emissions from deadcrootn_storage
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_storage_to_fire
    ! (gN/m2/s) fire N emissions from deadcrootn_xfer
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_xfer_to_fire
    ! (gN/m2/s) fire N emissions from retransn
    real(rk8) , pointer , dimension(:) :: m_retransn_to_fire
    ! (gN/m2/s) from leafn to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_leafn_to_litter_fire
    ! (gN/m2/s) from leafn_storage to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_leafn_storage_to_litter_fire
    ! (gN/m2/s) from leafn_xfer to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_leafn_xfer_to_litter_fire
    ! (gN/m2/s) from livestemn to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemn_to_litter_fire
    ! (gN/m2/s) from livestemn_storage to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemn_storage_to_litter_fire
    ! (gN/m2/s) from livestemn_xfer to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemn_xfer_to_litter_fire
    ! (gN/m2/s) from livestemn to deadstemn N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livestemn_to_deadstemn_fire
    ! (gN/m2/s) from deadstemn to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemn_to_litter_fire
    ! (gN/m2/s) from deadstemn_storage to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemn_storage_to_litter_fire
    ! (gN/m2/s) from deadstemn_xfer to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_deadstemn_xfer_to_litter_fire
    ! (gN/m2/s) from frootn to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_frootn_to_litter_fire
    ! (gN/m2/s) from frootn_storage to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_frootn_storage_to_litter_fire
    ! (gN/m2/s) from frootn_xfer to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_frootn_xfer_to_litter_fire
    ! (gN/m2/s) from livecrootn to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootn_to_litter_fire
    ! (gN/m2/s) from livecrootn_storage to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootn_storage_to_litter_fire
    ! (gN/m2/s) from livecrootn_xfer to litter N  due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootn_xfer_to_litter_fire
    ! (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire
    real(rk8) , pointer , dimension(:) :: m_livecrootn_to_deadcrootn_fire
    ! (gN/m2/s) from deadcrootn to deadcrootn due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_to_litter_fire
    ! (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_storage_to_litter_fire
    ! (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire
    real(rk8) , pointer , dimension(:) :: m_deadcrootn_xfer_to_litter_fire
    ! (gN/m2/s) from retransn to deadcrootn due to fire
    real(rk8) , pointer , dimension(:) :: m_retransn_to_litter_fire
    ! phenology fluxes from transfer pool
    ! grain N growth from storage for prognostic crop model (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: grainn_xfer_to_grainn
    ! leaf N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: leafn_xfer_to_leafn
    ! fine root N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: frootn_xfer_to_frootn
    ! live stem N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemn_xfer_to_livestemn
    ! dead stem N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: deadstemn_xfer_to_deadstemn
    ! live coarse root N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootn_xfer_to_livecrootn
    ! dead coarse root N growth from storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: deadcrootn_xfer_to_deadcrootn
    ! litterfall fluxes
    ! livestem N to litter (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemn_to_litter
    ! grain N to food for prognostic crop (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: grainn_to_food
    ! leaf N litterfall (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: leafn_to_litter
    ! leaf N to retranslocated N pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: leafn_to_retransn
    ! fine root N to retranslocated N pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: frootn_to_retransn
    ! fine root N litterfall (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: frootn_to_litter
    ! allocation fluxes
    ! deployment of retranslocated N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: retransn_to_npool
    ! deployment of soil mineral N uptake (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: sminn_to_npool
    ! allocation to grain N for prognostic crop (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_grainn
    ! allocation to grain N storage for prognostic crop (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_grainn_storage
    ! allocation to leaf N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_leafn
    ! allocation to leaf N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_leafn_storage
    ! allocation to fine root N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_frootn
    ! allocation to fine root N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_frootn_storage
    ! allocation to live stem N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_livestemn
    ! allocation to live stem N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_livestemn_storage
    ! allocation to dead stem N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_deadstemn
    ! allocation to dead stem N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_deadstemn_storage
    ! allocation to live coarse root N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_livecrootn
    ! allocation to live coarse root N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_livecrootn_storage
    ! allocation to dead coarse root N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_deadcrootn
    ! allocation to dead coarse root N storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: npool_to_deadcrootn_storage
    ! annual turnover of storage to transfer pools
    ! grain N shift storage to transfer for prognostic crop (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: grainn_storage_to_xfer
    ! leaf N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: leafn_storage_to_xfer
    ! fine root N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: frootn_storage_to_xfer
    ! live stem N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemn_storage_to_xfer
    ! dead stem N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: deadstemn_storage_to_xfer
    ! live coarse root N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootn_storage_to_xfer
    ! dead coarse root N shift storage to transfer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: deadcrootn_storage_to_xfer
    ! applied fertilizer (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: fert
    ! soybean fixed N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: soyfixn
    ! turnover of livewood to deadwood, with retranslocation
    ! live stem N turnover (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemn_to_deadstemn
    ! live stem N to retranslocated N pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livestemn_to_retransn
    ! live coarse root N turnover (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootn_to_deadcrootn
    ! live coarse root N to retranslocated N pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: livecrootn_to_retransn
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! total N deployed to growth and storage (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: ndeploy
    ! total N inputs to pft-level (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: pft_ninputs
    ! total N outputs from pft-level (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: pft_noutputs
    ! total N losses to wood product pools (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: wood_harvestn
    ! new variables for fire code
    ! total pft-level fire N loss (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: pft_fire_nloss
  end type pft_nflux_type

  !----------------------------------------------------
  ! pft VOC fluxes structure for history output
  !----------------------------------------------------
  type , public :: megan_out_type
    !(n_megan_comps) MEGAN flux [ug C m-2 h-1]
    real(rk8) , pointer , dimension(:) :: flux_out
  end type megan_out_type

  !----------------------------------------------------
  ! pft VOC flux variables structure
  !----------------------------------------------------
  type , public :: pft_vflux_type
    !total VOC flux into atmosphere [moles/m2/sec]
    real(rk8) , pointer , dimension(:) :: vocflx_tot
    !(num_mech_comps) MEGAN flux [moles/m2/sec]
    real(rk8) , pointer , dimension(:,:) :: vocflx
    real(rk8) , pointer , dimension(:) :: Eopt_out   !Eopt coefficient
    real(rk8) , pointer , dimension(:) :: topt_out   !topt coefficient
    real(rk8) , pointer , dimension(:) :: alpha_out  !alpha coefficient
    real(rk8) , pointer , dimension(:) :: cp_out     !cp coefficient
    real(rk8) , pointer , dimension(:) :: paru_out
    real(rk8) , pointer , dimension(:) :: par24u_out
    real(rk8) , pointer , dimension(:) :: par240u_out
    real(rk8) , pointer , dimension(:) :: para_out
    real(rk8) , pointer , dimension(:) :: par24a_out
    real(rk8) , pointer , dimension(:) :: par240a_out
    real(rk8) , pointer , dimension(:) :: gamma_out
    real(rk8) , pointer , dimension(:) :: gammaL_out
    real(rk8) , pointer , dimension(:) :: gammaT_out
    real(rk8) , pointer , dimension(:) :: gammaP_out
    real(rk8) , pointer , dimension(:) :: gammaA_out
    real(rk8) , pointer , dimension(:) :: gammaS_out
    real(rk8) , pointer , dimension(:) :: gammaC_out
    ! points to output fluxes
    type(megan_out_type) , pointer , dimension(:) :: meg
  end type pft_vflux_type

  !----------------------------------------------------
  ! pft dry dep velocity variables structure
  !----------------------------------------------------
  type , public :: pft_depvd_type
    real(rk8) , pointer , dimension(:,:) :: drydepvel
  end type pft_depvd_type

  !----------------------------------------------------
  ! pft dust flux variables structure
  !----------------------------------------------------
  type , public :: pft_dflux_type
    !(ndst)  !surface dust emission (kg/m**2/s) [ + = to atm]
    real(rk8) , pointer , dimension(:,:) :: flx_mss_vrt_dst
    !total dust flux into atmosphere
    real(rk8) , pointer , dimension(:) :: flx_mss_vrt_dst_tot
    !(ndst) turbulent deposition velocity (m/s)
    real(rk8) , pointer , dimension(:,:) :: vlc_trb
    !turbulent deposition velocity 1(m/s)
    real(rk8) , pointer , dimension(:) :: vlc_trb_1
    !turbulent deposition velocity 2(m/s)
    real(rk8) , pointer , dimension(:) :: vlc_trb_2
    !turbulent deposition velocity 3(m/s)
    real(rk8) , pointer , dimension(:) :: vlc_trb_3
    !turbulent deposition velocity 4(m/s)
    real(rk8) , pointer , dimension(:) :: vlc_trb_4
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
  type , public :: column_pstate_type
    !pft-level pstate variables averaged to the column
    type(pft_pstate_type) :: pps_a
    !number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    !soil color class
    integer(ik4) , pointer , dimension(:) :: isoicol
    !F. Li and S. Levis
    ! global real gdp data (k US$/capita)
    real(rk8) , pointer , dimension(:) :: gdp_lf
    ! global peatland fraction data (0-1)
    real(rk8) , pointer , dimension(:) :: peatf_lf
    ! global peak month of crop fire emissions
    integer(ik4) , pointer , dimension(:) :: abm_lf
    !gdp limitation factor for fire occurrence (0-1)
    real(rk8) , pointer , dimension(:) :: lgdp_col
    !gdp limitation factor for fire spreading (0-1)
    real(rk8) , pointer , dimension(:) :: lgdp1_col
    !pop limitation factor for fire spreading (0-1)
    real(rk8) , pointer , dimension(:) :: lpop_col
    !Clapp and Hornberger "b" (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: bsw
    !volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: watsat
    !btran parameter for btran=0
    real(rk8) , pointer , dimension(:,:) :: watdry
    !btran parameter for btran = 1
    real(rk8) , pointer , dimension(:,:) :: watopt
    !hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: hksat
    !mineral hksat
    real(rk8) , pointer , dimension(:,:) :: hksat_min
    !thermal conductivity
    real(rk8) , pointer , dimension(:,:) :: tk_hist
    !heat capacity
    real(rk8) , pointer , dimension(:,:) :: cv_hist
    !minimum soil suction (mm) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: sucsat
    !decay factor (m)
    real(rk8) , pointer , dimension(:) :: hkdepth
    !maximum saturated fraction for a gridcell
    real(rk8) , pointer , dimension(:) :: wtfact
    !fractional impermeability (-)
    real(rk8) , pointer , dimension(:,:) :: fracice
    !heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: csol
    !thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tkmg
    !thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tkdry
    !thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tksatu
    !restriction for min of soil potential (mm) (new)
    real(rk8) , pointer , dimension(:) :: smpmin
    !threshold soil moisture based on clay content
    real(rk8) , pointer , dimension(:) :: gwc_thr
    ![frc] Mass fraction clay limited to 0.20
    real(rk8) , pointer , dimension(:) :: mss_frc_cly_vld
    !basin factor
    real(rk8) , pointer , dimension(:) :: mbl_bsn_fct
    !true => do snow capping
    logical , pointer , dimension(:) :: do_capsnow
    !snow height of snow covered area (m)
    real(rk8) , pointer , dimension(:) :: snow_depth
    ! gridcell averaged snow height (m)
    real(rk8) , pointer , dimension(:) :: snowdp
    !fraction of ground covered by snow (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_sno
    !fraction of ground covered by snow (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_sno_eff
    !fractional area with surface water greater than zero
    real(rk8) , pointer , dimension(:) :: frac_h2osfc
    !temporay fractional area with surface water greater than zero
    real(rk8) , pointer , dimension(:) :: frac_h2osfc_temp
    !gridcell topographic standard deviation (m)
    real(rk8) , pointer , dimension(:) :: topo_std
    !gridcell topographic index
    real(rk8) , pointer , dimension(:) :: topo_ndx
    !gridcell topographic slope
    real(rk8) , pointer , dimension(:) :: topo_slope
    ! microtopography pdf sigma (m)
    real(rk8) , pointer , dimension(:) :: micro_sigma
    ! level at which h2osfc "percolates"
    real(rk8) , pointer , dimension(:) :: h2osfc_thresh
    ! SCA shape parameter
    real(rk8) , pointer , dimension(:) :: n_melt
    !interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: zi
    !layer thickness (m)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: dz
    !layer depth (m) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: z
    !fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: frac_iceold
    !flag for melting (=1), freezing (=2), Not=0 (new) (-nlevsno+1:nlevgrnd)
    integer(ik4) , pointer , dimension(:,:) :: imelt
    !effective porosity = porosity - vol_ice (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: eff_porosity
    !ground emissivity
    real(rk8) , pointer , dimension(:) :: emg
    !roughness length over ground, momentum [m]
    real(rk8) , pointer , dimension(:) :: z0mg
    !roughness length over ground, sensible heat [m]
    real(rk8) , pointer , dimension(:) :: z0hg
    !roughness length over ground, latent heat [m]
    real(rk8) , pointer , dimension(:) :: z0qg
    !latent heat of vapor of water (or sublimation) [j/kg]
    real(rk8) , pointer , dimension(:) :: htvp
    !coefficient of convective velocity [-]
    real(rk8) , pointer , dimension(:) :: beta
    !convective boundary height [m]
    real(rk8) , pointer , dimension(:) :: zii
    !ground albedo (direct) (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgrd
    !ground albedo (diffuse) (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgri
    !effective fraction of roots in each soil layer (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: rootr_column
    !fraction of roots in each soil layer for urban pervious road
    real(rk8) , pointer , dimension(:,:) :: rootfr_road_perv
    !effective fraction of roots in each soil layer of urban pervious road
    real(rk8) , pointer , dimension(:,:) :: rootr_road_perv
    !soil water as frac. of whc for top 0.05 m (0-1)
    ! (only comment changed by F. Li and S. Levis)
    real(rk8) , pointer , dimension(:) :: wf
    !soil water as frac. of whc for top 0.17 m (0-1)
    ! added by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: wf2
    !irrigation rate
!   real(rk8) , pointer , dimension(:) :: xirrig
    !maximum daylength for this column (s)
    real(rk8) , pointer , dimension(:) :: max_dayl
#if (defined VICHYDRO)
    !b infiltration parameter
    real(rk8) , pointer , dimension(:) :: b_infil
    !fracton of Dsmax where non-linear baseflow begins
    real(rk8) , pointer , dimension(:) :: ds
    !max. velocity of baseflow (mm/day)
    real(rk8) , pointer , dimension(:) :: dsmax
    !fraction of maximum soil moisutre where non-liear base flow occurs
    real(rk8) , pointer , dimension(:) :: wsvic
    !baseflow exponent (Qb)
    real(rk8) , pointer , dimension(:) :: c_param
    !pore-size distribution related paramter(Q12)
    real(rk8) , pointer , dimension(:,:) :: expt
    !Saturated hydrologic conductivity
    real(rk8) , pointer , dimension(:,:) :: ksat
    !soil moisture dissusion parameter
    real(rk8) , pointer , dimension(:,:) :: phi_s
    !layer depth of upper layer
    real(rk8) , pointer , dimension(:,:) :: depth
    !soil porisity (1-bulk_density/soil_density)
    real(rk8) , pointer , dimension(:,:) :: porosity
    !max layer moist + ice (mm)
    real(rk8) , pointer , dimension(:,:) :: max_moist
    !fraction of VIC layers in CLM layers
    real(rk8) , pointer , dimension(:,:,:) :: vic_clm_fract
#endif
    ! new variables for CN code
    ! solar declination angle (radians)
    real(rk8) , pointer , dimension(:) :: decl
    ! cosine of solar zenith angle
    real(rk8) , pointer , dimension(:) :: coszen
    ! soil water potential in each soil layer (MPa)
    real(rk8) , pointer , dimension(:,:) :: soilpsi
    ! bulk density of dry soil material [kg/m^3]
    real(rk8) , pointer , dimension(:,:) :: bd
    ! fraction of potential immobilization (no units)
    real(rk8) , pointer , dimension(:,:) :: fpi_vr
    ! fraction of potential immobilization (no units)
    real(rk8) , pointer , dimension(:) :: fpi
    ! respired fraction in decomposition step (frac)
    real(rk8) , pointer , dimension(:,:,:) :: rf_decomp_cascade
    ! what fraction of C leaving a given pool passes through a given
    ! transition (frac)
    real(rk8) , pointer , dimension(:,:,:) :: pathfrac_decomp_cascade
    ! (1/m) profile for N fixation additions
    real(rk8) , pointer , dimension(:,:) :: nfixation_prof
    ! (1/m) profile for N fixation additions
    real(rk8) , pointer , dimension(:,:) :: ndep_prof
    ! current depth of thaw
    real(rk8) , pointer , dimension(:) :: alt
    ! maximum annual depth of thaw
    real(rk8) , pointer , dimension(:) :: altmax
    ! prior year maximum annual depth of thaw
    real(rk8) , pointer , dimension(:) :: altmax_lastyear
    ! current depth of thaw
    integer(ik4) , pointer , dimension(:) :: alt_indx
    ! maximum annual depth of thaw
    integer(ik4) , pointer , dimension(:) :: altmax_indx
    ! prior year maximum annual depth of thaw
    integer(ik4) , pointer , dimension(:) :: altmax_lastyear_indx
    ! SOM advective flux (m/s)
    real(rk8) , pointer , dimension(:,:) :: som_adv_coef
    ! SOM diffusivity due to bio/cryo-turbation (m2/s)
    real(rk8) , pointer , dimension(:,:) :: som_diffus_coef
#ifdef NITRIF_DENITRIF
    ! maximumn monthly-mean soil temperature
    real(rk8) , pointer , dimension(:,:) :: tmean_monthly_max_vr
    ! monthly-mean soil temperature
    real(rk8) , pointer , dimension(:,:) :: tmean_monthly_vr
#endif
    !fraction of potential gpp (no units)
    real(rk8) , pointer , dimension(:) :: fpg
    !seconds since last annual accumulator turnover
    real(rk8) , pointer , dimension(:) :: annsum_counter
    !annual sum of NPP, averaged from pft-level (gC/m2/yr)
    real(rk8) , pointer , dimension(:) :: cannsum_npp
    ! (gC/m2/s) lagged net primary production
    real(rk8) , pointer , dimension(:) :: col_lag_npp
    !annual average of 2m air temperature, averaged from pft-level (K)
    real(rk8) , pointer , dimension(:) :: cannavg_t2m
    !volumetric soil water at field capacity (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: watfc
    ! F. Li and S. Levis
    ! fire counts (count/km2/timestep), valid only in Reg. C
    real(rk8) , pointer , dimension(:) :: nfire
    ! fire spread rate in pft level (m/s)
    real(rk8) , pointer , dimension(:) :: fsr_pft
    ! fire spread rate at column level (m/s)
    real(rk8) , pointer , dimension(:) :: fsr_col
    ! fire duration at column level (hr)
    real(rk8) , pointer , dimension(:) :: fd_col
    ! fire duration in pft level    (hr)
    real(rk8) , pointer , dimension(:) :: fd_pft
    !60-day running mean of tot. precipitation (mm/s)
    real(rk8) , pointer , dimension(:) :: prec60_col
    !10-day running mean of tot. precipitation (mm/s)
    real(rk8) , pointer , dimension(:) :: prec10_col
    ! conversion area fraction of BET and BDT that haven't burned before (0-1)
    real(rk8) , pointer , dimension(:) :: lfc
    ! conversion area fraction of BET and BDT that burned in this
    ! timestep ((timestep)-1)
    real(rk8) , pointer , dimension(:) :: lfc2
    ! annual decreased fraction coverage of BET on the gridcell (0-1)
    real(rk8) , pointer , dimension(:) :: dtrotr_col
    ! pft weight of BET and BDT on the gridcell(0-1)
    real(rk8) , pointer , dimension(:) :: trotr1_col
    ! pft weight of BDT on the gridcell (0-1)
    real(rk8) , pointer , dimension(:) :: trotr2_col
    ! crop fraction in veg column (0-1)
    real(rk8) , pointer , dimension(:) :: cropf_col
    ! baf for cropland per time step(0-1)
    real(rk8) , pointer , dimension(:) :: baf_crop
    ! baf for peatland per time step (0-1)
    real(rk8) , pointer , dimension(:) :: baf_peatf
    ! total burned area out of conversion (0-1)
    real(rk8) , pointer , dimension(:) :: fbac
    ! burned area out of conversion region due to land use fire (0-1)
    real(rk8) , pointer , dimension(:) :: fbac1
    ! btran2 at column level (0-1)
    real(rk8) , pointer , dimension(:) :: btran_col
    ! fractional coverage of non-crop PFTs (0-1)
    real(rk8) , pointer , dimension(:) :: wtlf
    ! fractional coverage of non-crop and non-bare-soil PFTs (0-1)
    real(rk8) , pointer , dimension(:) :: lfwt
    !timestep fractional area burned (0-1)
    real(rk8) , pointer , dimension(:) :: farea_burned
    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsnd_hst
    ! snow albedo, diffuse, for history files (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsni_hst
    ! soil albedo: direct (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsod
    ! soil albedo: diffuse (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsoi
    ! absorbed flux per unit incident direct flux: VIS (col,lyr) [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absdv
    ! absorbed flux per unit incident direct flux: NIR (col,lyr) [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absdn
    ! absorbed flux per unit incident diffuse flux: VIS (col,lyr) [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absiv
    ! absorbed flux per unit incident diffuse flux: NIR (col,lyr) [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absin
    ! snow grain radius (col,lyr) [m^-6, microns]
    real(rk8) , pointer , dimension(:,:) :: snw_rds
    ! snow grain radius, top layer (col) [m^-6, microns]
    real(rk8) , pointer , dimension(:) :: snw_rds_top
    ! snow liquid water fraction (mass), top layer (col) [fraction]
    real(rk8) , pointer , dimension(:) :: sno_liq_top
    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcpho
    ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcphi
    ! total mass of BC in snow (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bctot
    ! column-integrated mass of total BC (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_col
    ! top-layer mass of total BC (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_top
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocpho
    ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocphi
    ! total mass of OC in snow (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_octot
    ! column-integrated mass of total OC (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_col
    ! top-layer mass of total OC (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_top
    ! mass of dust species 1 in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_dst1
    ! mass of dust species 2 in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_dst2
    ! mass of dust species 3 in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_dst3
    ! mass of dust species 4 in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_dst4
    ! total mass of dust in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_dsttot
    ! column-integrated mass of dust in snow (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_col
    ! top-layer mass of dust in snow (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_top
    ! top-layer mass of snow (col) [kg]
    real(rk8) , pointer , dimension(:) :: h2osno_top
    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcphi
    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcpho
    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocphi
    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocpho
    ! mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst1
    ! mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst2
    ! mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst3
    ! mass concentration of dust species 4 in snow (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst4
    ! pure snow ground direct albedo (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgrd_pur
    ! pure snow ground diffuse albedo (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgri_pur
    ! ground direct albedo without BC  (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgrd_bc
    ! ground diffuse albedo without BC (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgri_bc
    ! ground direct albedo without OC  (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgrd_oc
    ! ground diffuse albedo without OC (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgri_oc
    ! ground direct albedo without dust  (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgrd_dst
    ! ground diffuse albedo without dust (numrad)
    real(rk8) , pointer , dimension(:,:) :: albgri_dst
    ! temperature gradient in top layer  [K m-1]
    real(rk8) , pointer , dimension(:) :: dTdz_top
    ! temperature of top snow layer [K]
    real(rk8) , pointer , dimension(:) :: snot_top
    ! new variables for S Lake code
    ! surface friction velocity (m/s)
    real(rk8) , pointer , dimension(:) :: ws
    ! coefficient for calculation of decay of eddy diffusivity with depth
    real(rk8) , pointer , dimension(:) :: ks
    ! lake layer thickness (m)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: dz_lake
    ! layer depth for lake (m)
    real(rk8) , pointer , dimension(:,:) :: z_lake
    ! top level eddy conductivity from previous timestep (W/mK)
    real(rk8) , pointer , dimension(:) :: savedtke1
    ! sand value for gridcell containing column (1:nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: cellsand
    ! clay value for gridcell containing column (1:nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: cellclay
    ! organic matter for gridcell containing column (1:nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: cellorg
    ! variable lake depth (m)
    real(rk8) , pointer , dimension(:) :: lakedepth
    ! lake extinction coefficient from surface data (1/m)
    real(rk8) , pointer , dimension(:) :: etal
    ! lake fetch from surface data (m)
    real(rk8) , pointer , dimension(:) :: lakefetch
    ! friction velocity (m/s)
    real(rk8) , pointer , dimension(:) :: ust_lake
    ! end new variables for S Lake code
    ! New variables for finundated in methane code
#ifdef LCH4
    ! coefficient for determining finundated (m)
    real(rk8) , pointer , dimension(:) :: zwt0
    ! maximum inundated fraction for a gridcell (for methane code)
    real(rk8) , pointer , dimension(:) :: f0
    ! coefficient for determining finundated (m)
    real(rk8) , pointer , dimension(:) :: p3
    ! added by Lei Meng for pH effects of methane production
    ! pH values
    real(rk8) , pointer , dimension(:) :: pH
#endif
    ! End New variables for methane code
    ! current irrigation rate [mm/s]
    real(rk8) , pointer , dimension(:) :: irrig_rate
    ! number of time steps for which we still need to irrigate today
    ! (if 0, ignore irrig_rate)
    integer(ik4) , pointer  , dimension(:) :: n_irrig_steps_left
    ! surface atm pressure, downscaled to column (Pa)
    real(rk8) , pointer , dimension(:) :: forc_pbot
    ! surface air density, downscaled to column (kg/m^3)
    real(rk8) , pointer , dimension(:) :: forc_rho
    ! ice fractional area
    real(rk8) , pointer , dimension(:) :: glc_frac
    ! surface elevation (m)
    real(rk8) , pointer , dimension(:) :: glc_topo
  end type column_pstate_type

  !----------------------------------------------------
  ! column energy state variables structure
  !----------------------------------------------------
  type , public :: column_estate_type
    !pft-level energy state variables averaged to the column
    type(pft_estate_type) :: pes_a
    !ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    !Urban ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd_u
    !Rural ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd_r
    !change in t_grnd, last iteration (Kelvin)
    real(rk8) , pointer , dimension(:) :: dt_grnd
    !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    !soil temperature in top 10cm of soil (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_soi_10cm
    !soil temperature in top 17cm of soil (Kelvin) by F. Li and S. Levis
    real(rk8) , pointer , dimension(:) :: tsoi17
    !lake temperature (Kelvin)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: t_lake
    !soil/snow temperature before update (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: tssbef
    !virtual potential temperature (kelvin)
    real(rk8) , pointer , dimension(:) :: thv
    !soil heat content (MJ/m2)
    real(rk8) , pointer , dimension(:) :: hc_soi
    !soil plus snow heat content (MJ/m2)
    real(rk8) , pointer , dimension(:) :: hc_soisno
    !atm temperature, downscaled to column (Kelvin)
    real(rk8) , pointer , dimension(:) :: forc_t
    !atm potl temperature, downscaled to column (Kelvin)
    real(rk8) , pointer , dimension(:) :: forc_th
    !surface water temperature
    real(rk8) , pointer , dimension(:) :: t_h2osfc
    !surface water temperature from time-step before
    real(rk8) , pointer , dimension(:) :: t_h2osfc_bef
  end type column_estate_type

  !----------------------------------------------------
  ! column water state variables structure
  !----------------------------------------------------
  type , public :: column_wstate_type
    !pft-level water state variables averaged to the column
    type(pft_wstate_type) :: pws_a
    !surface water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osfc
    !ground specific humidity [kg/kg]
    real(rk8) , pointer , dimension(:) :: qg_snow
    !ground specific humidity [kg/kg]
    real(rk8) , pointer , dimension(:) :: qg_soil
    !ground specific humidity [kg/kg]
    real(rk8) , pointer , dimension(:) :: qg_h2osfc
    !initial snow water
    real(rk8) , pointer , dimension(:,:) :: swe_old
    !snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    !imbalance in snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: errh2osno
    !snow sources (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: snow_sources
    !snow sinks (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: snow_sinks
    !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    !liquid water + ice lens in top 10cm of soil (kg/m2)
    real(rk8) , pointer , dimension(:) :: h2osoi_liqice_10cm
    !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    !snow mass for previous time step (kg/m2) (new)
    real(rk8) , pointer , dimension(:) :: h2osno_old
    !ground specific humidity [kg/kg]
    real(rk8) , pointer , dimension(:) :: qg
    !d(qg)/dT
    real(rk8) , pointer , dimension(:) :: dqgdT
    !average snow ice lens
    real(rk8) , pointer , dimension(:) :: snowice
    !average snow liquid water
    real(rk8) , pointer , dimension(:) :: snowliq
    !factor that reduces ground saturated specific humidity (-)
    real(rk8) , pointer , dimension(:) :: soilalpha
    !factor that reduces ground evaporation L&P1992(-)
    real(rk8) , pointer , dimension(:) :: soilbeta
    !urban factor that reduces ground saturated specific humidity (-)
    real(rk8) , pointer , dimension(:) :: soilalpha_u
    !water table depth
    real(rk8) , pointer , dimension(:) :: zwt
    !frost table depth
    real(rk8) , pointer , dimension(:) :: frost_table
    !perched water table depth
    real(rk8) , pointer , dimension(:) :: zwt_perched
    !integrated snowfall (mm H2O)
    real(rk8) , pointer , dimension(:) :: int_snow
    !fractional impermeable area
    real(rk8) , pointer , dimension(:) :: fcov
    !water in the unconfined aquifer (mm)
    real(rk8) , pointer , dimension(:) :: wa
    !aquifer recharge rate (mm/s)
    real(rk8) , pointer , dimension(:) :: qcharge
    !soil matric potential (mm)
    real(rk8) , pointer , dimension(:,:) :: smp_l
    !hydraulic conductivity (mm/s)
    real(rk8) , pointer , dimension(:,:) :: hk_l
    !fractional area with water table at surface
    real(rk8) , pointer , dimension(:) :: fsat
    !atm specific humidity, downscaled to column (kg/kg)
    real(rk8) , pointer , dimension(:) :: forc_q
#if (defined VICHYDRO)
    !soil moisture (kg/m2) for VIC soil layers
    real(rk8) , pointer , dimension(:,:) :: moist
    !soil ice (kg/m2) for VIC soil layers
    real(rk8) , pointer , dimension(:,:) :: ice
    !volumetric soil moisture for VIC soil layers
    real(rk8) , pointer , dimension(:,:) :: moist_vol
    !maximum infiltration rate calculated by VIC
    real(rk8) , pointer , dimension(:) :: max_infil
    !average saturation in top soil layers in VIC
    real(rk8) , pointer , dimension(:) :: i_0
#endif
#ifdef LCH4
    !fractional inundated area (excluding dedicated wetland columns)
    real(rk8) , pointer , dimension(:) :: finundated
#endif
    ! new variables for S Lake code
    ! mass fraction of lake layer that is frozen
    real(rk8) , pointer , dimension(:,:) :: lake_icefrac
    ! ice thickness (m) (integrated if lakepuddling)
    real(rk8) , pointer , dimension(:) :: lake_icethick
    ! end new variables for S Lake code
  end type column_wstate_type

  !----------------------------------------------------
  ! column carbon state variables structure
  !----------------------------------------------------
  type , public :: column_cstate_type
    !pft-level carbon state variables averaged to the column
    type(pft_cstate_type) :: pcs_a
    ! NOTE: the soilc variable is used by the original CLM C-cycle code,
    ! and is not used by the CN code
    !soil carbon (kg C /m**2)
    real(rk8) , pointer , dimension(:) :: soilc
    ! all c pools involved in decomposition
    ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cpools_vr
    ! (gC/m3) vertically-resolved column-level sink for C truncation
    real(rk8) , pointer , dimension(:,:) :: col_ctrunc_vr
    !fire-related variables added by F. Li and S. Levis
    !root carbon at column level (gC/m2)
    real(rk8) , pointer , dimension(:) :: rootc_col
    !column-level totvegc (gC/m2)
    real(rk8) , pointer , dimension(:) :: totvegc_col
    !column-level leafc (gC/m2)
    real(rk8) , pointer , dimension(:) :: leafc_col
    ! fuel avalability factor for Reg.C (0-1)
    real(rk8) , pointer , dimension(:) :: fuelc
    ! fuel avalability factor for Reg.A (0-1)
    real(rk8) , pointer , dimension(:) :: fuelc_crop
    ! pools for dynamic landcover
    ! (gC/m2) column-level pool for seeding new PFTs
    real(rk8) , pointer , dimension(:) :: seedc
    ! (gC/m2) wood product C pool, 10-year lifespan
    real(rk8) , pointer , dimension(:) :: prod10c
    ! (gC/m2) wood product C pool, 100-year lifespan
    real(rk8) , pointer , dimension(:) :: prod100c
    ! (gC/m2) total wood product C
    real(rk8) , pointer , dimension(:) :: totprodc
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(rk8) , pointer , dimension(:,:) :: decomp_cpools
    ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(rk8) , pointer , dimension(:,:) :: decomp_cpools_1m
    ! (gC/m2) Diagnostic: coarse woody debris C
    real(rk8) , pointer , dimension(:) :: cwdc
    ! (gC/m2) column-level sink for C truncation
    real(rk8) , pointer , dimension(:) :: col_ctrunc
    ! (gC/m2) total litter carbon
    real(rk8) , pointer , dimension(:) :: totlitc
    ! (gC/m2) total soil organic matter carbon
    real(rk8) , pointer , dimension(:) :: totsomc
    ! (gC/m2) total litter carbon to 1 meter
    real(rk8) , pointer , dimension(:) :: totlitc_1m
    ! (gC/m2) total soil organic matter carbon to 1 meter
    real(rk8) , pointer , dimension(:) :: totsomc_1m
    ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
    real(rk8) , pointer , dimension(:) :: totecosysc
    ! (gC/m2) total column carbon, incl veg and cpool
    real(rk8) , pointer , dimension(:) :: totcolc
  end type column_cstate_type

#ifdef LCH4
  !----------------------------------------------------
  ! column methane variables structure
  !----------------------------------------------------
  type , public :: column_ch4_type
    ! new variables for CH4 code
    ! column-level methane fluxes
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_sat
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_unsat
    ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_prod_depth_lake
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_sat
    !CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_unsat
    ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_oxid_depth_lake
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_aere_depth_sat
    ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_aere_depth_unsat
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_tran_depth_sat
    ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_tran_depth_unsat
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_ebul_depth_sat
    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: ch4_ebul_depth_unsat
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_ebul_total_sat
    ! Total column CH4 ebullition (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_ebul_total_unsat
    ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_aere_sat
    ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_aere_unsat
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_sat
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_unsat
    ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_ebul_lake
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_aere_depth_sat
    ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_aere_depth_unsat
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_oxid_depth_sat
    ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_oxid_depth_unsat
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_aere_depth_sat
    ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: o2_aere_depth_unsat
    !O2 consumption during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: o2_decomp_depth_sat
    !O2 consumption during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: o2_decomp_depth_unsat
    ! CO2 production during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: co2_decomp_depth_sat
    ! CO2 production during decomposition in each soil layer (nlevsoi)(mol/m3/s)
    real(rk8) , pointer , dimension(:,:) :: co2_decomp_depth_unsat
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_oxid_depth_sat
    ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: co2_oxid_depth_unsat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_sat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_unsat
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_o2_lake
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_sat
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_unsat
    ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: conc_ch4_lake
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_sat
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_unsat
    ! CH4 surface flux (mol/m2/s)
    real(rk8) , pointer , dimension(:) :: ch4_surf_diff_lake
    ! CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
    real(rk8) , pointer , dimension(:) :: ch4_dfsat_flux
    ! Other variables
    ! depth of water table for unsaturated fraction (m)
    real(rk8) , pointer , dimension(:) :: zwt_ch4_unsat
    !fsat from previous timestep
    real(rk8) , pointer , dimension(:) :: fsat_bef
    ! total soil organic matter found in level (g C / m^3) (nlevsoi)
    real(rk8) , pointer , dimension(:,:) :: lake_soilc
    !aerodynamic resistance for moisture (s/m)
    real(rk8) , pointer , dimension(:) :: lake_raw
    ! total methane found in soil column (g C / m^2)
    real(rk8) , pointer , dimension(:) :: totcolch4
    ! fraction of potential heterotrophic respiration
    real(rk8) , pointer , dimension(:,:) :: fphr
    ! seconds since last annual accumulator turnover
    real(rk8) , pointer , dimension(:) :: annsum_counter
    ! temporary average SOM heterotrophic resp. (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: tempavg_somhr
    ! annual average SOM heterotrophic resp. (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: annavg_somhr
    ! respiration-weighted annual average of finundated
    real(rk8) , pointer , dimension(:) :: tempavg_finrw
    ! respiration-weighted annual average of finundated
    real(rk8) , pointer , dimension(:) :: annavg_finrw
    ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
    real(rk8) , pointer , dimension(:) :: sif
    ! Ratio of oxygen available to that demanded by roots, aerobes,
    !  & methanotrophs (nlevsoi)
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
    ! time-lagged surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_surf_lag
    ! time-lagged fractional inundated area
    real(rk8) , pointer , dimension(:) :: finundated_lag
    ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
    real(rk8) , pointer , dimension(:,:) :: layer_sat_lag
  end type column_ch4_type
#endif


  !----------------------------------------------------
  ! column nitrogen state variables structure
  !----------------------------------------------------
  type , public :: column_nstate_type
    !pft-level nitrogen state variables averaged to the column
    type(pft_nstate_type) :: pns_a
    ! all n pools involved in decomposition
    ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8) , pointer , dimension(:,:,:) :: decomp_npools_vr
    ! (gN/m3) vertically-resolved soil mineral N
    real(rk8) , pointer , dimension(:,:) :: sminn_vr
    ! (gN/m3) vertically-resolved column-level sink for N truncation
    real(rk8) , pointer , dimension(:,:) :: col_ntrunc_vr
#ifdef NITRIF_DENITRIF
    ! (gN/m3) vertically-resolved soil mineral NO3
    real(rk8) , pointer , dimension(:,:) :: smin_no3_vr
    ! (gN/m2) soil mineral NO3 pool
    real(rk8) , pointer , dimension(:) :: smin_no3
    ! (gN/m3) vertically-resolved soil mineral NH4
    real(rk8) , pointer , dimension(:,:) :: smin_nh4_vr
    ! (gN/m2) soil mineral NH4 pool
    real(rk8) , pointer , dimension(:) :: smin_nh4
#endif
    ! wood product pools, for dynamic landcover
    ! (gN/m2) column-level pool for seeding new PFTs
    real(rk8) , pointer , dimension(:) :: seedn
    ! (gN/m2) wood product N pool, 10-year lifespan
    real(rk8) , pointer , dimension(:) :: prod10n
    ! (gN/m2) wood product N pool, 100-year lifespan
    real(rk8) , pointer , dimension(:) :: prod100n
    ! (gN/m2) total wood product N
    real(rk8) , pointer , dimension(:) :: totprodn
    ! summary (diagnostic) state variables, not involved in mass balance
    ! (gN/m2)  decomposing (litter, cwd, soil) N pools
    real(rk8) , pointer , dimension(:,:) :: decomp_npools
    ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
    real(rk8) , pointer , dimension(:,:) :: decomp_npools_1m
    ! (gN/m2) soil mineral N
    real(rk8) , pointer , dimension(:) :: sminn
    ! (gN/m2) column-level sink for N truncation
    real(rk8) , pointer , dimension(:) :: col_ntrunc
    ! (gN/m2) Diagnostic: coarse woody debris N
    real(rk8) , pointer , dimension(:) :: cwdn
    ! (gN/m2) total litter nitrogen
    real(rk8) , pointer , dimension(:) :: totlitn
    ! (gN/m2) total soil organic matter nitrogen
    real(rk8) , pointer , dimension(:) :: totsomn
    ! (gN/m2) total litter nitrogen to 1 meter
    real(rk8) , pointer , dimension(:) :: totlitn_1m
    ! (gN/m2) total soil organic matter nitrogen to 1 meter
    real(rk8) , pointer , dimension(:) :: totsomn_1m
    ! (gN/m2) total ecosystem nitrogen, incl veg
    real(rk8) , pointer , dimension(:) :: totecosysn
    ! (gN/m2) total column nitrogen, incl veg
    real(rk8) , pointer , dimension(:) :: totcoln
  end type column_nstate_type

  !----------------------------------------------------
  ! column VOC state variables structure
  !----------------------------------------------------
  type , public :: column_vstate_type
    !pft-level VOC state variables averaged to the column
    type(pft_vstate_type) :: pvs_a
  end type column_vstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! column DGVM state variables structure
  !----------------------------------------------------
  type , public :: column_dgvstate_type
    type(pft_dgvstate_type) :: pdgvs_a
  end type column_dgvstate_type
#endif

  !----------------------------------------------------
  ! column dust state variables structure
  !----------------------------------------------------
  type , public :: column_dstate_type
    real(rk8) , pointer , dimension(:) :: dummy_entry(:)
  end type column_dstate_type

  !----------------------------------------------------
  ! column energy flux variables structure
  !----------------------------------------------------
  type , public :: column_eflux_type
    ! pft-level energy flux variables averaged to the column
    type(pft_eflux_type) :: pef_a
    ! snow melt heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_snomelt
    ! urban snow melt heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_snomelt_u
    ! rural snow melt heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_snomelt_r
    ! implicit evaporation for soil temperature equation
    real(rk8) , pointer , dimension(:) :: eflx_impsoil
    ! ground heat flux between soil layers 1 and 2 (W/m2)
    real(rk8) , pointer , dimension(:) :: eflx_fgr12
    ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: eflx_fgr
    ! Urban variable
    ! heat flux from urban building interior to urban walls, roof (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_building_heat
    ! urban air conditioning flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_urban_ac
    ! urban heating flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_urban_heat
    ! heat flux from beneath the soil or ice column (W/m**2)
    ! positive upward; usually eflx_bot >= 0
    real(rk8) , pointer , dimension(:) :: eflx_bot
  end type column_eflux_type

  !----------------------------------------------------
  ! column momentum flux variables structure
  !----------------------------------------------------
  type , public :: column_mflux_type
    ! pft-level momentum flux variables averaged to the column
    type(pft_mflux_type) ::  pmf_a
  end type column_mflux_type

  !----------------------------------------------------
  ! column water flux variables structure
  !----------------------------------------------------
  type , public :: column_wflux_type
    ! pft-level water flux variables averaged to the column
    type(pft_wflux_type) :: pwf_a
    ! infiltration (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_infl
    ! surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_surf
    ! sub-surface runoff (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_drain
    ! net water input into soil from top (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_top_soil
    ! conversion of h2osfc to ice
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_to_ice
    !surface water runoff
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_surf
    !snow falling on surface water
    real(rk8) , pointer , dimension(:) :: qflx_snow_h2osfc
    ! sub-surface runoff from perched wt (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_drain_perched
    ! flood water flux at column level
    real(rk8) , pointer , dimension(:) :: qflx_floodc
    ! liquid water + ice from layer above soil to top soil layer or sent
    ! to qflx_qrgwl (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_sl_top_soil
    ! snow melt (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_snomelt
    ! snow melt (net)
    real(rk8) , pointer , dimension(:) :: qflx_snow_melt
    ! qflx_surf at glaciers, wetlands, lakes
    real(rk8) , pointer , dimension(:) :: qflx_qrgwl
    ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_runoff
    ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_runoff_u
    ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: qflx_runoff_r
    ! snow melt [mm/s]
    real(rk8) , pointer , dimension(:) :: qmelt
    ! mass balance correction term for dynamic weights
    real(rk8) , pointer , dimension(:) :: h2ocan_loss
    ! soil saturation excess [mm/s]
    real(rk8) , pointer , dimension(:) :: qflx_rsub_sat
    ! dry (BCPHO+BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_bc_dep_dry
    ! wet (BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_bc_dep_wet
    ! hydrophobic BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_bc_dep_pho
    ! hydrophillic BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_bc_dep_phi
    ! total (dry+wet) BC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_bc_dep
    ! dry (OCPHO+OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_oc_dep_dry
    ! wet (OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_oc_dep_wet
    ! hydrophobic OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_oc_dep_pho
    ! hydrophillic OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_oc_dep_phi
    ! total (dry+wet) OC deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_oc_dep
    ! dust species 1 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_dry1
    ! dust species 1 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_wet1
    ! dust species 2 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_dry2
    ! dust species 2 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_wet2
    ! dust species 3 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_dry3
    ! dust species 3 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_wet3
    ! dust species 4 dry deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_dry4
    ! dust species 4 wet deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep_wet4
    ! total (dry+wet) dust deposition on ground (positive definite) (col) [kg/s]
    real(rk8) , pointer , dimension(:) :: flx_dst_dep
    ! snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(rk8) , pointer , dimension(:,:) :: qflx_snofrz_lyr
    ! column-integrated snow freezing rate (positive definite) (col)
    ! [kg m-2 s-1]
    real(rk8) , pointer , dimension(:) :: qflx_snofrz_col
    !irrigation flux (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_irrig
    ! net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
    real(rk8) , pointer , dimension(:) :: qflx_glcice
    ! ice growth (positive definite) (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_glcice_frz
    ! ice melt (positive definite) (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_glcice_melt
    ! ice runoff passed from GLC to CLM (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: glc_rofi
    ! liquid runoff passed from GLC to CLM (mm H2O /s)
    real(rk8) , pointer , dimension(:) :: glc_rofl
  end type column_wflux_type

  !----------------------------------------------------
  ! column carbon flux variables structure
  !----------------------------------------------------
  type , public :: column_cflux_type
    ! pft-level carbon flux variables averaged to the column
    type(pft_cflux_type) :: pcf_a
    ! phenology: litterfall and crop fluxes
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! metabolic pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_c_to_litr_met_c
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! cellulose pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_c_to_litr_cel_c
    ! C fluxes associated with phenology (litterfall and crop) to litter
    ! lignin pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_c_to_litr_lig_c
    ! gap mortality
    ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_c_to_litr_met_c
    ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_c_to_litr_cel_c
    ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_c_to_litr_lig_c
    ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_c_to_cwdc
    ! fire
    ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: fire_mortality_c_to_cwdc
    ! harvest
    ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_c_to_litr_met_c
    ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_c_to_litr_cel_c
    ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_c_to_litr_lig_c
    ! C fluxes associated with harvest to CWD pool (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_c_to_cwdc
    ! new variables for CN code
    ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_to_prod10c
    ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemc_to_prod100c
    ! column-level fire fluxes
    ! vertically-resolved decomposing C fire loss (gC/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: m_decomp_cpools_to_fire_vr
    ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
    real(rk8) , pointer , dimension(:,:) :: m_decomp_cpools_to_fire
    ! C from leaf, froot, xfer and storage C to litter labile C by fire
    ! (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_c_to_litr_met_fire
    ! C from leaf, froot, xfer and storage C to litter cellulose C by fire
    ! (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_c_to_litr_cel_fire
    ! C from leaf, froot, xfer and storage C to litter lignin C by fire
    ! (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_c_to_litr_lig_fire
    ! (gC/m2/s) conversion C flux due to BET and BDT area decreasing
    ! (immediate loss to atm)
    real(rk8) , pointer , dimension(:) :: lf_conv_cflux
    ! (gC/m2/s) carbon emissions due to peat burning
    real(rk8) , pointer , dimension(:) :: somc_fire

    ! decomposition fluxes
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cascade_hr_vr
    ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_cascade_hr
    ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cascade_ctransfer_vr
    ! vertically-integrated (diagnostic) C transferred along deomposition cascade (gC/m2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_cascade_ctransfer
    ! (gC/m3/timestep)  change in decomposing c pools.
    !  Used to update concentrations concurrently with vertical transport
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cpools_sourcesink
    ! rate constant for decomposition (1./sec)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_k
    ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: hr_vr
    ! fraction by which decomposition is limited by anoxia
    real(rk8) , pointer , dimension(:,:) :: o_scalar
    ! fraction by which decomposition is limited by moisture availability
    real(rk8) , pointer , dimension(:,:) :: w_scalar
    ! fraction by which decomposition is limited by temperature
    real(rk8) , pointer , dimension(:,:) :: t_scalar
    ! total SOM C loss from vertical transport (gC/m^2/s)
    real(rk8) , pointer , dimension(:) :: som_c_leached
    ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_cpools_leached
    ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cpools_transport_tendency
#ifdef NITRIF_DENITRIF
    ! potential hr (not N-limited) (gC/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phr_vr
#endif
    ! dynamic landcover fluxes
#ifdef CN
    ! (gC/m2/s) seed source to PFT-level
    real(rk8) , pointer , dimension(:) :: dwt_seedc_to_leaf
    ! (gC/m2/s) seed source to PFT-level
    real(rk8) , pointer , dimension(:) :: dwt_seedc_to_deadstem
    ! (gC/m2/s) conversion C flux (immediate loss to atm)
    real(rk8) , pointer , dimension(:) :: dwt_conv_cflux
    ! (gC/m2/s) addition to 10-yr wood product pool
    real(rk8) , pointer , dimension(:) :: dwt_prod10c_gain
    ! (gC/m2/s) addition to 100-yr wood product pool
    real(rk8) , pointer , dimension(:) :: dwt_prod100c_gain
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootc_to_litr_met_c
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootc_to_litr_cel_c
    ! (gC/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootc_to_litr_lig_c
    ! (gC/m3/s) live coarse root to CWD due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_livecrootc_to_cwdc
    ! (gC/m3/s) dead coarse root to CWD due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_deadcrootc_to_cwdc
    ! (gC/m2/s) total carbon loss from product pools and conversion
    real(rk8) , pointer , dimension(:) :: dwt_closs
    ! (gC/m2/s) dwt_closs+product_closs
    real(rk8) , pointer , dimension(:) :: landuseflux
    ! (gC/m2/s) nee-landuseflux
    real(rk8) , pointer , dimension(:) :: landuptake
    ! wood product pool loss fluxes
    ! (gC/m2/s) decomposition loss from 10-yr wood product pool
    real(rk8) , pointer , dimension(:) :: prod10c_loss
    ! (gC/m2/s) decomposition loss from 100-yr wood product pool
    real(rk8) , pointer , dimension(:) :: prod100c_loss
    ! (gC/m2/s) total wood product carbon loss
    real(rk8) , pointer , dimension(:) :: product_closs
#endif
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) litter heterotrophic respiration
    real(rk8) , pointer , dimension(:) :: lithr
    ! (gC/m2/s) soil organic matter heterotrophic respiration
    real(rk8) , pointer , dimension(:) :: somhr
    ! (gC/m2/s) total heterotrophic respiration
    real(rk8) , pointer , dimension(:) :: hr
    ! (gC/m2/s) total soil respiration (HR + root resp)
    real(rk8) , pointer , dimension(:) :: sr
    ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
    real(rk8) , pointer , dimension(:) :: er
    ! (gC/m2/s) litter fire losses
    real(rk8) , pointer , dimension(:) :: litfire
    ! (gC/m2/s) soil organic matter fire losses
    real(rk8) , pointer , dimension(:) :: somfire
    ! (gC/m2/s) total ecosystem fire losses
    real(rk8) , pointer , dimension(:) :: totfire
    ! (gC/m2/s) net ecosystem production, excludes fire, landuse,
    ! and harvest flux, positive for sink
    real(rk8) , pointer , dimension(:) :: nep
    ! (gC/m2/s) net biome production, includes fire, landuse,
    ! and harvest flux, positive for sink
    real(rk8) , pointer , dimension(:) :: nbp
    ! (gC/m2/s) net ecosystem exchange of carbon, includes fire,
    ! landuse, harvest, and hrv_xsmrpool flux, positive for source
    real(rk8) , pointer , dimension(:) :: nee
    ! (gC/m2/s) total column-level carbon inputs (for balance check)
    real(rk8) , pointer , dimension(:) :: col_cinputs
    ! (gC/m2/s) total column-level carbon outputs (for balance check)
    real(rk8) , pointer , dimension(:) :: col_coutputs

#if (defined CN)
    ! CLAMP summary (diagnostic) flux variables, not involved in mass balance
    ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
    real(rk8) , pointer , dimension(:) :: cwdc_hr
    ! (gC/m2/s) col-level coarse woody debris C loss
    real(rk8) , pointer , dimension(:) :: cwdc_loss
    ! (gC/m2/s) col-level litter C loss
    real(rk8) , pointer , dimension(:) :: litterc_loss
#endif

    ! new variables for fire
    ! (gC/m2/s) total column-level fire C loss
    real(rk8) , pointer , dimension(:) :: col_fire_closs
  end type column_cflux_type

  !----------------------------------------------------
  ! column nitrogen flux variables structure
  !----------------------------------------------------
  type , public :: column_nflux_type
    !pft-level nitrogen flux variables averaged to the column
    type(pft_nflux_type) :: pnf_a
    ! new variables for CN code
    ! deposition fluxes
    ! atmospheric N deposition to soil mineral N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: ndep_to_sminn
    ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: nfix_to_sminn
    ! fertilizer N to soil mineral N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: fert_to_sminn
    ! soybean fixation to soil mineral N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: soyfixn_to_sminn
    ! phenology: litterfall and crop fluxes
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! metabolic pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_n_to_litr_met_n
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! cellulose pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_n_to_litr_cel_n
    ! N fluxes associated with phenology (litterfall and crop) to litter
    ! lignin pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: phenology_n_to_litr_lig_n
    ! gap mortality
    ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_n_to_litr_met_n
    ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_n_to_litr_cel_n
    ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_n_to_litr_lig_n
    ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gap_mortality_n_to_cwdn
    ! fire
    ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: fire_mortality_n_to_cwdn
    ! harvest
    ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_n_to_litr_met_n
    ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_n_to_litr_cel_n
    ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_n_to_litr_lig_n
    ! N fluxes associated with harvest to CWD pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: harvest_n_to_cwdn
    !
    ! dead stem N harvest mortality to 10-year product pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_to_prod10n
    ! dead stem N harvest mortality to 100-year product pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: hrv_deadstemn_to_prod100n
    ! vertically-resolved decomposing N fire loss (gN/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: m_decomp_npools_to_fire_vr
    ! vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
    real(rk8) , pointer , dimension(:,:) :: m_decomp_npools_to_fire
    ! column-level fire N fluxes added by F. Li and S. Levis
    ! N from leaf, froot, xfer and storage N to litter labile N
    ! by fire (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_n_to_litr_met_fire
    ! N from leaf, froot, xfer and storage N to litter cellulose N
    ! by fire (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_n_to_litr_cel_fire
    ! N from leaf, froot, xfer and storage N to litter lignin N
    ! by fire (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: m_n_to_litr_lig_fire

    ! decomposition fluxes
    ! vert-res transfer of N from donor to receiver pool along decomp.
    ! cascade (gN/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cascade_ntransfer_vr
    ! vert-int (diagnostic) transfer of N from donor to receiver pool
    ! along decomp. cascade (gN/m2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_cascade_ntransfer
    ! vert-res mineral N flux for transition along decomposition
    ! cascade (gN/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_cascade_sminn_flux_vr
    ! vert-int (diagnostic) mineral N flux for transition along
    ! decomposition cascade (gN/m2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_cascade_sminn_flux
    ! (gN/m3)  change in decomposing n pools (sum of all additions and
    ! subtractions from stateupdate1).  Used to update concentrations
    ! concurrently with vertical transport
    real(rk8) , pointer , dimension(:,:,:) :: decomp_npools_sourcesink
    ! vertically-resolved immobilization fluxes
    ! vertically-resolved potential N immobilization (gN/m3/s) at each level
    real(rk8) , pointer , dimension(:,:) :: potential_immob_vr
    ! vert-int (diagnostic) potential N immobilization (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: potential_immob
    ! vertically-resolved actual N immobilization (gN/m3/s) at each level
    real(rk8) , pointer , dimension(:,:) :: actual_immob_vr
    ! vert-int (diagnostic) actual N immobilization (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: actual_immob
    ! vertically-resolved plant uptake of soil mineral N (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: sminn_to_plant_vr
    ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: sminn_to_plant
    ! vertically-resolved supplemental N supply (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: supplement_to_sminn_vr
    ! vert-int (diagnostic) supplemental N supply (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: supplement_to_sminn
    ! vertically-resolved gross rate of N mineralization (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: gross_nmin_vr
    ! vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: gross_nmin
    ! vertically-resolved net rate of N mineralization (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: net_nmin_vr
    ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: net_nmin
#ifdef NITRIF_DENITRIF
    ! nitrification / denitrification fluxes
    ! (gN/m3/s) soil nitrification flux
    real(rk8) , pointer , dimension(:,:) :: f_nit_vr
    ! (gN/m3/s) soil denitrification flux
    real(rk8) , pointer , dimension(:,:) :: f_denit_vr
    ! (gN/m2/s) soil nitrification flux
    real(rk8) , pointer , dimension(:) :: f_nit
    ! (gN/m2/s) soil denitrification flux
    real(rk8) , pointer , dimension(:) :: f_denit

    ! (gN/m3/s) potential soil nitrification flux
    real(rk8) , pointer , dimension(:,:) :: pot_f_nit_vr
    ! (gN/m3/s) potential soil denitrification flux
    real(rk8) , pointer , dimension(:,:) :: pot_f_denit_vr
    ! (gN/m2/s) potential soil nitrification flux
    real(rk8) , pointer , dimension(:) :: pot_f_nit
    ! (gN/m2/s) potential soil denitrification flux
    real(rk8) , pointer , dimension(:) :: pot_f_denit
    ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(rk8) , pointer , dimension(:,:) :: n2_n2o_ratio_denit_vr
    ! flux of N2o from denitrification [gN/m^3/s]
    real(rk8) , pointer , dimension(:,:) :: f_n2o_denit_vr
    ! flux of N2o from denitrification [gN/m^2/s]
    real(rk8) , pointer , dimension(:) :: f_n2o_denit
    ! flux of N2o from nitrification [gN/m^3/s]
    real(rk8) , pointer , dimension(:,:) :: f_n2o_nit_vr
    ! flux of N2o from nitrification [gN/m^2/s]
    real(rk8) , pointer , dimension(:) :: f_n2o_nit

    ! immobilization / uptake fluxes
    ! vertically-resolved actual immobilization of NO3 (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: actual_immob_no3_vr
    ! vertically-resolved actual immobilization of NH4 (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: actual_immob_nh4_vr
    ! vertically-resolved plant uptake of soil NO3 (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: smin_no3_to_plant_vr
    ! vertically-resolved plant uptake of soil NH4 (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: smin_nh4_to_plant_vr
    ! actual immobilization of NO3 (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: actual_immob_no3
    ! actual immobilization of NH4 (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: actual_immob_nh4
    ! plant uptake of soil NO3 (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: smin_no3_to_plant
    ! plant uptake of soil Nh4 (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: smin_nh4_to_plant
    ! leaching fluxes
    ! vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: smin_no3_leached_vr
    ! soil mineral NO3 pool loss to leaching (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: smin_no3_leached
    ! vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: smin_no3_runoff_vr
    ! soil mineral NO3 pool loss to runoff (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: smin_no3_runoff

    ! NITRIF_DENITRIF diagnostic quantities
    ! (ugN / g soil) soil nitrate concentration
    real(rk8) , pointer , dimension(:,:) :: smin_no3_massdens_vr
    ! (kg soil / m3) bulk density of soil
    real(rk8) , pointer , dimension(:,:) :: soil_bulkdensity
    real(rk8) , pointer , dimension(:,:) :: k_nitr_t_vr
    real(rk8) , pointer , dimension(:,:) :: k_nitr_ph_vr
    real(rk8) , pointer , dimension(:,:) :: k_nitr_h2o_vr
    real(rk8) , pointer , dimension(:,:) :: k_nitr_vr
    real(rk8) , pointer , dimension(:,:) :: wfps_vr
    real(rk8) , pointer , dimension(:,:) :: fmax_denit_carbonsubstrate_vr
    real(rk8) , pointer , dimension(:,:) :: fmax_denit_nitrate_vr
    ! nitrification and denitrification fluxes
    real(rk8) , pointer , dimension(:,:) :: f_denit_base_vr
    real(rk8) , pointer , dimension(:,:) :: diffus ! diffusivity (m2/s)
    real(rk8) , pointer , dimension(:,:) :: ratio_k1
    real(rk8) , pointer , dimension(:,:) :: ratio_no3_co2
    real(rk8) , pointer , dimension(:,:) :: soil_co2_prod
    real(rk8) , pointer , dimension(:,:) :: fr_WFPS

    real(rk8) , pointer , dimension(:,:) :: r_psi
    real(rk8) , pointer , dimension(:,:) :: anaerobic_frac
#else
    ! denitrification fluxes
    ! vertically-resolved denitrification along decomp cascade (gN/m3/s)
    real(rk8) , pointer , dimension(:,:,:) :: sminn_to_denit_decomp_cascade_vr
    ! vertically-integrated (diagnostic) denitrification along decomp
    ! cascade (gN/m2/s)
    real(rk8) , pointer , dimension(:,:) :: sminn_to_denit_decomp_cascade
    ! vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: sminn_to_denit_excess_vr
    ! vertically-integrated (diagnostic) denitrification from excess
    ! mineral N pool (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: sminn_to_denit_excess
    ! leaching fluxes
    ! vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
    real(rk8) , pointer , dimension(:,:) :: sminn_leached_vr
    ! soil mineral N pool loss to leaching (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: sminn_leached
#endif
    ! dynamic landcover fluxes
    ! (gN/m2/s) seed source to PFT-level
    real(rk8) , pointer , dimension(:) :: dwt_seedn_to_leaf
    ! (gN/m2/s) seed source to PFT-level
    real(rk8) , pointer , dimension(:) :: dwt_seedn_to_deadstem
    ! (gN/m2/s) conversion N flux (immediate loss to atm)
    real(rk8) , pointer , dimension(:) :: dwt_conv_nflux
    ! (gN/m2/s) addition to 10-yr wood product pool
    real(rk8) , pointer , dimension(:) :: dwt_prod10n_gain
    ! (gN/m2/s) addition to 100-yr wood product pool
    real(rk8) , pointer , dimension(:) :: dwt_prod100n_gain
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootn_to_litr_met_n
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootn_to_litr_cel_n
    ! (gN/m3/s) fine root to litter due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_frootn_to_litr_lig_n
    ! (gN/m3/s) live coarse root to CWD due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_livecrootn_to_cwdn
    ! (gN/m3/s) dead coarse root to CWD due to landcover change
    real(rk8) , pointer , dimension(:,:) :: dwt_deadcrootn_to_cwdn
    ! (gN/m2/s) total nitrogen loss from product pools and conversion
    real(rk8) , pointer , dimension(:) :: dwt_nloss
    ! wood product pool loss fluxes
    ! (gN/m2/s) decomposition loss from 10-yr wood product pool
    real(rk8) , pointer , dimension(:) :: prod10n_loss
    ! (gN/m2/s) decomposition loss from 100-yr wood product pool
    real(rk8) , pointer , dimension(:) :: prod100n_loss
    ! (gN/m2/s) total wood product nitrogen loss
    real(rk8) , pointer , dimension(:) :: product_nloss
    ! summary (diagnostic) flux variables, not involved in mass balance
    ! total rate of denitrification (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: denit
    ! column-level N inputs (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: col_ninputs
    ! column-level N outputs (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: col_noutputs
    ! new variables for fire
    ! total column-level fire N loss (gN/m2/s)
    real(rk8) , pointer , dimension(:) :: col_fire_nloss

    ! total SOM N loss from vertical transport (gN/m^2/s)
    real(rk8) , pointer , dimension(:) :: som_n_leached
    ! N loss from vertical transport from each decomposing N pool (gN/m^2/s)
    real(rk8) , pointer , dimension(:,:) :: decomp_npools_leached
    ! N tendency due to vertical transport in decomposing N pools (gN/m^3/s)
    real(rk8) , pointer , dimension(:,:,:) :: decomp_npools_transport_tendency
  end type column_nflux_type

  !----------------------------------------------------
  ! column dust flux variables structure
  !----------------------------------------------------
  type , public :: column_dflux_type
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
  type , public :: landunit_pstate_type
    !column-level physical state variables averaged to landunit
    type(column_pstate_type) :: cps_a
    ! Urban variables
    ! internal building temperature (K)
    real(rk8) , pointer , dimension(:) :: t_building
    ! maximum internal building temperature (K)
    real(rk8) , pointer , dimension(:) :: t_building_max
    ! minimum internal building temperature (K)
    real(rk8) , pointer , dimension(:) :: t_building_min
    ! thermal conductivity of urban wall (W/m/K)
    real(rk8) , pointer , dimension(:,:) :: tk_wall
    ! thermal conductivity of urban roof (W/m/K)
    real(rk8) , pointer , dimension(:,:) :: tk_roof
    ! thermal conductivity of urban impervious road (W/m/K)
    real(rk8) , pointer , dimension(:,:) :: tk_improad
    ! heat capacity of urban wall (J/m^3/K)
    real(rk8) , pointer , dimension(:,:) :: cv_wall
    ! heat capacity of urban roof (J/m^3/K)
    real(rk8) , pointer , dimension(:,:) :: cv_roof
    ! heat capacity of urban impervious road (J/m^3/K)
    real(rk8) , pointer , dimension(:,:) :: cv_improad
    ! total thickness of urban wall (m)
    real(rk8) , pointer , dimension(:) :: thick_wall
    ! total thickness of urban roof (m)
    real(rk8) , pointer , dimension(:) :: thick_roof
    ! number of impervious road layers (-)
    integer(ik4) , pointer , dimension(:) :: nlev_improad
    ! view factor of sky for road
    real(rk8) , pointer , dimension(:) :: vf_sr
    ! view factor of one wall for road
    real(rk8) , pointer , dimension(:) :: vf_wr
    ! view factor of sky for one wall
    real(rk8) , pointer , dimension(:) :: vf_sw
    ! view factor of road for one wall
    real(rk8) , pointer , dimension(:) :: vf_rw
    ! view factor of opposing wall for one wall
    real(rk8) , pointer , dimension(:) :: vf_ww
    ! urban canopy air temperature (K)
    real(rk8) , pointer , dimension(:) :: taf
    ! urban canopy air specific humidity (kg/kg)
    real(rk8) , pointer , dimension(:) :: qaf
    ! direct solar absorbed by roof per unit ground area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_roof_dir
    ! diffuse solar absorbed by roof per unit ground area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_roof_dif
    ! direct  solar absorbed by sunwall per unit wall area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_sunwall_dir
    ! diffuse solar absorbed by sunwall per unit wall area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_sunwall_dif
    ! direct  solar absorbed by shadewall per unit wall area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_shadewall_dir
    ! diffuse solar absorbed by shadewall per unit wall area per unit
    ! incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_shadewall_dif
    ! direct  solar absorbed by impervious road per unit ground area
    ! per unit incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_improad_dir
    ! diffuse solar absorbed by impervious road per unit ground area
    ! per unit incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_improad_dif
    ! direct  solar absorbed by pervious road per unit ground area
    ! per unit incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_perroad_dir
    ! diffuse solar absorbed by pervious road per unit ground area
    ! per unit incident flux
    real(rk8) , pointer , dimension(:,:) :: sabs_perroad_dif
  end type landunit_pstate_type

  !----------------------------------------------------
  ! landunit energy flux variables structure
  !----------------------------------------------------
  type , public :: landunit_eflux_type
    ! column-level energy flux variables averaged to landunit
    type(column_eflux_type) ::  cef_a
    ! Urban variables
    ! multiplicative traffic factor for sensible heat flux
    ! from urban traffic (-)
    real(rk8) , pointer , dimension(:) :: eflx_traffic_factor
    ! traffic sensible heat flux (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_traffic
    ! sensible heat flux from domestic heating/cooling sources of waste
    ! heat (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_wasteheat
    ! sensible heat flux to be put back into canyon due to removal
    ! by AC (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_heat_from_ac
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
  type , public :: gridcell_pstate_type
    !column-level physical state variables averaged to gridcell
    type(column_pstate_type) :: cps_a
  end type gridcell_pstate_type

  !----------------------------------------------------
  ! gridcell energy state variables structure
  !----------------------------------------------------
  type , public :: gridcell_estate_type
    !column-level energy state variables averaged to gridcell
    type(column_estate_type) :: ces_a
    ! initial gridcell total heat content
    real(rk8) , pointer , dimension(:) :: gc_heat1
    ! post land cover change total heat content
    real(rk8) , pointer , dimension(:) :: gc_heat2
  end type gridcell_estate_type

  !----------------------------------------------------
  ! gridcell water state variables structure
  !----------------------------------------------------
  type , public :: gridcell_wstate_type
    !column-level water state variables averaged to gridcell
    type(column_wstate_type) :: cws_a
    ! initial gridcell total h2o liq content
    real(rk8) , pointer , dimension(:) :: gc_liq1
    ! post land cover change total liq content
    real(rk8) , pointer , dimension(:) :: gc_liq2
    ! initial gridcell total h2o liq content
    real(rk8) , pointer , dimension(:) :: gc_ice1
    ! post land cover change total ice content
    real(rk8) , pointer , dimension(:) :: gc_ice2
  end type gridcell_wstate_type

  !----------------------------------------------------
  ! gridcell carbon state variables structure
  !----------------------------------------------------
  type , public :: gridcell_cstate_type
    !column-level carbon state variables averaged to gridcell
    type(column_cstate_type) :: ccs_a
  end type gridcell_cstate_type

#ifdef LCH4
  !----------------------------------------------------
  ! gridcell CH4 flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_ch4_type
    ! Atmospheric conc of CH4, O2, CO2 (mol/m3)
    real(rk8) , pointer , dimension(:,:) :: c_atm
    ! gridcell CO2 production from CH4 oxidation (g C/m**2/s)
    real(rk8) , pointer , dimension(:) :: ch4co2f
    !gridcell average CH4 production (g C/m^2/s)
    real(rk8) , pointer , dimension(:) :: ch4prodg
    !gridcell average net methane correction to CO2 flux (g C/m^2/s)
    real(rk8) , pointer , dimension(:) :: nem
  end type gridcell_ch4_type
#endif

  !----------------------------------------------------
  ! gridcell nitrogen state variables structure
  !----------------------------------------------------
  type , public :: gridcell_nstate_type
    !column-level nitrogen state variables averaged to gridcell
    type(column_nstate_type) :: cns_a
  end type gridcell_nstate_type

  !----------------------------------------------------
  ! gridcell VOC state variables structure
  !----------------------------------------------------
  type , public :: gridcell_vstate_type
    !column-level VOC state variables averaged to gridcell
    type(column_vstate_type):: cvs_a
  end type gridcell_vstate_type

  !----------------------------------------------------
  ! gridcell VOC emission factor variables structure (heald)
  !----------------------------------------------------
  type , public :: gridcell_efstate_type
    ! isoprene emission factors
    real(rk8) , pointer , dimension(:,:) :: efisop
  end type gridcell_efstate_type

  !----------------------------------------------------
  ! gridcell dust state variables structure
  !----------------------------------------------------
  type , public :: gridcell_dstate_type
    !column-level dust state variables averaged to gridcell
    type(column_dstate_type) :: cds_a
  end type gridcell_dstate_type

#if (defined CNDV)
  !----------------------------------------------------
  ! gridcell DGVM state variables structure
  !----------------------------------------------------
  type , public :: gridcell_dgvstate_type
    !20-yr running mean of agdd
    real(rk8) , pointer , dimension(:) :: agdd20
    !20-yr running mean of tmomin
    real(rk8) , pointer , dimension(:) :: tmomin20
    !ann minimum of 10-day running mean (K)
    real(rk8) , pointer , dimension(:) :: t10min
  end type gridcell_dgvstate_type
#endif

  !----------------------------------------------------
  ! gridcell energy flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_eflux_type
    !column-level energy flux variables averaged to gridcell
    type(column_eflux_type) :: cef_a
    ! total grid-level sensible heat flux
    real(rk8) , pointer , dimension(:) :: eflx_sh_totg
    ! dynamic land cover change conversion energy flux
    real(rk8) , pointer , dimension(:) :: eflx_dynbal
  end type gridcell_eflux_type

  !----------------------------------------------------
  ! gridcell momentum flux variables structure
  !-- -------------------------------------------------
  type , public :: gridcell_mflux_type
    !pft-level momentum flux variables averaged to gridcell
    type(pft_mflux_type) :: pmf_a
  end type gridcell_mflux_type

  !----------------------------------------------------
  ! gridcell water flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_wflux_type
    ! column-level water flux variables averaged to gridcell
    type(column_wflux_type) :: cwf_a
    ! total grid-level liq runoff
    real(rk8) , pointer , dimension(:) :: qflx_runoffg
    ! total grid-level ice runoff
    real(rk8) , pointer , dimension(:) :: qflx_snwcp_iceg
    ! liq dynamic land cover change conversion runoff flux
    real(rk8) , pointer , dimension(:) :: qflx_liq_dynbal
    ! ice dynamic land cover change conversion runoff flux
    real(rk8) , pointer , dimension(:) :: qflx_ice_dynbal
    ! total grid-level flood water flux
    real(rk8) , pointer , dimension(:) :: qflx_floodg
  end type gridcell_wflux_type

  !----------------------------------------------------
  ! gridcell carbon flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_cflux_type
    !column-level carbon flux variables averaged to gridcell
    type(column_cflux_type) :: ccf_a
  end type gridcell_cflux_type

  !----------------------------------------------------
  ! gridcell nitrogen flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_nflux_type
    !column-level nitrogen flux variables averaged to gridcell
    type(column_nflux_type) :: cnf_a
  end type gridcell_nflux_type

  !----------------------------------------------------
  ! gridcell dust flux variables structure
  !----------------------------------------------------
  type , public :: gridcell_dflux_type
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

  type , public :: pft_type
    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into column level quantities
    integer(ik4) , pointer , dimension(:) :: column
    !weight (relative to column)
    real(rk8) , pointer , dimension(:) :: wtcol
    !index into landunit level quantities
    integer(ik4) , pointer , dimension(:) :: landunit
    !weight (relative to landunit)
    real(rk8) , pointer , dimension(:) :: wtlunit
    !index into gridcell level quantities
    integer(ik4) , pointer , dimension(:) :: gridcell
    !weight (relative to gridcell)
    real(rk8) , pointer , dimension(:) :: wtgcell

    ! topological mapping functionality
    !pft vegetation
    integer(ik4) , pointer , dimension(:) :: itype
    !m index for laixy(i,j,m),etc.
    integer(ik4) , pointer , dimension(:) :: mxy
    !true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: active

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

  type , public :: column_type
    !plant functional type (pft) data structure
    type(pft_type)   :: p

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into landunit level quantities
    integer(ik4) , pointer , dimension(:) :: landunit
    !weight (relative to landunit)
    real(rk8) , pointer , dimension(:) :: wtlunit
    !index into gridcell level quantities
    integer(ik4) , pointer , dimension(:) :: gridcell
    !weight (relative to gridcell)
    real(rk8) , pointer , dimension(:) :: wtgcell
    !beginning pft index for each column
    integer(ik4) , pointer , dimension(:) :: pfti
    !ending pft index for each column
    integer(ik4) , pointer , dimension(:) :: pftf
    !number of pfts for each column
    integer(ik4) , pointer , dimension(:) :: npfts

    ! topological mapping functionality
    !column type
    integer(ik4) , pointer , dimension(:) :: itype
    !true=>do computations on this column (see reweightMod for details)
    logical , pointer , dimension(:) :: active

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

  type , public :: landunit_type
    ! column data structure (soil/snow/canopy columns)
    type(column_type) :: c

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    !index into gridcell level quantities
    integer(ik4) , pointer , dimension(:) :: gridcell
    !weight (relative to gridcell)
    real(rk8) , pointer , dimension(:) :: wtgcell
    !beginning column index per landunit
    integer(ik4) , pointer , dimension(:) :: coli
    !ending column index for each landunit
    integer(ik4) , pointer , dimension(:) :: colf
    !number of columns for each landunit
    integer(ik4) , pointer , dimension(:) :: ncolumns
    !beginning pft index for each landunit
    integer(ik4) , pointer , dimension(:) :: pfti
    !ending pft index for each landunit
    integer(ik4) , pointer , dimension(:) :: pftf
    !number of pfts for each landunit
    integer(ik4) , pointer , dimension(:) :: npfts

    ! Urban canyon related properties
    ! urban landunit canyon height to width ratio (-)
    real(rk8) , pointer , dimension(:) :: canyon_hwr
    ! urban landunit weight of pervious road column to total road (-)
    real(rk8) , pointer , dimension(:) :: wtroad_perv
    ! weight of roof with respect to urban landunit (-)
    real(rk8) , pointer , dimension(:) :: wtlunit_roof

    ! Urban related info MV - this should be moved to land physical state - MV
    ! height of urban roof (m)
    real(rk8) , pointer , dimension(:) :: ht_roof
    ! height above road at which wind in canyon is to be computed (m)
    real(rk8) , pointer , dimension(:) :: wind_hgt_canyon
    ! urban landunit momentum roughness length (m)
    real(rk8) , pointer , dimension(:) :: z_0_town
    ! urban landunit displacement height (m)
    real(rk8) , pointer , dimension(:) :: z_d_town

    ! topological mapping functionality
    !landunit type
    integer(ik4) , pointer , dimension(:) :: itype
    !true=>landunit is not vegetated
    logical , pointer , dimension(:) :: ifspecial
    !true=>lake point
    logical , pointer , dimension(:) :: lakpoi
    !true=>urban point
    logical , pointer , dimension(:) :: urbpoi
    !true=>glacier_mec point
    logical , pointer , dimension(:) :: glcmecpoi
    !urban density type
    integer(ik4) , pointer , dimension(:) :: udenstype
    !true=>do computations on this landunit (see reweightMod for details)
    logical , pointer , dimension(:) :: active

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

  type , public :: gridcell_type
    !geomorphological landunits
    type(landunit_type) :: l

    ! g/l/c/p hierarchy, local g/l/c/p cells only
    ! Beginning landunit index
    integer(ik4) , pointer , dimension(:) :: luni
    ! Ending landunit index
    integer(ik4) , pointer , dimension(:) :: lunf
    ! Number of landunit for each gridcell
    integer(ik4) , pointer , dimension(:) :: nlandunits
    ! Beginning column index
    integer(ik4) , pointer , dimension(:) :: coli
    ! Ending column index
    integer(ik4) , pointer , dimension(:) :: colf
    ! Number of columns for each gridcell
    integer(ik4) , pointer , dimension(:) :: ncolumns
    ! Beginning pft index
    integer(ik4) , pointer , dimension(:) :: pfti
    ! Ending pft index
    integer(ik4) , pointer , dimension(:) :: pftf
    ! Number of pfts for each gridcell
    integer(ik4) , pointer , dimension(:) :: npfts

    ! Topological mapping functionality, local 1d gdc arrays

    ! Global index
    integer(ik4) , pointer , dimension(:) :: gindex
    ! Total land area, gridcell (km^2)
    real(rk8) , pointer , dimension(:) :: area
    ! Latitude (radians)
    real(rk8) , pointer , dimension(:) :: lat
    ! Longitude (radians)
    real(rk8) , pointer , dimension(:) :: lon
    ! Latitude (degrees)
    real(rk8) , pointer , dimension(:) :: latdeg
    ! Longitude (degrees)
    real(rk8) , pointer , dimension(:) :: londeg
    ! "atm" global index
    integer(ik4) , pointer , dimension(:) :: gindex_a
    ! "atm" latitude (radians) for albedo
    real(rk8) , pointer , dimension(:) :: lat_a
    ! "atm" longitude (radians) for albedo
    real(rk8) , pointer , dimension(:) :: lon_a
    ! "atm" latitude (degrees) for albedo
    real(rk8) , pointer , dimension(:) :: latdeg_a
    ! "atm" longitude (degrees) for albedo
    real(rk8) , pointer , dimension(:) :: londeg_a

    ! Greenland ice sheet mask
    real(rk8) , pointer , dimension(:) :: gris_mask
    ! Greenland ice-covered area per gridcell (km^2)
    real(rk8) , pointer , dimension(:) :: gris_area
    ! Antarctic ice sheet mask
    real(rk8) , pointer , dimension(:) :: aais_mask
    ! Antarctic ice-covered area per gridcell (km^2)
    real(rk8) , pointer , dimension(:) :: aais_area
    ! Total water storage (mm H2O)
    real(rk8) , pointer , dimension(:) :: tws

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

  type , public :: model_type
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
  type(model_type) , public , target , save :: clm3

  !----------------------------------------------------
  ! Declare single instance of array of ecophysiological constant types
  !----------------------------------------------------
  type(pft_epc_type) , public , target , save :: pftcon

  !----------------------------------------------------
  ! Declare single instance of array of decomposition cascade constant types
  !----------------------------------------------------
  type(decomp_cascade_type) , public , target , save :: decomp_cascade_con

#if (defined CNDV)
  !----------------------------------------------------
  ! Declare single instance of array of dgvm ecophysiological constant types
  !----------------------------------------------------
  type(pft_dgvepc_type) , public , target , save :: dgv_pftcon
#endif

  ! name of lndgrid
  character(len=16) , parameter , public :: grlnd  = 'lndgrid'
  ! name of atmgrid
  character(len=16) , parameter , public :: namea  = 'gridcellatm'
  ! name of gridcells
  character(len=16) , parameter , public :: nameg  = 'gridcell'
  ! name of landunits
  character(len=16) , parameter , public :: namel  = 'landunit'
  ! name of columns
  character(len=16) , parameter , public :: namec  = 'column'
  ! name of pfts
  character(len=16) , parameter , public :: namep  = 'pft'

end module mod_clm_type
