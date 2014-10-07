module mod_clm_slakecon
  !
  ! Module containing constants and parameters for the SLake code
  !   (CLM4-LISSS, documented in Subin et al. 2011, JAMES)
  use mod_intkinds
  use mod_realkinds
  use mod_clm_varpar , only : numrad

  implicit none

  private

  save

  integer(ik4) :: i

  !------------------------------------------------------------------
  ! Lake Model Constants
  !------------------------------------------------------------------

  ! Non-tuneable constants for the lake model
  ! temperature of maximum water density (K)
  ! This is from Hostetler and Bartlein (1990); more updated sources suggest
  ! 277.13 K.
  real(rk8) , parameter , public :: tdmax = 277.D0

  ! Tuneable constants for the lake model
  ! lake emissivity. This
  ! is used for both frozen and unfrozen lakes. This is pulled in from CLM4
  ! and the reference is unclear.
  real(rk8) , public , parameter :: emg_lake = 0.97D0

  ! albedo frozen lakes by waveband (1=vis, 2=nir)
  ! Also unclear what the reference is for this.
  real(rk8) , public , parameter :: alblak(numrad) = (/0.60D0, 0.40D0/)
  ! albedo of melting lakes due to puddling, open water, or white ice
  ! From D. Mironov (2010) Boreal Env. Research
  ! To revert albedo of melting lakes to the cold snow-free value, set
  ! lake_melt_icealb namelist to 0.60, 0.40 like alblak above.
  real(rk8) , public :: alblakwi(numrad)

  ! Coefficient for calculating ice "fraction" for lake surface albedo
  ! From D. Mironov (2010) Boreal Env. Research
  real(rk8) , public , parameter :: calb = 95.6D0
  ! The fraction of the visible (e.g. vis not nir from atm) sunlight
  ! absorbed in ~1 m of water (the surface layer za_lake).
  ! This is roughly the fraction over 700 nm but may depend on the details
  ! of atmospheric radiative transfer.
  ! As long as NIR = 700 nm and up, this can be zero.
  real(rk8) , public :: betavis = 0.0D0
  ! Momentum Roughness length over frozen lakes without snow  (m)
  ! Typical value found in the literature, and consistent with Mironov
  ! expressions. See e.g. Morris EM 1989, Andreas EL 1987,
  ! Guest & Davidson 1991 (as cited in Vavrus 1996)
  real(rk8) , public , parameter :: z0frzlake = 0.001D0
  ! Base of surface light absorption layer for lakes (m)
  real(rk8) , public , parameter :: za_lake = 0.6D0

  ! For calculating prognostic roughness length
  ! min. Charnock parameter
  real(rk8) , public , parameter :: cur0    = 0.01D0
  ! empirical constant for roughness under smooth flow
  real(rk8) , public , parameter :: cus     = 0.1D0
  ! maximum Charnock parameter
  real(rk8) , public , parameter :: curm    = 0.1D0

  !! The following will be set in initSLake based on namelists.
  ! critical dimensionless fetch for Charnock parameter.
  real(rk8) , public :: fcrit
  ! (m) Minimum allowed roughness length for unfrozen lakes.
  real(rk8) , public :: minz0lake

  ! For calculating enhanced diffusivity
  ! (s^-2) (yields diffusivity about 6 times km) ! Fang & Stefan 1996
  real(rk8) , public , parameter :: n2min = 7.5D-5

  real(rk8) , public :: lsadz = 0.03D0 ! m

  ! Note, this will be adjusted in initSLake if the timestep is not 1800 s.
  ! Lake top numerics can oscillate with 0.01m top layer and 1800 s timestep.
  ! The problem is that the surface flux is fixed during the calculation of
  ! the top layer temperature in the diffusion and not corrected for the
  ! tendency of the top layer.
  ! This thickness will be added to all minimum and maximum snow layer
  ! thicknesses compared to that used over non-lakes.
  ! Analysis of the CFL condition suggests that the minimum snow layer
  ! thickness for 1800 s needs to be at least ~1.2 cm for the bulk snow
  ! values of conductivity and heat capacity and as much as 2.3 cm for
  ! pure ice.
  ! Alternatively, a check could be done in SLakeTemperature in case
  ! t_grnd(c) - t_soisno(c,snl(c)+1) changed sign after the Crank-Nicholson
  ! step.
  ! Such an approach, while perhaps allowing additional snow layer resolution,
  ! has not been tested.
  ! The approach used over non-lakes is to have a first-order surface flux
  ! correction.
  ! We choose not to do that here because t_grnd can vary independently of
  ! the top model layer temperature, while it is fixed to the top layer
  ! temperature if tbot > tfrz and the lake is frozen, or if there is an
  ! unstable density gradient in the top unfrozen lake layer.

  !! The following will be set in initSLake based on namelists.
  ! (m) Optional minimum total ice thickness required to allow lake puddling.
  ! Currently used for sensitivity tests only.
  real(rk8) , public :: pudz
  ! (m) Depth beneath which to increase mixing.
  ! See discussion in Subin et al. 2011
  real(rk8) , public :: depthcrit
  real(rk8) , public :: mixfact ! Mixing increase factor.

  !!!!!!!!!!!
  ! Namelists (some of these have not been extensively tested and are
  !            hardwired to default values currently).
  !!!!!!!!!!!

  ! used in SLakeFluxes
  ! true => use old fcrit & minz0 as per Subin et al 2011 form
  ! See initSLakeMod for details. Difference is very small for
  ! small lakes and negligible for large lakes.
  ! Currently hardwired off.
  logical ,  public :: lake_use_old_fcrit_minz0 = .false.

  ! used in SLakeTemperature
  ! Increase mixing by a large factor for deep lakes
  ! Crude but enhanced performance at all 4 deep lakes tested.
  ! See Subin et al 2011 (JAMES) for details
  ! (m) minimum lake depth to invoke deepmixing
  real(rk8) , public :: deepmixing_depthcrit = 25.D0
  ! factor to increase mixing by
  real(rk8) , public :: deepmixing_mixfact   = 10.D0
  ! true => Suppress enhanced diffusion. Small differences.
  ! Currently hardwired .false.
  ! See Subin et al 2011 for details.
  ! Enhanced diffusion is intended for under ice and at large depths.
  ! It is a much smaller change on its own than the "deepmixing"
  ! above, but it increases the effect of deepmixing under ice and for
  ! large depths.
  logical , public :: lake_no_ed = .false.

  ! puddling (not extensively tested and currently hardwired off)
  ! used in SLakeTemperature and SurfaceAlbedo
  ! true => suppress convection when greater than minimum amount
  ! of ice is present. This also effectively sets lake_no_melt_icealb.
  logical , public :: lakepuddling = .false.
  ! (m) minimum amount of total ice nominal thickness before
  ! convection is suppressed
  real(rk8) , public :: lake_puddle_thick = 0.2D0

  ! alblakwi used in SurfaceAlbedo. Will be set by lake_melt_icealb in
  ! initSLake.
  ! Namelist for inputting alblakwi
  real(rk8) , public :: lake_melt_icealb(numrad) = (/ 0.10D0, 0.10D0/)

end module mod_clm_slakecon
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
