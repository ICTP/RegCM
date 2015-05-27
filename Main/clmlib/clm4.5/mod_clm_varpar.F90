module mod_clm_varpar
  !
  ! Module containing CLM parameters
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam , only : num_soil_layers

  implicit none

  private

  save

  ! -------------------------------------------------------
  ! Module Parameters
  ! -------------------------------------------------------

  ! true => run with more vertical soil layers
  logical , public :: more_vertlayers = .false.


  ! Note - model resolution is read in from the surface dataset

  integer(ik4) , public , parameter :: nlev_equalspace   = 15
  integer(ik4) , public , parameter :: toplev_equalspace =  6

  ! Number of hydrologically active soil layers
  integer(ik4) , public :: nlevsoi
  ! Number of soil layers on input file
  integer(ik4) , public :: nlevsoifl
  ! Number of ground layers (includes lower layers that are hydrologically
  ! inactive)
  integer(ik4) , public :: nlevgrnd
  ! Number of urban layers
  integer(ik4) , public :: nlevurb

#ifndef EXTRALAKELAYERS
  ! Number of lake layers
  integer(ik4) , public , parameter :: nlevlak =  10
#else
  ! Yields better results for site simulations. Number of lake layers
  integer(ik4) , public , parameter :: nlevlak =  25
#endif
  ! Maximum number of snow layers
  integer(ik4) , public , parameter :: nlevsno = 5
  ! For CH4 code
  integer(ik4) , public , parameter :: ngases = 3 ! CH4, O2, & CO2

  ! Number of leaf layers in canopy
  integer(ik4) , public , parameter :: nlevcan = 1
  ! Number of water types (soil, ice, 2 lakes, wetland)
  integer(ik4) , public , parameter :: numwat = 5
  ! Number of solar radiation bands: vis, nir
  integer(ik4) , public , parameter :: numrad = 2
  ! Index for visible band
  integer(ik4) , public , parameter :: ivis   = 1
  ! Index for near-infrared band
  integer(ik4) , public , parameter :: inir   = 2
  ! Number of solar type bands: direct, diffuse
  integer(ik4) , public , parameter :: numsolar = 2
  ! Number of dust size classes (BGC only)
  integer(ik4) , public , parameter :: ndst = 4
  ! Number of size distns in src soil (BGC only)
  integer(ik4) , public , parameter :: dst_src_nbr = 3
  ! Number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer(ik4) , public , parameter :: sz_nbr = 200
  ! Maximum number of PFT's for any mode; might we set some of these
  ! automatically from reading pft-physiology?
  integer(ik4) , public , parameter :: mxpft = 24
  ! ADD_MORE_CROP_PFT
  !integer(ik4) , public , parameter :: numveg = 20
  integer(ik4) , public , parameter :: numveg = 16
  ! Number of veg types (without specific crop)
#if (defined VICHYDRO)
  ! Number of VIC soil layer -- Added by AWang
  integer(ik4) , public , parameter :: nlayer = 3
  ! Number of VIC soil layer + 3 lower thermal layers
  integer(ik4) , public , parameter :: nlayert = 8
#endif
#if (defined CROP)
  ! Actual # of pfts (without bare)
  integer(ik4) , public , parameter :: numpft = mxpft
  ! Actual # of crops
  integer(ik4) , public , parameter :: numcft = 10
#else
  ! Actual # of pfts (without bare)
  integer(ik4) , public , parameter :: numpft = numveg
  ! Actual # of crops
  ! ADD_MORE_CROP_PFT
  !integer(ik4) , public , parameter :: numcft = 6
  integer(ik4) , public , parameter :: numcft = 2
#endif
  ! Max number of plant functional types in naturally vegetated landunit
  integer(ik4) , public , parameter :: maxpatch_pft = numpft+1
  ! Number of urban landunits
  integer(ik4) , public , parameter :: numurbl = 3

  ! Number of biogeochemically active soil layers
  integer(ik4) , public :: nlevdecomp
  ! Number of biogeochemical layers (includes lower layers that are
  ! biogeochemically inactive)
  integer(ik4) , public :: nlevdecomp_full

  real(rk8) , public :: outfrq

  ! -------------------------------------------------------
  ! Module Varaibles (initialized in clm_varpar_init)
  ! -------------------------------------------------------

#ifndef CENTURY_DECOMP
  ! parameters for decomposition cascade
  integer(ik4) , public , parameter :: ndecomp_pools = 8
  integer(ik4) , public , parameter :: ndecomp_cascade_transitions = 9
  integer(ik4) , public , parameter :: i_met_lit = 1
  integer(ik4) , public , parameter :: i_cel_lit = 2
  integer(ik4) , public , parameter :: i_lig_lit = 3
  integer(ik4) , public , parameter :: i_cwd = 4
  integer(ik4) , public , parameter :: nsompools = 4
#else
  ! parameters for decomposition cascade
  integer(ik4) , public , parameter :: ndecomp_pools = 7
  integer(ik4) , public , parameter :: ndecomp_cascade_transitions = 10
  integer(ik4) , public , parameter :: i_met_lit = 1
  integer(ik4) , public , parameter :: i_cel_lit = 2
  integer(ik4) , public , parameter :: i_lig_lit = 3
  integer(ik4) , public , parameter :: i_cwd = 4
  integer(ik4) , public , parameter :: nsompools = 3
#endif

  ! Indices used in surface file read and set in clm_varpar_init

  ! Max number of patches
  integer(ik4) , public :: maxpatch
  ! Max number of urban pfts (columns) in urban landunit
  integer(ik4) , public :: maxpatch_urb
  ! Number of urban pfts (columns) in urban landunit
  integer(ik4) , public :: npatch_urban
  ! Number of lake pfts (columns) in lake landunit
  integer(ik4) , public :: npatch_lake
  ! Number of wetland pfts (columns) in wetland landunit
  integer(ik4) , public :: npatch_wet
  ! Number of glacier pfts (columns) in glacier landunit
  integer(ik4) , public :: npatch_glacier
  integer(ik4) , public :: max_pft_per_gcell
  integer(ik4) , public :: max_pft_per_lu
  integer(ik4) , public :: max_pft_per_col
  integer(ik4) , public :: npatch_urban_tbd
  integer(ik4) , public :: npatch_urban_hd
  integer(ik4) , public :: npatch_urban_md

  ! Machine epsilon
  real(rk8) , public :: mach_eps

  public :: clm_varpar_init          ! set parameters

  contains
  !
  ! This subroutine initializes parameters in clm_varpar
  !
  subroutine clm_varpar_init()
    implicit none
    maxpatch_urb   = 5
    npatch_urban_tbd = maxpatch_pft + 1
    npatch_urban_hd  = npatch_urban_tbd + maxpatch_urb
    npatch_urban_md  = npatch_urban_hd + maxpatch_urb
    npatch_lake      = npatch_urban_md + maxpatch_urb
    npatch_wet     = npatch_lake  + 1
    npatch_glacier = npatch_wet   + 1
    maxpatch       = npatch_glacier
    mach_eps       = epsilon(1.0D0)

    max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb*numurbl
#if (defined CROP)
    max_pft_per_gcell = max_pft_per_gcell +  numcft
#endif
    max_pft_per_lu    = max(numpft+1, numcft, maxpatch_urb)
    max_pft_per_col   = max(numpft+1, numcft, maxpatch_urb)

    nlevsoifl   =  num_soil_layers
    nlevurb     =  5
    if ( .not. more_vertlayers ) then
      nlevsoi     =  nlevsoifl
      nlevgrnd    =  15
    else
      nlevsoi     =  8  + nlev_equalspace
      nlevgrnd    =  15 + nlev_equalspace
    end if

    ! Here is a switch to set the number of soil levels for the biogeochemistry
    ! calculations. Currently it works on either a single level or on nlevsoi
    ! and nlevgrnd levels
#ifdef VERTSOILC
    nlevdecomp      = nlevsoi
    nlevdecomp_full = nlevgrnd
#else
    nlevdecomp      = 1
    nlevdecomp_full = 1
#endif
  end subroutine clm_varpar_init

end module mod_clm_varpar
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
