module mod_clm_surfacealbedo
  !
  ! Performs surface albedo calculations
  !
  use mod_realkinds
  use mod_intkinds
  use mod_mpmessage
  use mod_stdio
  use mod_dynparam
  use mod_constants , only : mathpi
  use mod_clm_varcon , only : istsoil
  use mod_clm_varpar , only : numrad, nlevcan
  use mod_clm_varcon , only : istcrop
  use mod_clm_varpar , only : nlevsno
  use mod_clm_snicar , only : sno_nbr_aer, SNICAR_RT, DO_SNO_AER, DO_SNO_OC

  implicit none

  ! undefined real
  real(rk8) , public , parameter :: SHR_ORB_UNDEF_REAL = 1.D36
  ! undefined int
  integer(ik4) , public , parameter :: SHR_ORB_UNDEF_INT  = 2000000000

  !----------------------------------------------------------------------------
  ! PRIVATE: by default everything else is private to this module
  !----------------------------------------------------------------------------
  private

  real(rk8) , parameter :: pi                 = mathpi
  real(rk8) , parameter :: SHR_ORB_ECCEN_MIN  =   0.0D0 ! min value for eccen
  real(rk8) , parameter :: SHR_ORB_ECCEN_MAX  =   0.1D0 ! max value for eccen
  real(rk8) , parameter :: SHR_ORB_OBLIQ_MIN  = -90.0D0 ! min value for obliq
  real(rk8) , parameter :: SHR_ORB_OBLIQ_MAX  = +90.0D0 ! max value for obliq
  real(rk8) , parameter :: SHR_ORB_MVELP_MIN  =   0.0D0 ! min value for mvelp
  real(rk8) , parameter :: SHR_ORB_MVELP_MAX  = 360.0D0 ! max value for mvelp

  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
  public :: shr_orb_cosz , shr_orb_decl
  !
  ! The CLM default albice values are too high.
  ! Full-spectral albedo for land ice is ~0.5
  ! (Paterson, Physics of Glaciers, 1994, p. 59)
  ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(rk8) , public :: albice(numrad) = (/ 0.80D0, 0.55D0 /)

  private :: SoilAlbedo    ! Determine ground surface albedo
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer

  contains
  !
  ! Surface albedo and two-stream fluxes
  ! Surface albedos. Also fluxes (per unit incoming direct and diffuse
  ! radiation) reflected, transmitted, and absorbed by vegetation.
  ! Calculate sunlit and shaded fluxes as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy to calculate APAR profile
  ! The calling sequence is:
  ! -> SurfaceAlbedo:   albedos for next time step
  !    -> SoilAlbedo:   soil/lake/glacier/wetland albedos
  !    -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
  !    -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
  !    -> TwoStream:    absorbed, reflected, transmitted solar fluxes
  !                      (vis dir,vis dif, nir dir, nir dif)
  !
  subroutine SurfaceAlbedo(lbg, ubg, lbc, ubc, lbp, ubp, &
                           num_nourbanc, filter_nourbanc, &
                           num_nourbanp, filter_nourbanp, &
                           nextsw_cday, declinp1)
    use mod_clm_type
    use mod_clm_varctl , only : subgridflag
    implicit none
    integer(ik4) , intent(in) :: lbg , ubg  ! gridcell bounds
    integer(ik4) , intent(in) :: lbc , ubc  ! column bounds
    integer(ik4) , intent(in) :: lbp , ubp  ! pft bounds
    ! number of columns in non-urban filter
    integer(ik4) , intent(in) :: num_nourbanc
    ! column filter for non-urban points
    integer(ik4) , intent(in) , dimension(ubc-lbc+1) :: filter_nourbanc
    ! number of pfts in non-urban filter
    integer(ik4) , intent(in) :: num_nourbanp
    ! pft filter for non-urban points
    integer(ik4) , intent(in) , dimension(ubp-lbp+1) :: filter_nourbanp
    ! calendar day at Greenwich (1.00, ..., days/year)
    real(rk8) , intent(in) :: nextsw_cday
    ! declination angle (radians) for next time step
    real(rk8) , intent(in) :: declinp1
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer , dimension(:) :: pactive
    ! gridcell of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pgridcell
    ! index into landunit level quantities
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: itypelun
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! gridcell of corresponding column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! gridcell latitude (radians)
    real(rk8) , pointer , dimension(:) :: lat
    ! gridcell longitude (radians)
    real(rk8) , pointer , dimension(:) :: lon
    ! one-sided leaf area index, no burying by snow
    real(rk8) , pointer , dimension(:) :: tlai
    ! one-sided stem area index, no burying by snow
    real(rk8) , pointer , dimension(:) :: tsai
    ! one-sided leaf area index with burying by snow
    real(rk8) , pointer , dimension(:) :: elai
    ! one-sided stem area index with burying by snow
    real(rk8) , pointer , dimension(:) :: esai
    ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    ! leaf reflectance: 1=vis, 2=nir
    real(rk8) , pointer , dimension(:,:) :: rhol
    ! stem reflectance: 1=vis, 2=nir
    real(rk8) , pointer , dimension(:,:) :: rhos
    ! leaf transmittance: 1=vis, 2=nir
    real(rk8) , pointer , dimension(:,:) :: taul
    ! stem transmittance: 1=vis, 2=nir
    real(rk8) , pointer , dimension(:,:) :: taus
    ! pft vegetation type
    integer(ik4) , pointer , dimension(:) :: ivt
    !
    ! local pointers toimplicit out arguments
    !
    ! cosine of solar zenith angle
    real(rk8) , pointer , dimension(:) :: coszen
    ! ground albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd
    ! ground albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri
    ! surface albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albd
    ! surface albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albi
    ! flux absorbed by canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd
    ! flux absorbed by sunlit canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd_sun
    ! flux absorbed by shaded canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd_sha
    ! flux absorbed by canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi
    ! flux absorbed by sunlit canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi_sun
    ! flux absorbed by shaded canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi_sha
    ! down direct flux below canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: ftdd
    ! down diffuse flux below canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: ftid
    ! down diffuse flux below canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: ftii
    ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsun
    ! leaf to canopy scaling coefficient, shaded leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsha
    ! number of canopy layers
    integer(ik4) , pointer , dimension(:) :: ncan
    ! number of canopy layers, above snow for radiative transfer
    integer(ik4) , pointer , dimension(:) :: nrad
    ! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sun_z
    ! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sha_z
    ! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sun_z
    ! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sha_z
    ! sunlit fraction of canopy layer
    real(rk8) , pointer , dimension(:,:) :: fsun_z
    ! tlai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tlai_z
    ! tsai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tsai_z
    ! solar declination angle (radians)
    real(rk8) , pointer , dimension(:) :: decl
    ! fraction of ground covered by snow (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_sno
    ! liquid water content (col,lyr) [kg/m2]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    ! ice lens content (col,lyr) [kg/m2]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    ! mass concentration of hydrophilic BC (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcphi
    ! mass concentration of hydrophobic BC (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcpho
    ! mass concentration of hydrophilic OC (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocphi
    ! mass concentration of hydrophobic OC (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocpho
    ! mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst1
    ! mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst2
    ! mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst3
    ! mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst4
    ! direct-beam soil albedo (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsod
    ! diffuse soil albedo (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsoi
    ! direct flux absorption factor (col,lyr): VIS [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absdv
    ! direct flux absorption factor (col,lyr): NIR [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absdn
    ! diffuse flux absorption factor (col,lyr): VIS [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absiv
    ! diffuse flux absorption factor (col,lyr): NIR [frc]
    real(rk8) , pointer , dimension(:,:) :: flx_absin
    ! snow grain radius (col,lyr) [microns]
    real(rk8) , pointer , dimension(:,:) :: snw_rds
    ! pure snow ground albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd_pur
    ! pure snow ground albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri_pur
    ! ground albedo without BC (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd_bc
    ! ground albedo without BC (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri_bc
    ! ground albedo without OC (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd_oc
    ! ground albedo without OC (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri_oc
    ! ground albedo without dust (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd_dst
    ! ground albedo without dust (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri_dst
    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsnd_hst
    ! snow ground albedo, diffuse, for history files (col,bnd) [frc]
    real(rk8) , pointer , dimension(:,:) :: albsni_hst

    ! prevents overflow for division by zero
    real(rk8) , parameter :: mpe = 1.D-06
    real(rk8) :: extkn   ! nitrogen allocation coefficient
    integer(ik4) :: fp , fc , g , c , p , iv  ! indices
    integer(ik4) :: ib   ! band index
    integer(ik4) :: ic   ! 0=unit incoming direct; 1=unit incoming diffuse
    real(rk8) , dimension(lbp:ubp) :: wl ! fraction of LAI+SAI that is LAI
    real(rk8) , dimension(lbp:ubp) :: ws ! fraction of LAI+SAI that is SAI
    real(rk8) :: dinc       ! lai+sai increment for canopy layer
    real(rk8) :: dincmax    ! maximum lai+sai increment for canopy layer
    ! cumulative sum of maximum lai+sai increment for canopy layer
    real(rk8) :: dincmax_sum
    ! sum of canopy layer lai for error check
    real(rk8) :: laisum
    real(rk8) :: saisum ! sum of canopy layer sai for error check
    ! lai buried by snow: tlai - elai
    real(rk8) , dimension(lbp:ubp) :: blai
    ! sai buried by snow: tsai - esai
    real(rk8) , dimension(lbp:ubp) :: bsai
    ! leaf/stem refl weighted by fraction LAI and SAI
    real(rk8) , dimension(lbp:ubp,numrad) :: rho
    ! leaf/stem tran weighted by fraction LAI and SAI
    real(rk8) , dimension(lbp:ubp,numrad) :: tau
    real(rk8) , dimension(lbc:ubc,numrad) :: albsnd    ! snow albedo (direct)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsni    ! snow albedo (diffuse)
    ! cosine solar zenith angle for next time step (gridcell level)
    real(rk8) , dimension(lbg:ubg) :: coszen_gcell
    ! cosine solar zenith angle for next time step (pft level)
    real(rk8) , dimension(lbc:ubc) :: coszen_col
    ! cosine solar zenith angle for next time step (pft level)
    real(rk8) , dimension(lbp:ubp) :: coszen_pft
    ! number of vegetated pfts where coszen>0
    integer(ik4) :: num_vegsol
    ! pft filter where vegetated and coszen>0
    integer(ik4) , dimension(ubp-lbp+1) :: filter_vegsol
    ! number of vegetated pfts where coszen>0
    integer(ik4) :: num_novegsol
    ! pft filter where vegetated and coszen>0
    integer(ik4) , dimension(ubp-lbp+1) :: filter_novegsol
    ! number of solar radiation waveband classes
    integer(ik4) , parameter :: nband = numrad
    ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer(ik4) :: flg_slr
    ! flag for SNICAR (=1 when called from CLM, =2 when called from sea-ice)
    integer(ik4) :: flg_snw_ice
    ! direct pure snow albedo (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsnd_pur
    ! diffuse pure snow albedo (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsni_pur
    ! direct snow albedo without BC (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsnd_bc
    ! diffuse snow albedo without BC (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsni_bc
    ! direct snow albedo without OC (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsnd_oc
    ! diffuse snow albedo without OC (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsni_oc
    ! direct snow albedo without dust (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsnd_dst
    ! diffuse snow albedo without dust (radiative forcing)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsni_dst
    integer(ik4) :: i ! index for layers [idx]
    ! flux absorption factor for just snow (direct) [frc]
    real(rk8) , dimension(lbc:ubc,-nlevsno+1:1,numrad) :: flx_absd_snw
    ! flux absorption factor for just snow (diffuse) [frc]
    real(rk8) , dimension(lbc:ubc,-nlevsno+1:1,numrad) :: flx_absi_snw
    ! dummy array for forcing calls
    real(rk8) , dimension(lbc:ubc,-nlevsno+1:1,numrad) :: foo_snw
    ! albedo of surface underneath snow (col,bnd)
    real(rk8) , dimension(lbc:ubc,numrad) :: albsfc
    ! liquid snow content (col,lyr) [kg m-2]
    real(rk8) , dimension(lbc:ubc,-nlevsno+1:0) :: h2osno_liq
    ! ice content in snow (col,lyr) [kg m-2]
    real(rk8) , dimension(lbc:ubc,-nlevsno+1:0) :: h2osno_ice
    ! snow grain size sent to SNICAR (col,lyr) [microns]
    integer(ik4) , dimension(lbc:ubc,-nlevsno+1:0) :: snw_rds_in

    ! mass concentration of aerosol species for forcing calculation
    ! (zero) (col,lyr,aer) [kg kg-1]
    real(rk8) :: mss_cnc_aer_in_frc_pur(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)
    ! mass concentration of aerosol species for BC forcing
    ! (col,lyr,aer) [kg kg-1]
    real(rk8) :: mss_cnc_aer_in_frc_bc(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)
    ! mass concentration of aerosol species for OC forcing
    ! (col,lyr,aer) [kg kg-1]
    real(rk8) :: mss_cnc_aer_in_frc_oc(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)
    ! mass concentration of aerosol species for dust forcing
    ! (col,lyr,aer) [kg kg-1]
    real(rk8) :: mss_cnc_aer_in_frc_dst(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)
    ! mass concentration of all aerosol species for feedback calculation
    ! (col,lyr,aer) [kg kg-1]
    real(rk8) :: mss_cnc_aer_in_fdb(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    lat => clm3%g%lat
    lon => clm3%g%lon

    ! Assign local pointers to derived subtypes components (landunit level)

    itypelun       => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    cgridcell      => clm3%g%l%c%gridcell
    h2osno         => clm3%g%l%c%cws%h2osno
    albgrd         => clm3%g%l%c%cps%albgrd
    albgri         => clm3%g%l%c%cps%albgri
    decl           => clm3%g%l%c%cps%decl
    coszen         => clm3%g%l%c%cps%coszen
    albsod         => clm3%g%l%c%cps%albsod
    albsoi         => clm3%g%l%c%cps%albsoi
    frac_sno       => clm3%g%l%c%cps%frac_sno
    flx_absdv      => clm3%g%l%c%cps%flx_absdv
    flx_absdn      => clm3%g%l%c%cps%flx_absdn
    flx_absiv      => clm3%g%l%c%cps%flx_absiv
    flx_absin      => clm3%g%l%c%cps%flx_absin
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    snw_rds        => clm3%g%l%c%cps%snw_rds
    albgrd_pur     => clm3%g%l%c%cps%albgrd_pur
    albgri_pur     => clm3%g%l%c%cps%albgri_pur
    albgrd_bc      => clm3%g%l%c%cps%albgrd_bc
    albgri_bc      => clm3%g%l%c%cps%albgri_bc
    albgrd_oc      => clm3%g%l%c%cps%albgrd_oc
    albgri_oc      => clm3%g%l%c%cps%albgri_oc
    albgrd_dst     => clm3%g%l%c%cps%albgrd_dst
    albgri_dst     => clm3%g%l%c%cps%albgri_dst
    mss_cnc_bcphi  => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho  => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_cnc_ocphi  => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho  => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_cnc_dst1   => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2   => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3   => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4   => clm3%g%l%c%cps%mss_cnc_dst4
    albsnd_hst     => clm3%g%l%c%cps%albsnd_hst
    albsni_hst     => clm3%g%l%c%cps%albsni_hst

    ! Assign local pointers to derived subtypes components (pft-level)

    pactive   => clm3%g%l%c%p%active
    plandunit => clm3%g%l%c%p%landunit
    pgridcell => clm3%g%l%c%p%gridcell
    pcolumn   => clm3%g%l%c%p%column
    albd      => clm3%g%l%c%p%pps%albd
    albi      => clm3%g%l%c%p%pps%albi
    fabd      => clm3%g%l%c%p%pps%fabd
    fabd_sun  => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha  => clm3%g%l%c%p%pps%fabd_sha
    fabi      => clm3%g%l%c%p%pps%fabi
    fabi_sun  => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha  => clm3%g%l%c%p%pps%fabi_sha
    ftdd      => clm3%g%l%c%p%pps%ftdd
    ftid      => clm3%g%l%c%p%pps%ftid
    ftii      => clm3%g%l%c%p%pps%ftii
    vcmaxcintsun => clm3%g%l%c%p%pps%vcmaxcintsun
    vcmaxcintsha => clm3%g%l%c%p%pps%vcmaxcintsha
    ncan      => clm3%g%l%c%p%pps%ncan
    nrad      => clm3%g%l%c%p%pps%nrad
    fabd_sun_z  => clm3%g%l%c%p%pps%fabd_sun_z
    fabd_sha_z  => clm3%g%l%c%p%pps%fabd_sha_z
    fabi_sun_z  => clm3%g%l%c%p%pps%fabi_sun_z
    fabi_sha_z  => clm3%g%l%c%p%pps%fabi_sha_z
    fsun_z      => clm3%g%l%c%p%pps%fsun_z
    tlai_z    => clm3%g%l%c%p%pps%tlai_z
    tsai_z    => clm3%g%l%c%p%pps%tsai_z
    tlai      => clm3%g%l%c%p%pps%tlai
    tsai      => clm3%g%l%c%p%pps%tsai
    elai      => clm3%g%l%c%p%pps%elai
    esai      => clm3%g%l%c%p%pps%esai
    ivt       => clm3%g%l%c%p%itype
    rhol      => pftcon%rhol
    rhos      => pftcon%rhos
    taul      => pftcon%taul
    taus      => pftcon%taus

    ! Cosine solar zenith angle for next time step

    do g = lbg , ubg
      coszen_gcell(g) = shr_orb_cosz (nextsw_cday, lat(g), lon(g), declinp1)
    end do

    ! Save coszen and declination values to  clm3 data structures for
    ! use in other places in the CN and urban code

    do c = lbc , ubc
      g = cgridcell(c)
      coszen_col(c) = coszen_gcell(g)
      coszen(c) = coszen_col(c)
      decl(c) = declinp1
    end do

    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)
      ! if (pwtgcell(p)>0.D0) then ! "if" added due to chg in filter definition
      g = pgridcell(p)
      coszen_pft(p) = coszen_gcell(g)
      ! end if ! then removed for CNDV (and dyn. landuse?) cases to work
    end do

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1 , numrad
      do fc = 1 , num_nourbanc
        c = filter_nourbanc(fc)
        albsod(c,ib)     = 0.D0
        albsoi(c,ib)     = 0.D0
        albgrd(c,ib)     = 0.D0
        albgri(c,ib)     = 0.D0
        albgrd_pur(c,ib) = 0.D0
        albgri_pur(c,ib) = 0.D0
        albgrd_bc(c,ib)  = 0.D0
        albgri_bc(c,ib)  = 0.D0
        albgrd_oc(c,ib)  = 0.D0
        albgri_oc(c,ib)  = 0.D0
        albgrd_dst(c,ib) = 0.D0
        albgri_dst(c,ib) = 0.D0
        do i = -nlevsno+1 , 1 , 1
          flx_absdv(c,i) = 0.D0
          flx_absdn(c,i) = 0.D0
          flx_absiv(c,i) = 0.D0
          flx_absin(c,i) = 0.D0
        end do
      end do
      do fp = 1 , num_nourbanp
        p = filter_nourbanp(fp)
!       if (pwtgcell(p)>0.D0) then ! "if" added due to chg in filter definition
        albd(p,ib) = 1.D0
        albi(p,ib) = 1.D0
        fabd(p,ib) = 0.D0
        fabd_sun(p,ib) = 0.D0
        fabd_sha(p,ib) = 0.D0
        fabi(p,ib) = 0.D0
        fabi_sun(p,ib) = 0.D0
        fabi_sha(p,ib) = 0.D0
        ftdd(p,ib) = 0.D0
        ftid(p,ib) = 0.D0
        ftii(p,ib) = 0.D0
!       if (ib==1) then
!         ncan(p) = 0
!         nrad(p) = 0
!         do iv = 1, nlevcan
!           fabd_sun_z(p,iv) = 0.D0
!           fabd_sha_z(p,iv) = 0.D0
!           fabi_sun_z(p,iv) = 0.D0
!           fabi_sha_z(p,iv) = 0.D0
!           fsun_z(p,iv) = 0.D0
!           tlai_z(p,iv) = 0.D0
!           tsai_z(p,iv) = 0.D0
!         end do
!       end if
!       end if ! then removed for CNDV (and dyn. landuse?) cases to work
      end do
    end do

    ! SoilAlbedo called before SNICAR_RT
    ! so that reflectance of soil beneath snow column is known
    ! ahead of time for snow RT calculation.

    ! Snow albedos
    ! Note that snow albedo routine will only compute nonzero snow albedos
    ! where h2osno> 0 and coszen > 0

    ! Ground surface albedos
    ! Note that ground albedo routine will only compute nonzero snow albedos
    ! where coszen > 0

    call SoilAlbedo(lbc, ubc, num_nourbanc, filter_nourbanc, &
                    coszen_col, albsnd, albsni)

    ! set variables to pass to SNICAR.

    flg_snw_ice = 1   ! calling from CLM, not CSIM
    do c = lbc , ubc
      albsfc(c,:)     = albsoi(c,:)
      h2osno_liq(c,:) = h2osoi_liq(c,-nlevsno+1:0)
      h2osno_ice(c,:) = h2osoi_ice(c,-nlevsno+1:0)
      snw_rds_in(c,:) = nint(snw_rds(c,:))

      ! zero aerosol input arrays
      mss_cnc_aer_in_frc_pur(c,:,:) = 0.D0
      mss_cnc_aer_in_frc_bc(c,:,:)  = 0.D0
      mss_cnc_aer_in_frc_oc(c,:,:)  = 0.D0
      mss_cnc_aer_in_frc_dst(c,:,:) = 0.D0
      mss_cnc_aer_in_fdb(c,:,:)     = 0.D0
    end do

    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if ( DO_SNO_AER ) then
      mss_cnc_aer_in_fdb(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
      mss_cnc_aer_in_fdb(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)

      ! DO_SNO_OC is set in SNICAR_varpar.
      ! Default case is to ignore OC concentrations because:
      !  1) Knowledge of their optical properties is primitive
      !  2) When 'water-soluble' OPAC optical properties are applied to OC
      !     in snow, it has a negligible darkening effect.
      if ( DO_SNO_OC ) then
        mss_cnc_aer_in_fdb(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
        mss_cnc_aer_in_fdb(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
      end if

      mss_cnc_aer_in_fdb(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
      mss_cnc_aer_in_fdb(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
      mss_cnc_aer_in_fdb(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
      mss_cnc_aer_in_fdb(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)
    end if

    ! If radiative forcing is being calculated, first estimate
    ! clean-snow albedo
#if (defined SNICAR_FRC)
    ! 1. BC input array:
    !  set dust and (optionally) OC concentrations, so
    !       BC_FRC=[(BC+OC+dust)-(OC+dust)]
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)
    if ( DO_SNO_OC ) then
      mss_cnc_aer_in_frc_bc(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
      mss_cnc_aer_in_frc_bc(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
    end if

    ! BC FORCING CALCULATIONS
    flg_slr = 1 ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_bc, albsfc, albsnd_bc, foo_snw)

    flg_slr = 2 ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_bc, albsfc, albsni_bc, foo_snw)

    ! 2. OC input array:
    !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
    if ( DO_SNO_OC ) then
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
      mss_cnc_aer_in_frc_oc(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)

      ! OC FORCING CALCULATIONS
      flg_slr = 1 ! direct-beam
      call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                     coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                     mss_cnc_aer_in_frc_oc, albsfc, albsnd_oc, foo_snw)

      flg_slr = 2 ! diffuse
      call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                     coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                     mss_cnc_aer_in_frc_oc, albsfc, albsni_oc, foo_snw)
    end if

    ! 3. DUST input array:
    ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
    mss_cnc_aer_in_frc_dst(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
    mss_cnc_aer_in_frc_dst(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)
    if ( DO_SNO_OC ) then
      mss_cnc_aer_in_frc_dst(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
      mss_cnc_aer_in_frc_dst(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
    end if

    ! DUST FORCING CALCULATIONS
    flg_slr = 1 ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_dst, albsfc, albsnd_dst, foo_snw)

    flg_slr = 2 ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_dst, albsfc, albsni_dst, foo_snw)

    ! 4. ALL AEROSOL FORCING CALCULATION
    ! (pure snow albedo)
    flg_slr = 1 ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_pur, albsfc, albsnd_pur, foo_snw)

    flg_slr = 2 ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_pur, albsfc, albsni_pur, foo_snw)

#endif

    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
    flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_fdb, albsfc, albsnd, flx_absd_snw)

    flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_fdb, albsfc, albsni, flx_absi_snw)

    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1 , nband
      do fc = 1 , num_nourbanc
        c = filter_nourbanc(fc)
        if (coszen(c) > 0.D0) then
          ! ground albedo was originally computed in SoilAlbedo, but is now
          ! computed here because the order of SoilAlbedo and SNICAR_RT was
          ! switched for SNICAR.
          albgrd(c,ib) = albsod(c,ib)*(1.D0-frac_sno(c)) + &
                  albsnd(c,ib)*frac_sno(c)
          albgri(c,ib) = albsoi(c,ib)*(1.D0-frac_sno(c)) + &
                  albsni(c,ib)*frac_sno(c)

          ! albedos for radiative forcing calculations:
#if (defined SNICAR_FRC)
          ! BC forcing albedo
          albgrd_bc(c,ib) = albsod(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsnd_bc(c,ib)*frac_sno(c)
          albgri_bc(c,ib) = albsoi(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsni_bc(c,ib)*frac_sno(c)
          if ( DO_SNO_OC ) then
            ! OC forcing albedo
            albgrd_oc(c,ib) = albsod(c,ib)*(1.0D0-frac_sno(c)) + &
                    albsnd_oc(c,ib)*frac_sno(c)
            albgri_oc(c,ib) = albsoi(c,ib)*(1.0D0-frac_sno(c)) + &
                    albsni_oc(c,ib)*frac_sno(c)
          end if

          ! dust forcing albedo
          albgrd_dst(c,ib) = albsod(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsnd_dst(c,ib)*frac_sno(c)
          albgri_dst(c,ib) = albsoi(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsni_dst(c,ib)*frac_sno(c)

          ! pure snow albedo for all-aerosol radiative forcing
          albgrd_pur(c,ib) = albsod(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsnd_pur(c,ib)*frac_sno(c)
          albgri_pur(c,ib) = albsoi(c,ib)*(1.0D0-frac_sno(c)) + &
                  albsni_pur(c,ib)*frac_sno(c)
#endif

          ! also in this loop (but optionally in a different loop for
          ! vectorized code) weight snow layer radiative absorption factors
          ! based on snow fraction and soil albedo
          !  (NEEDED FOR ENERGY CONSERVATION)
          do i = -nlevsno+1 , 1 , 1
            if ( subgridflag == 0 ) then
              if ( ib == 1 ) then
                flx_absdv(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.0D0-frac_sno(c))*(1.0D0-albsod(c,ib)) * &
                        (flx_absd_snw(c,i,ib)/(1.0D0-albsnd(c,ib))))
                flx_absiv(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.0D0-frac_sno(c))*(1.0D0-albsoi(c,ib)) * &
                        (flx_absi_snw(c,i,ib)/(1.0D0-albsni(c,ib))))
              else if ( ib == 2 ) then
                flx_absdn(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.0D0-frac_sno(c))*(1.0D0-albsod(c,ib)) * &
                        (flx_absd_snw(c,i,ib)/(1.0D0-albsnd(c,ib))))
                flx_absin(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.0D0-frac_sno(c))*(1.0D0-albsoi(c,ib)) * &
                        (flx_absi_snw(c,i,ib)/(1.0D0-albsni(c,ib))))
              end if
            else
              if ( ib == 1 ) then
                flx_absdv(c,i) = flx_absd_snw(c,i,ib)*(1.0D0-albsnd(c,ib))
                flx_absiv(c,i) = flx_absi_snw(c,i,ib)*(1.0D0-albsni(c,ib))
              else if ( ib == 2 ) then
                flx_absdn(c,i) = flx_absd_snw(c,i,ib)*(1.0D0-albsnd(c,ib))
                flx_absin(c,i) = flx_absi_snw(c,i,ib)*(1.0D0-albsni(c,ib))
              end if
            end if
          end do
        end if
      end do
    end do

    ! for diagnostics, set snow albedo to spval over non-snow points
    ! so that it is not averaged in history buffer
    ! (OPTIONAL)
    do ib = 1 , nband
      do fc = 1 , num_nourbanc
        c = filter_nourbanc(fc)
        if ( (coszen(c) > 0.D0) .and. (h2osno(c) > 0.D0) ) then
          albsnd_hst(c,ib) = albsnd(c,ib)
          albsni_hst(c,ib) = albsni(c,ib)
        else
          albsnd_hst(c,ib) = 0.D0
          albsni_hst(c,ib) = 0.D0
        end if
      end do
    end do

    ! Create solar-vegetated filter for the following calculations

    num_vegsol = 0
    num_novegsol = 0
    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)
      if ( coszen_pft(p) > 0.D0 ) then
        if ( (itypelun(plandunit(p)) == istsoil .or.        &
              itypelun(plandunit(p)) == istcrop     ) .and. &
              (elai(p) + esai(p)) > 0.D0 .and. pactive(p) ) then
          num_vegsol = num_vegsol + 1
          filter_vegsol(num_vegsol) = p
        else
          num_novegsol = num_novegsol + 1
          filter_novegsol(num_novegsol) = p
        end if
      end if
    end do

    ! Weight reflectance/transmittance by lai and sai
    ! Only perform on vegetated pfts where coszen > 0

    do fp = 1 , num_vegsol
      p = filter_vegsol(fp)
      wl(p) = elai(p) / max( elai(p)+esai(p), mpe )
      ws(p) = esai(p) / max( elai(p)+esai(p), mpe )
    end do

    do ib = 1 , numrad
      do fp = 1 , num_vegsol
        p = filter_vegsol(fp)
        rho(p,ib) = max( rhol(ivt(p),ib)*wl(p) + rhos(ivt(p),ib)*ws(p), mpe )
        tau(p,ib) = max( taul(ivt(p),ib)*wl(p) + taus(ivt(p),ib)*ws(p), mpe )
      end do
    end do

    ! Diagnose number of canopy layers for radiative transfer, in increments
    ! of dincmax.
    ! Add to number of layers so long as cumulative leaf+stem area does not
    ! exceed total leaf+stem area. Then add any remaining leaf+stem area to
    ! next layer and exit the loop.
    ! Do this first for elai and esai (not buried by snow) and then for the
    ! part of the canopy that is buried by snow.
    ! ------------------
    ! tlai_z = leaf area increment for a layer
    ! tsai_z = stem area increment for a layer
    ! nrad   = number of canopy layers above snow
    ! ncan   = total number of canopy layers
    !
    ! tlai_z summed from 1 to nrad = elai
    ! tlai_z summed from 1 to ncan = tlai

    ! tsai_z summed from 1 to nrad = esai
    ! tsai_z summed from 1 to ncan = tsai
    ! ------------------
    !
    ! Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol"
    ! because layering is needed for all time steps regardless of radiation
    !
    ! Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered
    ! by nlevcan = 1

    dincmax = 0.25D0
    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)

      if ( nlevcan == 1 ) then
        nrad(p) = 1
        ncan(p) = 1
        tlai_z(p,1) = elai(p)
        tsai_z(p,1) = esai(p)
      else if (nlevcan > 1) then
        if ( elai(p)+esai(p) == 0.D0 ) then
          nrad(p) = 0
        else
          dincmax_sum = 0.D0
          do iv = 1 , nlevcan
            dincmax_sum = dincmax_sum + dincmax
            if ( ((elai(p)+esai(p))-dincmax_sum) > 1.D-06 ) then
              nrad(p) = iv
              dinc = dincmax
              tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
              tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
            else
              nrad(p) = iv
              dinc = dincmax - (dincmax_sum - (elai(p)+esai(p)))
              tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
              tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
              exit
            end if
          end do

          ! Mimumum of 4 canopy layers

          if ( nrad(p) < 4 ) then
            nrad(p) = 4
            do iv = 1 , nrad(p)
              tlai_z(p,iv) = elai(p) / nrad(p)
              tsai_z(p,iv) = esai(p) / nrad(p)
            end do
          end if
        end if
      end if

      ! Error check: make sure cumulative of increments does not exceed total

      laisum = 0.D0
      saisum = 0.D0
      do iv = 1 , nrad(p)
        laisum = laisum + tlai_z(p,iv)
        saisum = saisum + tsai_z(p,iv)
      end do
      if ( abs(laisum-elai(p)) > 1.D-06 .or. &
           abs(saisum-esai(p)) > 1.D-06) then
        write (stderr,*) &
                'multi-layer canopy error 01 in SurfaceAlbedo: ', &
                nrad(p),elai(p),laisum,esai(p),saisum
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if

      ! Repeat to find canopy layers buried by snow

      if ( nlevcan > 1 ) then
        blai(p) = tlai(p) - elai(p)
        bsai(p) = tsai(p) - esai(p)
        if ( blai(p)+bsai(p) == 0.D0 ) then
          ncan(p) = nrad(p)
        else
          dincmax_sum = 0.D0
          do iv = nrad(p)+1 , nlevcan
            dincmax_sum = dincmax_sum + dincmax
            if ( ((blai(p)+bsai(p))-dincmax_sum) > 1.D-06 ) then
              ncan(p) = iv
              dinc = dincmax
              tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
              tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
            else
              ncan(p) = iv
              dinc = dincmax - (dincmax_sum - (blai(p)+bsai(p)))
              tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
              tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
              exit
            end if
          end do
        end if

        ! Error check: make sure cumulative of increments does not exceed total

        laisum = 0.D0
        saisum = 0.D0
        do iv = 1 , ncan(p)
          laisum = laisum + tlai_z(p,iv)
          saisum = saisum + tsai_z(p,iv)
        end do
        if ( abs(laisum-tlai(p)) > 1.D-06 .or. &
             abs(saisum-tsai(p)) > 1.D-06) then
          write (stderr,*) 'multi-layer canopy error 02 in SurfaceAlbedo: ', &
                  nrad(p),ncan(p)
          write (stderr,*) tlai(p),elai(p),blai(p),laisum,tsai(p), &
                  esai(p),bsai(p),saisum
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      end if
    end do

    ! Zero fluxes for active canopy layers

    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)
      do iv = 1 , nrad(p)
        fabd_sun_z(p,iv) = 0.D0
        fabd_sha_z(p,iv) = 0.D0
        fabi_sun_z(p,iv) = 0.D0
        fabi_sha_z(p,iv) = 0.D0
        fsun_z(p,iv) = 0.D0
      end do
    end do

    ! Default leaf to canopy scaling coefficients, used when coszen <= 0.
    ! This is the leaf nitrogen profile integrated over the full canopy.
    ! Integrate exp(-kn*x) over x=0 to x=elai and assign to shaded canopy,
    ! because sunlit fraction is 0. Canopy scaling coefficients are set in
    ! TwoStream for coszen > 0. So kn must be set here and in TwoStream.

    extkn = 0.30D0
    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)
      if ( nlevcan == 1 ) then
        vcmaxcintsun(p) = 0.D0
        vcmaxcintsha(p) = (1.D0 - exp(-extkn*elai(p))) / extkn
        if ( elai(p) > 0.D0 ) then
          vcmaxcintsha(p) = vcmaxcintsha(p) / elai(p)
        else
          vcmaxcintsha(p) = 0.D0
        end if
      else if ( nlevcan > 1 ) then
        vcmaxcintsun(p) = 0.D0
        vcmaxcintsha(p) = 0.D0
      end if
    end do

    ! Calculate surface albedos and fluxes
    ! Only perform on vegetated pfts where coszen > 0

    call TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, &
            num_vegsol, coszen_pft, rho, tau)

    ! Determine values for non-vegetated pfts where coszen > 0

    do ib = 1 , numrad
      do fp = 1 , num_novegsol
        p = filter_novegsol(fp)
        c = pcolumn(p)
        fabd(p,ib) = 0.D0
        fabd_sun(p,ib) = 0.D0
        fabd_sha(p,ib) = 0.D0
        fabi(p,ib) = 0.D0
        fabi_sun(p,ib) = 0.D0
        fabi_sha(p,ib) = 0.D0
        ftdd(p,ib) = 1.D0
        ftid(p,ib) = 0.D0
        ftii(p,ib) = 1.D0
        albd(p,ib) = albgrd(c,ib)
        albi(p,ib) = albgri(c,ib)
!       if (ib == 1) then
!         do iv = 1, nrad(p)
!           fabd_sun_z(p,iv) = 0.D0
!           fabd_sha_z(p,iv) = 0.D0
!           fabi_sun_z(p,iv) = 0.D0
!           fabi_sha_z(p,iv) = 0.D0
!           fsun_z(p,iv) = 0.D0
!         end do
!       end if
      end do
    end do
  end subroutine SurfaceAlbedo
  !
  ! Determine ground surface albedo, accounting for snow
  !
  subroutine SoilAlbedo (lbc, ubc, num_nourbanc, filter_nourbanc, &
                  coszen, albsnd, albsni)
    use mod_clm_type
    use mod_clm_varpar , only : numrad
    use mod_clm_varcon , only : albsat , albdry , tfrz , istice , istice_mec
    use mod_clm_varcon , only : istdlak
    use mod_clm_slakecon , only : alblak
    use mod_clm_slakecon , only : alblakwi , calb , lakepuddling
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! column bounds
    ! number of columns in non-urban points in column filter
    integer(ik4) , intent(in) :: num_nourbanc
    ! column filter for non-urban points
    integer(ik4) , dimension(ubc-lbc+1) , intent(in) :: filter_nourbanc
    ! cos solar zenith angle next time step (column-level)
    real(rk8) , dimension(lbc:ubc) , intent(in) :: coszen
    ! snow albedo (direct)
    real(rk8) , dimension(lbc:ubc,numrad) , intent(in) :: albsnd
    ! snow albedo (diffuse)
    real(rk8) , dimension(lbc:ubc,numrad) , intent(in) :: albsni
    !
    ! local pointers to original implicit in arguments
    !
    ! landunit of corresponding column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! soil color class
    integer(ik4) , pointer , dimension(:) :: isoicol
    ! ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    ! volumetric soil water [m3/m3]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    ! mass fraction of lake layer that is frozen
    real(rk8) , pointer , dimension(:,:) :: lake_icefrac
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    !
    ! local pointers to original implicit out arguments
    !
    ! ground albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albgrd
    ! ground albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albgri
    ! albsod and albsoi are now clm_type variables so they can
    ! be used by SNICAR.
    ! soil albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albsod
    ! soil albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albsoi
    ! number of solar radiation waveband classes
    integer(ik4) , parameter :: nband = numrad
    integer(ik4) :: fc         ! non-urban filter column index
    integer(ik4) :: c , l      ! indices
    integer(ik4) :: ib         ! waveband number (1=vis, 2=nir)
    real(rk8) :: inc           ! soil water correction factor for soil albedo
    ! albsod and albsoi are now clm_type variables so they can be
    ! used by SNICAR.
    !real(rk8) :: albsod       ! soil albedo (direct)
    !real(rk8) :: albsoi       ! soil albedo (diffuse)
    integer(ik4) :: soilcol    ! soilcolor
    real(rk8) :: sicefr        ! Lake surface ice fraction
                               ! (based on D. Mironov 2010)

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit  => clm3%g%l%c%landunit
    isoicol    => clm3%g%l%c%cps%isoicol
    t_grnd     => clm3%g%l%c%ces%t_grnd
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
    albgrd     => clm3%g%l%c%cps%albgrd
    albgri     => clm3%g%l%c%cps%albgri
    albsod     => clm3%g%l%c%cps%albsod
    albsoi     => clm3%g%l%c%cps%albsoi
    snl        => clm3%g%l%c%cps%snl
    lake_icefrac => clm3%g%l%c%cws%lake_icefrac

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype

    ! Compute soil albedos

    do ib = 1 , nband
      do fc = 1 , num_nourbanc
        c = filter_nourbanc(fc)
        if ( coszen(c) > 0.D0 ) then
          l = clandunit(c)
          if ( ltype(l) == istsoil .or. ltype(l) == istcrop )  then ! soil
            inc    = max(0.11D0-0.40D0*h2osoi_vol(c,1), 0.D0)
            soilcol = isoicol(c)
            ! changed from local variable to clm_type:
            !albsod = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
            !albsoi = albsod
            albsod(c,ib) = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
            albsoi(c,ib) = albsod(c,ib)
          else if ( ltype(l) == istice .or. &
                    ltype(l) == istice_mec)  then  ! land ice
            ! changed from local variable to clm_type:
            !albsod = albice(ib)
            !albsoi = albsod
            albsod(c,ib) = albice(ib)
            albsoi(c,ib) = albsod(c,ib)
            ! unfrozen lake, wetland
          else if ( t_grnd(c) > tfrz .or. &
                  (lakepuddling .and. ltype(l) == istdlak .and. &
                   t_grnd(c) == tfrz .and. lake_icefrac(c,1) < 1.D0 .and. &
                   lake_icefrac(c,2) > 0.D0) ) then
            albsod(c,ib) = 0.05D0/(max(0.001D0,coszen(c)) + 0.15D0)
            ! This expression is apparently from BATS according to Yongjiu Dai.
            ! The diffuse albedo should be an average over the whole sky of
            ! an angular-dependent direct expression.
            ! The expression above may have been derived to encompass both
            ! (e.g. Henderson-Sellers 1986),
            ! but I'll assume it applies more appropriately to the direct
            ! form for now.
            ! ZMS: Attn EK, currently restoring this for wetlands even though
            ! it is wrong in order to try to get
            ! bfb baseline comparison when no lakes are present. I'm assuming
            ! wetlands will be phased out anyway.
            if ( ltype(l) == istdlak ) then
              albsoi(c,ib) = 0.10D0
            else
              albsoi(c,ib) = albsod(c,ib)
            end if
          else        ! frozen lake, wetland
            ! Introduce crude surface frozen fraction according to
            ! D. Mironov (2010)
            ! Attn EK: This formulation is probably just as good for
            ! "wetlands" if they are not phased out.
            ! Tenatively I'm restricting this to lakes because I haven't
            ! tested it for wetlands. But if anything the albedo should be
            ! lower when melting over frozen ground than a solid frozen lake.
            !
            if ( ltype(l) == istdlak .and. &
                 .not. lakepuddling .and. snl(c) == 0 ) then
              ! Need to reference snow layers here because t_grnd could be
              ! over snow or ice but we really want the ice surface temperature
              ! with no snow
              sicefr = 1.D0 - exp(-calb * (tfrz - t_grnd(c))/tfrz)
              albsod(c,ib) = sicefr*alblak(ib) + &
                      (1.D0-sicefr)*max(alblakwi(ib), &
                       0.05D0/(max(0.001D0,coszen(c)) + 0.15D0))
              albsoi(c,ib) = sicefr*alblak(ib) + &
                      (1.D0-sicefr)*max(alblakwi(ib), 0.10D0)
              ! Make sure this is no less than the open water albedo above.
              ! Setting lake_melt_icealb(:) = alblak(:) in namelist reverts
              ! the melting albedo to the cold snow-free value.
            else
              albsod(c,ib) = alblak(ib)
              albsoi(c,ib) = albsod(c,ib)
            end if
          end if

          ! Weighting is done in SurfaceAlbedo, after the call to SNICAR_RT
          ! This had to be done, because SoilAlbedo is called before
          ! SNICAR_RT, so at this point, snow albedo is not yet known.
        end if
      end do
    end do
  end subroutine SoilAlbedo
  !
  ! Two-stream fluxes for canopy radiative transfer
  ! Use two-stream approximation of Dickinson (1983) Adv Geophysics
  ! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
  ! to calculate fluxes absorbed by vegetation, reflected by vegetation,
  ! and transmitted through vegetation for unit incoming direct or diffuse
  ! flux given an underlying surface with known albedo.
  ! Calculate sunlit and shaded fluxes as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy to calculate APAR profile
  !
  subroutine TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, &
                  num_vegsol, coszen, rho, tau)
    use mod_clm_type
    use mod_clm_varpar , only : numrad , nlevcan
    use mod_clm_varcon , only : omegas , tfrz , betads , betais
    implicit none
    integer(ik4) , intent(in)  :: lbc , ubc ! column bounds
    integer(ik4) , intent(in)  :: lbp , ubp ! pft bounds
    ! filter for vegetated pfts with coszen>0
    integer(ik4) , dimension(ubp-lbp+1) , intent(in)  :: filter_vegsol
    ! number of vegetated pfts where coszen>0
    integer(ik4) , intent(in) :: num_vegsol
    ! cosine solar zenith angle for next time step
    real(rk8) , dimension(lbp:ubp) , intent(in) :: coszen
    ! leaf/stem refl weighted by fraction LAI and SAI
    real(rk8) , dimension(lbp:ubp,numrad) , intent(in) :: rho
    ! leaf/stem tran weighted by fraction LAI and SAI
    real(rk8) , dimension(lbp:ubp,numrad) , intent(in) :: tau
    ! column of corresponding pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! ground albedo (direct) (column-level)
    real(rk8) , pointer , dimension(:,:) :: albgrd
    ! ground albedo (diffuse)(column-level)
    real(rk8) , pointer , dimension(:,:) :: albgri
    ! vegetation temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_veg
    ! fraction of canopy that is wet (0 to 1)
    real(rk8) , pointer , dimension(:) :: fwet
    ! pft vegetation type
    integer(ik4) , pointer , dimension(:) :: ivt
    ! ecophys const - leaf/stem orientation index
    real(rk8) , pointer , dimension(:) :: xl
    ! one-sided leaf area index with burying by snow
    real(rk8) , pointer , dimension(:) :: elai
    ! one-sided stem area index with burying by snow
    real(rk8) , pointer , dimension(:) :: esai
    ! number of canopy layers, above snow for radiative transfer
    integer(ik4) , pointer , dimension(:) :: nrad
    ! tlai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tlai_z
    ! tsai increment for canopy layer
    real(rk8) , pointer , dimension(:,:) :: tsai_z
    ! surface albedo (direct)
    real(rk8) , pointer , dimension(:,:) :: albd
    ! surface albedo (diffuse)
    real(rk8) , pointer , dimension(:,:) :: albi
    ! flux absorbed by canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd
    ! flux absorbed by sunlit canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd_sun
    ! flux absorbed by shaded canopy per unit direct flux
    real(rk8) , pointer , dimension(:,:) :: fabd_sha
    ! flux absorbed by canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi
    ! flux absorbed by sunlit canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi_sun
    ! flux absorbed by shaded canopy per unit diffuse flux
    real(rk8) , pointer , dimension(:,:) :: fabi_sha
    ! down direct flux below canopy per unit direct flx
    real(rk8) , pointer , dimension(:,:) :: ftdd
    ! down diffuse flux below canopy per unit direct flx
    real(rk8) , pointer , dimension(:,:) :: ftid
    ! down diffuse flux below canopy per unit diffuse flx
    real(rk8) , pointer , dimension(:,:) :: ftii
    ! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sun_z
    ! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabd_sha_z
    ! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sun_z
    ! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8) , pointer , dimension(:,:) :: fabi_sha_z
    ! sunlit fraction of canopy layer
    real(rk8) , pointer , dimension(:,:) :: fsun_z
    ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsun
    ! leaf to canopy scaling coefficient, shaded leaf vcmax
    real(rk8) , pointer , dimension(:) :: vcmaxcintsha
    integer(ik4) :: fp , p , c , iv  ! array indices
    integer(ik4) :: ib               ! waveband number
    real(rk8) :: cosz ! 0.001 <= coszen <= 1.000
    real(rk8) :: asu  ! single scattering albedo
    real(rk8) , dimension(lbp:ubp) :: chil    ! -0.4 <= xl <= 0.6
    ! leaf projection in solar direction (0 to 1)
    real(rk8) , dimension(lbp:ubp) :: gdir
    ! optical depth of direct beam per unit leaf area
    real(rk8) , dimension(lbp:ubp) :: twostext
    ! average diffuse optical depth
    real(rk8) , dimension(lbp:ubp) :: avmu
    ! fraction of intercepted radiation that is scattered (0 to 1)
    real(rk8) , dimension(lbp:ubp,numrad) :: omega
    real(rk8) :: omegal ! omega for leaves
    real(rk8) :: betai  ! upscatter parameter for diffuse radiation
    real(rk8) :: betail ! betai for leaves
    real(rk8) :: betad  ! upscatter parameter for direct beam radiation
    real(rk8) :: betadl ! betad for leaves
    real(rk8) :: tmp0 , tmp1 , tmp2 , tmp3 , tmp4 , tmp5 , tmp6
    real(rk8) :: tmp7 , tmp8 , tmp9 ! temporary
    real(rk8) :: p1 , p2 , p3 , p4 , s1 , s2 , u1 , u2 , u3 ! temporary
    real(rk8) :: b , c1 , d , d1 , d2 , f , h
    real(rk8) :: h1 , h2 , h3 , h4 , h5 , h6 , h7 , h8 , h9 , h10   ! temporary
    real(rk8) :: phi1 , phi2 , sigma , temp1   ! temporary
    real(rk8) , dimension(lbp:ubp) :: temp0 , temp2       ! temporary
    real(rk8) :: t1  ! temporary
    ! parameter for sunlit/shaded leaf radiation absorption
    real(rk8) :: a1 , a2
    ! temporary for flux derivatives
    real(rk8) :: v , dv , u , du
    ! temporary for flux derivatives
    real(rk8) :: dh2 , dh3 , dh5 , dh6 , dh7 , dh8 , dh9 , dh10
    ! temporary for flux derivatives
    real(rk8) :: da1 , da2
    ! ftid, ftii derivative with respect to lai+sai
    real(rk8) :: d_ftid , d_ftii
    ! fabd, fabi derivative with respect to lai+sai
    real(rk8) :: d_fabd , d_fabi
    ! fabd_sun, fabd_sha derivative with respect to lai+sai
    real(rk8) :: d_fabd_sun , d_fabd_sha
    ! fabi_sun, fabi_sha derivative with respect to lai+sai
    real(rk8) :: d_fabi_sun , d_fabi_sha
    ! cumulative lai+sai for canopy layer (at middle of layer)
    real(rk8) :: laisum
    ! direct beam extinction coefficient
    real(rk8) :: extkb
    ! nitrogen allocation coefficient
    real(rk8) :: extkn

    ! Assign local pointers to derived subtypes components (column-level)

    albgrd  => clm3%g%l%c%cps%albgrd
    albgri  => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn  => clm3%g%l%c%p%column
    fwet     => clm3%g%l%c%p%pps%fwet
    t_veg    => clm3%g%l%c%p%pes%t_veg
    ivt      => clm3%g%l%c%p%itype
    elai     => clm3%g%l%c%p%pps%elai
    esai     => clm3%g%l%c%p%pps%esai
    albd     => clm3%g%l%c%p%pps%albd
    albi     => clm3%g%l%c%p%pps%albi
    fabd     => clm3%g%l%c%p%pps%fabd
    fabd_sun => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha => clm3%g%l%c%p%pps%fabd_sha
    fabi     => clm3%g%l%c%p%pps%fabi
    fabi_sun => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha => clm3%g%l%c%p%pps%fabi_sha
    ftdd     => clm3%g%l%c%p%pps%ftdd
    ftid     => clm3%g%l%c%p%pps%ftid
    ftii     => clm3%g%l%c%p%pps%ftii
    xl       => pftcon%xl
    nrad     => clm3%g%l%c%p%pps%nrad
    fabd_sun_z  => clm3%g%l%c%p%pps%fabd_sun_z
    fabd_sha_z  => clm3%g%l%c%p%pps%fabd_sha_z
    fabi_sun_z  => clm3%g%l%c%p%pps%fabi_sun_z
    fabi_sha_z  => clm3%g%l%c%p%pps%fabi_sha_z
    fsun_z      => clm3%g%l%c%p%pps%fsun_z
    tlai_z   => clm3%g%l%c%p%pps%tlai_z
    tsai_z   => clm3%g%l%c%p%pps%tsai_z
    vcmaxcintsun => clm3%g%l%c%p%pps%vcmaxcintsun
    vcmaxcintsha => clm3%g%l%c%p%pps%vcmaxcintsha

    ! Calculate two-stream parameters that are independent of waveband:
    ! chil, gdir, twostext, avmu, and temp0 and temp2 (used for asu)

    do fp = 1 , num_vegsol
      p = filter_vegsol(fp)

      ! note that the following limit only acts on cosz values > 0 and
      ! less than 0.001, not on values cosz = 0, since these zero have
      ! already been filtered out in filter_vegsol
      cosz = max(0.001D0, coszen(p))

      chil(p) = min( max(xl(ivt(p)), -0.4D0), 0.6D0 )
      if (abs(chil(p)) <= 0.01D0) chil(p) = 0.01D0
      phi1 = 0.5D0 - 0.633D0*chil(p) - 0.330D0*chil(p)*chil(p)
      phi2 = 0.877D0 * (1.D0-2.D0*phi1)
      gdir(p) = phi1 + phi2*cosz
      twostext(p) = gdir(p)/cosz
      avmu(p) = ( 1.D0 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
      temp0(p) = gdir(p) + phi2*cosz
      temp1 = phi1*cosz
      temp2(p) = ( 1.D0 - temp1/temp0(p) * log((temp1+temp0(p))/temp1) )
    end do
    !
    ! Loop over all wavebands to calculate for the full canopy the scattered
    ! fluxes reflected upward and transmitted downward by the canopy and the
    ! flux absorbed by the canopy for a unit incoming direct beam and diffuse
    ! flux at the top of the canopy given an underlying surface of known albedo.
    !
    ! Output:
    ! ------------------
    ! Direct beam fluxes
    ! ------------------
    ! albd       - Upward scattered flux above canopy
    !               (per unit direct beam flux)
    ! ftid       - Downward scattered flux below canopy
    !               (per unit direct beam flux)
    ! ftdd       - Transmitted direct beam flux below canopy
    !               (per unit direct beam flux)
    ! fabd       - Flux absorbed by canopy
    !               (per unit direct beam flux)
    ! fabd_sun   - Sunlit portion of fabd
    ! fabd_sha   - Shaded portion of fabd
    ! fabd_sun_z - absorbed sunlit leaf direct PAR
    !               (per unit sunlit lai+sai) for each canopy layer
    ! fabd_sha_z - absorbed shaded leaf direct PAR
    !               (per unit shaded lai+sai) for each canopy layer
    ! ------------------
    ! Diffuse fluxes
    ! ------------------
    ! albi       - Upward scattered flux above canopy (per unit diffuse flux)
    ! ftii       - Downward scattered flux below canopy (per unit diffuse flux)
    ! fabi       - Flux absorbed by canopy (per unit diffuse flux)
    ! fabi_sun   - Sunlit portion of fabi
    ! fabi_sha   - Shaded portion of fabi
    ! fabi_sun_z - absorbed sunlit leaf diffuse PAR
    !               (per unit sunlit lai+sai) for each canopy layer
    ! fabi_sha_z - absorbed shaded leaf diffuse PAR
    !               (per unit shaded lai+sai) for each canopy layer
    !
    do ib = 1 , numrad
      do fp = 1 , num_vegsol
        p = filter_vegsol(fp)
        c = pcolumn(p)
        ! Calculate two-stream parameters omega, betad, and betai.
        ! Omega, betad, betai are adjusted for snow. Values for omega*betad
        ! and omega*betai are calculated and then divided by the new omega
        ! because the product omega*betai, omega*betad is used in solution.
        ! Also, the transmittances and reflectances (tau, rho) are linear
        ! weights of leaf and stem values.
        omegal = rho(p,ib) + tau(p,ib)
        asu = 0.5D0*omegal*gdir(p)/temp0(p) *temp2(p)
        betadl = (1.D0+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu
        betail = 0.5D0 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) * &
                 ((1.D0+chil(p))/2.D0)**2) / omegal
        ! Adjust omega, betad, and betai for intercepted snow
        if (t_veg(p) > tfrz) then                             !no snow
          tmp0 = omegal
          tmp1 = betadl
          tmp2 = betail
        else
          tmp0 =   (1.D0-fwet(p))*omegal        + &
                  fwet(p)*omegas(ib)
          tmp1 = ( (1.D0-fwet(p))*omegal*betadl + &
                  fwet(p)*omegas(ib)*betads ) / tmp0
          tmp2 = ( (1.D0-fwet(p))*omegal*betail + &
                  fwet(p)*omegas(ib)*betais ) / tmp0
        end if
        omega(p,ib) = tmp0
        betad = tmp1
        betai = tmp2

        ! Common terms

        b = 1.D0 - omega(p,ib) + omega(p,ib)*betai
        c1 = omega(p,ib)*betai
        tmp0 = avmu(p)*twostext(p)
        d = tmp0 * omega(p,ib)*betad
        f = tmp0 * omega(p,ib)*(1.D0-betad)
        tmp1 = b*b - c1*c1
        h = sqrt(tmp1) / avmu(p)
        sigma = tmp0*tmp0 - tmp1
        p1 = b + avmu(p)*h
        p2 = b - avmu(p)*h
        p3 = b + tmp0
        p4 = b - tmp0

        ! Absorbed, reflected, transmitted fluxes per unit incoming radiation
        ! for full canopy

        t1 = min(h*(elai(p)+esai(p)), 40.D0)
        s1 = exp(-t1)
        t1 = min(twostext(p)*(elai(p)+esai(p)), 40.D0)
        s2 = exp(-t1)

        ! Direct beam

        u1 = b - c1/albgrd(c,ib)
        u2 = b - c1*albgrd(c,ib)
        u3 = f + c1*albgrd(c,ib)
        tmp2 = u1 - avmu(p)*h
        tmp3 = u1 + avmu(p)*h
        d1 = p1*tmp2/s1 - p2*tmp3*s1
        tmp4 = u2 + avmu(p)*h
        tmp5 = u2 - avmu(p)*h
        d2 = tmp4/s1 - tmp5*s1
        h1 = -d*p4 - c1*f
        tmp6 = d - h1*p3/sigma
        tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
        h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
        h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
        h4 = -f*p3 - c1*d
        tmp8 = h4/sigma
        tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
        h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
        h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

        albd(p,ib) = h1/sigma + h2 + h3
        ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1
        ftdd(p,ib) = s2
        fabd(p,ib) = 1.D0 - albd(p,ib) - &
                (1.D0-albgrd(c,ib))*ftdd(p,ib) - (1.D0-albgri(c,ib))*ftid(p,ib)

        a1 = h1 / sigma * (1.D0 - s2*s2) / (2.D0 * twostext(p)) &
             + h2         * (1.D0 - s2*s1) / (twostext(p) + h) &
             + h3         * (1.D0 - s2/s1) / (twostext(p) - h)

        a2 = h4 / sigma * (1.D0 - s2*s2) / (2.D0 * twostext(p)) &
             + h5         * (1.D0 - s2*s1) / (twostext(p) + h) &
             + h6         * (1.D0 - s2/s1) / (twostext(p) - h)

        fabd_sun(p,ib) = (1.D0 - omega(p,ib)) * &
                ( 1.D0 - s2 + 1.D0 / avmu(p) * (a1 + a2) )
        fabd_sha(p,ib) = fabd(p,ib) - fabd_sun(p,ib)

        ! Diffuse

        u1 = b - c1/albgri(c,ib)
        u2 = b - c1*albgri(c,ib)
        tmp2 = u1 - avmu(p)*h
        tmp3 = u1 + avmu(p)*h
        d1 = p1*tmp2/s1 - p2*tmp3*s1
        tmp4 = u2 + avmu(p)*h
        tmp5 = u2 - avmu(p)*h
        d2 = tmp4/s1 - tmp5*s1
        h7 = (c1*tmp2) / (d1*s1)
        h8 = (-c1*tmp3*s1) / d1
        h9 = tmp4 / (d2*s1)
        h10 = (-tmp5*s1) / d2

        albi(p,ib) = h7 + h8
        ftii(p,ib) = h9*s1 + h10/s1
        fabi(p,ib) = 1.D0 - albi(p,ib) - (1.D0-albgri(c,ib))*ftii(p,ib)

        a1 = h7 * (1.D0 - s2*s1) / (twostext(p) + h) +  &
                h8 * (1.D0 - s2/s1) / (twostext(p) - h)
        a2 = h9 * (1.D0 - s2*s1) / (twostext(p) + h) + &
                h10 * (1.D0 - s2/s1) / (twostext(p) - h)

        fabi_sun(p,ib) = (1.D0 - omega(p,ib)) / avmu(p) * (a1 + a2)
        fabi_sha(p,ib) = fabi(p,ib) - fabi_sun(p,ib)

        ! Repeat two-stream calculations for each canopy layer to
        ! calculate derivatives.
        ! tlai_z and tsai_z are the leaf+stem area increment for a layer.
        ! Derivatives are calculated at the center of the layer.
        ! Derivatives are needed only for the visible waveband to calculate
        ! absorbed PAR (per unit lai+sai) for each canopy layer.
        ! Derivatives are calculated first per unit lai+sai and then
        ! normalized for sunlit or shaded fraction of canopy layer.

        ! Sun/shade big leaf code uses only one layer, with canopy
        ! integrated values from above and also canopy-integrated
        ! scaling coefficients

        if ( ib == 1 ) then
          if ( nlevcan == 1 ) then

            ! sunlit fraction of canopy
            fsun_z(p,1) = (1.D0 - s2) / t1

            ! absorbed PAR (per unit sun/shade lai+sai)
            laisum = elai(p)+esai(p)
            fabd_sun_z(p,1) = fabd_sun(p,ib) / (fsun_z(p,1)*laisum)
            fabi_sun_z(p,1) = fabi_sun(p,ib) / (fsun_z(p,1)*laisum)
            fabd_sha_z(p,1) = fabd_sha(p,ib) / ((1.D0 - fsun_z(p,1))*laisum)
            fabi_sha_z(p,1) = fabi_sha(p,ib) / ((1.D0 - fsun_z(p,1))*laisum)

            ! leaf to canopy scaling coefficients
            extkn = 0.30D0
            extkb = twostext(p)
            vcmaxcintsun(p) = (1.D0 - exp(-(extkn+extkb)*elai(p))) / &
                    (extkn + extkb)
            vcmaxcintsha(p) = (1.D0 - exp(-extkn*elai(p))) / &
                    extkn - vcmaxcintsun(p)
            if ( elai(p) .gt. 0.D0 ) then
              vcmaxcintsun(p) = vcmaxcintsun(p) / (fsun_z(p,1)*elai(p))
              vcmaxcintsha(p) = vcmaxcintsha(p) / ((1.D0 - fsun_z(p,1))*elai(p))
            else
              vcmaxcintsun(p) = 0.D0
              vcmaxcintsha(p) = 0.D0
            end if

          else if ( nlevcan > 1 ) then
            do iv = 1 , nrad(p)

              ! Cumulative lai+sai at center of layer

              if ( iv == 1 ) then
                laisum = 0.5D0 * (tlai_z(p,iv)+tsai_z(p,iv))
              else
                laisum = laisum + 0.5D0 * &
                     ((tlai_z(p,iv-1)+tsai_z(p,iv-1)) + &
                      (tlai_z(p,iv)+tsai_z(p,iv)))
              end if

              ! Coefficients s1 and s2 depend on cumulative lai+sai.
              ! s2 is the sunlit fraction

              t1 = min(h*laisum, 40.D0)
              s1 = exp(-t1)
              t1 = min(twostext(p)*laisum, 40.D0)
              s2 = exp(-t1)
              fsun_z(p,iv) = s2

              ! ===============
              ! Direct beam
              ! ===============

              ! Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai

              u1 = b - c1/albgrd(c,ib)
              u2 = b - c1*albgrd(c,ib)
              u3 = f + c1*albgrd(c,ib)
              tmp2 = u1 - avmu(p)*h
              tmp3 = u1 + avmu(p)*h
              d1 = p1*tmp2/s1 - p2*tmp3*s1
              tmp4 = u2 + avmu(p)*h
              tmp5 = u2 - avmu(p)*h
              d2 = tmp4/s1 - tmp5*s1
              h1 = -d*p4 - c1*f
              tmp6 = d - h1*p3/sigma
              tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
              h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
              h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
              h4 = -f*p3 - c1*d
              tmp8 = h4/sigma
              tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
              h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
              h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

              a1 = h1 / sigma * (1.D0 - s2*s2) / (2.D0 * twostext(p)) &
                   + h2         * (1.D0 - s2*s1) / (twostext(p) + h) &
                   + h3         * (1.D0 - s2/s1) / (twostext(p) - h)

              a2 = h4 / sigma * (1.D0 - s2*s2) / (2.D0 * twostext(p)) &
                   + h5         * (1.D0 - s2*s1) / (twostext(p) + h) &
                   + h6         * (1.D0 - s2/s1) / (twostext(p) - h)

              ! Derivatives for h2, h3, h5, h6 and a1, a2

              v = d1
              dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

              u = tmp6 * tmp2 / s1 - p2 * tmp7
              du = h * tmp6 * tmp2 / s1 + twostext(p) * p2 * tmp7
              dh2 = (v * du - u * dv) / (v * v)

              u = -tmp6 * tmp3 * s1 + p1 * tmp7
              du = h * tmp6 * tmp3 * s1 - twostext(p) * p1 * tmp7
              dh3 = (v * du - u * dv) / (v * v)

              v = d2
              dv = h * tmp4 / s1 + h * tmp5 * s1

              u = -h4/sigma * tmp4 / s1 - tmp9
              du = -h * h4/sigma * tmp4 / s1 + twostext(p) * tmp9
              dh5 = (v * du - u * dv) / (v * v)

              u = h4/sigma * tmp5 * s1 + tmp9
              du = -h * h4/sigma * tmp5 * s1 - twostext(p) * tmp9
              dh6 = (v * du - u * dv) / (v * v)

              da1 = h1/sigma * s2*s2 + h2 * s2*s1 + h3 * s2/s1 &
                    + (1.D0 - s2*s1) / (twostext(p) + h) * dh2 &
                    + (1.D0 - s2/s1) / (twostext(p) - h) * dh3
              da2 = h4/sigma * s2*s2 + h5 * s2*s1 + h6 * s2/s1 &
                    + (1.D0 - s2*s1) / (twostext(p) + h) * dh5 &
                    + (1.D0 - s2/s1) / (twostext(p) - h) * dh6

              ! Flux derivatives

              d_ftid = -twostext(p)*h4/sigma*s2 - &
                      h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1
              d_fabd = -(dh2+dh3) + &
                      (1.D0-albgrd(c,ib))*twostext(p)*s2 - &
                      (1.D0-albgri(c,ib))*d_ftid
              d_fabd_sun = (1.D0 - omega(p,ib)) * &
                      (twostext(p)*s2 + 1.D0 / avmu(p) * (da1 + da2))
              d_fabd_sha = d_fabd - d_fabd_sun

              fabd_sun_z(p,iv) = max(d_fabd_sun, 0.D0)
              fabd_sha_z(p,iv) = max(d_fabd_sha, 0.D0)

              ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI).
              ! Need to normalize derivatives by sunlit or shaded fraction
              ! to get APARsun per unit (LAI+SAI)sun and APARsha per unit
              ! (LAI+SAI)sha

              fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / fsun_z(p,iv)
              fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / (1.D0 - fsun_z(p,iv))

              ! ===============
              ! Diffuse
              ! ===============

              ! Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai

              u1 = b - c1/albgri(c,ib)
              u2 = b - c1*albgri(c,ib)
              tmp2 = u1 - avmu(p)*h
              tmp3 = u1 + avmu(p)*h
              d1 = p1*tmp2/s1 - p2*tmp3*s1
              tmp4 = u2 + avmu(p)*h
              tmp5 = u2 - avmu(p)*h
              d2 = tmp4/s1 - tmp5*s1
              h7 = (c1*tmp2) / (d1*s1)
              h8 = (-c1*tmp3*s1) / d1
              h9 = tmp4 / (d2*s1)
              h10 = (-tmp5*s1) / d2

              a1 = h7 * (1.D0 - s2*s1) / (twostext(p) + h) +  &
                      h8 * (1.D0 - s2/s1) / (twostext(p) - h)
              a2 = h9 * (1.D0 - s2*s1) / (twostext(p) + h) + &
                      h10 * (1.D0 - s2/s1) / (twostext(p) - h)

              ! Derivatives for h7, h8, h9, h10 and a1, a2

              v = d1
              dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

              u = c1 * tmp2 / s1
              du = h * c1 * tmp2 / s1
              dh7 = (v * du - u * dv) / (v * v)

              u = -c1 * tmp3 * s1
              du = h * c1 * tmp3 * s1
              dh8 = (v * du - u * dv) / (v * v)

              v = d2
              dv = h * tmp4 / s1 + h * tmp5 * s1

              u = tmp4 / s1
              du = h * tmp4 / s1
              dh9 = (v * du - u * dv) / (v * v)

              u = -tmp5 * s1
              du = h * tmp5 * s1
              dh10 = (v * du - u * dv) / (v * v)

              da1 = h7*s2*s1 +  &
                      h8*s2/s1 + (1.D0-s2*s1)/(twostext(p)+h)*dh7 + &
                      (1.D0-s2/s1)/(twostext(p)-h)*dh8
              da2 = h9*s2*s1 + h10*s2/s1 + &
                      (1.D0-s2*s1)/(twostext(p)+h)*dh9 + &
                      (1.D0-s2/s1)/(twostext(p)-h)*dh10

              ! Flux derivatives

              d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1
              d_fabi = -(dh7+dh8) - (1.D0-albgri(c,ib))*d_ftii
              d_fabi_sun = (1.D0 - omega(p,ib)) / avmu(p) * (da1 + da2)
              d_fabi_sha = d_fabi - d_fabi_sun

              fabi_sun_z(p,iv) = max(d_fabi_sun, 0.D0)
              fabi_sha_z(p,iv) = max(d_fabi_sha, 0.D0)

              ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI).
              ! Need to normalize derivatives by sunlit or shaded fraction
              ! to get APARsun per unit (LAI+SAI)sun and APARsha per unit
              ! (LAI+SAI)sha

              fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / fsun_z(p,iv)
              fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / (1.D0 - fsun_z(p,iv))

            end do   ! end of canopy layer loop
          end if
        end if
      end do   ! end of pft loop
    end do   ! end of radiation band loop
  end subroutine TwoStream
  !
  ! FUNCTION to return the cosine of the solar zenith angle.
  ! Assumes 365.0 days/year.
  !
  real(rk8) function shr_orb_cosz(jday,lat,lon,declin)

    real(rk8) , intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real(rk8) , intent(in) :: lat    ! Centered latitude (radians)
    real(rk8) , intent(in) :: lon    ! Centered longitude (radians)
    real(rk8) , intent(in) :: declin ! Solar declination (radians)

    shr_orb_cosz = sin(lat)*sin(declin) - &
                   cos(lat)*cos(declin)*cos(jday*2.0D0*pi + lon)

  end function shr_orb_cosz
  !
  ! Calculate earths orbital parameters using Dave Threshers formula which
  ! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term
  ! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy
  ! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
  !
  !
  subroutine shr_orb_params( iyear_AD , eccen  , obliq , mvelp     ,     &
                            obliqr   , lambm0 , mvelpp, log_print )
    integer(ik4) , intent(in) :: iyear_AD  ! Year to calculate orbit for
    real(rk8) , intent(inout) :: eccen   ! orbital eccentricity
    real(rk8) , intent(inout) :: obliq   ! obliquity in degrees
    real(rk8) , intent(inout) :: mvelp   ! moving vernal equinox long
    real(rk8) , intent(out) :: obliqr    ! Earths obliquity in rad
    real(rk8) , intent(out) :: lambm0    ! Mean long of perihelion at
                                         ! vernal equinox (radians)
    real(rk8) , intent(out)   :: mvelpp  ! moving vernal equinox long
                                         ! of perihelion plus pi (rad)
    logical , intent(in) :: log_print    ! Flags print of status/error

    ! # of elements in series wrt obliquity
    integer(ik4) , parameter :: poblen = 47
    ! # of elements in series wrt eccentricity
    integer(ik4) , parameter :: pecclen = 19
    ! # of elements in series wrt vernal equinox
    integer(ik4) , parameter :: pmvelen = 78
    ! arc sec to deg conversion
    real(rk8) , parameter :: psecdeg = 1.0D0/3600.0D0

    real(rk8) :: degrad = pi/180.D0 ! degree to radian conversion factor
    real(rk8) :: yb4_1950AD         ! number of years before 1950 AD

    character(len=*) , parameter :: subname = '(shr_orb_params)'

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    ! amplitudes for obliquity cos series
    real(rk8) , parameter :: obamp(poblen) =  &
          (/   -2462.2214466D0, -857.3232075D0, -629.3231835D0,   &
                -414.2804924D0, -311.7632587D0,  308.9408604D0,   &
                -162.5533601D0, -116.1077911D0,  101.1189923D0,   &
                 -67.6856209D0,   24.9079067D0,   22.5811241D0,   &
                 -21.1648355D0,  -15.6549876D0,   15.3936813D0,   &
                  14.6660938D0,  -11.7273029D0,   10.2742696D0,   &
                   6.4914588D0,    5.8539148D0,   -5.4872205D0,   &
                  -5.4290191D0,    5.1609570D0,    5.0786314D0,   &
                  -4.0735782D0,    3.7227167D0,    3.3971932D0,   &
                  -2.8347004D0,   -2.6550721D0,   -2.5717867D0,   &
                  -2.4712188D0,    2.4625410D0,    2.2464112D0,   &
                  -2.0755511D0,   -1.9713669D0,   -1.8813061D0,   &
                  -1.8468785D0,    1.8186742D0,    1.7601888D0,   &
                  -1.5428851D0,    1.4738838D0,   -1.4593669D0,   &
                   1.4192259D0,   -1.1818980D0,    1.1756474D0,   &
                  -1.1316126D0,    1.0896928D0/)

    ! rates for obliquity cosine series
    real(rk8) , parameter :: obrate(poblen) = &
            (/  31.609974D0, 32.620504D0, 24.172203D0,   &
                31.983787D0, 44.828336D0, 30.973257D0,   &
                43.668246D0, 32.246691D0, 30.599444D0,   &
                42.681324D0, 43.836462D0, 47.439436D0,   &
                63.219948D0, 64.230478D0,  1.010530D0,   &
                 7.437771D0, 55.782177D0,  0.373813D0,   &
                13.218362D0, 62.583231D0, 63.593761D0,   &
                76.438310D0, 45.815258D0,  8.448301D0,   &
                56.792707D0, 49.747842D0, 12.058272D0,   &
                75.278220D0, 65.241008D0, 64.604291D0,   &
                 1.647247D0,  7.811584D0, 12.207832D0,   &
                63.856665D0, 56.155990D0, 77.448840D0,   &
                 6.801054D0, 62.209418D0, 20.656133D0,   &
                48.344406D0, 55.145460D0, 69.000539D0,   &
                11.071350D0, 74.291298D0, 11.047742D0,   &
                 0.636717D0, 12.844549D0/)

    ! phases for obliquity cosine series
    real(rk8) , parameter :: obphas(poblen) = &
          (/    251.9025D0, 280.8325D0, 128.3057D0,   &
                292.7252D0,  15.3747D0, 263.7951D0,   &
                308.4258D0, 240.0099D0, 222.9725D0,   &
                268.7809D0, 316.7998D0, 319.6024D0,   &
                143.8050D0, 172.7351D0,  28.9300D0,   &
                123.5968D0,  20.2082D0,  40.8226D0,   &
                123.4722D0, 155.6977D0, 184.6277D0,   &
                267.2772D0,  55.0196D0, 152.5268D0,   &
                 49.1382D0, 204.6609D0,  56.5233D0,   &
                200.3284D0, 201.6651D0, 213.5577D0,   &
                 17.0374D0, 164.4194D0,  94.5422D0,   &
                131.9124D0,  61.0309D0, 296.2073D0,   &
                135.4894D0, 114.8750D0, 247.0691D0,   &
                256.6114D0,  32.1008D0, 143.6804D0,   &
                 16.8784D0, 160.6835D0,  27.5932D0,   &
                348.1074D0,  82.6496D0/)

    ! Cosine/sine series data for computation of eccentricity and fixed vernal
    ! equinox longitude of perihelion (fvelp): amplitude,
    ! rate (arc seconds/year), phase (degrees).

    ! ampl for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecamp (pecclen) = &
          (/   0.01860798D0,  0.01627522D0, -0.01300660D0,   &
               0.00988829D0, -0.00336700D0,  0.00333077D0,   &
              -0.00235400D0,  0.00140015D0,  0.00100700D0,   &
               0.00085700D0,  0.00064990D0,  0.00059900D0,   &
               0.00037800D0, -0.00033700D0,  0.00027600D0,   &
               0.00018200D0, -0.00017400D0, -0.00012400D0,   &
               0.00001250D0/)

    ! rates for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecrate(pecclen) = &
          (/    4.2072050D0,  7.3460910D0, 17.8572630D0,  &
               17.2205460D0, 16.8467330D0,  5.1990790D0,  &
               18.2310760D0, 26.2167580D0,  6.3591690D0,  &
               16.2100160D0,  3.0651810D0, 16.5838290D0,  &
               18.4939800D0,  6.1909530D0, 18.8677930D0,  &
               17.4255670D0,  6.1860010D0, 18.4174410D0,  &
                0.6678630D0/)

    ! phases for eccen/fvelp cos/sin series
    real(rk8) , parameter :: ecphas(pecclen) = &
          (/    28.620089D0, 193.788772D0, 308.307024D0,  &
               320.199637D0, 279.376984D0,  87.195000D0,  &
               349.129677D0, 128.443387D0, 154.143880D0,  &
               291.269597D0, 114.860583D0, 332.092251D0,  &
               296.414411D0, 145.769910D0, 337.237063D0,  &
               152.092288D0, 126.839891D0, 210.667199D0,  &
                72.108838D0/)

    ! Sine series data for computation of moving vernal equinox longitude of
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    ! amplitudes for mvelp sine series
    real(rk8) , parameter :: mvamp (pmvelen) = &
          (/   7391.0225890D0, 2555.1526947D0, 2022.7629188D0,  &
              -1973.6517951D0, 1240.2321818D0,  953.8679112D0,  &
               -931.7537108D0,  872.3795383D0,  606.3544732D0,  &
               -496.0274038D0,  456.9608039D0,  346.9462320D0,  &
               -305.8412902D0,  249.6173246D0, -199.1027200D0,  &
                191.0560889D0, -175.2936572D0,  165.9068833D0,  &
                161.1285917D0,  139.7878093D0, -133.5228399D0,  &
                117.0673811D0,  104.6907281D0,   95.3227476D0,  &
                 86.7824524D0,   86.0857729D0,   70.5893698D0,  &
                -69.9719343D0,  -62.5817473D0,   61.5450059D0,  &
                -57.9364011D0,   57.1899832D0,  -57.0236109D0,  &
                -54.2119253D0,   53.2834147D0,   52.1223575D0,  &
                -49.0059908D0,  -48.3118757D0,  -45.4191685D0,  &
                -42.2357920D0,  -34.7971099D0,   34.4623613D0,  &
                -33.8356643D0,   33.6689362D0,  -31.2521586D0,  &
                -30.8798701D0,   28.4640769D0,  -27.1960802D0,  &
                 27.0860736D0,  -26.3437456D0,   24.7253740D0,  &
                 24.6732126D0,   24.4272733D0,   24.0127327D0,  &
                 21.7150294D0,  -21.5375347D0,   18.1148363D0,  &
                -16.9603104D0,  -16.1765215D0,   15.5567653D0,  &
                 15.4846529D0,   15.2150632D0,   14.5047426D0,  &
                -14.3873316D0,   13.1351419D0,   12.8776311D0,  &
                 11.9867234D0,   11.9385578D0,   11.7030822D0,  &
                 11.6018181D0,  -11.2617293D0,  -10.4664199D0,  &
                 10.4333970D0,  -10.2377466D0,   10.1934446D0,  &
                -10.1280191D0,   10.0289441D0,  -10.0034259D0/)

    ! rates for mvelp sine series
    real(rk8) , parameter :: mvrate(pmvelen) = &
          (/    31.609974D0, 32.620504D0, 24.172203D0,   &
                 0.636717D0, 31.983787D0,  3.138886D0,   &
                30.973257D0, 44.828336D0,  0.991874D0,   &
                 0.373813D0, 43.668246D0, 32.246691D0,   &
                30.599444D0,  2.147012D0, 10.511172D0,   &
                42.681324D0, 13.650058D0,  0.986922D0,   &
                 9.874455D0, 13.013341D0,  0.262904D0,   &
                 0.004952D0,  1.142024D0, 63.219948D0,   &
                 0.205021D0,  2.151964D0, 64.230478D0,   &
                43.836462D0, 47.439436D0,  1.384343D0,   &
                 7.437771D0, 18.829299D0,  9.500642D0,   &
                 0.431696D0,  1.160090D0, 55.782177D0,   &
                12.639528D0,  1.155138D0,  0.168216D0,   &
                 1.647247D0, 10.884985D0,  5.610937D0,   &
                12.658184D0,  1.010530D0,  1.983748D0,   &
                14.023871D0,  0.560178D0,  1.273434D0,   &
                12.021467D0, 62.583231D0, 63.593761D0,   &
                76.438310D0,  4.280910D0, 13.218362D0,   &
                17.818769D0,  8.359495D0, 56.792707D0,   &
                8.448301D0,  1.978796D0,  8.863925D0,   &
                 0.186365D0,  8.996212D0,  6.771027D0,   &
                45.815258D0, 12.002811D0, 75.278220D0,   &
                65.241008D0, 18.870667D0, 22.009553D0,   &
                64.604291D0, 11.498094D0,  0.578834D0,   &
                 9.237738D0, 49.747842D0,  2.147012D0,   &
                 1.196895D0,  2.133898D0,  0.173168D0/)

    ! phases for mvelp sine series
    real(rk8) , parameter :: mvphas(pmvelen) = &
          (/    251.9025D0, 280.8325D0, 128.3057D0,   &
                348.1074D0, 292.7252D0, 165.1686D0,   &
                263.7951D0,  15.3747D0,  58.5749D0,   &
                 40.8226D0, 308.4258D0, 240.0099D0,   &
                222.9725D0, 106.5937D0, 114.5182D0,   &
                268.7809D0, 279.6869D0,  39.6448D0,   &
                126.4108D0, 291.5795D0, 307.2848D0,   &
                 18.9300D0, 273.7596D0, 143.8050D0,   &
                191.8927D0, 125.5237D0, 172.7351D0,   &
                316.7998D0, 319.6024D0,  69.7526D0,   &
                123.5968D0, 217.6432D0,  85.5882D0,   &
                156.2147D0,  66.9489D0,  20.2082D0,   &
                250.7568D0,  48.0188D0,   8.3739D0,   &
                 17.0374D0, 155.3409D0,  94.1709D0,   &
                221.1120D0,  28.9300D0, 117.1498D0,   &
                320.5095D0, 262.3602D0, 336.2148D0,   &
                233.0046D0, 155.6977D0, 184.6277D0,   &
                267.2772D0,  78.9281D0, 123.4722D0,   &
                188.7132D0, 180.1364D0,  49.1382D0,   &
                152.5268D0,  98.2198D0,  97.4808D0,   &
                221.5376D0, 168.2438D0, 161.1199D0,   &
                 55.0196D0, 262.6495D0, 200.3284D0,   &
                201.6651D0, 294.6547D0,  99.8233D0,   &
                213.5577D0, 154.1631D0, 232.7153D0,   &
                138.3034D0, 204.6609D0, 106.5938D0,   &
                250.4676D0, 332.3345D0,  27.3039D0/)

    integer(ik4) :: i       ! Index for series summations
    real(rk8) :: obsum   ! Obliquity series summation
    real(rk8) :: cossum  ! Cos series summation for eccentricity/fvelp
    real(rk8) :: sinsum  ! Sin series summation for eccentricity/fvelp
    real(rk8) :: fvelp   ! Fixed vernal equinox long of perihelion
    real(rk8) :: mvsum   ! mvelp series summation
    real(rk8) :: beta    ! Intermediate argument for lambm0
    real(rk8) :: years   ! Years to time of interest ( pos <=> future)
    real(rk8) :: eccen2  ! eccentricity squared
    real(rk8) :: eccen3  ! eccentricity cubed

    character(*) , parameter :: svnID  = "SVN " // &
         "$Id: shr_orb_mod.F90 25434 2010-11-04 22:46:24Z tcraig $"
    character(*) , parameter :: svnURL = "SVN <unknown URL>"
    character(len=*) , parameter :: F00 = "('(shr_orb_params) ',4a)"
    character(len=*) , parameter :: F01 = "('(shr_orb_params) ',a,i9)"
    character(len=*) , parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
    character(len=*) , parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

    ! radinp and algorithms below will need a degree to radian conversion factor

    if ( log_print .and. debug_level > 0 ) then
      write(stdout,F00) 'Calculate characteristics of the orbit:'
      write(stdout,F00) svnID
!     write(stdout,F00) svnURL
    end if

    ! Check for flag to use input orbit parameters

    if ( iyear_AD == SHR_ORB_UNDEF_INT ) then

      ! Check input obliq, eccen, and mvelp to ensure reasonable

      if ( obliq == SHR_ORB_UNDEF_REAL ) then
        write(stdout,F00) trim(subname)//' Have to specify orbital parameters:'
        write(stdout,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
        write(stdout,F00) 'iyear_AD is the year to simulate orbit (ie. 1950): '
        write(stdout,F00) 'obliq, eccen, mvelp specify the orbit directly:'
        write(stdout,F00) 'The AMIP II settings (for a 1995 orbit) are: '
        write(stdout,F00) ' obliq =  23.4441'
        write(stdout,F00) ' eccen =   0.016715'
        write(stdout,F00) ' mvelp = 102.7'
        call fatal(__FILE__,__LINE__,subname//' ERROR: unreasonable obliq')
      else if ( log_print ) then
        write(stdout,F00) 'Use input orbital parameters: '
      end if
      if ( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
        write(stdout,F03) 'Input obliquity unreasonable: ', obliq
        call fatal(__FILE__,__LINE__,subname//' ERROR: unreasonable obliq')
      end if
      if ( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
        write(stdout,F03) 'Input eccentricity unreasonable: ', eccen
        call fatal(__FILE__,__LINE__,subname//' ERROR: unreasonable eccen')
      end if
      if ( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
        write(stdout,F03) 'Input mvelp unreasonable: ' , mvelp
        call fatal(__FILE__,__LINE__,subname//' ERROR: unreasonable mvelp')
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

    else  ! Otherwise calculate based on years before present

      if ( log_print .and. debug_level > 0) then
        write(stdout,F01) 'Calculate orbit for year: ' , iyear_AD
      end if
      yb4_1950AD = 1950.0D0 - real(iyear_AD,rk8)
      if ( abs(yb4_1950AD) .gt. 1000000.0D0 )then
        write(stdout,F00) 'orbit only valid for years+-1000000'
        write(stdout,F00) 'Relative to 1950 AD'
        write(stdout,F03) '# of years before 1950: ',yb4_1950AD
        write(stdout,F01) 'Year to simulate was  : ',iyear_AD
        call fatal(__FILE__,__LINE__,subname//' ERROR: unreasonable year')
      end if

      ! The following calculates the earths obliquity, orbital eccentricity
      ! (and various powers of it) and vernal equinox mean longitude of
      ! perihelion for years in the past (future = negative of years past),
      ! using constants (see parameter section) given in the program of:
      !
      ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
      ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
      ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
      !
      ! and formulas given in the paper (where less precise constants are also
      ! given):
      !
      ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
      ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
      !
      ! The algorithm is valid only to 1,000,000 years past or hence.
      ! For a solution valid to 5-10 million years past see the above author.
      ! Algorithm below is better for years closer to present than is the
      ! 5-10 million year solution.
      !
      ! Years to time of interest must be negative of years before present
      ! (1950) in formulas that follow.

      years = - yb4_1950AD

      ! In the summations below, cosine or sine arguments, which end up in
      ! degrees, must be converted to radians via multiplication by degrad.
      !
      ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
      ! degrees. Convert the amplitudes and rates, which are in arc secs, into
      ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
      ! factor).  For obliq, first term is Berger 1978 epsilon star; second
      ! term is series summation in degrees.

      obsum = 0.0D0
      do i = 1 , poblen
        obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
                obphas(i))*degrad)
      end do
      obliq = 23.320556D0 + obsum

      ! Summation of cosine and sine series for computation of eccentricity
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation
      ! of moving vernal equinox longitude of perihelion.  Convert the rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.

      cossum = 0.0D0
      do i = 1 , pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      sinsum = 0.0D0
      do i = 1 , pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      ! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen

      ! A series of cases for fvelp, which is in radians.

      if (abs(cossum) .le. 1.0D-8) then
        if (sinsum .eq. 0.0D0) then
          fvelp = 0.0D0
        else if (sinsum .lt. 0.0D0) then
          fvelp = 1.5D0*pi
        else if (sinsum .gt. 0.0D0) then
          fvelp = .5D0*pi
        end if
      else if (cossum .lt. 0.0D0) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0D0) then
        if (sinsum .lt. 0.0D0) then
          fvelp = atan(sinsum/cossum) + 2.0D0*pi
        else
          fvelp = atan(sinsum/cossum)
        end if
      end if

      ! Summation of sin series for computation of moving vernal equinox long
      ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
      ! first term is fvelp in degrees; second term is Berger 1978 psi bar
      ! times years and in degrees; third term is Berger 1978 zeta; fourth
      ! term is series summation in degrees.  Convert the amplitudes and rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
      ! Series summation plus second and third terms constitute Berger 1978
      ! psi, which is the general precession.

      mvsum = 0.0D0
      do i = 1 , pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
                mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273D0*psecdeg*years + 3.392506D0 + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp .lt. 0.0D0)
        mvelp = mvelp + 360.0D0
      end do
      do while (mvelp .ge. 360.0D0)
        mvelp = mvelp - 360.0D0
      end do

    end if  ! end of test on whether to calculate or use input orbital params

    ! Orbit needs the obliquity in radians

    obliqr = obliq*degrad

    ! 180 degrees must be added to mvelp since observations are made from the
    ! earth and the sun is considered (wrongly for the algorithm) to go around
    ! the earth. For a more graphic explanation see Appendix B in:
    !
    ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
    ! Periods.  J. of Geophysical Research 98:10,341-10,362.
    !
    ! Additionally, orbit will need this value in radians. So mvelp becomes
    ! mvelpp (mvelp plus pi)

    mvelpp = (mvelp + 180.D0)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = sqrt(1.D0 - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = 2.D0*((.5D0*eccen + .125D0*eccen3)*(1.D0 + beta)*sin(mvelpp)  &
           - .250D0*eccen2*(.5D0    + beta)*sin(2.D0*mvelpp)               &
           + .125D0*eccen3*(1.D0/3.D0 + beta)*sin(3.D0*mvelpp))

    if ( log_print ) then
      write(stdout,F03) '------ Computed Orbital Parameters ------'
      write(stdout,F03) 'Eccentricity      = ',eccen
      write(stdout,F03) 'Obliquity (deg)   = ',obliq
      write(stdout,F03) 'Obliquity (rad)   = ',obliqr
      write(stdout,F03) 'Long of perh(deg) = ',mvelp
      write(stdout,F03) 'Long of perh(rad) = ',mvelpp
      write(stdout,F03) 'Long at v.e.(rad) = ',lambm0
      write(stdout,F03) '-----------------------------------------'
    end if
  end subroutine shr_orb_params
  !
  ! Compute earth/orbit parameters using formula suggested by
  ! Duane Thresher.
  !
  subroutine shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)
    real(rk8) , intent(in) :: calday ! Calendar day, including fraction
    real(rk8) , intent(in) :: eccen  ! Eccentricity
    real(rk8) , intent(in) :: obliqr ! Earths obliquity in radians
    real(rk8) , intent(in) :: lambm0 ! Mean long of perihelion at the
                                     ! vernal equinox (radians)
    real(rk8) , intent(in) :: mvelpp ! moving vernal equinox longitude
                                     ! of perihelion plus pi (radians)
    real(rk8) , intent(out) :: delta ! Solar declination angle in rad
    real(rk8) , intent(out) :: eccf  ! Earth-sun distance factor (ie. (1/r)**2)

    real(rk8) , parameter :: dayspy = 365.0D0  ! days per year
    real(rk8) , parameter :: ve     = 80.5D0   ! Calday of vernal equinox
                                              ! assumes Jan 1 = calday 1

    real(rk8) :: lambm  ! Lambda m, mean long of perihelion (rad)
    real(rk8) :: lmm    ! Intermediate argument involving lambm
    real(rk8) :: lamb   ! Lambda, the earths long of perihelion
    real(rk8) :: invrho ! Inverse normalized sun/earth distance
    real(rk8) :: sinl   ! Sine of lmm

    ! Compute eccentricity factor and solar declination using
    ! day value where a round day (such as 213.0) refers to 0z at
    ! Greenwich longitude.
    !
    ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    ! 35:2362-2367.
    !
    ! To get the earths true longitude (position in orbit; lambda in Berger
    ! 1978) which is necessary to find the eccentricity factor and declination,
    ! must first calculate the mean longitude (lambda m in Berger 1978) at
    ! the present day.  This is done by adding to lambm0 (the mean longitude
    ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    ! an increment (delta lambda m in Berger 1978) that is the number of
    ! days past or before (a negative increment) the vernal equinox divided by
    ! the days in a model year times the 2*pi radians in a complete orbit.

    lambm = lambm0 + (calday - ve)*2.D0*pi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(2.D0*sinl + eccen*(1.25D0*sin(2.D0*lmm)  &
          + eccen*((13.0D0/12.0D0)*sin(3.D0*lmm) - 0.25D0*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually
    ! inverse rho will be used), and thus the eccentricity factor (eccf),
    ! can be calculated from formulas given in Berger 1978.

    invrho = (1.D0 + eccen*cos(lamb - mvelpp)) / (1.D0 - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta  = asin(sin(obliqr)*sin(lamb))
    eccf   = invrho*invrho
  end subroutine shr_orb_decl

end module mod_clm_surfacealbedo
