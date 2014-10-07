module mod_clm_surfaceradiation
  !
  ! Calculate solar fluxes absorbed by vegetation and ground surface
  !
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_mpmessage
  use mod_runparams

  implicit none

  private

  save

  public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface

  contains
  !
  ! Solar fluxes absorbed by vegetation and ground surface
  ! Note possible problem when land is on different grid than atmosphere.
  ! Land may have sun above the horizon (coszen > 0) but atmosphere may
  ! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
  ! because all fluxes (absorbed, reflected, transmitted) are multiplied
  ! by the incoming flux and all will equal zero.
  ! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0)
  ! but land may have sun below horizon. This is okay because fabd, fabi,
  ! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
  ! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words,
  ! all the radiation is reflected. NDVI should equal zero in this case.
  ! However, the way the code is currently implemented this is only true
  ! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
  ! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
  !
  subroutine SurfaceRadiation(lbp, ubp, num_nourbanp, filter_nourbanp)
    use mod_clm_type
    use mod_clm_atmlnd      , only : clm_a2l
    use mod_clm_varpar      , only : numrad
    use mod_clm_varcon      , only : spval, istsoil, degpsec, isecspday
    use mod_clm_varcon      , only : istcrop
    use mod_clm_varctl      , only : subgridflag
    use mod_clm_varpar      , only : nlevsno
    use mod_clm_snicar      , only : DO_SNO_OC
    implicit none
    integer(ik4), intent(in) :: lbp, ubp   ! pft upper and lower bounds
    ! number of pfts in non-urban points in pft filter
    integer(ik4), intent(in) :: num_nourbanp
    ! pft filter for non-urban points
    integer(ik4), intent(in) :: filter_nourbanp(ubp-lbp+1)
    real(rk8), pointer :: albsod(:,:) ! direct-beam soil albedo (col,bnd) [frc]
    real(rk8), pointer :: albsoi(:,:) ! diffuse soil albedo (col,bnd) [frc]
    ! solar radiation absorbed by soil (W/m**2)
    real(rk8), pointer :: sabg_soil(:)
    ! solar radiation absorbed by snow (W/m**2)
    real(rk8), pointer :: sabg_snow(:)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: ivt(:)           ! pft vegetation type
    integer(ik4) , pointer :: pcolumn(:)       ! pft's column index
    integer(ik4) , pointer :: pgridcell(:)     ! pft's gridcell index
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer :: elai(:)
    ! one-sided stem area index with burying by snow
    real(rk8), pointer :: esai(:)
    real(rk8), pointer :: londeg(:)        ! longitude (degrees)
    real(rk8), pointer :: latdeg(:)        ! latitude (degrees)
    real(rk8), pointer :: coszen(:)        ! cosine of solar zenith angle
    real(rk8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
    real(rk8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
    ! flux absorbed by canopy per unit direct flux
    real(rk8), pointer :: fabd(:,:)
    ! flux absorbed by sunlit canopy per unit direct flux
    real(rk8), pointer :: fabd_sun(:,:)
    ! flux absorbed by shaded canopy per unit direct flux
    real(rk8), pointer :: fabd_sha(:,:)
    ! flux absorbed by canopy per unit diffuse flux
    real(rk8), pointer :: fabi(:,:)
    ! flux absorbed by sunlit canopy per unit diffuse flux
    real(rk8), pointer :: fabi_sun(:,:)
    ! flux absorbed by shaded canopy per unit diffuse flux
    real(rk8), pointer :: fabi_sha(:,:)
    ! down direct flux below canopy per unit direct flux
    real(rk8), pointer :: ftdd(:,:)
    ! down diffuse flux below canopy per unit direct flux
    real(rk8), pointer :: ftid(:,:)
    ! down diffuse flux below canopy per unit diffuse flux
    real(rk8), pointer :: ftii(:,:)
    ! number of canopy layers, above snow for radiative transfer
    integer(ik4) , pointer :: nrad(:)
    ! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer :: fabd_sun_z(:,:)
    ! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer :: fabd_sha_z(:,:)
    ! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer :: fabi_sun_z(:,:)
    ! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(rk8), pointer :: fabi_sha_z(:,:)
    real(rk8), pointer :: fsun_z(:,:)  ! sunlit fraction of canopy layer
    real(rk8), pointer :: tlai_z(:,:)  ! tlai increment for canopy layer
    real(rk8), pointer :: tsai_z(:,:)  ! tsai increment for canopy layer
    real(rk8), pointer :: albgrd(:,:)  ! ground albedo (direct)
    real(rk8), pointer :: albgri(:,:)  ! ground albedo (diffuse)
    real(rk8), pointer :: albd(:,:)    ! surface albedo (direct)
    real(rk8), pointer :: albi(:,:)    ! surface albedo (diffuse)

    real(rk8), pointer :: fsun(:)      ! sunlit fraction of canopy
    real(rk8), pointer :: laisun(:)    ! sunlit leaf area
    real(rk8), pointer :: laisha(:)    ! shaded leaf area
    ! sunlit leaf area for canopy layer
    real(rk8), pointer :: laisun_z(:,:)
    ! shaded leaf area for canopy layer
    real(rk8), pointer :: laisha_z(:,:)
    ! absorbed PAR for sunlit leaves in canopy layer
    real(rk8), pointer :: parsun_z(:,:)
    ! absorbed PAR for shaded leaves in canopy layer
    real(rk8), pointer :: parsha_z(:,:)
    ! solar radiation absorbed by ground (W/m**2)
    real(rk8), pointer :: sabg(:)
    ! solar radiation absorbed by vegetation (W/m**2)
    real(rk8), pointer :: sabv(:)
    ! solar radiation absorbed (total) (W/m**2)
    real(rk8), pointer :: fsa(:)
    ! rural solar radiation absorbed (total) (W/m**2)
    real(rk8), pointer :: fsa_r(:)
    ! landunit type
    integer(ik4) , pointer :: ityplun(:)
    ! index into landunit level quantities
    integer(ik4) , pointer :: plandunit(:)
    ! solar radiation reflected (W/m**2)
    real(rk8), pointer :: fsr(:)
    ! incident direct beam vis solar radiation (W/m**2)
    real(rk8), pointer :: fsds_vis_d(:)
    ! incident direct beam nir solar radiation (W/m**2)
    real(rk8), pointer :: fsds_nir_d(:)
    ! incident diffuse vis solar radiation (W/m**2)
    real(rk8), pointer :: fsds_vis_i(:)
    ! incident diffuse nir solar radiation (W/m**2)
    real(rk8), pointer :: fsds_nir_i(:)
    ! reflected direct beam vis solar radiation (W/m**2)
    real(rk8), pointer :: fsr_vis_d(:)
    ! reflected direct beam nir solar radiation (W/m**2)
    real(rk8), pointer :: fsr_nir_d(:)
    ! reflected diffuse vis solar radiation (W/m**2)
    real(rk8), pointer :: fsr_vis_i(:)
    ! reflected diffuse nir solar radiation (W/m**2)
    real(rk8), pointer :: fsr_nir_i(:)
    ! incident direct beam vis solar rad at local noon (W/m**2)
    real(rk8), pointer :: fsds_vis_d_ln(:)
    ! incident direct beam nir solar rad at local noon (W/m**2)
    real(rk8), pointer :: fsds_nir_d_ln(:)
    ! reflected direct beam vis solar rad at local noon (W/m**2)
    real(rk8), pointer :: fsr_vis_d_ln(:)
    ! reflected direct beam nir solar rad at local noon (W/m**2)
    real(rk8), pointer :: fsr_nir_d_ln(:)
    ! incident diffuse beam vis solar rad at local noon (W/m**2)
    real(rk8), pointer :: fsds_vis_i_ln(:)
    ! absorbed par by vegetation at local noon (W/m**2)
    real(rk8), pointer :: parveg_ln(:)
    ! direct flux absorption factor (col,lyr): VIS [frc]
    real(rk8), pointer :: flx_absdv(:,:)
    ! direct flux absorption factor (col,lyr): NIR [frc]
    real(rk8), pointer :: flx_absdn(:,:)
    ! diffuse flux absorption factor (col,lyr): VIS [frc]
    real(rk8), pointer :: flx_absiv(:,:)
    ! diffuse flux absorption factor (col,lyr): NIR [frc]
    real(rk8), pointer :: flx_absin(:,:)
    ! negative number of snow layers [nbr]
    integer(ik4) , pointer :: snl(:)
    ! pure snow ground albedo (direct)
    real(rk8), pointer :: albgrd_pur(:,:)
    ! pure snow ground albedo (diffuse)
    real(rk8), pointer :: albgri_pur(:,:)
    ! ground albedo without BC (direct) (col,bnd)
    real(rk8), pointer :: albgrd_bc(:,:)
    ! ground albedo without BC (diffuse) (col,bnd)
    real(rk8), pointer :: albgri_bc(:,:)
    ! ground albedo without OC (direct) (col,bnd)
    real(rk8), pointer :: albgrd_oc(:,:)
    ! ground albedo without OC (diffuse) (col,bnd)
    real(rk8), pointer :: albgri_oc(:,:)
    ! ground albedo without dust (direct) (col,bnd)
    real(rk8), pointer :: albgrd_dst(:,:)
    ! ground albedo without dust (diffuse) (col,bnd)
    real(rk8), pointer :: albgri_dst(:,:)
    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(rk8), pointer :: albsnd_hst(:,:)
    ! snow ground albedo, diffuse, for history files (col,bnd
    real(rk8), pointer :: albsni_hst(:,:)
    ! absorbed radiative flux (pft,lyr) [W/m2]
    real(rk8), pointer :: sabg_lyr(:,:)
    ! (rural) shortwave radiation penetrating top soisno layer [W/m2]
    real(rk8), pointer :: sabg_pen(:)
    ! surface forcing of snow with all aerosols (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_aer(:)
    ! surface forcing of snow with BC (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_bc(:)
    ! surface forcing of snow with OC (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_oc(:)
    ! surface forcing of snow with dust (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_dst(:)
    ! surface forcing of snow with all aerosols, averaged only
    ! when snow is present (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_aer_sno(:)
    ! surface forcing of snow with BC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_bc_sno(:)
    ! surface forcing of snow with OC, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_oc_sno(:)
    ! surface forcing of snow with dust, averaged only when snow
    ! is present (pft) [W/m2]
    real(rk8), pointer :: sfc_frc_dst_sno(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno(:)
    ! reflected visible, direct radiation from snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsr_sno_vd(:)
    ! reflected near-IR, direct radiation from snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsr_sno_nd(:)
    ! reflected visible, diffuse radiation from snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsr_sno_vi(:)
    ! reflected near-IR, diffuse radiation from snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsr_sno_ni(:)
    ! incident visible, direct radiation on snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsds_sno_vd(:)
    ! incident near-IR, direct radiation on snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsds_sno_nd(:)
    ! incident visible, diffuse radiation on snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsds_sno_vi(:)
    ! incident near-IR, diffuse radiation on snow
    ! (for history files) (pft) [W/m2]
    real(rk8), pointer :: fsds_sno_ni(:)
    ! snow height (m)
    real(rk8), pointer :: snow_depth(:)

    ! number of solar radiation waveband classes
    integer(ik4) , parameter :: nband = numrad

    integer(ik4)  :: fp         ! non-urban filter pft index
    integer(ik4)  :: p          ! pft index
    integer(ik4)  :: c          ! column index
    integer(ik4)  :: l          ! landunit index
    integer(ik4)  :: g          ! grid cell index
    integer(ik4)  :: ib         ! waveband number (1=vis, 2=nir)
    integer(ik4)  :: iv         ! canopy layer
    real(rk8) :: absrad    ! absorbed solar radiation (W/m**2)
    real(rk8) :: rnir      ! reflected solar radiation [nir] (W/m**2)
    real(rk8) :: rvis      ! reflected solar radiation [vis] (W/m**2)
    ! transmitted solar radiation: direct (W/m**2)
    real(rk8) :: trd(lbp:ubp,numrad)
    ! transmitted solar radiation: diffuse (W/m**2)
    real(rk8) :: tri(lbp:ubp,numrad)
    ! direct beam absorbed by canopy (W/m**2)
    real(rk8) :: cad(lbp:ubp,numrad)
    ! diffuse radiation absorbed by canopy (W/m**2)
    real(rk8) :: cai(lbp:ubp,numrad)
    ! seconds into current date in local time
    integer(ik4)  :: local_secp1
    ! land model time step (sec)
    real(rk8) :: dtime
    ! calendar info for current time step
    integer(ik4)  :: year,month,day,secs
    integer(ik4)  :: i   ! layer index [idx]
    ! temporary, absorbed energy in all active snow layers [W/m2]
    real(rk8) :: sabg_snl_sum
#ifdef SNICAR_FRC
    ! temp: absorbed solar radiation by pure snow [W/m2]
    real(rk8) :: absrad_pur
    ! temp: absorbed solar radiation without BC [W/m2]
    real(rk8) :: absrad_bc
    ! temp: absorbed solar radiation without OC [W/m2]
    real(rk8) :: absrad_oc
    ! temp: absorbed solar radiation without dust [W/m2]
    real(rk8) :: absrad_dst
#endif
    ! solar radiation absorbed by ground with pure snow [W/m2]
    real(rk8) :: sabg_pur(lbp:ubp)
    ! solar radiation absorbed by ground without BC [W/m2]
    real(rk8) :: sabg_bc(lbp:ubp)
    ! solar radiation absorbed by ground without OC [W/m2]
    real(rk8) :: sabg_oc(lbp:ubp)
    ! solar radiation absorbed by ground without dust [W/m2]
    real(rk8) :: sabg_dst(lbp:ubp)
    real(rk8) :: parveg(lbp:ubp)  ! absorbed par by vegetation (W/m**2)

    ! Assign local pointers to multi-level derived type members (gridcell level)
    londeg        => clm3%g%londeg
    latdeg        => clm3%g%latdeg
    forc_solad    => clm_a2l%forc_solad
    forc_solai    => clm_a2l%forc_solai

    ! Assign local pointers to multi-level derived type members (landunit level)

    ityplun       => clm3%g%l%itype

    ! Assign local pointers to multi-level derived type members (column level)

    albgrd        => clm3%g%l%c%cps%albgrd
    albgri        => clm3%g%l%c%cps%albgri
    coszen        => clm3%g%l%c%cps%coszen

    ! Assign local pointers to derived type members (pft-level)

    albsod        => clm3%g%l%c%cps%albsod
    albsoi        => clm3%g%l%c%cps%albsoi
    sabg_soil     => clm3%g%l%c%p%pef%sabg_soil
    sabg_snow     => clm3%g%l%c%p%pef%sabg_snow
    pactive       => clm3%g%l%c%p%active
    plandunit     => clm3%g%l%c%p%landunit
    ivt           => clm3%g%l%c%p%itype
    pcolumn       => clm3%g%l%c%p%column
    pgridcell     => clm3%g%l%c%p%gridcell
    elai          => clm3%g%l%c%p%pps%elai
    esai          => clm3%g%l%c%p%pps%esai
    laisun        => clm3%g%l%c%p%pps%laisun
    laisha        => clm3%g%l%c%p%pps%laisha
    laisun_z      => clm3%g%l%c%p%pps%laisun_z
    laisha_z      => clm3%g%l%c%p%pps%laisha_z
    albd          => clm3%g%l%c%p%pps%albd
    albi          => clm3%g%l%c%p%pps%albi
    fabd          => clm3%g%l%c%p%pps%fabd
    fabd_sun      => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha      => clm3%g%l%c%p%pps%fabd_sha
    fabi          => clm3%g%l%c%p%pps%fabi
    fabi_sun      => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha      => clm3%g%l%c%p%pps%fabi_sha
    ftdd          => clm3%g%l%c%p%pps%ftdd
    ftid          => clm3%g%l%c%p%pps%ftid
    ftii          => clm3%g%l%c%p%pps%ftii
    nrad          => clm3%g%l%c%p%pps%nrad
    fabd_sun_z    => clm3%g%l%c%p%pps%fabd_sun_z
    fabd_sha_z    => clm3%g%l%c%p%pps%fabd_sha_z
    fabi_sun_z    => clm3%g%l%c%p%pps%fabi_sun_z
    fabi_sha_z    => clm3%g%l%c%p%pps%fabi_sha_z
    fsun_z        => clm3%g%l%c%p%pps%fsun_z
    tlai_z        => clm3%g%l%c%p%pps%tlai_z
    tsai_z        => clm3%g%l%c%p%pps%tsai_z
    fsun          => clm3%g%l%c%p%pps%fsun
    sabg          => clm3%g%l%c%p%pef%sabg
    sabv          => clm3%g%l%c%p%pef%sabv
    snow_depth    => clm3%g%l%c%cps%snow_depth
    fsa           => clm3%g%l%c%p%pef%fsa
    fsa_r         => clm3%g%l%c%p%pef%fsa_r
    fsr           => clm3%g%l%c%p%pef%fsr
    parsun_z      => clm3%g%l%c%p%pef%parsun_z
    parsha_z      => clm3%g%l%c%p%pef%parsha_z
    fsds_vis_d    => clm3%g%l%c%p%pef%fsds_vis_d
    fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
    fsds_vis_i    => clm3%g%l%c%p%pef%fsds_vis_i
    fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
    fsr_vis_d     => clm3%g%l%c%p%pef%fsr_vis_d
    fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
    fsr_vis_i     => clm3%g%l%c%p%pef%fsr_vis_i
    fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
    fsds_vis_d_ln => clm3%g%l%c%p%pef%fsds_vis_d_ln
    fsds_nir_d_ln => clm3%g%l%c%p%pef%fsds_nir_d_ln
    parveg_ln     => clm3%g%l%c%p%pef%parveg_ln
    fsds_vis_i_ln => clm3%g%l%c%p%pef%fsds_vis_i_ln
    fsr_vis_d_ln  => clm3%g%l%c%p%pef%fsr_vis_d_ln
    fsr_nir_d_ln  => clm3%g%l%c%p%pef%fsr_nir_d_ln

    ! Assign local pointers to derived type members (ecophysiological)

    frac_sno         => clm3%g%l%c%cps%frac_sno
    flx_absdv        => clm3%g%l%c%cps%flx_absdv
    flx_absdn        => clm3%g%l%c%cps%flx_absdn
    flx_absiv        => clm3%g%l%c%cps%flx_absiv
    flx_absin        => clm3%g%l%c%cps%flx_absin
    sabg_lyr         => clm3%g%l%c%p%pef%sabg_lyr
    sabg_pen         => clm3%g%l%c%p%pef%sabg_pen
    snl              => clm3%g%l%c%cps%snl
    sfc_frc_aer      => clm3%g%l%c%p%pef%sfc_frc_aer
    sfc_frc_aer_sno  => clm3%g%l%c%p%pef%sfc_frc_aer_sno
    albgrd_pur       => clm3%g%l%c%cps%albgrd_pur
    albgri_pur       => clm3%g%l%c%cps%albgri_pur
    sfc_frc_bc       => clm3%g%l%c%p%pef%sfc_frc_bc
    sfc_frc_bc_sno   => clm3%g%l%c%p%pef%sfc_frc_bc_sno
    albgrd_bc        => clm3%g%l%c%cps%albgrd_bc
    albgri_bc        => clm3%g%l%c%cps%albgri_bc
    sfc_frc_oc       => clm3%g%l%c%p%pef%sfc_frc_oc
    sfc_frc_oc_sno   => clm3%g%l%c%p%pef%sfc_frc_oc_sno
    albgrd_oc        => clm3%g%l%c%cps%albgrd_oc
    albgri_oc        => clm3%g%l%c%cps%albgri_oc
    sfc_frc_dst      => clm3%g%l%c%p%pef%sfc_frc_dst
    sfc_frc_dst_sno  => clm3%g%l%c%p%pef%sfc_frc_dst_sno
    albgrd_dst       => clm3%g%l%c%cps%albgrd_dst
    albgri_dst       => clm3%g%l%c%cps%albgri_dst
    albsnd_hst       => clm3%g%l%c%cps%albsnd_hst
    albsni_hst       => clm3%g%l%c%cps%albsni_hst
    fsr_sno_vd       => clm3%g%l%c%p%pef%fsr_sno_vd
    fsr_sno_nd       => clm3%g%l%c%p%pef%fsr_sno_nd
    fsr_sno_vi       => clm3%g%l%c%p%pef%fsr_sno_vi
    fsr_sno_ni       => clm3%g%l%c%p%pef%fsr_sno_ni
    fsds_sno_vd      => clm3%g%l%c%p%pef%fsds_sno_vd
    fsds_sno_nd      => clm3%g%l%c%p%pef%fsds_sno_nd
    fsds_sno_vi      => clm3%g%l%c%p%pef%fsds_sno_vi
    fsds_sno_ni      => clm3%g%l%c%p%pef%fsds_sno_ni

    ! Determine seconds off current time step

    dtime = int(dtsrf)
    call curr_date(idatex,year,month,day,secs)

    ! Initialize fluxes

    do fp = 1,num_nourbanp
      p = filter_nourbanp(fp)
      ! was redundant b/c filter already included wt>0;
      ! not redundant anymore with chg in filter definition
      if (pactive(p)) then
        sabg_soil(p)  = 0.D0
        sabg_snow(p)  = 0.D0
        sabg(p)       = 0.D0
        sabv(p)       = 0.D0
        fsa(p)        = 0.D0
        l = plandunit(p)
        if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
          fsa_r(p)      = 0.D0
        end if
        sabg_lyr(p,:) = 0.D0
        sabg_pur(p)   = 0.D0
        sabg_bc(p)    = 0.D0
        sabg_oc(p)    = 0.D0
        sabg_dst(p)   = 0.D0
        do iv = 1, nrad(p)
          parsun_z(p,iv) = 0.D0
          parsha_z(p,iv) = 0.D0
          laisun_z(p,iv) = 0.D0
          laisha_z(p,iv) = 0.D0
        end do
      end if
    end do

    ! Loop over pfts to calculate laisun_z and laisha_z for each layer.
    ! Derive canopy laisun, laisha, and fsun from layer sums.
    ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
    ! SurfaceAlbedo is canopy integrated so that layer value equals
    ! canopy value.

    do fp = 1,num_nourbanp
      p = filter_nourbanp(fp)
      if (pactive(p)) then

        laisun(p) = 0.D0
        laisha(p) = 0.D0
        do iv = 1, nrad(p)
          laisun_z(p,iv) = tlai_z(p,iv) * fsun_z(p,iv)
          laisha_z(p,iv) = tlai_z(p,iv) * (1.D0 - fsun_z(p,iv))
          laisun(p) = laisun(p) + laisun_z(p,iv)
          laisha(p) = laisha(p) + laisha_z(p,iv)
        end do
        if (elai(p) > 0.D0) then
          fsun(p) = laisun(p) / elai(p)
        else
          fsun(p) = 0.D0
        end if

      end if
    end do

    ! Loop over nband wavebands
    do ib = 1, nband
      do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        if (pactive(p)) then
          c = pcolumn(p)
          l = plandunit(p)
          g = pgridcell(p)

          ! Absorbed by canopy

          cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
          cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
          sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
          fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
          if (ib == 1) then
            parveg(p) = cad(p,ib) + cai(p,ib)
          end if
          if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
            fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
          end if

          ! Absorbed PAR profile through canopy
          ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
          ! are canopy integrated so that layer values equal big leaf values.

          if (ib == 1) then
            do iv = 1, nrad(p)
              parsun_z(p,iv) = forc_solad(g,ib)*fabd_sun_z(p,iv) + &
                       forc_solai(g,ib)*fabi_sun_z(p,iv)
              parsha_z(p,iv) = forc_solad(g,ib)*fabd_sha_z(p,iv) + &
                       forc_solai(g,ib)*fabi_sha_z(p,iv)
            end do
          end if

          ! Transmitted = solar fluxes incident on ground

          trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
          tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)

          ! Solar radiation absorbed by ground surface

          ! calculate absorbed solar by soil/snow separately
          absrad  = trd(p,ib)*(1.D0-albsod(c,ib)) + &
                  tri(p,ib)*(1.D0-albsoi(c,ib))
          sabg_soil(p) = sabg_soil(p) + absrad
          absrad  = trd(p,ib)*(1.D0-albsnd_hst(c,ib)) + &
                  tri(p,ib)*(1.D0-albsni_hst(c,ib))
          sabg_snow(p) = sabg_snow(p) + absrad
          absrad  = trd(p,ib)*(1.D0-albgrd(c,ib)) + &
                  tri(p,ib)*(1.D0-albgri(c,ib))
          sabg(p) = sabg(p) + absrad
          fsa(p)  = fsa(p)  + absrad
          if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
            fsa_r(p)  = fsa_r(p)  + absrad
          end if
          if (snl(c) == 0) then
            sabg_snow(p) = sabg(p)
            sabg_soil(p) = sabg(p)
          end if
          ! if no subgrid fluxes, make sure to set both components
          ! equal to weighted average
          if (subgridflag == 0) then
            sabg_snow(p) = sabg(p)
            sabg_soil(p) = sabg(p)
          end if

#if (defined SNICAR_FRC)
          ! Solar radiation absorbed by ground surface without BC
          absrad_bc = trd(p,ib)*(1.D0-albgrd_bc(c,ib)) + &
                   tri(p,ib)*(1.D0-albgri_bc(c,ib))
          sabg_bc(p) = sabg_bc(p) + absrad_bc

          ! Solar radiation absorbed by ground surface without OC
          absrad_oc = trd(p,ib)*(1.D0-albgrd_oc(c,ib)) + &
                   tri(p,ib)*(1.D0-albgri_oc(c,ib))
          sabg_oc(p) = sabg_oc(p) + absrad_oc

          ! Solar radiation absorbed by ground surface without dust
          absrad_dst = trd(p,ib)*(1.D0-albgrd_dst(c,ib)) + &
                   tri(p,ib)*(1.D0-albgri_dst(c,ib))
           sabg_dst(p) = sabg_dst(p) + absrad_dst

          ! Solar radiation absorbed by ground surface without any aerosols
          absrad_pur = trd(p,ib)*(1.D0-albgrd_pur(c,ib)) + &
                   tri(p,ib)*(1.D0-albgri_pur(c,ib))
          sabg_pur(p) = sabg_pur(p) + absrad_pur
#endif

        end if
      end do ! end of pft loop
    end do ! end nbands loop

    ! compute absorbed flux in each snow layer and top soil layer,
    ! based on flux factors computed in the radiative transfer portion
    ! of SNICAR.

    do fp = 1,num_nourbanp
      p = filter_nourbanp(fp)
      if (pactive(p)) then
        c = pcolumn(p)
        l = plandunit(p)
        sabg_snl_sum = 0.D0

        ! CASE1: No snow layers: all energy is absorbed in top soil layer
        if (snl(c) == 0) then
          sabg_lyr(p,:) = 0.D0
          sabg_lyr(p,1) = sabg(p)
          sabg_snl_sum  = sabg_lyr(p,1)

          ! CASE 2: Snow layers present: absorbed radiation is scaled
          ! according to flux factors computed by SNICAR
        else
          do i = -nlevsno+1,1,1
            sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + &
                            flx_absdn(c,i)*trd(p,2) + &
                            flx_absiv(c,i)*tri(p,1) + &
                            flx_absin(c,i)*tri(p,2)
            ! summed radiation in active snow layers:
            if (i >= snl(c)+1) then
              sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
            end if
          end do

          ! Error handling: The situation below can occur when solar
          ! radiation is NOT computed every timestep.
          ! When the number of snow layers has changed in between
          ! computations of the absorbed solar energy in each layer,
          ! we must redistribute the absorbed energy to avoid physically
          ! unrealistic conditions. The assumptions made below are
          ! somewhat arbitrary, but this situation does not arise very
          ! frequently.
          ! This error handling is implemented to accomodate any value of the
          ! radiation frequency.
          ! change condition to match sabg_snow isntead of sabg
          if (abs(sabg_snl_sum-sabg_snow(p)) > 0.00001D0) then
            if (snl(c) == 0) then
              sabg_lyr(p,-4:0) = 0.D0
              sabg_lyr(p,1) = sabg(p)
            else if (snl(c) == -1) then
              sabg_lyr(p,-4:-1) = 0.D0
              sabg_lyr(p,0) = sabg_snow(p)*0.6D0
              sabg_lyr(p,1) = sabg_snow(p)*0.4D0
            else
              sabg_lyr(p,:) = 0.D0
              sabg_lyr(p,snl(c)+1) = sabg_snow(p)*0.75D0
              sabg_lyr(p,snl(c)+2) = sabg_snow(p)*0.25D0
            end if
          end if

          ! If shallow snow depth, all solar radiation absorbed in top or
          ! top two snow layers to prevent unrealistic timestep soil warming
          if (subgridflag == 0) then
            if (snow_depth(c) < 0.10D0) then
              if (snl(c) == 0) then
                sabg_lyr(p,-4:0) = 0.D0
                sabg_lyr(p,1) = sabg(p)
              else if (snl(c) == -1) then
                sabg_lyr(p,-4:-1) = 0.D0
                sabg_lyr(p,0) = sabg(p)
                sabg_lyr(p,1) = 0.D0
              else
                sabg_lyr(p,:) = 0.D0
                sabg_lyr(p,snl(c)+1) = sabg(p)*0.75D0
                sabg_lyr(p,snl(c)+2) = sabg(p)*0.25D0
              end if
            end if
          end if
        end if

        ! This situation should not happen:
        if (abs(sum(sabg_lyr(p,:))-sabg_snow(p)) > 0.00001D0) then
          write(stderr,*) "SNICAR ERROR: Absorbed ground radiation "// &
                "not equal to summed snow layer radiation. pft = ",   &
                p," Col= ", c, " Diff= ", &
                sum(sabg_lyr(p,:))-sabg_snow(p), &
                " sabg_snow(p)= ", sabg_snow(p), &
                " sabg_sum(p)= ", &
                sum(sabg_lyr(p,:)), " snl(c)= ", snl(c)
          write(stderr,*) "flx_absdv1= ", trd(p,1)*(1.-albgrd(c,1)), &
                   "flx_absdv2= ", sum(flx_absdv(c,:))*trd(p,1)
          write(stderr,*) "flx_absiv1= ", tri(p,1)*(1.-albgri(c,1)), &
                   " flx_absiv2= ", sum(flx_absiv(c,:))*tri(p,1)
          write(stderr,*) "flx_absdn1= ", trd(p,2)*(1.-albgrd(c,2)), &
                   " flx_absdn2= ", sum(flx_absdn(c,:))*trd(p,2)
          write(stderr,*) "flx_absin1= ", tri(p,2)*(1.-albgri(c,2)), &
                   " flx_absin2= ", sum(flx_absin(c,:))*tri(p,2)
          write(stderr,*) "albgrd_nir= ", albgrd(c,2)
          write(stderr,*) "coszen= ", coszen(c)
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if

        ! Diagnostic: shortwave penetrating ground (e.g. top layer)
        if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
          sabg_pen(p) = sabg(p) - sabg_lyr(p, snl(c)+1)
        end if

#if (defined SNICAR_FRC)

        ! BC aerosol forcing (pft-level):
        sfc_frc_bc(p) = sabg(p) - sabg_bc(p)

        ! OC aerosol forcing (pft-level):
        if (DO_SNO_OC) then
          sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
        else
          sfc_frc_oc(p) = 0.D0
        end if

        ! dust aerosol forcing (pft-level):
        sfc_frc_dst(p) = sabg(p) - sabg_dst(p)

        ! all-aerosol forcing (pft-level):
        sfc_frc_aer(p) = sabg(p) - sabg_pur(p)

        ! forcings averaged only over snow:
        if (frac_sno(c) > 0.D0) then
          sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
          sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
          sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
          sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
        else
          sfc_frc_bc_sno(p)  = spval
          sfc_frc_oc_sno(p)  = spval
          sfc_frc_dst_sno(p) = spval
          sfc_frc_aer_sno(p) = spval
        end if
#endif
      end if
    end do

    ! Radiation diagnostics
    do fp = 1 , num_nourbanp
      p = filter_nourbanp(fp)
      if (pactive(p)) then
        g = pgridcell(p)

        ! NDVI and reflected solar radiation

        rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
        rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
        fsr(p) = rvis + rnir

        fsds_vis_d(p) = forc_solad(g,1)
        fsds_nir_d(p) = forc_solad(g,2)
        fsds_vis_i(p) = forc_solai(g,1)
        fsds_nir_i(p) = forc_solai(g,2)
        fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
        fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
        fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
        fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)

        local_secp1 = secs + nint((londeg(g)/degpsec)/dtime)*int(dtime)
        local_secp1 = mod(local_secp1,isecspday)
        if (local_secp1 == isecspday/2) then
          fsds_vis_d_ln(p) = forc_solad(g,1)
          fsds_nir_d_ln(p) = forc_solad(g,2)
          fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
          fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
          fsds_vis_i_ln(p) = forc_solai(g,1)
          parveg_ln(p)     = parveg(p)
        else
          fsds_vis_d_ln(p) = spval
          fsds_nir_d_ln(p) = spval
          fsr_vis_d_ln(p) = spval
          fsr_nir_d_ln(p) = spval
          fsds_vis_i_ln(p) = spval
          parveg_ln(p)     = spval
        end if

        ! diagnostic variables (downwelling and absorbed radiation
        ! partitioning) for history files (OPTIONAL)
        c = pcolumn(p)
        if (snl(c) < 0) then
          fsds_sno_vd(p) = forc_solad(g,1)
          fsds_sno_nd(p) = forc_solad(g,2)
          fsds_sno_vi(p) = forc_solai(g,1)
          fsds_sno_ni(p) = forc_solai(g,2)

          fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
          fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
          fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
          fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
        else
          fsds_sno_vd(p) = spval
          fsds_sno_nd(p) = spval
          fsds_sno_vi(p) = spval
          fsds_sno_ni(p) = spval

          fsr_sno_vd(p) = spval
          fsr_sno_nd(p) = spval
          fsr_sno_vi(p) = spval
          fsr_sno_ni(p) = spval
        end if
      end if
    end do
  end subroutine SurfaceRadiation

end module mod_clm_surfaceradiation
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
