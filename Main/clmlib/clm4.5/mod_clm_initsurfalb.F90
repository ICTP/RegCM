module mod_clm_initsurfalb
  !
  ! Computes initial surface albedo calculation -
  ! Initialization of ecosystem dynamics is needed for this
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants , only : degrad
  use mod_clm_type
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_filter , only : filter
  use mod_clm_varpar , only : nlevsoi , nlevsno , nlevlak , nlevgrnd
  use mod_clm_varcon , only : zlnd , istsoil , isturb , denice , denh2o , &
            icol_roof , icol_road_imperv , icol_road_perv , istcrop
  use mod_clm_fracwet , only : FracWet
  use mod_clm_surfacealbedo , only : SurfaceAlbedo
#if (defined CN)
  use mod_stdio
  use mod_mpmessage
  use mod_clm_cnecosystemdyn , only : CNEcosystemDyn
  use mod_clm_cnvegstructupdate , only : CNVegStructUpdate
#else
  use mod_clm_staticecosysdyn , only : EcosystemDyn , interpMonthlyVeg
#endif
  use mod_clm_urban , only : UrbanAlbedo

  implicit none

  private

  save

  logical, public :: do_initsurfalb

  public :: InitSurfAlb

contains
  !
  ! The variable, h2osoi_vol, is needed by the soil albedo routine
  ! this is not needed on restart since it is computed before the
  ! soil albedo computation is called.
  ! The remaining variables are initialized by calls to ecosystem
  ! dynamics and albedo subroutines.
  !
  subroutine initSurfalb(calday,declin,declinm1)
    implicit none
    real(rk8) , intent(in) :: calday  ! calendar day for declin
    real(rk8) , intent(in) :: declin  ! declination angle (radians) for calday
    ! declination angle (radians) for caldaym1
    real(rk8) , intent(in) , optional :: declinm1
    ! landunit index associated with each pft
    integer(ik4) , pointer :: plandunit(:)
    ! column type
    integer(ik4) , pointer :: ctype(:)
    ! landunit index associated with each column
    integer(ik4) , pointer :: clandunit(:)
    ! gridcell associated with each pft
    integer(ik4),  pointer :: pgridcell(:)
    ! landunit type
    integer(ik4) , pointer :: itypelun(:)
    ! true => landunit is a lake point
    logical , pointer :: lakpoi(:)
    ! layer thickness depth (m)
    real(rk8) , pointer :: dz(:,:)
    ! ice lens (kg/m2)
    real(rk8) , pointer :: h2osoi_ice(:,:)
    ! liquid water (kg/m2)
    real(rk8) , pointer :: h2osoi_liq(:,:)
    ! snow water (mm H2O)
    real(rk8) , pointer :: h2osno(:)
    ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer :: frac_veg_nosno_alb(:)
    ! daylength (seconds)
    real(rk8) , pointer :: dayl(:)
    ! latitude (degrees)
    real(rk8) , pointer :: latdeg(:)
    ! index into column level quantities
    integer(ik4) , pointer :: pcolumn(:)
    ! soil water potential in each soil layer (MPa)
    real(rk8) , pointer :: soilpsi(:,:)
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8) , pointer :: h2osoi_vol(:,:)
    ! snow height (m)
    real(rk8) , pointer :: snow_depth(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rk8) , pointer :: frac_sno(:)
    ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer :: frac_veg_nosno(:)
    ! fraction of canopy that is wet (0 to 1) (pft-level)
    real(rk8) , pointer :: fwet(:)
    ! solar declination angle (radians)
    real(rk8) , pointer :: decl(:)
    ! fraction of foliage that is green and dry [-] (new)
    real(rk8) , pointer :: fdry(:)
    ! one-sided leaf area index, no burying by snow
    real(rk8) , pointer :: tlai(:)
    ! one-sided stem area index, no burying by snow
    real(rk8) , pointer :: tsai(:)
    ! canopy top (m)
    real(rk8) , pointer :: htop(:)
    ! canopy bottom (m)
    real(rk8) , pointer :: hbot(:)
    ! one-sided leaf area index with burying by snow
    real(rk8) , pointer :: elai(:)
    ! one-sided stem area index with burying by snow
    real(rk8) , pointer :: esai(:)
    integer(ik4) :: j , l , c , p , fc ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending ldunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    real(rk8) :: lat    ! latitude (radians) for daylength calculation
    real(rk8) :: temp   ! temporary variable for daylength
    real(rk8) :: snowbd ! temporary calculation of snow bulk density (kg/m3)
    real(rk8) :: fmelt  ! snowbd/100

    ! Assign local pointers to derived subtypes components (landunit-level)

    lakpoi      => clm3%g%l%lakpoi
    itypelun    => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    dz          => clm3%g%l%c%cps%dz
    h2osoi_ice  => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol  => clm3%g%l%c%cws%h2osoi_vol
    snow_depth  => clm3%g%l%c%cps%snow_depth
    h2osno      => clm3%g%l%c%cws%h2osno
    frac_sno    => clm3%g%l%c%cps%frac_sno
    ctype       => clm3%g%l%c%itype
    clandunit   => clm3%g%l%c%landunit
    soilpsi     => clm3%g%l%c%cps%soilpsi

    ! Assign local pointers to derived subtypes components (pft-level)

    plandunit          => clm3%g%l%c%p%landunit
    frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno
    fwet               => clm3%g%l%c%p%pps%fwet

    ! Assign local pointers to derived subtypes components (pft-level)
    ! The folowing pointers will only be used for lake points in this routine

    htop => clm3%g%l%c%p%pps%htop
    hbot => clm3%g%l%c%p%pps%hbot
    tlai => clm3%g%l%c%p%pps%tlai
    tsai => clm3%g%l%c%p%pps%tsai
    elai => clm3%g%l%c%p%pps%elai
    esai => clm3%g%l%c%p%pps%esai
    fdry => clm3%g%l%c%p%pps%fdry

    decl      => clm3%g%l%c%cps%decl
    dayl      => clm3%g%l%c%p%pepv%dayl
    pcolumn   => clm3%g%l%c%p%column
    pgridcell => clm3%g%l%c%p%gridcell
    latdeg    => clm3%g%latdeg

    ! ========================================================================
    ! Determine surface albedo - initialized by calls to ecosystem dynamics and
    ! albedo subroutines. Note: elai, esai, frac_veg_nosno_alb are computed in
    ! Ecosysdyn and needed by routines FracWet and SurfaceAlbedo and
    ! frac_veg_nosno is needed by FracWet
    ! fwet is needed in routine TwoStream (called by SurfaceAlbedo)
    ! frac_sno is needed by SoilAlbedo (called by SurfaceAlbedo)
    ! ========================================================================

#if (!defined CN)
    ! the default mode uses prescribed vegetation structure
    ! Read monthly vegetation data for interpolation to daily values

    call interpMonthlyVeg()
#endif

    ! Determine proc bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Determine variables needed by SurfaceAlbedo for lake points

    do p = begp , endp
      l = plandunit(p)
      if ( lakpoi(l) ) then
        fwet(p) = 0.D0
        fdry(p) = 0.D0
        elai(p) = 0.D0
        esai(p) = 0.D0
        htop(p) = 0.D0
        hbot(p) = 0.D0
        tlai(p) = 0.D0
        tsai(p) = 0.D0
        frac_veg_nosno_alb(p) = 0
        frac_veg_nosno(p) = 0
      end if
    end do

    !
    ! Ecosystem dynamics: Uses CN, or static parameterizations
    !

#if (defined CN)
    do j = 1 , nlevgrnd
      do fc = 1 , filter%num_soilc
        c = filter%soilc(fc)
        soilpsi(c,j) = -15.0D0
      end do
    end do
#endif

    ! Determine variables needed for SurfaceAlbedo for non-lake points

    do c = begc , endc
      l = clandunit(c)
      if ( itypelun(l) == isturb ) then
        ! From Bonan 1996 (LSM technical note)
        frac_sno(c) = min( snow_depth(c)/0.05D0, 1.D0)
      else
        frac_sno(c) = 0.D0
        ! snow cover fraction as in Niu and Yang 2007
        if ( snow_depth(c) .gt. 0.0 ) then
          !bulk density of snow (kg/m3)
          snowbd = min(400.D0,h2osno(c)/snow_depth(c))
          fmelt    = (snowbd/100.)**1.
          ! 100 is the assumed fresh snow density
          !   1 is a melting factor that could be reconsidered
          ! optimal value of 1.5 in Niu et al., 2007
          frac_sno(c) = tanh( snow_depth(c) /(2.5 * zlnd * fmelt) )
        end if
      end if
    end do

#if defined (CN)
    ! CN initialization is done only on the soil landunits.

    if ( .not. present(declinm1) ) then
      write(stderr,*) &
              'declination for the previous timestep (declinm1) must be ',&
              ' present as argument in CN mode'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! it is necessary to initialize the solar declination for the previous
    ! timestep (caldaym1) so that the CNphenology routines know if this is
    ! before or after the summer solstice.

    ! declination for previous timestep
    do c = begc , endc
      l = clandunit(c)
      if ( itypelun(l) == istsoil .or. itypelun(l) == istcrop ) then
        decl(c) = declinm1
      end if
    end do

    ! daylength for previous timestep
    do p = begp , endp
      c = pcolumn(p)
      l = plandunit(p)
      if ( itypelun(l) == istsoil .or. itypelun(l) == istcrop ) then
        lat = latdeg(pgridcell(p)) * degrad
        temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
        temp = min(1.D0,max(-1.D0,temp))
        dayl(p) = 2.0D0 * 13750.9871D0 * acos(temp)
      end if
    end do

    ! declination for current timestep
    do c = begc , endc
      l = clandunit(c)
      if ( itypelun(l) == istsoil .or. itypelun(l) == istcrop ) then
        decl(c) = declin
      end if
    end do

    call CNEcosystemDyn(begc,endc,begp,endp, &
            filter%num_soilc, filter%soilc,  &
            filter%num_soilp, filter%soilp,  &
            filter%num_pcropp, filter%pcropp, doalb=.true.)
#else
    ! this is the default call if CN not set

    call EcosystemDyn(begp,endp,filter%num_nolakep,filter%nolakep,doalb=.true.)
#endif

    do p = begp , endp
      l = plandunit(p)
      if ( .not. lakpoi(l) ) then
        frac_veg_nosno(p) = frac_veg_nosno_alb(p)
        fwet(p) = 0.D0
      end if
    end do

    call FracWet(filter%num_nolakep, filter%nolakep)

    ! Compute Surface Albedo - all land points (including lake)
    ! other than urban. Needs as input fracion of soil covered by snow
    ! (Z.-L. Yang U. Texas)

    call SurfaceAlbedo(begg,endg,begc,endc,begp,endp, &
                       filter%num_nourbanc, filter%nourbanc, &
                       filter%num_nourbanp, filter%nourbanp, &
                       calday, declin)


    ! Determine albedos for urban landunits

    if ( filter%num_urbanl > 0 ) then
      call UrbanAlbedo(begl,endl,begc,endc,begp,endp, &
                       filter%num_urbanl, filter%urbanl, &
                       filter%num_urbanc, filter%urbanc, &
                       filter%num_urbanp, filter%urbanp )

    end if
  end subroutine initSurfalb

end module mod_clm_initsurfalb
