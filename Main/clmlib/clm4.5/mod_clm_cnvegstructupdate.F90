module mod_clm_cnvegstructupdate

#ifdef CN
  !
  ! Module for vegetation structure updates (LAI, SAI, htop, hbot)
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams, only : dtsrf

  implicit none

  save

  private

  logical :: modis_limited = .false.

  public :: CNVegStructUpdate

  contains
  !
  ! On the radiation time step, use C state variables and epc to diagnose
  ! vegetation structure (LAI, SAI, height)
  !
  subroutine CNVegStructUpdate(num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_pftvarcon, only : noveg, nc3crop, nc3irrig, &
      nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub
    use mod_clm_pftvarcon, only : ncorn, ncornirrig, npcropmin, ztopmx, laimx
    use mod_clm_varcon, only : rpi
    !use mod_clm_pftvarcon, only : nbrdlf_evr_trp_tree, &
     ! nc4_grass, nc3_nonarctic_grass
    implicit none
    ! number of column soil points in pft filter
    integer(ik4), intent(in) :: num_soilp
    ! pft filter for soil points
    integer(ik4), intent(in) :: filter_soilp(:)

#if (defined CNDV)
    real(rk8), pointer, contiguous :: allom2(:) ! ecophys const
    real(rk8), pointer, contiguous :: allom3(:) ! ecophys const
    real(rk8), pointer, contiguous :: nind(:)   ! number of individuals (#/m**2)
    ! fractional area of pft (pft area/nat veg area)
    real(rk8), pointer, contiguous :: fpcgrid(:)
#endif
    integer(ik4), pointer, contiguous :: ivt(:)     ! pft vegetation type
    integer(ik4), pointer, contiguous :: pcolumn(:) ! column index associated with each pft
    integer(ik4), pointer, contiguous :: pgridcell(:) ! pft's gridcell index
    real(rk8), pointer, contiguous :: snow_depth(:)    ! snow height (m)
    real(rk8), pointer, contiguous :: leafc(:)      ! (gC/m2) leaf C
    real(rk8), pointer, contiguous :: deadstemc(:)  ! (gC/m2) dead stem C
    !binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8), pointer, contiguous :: woody(:)
    !specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous :: slatop(:)
    !dSLA/dLAI, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous :: dsladlai(:)
    !ratio of momentum roughness length to canopy top height (-)
    real(rk8), pointer, contiguous :: z0mr(:)
    !ratio of displacement height to canopy top height (-)
    real(rk8), pointer, contiguous :: displar(:)
    ! observational height of wind at pft-level [m]
    real(rk8), pointer, contiguous :: forc_hgt_u_pft(:)
    real(rk8), pointer, contiguous :: dwood(:)      ! density of wood (gC/m^3)
    real(rk8), pointer, contiguous :: farea_burned(:) !F. Li and S. Levis

    ! frac of vegetation not covered by snow [-]
    integer(ik4), pointer, contiguous :: frac_veg_nosno_alb(:)
    !one-sided leaf area index, no burying by snow
    real(rk8), pointer, contiguous :: tlai(:)
    !one-sided stem area index, no burying by snow
    real(rk8), pointer, contiguous :: tsai(:)
    real(rk8), pointer, contiguous :: htop(:) !canopy top (m)
    real(rk8), pointer, contiguous :: hbot(:) !canopy bottom (m)
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer, contiguous :: elai(:)
    ! one-sided stem area index with burying by snow
    real(rk8), pointer, contiguous :: esai(:)
    ! max hgt attained by a crop during yr (m)
    real(rk8), pointer, contiguous :: htmx(:)
    integer(ik4), pointer, contiguous :: peaklai(:)  ! 1: max allowed lai; 0: not at max
    integer(ik4), pointer, contiguous :: harvdate(:) ! harvest date

    integer(ik4) :: p,c,g  !indices
    integer(ik4) :: fp     !lake filter indices
    ! ratio of height:radius_breast_height (tree allometry)
    real(rk8):: taper
    real(rk8):: stocking   ! #stems / ha (stocking density)
    real(rk8):: ol         ! thickness of canopy layer covered by snow (m)
    real(rk8):: fb         ! fraction of canopy layer covered by snow
    real(rk8) :: tlai_old  ! for use in Zeng tsai formula
    real(rk8) :: tsai_old  ! for use in Zeng tsai formula
    real(rk8) :: tsai_min  ! PFT derived minimum tsai
    real(rk8) :: tsai_alpha  ! monthly decay rate of tsai
    ! samy : minimum leaf area index for pft retrieved from BATS scheme
    real(rk8) :: tlai_min
    real(rk8) dt             ! radiation time step (sec)

    ! number of seconds in a 30 day month (60x60x24x30)
    real(rk8), parameter :: dtsmonth = 2592000._rk8
    !-------------------------------------------------------------------
    ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835
    !
    ! tsai(p) = max( tsai_alpha(ivt(p))*tsai_old + &
    !           max(tlai_old-tlai(p),0_rk8), tsai_min(ivt(p)) )
    ! notes:
    ! * RHS tsai & tlai are from previous timestep
    ! * should create tsai_alpha(ivt(p)) & tsai_min(ivt(p))
    !      in pftvarcon.F90 - slevis
    ! * all non-crop pfts use same values:
    !   crop    tsai_alpha,tsai_min = 0.0,0.1
    !   noncrop tsai_alpha,tsai_min = 0.5,1.0 (includes bare soil and urban)
    !--------------------------------------------------------------------

    ! assign local pointers to derived type arrays (in)
#if (defined CNDV)
    allom2   => dgv_pftcon%allom2
    allom3   => dgv_pftcon%allom3
    nind     => clm3%g%l%c%p%pdgvs%nind
    fpcgrid  => clm3%g%l%c%p%pdgvs%fpcgrid
#endif
    ivt          => clm3%g%l%c%p%itype
    pcolumn      => clm3%g%l%c%p%column
    pgridcell    => clm3%g%l%c%p%gridcell
    leafc        => clm3%g%l%c%p%pcs%leafc
    deadstemc    => clm3%g%l%c%p%pcs%deadstemc
    snow_depth   => clm3%g%l%c%cps%snow_depth
    woody        => pftcon%woody
    slatop       => pftcon%slatop
    dsladlai     => pftcon%dsladlai
    z0mr         => pftcon%z0mr
    displar      => pftcon%displar
    dwood        => pftcon%dwood
    farea_burned => clm3%g%l%c%cps%farea_burned

    ! assign local pointers to derived type arrays (out)
    tlai               => clm3%g%l%c%p%pps%tlai
    tsai               => clm3%g%l%c%p%pps%tsai
    htop               => clm3%g%l%c%p%pps%htop
    hbot               => clm3%g%l%c%p%pps%hbot
    elai               => clm3%g%l%c%p%pps%elai
    esai               => clm3%g%l%c%p%pps%esai
    frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    htmx               => clm3%g%l%c%p%pps%htmx
    peaklai            => clm3%g%l%c%p%pps%peaklai
    harvdate           => clm3%g%l%c%p%pps%harvdate
    forc_hgt_u_pft     => clm3%g%l%c%p%pps%forc_hgt_u_pft

    dt = dtsrf

    ! constant allometric parameters
    taper = 200._rk8
    stocking = 1000._rk8

    ! convert from stems/ha -> stems/m^2
    stocking = stocking / 10000._rk8

    ! pft loop
    do fp = 1, num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)
      g = pgridcell(p)

      if ( ivt(p) /= noveg ) then

        tlai_old = tlai(p) ! n-1 value
        tsai_old = tsai(p) ! n-1 value

        ! update the leaf area index based on leafC and SLA
        ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923.
        ! Samy
        ! This modification is done to obtain LAI in the range of MODIS
        ! satellite derived values
        ! To get reseanoable values for LAI or GPP for Savannah region
        ! (as a  case study) BIOME-BGC model was used with multiplication
        ! coefficient specific to PFT under study
        if ( modis_limited ) then
          if ( dsladlai(ivt(p)) > 0._rk8 ) then
            tlai(p) = (slatop(ivt(p)) * &
              (exp(leafc(p)*dsladlai(ivt(p))) - 1.0_rk8))/dsladlai(ivt(p))
          else
             if ( ivt(p) >= npcropmin ) then
               tlai(p) = slatop(ivt(p)) * leafc(p)
             else
               tlai(p) = 0.25_rk8 * slatop(ivt(p)) * leafc(p)
             end if
          end if
          tlai(p) = max(0._rk8, tlai(p))
        else
          if (dsladlai(ivt(p)) > 0.0_rk8) then
            tlai(p) = (slatop(ivt(p)) * &
              (exp(leafc(p)*dsladlai(ivt(p))) - 1.0_rk8))/dsladlai(ivt(p))
          else
            tlai(p) = slatop(ivt(p)) * leafc(p)
          end if
          tlai(p) = max(0.0_rk8, tlai(p))
        end if

        ! update the stem area index and height based on LAI,
        ! stem mass, and veg type.
        ! With the exception of htop for woody vegetation, this follows
        ! the DGVM logic.

        ! tsai formula from Zeng et. al. 2002,
        ! Journal of Climate, p1835 (see notes)
        ! Assumes doalb time step == CLM time step,
        ! SAI min and monthly decay factor
        ! alpha are set by PFT, and alpha is scaled to CLM
        ! time step by multiplying by dt and dividing by dtsmonth
        ! (seconds in average 30 day month)
        ! tsai_min scaled by 0.5 to match MODIS satellite derived values
        if (ivt(p) == nc3crop .or. ivt(p) == nc3irrig) then ! generic crops
          tsai_alpha = 1.0_rk8-1.0_rk8*dt/dtsmonth
          tsai_min = 0.1_rk8
        else
          tsai_alpha = 1.0_rk8-0.5_rk8*dt/dtsmonth
          tsai_min = 1.0_rk8
        end if
        tsai_min = tsai_min * 0.5_rk8
        tsai(p) = max(tsai_alpha*tsai_old+max(tlai_old-tlai(p),0._rk8),tsai_min)

        if ( abs(woody(ivt(p))-1._rk8) < epsilon(1.0) ) then

          ! trees and shrubs

          ! if shrubs have a squat taper
          if ( ivt(p) >= nbrdlf_evr_shrub .and. &
               ivt(p) <= nbrdlf_dcd_brl_shrub ) then
            taper = 10._rk8
            ! otherwise have a tall taper
          else
            taper = 200._rk8
          end if

          ! trees and shrubs for now have a very simple allometry,
          ! with hard-wired stem taper (height:radius) and hard-wired
          ! stocking density (#individuals/area)
#if (defined CNDV)
          if ( fpcgrid(p) > 0._rk8 .and. nind(p) > 0._rk8 ) then
            !#ind/m2 nat veg area -> #ind/m2 pft area
            stocking = nind(p)/fpcgrid(p)
            htop(p) = allom2(ivt(p)) * ( (24._rk8 * deadstemc(p) / &
                   (rpi * stocking * dwood(ivt(p)) * taper))** &
                   (1._rk8/3._rk8) )**allom3(ivt(p)) ! lpj's htop w/ cn's stemdiam
          else
            htop(p) = 0._rk8
          end if
#else
          htop(p) = ((3._rk8 * deadstemc(p) * taper * taper)/ &
                   (rpi * stocking * dwood(ivt(p))))**(1._rk8/3._rk8)
#endif

          ! Peter Thornton, 5/3/2004
          ! Adding test to keep htop from getting too close to
          ! forcing height for windspeed
          ! Also added for grass, below, although it is not likely
          ! to ever be an issue.
          htop(p) = min(htop(p), &
            (forc_hgt_u_pft(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._rk8)

          ! Peter Thornton, 8/11/2004
          ! Adding constraint to keep htop from going to 0.0.
          ! This becomes an issue when fire mortality is pushing deadstemc
          ! to 0.0.
          htop(p) = max(htop(p), 0.01_rk8)
          hbot(p) = max(0._rk8, min(3._rk8, htop(p)-1._rk8))

        else if ( ivt(p) >= npcropmin ) then ! prognostic crops

          if (tlai(p) >= laimx(ivt(p))) peaklai(p) = 1 ! used in CNAllocation

          if ( ivt(p) == ncorn .or. ivt(p) == ncornirrig ) then
            tsai(p) = 0.1_rk8 * tlai(p)
          else
            tsai(p) = 0.2_rk8 * tlai(p)
          end if

          ! "stubble" after harvest
          if ( harvdate(p) < 999 .and. tlai(p) == 0._rk8 ) then
            !changed by F. Li and S. Levis
            tsai(p) = 0.25_rk8*(1._rk8-farea_burned(c)*0.90_rk8)
            htmx(p) = 0._rk8
            peaklai(p) = 0
          end if

          ! canopy top and bottom heights
          htop(p) = ztopmx(ivt(p)) * &
            (min(tlai(p)/(laimx(ivt(p))-1._rk8),1._rk8))**2
          htmx(p) = max(htmx(p), htop(p))
          htop(p) = max(0.05_rk8, max(htmx(p),htop(p)))
          hbot(p) = 0.02_rk8
        else ! generic crops and ...
          ! grasses

          ! height for grasses depends only on LAI
          htop(p) = max(0.25_rk8, tlai(p) * 0.25_rk8)

          htop(p) = min(htop(p), &
            (forc_hgt_u_pft(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._rk8)

          ! Peter Thornton, 8/11/2004
          ! Adding constraint to keep htop from going to 0.0.
          htop(p) = max(htop(p), 0.01_rk8)
          hbot(p) = max(0.0_rk8, min(0.05_rk8, htop(p)-0.20_rk8))
        end if
      else
        tlai(p) = 0._rk8
        tsai(p) = 0._rk8
        htop(p) = 0._rk8
        hbot(p) = 0._rk8
      end if

      ! adjust lai and sai for burying by snow.

      ! snow burial fraction for short vegetation (e.g. grasses) as in
      ! Wang and Zeng, 2007.
      if ( ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
        ol = min( max(snow_depth(c)-hbot(p), 0._rk8), htop(p)-hbot(p))
        fb = 1._rk8 - ol / max(1.e-6_rk8, htop(p)-hbot(p))
      else
        ! 0.2m is assumed depth of snow required for
        ! complete burial of grasses
        fb = 1._rk8 - max(min(snow_depth(c),0.2_rk8),0._rk8)/0.2_rk8
      end if

      elai(p) = max(tlai(p)*fb, 0.0_rk8)
      esai(p) = max(tsai(p)*fb, 0.0_rk8)

      ! Fraction of vegetation free of snow
      if ( (elai(p) + esai(p)) > 0._rk8 ) then
        frac_veg_nosno_alb(p) = 1
      else
        frac_veg_nosno_alb(p) = 0
      end if
    end do

  end subroutine CNVegStructUpdate

#endif

end module mod_clm_cnvegstructupdate
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
