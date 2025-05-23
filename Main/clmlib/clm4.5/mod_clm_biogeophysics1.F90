module mod_clm_biogeophysics1
  !
  ! Performs calculation of leaf temperature and surface fluxes.
  ! Biogeophysics2.F90 then determines soil/snow and ground
  ! temperatures and updates the surface fluxes for the new ground
  ! temperature.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_clm_atmlnd, only : clm_a2l
  use mod_clm_varcon, only : denh2o, denice, roverg, hvap, hsub, &
                       istice, istwet, istsoil, isturb, istdlak,  &
                       zlnd, zsno, tfrz, icol_roof, icol_sunwall, &
                       icol_shadewall, icol_road_imperv,          &
                       icol_road_perv, tfrz, spval, istdlak
  use mod_clm_varcon, only : istcrop
  use mod_clm_varpar, only : nlevgrnd, nlevurb, nlevsno, nlevsoi
  use mod_clm_varpar, only : max_pft_per_gcell
  use mod_clm_qsat, only : QSat
  use mod_constants, only : mathpi

  implicit none

  private

  save

  public :: Biogeophysics1   ! Calculate leaf temperature and surface fluxes

  contains
  !
  ! This is the main subroutine to execute the calculation of leaf temperature
  ! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
  ! temperatures and updates the surface fluxes for the new ground
  ! temperature.
  !
  ! Calling sequence is:
  ! Biogeophysics1:           surface biogeophysics driver
  !  -> QSat:                 saturated vapor pressure, specific humidity, and
  !                           derivatives at ground surface and derivatives at
  !                           leaf surface using updated leaf temperature
  ! Leaf temperature
  ! Foliage energy conservation is given by the foliage energy budget
  ! equation:
  !                Rnet - Hf - LEf = 0
  ! The equation is solved by Newton-Raphson iteration, in which this
  ! iteration includes the calculation of the photosynthesis and
  ! stomatal resistance, and the integration of turbulent flux profiles.
  ! The sensible and latent heat transfer between foliage and atmosphere
  ! and ground is linked by the equations:
  !                Ha = Hf + Hg and Ea = Ef + Eg
  !
  subroutine Biogeophysics1(lbg, ubg, lbc, ubc, lbp, ubp, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep)
    implicit none
    integer(ik4), intent(in) :: lbg, ubg  ! gridcell-index bounds
    integer(ik4), intent(in) :: lbc, ubc  ! column-index bounds
    integer(ik4), intent(in) :: lbp, ubp  ! pft-index bounds
    ! number of column non-lake points in column filter
    integer(ik4), intent(in) :: num_nolakec
    ! column filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakec(ubc-lbc+1)
    ! number of column non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakep(ubp-lbp+1)

    ! eff. fraction of ground covered by snow (0 to 1)
    real(rk8), pointer, contiguous :: frac_sno_eff(:)
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8), pointer, contiguous :: frac_h2osfc(:)
    ! surface water (mm)
    real(rk8), pointer, contiguous :: h2osfc(:)
    ! surface water temperature
    real(rk8), pointer, contiguous :: t_h2osfc(:)
    ! saved surface water temperature
    real(rk8), pointer, contiguous :: t_h2osfc_bef(:)
    ! specific humidity at snow surface [kg/kg]
    real(rk8), pointer, contiguous :: qg_snow(:)
    ! specific humidity at soil surface [kg/kg]
    real(rk8), pointer, contiguous :: qg_soil(:)
    ! specific humidity at h2osfc surface [kg/kg]
    real(rk8), pointer, contiguous :: qg_h2osfc(:)
    !true=>do computations on this pft (see reweightMod for details)
    logical, pointer, contiguous :: pactive(:)
    integer(ik4), pointer, contiguous :: ivt(:)         !pft vegetation type
    integer(ik4), pointer, contiguous :: ityplun(:)     !landunit type
    integer(ik4), pointer, contiguous :: clandunit(:)   !column's landunit index
    integer(ik4), pointer, contiguous :: cgridcell(:)   !column's gridcell index
    integer(ik4), pointer, contiguous :: ctype(:)       !column type
    real(rk8), pointer, contiguous :: forc_pbot(:)  !atmospheric pressure (Pa)
    real(rk8), pointer, contiguous :: forc_q(:)     !atmospheric specific humidity (kg/kg)
    real(rk8), pointer, contiguous :: forc_t(:)     !atmospheric temperature (Kelvin)
    real(rk8), pointer, contiguous :: forc_hgt_t(:) !observational height of temperature [m]
    real(rk8), pointer, contiguous :: forc_hgt_u(:) !observational height of wind [m]
    !observational height of specific humidity [m]
    real(rk8), pointer, contiguous :: forc_hgt_q(:)
    integer(ik4), pointer, contiguous :: npfts(:)         !number of pfts on gridcell
    integer(ik4), pointer, contiguous :: pfti(:)          !initial pft on gridcell
    integer(ik4), pointer, contiguous :: plandunit(:)     !pft's landunit index
    !observational height of wind at pft level [m]
    real(rk8), pointer, contiguous :: forc_hgt_u_pft(:)
    !observational height of temperature at pft level [m]
    real(rk8), pointer, contiguous :: forc_hgt_t_pft(:)
    !observational height of specific humidity at pft level [m]
    real(rk8), pointer, contiguous :: forc_hgt_q_pft(:)
    !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4), pointer, contiguous :: frac_veg_nosno(:)
    integer(ik4), pointer, contiguous :: pgridcell(:)    !pft's gridcell index
    integer(ik4), pointer, contiguous :: pcolumn(:)      !pft's column index
    !momentum roughness length of urban landunit (m)
    real(rk8), pointer, contiguous :: z_0_town(:)
    !displacement height of urban landunit (m)
    real(rk8), pointer, contiguous :: z_d_town(:)
    !atmospheric potential temperature (Kelvin)
    real(rk8), pointer, contiguous :: forc_th(:)
    !atmospheric wind speed in east direction (m/s)
    real(rk8), pointer, contiguous :: forc_u(:)
    !atmospheric wind speed in north direction (m/s)
    real(rk8), pointer, contiguous :: forc_v(:)
    !restriction for min of soil potential (mm)
    real(rk8), pointer, contiguous :: smpmin(:)
    integer(ik4), pointer, contiguous :: snl(:)  !number of snow layers
    !fraction of ground covered by snow (0 to 1)
    real(rk8), pointer, contiguous :: frac_sno(:)
    real(rk8), pointer, contiguous :: h2osno(:) !snow water (mm H2O)
    !one-sided leaf area index with burying by snow
    real(rk8), pointer, contiguous :: elai(:)
    !one-sided stem area index with burying by snow
    real(rk8), pointer, contiguous :: esai(:)
    !ratio of momentum roughness length to canopy top height (-)
    real(rk8), pointer, contiguous :: z0mr(:)
    !ratio of displacement height to canopy top height (-)
    real(rk8), pointer, contiguous :: displar(:)
    real(rk8), pointer, contiguous :: htop(:)    !canopy top (m)
    real(rk8), pointer, contiguous :: dz(:,:)    !layer depth (m)
    real(rk8), pointer, contiguous :: t_soisno(:,:) !soil temperature (Kelvin)
    real(rk8), pointer, contiguous :: h2osoi_liq(:,:) !liquid water (kg/m2)
    real(rk8), pointer, contiguous :: h2osoi_ice(:,:) !ice lens (kg/m2)
    !volumetric soil water at saturation (porosity)
    real(rk8), pointer, contiguous :: watsat(:,:)
    real(rk8), pointer, contiguous :: sucsat(:,:)   !minimum soil suction (mm)
    real(rk8), pointer, contiguous :: bsw(:,:)      !Clapp and Hornberger "b"
    !volumetric soil water at field capacity
    real(rk8), pointer, contiguous :: watfc(:,:)
    !volumetric soil moisture corresponding to no restriction on ET
    ! from urban pervious surface
    real(rk8), pointer, contiguous :: watopt(:,:)
    !volumetric soil moisture corresponding to no restriction on ET
    ! from urban pervious surface
    real(rk8), pointer, contiguous :: watdry(:,:)
    !fraction of roots in each soil layer for urban pervious road
    real(rk8), pointer, contiguous :: rootfr_road_perv(:,:)
    !effective fraction of roots in each soil layer for urban pervious road
    real(rk8), pointer, contiguous :: rootr_road_perv(:,:)

    real(rk8), pointer, contiguous :: t_grnd(:) !ground temperature (Kelvin)
    real(rk8), pointer, contiguous :: qg(:)     !ground specific humidity [kg/kg]
    real(rk8), pointer, contiguous :: dqgdT(:)  !d(qg)/dT
    real(rk8), pointer, contiguous :: emg(:)    !ground emissivity
    !latent heat of vapor of water (or sublimation) [j/kg]
    real(rk8), pointer, contiguous :: htvp(:)
    real(rk8), pointer, contiguous :: beta(:)   !coefficient of convective velocity [-]
    real(rk8), pointer, contiguous :: zii(:)    !convective boundary height [m]
    !intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(rk8), pointer, contiguous :: thm(:)
    real(rk8), pointer, contiguous :: thv(:)    !virtual potential temperature (kelvin)
    !roughness length over ground, momentum [m]
    real(rk8), pointer, contiguous :: z0mg(:)
    !roughness length over ground, sensible heat [m]
    real(rk8), pointer, contiguous :: z0hg(:)
    !roughness length over ground, latent heat [m]
    real(rk8), pointer, contiguous :: z0qg(:)
    real(rk8), pointer, contiguous :: emv(:)    !vegetation emissivity
    real(rk8), pointer, contiguous :: z0m(:)    !momentum roughness length (m)
    real(rk8), pointer, contiguous :: displa(:) !displacement height (m)
    !roughness length over vegetation, momentum [m]
    real(rk8), pointer, contiguous :: z0mv(:)
    !roughness length over vegetation, sensible heat [m]
    real(rk8), pointer, contiguous :: z0hv(:)
    !roughness length over vegetation, latent heat [m]
    real(rk8), pointer, contiguous :: z0qv(:)
    !total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous :: eflx_sh_tot(:)
    !urban total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous :: eflx_sh_tot_u(:)
    !rural total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous :: eflx_sh_tot_r(:)
    !total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous :: eflx_lh_tot(:)
    !urban total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous :: eflx_lh_tot_u(:)
    !rural total latent heat flux (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous :: eflx_lh_tot_r(:)
    !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous :: eflx_sh_veg(:)
    !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8), pointer, contiguous :: qflx_evap_tot(:)
    !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous :: qflx_evap_veg(:)
    !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(rk8), pointer, contiguous :: qflx_tran_veg(:)
    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rk8), pointer, contiguous :: cgrnd(:)
    !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
    real(rk8), pointer, contiguous :: cgrnds(:)
    !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
    real(rk8), pointer, contiguous :: cgrndl(:)
    !soil/snow temperature before update
    real(rk8) ,pointer, contiguous :: tssbef(:,:)
    !factor that reduces ground saturated specific humidity (-)
    real(rk8) ,pointer, contiguous :: soilalpha(:)
    !factor that reduces ground evaporation
    real(rk8) ,pointer, contiguous :: soilbeta(:)
    !Urban factor that reduces ground saturated specific humidity (-)
    real(rk8) ,pointer, contiguous :: soilalpha_u(:)

    integer(ik4)  :: g,l,c,p !indices
    integer(ik4)  :: j       !soil/snow level index
    integer(ik4)  :: fp      !lake filter pft index
    integer(ik4)  :: fc      !lake filter column index
    real(rk8) :: qred    !soil surface relative humidity
    real(rk8) :: avmuir  !ir inverse optical depth per unit leaf area
    real(rk8) :: eg      !water vapor pressure at temperature T [pa]
    real(rk8) :: qsatg   !saturated humidity [kg/kg]
    real(rk8) :: degdT   !d(eg)/dT
    real(rk8) :: qsatgdT !d(qsatg)/dT
    real(rk8) :: fac     !soil wetness of surface layer
    real(rk8) :: psit    !negative potential of soil
    real(rk8) :: hr      !relative humidity
    real(rk8) :: hr_road_perv  !relative humidity for urban pervious road
    real(rk8) :: wx      !partial volume of ice and water of surface layer
    !soil wetness of surface layer relative to field capacity
    real(rk8) :: fac_fc
    real(rk8) :: eff_porosity  ! effective porosity in layer
    real(rk8) :: vol_ice       ! partial volume of ice lens in layer
    real(rk8) :: vol_liq       ! partial volume of liquid water in layer
    integer(ik4)  :: pi            !index

    ! Assign local pointers to derived type members (gridcell-level)

    frac_sno_eff   => clm3%g%l%c%cps%frac_sno_eff
    frac_h2osfc   => clm3%g%l%c%cps%frac_h2osfc
    h2osfc        => clm3%g%l%c%cws%h2osfc
    t_h2osfc      => clm3%g%l%c%ces%t_h2osfc
    t_h2osfc_bef  => clm3%g%l%c%ces%t_h2osfc_bef
    qg_snow       => clm3%g%l%c%cws%qg_snow
    qg_soil       => clm3%g%l%c%cws%qg_soil
    qg_h2osfc     => clm3%g%l%c%cws%qg_h2osfc
    forc_hgt_t    => clm_a2l%forc_hgt_t
    forc_u        => clm_a2l%forc_u
    forc_v        => clm_a2l%forc_v
    forc_hgt_u    => clm_a2l%forc_hgt_u
    forc_hgt_q    => clm_a2l%forc_hgt_q
    npfts         => clm3%g%npfts
    pfti          => clm3%g%pfti

    ! Assign local pointers to derived type members (landunit-level)

    ityplun       => clm3%g%l%itype
    z_0_town      => clm3%g%l%z_0_town
    z_d_town      => clm3%g%l%z_d_town

    ! Assign local pointers to derived type members (column-level)

    forc_pbot     => clm3%g%l%c%cps%forc_pbot
    forc_q        => clm3%g%l%c%cws%forc_q
    forc_t        => clm3%g%l%c%ces%forc_t
    forc_th       => clm3%g%l%c%ces%forc_th

    cgridcell     => clm3%g%l%c%gridcell
    clandunit     => clm3%g%l%c%landunit
    ctype         => clm3%g%l%c%itype
    beta          => clm3%g%l%c%cps%beta
    dqgdT         => clm3%g%l%c%cws%dqgdT
    emg           => clm3%g%l%c%cps%emg
    frac_sno      => clm3%g%l%c%cps%frac_sno
    h2osno        => clm3%g%l%c%cws%h2osno
    htvp          => clm3%g%l%c%cps%htvp
    qg            => clm3%g%l%c%cws%qg
    smpmin        => clm3%g%l%c%cps%smpmin
    snl           => clm3%g%l%c%cps%snl
    t_grnd        => clm3%g%l%c%ces%t_grnd
    thv           => clm3%g%l%c%ces%thv
    z0hg          => clm3%g%l%c%cps%z0hg
    z0mg          => clm3%g%l%c%cps%z0mg
    z0qg          => clm3%g%l%c%cps%z0qg
    zii           => clm3%g%l%c%cps%zii
    bsw           => clm3%g%l%c%cps%bsw
    dz            => clm3%g%l%c%cps%dz
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    soilalpha     => clm3%g%l%c%cws%soilalpha
    soilbeta      => clm3%g%l%c%cws%soilbeta
    soilalpha_u   => clm3%g%l%c%cws%soilalpha_u
    sucsat        => clm3%g%l%c%cps%sucsat
    t_soisno      => clm3%g%l%c%ces%t_soisno
    tssbef        => clm3%g%l%c%ces%tssbef
    watsat        => clm3%g%l%c%cps%watsat
    watfc         => clm3%g%l%c%cps%watfc
    watdry        => clm3%g%l%c%cps%watdry
    watopt        => clm3%g%l%c%cps%watopt
    rootfr_road_perv => clm3%g%l%c%cps%rootfr_road_perv
    rootr_road_perv  => clm3%g%l%c%cps%rootr_road_perv

    ! Assign local pointers to derived type members (pft-level)

    pactive       => clm3%g%l%c%p%active
    ivt           => clm3%g%l%c%p%itype
    elai          => clm3%g%l%c%p%pps%elai
    esai          => clm3%g%l%c%p%pps%esai
    htop          => clm3%g%l%c%p%pps%htop
    emv           => clm3%g%l%c%p%pps%emv
    z0m           => clm3%g%l%c%p%pps%z0m
    displa        => clm3%g%l%c%p%pps%displa
    z0mv          => clm3%g%l%c%p%pps%z0mv
    z0hv          => clm3%g%l%c%p%pps%z0hv
    z0qv          => clm3%g%l%c%p%pps%z0qv
    eflx_sh_tot   => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_sh_tot_u => clm3%g%l%c%p%pef%eflx_sh_tot_u
    eflx_sh_tot_r => clm3%g%l%c%p%pef%eflx_sh_tot_r
    eflx_lh_tot   => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_lh_tot_u => clm3%g%l%c%p%pef%eflx_lh_tot_u
    eflx_lh_tot_r => clm3%g%l%c%p%pef%eflx_lh_tot_r
    eflx_sh_veg   => clm3%g%l%c%p%pef%eflx_sh_veg
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_veg => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_tran_veg => clm3%g%l%c%p%pwf%qflx_tran_veg
    cgrnd         => clm3%g%l%c%p%pef%cgrnd
    cgrnds        => clm3%g%l%c%p%pef%cgrnds
    cgrndl        => clm3%g%l%c%p%pef%cgrndl
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
    forc_hgt_t_pft => clm3%g%l%c%p%pps%forc_hgt_t_pft
    forc_hgt_q_pft => clm3%g%l%c%p%pps%forc_hgt_q_pft
    plandunit      => clm3%g%l%c%p%landunit
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    thm            => clm3%g%l%c%p%pes%thm
    pgridcell      => clm3%g%l%c%p%gridcell
    pcolumn        => clm3%g%l%c%p%column

    ! Assign local pointers to derived type members (ecophysiological)

    z0mr          => pftcon%z0mr
    displar       => pftcon%displar

    do j = -nlevsno+1, nlevgrnd
      do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof) .and. j > nlevurb) then
          tssbef(c,j) = spval
        else
          tssbef(c,j) = t_soisno(c,j)
        end if
        ! record t_h2osfc prior to updating
        t_h2osfc_bef(c) = t_h2osfc(c)
      end do
    end do

    do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)

      if (ctype(c) == icol_road_perv) then
        hr_road_perv = 0._rk8
      end if

      ! begin calculations that relate only to the column level
      ! Ground and soil temperatures from previous time step

      ! ground temperature is weighted average of exposed soil,
      ! snow, and h2osfc
      if (snl(c) < 0) then
        t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) + &
               (1.0_rk8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) + &
               frac_h2osfc(c) * t_h2osfc(c)
      else
        t_grnd(c) = (1.0_rk8 - frac_h2osfc(c)) * t_soisno(c,1) + &
          frac_h2osfc(c) * t_h2osfc(c)
      end if
      ! Saturated vapor pressure, specific humidity and their derivatives
      ! at ground surface

      qred = 1._rk8
      if (ityplun(l)/=istwet .AND. ityplun(l)/=istice ) then
        if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
          wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)
          fac  = min(1._rk8, wx/watsat(c,1))
          fac  = max( fac, 0.01_rk8 )
          psit = -sucsat(c,1) * fac ** (-bsw(c,1))
          psit = max(smpmin(c), psit)
          ! modify qred to account for h2osfc
          hr   = exp(psit/roverg/t_soisno(c,1))
          qred = (1._rk8 - frac_sno(c) - frac_h2osfc(c))*hr + &
                  frac_sno(c) + frac_h2osfc(c)

          !! Lee and Pielke 1992 beta, added by K.Sakaguchi
          if (wx < watfc(c,1) ) then
            !when water content of ths top layer is less than that at F.C.
            !eqn5.66 but divided by theta at field capacity
            fac_fc  = min(1._rk8, wx/watfc(c,1))
            fac_fc  = max( fac_fc, 0.01_rk8 )
            ! modify soil beta by snow cover. soilbeta for snow surface is one
            soilbeta(c) = (1._rk8-frac_sno(c)-frac_h2osfc(c)) * &
                     0.25_rk8*(1._rk8 - cos(mathpi*fac_fc))**2._rk8 + &
                              frac_sno(c)+ frac_h2osfc(c)
          else   !when water content of ths top layer is more than that at F.C.
            soilbeta(c) = 1._rk8
          end if

          soilalpha(c) = qred
          ! Pervious road depends on water in total soil column
        else if (ctype(c) == icol_road_perv) then
          do j = 1, nlevsoi
            if (t_soisno(c,j) >= tfrz) then
              vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
              eff_porosity = watsat(c,j)-vol_ice
              vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
              fac = min( max(vol_liq-watdry(c,j),0._rk8) / &
                             (watopt(c,j)-watdry(c,j)), 1._rk8 )
            else
              fac = 0._rk8
            end if
            rootr_road_perv(c,j) = rootfr_road_perv(c,j)*fac
            hr_road_perv = hr_road_perv + rootr_road_perv(c,j)
          end do
          ! Allows for sublimation of snow or dew on snow
          qred = (1.0_rk8-frac_sno(c))*hr_road_perv + frac_sno(c)

          ! Normalize root resistances to get layer contribution to total ET
          if (hr_road_perv > 0._rk8) then
            do j = 1, nlevsoi
              rootr_road_perv(c,j) = rootr_road_perv(c,j)/hr_road_perv
            end do
          end if
          soilalpha_u(c) = qred
          soilbeta(c) = 0._rk8
        else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
          qred = 0._rk8
          soilbeta(c) = 0._rk8
          soilalpha_u(c) = spval
        else if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
          qred = 1._rk8
          soilbeta(c) = 0._rk8
          soilalpha_u(c) = spval
        end if
      else
        soilalpha(c) = spval
        soilbeta(c) =   1._rk8
      end if

      ! compute humidities individually for snow, soil, h2osfc for
      ! vegetated landunits
      if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then

        call QSat(t_soisno(c,snl(c)+1), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
        if (qsatg > forc_q(c) .and. forc_q(c) > qsatg) then
          qsatg = forc_q(c)
          qsatgdT = 0._rk8
        end if

        qg_snow(c) = qsatg
        dqgdT(c) = frac_sno(c)*qsatgdT

        call QSat(t_soisno(c,1), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
        if (qsatg > forc_q(c) .and. forc_q(c) > hr*qsatg) then
          qsatg = forc_q(c)
          qsatgdT = 0._rk8
        end if
        qg_soil(c) = hr*qsatg
        dqgdT(c) = dqgdT(c) + (1._rk8 - frac_sno(c) - frac_h2osfc(c))*hr*qsatgdT

        ! to be consistent with hs_top values in SoilTemp, set
        ! qg_snow to qg_soil for snl = 0 case
        ! this ensures hs_top_snow will equal hs_top_soil
        if (snl(c) >= 0) then
          qg_snow(c) = qg_soil(c)
          dqgdT(c) = (1._rk8 - frac_h2osfc(c))*hr*dqgdT(c)
        end if

        call QSat(t_h2osfc(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
        if (qsatg > forc_q(c) .and. forc_q(c) > qsatg) then
          qsatg = forc_q(c)
          qsatgdT = 0._rk8
        end if
        qg_h2osfc(c) = qsatg
        dqgdT(c) = dqgdT(c) + frac_h2osfc(c) * qsatgdT

        ! qg(c) = frac_sno(c)*qg_snow(c) + (1._rk8 - frac_sno(c) - &
        !           frac_h2osfc(c))*qg_soil(c) &
        qg(c) = frac_sno_eff(c)*qg_snow(c) + &
          (1._rk8 - frac_sno_eff(c) - frac_h2osfc(c))*qg_soil(c) &
               + frac_h2osfc(c) * qg_h2osfc(c)

      else
        call QSat(t_grnd(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
        qg(c) = qred*qsatg
        dqgdT(c) = qred*qsatgdT

        if (qsatg > forc_q(c) .and. forc_q(c) > qred*qsatg) then
          qg(c) = forc_q(c)
          dqgdT(c) = 0._rk8
        end if

        qg_snow(c) = qg(c)
        qg_soil(c) = qg(c)
        qg_h2osfc(c) = qg(c)
      end if

      ! Ground emissivity - only calculate for non-urban landunits
      ! Urban emissivities are currently read in from data file

      if ( ityplun(l) /= isturb ) then
        if ( ityplun(l) == istice ) then
          emg(c) = 0.97_rk8
        else
          emg(c) = (1._rk8-frac_sno(c))*0.96_rk8 + frac_sno(c)*0.97_rk8
        end if
      end if

      ! Latent heat. We arbitrarily assume that the sublimation occurs
      ! only as h2osoi_liq = 0

      htvp(c) = hvap
      if ( h2osoi_liq(c,snl(c)+1) <= 0._rk8 .and. &
           h2osoi_ice(c,snl(c)+1) > 0._rk8 ) htvp(c) = hsub

      ! Switch between vaporization and sublimation causes rapid solution
      ! separation in perturbation growth test

#if (defined PERGRO)
      htvp(c) = hvap
#endif

      ! Ground roughness lengths over non-lake columns (includes bare
      ! ground, ground underneath canopy, wetlands, etc.)

      if ( frac_sno(c) > 0._rk8 ) then
        z0mg(c) = zsno
      else
        z0mg(c) = zlnd
      end if
      z0hg(c) = z0mg(c)            ! initial set only
      z0qg(c) = z0mg(c)            ! initial set only

      ! Potential, virtual potential temperature, and wind speed at the
      ! reference height

      beta(c) = 1._rk8
      zii(c)  = 1000._rk8
      thv(c)  = forc_th(c)*(1._rk8+0.61_rk8*forc_q(c))

    end do ! (end of columns loop)

    ! Initialization

    do fp = 1, num_nolakep
      p = filter_nolakep(fp)

      ! Initial set (needed for history tape fields)

      eflx_sh_tot(p) = 0._rk8
      l = plandunit(p)
      if ( ityplun(l) == isturb ) then
        eflx_sh_tot_u(p) = 0._rk8
      else if ( ityplun(l) == istsoil .or. ityplun(l) == istcrop ) then
        eflx_sh_tot_r(p) = 0._rk8
      end if
      eflx_lh_tot(p) = 0._rk8
      if ( ityplun(l) == isturb ) then
        eflx_lh_tot_u(p) = 0._rk8
      else if ( ityplun(l) == istsoil .or. ityplun(l) == istcrop ) then
        eflx_lh_tot_r(p) = 0._rk8
      end if
      eflx_sh_veg(p) = 0._rk8
      qflx_evap_tot(p) = 0._rk8
      qflx_evap_veg(p) = 0._rk8
      qflx_tran_veg(p) = 0._rk8

      ! Initial set for calculation

      cgrnd(p)  = 0._rk8
      cgrnds(p) = 0._rk8
      cgrndl(p) = 0._rk8

      ! Vegetation Emissivity

      avmuir = 1._rk8
      if ( (elai(p)+esai(p))/avmuir > 25.0_rk8 ) then
        emv(p) = 1._rk8
      else
        emv(p) = 1._rk8-exp(-(elai(p)+esai(p))/avmuir)
      end if

      ! Roughness lengths over vegetation

      z0m(p)    = z0mr(ivt(p)) * htop(p)
      displa(p) = displar(ivt(p)) * htop(p)

      if ( ityplun(l) == isturb ) then
        z0mv(p) = z_0_town(l)
      else
        z0mv(p) = z0m(p)
      end if
      z0hv(p)   = z0mv(p)
      z0qv(p)   = z0mv(p)
    end do

    ! Make forcing height a pft-level quantity that is the atmospheric forcing
    ! height plus each pft's z0m+displa
    do pi = 1, max_pft_per_gcell
      do g = lbg, ubg
        if ( pi <= npfts(g) ) then
          p = pfti(g) + pi - 1
          if ( pactive(p) ) then
            l = plandunit(p)
            c = pcolumn(p)
            if ( ityplun(l) == istsoil .or. ityplun(l) == istcrop ) then
              if ( frac_veg_nosno(p) == 0 ) then
                forc_hgt_u_pft(p) = forc_hgt_u(g) + z0mg(c) + displa(p)
                forc_hgt_t_pft(p) = forc_hgt_t(g) + z0mg(c) + displa(p)
                forc_hgt_q_pft(p) = forc_hgt_q(g) + z0mg(c) + displa(p)
              else
                forc_hgt_u_pft(p) = forc_hgt_u(g) + z0m(p) + displa(p)
                forc_hgt_t_pft(p) = forc_hgt_t(g) + z0m(p) + displa(p)
                forc_hgt_q_pft(p) = forc_hgt_q(g) + z0m(p) + displa(p)
              end if
            else if ( ityplun(l) == istwet .or. ityplun(l) == istice ) then
              forc_hgt_u_pft(p) = forc_hgt_u(g) + z0mg(c)
              forc_hgt_t_pft(p) = forc_hgt_t(g) + z0mg(c)
              forc_hgt_q_pft(p) = forc_hgt_q(g) + z0mg(c)
              ! Appropriate momentum roughness length will be
              ! added in SLakeFLuxesMod.
            else if ( ityplun(l) == istdlak ) then
              forc_hgt_u_pft(p) = forc_hgt_u(g)
              forc_hgt_t_pft(p) = forc_hgt_t(g)
              forc_hgt_q_pft(p) = forc_hgt_q(g)
            else if ( ityplun(l) == isturb ) then
              forc_hgt_u_pft(p) = forc_hgt_u(g) + z_0_town(l) + z_d_town(l)
              forc_hgt_t_pft(p) = forc_hgt_t(g) + z_0_town(l) + z_d_town(l)
              forc_hgt_q_pft(p) = forc_hgt_q(g) + z_0_town(l) + z_d_town(l)
            end if
          end if
        end if
      end do
    end do

    do fp = 1, num_nolakep
      p = filter_nolakep(fp)
      c = pcolumn(p)
      thm(p) = forc_t(c) + 0.0098_rk8*forc_hgt_t_pft(p)
    end do
  end subroutine Biogeophysics1

end module mod_clm_biogeophysics1
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
