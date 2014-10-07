module mod_clm_vocemission
  !
  ! Volatile organic compound emission
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_clm_varpar,   only : numpft, nlevcan
  use mod_clm_pftvarcon ,   only : ndllf_evr_tmp_tree,  ndllf_evr_brl_tree,    &
                           ndllf_dcd_brl_tree,  nbrdlf_evr_trp_tree,   &
                           nbrdlf_evr_tmp_tree, nbrdlf_dcd_brl_shrub,  &
                           nbrdlf_dcd_trp_tree, nbrdlf_dcd_tmp_tree,   &
                           nbrdlf_dcd_brl_tree, nbrdlf_evr_shrub,      &
                           nc3_arctic_grass,    nc3crop,               &
                           nc4_grass,           noveg
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_type , only : clm3 , megan_out_type
  use mod_clm_domain , only : ldomain
  use mod_clm_varcon , only : spval
  use mod_clm_megan , only : shr_megan_megcomps_n
  use mod_clm_megan , only : shr_megan_megcomp_t
  use mod_clm_megan , only : shr_megan_linkedlist
  use mod_clm_megan , only : shr_megan_mechcomps_n
  use mod_clm_megan , only : shr_megan_mechcomps
  use mod_clm_megan , only : shr_megan_mapped_emisfctrs
  use mod_clm_meganfactors , only : Agro, Amat, Anew, Aold
  use mod_clm_meganfactors , only : betaT, ct1, ct2, LDF, Ceo

  implicit none

  private

  save

  logical , parameter :: debug = .false.

  public :: VOCEmission
  public :: VOCEmission_init

  contains
  !
  ! Volatile organic compound emission
  ! This code simulates volatile organic compound emissions following
  ! MEGAN (Model of Emissions of Gases and Aerosols from Nature) v2.1
  ! for 20 compound classes. The original description of this
  ! algorithm (for isoprene only) can be found in Guenther et al., 2006
  ! (we follow equations 2-9, 16-17, 20 for explicit canopy).
  ! The model scheme came be described as:
  !    E= epsilon * gamma * rho
  ! VOC flux (E) [ug m-2 h-1] is calculated from baseline emission
  ! factors (epsilon) [ug m-2 h-1] which are specified for each of the 16
  ! CLM PFTs (in input file) OR in the case of isoprene, from
  ! mapped EFs for each PFT which reflect species divergence of emissions,
  ! particularly in North America.
  ! The emission activity factor (gamma) [unitless] for includes
  ! dependence on PPFT, temperature, LAI, leaf age and soil moisture.
  ! For isoprene only we also include the effect of CO2 inhibition as
  ! described by Heald et al., 2009.
  ! The canopy environment constant was calculated offline for CLM+CAM at
  ! standard conditions.
  ! We assume that the escape efficiency (rho) here is unity following
  ! Guenther et al., 2006.
  ! A manuscript describing MEGAN 2.1 and the implementation in CLM is
  ! in preparation: Guenther, Heald et al., 2012
  ! Subroutine written to operate at the patch level.
  !
  ! Input: <filename> to be read in with EFs and some parameters.
  !        Currently these are set in procedure init_EF_params
  ! Output: vocflx(shr_megan_mechcomps_n) !VOC flux [moles/m2/sec]
  !
  subroutine VOCEmission (lbp, ubp, num_soilp, filter_soilp )
    implicit none
    integer(ik4), intent(in) :: lbp, ubp  ! pft bounds
    integer(ik4), intent(in) :: num_soilp ! number of columns in soil pft filter
    integer(ik4), intent(in) :: filter_soilp(num_soilp) ! pft filter for soil
    integer(ik4) , pointer :: pgridcell(:) ! gridcell index of corresponding pft
    integer(ik4) , pointer :: pcolumn(:)   ! column index of corresponding pft
    integer(ik4) , pointer :: ivt(:)       ! pft vegetation type for current
    integer(ik4) , pointer :: clandunit(:) ! gridcell of corresponding column
    integer(ik4) , pointer :: ltype(:)     ! landunit type
    real(rk8), pointer :: t_veg(:)   ! pft vegetation temperature (Kelvin)
    real(rk8), pointer :: fsun(:)    ! sunlit fraction of canopy
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer :: elai(:)
    real(rk8), pointer :: clayfrac(:)      ! fraction of soil that is clay
    real(rk8), pointer :: sandfrac(:)      ! fraction of soil that is sand
    ! direct beam radiation (visible only)
    real(rk8), pointer :: forc_solad(:,:)
    ! diffuse radiation     (visible only)
    real(rk8), pointer :: forc_solai(:,:)
    ! one-sided leaf area index from previous timestep
    real(rk8), pointer :: elai_p(:)
    ! avg pft vegetation temperature for last 24 hrs
    real(rk8), pointer :: t_veg24(:)
    ! avg pft vegetation temperature for last 240 hrs
    real(rk8), pointer :: t_veg240(:)
    ! sunlit fraction of canopy last 24 hrs
    real(rk8), pointer :: fsun24(:)
    ! sunlit fraction of canopy last 240 hrs
    real(rk8), pointer :: fsun240(:)
    ! direct beam radiation last 24hrs  (visible only)
    real(rk8), pointer :: forc_solad24(:)
    ! diffuse radiation  last 24hrs     (visible only)
    real(rk8), pointer :: forc_solai24(:)
    ! direct beam radiation last 240hrs (visible only)
    real(rk8), pointer :: forc_solad240(:)
    ! diffuse radiation  last 240hrs    (visible only)
    real(rk8), pointer :: forc_solai240(:)
    real(rk8), pointer :: h2osoi_vol(:,:) ! volumetric soil water (m3/m3)
    real(rk8), pointer :: h2osoi_ice(:,:) ! ice soil content (kg/m3)
    real(rk8), pointer :: dz(:,:)         ! depth of layer (m)
    real(rk8), pointer :: bsw(:,:)        ! Clapp and Hornberger "b" (nlevgrnd)
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8), pointer :: watsat(:,:)
    real(rk8), pointer :: sucsat(:,:)  ! minimum soil suction (mm) (nlevgrnd)
    real(rk8), pointer :: cisun_z(:,:) ! sunlit intracellular CO2 (Pa)
    real(rk8), pointer :: cisha_z(:,:) ! shaded intracellular CO2 (Pa)
    real(rk8), pointer :: forc_pbot(:) ! atmospheric pressure (Pa)

    real(rk8), pointer :: vocflx(:,:)   ! VOC flux [moles/m2/sec]
    real(rk8), pointer :: vocflx_tot(:) ! VOC flux [moles/m2/sec]

    type(megan_out_type), pointer :: meg_out(:) ! fluxes for CLM history

    real(rk8), pointer :: gamma_out(:)
    real(rk8), pointer :: gammaT_out(:)
    real(rk8), pointer :: gammaP_out(:)
    real(rk8), pointer :: gammaL_out(:)
    real(rk8), pointer :: gammaA_out(:)
    real(rk8), pointer :: gammaS_out(:)
    real(rk8), pointer :: gammaC_out(:)

    real(rk8), pointer :: Eopt_out(:)
    real(rk8), pointer :: topt_out(:)
    real(rk8), pointer :: alpha_out(:)
    real(rk8), pointer :: cp_out(:)
    real(rk8), pointer :: paru_out(:)
    real(rk8), pointer :: par24u_out(:)
    real(rk8), pointer :: par240u_out(:)
    real(rk8), pointer :: para_out(:)
    real(rk8), pointer :: par24a_out(:)
    real(rk8), pointer :: par240a_out(:)
    integer(ik4)  :: fp,p,g,c    ! indices
    real(rk8) :: xepsilon        ! emission factor [ug m-2 h-1]
    real(rk8) :: par_sun         ! temporary
    real(rk8) :: par24_sun       ! temporary
    real(rk8) :: par240_sun      ! temporary
    real(rk8) :: par_sha         ! temporary
    real(rk8) :: par24_sha       ! temporary
    real(rk8) :: par240_sha      ! temporary
    ! activity factor (accounting for light, T, age, LAI conditions)
    real(rk8) :: gamma_x
    real(rk8) :: gamma_p   ! activity factor for PPFD
    real(rk8) :: gamma_l   ! activity factor for PPFD & LAI
    real(rk8) :: gamma_t   ! activity factor for temperature
    real(rk8) :: gamma_a   ! activity factor for leaf age
    real(rk8) :: gamma_sm  ! activity factor for soil moisture
    real(rk8) :: gamma_c   ! activity factor for CO2 (only isoprene)

    integer(ik4) :: class_num, n_meg_comps, imech, imeg, ii

    real(rk8) :: vocflx_meg(shr_megan_megcomps_n)
    type(shr_megan_megcomp_t), pointer :: meg_cmp

    real(rk8) :: cp, alpha,  Eopt, topt  ! for history output

    ! factor used convert MEGAN units [micro-grams/m2/hr] to
    ! CAM srf emis units [g/m2/sec]
    real(rk8), parameter :: megemis_units_factor = 1.D0/3600.D0/1.D6

!    real(rk8) :: root_depth(0:numpft)    ! Root depth [m]
    character(len=32), parameter :: subname = "VOCEmission"
!
!    ! root depth (m) (defined based on Zeng et al., 2001, cf Guenther 2006)
!    ! bare-soil
!    root_depth(noveg) = 0.D0
!    ! evergreen tree
!    root_depth(ndllf_evr_tmp_tree:ndllf_evr_brl_tree) = 1.8D0
!    ! needleleaf deciduous boreal tree
!    root_depth(ndllf_dcd_brl_tree) = 2.0D0
!    ! broadleaf evergreen tree
!    root_depth(nbrdlf_evr_trp_tree:nbrdlf_evr_tmp_tree) = 3.0D0
!    ! broadleaf deciduous tree
!    root_depth(nbrdlf_dcd_trp_tree:nbrdlf_dcd_brl_tree) = 2.0D0
!    ! shrub
!    root_depth(nbrdlf_evr_shrub:nbrdlf_dcd_brl_shrub) = 2.5D0
!    ! grass/crop
!    root_depth(nc3_arctic_grass:numpft) = 1.5D0
!
    if ( shr_megan_mechcomps_n < 1) return

    clandunit  => clm3%g%l%c%landunit
    ltype      => clm3%g%l%itype
    ! Assign local pointers to derived type members (gridcell-level)
    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai
    forc_pbot  => clm_a2l%forc_pbot

    ! Assign local pointers to derived subtypes components (column-level)
    h2osoi_vol       => clm3%g%l%c%cws%h2osoi_vol
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    dz               => clm3%g%l%c%cps%dz
    bsw              => clm3%g%l%c%cps%bsw
    watsat           => clm3%g%l%c%cps%watsat
    sucsat           => clm3%g%l%c%cps%sucsat

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell        => clm3%g%l%c%p%gridcell
    pcolumn          => clm3%g%l%c%p%column
    ivt              => clm3%g%l%c%p%itype
    t_veg            => clm3%g%l%c%p%pes%t_veg
    fsun             => clm3%g%l%c%p%pps%fsun
    elai             => clm3%g%l%c%p%pps%elai
    clayfrac         => clm3%g%l%c%p%pps%clayfrac
    sandfrac         => clm3%g%l%c%p%pps%sandfrac

    cisun_z          => clm3%g%l%c%p%pcf%cisun_z
    cisha_z          => clm3%g%l%c%p%pcf%cisha_z
    if ( nlevcan /= 1 )then
      call fatal(__FILE__,__LINE__,&
              subname//' error: can NOT work without nlevcan == 1' )
    end if

    vocflx           => clm3%g%l%c%p%pvf%vocflx
    vocflx_tot       => clm3%g%l%c%p%pvf%vocflx_tot
    meg_out          => clm3%g%l%c%p%pvf%meg

    gammaL_out       => clm3%g%l%c%p%pvf%gammaL_out
    gammaT_out       => clm3%g%l%c%p%pvf%gammaT_out
    gammaP_out       => clm3%g%l%c%p%pvf%gammaP_out
    gammaA_out       => clm3%g%l%c%p%pvf%gammaA_out
    gammaS_out       => clm3%g%l%c%p%pvf%gammaS_out
    gammaC_out       => clm3%g%l%c%p%pvf%gammaC_out
    gamma_out        => clm3%g%l%c%p%pvf%gamma_out

    Eopt_out         => clm3%g%l%c%p%pvf%Eopt_out
    topt_out         => clm3%g%l%c%p%pvf%topt_out
    alpha_out        => clm3%g%l%c%p%pvf%alpha_out
    cp_out           => clm3%g%l%c%p%pvf%cp_out
    paru_out         => clm3%g%l%c%p%pvf%paru_out
    par24u_out       => clm3%g%l%c%p%pvf%par24u_out
    par240u_out      => clm3%g%l%c%p%pvf%par240u_out
    para_out         => clm3%g%l%c%p%pvf%para_out
    par24a_out       => clm3%g%l%c%p%pvf%par24a_out
    par240a_out      => clm3%g%l%c%p%pvf%par240a_out

    t_veg24          => clm3%g%l%c%p%pvs%t_veg24
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240
    forc_solad24     => clm3%g%l%c%p%pvs%fsd24
    forc_solad240    => clm3%g%l%c%p%pvs%fsd240
    forc_solai24     => clm3%g%l%c%p%pvs%fsi24
    forc_solai240    => clm3%g%l%c%p%pvs%fsi240
    fsun24           => clm3%g%l%c%p%pvs%fsun24
    fsun240          => clm3%g%l%c%p%pvs%fsun240
    elai_p           => clm3%g%l%c%p%pvs%elai_p

    ! initialize variables which get passed to the atmosphere
    vocflx(lbp:ubp,:) = 0.D0
    vocflx_tot(lbp:ubp) = 0.D0

    do imeg=1,shr_megan_megcomps_n
      meg_out(imeg)%flux_out(lbp:ubp) = 0.D0
    end do

    gamma_out(lbp:ubp) = spval
    gammaP_out(lbp:ubp) = spval
    gammaT_out(lbp:ubp) = spval
    gammaA_out(lbp:ubp) = spval
    gammaS_out(lbp:ubp) = spval
    gammaL_out(lbp:ubp) = spval
    gammaC_out(lbp:ubp) = spval

    paru_out(lbp:ubp) = spval
    par24u_out(lbp:ubp) = spval
    par240u_out(lbp:ubp) = spval

    para_out(lbp:ubp) = spval
    par24a_out(lbp:ubp) = spval
    par240a_out(lbp:ubp) = spval

    alpha_out(lbp:ubp) = spval
    cp_out(lbp:ubp) = spval

    topt_out(lbp:ubp) = spval
    Eopt_out(lbp:ubp) = spval

    ! initalize to zero since this might not alway get set
    vocflx_meg(:) = 0.D0

    ! Begin loop over points
    !
    do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      c = pcolumn(p)

      ! initialize EF
      xepsilon = 0.D0

      ! calculate VOC emissions for non-bare ground PFTs
      if (ivt(p) > 0) then
        gamma_x = 0.D0

        ! Calculate PAR: multiply w/m2 by 4.6 to get umol/m2/s for par
        !    (added 8/14/02)
        !------------------------
        ! SUN:
        par_sun = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6D0
        par24_sun = (forc_solad24(p) + fsun24(p) * forc_solai24(p)) * 4.6D0
        par240_sun = (forc_solad240(p) + fsun240(p) * forc_solai240(p)) * 4.6D0
        ! SHADE:
        par_sha = ((1.D0 - fsun(p)) * forc_solai(g,1)) * 4.6D0
        par24_sha = ((1.D0 - fsun24(p)) * forc_solai24(p)) * 4.6D0
        par240_sha = ((1.D0 - fsun240(p)) * forc_solai240(p)) * 4.6D0

        ! Activity factor for LAI (Guenther et al., 2006): all species
        gamma_l = get_gamma_L(fsun240(p), elai(p))

        ! Activity factor for soil moisture: all species
        !  (commented out for now)
!       gamma_sm = get_gamma_SM(clayfrac(p), sandfrac(p), &
!                  h2osoi_vol(c,:), h2osoi_ice(c,:), &
!                  dz(c,:), bsw(c,:), watsat(c,:), sucsat(c,:), &
!                  root_depth(ivt(p)))
        gamma_sm = 1.0D0
        !
        ! Loop through VOCs for light, temperature and leaf age
        ! activity factor & apply all final activity factors to baseline
        ! emission factors
        !

        ! loop over megan compounds
        meg_cmp => shr_megan_linkedlist
        meg_cmp_loop: do while(associated(meg_cmp))
          imeg = meg_cmp%index

          ! set emis factor
          ! if specified, set EF for isoprene with mapped values
          if ( trim(meg_cmp%name) == 'isoprene' .and. &
                  shr_megan_mapped_emisfctrs) then
            xepsilon = get_map_EF(ivt(p),g)
          else
            xepsilon = meg_cmp%emis_factors(ivt(p))
          end if

          class_num = meg_cmp%class_number

          ! Activity factor for PPFD
          gamma_p = get_gamma_P(par_sun, par24_sun, par240_sun, &
                  par_sha, par240_sha, fsun(p), fsun240(p), forc_solad240(p), &
                  forc_solai240(p), LDF(class_num), cp, alpha)
          ! Activity factor for T
          gamma_t = get_gamma_T(t_veg240(p), t_veg24(p),t_veg(p), &
                  ct1(class_num), ct2(class_num),&
                  betaT(class_num),LDF(class_num), Ceo(class_num), Eopt, topt)

          ! Activity factor for Leaf Age
          gamma_a = get_gamma_A(ivt(p), elai_p(p),elai(p),class_num)

          ! Activity factor for CO2 (only for isoprene)
          if (trim(meg_cmp%name) == 'isoprene') then
            gamma_c = get_gamma_C(cisun_z(p,1),cisha_z(p,1), &
                                  forc_pbot(g),fsun(p))
          else
            gamma_c = 1.D0
          end if

          ! Calculate total scaling factor
          gamma_x = gamma_l * gamma_sm * gamma_a * gamma_p * gamma_T * gamma_c

          if ( ( gamma_x >= 0.0D0 ) .and. ( gamma_x < 100.D0 ) ) then

            vocflx_meg(imeg) = xepsilon * gamma_x * &
                    megemis_units_factor / meg_cmp%molec_weight ! moles/m2/sec

            ! assign to arrays for history file output
            ! (not weighted by landfrac)
            meg_out(imeg)%flux_out(p) = meg_out(imeg)%flux_out(p) + &
                     xepsilon * gamma_x * megemis_units_factor*1.D-3 ! Kg/m2/sec

            if (imeg==1) then
              gamma_out(p) = gamma_x
              gammaP_out(p) = gamma_p
              gammaT_out(p) = gamma_t
              gammaA_out(p) = gamma_a
              gammaS_out(p) = gamma_sm
              gammaL_out(p) = gamma_l
              gammaC_out(p) = gamma_c

              paru_out(p) = par_sun
              par24u_out(p) = par24_sun
              par240u_out(p) = par240_sun

              para_out(p) = par_sha
              par24a_out(p) = par24_sha
              par240a_out(p) = par240_sha

              alpha_out(p) = alpha
              cp_out(p) = cp

              topt_out(p) = topt
              Eopt_out(p) = Eopt
            end if
          end if
          if (debug .and. gamma_x > 0.0D0) then
            write(stdout,*) 'MEGAN: n, megan name, epsilon, gamma, vocflx: ', &
                imeg, meg_cmp%name, xepsilon, gamma_x, vocflx_meg(imeg), &
                gamma_p,gamma_t,gamma_a,gamma_sm,gamma_l
          end if
          meg_cmp => meg_cmp%next_megcomp
        end do meg_cmp_loop
        ! sum up the megan compound fluxes for the fluxes of chem
        ! mechanism compounds
        do imech = 1,shr_megan_mechcomps_n
          n_meg_comps = shr_megan_mechcomps(imech)%n_megan_comps
          ! loop over number of megan compounds that make up the
          ! nth mechanism compoud
          do imeg = 1,n_meg_comps
            ii = shr_megan_mechcomps(imech)%megan_comps(imeg)%ptr%index
            vocflx(p,imech) = vocflx(p,imech) + vocflx_meg(ii)
          end do
          vocflx_tot(p) = vocflx_tot(p) + vocflx(p,imech) ! moles/m2/sec
        end do
      end if ! ivt(1:15 only)
    end do ! fp
  end subroutine VOCEmission
  !
  ! Interface to set all input parameters for 20 VOC compound classes.
  ! including EFs for 16(+1 bare ground) PFTs.
  ! For now set all specified values, in future to be replaced with
  ! values read in from file.
  ! (heald, 04/27/11)
  !
  subroutine VOCEmission_init(  )
    use mod_clm_megan , only : shr_megan_factors_file
    use mod_clm_meganfactors , only : megan_factors_init , megan_factors_get
    implicit none
    type(shr_megan_megcomp_t), pointer :: meg_cmp
    integer(ik4)  :: class_num
    real(rk8) :: factors(numpft)
    real(rk8) :: molec_wght

    if ( shr_megan_mechcomps_n < 1) return
    call megan_factors_init( shr_megan_factors_file )
    meg_cmp => shr_megan_linkedlist
    do while(associated(meg_cmp))
      allocate(meg_cmp%emis_factors(numpft))
      call megan_factors_get( trim(meg_cmp%name), factors, &
                             class_num, molec_wght )
      meg_cmp%emis_factors = factors
      meg_cmp%class_number = class_num
      meg_cmp%molec_weight = molec_wght
      meg_cmp => meg_cmp%next_megcomp
    end do
  end subroutine VOCEmission_init
  !
  ! Get mapped EF for isoprene
  ! Use gridded values for 6 PFTs specified by MEGAN following
  ! Guenther et al. (2006).  Map the numpft CLM PFTs to these 6.
  ! Units: [ug m-2 h-1]
  !
  function get_map_EF(ivt_in,g_in)
    use mod_clm_type
    implicit none
    integer(ik4) , intent(in) :: ivt_in
    integer(ik4) , intent(in) :: g_in
    real(rk8) :: get_map_EF

    ! emission factors for isoprene for each pft [ug m-2 h-1]
    real(rk8) , pointer :: efisop(:,:)

    ! assign local pointer
    efisop     => clm3%g%gve%efisop

    get_map_EF = 0.D0
    if ( ivt_in == ndllf_evr_tmp_tree .or. &
         ivt_in == ndllf_evr_brl_tree) then     !fineleaf evergreen
      get_map_EF = efisop(2,g_in)
    else if (ivt_in == ndllf_dcd_brl_tree) then     !fineleaf deciduous
      get_map_EF = efisop(3,g_in)
    else if (ivt_in >= nbrdlf_evr_trp_tree .and. &
             ivt_in <= nbrdlf_dcd_brl_tree) then    !broadleaf trees
      get_map_EF = efisop(1,g_in)
    else if (ivt_in >= nbrdlf_evr_shrub .and. &
             ivt_in <= nbrdlf_dcd_brl_shrub) then   !shrubs
      get_map_EF = efisop(4,g_in)
    else if (ivt_in >= nc3_arctic_grass .and. &
             ivt_in <= nc4_grass) then              !grass
      get_map_EF = efisop(5,g_in)
    else if (ivt_in >= nc3crop) then                !crops
      get_map_EF =efisop(6,g_in)
    end if
  end function get_map_EF
  !
  ! Activity factor for PPFD (Guenther et al., 2006):
  !               all light dependent species
  !-------------------------
  ! With distinction between sunlit and shaded leafs, weight scalings by
  ! fsun and fshade
  ! Scale total incident par by fraction of sunlit leaves (added on 1/2002)
  ! fvitt -- forc_solad240, forc_solai240 can be zero when CLM finidat is
  !          specified which will cause par240 to be zero and produce NaNs
  !          via log(par240)
  ! dml   -- fsun240 can be equal to or greater than one before 10 day
  !          averages are set on startup or if a new pft comes online
  !          during land cover change.
  ! Avoid this problem by only doing calculations with fsun240 when fsun240 is
  ! between 0 and 1
  !
  real(rk8) function get_gamma_P(par_sun_in,par24_sun_in,par240_sun_in, &
                  par_sha_in,par240_sha_in,fsun_in,fsun240_in, &
                  forc_solad240_in,forc_solai240_in,LDF_in,cp,alpha)
    implicit none
    real(rk8) , intent(in) :: par_sun_in
    real(rk8) , intent(in) :: par24_sun_in
    real(rk8) , intent(in) :: par240_sun_in
    real(rk8) , intent(in) :: par_sha_in
    real(rk8) , intent(in) :: par240_sha_in
    real(rk8) , intent(in) :: fsun_in
    real(rk8) , intent(in) :: fsun240_in
    real(rk8) , intent(in) :: forc_solad240_in
    real(rk8) , intent(in) :: forc_solai240_in
    real(rk8) , intent(in) :: LDF_in
    real(rk8) , intent(out) :: cp                      ! temporary
    real(rk8) , intent(out) :: alpha                   ! temporary
    real(rk8) :: gamma_p_LDF             ! activity factor for PPFD

    real(rk8), parameter :: ca1 = 0.004D0  ! empirical coefficent for alpha
    real(rk8), parameter :: ca2 = 0.0005D0 ! empirical coefficent for alpha
    real(rk8), parameter :: ca3 = 0.0468D0 ! empirical coefficent for cp
    ! std conditions for past 24 hrs [umol/m2/s]
    real(rk8), parameter :: par0_sun = 200.D0
    ! std conditions for past 24 hrs [umol/m2/s]
    real(rk8), parameter :: par0_shade = 50.D0
    real(rk8), parameter :: alpha_fix = 0.001D0 ! empirical coefficient
    real(rk8), parameter :: cp_fix = 1.21D0     ! empirical coefficient

    if ( (fsun240_in > 0.D0) .and. (fsun240_in < 1.D0) .and. &
         (forc_solad240_in > 0.D0) .and. (forc_solai240_in > 0.D0)) then
      ! With alpha and cp calculated based on eq 6 and 7:
      ! Note indexing for accumulated variables is all at pft level
      ! SUN:
      alpha = ca1 - ca2 * log(par240_sun_in)
      cp = ca3 * exp(ca2 * (par24_sun_in-par0_sun))*par240_sun_in**(0.6D0)
      gamma_p_LDF = fsun_in * ( cp * alpha * par_sun_in * &
              (1.D0 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5D0) )
      ! SHADE:
      alpha = ca1 - ca2 * log(par240_sha_in)
      cp = ca3 * exp(ca2 * (par_sha_in-par0_shade))*par240_sha_in**(0.6D0)
      gamma_p_LDF = gamma_p_LDF + (1.D0-fsun_in) * (cp*alpha*par_sha_in* &
              (1.D0 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5D0))
    else
      ! With fixed alpha and cp (from MEGAN User's Guide):
      ! SUN: direct + diffuse
      alpha = alpha_fix
      cp = cp_fix
      gamma_p_LDF = fsun_in * ( cp * alpha*par_sun_in * &
              (1.D0 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5D0) )
      ! SHADE: diffuse
      gamma_p_LDF = gamma_p_LDF + (1.D0-fsun_in) * (cp*alpha*par_sha_in * &
              (1.D0 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5D0))
    end if
    ! Calculate total activity factor for PPFD accounting for
    ! light-dependent fraction
    get_gamma_P = (1.D0 - LDF_in) + LDF_in * gamma_p_LDF
  end function get_gamma_P
  !
  ! Activity factor for LAI (Guenther et al., 2006): all species
  ! Guenther et al., 2006 eq 3
  !
  function get_gamma_L(fsun240_in,elai_in)
    use mod_clm_varcon , only : denice
    use mod_clm_varpar , only : nlevsoi
    implicit none
    real(rk8) , intent(in) :: fsun240_in
    real(rk8) , intent(in) :: elai_in
    real(rk8) :: get_gamma_L             ! return value
    ! parameters
    ! factor to set emissions to unity @ std
    real(rk8) , parameter :: cce = 0.30D0
    ! same as Cce but for non-accumulated vars
    real(rk8) , parameter :: cce1 = 0.24D0

    if ( (fsun240_in > 0.0D0) .and. (fsun240_in < 1.D30) ) then
      get_gamma_L = cce * elai_in
    else
      get_gamma_L = cce1 * elai_in
    end if
  end function get_gamma_L
  !
  ! Activity factor for soil moisture (Guenther et al., 2006): all species
  !----------------------------------
  ! Calculate the mean scaling factor throughout the root depth.
  ! wilting point potential is in units of matric potential (mm)
  ! (1 J/Kg = 0.001 MPa, approx = 0.1 m)
  ! convert to volumetric soil water using equation 7.118 of the
  ! CLM4 Technical Note
  !
  function get_gamma_SM(clayfrac_in,sandfrac_in,h2osoi_vol_in,h2osoi_ice_in, &
                  dz_in,bsw_in,watsat_in,sucsat_in,root_depth_in)
    use mod_clm_varcon , only : denice
    use mod_clm_varpar , only : nlevsoi
    implicit none
    real(rk8) , intent(in) :: clayfrac_in
    real(rk8) , intent(in) :: sandfrac_in
    real(rk8) , intent(in) :: h2osoi_vol_in(nlevsoi)
    real(rk8) , intent(in) :: h2osoi_ice_in(nlevsoi)
    real(rk8) , intent(in) :: dz_in(nlevsoi)
    real(rk8) , intent(in) :: bsw_in(nlevsoi)
    real(rk8) , intent(in) :: watsat_in(nlevsoi)
    real(rk8) , intent(in) :: sucsat_in(nlevsoi)
    real(rk8) , intent(in) :: root_depth_in
    real(rk8) :: get_gamma_SM
    ! local variables
    integer(ik4) :: j
    real(rk8) :: nl                      ! temporary number of soil levels
    real(rk8) :: theta_ice               ! water content in ice in m3/m3
    real(rk8) :: wilt                    ! wilting point in m3/m3
    real(rk8) :: theta1                  ! temporary

    ! parameters
    real(rk8), parameter :: deltheta1=0.06D0 ! empirical coefficient
    real(rk8), parameter :: smpmax = 2.57D5  ! maximum soil matrix potential

    if ((clayfrac_in > 0) .and. (sandfrac_in > 0)) then
      get_gamma_SM = 0.D0
      nl=0.D0
      do j = 1,nlevsoi
        if  (sum(dz_in(1:j)) < root_depth_in)  then
          theta_ice = h2osoi_ice_in(j)/(dz_in(j)*denice)
          wilt = ((smpmax/sucsat_in(j))**(-1.D0/bsw_in(j))) * &
                  (watsat_in(j) - theta_ice)
          theta1 = wilt + deltheta1
          if (h2osoi_vol_in(j) >= theta1) then
            get_gamma_SM = get_gamma_SM + 1.D0
          else if ( (h2osoi_vol_in(j) > wilt) .and. &
                    (h2osoi_vol_in(j) < theta1) ) then
            get_gamma_SM = get_gamma_SM + &
                    ( h2osoi_vol_in(j) - wilt ) / deltheta1
          else
            get_gamma_SM = get_gamma_SM + 0.D0
          end if
          nl = nl+1.D0
        end if
      end do
      if (nl > 0.D0) then
        get_gamma_SM = get_gamma_SM/nl
      end if
      if (get_gamma_SM > 1.0D0) then
        write(stdout,*) 'healdSM > 1: gamma_SM, nl', get_gamma_SM, nl
        get_gamma_SM=1.0D0
      end if
    else
      get_gamma_SM = 1.0D0
    end if
  end function get_gamma_SM
  !
  ! Activity factor for temperature
  !--------------------------------
  ! Calculate both a light-dependent fraction as in Guenther et al., 2006
  ! for isoprene of a max saturation type form. Also caculate a
  ! light-independent fraction of the form of an exponential.
  ! Final activity factor depends on light dependent fraction
  ! of compound type.
  !
  function get_gamma_T(t_veg240_in,t_veg24_in,t_veg_in,ct1_in,ct2_in, &
                       betaT_in,LDF_in,Ceo_in,Eopt,topt)
    implicit none
    ! varibles in
    real(rk8) , intent(in) :: t_veg240_in
    real(rk8) , intent(in) :: t_veg24_in
    real(rk8) , intent(in) :: t_veg_in
    real(rk8) , intent(in) :: ct1_in
    real(rk8) , intent(in) :: ct2_in
    real(rk8) , intent(in) :: betaT_in
    real(rk8) , intent(in) :: LDF_in
    real(rk8) , intent(in) :: Ceo_in
    real(rk8) , intent(out) :: Eopt                    ! temporary
    real(rk8) , intent(out) :: topt                    ! temporary

    ! local variables
    real(rk8) :: get_gamma_T
    real(rk8) :: gamma_t_LDF             ! activity factor for temperature
    real(rk8) :: gamma_t_LIF             ! activity factor for temperature
    real(rk8) :: x                       ! temporary

    ! parameters
    real(rk8), parameter :: co1 = 313.D0  ! empirical coefficient
    real(rk8), parameter :: co2 = 0.6D0   ! empirical coefficient
    real(rk8), parameter :: co4 = 0.05D0  ! empirical coefficient
    real(rk8), parameter :: tstd0 = 297D0 ! std temperature [K]
    real(rk8), parameter :: topt_fix = 317.D0 ! std temperature [K]
    real(rk8), parameter :: Eopt_fix = 2.26D0 ! empirical coefficient
    ! empirical coefficient (0.0083 in User's Guide)
    real(rk8), parameter :: ct3 = 0.00831D0
    real(rk8), parameter :: tstd = 303.15D0  ! std temperature [K]

    ! Light dependent fraction (Guenther et al., 2006)
    if ( (t_veg240_in > 0.0D0) .and. (t_veg240_in < 1.D30) ) then
      ! topt and Eopt from eq 8 and 9:
      topt = co1 + (co2 * (t_veg240_in-tstd0))
      Eopt = Ceo_in * exp (co4 * (t_veg24_in-tstd0)) * &
                      exp (co4 * (t_veg240_in -tstd0))
    else
      topt = topt_fix
      Eopt = Eopt_fix
    end if
    x = ( (1.D0/topt) - (1.D0/(t_veg_in)) ) / ct3
    gamma_t_LDF = Eopt * ( ct2_in * exp(ct1_in * x) / &
                          (ct2_in - ct1_in * (1.D0 - exp(ct2_in * x))) )
    ! Light independent fraction (of exp(beta T) form)
    gamma_t_LIF = exp(betaT_in * (t_veg_in - tstd))
    ! Calculate total activity factor for light as a function of
    ! light-dependent fraction
    get_gamma_T = (1-LDF_in)*gamma_T_LIF + LDF_in*gamma_T_LDF
  end function get_gamma_T
  !
  ! Activity factor for leaf age (Guenther et al., 2006)
  !-----------------------------
  ! If not CNDV elai is constant therefore gamma_a=1.0
  ! gamma_a set to unity for evergreens (PFTs 1, 2, 4, 5)
  ! Note that we assume here that the time step is shorter than the number of
  !days after budbreak required to induce isoprene emissions (ti=12 days) and
  ! the number of days after budbreak to reach peak emission (tm=28 days)
  !
  function get_gamma_A(ivt_in, elai_p_in,elai_in,nclass_in)
    implicit none
    ! varibles in
    integer(ik4) , intent(in)  :: ivt_in
    integer(ik4) , intent(in)  :: nclass_in
    real(rk8) , intent(in) :: elai_p_in
    real(rk8) , intent(in) :: elai_in
    real(rk8) :: get_gamma_A
    ! local variables
    real(rk8) :: elai_prev ! lai for previous timestep
    ! fractions of leaves at different phenological stages
    real(rk8) :: fnew , fgro , fmat , fold
    if ( (ivt_in == ndllf_dcd_brl_tree) .or. &
         (ivt_in >= nbrdlf_dcd_trp_tree) ) then  ! non-evergreen
      if ( (elai_p_in > 0.0D0) .and. (elai_p_in < 1.D30) )then
        ! have accumulated average lai over last timestep
        elai_prev = 2.D0*elai_p_in-elai_in
        if (elai_prev == elai_in) then
          fnew = 0.0D0
          fgro = 0.0D0
          fmat = 1.0D0
          fold = 0.0D0
        else if (elai_prev > elai_in) then
          fnew = 0.0D0
          fgro = 0.0D0
          fmat = 1.0D0 - (elai_prev - elai_in)/elai_prev
          fold = (elai_prev - elai_in)/elai_prev
        else if (elai_prev < elai_in) then
          fnew = 1 - (elai_prev / elai_in)
          fgro = 0.0D0
          fmat = (elai_prev / elai_in)
          fold = 0.0D0
        end if
        get_gamma_A = fnew*Anew(nclass_in) + fgro*Agro(nclass_in) + &
                fmat*Amat(nclass_in) + fold*Aold(nclass_in)
      else
        get_gamma_A = 1.0D0
      end if
    else
      get_gamma_A = 1.0D0
    end if
  end function get_gamma_A
  !
  ! Activity factor for instantaneous CO2 changes (Heald et al., 2009)
  !-------------------------
  ! With distinction between sunlit and shaded leafs, weight scalings by
  ! fsun and fshade
  !
  function get_gamma_C(cisun_in,cisha_in,forc_pbot_in,fsun_in)
    ! corresponds to CCSM_CO2_PPMV set in env_conf.xml
    use mod_clm_varctl , only : co2_ppmv
    implicit none
    real(rk8) , intent(in) :: cisun_in
    real(rk8) , intent(in) :: cisha_in
    real(rk8) , intent(in) :: forc_pbot_in
    real(rk8) , intent(in) :: fsun_in
    real(rk8) :: get_gamma_C
    real(rk8) :: IEmin   ! empirical coeff for CO2
    real(rk8) :: IEmax   ! empirical coeff for CO2
    real(rk8) :: ECi50   ! empirical coeff for CO2
    real(rk8) :: Cislope ! empirical coeff for CO2
    real(rk8) :: fint    ! interpolation fraction for CO2
    ! temporary sunlight/shade weighted cisun & cisha (umolCO2/mol)
    real(rk8) :: ci

    ! Determine long-term CO2 growth environment (ie. ambient CO2)
    ! and interpolate parameters
    if ( co2_ppmv < 400.D0 ) then
      IEmin   = 0.7301D0
      IEmax   = 1.0542D0
      ECi50   = 457.D0
      Cislope = 3.1665D0
    else if ( (co2_ppmv > 400.D0) .and. (co2_ppmv < 500.D0) ) then
      fint    = (co2_ppmv - 400.D0)/100.D0
      IEmin   = fint*0.7301D0 + (1.D0 - fint)*0.7034D0
      IEmax   = fint*1.0542D0 + (1.D0 - fint)*0.9897D0
      ECi50   = fint*457.D0 + (1.D0 - fint)*472.D0
      Cislope = fint*3.1665D0 + (1.D0 - fint)*3.0652D0
    else if ( (co2_ppmv > 500.D0) .and. (co2_ppmv < 600.D0) ) then
      fint = (co2_ppmv - 500.D0)/100.D0
      IEmin   = fint*0.7034D0 + (1.D0 - fint)*0.6768D0
      IEmax   = fint*0.9897D0 + (1.D0 - fint)*0.9253D0
      ECi50   = fint*472.D0 + (1.D0 - fint)*488.D0
      Cislope = fint*3.0652D0 + (1.D0 - fint)*2.9321D0
    else if ( (co2_ppmv > 600.D0) .and. (co2_ppmv < 700.D0) ) then
      fint = (co2_ppmv - 600.D0)/100.D0
      IEmin   = fint*0.6768D0 + (1.D0 - fint)*0.6500D0
      IEmax   = fint*0.9253D0 + (1.D0 - fint)*0.8611D0
      ECi50   = fint*488.D0 + (1.D0 - fint)*508.D0
      Cislope = fint*2.9321D0 + (1.D0 - fint)*2.7497D0
    else if ( (co2_ppmv > 700.D0) .and. (co2_ppmv < 800.D0) ) then
      fint = (co2_ppmv - 700.D0)/100.D0
      IEmin   = fint*0.6500D0 + (1.D0 - fint)*0.6063D0
      IEmax   = fint*0.8611D0 + (1.D0 - fint)*0.7976D0
      ECi50   = fint*508.D0 + (1.D0 - fint)*575.D0
      Cislope = fint*2.7497D0 + (1.D0 - fint)*2.3643D0
    else if ( co2_ppmv > 800.D0 ) then
      IEmin   = 0.6063D0
      IEmax   = 0.7976D0
      ECi50   = 575.D0
      Cislope = 2.3643D0
    end if

    ! Intracellular CO2 concentrations (ci) given in Pa, divide by atmos
    ! pressure to get mixing ratio (umolCO2/mol)
    if ( (cisun_in > 0.D0) .and. (cisha_in > 0.D0) .and. &
         (forc_pbot_in > 0.D0) .and. (fsun_in > 0.D0) ) then
      ci = ( fsun_in*cisun_in + (1.D0-fsun_in)*cisha_in )/forc_pbot_in * 1.D6
      get_gamma_C = IEmin + ( (IEmax-IEmin)/(1.D0+(ci/ECi50)**Cislope) )
    else if ( (cisha_in > 0.D0) .and. (forc_pbot_in > 0.D0) ) then
      ci = cisha_in/forc_pbot_in * 1.D6
      get_gamma_C = IEmin + ( (IEmax-IEmin)/(1.D0+(ci/ECi50)**Cislope) )
    else if ( (cisun_in > 0.D0) .and. (forc_pbot_in > 0.D0) ) then
      ci = cisun_in/forc_pbot_in * 1.D6
      get_gamma_C = IEmin + ( (IEmax-IEmin)/(1.D0+(ci/ECi50)**Cislope) )
    else
      get_gamma_C = 1.D0
    end if
  end function get_gamma_C

end module mod_clm_vocemission
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
