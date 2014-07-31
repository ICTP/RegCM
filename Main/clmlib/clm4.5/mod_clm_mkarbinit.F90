module mod_clm_mkarbinit
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mppparam
  use mod_dynparam
  use mod_clm_atmlnd
  use mod_clm_type
  use mod_clm_varpar , only : nlevsoi , nlevgrnd , nlevsno , nlevlak , nlevurb
  use mod_clm_varcon , only : bdsno , istice , istwet , istsoil , isturb , &
       denice , denh2o , spval , sb , icol_road_perv, icol_road_imperv ,   &
       icol_roof , icol_sunwall , icol_shadewall , tfrz
  use mod_clm_varcon , only : istcrop
  use mod_clm_varcon , only : h2osno_max
  use mod_clm_varctl , only : pertlim
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_snicar , only : snw_rds_min
  use mod_bats_param

  implicit none

  private

  save

  public mkarbinit   ! Make arbitrary initial conditions
  public mkregcminit ! Get initial conditions from the regcm model
  public perturbIC   ! Perturb the initial conditions by pertlim

  contains
  !
  ! Initializes the following time varying variables:
  ! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
  ! snow       : snow_depth, snl, dz, z, zi
  ! temperature: t_soisno, t_veg, t_grnd
  !
  subroutine mkarbinit()
    implicit none
    ! column index associated with each pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! column type
    integer(ik4) , pointer , dimension(:) :: ctype
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! landunit index associated with each pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! true => landunit is an urban point
    logical , pointer , dimension(:) :: urbpoi
    ! true => landunit is not vegetated
    logical , pointer , dimension(:) :: ifspecial
    ! layer thickness depth (m)
    real(rk8) , pointer , dimension(:,:) :: dz
    ! volumetric soil water at saturation (porosity)
    real(rk8) , pointer , dimension(:,:) :: watsat
    ! ice lens (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    ! liquid water (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    ! Clapp and Hornberger "b"
    real(rk8) , pointer , dimension(:,:) :: bsw
    ! minimum soil suction (mm)
    real(rk8) , pointer , dimension(:,:) :: sucsat
    ! interface level below a "z" level (m)
    real(rk8) , pointer , dimension(:,:) :: zi
    ! water in the unconfined aquifer (mm)
    real(rk8) , pointer , dimension(:) :: wa
    ! water table depth (m)
    real(rk8) , pointer , dimension(:) :: zwt
    ! surface water (mm)
    real(rk8) , pointer , dimension(:) :: h2osfc
    ! surface water temperature
    real(rk8) , pointer , dimension(:) :: t_h2osfc
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_h2osfc
    !surface water runoff (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_surf
    ! frost table depth (m)
    real(rk8) , pointer , dimension(:) :: frost_table
    ! perched water table depth (m)
    real(rk8) , pointer , dimension(:) :: zwt_perched
    ! integrated snowfall
    real(rk8) , pointer , dimension(:) :: int_snow
    ! snow melt (net)
    real(rk8) , pointer , dimension(:) :: qflx_snow_melt
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    ! lake temperature (Kelvin)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: t_lake
    ! ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    ! vegetation temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_veg
    ! 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m
    ! Urban 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_u
    ! Rural 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_r
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    ! canopy water (mm H2O) (column-level)
    real(rk8) , pointer , dimension(:) :: h2ocan_col
    ! canopy water (mm H2O) (pft-level)
    real(rk8) , pointer , dimension(:) :: h2ocan_pft
    ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    ! snow height (m)
    real(rk8) , pointer , dimension(:) :: snow_depth
    ! irrigation flux (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_irrig
    ! emitted infrared (longwave) radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_out
    ! soil water potential in each soil layer (MPa)
    real(rk8) , pointer , dimension(:,:) :: soilpsi
    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(rk8) , pointer , dimension(:,:) :: snw_rds
    ! snow grain size, top (col) [microns]
    real(rk8) , pointer , dimension(:) :: snw_rds_top
    ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(rk8) , pointer , dimension(:) :: sno_liq_top
    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcpho
    ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcphi
    ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bctot
    ! total mass of BC in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_col
    ! total mass of BC in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_top
    ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcphi
    ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcpho
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocpho
    ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocphi
    ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_octot
    ! total mass of OC in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_col
    ! total mass of OC in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_top
    ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocphi
    ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocpho
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
    ! total mass of dust in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_col
    ! total mass of dust in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_top
    ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst1
    ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst2
    ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst3
    ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst4
    ! current irrigation rate [mm/s]
    real(rk8) , pointer , dimension(:) :: irrig_rate
    ! soil T for top 0.17 m
    real(rk8) , pointer , dimension(:) :: tsoi17
    !fractional area with water table at surface
    real(rk8) , pointer , dimension(:) :: fsat
    ! number of time steps for which we still need to irrigate today
    ! (if 0, ignore irrig_rate)
    integer(ik4) ,  pointer , dimension(:) :: n_irrig_steps_left

    integer(ik4) :: j , l , c , p  ! indices
    integer(ik4) :: nlevs       ! number of levels
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    real(rk8) :: vwc , psi      ! for calculating soilpsi

    if ( myid == italk )then
      write(stdout,*) 'Setting initial data to non-spun up values'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    h2osfc           => clm3%g%l%c%cws%h2osfc
    t_h2osfc         => clm3%g%l%c%ces%t_h2osfc
    frac_h2osfc      => clm3%g%l%c%cps%frac_h2osfc
    qflx_h2osfc_surf => clm3%g%l%c%cwf%qflx_h2osfc_surf
    qflx_snow_melt   => clm3%g%l%c%cwf%qflx_snow_melt
    frost_table      => clm3%g%l%c%cws%frost_table
    zwt_perched      => clm3%g%l%c%cws%zwt_perched
    int_snow         => clm3%g%l%c%cws%int_snow
    ltype            => clm3%g%l%itype
    lakpoi           => clm3%g%l%lakpoi
    ifspecial        => clm3%g%l%ifspecial
    urbpoi           => clm3%g%l%urbpoi

    ! Assign local pointers to derived subtypes components (column-level)

    ctype            => clm3%g%l%c%itype
    clandunit        => clm3%g%l%c%landunit
    snl              => clm3%g%l%c%cps%snl
    dz               => clm3%g%l%c%cps%dz
    watsat           => clm3%g%l%c%cps%watsat
    sucsat           => clm3%g%l%c%cps%sucsat
    bsw              => clm3%g%l%c%cps%bsw
    soilpsi          => clm3%g%l%c%cps%soilpsi
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol       => clm3%g%l%c%cws%h2osoi_vol
    h2ocan_col       => clm3%g%l%c%cws%pws_a%h2ocan
    qflx_irrig       => clm3%g%l%c%cwf%qflx_irrig
    snow_depth       => clm3%g%l%c%cps%snow_depth
    h2osno           => clm3%g%l%c%cws%h2osno
    t_soisno         => clm3%g%l%c%ces%t_soisno
    t_lake           => clm3%g%l%c%ces%t_lake
    t_grnd           => clm3%g%l%c%ces%t_grnd
    tsoi17           => clm3%g%l%c%ces%tsoi17
    zi               => clm3%g%l%c%cps%zi
    wa               => clm3%g%l%c%cws%wa
    zwt              => clm3%g%l%c%cws%zwt
    fsat             => clm3%g%l%c%cws%fsat
    snw_rds          => clm3%g%l%c%cps%snw_rds
    snw_rds_top      => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top      => clm3%g%l%c%cps%sno_liq_top
    mss_bcpho        => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi        => clm3%g%l%c%cps%mss_bcphi
    mss_bctot        => clm3%g%l%c%cps%mss_bctot
    mss_bc_col       => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top       => clm3%g%l%c%cps%mss_bc_top
    mss_cnc_bcphi    => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho    => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_ocpho        => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi        => clm3%g%l%c%cps%mss_ocphi
    mss_octot        => clm3%g%l%c%cps%mss_octot
    mss_oc_col       => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top       => clm3%g%l%c%cps%mss_oc_top
    mss_cnc_ocphi    => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho    => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_dst1         => clm3%g%l%c%cps%mss_dst1
    mss_dst2         => clm3%g%l%c%cps%mss_dst2
    mss_dst3         => clm3%g%l%c%cps%mss_dst3
    mss_dst4         => clm3%g%l%c%cps%mss_dst4
    mss_dsttot       => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col      => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top      => clm3%g%l%c%cps%mss_dst_top
    mss_cnc_dst1     => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2     => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3     => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4     => clm3%g%l%c%cps%mss_cnc_dst4
    n_irrig_steps_left => clm3%g%l%c%cps%n_irrig_steps_left
    irrig_rate       => clm3%g%l%c%cps%irrig_rate

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    h2ocan_pft     => clm3%g%l%c%p%pws%h2ocan
    t_veg          => clm3%g%l%c%p%pes%t_veg
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    t_ref2m_u      => clm3%g%l%c%p%pes%t_ref2m_u
    t_ref2m_r      => clm3%g%l%c%p%pes%t_ref2m_r
    plandunit      => clm3%g%l%c%p%landunit
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, and snow_depth has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0.D0

       ! added for canopy water mass balance under dynamic pft weights
       !clm3%g%l%c%p%pps%tlai(p) = 0.D0
       !clm3%g%l%c%p%pps%tsai(p) = 0.D0
       !clm3%g%l%c%p%pps%elai(p) = 0.D0
       !clm3%g%l%c%p%pps%esai(p) = 0.D0
       !clm3%g%l%c%p%pps%htop(p) = 0.D0
       !clm3%g%l%c%p%pps%hbot(p) = 0.D0
       !clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = 0
    end do

    ! initialize h2osfc, frac_h2osfc, t_h2osfc, qflx_snow_melt
    do c = begc,endc
       h2osfc(c)           = 0.D0
       frac_h2osfc(c)      = 0.D0
       !       t_h2osfc(c) = spval
       t_h2osfc(c)         = 274.D0
       qflx_h2osfc_surf(c) = 0.D0
       qflx_snow_melt(c)   = 0.D0
    enddo

    do c = begc , endc
      ! canopy water (column level)
      h2ocan_col(c) = 0.D0

      ! snow water

      l = clandunit(c)

      if ( ltype(l) == istice ) then
        h2osno(c) = h2osno_max
      else
        h2osno(c) = 0.D0
      end if

      ! initialize int_snow, int_melt
      int_snow(c) = h2osno(c)
      ! snow depth

      snow_depth(c)  = h2osno(c) / bdsno

      ! Initialize Irrigation to zero
      if ( ltype(l) == istsoil ) then
        n_irrig_steps_left(c) = 0
        irrig_rate(c)         = 0.0D0
      end if
    end do

    ! Set snow layer number, depth and thickiness

    call snow_depth2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    ! NOTE: THESE MEMORY COPIES ARE INEFFICIENT --
    !     SINCE nlev LOOP IS NESTED FIRST!!!!
    do c = begc , endc

      t_soisno(c,-nlevsno+1:nlevgrnd) = spval
      t_lake(c,1:nlevlak) = spval

      l = clandunit(c)
      if ( .not. lakpoi(l) ) then  !not lake
        t_soisno(c,-nlevsno+1:0) = spval
        if ( snl(c) < 0 ) then    !snow layer temperatures
          do j = snl(c)+1 , 0
            t_soisno(c,j) = 250.D0
          end do
        end if
        if ( ltype(l)==istice ) then
          do j = 1 , nlevgrnd
            t_soisno(c,j) = 250.D0
          end do
        else if (ltype(l) == istwet) then
          do j = 1 , nlevgrnd
            t_soisno(c,j) = 277.D0
          end do
        else if (ltype(l) == isturb) then
#if (defined VANCOUVER)
          if (ctype(c) == icol_road_perv .or. &
              ctype(c) == icol_road_imperv) then
            ! Set road top layer to initial air temperature and interpolate
            ! other layers down to 20C in bottom layer
            do j = 1 , nlevgrnd
              t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1))
            end do
          ! Set wall and roof layers to initial air temperature
          else if (ctype(c) == icol_sunwall .or. &
                   ctype(c) == icol_shadewall .or. &
                   ctype(c) == icol_roof) then
            do j = 1 , nlevurb
              t_soisno(c,j) = 297.56
            end do
          else
            do j = 1 , nlevgrnd
              t_soisno(c,j) = 283.D0
            end do
          end if
#elif (defined MEXICOCITY)
          if (ctype(c) == icol_road_perv .or. &
              ctype(c) == icol_road_imperv) then
            ! Set road top layer to initial air temperature and interpolate
            ! other layers down to 22C in bottom layer
            do j = 1 , nlevgrnd
              t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1))
            end do
          ! Set wall and roof layers to initial air temperature
          else if (ctype(c) == icol_sunwall .or. &
                   ctype(c) == icol_shadewall .or. &
                   ctype(c) == icol_roof) then
            do j = 1 , nlevurb
              t_soisno(c,j) = 289.46
            end do
          else
            do j = 1 , nlevgrnd
              t_soisno(c,j) = 283.D0
            end do
          end if
#else
          if (ctype(c) == icol_road_perv .or. &
              ctype(c) == icol_road_imperv) then
            do j = 1 , nlevgrnd
              t_soisno(c,j) = 274.D0
            end do
          ! Set sunwall, shadewall, roof to fairly high temperature to
          ! avoid initialization shock from large heating/air conditioning flux
          else if (ctype(c) == icol_sunwall .or. &
                   ctype(c) == icol_shadewall .or. &
                   ctype(c) == icol_roof) then
            do j = 1 , nlevurb
              t_soisno(c,j) = 292.D0
            end do
          end if
#endif
        else
          do j = 1 , nlevgrnd
            t_soisno(c,j) = 274.D0
          end do
        end if
        t_grnd(c) = t_soisno(c,snl(c)+1)
      else                     !lake
        t_lake(c,1:nlevlak) = 277.D0
        t_grnd(c) = t_lake(c,1)
      end if
      tsoi17(c) = t_grnd(c)
    end do

    call perturbIC( clm3%g%l )

    do p = begp , endp
      c = pcolumn(p)
      l = plandunit(p)

      ! Initialize Irrigation to zero
      if (ltype(l)==istsoil) then
        qflx_irrig(c)      = 0.0D0
      end if

#if (defined VANCOUVER)
      t_veg(p) = 297.56
      t_ref2m(p) = 297.56
      if (urbpoi(l)) then
        t_ref2m_u(p) = 297.56
      else
        t_ref2m_u(p) = spval
      end if
      if (ifspecial(l)) then
        t_ref2m_r(p) = spval
      else
        t_ref2m_r(p) = 297.56
      end if
#elif (defined MEXICOCITY)
      t_veg(p) = 289.46
      t_ref2m(p) = 289.46
      if (urbpoi(l)) then
        t_ref2m_u(p) = 289.46
      else
        t_ref2m_u(p) = spval
      end if
      if (ifspecial(l)) then
        t_ref2m_r(p) = spval
      else
        t_ref2m_r(p) = 289.46
      end if
#else
      t_veg(p) = 283.D0
      t_ref2m(p) = 283.D0
      if (urbpoi(l)) then
        t_ref2m_u(p) = 283.D0
      else
        t_ref2m_u(p) = spval
      end if
      if (ifspecial(l)) then
        t_ref2m_r(p) = spval
      else
        t_ref2m_r(p) = 283.D0
      end if
#endif
      eflx_lwrad_out(p) = sb * (t_grnd(c))**4
    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are
    ! obtained NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have
    ! valid values over soil and urban pervious road
    ! (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:) = spval

    wa(begc:endc)  = 5000.D0
    zwt(begc:endc) = 0.D0

    do c = begc , endc
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then  !not lake
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_road_perv ) then
            wa(c)  = 4800.D0
            ! One meter below soil column
            zwt(c) = (25.D0 + zi(c,nlevsoi)) - wa(c)/0.2D0 /1000.D0
          else
            wa(c)  = spval
            zwt(c) = spval
          end if
          ! initialize frost_table, zwt_perched
          zwt_perched(c) = spval
          frost_table(c) = spval
        else
          wa(c)  = 4000.D0
          ! One meter below soil column
          ! initialize frost_table, zwt_perched to bottom of soil column
          zwt(c) = (25.D0 + zi(c,nlevsoi)) - wa(c)/0.2D0 /1000.D0
          zwt_perched(c) = zi(c,nlevsoi)
          frost_table(c) = zi(c,nlevsoi)
        end if
      end if
    end do

    do c = begc , endc
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then  !not lake

        ! volumetric water
        if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( j > nlevsoi ) then
              h2osoi_vol(c,j) = 0.0D0
            else
              h2osoi_vol(c,j) = 0.15D0
            end if
          end do
        else if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_road_perv ) then
            nlevs = nlevgrnd
            do j = 1 , nlevs
              if ( j <= nlevsoi ) then
                h2osoi_vol(c,j) = 0.3D0
              else
                h2osoi_vol(c,j) = 0.0D0
              end if
            end do
          else if ( ctype(c) == icol_road_imperv ) then
            nlevs = nlevgrnd
            do j = 1 , nlevs
              h2osoi_vol(c,j) = 0.0D0
            end do
          else
            nlevs = nlevurb
            do j = 1 , nlevs
              h2osoi_vol(c,j) = 0.0D0
            end do
          end if
        else if ( ltype(l) == istwet ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( j > nlevsoi ) then
              h2osoi_vol(c,j) = 0.0D0
            else
              h2osoi_vol(c,j) = 1.0D0
            end if
          end do
        else if ( ltype(l) == istice ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            h2osoi_vol(c,j) = 1.0D0
          end do
        endif
        do j = 1 , nlevs
          h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
          ! soil layers
          if ( t_soisno(c,j) <= tfrz ) then
            h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
            h2osoi_liq(c,j) = 0.D0
          else
            h2osoi_ice(c,j) = 0.D0
            h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
          end if
        end do

#if (defined CN)
        ! soil water potential (added 10/21/03, PET)
        ! required for CN code
        if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( h2osoi_liq(c,j) > 0.D0 ) then
              vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
              psi = sucsat(c,j) * (-9.8D-6) * &
                      (vwc/watsat(c,j))**(-bsw(c,j))  ! Mpa
              soilpsi(c,j) = max(psi, -15.0D0)
              soilpsi(c,j) = min(soilpsi(c,j),0.0D0)
            end if
          end do
        end if
        fsat(c)   = 0.0D0
#endif
      end if
    end do

    ! Set snow
    do j = -nlevsno+1 , 0
      do c = begc , endc
        l = clandunit(c)
        if ( .not. lakpoi(l) ) then  !not lake
          if ( j > snl(c) ) then
            h2osoi_ice(c,j) = dz(c,j)*250.D0
            h2osoi_liq(c,j) = 0.D0
          end if
        end if
      end do
    end do
    ! initialize SNICAR fields:
    do c = begc ,endc
      mss_bctot(c,:) = 0.D0
      mss_bcpho(c,:) = 0.D0
      mss_bcphi(c,:) = 0.D0
      mss_cnc_bcphi(c,:)=0.D0
      mss_cnc_bcpho(c,:)=0.D0

      mss_octot(c,:) = 0.D0
      mss_ocpho(c,:) = 0.D0
      mss_ocphi(c,:) = 0.D0
      mss_cnc_ocphi(c,:)=0.D0
      mss_cnc_ocpho(c,:)=0.D0

      mss_dst1(c,:) = 0.D0
      mss_dst2(c,:) = 0.D0
      mss_dst3(c,:) = 0.D0
      mss_dst4(c,:) = 0.D0
      mss_dsttot(c,:) = 0.D0
      mss_cnc_dst1(c,:)=0.D0
      mss_cnc_dst2(c,:)=0.D0
      mss_cnc_dst3(c,:)=0.D0
      mss_cnc_dst4(c,:)=0.D0

      if (snl(c) < 0) then
        snw_rds(c,snl(c)+1:0)        = snw_rds_min
        snw_rds(c,-nlevsno+1:snl(c)) = 0.D0
        snw_rds_top(c)               = snw_rds_min
        sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / &
                (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
      else if (h2osno(c) > 0.D0) then
        snw_rds(c,0)             = snw_rds_min
        snw_rds(c,-nlevsno+1:-1) = 0.D0
        snw_rds_top(c)           = spval
        sno_liq_top(c)           = spval
      else
        snw_rds(c,:)   = 0.D0
        snw_rds_top(c) = spval
        sno_liq_top(c) = spval
      end if
    end do
  end subroutine mkarbinit
  !
  ! Initializes the following time varying variables:
  ! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
  ! snow       : snow_depth, snl, dz, z, zi
  ! temperature: t_soisno, t_veg, t_grnd
  !
  subroutine mkregcminit( adomain )
    implicit none
    type(atm_domain) , intent(in) :: adomain
    ! column index associated with each pft
    integer(ik4) , pointer , dimension(:) :: pcolumn
    ! column type
    integer(ik4) , pointer , dimension(:) :: ctype
    ! gridcell index associated with each column
    integer(ik4) , pointer , dimension(:) :: cgridcell
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! landunit index associated with each pft
    integer(ik4) , pointer , dimension(:) :: plandunit
    ! true => landunit is an urban point
    logical , pointer , dimension(:) :: urbpoi
    ! true => landunit is not vegetated
    logical , pointer , dimension(:) :: ifspecial
    ! layer thickness depth (m)
    real(rk8) , pointer , dimension(:,:) :: dz
    ! volumetric soil water at saturation (porosity)
    real(rk8) , pointer , dimension(:,:) :: watsat
    ! ice lens (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_ice
    ! liquid water (kg/m2)
    real(rk8) , pointer , dimension(:,:) :: h2osoi_liq
    ! Clapp and Hornberger "b"
    real(rk8) , pointer , dimension(:,:) :: bsw
    ! minimum soil suction (mm)
    real(rk8) , pointer , dimension(:,:) :: sucsat
    ! interface level below a "z" level (m)
    real(rk8) , pointer , dimension(:,:) :: zi
    ! water in the unconfined aquifer (mm)
    real(rk8) , pointer , dimension(:) :: wa
    ! water table depth (m)
    real(rk8) , pointer , dimension(:) :: zwt
    ! surface water (mm)
    real(rk8) , pointer , dimension(:) :: h2osfc
    ! surface water temperature
    real(rk8) , pointer , dimension(:) :: t_h2osfc
    ! fraction of ground covered by surface water (0 to 1)
    real(rk8) , pointer , dimension(:) :: frac_h2osfc
    !surface water runoff (mm/s)
    real(rk8) , pointer , dimension(:) :: qflx_h2osfc_surf
    ! frost table depth (m)
    real(rk8) , pointer , dimension(:) :: frost_table
    ! perched water table depth (m)
    real(rk8) , pointer , dimension(:) :: zwt_perched
    ! integrated snowfall
    real(rk8) , pointer , dimension(:) :: int_snow
    ! snow melt (net)
    real(rk8) , pointer , dimension(:) :: qflx_snow_melt
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    ! lake temperature (Kelvin)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: t_lake
    ! ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    ! vegetation temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_veg
    ! 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m
    ! Urban 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_u
    ! Rural 2 m height surface air temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_ref2m_r
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(rk8) , pointer , dimension(:,:) :: h2osoi_vol
    ! canopy water (mm H2O) (column-level)
    real(rk8) , pointer , dimension(:) :: h2ocan_col
    ! canopy water (mm H2O) (pft-level)
    real(rk8) , pointer , dimension(:) :: h2ocan_pft
    ! snow water (mm H2O)
    real(rk8) , pointer , dimension(:) :: h2osno
    ! snow height (m)
    real(rk8) , pointer , dimension(:) :: snow_depth
    ! irrigation flux (mm H2O/s)
    real(rk8) , pointer , dimension(:) :: qflx_irrig
    ! emitted infrared (longwave) radiation (W/m**2)
    real(rk8) , pointer , dimension(:) :: eflx_lwrad_out
    ! soil water potential in each soil layer (MPa)
    real(rk8) , pointer , dimension(:,:) :: soilpsi
    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(rk8) , pointer , dimension(:,:) :: snw_rds
    ! snow grain size, top (col) [microns]
    real(rk8) , pointer , dimension(:) :: snw_rds_top
    ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(rk8) , pointer , dimension(:) :: sno_liq_top
    ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcpho
    ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bcphi
    ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_bctot
    ! total mass of BC in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_col
    ! total mass of BC in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_bc_top
    ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcphi
    ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_bcpho
    ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocpho
    ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_ocphi
    ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(rk8) , pointer , dimension(:,:) :: mss_octot
    ! total mass of OC in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_col
    ! total mass of OC in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_oc_top
    ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocphi
    ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_ocpho
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
    ! total mass of dust in snow column (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_col
    ! total mass of dust in top snow layer (col) [kg]
    real(rk8) , pointer , dimension(:) :: mss_dst_top
    ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst1
    ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst2
    ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst3
    ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    real(rk8) , pointer , dimension(:,:) :: mss_cnc_dst4
    ! current irrigation rate [mm/s]
    real(rk8) , pointer , dimension(:) :: irrig_rate
    ! soil T for top 0.17 m
    real(rk8) , pointer , dimension(:) :: tsoi17
    !fractional area with water table at surface
    real(rk8) , pointer , dimension(:) :: fsat
    ! number of time steps for which we still need to irrigate today
    ! (if 0, ignore irrig_rate)
    integer(ik4) ,  pointer , dimension(:) :: n_irrig_steps_left

    integer(ik4) :: j , l , c , p  , g ! indices
    integer(ik4) :: nlevs       ! number of levels
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    real(rk8) :: vwc , psi      ! for calculating soilpsi

    if ( myid == italk )then
      write(stdout,*) 'Setting initial data to non-spun up values'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    h2osfc           => clm3%g%l%c%cws%h2osfc
    t_h2osfc         => clm3%g%l%c%ces%t_h2osfc
    frac_h2osfc      => clm3%g%l%c%cps%frac_h2osfc
    qflx_h2osfc_surf => clm3%g%l%c%cwf%qflx_h2osfc_surf
    qflx_snow_melt   => clm3%g%l%c%cwf%qflx_snow_melt
    frost_table      => clm3%g%l%c%cws%frost_table
    zwt_perched      => clm3%g%l%c%cws%zwt_perched
    int_snow         => clm3%g%l%c%cws%int_snow
    ltype            => clm3%g%l%itype
    lakpoi           => clm3%g%l%lakpoi
    ifspecial        => clm3%g%l%ifspecial
    urbpoi           => clm3%g%l%urbpoi

    ! Assign local pointers to derived subtypes components (column-level)

    ctype            => clm3%g%l%c%itype
    clandunit        => clm3%g%l%c%landunit
    cgridcell        => clm3%g%l%c%gridcell
    snl              => clm3%g%l%c%cps%snl
    dz               => clm3%g%l%c%cps%dz
    watsat           => clm3%g%l%c%cps%watsat
    sucsat           => clm3%g%l%c%cps%sucsat
    bsw              => clm3%g%l%c%cps%bsw
    soilpsi          => clm3%g%l%c%cps%soilpsi
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol       => clm3%g%l%c%cws%h2osoi_vol
    h2ocan_col       => clm3%g%l%c%cws%pws_a%h2ocan
    qflx_irrig       => clm3%g%l%c%cwf%qflx_irrig
    snow_depth       => clm3%g%l%c%cps%snow_depth
    h2osno           => clm3%g%l%c%cws%h2osno
    t_soisno         => clm3%g%l%c%ces%t_soisno
    t_lake           => clm3%g%l%c%ces%t_lake
    t_grnd           => clm3%g%l%c%ces%t_grnd
    tsoi17           => clm3%g%l%c%ces%tsoi17
    zi               => clm3%g%l%c%cps%zi
    wa               => clm3%g%l%c%cws%wa
    zwt              => clm3%g%l%c%cws%zwt
    fsat             => clm3%g%l%c%cws%fsat
    snw_rds          => clm3%g%l%c%cps%snw_rds
    snw_rds_top      => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top      => clm3%g%l%c%cps%sno_liq_top
    mss_bcpho        => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi        => clm3%g%l%c%cps%mss_bcphi
    mss_bctot        => clm3%g%l%c%cps%mss_bctot
    mss_bc_col       => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top       => clm3%g%l%c%cps%mss_bc_top
    mss_cnc_bcphi    => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho    => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_ocpho        => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi        => clm3%g%l%c%cps%mss_ocphi
    mss_octot        => clm3%g%l%c%cps%mss_octot
    mss_oc_col       => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top       => clm3%g%l%c%cps%mss_oc_top
    mss_cnc_ocphi    => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho    => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_dst1         => clm3%g%l%c%cps%mss_dst1
    mss_dst2         => clm3%g%l%c%cps%mss_dst2
    mss_dst3         => clm3%g%l%c%cps%mss_dst3
    mss_dst4         => clm3%g%l%c%cps%mss_dst4
    mss_dsttot       => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col      => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top      => clm3%g%l%c%cps%mss_dst_top
    mss_cnc_dst1     => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2     => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3     => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4     => clm3%g%l%c%cps%mss_cnc_dst4
    n_irrig_steps_left => clm3%g%l%c%cps%n_irrig_steps_left
    irrig_rate       => clm3%g%l%c%cps%irrig_rate

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    h2ocan_pft     => clm3%g%l%c%p%pws%h2ocan
    t_veg          => clm3%g%l%c%p%pes%t_veg
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    t_ref2m_u      => clm3%g%l%c%p%pes%t_ref2m_u
    t_ref2m_r      => clm3%g%l%c%p%pes%t_ref2m_r
    plandunit      => clm3%g%l%c%p%landunit
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, and snow_depth has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0.D0

       ! added for canopy water mass balance under dynamic pft weights
       !clm3%g%l%c%p%pps%tlai(p) = 0.D0
       !clm3%g%l%c%p%pps%tsai(p) = 0.D0
       !clm3%g%l%c%p%pps%elai(p) = 0.D0
       !clm3%g%l%c%p%pps%esai(p) = 0.D0
       !clm3%g%l%c%p%pps%htop(p) = 0.D0
       !clm3%g%l%c%p%pps%hbot(p) = 0.D0
       !clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = 0
    end do

    ! initialize h2osfc, frac_h2osfc, t_h2osfc, qflx_snow_melt
    do c = begc , endc
      g = cgridcell(c)

      h2osfc(c)           = 0.D0
      frac_h2osfc(c)      = 0.D0
      t_h2osfc(c)         = adomain%tgrd(g)
      qflx_h2osfc_surf(c) = 0.D0
      qflx_snow_melt(c)   = 0.D0
    enddo

    do c = begc , endc
      g = cgridcell(c)
      l = clandunit(c)

      ! canopy water (column level)
      h2ocan_col(c) = 0.D0

      if ( .not. lakpoi(l) ) then  !not lake
        if ( ltype(l) == istice ) then
          if ( adomain%tgrd(g) < 268.0D0 ) then
            h2osno(c) = h2osno_max
          else
            h2osno(c) = 0.0D0
          end if
        else if ( ltype(l) /= isturb ) then
          h2osno(c) = adomain%snow(g)
        else
          h2osno(c) = 0.0D0
        end if

        if ( h2osno(c) < 0.1D0 ) then
          if ( ltype(l) /= isturb ) then
            ! Start with some snow on mountains and on cold regions
            if ( adomain%tgrd(g) < 263.0D0 ) then
              h2osno(c) = adomain%topo(g)/10.0
            end if
          end if
        end if
      end if

      ! initialize int_snow, int_melt
      int_snow(c) = h2osno(c)
      ! snow depth
      snow_depth(c) = h2osno(c) / bdsno
      ! Initialize Irrigation to zero
      if ( ltype(l) == istsoil ) then
        n_irrig_steps_left(c) = 0
        irrig_rate(c)         = 0.0D0
      end if
    end do

    ! Set snow layer number, depth and thickiness

    call snow_depth2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    ! NOTE: THESE MEMORY COPIES ARE INEFFICIENT --
    !     SINCE nlev LOOP IS NESTED FIRST!!!!
    do c = begc , endc

      t_soisno(c,-nlevsno+1:nlevgrnd) = spval
      t_lake(c,1:nlevlak) = spval

      g = cgridcell(c)
      l = clandunit(c)

      if ( .not. lakpoi(l) ) then  !not lake
        t_soisno(c,-nlevsno+1:0) = spval
        if ( snl(c) < 0 ) then    !snow layer temperatures
          do j = snl(c)+1 , 0
            t_soisno(c,j) = adomain%tgrd(g) - 1.0d0 * dble(snl(c)-j)
          end do
        end if
        if ( ltype(l) == istice ) then
          do j = 1 , nlevgrnd
            t_soisno(c,j) = 250.0D0
          end do
        else if (ltype(l) == istwet) then
          do j = 1 , nlevgrnd
            t_soisno(c,j) = adomain%tgrd(g)
          end do
        else if (ltype(l) == isturb) then
          if (ctype(c) == icol_road_perv .or. &
              ctype(c) == icol_road_imperv) then
            do j = 1 , nlevgrnd
              t_soisno(c,j) = adomain%tgrd(g) + 2.0D0
            end do
          ! Set sunwall, shadewall, roof to fairly high temperature to
          ! avoid initialization shock from large heating/air conditioning flux
          else if (ctype(c) == icol_sunwall .or. &
                   ctype(c) == icol_shadewall .or. &
                   ctype(c) == icol_roof) then
            do j = 1 , nlevurb
              t_soisno(c,j) = adomain%tgrd(g) + 5.0D0
            end do
          end if
        else
          do j = 1 , nlevgrnd
            t_soisno(c,j) = adomain%tgrd(g)
          end do
        end if
        t_grnd(c) = adomain%tgrd(g)
      else  ! lake
        t_lake(c,1:nlevlak) = adomain%tgrd(g)
        t_grnd(c) = t_lake(c,1)
      end if
      tsoi17(c) = t_grnd(c)
    end do

    do p = begp , endp
      c = pcolumn(p)
      l = plandunit(p)

      ! Initialize Irrigation to zero
      if ( ltype(l) == istsoil ) then
        qflx_irrig(c) = 0.0D0
      end if

      t_veg(p) = t_grnd(c)
      t_ref2m(p) = t_grnd(c)
      if ( urbpoi(l) ) then
        t_ref2m_u(p) = t_grnd(c)
      else
        t_ref2m_u(p) = spval
      end if
      if ( ifspecial(l) ) then
        t_ref2m_r(p) = spval
      else
        t_ref2m_r(p) = t_grnd(c)
      end if
      eflx_lwrad_out(p) = sb * (t_grnd(c))**4
    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are
    ! obtained NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have
    ! valid values over soil and urban pervious road
    ! (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:) = spval

    wa(begc:endc)  = 5000.D0
    zwt(begc:endc) = 0.D0

    do c = begc , endc
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then  !not lake
        if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_road_perv ) then
            wa(c)  = 4800.D0
            ! One meter below soil column
            zwt(c) = (25.D0 + zi(c,nlevsoi)) - wa(c)/0.2D0 /1000.D0
          else
            wa(c)  = spval
            zwt(c) = spval
          end if
          ! initialize frost_table, zwt_perched
          zwt_perched(c) = spval
          frost_table(c) = spval
        else
          wa(c)  = 4000.D0
          ! One meter below soil column
          ! initialize frost_table, zwt_perched to bottom of soil column
          zwt(c) = (25.D0 + zi(c,nlevsoi)) - wa(c)/0.2D0 /1000.D0
          zwt_perched(c) = zi(c,nlevsoi)
          frost_table(c) = zi(c,nlevsoi)
        end if
      end if
    end do

    do c = begc , endc
      g = cgridcell(c)
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then  !not lake

        ! volumetric water
        if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( j > nlevsoi ) then
              h2osoi_vol(c,j) = 0.0D0
            else
              ! h2osoi_vol(c,j) = 0.15D0
              ! h2osoi_vol(c,j) = watsat(c,j)*0.10D0
              h2osoi_vol(c,j) = slmo(adomain%luse(g)) * &
                            xmopor(iexsol(adomain%luse(g)))
            end if
          end do
        else if ( ltype(l) == isturb ) then
          if ( ctype(c) == icol_road_perv ) then
            nlevs = nlevgrnd
            do j = 1 , nlevs
              if ( j <= nlevsoi ) then
                h2osoi_vol(c,j) = watsat(c,j)*0.30D0
                !h2osoi_vol(c,j) = 0.3D0
              else
                h2osoi_vol(c,j) = 0.0D0
              end if
            end do
          else if ( ctype(c) == icol_road_imperv ) then
            nlevs = nlevgrnd
            do j = 1 , nlevs
              h2osoi_vol(c,j) = 0.0D0
            end do
          else
            nlevs = nlevurb
            do j = 1 , nlevs
              h2osoi_vol(c,j) = 0.0D0
            end do
          end if
        else if ( ltype(l) == istwet ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( j > nlevsoi ) then
              h2osoi_vol(c,j) = 0.0D0
            else
              h2osoi_vol(c,j) = 1.0D0
            end if
          end do
        else if ( ltype(l) == istice ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            h2osoi_vol(c,j) = 1.0D0
          end do
        endif
        do j = 1 , nlevs
          h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
          ! soil layers
          if ( t_soisno(c,j) <= tfrz ) then
            h2osoi_ice(c,j) = dz(c,j)*denice*h2osoi_vol(c,j)
            h2osoi_liq(c,j) = 0.D0
          else
            h2osoi_ice(c,j) = 0.D0
            h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
          end if
        end do

#if (defined CN)
        ! soil water potential (added 10/21/03, PET)
        ! required for CN code
        if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
          nlevs = nlevgrnd
          do j = 1 , nlevs
            if ( h2osoi_liq(c,j) > 0.D0 ) then
              vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
              psi = sucsat(c,j) * (-9.8D-6) * &
                      (vwc/watsat(c,j))**(-bsw(c,j))  ! Mpa
              soilpsi(c,j) = max(psi, -15.0D0)
              soilpsi(c,j) = min(soilpsi(c,j),0.0D0)
            end if
          end do
        end if
        fsat(c)   = 0.0D0
#endif
      end if
    end do

    ! Set snow
    do j = -nlevsno+1 , 0
      do c = begc , endc
        l = clandunit(c)
        if ( .not. lakpoi(l) ) then  !not lake
          if ( j > snl(c) ) then
            h2osoi_ice(c,j) = dz(c,j)*250.D0
            h2osoi_liq(c,j) = 0.D0
          end if
        end if
      end do
    end do
    ! initialize SNICAR fields:
    do c = begc ,endc
      mss_bctot(c,:) = 0.D0
      mss_bcpho(c,:) = 0.D0
      mss_bcphi(c,:) = 0.D0
      mss_cnc_bcphi(c,:)=0.D0
      mss_cnc_bcpho(c,:)=0.D0

      mss_octot(c,:) = 0.D0
      mss_ocpho(c,:) = 0.D0
      mss_ocphi(c,:) = 0.D0
      mss_cnc_ocphi(c,:) = 0.D0
      mss_cnc_ocpho(c,:) = 0.D0

      mss_dst1(c,:) = 0.D0
      mss_dst2(c,:) = 0.D0
      mss_dst3(c,:) = 0.D0
      mss_dst4(c,:) = 0.D0
      mss_dsttot(c,:) = 0.D0
      mss_cnc_dst1(c,:) = 0.D0
      mss_cnc_dst2(c,:) = 0.D0
      mss_cnc_dst3(c,:) = 0.D0
      mss_cnc_dst4(c,:) = 0.D0

      if (snl(c) < 0) then
        snw_rds(c,snl(c)+1:0)        = snw_rds_min
        snw_rds(c,-nlevsno+1:snl(c)) = 0.D0
        snw_rds_top(c)               = snw_rds_min
        sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / &
                (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
      else if (h2osno(c) > 0.D0) then
        snw_rds(c,0)             = snw_rds_min
        snw_rds(c,-nlevsno+1:-1) = 0.D0
        snw_rds_top(c)           = spval
        sno_liq_top(c)           = spval
      else
        snw_rds(c,:)   = 0.D0
        snw_rds_top(c) = spval
        sno_liq_top(c) = spval
      end if
    end do
  end subroutine mkregcminit
  !
  ! Perturbs initial conditions by the amount in the namelist variable pertlim.
  !
  subroutine perturbIC( landunit )
    implicit none
    type(landunit_type) , intent(inout) :: landunit
    integer(ik4) :: j , l , c   ! indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    real(rk8) :: pertval  ! for calculating temperature perturbation
    integer(ik4) :: nlevs ! number of levels
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: ltype
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    ! lake temperature (Kelvin)  (1:nlevlak)
    real(rk8) , pointer , dimension(:,:) :: t_lake
    ! ground temperature (Kelvin)
    real(rk8) , pointer , dimension(:) :: t_grnd
    integer(ik4) , pointer , dimension(:) :: ctype ! column type

    if ( pertlim /= 0.0D0 )then
      if ( myid == italk ) then
        write(stdout,*) 'Applying perturbation to initial soil temperature'
      end if

      clandunit  => landunit%c%landunit
      lakpoi     => landunit%lakpoi
      ltype      => landunit%itype
      t_soisno   => landunit%c%ces%t_soisno
      t_lake     => landunit%c%ces%t_lake
      t_grnd     => landunit%c%ces%t_grnd
      snl        => landunit%c%cps%snl
      ctype      => clm3%g%l%c%itype

      ! Determine subgrid bounds on this processor

      call get_proc_bounds(begc=begc,endc=endc)

      do c = begc , endc
        l = clandunit(c)
        if ( .not. lakpoi(l) ) then  !not lake
          if ( ltype(l) == isturb ) then
            if ( ctype(c) == icol_road_perv .or. &
                 ctype(c) == icol_road_imperv ) then
              nlevs = nlevgrnd
            else
              nlevs = nlevurb
            end if
          else
            nlevs = nlevgrnd
          end if
          ! Randomly perturb soil temperature
          do j = 1 , nlevs
            call random_number (pertval)
            pertval       = 2.D0*pertlim*(0.5D0 - pertval)
            t_soisno(c,j) = t_soisno(c,j)*(1.D0 + pertval)
          end do
        else                       !lake
          ! Randomly perturb lake temperature
          do j = 1 , nlevlak
            call random_number (pertval)
            pertval     = 2.D0*pertlim*(0.5D0 - pertval)
            t_lake(c,j) = t_lake(c,j)*(1.D0 + pertval)
          end do
        end if
        ! Randomly perturb surface ground temp
        call random_number (pertval)
        pertval   = 2.D0*pertlim*(0.5D0 - pertval)
        t_grnd(c) = t_grnd(c)*(1.D0 + pertval)
      end do
    end if
  end subroutine perturbIC
  !
  ! Create snow layers and interfaces given snow depth.
  ! Note that cps%zi(0) is set in routine iniTimeConst.
  !
  subroutine snow_depth2lev(lbc, ubc)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc ! column bounds
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! snow height (m)
    real(rk8) , pointer , dimension(:) :: snow_depth
    ! true => landunit is a lake point
    logical , pointer , dimension(:) :: lakpoi
    ! number of snow layers
    integer(ik4) , pointer , dimension(:) :: snl
    ! layer depth  (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: z
    ! layer thickness depth (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: dz
    ! interface depth (m) over snow only
    real(rk8) , pointer , dimension(:,:) :: zi
    integer(ik4) :: c , l , j  !indices

    ! Assign local pointers to derived subtypes components (landunit-level)

    lakpoi => clm3%g%l%lakpoi

    ! Assign local pointers to derived type members (column-level)

    clandunit => clm3%g%l%c%landunit
    snow_depth    => clm3%g%l%c%cps%snow_depth
    snl       => clm3%g%l%c%cps%snl
    zi        => clm3%g%l%c%cps%zi
    dz        => clm3%g%l%c%cps%dz
    z         => clm3%g%l%c%cps%z

    ! Initialize snow levels and interfaces (lake and non-lake points)

    do c = lbc , ubc
      dz(c,-nlevsno+1: 0) = 1.D36
      z (c,-nlevsno+1: 0) = 1.D36
      zi(c,-nlevsno  :-1) = 1.D36
    end do

    ! Determine snow levels and interfaces for non-lake points

    do c = lbc , ubc
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then
        if (snow_depth(c) < 0.01D0) then
          snl(c) = 0
          dz(c,-nlevsno+1:0) = 0.D0
          z (c,-nlevsno+1:0) = 0.D0
          zi(c,-nlevsno+0:0) = 0.D0
        else
          if ( (snow_depth(c) >= 0.01D0) .and. &
               (snow_depth(c) <= 0.03D0) ) then
            snl(c) = -1
            dz(c,0)  = snow_depth(c)
          else if ( (snow_depth(c) > 0.03D0) .and. &
                    (snow_depth(c) <= 0.04D0) ) then
            snl(c) = -2
            dz(c,-1) = snow_depth(c)/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.04D0) .and. &
                    (snow_depth(c) <= 0.07D0) ) then
            snl(c) = -2
            dz(c,-1) = 0.02D0
            dz(c, 0) = snow_depth(c) - dz(c,-1)
          else if ( (snow_depth(c) > 0.07D0) .and. &
                    (snow_depth(c) <= 0.12D0) ) then
            snl(c) = -3
            dz(c,-2) = 0.02D0
            dz(c,-1) = (snow_depth(c) - 0.02D0)/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.12D0) .and. &
                    (snow_depth(c) <= 0.18D0) ) then
            snl(c) = -3
            dz(c,-2) = 0.02D0
            dz(c,-1) = 0.05D0
            dz(c, 0) = snow_depth(c) - dz(c,-2) - dz(c,-1)
          else if ( (snow_depth(c) > 0.18D0) .and. &
                    (snow_depth(c) <= 0.29D0) ) then
            snl(c) = -4
            dz(c,-3) = 0.02D0
            dz(c,-2) = 0.05D0
            dz(c,-1) = (snow_depth(c) - dz(c,-3) - dz(c,-2))/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( (snow_depth(c) > 0.29D0) .and. &
                    (snow_depth(c) <= 0.41D0) ) then
            snl(c) = -4
            dz(c,-3) = 0.02D0
            dz(c,-2) = 0.05D0
            dz(c,-1) = 0.11D0
            dz(c, 0) = snow_depth(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
          else if ( (snow_depth(c) > 0.41D0) .and. &
                    (snow_depth(c) <= 0.64D0) ) then
            snl(c) = -5
            dz(c,-4) = 0.02D0
            dz(c,-3) = 0.05D0
            dz(c,-2) = 0.11D0
            dz(c,-1) = (snow_depth(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2.D0
            dz(c, 0) = dz(c,-1)
          else if ( snow_depth(c) > 0.64D0 ) then
            snl(c) = -5
            dz(c,-4) = 0.02D0
            dz(c,-3) = 0.05D0
            dz(c,-2) = 0.11D0
            dz(c,-1) = 0.23D0
            dz(c, 0)=snow_depth(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
          end if
        end if
      end if
    end do

    ! The following loop is currently not vectorized

    do c = lbc , ubc
      l = clandunit(c)
      if ( .not. lakpoi(l) ) then
        do j = 0 , snl(c)+1 , -1
          z(c,j)    = zi(c,j) - 0.5D0*dz(c,j)
          zi(c,j-1) = zi(c,j) - dz(c,j)
        end do
      end if
    end do

    ! Determine snow levels and interfaces for lake points

    do c = lbc , ubc
      l = clandunit(c)
      if ( lakpoi(l) ) then
        snl(c) = 0
        dz(c,-nlevsno+1:0) = 0.D0
        z (c,-nlevsno+1:0) = 0.D0
        zi(c,-nlevsno+0:0) = 0.D0
      end if
    end do
  end subroutine snow_depth2lev

end module mod_clm_mkarbinit
