module mod_clm_typeinit
  !
  ! Allocate clmtype components and initialize them to signaling NaN.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_intkinds
  use mod_mpmessage
  use mod_clm_type
  use mod_clm_varpar
#if (defined VICHYDRO)
  use mod_clm_varpar
#endif
  use mod_clm_varctl
  use mod_clm_decomp , only : get_proc_bounds , get_proc_global
  use mod_clm_varcon , only : spval , ispval
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_megan , only : shr_megan_megcomps_n
  use mod_clm_drydep , only : n_drydep , drydep_method , DD_XLND
  use mod_clm_varpar , only : ngases

  implicit none

  private

  public :: initClmtype

  contains
  !
  ! Initialize clmtype components to signaling nan
  ! The following clmtype components should NOT be initialized here
  ! since they are set in routine clm_map which is called before this
  ! routine is invoked
  !    *%area, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
  !    *%ifspecial, *%ityplun, *%itype
  !    *%pfti, *%pftf, *%pftn
  !    *%coli, *%colf, *%coln
  !    *%luni, *%lunf, *%lunn
  !
  subroutine initClmtype()
    implicit none
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numg        ! total number of gridcells across all proc
    integer(ik4) :: numl        ! total number of landunits across all proc
    integer(ik4) :: numc        ! total number of columns across all proc
    integer(ik4) :: nump        ! total number of pfts across all proc
    character(len=32) , parameter :: subname = "initClmtype"

    ! Determine necessary indices

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call get_proc_global(numg,numl,numc,nump)

    call init_pft_type     (begp, endp, clm3%g%l%c%p)
    call init_column_type  (begc, endc, clm3%g%l%c)
    call init_landunit_type(begl, endl, clm3%g%l)
    call init_gridcell_type(begg, endg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    call init_decomp_cascade_constants()

#if (defined CNDV)
    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()
#endif

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, clm3%g%l%c%p%pebal)
    call init_energy_balance_type(begc, endc, clm3%g%l%c%cebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, clm3%g%l%c%p%pwbal)
    call init_water_balance_type(begc, endc, clm3%g%l%c%cwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, clm3%g%l%c%p%pcbal)
    call init_carbon_balance_type(begc, endc, clm3%g%l%c%ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, clm3%g%l%c%p%pnbal)
    call init_nitrogen_balance_type(begc, endc, clm3%g%l%c%cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
    call init_pft_pstate_type(begc, endc, clm3%g%l%c%cps%pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)

    !pft photosynthesis relevant variables
    call init_pft_psynstate_type(begp, endp, clm3%g%l%c%p%ppsyns)
#if (defined CNDV)
    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, clm3%g%l%c%p%pdgvs)
#endif
#if (defined CNDV)
    call init_pft_pdgvstate_type(begc, endc, clm3%g%l%c%cdgvs%pdgvs_a)
#endif
    call init_pft_vstate_type(begp, endp, clm3%g%l%c%p%pvs)

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
    call init_pft_estate_type(begc, endc, clm3%g%l%c%ces%pes_a)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
    call init_pft_wstate_type(begc, endc, clm3%g%l%c%cws%pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%ccs%pcs_a)
    
    if ( use_c13 ) then       
      call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc13s)
      call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc13s%pcs_a)
#ifdef CROP
      call fatal(__FILE__,__LINE__, &
       trim(subname)//" ERROR:: CROP and C13 can NOT be on at the same time")
#endif
    end if

    if ( use_c14 ) then
      call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc14s)
      call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc14s%pcs_a)
#ifdef CROP
      call fatal(__FILE__,__LINE__,&
       trim(subname)//" ERROR:: CROP and C14 can NOT be on at the same time")
#endif
    end if

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
    call init_pft_nstate_type(begc, endc, clm3%g%l%c%cns%pns_a)

    ! pft energy flux variables at pft level and averaged to column

    call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
    call init_pft_eflux_type(begc, endc, clm3%g%l%c%cef%pef_a)

    ! pft momentum flux variables at pft level and averaged to the column

    call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
    call init_pft_mflux_type(begc, endc, clm3%g%l%c%cmf%pmf_a)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
    call init_pft_wflux_type(begc, endc, clm3%g%l%c%cwf%pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)
    
    if ( use_c13 ) then       
      call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc13f)
      call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc13f%pcf_a)
    end if
    
    if ( use_c14 ) then
      call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc14f)
      call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc14f%pcf_a)
    end if

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
    call init_pft_nflux_type(begc, endc, clm3%g%l%c%cnf%pnf_a)

    ! pft VOC flux variables at pft level

    call init_pft_vflux_type(begp, endp, clm3%g%l%c%p%pvf)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(begg, endg, clm3%g%gve)

    ! pft dust flux variables at pft level and averaged to column

    call init_pft_dflux_type(begp, endp, clm3%g%l%c%p%pdf)
    call init_pft_dflux_type(begc, endc, clm3%g%l%c%cdf%pdf_a)

    ! pft dry dep velocity variables at pft level and averaged to column

    call init_pft_depvd_type(begp, endp, clm3%g%l%c%p%pdd)

    ! column physical state variables at column level and averaged to
    ! the landunit

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)
    call init_column_pstate_type(begl, endl, clm3%g%l%lps%cps_a)

    ! column energy state variables at column level and averaged to
    ! the gridcell

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
    call init_column_estate_type(begg, endg, clm3%g%ges%ces_a)

    ! column water state variables at column level and averaged to
    ! the gridcell

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
    call init_column_wstate_type(begg, endg, clm3%g%gws%cws_a)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
    
    if ( use_c13 ) then       
      call init_column_cstate_type(begc, endc, clm3%g%l%c%cc13s)
    end if

    if ( use_c14 ) then       
      call init_column_cstate_type(begc, endc, clm3%g%l%c%cc14s)
    end if
    ! column nitrogen state variables at column level

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)

    ! column energy flux variables at column level and averaged to
    ! the landunit and gridcell

    call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
    call init_column_eflux_type(begl, endl, clm3%g%l%lef%cef_a)
    call init_column_eflux_type(begg, endg, clm3%g%gef%cef_a)

    ! column water flux variables at column level and averaged to
    ! gridcell

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    call init_column_wflux_type(begg, endg, clm3%g%gwf%cwf_a)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)
    
    if ( use_c13 ) then       
      call init_column_cflux_type(begc, endc, clm3%g%l%c%cc13f)
    end if
    
    if ( use_c14 ) then       
      call init_column_cflux_type(begc, endc, clm3%g%l%c%cc14f)
    end if

#if (defined LCH4)
    ! column CH4 flux variables at column level
    call init_column_ch4_type(begc, endc, clm3%g%l%c%cch4)
#endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, clm3%g%l%lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(begl, endl, clm3%g%l%lef)

#if (defined CNDV)
    ! gridcell DGVM variables

    call init_gridcell_dgvstate_type(begg, endg, clm3%g%gdgvs)
#endif

    ! gridcell physical state variables


    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, clm3%g%gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(begg, endg, clm3%g%gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(begg, endg, clm3%g%gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(begg, endg, clm3%g%ges)

#if (defined LCH4)
    ! gridcell: ch4 variables

    call init_gridcell_ch4_type(begg, endg, clm3%g%gch4)
#endif

  contains
  !
  ! Initialize components of pft_type structure
  !
  subroutine init_pft_type (ibeg,iend,p)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(pft_type) , intent(inout) :: p

    allocate(p%gridcell(ibeg:iend),p%wtgcell(ibeg:iend))
    allocate(p%landunit(ibeg:iend),p%wtlunit(ibeg:iend))
    allocate(p%column  (ibeg:iend),p%wtcol  (ibeg:iend))

    allocate(p%itype(ibeg:iend))
    allocate(p%mxy(ibeg:iend))
    allocate(p%active(ibeg:iend))
  end subroutine init_pft_type
  !
  ! Initialize components of column_type structure
  !
  subroutine init_column_type (ibeg,iend,c)
! !ARGUMENTS:
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(column_type) , intent(inout) :: c

   allocate(c%gridcell(ibeg:iend),c%wtgcell(ibeg:iend))
   allocate(c%landunit(ibeg:iend),c%wtlunit(ibeg:iend))

   allocate(c%pfti(ibeg:iend),c%pftf(ibeg:iend),c%npfts(ibeg:iend))

   allocate(c%itype(ibeg:iend))
   allocate(c%active(ibeg:iend))
  end subroutine init_column_type
  !
  ! Initialize components of landunit_type structure
  !
  subroutine init_landunit_type (ibeg,iend,l)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(landunit_type) , intent(inout) :: l

    allocate(l%gridcell(ibeg:iend),l%wtgcell(ibeg:iend))

    allocate(l%coli(ibeg:iend),l%colf(ibeg:iend),l%ncolumns(ibeg:iend))
    allocate(l%pfti(ibeg:iend),l%pftf(ibeg:iend),l%npfts   (ibeg:iend))

    allocate(l%itype(ibeg:iend))
    allocate(l%ifspecial(ibeg:iend))
    allocate(l%lakpoi(ibeg:iend))
    allocate(l%urbpoi(ibeg:iend))
    allocate(l%glcmecpoi(ibeg:iend))
    allocate(l%udenstype(ibeg:iend))
    allocate(l%active(ibeg:iend))

    ! MV - these should be moved to landunit physical state -MV
    allocate(l%canyon_hwr(ibeg:iend))
    allocate(l%wtroad_perv(ibeg:iend))
    allocate(l%ht_roof(ibeg:iend))
    allocate(l%wtlunit_roof(ibeg:iend))
    allocate(l%wind_hgt_canyon(ibeg:iend))
    allocate(l%z_0_town(ibeg:iend))
    allocate(l%z_d_town(ibeg:iend))

    l%canyon_hwr(ibeg:iend)  = nan
    l%wtroad_perv(ibeg:iend) = nan
    l%ht_roof(ibeg:iend) = nan
    l%wtlunit_roof(ibeg:iend) = nan
    l%wind_hgt_canyon(ibeg:iend) = nan
    l%z_0_town(ibeg:iend) = nan
    l%z_d_town(ibeg:iend) = nan

    l%glcmecpoi(ibeg:iend) = .false.
  end subroutine init_landunit_type
  !
  ! Initialize components of gridcell_type structure
  !
  subroutine init_gridcell_type (ibeg,iend,g)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(gridcell_type) , intent(inout) :: g

    allocate(g%luni(ibeg:iend),g%lunf(ibeg:iend),g%nlandunits(ibeg:iend))
    allocate(g%coli(ibeg:iend),g%colf(ibeg:iend),g%ncolumns  (ibeg:iend))
    allocate(g%pfti(ibeg:iend),g%pftf(ibeg:iend),g%npfts     (ibeg:iend))

    allocate(g%gindex(ibeg:iend))
    allocate(g%area(ibeg:iend))
    allocate(g%lat(ibeg:iend))
    allocate(g%lon(ibeg:iend))
    allocate(g%latdeg(ibeg:iend))
    allocate(g%londeg(ibeg:iend))
    allocate(g%gindex_a(ibeg:iend))
    allocate(g%lat_a(ibeg:iend))
    allocate(g%lon_a(ibeg:iend))
    allocate(g%latdeg_a(ibeg:iend))
    allocate(g%londeg_a(ibeg:iend))

    allocate(g%gris_mask(ibeg:iend))
    allocate(g%gris_area(ibeg:iend))
    allocate(g%aais_mask(ibeg:iend))
    allocate(g%aais_area(ibeg:iend))
    allocate(g%tws(ibeg:iend))
    g%gris_mask(ibeg:iend) = nan
    g%gris_area(ibeg:iend) = nan
    g%aais_mask(ibeg:iend) = nan
    g%aais_area(ibeg:iend) = nan
    g%tws(ibeg:iend) = nan
  end subroutine init_gridcell_type
  !
  ! Initialize energy balance variables
  !
  subroutine init_energy_balance_type(ibeg,iend,ebal)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(energy_balance_type) , intent(inout) :: ebal

    allocate(ebal%errsoi(ibeg:iend))
    allocate(ebal%errseb(ibeg:iend))
    allocate(ebal%errsol(ibeg:iend))
    allocate(ebal%errlon(ibeg:iend))

    ebal%errsoi(ibeg:iend) = nan
    ebal%errseb(ibeg:iend) = nan
    ebal%errsol(ibeg:iend) = nan
    ebal%errlon(ibeg:iend) = nan
  end subroutine init_energy_balance_type
  !
  ! Initialize water balance variables
  !
  subroutine init_water_balance_type(ibeg,iend,wbal)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(water_balance_type) , intent(inout) :: wbal

    allocate(wbal%begwb(ibeg:iend))
    allocate(wbal%endwb(ibeg:iend))
    allocate(wbal%errh2o(ibeg:iend))

    wbal%begwb(ibeg:iend) = nan
    wbal%endwb(ibeg:iend) = nan
    wbal%errh2o(ibeg:iend) = nan
  end subroutine init_water_balance_type
  !
  ! Initialize carbon balance variables
  !
  subroutine init_carbon_balance_type(ibeg,iend,cbal)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(carbon_balance_type) , intent(inout) :: cbal

    allocate(cbal%begcb(ibeg:iend))
    allocate(cbal%endcb(ibeg:iend))
    allocate(cbal%errcb(ibeg:iend))

    cbal%begcb(ibeg:iend) = nan
    cbal%endcb(ibeg:iend) = nan
    cbal%errcb(ibeg:iend) = nan
  end subroutine init_carbon_balance_type
  !
  ! Initialize nitrogen balance variables
  !
  subroutine init_nitrogen_balance_type(ibeg,iend,nbal)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type(nitrogen_balance_type) , intent(inout) :: nbal

    allocate(nbal%begnb(ibeg:iend))
    allocate(nbal%endnb(ibeg:iend))
    allocate(nbal%errnb(ibeg:iend))

    nbal%begnb(ibeg:iend) = nan
    nbal%endnb(ibeg:iend) = nan
    nbal%errnb(ibeg:iend) = nan
  end subroutine init_nitrogen_balance_type
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_ecophys_constants()
    implicit none

    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%tree(0:numpft))
    allocate(pftcon%smpso(0:numpft)) 
    allocate(pftcon%smpsc(0:numpft)) 
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%graincn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%dwood(0:numpft))
    allocate(pftcon%rootprof_beta(0:numpft))
    allocate(pftcon%fertnitro(0:numpft))
    allocate(pftcon%fleafcn(0:numpft))
    allocate(pftcon%ffrootcn(0:numpft))
    allocate(pftcon%fstemcn(0:numpft))

    pftcon%noveg(:) = bigint
    pftcon%tree(:) = bigint
    pftcon%smpso(:) = nan
    pftcon%smpsc(:) = nan
    pftcon%fnitr(:) = nan
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%graincn(:) = nan
    pftcon%froot_leaf(:) = nan
    pftcon%stem_leaf(:) = nan
    pftcon%croot_stem(:) = nan
    pftcon%flivewd(:) = nan
    pftcon%fcur(:) = nan
    pftcon%lf_flab(:) = nan
    pftcon%lf_fcel(:) = nan
    pftcon%lf_flig(:) = nan
    pftcon%fr_flab(:) = nan
    pftcon%fr_fcel(:) = nan
    pftcon%fr_flig(:) = nan
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%dwood(:) = nan
    pftcon%rootprof_beta(:) = nan
    pftcon%fertnitro(:) = nan
    pftcon%fleafcn(:)   = nan
    pftcon%ffrootcn(:)  = nan
    pftcon%fstemcn(:)   = nan
  end subroutine init_pft_ecophys_constants
  !
  ! Initialize decomposition cascade state
  !
  subroutine init_decomp_cascade_constants()
    implicit none
    integer(ik4) :: nct , np
    nct = ndecomp_cascade_transitions
    np = ndecomp_pools
    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:nct))
    allocate(decomp_cascade_con%cascade_donor_pool(1:nct))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:nct))
    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:np))
    allocate(decomp_cascade_con%is_litter(0:np))
    allocate(decomp_cascade_con%is_soil(0:np))
    allocate(decomp_cascade_con%is_cwd(0:np))
    allocate(decomp_cascade_con%initial_cn_ratio(0:np))
    allocate(decomp_cascade_con%initial_stock(0:np))
    allocate(decomp_cascade_con%is_metabolic(0:np))
    allocate(decomp_cascade_con%is_cellulose(0:np))
    allocate(decomp_cascade_con%is_lignin(0:np))
    allocate(decomp_cascade_con%spinup_factor(0:np))
    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:nct) = ''
    decomp_cascade_con%cascade_donor_pool(1:nct) = 0
    decomp_cascade_con%cascade_receiver_pool(1:nct) = 0
    !-- properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:np) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:np) = ''
    decomp_cascade_con%decomp_pool_name_restart(0:np) = ''
    decomp_cascade_con%decomp_pool_name_long(0:np) = ''
    decomp_cascade_con%decomp_pool_name_short(0:np) = ''
    decomp_cascade_con%is_litter(0:np) = .false.
    decomp_cascade_con%is_soil(0:np) = .false.
    decomp_cascade_con%is_cwd(0:np) = .false.
    decomp_cascade_con%initial_cn_ratio(0:np) = nan
    decomp_cascade_con%initial_stock(0:np) = nan
    decomp_cascade_con%is_metabolic(0:np) = .false.
    decomp_cascade_con%is_cellulose(0:np) = .false.
    decomp_cascade_con%is_lignin(0:np) = .false.
    decomp_cascade_con%spinup_factor(0:np) = nan
  end subroutine init_decomp_cascade_constants

#if (defined CNDV)
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_DGVMecophys_constants()
    implicit none

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan
  end subroutine init_pft_DGVMecophys_constants

#endif
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_pstate_type(ibeg,iend,pps)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_pstate_type) , intent(inout) :: pps

    allocate(pps%prec10(ibeg:iend)) !F. Li and S. Levis
    allocate(pps%prec60(ibeg:iend)) !F. Li and S. Levis
    allocate(pps%frac_veg_nosno(ibeg:iend))
    allocate(pps%frac_veg_nosno_alb(ibeg:iend))
    allocate(pps%emv(ibeg:iend))
    allocate(pps%z0mv(ibeg:iend))
    allocate(pps%z0hv(ibeg:iend))
    allocate(pps%z0qv(ibeg:iend))
    allocate(pps%rootfr(ibeg:iend,1:nlevgrnd))
    allocate(pps%rootr(ibeg:iend,1:nlevgrnd))
    allocate(pps%rresis(ibeg:iend,1:nlevgrnd))
    allocate(pps%dewmx(ibeg:iend))
    allocate(pps%rssun(ibeg:iend))
    allocate(pps%rssha(ibeg:iend))
    allocate(pps%rhal(ibeg:iend))
    allocate(pps%vpdal(ibeg:iend))
    allocate(pps%rssun_z(ibeg:iend,1:nlevcan))
    allocate(pps%rssha_z(ibeg:iend,1:nlevcan))
    allocate(pps%laisun(ibeg:iend))
    allocate(pps%laisha(ibeg:iend))
    allocate(pps%laisun_z(ibeg:iend,1:nlevcan))
    allocate(pps%laisha_z(ibeg:iend,1:nlevcan))
    allocate(pps%btran(ibeg:iend))
    allocate(pps%btran2(ibeg:iend))   ! F. Li and S. Levis
    allocate(pps%fsun(ibeg:iend))
    allocate(pps%tlai(ibeg:iend))
    allocate(pps%tsai(ibeg:iend))
    allocate(pps%elai(ibeg:iend))
    allocate(pps%esai(ibeg:iend))
    allocate(pps%fwet(ibeg:iend))
    allocate(pps%fdry(ibeg:iend))
    allocate(pps%dt_veg(ibeg:iend))
    allocate(pps%htop(ibeg:iend))
    allocate(pps%hbot(ibeg:iend))
    allocate(pps%z0m(ibeg:iend))
    allocate(pps%displa(ibeg:iend))
    allocate(pps%albd(ibeg:iend,1:numrad))
    allocate(pps%albi(ibeg:iend,1:numrad))
    allocate(pps%fabd(ibeg:iend,1:numrad))
    allocate(pps%fabd_sun(ibeg:iend,1:numrad))
    allocate(pps%fabd_sha(ibeg:iend,1:numrad))
    allocate(pps%fabi(ibeg:iend,1:numrad))
    allocate(pps%fabi_sun(ibeg:iend,1:numrad))
    allocate(pps%fabi_sha(ibeg:iend,1:numrad))
    allocate(pps%ftdd(ibeg:iend,1:numrad))
    allocate(pps%ftid(ibeg:iend,1:numrad))
    allocate(pps%ftii(ibeg:iend,1:numrad))
    allocate(pps%vcmaxcintsun(ibeg:iend))
    allocate(pps%vcmaxcintsha(ibeg:iend))
    allocate(pps%ncan(ibeg:iend))
    allocate(pps%nrad(ibeg:iend))
    allocate(pps%fabd_sun_z(ibeg:iend,1:nlevcan))
    allocate(pps%fabd_sha_z(ibeg:iend,1:nlevcan))
    allocate(pps%fabi_sun_z(ibeg:iend,1:nlevcan))
    allocate(pps%fabi_sha_z(ibeg:iend,1:nlevcan))
    allocate(pps%fsun_z(ibeg:iend,1:nlevcan))
    allocate(pps%tlai_z(ibeg:iend,1:nlevcan))
    allocate(pps%tsai_z(ibeg:iend,1:nlevcan))
    allocate(pps%u10(ibeg:iend))
    allocate(pps%u10_clm(ibeg:iend))
    allocate(pps%va(ibeg:iend))
    allocate(pps%fv(ibeg:iend))
    allocate(pps%ram1(ibeg:iend))
    allocate(pps%burndate(ibeg:iend))   ! F. Li and S. Levis
    allocate(pps%ram1_lake(ibeg:iend))
    allocate(pps%rh_leaf(ibeg:iend))
    allocate(pps%rhaf(ibeg:iend))
    if ( crop_prog ) then
      allocate(pps%hdidx(ibeg:iend))
      allocate(pps%cumvd(ibeg:iend))
      allocate(pps%htmx(ibeg:iend))
      allocate(pps%vf(ibeg:iend))
      allocate(pps%gddmaturity(ibeg:iend))
      allocate(pps%gdd0(ibeg:iend))
      allocate(pps%gdd8(ibeg:iend))
      allocate(pps%gdd10(ibeg:iend))
      allocate(pps%gdd020(ibeg:iend))
      allocate(pps%gdd820(ibeg:iend))
      allocate(pps%gdd1020(ibeg:iend))
      allocate(pps%gddplant(ibeg:iend))
      allocate(pps%gddtsoi(ibeg:iend))
      allocate(pps%huileaf(ibeg:iend))
      allocate(pps%huigrain(ibeg:iend))
      allocate(pps%aleafi(ibeg:iend))
      allocate(pps%astemi(ibeg:iend))
      allocate(pps%aleaf(ibeg:iend))
      allocate(pps%astem(ibeg:iend))
      allocate(pps%croplive(ibeg:iend))
      allocate(pps%cropplant(ibeg:iend)) !,numpft)) ! make 2-D if using
      allocate(pps%harvdate(ibeg:iend))  !,numpft)) ! crop rotation
      allocate(pps%idop(ibeg:iend))
      allocate(pps%peaklai(ibeg:iend))
    end if
    allocate(pps%vds(ibeg:iend))
    allocate(pps%forc_hgt_u_pft(ibeg:iend))
    allocate(pps%forc_hgt_t_pft(ibeg:iend))
    allocate(pps%forc_hgt_q_pft(ibeg:iend))
    allocate(pps%lfpftd(ibeg:iend))      !F. Li and S. Levis

    ! 4/14/05: PET
    ! Adding isotope code
    
    if ( use_c13 ) then       
      allocate(pps%alphapsnsun(ibeg:iend))
      allocate(pps%alphapsnsha(ibeg:iend))
    endif

    allocate(pps%sandfrac(ibeg:iend))
    allocate(pps%clayfrac(ibeg:iend))
    pps%sandfrac(ibeg:iend) = nan
    pps%clayfrac(ibeg:iend) = nan
    allocate(pps%mlaidiff(ibeg:iend))
    allocate(pps%rb1(ibeg:iend))
    allocate(pps%annlai(12,ibeg:iend))
    pps%mlaidiff(ibeg:iend) = nan
    pps%rb1(ibeg:iend) = nan
    pps%annlai(:,:) = nan
    
#if (defined LCH4)
    ! CH4 code
    allocate(pps%grnd_ch4_cond(ibeg:iend))
    allocate(pps%canopy_cond(ibeg:iend))
#endif
   ! and vertical profiles for calculating fluxes
    allocate(pps%leaf_prof(ibeg:iend,1:nlevdecomp_full))
    allocate(pps%froot_prof(ibeg:iend,1:nlevdecomp_full))
    allocate(pps%croot_prof(ibeg:iend,1:nlevdecomp_full))
    allocate(pps%stem_prof(ibeg:iend,1:nlevdecomp_full))
    pps%prec10(ibeg:iend) = nan   ! F. Li and S. Levis
    pps%prec60(ibeg:iend) = nan   ! F. Li and S. Levis
    pps%frac_veg_nosno(ibeg:iend) = bigint
    pps%frac_veg_nosno_alb(ibeg:iend) = 0
    pps%emv(ibeg:iend) = nan
    pps%z0mv(ibeg:iend) = nan
    pps%z0hv(ibeg:iend) = nan
    pps%z0qv(ibeg:iend) = nan
    pps%rootfr(ibeg:iend,:nlevgrnd) = spval
    pps%rootr (ibeg:iend,:nlevgrnd) = spval
    pps%rresis(ibeg:iend,:nlevgrnd) = spval
    pps%dewmx(ibeg:iend) = nan
    pps%rssun(ibeg:iend) = nan
    pps%rhal(ibeg:iend) = nan
    pps%vpdal(ibeg:iend) = nan    
    pps%rssha(ibeg:iend) = nan
    pps%rssun_z(ibeg:iend,:nlevcan) = nan
    pps%rssha_z(ibeg:iend,:nlevcan) = nan
    pps%laisun(ibeg:iend) = nan
    pps%laisha(ibeg:iend) = nan
    pps%laisun_z(ibeg:iend,:nlevcan) = nan
    pps%laisha_z(ibeg:iend,:nlevcan) = nan
    pps%btran(ibeg:iend) = spval
    pps%btran2(ibeg:iend) = spval       !F. Li and S. Levis
    pps%fsun(ibeg:iend) = spval
    pps%fsun_z(ibeg:iend,:nlevcan) = 0.D0
    pps%tlai(ibeg:iend) = 0.D0
    pps%tsai(ibeg:iend) = 0.D0
    pps%elai(ibeg:iend) = 0.D0
    pps%tlai_z(ibeg:iend,:nlevcan) = 0.D0
    pps%tsai_z(ibeg:iend,:nlevcan) = 0.D0
    pps%esai(ibeg:iend) = 0.D0
    pps%ncan(ibeg:iend) = 0
    pps%nrad(ibeg:iend) = 0
    pps%fwet(ibeg:iend) = nan
    pps%fdry(ibeg:iend) = nan
    pps%dt_veg(ibeg:iend) = nan
    pps%htop(ibeg:iend) = 0.D0
    pps%hbot(ibeg:iend) = 0.D0
    pps%z0m(ibeg:iend) = nan
    pps%displa(ibeg:iend) = nan
    pps%albd(ibeg:iend,:numrad) = nan
    pps%albi(ibeg:iend,:numrad) = nan
    pps%fabd(ibeg:iend,:numrad) = nan
    pps%fabd_sun(ibeg:iend,:numrad) = nan
    pps%fabd_sha(ibeg:iend,:numrad) = nan
    pps%fabi(ibeg:iend,:numrad) = nan
    pps%fabi_sun(ibeg:iend,:numrad) = nan
    pps%fabi_sha(ibeg:iend,:numrad) = nan
    pps%ftdd(ibeg:iend,:numrad) = nan
    pps%ftid(ibeg:iend,:numrad) = nan
    pps%ftii(ibeg:iend,:numrad) = nan
    pps%vcmaxcintsun(ibeg:iend) = nan
    pps%vcmaxcintsha(ibeg:iend) = nan
    pps%fabd_sun_z(ibeg:iend,:nlevcan) = 0.D0
    pps%fabd_sha_z(ibeg:iend,:nlevcan) = 0.D0
    pps%fabi_sun_z(ibeg:iend,:nlevcan) = 0.D0
    pps%fabi_sha_z(ibeg:iend,:nlevcan) = 0.D0
    pps%u10(ibeg:iend) = nan
    pps%u10_clm(ibeg:iend) = nan
    pps%va(ibeg:iend) = nan
    pps%fv(ibeg:iend) = nan
    pps%ram1(ibeg:iend) = nan
    pps%burndate(ibeg:iend)    = ispval   ! F. Li and S. Levis
    pps%ram1_lake(ibeg:iend) = nan
    pps%rh_leaf(ibeg:iend) = spval
    pps%rhaf(ibeg:iend)    = spval
    if ( crop_prog ) then
      pps%hdidx(ibeg:iend)       = nan
      pps%cumvd(ibeg:iend)       = nan
      pps%htmx(ibeg:iend)        = 0.0D0
      pps%vf(ibeg:iend)          = 0.0D0
      pps%gddmaturity(ibeg:iend) = spval
      pps%gdd0(ibeg:iend)        = spval
      pps%gdd8(ibeg:iend)        = spval
      pps%gdd10(ibeg:iend)       = spval
      pps%gdd020(ibeg:iend)      = spval
      pps%gdd820(ibeg:iend)      = spval
      pps%gdd1020(ibeg:iend)     = spval
      pps%gddplant(ibeg:iend)    = spval
      pps%gddtsoi(ibeg:iend)     = spval
      pps%huileaf(ibeg:iend)     = nan
      pps%huigrain(ibeg:iend)    = nan
      pps%aleafi(ibeg:iend)      = nan
      pps%astemi(ibeg:iend)      = nan
      pps%aleaf(ibeg:iend)       = nan
      pps%astem(ibeg:iend)       = nan
      pps%croplive(ibeg:iend)    = .false.
      pps%cropplant(ibeg:iend)   = .false.
      pps%harvdate(ibeg:iend)    = bigint
      pps%idop(ibeg:iend)        = bigint
      pps%peaklai(ibeg:iend)     = 0
    end if
    pps%vds(ibeg:iend) = nan
    pps%forc_hgt_u_pft(ibeg:iend) = nan
    pps%forc_hgt_t_pft(ibeg:iend) = nan
    pps%forc_hgt_q_pft(ibeg:iend) = nan
    ! 4/14/05: PET
    ! Adding isotope code    ! EBK Check this!
    !!!pps%cisun(ibeg:iend) = spval
    !!!pps%cisha(ibeg:iend) = spval
    
    if ( use_c13 ) then       
      pps%alphapsnsun(ibeg:iend) = spval
      pps%alphapsnsha(ibeg:iend) = spval
    endif

#if defined (LCH4)
    ! CH4 code
    pps%grnd_ch4_cond(ibeg:iend) = nan
    pps%canopy_cond(ibeg:iend) = nan
#endif
   ! and vertical profiles for calculating fluxes
    pps%leaf_prof(ibeg:iend,1:nlevdecomp_full) = spval
    pps%froot_prof(ibeg:iend,1:nlevdecomp_full) = spval
    pps%croot_prof(ibeg:iend,1:nlevdecomp_full) = spval
    pps%stem_prof(ibeg:iend,1:nlevdecomp_full) = spval

  end subroutine init_pft_pstate_type
  !
  ! Initialize pft ecophysiological variables
  !
  subroutine init_pft_epv_type(ibeg,iend,pepv)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_epv_type) , intent(inout) :: pepv

    allocate(pepv%dormant_flag(ibeg:iend))
    allocate(pepv%days_active(ibeg:iend))
    allocate(pepv%onset_flag(ibeg:iend))
    allocate(pepv%onset_counter(ibeg:iend))
    allocate(pepv%onset_gddflag(ibeg:iend))
    allocate(pepv%onset_fdd(ibeg:iend))
    allocate(pepv%onset_gdd(ibeg:iend))
    allocate(pepv%onset_swi(ibeg:iend))
    allocate(pepv%offset_flag(ibeg:iend))
    allocate(pepv%offset_counter(ibeg:iend))
    allocate(pepv%offset_fdd(ibeg:iend))
    allocate(pepv%offset_swi(ibeg:iend))
    allocate(pepv%fert_counter(ibeg:iend))
    allocate(pepv%grain_flag(ibeg:iend))
    allocate(pepv%lgsf(ibeg:iend))
    allocate(pepv%bglfr(ibeg:iend))
    allocate(pepv%bgtr(ibeg:iend))
    allocate(pepv%dayl(ibeg:iend))
    allocate(pepv%prev_dayl(ibeg:iend))
    allocate(pepv%annavg_t2m(ibeg:iend))
    allocate(pepv%tempavg_t2m(ibeg:iend))
    allocate(pepv%gpp(ibeg:iend))
    allocate(pepv%availc(ibeg:iend))
    allocate(pepv%xsmrpool_recover(ibeg:iend))
    if ( use_c13 ) then
      allocate(pepv%xsmrpool_c13ratio(ibeg:iend))
    endif
    allocate(pepv%alloc_pnow(ibeg:iend))
    allocate(pepv%c_allometry(ibeg:iend))
    allocate(pepv%n_allometry(ibeg:iend))
    allocate(pepv%plant_ndemand(ibeg:iend))
    allocate(pepv%tempsum_potential_gpp(ibeg:iend))
    allocate(pepv%annsum_potential_gpp(ibeg:iend))
    allocate(pepv%tempmax_retransn(ibeg:iend))
    allocate(pepv%annmax_retransn(ibeg:iend))
    allocate(pepv%avail_retransn(ibeg:iend))
    allocate(pepv%plant_nalloc(ibeg:iend))
    allocate(pepv%plant_calloc(ibeg:iend))
    allocate(pepv%excess_cflux(ibeg:iend))
    allocate(pepv%downreg(ibeg:iend))
    allocate(pepv%prev_leafc_to_litter(ibeg:iend))
    allocate(pepv%prev_frootc_to_litter(ibeg:iend))
    allocate(pepv%tempsum_npp(ibeg:iend))
    allocate(pepv%annsum_npp(ibeg:iend))
#if (defined CNDV)
    allocate(pepv%tempsum_litfall(ibeg:iend))
    allocate(pepv%annsum_litfall(ibeg:iend))
#endif
    if ( use_c13 ) then
      allocate(pepv%rc13_canair(ibeg:iend))
      allocate(pepv%rc13_psnsun(ibeg:iend))
      allocate(pepv%rc13_psnsha(ibeg:iend))
    endif
    
    if ( use_c14 ) then
      allocate(pepv%rc14_atm(ibeg:iend))
      ! allocate(pepv%rc14_canair(ibeg:iend))
      ! allocate(pepv%rc14_psnsun(ibeg:iend))
      ! allocate(pepv%rc14_psnsha(ibeg:iend))
    endif

    pepv%dormant_flag(ibeg:iend) = nan
    pepv%days_active(ibeg:iend) = nan
    pepv%onset_flag(ibeg:iend) = nan
    pepv%onset_counter(ibeg:iend) = nan
    pepv%onset_gddflag(ibeg:iend) = nan
    pepv%onset_fdd(ibeg:iend) = nan
    pepv%onset_gdd(ibeg:iend) = nan
    pepv%onset_swi(ibeg:iend) = nan
    pepv%offset_flag(ibeg:iend) = nan
    pepv%offset_counter(ibeg:iend) = nan
    pepv%offset_fdd(ibeg:iend) = nan
    pepv%offset_swi(ibeg:iend) = nan
    pepv%fert_counter(ibeg:iend) = nan
    pepv%grain_flag(ibeg:iend) = nan
    pepv%lgsf(ibeg:iend) = nan
    pepv%bglfr(ibeg:iend) = nan
    pepv%bgtr(ibeg:iend) = nan
    pepv%dayl(ibeg:iend) = nan
    pepv%prev_dayl(ibeg:iend) = nan
    pepv%annavg_t2m(ibeg:iend) = nan
    pepv%tempavg_t2m(ibeg:iend) = nan
    pepv%gpp(ibeg:iend) = nan
    pepv%availc(ibeg:iend) = nan
    pepv%xsmrpool_recover(ibeg:iend) = nan
    if ( use_c13 ) then
      pepv%xsmrpool_c13ratio(ibeg:iend) = nan
    endif
    pepv%alloc_pnow(ibeg:iend) = nan
    pepv%c_allometry(ibeg:iend) = nan
    pepv%n_allometry(ibeg:iend) = nan
    pepv%plant_ndemand(ibeg:iend) = nan
    pepv%tempsum_potential_gpp(ibeg:iend) = nan
    pepv%annsum_potential_gpp(ibeg:iend) = nan
    pepv%tempmax_retransn(ibeg:iend) = nan
    pepv%annmax_retransn(ibeg:iend) = nan
    pepv%avail_retransn(ibeg:iend) = nan
    pepv%plant_nalloc(ibeg:iend) = nan
    pepv%plant_calloc(ibeg:iend) = nan
    pepv%excess_cflux(ibeg:iend) = nan
    pepv%downreg(ibeg:iend) = nan
    pepv%prev_leafc_to_litter(ibeg:iend) = nan
    pepv%prev_frootc_to_litter(ibeg:iend) = nan
    pepv%tempsum_npp(ibeg:iend) = nan
    pepv%annsum_npp(ibeg:iend) = nan
#if (defined CNDV)
    pepv%tempsum_litfall(ibeg:iend) = nan
    pepv%annsum_litfall(ibeg:iend) = nan
#endif
    if ( use_c13 ) then
      pepv%rc13_canair(ibeg:iend) = spval
      pepv%rc13_psnsun(ibeg:iend) = spval
      pepv%rc13_psnsha(ibeg:iend) = spval
    endif

    if ( use_c14 ) then
      pepv%rc14_atm(ibeg:iend) = nan
      ! pepv%rc14_canair(ibeg:iend) = nan
      ! pepv%rc14_psnsun(ibeg:iend) = nan
      ! pepv%rc14_psnsha(ibeg:iend) = nan
    endif
  end subroutine init_pft_epv_type

#if (defined CNDV)
  !
  ! Initialize pft DGVM state variables
  !
  subroutine init_pft_pdgvstate_type(ibeg,iend,pdgvs)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_dgvstate_type) , intent(inout) :: pdgvs

    allocate(pdgvs%agddtw(ibeg:iend))
    allocate(pdgvs%agdd(ibeg:iend))
    allocate(pdgvs%t_mo(ibeg:iend))
    allocate(pdgvs%t_mo_min(ibeg:iend))
    allocate(pdgvs%prec365(ibeg:iend))
    allocate(pdgvs%present(ibeg:iend))
    allocate(pdgvs%pftmayexist(ibeg:iend))
    allocate(pdgvs%nind(ibeg:iend))
    allocate(pdgvs%lm_ind(ibeg:iend))
    allocate(pdgvs%lai_ind(ibeg:iend))
    allocate(pdgvs%fpcinc(ibeg:iend))
    allocate(pdgvs%fpcgrid(ibeg:iend))
    allocate(pdgvs%fpcgridold(ibeg:iend))
    allocate(pdgvs%crownarea(ibeg:iend))
    allocate(pdgvs%greffic(ibeg:iend))
    allocate(pdgvs%heatstress(ibeg:iend))

    pdgvs%agddtw(ibeg:iend)           = nan
    pdgvs%agdd(ibeg:iend)             = nan
    pdgvs%t_mo(ibeg:iend)             = nan
    pdgvs%t_mo_min(ibeg:iend)         = nan
    pdgvs%prec365(ibeg:iend)          = nan
    pdgvs%present(ibeg:iend)          = .false.
    pdgvs%pftmayexist(ibeg:iend)      = .true.
    pdgvs%nind(ibeg:iend)             = nan
    pdgvs%lm_ind(ibeg:iend)           = nan
    pdgvs%lai_ind(ibeg:iend)          = nan
    pdgvs%fpcinc(ibeg:iend)           = nan
    pdgvs%fpcgrid(ibeg:iend)          = nan
    pdgvs%fpcgridold(ibeg:iend)       = nan
    pdgvs%crownarea(ibeg:iend)        = nan
    pdgvs%greffic(ibeg:iend)          = nan
    pdgvs%heatstress(ibeg:iend)       = nan
  end subroutine init_pft_pdgvstate_type
#endif
  !
  ! Initialize pft VOC variables
  !
  subroutine init_pft_vstate_type(ibeg,iend,pvs)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_vstate_type) , intent(inout) :: pvs

    allocate(pvs%t_veg24 (ibeg:iend))
    allocate(pvs%t_veg240(ibeg:iend))
    allocate(pvs%fsd24   (ibeg:iend))
    allocate(pvs%fsd240  (ibeg:iend))
    allocate(pvs%fsi24   (ibeg:iend))
    allocate(pvs%fsi240  (ibeg:iend))
    allocate(pvs%fsun24  (ibeg:iend))
    allocate(pvs%fsun240 (ibeg:iend))
    allocate(pvs%elai_p  (ibeg:iend))

    pvs%t_veg24 (ibeg:iend)   = spval
    pvs%t_veg240(ibeg:iend)   = spval
    pvs%fsd24   (ibeg:iend)   = spval
    pvs%fsd240  (ibeg:iend)   = spval
    pvs%fsi24   (ibeg:iend)   = spval
    pvs%fsi240  (ibeg:iend)   = spval
    pvs%fsun24  (ibeg:iend)   = spval
    pvs%fsun240 (ibeg:iend)   = spval
    pvs%elai_p  (ibeg:iend)   = spval
  end subroutine init_pft_vstate_type
  !
  ! Initialize pft energy state
  !
  subroutine init_pft_psynstate_type(ibeg,iend,ppsyns)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_psynstate_type) , intent(inout) :: ppsyns

    allocate(ppsyns%c3flag(ibeg:iend))
    allocate(ppsyns%ac(ibeg:iend,1:nlevcan))
    allocate(ppsyns%aj(ibeg:iend,1:nlevcan))
    allocate(ppsyns%ap(ibeg:iend,1:nlevcan))
    allocate(ppsyns%ag(ibeg:iend,1:nlevcan))
    allocate(ppsyns%an(ibeg:iend,1:nlevcan))
    allocate(ppsyns%vcmax_z(ibeg:iend,1:nlevcan))
    allocate(ppsyns%cp(ibeg:iend))
    allocate(ppsyns%kc(ibeg:iend))
    allocate(ppsyns%ko(ibeg:iend))
    allocate(ppsyns%qe(ibeg:iend))
    allocate(ppsyns%tpu_z(ibeg:iend,1:nlevcan))
    allocate(ppsyns%kp_z(ibeg:iend,1:nlevcan))   
    allocate(ppsyns%theta_cj(ibeg:iend))
    allocate(ppsyns%bbb(ibeg:iend))
    allocate(ppsyns%mbb(ibeg:iend))
    allocate(ppsyns%gb_mol(ibeg:iend))
    allocate(ppsyns%gs_mol(ibeg:iend,1:nlevcan))

    ppsyns%c3flag = .false.
    ppsyns%ac(ibeg:iend,1:nlevcan) = nan
    ppsyns%aj(ibeg:iend,1:nlevcan) = nan
    ppsyns%ap(ibeg:iend,1:nlevcan) = nan
    ppsyns%ag(ibeg:iend,1:nlevcan) = nan
    ppsyns%an(ibeg:iend,1:nlevcan) = nan
    ppsyns%vcmax_z(ibeg:iend,1:nlevcan) = nan
    ppsyns%cp(ibeg:iend) = nan
    ppsyns%kc(ibeg:iend) = nan
    ppsyns%ko(ibeg:iend) = nan
    ppsyns%qe(ibeg:iend) = nan
    ppsyns%tpu_z(ibeg:iend,1:nlevcan) = nan
    ppsyns%kp_z(ibeg:iend,1:nlevcan) = nan
    ppsyns%theta_cj(ibeg:iend) = nan
    ppsyns%bbb(ibeg:iend) = nan
    ppsyns%mbb(ibeg:iend) = nan
    ppsyns%gb_mol(ibeg:iend) = nan
    ppsyns%gs_mol(ibeg:iend,1:nlevcan) = nan
  end subroutine init_pft_psynstate_type
  !
  ! Initialize pft energy state
  !
  subroutine init_pft_estate_type(ibeg,iend,pes)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_estate_type) , intent(inout) :: pes

    allocate(pes%t_ref2m(ibeg:iend))
    allocate(pes%t_ref2m_min(ibeg:iend))
    allocate(pes%t_ref2m_max(ibeg:iend))
    allocate(pes%t_ref2m_min_inst(ibeg:iend))
    allocate(pes%t_ref2m_max_inst(ibeg:iend))
    allocate(pes%q_ref2m(ibeg:iend))
    allocate(pes%t_ref2m_u(ibeg:iend))
    allocate(pes%t_ref2m_r(ibeg:iend))
    allocate(pes%t_ref2m_min_u(ibeg:iend))
    allocate(pes%t_ref2m_min_r(ibeg:iend))
    allocate(pes%t_ref2m_max_u(ibeg:iend))
    allocate(pes%t_ref2m_max_r(ibeg:iend))
    allocate(pes%t_ref2m_min_inst_u(ibeg:iend))
    allocate(pes%t_ref2m_min_inst_r(ibeg:iend))
    allocate(pes%t_ref2m_max_inst_u(ibeg:iend))
    allocate(pes%t_ref2m_max_inst_r(ibeg:iend))
    allocate(pes%t10(ibeg:iend))
    if ( crop_prog )then
      allocate(pes%a10tmin(ibeg:iend))
      allocate(pes%a5tmin(ibeg:iend))
    end if
    allocate(pes%rh_ref2m(ibeg:iend))
    allocate(pes%rh_ref2m_u(ibeg:iend))
    allocate(pes%rh_ref2m_r(ibeg:iend))
    allocate(pes%t_veg(ibeg:iend))
    allocate(pes%thm(ibeg:iend))

    pes%t_ref2m(ibeg:iend) = nan
    pes%t_ref2m_min(ibeg:iend) = nan
    pes%t_ref2m_max(ibeg:iend) = nan
    pes%t_ref2m_min_inst(ibeg:iend) = nan
    pes%t_ref2m_max_inst(ibeg:iend) = nan
    pes%q_ref2m(ibeg:iend) = nan
    pes%t_ref2m_u(ibeg:iend) = nan
    pes%t_ref2m_r(ibeg:iend) = nan
    pes%t_ref2m_min_u(ibeg:iend) = nan
    pes%t_ref2m_min_r(ibeg:iend) = nan
    pes%t_ref2m_max_u(ibeg:iend) = nan
    pes%t_ref2m_max_r(ibeg:iend) = nan
    pes%t_ref2m_min_inst_u(ibeg:iend) = nan
    pes%t_ref2m_min_inst_r(ibeg:iend) = nan
    pes%t_ref2m_max_inst_u(ibeg:iend) = nan
    pes%t_ref2m_max_inst_r(ibeg:iend) = nan
    pes%t10(ibeg:iend)                = spval
    if ( crop_prog )then
      pes%a10tmin(ibeg:iend)     = spval
      pes%a5tmin(ibeg:iend)      = spval
    end if
    pes%rh_ref2m(ibeg:iend) = nan
    pes%rh_ref2m_u(ibeg:iend) = nan
    pes%rh_ref2m_r(ibeg:iend) = nan
    pes%t_veg(ibeg:iend) = nan
    pes%thm(ibeg:iend) = nan
  end subroutine init_pft_estate_type
  !
  ! Initialize pft water state
  !
  subroutine init_pft_wstate_type(ibeg,iend,pws)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_wstate_type) , intent(inout) :: pws !pft water state

    allocate(pws%h2ocan(ibeg:iend))
    pws%h2ocan(ibeg:iend) = spval
  end subroutine init_pft_wstate_type
  !
  ! Initialize pft carbon state
  !
  subroutine init_pft_cstate_type(ibeg,iend,pcs)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_cstate_type) , intent(inout) :: pcs !pft carbon state

    allocate(pcs%leafc(ibeg:iend))
    allocate(pcs%leafc_storage(ibeg:iend))
    allocate(pcs%leafc_xfer(ibeg:iend))
    allocate(pcs%frootc(ibeg:iend))
    allocate(pcs%frootc_storage(ibeg:iend))
    allocate(pcs%frootc_xfer(ibeg:iend))
    allocate(pcs%livestemc(ibeg:iend))
    allocate(pcs%livestemc_storage(ibeg:iend))
    allocate(pcs%livestemc_xfer(ibeg:iend))
    allocate(pcs%deadstemc(ibeg:iend))
    allocate(pcs%deadstemc_storage(ibeg:iend))
    allocate(pcs%deadstemc_xfer(ibeg:iend))
    allocate(pcs%livecrootc(ibeg:iend))
    allocate(pcs%livecrootc_storage(ibeg:iend))
    allocate(pcs%livecrootc_xfer(ibeg:iend))
    allocate(pcs%deadcrootc(ibeg:iend))
    allocate(pcs%deadcrootc_storage(ibeg:iend))
    allocate(pcs%deadcrootc_xfer(ibeg:iend))
    allocate(pcs%gresp_storage(ibeg:iend))
    allocate(pcs%gresp_xfer(ibeg:iend))
    allocate(pcs%cpool(ibeg:iend))
    allocate(pcs%xsmrpool(ibeg:iend))
    allocate(pcs%pft_ctrunc(ibeg:iend))
    allocate(pcs%dispvegc(ibeg:iend))
    allocate(pcs%storvegc(ibeg:iend))
    allocate(pcs%totvegc(ibeg:iend))
    allocate(pcs%totpftc(ibeg:iend))
    allocate(pcs%leafcmax(ibeg:iend))
    if ( crop_prog )then
      allocate(pcs%grainc(ibeg:iend))
      allocate(pcs%grainc_storage(ibeg:iend))
      allocate(pcs%grainc_xfer(ibeg:iend))
    end if
#ifdef CN
    allocate(pcs%woodc(ibeg:iend))
#endif

    pcs%leafc(ibeg:iend) = nan
    pcs%leafc_storage(ibeg:iend) = nan
    pcs%leafc_xfer(ibeg:iend) = nan
    pcs%frootc(ibeg:iend) = nan
    pcs%frootc_storage(ibeg:iend) = nan
    pcs%frootc_xfer(ibeg:iend) = nan
    pcs%livestemc(ibeg:iend) = nan
    pcs%livestemc_storage(ibeg:iend) = nan
    pcs%livestemc_xfer(ibeg:iend) = nan
    pcs%deadstemc(ibeg:iend) = nan
    pcs%deadstemc_storage(ibeg:iend) = nan
    pcs%deadstemc_xfer(ibeg:iend) = nan
    pcs%livecrootc(ibeg:iend) = nan
    pcs%livecrootc_storage(ibeg:iend) = nan
    pcs%livecrootc_xfer(ibeg:iend) = nan
    pcs%deadcrootc(ibeg:iend) = nan
    pcs%deadcrootc_storage(ibeg:iend) = nan
    pcs%deadcrootc_xfer(ibeg:iend) = nan
    pcs%gresp_storage(ibeg:iend) = nan
    pcs%gresp_xfer(ibeg:iend) = nan
    pcs%cpool(ibeg:iend) = nan
    pcs%xsmrpool(ibeg:iend) = nan
    pcs%pft_ctrunc(ibeg:iend) = nan
    pcs%dispvegc(ibeg:iend) = nan
    pcs%storvegc(ibeg:iend) = nan
    pcs%totvegc(ibeg:iend) = nan
    pcs%totpftc(ibeg:iend) = nan
    pcs%leafcmax(ibeg:iend) = nan
    if ( crop_prog )then
      pcs%grainc(ibeg:iend)         = nan
      pcs%grainc_storage(ibeg:iend) = nan
      pcs%grainc_xfer(ibeg:iend)    = nan
    end if
#ifdef CN
    pcs%woodc(ibeg:iend) = nan
#endif
  end subroutine init_pft_cstate_type
  !
  ! Initialize pft nitrogen state
  !
  subroutine init_pft_nstate_type(ibeg,iend,pns)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_nstate_type) , intent(inout) :: pns !pft nitrogen state

    if ( crop_prog )then
      allocate(pns%grainn(ibeg:iend))
      allocate(pns%grainn_storage(ibeg:iend))
      allocate(pns%grainn_xfer(ibeg:iend))
    end if
    allocate(pns%leafn(ibeg:iend))
    allocate(pns%leafn_storage(ibeg:iend))
    allocate(pns%leafn_xfer(ibeg:iend))
    allocate(pns%frootn(ibeg:iend))
    allocate(pns%frootn_storage(ibeg:iend))
    allocate(pns%frootn_xfer(ibeg:iend))
    allocate(pns%livestemn(ibeg:iend))
    allocate(pns%livestemn_storage(ibeg:iend))
    allocate(pns%livestemn_xfer(ibeg:iend))
    allocate(pns%deadstemn(ibeg:iend))
    allocate(pns%deadstemn_storage(ibeg:iend))
    allocate(pns%deadstemn_xfer(ibeg:iend))
    allocate(pns%livecrootn(ibeg:iend))
    allocate(pns%livecrootn_storage(ibeg:iend))
    allocate(pns%livecrootn_xfer(ibeg:iend))
    allocate(pns%deadcrootn(ibeg:iend))
    allocate(pns%deadcrootn_storage(ibeg:iend))
    allocate(pns%deadcrootn_xfer(ibeg:iend))
    allocate(pns%retransn(ibeg:iend))
    allocate(pns%npool(ibeg:iend))
    allocate(pns%pft_ntrunc(ibeg:iend))
    allocate(pns%dispvegn(ibeg:iend))
    allocate(pns%storvegn(ibeg:iend))
    allocate(pns%totvegn(ibeg:iend))
    allocate(pns%totpftn(ibeg:iend))

    if ( crop_prog )then
      pns%grainn(ibeg:iend)         = nan
      pns%grainn_storage(ibeg:iend) = nan
      pns%grainn_xfer(ibeg:iend)    = nan
    end if
    pns%leafn(ibeg:iend) = nan
    pns%leafn_storage(ibeg:iend) = nan
    pns%leafn_xfer(ibeg:iend) = nan
    pns%frootn(ibeg:iend) = nan
    pns%frootn_storage(ibeg:iend) = nan
    pns%frootn_xfer(ibeg:iend) = nan
    pns%livestemn(ibeg:iend) = nan
    pns%livestemn_storage(ibeg:iend) = nan
    pns%livestemn_xfer(ibeg:iend) = nan
    pns%deadstemn(ibeg:iend) = nan
    pns%deadstemn_storage(ibeg:iend) = nan
    pns%deadstemn_xfer(ibeg:iend) = nan
    pns%livecrootn(ibeg:iend) = nan
    pns%livecrootn_storage(ibeg:iend) = nan
    pns%livecrootn_xfer(ibeg:iend) = nan
    pns%deadcrootn(ibeg:iend) = nan
    pns%deadcrootn_storage(ibeg:iend) = nan
    pns%deadcrootn_xfer(ibeg:iend) = nan
    pns%retransn(ibeg:iend) = nan
    pns%npool(ibeg:iend) = nan
    pns%pft_ntrunc(ibeg:iend) = nan
    pns%dispvegn(ibeg:iend) = nan
    pns%storvegn(ibeg:iend) = nan
    pns%totvegn(ibeg:iend) = nan
    pns%totpftn(ibeg:iend) = nan
  end subroutine init_pft_nstate_type
  !
  ! Initialize pft energy flux variables
  !
  subroutine init_pft_eflux_type(ibeg,iend,pef)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_eflux_type) , intent(inout) :: pef

    allocate(pef%sabg(ibeg:iend))
    allocate(pef%sabv(ibeg:iend))
    allocate(pef%fsa(ibeg:iend))
    allocate(pef%fsa_u(ibeg:iend))
    allocate(pef%fsa_r(ibeg:iend))
    allocate(pef%fsr(ibeg:iend))
    allocate(pef%parsun_z(ibeg:iend,1:nlevcan))
    allocate(pef%parsha_z(ibeg:iend,1:nlevcan))
    allocate(pef%dlrad(ibeg:iend))
    allocate(pef%ulrad(ibeg:iend))
    allocate(pef%eflx_lh_tot(ibeg:iend))
    allocate(pef%eflx_lh_tot_u(ibeg:iend))
    allocate(pef%eflx_lh_tot_r(ibeg:iend))
    allocate(pef%eflx_lh_grnd(ibeg:iend))
    allocate(pef%eflx_soil_grnd(ibeg:iend))
    allocate(pef%eflx_soil_grnd_u(ibeg:iend))
    allocate(pef%eflx_soil_grnd_r(ibeg:iend))
    allocate(pef%eflx_sh_tot(ibeg:iend))
    allocate(pef%eflx_sh_tot_u(ibeg:iend))
    allocate(pef%eflx_sh_tot_r(ibeg:iend))
    allocate(pef%eflx_sh_grnd(ibeg:iend))
    allocate(pef%eflx_sh_veg(ibeg:iend))
    allocate(pef%eflx_lh_vege(ibeg:iend))
    allocate(pef%eflx_lh_vegt(ibeg:iend))
    allocate(pef%eflx_wasteheat_pft(ibeg:iend))
    allocate(pef%eflx_heat_from_ac_pft(ibeg:iend))
    allocate(pef%eflx_traffic_pft(ibeg:iend))
    allocate(pef%eflx_anthro(ibeg:iend))
    allocate(pef%cgrnd(ibeg:iend))
    allocate(pef%cgrndl(ibeg:iend))
    allocate(pef%cgrnds(ibeg:iend))
    allocate(pef%eflx_gnet(ibeg:iend))
    allocate(pef%eflx_grnd_lake(ibeg:iend))
    allocate(pef%dgnetdT(ibeg:iend))
    allocate(pef%eflx_lwrad_out(ibeg:iend))
    allocate(pef%eflx_lwrad_net(ibeg:iend))
    allocate(pef%eflx_lwrad_net_u(ibeg:iend))
    allocate(pef%eflx_lwrad_net_r(ibeg:iend))
    allocate(pef%netrad(ibeg:iend))
    allocate(pef%fsds_vis_d(ibeg:iend))
    allocate(pef%fsds_nir_d(ibeg:iend))
    allocate(pef%fsds_vis_i(ibeg:iend))
    allocate(pef%fsds_nir_i(ibeg:iend))
    allocate(pef%fsr_vis_d(ibeg:iend))
    allocate(pef%fsr_nir_d(ibeg:iend))
    allocate(pef%fsr_vis_i(ibeg:iend))
    allocate(pef%fsr_nir_i(ibeg:iend))
    allocate(pef%fsds_vis_d_ln(ibeg:iend))
    allocate(pef%fsds_vis_i_ln(ibeg:iend))
    allocate(pef%parveg_ln(ibeg:iend))
    allocate(pef%fsds_nir_d_ln(ibeg:iend))
    allocate(pef%fsr_vis_d_ln(ibeg:iend))
    allocate(pef%fsr_nir_d_ln(ibeg:iend))
    allocate(pef%sabg_lyr(ibeg:iend,-nlevsno+1:1))
    allocate(pef%sabg_pen(ibeg:iend))
    allocate(pef%sfc_frc_aer(ibeg:iend))
    allocate(pef%sfc_frc_bc(ibeg:iend))
    allocate(pef%sfc_frc_oc(ibeg:iend))
    allocate(pef%sfc_frc_dst(ibeg:iend))
    allocate(pef%sfc_frc_aer_sno(ibeg:iend))
    allocate(pef%sfc_frc_bc_sno(ibeg:iend))
    allocate(pef%sfc_frc_oc_sno(ibeg:iend))
    allocate(pef%sfc_frc_dst_sno(ibeg:iend))
    allocate(pef%fsr_sno_vd(ibeg:iend))
    allocate(pef%fsr_sno_nd(ibeg:iend))
    allocate(pef%fsr_sno_vi(ibeg:iend))
    allocate(pef%fsr_sno_ni(ibeg:iend))
    allocate(pef%fsds_sno_vd(ibeg:iend))
    allocate(pef%fsds_sno_nd(ibeg:iend))
    allocate(pef%fsds_sno_vi(ibeg:iend))
    allocate(pef%fsds_sno_ni(ibeg:iend))

    pef%sabg(ibeg:iend) = nan
    pef%sabv(ibeg:iend) = nan
    pef%fsa(ibeg:iend) = nan
    pef%fsa_u(ibeg:iend) = nan
    pef%fsa_r(ibeg:iend) = nan
    pef%fsr(ibeg:iend) = nan
    pef%parsun_z(ibeg:iend,:nlevcan) = nan
    pef%parsha_z(ibeg:iend,:nlevcan) = nan
    pef%dlrad(ibeg:iend) = nan
    pef%ulrad(ibeg:iend) = nan
    pef%eflx_lh_tot(ibeg:iend) = nan
    pef%eflx_lh_tot_u(ibeg:iend) = nan
    pef%eflx_lh_tot_r(ibeg:iend) = nan
    pef%eflx_lh_grnd(ibeg:iend) = nan
    pef%eflx_soil_grnd(ibeg:iend) = nan
    pef%eflx_soil_grnd_u(ibeg:iend) = nan
    pef%eflx_soil_grnd_r(ibeg:iend) = nan
    pef%eflx_sh_tot(ibeg:iend) = nan
    pef%eflx_sh_tot_u(ibeg:iend) = nan
    pef%eflx_sh_tot_r(ibeg:iend) = nan
    pef%eflx_sh_grnd(ibeg:iend) = nan
    pef%eflx_sh_veg(ibeg:iend) = nan
    pef%eflx_lh_vege(ibeg:iend) = nan
    pef%eflx_lh_vegt(ibeg:iend) = nan
    pef%eflx_wasteheat_pft(ibeg:iend) = nan
    pef%eflx_heat_from_ac_pft(ibeg:iend) = nan
    pef%eflx_traffic_pft(ibeg:iend) = nan
    pef%eflx_anthro(ibeg:iend) = nan
    pef%cgrnd(ibeg:iend) = nan
    pef%cgrndl(ibeg:iend) = nan
    pef%cgrnds(ibeg:iend) = nan
    pef%eflx_gnet(ibeg:iend) = nan
    pef%eflx_grnd_lake(ibeg:iend) = nan
    pef%dgnetdT(ibeg:iend) = nan
    pef%eflx_lwrad_out(ibeg:iend) = nan
    pef%eflx_lwrad_net(ibeg:iend) = nan
    pef%eflx_lwrad_net_u(ibeg:iend) = nan
    pef%eflx_lwrad_net_r(ibeg:iend) = nan
    pef%netrad(ibeg:iend) = nan
    pef%fsds_vis_d(ibeg:iend) = nan
    pef%fsds_vis_i_ln(ibeg:iend) = nan
    pef%parveg_ln(ibeg:iend) = nan
    pef%fsds_nir_d(ibeg:iend) = nan
    pef%fsds_vis_i(ibeg:iend) = nan
    pef%fsds_nir_i(ibeg:iend) = nan
    pef%fsr_vis_d(ibeg:iend) = nan
    pef%fsr_nir_d(ibeg:iend) = nan
    pef%fsr_vis_i(ibeg:iend) = nan
    pef%fsr_nir_i(ibeg:iend) = nan
    pef%fsds_vis_d_ln(ibeg:iend) = nan
    pef%fsds_nir_d_ln(ibeg:iend) = nan
    pef%fsr_vis_d_ln(ibeg:iend) = nan
    pef%fsr_nir_d_ln(ibeg:iend) = nan
    pef%sabg_lyr(ibeg:iend,-nlevsno+1:1) = nan
    pef%sabg_pen(ibeg:iend) = nan
    pef%sfc_frc_aer(ibeg:iend) = nan
    pef%sfc_frc_bc(ibeg:iend) = nan
    pef%sfc_frc_oc(ibeg:iend) = nan
    pef%sfc_frc_dst(ibeg:iend) = nan
    pef%sfc_frc_aer_sno(ibeg:iend) = nan
    pef%sfc_frc_bc_sno(ibeg:iend) = nan
    pef%sfc_frc_oc_sno(ibeg:iend) = nan
    pef%sfc_frc_dst_sno(ibeg:iend) = nan
    pef%fsr_sno_vd(ibeg:iend) = nan
    pef%fsr_sno_nd(ibeg:iend) = nan
    pef%fsr_sno_vi(ibeg:iend) = nan
    pef%fsr_sno_ni(ibeg:iend) = nan
    pef%fsds_sno_vd(ibeg:iend) = nan
    pef%fsds_sno_nd(ibeg:iend) = nan
    pef%fsds_sno_vi(ibeg:iend) = nan
    pef%fsds_sno_ni(ibeg:iend) = nan
    allocate(pef%eflx_sh_snow(ibeg:iend))
    allocate(pef%eflx_sh_soil(ibeg:iend))
    allocate(pef%eflx_sh_h2osfc(ibeg:iend))
    pef%eflx_sh_snow(ibeg:iend) = nan
    pef%eflx_sh_soil(ibeg:iend) = nan
    pef%eflx_sh_h2osfc(ibeg:iend) = nan

    allocate(pef%sabg_soil(ibeg:iend))
    allocate(pef%sabg_snow(ibeg:iend))
    allocate(pef%sabg_chk(ibeg:iend))
    pef%sabg_soil(ibeg:iend) = nan
    pef%sabg_snow(ibeg:iend) = nan
    pef%sabg_chk(ibeg:iend) = nan
  end subroutine init_pft_eflux_type
  !
  ! Initialize pft momentum flux variables
  !
  subroutine init_pft_mflux_type(ibeg,iend,pmf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_mflux_type) , intent(inout) :: pmf

    allocate(pmf%taux(ibeg:iend))
    allocate(pmf%tauy(ibeg:iend))

    pmf%taux(ibeg:iend) = nan
    pmf%tauy(ibeg:iend) = nan
  end subroutine init_pft_mflux_type
  !
  ! Initialize pft water flux variables
  !
  subroutine init_pft_wflux_type(ibeg,iend,pwf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_wflux_type) , intent(inout) :: pwf

    allocate(pwf%qflx_prec_intr(ibeg:iend))
    allocate(pwf%qflx_prec_grnd(ibeg:iend))
    allocate(pwf%qflx_rain_grnd(ibeg:iend))
    allocate(pwf%qflx_snow_grnd(ibeg:iend))
    allocate(pwf%qflx_snwcp_liq(ibeg:iend))
    allocate(pwf%qflx_snwcp_ice(ibeg:iend))
    allocate(pwf%qflx_evap_veg(ibeg:iend))
    allocate(pwf%qflx_tran_veg(ibeg:iend))
    allocate(pwf%qflx_evap_can(ibeg:iend))
    allocate(pwf%qflx_evap_soi(ibeg:iend))
    allocate(pwf%qflx_evap_tot(ibeg:iend))
    allocate(pwf%qflx_evap_grnd(ibeg:iend))
    allocate(pwf%qflx_dew_grnd(ibeg:iend))
    allocate(pwf%qflx_sub_snow(ibeg:iend))
    allocate(pwf%qflx_dew_snow(ibeg:iend))

    pwf%qflx_prec_intr(ibeg:iend) = nan
    pwf%qflx_prec_grnd(ibeg:iend) = nan
    pwf%qflx_rain_grnd(ibeg:iend) = nan
    pwf%qflx_snow_grnd(ibeg:iend) = nan
    pwf%qflx_snwcp_liq(ibeg:iend) = nan
    pwf%qflx_snwcp_ice(ibeg:iend) = nan
    pwf%qflx_evap_veg(ibeg:iend) = nan
    pwf%qflx_tran_veg(ibeg:iend) = nan
    pwf%qflx_evap_can(ibeg:iend) = nan
    pwf%qflx_evap_soi(ibeg:iend) = nan
    pwf%qflx_evap_tot(ibeg:iend) = nan
    pwf%qflx_evap_grnd(ibeg:iend) = 0.0D0
    pwf%qflx_dew_grnd(ibeg:iend) = 0.0D0
    pwf%qflx_sub_snow(ibeg:iend) = 0.0D0
    pwf%qflx_dew_snow(ibeg:iend) = 0.0D0

    allocate(pwf%qflx_ev_snow(ibeg:iend))
    allocate(pwf%qflx_ev_soil(ibeg:iend))
    allocate(pwf%qflx_ev_h2osfc(ibeg:iend))
    pwf%qflx_ev_snow(ibeg:iend) = nan
    pwf%qflx_ev_soil(ibeg:iend) = nan
    pwf%qflx_ev_h2osfc(ibeg:iend) = nan
  end subroutine init_pft_wflux_type
  !
  ! Initialize pft carbon flux variables
  !
  subroutine init_pft_cflux_type(ibeg,iend,pcf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_cflux_type) , intent(inout) :: pcf

    allocate(pcf%psnsun(ibeg:iend))
    allocate(pcf%psnsha(ibeg:iend))
    allocate(pcf%psnsun_z(ibeg:iend,1:nlevcan))
    allocate(pcf%psnsha_z(ibeg:iend,1:nlevcan))
    allocate(pcf%cisun_z(ibeg:iend,1:nlevcan))
    allocate(pcf%cisha_z(ibeg:iend,1:nlevcan))
    allocate(pcf%lmrsun(ibeg:iend))
    allocate(pcf%lmrsha(ibeg:iend))
    allocate(pcf%lmrsun_z(ibeg:iend,1:nlevcan))
    allocate(pcf%lmrsha_z(ibeg:iend,1:nlevcan))
    allocate(pcf%fpsn(ibeg:iend))
    allocate(pcf%fco2(ibeg:iend))
    allocate(pcf%psnsun_wc(ibeg:iend))
    allocate(pcf%psnsha_wc(ibeg:iend))
    allocate(pcf%fpsn_wc(ibeg:iend))
    allocate(pcf%psnsun_wj(ibeg:iend))
    allocate(pcf%psnsha_wj(ibeg:iend))
    allocate(pcf%fpsn_wj(ibeg:iend))
    allocate(pcf%psnsun_wp(ibeg:iend))
    allocate(pcf%psnsha_wp(ibeg:iend))
    allocate(pcf%fpsn_wp(ibeg:iend))

    allocate(pcf%m_leafc_to_litter(ibeg:iend))
    allocate(pcf%m_frootc_to_litter(ibeg:iend))
    allocate(pcf%m_leafc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_frootc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_livestemc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_deadstemc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_livecrootc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_deadcrootc_storage_to_litter(ibeg:iend))
    allocate(pcf%m_leafc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_frootc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_livestemc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_deadstemc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_livecrootc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_deadcrootc_xfer_to_litter(ibeg:iend))
    allocate(pcf%m_livestemc_to_litter(ibeg:iend))
    allocate(pcf%m_deadstemc_to_litter(ibeg:iend))
    allocate(pcf%m_livecrootc_to_litter(ibeg:iend))
    allocate(pcf%m_deadcrootc_to_litter(ibeg:iend))
    allocate(pcf%m_gresp_storage_to_litter(ibeg:iend))
    allocate(pcf%m_gresp_xfer_to_litter(ibeg:iend))
    allocate(pcf%hrv_leafc_to_litter(ibeg:iend))             
    allocate(pcf%hrv_leafc_storage_to_litter(ibeg:iend))     
    allocate(pcf%hrv_leafc_xfer_to_litter(ibeg:iend))        
    allocate(pcf%hrv_frootc_to_litter(ibeg:iend))            
    allocate(pcf%hrv_frootc_storage_to_litter(ibeg:iend))    
    allocate(pcf%hrv_frootc_xfer_to_litter(ibeg:iend))       
    allocate(pcf%hrv_livestemc_to_litter(ibeg:iend))         
    allocate(pcf%hrv_livestemc_storage_to_litter(ibeg:iend)) 
    allocate(pcf%hrv_livestemc_xfer_to_litter(ibeg:iend))    
    allocate(pcf%hrv_deadstemc_to_prod10c(ibeg:iend))        
    allocate(pcf%hrv_deadstemc_to_prod100c(ibeg:iend))       
    allocate(pcf%hrv_deadstemc_storage_to_litter(ibeg:iend)) 
    allocate(pcf%hrv_deadstemc_xfer_to_litter(ibeg:iend))    
    allocate(pcf%hrv_livecrootc_to_litter(ibeg:iend))        
    allocate(pcf%hrv_livecrootc_storage_to_litter(ibeg:iend))
    allocate(pcf%hrv_livecrootc_xfer_to_litter(ibeg:iend))   
    allocate(pcf%hrv_deadcrootc_to_litter(ibeg:iend))        
    allocate(pcf%hrv_deadcrootc_storage_to_litter(ibeg:iend))
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(ibeg:iend))   
    allocate(pcf%hrv_gresp_storage_to_litter(ibeg:iend))     
    allocate(pcf%hrv_gresp_xfer_to_litter(ibeg:iend))        
    allocate(pcf%hrv_xsmrpool_to_atm(ibeg:iend))    
             
    ! fire related variables changed by F. Li and S. Levis           
    allocate(pcf%m_leafc_to_fire(ibeg:iend))
    allocate(pcf%m_leafc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_leafc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_livestemc_to_fire(ibeg:iend))
    allocate(pcf%m_livestemc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_livestemc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_to_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_frootc_to_fire(ibeg:iend))
    allocate(pcf%m_frootc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_frootc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_to_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_to_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_gresp_storage_to_fire(ibeg:iend))
    allocate(pcf%m_gresp_xfer_to_fire(ibeg:iend))
    allocate(pcf%m_leafc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_leafc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_leafc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livestemc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livestemc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livestemc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livestemc_to_deadstemc_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_deadstemc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_frootc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_frootc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_frootc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_livecrootc_to_deadcrootc_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_to_litter_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_deadcrootc_xfer_to_litter_fire(ibeg:iend))
    allocate(pcf%m_gresp_storage_to_litter_fire(ibeg:iend))
    allocate(pcf%m_gresp_xfer_to_litter_fire(ibeg:iend))

    allocate(pcf%leafc_xfer_to_leafc(ibeg:iend))
    allocate(pcf%frootc_xfer_to_frootc(ibeg:iend))
    allocate(pcf%livestemc_xfer_to_livestemc(ibeg:iend))
    allocate(pcf%deadstemc_xfer_to_deadstemc(ibeg:iend))
    allocate(pcf%livecrootc_xfer_to_livecrootc(ibeg:iend))
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(ibeg:iend))
    allocate(pcf%leafc_to_litter(ibeg:iend))
    allocate(pcf%frootc_to_litter(ibeg:iend))
    allocate(pcf%leaf_mr(ibeg:iend))
    allocate(pcf%froot_mr(ibeg:iend))
    allocate(pcf%livestem_mr(ibeg:iend))
    allocate(pcf%livecroot_mr(ibeg:iend))
    allocate(pcf%grain_mr(ibeg:iend))
    allocate(pcf%leaf_curmr(ibeg:iend))
    allocate(pcf%froot_curmr(ibeg:iend))
    allocate(pcf%livestem_curmr(ibeg:iend))
    allocate(pcf%livecroot_curmr(ibeg:iend))
    allocate(pcf%grain_curmr(ibeg:iend))
    allocate(pcf%leaf_xsmr(ibeg:iend))
    allocate(pcf%froot_xsmr(ibeg:iend))
    allocate(pcf%livestem_xsmr(ibeg:iend))
    allocate(pcf%livecroot_xsmr(ibeg:iend))
    allocate(pcf%grain_xsmr(ibeg:iend))
    allocate(pcf%psnsun_to_cpool(ibeg:iend))
    allocate(pcf%psnshade_to_cpool(ibeg:iend))
    allocate(pcf%cpool_to_xsmrpool(ibeg:iend))
    allocate(pcf%cpool_to_leafc(ibeg:iend))
    allocate(pcf%cpool_to_leafc_storage(ibeg:iend))
    allocate(pcf%cpool_to_frootc(ibeg:iend))
    allocate(pcf%cpool_to_frootc_storage(ibeg:iend))
    allocate(pcf%cpool_to_livestemc(ibeg:iend))
    allocate(pcf%cpool_to_livestemc_storage(ibeg:iend))
    allocate(pcf%cpool_to_deadstemc(ibeg:iend))
    allocate(pcf%cpool_to_deadstemc_storage(ibeg:iend))
    allocate(pcf%cpool_to_livecrootc(ibeg:iend))
    allocate(pcf%cpool_to_livecrootc_storage(ibeg:iend))
    allocate(pcf%cpool_to_deadcrootc(ibeg:iend))
    allocate(pcf%cpool_to_deadcrootc_storage(ibeg:iend))
    allocate(pcf%cpool_to_gresp_storage(ibeg:iend))
    allocate(pcf%cpool_leaf_gr(ibeg:iend))
    allocate(pcf%cpool_leaf_storage_gr(ibeg:iend))
    allocate(pcf%transfer_leaf_gr(ibeg:iend))
    allocate(pcf%cpool_froot_gr(ibeg:iend))
    allocate(pcf%cpool_froot_storage_gr(ibeg:iend))
    allocate(pcf%transfer_froot_gr(ibeg:iend))
    allocate(pcf%cpool_livestem_gr(ibeg:iend))
    allocate(pcf%cpool_livestem_storage_gr(ibeg:iend))
    allocate(pcf%transfer_livestem_gr(ibeg:iend))
    allocate(pcf%cpool_deadstem_gr(ibeg:iend))
    allocate(pcf%cpool_deadstem_storage_gr(ibeg:iend))
    allocate(pcf%transfer_deadstem_gr(ibeg:iend))
    allocate(pcf%cpool_livecroot_gr(ibeg:iend))
    allocate(pcf%cpool_livecroot_storage_gr(ibeg:iend))
    allocate(pcf%transfer_livecroot_gr(ibeg:iend))
    allocate(pcf%cpool_deadcroot_gr(ibeg:iend))
    allocate(pcf%cpool_deadcroot_storage_gr(ibeg:iend))
    allocate(pcf%transfer_deadcroot_gr(ibeg:iend))
    allocate(pcf%leafc_storage_to_xfer(ibeg:iend))
    allocate(pcf%frootc_storage_to_xfer(ibeg:iend))
    allocate(pcf%livestemc_storage_to_xfer(ibeg:iend))
    allocate(pcf%deadstemc_storage_to_xfer(ibeg:iend))
    allocate(pcf%livecrootc_storage_to_xfer(ibeg:iend))
    allocate(pcf%deadcrootc_storage_to_xfer(ibeg:iend))
    allocate(pcf%gresp_storage_to_xfer(ibeg:iend))
    allocate(pcf%livestemc_to_deadstemc(ibeg:iend))
    allocate(pcf%livecrootc_to_deadcrootc(ibeg:iend))
    allocate(pcf%gpp(ibeg:iend))
    allocate(pcf%mr(ibeg:iend))
    allocate(pcf%current_gr(ibeg:iend))
    allocate(pcf%transfer_gr(ibeg:iend))
    allocate(pcf%storage_gr(ibeg:iend))
    allocate(pcf%gr(ibeg:iend))
    allocate(pcf%ar(ibeg:iend))
    allocate(pcf%rr(ibeg:iend))
    allocate(pcf%npp(ibeg:iend))
    allocate(pcf%agnpp(ibeg:iend))
    allocate(pcf%bgnpp(ibeg:iend))
    allocate(pcf%litfall(ibeg:iend))
    allocate(pcf%vegfire(ibeg:iend))
    allocate(pcf%wood_harvestc(ibeg:iend))
    allocate(pcf%pft_cinputs(ibeg:iend))
    allocate(pcf%pft_coutputs(ibeg:iend))
    allocate(pcf%pft_fire_closs(ibeg:iend))
    if ( crop_prog )then
      allocate(pcf%xsmrpool_to_atm(ibeg:iend))
      allocate(pcf%grainc_xfer_to_grainc(ibeg:iend))
      allocate(pcf%livestemc_to_litter(ibeg:iend))
      allocate(pcf%grainc_to_food(ibeg:iend))
      allocate(pcf%cpool_to_grainc(ibeg:iend))
      allocate(pcf%cpool_to_grainc_storage(ibeg:iend))
      allocate(pcf%cpool_grain_gr(ibeg:iend))
      allocate(pcf%cpool_grain_storage_gr(ibeg:iend))
      allocate(pcf%transfer_grain_gr(ibeg:iend))
      allocate(pcf%grainc_storage_to_xfer(ibeg:iend))
    end if
#ifdef CN
    allocate(pcf%frootc_alloc(ibeg:iend))
    allocate(pcf%frootc_loss(ibeg:iend))
    allocate(pcf%leafc_alloc(ibeg:iend))
    allocate(pcf%leafc_loss(ibeg:iend))
    allocate(pcf%woodc_alloc(ibeg:iend))
    allocate(pcf%woodc_loss(ibeg:iend))
#endif
#ifdef LCH4
    allocate(pcf%tempavg_agnpp(ibeg:iend))
    allocate(pcf%tempavg_bgnpp(ibeg:iend))
    allocate(pcf%annavg_agnpp(ibeg:iend))
    allocate(pcf%annavg_bgnpp(ibeg:iend))
#endif

    pcf%psnsun(ibeg:iend) = nan
    pcf%psnsha(ibeg:iend) = nan
    pcf%psnsun_z(ibeg:iend,:nlevcan) = nan
    pcf%psnsha_z(ibeg:iend,:nlevcan) = nan
    pcf%cisun_z(ibeg:iend,:nlevcan) = nan
    pcf%cisha_z(ibeg:iend,:nlevcan) = nan
    pcf%lmrsun(ibeg:iend) = nan
    pcf%lmrsha(ibeg:iend) = nan
    pcf%lmrsun_z(ibeg:iend,:nlevcan) = nan
    pcf%lmrsha_z(ibeg:iend,:nlevcan) = nan
    pcf%fpsn(ibeg:iend) = spval
    pcf%fco2(ibeg:iend) = 0.D0
    pcf%psnsun_wc(ibeg:iend) = nan
    pcf%psnsha_wc(ibeg:iend) = nan
    pcf%fpsn_wc(ibeg:iend) = nan
    pcf%psnsun_wj(ibeg:iend) = nan
    pcf%psnsha_wj(ibeg:iend) = nan
    pcf%fpsn_wj(ibeg:iend) = nan
    pcf%psnsun_wp(ibeg:iend) = nan
    pcf%psnsha_wp(ibeg:iend) = nan
    pcf%fpsn_wp(ibeg:iend) = nan

    pcf%m_leafc_to_litter(ibeg:iend) = nan
    pcf%m_frootc_to_litter(ibeg:iend) = nan
    pcf%m_leafc_storage_to_litter(ibeg:iend) = nan
    pcf%m_frootc_storage_to_litter(ibeg:iend) = nan
    pcf%m_livestemc_storage_to_litter(ibeg:iend) = nan
    pcf%m_deadstemc_storage_to_litter(ibeg:iend) = nan
    pcf%m_livecrootc_storage_to_litter(ibeg:iend) = nan
    pcf%m_deadcrootc_storage_to_litter(ibeg:iend) = nan
    pcf%m_leafc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_frootc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_livestemc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_deadstemc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_livecrootc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_deadcrootc_xfer_to_litter(ibeg:iend) = nan
    pcf%m_livestemc_to_litter(ibeg:iend) = nan
    pcf%m_deadstemc_to_litter(ibeg:iend) = nan
    pcf%m_livecrootc_to_litter(ibeg:iend) = nan
    pcf%m_deadcrootc_to_litter(ibeg:iend) = nan
    pcf%m_gresp_storage_to_litter(ibeg:iend) = nan
    pcf%m_gresp_xfer_to_litter(ibeg:iend) = nan
    pcf%hrv_leafc_to_litter(ibeg:iend) = nan             
    pcf%hrv_leafc_storage_to_litter(ibeg:iend) = nan     
    pcf%hrv_leafc_xfer_to_litter(ibeg:iend) = nan        
    pcf%hrv_frootc_to_litter(ibeg:iend) = nan            
    pcf%hrv_frootc_storage_to_litter(ibeg:iend) = nan    
    pcf%hrv_frootc_xfer_to_litter(ibeg:iend) = nan       
    pcf%hrv_livestemc_to_litter(ibeg:iend) = nan         
    pcf%hrv_livestemc_storage_to_litter(ibeg:iend) = nan 
    pcf%hrv_livestemc_xfer_to_litter(ibeg:iend) = nan    
    pcf%hrv_deadstemc_to_prod10c(ibeg:iend) = nan        
    pcf%hrv_deadstemc_to_prod100c(ibeg:iend) = nan       
    pcf%hrv_deadstemc_storage_to_litter(ibeg:iend) = nan 
    pcf%hrv_deadstemc_xfer_to_litter(ibeg:iend) = nan    
    pcf%hrv_livecrootc_to_litter(ibeg:iend) = nan        
    pcf%hrv_livecrootc_storage_to_litter(ibeg:iend) = nan
    pcf%hrv_livecrootc_xfer_to_litter(ibeg:iend) = nan   
    pcf%hrv_deadcrootc_to_litter(ibeg:iend) = nan        
    pcf%hrv_deadcrootc_storage_to_litter(ibeg:iend) = nan
    pcf%hrv_deadcrootc_xfer_to_litter(ibeg:iend) = nan   
    pcf%hrv_gresp_storage_to_litter(ibeg:iend) = nan     
    pcf%hrv_gresp_xfer_to_litter(ibeg:iend) = nan        
    pcf%hrv_xsmrpool_to_atm(ibeg:iend) = nan            
     
    ! fire variable changed by F. Li and S. Levis
    pcf%m_leafc_to_fire(ibeg:iend) = nan
    pcf%m_leafc_storage_to_fire(ibeg:iend) = nan
    pcf%m_leafc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_livestemc_to_fire(ibeg:iend) = nan
    pcf%m_livestemc_storage_to_fire(ibeg:iend) = nan
    pcf%m_livestemc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_deadstemc_to_fire(ibeg:iend) = nan
    pcf%m_deadstemc_storage_to_fire(ibeg:iend) = nan
    pcf%m_deadstemc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_frootc_to_fire(ibeg:iend) = nan
    pcf%m_frootc_storage_to_fire(ibeg:iend) = nan
    pcf%m_frootc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_livecrootc_to_fire(ibeg:iend) = nan
    pcf%m_livecrootc_storage_to_fire(ibeg:iend) = nan
    pcf%m_livecrootc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_to_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_storage_to_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_xfer_to_fire(ibeg:iend) = nan
    pcf%m_gresp_storage_to_fire(ibeg:iend) = nan
    pcf%m_gresp_xfer_to_fire(ibeg:iend) = nan
    
    pcf%m_leafc_to_litter_fire(ibeg:iend) = nan
    pcf%m_leafc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_leafc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_livestemc_to_litter_fire(ibeg:iend) = nan
    pcf%m_livestemc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_livestemc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_livestemc_to_deadstemc_fire(ibeg:iend) = nan
    pcf%m_deadstemc_to_litter_fire(ibeg:iend) = nan
    pcf%m_deadstemc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_deadstemc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_frootc_to_litter_fire(ibeg:iend) = nan
    pcf%m_frootc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_frootc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_livecrootc_to_litter_fire(ibeg:iend) = nan
    pcf%m_livecrootc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_livecrootc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_livecrootc_to_deadcrootc_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_to_litter_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_deadcrootc_xfer_to_litter_fire(ibeg:iend) = nan
    pcf%m_gresp_storage_to_litter_fire(ibeg:iend) = nan
    pcf%m_gresp_xfer_to_litter_fire(ibeg:iend) = nan

    pcf%leafc_xfer_to_leafc(ibeg:iend) = nan
    pcf%frootc_xfer_to_frootc(ibeg:iend) = nan
    pcf%livestemc_xfer_to_livestemc(ibeg:iend) = nan
    pcf%deadstemc_xfer_to_deadstemc(ibeg:iend) = nan
    pcf%livecrootc_xfer_to_livecrootc(ibeg:iend) = nan
    pcf%deadcrootc_xfer_to_deadcrootc(ibeg:iend) = nan
    pcf%leafc_to_litter(ibeg:iend) = nan
    pcf%frootc_to_litter(ibeg:iend) = nan
    pcf%leaf_mr(ibeg:iend) = nan
    pcf%froot_mr(ibeg:iend) = nan
    pcf%livestem_mr(ibeg:iend) = nan
    pcf%livecroot_mr(ibeg:iend) = nan
    pcf%grain_mr(ibeg:iend) = nan
    pcf%leaf_curmr(ibeg:iend) = nan
    pcf%froot_curmr(ibeg:iend) = nan
    pcf%livestem_curmr(ibeg:iend) = nan
    pcf%livecroot_curmr(ibeg:iend) = nan
    pcf%grain_curmr(ibeg:iend) = nan
    pcf%leaf_xsmr(ibeg:iend) = nan
    pcf%froot_xsmr(ibeg:iend) = nan
    pcf%livestem_xsmr(ibeg:iend) = nan
    pcf%livecroot_xsmr(ibeg:iend) = nan
    pcf%grain_xsmr(ibeg:iend) = nan
    pcf%psnsun_to_cpool(ibeg:iend) = nan
    pcf%psnshade_to_cpool(ibeg:iend) = nan
    pcf%cpool_to_xsmrpool(ibeg:iend) = nan
    pcf%cpool_to_leafc(ibeg:iend) = nan
    pcf%cpool_to_leafc_storage(ibeg:iend) = nan
    pcf%cpool_to_frootc(ibeg:iend) = nan
    pcf%cpool_to_frootc_storage(ibeg:iend) = nan
    pcf%cpool_to_livestemc(ibeg:iend) = nan
    pcf%cpool_to_livestemc_storage(ibeg:iend) = nan
    pcf%cpool_to_deadstemc(ibeg:iend) = nan
    pcf%cpool_to_deadstemc_storage(ibeg:iend) = nan
    pcf%cpool_to_livecrootc(ibeg:iend) = nan
    pcf%cpool_to_livecrootc_storage(ibeg:iend) = nan
    pcf%cpool_to_deadcrootc(ibeg:iend) = nan
    pcf%cpool_to_deadcrootc_storage(ibeg:iend) = nan
    pcf%cpool_to_gresp_storage(ibeg:iend) = nan
    pcf%cpool_leaf_gr(ibeg:iend) = nan
    pcf%cpool_leaf_storage_gr(ibeg:iend) = nan
    pcf%transfer_leaf_gr(ibeg:iend) = nan
    pcf%cpool_froot_gr(ibeg:iend) = nan
    pcf%cpool_froot_storage_gr(ibeg:iend) = nan
    pcf%transfer_froot_gr(ibeg:iend) = nan
    pcf%cpool_livestem_gr(ibeg:iend) = nan
    pcf%cpool_livestem_storage_gr(ibeg:iend) = nan
    pcf%transfer_livestem_gr(ibeg:iend) = nan
    pcf%cpool_deadstem_gr(ibeg:iend) = nan
    pcf%cpool_deadstem_storage_gr(ibeg:iend) = nan
    pcf%transfer_deadstem_gr(ibeg:iend) = nan
    pcf%cpool_livecroot_gr(ibeg:iend) = nan
    pcf%cpool_livecroot_storage_gr(ibeg:iend) = nan
    pcf%transfer_livecroot_gr(ibeg:iend) = nan
    pcf%cpool_deadcroot_gr(ibeg:iend) = nan
    pcf%cpool_deadcroot_storage_gr(ibeg:iend) = nan
    pcf%transfer_deadcroot_gr(ibeg:iend) = nan
    pcf%leafc_storage_to_xfer(ibeg:iend) = nan
    pcf%frootc_storage_to_xfer(ibeg:iend) = nan
    pcf%livestemc_storage_to_xfer(ibeg:iend) = nan
    pcf%deadstemc_storage_to_xfer(ibeg:iend) = nan
    pcf%livecrootc_storage_to_xfer(ibeg:iend) = nan
    pcf%deadcrootc_storage_to_xfer(ibeg:iend) = nan
    pcf%gresp_storage_to_xfer(ibeg:iend) = nan
    pcf%livestemc_to_deadstemc(ibeg:iend) = nan
    pcf%livecrootc_to_deadcrootc(ibeg:iend) = nan
    pcf%gpp(ibeg:iend) = nan
    pcf%mr(ibeg:iend) = nan
    pcf%current_gr(ibeg:iend) = nan
    pcf%transfer_gr(ibeg:iend) = nan
    pcf%storage_gr(ibeg:iend) = nan
    pcf%gr(ibeg:iend) = nan
    pcf%ar(ibeg:iend) = nan
    pcf%rr(ibeg:iend) = nan
    pcf%npp(ibeg:iend) = nan
    pcf%agnpp(ibeg:iend) = nan
    pcf%bgnpp(ibeg:iend) = nan
    pcf%litfall(ibeg:iend) = nan
    pcf%vegfire(ibeg:iend) = nan
    pcf%wood_harvestc(ibeg:iend) = nan
    pcf%pft_cinputs(ibeg:iend) = nan
    pcf%pft_coutputs(ibeg:iend) = nan
    pcf%pft_fire_closs(ibeg:iend) = nan
    if ( crop_prog )then
      pcf%xsmrpool_to_atm(ibeg:iend)         = nan
      pcf%grainc_xfer_to_grainc(ibeg:iend)   = nan
      pcf%livestemc_to_litter(ibeg:iend)     = nan
      pcf%grainc_to_food(ibeg:iend)          = nan
      pcf%cpool_to_grainc(ibeg:iend)         = nan
      pcf%cpool_to_grainc_storage(ibeg:iend) = nan
      pcf%cpool_grain_gr(ibeg:iend)          = nan
      pcf%cpool_grain_storage_gr(ibeg:iend)  = nan
      pcf%transfer_grain_gr(ibeg:iend)       = nan
      pcf%grainc_storage_to_xfer(ibeg:iend)  = nan
    end if
#if (defined CN)
    pcf%frootc_alloc(ibeg:iend) = nan
    pcf%frootc_loss(ibeg:iend) = nan
    pcf%leafc_alloc(ibeg:iend) = nan
    pcf%leafc_loss(ibeg:iend) = nan
    pcf%woodc_alloc(ibeg:iend) = nan
    pcf%woodc_loss(ibeg:iend) = nan
#endif
#if (defined LCH4)
    pcf%tempavg_agnpp(ibeg:iend) = spval ! For back-compatibility
    pcf%tempavg_bgnpp(ibeg:iend) = spval ! For back-compatibility
    pcf%annavg_agnpp(ibeg:iend) = spval ! To detect first year
    pcf%annavg_bgnpp(ibeg:iend) = spval ! To detect first year
#endif
  end subroutine init_pft_cflux_type
  !
  ! Initialize pft nitrogen flux variables
  !
  subroutine init_pft_nflux_type(ibeg,iend,pnf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_nflux_type) , intent(inout) :: pnf

    allocate(pnf%m_leafn_to_litter(ibeg:iend))
    allocate(pnf%m_frootn_to_litter(ibeg:iend))
    allocate(pnf%m_leafn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_frootn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_livestemn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_deadstemn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_livecrootn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_deadcrootn_storage_to_litter(ibeg:iend))
    allocate(pnf%m_leafn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_frootn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_livestemn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_deadstemn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_livecrootn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_deadcrootn_xfer_to_litter(ibeg:iend))
    allocate(pnf%m_livestemn_to_litter(ibeg:iend))
    allocate(pnf%m_deadstemn_to_litter(ibeg:iend))
    allocate(pnf%m_livecrootn_to_litter(ibeg:iend))
    allocate(pnf%m_deadcrootn_to_litter(ibeg:iend))
    allocate(pnf%m_retransn_to_litter(ibeg:iend))
    allocate(pnf%hrv_leafn_to_litter(ibeg:iend))             
    allocate(pnf%hrv_frootn_to_litter(ibeg:iend))            
    allocate(pnf%hrv_leafn_storage_to_litter(ibeg:iend))     
    allocate(pnf%hrv_frootn_storage_to_litter(ibeg:iend))    
    allocate(pnf%hrv_livestemn_storage_to_litter(ibeg:iend)) 
    allocate(pnf%hrv_deadstemn_storage_to_litter(ibeg:iend)) 
    allocate(pnf%hrv_livecrootn_storage_to_litter(ibeg:iend))
    allocate(pnf%hrv_deadcrootn_storage_to_litter(ibeg:iend))
    allocate(pnf%hrv_leafn_xfer_to_litter(ibeg:iend))        
    allocate(pnf%hrv_frootn_xfer_to_litter(ibeg:iend))       
    allocate(pnf%hrv_livestemn_xfer_to_litter(ibeg:iend))    
    allocate(pnf%hrv_deadstemn_xfer_to_litter(ibeg:iend))    
    allocate(pnf%hrv_livecrootn_xfer_to_litter(ibeg:iend))   
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(ibeg:iend))   
    allocate(pnf%hrv_livestemn_to_litter(ibeg:iend))         
    allocate(pnf%hrv_deadstemn_to_prod10n(ibeg:iend))        
    allocate(pnf%hrv_deadstemn_to_prod100n(ibeg:iend))       
    allocate(pnf%hrv_livecrootn_to_litter(ibeg:iend))        
    allocate(pnf%hrv_deadcrootn_to_litter(ibeg:iend))        
    allocate(pnf%hrv_retransn_to_litter(ibeg:iend))   
           
    ! fire variables changed by F. Li and S. Levis            
    allocate(pnf%m_leafn_to_fire(ibeg:iend))
    allocate(pnf%m_leafn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_leafn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_livestemn_to_fire(ibeg:iend))
    allocate(pnf%m_livestemn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_livestemn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_to_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_frootn_to_fire(ibeg:iend))
    allocate(pnf%m_frootn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_frootn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_to_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_to_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_storage_to_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_xfer_to_fire(ibeg:iend))
    allocate(pnf%m_retransn_to_fire(ibeg:iend))

    allocate(pnf%m_leafn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_leafn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_leafn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livestemn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livestemn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livestemn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livestemn_to_deadstemn_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_deadstemn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_frootn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_frootn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_frootn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_livecrootn_to_deadcrootn_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_to_litter_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_storage_to_litter_fire(ibeg:iend))
    allocate(pnf%m_deadcrootn_xfer_to_litter_fire(ibeg:iend))
    allocate(pnf%m_retransn_to_litter_fire(ibeg:iend))

    allocate(pnf%leafn_xfer_to_leafn(ibeg:iend))
    allocate(pnf%frootn_xfer_to_frootn(ibeg:iend))
    allocate(pnf%livestemn_xfer_to_livestemn(ibeg:iend))
    allocate(pnf%deadstemn_xfer_to_deadstemn(ibeg:iend))
    allocate(pnf%livecrootn_xfer_to_livecrootn(ibeg:iend))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(ibeg:iend))
    allocate(pnf%leafn_to_litter(ibeg:iend))
    allocate(pnf%leafn_to_retransn(ibeg:iend))
    allocate(pnf%frootn_to_retransn(ibeg:iend))
    allocate(pnf%frootn_to_litter(ibeg:iend))
    allocate(pnf%retransn_to_npool(ibeg:iend))
    allocate(pnf%sminn_to_npool(ibeg:iend))
    allocate(pnf%npool_to_leafn(ibeg:iend))
    allocate(pnf%npool_to_leafn_storage(ibeg:iend))
    allocate(pnf%npool_to_frootn(ibeg:iend))
    allocate(pnf%npool_to_frootn_storage(ibeg:iend))
    allocate(pnf%npool_to_livestemn(ibeg:iend))
    allocate(pnf%npool_to_livestemn_storage(ibeg:iend))
    allocate(pnf%npool_to_deadstemn(ibeg:iend))
    allocate(pnf%npool_to_deadstemn_storage(ibeg:iend))
    allocate(pnf%npool_to_livecrootn(ibeg:iend))
    allocate(pnf%npool_to_livecrootn_storage(ibeg:iend))
    allocate(pnf%npool_to_deadcrootn(ibeg:iend))
    allocate(pnf%npool_to_deadcrootn_storage(ibeg:iend))
    allocate(pnf%leafn_storage_to_xfer(ibeg:iend))
    allocate(pnf%frootn_storage_to_xfer(ibeg:iend))
    allocate(pnf%livestemn_storage_to_xfer(ibeg:iend))
    allocate(pnf%deadstemn_storage_to_xfer(ibeg:iend))
    allocate(pnf%livecrootn_storage_to_xfer(ibeg:iend))
    allocate(pnf%deadcrootn_storage_to_xfer(ibeg:iend))
    allocate(pnf%livestemn_to_deadstemn(ibeg:iend))
    allocate(pnf%livestemn_to_retransn(ibeg:iend))
    allocate(pnf%livecrootn_to_deadcrootn(ibeg:iend))
    allocate(pnf%livecrootn_to_retransn(ibeg:iend))
    allocate(pnf%ndeploy(ibeg:iend))
    allocate(pnf%pft_ninputs(ibeg:iend))
    allocate(pnf%pft_noutputs(ibeg:iend))
    allocate(pnf%wood_harvestn(ibeg:iend))
    allocate(pnf%pft_fire_nloss(ibeg:iend))
    if ( crop_prog )then
      allocate(pnf%grainn_xfer_to_grainn(ibeg:iend))
      allocate(pnf%livestemn_to_litter(ibeg:iend))
      allocate(pnf%grainn_to_food(ibeg:iend))
      allocate(pnf%npool_to_grainn(ibeg:iend))
      allocate(pnf%npool_to_grainn_storage(ibeg:iend))
      allocate(pnf%grainn_storage_to_xfer(ibeg:iend))
      allocate(pnf%fert(ibeg:iend))
      allocate(pnf%soyfixn(ibeg:iend))
    end if

    pnf%m_leafn_to_litter(ibeg:iend) = nan
    pnf%m_frootn_to_litter(ibeg:iend) = nan
    pnf%m_leafn_storage_to_litter(ibeg:iend) = nan
    pnf%m_frootn_storage_to_litter(ibeg:iend) = nan
    pnf%m_livestemn_storage_to_litter(ibeg:iend) = nan
    pnf%m_deadstemn_storage_to_litter(ibeg:iend) = nan
    pnf%m_livecrootn_storage_to_litter(ibeg:iend) = nan
    pnf%m_deadcrootn_storage_to_litter(ibeg:iend) = nan
    pnf%m_leafn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_frootn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_livestemn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_deadstemn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_livecrootn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_deadcrootn_xfer_to_litter(ibeg:iend) = nan
    pnf%m_livestemn_to_litter(ibeg:iend) = nan
    pnf%m_deadstemn_to_litter(ibeg:iend) = nan
    pnf%m_livecrootn_to_litter(ibeg:iend) = nan
    pnf%m_deadcrootn_to_litter(ibeg:iend) = nan
    pnf%m_retransn_to_litter(ibeg:iend) = nan
    pnf%hrv_leafn_to_litter(ibeg:iend) = nan             
    pnf%hrv_frootn_to_litter(ibeg:iend) = nan            
    pnf%hrv_leafn_storage_to_litter(ibeg:iend) = nan     
    pnf%hrv_frootn_storage_to_litter(ibeg:iend) = nan    
    pnf%hrv_livestemn_storage_to_litter(ibeg:iend) = nan 
    pnf%hrv_deadstemn_storage_to_litter(ibeg:iend) = nan 
    pnf%hrv_livecrootn_storage_to_litter(ibeg:iend) = nan
    pnf%hrv_deadcrootn_storage_to_litter(ibeg:iend) = nan
    pnf%hrv_leafn_xfer_to_litter(ibeg:iend) = nan        
    pnf%hrv_frootn_xfer_to_litter(ibeg:iend) = nan       
    pnf%hrv_livestemn_xfer_to_litter(ibeg:iend) = nan    
    pnf%hrv_deadstemn_xfer_to_litter(ibeg:iend) = nan    
    pnf%hrv_livecrootn_xfer_to_litter(ibeg:iend) = nan   
    pnf%hrv_deadcrootn_xfer_to_litter(ibeg:iend) = nan   
    pnf%hrv_livestemn_to_litter(ibeg:iend) = nan         
    pnf%hrv_deadstemn_to_prod10n(ibeg:iend) = nan        
    pnf%hrv_deadstemn_to_prod100n(ibeg:iend) = nan       
    pnf%hrv_livecrootn_to_litter(ibeg:iend) = nan        
    pnf%hrv_deadcrootn_to_litter(ibeg:iend) = nan        
    pnf%hrv_retransn_to_litter(ibeg:iend) = nan   
        
    ! fire varibles changed by F. Li and S. Levis         
    pnf%m_leafn_to_fire(ibeg:iend) = nan
    pnf%m_leafn_storage_to_fire(ibeg:iend) = nan
    pnf%m_leafn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_livestemn_to_fire(ibeg:iend) = nan
    pnf%m_livestemn_storage_to_fire(ibeg:iend) = nan
    pnf%m_livestemn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_deadstemn_to_fire(ibeg:iend) = nan
    pnf%m_deadstemn_storage_to_fire(ibeg:iend) = nan
    pnf%m_deadstemn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_frootn_to_fire(ibeg:iend) = nan
    pnf%m_frootn_storage_to_fire(ibeg:iend) = nan
    pnf%m_frootn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_livestemn_to_fire(ibeg:iend) = nan
    pnf%m_livecrootn_storage_to_fire(ibeg:iend) = nan
    pnf%m_livecrootn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_to_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_storage_to_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_xfer_to_fire(ibeg:iend) = nan
    pnf%m_retransn_to_fire(ibeg:iend) = nan

    pnf%m_leafn_to_litter_fire(ibeg:iend) = nan
    pnf%m_leafn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_leafn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_livestemn_to_litter_fire(ibeg:iend) = nan
    pnf%m_livestemn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_livestemn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_livestemn_to_deadstemn_fire(ibeg:iend) = nan
    pnf%m_deadstemn_to_litter_fire(ibeg:iend) = nan
    pnf%m_deadstemn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_deadstemn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_frootn_to_litter_fire(ibeg:iend) = nan
    pnf%m_frootn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_frootn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_livecrootn_to_litter_fire(ibeg:iend) = nan
    pnf%m_livecrootn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_livecrootn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_livecrootn_to_deadcrootn_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_to_litter_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_storage_to_litter_fire(ibeg:iend) = nan
    pnf%m_deadcrootn_xfer_to_litter_fire(ibeg:iend) = nan
    pnf%m_retransn_to_litter_fire(ibeg:iend) = nan

    pnf%leafn_xfer_to_leafn(ibeg:iend) = nan
    pnf%frootn_xfer_to_frootn(ibeg:iend) = nan
    pnf%livestemn_xfer_to_livestemn(ibeg:iend) = nan
    pnf%deadstemn_xfer_to_deadstemn(ibeg:iend) = nan
    pnf%livecrootn_xfer_to_livecrootn(ibeg:iend) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(ibeg:iend) = nan
    pnf%leafn_to_litter(ibeg:iend) = nan
    pnf%leafn_to_retransn(ibeg:iend) = nan
    pnf%frootn_to_retransn(ibeg:iend) = nan
    pnf%frootn_to_litter(ibeg:iend) = nan
    pnf%retransn_to_npool(ibeg:iend) = nan
    pnf%sminn_to_npool(ibeg:iend) = nan
    pnf%npool_to_leafn(ibeg:iend) = nan
    pnf%npool_to_leafn_storage(ibeg:iend) = nan
    pnf%npool_to_frootn(ibeg:iend) = nan
    pnf%npool_to_frootn_storage(ibeg:iend) = nan
    pnf%npool_to_livestemn(ibeg:iend) = nan
    pnf%npool_to_livestemn_storage(ibeg:iend) = nan
    pnf%npool_to_deadstemn(ibeg:iend) = nan
    pnf%npool_to_deadstemn_storage(ibeg:iend) = nan
    pnf%npool_to_livecrootn(ibeg:iend) = nan
    pnf%npool_to_livecrootn_storage(ibeg:iend) = nan
    pnf%npool_to_deadcrootn(ibeg:iend) = nan
    pnf%npool_to_deadcrootn_storage(ibeg:iend) = nan
    pnf%leafn_storage_to_xfer(ibeg:iend) = nan
    pnf%frootn_storage_to_xfer(ibeg:iend) = nan
    pnf%livestemn_storage_to_xfer(ibeg:iend) = nan
    pnf%deadstemn_storage_to_xfer(ibeg:iend) = nan
    pnf%livecrootn_storage_to_xfer(ibeg:iend) = nan
    pnf%deadcrootn_storage_to_xfer(ibeg:iend) = nan
    pnf%livestemn_to_deadstemn(ibeg:iend) = nan
    pnf%livestemn_to_retransn(ibeg:iend) = nan
    pnf%livecrootn_to_deadcrootn(ibeg:iend) = nan
    pnf%livecrootn_to_retransn(ibeg:iend) = nan
    pnf%ndeploy(ibeg:iend) = nan
    pnf%pft_ninputs(ibeg:iend) = nan
    pnf%pft_noutputs(ibeg:iend) = nan
    pnf%wood_harvestn(ibeg:iend) = nan
    pnf%pft_fire_nloss(ibeg:iend) = nan
    if ( crop_prog )then
      pnf%grainn_xfer_to_grainn(ibeg:iend)   = nan
      pnf%livestemn_to_litter(ibeg:iend)     = nan
      pnf%grainn_to_food(ibeg:iend)          = nan
      pnf%npool_to_grainn(ibeg:iend)         = nan
      pnf%npool_to_grainn_storage(ibeg:iend) = nan
      pnf%grainn_storage_to_xfer(ibeg:iend)  = nan
      pnf%fert(ibeg:iend)                    = nan
      pnf%soyfixn(ibeg:iend)                 = nan
    end if
  end subroutine init_pft_nflux_type
  !
  ! Initialize pft VOC flux variables
  !
  subroutine init_pft_vflux_type(ibeg,iend,pvf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_vflux_type) , intent(inout) :: pvf
    integer(ik4) :: i

    if ( shr_megan_megcomps_n < 1 ) return

    allocate(pvf%vocflx_tot(ibeg:iend))
    allocate(pvf%vocflx(ibeg:iend,1:shr_megan_megcomps_n))
    allocate(pvf%Eopt_out(ibeg:iend))
    allocate(pvf%topt_out(ibeg:iend))
    allocate(pvf%alpha_out(ibeg:iend))
    allocate(pvf%cp_out(ibeg:iend))
    allocate(pvf%para_out(ibeg:iend))
    allocate(pvf%par24a_out(ibeg:iend))
    allocate(pvf%par240a_out(ibeg:iend))
    allocate(pvf%paru_out(ibeg:iend))
    allocate(pvf%par24u_out(ibeg:iend))
    allocate(pvf%par240u_out(ibeg:iend))
    allocate(pvf%gamma_out(ibeg:iend))
    allocate(pvf%gammaL_out(ibeg:iend))
    allocate(pvf%gammaT_out(ibeg:iend))
    allocate(pvf%gammaP_out(ibeg:iend))
    allocate(pvf%gammaA_out(ibeg:iend))
    allocate(pvf%gammaS_out(ibeg:iend))
    allocate(pvf%gammaC_out(ibeg:iend))

    pvf%vocflx_tot(ibeg:iend) = nan
    pvf%vocflx(ibeg:iend,1:shr_megan_megcomps_n) = nan
    pvf%Eopt_out(ibeg:iend) = nan
    pvf%topt_out(ibeg:iend) = nan
    pvf%alpha_out(ibeg:iend) = nan
    pvf%cp_out(ibeg:iend) = nan
    pvf%para_out(ibeg:iend) = nan
    pvf%par24a_out(ibeg:iend) = nan
    pvf%par240a_out(ibeg:iend) = nan
    pvf%paru_out(ibeg:iend) = nan
    pvf%par24u_out(ibeg:iend) = nan
    pvf%par240u_out(ibeg:iend) = nan
    pvf%gamma_out(ibeg:iend) = nan
    pvf%gammaL_out(ibeg:iend) = nan
    pvf%gammaT_out(ibeg:iend) = nan
    pvf%gammaP_out(ibeg:iend) = nan
    pvf%gammaA_out(ibeg:iend) = nan
    pvf%gammaS_out(ibeg:iend) = nan
    pvf%gammaC_out(ibeg:iend) = nan

    allocate(pvf%meg(shr_megan_megcomps_n))

    do i = 1 , shr_megan_megcomps_n
      allocate(pvf%meg(i)%flux_out(ibeg:iend))
      pvf%meg(i)%flux_out(ibeg:iend) = nan
    enddo
  end subroutine init_pft_vflux_type
  !
  ! Initialize pft dust flux variables
  !
  subroutine init_pft_dflux_type(ibeg,iend,pdf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_dflux_type) , intent(inout) :: pdf

    allocate(pdf%flx_mss_vrt_dst(ibeg:iend,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(ibeg:iend))
    allocate(pdf%vlc_trb(ibeg:iend,1:ndst))
    allocate(pdf%vlc_trb_1(ibeg:iend))
    allocate(pdf%vlc_trb_2(ibeg:iend))
    allocate(pdf%vlc_trb_3(ibeg:iend))
    allocate(pdf%vlc_trb_4(ibeg:iend))

    pdf%flx_mss_vrt_dst(ibeg:iend,1:ndst) = nan
    pdf%flx_mss_vrt_dst_tot(ibeg:iend) = nan
    pdf%vlc_trb(ibeg:iend,1:ndst) = nan
    pdf%vlc_trb_1(ibeg:iend) = nan
    pdf%vlc_trb_2(ibeg:iend) = nan
    pdf%vlc_trb_3(ibeg:iend) = nan
    pdf%vlc_trb_4(ibeg:iend) = nan
  end subroutine init_pft_dflux_type
  !
  ! Initialize pft dep velocity variables
  !
  subroutine init_pft_depvd_type(ibeg,iend,pdd)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (pft_depvd_type) , intent(inout) :: pdd
    integer(ik4) :: i

    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
      allocate(pdd%drydepvel(ibeg:iend,n_drydep))
      pdd%drydepvel = nan
    end if
  end subroutine init_pft_depvd_type
  !
  ! Initialize column physical state variables
  !
  subroutine init_column_pstate_type(ibeg,iend,cps)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_pstate_type) , intent(inout) :: cps

    allocate(cps%snl(ibeg:iend))      !* cannot be averaged up
    allocate(cps%isoicol(ibeg:iend))  !* cannot be averaged up
   
    ! F. Li and S. Levis
    allocate(cps%gdp_lf(ibeg:iend))  
    allocate(cps%peatf_lf(ibeg:iend))  
    allocate(cps%abm_lf(ibeg:iend))  
    allocate(cps%lgdp_col(ibeg:iend))  
    allocate(cps%lgdp1_col(ibeg:iend))
    allocate(cps%lpop_col(ibeg:iend))  

    allocate(cps%bsw(ibeg:iend,nlevgrnd))
    allocate(cps%watsat(ibeg:iend,nlevgrnd))
    allocate(cps%watfc(ibeg:iend,nlevgrnd))
    allocate(cps%watdry(ibeg:iend,nlevgrnd)) 
    allocate(cps%watopt(ibeg:iend,nlevgrnd)) 
    allocate(cps%hksat(ibeg:iend,nlevgrnd))
    allocate(cps%sucsat(ibeg:iend,nlevgrnd))
    allocate(cps%csol(ibeg:iend,nlevgrnd))
    allocate(cps%tkmg(ibeg:iend,nlevgrnd))
    allocate(cps%tkdry(ibeg:iend,nlevgrnd))
    allocate(cps%tksatu(ibeg:iend,nlevgrnd))
    allocate(cps%smpmin(ibeg:iend))
    allocate(cps%hkdepth(ibeg:iend))
    allocate(cps%wtfact(ibeg:iend))
    allocate(cps%fracice(ibeg:iend,nlevgrnd))
    allocate(cps%gwc_thr(ibeg:iend))
    allocate(cps%mss_frc_cly_vld(ibeg:iend))
    allocate(cps%mbl_bsn_fct(ibeg:iend))
    allocate(cps%do_capsnow(ibeg:iend))
    allocate(cps%snow_depth(ibeg:iend))
    allocate(cps%snowdp(ibeg:iend))
    allocate(cps%frac_sno (ibeg:iend))
    allocate(cps%zi(ibeg:iend,-nlevsno+0:nlevgrnd))
    allocate(cps%dz(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cps%z (ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cps%frac_iceold(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cps%imelt(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cps%eff_porosity(ibeg:iend,nlevgrnd))
    allocate(cps%emg(ibeg:iend))
    allocate(cps%z0mg(ibeg:iend))
    allocate(cps%z0hg(ibeg:iend))
    allocate(cps%z0qg(ibeg:iend))
    allocate(cps%htvp(ibeg:iend))
    allocate(cps%beta(ibeg:iend))
    allocate(cps%zii(ibeg:iend))
    allocate(cps%albgrd(ibeg:iend,numrad))
    allocate(cps%albgri(ibeg:iend,numrad))
    allocate(cps%rootr_column(ibeg:iend,nlevgrnd))
    allocate(cps%rootfr_road_perv(ibeg:iend,nlevgrnd))
    allocate(cps%rootr_road_perv(ibeg:iend,nlevgrnd))
    allocate(cps%wf(ibeg:iend))
    allocate(cps%wf2(ibeg:iend))
!   allocate(cps%xirrig(ibeg:iend))
    allocate(cps%max_dayl(ibeg:iend))
    allocate(cps%soilpsi(ibeg:iend,nlevgrnd))
    allocate(cps%decl(ibeg:iend))
    allocate(cps%coszen(ibeg:iend))
    allocate(cps%bd(ibeg:iend,nlevgrnd))
    allocate(cps%fpi(ibeg:iend))
    allocate(cps%fpi_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cps%fpg(ibeg:iend))
    allocate(cps%annsum_counter(ibeg:iend))
    allocate(cps%cannsum_npp(ibeg:iend))
    allocate(cps%col_lag_npp(ibeg:iend))
    allocate(cps%cannavg_t2m(ibeg:iend))

    ! fire-related variables changed by F. Li and S. Levis
    allocate(cps%nfire(ibeg:iend))
    allocate(cps%farea_burned(ibeg:iend))
    allocate(cps%fsr_col(ibeg:iend))
    allocate(cps%fd_col(ibeg:iend))
    allocate(cps%cropf_col(ibeg:iend))
    allocate(cps%prec10_col(ibeg:iend))
    allocate(cps%prec60_col(ibeg:iend))
    allocate(cps%lfc(ibeg:iend))
    allocate(cps%lfc2(ibeg:iend))
    allocate(cps%trotr1_col(ibeg:iend))
    allocate(cps%trotr2_col(ibeg:iend))
    allocate(cps%dtrotr_col(ibeg:iend))
    allocate(cps%baf_crop(ibeg:iend))
    allocate(cps%baf_peatf(ibeg:iend))
    allocate(cps%fbac(ibeg:iend))
    allocate(cps%fbac1(ibeg:iend))
    allocate(cps%btran_col(ibeg:iend))
    allocate(cps%wtlf(ibeg:iend))
    allocate(cps%lfwt(ibeg:iend))

    allocate(cps%albsnd_hst(ibeg:iend,numrad))
    allocate(cps%albsni_hst(ibeg:iend,numrad))
    allocate(cps%albsod(ibeg:iend,numrad))
    allocate(cps%albsoi(ibeg:iend,numrad))
    allocate(cps%flx_absdv(ibeg:iend,-nlevsno+1:1))
    allocate(cps%flx_absdn(ibeg:iend,-nlevsno+1:1))
    allocate(cps%flx_absiv(ibeg:iend,-nlevsno+1:1))
    allocate(cps%flx_absin(ibeg:iend,-nlevsno+1:1))
    allocate(cps%snw_rds(ibeg:iend,-nlevsno+1:0))
    allocate(cps%snw_rds_top(ibeg:iend))
    allocate(cps%sno_liq_top(ibeg:iend))
    allocate(cps%mss_bcpho(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_bcphi(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_bctot(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_bc_col(ibeg:iend))
    allocate(cps%mss_bc_top(ibeg:iend))
    allocate(cps%mss_ocpho(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_ocphi(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_octot(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_oc_col(ibeg:iend))
    allocate(cps%mss_oc_top(ibeg:iend))
    allocate(cps%mss_dst1(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_dst2(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_dst3(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_dst4(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_dsttot(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_dst_col(ibeg:iend))
    allocate(cps%mss_dst_top(ibeg:iend))
    allocate(cps%h2osno_top(ibeg:iend))
    allocate(cps%mss_cnc_bcphi(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_bcpho(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocphi(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocpho(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst1(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst2(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst3(ibeg:iend,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst4(ibeg:iend,-nlevsno+1:0))
    allocate(cps%albgrd_pur(ibeg:iend,numrad))
    allocate(cps%albgri_pur(ibeg:iend,numrad))
    allocate(cps%albgrd_bc(ibeg:iend,numrad))
    allocate(cps%albgri_bc(ibeg:iend,numrad))
    allocate(cps%albgrd_oc(ibeg:iend,numrad))
    allocate(cps%albgri_oc(ibeg:iend,numrad))
    allocate(cps%albgrd_dst(ibeg:iend,numrad))
    allocate(cps%albgri_dst(ibeg:iend,numrad))
    allocate(cps%dTdz_top(ibeg:iend))
    allocate(cps%snot_top(ibeg:iend))
    ! New variables for "S" Lakes
    allocate(cps%ws(ibeg:iend))
    allocate(cps%ks(ibeg:iend))
    allocate(cps%dz_lake(ibeg:iend,nlevlak))
    allocate(cps%z_lake(ibeg:iend,nlevlak))
    allocate(cps%savedtke1(ibeg:iend))
    allocate(cps%cellsand(ibeg:iend,nlevsoi))
    allocate(cps%cellclay(ibeg:iend,nlevsoi))
    allocate(cps%cellorg(ibeg:iend,nlevsoi))
    allocate(cps%lakedepth(ibeg:iend))
    allocate(cps%etal(ibeg:iend))
    allocate(cps%lakefetch(ibeg:iend))
    allocate(cps%ust_lake(ibeg:iend))
    ! End new variables for S Lakes
#if (defined VICHYDRO)
    ! new variables for VIC hydrology 
    allocate(cps%b_infil(ibeg:iend))
    allocate(cps%dsmax(ibeg:iend))
    allocate(cps%ds(ibeg:iend))
    allocate(cps%Wsvic(ibeg:iend))
    allocate(cps%c_param(ibeg:iend))
    allocate(cps%expt(ibeg:iend, nlayer))
    allocate(cps%ksat(ibeg:iend, nlayer))
    allocate(cps%phi_s(ibeg:iend, nlayer))
    allocate(cps%depth(ibeg:iend, nlayert))
    allocate(cps%porosity(ibeg:iend, nlayer))
    allocate(cps%max_moist(ibeg:iend, nlayer))
    allocate(cps%vic_clm_fract(ibeg:iend, nlayer, nlevsoi))
#endif
#ifdef LCH4
    ! New variable for finundated parameterization
    allocate(cps%zwt0(ibeg:iend))
    allocate(cps%f0(ibeg:iend))
    allocate(cps%p3(ibeg:iend))
    ! New variable for methane
    allocate(cps%pH(ibeg:iend))
#endif
    allocate(cps%irrig_rate(ibeg:iend))
    allocate(cps%n_irrig_steps_left(ibeg:iend))
    allocate(cps%forc_pbot(ibeg:iend))
    allocate(cps%forc_rho(ibeg:iend))
    allocate(cps%glc_topo(ibeg:iend))

    allocate(cps%rf_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
                                   1:ndecomp_cascade_transitions))
    allocate(cps%pathfrac_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
                                   1:ndecomp_cascade_transitions))
    allocate(cps%nfixation_prof(ibeg:iend,1:nlevdecomp_full))
    allocate(cps%ndep_prof(ibeg:iend,1:nlevdecomp_full))
    allocate(cps%alt(ibeg:iend))
    allocate(cps%altmax(ibeg:iend))
    allocate(cps%altmax_lastyear(ibeg:iend))
    allocate(cps%alt_indx(ibeg:iend))
    allocate(cps%altmax_indx(ibeg:iend))
    allocate(cps%altmax_lastyear_indx(ibeg:iend))
    allocate(cps%som_adv_coef(ibeg:iend,1:nlevdecomp_full))
    allocate(cps%som_diffus_coef(ibeg:iend,1:nlevdecomp_full))

    cps%isoicol(ibeg:iend) = bigint

    !F. Li and S. Levis
    cps%gdp_lf(ibeg:iend) = nan
    cps%peatf_lf(ibeg:iend) = nan
    cps%abm_lf(ibeg:iend) = 13 
    cps%lgdp_col(ibeg:iend) = nan
    cps%lgdp1_col(ibeg:iend) = nan
    cps%lpop_col(ibeg:iend) = nan

    cps%bsw(ibeg:iend,1:nlevgrnd) = nan
    cps%watsat(ibeg:iend,1:nlevgrnd) = nan
    cps%watfc(ibeg:iend,1:nlevgrnd) = nan
    cps%watdry(ibeg:iend,1:nlevgrnd) = nan
    cps%watopt(ibeg:iend,1:nlevgrnd) = nan
    cps%hksat(ibeg:iend,1:nlevgrnd) = nan
    cps%sucsat(ibeg:iend,1:nlevgrnd) = nan
    cps%csol(ibeg:iend,1:nlevgrnd) = nan
    cps%tkmg(ibeg:iend,1:nlevgrnd) = nan
    cps%tkdry(ibeg:iend,1:nlevgrnd) = nan
    cps%tksatu(ibeg:iend,1:nlevgrnd) = nan
    cps%smpmin(ibeg:iend) = nan
    cps%hkdepth(ibeg:iend) = nan
    cps%wtfact(ibeg:iend) = nan
    cps%fracice(ibeg:iend,1:nlevgrnd) = nan
    cps%gwc_thr(ibeg:iend) = nan
    cps%mss_frc_cly_vld(ibeg:iend) = nan
    cps%mbl_bsn_fct(ibeg:iend) = nan
    cps%do_capsnow (ibeg:iend)= .false.
    cps%snow_depth(ibeg:iend) = nan
    cps%snowdp(ibeg:iend) = nan
    cps%frac_sno(ibeg:iend) = nan
    cps%zi(ibeg:iend,-nlevsno+0:nlevgrnd) = nan
    cps%dz(ibeg:iend,-nlevsno+1:nlevgrnd) = nan
    cps%z (ibeg:iend,-nlevsno+1:nlevgrnd) = nan
    cps%frac_iceold(ibeg:iend,-nlevsno+1:nlevgrnd) = spval
    cps%imelt(ibeg:iend,-nlevsno+1:nlevgrnd) = bigint
    cps%eff_porosity(ibeg:iend,1:nlevgrnd) = spval
    cps%emg(ibeg:iend) = nan
    cps%z0mg(ibeg:iend) = nan
    cps%z0hg(ibeg:iend) = nan
    cps%z0qg(ibeg:iend) = nan
    cps%htvp(ibeg:iend) = nan
    cps%beta(ibeg:iend) = nan
    cps%zii(ibeg:iend) = nan
    cps%albgrd(ibeg:iend,:numrad) = nan
    cps%albgri(ibeg:iend,:numrad) = nan
    cps%rootr_column(ibeg:iend,1:nlevgrnd) = spval
    cps%rootfr_road_perv(ibeg:iend,1:nlevgrnd) = nan
    cps%rootr_road_perv(ibeg:iend,1:nlevgrnd) = nan
    cps%wf(ibeg:iend) = nan
!   cps%xirrig(ibeg:iend) = 0.D0
    cps%soilpsi(ibeg:iend,1:nlevgrnd) = spval
    cps%decl(ibeg:iend) = nan
    cps%coszen(ibeg:iend) = nan
    cps%bd(ibeg:iend,1:nlevgrnd) = spval
    cps%fpi(ibeg:iend) = nan
    cps%fpi_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cps%fpg(ibeg:iend) = nan
    cps%annsum_counter(ibeg:iend) = nan
    cps%cannsum_npp(ibeg:iend) = nan
    cps%col_lag_npp(ibeg:iend) = spval
    cps%cannavg_t2m(ibeg:iend) = nan

    ! fire-related varibles changed by F. Li and S. Levis
    cps%nfire(ibeg:iend) = spval
    cps%farea_burned(ibeg:iend) = nan
    cps%btran_col(ibeg:iend) = nan
    cps%wtlf(ibeg:iend) = nan
    cps%lfwt(ibeg:iend) = nan
    cps%fsr_col(ibeg:iend) = nan
    cps%fd_col(ibeg:iend) = nan
    cps%cropf_col(ibeg:iend) = nan
    cps%baf_crop(ibeg:iend) = nan
    cps%baf_peatf(ibeg:iend) = nan
    cps%fbac(ibeg:iend) = nan
    cps%fbac1(ibeg:iend) = nan
    cps%trotr1_col(ibeg:iend) = 0.D0
    cps%trotr2_col(ibeg:iend) = 0.D0
    cps%dtrotr_col(ibeg:iend) = 0.D0
    cps%prec10_col(ibeg:iend) = nan
    cps%prec60_col(ibeg:iend) = nan
    cps%lfc(ibeg:iend) = spval
    cps%lfc2(ibeg:iend) = 0.D0

    cps%albsnd_hst(ibeg:iend,:numrad) = spval
    cps%albsni_hst(ibeg:iend,:numrad) = spval
    cps%albsod(ibeg:iend,:numrad)     = spval
    cps%albsoi(ibeg:iend,:numrad)     = spval
    cps%flx_absdv(ibeg:iend,-nlevsno+1:1) = spval
    cps%flx_absdn(ibeg:iend,-nlevsno+1:1) = spval
    cps%flx_absiv(ibeg:iend,-nlevsno+1:1) = spval
    cps%flx_absin(ibeg:iend,-nlevsno+1:1) = spval
    cps%snw_rds(ibeg:iend,-nlevsno+1:0) = nan
    cps%snw_rds_top(ibeg:iend) = nan
    cps%sno_liq_top(ibeg:iend) = nan
    cps%mss_bcpho(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_bcphi(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_bctot(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_bc_col(ibeg:iend) = nan
    cps%mss_bc_top(ibeg:iend) = nan
    cps%mss_ocpho(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_ocphi(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_octot(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_oc_col(ibeg:iend) = nan
    cps%mss_oc_top(ibeg:iend) = nan
    cps%mss_dst1(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_dst2(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_dst3(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_dst4(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_dsttot(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_dst_col(ibeg:iend) = nan
    cps%mss_dst_top(ibeg:iend) = nan
    cps%h2osno_top(ibeg:iend) = nan
    cps%mss_cnc_bcphi(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_bcpho(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_ocphi(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_ocpho(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_dst1(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_dst2(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_dst3(ibeg:iend,-nlevsno+1:0) = nan
    cps%mss_cnc_dst4(ibeg:iend,-nlevsno+1:0) = nan
    cps%albgrd_pur(ibeg:iend,:numrad) = nan
    cps%albgri_pur(ibeg:iend,:numrad) = nan
    cps%albgrd_bc(ibeg:iend,:numrad) = nan
    cps%albgri_bc(ibeg:iend,:numrad) = nan
    cps%albgrd_oc(ibeg:iend,:numrad) = nan
    cps%albgri_oc(ibeg:iend,:numrad) = nan 
    cps%albgrd_dst(ibeg:iend,:numrad) = nan
    cps%albgri_dst(ibeg:iend,:numrad) = nan
    cps%dTdz_top(ibeg:iend) = nan
    cps%snot_top(ibeg:iend) = nan
    ! New variables for "S" Lakes
    cps%ws(ibeg:iend) = nan
    cps%ks(ibeg:iend) = nan
    cps%dz_lake(ibeg:iend,1:nlevlak) = nan
    cps%z_lake(ibeg:iend,1:nlevlak) = nan
    ! Initialize to spval so that c->g averaging will be done properly
    cps%savedtke1(ibeg:iend) = spval
    cps%cellsand(ibeg:iend,1:nlevsoi) = nan
    cps%cellclay(ibeg:iend,1:nlevsoi) = nan
    cps%cellorg(ibeg:iend,1:nlevsoi) = nan
    ! Initialize to spval so that it can be a placeholder
    ! for future file input
    cps%lakedepth(ibeg:iend) = spval
    cps%etal(ibeg:iend) = nan
    cps%lakefetch(ibeg:iend) = nan
    ! Initial to spval to detect input from restart file if not arbinit
    cps%ust_lake(ibeg:iend) = spval
    ! End new variables for S Lakes
#ifdef LCH4
    cps%zwt0(ibeg:iend) = nan
    cps%f0(ibeg:iend)   = nan
    cps%p3(ibeg:iend)   = nan
    ! New variable for methane
    cps%pH(ibeg:iend)   = nan
#endif

    cps%irrig_rate(ibeg:iend) = nan
    cps%n_irrig_steps_left(ibeg:iend) = 0
    cps%forc_pbot(ibeg:iend) = nan
    cps%forc_rho(ibeg:iend) = nan
    cps%glc_topo(ibeg:iend) = nan

    cps%rf_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions) = nan
    cps%pathfrac_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions) = nan
    cps%nfixation_prof(ibeg:iend,1:nlevdecomp_full) = spval
    cps%ndep_prof(ibeg:iend,1:nlevdecomp_full) = spval
    cps%alt(ibeg:iend) = spval
    cps%altmax(ibeg:iend) = spval
    cps%altmax_lastyear(ibeg:iend) = spval
    cps%alt_indx(ibeg:iend) = bigint
    cps%altmax_indx(ibeg:iend) = bigint
    cps%altmax_lastyear_indx(ibeg:iend) = bigint
    cps%som_adv_coef(ibeg:iend,1:nlevdecomp_full) = spval
    cps%som_diffus_coef(ibeg:iend,1:nlevdecomp_full) = spval

    allocate(cps%frac_sno_eff(ibeg:iend))
    allocate(cps%topo_std(ibeg:iend))
    allocate(cps%topo_ndx(ibeg:iend))
    allocate(cps%topo_slope(ibeg:iend))
    allocate(cps%hksat_min(ibeg:iend,nlevgrnd))
    allocate(cps%frac_h2osfc(ibeg:iend))
    allocate(cps%micro_sigma(ibeg:iend))
    allocate(cps%h2osfc_thresh(ibeg:iend))
    allocate(cps%frac_h2osfc_temp(ibeg:iend))
    allocate(cps%n_melt(ibeg:iend))

    cps%frac_sno_eff(ibeg:iend) = spval
    cps%topo_std(ibeg:iend) = nan
    cps%topo_ndx(ibeg:iend) = nan
    cps%topo_slope(ibeg:iend) = nan
    cps%hksat_min(ibeg:iend,1:nlevgrnd) = nan
    cps%frac_h2osfc(ibeg:iend) = spval
    cps%micro_sigma(ibeg:iend) = nan
    cps%h2osfc_thresh(ibeg:iend) = nan
    cps%frac_h2osfc_temp(ibeg:iend) = 0.0D0
    cps%n_melt(ibeg:iend) = nan
#if (defined VICHYDRO)
    ! new variables for VIC hydrology
    cps%b_infil(ibeg:iend)  = nan
    cps%dsmax(ibeg:iend)    = nan
    cps%ds(ibeg:iend)       = nan
    cps%Wsvic(ibeg:iend)    = nan
    cps%c_param(ibeg:iend)  = nan
    cps%expt(ibeg:iend, 1:nlayer)     = nan
    cps%ksat(ibeg:iend, 1:nlayer)     = nan
    cps%phi_s(ibeg:iend, 1:nlayer)    = nan
    cps%depth(ibeg:iend, 1:nlayert)    = nan
    cps%porosity(ibeg:iend, 1:nlayer)   = nan
    cps%max_moist(ibeg:iend, 1:nlayer)   = nan
    cps%vic_clm_fract(ibeg:iend,1:nlayer,1:nlevsoi) = nan
#endif
  end subroutine init_column_pstate_type
  !
  ! Initialize column energy state variables
  !
  subroutine init_column_estate_type(ibeg,iend,ces)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_estate_type) , intent(inout) :: ces

    allocate(ces%t_grnd(ibeg:iend))
    allocate(ces%t_grnd_u(ibeg:iend))
    allocate(ces%t_grnd_r(ibeg:iend))
    allocate(ces%dt_grnd(ibeg:iend))
    allocate(ces%t_soisno(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(ces%t_soi_10cm(ibeg:iend))
    allocate(ces%tsoi17(ibeg:iend))
    allocate(ces%t_lake(ibeg:iend,1:nlevlak))
    allocate(ces%tssbef(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(ces%thv(ibeg:iend))
    allocate(ces%hc_soi(ibeg:iend))
    allocate(ces%hc_soisno(ibeg:iend))
    allocate(ces%forc_t(ibeg:iend))
    allocate(ces%forc_th(ibeg:iend))

    ces%t_grnd(ibeg:iend)    = nan
    ces%t_grnd_u(ibeg:iend)  = nan
    ces%t_grnd_r(ibeg:iend)  = nan
    ces%dt_grnd(ibeg:iend)   = nan
    ces%t_soisno(ibeg:iend,-nlevsno+1:nlevgrnd) = spval
    ces%t_soi_10cm(ibeg:iend) = spval
    ces%tsoi17(ibeg:iend) = spval
    ces%t_lake(ibeg:iend,1:nlevlak)            = nan
    ces%tssbef(ibeg:iend,-nlevsno+1:nlevgrnd)   = nan
    ces%thv(ibeg:iend)       = nan
    ces%hc_soi(ibeg:iend)    = nan
    ces%hc_soisno(ibeg:iend) = nan
    ces%forc_t(ibeg:iend) = nan
    ces%forc_th(ibeg:iend) = nan

    allocate(ces%t_h2osfc(ibeg:iend))
    allocate(ces%t_h2osfc_bef(ibeg:iend))

    ces%t_h2osfc(ibeg:iend) = spval
    ces%t_h2osfc_bef(ibeg:iend)   = nan
  end subroutine init_column_estate_type
  !
  ! Initialize column water state variables
  !
  subroutine init_column_wstate_type(ibeg,iend,cws)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_wstate_type) , intent(inout) :: cws !column water state

    allocate(cws%h2osno(ibeg:iend))
    allocate(cws%errh2osno(ibeg:iend))
    allocate(cws%snow_sources(ibeg:iend))
    allocate(cws%snow_sinks(ibeg:iend))
    allocate(cws%h2osoi_liq(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_ice(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_liqice_10cm(ibeg:iend))
    allocate(cws%h2osoi_vol(ibeg:iend,1:nlevgrnd))
    allocate(cws%h2osno_old(ibeg:iend))
    allocate(cws%qg(ibeg:iend))
    allocate(cws%dqgdT(ibeg:iend))
    allocate(cws%snowice(ibeg:iend))
    allocate(cws%snowliq(ibeg:iend))
    allocate(cws%soilalpha(ibeg:iend))
    allocate(cws%soilbeta(ibeg:iend))
    allocate(cws%soilalpha_u(ibeg:iend))
    allocate(cws%zwt(ibeg:iend))
    allocate(cws%fcov(ibeg:iend))
    allocate(cws%fsat(ibeg:iend))
    !New variable for methane code
#ifdef LCH4
    allocate(cws%finundated(ibeg:iend))
#endif
    allocate(cws%wa(ibeg:iend))
    allocate(cws%qcharge(ibeg:iend))
    allocate(cws%smp_l(ibeg:iend,1:nlevgrnd))
    allocate(cws%hk_l(ibeg:iend,1:nlevgrnd))
    ! New variables for "S" lakes
    allocate(cws%lake_icefrac(ibeg:iend,1:nlevlak))
    allocate(cws%lake_icethick(ibeg:iend))
    ! End new variables for S lakes
    allocate(cws%forc_q(ibeg:iend))
#if (defined VICHYDRO)
    allocate(cws%moist(ibeg:iend,1:nlayert))
    allocate(cws%ice(ibeg:iend,1:nlayert))
    allocate(cws%moist_vol(ibeg:iend,1:nlayert))
    allocate(cws%max_infil(ibeg:iend))
    allocate(cws%i_0(ibeg:iend))
#endif

    cws%h2osno(ibeg:iend) = nan
    cws%errh2osno(ibeg:iend) = nan
    cws%snow_sources(ibeg:iend) = nan
    cws%snow_sinks(ibeg:iend) = nan
    cws%h2osoi_liq(ibeg:iend,-nlevsno+1:nlevgrnd)= spval
    cws%h2osoi_ice(ibeg:iend,-nlevsno+1:nlevgrnd) = spval
    cws%h2osoi_liqice_10cm(ibeg:iend) = spval
    cws%h2osoi_vol(ibeg:iend,1:nlevgrnd) = spval
    cws%h2osno_old(ibeg:iend) = nan
    cws%qg(ibeg:iend) = nan
    cws%dqgdT(ibeg:iend) = nan
    cws%snowice(ibeg:iend) = nan
    cws%snowliq(ibeg:iend) = nan
    cws%soilalpha(ibeg:iend) = nan
    cws%soilbeta(ibeg:iend) = nan
    cws%soilalpha_u(ibeg:iend) = nan
    cws%zwt(ibeg:iend) = nan
    cws%fcov(ibeg:iend) = nan
    cws%fsat(ibeg:iend) = nan
#ifdef LCH4
    cws%finundated(ibeg:iend) = nan
#endif
    cws%wa(ibeg:iend) = spval
    cws%qcharge(ibeg:iend) = nan
    cws%smp_l(ibeg:iend,1:nlevgrnd) = spval
    cws%hk_l(ibeg:iend,1:nlevgrnd) = spval
    ! New variables for "S" lakes
    ! Initialize to spval for detection of whether
    ! it has not been initialized by file input
    cws%lake_icefrac(ibeg:iend,1:nlevlak) = spval
    ! and so c->g averaging will be done properly
    cws%lake_icethick(ibeg:iend) = nan
    ! End new variables for S lakes
    cws%forc_q(ibeg:iend) = nan

    allocate(cws%h2osfc(ibeg:iend))
    allocate(cws%qg_snow(ibeg:iend))
    allocate(cws%qg_soil(ibeg:iend))
    allocate(cws%qg_h2osfc(ibeg:iend))
    allocate(cws%frost_table(ibeg:iend))
    allocate(cws%zwt_perched(ibeg:iend))
    allocate(cws%int_snow(ibeg:iend))
    allocate(cws%swe_old(ibeg:iend,-nlevsno+1:0))

    cws%h2osfc(ibeg:iend)      = spval
    cws%qg_snow(ibeg:iend)     = nan
    cws%qg_soil(ibeg:iend)     = nan
    cws%qg_h2osfc(ibeg:iend)   = nan
    cws%frost_table(ibeg:iend) = spval
    cws%zwt_perched(ibeg:iend) = spval
    cws%int_snow(ibeg:iend)    = spval
    cws%swe_old(ibeg:iend,-nlevsno+1:0)= nan
#if (defined VICHYDRO)
    cws%moist(ibeg:iend,1:nlayert) = spval
    cws%ice(ibeg:iend,1:nlayert) = spval
    cws%moist_vol(ibeg:iend,1:nlayert) = spval
    cws%max_infil(ibeg:iend) = spval
    cws%i_0(ibeg:iend) = spval
#endif
  end subroutine init_column_wstate_type
  !
  ! Initialize column carbon state variables
  !
  subroutine init_column_cstate_type(ibeg,iend,ccs)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_cstate_type) , intent(inout) :: ccs

    allocate(ccs%soilc(ibeg:iend))
    allocate(ccs%cwdc(ibeg:iend))
    allocate(ccs%col_ctrunc(ibeg:iend))
    allocate(ccs%decomp_cpools_vr(ibeg:iend,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccs%decomp_cpools(ibeg:iend,1:ndecomp_pools))
    allocate(ccs%decomp_cpools_1m(ibeg:iend,1:ndecomp_pools))
    allocate(ccs%col_ctrunc_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(ccs%seedc(ibeg:iend))
    allocate(ccs%prod10c(ibeg:iend))
    allocate(ccs%prod100c(ibeg:iend))
    allocate(ccs%totprodc(ibeg:iend))
    allocate(ccs%totlitc(ibeg:iend))
    allocate(ccs%totsomc(ibeg:iend))
    allocate(ccs%totlitc_1m(ibeg:iend))
    allocate(ccs%totsomc_1m(ibeg:iend))
    allocate(ccs%totecosysc(ibeg:iend))
    allocate(ccs%totcolc(ibeg:iend))

    !F. Li and S. Levis
    allocate(ccs%rootc_col(ibeg:iend))
    allocate(ccs%totvegc_col(ibeg:iend))
    allocate(ccs%leafc_col(ibeg:iend))
    allocate(ccs%fuelc(ibeg:iend))
    allocate(ccs%fuelc_crop(ibeg:iend))

    ccs%soilc(ibeg:iend) = nan
    ccs%cwdc(ibeg:iend) = nan
    ccs%decomp_cpools_vr(ibeg:iend,1:nlevdecomp_full,1:ndecomp_pools) = nan
    ccs%decomp_cpools(ibeg:iend,1:ndecomp_pools) = nan
    ccs%decomp_cpools_1m(ibeg:iend,1:ndecomp_pools) = nan
    ccs%col_ctrunc(ibeg:iend) = nan
    ccs%col_ctrunc_vr(ibeg:iend,1:nlevdecomp_full) = nan
    ccs%seedc(ibeg:iend) = nan
    ccs%prod10c(ibeg:iend) = nan
    ccs%prod100c(ibeg:iend) = nan
    ccs%totprodc(ibeg:iend) = nan
    ccs%totlitc(ibeg:iend) = nan
    ccs%totsomc(ibeg:iend) = nan
    ccs%totlitc_1m(ibeg:iend) = nan
    ccs%totsomc_1m(ibeg:iend) = nan
    ccs%totecosysc(ibeg:iend) = nan
    ccs%totcolc(ibeg:iend) = nan

    ccs%rootc_col(ibeg:iend) = nan
    ccs%totvegc_col(ibeg:iend) = nan
    ccs%leafc_col(ibeg:iend) = nan
    ccs%fuelc(ibeg:iend) = spval
    ccs%fuelc_crop(ibeg:iend) = nan
  end subroutine init_column_cstate_type
  !
  ! Initialize column nitrogen state variables
  !
  subroutine init_column_nstate_type(ibeg,iend,cns)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_nstate_type) , intent(inout) :: cns

    allocate(cns%decomp_npools(ibeg:iend,1:ndecomp_pools))
    allocate(cns%decomp_npools_1m(ibeg:iend,1:ndecomp_pools))
    allocate(cns%decomp_npools_vr(ibeg:iend,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cns%sminn_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cns%col_ntrunc_vr(ibeg:iend,1:nlevdecomp_full))
#ifdef NITRIF_DENITRIF
    allocate(cns%smin_no3_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cns%smin_nh4_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cns%smin_no3(ibeg:iend))
    allocate(cns%smin_nh4(ibeg:iend))
#endif
    allocate(cns%cwdn(ibeg:iend))
    allocate(cns%sminn(ibeg:iend))
    allocate(cns%col_ntrunc(ibeg:iend))
    allocate(cns%seedn(ibeg:iend))
    allocate(cns%prod10n(ibeg:iend))
    allocate(cns%prod100n(ibeg:iend))
    allocate(cns%totprodn(ibeg:iend))
    allocate(cns%totlitn(ibeg:iend))
    allocate(cns%totsomn(ibeg:iend))
    allocate(cns%totlitn_1m(ibeg:iend))
    allocate(cns%totsomn_1m(ibeg:iend))
    allocate(cns%totecosysn(ibeg:iend))
    allocate(cns%totcoln(ibeg:iend))

    cns%decomp_npools(ibeg:iend,1:ndecomp_pools) = nan
    cns%decomp_npools_1m(ibeg:iend,1:ndecomp_pools) = nan
    cns%decomp_npools_vr(ibeg:iend,1:nlevdecomp_full,1:ndecomp_pools) = nan
    cns%sminn_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cns%col_ntrunc_vr(ibeg:iend,1:nlevdecomp_full) = nan
#ifdef NITRIF_DENITRIF
    cns%smin_no3_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cns%smin_nh4_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cns%smin_no3(ibeg:iend) = nan
    cns%smin_nh4(ibeg:iend) = nan
#endif
    cns%cwdn(ibeg:iend) = nan
    cns%sminn(ibeg:iend) = nan
    cns%col_ntrunc(ibeg:iend) = nan
    cns%seedn(ibeg:iend) = nan
    cns%prod10n(ibeg:iend) = nan
    cns%prod100n(ibeg:iend) = nan
    cns%totprodn(ibeg:iend) = nan
    cns%totlitn(ibeg:iend) = nan
    cns%totsomn(ibeg:iend) = nan
    cns%totlitn_1m(ibeg:iend) = nan
    cns%totsomn_1m(ibeg:iend) = nan
    cns%totecosysn(ibeg:iend) = nan
    cns%totcoln(ibeg:iend) = nan
  end subroutine init_column_nstate_type
  !
  ! Initialize column energy flux variables
  !
  subroutine init_column_eflux_type(ibeg,iend,cef)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_eflux_type) , intent(inout) :: cef

    allocate(cef%eflx_snomelt(ibeg:iend))
    allocate(cef%eflx_snomelt_u(ibeg:iend))
    allocate(cef%eflx_snomelt_r(ibeg:iend))
    allocate(cef%eflx_impsoil(ibeg:iend))
    allocate(cef%eflx_fgr12(ibeg:iend))
    allocate(cef%eflx_fgr(ibeg:iend, 1:nlevgrnd))
    allocate(cef%eflx_building_heat(ibeg:iend))
    allocate(cef%eflx_urban_ac(ibeg:iend))
    allocate(cef%eflx_urban_heat(ibeg:iend))
    allocate(cef%eflx_bot(ibeg:iend))

    cef%eflx_snomelt(ibeg:iend)       = spval
    cef%eflx_snomelt_u(ibeg:iend)       = spval
    cef%eflx_snomelt_r(ibeg:iend)       = spval
    cef%eflx_impsoil(ibeg:iend)       = nan
    cef%eflx_fgr12(ibeg:iend)         = nan
    cef%eflx_fgr(ibeg:iend, 1:nlevgrnd) = nan
    cef%eflx_building_heat(ibeg:iend) = nan
    cef%eflx_urban_ac(ibeg:iend) = nan
    cef%eflx_urban_heat(ibeg:iend) = nan
    cef%eflx_bot(ibeg:iend) = nan
  end subroutine init_column_eflux_type
  !
  ! Initialize column water flux variables
  !
  subroutine init_column_wflux_type(ibeg,iend,cwf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_wflux_type) , intent(inout) :: cwf

    allocate(cwf%qflx_infl(ibeg:iend))
    allocate(cwf%qflx_surf(ibeg:iend))
    allocate(cwf%qflx_drain(ibeg:iend))
    allocate(cwf%qflx_top_soil(ibeg:iend))
    allocate(cwf%qflx_sl_top_soil(ibeg:iend))
    allocate(cwf%qflx_snomelt(ibeg:iend))
    allocate(cwf%qflx_qrgwl(ibeg:iend))
    allocate(cwf%qflx_runoff(ibeg:iend))
    allocate(cwf%qflx_runoff_u(ibeg:iend))
    allocate(cwf%qflx_runoff_r(ibeg:iend))
    allocate(cwf%qmelt(ibeg:iend))
    allocate(cwf%h2ocan_loss(ibeg:iend))
    allocate(cwf%qflx_rsub_sat(ibeg:iend))
    allocate(cwf%flx_bc_dep_dry(ibeg:iend))
    allocate(cwf%flx_bc_dep_wet(ibeg:iend))
    allocate(cwf%flx_bc_dep_pho(ibeg:iend))
    allocate(cwf%flx_bc_dep_phi(ibeg:iend))
    allocate(cwf%flx_bc_dep(ibeg:iend))
    allocate(cwf%flx_oc_dep_dry(ibeg:iend))
    allocate(cwf%flx_oc_dep_wet(ibeg:iend))
    allocate(cwf%flx_oc_dep_pho(ibeg:iend))
    allocate(cwf%flx_oc_dep_phi(ibeg:iend))
    allocate(cwf%flx_oc_dep(ibeg:iend))
    allocate(cwf%flx_dst_dep_dry1(ibeg:iend))
    allocate(cwf%flx_dst_dep_wet1(ibeg:iend))
    allocate(cwf%flx_dst_dep_dry2(ibeg:iend))
    allocate(cwf%flx_dst_dep_wet2(ibeg:iend))
    allocate(cwf%flx_dst_dep_dry3(ibeg:iend))
    allocate(cwf%flx_dst_dep_wet3(ibeg:iend))
    allocate(cwf%flx_dst_dep_dry4(ibeg:iend))
    allocate(cwf%flx_dst_dep_wet4(ibeg:iend))
    allocate(cwf%flx_dst_dep(ibeg:iend))
    allocate(cwf%qflx_snofrz_lyr(ibeg:iend,-nlevsno+1:0))
    allocate(cwf%qflx_snofrz_col(ibeg:iend))
    allocate(cwf%qflx_irrig(ibeg:iend))
    allocate(cwf%qflx_glcice(ibeg:iend))
    allocate(cwf%qflx_glcice_frz(ibeg:iend))
    allocate(cwf%qflx_glcice_melt(ibeg:iend))
    allocate(cwf%glc_rofi(ibeg:iend))
    allocate(cwf%glc_rofl(ibeg:iend))
    allocate(cwf%qflx_floodc(ibeg:iend))

    cwf%qflx_infl(ibeg:iend) = nan
    cwf%qflx_surf(ibeg:iend) = nan
    cwf%qflx_drain(ibeg:iend) = nan
    cwf%qflx_top_soil(ibeg:iend) = spval
    cwf%qflx_sl_top_soil(ibeg:iend) = nan
    cwf%qflx_snomelt(ibeg:iend) = nan
    cwf%qflx_qrgwl(ibeg:iend) = nan
    cwf%qflx_runoff(ibeg:iend) = nan
    cwf%qflx_runoff_u(ibeg:iend) = nan
    cwf%qflx_runoff_r(ibeg:iend) = nan
    cwf%qmelt(ibeg:iend) = nan
    cwf%h2ocan_loss(ibeg:iend) = nan
    cwf%qflx_rsub_sat(ibeg:iend) = spval
    cwf%flx_bc_dep_dry(ibeg:iend) = nan
    cwf%flx_bc_dep_wet(ibeg:iend) = nan
    cwf%flx_bc_dep_pho(ibeg:iend) = nan
    cwf%flx_bc_dep_phi(ibeg:iend) = nan
    cwf%flx_bc_dep(ibeg:iend) = nan
    cwf%flx_oc_dep_dry(ibeg:iend) = nan
    cwf%flx_oc_dep_wet(ibeg:iend) = nan
    cwf%flx_oc_dep_pho(ibeg:iend) = nan
    cwf%flx_oc_dep_phi(ibeg:iend) = nan
    cwf%flx_oc_dep(ibeg:iend) = nan
    cwf%flx_dst_dep_dry1(ibeg:iend) = nan
    cwf%flx_dst_dep_wet1(ibeg:iend) = nan
    cwf%flx_dst_dep_dry2(ibeg:iend) = nan
    cwf%flx_dst_dep_wet2(ibeg:iend) = nan
    cwf%flx_dst_dep_dry3(ibeg:iend) = nan
    cwf%flx_dst_dep_wet3(ibeg:iend) = nan
    cwf%flx_dst_dep_dry4(ibeg:iend) = nan
    cwf%flx_dst_dep_wet4(ibeg:iend) = nan
    cwf%flx_dst_dep(ibeg:iend) = nan
    cwf%qflx_snofrz_lyr(ibeg:iend,-nlevsno+1:0) = spval
    cwf%qflx_snofrz_col(ibeg:iend) = nan
    cwf%qflx_irrig(ibeg:iend)  = spval
    cwf%qflx_glcice(ibeg:iend) = nan
    cwf%qflx_glcice_frz(ibeg:iend) = nan
    cwf%qflx_glcice_melt(ibeg:iend) = spval
    cwf%glc_rofi(ibeg:iend)    = nan
    cwf%glc_rofl(ibeg:iend)    = nan
    cwf%qflx_floodc(ibeg:iend) = spval

    allocate(cwf%qflx_h2osfc_to_ice(ibeg:iend))
    allocate(cwf%qflx_h2osfc_surf(ibeg:iend))
    allocate(cwf%qflx_snow_h2osfc(ibeg:iend))
    allocate(cwf%qflx_drain_perched(ibeg:iend))
    allocate(cwf%qflx_floodc(ibeg:iend))
    allocate(cwf%qflx_snow_melt(ibeg:iend))

    cwf%qflx_h2osfc_to_ice(ibeg:iend) = spval
    cwf%qflx_h2osfc_surf(ibeg:iend) = spval
    cwf%qflx_snow_h2osfc(ibeg:iend) = nan
    cwf%qflx_drain_perched(ibeg:iend) = spval
    cwf%qflx_floodc(ibeg:iend) = spval
    cwf%qflx_snow_melt(ibeg:iend) = spval
  end subroutine init_column_wflux_type
  !
  ! Initialize column carbon flux variables
  !
  subroutine init_column_cflux_type(ibeg,iend,ccf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_cflux_type) , intent(inout) :: ccf

    allocate(ccf%hrv_deadstemc_to_prod10c(ibeg:iend))        
    allocate(ccf%hrv_deadstemc_to_prod100c(ibeg:iend))       
    allocate(ccf%m_decomp_cpools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%m_decomp_cpools_to_fire(ibeg:iend,1:ndecomp_pools))
    allocate(ccf%decomp_cascade_hr_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_hr(ibeg:iend,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cpools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%decomp_k(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions))
    allocate(ccf%t_scalar(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%w_scalar(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%hr_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%o_scalar(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%som_c_leached(ibeg:iend))
    allocate(ccf%decomp_cpools_leached(ibeg:iend,1:ndecomp_pools))
    allocate(ccf%decomp_cpools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))

    allocate(ccf%phenology_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%fire_mortality_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full))

#ifdef NITRIF_DENITRIF
    allocate(ccf%phr_vr(ibeg:iend,1:nlevdecomp_full))
#endif

#ifdef CN
    !F. Li and S. Levis
    allocate(ccf%somc_fire(ibeg:iend))
    allocate(ccf%lf_conv_cflux(ibeg:iend))
    allocate(ccf%dwt_seedc_to_leaf(ibeg:iend))
    allocate(ccf%dwt_seedc_to_deadstem(ibeg:iend))
    allocate(ccf%dwt_conv_cflux(ibeg:iend))
    allocate(ccf%dwt_prod10c_gain(ibeg:iend))
    allocate(ccf%dwt_prod100c_gain(ibeg:iend))
    allocate(ccf%dwt_frootc_to_litr_met_c(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_cel_c(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_lig_c(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%dwt_livecrootc_to_cwdc(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%dwt_deadcrootc_to_cwdc(ibeg:iend,1:nlevdecomp_full))
    allocate(ccf%dwt_closs(ibeg:iend))
    allocate(ccf%landuseflux(ibeg:iend))
    allocate(ccf%landuptake(ibeg:iend))
    allocate(ccf%prod10c_loss(ibeg:iend))
    allocate(ccf%prod100c_loss(ibeg:iend))
    allocate(ccf%product_closs(ibeg:iend))
#endif
    allocate(ccf%lithr(ibeg:iend))
    allocate(ccf%somhr(ibeg:iend))
    allocate(ccf%hr(ibeg:iend))
    allocate(ccf%sr(ibeg:iend))
    allocate(ccf%er(ibeg:iend))
    allocate(ccf%litfire(ibeg:iend))
    allocate(ccf%somfire(ibeg:iend))
    allocate(ccf%totfire(ibeg:iend))
    allocate(ccf%nep(ibeg:iend))
    allocate(ccf%nbp(ibeg:iend))
    allocate(ccf%nee(ibeg:iend))
    allocate(ccf%col_cinputs(ibeg:iend))
    allocate(ccf%col_coutputs(ibeg:iend))
    allocate(ccf%col_fire_closs(ibeg:iend))

#if (defined CN)
    allocate(ccf%cwdc_hr(ibeg:iend))
    allocate(ccf%cwdc_loss(ibeg:iend))
    allocate(ccf%litterc_loss(ibeg:iend))
#endif

    ccf%m_c_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full) = nan
    ccf%m_c_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full) = nan
    ccf%m_c_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full) = nan
    ccf%hrv_deadstemc_to_prod10c(ibeg:iend)               = nan        
    ccf%hrv_deadstemc_to_prod100c(ibeg:iend)              = nan       
    ccf%m_decomp_cpools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools)                = nan
    ccf%m_decomp_cpools_to_fire(ibeg:iend,1:ndecomp_pools)    = nan
    ccf%decomp_cascade_hr_vr(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions)                    = nan
    ccf%decomp_cascade_hr(ibeg:iend,1:ndecomp_cascade_transitions) = nan
    ccf%decomp_cascade_ctransfer_vr(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions)                    = nan
    ccf%decomp_cascade_ctransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions)                    = nan
    ccf%decomp_cpools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools)                = nan
    ccf%decomp_k(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions)                    = spval
    ! Initialize these four below to spval to allow history to not average
    ! over inactive points.
    ccf%t_scalar(ibeg:iend,1:nlevdecomp_full)  = spval
    ccf%w_scalar(ibeg:iend,1:nlevdecomp_full)  = spval
    ccf%hr_vr(ibeg:iend, 1:nlevdecomp_full)    = spval
    ccf%o_scalar(ibeg:iend, 1:nlevdecomp_full) = spval
    ccf%som_c_leached(ibeg:iend) = nan 
    ccf%decomp_cpools_leached(ibeg:iend,1:ndecomp_pools) = nan
    ccf%decomp_cpools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools) = nan

    ccf%phenology_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%phenology_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%phenology_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%gap_mortality_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%gap_mortality_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%gap_mortality_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%gap_mortality_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full)  = nan
    ccf%fire_mortality_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full) = nan
    ccf%harvest_c_to_litr_met_c(ibeg:iend, 1:nlevdecomp_full)  = nan
    ccf%harvest_c_to_litr_cel_c(ibeg:iend, 1:nlevdecomp_full)  = nan
    ccf%harvest_c_to_litr_lig_c(ibeg:iend, 1:nlevdecomp_full)  = nan
    ccf%harvest_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full)        = nan

#ifdef NITRIF_DENITRIF
    ccf%phr_vr(ibeg:iend,1:nlevdecomp_full) = nan
#endif
#if (defined CN)
    !F. Li and S. Levis
    ccf%somc_fire(ibeg:iend)              = nan
    ccf%lf_conv_cflux(ibeg:iend)          = nan
    ccf%dwt_seedc_to_leaf(ibeg:iend)      = nan
    ccf%dwt_seedc_to_deadstem(ibeg:iend)  = nan
    ccf%dwt_conv_cflux(ibeg:iend)         = nan
    ccf%dwt_prod10c_gain(ibeg:iend)       = nan
    ccf%dwt_prod100c_gain(ibeg:iend)      = nan
    ccf%dwt_frootc_to_litr_met_c(ibeg:iend,1:nlevdecomp_full)  = nan
    ccf%dwt_frootc_to_litr_cel_c(ibeg:iend,1:nlevdecomp_full)  = nan
    ccf%dwt_frootc_to_litr_lig_c(ibeg:iend,1:nlevdecomp_full)  = nan
    ccf%dwt_livecrootc_to_cwdc(ibeg:iend,1:nlevdecomp_full)    = nan
    ccf%dwt_deadcrootc_to_cwdc(ibeg:iend,1:nlevdecomp_full)    = nan
    ccf%dwt_closs(ibeg:iend)      = nan
    ccf%landuseflux(ibeg:iend)    = nan
    ccf%landuptake(ibeg:iend)     = nan
    ccf%prod10c_loss(ibeg:iend)   = nan
    ccf%prod100c_loss(ibeg:iend)  = nan
    ccf%product_closs(ibeg:iend)  = nan
#endif
    ccf%lithr(ibeg:iend)          = nan
    ccf%somhr(ibeg:iend)          = nan
    ccf%hr(ibeg:iend)             = nan
    ccf%sr(ibeg:iend)             = nan
    ccf%er(ibeg:iend)             = nan
    ccf%litfire(ibeg:iend)        = nan
    ccf%somfire(ibeg:iend)        = nan
    ccf%totfire(ibeg:iend)        = nan
    ccf%nep(ibeg:iend)            = nan
    ccf%nbp(ibeg:iend)            = nan
    ccf%nee(ibeg:iend)            = nan
    ccf%col_cinputs(ibeg:iend)    = nan
    ccf%col_coutputs(ibeg:iend)   = nan
    ccf%col_fire_closs(ibeg:iend) = nan

#if (defined CN)
    ccf%cwdc_hr(ibeg:iend)        = nan
    ccf%cwdc_loss(ibeg:iend)      = nan
    ccf%litterc_loss(ibeg:iend)   = nan
#endif
  end subroutine init_column_cflux_type
#ifdef LCH4
  !
  ! Initialize column methane flux variables
  !
  subroutine init_column_ch4_type(ibeg,iend,cch4)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_ch4_type) , intent(inout) :: cch4

    allocate(cch4%ch4_prod_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_lake(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_lake(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_surf_aere_sat(ibeg:iend))
    allocate(cch4%ch4_surf_aere_unsat(ibeg:iend))
    allocate(cch4%ch4_ebul_depth_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_ebul_depth_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_ebul_total_sat(ibeg:iend))
    allocate(cch4%ch4_ebul_total_unsat(ibeg:iend))
    allocate(cch4%ch4_surf_ebul_sat(ibeg:iend))
    allocate(cch4%ch4_surf_ebul_unsat(ibeg:iend))
    allocate(cch4%ch4_surf_ebul_lake(ibeg:iend))
    allocate(cch4%conc_ch4_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%conc_ch4_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%conc_ch4_lake(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_surf_diff_sat(ibeg:iend))
    allocate(cch4%ch4_surf_diff_unsat(ibeg:iend))
    allocate(cch4%ch4_surf_diff_lake(ibeg:iend))
    allocate(cch4%conc_o2_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%conc_o2_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%conc_o2_lake(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4_dfsat_flux(ibeg:iend))
    allocate(cch4%zwt_ch4_unsat(ibeg:iend))
    allocate(cch4%fsat_bef(ibeg:iend))
    allocate(cch4%lake_soilc(ibeg:iend,1:nlevgrnd))
    allocate(cch4%lake_raw(ibeg:iend))
    allocate(cch4%totcolch4(ibeg:iend))
    allocate(cch4%fphr(ibeg:iend,1:nlevgrnd))
    allocate(cch4%annsum_counter(ibeg:iend))
    allocate(cch4%tempavg_somhr(ibeg:iend))
    allocate(cch4%annavg_somhr(ibeg:iend))
    allocate(cch4%tempavg_finrw(ibeg:iend))
    allocate(cch4%annavg_finrw(ibeg:iend))
    allocate(cch4%sif(ibeg:iend))
    allocate(cch4%o2stress_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%o2stress_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4stress_unsat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%ch4stress_sat(ibeg:iend,1:nlevgrnd))
    allocate(cch4%qflx_surf_lag(ibeg:iend))
    allocate(cch4%finundated_lag(ibeg:iend))
    allocate(cch4%layer_sat_lag(ibeg:iend,1:nlevgrnd))

    cch4%ch4_prod_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_lake(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_lake(ibeg:iend,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%o2_decomp_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    ! To detect first time-step for denitrification code
    cch4%o2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd) = spval
    cch4%o2_aere_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%o2_aere_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_aere_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%co2_aere_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_surf_aere_sat(ibeg:iend) = nan
    cch4%ch4_surf_aere_unsat(ibeg:iend) = nan
    cch4%ch4_ebul_depth_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_ebul_depth_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_ebul_total_sat(ibeg:iend) = nan
    cch4%ch4_ebul_total_unsat(ibeg:iend) = nan
    cch4%ch4_surf_ebul_sat(ibeg:iend) = nan
    cch4%ch4_surf_ebul_unsat(ibeg:iend) = nan
    cch4%ch4_surf_ebul_lake(ibeg:iend) = nan
    cch4%conc_ch4_sat(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_ch4_unsat(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    ! Just a diagnostic, so nan is fine
    cch4%conc_ch4_lake(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_surf_diff_sat(ibeg:iend) = nan
    cch4%ch4_surf_diff_unsat(ibeg:iend) = nan
    cch4%ch4_surf_diff_lake(ibeg:iend) = nan
    cch4%conc_o2_sat(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    ! To detect file input and detect first time-step for denitrification code
    cch4%conc_o2_unsat(ibeg:iend,1:nlevgrnd) = spval
    ! Just a diagnostic, so nan is fine
    cch4%conc_o2_lake(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4_dfsat_flux(ibeg:iend) = nan
    cch4%zwt_ch4_unsat(ibeg:iend) = nan
    cch4%fsat_bef(ibeg:iend) = spval ! To detect first time-step
    cch4%lake_soilc(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    cch4%lake_raw(ibeg:iend) = nan
    cch4%totcolch4(ibeg:iend) = spval ! To detect first time-step
    cch4%fphr(ibeg:iend,1:nlevgrnd) = nan
    cch4%annsum_counter(ibeg:iend) = spval ! To detect first time-step
    cch4%tempavg_somhr(ibeg:iend) = nan
    cch4%annavg_somhr(ibeg:iend) = spval ! To detect first year
    cch4%tempavg_finrw(ibeg:iend) = nan
    cch4%annavg_finrw(ibeg:iend) = spval ! To detect first year
    cch4%sif(ibeg:iend) = nan
    cch4%o2stress_unsat(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    cch4%o2stress_sat(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    cch4%ch4stress_unsat(ibeg:iend,1:nlevgrnd) = nan
    cch4%ch4stress_sat(ibeg:iend,1:nlevgrnd) = nan
    cch4%qflx_surf_lag(ibeg:iend) = spval ! To detect file input
    cch4%finundated_lag(ibeg:iend) = spval ! To detect file input
    cch4%layer_sat_lag(ibeg:iend,1:nlevgrnd) = spval ! To detect file input
    ! def CH4
  end subroutine init_column_ch4_type
#endif
  !
  ! Initialize column nitrogen flux variables
  !
  subroutine init_column_nflux_type(ibeg,iend,cnf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (column_nflux_type) , intent(inout) :: cnf

    allocate(cnf%ndep_to_sminn(ibeg:iend))
    allocate(cnf%nfix_to_sminn(ibeg:iend))
    allocate(cnf%fert_to_sminn(ibeg:iend))
    allocate(cnf%soyfixn_to_sminn(ibeg:iend))    
    allocate(cnf%hrv_deadstemn_to_prod10n(ibeg:iend))        
    allocate(cnf%hrv_deadstemn_to_prod100n(ibeg:iend))       

    allocate(cnf%m_n_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant(ibeg:iend))
    allocate(cnf%potential_immob(ibeg:iend))
    allocate(cnf%actual_immob(ibeg:iend))
    allocate(cnf%gross_nmin(ibeg:iend))
    allocate(cnf%net_nmin(ibeg:iend))
    allocate(cnf%denit(ibeg:iend))
    allocate(cnf%supplement_to_sminn(ibeg:iend))
    allocate(cnf%m_decomp_npools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cnf%m_decomp_npools_to_fire(ibeg:iend,1:ndecomp_pools))
    allocate(cnf%decomp_cascade_ntransfer_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_ntransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux(ibeg:iend, &
            1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_npools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))

    allocate(cnf%phenology_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%fire_mortality_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full))

#ifndef NITRIF_DENITRIF
    allocate(cnf%sminn_to_denit_decomp_cascade_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_decomp_cascade(ibeg:iend, &
            1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_excess_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%sminn_to_denit_excess(ibeg:iend))
    allocate(cnf%sminn_leached_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%sminn_leached(ibeg:iend))
#else
    allocate(cnf%f_nit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%f_denit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached(ibeg:iend))
    allocate(cnf%smin_no3_runoff_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%smin_no3_runoff(ibeg:iend))
    allocate(cnf%pot_f_nit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%pot_f_nit(ibeg:iend))
    allocate(cnf%pot_f_denit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%pot_f_denit(ibeg:iend))
    allocate(cnf%actual_immob_no3_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%actual_immob_nh4_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%smin_no3_to_plant_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%smin_nh4_to_plant_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%f_nit(ibeg:iend))
    allocate(cnf%f_denit(ibeg:iend))
    allocate(cnf%n2_n2o_ratio_denit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%f_n2o_denit(ibeg:iend))
    allocate(cnf%f_n2o_denit_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%f_n2o_nit(ibeg:iend))
    allocate(cnf%f_n2o_nit_vr(ibeg:iend,1:nlevdecomp_full))

    allocate(cnf%smin_no3_massdens_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%soil_bulkdensity(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%k_nitr_t_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%k_nitr_ph_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%k_nitr_h2o_vr(ibeg:iend,1:nlevdecomp_full)) 
    allocate(cnf%k_nitr_vr(ibeg:iend,1:nlevdecomp_full)) 
    allocate(cnf%wfps_vr(ibeg:iend,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_carbonsubstrate_vr(ibeg:iend,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_nitrate_vr(ibeg:iend,1:nlevdecomp_full)) 
    allocate(cnf%f_denit_base_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%diffus(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%ratio_k1(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%ratio_no3_co2(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%soil_co2_prod(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%fr_WFPS(ibeg:iend,1:nlevdecomp_full))

    allocate(cnf%r_psi(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%anaerobic_frac(ibeg:iend,1:nlevdecomp_full))
#endif
    allocate(cnf%potential_immob_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%actual_immob_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%supplement_to_sminn_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%gross_nmin_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%net_nmin_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_seedn_to_leaf(ibeg:iend))
    allocate(cnf%dwt_seedn_to_deadstem(ibeg:iend))
    allocate(cnf%dwt_conv_nflux(ibeg:iend))
    allocate(cnf%dwt_prod10n_gain(ibeg:iend))
    allocate(cnf%dwt_prod100n_gain(ibeg:iend))
    allocate(cnf%dwt_frootn_to_litr_met_n(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_cel_n(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_lig_n(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_livecrootn_to_cwdn(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_deadcrootn_to_cwdn(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%dwt_nloss(ibeg:iend))
    allocate(cnf%prod10n_loss(ibeg:iend))
    allocate(cnf%prod100n_loss(ibeg:iend))
    allocate(cnf%product_nloss(ibeg:iend))
    allocate(cnf%col_ninputs(ibeg:iend))
    allocate(cnf%col_noutputs(ibeg:iend))
    allocate(cnf%col_fire_nloss(ibeg:iend))
    allocate(cnf%som_n_leached(ibeg:iend))
    allocate(cnf%decomp_npools_leached(ibeg:iend,1:ndecomp_pools))
    allocate(cnf%decomp_npools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools))
    
    cnf%ndep_to_sminn(ibeg:iend) = nan
    cnf%nfix_to_sminn(ibeg:iend) = nan
    cnf%fert_to_sminn(ibeg:iend) = nan
    cnf%soyfixn_to_sminn(ibeg:iend) = nan
    cnf%hrv_deadstemn_to_prod10n(ibeg:iend) = nan        
    cnf%hrv_deadstemn_to_prod100n(ibeg:iend) = nan       
    cnf%m_n_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant(ibeg:iend) = nan
    cnf%potential_immob(ibeg:iend) = nan
    cnf%actual_immob(ibeg:iend) = nan
    cnf%gross_nmin(ibeg:iend) = nan
    cnf%net_nmin(ibeg:iend) = nan
    cnf%denit(ibeg:iend) = nan
    cnf%supplement_to_sminn(ibeg:iend) = nan
    cnf%m_decomp_npools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools)                 = nan
    cnf%m_decomp_npools_to_fire(ibeg:iend,1:ndecomp_pools)     = nan
    cnf%decomp_cascade_ntransfer_vr(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions)  = nan
    cnf%decomp_cascade_ntransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions)  = nan
    cnf%decomp_cascade_sminn_flux_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%decomp_cascade_sminn_flux(ibeg:iend, &
            1:ndecomp_cascade_transitions) = nan
    cnf%decomp_npools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools) = nan
    
    cnf%phenology_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%phenology_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%phenology_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%gap_mortality_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%gap_mortality_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%gap_mortality_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%gap_mortality_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full)  = nan
    cnf%fire_mortality_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full) = nan
    cnf%harvest_n_to_litr_met_n(ibeg:iend, 1:nlevdecomp_full)  = nan
    cnf%harvest_n_to_litr_cel_n(ibeg:iend, 1:nlevdecomp_full)  = nan
    cnf%harvest_n_to_litr_lig_n(ibeg:iend, 1:nlevdecomp_full)  = nan
    cnf%harvest_n_to_cwdn(ibeg:iend, 1:nlevdecomp_full)        = nan

#ifndef NITRIF_DENITRIF
    cnf%sminn_to_denit_decomp_cascade_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_decomp_cascade(ibeg:iend, &
            1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_excess_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%sminn_to_denit_excess(ibeg:iend) = nan
    cnf%sminn_leached_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%sminn_leached(ibeg:iend) = nan
#else
    cnf%f_nit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%f_denit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached(ibeg:iend) = nan
    cnf%smin_no3_runoff_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%smin_no3_runoff(ibeg:iend) = nan
    cnf%pot_f_nit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%pot_f_nit(ibeg:iend) = nan
    cnf%pot_f_denit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%pot_f_denit(ibeg:iend) = nan
    cnf%actual_immob_no3_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%actual_immob_nh4_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%smin_no3_to_plant_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%smin_nh4_to_plant_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%f_nit(ibeg:iend) = nan
    cnf%f_denit(ibeg:iend) = nan
    cnf%n2_n2o_ratio_denit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%f_n2o_denit(ibeg:iend) = nan
    cnf%f_n2o_denit_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%f_n2o_nit(ibeg:iend) = nan
    cnf%f_n2o_nit_vr(ibeg:iend,1:nlevdecomp_full) = nan

    cnf%smin_no3_massdens_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%soil_bulkdensity(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%k_nitr_t_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%k_nitr_ph_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%k_nitr_h2o_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%k_nitr_vr(ibeg:iend,1:nlevdecomp_full) = nan 
    cnf%wfps_vr(ibeg:iend,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_carbonsubstrate_vr(ibeg:iend,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_nitrate_vr(ibeg:iend,1:nlevdecomp_full) = nan 
    cnf%f_denit_base_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%diffus(ibeg:iend,1:nlevdecomp_full) = spval
    cnf%ratio_k1(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%ratio_no3_co2(ibeg:iend,1:nlevdecomp_full) = spval
    cnf%soil_co2_prod(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%fr_WFPS(ibeg:iend,1:nlevdecomp_full) = spval

    cnf%r_psi(ibeg:iend,1:nlevdecomp_full) = spval
    cnf%anaerobic_frac(ibeg:iend,1:nlevdecomp_full) = spval
#endif
    cnf%potential_immob_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%actual_immob_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%supplement_to_sminn_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%gross_nmin_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%net_nmin_vr(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_seedn_to_leaf(ibeg:iend) = nan
    cnf%dwt_seedn_to_deadstem(ibeg:iend) = nan
    cnf%dwt_conv_nflux(ibeg:iend) = nan
    cnf%dwt_prod10n_gain(ibeg:iend) = nan
    cnf%dwt_prod100n_gain(ibeg:iend) = nan
    cnf%dwt_frootn_to_litr_met_n(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_cel_n(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_lig_n(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_livecrootn_to_cwdn(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_deadcrootn_to_cwdn(ibeg:iend,1:nlevdecomp_full) = nan
    cnf%dwt_nloss(ibeg:iend) = nan
    cnf%prod10n_loss(ibeg:iend) = nan
    cnf%prod100n_loss(ibeg:iend) = nan
    cnf%product_nloss(ibeg:iend) = nan
    cnf%col_ninputs(ibeg:iend) = nan
    cnf%col_noutputs(ibeg:iend) = nan
    cnf%col_fire_nloss(ibeg:iend) = nan
    cnf%som_n_leached(ibeg:iend)  = nan 
    cnf%decomp_npools_leached(ibeg:iend,1:ndecomp_pools) = nan
    cnf%decomp_npools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools) = nan
  end subroutine init_column_nflux_type
  !
  ! Initialize landunit physical state variables
  !
  subroutine init_landunit_pstate_type(ibeg,iend,lps)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (landunit_pstate_type) , intent(inout) :: lps

    allocate(lps%t_building(ibeg:iend))
    allocate(lps%t_building_max(ibeg:iend))
    allocate(lps%t_building_min(ibeg:iend))
    if ( nlevurb > 0 )then
      allocate(lps%tk_wall(ibeg:iend,nlevurb))
      allocate(lps%tk_roof(ibeg:iend,nlevurb))
      allocate(lps%cv_wall(ibeg:iend,nlevurb))
      allocate(lps%cv_roof(ibeg:iend,nlevurb))
    end if
    allocate(lps%tk_improad(ibeg:iend,nlevurb))
    allocate(lps%cv_improad(ibeg:iend,nlevurb))
    allocate(lps%thick_wall(ibeg:iend))
    allocate(lps%thick_roof(ibeg:iend))
    allocate(lps%nlev_improad(ibeg:iend))
    allocate(lps%vf_sr(ibeg:iend))
    allocate(lps%vf_wr(ibeg:iend))
    allocate(lps%vf_sw(ibeg:iend))
    allocate(lps%vf_rw(ibeg:iend))
    allocate(lps%vf_ww(ibeg:iend))
    allocate(lps%taf(ibeg:iend))
    allocate(lps%qaf(ibeg:iend))
    allocate(lps%sabs_roof_dir(ibeg:iend,1:numrad))
    allocate(lps%sabs_roof_dif(ibeg:iend,1:numrad))
    allocate(lps%sabs_sunwall_dir(ibeg:iend,1:numrad))
    allocate(lps%sabs_sunwall_dif(ibeg:iend,1:numrad))
    allocate(lps%sabs_shadewall_dir(ibeg:iend,1:numrad))
    allocate(lps%sabs_shadewall_dif(ibeg:iend,1:numrad))
    allocate(lps%sabs_improad_dir(ibeg:iend,1:numrad))
    allocate(lps%sabs_improad_dif(ibeg:iend,1:numrad))
    allocate(lps%sabs_perroad_dir(ibeg:iend,1:numrad))
    allocate(lps%sabs_perroad_dif(ibeg:iend,1:numrad))

    lps%t_building(ibeg:iend) = nan
    lps%t_building_max(ibeg:iend) = nan
    lps%t_building_min(ibeg:iend) = nan
    if ( nlevurb > 0 )then
      lps%tk_wall(ibeg:iend,1:nlevurb) = nan
      lps%tk_roof(ibeg:iend,1:nlevurb) = nan
      lps%cv_wall(ibeg:iend,1:nlevurb) = nan
      lps%cv_roof(ibeg:iend,1:nlevurb) = nan
    end if
    lps%tk_improad(ibeg:iend,1:nlevurb) = nan
    lps%cv_improad(ibeg:iend,1:nlevurb) = nan
    lps%thick_wall(ibeg:iend) = nan
    lps%thick_roof(ibeg:iend) = nan
    lps%nlev_improad(ibeg:iend) = bigint
    lps%vf_sr(ibeg:iend) = nan
    lps%vf_wr(ibeg:iend) = nan
    lps%vf_sw(ibeg:iend) = nan
    lps%vf_rw(ibeg:iend) = nan
    lps%vf_ww(ibeg:iend) = nan
    lps%taf(ibeg:iend) = nan
    lps%qaf(ibeg:iend) = nan
    lps%sabs_roof_dir(ibeg:iend,1:numrad) = nan
    lps%sabs_roof_dif(ibeg:iend,1:numrad) = nan
    lps%sabs_sunwall_dir(ibeg:iend,1:numrad) = nan
    lps%sabs_sunwall_dif(ibeg:iend,1:numrad) = nan
    lps%sabs_shadewall_dir(ibeg:iend,1:numrad) = nan
    lps%sabs_shadewall_dif(ibeg:iend,1:numrad) = nan
    lps%sabs_improad_dir(ibeg:iend,1:numrad) = nan
    lps%sabs_improad_dif(ibeg:iend,1:numrad) = nan
    lps%sabs_perroad_dir(ibeg:iend,1:numrad) = nan
    lps%sabs_perroad_dif(ibeg:iend,1:numrad) = nan
  end subroutine init_landunit_pstate_type
  !
  ! Initialize landunit energy flux variables
  !
  subroutine init_landunit_eflux_type(ibeg,iend,lef)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend 
    type (landunit_eflux_type) , intent(inout) :: lef 

    allocate(lef%eflx_traffic(ibeg:iend))
    allocate(lef%eflx_traffic_factor(ibeg:iend))
    allocate(lef%eflx_wasteheat(ibeg:iend))
    allocate(lef%eflx_heat_from_ac(ibeg:iend))

    lef%eflx_traffic(ibeg:iend) = nan
    lef%eflx_traffic_factor(ibeg:iend) = nan
    lef%eflx_wasteheat(ibeg:iend) = nan
    lef%eflx_heat_from_ac(ibeg:iend) = nan
  end subroutine init_landunit_eflux_type

#if (defined CNDV)
  !
  ! Initialize gridcell DGVM variables
  !
  subroutine init_gridcell_dgvstate_type(ibeg,iend,gps)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_dgvstate_type) , intent(inout) :: gps

    allocate(gps%agdd20(ibeg:iend))
    allocate(gps%tmomin20(ibeg:iend))
    allocate(gps%t10min(ibeg:iend))
    gps%agdd20(ibeg:iend) = nan
    gps%tmomin20(ibeg:iend) = nan
    gps%t10min(ibeg:iend) = nan
  end subroutine init_gridcell_dgvstate_type
#endif
  !
  ! Initialize gridcell isoprene emission factor variables
  !
  subroutine init_gridcell_efstate_type(ibeg,iend,gve)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_efstate_type) , intent(inout) :: gve

    allocate(gve%efisop(6,ibeg:iend))
    gve%efisop(:,ibeg:iend) = nan
  end subroutine init_gridcell_efstate_type
  !
  ! Initialize gridcell water flux variables
  !
  subroutine init_gridcell_wflux_type(ibeg,iend,gwf)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_wflux_type) , intent(inout) :: gwf

    allocate(gwf%qflx_runoffg(ibeg:iend))
    allocate(gwf%qflx_snwcp_iceg(ibeg:iend))
    allocate(gwf%qflx_liq_dynbal(ibeg:iend))
    allocate(gwf%qflx_ice_dynbal(ibeg:iend))

    gwf%qflx_runoffg(ibeg:iend) = 0.D0
    gwf%qflx_snwcp_iceg(ibeg:iend) = 0.D0
    gwf%qflx_liq_dynbal(ibeg:iend) = nan
    gwf%qflx_ice_dynbal(ibeg:iend) = nan
  end subroutine init_gridcell_wflux_type
  !
  ! Initialize gridcell energy flux variables
  !
  subroutine init_gridcell_eflux_type(ibeg,iend,gef)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_eflux_type) , intent(inout) :: gef

    allocate(gef%eflx_sh_totg(ibeg:iend))
    allocate(gef%eflx_dynbal(ibeg:iend))

    gef%eflx_sh_totg(ibeg:iend) = nan
    gef%eflx_dynbal(ibeg:iend) = nan
  end subroutine init_gridcell_eflux_type
  !
  ! Initialize gridcell water state variables
  !
  subroutine init_gridcell_wstate_type(ibeg,iend,gws)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_wstate_type) , intent(inout) :: gws

    allocate(gws%gc_liq1(ibeg:iend))
    allocate(gws%gc_liq2(ibeg:iend))
    allocate(gws%gc_ice1(ibeg:iend))     
    allocate(gws%gc_ice2(ibeg:iend))    

    gws%gc_liq1(ibeg:iend) = nan
    gws%gc_liq2(ibeg:iend) = nan
    gws%gc_ice1(ibeg:iend) = nan     
    gws%gc_ice2(ibeg:iend) = nan    
  end subroutine init_gridcell_wstate_type
  !
  ! Initialize gridcell energy state variables     
  !
  subroutine init_gridcell_estate_type(ibeg,iend,ges)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_estate_type) , intent(inout) :: ges

    allocate(ges%gc_heat1(ibeg:iend))     
    allocate(ges%gc_heat2(ibeg:iend))    

    ges%gc_heat1(ibeg:iend) = nan     
    ges%gc_heat2(ibeg:iend) = nan    
  end subroutine init_gridcell_estate_type
#ifdef LCH4
  !
  ! Initialize gridcell ch4 variables
  !
  subroutine init_gridcell_ch4_type(ibeg,iend,gch4)
    implicit none
    integer(ik4) , intent(in) :: ibeg , iend
    type (gridcell_ch4_type) , intent(inout) :: gch4

    allocate(gch4%c_atm(ibeg:iend,1:ngases))
    allocate(gch4%ch4co2f(ibeg:iend))
    allocate(gch4%ch4prodg(ibeg:iend))
    allocate(gch4%nem(ibeg:iend))

    gch4%c_atm(ibeg:iend,1:ngases) = nan
    gch4%ch4co2f(ibeg:iend) = nan
    gch4%ch4prodg(ibeg:iend) = nan
    gch4%nem(ibeg:iend) = nan
  end subroutine init_gridcell_ch4_type
#endif

  end subroutine initClmtype

end module mod_clm_typeinit
