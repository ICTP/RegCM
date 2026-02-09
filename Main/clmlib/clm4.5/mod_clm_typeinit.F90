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
  use mod_clm_decomp, only : get_proc_bounds, get_proc_global
  use mod_clm_varcon, only : spval, ispval
  use mod_clm_surfrd, only : crop_prog
  use mod_clm_megan, only : shr_megan_megcomps_n
  use mod_clm_drydep, only : n_drydep, drydep_method, DD_XLND
  use mod_clm_varpar

  implicit none

  private

  save

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
    integer(ik4) :: begp, endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numg        ! total number of gridcells across all proc
    integer(ik4) :: numl        ! total number of landunits across all proc
    integer(ik4) :: numc        ! total number of columns across all proc
    integer(ik4) :: nump        ! total number of pfts across all proc

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
       "initClmtype ERROR:: CROP and C13 can NOT be on at the same time")
#endif
    end if

    if ( use_c14 ) then
      call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc14s)
      call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc14s%pcs_a)
#ifdef CROP
      call fatal(__FILE__,__LINE__,&
       "initClmtype ERROR:: CROP and C14 can NOT be on at the same time")
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
    integer(ik4), intent(in) :: ibeg, iend
    type(pft_type), intent(inout) :: p

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
    integer(ik4), intent(in) :: ibeg, iend
    type(column_type), intent(inout) :: c

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
    integer(ik4), intent(in) :: ibeg, iend
    type(landunit_type), intent(inout) :: l

    allocate(l%gridcell(ibeg:iend),l%wtgcell(ibeg:iend))

    allocate(l%coli(ibeg:iend),l%colf(ibeg:iend),l%ncolumns(ibeg:iend))
    allocate(l%pfti(ibeg:iend),l%pftf(ibeg:iend),l%npfts   (ibeg:iend))

    allocate(l%itype(ibeg:iend))
    allocate(l%ifspecial(ibeg:iend))
    allocate(l%lakpoi(ibeg:iend))
    allocate(l%urbpoi(ibeg:iend))
    allocate(l%udenstype(ibeg:iend))
    allocate(l%active(ibeg:iend))

    ! MV - these should be moved to landunit physical state -MV
    allocate(l%canyon_hwr(ibeg:iend), source=nan)
    allocate(l%wtroad_perv(ibeg:iend), source=nan)
    allocate(l%ht_roof(ibeg:iend), source=nan)
    allocate(l%wtlunit_roof(ibeg:iend), source=nan)
    allocate(l%wind_hgt_canyon(ibeg:iend), source=nan)
    allocate(l%z_0_town(ibeg:iend), source=nan)
    allocate(l%z_d_town(ibeg:iend), source=nan)
  end subroutine init_landunit_type
  !
  ! Initialize components of gridcell_type structure
  !
  subroutine init_gridcell_type (ibeg,iend,g)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type(gridcell_type), intent(inout) :: g

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

    allocate(g%gris_mask(ibeg:iend), source=nan)
    allocate(g%gris_area(ibeg:iend), source=nan)
    allocate(g%aais_mask(ibeg:iend), source=nan)
    allocate(g%aais_area(ibeg:iend), source=nan)
    allocate(g%tws(ibeg:iend), source=nan)
  end subroutine init_gridcell_type
  !
  ! Initialize energy balance variables
  !
  subroutine init_energy_balance_type(ibeg,iend,ebal)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type(energy_balance_type), intent(inout) :: ebal

    allocate(ebal%errsoi(ibeg:iend), source=nan)
    allocate(ebal%errseb(ibeg:iend), source=nan)
    allocate(ebal%errsol(ibeg:iend), source=nan)
    allocate(ebal%errlon(ibeg:iend), source=nan)
  end subroutine init_energy_balance_type
  !
  ! Initialize water balance variables
  !
  subroutine init_water_balance_type(ibeg,iend,wbal)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type(water_balance_type), intent(inout) :: wbal

    allocate(wbal%begwb(ibeg:iend), source=nan)
    allocate(wbal%endwb(ibeg:iend), source=nan)
    allocate(wbal%errh2o(ibeg:iend), source=nan)
  end subroutine init_water_balance_type
  !
  ! Initialize carbon balance variables
  !
  subroutine init_carbon_balance_type(ibeg,iend,cbal)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type(carbon_balance_type), intent(inout) :: cbal

    allocate(cbal%begcb(ibeg:iend), source=nan)
    allocate(cbal%endcb(ibeg:iend), source=nan)
    allocate(cbal%errcb(ibeg:iend), source=nan)
  end subroutine init_carbon_balance_type
  !
  ! Initialize nitrogen balance variables
  !
  subroutine init_nitrogen_balance_type(ibeg,iend,nbal)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type(nitrogen_balance_type), intent(inout) :: nbal

    allocate(nbal%begnb(ibeg:iend), source=nan)
    allocate(nbal%endnb(ibeg:iend), source=nan)
    allocate(nbal%errnb(ibeg:iend), source=nan)
  end subroutine init_nitrogen_balance_type
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_ecophys_constants()
    implicit none

    allocate(pftcon%noveg(0:numpft), source=bigint)
    allocate(pftcon%tree(0:numpft), source=bigint)
    allocate(pftcon%smpso(0:numpft), source=nan)
    allocate(pftcon%smpsc(0:numpft), source=nan)
    allocate(pftcon%fnitr(0:numpft), source=nan)
    allocate(pftcon%foln(0:numpft), source=nan)
    allocate(pftcon%dleaf(0:numpft), source=nan)
    allocate(pftcon%c3psn(0:numpft), source=nan)
    allocate(pftcon%xl(0:numpft), source=nan)
    allocate(pftcon%rhol(0:numpft,numrad), source=nan)
    allocate(pftcon%rhos(0:numpft,numrad), source=nan)
    allocate(pftcon%taul(0:numpft,numrad), source=nan)
    allocate(pftcon%taus(0:numpft,numrad), source=nan)
    allocate(pftcon%z0mr(0:numpft), source=nan)
    allocate(pftcon%displar(0:numpft), source=nan)
    allocate(pftcon%roota_par(0:numpft), source=nan)
    allocate(pftcon%rootb_par(0:numpft), source=nan)
    allocate(pftcon%slatop(0:numpft), source=nan)
    allocate(pftcon%dsladlai(0:numpft), source=nan)
    allocate(pftcon%leafcn(0:numpft), source=nan)
    allocate(pftcon%flnr(0:numpft), source=nan)
    allocate(pftcon%woody(0:numpft), source=nan)
    allocate(pftcon%lflitcn(0:numpft), source=nan)
    allocate(pftcon%frootcn(0:numpft), source=nan)
    allocate(pftcon%livewdcn(0:numpft), source=nan)
    allocate(pftcon%deadwdcn(0:numpft), source=nan)
    allocate(pftcon%graincn(0:numpft), source=nan)
    allocate(pftcon%froot_leaf(0:numpft), source=nan)
    allocate(pftcon%stem_leaf(0:numpft), source=nan)
    allocate(pftcon%croot_stem(0:numpft), source=nan)
    allocate(pftcon%flivewd(0:numpft), source=nan)
    allocate(pftcon%fcur(0:numpft), source=nan)
    allocate(pftcon%lf_flab(0:numpft), source=nan)
    allocate(pftcon%lf_fcel(0:numpft), source=nan)
    allocate(pftcon%lf_flig(0:numpft), source=nan)
    allocate(pftcon%fr_flab(0:numpft), source=nan)
    allocate(pftcon%fr_fcel(0:numpft), source=nan)
    allocate(pftcon%fr_flig(0:numpft), source=nan)
    allocate(pftcon%leaf_long(0:numpft), source=nan)
    allocate(pftcon%evergreen(0:numpft), source=nan)
    allocate(pftcon%stress_decid(0:numpft), source=nan)
    allocate(pftcon%season_decid(0:numpft), source=nan)
    allocate(pftcon%dwood(0:numpft), source=nan)
    allocate(pftcon%rootprof_beta(0:numpft), source=nan)
    allocate(pftcon%fertnitro(0:numpft), source=nan)
    allocate(pftcon%fleafcn(0:numpft), source=nan)
    allocate(pftcon%ffrootcn(0:numpft), source=nan)
    allocate(pftcon%fstemcn(0:numpft), source=nan)
  end subroutine init_pft_ecophys_constants
  !
  ! Initialize decomposition cascade state
  !
  subroutine init_decomp_cascade_constants()
    implicit none
    integer(ik4) :: nct, np
    nct = ndecomp_cascade_transitions
    np = ndecomp_pools
    !-- properties of each pathway along decomposition cascade
    allocate(decomp_cascade_con%cascade_step_name(1:nct))
    allocate(decomp_cascade_con%cascade_donor_pool(1:nct), source=0)
    allocate(decomp_cascade_con%cascade_receiver_pool(1:nct), source=0)
    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:np), &
            source=.false.)
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:np))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:np))
    allocate(decomp_cascade_con%is_litter(0:np), source=.false.)
    allocate(decomp_cascade_con%is_soil(0:np), source=.false.)
    allocate(decomp_cascade_con%is_cwd(0:np), source=.false.)
    allocate(decomp_cascade_con%initial_cn_ratio(0:np), source=nan)
    allocate(decomp_cascade_con%initial_stock(0:np), source=nan)
    allocate(decomp_cascade_con%is_metabolic(0:np), source=.false.)
    allocate(decomp_cascade_con%is_cellulose(0:np), source=.false.)
    allocate(decomp_cascade_con%is_lignin(0:np), source=.false.)
    allocate(decomp_cascade_con%spinup_factor(0:np), source=nan)
    !-- properties of each pathway along decomposition cascade
    decomp_cascade_con%cascade_step_name(1:nct) = ''
    !-- properties of each decomposing pool
    decomp_cascade_con%decomp_pool_name_history(0:np) = ''
    decomp_cascade_con%decomp_pool_name_restart(0:np) = ''
    decomp_cascade_con%decomp_pool_name_long(0:np) = ''
    decomp_cascade_con%decomp_pool_name_short(0:np) = ''
  end subroutine init_decomp_cascade_constants

#if (defined CNDV)
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_DGVMecophys_constants()
    implicit none

    allocate(dgv_pftcon%crownarea_max(0:numpft), source=nan)
    allocate(dgv_pftcon%tcmin(0:numpft), source=nan)
    allocate(dgv_pftcon%tcmax(0:numpft), source=nan)
    allocate(dgv_pftcon%gddmin(0:numpft), source=nan)
    allocate(dgv_pftcon%twmax(0:numpft), source=nan)
    allocate(dgv_pftcon%reinickerp(0:numpft), source=nan)
    allocate(dgv_pftcon%allom1(0:numpft), source=nan)
    allocate(dgv_pftcon%allom2(0:numpft), source=nan)
    allocate(dgv_pftcon%allom3(0:numpft), source=nan)
  end subroutine init_pft_DGVMecophys_constants

#endif
  !
  ! Initialize pft physical state
  !
  subroutine init_pft_pstate_type(ibeg,iend,pps)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_pstate_type), intent(inout) :: pps

    allocate(pps%prec10(ibeg:iend), source=nan) !F. Li and S. Levis
    allocate(pps%prec60(ibeg:iend), source=nan) !F. Li and S. Levis
    allocate(pps%frac_veg_nosno(ibeg:iend), source=bigint)
    allocate(pps%frac_veg_nosno_alb(ibeg:iend), source=0)
    allocate(pps%emv(ibeg:iend), source=nan)
    allocate(pps%z0mv(ibeg:iend), source=0.0_rk8)
    allocate(pps%z0hv(ibeg:iend), source=0.0_rk8)
    allocate(pps%z0qv(ibeg:iend), source=0.0_rk8)
    allocate(pps%rootfr(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(pps%rootr(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(pps%rresis(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(pps%dewmx(ibeg:iend), source=nan)
    allocate(pps%rssun(ibeg:iend), source=nan)
    allocate(pps%rssha(ibeg:iend), source=nan)
    allocate(pps%rhal(ibeg:iend), source=nan)
    allocate(pps%vpdal(ibeg:iend), source=nan)
    allocate(pps%rssun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pps%rssha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pps%laisun(ibeg:iend), source=nan)
    allocate(pps%laisha(ibeg:iend), source=nan)
    allocate(pps%laisun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pps%laisha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pps%btran(ibeg:iend), source=spval)
    allocate(pps%btran2(ibeg:iend), source=spval)   ! F. Li and S. Levis
    allocate(pps%fsun(ibeg:iend), source=spval)
    allocate(pps%tlai(ibeg:iend), source=0.0_rk8)
    allocate(pps%tsai(ibeg:iend), source=0.0_rk8)
    allocate(pps%elai(ibeg:iend), source=0.0_rk8)
    allocate(pps%esai(ibeg:iend), source=0.0_rk8)
    allocate(pps%fwet(ibeg:iend), source=nan)
    allocate(pps%fdry(ibeg:iend), source=nan)
    allocate(pps%dt_veg(ibeg:iend), source=nan)
    allocate(pps%htop(ibeg:iend), source=0.0_rk8)
    allocate(pps%hbot(ibeg:iend), source=0.0_rk8)
    allocate(pps%z0m(ibeg:iend), source=nan)
    allocate(pps%displa(ibeg:iend), source=nan)
    allocate(pps%albd(ibeg:iend,1:numrad), source=nan)
    allocate(pps%albi(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabd(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabd_sun(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabd_sha(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabi(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabi_sun(ibeg:iend,1:numrad), source=nan)
    allocate(pps%fabi_sha(ibeg:iend,1:numrad), source=nan)
    allocate(pps%ftdd(ibeg:iend,1:numrad), source=nan)
    allocate(pps%ftid(ibeg:iend,1:numrad), source=nan)
    allocate(pps%ftii(ibeg:iend,1:numrad), source=nan)
    allocate(pps%vcmaxcintsun(ibeg:iend), source=nan)
    allocate(pps%vcmaxcintsha(ibeg:iend), source=nan)
    allocate(pps%ncan(ibeg:iend), source=0)
    allocate(pps%nrad(ibeg:iend), source=0)
    allocate(pps%fabd_sun_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%fabd_sha_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%fabi_sun_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%fabi_sha_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%fsun_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%tlai_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%tsai_z(ibeg:iend,1:nlevcan), source=0.0_rk8)
    allocate(pps%u10(ibeg:iend), source=nan)
    allocate(pps%u10_clm(ibeg:iend), source=nan)
    allocate(pps%va(ibeg:iend), source=nan)
    allocate(pps%fv(ibeg:iend), source=nan)
    allocate(pps%ram1(ibeg:iend), source=nan)
    allocate(pps%rah1(ibeg:iend), source=nan)
    allocate(pps%br1(ibeg:iend), source=nan)
    allocate(pps%burndate(ibeg:iend), source=ispval)   ! F. Li and S. Levis
    allocate(pps%ram1_lake(ibeg:iend), source=nan)
    allocate(pps%rh_leaf(ibeg:iend), source=spval)
    allocate(pps%rhaf(ibeg:iend), source=spval)
    if ( crop_prog ) then
      allocate(pps%hdidx(ibeg:iend), source=nan)
      allocate(pps%cumvd(ibeg:iend), source=nan)
      allocate(pps%htmx(ibeg:iend), source=0.0_rk8)
      allocate(pps%vf(ibeg:iend), source=0.0_rk8)
      allocate(pps%gddmaturity(ibeg:iend), source=spval)
      allocate(pps%gdd0(ibeg:iend), source=spval)
      allocate(pps%gdd8(ibeg:iend), source=spval)
      allocate(pps%gdd10(ibeg:iend), source=spval)
      allocate(pps%gdd020(ibeg:iend), source=spval)
      allocate(pps%gdd820(ibeg:iend), source=spval)
      allocate(pps%gdd1020(ibeg:iend), source=spval)
      allocate(pps%gddplant(ibeg:iend), source=spval)
      allocate(pps%gddtsoi(ibeg:iend), source=spval)
      allocate(pps%huileaf(ibeg:iend), source=nan)
      allocate(pps%huigrain(ibeg:iend), source=nan)
      allocate(pps%aleafi(ibeg:iend), source=nan)
      allocate(pps%astemi(ibeg:iend), source=nan)
      allocate(pps%aleaf(ibeg:iend), source=nan)
      allocate(pps%astem(ibeg:iend), source=nan)
      allocate(pps%croplive(ibeg:iend), source=.false.)
      ! make 2-D if using crop rotation
      allocate(pps%cropplant(ibeg:iend), source=.false.) !,numpft))
      allocate(pps%harvdate(ibeg:iend), source=bigint)  !,numpft))
      allocate(pps%idop(ibeg:iend), source=bigint)
      allocate(pps%peaklai(ibeg:iend), source=0)
    end if
    allocate(pps%vds(ibeg:iend), source=nan)
    allocate(pps%forc_hgt_u_pft(ibeg:iend), source=nan)
    allocate(pps%forc_hgt_t_pft(ibeg:iend), source=nan)
    allocate(pps%forc_hgt_q_pft(ibeg:iend), source=nan)
    allocate(pps%lfpftd(ibeg:iend))      !F. Li and S. Levis

    ! 4/14/05: PET
    ! Adding isotope code

    if ( use_c13 ) then
      allocate(pps%alphapsnsun(ibeg:iend), source=spval)
      allocate(pps%alphapsnsha(ibeg:iend), source=spval)
    endif

    allocate(pps%sandfrac(ibeg:iend), source=nan)
    allocate(pps%clayfrac(ibeg:iend), source=nan)
    allocate(pps%mlaidiff(ibeg:iend), source=nan)
    allocate(pps%rb1(ibeg:iend), source=nan)
    allocate(pps%annlai(12,ibeg:iend), source=nan)

#if (defined LCH4)
    ! CH4 code
    allocate(pps%grnd_ch4_cond(ibeg:iend), source=nan)
    allocate(pps%canopy_cond(ibeg:iend), source=nan)
#endif
   ! and vertical profiles for calculating fluxes
    allocate(pps%leaf_prof(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(pps%froot_prof(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(pps%croot_prof(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(pps%stem_prof(ibeg:iend,1:nlevdecomp_full), source=spval)

    ! 4/14/05: PET
    ! Adding isotope code    ! EBK Check this!
    !!!pps%cisun(ibeg:iend) = spval
    !!!pps%cisha(ibeg:iend) = spval
  end subroutine init_pft_pstate_type
  !
  ! Initialize pft ecophysiological variables
  !
  subroutine init_pft_epv_type(ibeg,iend,pepv)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_epv_type), intent(inout) :: pepv

    allocate(pepv%dormant_flag(ibeg:iend), source=nan)
    allocate(pepv%days_active(ibeg:iend), source=nan)
    allocate(pepv%onset_flag(ibeg:iend), source=nan)
    allocate(pepv%onset_counter(ibeg:iend), source=nan)
    allocate(pepv%onset_gddflag(ibeg:iend), source=nan)
    allocate(pepv%onset_fdd(ibeg:iend), source=nan)
    allocate(pepv%onset_gdd(ibeg:iend), source=nan)
    allocate(pepv%onset_swi(ibeg:iend), source=nan)
    allocate(pepv%offset_flag(ibeg:iend), source=nan)
    allocate(pepv%offset_counter(ibeg:iend), source=nan)
    allocate(pepv%offset_fdd(ibeg:iend), source=nan)
    allocate(pepv%offset_swi(ibeg:iend), source=nan)
    allocate(pepv%fert_counter(ibeg:iend), source=nan)
    allocate(pepv%grain_flag(ibeg:iend), source=nan)
    allocate(pepv%lgsf(ibeg:iend), source=nan)
    allocate(pepv%bglfr(ibeg:iend), source=nan)
    allocate(pepv%bgtr(ibeg:iend), source=nan)
    allocate(pepv%dayl(ibeg:iend), source=nan)
    allocate(pepv%prev_dayl(ibeg:iend), source=nan)
    allocate(pepv%annavg_t2m(ibeg:iend), source=nan)
    allocate(pepv%tempavg_t2m(ibeg:iend), source=nan)
    allocate(pepv%gpp(ibeg:iend), source=nan)
    allocate(pepv%availc(ibeg:iend), source=nan)
    allocate(pepv%xsmrpool_recover(ibeg:iend), source=nan)
    if ( use_c13 ) then
      allocate(pepv%xsmrpool_c13ratio(ibeg:iend), source=nan)
    endif
    allocate(pepv%alloc_pnow(ibeg:iend), source=nan)
    allocate(pepv%c_allometry(ibeg:iend), source=nan)
    allocate(pepv%n_allometry(ibeg:iend), source=nan)
    allocate(pepv%plant_ndemand(ibeg:iend), source=nan)
    allocate(pepv%tempsum_potential_gpp(ibeg:iend), source=nan)
    allocate(pepv%annsum_potential_gpp(ibeg:iend), source=nan)
    allocate(pepv%tempmax_retransn(ibeg:iend), source=nan)
    allocate(pepv%annmax_retransn(ibeg:iend), source=nan)
    allocate(pepv%avail_retransn(ibeg:iend), source=nan)
    allocate(pepv%plant_nalloc(ibeg:iend), source=nan)
    allocate(pepv%plant_calloc(ibeg:iend), source=nan)
    allocate(pepv%excess_cflux(ibeg:iend), source=nan)
    allocate(pepv%downreg(ibeg:iend), source=nan)
    allocate(pepv%prev_leafc_to_litter(ibeg:iend), source=nan)
    allocate(pepv%prev_frootc_to_litter(ibeg:iend), source=nan)
    allocate(pepv%tempsum_npp(ibeg:iend), source=nan)
    allocate(pepv%annsum_npp(ibeg:iend), source=nan)
#if (defined CNDV)
    allocate(pepv%tempsum_litfall(ibeg:iend), source=nan)
    allocate(pepv%annsum_litfall(ibeg:iend), source=nan)
#endif
    if ( use_c13 ) then
      allocate(pepv%rc13_canair(ibeg:iend), source=spval)
      allocate(pepv%rc13_psnsun(ibeg:iend), source=spval)
      allocate(pepv%rc13_psnsha(ibeg:iend), source=spval)
    endif

    if ( use_c14 ) then
      allocate(pepv%rc14_atm(ibeg:iend), source=nan)
      ! allocate(pepv%rc14_canair(ibeg:iend), source=nan)
      ! allocate(pepv%rc14_psnsun(ibeg:iend), source=nan)
      ! allocate(pepv%rc14_psnsha(ibeg:iend), source=nan)
    endif
  end subroutine init_pft_epv_type

#if (defined CNDV)
  !
  ! Initialize pft DGVM state variables
  !
  subroutine init_pft_pdgvstate_type(ibeg,iend,pdgvs)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_dgvstate_type), intent(inout) :: pdgvs

    allocate(pdgvs%agddtw(ibeg:iend), source=nan)
    allocate(pdgvs%agdd(ibeg:iend), source=nan)
    allocate(pdgvs%t_mo(ibeg:iend), source=nan)
    allocate(pdgvs%t_mo_min(ibeg:iend), source=nan)
    allocate(pdgvs%prec365(ibeg:iend), source=nan)
    allocate(pdgvs%present(ibeg:iend), source=.false.)
    allocate(pdgvs%pftmayexist(ibeg:iend), source=.true.)
    allocate(pdgvs%nind(ibeg:iend), source=nan)
    allocate(pdgvs%lm_ind(ibeg:iend), source=nan)
    allocate(pdgvs%lai_ind(ibeg:iend), source=nan)
    allocate(pdgvs%fpcinc(ibeg:iend), source=nan)
    allocate(pdgvs%fpcgrid(ibeg:iend), source=nan)
    allocate(pdgvs%fpcgridold(ibeg:iend), source=nan)
    allocate(pdgvs%crownarea(ibeg:iend), source=nan)
    allocate(pdgvs%greffic(ibeg:iend), source=nan)
    allocate(pdgvs%heatstress(ibeg:iend), source=nan)
  end subroutine init_pft_pdgvstate_type
#endif
  !
  ! Initialize pft VOC variables
  !
  subroutine init_pft_vstate_type(ibeg,iend,pvs)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_vstate_type), intent(inout) :: pvs

    allocate(pvs%t_veg24 (ibeg:iend), source=spval)
    allocate(pvs%t_veg240(ibeg:iend), source=spval)
    allocate(pvs%fsd24   (ibeg:iend), source=spval)
    allocate(pvs%fsd240  (ibeg:iend), source=spval)
    allocate(pvs%fsi24   (ibeg:iend), source=spval)
    allocate(pvs%fsi240  (ibeg:iend), source=spval)
    allocate(pvs%fsun24  (ibeg:iend), source=spval)
    allocate(pvs%fsun240 (ibeg:iend), source=spval)
    allocate(pvs%elai_p  (ibeg:iend), source=spval)
  end subroutine init_pft_vstate_type
  !
  ! Initialize pft energy state
  !
  subroutine init_pft_psynstate_type(ibeg,iend,ppsyns)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_psynstate_type), intent(inout) :: ppsyns

    allocate(ppsyns%c3flag(ibeg:iend), source=.false.)
    allocate(ppsyns%ac(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%aj(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%ap(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%ag(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%an(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%vcmax_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%cp(ibeg:iend), source=nan)
    allocate(ppsyns%kc(ibeg:iend), source=nan)
    allocate(ppsyns%ko(ibeg:iend), source=nan)
    allocate(ppsyns%qe(ibeg:iend), source=nan)
    allocate(ppsyns%tpu_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%kp_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(ppsyns%theta_cj(ibeg:iend), source=nan)
    allocate(ppsyns%bbb(ibeg:iend), source=nan)
    allocate(ppsyns%mbb(ibeg:iend), source=nan)
    allocate(ppsyns%gb_mol(ibeg:iend), source=nan)
    allocate(ppsyns%gs_mol(ibeg:iend,1:nlevcan), source=nan)
  end subroutine init_pft_psynstate_type
  !
  ! Initialize pft energy state
  !
  subroutine init_pft_estate_type(ibeg,iend,pes)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_estate_type), intent(inout) :: pes

    allocate(pes%t_ref2m(ibeg:iend), source=nan)
    allocate(pes%q10m(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min_inst(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max_inst(ibeg:iend), source=nan)
    allocate(pes%q_ref2m(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_u(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_r(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min_u(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min_r(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max_u(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max_r(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min_inst_u(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_min_inst_r(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max_inst_u(ibeg:iend), source=nan)
    allocate(pes%t_ref2m_max_inst_r(ibeg:iend), source=nan)
    allocate(pes%t10(ibeg:iend), source=spval)
    if ( crop_prog )then
      allocate(pes%a10tmin(ibeg:iend), source=spval)
      allocate(pes%a5tmin(ibeg:iend), source=spval)
    end if
    allocate(pes%rh_ref2m(ibeg:iend), source=nan)
    allocate(pes%rh_ref2m_u(ibeg:iend), source=nan)
    allocate(pes%rh_ref2m_r(ibeg:iend), source=nan)
    allocate(pes%t_veg(ibeg:iend), source=nan)
    allocate(pes%thm(ibeg:iend), source=nan)
  end subroutine init_pft_estate_type
  !
  ! Initialize pft water state
  !
  subroutine init_pft_wstate_type(ibeg,iend,pws)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_wstate_type), intent(inout) :: pws !pft water state

    allocate(pws%h2ocan(ibeg:iend), source=spval)
  end subroutine init_pft_wstate_type
  !
  ! Initialize pft carbon state
  !
  subroutine init_pft_cstate_type(ibeg,iend,pcs)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_cstate_type), intent(inout) :: pcs !pft carbon state

    allocate(pcs%leafc(ibeg:iend), source=nan)
    allocate(pcs%leafc_storage(ibeg:iend), source=nan)
    allocate(pcs%leafc_xfer(ibeg:iend), source=nan)
    allocate(pcs%frootc(ibeg:iend), source=nan)
    allocate(pcs%frootc_storage(ibeg:iend), source=nan)
    allocate(pcs%frootc_xfer(ibeg:iend), source=nan)
    allocate(pcs%livestemc(ibeg:iend), source=nan)
    allocate(pcs%livestemc_storage(ibeg:iend), source=nan)
    allocate(pcs%livestemc_xfer(ibeg:iend), source=nan)
    allocate(pcs%deadstemc(ibeg:iend), source=nan)
    allocate(pcs%deadstemc_storage(ibeg:iend), source=nan)
    allocate(pcs%deadstemc_xfer(ibeg:iend), source=nan)
    allocate(pcs%livecrootc(ibeg:iend), source=nan)
    allocate(pcs%livecrootc_storage(ibeg:iend), source=nan)
    allocate(pcs%livecrootc_xfer(ibeg:iend), source=nan)
    allocate(pcs%deadcrootc(ibeg:iend), source=nan)
    allocate(pcs%deadcrootc_storage(ibeg:iend), source=nan)
    allocate(pcs%deadcrootc_xfer(ibeg:iend), source=nan)
    allocate(pcs%gresp_storage(ibeg:iend), source=nan)
    allocate(pcs%gresp_xfer(ibeg:iend), source=nan)
    allocate(pcs%cpool(ibeg:iend), source=nan)
    allocate(pcs%xsmrpool(ibeg:iend), source=nan)
    allocate(pcs%pft_ctrunc(ibeg:iend), source=nan)
    allocate(pcs%dispvegc(ibeg:iend), source=nan)
    allocate(pcs%storvegc(ibeg:iend), source=nan)
    allocate(pcs%totvegc(ibeg:iend), source=nan)
    allocate(pcs%totpftc(ibeg:iend), source=nan)
    allocate(pcs%leafcmax(ibeg:iend), source=nan)
    if ( crop_prog )then
      allocate(pcs%grainc(ibeg:iend), source=nan)
      allocate(pcs%grainc_storage(ibeg:iend), source=nan)
      allocate(pcs%grainc_xfer(ibeg:iend), source=nan)
    end if
#ifdef CN
    allocate(pcs%woodc(ibeg:iend), source=nan)
#endif
  end subroutine init_pft_cstate_type
  !
  ! Initialize pft nitrogen state
  !
  subroutine init_pft_nstate_type(ibeg,iend,pns)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_nstate_type), intent(inout) :: pns !pft nitrogen state

    if ( crop_prog )then
      allocate(pns%grainn(ibeg:iend), source=nan)
      allocate(pns%grainn_storage(ibeg:iend), source=nan)
      allocate(pns%grainn_xfer(ibeg:iend), source=nan)
    end if
    allocate(pns%leafn(ibeg:iend), source=nan)
    allocate(pns%leafn_storage(ibeg:iend), source=nan)
    allocate(pns%leafn_xfer(ibeg:iend), source=nan)
    allocate(pns%frootn(ibeg:iend), source=nan)
    allocate(pns%frootn_storage(ibeg:iend), source=nan)
    allocate(pns%frootn_xfer(ibeg:iend), source=nan)
    allocate(pns%livestemn(ibeg:iend), source=nan)
    allocate(pns%livestemn_storage(ibeg:iend), source=nan)
    allocate(pns%livestemn_xfer(ibeg:iend), source=nan)
    allocate(pns%deadstemn(ibeg:iend), source=nan)
    allocate(pns%deadstemn_storage(ibeg:iend), source=nan)
    allocate(pns%deadstemn_xfer(ibeg:iend), source=nan)
    allocate(pns%livecrootn(ibeg:iend), source=nan)
    allocate(pns%livecrootn_storage(ibeg:iend), source=nan)
    allocate(pns%livecrootn_xfer(ibeg:iend), source=nan)
    allocate(pns%deadcrootn(ibeg:iend), source=nan)
    allocate(pns%deadcrootn_storage(ibeg:iend), source=nan)
    allocate(pns%deadcrootn_xfer(ibeg:iend), source=nan)
    allocate(pns%retransn(ibeg:iend), source=nan)
    allocate(pns%npool(ibeg:iend), source=nan)
    allocate(pns%pft_ntrunc(ibeg:iend), source=nan)
    allocate(pns%dispvegn(ibeg:iend), source=nan)
    allocate(pns%storvegn(ibeg:iend), source=nan)
    allocate(pns%totvegn(ibeg:iend), source=nan)
    allocate(pns%totpftn(ibeg:iend), source=nan)
  end subroutine init_pft_nstate_type
  !
  ! Initialize pft energy flux variables
  !
  subroutine init_pft_eflux_type(ibeg,iend,pef)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_eflux_type), intent(inout) :: pef

    allocate(pef%sabg(ibeg:iend), source=nan)
    allocate(pef%sabv(ibeg:iend), source=nan)
    allocate(pef%fsa(ibeg:iend), source=nan)
    allocate(pef%fsa_u(ibeg:iend), source=nan)
    allocate(pef%fsa_r(ibeg:iend), source=nan)
    allocate(pef%fsr(ibeg:iend), source=nan)
    allocate(pef%parsun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pef%parsha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pef%dlrad(ibeg:iend), source=nan)
    allocate(pef%ulrad(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_tot(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_tot_u(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_tot_r(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_grnd(ibeg:iend), source=nan)
    allocate(pef%eflx_soil_grnd(ibeg:iend), source=nan)
    allocate(pef%eflx_soil_grnd_u(ibeg:iend), source=nan)
    allocate(pef%eflx_soil_grnd_r(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_tot(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_tot_u(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_tot_r(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_grnd(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_veg(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_vege(ibeg:iend), source=nan)
    allocate(pef%eflx_lh_vegt(ibeg:iend), source=nan)
    allocate(pef%eflx_wasteheat_pft(ibeg:iend), source=nan)
    allocate(pef%eflx_heat_from_ac_pft(ibeg:iend), source=nan)
    allocate(pef%eflx_traffic_pft(ibeg:iend), source=nan)
    allocate(pef%eflx_anthro(ibeg:iend), source=nan)
    allocate(pef%cgrnd(ibeg:iend), source=nan)
    allocate(pef%cgrndl(ibeg:iend), source=nan)
    allocate(pef%cgrnds(ibeg:iend), source=nan)
    allocate(pef%eflx_gnet(ibeg:iend), source=nan)
    allocate(pef%eflx_grnd_lake(ibeg:iend), source=nan)
    allocate(pef%dgnetdT(ibeg:iend), source=nan)
    allocate(pef%eflx_lwrad_out(ibeg:iend), source=nan)
    allocate(pef%eflx_lwrad_net(ibeg:iend), source=nan)
    allocate(pef%eflx_lwrad_net_u(ibeg:iend), source=nan)
    allocate(pef%eflx_lwrad_net_r(ibeg:iend), source=nan)
    allocate(pef%netrad(ibeg:iend), source=nan)
    allocate(pef%fsds_vis_d(ibeg:iend), source=nan)
    allocate(pef%fsds_nir_d(ibeg:iend), source=nan)
    allocate(pef%fsds_vis_i(ibeg:iend), source=nan)
    allocate(pef%fsds_nir_i(ibeg:iend), source=nan)
    allocate(pef%fsr_vis_d(ibeg:iend), source=nan)
    allocate(pef%fsr_nir_d(ibeg:iend), source=nan)
    allocate(pef%fsr_vis_i(ibeg:iend), source=nan)
    allocate(pef%fsr_nir_i(ibeg:iend), source=nan)
    allocate(pef%fsds_vis_d_ln(ibeg:iend), source=nan)
    allocate(pef%fsds_vis_i_ln(ibeg:iend), source=nan)
    allocate(pef%parveg_ln(ibeg:iend), source=nan)
    allocate(pef%fsds_nir_d_ln(ibeg:iend), source=nan)
    allocate(pef%fsr_vis_d_ln(ibeg:iend), source=nan)
    allocate(pef%fsr_nir_d_ln(ibeg:iend), source=nan)
    allocate(pef%sabg_lyr(ibeg:iend,-nlevsno+1:1), source=nan)
    allocate(pef%sabg_pen(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_aer(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_bc(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_oc(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_dst(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_aer_sno(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_bc_sno(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_oc_sno(ibeg:iend), source=nan)
    allocate(pef%sfc_frc_dst_sno(ibeg:iend), source=nan)
    allocate(pef%fsr_sno_vd(ibeg:iend), source=nan)
    allocate(pef%fsr_sno_nd(ibeg:iend), source=nan)
    allocate(pef%fsr_sno_vi(ibeg:iend), source=nan)
    allocate(pef%fsr_sno_ni(ibeg:iend), source=nan)
    allocate(pef%fsds_sno_vd(ibeg:iend), source=nan)
    allocate(pef%fsds_sno_nd(ibeg:iend), source=nan)
    allocate(pef%fsds_sno_vi(ibeg:iend), source=nan)
    allocate(pef%fsds_sno_ni(ibeg:iend), source=nan)

    allocate(pef%eflx_sh_snow(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_soil(ibeg:iend), source=nan)
    allocate(pef%eflx_sh_h2osfc(ibeg:iend), source=nan)

    allocate(pef%sabg_soil(ibeg:iend), source=nan)
    allocate(pef%sabg_snow(ibeg:iend), source=nan)
    allocate(pef%sabg_chk(ibeg:iend), source=nan)
  end subroutine init_pft_eflux_type
  !
  ! Initialize pft momentum flux variables
  !
  subroutine init_pft_mflux_type(ibeg,iend,pmf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_mflux_type), intent(inout) :: pmf

    allocate(pmf%taux(ibeg:iend), source=nan)
    allocate(pmf%tauy(ibeg:iend), source=nan)
  end subroutine init_pft_mflux_type
  !
  ! Initialize pft water flux variables
  !
  subroutine init_pft_wflux_type(ibeg,iend,pwf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_wflux_type), intent(inout) :: pwf

    allocate(pwf%qflx_prec_intr(ibeg:iend), source=nan)
    allocate(pwf%qflx_prec_grnd(ibeg:iend), source=nan)
    allocate(pwf%qflx_rain_grnd(ibeg:iend), source=nan)
    allocate(pwf%qflx_snow_grnd(ibeg:iend), source=nan)
    allocate(pwf%qflx_snwcp_liq(ibeg:iend), source=nan)
    allocate(pwf%qflx_snwcp_ice(ibeg:iend), source=nan)
    allocate(pwf%qflx_evap_veg(ibeg:iend), source=nan)
    allocate(pwf%qflx_tran_veg(ibeg:iend), source=nan)
    allocate(pwf%qflx_evap_can(ibeg:iend), source=nan)
    allocate(pwf%qflx_evap_soi(ibeg:iend), source=nan)
    allocate(pwf%qflx_evap_tot(ibeg:iend), source=nan)
    allocate(pwf%qflx_evap_grnd(ibeg:iend), source=0.0_rk8)
    allocate(pwf%qflx_dew_grnd(ibeg:iend), source=0.0_rk8)
    allocate(pwf%qflx_sub_snow(ibeg:iend), source=0.0_rk8)
    allocate(pwf%qflx_dew_snow(ibeg:iend), source=0.0_rk8)

    allocate(pwf%qflx_ev_snow(ibeg:iend), source=nan)
    allocate(pwf%qflx_ev_soil(ibeg:iend), source=nan)
    allocate(pwf%qflx_ev_h2osfc(ibeg:iend), source=nan)
  end subroutine init_pft_wflux_type
  !
  ! Initialize pft carbon flux variables
  !
  subroutine init_pft_cflux_type(ibeg,iend,pcf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_cflux_type), intent(inout) :: pcf

    allocate(pcf%psnsun(ibeg:iend), source=nan)
    allocate(pcf%psnsha(ibeg:iend), source=nan)
    allocate(pcf%psnsun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%psnsha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%cisun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%cisha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%lmrsun(ibeg:iend), source=nan)
    allocate(pcf%lmrsha(ibeg:iend), source=nan)
    allocate(pcf%lmrsun_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%lmrsha_z(ibeg:iend,1:nlevcan), source=nan)
    allocate(pcf%fpsn(ibeg:iend), source=spval)
    allocate(pcf%fco2(ibeg:iend), source=0.0_rk8)
    allocate(pcf%psnsun_wc(ibeg:iend), source=nan)
    allocate(pcf%psnsha_wc(ibeg:iend), source=nan)
    allocate(pcf%fpsn_wc(ibeg:iend), source=nan)
    allocate(pcf%psnsun_wj(ibeg:iend), source=nan)
    allocate(pcf%psnsha_wj(ibeg:iend), source=nan)
    allocate(pcf%fpsn_wj(ibeg:iend), source=nan)
    allocate(pcf%psnsun_wp(ibeg:iend), source=nan)
    allocate(pcf%psnsha_wp(ibeg:iend), source=nan)
    allocate(pcf%fpsn_wp(ibeg:iend), source=nan)

    allocate(pcf%m_leafc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_leafc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_leafc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_leafc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_frootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_frootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_frootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livestemc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livestemc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livestemc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadstemc_to_prod10c(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadstemc_to_prod100c(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadstemc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadstemc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livecrootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livecrootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_livecrootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadcrootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadcrootc_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_gresp_storage_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_gresp_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pcf%hrv_xsmrpool_to_atm(ibeg:iend), source=nan)

    ! fire related variables changed by F. Li and S. Levis
    allocate(pcf%m_leafc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_storage_to_fire(ibeg:iend))
    allocate(pcf%m_livestemc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_storage_to_fire(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_xfer_to_fire(ibeg:iend), source=nan)

    allocate(pcf%m_leafc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_leafc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livestemc_to_deadstemc_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadstemc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_frootc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_livecrootc_to_deadcrootc_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_deadcrootc_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pcf%m_gresp_xfer_to_litter_fire(ibeg:iend), source=nan)

    allocate(pcf%leafc_xfer_to_leafc(ibeg:iend), source=nan)
    allocate(pcf%frootc_xfer_to_frootc(ibeg:iend), source=nan)
    allocate(pcf%livestemc_xfer_to_livestemc(ibeg:iend), source=nan)
    allocate(pcf%deadstemc_xfer_to_deadstemc(ibeg:iend), source=nan)
    allocate(pcf%livecrootc_xfer_to_livecrootc(ibeg:iend), source=nan)
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(ibeg:iend), source=nan)
    allocate(pcf%leafc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%frootc_to_litter(ibeg:iend), source=nan)
    allocate(pcf%leaf_mr(ibeg:iend), source=nan)
    allocate(pcf%froot_mr(ibeg:iend), source=nan)
    allocate(pcf%livestem_mr(ibeg:iend), source=nan)
    allocate(pcf%livecroot_mr(ibeg:iend), source=nan)
    allocate(pcf%grain_mr(ibeg:iend), source=nan)
    allocate(pcf%leaf_curmr(ibeg:iend), source=nan)
    allocate(pcf%froot_curmr(ibeg:iend), source=nan)
    allocate(pcf%livestem_curmr(ibeg:iend), source=nan)
    allocate(pcf%livecroot_curmr(ibeg:iend), source=nan)
    allocate(pcf%grain_curmr(ibeg:iend), source=nan)
    allocate(pcf%leaf_xsmr(ibeg:iend), source=nan)
    allocate(pcf%froot_xsmr(ibeg:iend), source=nan)
    allocate(pcf%livestem_xsmr(ibeg:iend), source=nan)
    allocate(pcf%livecroot_xsmr(ibeg:iend), source=nan)

    allocate(pcf%grain_xsmr(ibeg:iend), source=nan)
    allocate(pcf%psnsun_to_cpool(ibeg:iend), source=nan)
    allocate(pcf%psnshade_to_cpool(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_xsmrpool(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_leafc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_leafc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_frootc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_frootc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_livestemc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_livestemc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_deadstemc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_deadstemc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_livecrootc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_livecrootc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_deadcrootc(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_deadcrootc_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_to_gresp_storage(ibeg:iend), source=nan)
    allocate(pcf%cpool_leaf_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_leaf_storage_gr(ibeg:iend), source=nan)

    allocate(pcf%transfer_leaf_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_froot_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_froot_storage_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_froot_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_livestem_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_livestem_storage_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_livestem_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_deadstem_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_deadstem_storage_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_deadstem_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_livecroot_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_livecroot_storage_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_livecroot_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_deadcroot_gr(ibeg:iend), source=nan)
    allocate(pcf%cpool_deadcroot_storage_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_deadcroot_gr(ibeg:iend), source=nan)
    allocate(pcf%leafc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%frootc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%livestemc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%deadstemc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%livecrootc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%deadcrootc_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%gresp_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pcf%livestemc_to_deadstemc(ibeg:iend), source=nan)
    allocate(pcf%livecrootc_to_deadcrootc(ibeg:iend), source=nan)

    allocate(pcf%gpp(ibeg:iend), source=nan)
    allocate(pcf%mr(ibeg:iend), source=nan)
    allocate(pcf%current_gr(ibeg:iend), source=nan)
    allocate(pcf%transfer_gr(ibeg:iend), source=nan)
    allocate(pcf%storage_gr(ibeg:iend), source=nan)
    allocate(pcf%gr(ibeg:iend), source=nan)
    allocate(pcf%ar(ibeg:iend), source=nan)
    allocate(pcf%rr(ibeg:iend), source=nan)
    allocate(pcf%npp(ibeg:iend), source=nan)
    allocate(pcf%agnpp(ibeg:iend), source=nan)
    allocate(pcf%bgnpp(ibeg:iend), source=nan)
    allocate(pcf%litfall(ibeg:iend), source=nan)
    allocate(pcf%vegfire(ibeg:iend), source=nan)
    allocate(pcf%wood_harvestc(ibeg:iend), source=nan)
    allocate(pcf%pft_cinputs(ibeg:iend), source=nan)
    allocate(pcf%pft_coutputs(ibeg:iend), source=nan)
    allocate(pcf%pft_fire_closs(ibeg:iend), source=nan)

    if ( crop_prog )then
      allocate(pcf%xsmrpool_to_atm(ibeg:iend), source=nan)
      allocate(pcf%grainc_xfer_to_grainc(ibeg:iend), source=nan)
      allocate(pcf%livestemc_to_litter(ibeg:iend), source=nan)
      allocate(pcf%grainc_to_food(ibeg:iend), source=nan)
      allocate(pcf%cpool_to_grainc(ibeg:iend), source=nan)
      allocate(pcf%cpool_to_grainc_storage(ibeg:iend), source=nan)
      allocate(pcf%cpool_grain_gr(ibeg:iend), source=nan)
      allocate(pcf%cpool_grain_storage_gr(ibeg:iend), source=nan)
      allocate(pcf%transfer_grain_gr(ibeg:iend), source=nan)
      allocate(pcf%grainc_storage_to_xfer(ibeg:iend), source=nan)
    end if
#ifdef CN
    allocate(pcf%frootc_alloc(ibeg:iend), source=nan)
    allocate(pcf%frootc_loss(ibeg:iend), source=nan)
    allocate(pcf%leafc_alloc(ibeg:iend), source=nan)
    allocate(pcf%leafc_loss(ibeg:iend), source=nan)
    allocate(pcf%woodc_alloc(ibeg:iend), source=nan)
    allocate(pcf%woodc_loss(ibeg:iend), source=nan)
#endif
#ifdef LCH4
    allocate(pcf%tempavg_agnpp(ibeg:iend), source=spval)
    allocate(pcf%tempavg_bgnpp(ibeg:iend), source=spval)
    allocate(pcf%annavg_agnpp(ibeg:iend), source=spval)
    allocate(pcf%annavg_bgnpp(ibeg:iend), source=spval)
#endif
  end subroutine init_pft_cflux_type
  !
  ! Initialize pft nitrogen flux variables
  !
  subroutine init_pft_nflux_type(ibeg,iend,pnf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_nflux_type), intent(inout) :: pnf

    allocate(pnf%m_leafn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%m_retransn_to_litter(ibeg:iend), source=nan)

    allocate(pnf%hrv_leafn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_frootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_leafn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_frootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_livestemn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadstemn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_livecrootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadcrootn_storage_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_leafn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_frootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_livestemn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadstemn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_livecrootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_livestemn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadstemn_to_prod10n(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadstemn_to_prod100n(ibeg:iend), source=nan)
    allocate(pnf%hrv_livecrootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_deadcrootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%hrv_retransn_to_litter(ibeg:iend), source=nan)

    ! fire varibles changed by F. Li and S. Levis
    allocate(pnf%m_leafn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_storage_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_xfer_to_fire(ibeg:iend), source=nan)
    allocate(pnf%m_retransn_to_fire(ibeg:iend), source=nan)

    allocate(pnf%m_leafn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_leafn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livestemn_to_deadstemn_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadstemn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_frootn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_livecrootn_to_deadcrootn_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_storage_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_deadcrootn_xfer_to_litter_fire(ibeg:iend), source=nan)
    allocate(pnf%m_retransn_to_litter_fire(ibeg:iend), source=nan)

    allocate(pnf%leafn_xfer_to_leafn(ibeg:iend), source=nan)
    allocate(pnf%frootn_xfer_to_frootn(ibeg:iend), source=nan)
    allocate(pnf%livestemn_xfer_to_livestemn(ibeg:iend), source=nan)
    allocate(pnf%deadstemn_xfer_to_deadstemn(ibeg:iend), source=nan)
    allocate(pnf%livecrootn_xfer_to_livecrootn(ibeg:iend), source=nan)
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(ibeg:iend), source=nan)
    allocate(pnf%leafn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%leafn_to_retransn(ibeg:iend), source=nan)
    allocate(pnf%frootn_to_retransn(ibeg:iend), source=nan)
    allocate(pnf%frootn_to_litter(ibeg:iend), source=nan)
    allocate(pnf%retransn_to_npool(ibeg:iend), source=nan)
    allocate(pnf%sminn_to_npool(ibeg:iend), source=nan)

    allocate(pnf%npool_to_leafn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_leafn_storage(ibeg:iend), source=nan)
    allocate(pnf%npool_to_frootn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_frootn_storage(ibeg:iend), source=nan)
    allocate(pnf%npool_to_livestemn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_livestemn_storage(ibeg:iend), source=nan)
    allocate(pnf%npool_to_deadstemn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_deadstemn_storage(ibeg:iend), source=nan)
    allocate(pnf%npool_to_livecrootn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_livecrootn_storage(ibeg:iend), source=nan)
    allocate(pnf%npool_to_deadcrootn(ibeg:iend), source=nan)
    allocate(pnf%npool_to_deadcrootn_storage(ibeg:iend), source=nan)

    allocate(pnf%leafn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%frootn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%livestemn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%deadstemn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%livecrootn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%deadcrootn_storage_to_xfer(ibeg:iend), source=nan)
    allocate(pnf%livestemn_to_deadstemn(ibeg:iend), source=nan)
    allocate(pnf%livestemn_to_retransn(ibeg:iend), source=nan)
    allocate(pnf%livecrootn_to_deadcrootn(ibeg:iend), source=nan)
    allocate(pnf%livecrootn_to_retransn(ibeg:iend), source=nan)
    allocate(pnf%ndeploy(ibeg:iend), source=nan)
    allocate(pnf%pft_ninputs(ibeg:iend), source=nan)
    allocate(pnf%pft_noutputs(ibeg:iend), source=nan)
    allocate(pnf%wood_harvestn(ibeg:iend), source=nan)
    allocate(pnf%pft_fire_nloss(ibeg:iend), source=nan)
    if ( crop_prog )then
      allocate(pnf%grainn_xfer_to_grainn(ibeg:iend), source=nan)
      allocate(pnf%livestemn_to_litter(ibeg:iend), source=nan)
      allocate(pnf%grainn_to_food(ibeg:iend), source=nan)
      allocate(pnf%npool_to_grainn(ibeg:iend), source=nan)
      allocate(pnf%npool_to_grainn_storage(ibeg:iend), source=nan)
      allocate(pnf%grainn_storage_to_xfer(ibeg:iend), source=nan)
      allocate(pnf%fert(ibeg:iend), source=nan)
      allocate(pnf%soyfixn(ibeg:iend), source=nan)
    end if
  end subroutine init_pft_nflux_type
  !
  ! Initialize pft VOC flux variables
  !
  subroutine init_pft_vflux_type(ibeg,iend,pvf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_vflux_type), intent(inout) :: pvf
    integer(ik4) :: i

    if ( shr_megan_megcomps_n < 1 ) return

    allocate(pvf%vocflx_tot(ibeg:iend), source=nan)
    allocate(pvf%vocflx(ibeg:iend,1:shr_megan_megcomps_n), source=nan)
    allocate(pvf%Eopt_out(ibeg:iend), source=nan)
    allocate(pvf%topt_out(ibeg:iend), source=nan)
    allocate(pvf%alpha_out(ibeg:iend), source=nan)
    allocate(pvf%cp_out(ibeg:iend), source=nan)
    allocate(pvf%para_out(ibeg:iend), source=nan)
    allocate(pvf%par24a_out(ibeg:iend), source=nan)
    allocate(pvf%par240a_out(ibeg:iend), source=nan)
    allocate(pvf%paru_out(ibeg:iend), source=nan)
    allocate(pvf%par24u_out(ibeg:iend), source=nan)
    allocate(pvf%par240u_out(ibeg:iend), source=nan)
    allocate(pvf%gamma_out(ibeg:iend), source=nan)
    allocate(pvf%gammaL_out(ibeg:iend), source=nan)
    allocate(pvf%gammaT_out(ibeg:iend), source=nan)
    allocate(pvf%gammaP_out(ibeg:iend), source=nan)
    allocate(pvf%gammaA_out(ibeg:iend), source=nan)
    allocate(pvf%gammaS_out(ibeg:iend), source=nan)
    allocate(pvf%gammaC_out(ibeg:iend), source=nan)

    allocate(pvf%meg(shr_megan_megcomps_n))

    do i = 1, shr_megan_megcomps_n
      allocate(pvf%meg(i)%flux_out(ibeg:iend), source=nan)
    enddo
  end subroutine init_pft_vflux_type
  !
  ! Initialize pft dust flux variables
  !
  subroutine init_pft_dflux_type(ibeg,iend,pdf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_dflux_type), intent(inout) :: pdf

    allocate(pdf%flx_mss_vrt_dst(ibeg:iend,1:ndst), source=nan)
    allocate(pdf%flx_mss_vrt_dst_tot(ibeg:iend), source=nan)
    allocate(pdf%vlc_trb(ibeg:iend,1:ndst), source=nan)
    allocate(pdf%lnd_frc_mbl_dst(ibeg:iend), source=nan)
    allocate(pdf%vlc_trb_1(ibeg:iend), source=nan)
    allocate(pdf%vlc_trb_2(ibeg:iend), source=nan)
    allocate(pdf%vlc_trb_3(ibeg:iend), source=nan)
    allocate(pdf%vlc_trb_4(ibeg:iend), source=nan)
  end subroutine init_pft_dflux_type
  !
  ! Initialize pft dep velocity variables
  !
  subroutine init_pft_depvd_type(ibeg,iend,pdd)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (pft_depvd_type), intent(inout) :: pdd

    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
      allocate(pdd%drydepvel(ibeg:iend,n_drydep), source=nan)
    end if
  end subroutine init_pft_depvd_type
  !
  ! Initialize column physical state variables
  !
  subroutine init_column_pstate_type(ibeg,iend,cps)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_pstate_type), intent(inout) :: cps

    allocate(cps%snl(ibeg:iend))      !* cannot be averaged up
    allocate(cps%isoicol(ibeg:iend), source=bigint)  !* cannot be averaged up
#ifdef HAMSTER_ALBEDO
    allocate(cps%hamster_alb(ibeg:iend,numrad), source=nan)
#endif
    ! F. Li and S. Levis
    allocate(cps%gdp_lf(ibeg:iend), source=nan)
    allocate(cps%peatf_lf(ibeg:iend), source=nan)
    allocate(cps%abm_lf(ibeg:iend), source=13)
    allocate(cps%lgdp_col(ibeg:iend), source=nan)
    allocate(cps%lgdp1_col(ibeg:iend), source=nan)
    allocate(cps%lpop_col(ibeg:iend), source=nan)

    allocate(cps%bsw(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%watsat(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%watfc(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%watdry(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%watopt(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%hksat(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%sucsat(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%csol(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%tkmg(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%tkdry(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%tksatu(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%smpmin(ibeg:iend), source=nan)
    allocate(cps%hkdepth(ibeg:iend), source=nan)
    allocate(cps%wtfact(ibeg:iend), source=nan)
    allocate(cps%fracice(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%gwc_thr(ibeg:iend), source=nan)
    allocate(cps%mss_frc_cly_vld(ibeg:iend), source=nan)
    allocate(cps%mbl_bsn_fct(ibeg:iend), source=nan)
    allocate(cps%do_capsnow(ibeg:iend), source=.false.)
    allocate(cps%snow_depth(ibeg:iend), source=nan)
    allocate(cps%snowdp(ibeg:iend), source=nan)
    allocate(cps%frac_sno(ibeg:iend), source=nan)
    allocate(cps%zi(ibeg:iend,-nlevsno+0:nlevgrnd), source=nan)
    allocate(cps%dz(ibeg:iend,-nlevsno+1:nlevgrnd), source=nan)
    allocate(cps%z (ibeg:iend,-nlevsno+1:nlevgrnd), source=nan)
    allocate(cps%frac_iceold(ibeg:iend,-nlevsno+1:nlevgrnd))
    allocate(cps%imelt(ibeg:iend,-nlevsno+1:nlevgrnd), source=bigint)
    allocate(cps%eff_porosity(ibeg:iend,nlevgrnd), source=spval)
    allocate(cps%emg(ibeg:iend), source=nan)
    allocate(cps%z0mg(ibeg:iend), source=nan)
    allocate(cps%z0hg(ibeg:iend), source=nan)
    allocate(cps%z0qg(ibeg:iend), source=nan)
    allocate(cps%htvp(ibeg:iend), source=nan)
    allocate(cps%beta(ibeg:iend), source=nan)
    allocate(cps%zii(ibeg:iend), source=nan)
    allocate(cps%albgrd(ibeg:iend,numrad), source=nan)
    allocate(cps%albgri(ibeg:iend,numrad), source=nan)
    allocate(cps%rootr_column(ibeg:iend,nlevgrnd), source=spval)
    allocate(cps%rootfr_road_perv(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%rootr_road_perv(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%wf(ibeg:iend), source=nan)
    allocate(cps%wf2(ibeg:iend))
!   allocate(cps%xirrig(ibeg:iend))
    allocate(cps%max_dayl(ibeg:iend))
    allocate(cps%soilpsi(ibeg:iend,nlevgrnd), source=spval)
    allocate(cps%decl(ibeg:iend), source=nan)
    allocate(cps%coszen(ibeg:iend), source=nan)
    allocate(cps%bd(ibeg:iend,nlevgrnd), source=spval)
    allocate(cps%fpi(ibeg:iend), source=nan)
    allocate(cps%fpi_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cps%fpg(ibeg:iend), source=nan)
    allocate(cps%annsum_counter(ibeg:iend), source=nan)
    allocate(cps%cannsum_npp(ibeg:iend), source=nan)
    allocate(cps%col_lag_npp(ibeg:iend), source=spval)
    allocate(cps%cannavg_t2m(ibeg:iend), source=nan)

    ! fire-related variables changed by F. Li and S. Levis
    allocate(cps%nfire(ibeg:iend), source=spval)
    allocate(cps%farea_burned(ibeg:iend), source=nan)
    allocate(cps%fsr_col(ibeg:iend), source=nan)
    allocate(cps%fd_col(ibeg:iend), source=nan)
    allocate(cps%cropf_col(ibeg:iend), source=nan)
    allocate(cps%prec10_col(ibeg:iend), source=nan)
    allocate(cps%prec60_col(ibeg:iend), source=nan)
    allocate(cps%lfc(ibeg:iend), source=spval)
    allocate(cps%lfc2(ibeg:iend), source=0.0_rk8)
    allocate(cps%trotr1_col(ibeg:iend), source=0.0_rk8)
    allocate(cps%trotr2_col(ibeg:iend), source=0.0_rk8)
    allocate(cps%dtrotr_col(ibeg:iend), source=0.0_rk8)
    allocate(cps%baf_crop(ibeg:iend), source=nan)
    allocate(cps%baf_peatf(ibeg:iend), source=nan)
    allocate(cps%fbac(ibeg:iend), source=nan)
    allocate(cps%fbac1(ibeg:iend), source=nan)
    allocate(cps%btran_col(ibeg:iend), source=nan)
    allocate(cps%wtlf(ibeg:iend), source=nan)
    allocate(cps%lfwt(ibeg:iend), source=nan)
    allocate(cps%albsnd_hst(ibeg:iend,numrad), source=spval)
    allocate(cps%albsni_hst(ibeg:iend,numrad), source=spval)
    allocate(cps%albsod(ibeg:iend,numrad), source=spval)
    allocate(cps%albsoi(ibeg:iend,numrad), source=spval)
    allocate(cps%flx_absdv(ibeg:iend,-nlevsno+1:1), source=spval)
    allocate(cps%flx_absdn(ibeg:iend,-nlevsno+1:1), source=spval)
    allocate(cps%flx_absiv(ibeg:iend,-nlevsno+1:1), source=spval)
    allocate(cps%flx_absin(ibeg:iend,-nlevsno+1:1), source=spval)
    allocate(cps%snw_rds(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%snw_rds_top(ibeg:iend), source=nan)
    allocate(cps%sno_liq_top(ibeg:iend), source=nan)
    allocate(cps%mss_bcpho(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_bcphi(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_bctot(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_bc_col(ibeg:iend), source=nan)
    allocate(cps%mss_bc_top(ibeg:iend), source=nan)
    allocate(cps%mss_ocpho(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_ocphi(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_octot(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_oc_col(ibeg:iend), source=nan)
    allocate(cps%mss_oc_top(ibeg:iend), source=nan)
    allocate(cps%mss_dst1(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_dst2(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_dst3(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_dst4(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_dsttot(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_dst_col(ibeg:iend), source=nan)
    allocate(cps%mss_dst_top(ibeg:iend), source=nan)
    allocate(cps%h2osno_top(ibeg:iend), source=nan)
    allocate(cps%mss_cnc_bcphi(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_bcpho(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_ocphi(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_ocpho(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_dst1(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_dst2(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_dst3(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%mss_cnc_dst4(ibeg:iend,-nlevsno+1:0), source=nan)
    allocate(cps%albgrd_pur(ibeg:iend,numrad), source=nan)
    allocate(cps%albgri_pur(ibeg:iend,numrad), source=nan)
    allocate(cps%albgrd_bc(ibeg:iend,numrad), source=nan)
    allocate(cps%albgri_bc(ibeg:iend,numrad), source=nan)
    allocate(cps%albgrd_oc(ibeg:iend,numrad), source=nan)
    allocate(cps%albgri_oc(ibeg:iend,numrad), source=nan)
    allocate(cps%albgrd_dst(ibeg:iend,numrad), source=nan)
    allocate(cps%albgri_dst(ibeg:iend,numrad), source=nan)
    allocate(cps%dTdz_top(ibeg:iend), source=nan)
    allocate(cps%snot_top(ibeg:iend), source=nan)
    ! New variables for "S" Lakes
    allocate(cps%ws(ibeg:iend), source=nan)
    allocate(cps%ks(ibeg:iend), source=nan)
    allocate(cps%dz_lake(ibeg:iend,nlevlak), source=nan)
    allocate(cps%z_lake(ibeg:iend,nlevlak), source=nan)
    ! Initialize to spval so that c->g averaging will be done properly
    allocate(cps%savedtke1(ibeg:iend), source=spval)
    allocate(cps%cellsand(ibeg:iend,nlevsoi), source=nan)
    allocate(cps%cellclay(ibeg:iend,nlevsoi), source=nan)
    allocate(cps%cellorg(ibeg:iend,nlevsoi), source=nan)
    ! Initialize to spval so that it can be a placeholder
    ! for future file input
    allocate(cps%lakedepth(ibeg:iend), source=spval)
    allocate(cps%etal(ibeg:iend), source=nan)
    allocate(cps%lakefetch(ibeg:iend), source=nan)
    ! Initial to spval to detect input from restart file if not arbinit
    allocate(cps%ust_lake(ibeg:iend), source=spval)
    ! End new variables for S Lakes
#if (defined VICHYDRO)
    ! new variables for VIC hydrology
    allocate(cps%b_infil(ibeg:iend), source=nan)
    allocate(cps%dsmax(ibeg:iend), source=nan)
    allocate(cps%ds(ibeg:iend), source=nan)
    allocate(cps%Wsvic(ibeg:iend), source=nan)
    allocate(cps%c_param(ibeg:iend), source=nan)
    allocate(cps%expt(ibeg:iend, nlayer), source=nan)
    allocate(cps%ksat(ibeg:iend, nlayer), source=nan)
    allocate(cps%phi_s(ibeg:iend, nlayer), source=nan)
    allocate(cps%depth(ibeg:iend, nlayert), source=nan)
    allocate(cps%porosity(ibeg:iend, nlayer), source=nan)
    allocate(cps%max_moist(ibeg:iend, nlayer), source=nan)
    allocate(cps%vic_clm_fract(ibeg:iend, nlayer, nlevsoi), source=nan)
#endif
#ifdef LCH4
    ! New variable for finundated parameterization
    !allocate(cps%zwt0(ibeg:iend), source=nan)
    !allocate(cps%f0(ibeg:iend), source=nan)
    !allocate(cps%p3(ibeg:iend), source=nan)
    allocate(cps%k(ibeg:iend), source=nan)
    allocate(cps%q(ibeg:iend), source=nan)
    allocate(cps%v(ibeg:iend), source=nan)
    allocate(cps%maxf(ibeg:iend), source=nan)
    ! New variable for methane
    allocate(cps%pH(ibeg:iend), source=nan)
#endif

#ifdef CN
    allocate(cps%q10(ibeg:iend), source=nan)
    allocate(cps%ndep(ibeg:iend), source=nan)
#endif

    allocate(cps%irrig_rate(ibeg:iend), source=nan)
    allocate(cps%n_irrig_steps_left(ibeg:iend), source=0)
    allocate(cps%forc_pbot(ibeg:iend), source=nan)
    allocate(cps%forc_rho(ibeg:iend), source=nan)

    allocate(cps%rf_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
                     1:ndecomp_cascade_transitions), source=nan)
    allocate(cps%pathfrac_decomp_cascade(ibeg:iend,1:nlevdecomp_full, &
                     1:ndecomp_cascade_transitions), source=nan)
    allocate(cps%nfixation_prof(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cps%ndep_prof(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cps%alt(ibeg:iend), source=spval)
    allocate(cps%altmax(ibeg:iend), source=spval)
    allocate(cps%altmax_lastyear(ibeg:iend), source=spval)
    allocate(cps%alt_indx(ibeg:iend), source=bigint)
    allocate(cps%altmax_indx(ibeg:iend), source=bigint)
    allocate(cps%altmax_lastyear_indx(ibeg:iend), source=bigint)
    allocate(cps%som_adv_coef(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cps%som_diffus_coef(ibeg:iend,1:nlevdecomp_full), source=spval)

    allocate(cps%frac_sno_eff(ibeg:iend), source=spval)
    allocate(cps%topo_std(ibeg:iend), source=nan)
    allocate(cps%topo_ndx(ibeg:iend), source=nan)
    allocate(cps%topo_slope(ibeg:iend), source=nan)
    allocate(cps%hksat_min(ibeg:iend,nlevgrnd), source=nan)
    allocate(cps%frac_h2osfc(ibeg:iend), source=spval)
    allocate(cps%micro_sigma(ibeg:iend), source=nan)
    allocate(cps%h2osfc_thresh(ibeg:iend), source=nan)
    allocate(cps%frac_h2osfc_temp(ibeg:iend), source=0.0_rk8)
    allocate(cps%n_melt(ibeg:iend), source=nan)
  end subroutine init_column_pstate_type
  !
  ! Initialize column energy state variables
  !
  subroutine init_column_estate_type(ibeg,iend,ces)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_estate_type), intent(inout) :: ces

    allocate(ces%t_grnd(ibeg:iend), source=nan)
    allocate(ces%t_grnd_u(ibeg:iend), source=nan)
    allocate(ces%t_grnd_r(ibeg:iend), source=nan)
    allocate(ces%dt_grnd(ibeg:iend), source=nan)
    allocate(ces%t_soisno(ibeg:iend,-nlevsno+1:nlevgrnd), source=spval)
    allocate(ces%t_soi_10cm(ibeg:iend), source=spval)
    allocate(ces%tsoi17(ibeg:iend), source=spval)
    allocate(ces%t_lake(ibeg:iend,1:nlevlak), source=nan)
    allocate(ces%tssbef(ibeg:iend,-nlevsno+1:nlevgrnd), source=nan)
    allocate(ces%thv(ibeg:iend), source=nan)
    allocate(ces%hc_soi(ibeg:iend), source=nan)
    allocate(ces%hc_soisno(ibeg:iend), source=nan)
    allocate(ces%forc_t(ibeg:iend), source=nan)
    allocate(ces%forc_th(ibeg:iend), source=nan)

    allocate(ces%t_h2osfc(ibeg:iend), source=spval)
    allocate(ces%t_h2osfc_bef(ibeg:iend), source=nan)
  end subroutine init_column_estate_type
  !
  ! Initialize column water state variables
  !
  subroutine init_column_wstate_type(ibeg,iend,cws)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_wstate_type), intent(inout) :: cws !column water state

    allocate(cws%h2osno(ibeg:iend), source=nan)
    allocate(cws%errh2osno(ibeg:iend), source=nan)
    allocate(cws%snow_sources(ibeg:iend), source=nan)
    allocate(cws%snow_sinks(ibeg:iend), source=nan)
    allocate(cws%h2osoi_liq(ibeg:iend,-nlevsno+1:nlevgrnd), source=spval)
    allocate(cws%h2osoi_ice(ibeg:iend,-nlevsno+1:nlevgrnd), source=spval)
    allocate(cws%h2osoi_liqice_10cm(ibeg:iend), source=spval)
    allocate(cws%h2osoi_vol(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cws%h2osno_old(ibeg:iend), source=nan)
    allocate(cws%qg(ibeg:iend), source=nan)
    allocate(cws%dqgdT(ibeg:iend), source=nan)
    allocate(cws%snowice(ibeg:iend), source=nan)
    allocate(cws%snowliq(ibeg:iend), source=nan)
    allocate(cws%soilalpha(ibeg:iend), source=nan)
    allocate(cws%soilbeta(ibeg:iend), source=nan)
    allocate(cws%soilalpha_u(ibeg:iend), source=nan)
    allocate(cws%zwt(ibeg:iend), source=nan)
    allocate(cws%fcov(ibeg:iend), source=nan)
    allocate(cws%fsat(ibeg:iend), source=nan)
    !New variable for methane code
#ifdef LCH4
    allocate(cws%finundated(ibeg:iend), source=nan)
#endif
    allocate(cws%wa(ibeg:iend), source=spval)
    allocate(cws%qcharge(ibeg:iend), source=nan)
    allocate(cws%smp_l(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cws%hk_l(ibeg:iend,1:nlevgrnd), source=spval)
    ! New variables for "S" lakes
    ! Initialize to spval for detection of whether
    ! it has not been initialized by file input
    allocate(cws%lake_icefrac(ibeg:iend,1:nlevlak), source=spval)
    ! and so c->g averaging will be done properly
    allocate(cws%lake_icethick(ibeg:iend), source=nan)
    ! End new variables for S lakes
    allocate(cws%forc_q(ibeg:iend), source=nan)
#if (defined VICHYDRO)
    allocate(cws%moist(ibeg:iend,1:nlayert), source=spval)
    allocate(cws%ice(ibeg:iend,1:nlayert), source=spval)
    allocate(cws%moist_vol(ibeg:iend,1:nlayert), source=spval)
    allocate(cws%max_infil(ibeg:iend), source=spval)
    allocate(cws%i_0(ibeg:iend), source=spval)
#endif

    allocate(cws%h2osfc(ibeg:iend), source=spval)
    allocate(cws%qg_snow(ibeg:iend), source=nan)
    allocate(cws%qg_soil(ibeg:iend), source=nan)
    allocate(cws%qg_h2osfc(ibeg:iend), source=nan)
    allocate(cws%frost_table(ibeg:iend), source=spval)
    allocate(cws%zwt_perched(ibeg:iend), source=spval)
    allocate(cws%int_snow(ibeg:iend), source=spval)
    allocate(cws%swe_old(ibeg:iend,-nlevsno+1:0), source=nan)
  end subroutine init_column_wstate_type
  !
  ! Initialize column carbon state variables
  !
  subroutine init_column_cstate_type(ibeg,iend,ccs)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_cstate_type), intent(inout) :: ccs

    allocate(ccs%soilc(ibeg:iend), source=nan)
    allocate(ccs%cwdc(ibeg:iend), source=nan)
    allocate(ccs%col_ctrunc(ibeg:iend), source=nan)
    allocate(ccs%decomp_cpools_vr(ibeg:iend, &
             1:nlevdecomp_full,1:ndecomp_pools), source=nan)
    allocate(ccs%decomp_cpools(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(ccs%decomp_cpools_1m(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(ccs%col_ctrunc_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(ccs%seedc(ibeg:iend), source=nan)
    allocate(ccs%prod10c(ibeg:iend), source=nan)
    allocate(ccs%prod100c(ibeg:iend), source=nan)
    allocate(ccs%totprodc(ibeg:iend), source=nan)
    allocate(ccs%totlitc(ibeg:iend), source=nan)
    allocate(ccs%totsomc(ibeg:iend), source=nan)
    allocate(ccs%totlitc_1m(ibeg:iend), source=nan)
    allocate(ccs%totsomc_1m(ibeg:iend), source=nan)
    allocate(ccs%totecosysc(ibeg:iend), source=nan)
    allocate(ccs%totcolc(ibeg:iend), source=nan)

    !F. Li and S. Levis
    allocate(ccs%rootc_col(ibeg:iend), source=nan)
    allocate(ccs%totvegc_col(ibeg:iend), source=nan)
    allocate(ccs%leafc_col(ibeg:iend), source=nan)
    allocate(ccs%fuelc(ibeg:iend), source=spval)
    allocate(ccs%fuelc_crop(ibeg:iend), source=nan)
  end subroutine init_column_cstate_type
  !
  ! Initialize column nitrogen state variables
  !
  subroutine init_column_nstate_type(ibeg,iend,cns)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_nstate_type), intent(inout) :: cns

    allocate(cns%decomp_npools(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(cns%decomp_npools_1m(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(cns%decomp_npools_vr(ibeg:iend, &
             1:nlevdecomp_full,1:ndecomp_pools), source=nan)
    allocate(cns%sminn_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cns%col_ntrunc_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
#ifdef NITRIF_DENITRIF
    allocate(cns%smin_no3_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cns%smin_nh4_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cns%smin_no3(ibeg:iend), source=nan)
    allocate(cns%smin_nh4(ibeg:iend), source=nan)
#endif
    allocate(cns%cwdn(ibeg:iend), source=nan)
    allocate(cns%sminn(ibeg:iend), source=nan)
    allocate(cns%col_ntrunc(ibeg:iend), source=nan)
    allocate(cns%seedn(ibeg:iend), source=nan)
    allocate(cns%prod10n(ibeg:iend), source=nan)
    allocate(cns%prod100n(ibeg:iend), source=nan)
    allocate(cns%totprodn(ibeg:iend), source=nan)
    allocate(cns%totlitn(ibeg:iend), source=nan)
    allocate(cns%totsomn(ibeg:iend), source=nan)
    allocate(cns%totlitn_1m(ibeg:iend), source=nan)
    allocate(cns%totsomn_1m(ibeg:iend), source=nan)
    allocate(cns%totecosysn(ibeg:iend), source=nan)
    allocate(cns%totcoln(ibeg:iend), source=nan)
  end subroutine init_column_nstate_type
  !
  ! Initialize column energy flux variables
  !
  subroutine init_column_eflux_type(ibeg,iend,cef)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_eflux_type), intent(inout) :: cef

    allocate(cef%eflx_snomelt(ibeg:iend), source=spval)
    allocate(cef%eflx_snomelt_u(ibeg:iend), source=spval)
    allocate(cef%eflx_snomelt_r(ibeg:iend), source=spval)
    allocate(cef%eflx_impsoil(ibeg:iend), source=nan)
    allocate(cef%eflx_fgr12(ibeg:iend), source=nan)
    allocate(cef%eflx_fgr(ibeg:iend, 1:nlevgrnd), source=nan)
    allocate(cef%eflx_building_heat(ibeg:iend), source=nan)
    allocate(cef%eflx_urban_ac(ibeg:iend), source=nan)
    allocate(cef%eflx_urban_heat(ibeg:iend), source=nan)
    allocate(cef%eflx_bot(ibeg:iend), source=nan)
  end subroutine init_column_eflux_type
  !
  ! Initialize column water flux variables
  !
  subroutine init_column_wflux_type(ibeg,iend,cwf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_wflux_type), intent(inout) :: cwf

    allocate(cwf%qflx_infl(ibeg:iend), source=nan)
    allocate(cwf%qflx_surf(ibeg:iend), source=nan)
    allocate(cwf%qflx_drain(ibeg:iend), source=nan)
    allocate(cwf%qflx_top_soil(ibeg:iend), source=spval)
    allocate(cwf%qflx_sl_top_soil(ibeg:iend), source=nan)
    allocate(cwf%qflx_snomelt(ibeg:iend), source=nan)
    allocate(cwf%qflx_qrgwl(ibeg:iend), source=nan)
    allocate(cwf%qflx_runoff(ibeg:iend), source=nan)
    allocate(cwf%qflx_runoff_u(ibeg:iend), source=nan)
    allocate(cwf%qflx_runoff_r(ibeg:iend), source=nan)
    allocate(cwf%qmelt(ibeg:iend), source=nan)
    allocate(cwf%h2ocan_loss(ibeg:iend), source=nan)
    allocate(cwf%qflx_rsub_sat(ibeg:iend), source=spval)
    allocate(cwf%flx_bc_dep_dry(ibeg:iend), source=nan)
    allocate(cwf%flx_bc_dep_wet(ibeg:iend), source=nan)
    allocate(cwf%flx_bc_dep_pho(ibeg:iend), source=nan)
    allocate(cwf%flx_bc_dep_phi(ibeg:iend), source=nan)
    allocate(cwf%flx_bc_dep(ibeg:iend), source=nan)
    allocate(cwf%flx_oc_dep_dry(ibeg:iend), source=nan)
    allocate(cwf%flx_oc_dep_wet(ibeg:iend), source=nan)
    allocate(cwf%flx_oc_dep_pho(ibeg:iend), source=nan)
    allocate(cwf%flx_oc_dep_phi(ibeg:iend), source=nan)
    allocate(cwf%flx_oc_dep(ibeg:iend), source=nan)

    allocate(cwf%flx_dst_dep_dry1(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_wet1(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_dry2(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_wet2(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_dry3(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_wet3(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_dry4(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep_wet4(ibeg:iend), source=nan)
    allocate(cwf%flx_dst_dep(ibeg:iend), source=nan)
    allocate(cwf%qflx_snofrz_lyr(ibeg:iend,-nlevsno+1:0), source=spval)
    allocate(cwf%qflx_snofrz_col(ibeg:iend), source=nan)
    allocate(cwf%qflx_irrig(ibeg:iend), source=spval)
    allocate(cwf%qflx_floodc(ibeg:iend), source=spval)

    allocate(cwf%qflx_h2osfc_to_ice(ibeg:iend), source=spval)
    allocate(cwf%qflx_h2osfc_surf(ibeg:iend), source=spval)
    allocate(cwf%qflx_snow_h2osfc(ibeg:iend), source=nan)
    allocate(cwf%qflx_drain_perched(ibeg:iend), source=spval)
    allocate(cwf%qflx_floodc(ibeg:iend), source=spval)
    allocate(cwf%qflx_snow_melt(ibeg:iend), source=spval)
  end subroutine init_column_wflux_type
  !
  ! Initialize column carbon flux variables
  !
  subroutine init_column_cflux_type(ibeg,iend,ccf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_cflux_type), intent(inout) :: ccf

    allocate(ccf%hrv_deadstemc_to_prod10c(ibeg:iend), source=nan)
    allocate(ccf%hrv_deadstemc_to_prod100c(ibeg:iend), source=nan)
    allocate(ccf%m_decomp_cpools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)
    allocate(ccf%m_decomp_cpools_to_fire(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(ccf%decomp_cascade_hr_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions), source=nan)
    allocate(ccf%decomp_cascade_hr(ibeg:iend, &
            1:ndecomp_cascade_transitions), source=nan)
    allocate(ccf%decomp_cascade_ctransfer_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions), source=nan)
    allocate(ccf%decomp_cascade_ctransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions), source=nan)
    allocate(ccf%decomp_cpools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)
    allocate(ccf%decomp_k(ibeg:iend,1:nlevdecomp_full, &
            1:ndecomp_cascade_transitions), source=spval)
    ! Initialize these four below to spval to allow history to not average
    ! over inactive points.
    allocate(ccf%t_scalar(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(ccf%w_scalar(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(ccf%hr_vr(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(ccf%o_scalar(ibeg:iend,1:nlevdecomp_full), source=spval)

    allocate(ccf%som_c_leached(ibeg:iend), source=nan)
    allocate(ccf%decomp_cpools_leached(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(ccf%decomp_cpools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)

    allocate(ccf%phenology_c_to_litr_met_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%phenology_c_to_litr_cel_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%phenology_c_to_litr_lig_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%gap_mortality_c_to_litr_met_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%gap_mortality_c_to_litr_cel_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%gap_mortality_c_to_litr_lig_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%gap_mortality_c_to_cwdc(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%fire_mortality_c_to_cwdc(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%m_c_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(ccf%m_c_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(ccf%m_c_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(ccf%harvest_c_to_litr_met_c(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(ccf%harvest_c_to_litr_cel_c(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(ccf%harvest_c_to_litr_lig_c(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(ccf%harvest_c_to_cwdc(ibeg:iend, 1:nlevdecomp_full), source=nan)

#ifdef NITRIF_DENITRIF
    allocate(ccf%phr_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
#endif

#ifdef CN
    !F. Li and S. Levis
    allocate(ccf%somc_fire(ibeg:iend), source=nan)
    allocate(ccf%lf_conv_cflux(ibeg:iend), source=nan)
    allocate(ccf%dwt_seedc_to_leaf(ibeg:iend), source=nan)
    allocate(ccf%dwt_seedc_to_deadstem(ibeg:iend), source=nan)
    allocate(ccf%dwt_conv_cflux(ibeg:iend), source=nan)
    allocate(ccf%dwt_prod10c_gain(ibeg:iend), source=nan)
    allocate(ccf%dwt_prod100c_gain(ibeg:iend), source=nan)
    allocate(ccf%dwt_frootc_to_litr_met_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%dwt_frootc_to_litr_cel_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%dwt_frootc_to_litr_lig_c(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%dwt_livecrootc_to_cwdc(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%dwt_deadcrootc_to_cwdc(ibeg:iend, &
            1:nlevdecomp_full), source=nan)
    allocate(ccf%dwt_closs(ibeg:iend), source=nan)
    allocate(ccf%landuseflux(ibeg:iend), source=nan)
    allocate(ccf%landuptake(ibeg:iend), source=nan)
    allocate(ccf%prod10c_loss(ibeg:iend), source=nan)
    allocate(ccf%prod100c_loss(ibeg:iend), source=nan)
    allocate(ccf%product_closs(ibeg:iend), source=nan)
#endif
    allocate(ccf%lithr(ibeg:iend), source=nan)
    allocate(ccf%somhr(ibeg:iend), source=nan)
    allocate(ccf%hr(ibeg:iend), source=nan)
    allocate(ccf%sr(ibeg:iend), source=nan)
    allocate(ccf%er(ibeg:iend), source=nan)
    allocate(ccf%litfire(ibeg:iend), source=nan)
    allocate(ccf%somfire(ibeg:iend), source=nan)
    allocate(ccf%totfire(ibeg:iend), source=nan)
    allocate(ccf%nep(ibeg:iend), source=nan)
    allocate(ccf%nbp(ibeg:iend), source=nan)
    allocate(ccf%nee(ibeg:iend), source=nan)
    allocate(ccf%col_cinputs(ibeg:iend), source=nan)
    allocate(ccf%col_coutputs(ibeg:iend), source=nan)
    allocate(ccf%col_fire_closs(ibeg:iend), source=nan)

#if (defined CN)
    allocate(ccf%cwdc_hr(ibeg:iend), source=nan)
    allocate(ccf%cwdc_loss(ibeg:iend), source=nan)
    allocate(ccf%litterc_loss(ibeg:iend), source=nan)
#endif
  end subroutine init_column_cflux_type
#ifdef LCH4
  !
  ! Initialize column methane flux variables
  !
  subroutine init_column_ch4_type(ibeg,iend,cch4)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_ch4_type), intent(inout) :: cch4

    allocate(cch4%ch4_prod_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_prod_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_prod_depth_lake(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_oxid_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_oxid_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_oxid_depth_lake(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%o2_oxid_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%o2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%o2_decomp_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    ! To detect first time-step for denitrification code
    allocate(cch4%o2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%o2_aere_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%o2_aere_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_decomp_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_decomp_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_oxid_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_oxid_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_aere_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_aere_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_tran_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_tran_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_aere_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%co2_aere_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_surf_aere_sat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_aere_unsat(ibeg:iend), source=nan)
    allocate(cch4%ch4_ebul_depth_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_ebul_depth_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_ebul_total_sat(ibeg:iend), source=nan)
    allocate(cch4%ch4_ebul_total_unsat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_ebul_sat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_ebul_unsat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_ebul_lake(ibeg:iend), source=nan)
    allocate(cch4%conc_ch4_sat(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%conc_ch4_unsat(ibeg:iend,1:nlevgrnd), source=spval)

    ! Just a diagnostic, so nan is fine
    allocate(cch4%conc_ch4_lake(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_surf_diff_sat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_diff_unsat(ibeg:iend), source=nan)
    allocate(cch4%ch4_surf_diff_lake(ibeg:iend), source=nan)
    allocate(cch4%conc_o2_sat(ibeg:iend,1:nlevgrnd), source=spval)

    ! To detect file input and detect first time-step for denitrification code
    allocate(cch4%conc_o2_unsat(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%conc_o2_lake(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4_dfsat_flux(ibeg:iend), source=nan)
    allocate(cch4%zwt_ch4_unsat(ibeg:iend), source=nan)
    allocate(cch4%fsat_bef(ibeg:iend), source=spval)
    allocate(cch4%lake_soilc(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%lake_raw(ibeg:iend), source=nan)
    allocate(cch4%totcolch4(ibeg:iend), source=spval)

    allocate(cch4%fphr(ibeg:iend,1:nlevgrnd))
    allocate(cch4%annsum_counter(ibeg:iend), source=spval)
    allocate(cch4%tempavg_somhr(ibeg:iend), source=nan)
    allocate(cch4%annavg_somhr(ibeg:iend), source=spval)
    allocate(cch4%tempavg_finrw(ibeg:iend), source=nan)
    allocate(cch4%annavg_finrw(ibeg:iend), source=spval)
    allocate(cch4%sif(ibeg:iend), source=nan)
    allocate(cch4%o2stress_unsat(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%o2stress_sat(ibeg:iend,1:nlevgrnd), source=spval)
    allocate(cch4%ch4stress_unsat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%ch4stress_sat(ibeg:iend,1:nlevgrnd), source=nan)
    allocate(cch4%qflx_surf_lag(ibeg:iend), source=spval)
    allocate(cch4%finundated_lag(ibeg:iend), source=spval)
    allocate(cch4%layer_sat_lag(ibeg:iend,1:nlevgrnd), source=spval)
  end subroutine init_column_ch4_type
#endif
  !
  ! Initialize column nitrogen flux variables
  !
  subroutine init_column_nflux_type(ibeg,iend,cnf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (column_nflux_type), intent(inout) :: cnf

    allocate(cnf%ndep_to_sminn(ibeg:iend), source=nan)
    allocate(cnf%nfix_to_sminn(ibeg:iend), source=nan)
    allocate(cnf%fert_to_sminn(ibeg:iend), source=nan)
    allocate(cnf%soyfixn_to_sminn(ibeg:iend), source=nan)
    allocate(cnf%hrv_deadstemn_to_prod10n(ibeg:iend), source=nan)
    allocate(cnf%hrv_deadstemn_to_prod100n(ibeg:iend), source=nan)

    allocate(cnf%m_n_to_litr_met_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%m_n_to_litr_cel_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%m_n_to_litr_lig_fire(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%sminn_to_plant(ibeg:iend), source=nan)
    allocate(cnf%potential_immob(ibeg:iend), source=nan)
    allocate(cnf%actual_immob(ibeg:iend), source=nan)
    allocate(cnf%gross_nmin(ibeg:iend), source=nan)
    allocate(cnf%net_nmin(ibeg:iend), source=nan)
    allocate(cnf%denit(ibeg:iend), source=nan)
    allocate(cnf%supplement_to_sminn(ibeg:iend), source=nan)

    allocate(cnf%m_decomp_npools_to_fire_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)
    allocate(cnf%m_decomp_npools_to_fire(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(cnf%decomp_cascade_ntransfer_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%decomp_cascade_ntransfer(ibeg:iend, &
            1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%decomp_cascade_sminn_flux_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%decomp_cascade_sminn_flux(ibeg:iend, &
            1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%decomp_npools_sourcesink(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)

    allocate(cnf%phenology_n_to_litr_met_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%phenology_n_to_litr_cel_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%phenology_n_to_litr_lig_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%gap_mortality_n_to_litr_met_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%gap_mortality_n_to_litr_cel_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%gap_mortality_n_to_litr_lig_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%gap_mortality_n_to_cwdn(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%fire_mortality_n_to_cwdn(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%harvest_n_to_litr_met_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%harvest_n_to_litr_cel_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%harvest_n_to_litr_lig_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%harvest_n_to_cwdn(ibeg:iend, &
             1:nlevdecomp_full), source=nan)

#ifndef NITRIF_DENITRIF
    allocate(cnf%sminn_to_denit_decomp_cascade_vr(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%sminn_to_denit_decomp_cascade(ibeg:iend, &
            1:ndecomp_cascade_transitions), source=nan)
    allocate(cnf%sminn_to_denit_excess_vr(ibeg:iend,1:nlevdecomp_full))
    allocate(cnf%sminn_to_denit_excess(ibeg:iend), source=nan)
    allocate(cnf%sminn_leached_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%sminn_leached(ibeg:iend), source=nan)
#else
    allocate(cnf%f_nit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_denit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%smin_no3_leached_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%smin_no3_leached(ibeg:iend), source=nan)
    allocate(cnf%smin_no3_runoff_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%smin_no3_runoff(ibeg:iend), source=nan)
    allocate(cnf%pot_f_nit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%pot_f_nit(ibeg:iend), source=nan)
    allocate(cnf%pot_f_denit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%pot_f_denit(ibeg:iend), source=nan)
    allocate(cnf%actual_immob_no3_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%actual_immob_nh4_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%smin_no3_to_plant_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%smin_nh4_to_plant_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_nit(ibeg:iend), source=nan)
    allocate(cnf%f_denit(ibeg:iend), source=nan)
    allocate(cnf%n2_n2o_ratio_denit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_n2o_denit(ibeg:iend), source=nan)
    allocate(cnf%f_n2o_denit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_n2o_nit(ibeg:iend), source=nan)
    allocate(cnf%f_n2o_nit_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_n2o_tot(ibeg:iend), source=nan)
    allocate(cnf%f_n2o_tot_vr(ibeg:iend,1:nlevdecomp_full), source=nan)

    allocate(cnf%smin_no3_massdens_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%soil_bulkdensity(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%k_nitr_t_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%k_nitr_ph_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%k_nitr_h2o_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%k_nitr_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%wfps_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%fmax_denit_carbonsubstrate_vr(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%fmax_denit_nitrate_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%f_denit_base_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%diffus(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cnf%ratio_k1(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%ratio_no3_co2(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cnf%soil_co2_prod(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%fr_WFPS(ibeg:iend,1:nlevdecomp_full), source=spval)

    allocate(cnf%r_psi(ibeg:iend,1:nlevdecomp_full), source=spval)
    allocate(cnf%anaerobic_frac(ibeg:iend,1:nlevdecomp_full), source=spval)
#endif
    allocate(cnf%potential_immob_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%actual_immob_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%sminn_to_plant_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%supplement_to_sminn_vr(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%gross_nmin_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%net_nmin_vr(ibeg:iend,1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_seedn_to_leaf(ibeg:iend), source=nan)
    allocate(cnf%dwt_seedn_to_deadstem(ibeg:iend), source=nan)
    allocate(cnf%dwt_conv_nflux(ibeg:iend), source=nan)
    allocate(cnf%dwt_prod10n_gain(ibeg:iend), source=nan)
    allocate(cnf%dwt_prod100n_gain(ibeg:iend), source=nan)
    allocate(cnf%dwt_frootn_to_litr_met_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_frootn_to_litr_cel_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_frootn_to_litr_lig_n(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_livecrootn_to_cwdn(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_deadcrootn_to_cwdn(ibeg:iend, &
             1:nlevdecomp_full), source=nan)
    allocate(cnf%dwt_nloss(ibeg:iend), source=nan)
    allocate(cnf%prod10n_loss(ibeg:iend), source=nan)
    allocate(cnf%prod100n_loss(ibeg:iend), source=nan)
    allocate(cnf%product_nloss(ibeg:iend), source=nan)
    allocate(cnf%col_ninputs(ibeg:iend), source=nan)
    allocate(cnf%col_noutputs(ibeg:iend), source=nan)
    allocate(cnf%col_fire_nloss(ibeg:iend), source=nan)
    allocate(cnf%som_n_leached(ibeg:iend), source=nan)
    allocate(cnf%decomp_npools_leached(ibeg:iend,1:ndecomp_pools), source=nan)
    allocate(cnf%decomp_npools_transport_tendency(ibeg:iend, &
            1:nlevdecomp_full,1:ndecomp_pools), source=nan)
  end subroutine init_column_nflux_type
  !
  ! Initialize landunit physical state variables
  !
  subroutine init_landunit_pstate_type(ibeg,iend,lps)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (landunit_pstate_type), intent(inout) :: lps

    allocate(lps%t_building(ibeg:iend), source=nan)
    allocate(lps%t_building_max(ibeg:iend), source=nan)
    allocate(lps%t_building_min(ibeg:iend), source=nan)
    if ( nlevurb > 0 )then
      allocate(lps%tk_wall(ibeg:iend,nlevurb), source=nan)
      allocate(lps%tk_roof(ibeg:iend,nlevurb), source=nan)
      allocate(lps%cv_wall(ibeg:iend,nlevurb), source=nan)
      allocate(lps%cv_roof(ibeg:iend,nlevurb), source=nan)
    end if
    allocate(lps%tk_improad(ibeg:iend,nlevurb), source=nan)
    allocate(lps%cv_improad(ibeg:iend,nlevurb), source=nan)
    allocate(lps%thick_wall(ibeg:iend), source=nan)
    allocate(lps%thick_roof(ibeg:iend), source=nan)
    allocate(lps%nlev_improad(ibeg:iend), source=bigint)
    allocate(lps%vf_sr(ibeg:iend), source=nan)
    allocate(lps%vf_wr(ibeg:iend), source=nan)
    allocate(lps%vf_sw(ibeg:iend), source=nan)
    allocate(lps%vf_rw(ibeg:iend), source=nan)
    allocate(lps%vf_ww(ibeg:iend), source=nan)
    allocate(lps%taf(ibeg:iend), source=nan)
    allocate(lps%qaf(ibeg:iend), source=nan)
    allocate(lps%sabs_roof_dir(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_roof_dif(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_sunwall_dir(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_sunwall_dif(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_shadewall_dir(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_shadewall_dif(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_improad_dir(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_improad_dif(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_perroad_dir(ibeg:iend,1:numrad), source=nan)
    allocate(lps%sabs_perroad_dif(ibeg:iend,1:numrad), source=nan)
  end subroutine init_landunit_pstate_type
  !
  ! Initialize landunit energy flux variables
  !
  subroutine init_landunit_eflux_type(ibeg,iend,lef)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (landunit_eflux_type), intent(inout) :: lef

    allocate(lef%eflx_traffic(ibeg:iend), source=nan)
    allocate(lef%eflx_traffic_factor(ibeg:iend), source=nan)
    allocate(lef%eflx_wasteheat(ibeg:iend), source=nan)
    allocate(lef%eflx_heat_from_ac(ibeg:iend), source=nan)
  end subroutine init_landunit_eflux_type

#if (defined CNDV)
  !
  ! Initialize gridcell DGVM variables
  !
  subroutine init_gridcell_dgvstate_type(ibeg,iend,gps)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_dgvstate_type), intent(inout) :: gps

    allocate(gps%agdd20(ibeg:iend), source=nan)
    allocate(gps%tmomin20(ibeg:iend), source=nan)
    allocate(gps%t10min(ibeg:iend), source=nan)
  end subroutine init_gridcell_dgvstate_type
#endif
  !
  ! Initialize gridcell isoprene emission factor variables
  !
  subroutine init_gridcell_efstate_type(ibeg,iend,gve)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_efstate_type), intent(inout) :: gve

    allocate(gve%efisop(6,ibeg:iend), source=nan)
  end subroutine init_gridcell_efstate_type
  !
  ! Initialize gridcell water flux variables
  !
  subroutine init_gridcell_wflux_type(ibeg,iend,gwf)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_wflux_type), intent(inout) :: gwf

    allocate(gwf%qflx_runoffg(ibeg:iend), source=0.0_rk8)
    allocate(gwf%qflx_snwcp_iceg(ibeg:iend), source=0.0_rk8)
    allocate(gwf%qflx_liq_dynbal(ibeg:iend), source=nan)
    allocate(gwf%qflx_ice_dynbal(ibeg:iend), source=nan)
  end subroutine init_gridcell_wflux_type
  !
  ! Initialize gridcell energy flux variables
  !
  subroutine init_gridcell_eflux_type(ibeg,iend,gef)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_eflux_type), intent(inout) :: gef

    allocate(gef%eflx_sh_totg(ibeg:iend), source=nan)
    allocate(gef%eflx_dynbal(ibeg:iend), source=nan)
  end subroutine init_gridcell_eflux_type
  !
  ! Initialize gridcell water state variables
  !
  subroutine init_gridcell_wstate_type(ibeg,iend,gws)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_wstate_type), intent(inout) :: gws

    allocate(gws%gc_liq1(ibeg:iend), source=nan)
    allocate(gws%gc_liq2(ibeg:iend), source=nan)
    allocate(gws%gc_ice1(ibeg:iend), source=nan)
    allocate(gws%gc_ice2(ibeg:iend), source=nan)
  end subroutine init_gridcell_wstate_type
  !
  ! Initialize gridcell energy state variables
  !
  subroutine init_gridcell_estate_type(ibeg,iend,ges)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_estate_type), intent(inout) :: ges

    allocate(ges%gc_heat1(ibeg:iend), source=nan)
    allocate(ges%gc_heat2(ibeg:iend), source=nan)
  end subroutine init_gridcell_estate_type
#ifdef LCH4
  !
  ! Initialize gridcell ch4 variables
  !
  subroutine init_gridcell_ch4_type(ibeg,iend,gch4)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (gridcell_ch4_type), intent(inout) :: gch4

    allocate(gch4%c_atm(ibeg:iend,1:ngases), source=nan)
    allocate(gch4%ch4co2f(ibeg:iend), source=nan)
    allocate(gch4%ch4prodg(ibeg:iend), source=nan)
    allocate(gch4%nem(ibeg:iend), source=nan)
  end subroutine init_gridcell_ch4_type
#endif

  end subroutine initClmtype

end module mod_clm_typeinit
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
