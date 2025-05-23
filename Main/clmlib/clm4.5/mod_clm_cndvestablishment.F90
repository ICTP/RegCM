
module mod_clm_cndvestablishment

#if (defined CNDV)
  !
  ! Calculates establishment of new pfts
  ! Called once per year
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio

  implicit none

  private

  save

  public :: Establishment

  contains
  !
  ! Calculates establishment of new pfts
  ! Called once per year
  !
  subroutine Establishment(lbg, ubg, lbp, ubp)
    use mod_clm_type
    use mod_clm_varpar  , only : numpft
    use mod_clm_varcon  , only : istsoil
    use mod_clm_pftvarcon   , only : noveg, nc3_arctic_grass
    use mod_clm_varcon, only : secspday, rpi, tfrz
    implicit none
    integer(ik4), intent(in) :: lbg, ubg         ! gridcell bounds
    integer(ik4), intent(in) :: lbp, ubp         ! pft bounds

    ! exclude seasonal decid pfts from tropics [1=true, 0=false]
    logical, pointer, contiguous :: pftmayexist(:)
    integer(ik4), pointer, contiguous :: plandunit(:)  ! landunit of corresponding pft
    integer(ik4), pointer, contiguous :: pgridcell(:)  ! gridcell of corresponding pft
    integer(ik4), pointer, contiguous :: ltype(:) ! landunit type for corresponding pft
    real(rk8), pointer, contiguous :: tmomin20(:)  ! 20-yr running mean of tmomin
    real(rk8), pointer, contiguous :: agdd20(:)    ! 20-yr running mean of agdd
    ! accumulated growing degree days above twmax
    real(rk8), pointer, contiguous :: agddtw(:)
    ! 365-day running mean of tot. precipitation
    real(rk8), pointer, contiguous :: prec365(:)
    !specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous :: slatop(:)
    !dSLA/dLAI, projected area basis [m^2/gC]
    real(rk8), pointer, contiguous :: dsladlai(:)
    ! ecophys const - woody pft or not
    real(rk8), pointer, contiguous :: woody(:)
    ! ecophys const - tree maximum crown area [m2]
    real(rk8), pointer, contiguous :: crownarea_max(:)
    ! ecophys const - upper limit of temperature of the warmest month
    real(rk8), pointer, contiguous :: twmax(:)
    ! ecophys const - parameter in allometric equation
    real(rk8), pointer, contiguous :: reinickerp(:)
    ! ecophys const - wood density (gC/m3)
    real(rk8), pointer, contiguous :: dwood(:)
    ! ecophys const - parameter in allometric
    real(rk8), pointer, contiguous :: allom1(:)
    ! ecophys const - minimum coldest monthly mean temperature
    real(rk8), pointer, contiguous :: tcmin(:)
    ! ecophys const - maximum coldest monthly mean temperature
    real(rk8), pointer, contiguous :: tcmax(:)
    ! ecophys const - minimum growing degree days (at or above 5 C)
    real(rk8), pointer, contiguous :: gddmin(:)
    real(rk8), pointer, contiguous :: leafcmax(:)        ! (gC/m2) ann max leaf C
    real(rk8), pointer, contiguous :: deadstemc(:)       ! (gC/m2) dead stem C
    real(rk8), pointer, contiguous :: annsum_npp(:)      ! annual sum NPP (gC/m2/yr)
    real(rk8), pointer, contiguous :: annsum_litfall(:)  ! annual sum litfall (gC/m2/yr)

    integer(ik4), pointer, contiguous :: ivt(:)  ! vegetation type for this pft
    logical, pointer, contiguous :: present(:)   ! true=> PFT present in patch
    real(rk8), pointer, contiguous :: nind(:)     ! number of individuals (#/m**2)

    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer, contiguous :: fpcgrid(:)
    ! area that each individual tree takes up (m^2)
    real(rk8), pointer, contiguous :: crownarea(:)
    ! lpj's growth efficiency
    real(rk8), pointer, contiguous :: greffic(:)
    real(rk8), pointer, contiguous :: heatstress(:)

    integer(ik4)  :: g,l,p,m   ! indices
    ! local gridcell filter for error check
    integer(ik4)  :: fn, filterg(ubg-lbg+1)
    !
    ! gridcell level variables
    !
    integer(ik4)  :: ngrass(lbg:ubg)      ! counter
    integer(ik4)  :: npft_estab(lbg:ubg)  ! counter
    ! total fractional cover of trees in vegetated portion of gridcell
    real(rk8) :: fpc_tree_total(lbg:ubg)
    ! old-total fractional vegetated portion of gridcell (without bare ground)
    real(rk8) :: fpc_total(lbg:ubg)
    ! new-total fractional vegetated portion of gridcell (without bare ground)
    real(rk8) :: fpc_total_new(lbg:ubg)
    !
    ! pft level variables
    !
    logical  :: survive(lbp:ubp)    ! true=>pft survives
    logical  :: estab(lbp:ubp)      ! true=>pft is established
    logical  :: iswood(lbp:ubp)     ! true=>pft is woody
    real(rk8) :: dstemc(lbp:ubp)    ! local copy of deadstemc

    real(rk8) :: taper ! ratio of height:radius_breast_height (tree allometry)
    real(rk8) :: estab_rate        !establishment rate
    real(rk8) :: estab_grid        !establishment rate on grid cell
    real(rk8) :: fpcgridtemp       ! temporary
    real(rk8) :: stemdiam          ! stem diameter
    real(rk8) :: stocking          ! #stems / ha (stocking density)
    real(rk8) :: lai_ind           ! LAI per individual
    real(rk8) :: lm_ind            !leaf carbon (gC/ind)
    real(rk8) :: fpc_ind           !individual foliage projective cover
    real(rk8):: bm_delta

    real(rk8), parameter :: ramp_agddtw = 300.0
    !
    ! minimum individual density for persistence of PFT (indiv/m2)
    !
    real(rk8), parameter :: nind_min = 1.0e-10_rk8
    !
    ! minimum precip. for establishment (mm/s)
    !
    real(rk8), parameter :: prec_min_estab = 100.0_rk8/(365.0_rk8*secspday)
    !
    ! maximum sapling establishment rate (indiv/m2)
    !
    real(rk8), parameter :: estab_max = 0.24_rk8

    ! Assign local pointers to derived type members (gridcell-level)

    agdd20        => clm3%g%gdgvs%agdd20
    tmomin20      => clm3%g%gdgvs%tmomin20

    ! Assign local pointers to derived type members (landunit-level)

    ltype => clm3%g%l%itype

    ! Assign local pointers to derived type members (pft-level)

    ivt           => clm3%g%l%c%p%itype
    pgridcell     => clm3%g%l%c%p%gridcell
    plandunit     => clm3%g%l%c%p%landunit
    present       => clm3%g%l%c%p%pdgvs%present
    nind          => clm3%g%l%c%p%pdgvs%nind
    fpcgrid       => clm3%g%l%c%p%pdgvs%fpcgrid
    crownarea     => clm3%g%l%c%p%pdgvs%crownarea
    greffic       => clm3%g%l%c%p%pdgvs%greffic
    heatstress    => clm3%g%l%c%p%pdgvs%heatstress
    annsum_npp     => clm3%g%l%c%p%pepv%annsum_npp
    annsum_litfall => clm3%g%l%c%p%pepv%annsum_litfall
    prec365       => clm3%g%l%c%p%pdgvs%prec365
    agddtw        => clm3%g%l%c%p%pdgvs%agddtw
    pftmayexist   => clm3%g%l%c%p%pdgvs%pftmayexist

    ! Assign local pointers to derived type members (vegetation types)

    crownarea_max => dgv_pftcon%crownarea_max
    twmax         => dgv_pftcon%twmax
    reinickerp    => dgv_pftcon%reinickerp
    allom1        => dgv_pftcon%allom1
    tcmax         => dgv_pftcon%tcmax
    tcmin         => dgv_pftcon%tcmin
    gddmin        => dgv_pftcon%gddmin
    leafcmax      => clm3%g%l%c%p%pcs%leafcmax
    deadstemc     => clm3%g%l%c%p%pcs%deadstemc
    slatop        => pftcon%slatop
    dsladlai      => pftcon%dsladlai
    dwood         => pftcon%dwood
    woody         => pftcon%woody

    ! **********************************************************************
    ! Slevis version of LPJ's subr. bioclim
    ! Limits based on 20-year running averages of coldest-month mean
    ! temperature and growing degree days (5 degree base).
    ! For SURVIVAL, coldest month temperature and GDD should be
    ! at least as high as PFT-specific limits.
    ! For REGENERATION, PFT must be able to survive AND coldest month
    ! temperature should be no higher than a PFT-specific limit.
    ! **********************************************************************

    taper = 200.0_rk8 ! make a global constant as with dwood (lpj's wooddens)

    ! Initialize gridcell-level metrics

    do g = lbg, ubg
      ngrass(g) = 0
      npft_estab(g) = 0
      fpc_tree_total(g) = 0.0_rk8
      fpc_total(g) = 0.0_rk8
      fpc_total_new(g) = 0.0_rk8
    end do

    do p = lbp, ubp
      g = pgridcell(p)

      ! Set the presence of pft for this gridcell

      if ( nind(p) < 1.0e-10_rk8 ) present(p) = .false.
      if (.not. present(p)) then
        nind(p) = 0.0_rk8
        fpcgrid(p) = 0.0_rk8
      end if
      survive(p) = .false.
      estab(p)   = .false.
      dstemc(p)  = deadstemc(p)
      iswood(p)  = abs(woody(ivt(p))-1.0_rk8) < epsilon(1.0)
      if ( fpcgrid(p) < 1.0e-10_rk8 ) then
        nind(p) = 0.0_rk8
        fpcgrid(p) = 0.0_rk8
        present(p) = .false.
      end if
    end do

    ! Must go thru all 16 pfts and decide which can/cannot establish or survive
    ! Determine present, survive, estab.  Note: Even if tmomin20>tcmax, crops
    ! and 2nd boreal summergreen tree cannot exist (see
    ! EcosystemDynini) because this model cannot simulate such pfts, yet.
    ! Note - agddtw is only defined at the pft level and has now been moved
    ! to an if-statement below to determine establishment of boreal trees

    do p = lbp, ubp
      g = pgridcell(p)
      if (tmomin20(g) >= tcmin(ivt(p)) + tfrz ) then
        if ( tmomin20(g) <= tcmax(ivt(p)) + tfrz  .and. &
             agdd20(g) >= gddmin(ivt(p))) then
          estab(p) = .true.
        end if
        survive(p) = .true.
        ! seasonal decid. pfts that would have occurred in regions without
        ! short winter day lengths (see CNPhenology)
        if ( .not. pftmayexist(p) ) then
          survive(p) = .false.
          estab(p) = .false.
          pftmayexist(p) = .true.
        end if
      end if
    end do

    do p = lbp, ubp
      g = pgridcell(p)
      l = plandunit(p)

      ! Case 1 -- pft ceases to exist -kill pfts not adapted to current climate

      if ( present(p) .and. (.not. survive(p) .or. nind(p) < nind_min) ) then
        present(p) = .false.
        fpcgrid(p) = 0.0_rk8
        nind(p) = 0.0_rk8
      end if

      ! Case 2 -- pft begins to exist - introduce newly "adapted" pfts

      if ( ltype(l) == istsoil ) then
        if ( .not. present(p) .and. &
             prec365(p) >= prec_min_estab .and. estab(p) ) then
          if ( twmax(ivt(p)) > 999.0_rk8 .or. agddtw(p) == 0.0_rk8 ) then

            present(p) = .true.
            nind(p) = 0.0_rk8
            ! lpj starts with fpcgrid=0 and calculates
            ! seed fpcgrid from the carbon of saplings;
            ! with CN we need the seed fpcgrid up front
            ! to scale seed leafc to lm_ind to get fpcgrid;
            ! sounds circular; also seed fpcgrid depends on sla,
            ! so theoretically need diff value for each pft;slevis
            fpcgrid(p) = 0.000844_rk8
            if ( iswood(p) ) fpcgrid(p) = 0.05_rk8

            ! Seed carbon for newly established pfts
            ! Equiv. to pleaf=1 & pstor=1 set in subr pftwt_cnbal (slevis)
            ! ***Dangerous*** to hardwire leafcmax here; find alternative!
            ! Consider just assigning nind and fpcgrid for newly
            ! established pfts instead of entering the circular procedure
            ! outlined in the paragraph above
            leafcmax(p) = 1.0_rk8
            if ( dstemc(p) <= 0.0_rk8 ) dstemc(p) = 0.1_rk8
          end if   ! conditions required for establishment
        end if   ! conditions required for establishment
      end if   ! if soil

      ! Case 3 -- some pfts continue to exist (no change) and some pfts
      ! continue to not exist (no change). Do nothing for this case.

    end do

    ! Sapling and grass establishment
    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)
    ! Calculate total woody FPC and number of woody PFTs present and
    ! able to establish

    do p = lbp, ubp
      g = pgridcell(p)
      if ( present(p) ) then
        if ( iswood(p) ) then
          fpc_tree_total(g) = fpc_tree_total(g) + fpcgrid(p)
          if ( estab(p) ) npft_estab(g) = npft_estab(g) + 1
        else
          if ( ivt(p) > noveg ) then !grass
            ngrass(g) = ngrass(g) + 1
          end if
        end if
      end if
    end do

    ! Above grid-level establishment counters are required for the next steps.

    do p = lbp, ubp
      g = pgridcell(p)

      if ( present(p) .and. estab(p) .and. iswood(p) ) then

        ! Calculate establishment rate over available space, per tree PFT
        ! Max establishment rate reduced by shading as tree FPC approaches 1
        ! Total establishment rate partitioned equally among
        ! regenerating woody PFTs

        estab_rate = estab_max*(1.0_rk8-exp(5.0_rk8 * &
            (fpc_tree_total(g)-1.0_rk8))) / dble(npft_estab(g))

        ! Calculate grid-level establishment rate per woody PFT
        ! Space available for woody PFT establishment is fraction of grid cell
        ! not currently occupied by woody PFTs

        estab_grid = estab_rate * (1.0_rk8-fpc_tree_total(g))

        ! Add new saplings to current population

        nind(p) = nind(p) + estab_grid

        !slevis: lpj's lm_ind was the max leaf mass for the year;
        !now lm_ind is the max leaf mass for the year calculated in CNFire
        !except when a pft is newly established (nind==0); then lm_ind
        !is assigned a leafcmax above

        lm_ind = leafcmax(p) * fpcgrid(p) / nind(p) ! nind>0 for sure
        if ( fpcgrid(p) > 0.0_rk8 .and. nind(p) > 0.0_rk8 ) then
          !#ind/m2 nat veg area -> #ind/m2 pft area
          stocking = nind(p)/fpcgrid(p)
          ! stemdiam derived here from cn's formula for htop found in
          ! CNVegStructUpdate and cn's assumption stemdiam=2*htop/taper
          ! this derivation neglects upper htop limit enforced elsewhere
          stemdiam = (24.0_rk8 * dstemc(p) / &
                  (rpi * stocking * dwood(ivt(p)) * taper))**(1.0_rk8/3.0_rk8)
        else
          stemdiam = 0.0_rk8
        end if
        ! Eqn D (now also in Light; need here for 1st yr when pfts
        ! haven't established, yet)
        crownarea(p) = min(crownarea_max(ivt(p)), allom1(ivt(p)) * &
                stemdiam**reinickerp(ivt(p)))

        ! Update LAI and FPC

        if (crownarea(p) > 0.0_rk8) then
          if (dsladlai(ivt(p)) > 0.0_rk8) then
            ! make lai_ind >= 0.001 to avoid killing plants at this stage
            lai_ind = max(0.001_rk8,((exp(lm_ind*dsladlai(ivt(p)) + &
                    log(slatop(ivt(p)))) - slatop(ivt(p))) / &
                    dsladlai(ivt(p))) / crownarea(p))
          else ! currently redundant because dsladlai=0 for grasses only
            lai_ind = lm_ind * slatop(ivt(p)) / crownarea(p) ! lpj's formula
          end if
        else
          lai_ind = 0.0_rk8
        end if
        fpc_ind = 1.0_rk8 - exp(-0.5_rk8*lai_ind)
        fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind

      end if   ! add new saplings block

      if ( present(p) .and. iswood(p) ) then
        fpc_total_new(g) = fpc_total_new(g) + fpcgrid(p)
      end if

    end do   ! close loop to update fpc_total_new

    ! Adjustments- don't allow trees to exceed 95% of vegetated landunit

    do p = lbp, ubp
      g = pgridcell(p)
      if ( fpc_total_new(g) > 0.95_rk8 ) then
        if ( present(p) .and. iswood(p) ) then
          nind(p) = nind(p) * 0.95_rk8 / fpc_total_new(g)
          fpcgrid(p) = fpcgrid(p) * 0.95_rk8 / fpc_total_new(g)
        end if
        fpc_total(g) = 0.95_rk8
      else
        fpc_total(g) = fpc_total_new(g)
      end if
    end do

    ! Section for grasses. Grasses can establish in non-vegetated areas

    do p = lbp, ubp
      g = pgridcell(p)
      if ( present(p) .and. .not. iswood(p) ) then
        if ( leafcmax(p) <= 0.0_rk8 .or. fpcgrid(p) <= 0.0_rk8 ) then
          present(p) = .false.
          nind(p) = 0.0_rk8
        else
          nind(p) = 1.0_rk8 ! in case these grasses just established
          crownarea(p) = 1.0_rk8
          lm_ind = leafcmax(p) * fpcgrid(p) / nind(p)
          if (dsladlai(ivt(p)) > 0.0_rk8) then
            lai_ind = max(0.001_rk8,((exp(lm_ind*dsladlai(ivt(p)) + &
                    log(slatop(ivt(p)))) - slatop(ivt(p))) / &
                    dsladlai(ivt(p))) / crownarea(p))
          else ! 'if' is currently redundant b/c dsladlai=0 for grasses only
            lai_ind = lm_ind * slatop(ivt(p)) / crownarea(p)
          end if
          fpc_ind = 1.0_rk8 - exp(-0.5_rk8*lai_ind)
          fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind
          fpc_total(g) = fpc_total(g) + fpcgrid(p)
        end if
      end if
    end do   ! end of pft-loop

    ! Adjustment of fpc_total > 1 due to grasses (ivt >= nc3_arctic_grass)

    do p = lbp, ubp
      g = pgridcell(p)

      if ( fpc_total(g) > 1.0_rk8 ) then
        if ( ivt(p) >= nc3_arctic_grass .and. fpcgrid(p) > 0.0_rk8 ) then
          fpcgridtemp = fpcgrid(p)
          fpcgrid(p) = max(0.0_rk8, fpcgrid(p) - (fpc_total(g)-1.0_rk8))
          fpc_total(g) = fpc_total(g) - fpcgridtemp + fpcgrid(p)
        end if
      end if

      ! Remove tiny fpcgrid amounts

      if ( fpcgrid(p) < 1.e-15_rk8 ) then
        fpc_total(g) = fpc_total(g) - fpcgrid(p)
        fpcgrid(p) = 0.0_rk8
        present(p) = .false.
        nind(p) = 0.0_rk8
      end if

    end do

    ! Set the fpcgrid for bare ground if there is bare ground in
    ! vegetated landunit and pft is bare ground so that everything
    ! can add up to one.

    do p = lbp, ubp

      g = pgridcell(p)
      if ( ivt(p) == noveg ) then
        if ( fpc_total(g) < 1.0_rk8 ) then
          fpcgrid(p) = 1.0_rk8 - fpc_total(g)
        else
          fpcgrid(p) = 0.0_rk8
        end if
        fpc_total(g) = fpc_total(g) + fpcgrid(p)
      end if

    end do

    ! Annual calculations used hourly in GapMortality
    ! Ultimately may wish to place in separate subroutine...

    do p = lbp, ubp
      g = pgridcell(p)

      ! Stress mortality from lpj's subr Mortality

      if ( iswood(p) .and. nind(p) > 0.0_rk8 .and. &
           leafcmax(p) > 0.0_rk8 .and. fpcgrid(p) > 0.0_rk8 ) then
        if ( twmax(ivt(p)) < 999.0_rk8 ) then
          heatstress(p) = max(0.0_rk8, min(1.0_rk8, agddtw(p) / ramp_agddtw))
        else
          heatstress(p) = 0.0_rk8
        end if

        ! Net individual living biomass increment
        ! NB: lpj's turnover not exactly same as cn's litfall:
        ! lpj's sap->heartwood turnover not included in litfall (slevis)

        bm_delta = max(0.0_rk8, annsum_npp(p) - annsum_litfall(p))
        lm_ind = leafcmax(p) * fpcgrid(p) / nind(p)

        ! Growth efficiency (net biomass increment per unit leaf area)
        if ( dsladlai(ivt(p)) > 0.0_rk8 ) then
          greffic(p) = bm_delta / (max(0.001_rk8,                   &
             ( ( exp(lm_ind*dsladlai(ivt(p)) + log(slatop(ivt(p)))) &
                 - slatop(ivt(p)) ) / dsladlai(ivt(p)) )))
        else ! currently redundant because dsladlai=0 for grasses only
          greffic(p) = bm_delta / (lm_ind * slatop(ivt(p)))
        end if
      else
        greffic(p) = 0.0_rk8
        heatstress(p) = 0.0_rk8
      end if
    end do

    ! Check for error in establishment
    do p = lbp, ubp
      g = pgridcell(p)
      if ( is_nan(fpcgrid(p)) ) then
        write(stderr,*) 'Error in Establishment: fpcgrid =', &
                fpcgrid(p), ' at pft ',p
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do

    fn = 0
    do g = lbg, ubg
      if ( abs(fpc_total(g) - 1.0_rk8) > 1.e-6_rk8 ) then
        fn = fn + 1
        filterg(fn) = g
      end if
    end do
    ! Just print out the first error
    if ( fn > 0 ) then
      g = filterg(1)
      write(stderr,*) 'Error in Establishment: fpc_total =', &
              fpc_total(g), ' at gridcell ',g
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

  end subroutine Establishment

#endif

end module mod_clm_cndvestablishment
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
