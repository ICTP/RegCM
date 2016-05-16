
module mod_clm_cndvlight

#if (defined CNDV)
  !
  ! Calculate light competition
  ! Update fpc for establishment routine
  ! Called once per year
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants , only : mathpi

  implicit none

  private

  save

  public :: Light

  contains
  !
  ! Calculate light competition
  ! Update fpc for establishment routine
  ! Called once per year
  subroutine Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)
    use mod_clm_type
    implicit none
    integer(ik4), intent(in) :: lbg, ubg     ! gridcell bounds
    integer(ik4), intent(in) :: lbp, ubp     ! pft bounds
    ! number of naturally-vegetated pfts in filter
    integer(ik4), intent(in) :: num_natvegp
    ! pft filter for naturally-vegetated points
    integer(ik4), intent(in) :: filter_natvegp(ubp-lbp+1)

    integer(ik4) , pointer :: ivt(:)       ! pft vegetation type
    integer(ik4) , pointer :: pgridcell(:) ! gridcell index of corresponding pft
    integer(ik4) , pointer :: tree(:)      ! ecophys const - tree pft or not
    !specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rk8), pointer :: slatop(:)
    real(rk8), pointer :: dsladlai(:) !dSLA/dLAI, projected area basis [m^2/gC]
    real(rk8), pointer :: woody(:)    ! ecophys const - woody pft or not
    real(rk8), pointer :: leafcmax(:) ! (gC/m2) leaf C storage
    real(rk8), pointer :: deadstemc(:) ! (gC/m2) dead stem C
    real(rk8), pointer :: dwood(:)     ! ecophys const - wood density (gC/m3)
    real(rk8), pointer :: reinickerp(:) ! ecophys const - parameter in allomet
    real(rk8), pointer :: crownarea_max(:) ! ecophys const - tree max crown a
    real(rk8), pointer :: allom1(:)        ! ecophys const - param in allomet

    ! area that each individual tree takes up (m^2)
    real(rk8), pointer :: crownarea(:)
    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer :: fpcgrid(:)
    real(rk8), pointer :: nind(:) ! number of individuals

    real(rk8), parameter :: fpc_tree_max = 0.95_rk8 !maximum total tree FPC
    integer(ik4)  :: p,fp, g   ! indices
    real(rk8) :: fpc_tree_total(lbg:ubg)
    real(rk8) :: fpc_inc_tree(lbg:ubg)
    real(rk8) :: fpc_inc(lbp:ubp) ! foliar projective cover increment (fraction)
    real(rk8) :: fpc_grass_total(lbg:ubg)
    real(rk8) :: fpc_shrub_total(lbg:ubg)
    real(rk8) :: fpc_grass_max(lbg:ubg)
    real(rk8) :: fpc_shrub_max(lbg:ubg)
    integer(ik4)  :: numtrees(lbg:ubg)
    real(rk8) :: excess
    real(rk8) :: nind_kill
    real(rk8) :: lai_ind
    real(rk8) :: fpc_ind
    real(rk8) :: fpcgrid_old
    real(rk8) :: lm_ind           !leaf carbon (gC/individual)
    real(rk8) :: stemdiam         ! stem diameter
    real(rk8) :: stocking         ! #stems / ha (stocking density)
    real(rk8) :: taper  ! ratio of height:radius_breast_height (tree allometry)

    ! Assign local pointers to derived type scalar members

    ivt           => clm3%g%l%c%p%itype
    pgridcell     => clm3%g%l%c%p%gridcell
    nind          => clm3%g%l%c%p%pdgvs%nind
    fpcgrid       => clm3%g%l%c%p%pdgvs%fpcgrid
    leafcmax      => clm3%g%l%c%p%pcs%leafcmax
    deadstemc     => clm3%g%l%c%p%pcs%deadstemc
    crownarea     => clm3%g%l%c%p%pdgvs%crownarea
    crownarea_max => dgv_pftcon%crownarea_max
    reinickerp    => dgv_pftcon%reinickerp
    allom1        => dgv_pftcon%allom1
    dwood         => pftcon%dwood
    slatop        => pftcon%slatop
    dsladlai      => pftcon%dsladlai
    woody         => pftcon%woody
    tree          => pftcon%tree

    taper = 200._rk8 ! make a global constant; used in Establishment + ?

    ! Initialize gridcell-level metrics

    do g = lbg, ubg
      fpc_tree_total(g) = 0._rk8
      fpc_inc_tree(g) = 0._rk8
      fpc_grass_total(g) = 0._rk8
      fpc_shrub_total(g) = 0._rk8
      numtrees(g) = 0
    end do

    do fp = 1,num_natvegp
      p = filter_natvegp(fp)
      g = pgridcell(p)

      ! Update LAI and FPC as in the last lines of DGVMAllocation

      if ( woody(ivt(p)) == 1._rk8 ) then
        if (fpcgrid(p) > 0._rk8 .and. nind(p) > 0._rk8) then
          !#ind/m2 nat veg area -> #ind/m2 pft area
          stocking = nind(p)/fpcgrid(p)
          ! stemdiam derived here from cn's formula for htop found in
          ! CNVegStructUpdate and cn's assumption stemdiam=2*htop/taper
          ! this derivation neglects upper htop limit enforced elsewhere
          stemdiam = (24._rk8 * deadstemc(p) / &
                  (mathpi * stocking * dwood(ivt(p)) * taper))**(1._rk8/3._rk8)
        else
          stemdiam = 0._rk8
        end if
        crownarea(p) = min(crownarea_max(ivt(p)), allom1(ivt(p)) * &
                stemdiam**reinickerp(ivt(p))) ! Eqn D (from Establishment)
!     else ! crownarea is 1 and does not need updating
      end if

      if ( crownarea(p) > 0._rk8 .and. nind(p) > 0._rk8 ) then
        lm_ind  = leafcmax(p) * fpcgrid(p) / nind(p)
        if ( dsladlai(ivt(p) ) > 0._rk8) then
          lai_ind = max(0.001_rk8,((exp(lm_ind*dsladlai(ivt(p)) + &
                  log(slatop(ivt(p)))) - &
                  slatop(ivt(p)))/dsladlai(ivt(p))) / crownarea(p))
        else
          lai_ind = lm_ind * slatop(ivt(p)) / crownarea(p)
        end if
      else
        lai_ind = 0._rk8
      end if

      fpc_ind = 1._rk8 - exp(-0.5_rk8*lai_ind)
      fpcgrid_old = fpcgrid(p)
      fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind
      fpc_inc(p) = max(0._rk8, fpcgrid(p) - fpcgrid_old)

      if ( woody(ivt(p)) == 1._rk8 ) then
        if ( tree(ivt(p)) == 1 ) then
          numtrees(g) = numtrees(g) + 1
          fpc_tree_total(g) = fpc_tree_total(g) + fpcgrid(p)
          fpc_inc_tree(g) = fpc_inc_tree(g) + fpc_inc(p)
        else ! if shrubs
          fpc_shrub_total(g) = fpc_shrub_total(g) + fpcgrid(p)
        end if
      else    ! if grass
        fpc_grass_total(g) = fpc_grass_total(g) + fpcgrid(p)
      end if
    end do

    do g = lbg, ubg
      fpc_grass_max(g) = 1._rk8 - min(fpc_tree_total(g), fpc_tree_max)
      fpc_shrub_max(g) = max(0._rk8, fpc_grass_max(g) - fpc_grass_total(g))
    end do

    ! The gridcell level metrics are now in place; continue...
    ! slevis replaced the previous code that updated pfpcgrid
    ! with a simpler way of doing so:
    ! fpcgrid(p) = fpcgrid(p) - excess
    ! Later we may wish to update this subroutine
    ! according to Strassmann's recommendations (see relevant pdf)

    do fp = 1,num_natvegp
      p = filter_natvegp(fp)
      g = pgridcell(p)

      ! light competition

      if ( woody(ivt(p)) == 1._rk8 .and. tree(ivt(p)) == 1._rk8 ) then
        if ( fpc_tree_total(g) > fpc_tree_max ) then
          if ( fpc_inc_tree(g) > 0._rk8 ) then
            excess = (fpc_tree_total(g) - fpc_tree_max) * &
                     fpc_inc(p) / fpc_inc_tree(g)
          else
            excess = (fpc_tree_total(g) - fpc_tree_max) / &
                     dble(numtrees(g))
          end if

          ! Reduce individual density (and thereby gridcell-level biomass)
          ! so that total tree FPC reduced to 'fpc_tree_max'

          if ( fpcgrid(p) > 0._rk8 ) then
            nind_kill = nind(p) * excess / fpcgrid(p)
            nind(p) = max(0._rk8, nind(p) - nind_kill)
            fpcgrid(p) = max(0._rk8, fpcgrid(p) - excess)
          else
            nind(p) = 0._rk8
            fpcgrid(p) = 0._rk8
          end if

          ! Transfer lost biomass to litter

        end if ! if tree cover exceeds max allowed
      else if ( woody(ivt(p)) == 0._rk8 ) then ! grass
        if ( fpc_grass_total(g) > fpc_grass_max(g) ) then
          ! grass competes with itself if total fpc exceeds 1
          excess = (fpc_grass_total(g) - fpc_grass_max(g)) * &
                  fpcgrid(p) / fpc_grass_total(g)
          fpcgrid(p) = max(0._rk8, fpcgrid(p) - excess)
        end if
      else if ( woody(ivt(p)) == 1._rk8 .and. tree(ivt(p)) == 0._rk8 ) then ! shrub
        if ( fpc_shrub_total(g) > fpc_shrub_max(g) ) then
          excess = 1._rk8 - fpc_shrub_max(g) / fpc_shrub_total(g)
          ! Reduce individual density (and thereby gridcell-level biomass)
          ! so that total shrub FPC reduced to fpc_shrub_max(g)
          if ( fpcgrid(p) > 0._rk8 ) then
            nind_kill = nind(p) * excess / fpcgrid(p)
            nind(p) = max(0._rk8, nind(p) - nind_kill)
            fpcgrid(p) = max(0._rk8, fpcgrid(p) - excess)
          else
            nind(p) = 0._rk8
            fpcgrid(p) = 0._rk8
          end if
        end if
      end if   ! end of if-tree
    end do
  end subroutine Light
#endif

end module mod_clm_cndvlight
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
