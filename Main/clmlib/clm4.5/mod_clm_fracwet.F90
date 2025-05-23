module mod_clm_fracwet
  !
  ! Determine fraction of vegetated surfaces which are wet and
  ! fraction of elai which is dry.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type

  implicit none

  private

  save

  public :: FracWet

  contains
  !
  ! Determine fraction of vegetated surfaces which are wet and
  ! fraction of elai which is dry. The variable ``fwet'' is the
  ! fraction of all vegetation surfaces which are wet including
  ! stem area which contribute to evaporation. The variable ``fdry''
  ! is the fraction of elai which is dry because only leaves
  ! can transpire.  Adjusted for stem area which does not transpire.
  !
  subroutine FracWet(numf, filter)
    implicit none
    ! number of filter non-lake points
    integer(ik4), intent(in) :: numf
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter(numf)
    ! fraction of veg not covered by snow (0/1 now) [-]
    integer(ik4), pointer, contiguous :: frac_veg_nosno(:)
    ! Maximum allowed dew [mm]
    real(rk8), pointer, contiguous :: dewmx(:)
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer, contiguous :: elai(:)
    ! one-sided stem area index with burying by snow
    real(rk8), pointer, contiguous :: esai(:)
    ! total canopy water (mm H2O)
    real(rk8), pointer, contiguous :: h2ocan(:)
    ! fraction of canopy that is wet (0 to 1)
    real(rk8), pointer, contiguous :: fwet(:)
    ! fraction of foliage that is green and dry [-] (new)
    real(rk8), pointer, contiguous :: fdry(:)
    integer(ik4) :: fp, p   ! indices
    real(rk8) :: vegt        ! frac_veg_nosno*lsai
    real(rk8) :: dewmxi      ! inverse of maximum allowed dew [1/mm]

    ! Assign local pointers to derived subtypes components (pft-level)

    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    dewmx => clm3%g%l%c%p%pps%dewmx
    elai => clm3%g%l%c%p%pps%elai
    esai => clm3%g%l%c%p%pps%esai
    h2ocan => clm3%g%l%c%p%pws%h2ocan
    fwet => clm3%g%l%c%p%pps%fwet
    fdry => clm3%g%l%c%p%pps%fdry

    ! Compute fraction of canopy that is wet and dry

    do fp = 1, numf
      p = filter(fp)
      if ( frac_veg_nosno(p) == 1 ) then
        if ( h2ocan(p) > 0._rk8 ) then
          vegt    = frac_veg_nosno(p)*(elai(p) + esai(p))
          dewmxi  = 1.0_rk8/dewmx(p)
          fwet(p) = ((dewmxi/vegt)*h2ocan(p))**0.666666666666_rk8
          fwet(p) = min (fwet(p),1.0_rk8)   ! Check for maximum limit of fwet
        else
          fwet(p) = 0._rk8
        end if
        fdry(p) = (1._rk8-fwet(p))*elai(p)/(elai(p)+esai(p))
#if (defined PERGRO)
        fwet(p) = 0._rk8
        fdry(p) = elai(p)/(elai(p)+esai(p))
#endif
      else
        fwet(p) = 0._rk8
        fdry(p) = 0._rk8
      end if
    end do
  end subroutine FracWet

end module mod_clm_fracwet

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
