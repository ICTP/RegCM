module mod_clm_fractwet

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FracWetMod
!
! !DESCRIPTION:
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FracWet
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FracWet
!
! !INTERFACE:
  subroutine FracWet(numf, filter)
!
! !DESCRIPTION:
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry. The variable ``fwet'' is the
! fraction of all vegetation surfaces which are wet including
! stem area which contribute to evaporation. The variable ``fdry''
! is the fraction of elai which is dry because only leaves
! can transpire.  Adjusted for stem area which does not transpire.
!
! !USES:
    use mod_realkinds
    use mod_clm_type
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: numf                  ! number of filter non-lake points
    integer, intent(in) :: filter(numf)          ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Hydrology1 in module Hydrology1Mod
!
! !REVISION HISTORY:
! Created by Keith Oleson and M. Vertenstein
! 03/08/29 Mariana Vertenstein : Migrated to vectorized code
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: frac_veg_nosno(:) ! fraction of veg not covered by snow (0/1 now) [-]
    real(rk8), pointer :: dewmx(:)          ! Maximum allowed dew [mm]
    real(rk8), pointer :: elai(:)           ! one-sided leaf area index with burying by snow
    real(rk8), pointer :: esai(:)           ! one-sided stem area index with burying by snow
    real(rk8), pointer :: h2ocan(:)         ! total canopy water (mm H2O)
!
! local pointers to implicit out arguments
!
    real(rk8), pointer :: fwet(:)           ! fraction of canopy that is wet (0 to 1)
    real(rk8), pointer :: fdry(:)           ! fraction of foliage that is green and dry [-] (new)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: fp,p             ! indices
    real(rk8) :: vegt             ! frac_veg_nosno*lsai
    real(rk8) :: dewmxi           ! inverse of maximum allowed dew [1/mm]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (pft-level)

    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    dewmx => clm3%g%l%c%p%pps%dewmx
    elai => clm3%g%l%c%p%pps%elai
    esai => clm3%g%l%c%p%pps%esai
    h2ocan => clm3%g%l%c%p%pws%h2ocan
    fwet => clm3%g%l%c%p%pps%fwet
    fdry => clm3%g%l%c%p%pps%fdry

    ! Compute fraction of canopy that is wet and dry

    do fp = 1,numf
       p = filter(fp)
       if (frac_veg_nosno(p) == 1) then
          if (h2ocan(p) > 0.D0) then
             vegt    = frac_veg_nosno(p)*(elai(p) + esai(p))
             dewmxi  = 1.0D0/dewmx(p)
             fwet(p) = ((dewmxi/vegt)*h2ocan(p))**0.666666666666D0
             fwet(p) = min (fwet(p),1.0D0)   ! Check for maximum limit of fwet
          else
             fwet(p) = 0.D0
          end if
          fdry(p) = (1.D0-fwet(p))*elai(p)/(elai(p)+esai(p))
#if (defined PERGRO)
          fwet(p) = 0.D0
          fdry(p) = elai(p)/(elai(p)+esai(p))
#endif
       else
          fwet(p) = 0.D0
          fdry(p) = 0.D0
       end if
    end do

  end subroutine FracWet

end module mod_clm_fractwet
