
module mod_clm_cndvecosystemdynini

#if (defined CNDV)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDVEcosystemDyniniMod
!
! !DESCRIPTION:
!
! !USES:
  use mod_realkinds
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: CNDVEcosystemDynini ! CNDV related initializations
!
! !REVISION HISTORY:
! Created by Sam Levis following DGVMEcosystemDynMod by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNDVEcosystemDynini
!
! !INTERFACE:
  subroutine CNDVEcosystemDynini()
!
! !DESCRIPTION:
! CNDV related initializations
!
! !USES:
    use mod_clm_type
    use mod_clm_decomp , only : get_proc_bounds, get_proc_global
    use mod_constants , only : tzero
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from LPJ initialization subroutines)
!         Sam Levis (adapted for CNDV coupling; eliminated redunant parameters)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p,n           ! indices
    integer  :: begp, endp      ! per-proc beginning and ending pft indices
    integer  :: begc, endc      !              "                column indices
    integer  :: begl, endl      !              "                landunit indices
    integer  :: begg, endg      !              "                gridcell indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    pptr => clm3%g%l%c%p

    ! ---------------------------------------------------------------
    ! Some of the following came from LPJ subroutine initgrid
    ! ---------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do p = begp,endp
       pptr%pdgvs%present(p)   = .false.
       pptr%pdgvs%crownarea(p) = 0.D0
       pptr%pdgvs%nind(p)      = 0.D0
       pptr%pcs%leafcmax(p)    = 0.D0
       pptr%pdgvs%t_mo_min(p)  = 1.0D+36
    end do

    do g = begg,endg
       gptr%gdgvs%agdd20(g)   = 0.D0
       gptr%gdgvs%tmomin20(g) = tzero - 5.D0 !initialize this way for Phenology code
    end do

  end subroutine CNDVEcosystemDynini

#endif

end module mod_clm_cndvecosystemdynini
