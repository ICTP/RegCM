!===============================================================================
! SVN $Id$
! SVN $URL$
!===============================================================================

!BOP ===========================================================================
!
! !MODULE: eshr_rc_mod --- Functions to deal with the return codes from ESMF.
!
! !DESCRIPTION:
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2006-Aug-24 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

module eshr_rc_mod

! !USES:
   use ESMF_Mod

   implicit none

   private    ! default private

! ! PUBLIC TYPES:

  ! None

! ! PUBLIC MEMBER FUNCTIONS

  public :: eshr_rc_set   ! Set conditions for aborting on error
  public :: eshr_rc_check ! Check the return code and abort if not success

! !PUBLIC DATA MEMBERS:


  ! None


!EOP

   ! Private member functions:

   ! None

   ! Private data members:

   logical, save :: eshr_rc_continueOnFailure = .false.  ! Continue if hit failure

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_rc_set -- Set conditions for aborting
!   
! !DESCRIPTION:
!   
!  Set the conditions for aborting when checking the return code.
!  By default if this routine is NOT set will abort if rc /= ESMF_SUCCESS.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_rc_set( continueOnFailure )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    logical, intent(IN) :: continueOnFailure ! Continue even on ESMF failure
!EOP
    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    eshr_rc_continueOnFailure = continueOnFailure
    if ( eshr_rc_continueOnFailure )then
       write(6,'(" (eshr_rc_set) Continuing even on ESMF failure" )') 
    end if

END SUBROUTINE eshr_rc_set

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_rc_check -- Check the return code from ESMF
!   
! !DESCRIPTION:
!   
!  Check the return code from an ESMF function or subroutine call and abort
!  and give the input message if it was NOT succesful.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_rc_check( rc, msg )
    use shr_sys_mod, only: shr_sys_flush, shr_sys_abort
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    integer,          intent(IN)           :: rc  ! Return code from ESMF routine
    character(len=*), intent(IN), optional :: msg ! Error message to write out
!EOP
    !----- local -----
   character(len=*), parameter :: FA  = "( 'eshr_rc_check', A)"
   character(len=*), parameter :: FI  = "( 'eshr_rc_check', I5)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
   if ( rc == ESMF_SUCCESS ) return
   if ( present(msg) ) write(6,FA) trim(msg)
   call shr_sys_flush(6)
   if ( .not. eshr_rc_continueOnFailure )then
      call shr_sys_abort('eshr_rc_check')
   else
       write(6,FI) rc
   end if

END SUBROUTINE eshr_rc_check

!===============================================================================

end module eshr_rc_mod
