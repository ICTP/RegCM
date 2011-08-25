!===============================================================================
! SVN $Id: cpl_const_mod.F90 2404 2006-11-07 22:05:45Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_const_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_const_mod - defines/provides common constants.
!
! !DESCRIPTION:
!    Defines/provides common constants.
!
! !REVISION HISTORY:
!     2002-jun-10 - B. Kauffman - created module
!     2002-dec-5  - T. Craig    - names consistent with convention, cpl_const_*
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_const_mod

! !USES:

   use cpl_kind_mod   ! kinds
   use shr_const_mod  ! shared physical constants

   implicit none

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

  ! none

! !PUBLIC DATA MEMBERS:

   public

   !----------------------------------------------------------------------------
   ! physical constants
   !----------------------------------------------------------------------------
   real(R8),parameter :: cpl_const_pi      = SHR_CONST_PI     ! pi
   real(R8),parameter :: cpl_const_rearth  = SHR_CONST_REARTH ! radius of earth ~ m
   real(R8),parameter :: cpl_const_rearth2 = SHR_CONST_REARTH*SHR_CONST_REARTH ! rad**2
   real(R8),parameter :: cpl_const_deg2rad = cpl_const_pi/180.0_R8  ! deg to rads
   real(R8),parameter :: cpl_const_rad2deg = 180.0_R8/cpl_const_pi  ! rad to degs

   real(R8),parameter :: cpl_const_latice  = SHR_CONST_LATICE       ! latent heat of fusion
   real(R8),parameter :: cpl_const_ocn_ref_sal = SHR_CONST_OCN_REF_SAL ! ocn ref salt
   real(R8),parameter :: cpl_const_ice_ref_sal = SHR_CONST_ICE_REF_SAL ! ice ref salt

   real(R8),parameter :: cpl_const_spval       = SHR_CONST_SPVAL       ! special value

!EOP

end module cpl_const_mod
