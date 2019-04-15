!===============================================================================
! SVN $Id: cpl_kind_mod.F90 238 2006-02-08 18:13:46Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_kind_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_kind_mod -- F90 kind declarations
!
! !DESCRIPTION:
!   F90 kind declarations.
!
! !REVISION HISTORY:
!     2002-Nov-04 - B. Kauffman - created initial version
!
! !REMARKS:
!   This module does not use the standard cpl6 module variable naming convention
!   because this would results in excessively long variable declarations.
!   ie. we want to see real(R8) and not real(cpl_kind_r8)
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_kind_mod

! !USES:

   use shr_kind_mod  !  shared kind declaration

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
  ! none

! !PUBLIC MEMBER FUNCTIONS:

  ! none

! !PUBLIC DATA MEMBERS:

  integer,parameter,public :: R8 = SHR_KIND_R8  ! 8 byte real
  integer,parameter,public :: R4 = SHR_KIND_R4  ! 4 byte real
  integer,parameter,public :: RN = SHR_KIND_RN  ! native/default real
  integer,parameter,public :: I8 = SHR_KIND_I8  ! 8 byte integer
  integer,parameter,public :: I4 = SHR_KIND_I4  ! 4 byte integer
  integer,parameter,public :: IN = SHR_KIND_IN  ! native/default integer

  integer,parameter,public :: CL = SHR_KIND_CL  ! generic "long" char string

!EOP

end module cpl_kind_mod
