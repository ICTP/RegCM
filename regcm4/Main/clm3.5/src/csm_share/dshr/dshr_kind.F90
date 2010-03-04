!===============================================================================
! SVN $Id: dshr_kind.F90 1048 2006-05-25 22:33:10Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_kind.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_kind -- F90 kind declarations
!
! !DESCRIPTION:
!   F90 kind declarations.
!
! !REVISION HISTORY:
!     2004-Dec-10 - B. Kauffman - created initial version
!
! !REMARKS:
!   This module does not use the standard dshr module variable naming convention
!   because this would result in excessively long variable declarations.
!   ie. we want to see real(R8) and not real(dshr_kind_r8)
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_kind

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
  integer,parameter,public :: CS = SHR_KIND_CS  ! short      char string
  integer,parameter,public :: CL = SHR_KIND_CL  ! long       char string
  integer,parameter,public :: CX = SHR_KIND_CX  ! extra-long char string

!EOP

!===============================================================================
end module dshr_kind
