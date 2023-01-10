MODULE mo_kind

  ! L. Kornblueh, MPI, August 2001, added working precision and comments 

  !FAB : 
!  use  mod_realkinds

  IMPLICIT NONE

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.

  ! Floating point section: 

  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: pi4 = 9
  INTEGER, PARAMETER :: pi8 = 14

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)
! Integer section

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)

CONTAINS

END MODULE mo_kind
