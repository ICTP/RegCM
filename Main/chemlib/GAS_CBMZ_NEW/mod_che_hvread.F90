module mod_che_hvread
  use mod_cbmz_jval1

  contains

  subroutine hvread

! This reads the MADRONICH LOOKUP TABLE (2002 VERSION).
!   Input file:  TUVGRID2  (kept in dhvmad/TUVcode4.1a)
!    (see TUVINFO in /l/kudzu/k-1/sillman/dhvmad)
!                 (also tuvtab2.f, tested in tuvtest2.f)
!
!   OPTION:  Link to jval2.f for altitude in km (standard)
!            or to   jval1.f for altitude in kPa
!
! Result: creates the hv input arrays:
!     c_hvin,c_nhv,c_hvmat, c_hvmatb, c_jarray
!
! Called by:  boxmain (as part of chemical setup)
! Calls to:   readhv (in file jval2.f (KM) or jval1.f (kPA).)
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
    IMPLICIT NONE

    integer c_hvin
    integer c_nhv(22)
    real(kind=8) c_hvmat(22,40)
    real(kind=8) c_hvmatb(22)
    real(kind=8) c_jarray(80,510,56)
    common/HVRVARS/c_hvmat, c_hvmatb, c_jarray
    common/HVINT/ c_hvin,c_nhv
    c_hvin = 26

    call readhv(c_hvin,c_nhv,c_hvmat, c_hvmatb, c_jarray)

  end subroutine hvread

end module mod_che_hvread
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
