module mod_cbmz_hvread

  use mod_cbmz_precision
  use mod_cbmz_jval1

  integer , public :: c_hvin
  integer , public , dimension(22) :: c_nhv
  real(kind=dp) , public , dimension(22) :: c_hvmatb
  real(kind=dp) , public , dimension(22,40) :: c_hvmat
  real(kind=dp) , public , dimension(80,510,56) :: c_jarray

  contains

  ! This reads the MADRONICH LOOKUP TABLE (2002 VERSION).
  !   Input file:  TUVGRID2
  !
  ! Result: creates the hv input arrays:
  !     c_hvin, c_nhv, c_hvmat, c_hvmatb, c_jarray
  !
  subroutine hvread
    implicit none
    c_hvin = 26
    call readhv(c_hvin,c_nhv,c_hvmat,c_hvmatb,c_jarray)
  end subroutine hvread

end module mod_cbmz_hvread

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
