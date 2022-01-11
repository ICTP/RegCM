module parkind

  use mod_realkinds
  use mod_intkinds

  implicit none

  !------------------------------------------------------------------
  ! rrtmg kinds
  ! Define integer and real kinds for various types.
  !------------------------------------------------------------------
  !
  ! integer kinds
  ! -------------
  !
  integer, parameter :: kind_ib = ik8
  integer, parameter :: kind_im = ik4
  integer, parameter :: kind_in = kind(1) ! native integer
  !
  ! real kinds
  ! ----------
  !
#ifdef SINGLE_PRECISION_REAL
  integer, parameter :: kind_rb = rk4
  real(kind_rb) , parameter :: almostzero = 1.e-10_kind_rb
#else
  integer, parameter :: kind_rb = rk8
  real(kind_rb) , parameter :: almostzero = 1.e-20_kind_rb
#endif
  integer, parameter :: kind_rm = rk4
  integer, parameter :: kind_rn = kind(1.0) ! native real

end module parkind
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
