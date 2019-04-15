module mod_cbmz_precision

  use mod_realkinds , only : rk4 , rk8 , rk16
  use mod_intkinds , only : ik4
  !
  ! definition of different levels of accuracy
  ! for real variables using kind parameterization
  !
  ! kpp sp - single precision kind
  integer, parameter :: sp = rk4
  ! kpp dp - double precision kind
  integer, parameter :: dp = rk8
  ! kpp qp - quadruple precision kind
#ifndef QUAD_PRECISION
  integer, parameter :: qp = rk8
#else
  integer, parameter :: qp = rk16
#endif

end module mod_cbmz_precision

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
