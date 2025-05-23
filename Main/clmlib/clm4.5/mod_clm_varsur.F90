module mod_clm_varsur
  !
  ! Module containing 2-d surface boundary data information
  !
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save
  !
  ! surface boundary data, these are all "gdc" local
  !
  ! percent of spec lunits wrt gcell
  real(rk8), public, allocatable, dimension(:) :: pctspec
  ! vegetation type
  integer(ik4), public, allocatable, dimension(:,:) :: vegxy
  ! subgrid weights
  real(rk8), public, allocatable, dimension(:,:), target :: wtxy

end module mod_clm_varsur
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
