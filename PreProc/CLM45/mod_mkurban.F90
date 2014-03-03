module mod_mkurban
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_getwindow
  use mod_bilinear
  use mod_nchelper
  use mod_memutil
  use netcdf

  implicit none

  private

  public :: mkurban

  real :: vmin = 0.0
  real :: vmisdat = -9999.0

  contains

  subroutine mkurban(urbanfile,urban)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rk4) , dimension(:,:) , intent(out) :: urban
    urban = 0.0D0
  end subroutine mkurban

end module mod_mkurban
