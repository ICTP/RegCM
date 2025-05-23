!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_maputils

  use mod_constants
  use mod_intkinds
  use mod_realkinds
  use mod_projections

  private

  public :: getcoord, corpar, mappar, anyprojparams, regcm_projection

  contains

  subroutine getcoord(pjpara,lon,lat,pj,jx,iy)
    implicit none
    type(anyprojparams), intent(in) :: pjpara
    type(regcm_projection), intent(out) :: pj
    integer(ik4), intent(in) :: jx, iy
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: lat, lon
    integer :: i, j

    call pj%initialize(pjpara)
    do i = 1, iy
      do j = 1, jx
        call pj%ijll(real(j,rkx),real(i,rkx),lat(j,i),lon(j,i))
      end do
    end do
  end subroutine getcoord

  subroutine corpar(lat,coriol)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: lat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: coriol
    coriol = real(eomeg2*sin(lat*degrad),rkx)
  end subroutine corpar

  subroutine mappar(pj,xlat,xlon,mapf)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: mapf
    call pj%mapfac(xlat,xlon,mapf)
  end subroutine mappar

end module mod_maputils

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
