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
  implicit none

  private

  public :: getcoord, get_coord_bounds
  public :: corpar, mappar, anyprojparams, regcm_projection

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

  ! Compute CF-convention cell-corner coordinates for a grid.
  !
  ! The projection pj must already be initialised by getcoord for the
  ! desired stagger (cross, dot, u or v).  For each cell centre (j,i) the
  ! four corners are evaluated by offsetting the real grid indices by ±0.5
  ! and stored in counterclockwise order:
  !
  !   nv=1 : SW  (j-0.5, i-0.5)
  !   nv=2 : SE  (j+0.5, i-0.5)
  !   nv=3 : NE  (j+0.5, i+0.5)
  !   nv=4 : NW  (j-0.5, i+0.5)
  !
  ! Output arrays have shape (4, jx, iy), i.e. CF convention (nv, jx, iy).
  !
  subroutine get_coord_bounds(pj,lon_bnds,lat_bnds,jx,iy)
    implicit none
    type(regcm_projection), intent(in) :: pj
    integer(ik4), intent(in) :: jx, iy
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: lat_bnds
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: lon_bnds
    integer(ik4) :: i, j
    real(rkx) :: clat, clon
    real(rkx), parameter :: half = 0.5_rkx

    do i = 1, iy
      do j = 1, jx
        ! SW corner
        call pj%ijll(real(j,rkx)-half, real(i,rkx)-half, clat, clon)
        lat_bnds(1,j,i) = clat
        lon_bnds(1,j,i) = clon
        ! SE corner
        call pj%ijll(real(j,rkx)+half, real(i,rkx)-half, clat, clon)
        lat_bnds(2,j,i) = clat
        lon_bnds(2,j,i) = clon
        ! NE corner
        call pj%ijll(real(j,rkx)+half, real(i,rkx)+half, clat, clon)
        lat_bnds(3,j,i) = clat
        lon_bnds(3,j,i) = clon
        ! NW corner
        call pj%ijll(real(j,rkx)-half, real(i,rkx)+half, clat, clon)
        lat_bnds(4,j,i) = clat
        lon_bnds(4,j,i) = clon
      end do
    end do
  end subroutine get_coord_bounds

end module mod_maputils

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
