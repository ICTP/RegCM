!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef OASIS_CPL

module mod_oasis
  !
  ! Placeholder for OASIS coupling
  ! ------------------------------
  !
  ! RegCM grid is like this:
  !
  ! O   O   O   O   O      This is an example jx = 5 , iy = 4 GRID
  !   X   X   X   X        The J is the index on WE (Left -> Right)
  ! O   O   O   O   O      The I is the index on SN (Bottom -> Top)
  !   X   X   X   X        On dot (O) point we have VECTORS (U,V)
  ! O   O   O   O   O      On cross (X) point we have SCALARS (T,Q,P,..)
  !   X   X   X   X        Note that the jx,iy refer to the DOT GRID !!
  ! O   O   O   O   O      The indexes to move on grid are:
  !                            jc,ic for cross grid
  !                            jd,id for dot grid
  !                        Modifiers :
  !                            e stands for EXTERNAL GRID (WITH BOUNDARY)
  !                            i stands for INTERNAL GRID (NO BOUNDARY)
  !                        Limit :
  !                            1 is LOW limit
  !                            2 is UP limit
  !  Therefore:
  !
  !               ici1 , ici2
  !
  !      Means from bottom to top for this processor without external
  !      boundary points on CROSS points
  !
  !               jde1 , jde2
  !
  !      Meeans from left to right for this processor with exernal
  !      boundary points on DOT points
  !
  ! The External Boundary Points are updated only in mod_bdycod and used
  ! in advection computation or to go from DOT to CROSS and viceversa.
  ! All model physics is on INTERNAL CROSS POINTS.
  !
  use mod_realkinds

  implicit none

  private

  real(rk8) , pointer , dimension(:,:) :: dot_point_grid_lat
  real(rk8) , pointer , dimension(:,:) :: dot_point_grid_lon
  real(rk8) , pointer , dimension(:,:) :: cross_point_grid_lat
  real(rk8) , pointer , dimension(:,:) :: cross_point_grid_lon

  real(rk8) , pointer , dimension(:,:) :: global_cross_lats
  real(rk8) , pointer , dimension(:,:) :: global_cross_lons

  public :: oasis_init

  contains

  subroutine oasis_init
    use mod_atm_interface , only : mddom
    implicit none

    ! MDDOM type contains per processor grid information

    call assignpnt(mddom%xlat,cross_point_grid_lat)
    call assignpnt(mddom%xlon,cross_point_grid_lon)
    call assignpnt(mddom%dlat,dot_point_grid_lat)
    call assignpnt(mddom%dlon,dot_point_grid_lon)

    ! In cross_point_grid_lat , cross_point_grid_lon we have
    ! processor specific coordinates on CROSS points
    ! In dot_point_grid_lat , dot_point_grid_lon we have
    ! processor specific coordinates on DOT points

    ! Global data can be collected like this:

    call getmem2d(global_cross_lats, &
                  jcross1,jcross2,icross1,icross2,'GLOBAL_CROSS_LAT')
    call getmem2d(global_cross_lons, &
                  jcross1,jcross2,icross1,icross2,'GLOBAL_CROSS_LON')

    call grid_collect(cross_point_grid_lat,global_cross_lats, &
                      jce1,jce2,ice1,ice2)
    call grid_collect(cross_point_grid_lon,global_cross_lons, &
                      jce1,jce2,ice1,ice2)

    ! global_cross_lats , global_cross_lons have now the
    ! GLOBAL CROSS GRID COORDINATES
  end subroutine oasis_init

end module mod_oasis

#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
