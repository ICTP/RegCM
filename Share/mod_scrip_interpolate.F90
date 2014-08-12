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

module mod_scrip_interpolate

  use mod_realkinds
  use mod_intkinds
  use mod_scrip_grids
  use mod_scrip_remap
  use mod_scrip_remap_vars
  use mod_scrip_remap_bicubic

  implicit none

  private

  interface interpolate
    module procedure interpolate_2d
  end interface interpolate

  public :: init_scrip_library
  public :: interpolate
  public :: release_scrip_library

  real(rk8) , pointer , dimension(:) :: inpgrid1d
  real(rk8) , pointer , dimension(:) :: outgrid1d
  integer , dimension(1) :: inpgridsize
  integer , dimension(1) :: outgridsize

  contains

    subroutine init_scrip_library(map,g1sn,g1we,g2sn,g2we, &
                                  g1lat,g1lon,g1mask,g2lat,g2lon,g2mask)
      implicit none
      integer(ik4), intent(in) :: map
      integer(ik4), intent(in) :: g1sn , g1we
      integer(ik4), intent(in) :: g2sn , g2we
      real(rk8) , pointer , dimension(:,:) , intent(in) :: g1lat , g1lon
      real(rk8) , pointer , dimension(:,:) , intent(in) :: g2lat , g2lon
      real(rk8) , pointer , dimension(:,:) , intent(in) :: g1mask , g2mask

      real(rk8) , allocatable , dimension(:) :: a1lat , a1lon
      real(rk8) , allocatable , dimension(:) :: a2lat , a2lon
      real(rk8) , allocatable , dimension(:,:) :: a1dlat , a1dlon
      real(rk8) , allocatable , dimension(:,:) :: a2dlat , a2dlon
      integer(ik4) , allocatable , dimension(:) :: a1mask , a2mask

      call scrip_dims_init(g1sn*g1we,2,(/g1we,g1sn/),4, &
                           g2sn*g2we,2,(/g2we,g2sn/),4)

    end subroutine init_scrip_library

    subroutine interpolate_2d(inpgrid,outgrid)
      implicit none
      real(rk8), pointer , dimension(:,:) :: inpgrid
      real(rk8), pointer , dimension(:,:) :: outgrid

      inpgrid1d = reshape(inpgrid, inpgridsize)

      select case (map_type)
        case(map_type_conserv)
          call remap(outgrid1d, &
                     wts_map1,grid2_add_map1,grid1_add_map1, &
                     inpgrid1d)
        case(map_type_bilinear)
          call remap(outgrid1d, &
                     wts_map1,grid2_add_map1,grid1_add_map1, &
                     inpgrid1d)
        case(map_type_bicubic)
          call remap(outgrid1d, &
                     wts_map1,grid2_add_map1,grid1_add_map1, &
                     inpgrid1d, &
                     src_grad1=grad1_lat, &
                     src_grad2=grad1_lon, &
                     src_grad3=grad1_latlon)
        case(map_type_distwgt)
          call remap(outgrid1d, &
                     wts_map1,grid2_add_map1,grid1_add_map1, &
                     inpgrid1d)
      end select
    end subroutine interpolate_2d

    subroutine release_scrip_library
      implicit none
    end subroutine release_scrip_library

end module mod_scrip_interpolate
