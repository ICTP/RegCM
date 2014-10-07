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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This code in its original form is coming from SCRIP v 1.4 library
!
!   THIS IS A DERIVATIVE WORK AND IS MARKED AS SUCH AS REQUESTED BY
!   THE ORIGINAL LICENSE BELOW.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Copyright (c) 1997, 1998 the Regents of the University of
!       California.
!
!     This software and ancillary information (herein called software)
!     called SCRIP is made available under the terms described here.
!     The software has been approved for release with associated
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with
!     the version available from Los Alamos National Laboratory.
!
!********************************************************************

module mod_scrip_grids

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_stdio
  use mod_message

  implicit none

  private

  logical :: init_sizes = .false.
  integer(ik4) :: grid1_size      ! total points of the grid
  integer(ik4) :: grid2_size      ! total points of the grid
  integer(ik4) :: grid1_rank      ! rank of the grid
  integer(ik4) :: grid2_rank      ! rank of the grid
  integer(ik4) :: grid1_corners   ! number of corners
  integer(ik4) :: grid2_corners   ! number of corners
  integer(ik4) , pointer , dimension(:) :: grid1_dims ! size of grid dimensions
  integer(ik4) , pointer , dimension(:) :: grid2_dims ! size of grid dimensions

  logical , pointer , dimension(:) :: grid1_mask ! flag which cells participate
  logical , pointer , dimension(:) :: grid2_mask ! flag which cells participate
  real(rk8) , pointer , dimension(:) :: grid1_center_lat  ! coordinates for
  real(rk8) , pointer , dimension(:) :: grid1_center_lon  ! each grid center
  real(rk8) , pointer , dimension(:) :: grid2_center_lat  ! lon coordinates for
  real(rk8) , pointer , dimension(:) :: grid2_center_lon  ! each grid center
  real(rk8) , pointer , dimension(:) :: grid1_area  ! tot area of each cell
  real(rk8) , pointer , dimension(:) :: grid2_area  ! tot area of each cell
  real(rk8) , pointer , dimension(:) :: grid1_frac  ! fractional area of grid
  real(rk8) , pointer , dimension(:) :: grid2_frac  ! fractional area of grid
  real(rk8) , pointer , dimension(:,:) :: grid1_corner_lat  ! coordinates for
  real(rk8) , pointer , dimension(:,:) :: grid1_corner_lon  ! each grid corner
  real(rk8) , pointer , dimension(:,:) :: grid2_corner_lat  ! coordinates for
  real(rk8) , pointer , dimension(:,:) :: grid2_corner_lon  ! each grid corner
  real(rk8) , pointer , dimension(:,:) :: grid1_bound_box ! lat/lon bounding box
  real(rk8) , pointer , dimension(:,:) :: grid2_bound_box ! lat/lon bounding box

  public :: grid1_size , grid1_rank , grid1_corners
  public :: grid1_dims , grid1_mask
  public :: grid1_center_lat , grid1_center_lon
  public :: grid1_area , grid1_frac
  public :: grid1_corner_lat , grid1_corner_lon
  public :: grid1_bound_box
  public :: grid2_size , grid2_rank , grid2_corners
  public :: grid2_dims , grid2_mask
  public :: grid2_center_lat , grid2_center_lon
  public :: grid2_area , grid2_frac
  public :: grid2_corner_lat , grid2_corner_lon
  public :: grid2_bound_box

  ! min,max adds for grid1 cells in this lat bin
  integer(ik4) , public , pointer , dimension(:,:) :: bin_addr1
  integer(ik4) , public , pointer , dimension(:,:) :: bin_addr2
  real(rk8) , public , pointer , dimension(:,:) :: bin_lats
  real(rk8) , public , pointer , dimension(:,:) :: bin_lons

  character(len=80) , public :: restrict_type = 'latlon'  ! type of bins to use
  integer(ik4) , public :: num_srch_bins = 90 ! num of bins for restricted srch
  logical , public :: luse_grid_centers = .true.

  public :: scrip_dims_init , scrip_grid_init

  contains

    subroutine scrip_dims_init(gasize,garank,gadims,gacorners, &
                               gbsize,gbrank,gbdims,gbcorners)
      implicit none
      integer(ik4) , intent(in) :: gasize , garank , gacorners
      integer(ik4) , intent(in) :: gbsize , gbrank , gbcorners
      integer(ik4) , dimension(:) , intent(in) :: gadims , gbdims

      grid1_size = gasize
      grid1_rank = garank
      grid1_corners = gacorners
      grid2_size = gbsize
      grid2_rank = gbrank
      grid2_corners = gbcorners

      call getmem1d(grid1_dims,1,grid1_rank,'scrip:grid1:rank')
      call getmem1d(grid2_dims,1,grid2_rank,'scrip:grid2:rank')

      grid1_dims(1:grid1_rank) = gadims(1:grid1_rank)
      grid2_dims(1:grid2_rank) = gbdims(1:grid2_rank)
      init_sizes = .true.
    end subroutine scrip_dims_init

    subroutine scrip_grid_init(gaclat,gaclon,gadlat,gadlon,gamask, &
                               gbclat,gbclon,gbdlat,gbdlon,gbmask)
      implicit none
      real(rk8) , dimension(:) , intent(in) :: gaclat , gaclon
      real(rk8) , dimension(:,:) , intent(in) :: gadlat , gadlon
      integer(ik4) , dimension(:) , intent(in) :: gamask
      real(rk8) , dimension(:) , intent(in) :: gbclat , gbclon
      real(rk8) , dimension(:,:) , intent(in) :: gbdlat , gbdlon
      integer(ik4) , dimension(:) , intent(in) :: gbmask

      integer(ik4) :: nx , ny , i , j , ip1 , jp1 , n
      integer(ik4) :: n_add , e_add , ne_add , nele , nsbsq

      ! lat/lon intervals for search bins
      real(rk8) :: dlat , dlon
      ! temps for computing bounding boxes
      real(rk8) , dimension(4) :: tmp_lats , tmp_lons

      if ( .not. init_sizes ) then
        call die('scrip_grid_init','Grid dimensions unknown',1)
      end if

      call getmem1d(grid1_mask,1,grid1_size,'scrip:grid1:mask')
      call getmem1d(grid1_center_lat,1,grid1_size,'scrip:grid1:centerlat')
      call getmem1d(grid1_center_lon,1,grid1_size,'scrip:grid1:centerlon')
      call getmem1d(grid1_area,1,grid1_size,'scrip:grid1:area')
      call getmem1d(grid1_frac,1,grid1_size,'scrip:grid1:frac')
      call getmem1d(grid1_mask,1,grid1_size,'scrip:grid1:mask')

      call getmem2d(grid1_corner_lat,1,grid1_corners, &
                                    1,grid1_size,'scrip:grid1:cornerlat')
      call getmem2d(grid1_corner_lon,1,grid1_corners, &
                                    1,grid1_size,'scrip:grid1:cornerlon')
      call getmem2d(grid1_bound_box,1,4, &
                                    1,grid1_size,'scrip:grid1:bound_box')

      call getmem1d(grid2_mask,1,grid2_size,'scrip:grid2:mask')
      call getmem1d(grid2_center_lat,1,grid2_size,'scrip:grid2:centerlat')
      call getmem1d(grid2_center_lon,1,grid2_size,'scrip:grid2:centerlon')
      call getmem1d(grid2_area,1,grid2_size,'scrip:grid2:area')
      call getmem1d(grid2_frac,1,grid2_size,'scrip:grid2:frac')
      call getmem1d(grid2_mask,1,grid2_size,'scrip:grid2:mask')

      call getmem2d(grid2_corner_lat,1,grid2_corners, &
                                    1,grid2_size,'scrip:grid2:cornerlat')
      call getmem2d(grid2_corner_lon,1,grid2_corners, &
                                    1,grid2_size,'scrip:grid2:cornerlon')
      call getmem2d(grid2_bound_box,1,4, &
                                    1,grid2_size,'scrip:grid2:bound_box')

      grid1_area = d_zero
      grid1_frac = d_zero
      grid1_center_lat = gaclat
      grid1_center_lon = gaclon
      grid1_corner_lat = gadlat
      grid1_corner_lon = gadlon
      where (gamask == 1)
        grid1_mask = .true.
      elsewhere
        grid1_mask = .false.
      endwhere
      grid1_center_lat = grid1_center_lat*degrad
      grid1_center_lon = grid1_center_lon*degrad
      grid1_corner_lat = grid1_corner_lat*degrad
      grid1_corner_lon = grid1_corner_lon*degrad
      grid2_area = d_zero
      grid2_frac = d_zero
      grid2_center_lat = gbclat
      grid2_center_lon = gbclon
      grid2_corner_lat = gbdlat
      grid2_corner_lon = gbdlon
      where (gbmask == 1)
        grid2_mask = .true.
      elsewhere
        grid2_mask = .false.
      endwhere
      grid2_center_lat = grid2_center_lat*degrad
      grid2_center_lon = grid2_center_lon*degrad
      grid2_corner_lat = grid2_corner_lat*degrad
      grid2_corner_lon = grid2_corner_lon*degrad
      !
      !  convert longitudes to 0,2pi interval
      !
      where (grid1_center_lon > twopi) &
           grid1_center_lon = grid1_center_lon - twopi
      where (grid1_center_lon < d_zero)  &
           grid1_center_lon = grid1_center_lon + twopi
      where (grid2_center_lon > twopi) &
           grid2_center_lon = grid2_center_lon - twopi
      where (grid2_center_lon < d_zero)  &
           grid2_center_lon = grid2_center_lon + twopi
      where (grid1_corner_lon > twopi) &
           grid1_corner_lon = grid1_corner_lon - twopi
      where (grid1_corner_lon < d_zero)  &
           grid1_corner_lon = grid1_corner_lon + twopi
      where (grid2_corner_lon > twopi) &
           grid2_corner_lon = grid2_corner_lon - twopi
      where (grid2_corner_lon < d_zero)  &
           grid2_corner_lon = grid2_corner_lon + twopi
      !
      ! make sure input latitude range is within the machine values
      ! for +/- pi/2
      !
      where (grid1_center_lat >  halfpi) grid1_center_lat =  halfpi
      where (grid1_corner_lat >  halfpi) grid1_corner_lat =  halfpi
      where (grid1_center_lat < -halfpi) grid1_center_lat = -halfpi
      where (grid1_corner_lat < -halfpi) grid1_corner_lat = -halfpi
      where (grid2_center_lat >  halfpi) grid2_center_lat =  halfpi
      where (grid2_corner_lat >  halfpi) grid2_corner_lat =  halfpi
      where (grid2_center_lat < -halfpi) grid2_center_lat = -halfpi
      where (grid2_corner_lat < -halfpi) grid2_corner_lat = -halfpi
      !
      ! compute bounding boxes for restricting future grid searches
      !
      if ( .not. luse_grid_centers ) then
        grid1_bound_box(1,:) = minval(grid1_corner_lat, dim=1)
        grid1_bound_box(2,:) = maxval(grid1_corner_lat, dim=1)
        grid1_bound_box(3,:) = minval(grid1_corner_lon, dim=1)
        grid1_bound_box(4,:) = maxval(grid1_corner_lon, dim=1)
        grid2_bound_box(1,:) = minval(grid2_corner_lat, dim=1)
        grid2_bound_box(2,:) = maxval(grid2_corner_lat, dim=1)
        grid2_bound_box(3,:) = minval(grid2_corner_lon, dim=1)
        grid2_bound_box(4,:) = maxval(grid2_corner_lon, dim=1)
      else
        nx = grid1_dims(1)
        ny = grid1_dims(2)
        do n = 1 , grid1_size
          !
          ! find N,S and NE points to this grid point
          !
          j = (n - 1)/nx +1
          i = n - (j-1)*nx
          if ( i < nx ) then
            ip1 = i + 1
          else
            !
            ! assume cyclic
            !
            ip1 = 1
            !
            ! but if it is not, correct
            !
            e_add = (j - 1)*nx + ip1
            if ( abs(grid1_center_lat(e_add) - &
                     grid1_center_lat(n   )) > halfpi ) then
              ip1 = i
            end if
          end if
          if ( j < ny ) then
            jp1 = j+1
          else
            !
            ! assume cyclic
            !
            jp1 = 1
            !
            ! but if it is not, correct
            !
            n_add = (jp1 - 1)*nx + i
            if ( abs(grid1_center_lat(n_add) - &
                     grid1_center_lat(n   )) > halfpi ) then
              jp1 = j
            end if
          end if
          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1
          !
          ! find N,S and NE lat/lon coords and check bounding box
          !
          tmp_lats(1) = grid1_center_lat(n)
          tmp_lats(2) = grid1_center_lat(e_add)
          tmp_lats(3) = grid1_center_lat(ne_add)
          tmp_lats(4) = grid1_center_lat(n_add)
          tmp_lons(1) = grid1_center_lon(n)
          tmp_lons(2) = grid1_center_lon(e_add)
          tmp_lons(3) = grid1_center_lon(ne_add)
          tmp_lons(4) = grid1_center_lon(n_add)
          grid1_bound_box(1,n) = minval(tmp_lats)
          grid1_bound_box(2,n) = maxval(tmp_lats)
          grid1_bound_box(3,n) = minval(tmp_lons)
          grid1_bound_box(4,n) = maxval(tmp_lons)
        end do
        nx = grid2_dims(1)
        ny = grid2_dims(2)
        do n = 1 , grid2_size
          !
          ! find N,S and NE points to this grid point
          !
          j = (n - 1)/nx +1
          i = n - (j-1)*nx
          if ( i < nx ) then
            ip1 = i + 1
          else
            !
            ! assume cyclic
            !
            ip1 = 1
            !
            ! but if it is not, correct
            !
            e_add = (j - 1)*nx + ip1
            if ( abs(grid2_center_lat(e_add) - &
                     grid2_center_lat(n   )) > halfpi ) then
              ip1 = i
            end if
          end if
          if ( j < ny ) then
            jp1 = j+1
          else
            !
            ! assume cyclic
            !
            jp1 = 1
            !
            ! but if it is not, correct
            !
            n_add = (jp1 - 1)*nx + i
            if ( abs(grid2_center_lat(n_add) - &
                     grid2_center_lat(n   )) > halfpi ) then
              jp1 = j
            end if
          end if
          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1
          !
          ! find N,S and NE lat/lon coords and check bounding box
          !
          tmp_lats(1) = grid2_center_lat(n)
          tmp_lats(2) = grid2_center_lat(e_add)
          tmp_lats(3) = grid2_center_lat(ne_add)
          tmp_lats(4) = grid2_center_lat(n_add)
          tmp_lons(1) = grid2_center_lon(n)
          tmp_lons(2) = grid2_center_lon(e_add)
          tmp_lons(3) = grid2_center_lon(ne_add)
          tmp_lons(4) = grid2_center_lon(n_add)
          grid2_bound_box(1,n) = minval(tmp_lats)
          grid2_bound_box(2,n) = maxval(tmp_lats)
          grid2_bound_box(3,n) = minval(tmp_lons)
          grid2_bound_box(4,n) = maxval(tmp_lons)
        end do
      end if
      where ( abs(grid1_bound_box(4,:) - &
                  grid1_bound_box(3,:)) > mathpi)
        grid1_bound_box(3,:) = d_zero
        grid1_bound_box(4,:) = twopi
      end where
      where ( abs(grid2_bound_box(4,:) - &
                  grid2_bound_box(3,:)) > mathpi)
        grid2_bound_box(3,:) = d_zero
        grid2_bound_box(4,:) = twopi
      end where
      !
      ! try to check for cells that overlap poles
      !
      where (grid1_center_lat > grid1_bound_box(2,:)) &
        grid1_bound_box(2,:) = halfpi
      where (grid1_center_lat < grid1_bound_box(1,:)) &
        grid1_bound_box(1,:) = -halfpi
      where (grid2_center_lat > grid2_bound_box(2,:)) &
        grid2_bound_box(2,:) = halfpi
      where (grid2_center_lat < grid2_bound_box(1,:)) &
        grid2_bound_box(1,:) = -halfpi
      !
      ! set up and assign address ranges to search bins in order to
      ! further restrict later searches
      !
      select case ( restrict_type )
        case ('latitude')
          if ( num_srch_bins <= 0 ) then
            call die('scrip_grid_init','num_srch_bins <= 0 ?',1)
          end if
          call getmem2d(bin_addr1,1,2,1,num_srch_bins,'scrip:bin_addr1')
          call getmem2d(bin_addr2,1,2,1,num_srch_bins,'scrip:bin_addr2')
          call getmem2d(bin_lats,1,2,1,num_srch_bins,'scrip:bin_lats')
          call getmem2d(bin_lons,1,2,1,num_srch_bins,'scrip:bin_lons')

          dlat = mathpi/num_srch_bins
          do n = 1 , num_srch_bins
            bin_lats(1,n) = (n-1)*dlat - halfpi
            bin_lats(2,n) =     n*dlat - halfpi
            bin_lons(1,n) = d_zero
            bin_lons(2,n) = twopi
            bin_addr1(1,n) = grid1_size + 1
            bin_addr1(2,n) = 0
            bin_addr2(1,n) = grid2_size + 1
            bin_addr2(2,n) = 0
          end do
          do nele = 1 , grid1_size
            do n = 1 , num_srch_bins
              if ( grid1_bound_box(1,nele) <= bin_lats(2,n) .and. &
                   grid1_bound_box(2,nele) >= bin_lats(1,n) ) then
                bin_addr1(1,n) = min(nele,bin_addr1(1,n))
                bin_addr1(2,n) = max(nele,bin_addr1(2,n))
              end if
            end do
          end do
          do nele = 1 , grid2_size
            do n = 1 , num_srch_bins
              if ( grid2_bound_box(1,nele) <= bin_lats(2,n) .and. &
                   grid2_bound_box(2,nele) >= bin_lats(1,n) ) then
                bin_addr2(1,n) = min(nele,bin_addr2(1,n))
                bin_addr2(2,n) = max(nele,bin_addr2(2,n))
              end if
            end do
          end do
        case ('latlon')
          if ( num_srch_bins <= 0 ) then
            call die('scrip_grid_init','num_srch_bins <= 0 ?',1)
          end if
          dlat = mathpi/num_srch_bins
          dlon = twopi/num_srch_bins
          nsbsq = num_srch_bins*num_srch_bins
          call getmem2d(bin_addr1,1,2,1,nsbsq,'scrip:bin_addr1')
          call getmem2d(bin_addr2,1,2,1,nsbsq,'scrip:bin_addr2')
          call getmem2d(bin_lats,1,2,1,nsbsq,'scrip:bin_lats')
          call getmem2d(bin_lons,1,2,1,nsbsq,'scrip:bin_lons')
          n = 0
          do j = 1 , num_srch_bins
            do i = 1 , num_srch_bins
              n = n + 1
              bin_lats(1,n) = (j-1)*dlat - halfpi
              bin_lats(2,n) =     j*dlat - halfpi
              bin_lons(1,n) = (i-1)*dlon
              bin_lons(2,n) =     i*dlon
              bin_addr1(1,n) = grid1_size + 1
              bin_addr1(2,n) = 0
              bin_addr2(1,n) = grid2_size + 1
              bin_addr2(2,n) = 0
            end do
          end do

          num_srch_bins = num_srch_bins**2
          do nele = 1 , grid1_size
            do n = 1 , num_srch_bins
              if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and. &
                  grid1_bound_box(2,nele) >= bin_lats(1,n) .and. &
                  grid1_bound_box(3,nele) <= bin_lons(2,n) .and. &
                  grid1_bound_box(4,nele) >= bin_lons(1,n) ) then
                bin_addr1(1,n) = min(nele,bin_addr1(1,n))
                bin_addr1(2,n) = max(nele,bin_addr1(2,n))
              end if
            end do
          end do
          do nele = 1 , grid2_size
            do n = 1 , num_srch_bins
              if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and. &
                  grid2_bound_box(2,nele) >= bin_lats(1,n) .and. &
                  grid2_bound_box(3,nele) <= bin_lons(2,n) .and. &
                  grid2_bound_box(4,nele) >= bin_lons(1,n) ) then
                bin_addr2(1,n) = min(nele,bin_addr1(1,n))
                bin_addr2(2,n) = max(nele,bin_addr1(2,n))
              end if
            end do
          end do
        case default
          call die('scrip_grid_init','unknown search restriction method',1)
      end select
    end subroutine scrip_grid_init

end module mod_scrip_grids
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
