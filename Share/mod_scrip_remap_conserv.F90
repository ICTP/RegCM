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
!
! this module contains necessary routines for computing addresses
! and weights for a conservative interpolation  between any two
! grids on a sphere.  the weights are computed by performing line
! integrals around all overlap regions of the two grids.  see
! Dukowicz and Kodis, SIAM J. Sci. Stat. Comput. 8, 305 (1987) and
! Jones, P.W. Monthly Weather Review (submitted).
!
!
module mod_scrip_remap_conserv

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_scrip_grids
  use mod_scrip_remap_vars

  implicit none

  private

  ! num cells in restricted search arrays
  integer(ik4) , public :: num_srch_cells

  ! thresholds for coord transf.
  real(rk8) , parameter :: north_thresh =  1.45D0
  real(rk8) , parameter :: south_thresh = -2.00D0

  ! global address of cells in srch arrays
  integer(ik4) , dimension(:) , allocatable , save :: srch_add

  ! coordinates of each corner of srch cells
  real(rk8) , dimension(:,:) , allocatable , save :: srch_corner_lat
  real(rk8) , dimension(:,:) , allocatable , save :: srch_corner_lon

  public :: remap_conserv

  contains

    subroutine remap_conserv
      implicit none
      !
      ! this routine traces the perimeters of every grid cell on each
      ! grid checking for intersections with the other grid and computing
      ! line integrals for each subsegment.
      !
      ! max number of subsegments per segment to prevent infinite loop
      integer(ik4) , parameter :: max_subseg = 10000

      integer(ik4) :: grid1_add  ! current linear address for grid1 cell
      integer(ik4) :: grid2_add  ! current linear address for grid2 cell
      integer(ik4) :: min_add    ! addresses for restricting search of
      integer(ik4) :: max_add    !   destination grid
      integer(ik4) :: n , nwgt   ! generic counters
      integer(ik4) :: corner     ! corner of cell that segment starts from
      integer(ik4) :: next_corn  ! corner of cell that segment ends on
      integer(ik4) :: num_subseg ! number of subsegments
      logical :: lcoinc   ! flag for coincident segments
      logical :: lrevers  ! flag for reversing direction of segment
      logical :: lbegin   ! flag for first integration of a segment
      ! mask for restricting searches
      logical , dimension(:) , allocatable :: srch_mask
      real(rk8) :: intrsct_lat , intrsct_lon ! lat/lon of next intersect
      real(rk8) :: beglat , endlat           ! endpoints of current seg.
      real(rk8) :: beglon , endlon           ! endpoints of current seg.
      real(rk8) :: norm_factor               ! factor for normalizing wts
      ! centroid coords on each grid
      real(rk8) , dimension(:) , allocatable :: grid1_centroid_lat
      real(rk8) , dimension(:) , allocatable :: grid1_centroid_lon
      real(rk8) , dimension(:) , allocatable :: grid2_centroid_lat
      real(rk8) , dimension(:) , allocatable :: grid2_centroid_lon
      real(rk8) , dimension(2) :: begseg ! begin lat/lon for
                                         ! full segment
      real(rk8) , dimension(6) :: weights ! local wgt array
      !
      ! initialize centroid arrays
      !
      allocate(grid1_centroid_lat(grid1_size), &
               grid1_centroid_lon(grid1_size), &
               grid2_centroid_lat(grid2_size), &
               grid2_centroid_lon(grid2_size))
      grid1_centroid_lat = d_zero
      grid1_centroid_lon = d_zero
      grid2_centroid_lat = d_zero
      grid2_centroid_lon = d_zero
      !
      ! integrate around each cell on grid1
      !
      allocate(srch_mask(grid2_size))
      do grid1_add = 1 , grid1_size
        !
        ! restrict searches first using search bins
        !
        ! call timer_start(1)
        min_add = grid2_size
        max_add = 1
        do n = 1 , num_srch_bins
          if ( grid1_add >= bin_addr1(1,n) .and. &
               grid1_add <= bin_addr1(2,n)) then
            min_add = min(min_add, bin_addr2(1,n))
            max_add = max(max_add, bin_addr2(2,n))
          end if
        end do
        !
        ! further restrict searches using bounding boxes
        !
        num_srch_cells = 0
        do grid2_add = min_add , max_add
          srch_mask(grid2_add) = (grid2_bound_box(1,grid2_add) <=     &
                                  grid1_bound_box(2,grid1_add)) .and. &
                                 (grid2_bound_box(2,grid2_add) >=     &
                                  grid1_bound_box(1,grid1_add)) .and. &
                                 (grid2_bound_box(3,grid2_add) <=     &
                                  grid1_bound_box(4,grid1_add)) .and. &
                                 (grid2_bound_box(4,grid2_add) >=     &
                                  grid1_bound_box(3,grid1_add))
          if ( srch_mask(grid2_add) ) num_srch_cells = num_srch_cells+1
        end do
        !
        ! create search arrays
        !
        allocate(srch_add(num_srch_cells),                      &
                 srch_corner_lat(grid2_corners,num_srch_cells), &
                 srch_corner_lon(grid2_corners,num_srch_cells))
        n = 0
        gather1: &
        do grid2_add = min_add , max_add
          if ( srch_mask(grid2_add) ) then
            n = n+1
            srch_add(n) = grid2_add
            srch_corner_lat(:,n) = grid2_corner_lat(:,grid2_add)
            srch_corner_lon(:,n) = grid2_corner_lon(:,grid2_add)
          end if
        end do gather1
        ! call timer_stop(1)
        !
        ! integrate around this cell
        !
        do corner = 1 , grid1_corners
          next_corn = mod(corner,grid1_corners) + 1
          !
          ! define endpoints of the current segment
          !
          beglat = grid1_corner_lat(corner,grid1_add)
          beglon = grid1_corner_lon(corner,grid1_add)
          endlat = grid1_corner_lat(next_corn,grid1_add)
          endlon = grid1_corner_lon(next_corn,grid1_add)
          lrevers = .false.
          !
          ! to ensure exact path taken during both
          ! sweeps, always integrate segments in the same
          ! direction (SW to NE)
          !
          if ( (endlat < beglat) .or. &
               (endlat == beglat .and. endlon < beglon) ) then
            beglat = grid1_corner_lat(next_corn,grid1_add)
            beglon = grid1_corner_lon(next_corn,grid1_add)
            endlat = grid1_corner_lat(corner,grid1_add)
            endlon = grid1_corner_lon(corner,grid1_add)
            lrevers = .true.
          end if
          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.
          num_subseg = 0
          !
          ! if this is a constant-longitude segment, skip the rest
          ! since the line integral contribution will be zero.
          !
          if ( endlon /= beglon ) then
            !
            ! integrate along this segment, detecting intersections
            ! and computing the line integral for each sub-segment
            !
            do while ( beglat /= endlat .or. beglon /= endlon )
              !
              ! prevent infinite loops if integration gets stuck
              ! near cell or threshold boundary
              !
              num_subseg = num_subseg + 1
              if ( num_subseg > max_subseg ) then
                call die('remap_conserv', &
                  'integration stalled: num_subseg exceeded limit',1)
              end if
              !
              ! find next intersection of this segment with a grid
              ! line on grid 2
              !
              ! call timer_start(2)
              call intersection(grid2_add,intrsct_lat,intrsct_lon,lcoinc, &
                                beglat, beglon, endlat, endlon, begseg,   &
                                lbegin, lrevers)
              ! call timer_stop(2)
              lbegin = .false.
              !
              ! compute line integral for this subsegment
              !
              ! call timer_start(3)
              if ( grid2_add /= 0 ) then
                call line_integral(weights, num_wts,                         &
                                   beglon, intrsct_lon, beglat, intrsct_lat, &
                                   grid1_center_lat(grid1_add),              &
                                   grid1_center_lon(grid1_add),              &
                                   grid2_center_lat(grid2_add),              &
                                   grid2_center_lon(grid2_add))
              else
                call line_integral(weights, num_wts,                         &
                                   beglon, intrsct_lon, beglat, intrsct_lat, &
                                   grid1_center_lat(grid1_add),              &
                                   grid1_center_lon(grid1_add),              &
                                   grid1_center_lat(grid1_add),              &
                                   grid1_center_lon(grid1_add))
              end if
              !call timer_stop(3)
              !
              ! if integrating in reverse order, change
              ! sign of weights
              !
              if ( lrevers ) then
                weights = -weights
              end if
              !
              ! store the appropriate addresses and weights.
              ! also add contributions to cell areas and centroids.
              !
              !if ( grid1_add == 119247 ) then
              !  print *,grid1_add,grid2_add,corner,weights(1)
              !  print *,grid1_corner_lat(:,grid1_add)
              !  print *,grid1_corner_lon(:,grid1_add)
              !  print *,grid2_corner_lat(:,grid2_add)
              !  print *,grid2_corner_lon(:,grid2_add)
              !  print *,beglat,beglon,intrsct_lat,intrsct_lon
              !end if
              !
              if ( grid2_add /= 0 ) then
                if ( grid1_mask(grid1_add) ) then
                  ! call timer_start(4)
                  call store_link_cnsrv(grid1_add, grid2_add, weights)
                  ! call timer_stop(4)
                  grid1_frac(grid1_add) = grid1_frac(grid1_add) + &
                                          weights(1)
                  grid2_frac(grid2_add) = grid2_frac(grid2_add) + &
                                          weights(num_wts+1)
                end if
              end if

              grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
              grid1_centroid_lat(grid1_add) = &
                grid1_centroid_lat(grid1_add) + weights(2)
              grid1_centroid_lon(grid1_add) = &
                grid1_centroid_lon(grid1_add) + weights(3)
              !
              ! reset beglat and beglon for next subsegment.
              !
              beglat = intrsct_lat
              beglon = intrsct_lon
            end do
          end if
          !
          ! end of segment
          !
        end do
        !
        ! finished with this cell: deallocate search array and
        ! start on next cell
        !
        deallocate(srch_add, srch_corner_lat, srch_corner_lon)
      end do
      deallocate(srch_mask)
      !
      ! integrate around each cell on grid2
      !
      allocate(srch_mask(grid1_size))
      do grid2_add = 1 , grid2_size
        !
        ! restrict searches first using search bins
        !
        !call timer_start(5)
        min_add = grid1_size
        max_add = 1
        do n = 1 , num_srch_bins
          if ( grid2_add >= bin_addr2(1,n) .and. &
               grid2_add <= bin_addr2(2,n) ) then
            min_add = min(min_add, bin_addr1(1,n))
            max_add = max(max_add, bin_addr1(2,n))
          end if
        end do
        !
        ! further restrict searches using bounding boxes
        !
        num_srch_cells = 0
        do grid1_add = min_add , max_add
          srch_mask(grid1_add) = (grid1_bound_box(1,grid1_add) <=     &
                                  grid2_bound_box(2,grid2_add)) .and. &
                                 (grid1_bound_box(2,grid1_add) >=     &
                                  grid2_bound_box(1,grid2_add)) .and. &
                                 (grid1_bound_box(3,grid1_add) <=     &
                                  grid2_bound_box(4,grid2_add)) .and. &
                                 (grid1_bound_box(4,grid1_add) >=     &
                                  grid2_bound_box(3,grid2_add))
          if ( srch_mask(grid1_add) ) num_srch_cells = num_srch_cells+1
        end do
        allocate(srch_add(num_srch_cells),                      &
                 srch_corner_lat(grid1_corners,num_srch_cells), &
                 srch_corner_lon(grid1_corners,num_srch_cells))
        n = 0
        gather2: &
        do grid1_add = min_add , max_add
          if ( srch_mask(grid1_add) ) then
            n = n+1
            srch_add(n) = grid1_add
            srch_corner_lat(:,n) = grid1_corner_lat(:,grid1_add)
            srch_corner_lon(:,n) = grid1_corner_lon(:,grid1_add)
          end if
        end do gather2
        ! call timer_stop(5)
        !
        ! integrate around this cell
        !
        do corner = 1 , grid2_corners
          next_corn = mod(corner,grid2_corners) + 1
          beglat = grid2_corner_lat(corner,grid2_add)
          beglon = grid2_corner_lon(corner,grid2_add)
          endlat = grid2_corner_lat(next_corn,grid2_add)
          endlon = grid2_corner_lon(next_corn,grid2_add)
          lrevers = .false.
          !
          ! to ensure exact path taken during both
          ! sweeps, always integrate in the same direction
          !
          if ( (endlat < beglat) .or. &
               (endlat == beglat .and. endlon < beglon) ) then
            beglat = grid2_corner_lat(next_corn,grid2_add)
            beglon = grid2_corner_lon(next_corn,grid2_add)
            endlat = grid2_corner_lat(corner,grid2_add)
            endlon = grid2_corner_lon(corner,grid2_add)
            lrevers = .true.
          end if
          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.
          !
          ! if this is a constant-longitude segment, skip the rest
          ! since the line integral contribution will be zero.
          !
          if ( endlon /= beglon ) then
            num_subseg = 0
            !
            ! integrate along this segment, detecting intersections
            ! and computing the line integral for each sub-segment
            !
            do while ( beglat /= endlat .or. beglon /= endlon )
              !
              ! prevent infinite loops if integration gets stuck
              ! near cell or threshold boundary
              !
              num_subseg = num_subseg + 1
              if ( num_subseg > max_subseg ) then
                call die('remap_conserv', &
                  'integration stalled: num_subseg exceeded limit',1)
              end if
              !
              ! find next intersection of this segment with a line
              ! on grid 2.
              !
              ! call timer_start(6)
              call intersection(grid1_add,intrsct_lat,intrsct_lon,lcoinc, &
                                beglat, beglon, endlat, endlon, begseg,   &
                                lbegin, lrevers)
              ! call timer_stop(6)
              lbegin = .false.
              !
              ! compute line integral for this subsegment.
              !
              ! call timer_start(7)
              if ( grid1_add /= 0 ) then
                call line_integral(weights, num_wts,                         &
                                   beglon, intrsct_lon, beglat, intrsct_lat, &
                                   grid1_center_lat(grid1_add),              &
                                   grid1_center_lon(grid1_add),              &
                                   grid2_center_lat(grid2_add),              &
                                   grid2_center_lon(grid2_add))
              else
                call line_integral(weights, num_wts,                         &
                                   beglon, intrsct_lon, beglat, intrsct_lat, &
                                   grid2_center_lat(grid2_add),              &
                                   grid2_center_lon(grid2_add),              &
                                   grid2_center_lat(grid2_add),              &
                                   grid2_center_lon(grid2_add))
              end if
              ! call timer_stop(7)
              if (lrevers) then
                weights = -weights
              end if
              !
              ! store the appropriate addresses and weights.
              ! also add contributions to cell areas and centroids.
              ! if there is a coincidence, do not store weights
              ! because they have been captured in the previous loop.
              ! the grid1 mask is the master mask
              !
              !if ( grid1_add == 119247 ) then
              !  print *,grid1_add,grid2_add,corner,weights(1)
              !  print *,grid1_corner_lat(:,grid1_add)
              !  print *,grid1_corner_lon(:,grid1_add)
              !  print *,grid2_corner_lat(:,grid2_add)
              !  print *,grid2_corner_lon(:,grid2_add)
              !  print *,beglat,beglon,intrsct_lat,intrsct_lon
              !end if
              !
              if ( .not. lcoinc .and. grid1_add /= 0 ) then
                if ( grid1_mask(grid1_add) ) then
                  ! call timer_start(8)
                  call store_link_cnsrv(grid1_add, grid2_add, weights)
                  ! call timer_stop(8)
                  grid1_frac(grid1_add) = grid1_frac(grid1_add) + &
                                          weights(1)
                  grid2_frac(grid2_add) = grid2_frac(grid2_add) + &
                                          weights(num_wts+1)
                end if
              end if
              grid2_area(grid2_add) = grid2_area(grid2_add) + &
                                      weights(num_wts+1)
              grid2_centroid_lat(grid2_add) = &
                grid2_centroid_lat(grid2_add) + weights(num_wts+2)
              grid2_centroid_lon(grid2_add) = &
                grid2_centroid_lon(grid2_add) + weights(num_wts+3)
              !
              ! reset beglat and beglon for next subsegment.
              !
              beglat = intrsct_lat
              beglon = intrsct_lon
            end do
          end if
          !
          ! end of segment
          !
        end do
        !
        ! finished with this cell: deallocate search array and
        ! start on next cell
        !
        deallocate(srch_add, srch_corner_lat, srch_corner_lon)
      end do
      deallocate(srch_mask)
      !
      ! correct for situations where N/S pole not explicitly included in
      ! grid (i.e. as a grid corner point). if pole is missing from only
      ! one grid, need to correct only the area and centroid of that
      ! grid.  if missing from both, do complete weight calculation.
      !
      ! North Pole
      !
      weights(1) =  twopi
      weights(2) =  mathpi*mathpi
      weights(3) =  d_zero
      weights(4) =  twopi
      weights(5) =  mathpi*mathpi
      weights(6) =  d_zero
      grid1_add = 0
      pole_loop1: &
      do n = 1 , grid1_size
        if ( grid1_area(n) < -d_three*halfpi .and. &
             grid1_center_lat(n) > d_zero ) then
          grid1_add = n
          exit pole_loop1
        end if
      end do pole_loop1
      grid2_add = 0
      pole_loop2: &
      do n = 1 , grid2_size
        if ( grid2_area(n) < -d_three*halfpi .and. &
             grid2_center_lat(n) > d_zero ) then
          grid2_add = n
          exit pole_loop2
        end if
      end do pole_loop2
      if ( grid1_add /=0 ) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) = &
          grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) = &
          grid1_centroid_lon(grid1_add) + weights(3)
      end if
      if ( grid2_add /=0 ) then
        grid2_area(grid2_add) = grid2_area(grid2_add) + weights(num_wts+1)
        grid2_centroid_lat(grid2_add) =  &
          grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) = &
          grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      end if
      if ( grid1_add /= 0 .and. grid2_add /=0 ) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)
        grid1_frac(grid1_add) = grid1_frac(grid1_add) + weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + weights(num_wts+1)
      end if
      !
      ! South Pole
      !
      weights(1) =  twopi
      weights(2) = -mathpi*mathpi
      weights(3) =  d_zero
      weights(4) =  twopi
      weights(5) = -mathpi*mathpi
      weights(6) =  d_zero
      grid1_add = 0
      pole_loop3: &
      do n = 1 , grid1_size
        if ( grid1_area(n) < -d_three*halfpi .and. &
             grid1_center_lat(n) < d_zero ) then
          grid1_add = n
          exit pole_loop3
        end if
      end do pole_loop3
      grid2_add = 0
      pole_loop4: &
      do n = 1 , grid2_size
        if ( grid2_area(n) < -d_three*halfpi .and. &
             grid2_center_lat(n) < d_zero ) then
          grid2_add = n
          exit pole_loop4
        end if
      end do pole_loop4
      if ( grid1_add /=0 ) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) = &
          grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) = &
          grid1_centroid_lon(grid1_add) + weights(3)
      end if
      if ( grid2_add /=0 ) then
        grid2_area(grid2_add) = grid2_area(grid2_add) + weights(num_wts+1)
        grid2_centroid_lat(grid2_add) = &
          grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) = &
          grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      end if
      if ( grid1_add /= 0 .and. grid2_add /=0 ) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)
        grid1_frac(grid1_add) = grid1_frac(grid1_add) + weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + weights(num_wts+1)
      end if
      !
      ! finish centroid computation
      !
      where ( grid1_area /= d_zero )
        grid1_centroid_lat = grid1_centroid_lat/grid1_area
        grid1_centroid_lon = grid1_centroid_lon/grid1_area
      end where
      where ( grid2_area /= d_zero )
        grid2_centroid_lat = grid2_centroid_lat/grid2_area
        grid2_centroid_lon = grid2_centroid_lon/grid2_area
      end where
      !
      ! include centroids in weights and normalize using destination
      ! area if requested
      !
      do n = 1 , num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        do nwgt = 1 , num_wts
          weights(        nwgt) = wts_map1(nwgt,n)
          if ( num_maps > 1 ) then
            weights(num_wts+nwgt) = wts_map2(nwgt,n)
          end if
        end do
        select case(norm_opt)
          case (norm_opt_dstarea)
            if ( grid2_area(grid2_add) /= d_zero ) then
              norm_factor = d_one/grid2_area(grid2_add)
            else
              norm_factor = d_zero
            end if
          case (norm_opt_frcarea)
            if ( grid2_frac(grid2_add ) /= d_zero) then
              norm_factor = d_one/grid2_frac(grid2_add)
            else
              norm_factor = d_zero
            end if
          case (norm_opt_none)
            norm_factor = d_one
        end select

        wts_map1(1,n) =  weights(1)*norm_factor
        wts_map1(2,n) = (weights(2) - weights(1)* &
                         grid1_centroid_lat(grid1_add))*norm_factor
        wts_map1(3,n) = (weights(3) - weights(1)* &
                         grid1_centroid_lon(grid1_add))*norm_factor
        if ( num_maps > 1 ) then
          select case(norm_opt)
            case (norm_opt_dstarea)
              if ( grid1_area(grid1_add) /= d_zero ) then
                norm_factor = d_one/grid1_area(grid1_add)
              else
                norm_factor = d_zero
              end if
            case (norm_opt_frcarea)
              if ( grid1_frac(grid1_add) /= d_zero ) then
                norm_factor = d_one/grid1_frac(grid1_add)
              else
                norm_factor = d_zero
              end if
            case (norm_opt_none)
              norm_factor = d_one
          end select
          wts_map2(1,n) =  weights(num_wts+1)*norm_factor
          wts_map2(2,n) = (weights(num_wts+2) - weights(num_wts+1)* &
                           grid2_centroid_lat(grid2_add))*norm_factor
          wts_map2(3,n) = (weights(num_wts+3) - weights(num_wts+1)* &
                           grid2_centroid_lon(grid2_add))*norm_factor
        end if
      end do
      write(stdout,*) 'Total number of links = ',num_links_map1
      where ( grid1_area /= d_zero ) grid1_frac = grid1_frac/grid1_area
      where ( grid2_area /= d_zero ) grid2_frac = grid2_frac/grid2_area
      !
      ! perform some error checking on final weights
      !
      grid2_centroid_lat = d_zero
      grid2_centroid_lon = d_zero
      do n = 1 , grid1_size
        if ( grid1_area(n) < -.01D0 ) then
          write(stderr,*) 'Grid 1 area error: ',n,grid1_area(n)
        end if
        if ( grid1_centroid_lat(n) < -halfpi-.01D0 .or. &
             grid1_centroid_lat(n) >  halfpi+.01D0 ) then
          write(stderr,*) 'Grid 1 centroid lat error: ',n,grid1_centroid_lat(n)
        end if
        grid1_centroid_lat(n) = d_zero
        grid1_centroid_lon(n) = d_zero
      end do
      do n = 1 , grid2_size
        if ( grid2_area(n) < -.01D0 ) then
          write(stderr,*) 'Grid 2 area error: ',n,grid2_area(n)
        end if
        if ( grid2_centroid_lat(n) < -halfpi-.01D0 .or. &
             grid2_centroid_lat(n) >  halfpi+.01D0 ) then
          write(stderr,*) 'Grid 2 centroid lat error: ',n,grid2_centroid_lat(n)
        end if
        grid2_centroid_lat(n) = d_zero
        grid2_centroid_lon(n) = d_zero
      end do
      do n = 1 , num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        if ( wts_map1(1,n) < -.01D0 ) then
          write(stderr,*) 'Map 1 weight < 0 ',grid1_add,grid2_add,wts_map1(1,n)
        end if
        if ( norm_opt /= norm_opt_none .and. wts_map1(1,n) > 1.01D0 ) then
          write(stderr,*) 'Map 1 weight > 1 ',grid1_add,grid2_add,wts_map1(1,n)
        end if
        grid2_centroid_lat(grid2_add) = &
          grid2_centroid_lat(grid2_add) + wts_map1(1,n)
        if ( num_maps > 1 ) then
          if ( wts_map2(1,n) < -.01D0 ) then
            write(stderr,*) &
              'Map 2 weight < 0 ',grid1_add,grid2_add,wts_map2(1,n)
          end if
          if ( norm_opt /= norm_opt_none .and. wts_map2(1,n) > 1.01D0 ) then
            write(stderr,*) &
              'Map 2 weight < 0 ',grid1_add,grid2_add,wts_map2(1,n)
          end if
          grid1_centroid_lat(grid1_add) = &
            grid1_centroid_lat(grid1_add) + wts_map2(1,n)
        end if
      end do
      do n = 1 , grid2_size
        select case(norm_opt)
          case (norm_opt_dstarea)
            norm_factor = grid2_frac(grid2_add)
          case (norm_opt_frcarea)
            norm_factor = d_one
          case (norm_opt_none)
            norm_factor = grid2_area(grid2_add)
        end select
        if ( abs(grid2_centroid_lat(grid2_add)-norm_factor) > .01D0 ) then
          write(stderr,*) 'Error: sum of wts for map1 ',grid2_add, &
                  grid2_centroid_lat(grid2_add),norm_factor
        end if
      end do
      if ( num_maps > 1 ) then
        do n = 1 , grid1_size
          select case(norm_opt)
            case (norm_opt_dstarea)
              norm_factor = grid1_frac(grid1_add)
            case (norm_opt_frcarea)
              norm_factor = d_one
            case (norm_opt_none)
              norm_factor = grid1_area(grid1_add)
          end select
          if ( abs(grid1_centroid_lat(grid1_add)-norm_factor) > .01D0 ) then
            write(stderr,*) 'Error: sum of wts for map2 ',grid1_add, &
                    grid1_centroid_lat(grid1_add),norm_factor
          end if
        end do
      end if
    end subroutine remap_conserv

    subroutine intersection(location,intrsct_lat,intrsct_lon,lcoinc, &
                            beglat, beglon, endlat, endlon, begseg,  &
                            lbegin, lrevers)
      implicit none
      !
      ! this routine finds the next intersection of a destination grid
      ! line with the line segment given by beglon, endlon, etc.
      ! a coincidence flag is returned if the segment is entirely
      ! coincident with an ocean grid line.  the cells in which to search
      ! for an intersection must have already been restricted in the
      ! calling routine.
      !
      ! flag for first integration along this segment
      logical , intent(in) :: lbegin
      ! flag whether segment integrated in reverse
      logical , intent(in) :: lrevers
      real(rk8) , intent(in) :: beglat , beglon ! beginning lat/lon for segment
      real(rk8) , intent(in) :: endlat , endlon ! ending lat/lon for segment
      ! begin lat/lon of full segment
      real(rk8) , dimension(2) , intent(inout) :: begseg

      ! address in destination array containing this segment
      integer(ik4) , intent(out) :: location
      ! flag segments which are entirely coincident with a grid line
      logical , intent(out) :: lcoinc
      ! lat/lon coords of next intersect.
      real(rk8) , intent(out) :: intrsct_lat , intrsct_lon

      integer(ik4) :: n , next_n , cell , srch_corners
      integer(ik4) , save :: last_loc ! save location when crossing threshold
      logical :: loutside             ! flags points outside grid
      ! flags segments crossing threshold bndy
      logical , save :: lthresh = .false.
      real(rk8) :: lon1 , lon2         ! local longitude variables for segment
      real(rk8) :: lat1 , lat2         ! local latitude  variables for segment
      real(rk8) :: grdlon1 , grdlon2   ! local longitude variables for grid cell
      real(rk8) :: grdlat1 , grdlat2   ! local latitude  variables for grid cell
      real(rk8) :: vec1_lat , vec1_lon ! vectors and cross products used
      real(rk8) :: vec2_lat , vec2_lon ! during grid search
      real(rk8) :: cross_product
      real(rk8) :: eps , offset        ! small offset away from intersect
      real(rk8) :: s1 , s2 , determ    ! variables used for linear solve to
      real(rk8) :: mat1 , mat2 , mat3
      real(rk8) :: mat4 , rhs1 , rhs2  ! find intersection
      ! lat/lon coords offset for next search
      real(rk8) , save :: intrsct_lat_off , intrsct_lon_off
      !
      ! initialize defaults, flags, etc.
      !
      location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon
      mat1 = d_zero
      mat2 = d_zero
      mat3 = d_zero
      if ( num_srch_cells == 0 ) return
      if ( beglat > north_thresh .or. beglat < south_thresh ) then
        if ( lthresh ) location = last_loc
        call pole_intersection(location,                               &
                               intrsct_lat,intrsct_lon,lcoinc,lthresh, &
                               beglat, beglon, endlat, endlon, begseg, lrevers)
        if ( lthresh ) then
          last_loc = location
          intrsct_lat_off = intrsct_lat
          intrsct_lon_off = intrsct_lon
        end if
        return
      end if
      loutside = .false.
      if ( lbegin ) then
        lat1 = beglat
        lon1 = beglon
      else
        lat1 = intrsct_lat_off
        lon1 = intrsct_lon_off
      end if
      lat2 = endlat
      lon2 = endlon
      if ( (lon2-lon1) > d_three*halfpi ) then
        lon2 = lon2 - twopi
      else if ( (lon2-lon1) < -d_three*halfpi ) then
        lon2 = lon2 + twopi
      end if
      s1 = d_zero
      !
      ! search for location of this segment in ocean grid using cross
      ! product method to determine whether a point is enclosed by a cell
      !
      ! call timer_start(12)
      srch_corners = size(srch_corner_lat,dim=1)
      srch_loop: &
      do
        !
        ! if last segment crossed threshold, use that location
        !
        if ( lthresh ) then
          do cell = 1 , num_srch_cells
            if ( srch_add(cell) == last_loc ) then
              location = last_loc
              eps = dlowval
              exit srch_loop
            end if
          end do
        end if
        !
        ! otherwise normal search algorithm
        !
        cell_loop: &
        do cell = 1 , num_srch_cells
          corner_loop: &
          do n = 1 , srch_corners
            next_n = mod(n,srch_corners) + 1
            !
            ! here we take the cross product of the vector making
            ! up each cell side with the vector formed by the vertex
            ! and search point.  if all the cross products are
            ! positive, the point is contained in the cell.
            !
            vec1_lat = srch_corner_lat(next_n,cell) - &
                       srch_corner_lat(n     ,cell)
            vec1_lon = srch_corner_lon(next_n,cell) - &
                       srch_corner_lon(n     ,cell)
            vec2_lat = lat1 - srch_corner_lat(n,cell)
            vec2_lon = lon1 - srch_corner_lon(n,cell)
            !
            ! if endpoint coincident with vertex, offset
            ! the endpoint
            !
            if ( vec2_lat == 0 .and. vec2_lon == 0 ) then
              lat1 = lat1 + 1.D-10*(lat2-lat1)
              lon1 = lon1 + 1.D-10*(lon2-lon1)
              vec2_lat = lat1 - srch_corner_lat(n,cell)
              vec2_lon = lon1 - srch_corner_lon(n,cell)
            end if
            !
            ! check for 0,2pi crossings
            !
            if ( vec1_lon >  mathpi ) then
              vec1_lon = vec1_lon - twopi
            else if ( vec1_lon < -mathpi ) then
              vec1_lon = vec1_lon + twopi
            end if
            if ( vec2_lon >  mathpi ) then
              vec2_lon = vec2_lon - twopi
            else if ( vec2_lon < -mathpi ) then
              vec2_lon = vec2_lon + twopi
            end if
            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
            !
            ! if the cross product for a side is zero, the point
            !   lies exactly on the side or the side is degenerate
            !   (zero length).  if degenerate, set the cross
            !   product to a positive number.  otherwise perform
            !   another cross product between the side and the
            !   segment itself.
            ! if this cross product is also zero, the line is
            !   coincident with the cell boundary - perform the
            !   dot product and only choose the cell if the dot
            !   product is positive (parallel vs anti-parallel).
            !
            if ( dabs(cross_product) < dlowval ) then
              if ( vec1_lat /= d_zero .or. vec1_lon /= d_zero ) then
                vec2_lat = lat2 - lat1
                vec2_lon = lon2 - lon1
                if (vec2_lon >  mathpi) then
                  vec2_lon = vec2_lon - twopi
                else if (vec2_lon < -mathpi) then
                  vec2_lon = vec2_lon + twopi
                end if
                cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
              else
                cross_product = d_one
              end if
              if ( dabs(cross_product) < dlowval ) then
                lcoinc = .true.
                cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat
                if ( lrevers ) cross_product = -cross_product
              end if
            end if
            !
            ! if cross product is less than zero, this cell
            ! doesn't work
            !
            if ( cross_product < d_zero ) exit corner_loop
          end do corner_loop
          !
          ! if cross products all positive, we found the location
          !
          if ( n > srch_corners ) then
            location = srch_add(cell)
            !
            ! if the beginning of this segment was outside the
            ! grid, invert the segment so the intersection found
            ! will be the first intersection with the grid
            !
            if ( loutside ) then
              lat2 = beglat
              lon2 = beglon
              location = 0
              eps  = -dlowval
            else
              eps  = dlowval
            end if
            exit srch_loop
          end if
          !
          ! otherwise move on to next cell
          !
        end do cell_loop
        !
        ! if still no cell found, the point lies outside the grid.
        !   take some baby steps along the segment to see if any
        !   part of the segment lies inside the grid.
        !
        loutside = .true.
        s1 = s1 + 0.001D0
        lat1 = beglat + s1*(endlat - beglat)
        lon1 = beglon + s1*(lon2   - beglon)
        !
        ! reached the end of the segment and still outside the grid
        ! return no intersection
        !
        if ( s1 >= d_one ) return
      end do srch_loop
      ! call timer_stop(12)
      !
      ! now that a cell is found, search for the next intersection.
      ! loop over sides of the cell to find intersection with side
      ! must check all sides for coincidences or intersections
      !
      ! call timer_start(13)
      intrsct_loop: &
      do n = 1 , srch_corners
        next_n = mod(n,srch_corners) + 1
        grdlon1 = srch_corner_lon(n     ,cell)
        grdlon2 = srch_corner_lon(next_n,cell)
        grdlat1 = srch_corner_lat(n     ,cell)
        grdlat2 = srch_corner_lat(next_n,cell)
        !
        ! set up linear system to solve for intersection
        !
        mat1 = lat2 - lat1
        mat2 = grdlat1 - grdlat2
        mat3 = lon2 - lon1
        mat4 = grdlon1 - grdlon2
        rhs1 = grdlat1 - lat1
        rhs2 = grdlon1 - lon1
        if (mat3 >  mathpi) then
          mat3 = mat3 - twopi
        else if (mat3 < -mathpi) then
          mat3 = mat3 + twopi
        end if
        if (mat4 >  mathpi) then
          mat4 = mat4 - twopi
        else if (mat4 < -mathpi) then
          mat4 = mat4 + twopi
        end if
        if (rhs2 >  mathpi) then
          rhs2 = rhs2 - twopi
        else if (rhs2 < -mathpi) then
          rhs2 = rhs2 + twopi
        end if
        determ = mat1*mat4 - mat2*mat3
        !
        ! if the determinant is zero, the segments are either
        !   parallel or coincident.  coincidences were detected
        !   above so do nothing.
        ! if the determinant is non-zero, solve for the linear
        !   parameters s for the intersection point on each line
        !   segment.
        ! if 0<s1,s2<1 then the segment intersects with this side.
        !   return the point of intersection (adding a small
        !   number so the intersection is off the grid line).
        !
        if ( abs(determ) > 1.D-30 ) then
          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ
          if ( s2 >= d_zero .and. s2 <= d_one .and. &
               s1 >  d_zero .and. s1 <= d_one) then
            !
            ! recompute intersection based on full segment
            ! so intersections are consistent for both sweeps
            !
            if ( .not. loutside ) then
              mat1 = lat2 - begseg(1)
              mat3 = lon2 - begseg(2)
              rhs1 = grdlat1 - begseg(1)
              rhs2 = grdlon1 - begseg(2)
            else
              mat1 = begseg(1) - endlat
              mat3 = begseg(2) - endlon
              rhs1 = grdlat1 - endlat
              rhs2 = grdlon1 - endlon
            end if
            if ( mat3 >  mathpi ) then
              mat3 = mat3 - twopi
            else if ( mat3 < -mathpi ) then
              mat3 = mat3 + twopi
            end if
            if ( rhs2 >  mathpi ) then
              rhs2 = rhs2 - twopi
            else if ( rhs2 < -mathpi ) then
              rhs2 = rhs2 + twopi
            end if
            determ = mat1*mat4 - mat2*mat3
            !
            ! sometimes due to roundoff, the previous
            ! determinant is non-zero, but the lines
            ! are actually coincident.  if this is the
            ! case, skip the rest.
            !
            if ( determ /= d_zero ) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ
              offset = s1 + eps/determ
              if (offset > d_one) offset = d_one
              if ( .not. loutside ) then
                intrsct_lat = begseg(1) + mat1*s1
                intrsct_lon = begseg(2) + mat3*s1
                intrsct_lat_off = begseg(1) + mat1*offset
                intrsct_lon_off = begseg(2) + mat3*offset
              else
                intrsct_lat = endlat + mat1*s1
                intrsct_lon = endlon + mat3*s1
                intrsct_lat_off = endlat + mat1*offset
                intrsct_lon_off = endlon + mat3*offset
              end if
              exit intrsct_loop
            end if
          end if
        end if
        !
        ! no intersection this side, move on to next side
        !
      end do intrsct_loop
      ! call timer_stop(13)
      !
      ! if the segment crosses a pole threshold, reset the intersection
      ! to be the threshold latitude.  only check if this was not a
      ! threshold segment since sometimes coordinate transform can end
      ! up on other side of threshold again.
      !
      if ( lthresh ) then
        if ( intrsct_lat < north_thresh .or. intrsct_lat > south_thresh ) &
          lthresh = .false.
      else if ( lat1 > d_zero .and. intrsct_lat > north_thresh ) then
        intrsct_lat = north_thresh + dlowval
        intrsct_lat_off = north_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      else if ( lat1 < d_zero .and. intrsct_lat < south_thresh ) then
        intrsct_lat = south_thresh - dlowval
        intrsct_lat_off = south_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      end if
    end subroutine intersection

    subroutine pole_intersection(location,                               &
                                 intrsct_lat,intrsct_lon,lcoinc,lthresh, &
                                 beglat, beglon, endlat, endlon, begseg, &
                                 lrevers)
      implicit none
      !
      ! this routine is identical to the intersection routine except
      ! that a coordinate transformation (using a Lambert azimuthal
      ! equivalent projection) is performed to treat polar cells more
      ! accurately.
      !
      ! beginning lat/lon endpoints for segment
      real(rk8) , intent(in) :: beglat , beglon
      ! ending    lat/lon endpoints for segment
      real(rk8) , intent(in) :: endlat , endlon
      ! begin lat/lon of full segment
      real(rk8) , dimension(2) , intent(inout) :: begseg
      ! flag true if segment integrated in reverse
      logical , intent(in) :: lrevers

      ! address in destination array containing this
      ! segment -- also may contain last location on entry
      integer(ik4) , intent(inout) :: location
      ! flag segment coincident with grid line
      logical , intent(out) :: lcoinc
      ! flag segment crossing threshold boundary
      logical , intent(inout) :: lthresh
      ! lat/lon coords of next intersect.
      real(rk8) , intent(out) :: intrsct_lat , intrsct_lon

      integer(ik4) :: n , next_n , cell , srch_corners
      logical :: loutside ! flags points outside grid
      real(rk8) :: pi4 , rns, & ! north/south conversion
           x1, x2,       & ! local x variables for segment
           y1, y2,       & ! local y variables for segment
           begx, begy,   & ! beginning x,y variables for segment
           endx, endy,   & ! beginning x,y variables for segment
           begsegx, begsegy,   & ! beginning x,y variables for segment
           grdx1, grdx2, & ! local x variables for grid cell
           grdy1, grdy2, & ! local y variables for grid cell
           vec1_y, vec1_x, & ! vectors and cross products used
           vec2_y, vec2_x, & ! during grid search
           cross_product, eps, & ! eps=small offset away from intersect
           s1, s2, determ,     & ! variables used for linear solve to
           mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection
      ! x,y of each corner of srch cells
      real(rk8) , dimension(:,:) , allocatable :: srch_corner_x
      real(rk8) , dimension(:,:) , allocatable :: srch_corner_y
      !
      ! save last intersection to avoid roundoff during coord
      ! transformation
      !
      logical , save :: luse_last = .false.
      real(rk8) , save :: intrsct_x, intrsct_y  ! x,y for intersection
      !
      ! variables necessary if segment manages to hit pole
      !
      ! count attempts to avoid pole
      integer(ik4) , save :: avoid_pole_count = 0
      ! endpoint offset to avoid pole
      real(rk8) , save :: avoid_pole_offset = dlowval
      !
      ! initialize defaults, flags, etc.
      !
      if ( .not. lthresh ) location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon
      loutside = .false.
      s1 = d_zero
      !
      ! convert coordinates
      !
      allocate(srch_corner_x(size(srch_corner_lat,DIM=1),  &
                             size(srch_corner_lat,DIM=2)), &
               srch_corner_y(size(srch_corner_lat,DIM=1),  &
                             size(srch_corner_lat,DIM=2)))
      if (beglat > d_zero) then
        pi4 = d_rfour*mathpi
        rns = d_one
      else
        pi4 = -d_rfour*mathpi
        rns = -d_one
      end if
      if ( luse_last ) then
        x1 = intrsct_x
        y1 = intrsct_y
      else
        x1 = rns*d_two*sin(pi4 - d_half*beglat)*cos(beglon)
        y1 =     d_two*sin(pi4 - d_half*beglat)*sin(beglon)
        luse_last = .true.
      end if
      x2 = rns*d_two*sin(pi4 - d_half*endlat)*cos(endlon)
      y2 =     d_two*sin(pi4 - d_half*endlat)*sin(endlon)
      srch_corner_x = rns*d_two*sin(pi4 - d_half*srch_corner_lat)* &
                                cos(srch_corner_lon)
      srch_corner_y =     d_two*sin(pi4 - d_half*srch_corner_lat)* &
                                sin(srch_corner_lon)
      begx = x1
      begy = y1
      endx = x2
      endy = y2
      begsegx = rns*d_two*sin(pi4 - d_half*begseg(1))*cos(begseg(2))
      begsegy =     d_two*sin(pi4 - d_half*begseg(1))*sin(begseg(2))
      intrsct_x = endx
      intrsct_y = endy
      !
      ! search for location of this segment in ocean grid using cross
      ! product method to determine whether a point is enclosed by a cell
      !
      ! call timer_start(12)
      srch_corners = size(srch_corner_lat,DIM=1)
      srch_loop: &
      do
        !
        ! if last segment crossed threshold, use that location
        !
        if ( lthresh ) then
          do cell = 1 , num_srch_cells
            if ( srch_add(cell) == location ) then
              eps = dlowval
              exit srch_loop
            end if
          end do
        end if
        !
        ! otherwise normal search algorithm
        !
        cell_loop: &
        do cell = 1 , num_srch_cells
          corner_loop: &
          do n = 1 , srch_corners
            next_n = mod(n,srch_corners) + 1
            !
            ! here we take the cross product of the vector making
            ! up each cell side with the vector formed by the vertex
            ! and search point.  if all the cross products are
            ! positive, the point is contained in the cell.
            !
            vec1_x = srch_corner_x(next_n,cell) - &
                     srch_corner_x(n     ,cell)
            vec1_y = srch_corner_y(next_n,cell) - &
                     srch_corner_y(n     ,cell)
            vec2_x = x1 - srch_corner_x(n,cell)
            vec2_y = y1 - srch_corner_y(n,cell)
            !
            ! if endpoint coincident with vertex, offset
            ! the endpoint
            !
            if ( vec2_x == 0 .and. vec2_y == 0 ) then
              x1 = x1 + 1.D-10*(x2-x1)
              y1 = y1 + 1.D-10*(y2-y1)
              vec2_x = x1 - srch_corner_x(n,cell)
              vec2_y = y1 - srch_corner_y(n,cell)
            end if
            cross_product = vec1_x*vec2_y - vec2_x*vec1_y
            !
            ! if the cross product for a side is zero, the point
            !   lies exactly on the side or the length of a side
            !   is zero.  if the length is zero set det > 0.
            !   otherwise, perform another cross
            !   product between the side and the segment itself.
            ! if this cross product is also zero, the line is
            !   coincident with the cell boundary - perform the
            !   dot product and only choose the cell if the dot
            !   product is positive (parallel vs anti-parallel).
            !

            if ( dabs(cross_product) < dlowval ) then
              if ( vec1_x /= d_zero .or. vec1_y /= 0 ) then
                vec2_x = x2 - x1
                vec2_y = y2 - y1
                cross_product = vec1_x*vec2_y - vec2_x*vec1_y
              else
                cross_product = d_one
              end if

              if ( dabs(cross_product) < dlowval ) then
                lcoinc = .true.
                cross_product = vec1_x*vec2_x + vec1_y*vec2_y
                if ( lrevers ) cross_product = -cross_product
              end if
            end if
            !
            ! if cross product is less than zero, this cell
            ! doesn't work
            !
            if ( cross_product < d_zero ) exit corner_loop
          end do corner_loop
          !
          ! if cross products all positive, we found the location
          !
          if ( n > srch_corners ) then
            location = srch_add(cell)
            !
            ! if the beginning of this segment was outside the
            ! grid, invert the segment so the intersection found
            ! will be the first intersection with the grid
            !
            if ( loutside ) then
              x2 = begx
              y2 = begy
              location = 0
              eps  = -dlowval
            else
              eps  = dlowval
            end if
            exit srch_loop
          end if
          !
          ! otherwise move on to next cell
          !
        end do cell_loop
        !
        ! if no cell found, the point lies outside the grid.
        !   take some baby steps along the segment to see if any
        !   part of the segment lies inside the grid.
        !
        loutside = .true.
        s1 = s1 + 0.001D0
        x1 = begx + s1*(x2 - begx)
        y1 = begy + s1*(y2 - begy)
        !
        ! reached the end of the segment and still outside the grid
        ! return no intersection
        !
        if ( s1 >= d_one ) then
          deallocate(srch_corner_x, srch_corner_y)
          luse_last = .false.
          return
        end if
      end do srch_loop
      !call timer_stop(12)
      !
      ! now that a cell is found, search for the next intersection.
      ! loop over sides of the cell to find intersection with side
      ! must check all sides for coincidences or intersections
      !
      ! call timer_start(13)
      intrsct_loop: &
      do n=1,srch_corners
        next_n = mod(n,srch_corners) + 1
        grdy1 = srch_corner_y(n     ,cell)
        grdy2 = srch_corner_y(next_n,cell)
        grdx1 = srch_corner_x(n     ,cell)
        grdx2 = srch_corner_x(next_n,cell)
        !
        ! set up linear system to solve for intersection
        !
        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1
        determ = mat1*mat4 - mat2*mat3
        !
        ! if the determinant is zero, the segments are either
        !   parallel or coincident or one segment has zero length.
        !   coincidences were detected above so do nothing.
        ! if the determinant is non-zero, solve for the linear
        !   parameters s for the intersection point on each line
        !   segment.
        ! if 0<s1,s2<1 then the segment intersects with this side.
        !   return the point of intersection (adding a small
        !   number so the intersection is off the grid line).
        !
        if ( abs(determ) > 1.D-30 ) then
          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ
          if ( s2 >= d_zero .and. s2 <= d_one .and. &
               s1 >  d_zero .and. s1 <= d_one ) then
            !
            ! recompute intersection using entire segment
            ! for consistency between sweeps
            !
            if ( .not. loutside ) then
              mat1 = x2 - begsegx
              mat3 = y2 - begsegy
              rhs1 = grdx1 - begsegx
              rhs2 = grdy1 - begsegy
            else
              mat1 = x2 - endx
              mat3 = y2 - endy
              rhs1 = grdx1 - endx
              rhs2 = grdy1 - endy
            end if
            determ = mat1*mat4 - mat2*mat3
            !
            ! sometimes due to roundoff, the previous
            ! determinant is non-zero, but the lines
            ! are actually coincident.  if this is the
            ! case, skip the rest.
            !
            if ( determ /= d_zero ) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ
              if ( .not. loutside ) then
                intrsct_x = begsegx + s1*mat1
                intrsct_y = begsegy + s1*mat3
              else
                intrsct_x = endx + s1*mat1
                intrsct_y = endy + s1*mat3
              end if
              !
              ! convert back to lat/lon coordinates
              !
              intrsct_lon = rns*atan2(intrsct_y,intrsct_x)
              if ( intrsct_lon < d_zero ) &
                intrsct_lon = intrsct_lon + twopi
              if ( abs(intrsct_x) > 1.D-10 ) then
                intrsct_lat = (pi4 - &
                  asin(rns*d_half*intrsct_x/cos(intrsct_lon)))*d_two
              else if ( abs(intrsct_y) > 1.D-10 ) then
                intrsct_lat = (pi4 - &
                  asin(d_half*intrsct_y/sin(intrsct_lon)))*d_two
              else
                intrsct_lat = d_two*pi4
              end if
              !
              ! add offset in transformed space for next pass.
              !
              if ( s1 - eps/determ < d_one ) then
                intrsct_x = intrsct_x - mat1*(eps/determ)
                intrsct_y = intrsct_y - mat3*(eps/determ)
              else
                if ( .not. loutside ) then
                  intrsct_x = endx
                  intrsct_y = endy
                  intrsct_lat = endlat
                  intrsct_lon = endlon
                else
                  intrsct_x = begsegx
                  intrsct_y = begsegy
                  intrsct_lat = begseg(1)
                  intrsct_lon = begseg(2)
                end if
              end if
              exit intrsct_loop
            end if
          end if
        end if
        !
        ! no intersection this side, move on to next side
        !
      end do intrsct_loop
      ! call timer_stop(13)
      deallocate(srch_corner_x, srch_corner_y)
      !
      ! if segment manages to cross over pole, shift the beginning
      ! endpoint in order to avoid hitting pole directly
      ! (it is ok for endpoint to be pole point)
      !
      if ( abs(intrsct_x) < 1.D-10 .and. abs(intrsct_y) < 1.D-10 .and. &
          (endx /= d_zero .and. endy /=0) ) then
        if ( avoid_pole_count > 2 ) then
          avoid_pole_count = 0
          avoid_pole_offset = 10.0D0*avoid_pole_offset
        end if
        cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx)
        intrsct_lat = begseg(1)
        if ( cross_product*intrsct_lat > d_zero ) then
          intrsct_lon = beglon    + avoid_pole_offset
          begseg(2)   = begseg(2) + avoid_pole_offset
        else
          intrsct_lon = beglon    - avoid_pole_offset
          begseg(2)   = begseg(2) - avoid_pole_offset
        end if
        avoid_pole_count = avoid_pole_count + 1
        luse_last = .false.
      else
        avoid_pole_count = 0
        avoid_pole_offset = dlowval
      end if
      !
      ! if the segment crosses a pole threshold, reset the intersection
      ! to be the threshold latitude and do not reuse x,y intersect
      ! on next entry.  only check if did not cross threshold last
      ! time - sometimes the coordinate transformation can place a
      ! segment on the other side of the threshold again
      !
      if ( lthresh ) then
        if ( intrsct_lat > north_thresh .or. intrsct_lat < south_thresh ) &
          lthresh = .false.
      else if ( beglat > d_zero .and. intrsct_lat < north_thresh ) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  mathpi) mat3 = mat3 - twopi
        if (mat3 < -mathpi) mat3 = mat3 + twopi
        intrsct_lat = north_thresh - dlowval
        s1 = (north_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      else if ( beglat < d_zero .and. intrsct_lat > south_thresh ) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  mathpi) mat3 = mat3 - twopi
        if (mat3 < -mathpi) mat3 = mat3 + twopi
        intrsct_lat = south_thresh + dlowval
        s1 = (south_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      end if
      !
      ! if reached end of segment, do not use x,y intersect
      ! on next entry
      !
      if ( intrsct_lat == endlat .and. intrsct_lon == endlon ) then
        luse_last = .false.
      end if
    end subroutine pole_intersection

    subroutine line_integral(weights, num_wts,                 &
                             in_phi1, in_phi2, theta1, theta2, &
                             grid1_lat, grid1_lon, grid2_lat, grid2_lon)
      implicit none
      !
      ! this routine computes the line integral of the flux function
      ! that results in the interpolation weights.  the line is defined
      ! by the input lat/lon of the endpoints.
      !
      integer(ik4) , intent(in) :: num_wts  ! number of weights to compute
      ! longitude endpoints for the segment
      real(rk8) , intent(in) :: in_phi1, in_phi2
      ! latitude  endpoints for the segment
      real(rk8) , intent(in) :: theta1 , theta2
      ! reference coordinates for each grid (to ensure correct 0,2pi interv.)
      real(rk8) , intent(in) :: grid1_lat , grid1_lon
      real(rk8) , intent(in) :: grid2_lat , grid2_lon

      ! line integral contribution to weights
      real(rk8) , dimension(2*num_wts) , intent(out) :: weights

      real(rk8) :: dphi , sinth1 , sinth2 , costh1 , costh2 , fac , &
                   phi1 , phi2
      real(rk8) :: f1 , f2 , fint
      !
      ! weights for the general case based on a trapezoidal approx to
      ! the integrals.
      !
      sinth1 = sin(theta1)
      sinth2 = sin(theta2)
      costh1 = cos(theta1)
      costh2 = cos(theta2)
      dphi = in_phi1 - in_phi2
      if (dphi >  mathpi) then
        dphi = dphi - twopi
      else if (dphi < -mathpi) then
        dphi = dphi + twopi
      end if
      dphi = d_half*dphi
      !
      ! the first weight is the area overlap integral. the second and
      ! fourth are second-order latitude gradient weights.
      !
      weights(        1) = dphi*(sinth1 + sinth2)
      weights(num_wts+1) = dphi*(sinth1 + sinth2)
      weights(        2) = dphi*(costh1 + costh2 + (theta1*sinth1 + &
                                                    theta2*sinth2))
      weights(num_wts+2) = dphi*(costh1 + costh2 + (theta1*sinth1 + &
                                                    theta2*sinth2))
      !
      ! the third and fifth weights are for the second-order phi gradient
      ! component.  must be careful of longitude range.
      !
      f1 = d_half*(costh1*sinth1 + theta1)
      f2 = d_half*(costh2*sinth2 + theta2)
      phi1 = in_phi1 - grid1_lon
      if (phi1 >  mathpi) then
        phi1 = phi1 - twopi
      else if (phi1 < -mathpi) then
        phi1 = phi1 + twopi
      end if
      phi2 = in_phi2 - grid1_lon
      if (phi2 >  mathpi) then
        phi2 = phi2 - twopi
      else if (phi2 < -mathpi) then
        phi2 = phi2 + twopi
      end if
      if ((phi2-phi1) <  mathpi .and. (phi2-phi1) > -mathpi) then
        weights(3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > d_zero) then
          fac = mathpi
        else
          fac = -mathpi
        end if
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(3) = d_half*phi1*(phi1-fac)*f1 - &
                     d_half*phi2*(phi2+fac)*f2 + &
                     d_half*fac*(phi1+phi2)*fint
      end if
      phi1 = in_phi1 - grid2_lon
      if (phi1 >  mathpi) then
        phi1 = phi1 - twopi
      else if (phi1 < -mathpi) then
        phi1 = phi1 + twopi
      end if
      phi2 = in_phi2 - grid2_lon
      if (phi2 >  mathpi) then
        phi2 = phi2 - twopi
      else if (phi2 < -mathpi) then
        phi2 = phi2 + twopi
      end if
      if ((phi2-phi1) <  mathpi .and. (phi2-phi1) > -mathpi) then
        weights(num_wts+3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > d_zero) then
          fac = mathpi
        else
          fac = -mathpi
        end if
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(num_wts+3) = d_half*phi1*(phi1-fac)*f1 - &
                             d_half*phi2*(phi2+fac)*f2 + &
                             d_half*fac*(phi1+phi2)*fint
      end if
    end subroutine line_integral

    subroutine store_link_cnsrv(add1, add2, weights)
      implicit none
      !
      ! this routine stores the address and weight for this link in
      ! the appropriate address and weight arrays and resizes those
      ! arrays if necessary.
      !
      integer(ik4) , intent(in) :: add1 , add2
      ! array of remapping weights for this link
      real(rk8) , dimension(:) , intent(in) :: weights
      integer(ik4) :: nlink , min_link , max_link ! link index
      ! min,max link add to restrict search
      integer(ik4) , dimension(:,:) , allocatable , save :: link_add1
      integer(ik4) , dimension(:,:) , allocatable , save :: link_add2
      logical , save :: first_call = .true.
      !
      ! if all weights are zero, do not bother storing the link
      !
      if ( all(dabs(weights) < dlowval) ) return
      !
      ! restrict the range of links to search for existing links
      !
      if ( first_call ) then
        allocate(link_add1(2,grid1_size), link_add2(2,grid2_size))
        link_add1 = 0
        link_add2 = 0
        first_call = .false.
        min_link = 1
        max_link = 0
      else
        min_link = min(link_add1(1,add1),link_add2(1,add2))
        max_link = max(link_add1(2,add1),link_add2(2,add2))
        if (min_link == 0) then
          min_link = 1
          max_link = 0
        end if
      end if
      !
      ! if the link already exists, add the weight to the current weight
      ! arrays
      !
      do nlink = min_link , max_link
        if ( add1 == grid1_add_map1(nlink) ) then
          if ( add2 == grid2_add_map1(nlink) ) then
            wts_map1(:,nlink) = wts_map1(:,nlink) + weights(1:num_wts)
            if ( num_maps == 2 ) then
              wts_map2(:,nlink) = wts_map2(:,nlink) + &
                                  weights(num_wts+1:2*num_wts)
            end if
            return
          end if
        end if
      end do
      !
      ! if the link does not yet exist, increment number of links and
      ! check to see if remap arrays need to be increased to accomodate
      ! the new link.  then store the link.
      !
      num_links_map1 = num_links_map1 + 1
      if ( num_links_map1 > max_links_map1 ) &
        call resize_remap_vars(1,resize_increment)
      grid1_add_map1(num_links_map1) = add1
      grid2_add_map1(num_links_map1) = add2
      wts_map1    (:,num_links_map1) = weights(1:num_wts)
      if ( num_maps > 1 ) then
        num_links_map2  = num_links_map2 + 1
        if ( num_links_map2 > max_links_map2 ) &
          call resize_remap_vars(2,resize_increment)
        grid1_add_map2(num_links_map2) = add1
        grid2_add_map2(num_links_map2) = add2
        wts_map2    (:,num_links_map2) = weights(num_wts+1:2*num_wts)
      end if
      if ( link_add1(1,add1) == 0 ) link_add1(1,add1) = num_links_map1
      if ( link_add2(1,add2) == 0 ) link_add2(1,add2) = num_links_map1
      link_add1(2,add1) = num_links_map1
      link_add2(2,add2) = num_links_map1
    end subroutine store_link_cnsrv

end module mod_scrip_remap_conserv
