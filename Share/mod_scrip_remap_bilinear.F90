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
!***********************************************************************

module mod_scrip_remap_bilinear

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_scrip_grids
  use mod_scrip_remap_vars

  implicit none

  private

  ! max iteration count for i,j iteration
  integer(ik4) , parameter :: max_iter = 100

  real(rk8) , parameter :: converge = 1.D-10 ! convergence criterion

  public :: remap_bilin

  contains

    subroutine remap_bilin
      implicit none
      !
      ! this routine computes the weights for a bilinear interpolation.
      !
      integer(ik4) :: n , icount
      integer(ik4) :: dst_add ! destination address
      integer(ik4) :: iter    ! iteration counter
      integer(ik4) :: nmap    ! index of current map being computed

      ! address for the four source points
      integer(ik4) , dimension(4) :: src_add

      ! coordinates of four bilinear corners
      real(rk8) , dimension(4) :: src_lats , src_lons
      ! bilinear weights for four corners
      real(rk8) , dimension(4) :: wgts

      real(rk8) :: plat , plon        ! lat/lon coords of destination point
      real(rk8) :: iguess , jguess    ! current guess for bilinear coordinate
      real(rk8) :: deli , delj        ! corrections to i,j
      real(rk8) :: dth1 , dth2 , dth3 ! some latitude  differences
      real(rk8) :: dph1 , dph2 , dph3 ! some longitude differences
      real(rk8) :: dthp , dphp        ! difference between point and sw corner
      real(rk8) :: mat1 , mat2 , mat3 , mat4 ! matrix elements
      real(rk8) :: determinant , sum_wgts    ! matrix determinant , sum of
      !
      ! compute mappings from grid1 to grid2
      !
      nmap = 1
      if ( grid1_rank /= 2 ) then
        call die('remap_bilin', &
          'Can not do bilinear interpolation when grid1_rank /= 2',1)
      end if
      !
      ! loop over destination grid
      !
      grid_loop1: &
      do dst_add = 1 , grid2_size
        if ( .not. grid2_mask(dst_add) ) cycle grid_loop1
        plat = grid2_center_lat(dst_add)
        plon = grid2_center_lon(dst_add)
        !
        ! find nearest square of grid points on source grid
        !
        call grid_search_bilin(src_add, src_lats, src_lons,      &
                               plat, plon, grid1_dims,           &
                               grid1_center_lat, grid1_center_lon, &
                               grid1_bound_box, bin_addr1, bin_addr2)
        !
        ! check to see if points are land points
        !
        if ( src_add(1) > 0 ) then
          do n = 1 , 4
            if ( .not. grid1_mask(src_add(n)) ) src_add(1) = 0
          end do
        end if
        !
        ! if point found, find local i,j coordinates for weights
        !
        if ( src_add(1) > 0 ) then
          grid2_frac(dst_add) = d_one
          !
          ! iterate to find i,j for bilinear approximation
          !
          dth1 = src_lats(2) - src_lats(1)
          dth2 = src_lats(4) - src_lats(1)
          dth3 = src_lats(3) - src_lats(2) - dth2

          dph1 = src_lons(2) - src_lons(1)
          dph2 = src_lons(4) - src_lons(1)
          dph3 = src_lons(3) - src_lons(2)

          if ( dph1 >  d_three*halfpi ) dph1 = dph1 - twopi
          if ( dph2 >  d_three*halfpi ) dph2 = dph2 - twopi
          if ( dph3 >  d_three*halfpi ) dph3 = dph3 - twopi
          if ( dph1 < -d_three*halfpi ) dph1 = dph1 + twopi
          if ( dph2 < -d_three*halfpi ) dph2 = dph2 + twopi
          if ( dph3 < -d_three*halfpi ) dph3 = dph3 + twopi

          dph3 = dph3 - dph2

          iguess = d_half
          jguess = d_half

          iter_loop1: &
          do iter = 1 , max_iter
            dthp = plat - src_lats(1) - dth1*iguess - &
                   dth2*jguess - dth3*iguess*jguess
            dphp = plon - src_lons(1)

            if ( dphp >  d_three*halfpi ) dphp = dphp - twopi
            if ( dphp < -d_three*halfpi ) dphp = dphp + twopi

            dphp = dphp - dph1*iguess - dph2*jguess - &
                   dph3*iguess*jguess

            mat1 = dth1 + dth3*jguess
            mat2 = dth2 + dth3*iguess
            mat3 = dph1 + dph3*jguess
            mat4 = dph2 + dph3*iguess

            determinant = mat1*mat4 - mat2*mat3

            deli = (dthp*mat4 - mat2*dphp)/determinant
            delj = (mat1*dphp - dthp*mat3)/determinant

            if ( abs(deli) < converge .and. &
                 abs(delj) < converge ) exit iter_loop1

            iguess = iguess + deli
            jguess = jguess + delj
          end do iter_loop1
          if ( iter <= max_iter ) then
            !
            ! successfully found i,j - compute weights
            !
            wgts(1) = (d_one-iguess)*(d_one-jguess)
            wgts(2) = iguess*(d_one-jguess)
            wgts(3) = iguess*jguess
            wgts(4) = (d_one-iguess)*jguess
            call store_link_bilin(dst_add, src_add, wgts, nmap)
          else
            write(stderr,*) 'Point coords: ',plat,plon
            write(stderr,*) 'Dest grid lats: ',src_lats
            write(stderr,*) 'Dest grid lons: ',src_lons
            write(stderr,*) 'Dest grid addresses: ',src_add
            write(stderr,*) 'Current i,j : ',iguess, jguess
            call die('remap_bilin', &
              'Iteration for i,j exceed max iteration count',1)
          end if
        !
        ! search for bilinear failed - use a distance-weighted
        ! average instead (this is typically near the pole)
        !
        else if ( src_add(1) < 0 ) then
          src_add = abs(src_add)
          icount = 0
          do n = 1 , 4
            if ( grid1_mask(src_add(n)) ) then
              icount = icount + 1
            else
              src_lats(n) = d_zero
            end if
          end do
          if ( icount > 0 ) then
            !
            ! renormalize weights
            !
            sum_wgts = sum(src_lats)
            wgts(1) = src_lats(1)/sum_wgts
            wgts(2) = src_lats(2)/sum_wgts
            wgts(3) = src_lats(3)/sum_wgts
            wgts(4) = src_lats(4)/sum_wgts

            grid2_frac(dst_add) = d_one
            call store_link_bilin(dst_add, src_add, wgts, nmap)
          end if
        end if
      end do grid_loop1
      !
      ! compute mappings from grid2 to grid1 if necessary
      !
      if ( num_maps > 1 ) then
        nmap = 2
        if ( grid2_rank /= 2 ) then
          call die('remap_bilin', &
            'Can not do bilinear interpolation when grid2_rank /= 2',1)
        end if
        !
        ! loop over destination grid
        !
        grid_loop2: &
        do dst_add = 1 , grid1_size
          if ( .not. grid1_mask(dst_add) ) cycle grid_loop2
          plat = grid1_center_lat(dst_add)
          plon = grid1_center_lon(dst_add)
          !
          ! find nearest square of grid points on source grid
          !
          call grid_search_bilin(src_add, src_lats, src_lons,      &
                                 plat, plon, grid2_dims,           &
                                 grid2_center_lat, grid2_center_lon, &
                                 grid2_bound_box, bin_addr2, bin_addr1)
          !
          ! check to see if points are land points
          !
          if ( src_add(1) > 0 ) then
            do n = 1 , 4
              if ( .not. grid2_mask(src_add(n)) ) src_add(1) = 0
            end do
          end if
          !
          ! if point found, find i,j coordinates for weights
          !
          if ( src_add(1) > 0 ) then
            grid1_frac(dst_add) = d_one
            !
            ! iterate to find i,j for bilinear approximation
            !
            dth1 = src_lats(2) - src_lats(1)
            dth2 = src_lats(4) - src_lats(1)
            dth3 = src_lats(3) - src_lats(2) - dth2

            dph1 = src_lons(2) - src_lons(1)
            dph2 = src_lons(4) - src_lons(1)
            dph3 = src_lons(3) - src_lons(2)

            if ( dph1 >  mathpi ) dph1 = dph1 - twopi
            if ( dph2 >  mathpi ) dph2 = dph2 - twopi
            if ( dph3 >  mathpi ) dph3 = dph3 - twopi
            if ( dph1 < -mathpi ) dph1 = dph1 + twopi
            if ( dph2 < -mathpi ) dph2 = dph2 + twopi
            if ( dph3 < -mathpi ) dph3 = dph3 + twopi

            dph3 = dph3 - dph2

            iguess = d_zero
            jguess = d_zero

            iter_loop2: &
            do iter = 1 , max_iter
              dthp = plat - src_lats(1) - dth1*iguess - &
                     dth2*jguess - dth3*iguess*jguess
              dphp = plon - src_lons(1)

              if ( dphp >  mathpi ) dphp = dphp - twopi
              if ( dphp < -mathpi ) dphp = dphp + twopi

              dphp = dphp - dph1*iguess - dph2*jguess - &
                     dph3*iguess*jguess

              mat1 = dth1 + dth3*jguess
              mat2 = dth2 + dth3*iguess
              mat3 = dph1 + dph3*jguess
              mat4 = dph2 + dph3*iguess

              determinant = mat1*mat4 - mat2*mat3

              deli = (dthp*mat4 - mat2*dphp)/determinant
              delj = (mat1*dphp - dthp*mat3)/determinant

              if ( abs(deli) < converge .and. &
                   abs(delj) < converge) exit iter_loop2

              iguess = iguess + deli
              jguess = jguess + delj
            end do iter_loop2

            if ( iter <= max_iter ) then
              !
              ! successfully found i,j - compute weights
              !
              wgts(1) = (d_one-iguess)*(d_one-jguess)
              wgts(2) = iguess*(d_one-jguess)
              wgts(3) = iguess*jguess
              wgts(4) = (d_one-iguess)*jguess
              call store_link_bilin(dst_add, src_add, wgts, nmap)
            else
              write(stderr,*) 'Point coords: ',plat,plon
              write(stderr,*) 'Dest grid lats: ',src_lats
              write(stderr,*) 'Dest grid lons: ',src_lons
              write(stderr,*) 'Dest grid addresses: ',src_add
              write(stderr,*) 'Current i,j : ',iguess, jguess
              call die('remap_bilin', &
                'Iteration for i,j exceed max iteration count',1)
            end if
            !
            ! search for bilinear failed - us a distance-weighted
            ! average instead
            !
          else if ( src_add(1) < 0 ) then
            src_add = abs(src_add)
            icount = 0
            do n = 1 , 4
              if ( grid2_mask(src_add(n)) ) then
                icount = icount + 1
              else
                src_lats(n) = d_zero
              end if
            end do
            if ( icount > 0 ) then
              !
              ! renormalize weights
              !
              sum_wgts = sum(src_lats)
              wgts(1) = src_lats(1)/sum_wgts
              wgts(2) = src_lats(2)/sum_wgts
              wgts(3) = src_lats(3)/sum_wgts
              wgts(4) = src_lats(4)/sum_wgts
              grid1_frac(dst_add) = d_one
              call store_link_bilin(dst_add, src_add, wgts, nmap)
            end if
          end if
        end do grid_loop2
      end if ! nmap=2
    end subroutine remap_bilin

    subroutine grid_search_bilin(src_add, src_lats, src_lons,    &
                                 plat, plon, src_grid_dims,      &
                                 src_center_lat, src_center_lon, &
                                 src_grid_bound_box,             &
                                 src_bin_add, dst_bin_add)
      implicit none
      !
      ! this routine finds the location of the search point plat, plon
      ! in the source grid and returns the corners needed for a bilinear
      ! interpolation.
      !
      ! address of each corner point enclosing P
      integer(ik4) , dimension(4) , intent(out) :: src_add
      ! coordinates of the four corner points
      real(rk8) , dimension(4) , intent(out) :: src_lats , src_lons

      ! coordinates of the search point
      real(rk8) , intent(in) :: plat , plon
      ! size of each src grid dimension
      integer(ik4) , dimension(2) , intent(in) :: src_grid_dims
      ! coordinates of each src grid center
      real(rk8) , dimension(:) , intent(in) :: src_center_lat , src_center_lon
      ! bound box for source grid
      real(rk8) , dimension(:,:) , intent(in) :: src_grid_bound_box
      ! latitude bins for restricting searches
      integer(ik4) , dimension(:,:) , intent(in) :: src_bin_add , dst_bin_add

      integer(ik4) :: n , next_n , srch_add  ! dummy indices
      integer(ik4) :: nx , ny                ! dimensions of src grid
      integer(ik4) :: min_add , max_add      ! addresses for restricting search
      integer(ik4) :: i , j , jp1 , ip1
      integer(ik4) :: n_add , e_add , ne_add ! addresses

      ! vectors for cross-product check
      real(rk8) :: vec1_lat , vec1_lon , vec2_lat , vec2_lon
      real(rk8) :: cross_product , cross_product_last
      real(rk8) :: coslat_dst , sinlat_dst , coslon_dst , sinlon_dst
      real(rk8) :: dist_min , distance ! for computing dist-weighted avg
      !
      ! restrict search first using bins
      !
      src_add = 0
      min_add = size(src_center_lat)
      max_add = 1
      do n = 1 , num_srch_bins
        if ( plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n) .and. &
             plon >= bin_lons(1,n) .and. plon <= bin_lons(2,n) ) then
          min_add = min(min_add, src_bin_add(1,n))
          max_add = max(max_add, src_bin_add(2,n))
        end if
      end do
      !
      ! now perform a more detailed search
      !
      nx = src_grid_dims(1)
      ny = src_grid_dims(2)

      srch_loop: &
      do srch_add = min_add , max_add
        !
        ! first check bounding box
        !
        if ( plat <= src_grid_bound_box(2,srch_add) .and. &
             plat >= src_grid_bound_box(1,srch_add) .and. &
             plon <= src_grid_bound_box(4,srch_add) .and. &
             plon >= src_grid_bound_box(3,srch_add) ) then
          !
          ! we are within bounding box so get really serious
          !
          ! determine neighbor addresses
          !
          j = (srch_add - 1)/nx +1
          i = srch_add - (j-1)*nx
          if ( i < nx ) then
            ip1 = i + 1
          else
            ip1 = 1
          end if
          if ( j < ny ) then
            jp1 = j+1
          else
            jp1 = 1
          end if
          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1
          src_lats(1) = src_center_lat(srch_add)
          src_lats(2) = src_center_lat(e_add)
          src_lats(3) = src_center_lat(ne_add)
          src_lats(4) = src_center_lat(n_add)
          src_lons(1) = src_center_lon(srch_add)
          src_lons(2) = src_center_lon(e_add)
          src_lons(3) = src_center_lon(ne_add)
          src_lons(4) = src_center_lon(n_add)
          !
          ! for consistency, we must make sure all lons are in
          ! same 2pi interval
          !
          vec1_lon = src_lons(1) - plon
          if ( vec1_lon >  mathpi ) then
            src_lons(1) = src_lons(1) - twopi
          else if ( vec1_lon < -mathpi ) then
            src_lons(1) = src_lons(1) + twopi
          end if
          do n = 2 , 4
            vec1_lon = src_lons(n) - src_lons(1)
            if ( vec1_lon >  mathpi ) then
              src_lons(n) = src_lons(n) - twopi
            else if ( vec1_lon < -mathpi ) then
              src_lons(n) = src_lons(n) + twopi
            end if
          end do
          corner_loop: &
          do n = 1 , 4
            next_n = mod(n,4) + 1
            !
            ! here we take the cross product of the vector making
            ! up each box side with the vector formed by the vertex
            ! and search point.  if all the cross products are
            ! positive, the point is contained in the box.
            !
            vec1_lat = src_lats(next_n) - src_lats(n)
            vec1_lon = src_lons(next_n) - src_lons(n)
            vec2_lat = plat - src_lats(n)
            vec2_lon = plon - src_lons(n)
            !
            ! check for 0,2pi crossings
            !
            if ( vec1_lon >  d_three*halfpi ) then
              vec1_lon = vec1_lon - twopi
            else if ( vec1_lon < -d_three*halfpi ) then
              vec1_lon = vec1_lon + twopi
            end if
            if ( vec2_lon >  d_three*halfpi ) then
              vec2_lon = vec2_lon - twopi
            else if ( vec2_lon < -d_three*halfpi ) then
              vec2_lon = vec2_lon + twopi
            end if

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
            !
            ! if cross product is less than zero, this cell
            ! doesn't work
            !
            if ( n == 1 ) cross_product_last = cross_product
            if ( cross_product*cross_product_last < d_zero ) &
              exit corner_loop
            cross_product_last = cross_product
          end do corner_loop
          !
          ! if cross products all same sign, we found the location
          !
          if ( n > 4 ) then
            src_add(1) = srch_add
            src_add(2) = e_add
            src_add(3) = ne_add
            src_add(4) = n_add
            return
          end if
          !
          ! otherwise move on to next cell
          !
        end if !bounding box check
      end do srch_loop
      !
      ! if no cell found, point is likely either in a box that
      ! straddles either pole or is outside the grid.  fall back
      ! to a distance-weighted average of the four closest
      ! points.  go ahead and compute weights here, but store
      ! in src_lats and return -add to prevent the parent
      ! routine from computing bilinear weights
      !
      !print *,'Could not find location for ',plat,plon
      !print *,'Using nearest-neighbor average for this point'
      !
      coslat_dst = cos(plat)
      sinlat_dst = sin(plat)
      coslon_dst = cos(plon)
      sinlon_dst = sin(plon)

      dist_min = dhival
      src_lats = dhival
      do srch_add = min_add , max_add
        distance = acos(coslat_dst*cos(src_center_lat(srch_add))*   &
                        (coslon_dst*cos(src_center_lon(srch_add)) + &
                         sinlon_dst*sin(src_center_lon(srch_add)))+ &
                         sinlat_dst*sin(src_center_lat(srch_add)))
        if ( distance < dist_min ) then
          sort_loop: &
          do n = 1 , 4
            if ( distance < src_lats(n) ) then
              do i = 4 , n+1 , -1
                src_add (i) = src_add (i-1)
                src_lats(i) = src_lats(i-1)
              end do
              src_add (n) = -srch_add
              src_lats(n) = distance
              dist_min = src_lats(4)
              exit sort_loop
            end if
          end do sort_loop
        end if
      end do

      src_lons = d_one/(src_lats + dlowval)
      distance = sum(src_lons)
      src_lats = src_lons/distance
    end subroutine grid_search_bilin

    subroutine store_link_bilin(dst_add, src_add, weights, nmap)
      implicit none
      !
      ! this routine stores the address and weight for four links
      ! associated with one destination point in the appropriate address
      ! and weight arrays and resizes those arrays if necessary.
      !
      ! address on destination grid
      integer(ik4) , intent(in) :: dst_add
      ! identifies which direction for mapping
      integer(ik4) , intent(in) :: nmap
      ! addresses on source grid
      integer(ik4) , dimension(4) , intent(in) :: src_add
      ! array of remapping weights for these links
      real(rk8) , dimension(4) , intent(in) :: weights

      integer(ik4) :: n , num_links_old ! placeholder for old link number
      !
      ! increment number of links and check to see if remap arrays need
      ! to be increased to accomodate the new link.  then store the
      ! link.
      !
      select case (nmap)
        case(1)
          num_links_old  = num_links_map1
          num_links_map1 = num_links_old + 4
          if ( num_links_map1 > max_links_map1 )  &
            call resize_remap_vars(1,resize_increment)
          do n = 1 , 4
            grid1_add_map1(num_links_old+n) = src_add(n)
            grid2_add_map1(num_links_old+n) = dst_add
            wts_map1    (1,num_links_old+n) = weights(n)
          end do
        case(2)
          num_links_old  = num_links_map2
          num_links_map2 = num_links_old + 4
          if ( num_links_map2 > max_links_map2 ) &
            call resize_remap_vars(2,resize_increment)
          do n = 1 , 4
            grid1_add_map2(num_links_old+n) = dst_add
            grid2_add_map2(num_links_old+n) = src_add(n)
            wts_map2    (1,num_links_old+n) = weights(n)
          end do
      end select
    end subroutine store_link_bilin

end module mod_scrip_remap_bilinear
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
