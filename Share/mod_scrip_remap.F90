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

module mod_scrip_remap

  use mod_realkinds
  use mod_intkinds
  use mod_constants

  private

  public :: remap

  contains

    !-----------------------------------------------------------
    !
    ! Performs the remapping based on weights computed elsewhere
    !
    !-----------------------------------------------------------

    subroutine remap(dst_array, map_wts, dst_add, src_add, &
                     src_array, src_grad1, src_grad2, src_grad3)
      implicit none
      integer(ik4), dimension(:), intent(in) :: &
           dst_add, &   ! destination address for each link
           src_add      ! source      address for each link
      real(rk8), dimension(:,:), intent(in) :: &
           map_wts      ! remapping weights for each link
      real(rk8), dimension(:), intent(in) ::   &
           src_array    ! array with source field to be remapped
      real(rk8), dimension(:), intent(in), optional :: &
           src_grad1, & ! gradient arrays on source grid necessary for
           src_grad2, & ! higher-order remappings
           src_grad3
      real(rk8), dimension(:), intent(inout) :: &
           dst_array    ! array for remapped field on destination grid
      integer(ik4) :: n , iorder
      !
      ! Check interpolation order
      !
      if ( present(src_grad1) ) then
        iorder = 2
      else
        iorder = 1
      endif
      !
      ! first or second order remapping ?
      !
      dst_array = d_zero
      select case (iorder)
        case(1)
          do n = 1 , size(dst_add)
            dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
                                    src_array(src_add(n))*map_wts(1,n)
          end do
        case(2)
          if ( size(map_wts,dim=1) == 3 ) then
            do n = 1 , size(dst_add)
              dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
                                      src_array(src_add(n))*map_wts(1,n) + &
                                      src_grad1(src_add(n))*map_wts(2,n) + &
                                      src_grad2(src_add(n))*map_wts(3,n)
            end do
          else if ( size(map_wts,dim=1) == 4 ) then
            do n = 1 , size(dst_add)
              dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
                                      src_array(src_add(n))*map_wts(1,n) + &
                                      src_grad1(src_add(n))*map_wts(2,n) + &
                                      src_grad2(src_add(n))*map_wts(3,n) + &
                                      src_grad3(src_add(n))*map_wts(4,n)
            end do
          endif
      end select
    end subroutine remap

end module mod_scrip_remap
