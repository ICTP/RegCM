!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! this for AQCWA test 
!
      module mod_regcm_param
      implicit none
!
! PARAMETER definitions
!

! Point in X (longitude) direction

      integer , parameter :: ix = 120

! Point in Y (latitude) direction

      integer , parameter :: jx = 128

! Point in vertical

      integer , parameter :: kx = 18

#ifdef MPP1

! Number of processor used

      integer , parameter :: nproc = 16

! Point in Y (latitude) direction

      integer , parameter :: jxp = jx/nproc

#endif

! Sub grid decomposition

      integer , parameter :: nsg = 5

! Number of bytes in reclen. Usually 4

      integer , parameter :: ibyte = 4

! Set amount of printout

      integer , parameter :: debug_level = 1

! Buffer Zone Depth
! nspgx-1,nspgd-1 represent the number of cross/dot point slices
! on the boundary sponge or relaxation boundary conditions.

      integer , parameter :: nspgx = 12
      integer , parameter :: nspgd = 12

! Number od split exp modes

      integer , parameter :: nsplit = 2

! Number of lake points for lake model

      integer , parameter :: lkpts = 10

! BATS parameters

      character(5) , parameter :: dattyp = 'FNEST'

      logical , parameter :: ehso4 = .false.

      character(4) , parameter :: lsmtyp = 'BATS'

      character(7) , parameter :: aertyp = 'AER00D0'

! Tracer parameters: number of tracers and bins number for dust

      integer , parameter :: ntr = 10
      integer , parameter :: nbin = 4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of configureation. Below this point things are
!    calculated from above
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer , parameter :: ixm1 = ix - 1
      integer , parameter :: ixm2 = ix - 2
      integer , parameter :: ixm3 = ix - 3

      integer , parameter :: jxp1 = jx + 1
      integer , parameter :: jxm1 = jx - 1
      integer , parameter :: jxm2 = jx - 2

      integer , parameter :: kxm1 = kx - 1
      integer , parameter :: kxm2 = kx - 2
      integer , parameter :: kxp1 = kx + 1
      integer , parameter :: kxp2 = kx + 2
      integer , parameter :: kxp3 = kx + 3
      integer , parameter :: kxp4 = kx + 4

      integer , parameter :: nnsg = nsg*nsg
      integer , parameter :: nspgv = (nspgd+nspgx)*8 + 8
      integer , parameter :: nspgp = nspgx*4

#ifdef MPP1
      integer :: myid
      integer :: iwest , ieast , isouth , inorth
      integer :: jbegin , ibegin
      integer :: jendl , iendl
      integer :: jendx , iendx
      integer :: jendm , iendm
#endif

      integer , dimension(289276) :: mdatez

      integer , parameter :: numbat = 21 + 6
      integer , parameter :: numsub = 16
      integer , parameter :: nrad2d = 21
      integer , parameter :: nrad3d = 5

      end module mod_regcm_param
