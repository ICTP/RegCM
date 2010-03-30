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
!
      module mod_regcm_param

      implicit none
!
! PARAMETER definitions
!

!
!################### GRID DIMENSION ####################################
!

! Point in X (longitude) direction

      integer , parameter :: ix = 34

! Point in Y (latitude) direction

      integer , parameter :: jx = 48

! Point in vertical

      integer , parameter :: kx = 18

!###################### I/O control flag ###############################

! Create GrADS CTL files

      integer , parameter :: igrads = 1

! Machine endianess. LEAVE IT UNTOUCHED IF WANT TO EXCHANGE I/O FILES

      integer , parameter :: ibigend = 1

! Number of bytes in reclen. Usually 4

      integer , parameter :: ibyte = 4

!####################### MPI parameters ################################

#ifdef MPP1

! Number of processor used

      integer , parameter :: nproc = 16

! Point in Y (latitude) direction

      integer , parameter :: jxp = jx/nproc

#endif

! Sub grid decomposition

      integer , parameter :: nsg = 1

! Set amount of printout (still unused, sorry)

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

! Type of global analysis datasets used in Pre processing
!
! One in: ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
!         NRP2W,GFS11,FVGCM,FNEST,EH5OM
!

      character(5) , parameter :: dattyp = 'EIN15'

! SO4 Control Flag

      logical , parameter :: ehso4 = .false.

! Land Surface Legend type
!
! WARNING : SET TOGETHER lsmtyp and nveg. THEY MUST BE CONSISTENT.
!
! lsmtyp -> One in : BATS,USGS
! nveg   -> Set to 20 for BATS, 25 for USGS

      character(4) , parameter :: lsmtyp = 'BATS'
      integer , parameter :: nveg = 20

! Aerosol dataset used
!
! One in : AER00D0 -> Neither aerosol, nor dust used
!          AER01D0 -> Biomass, SO2 + BC + OC, no dust
!          AER10D0 -> Anthropogenic, SO2 + BC + OC, no dust
!          AER11D0 -> Anthropogenic+Biomass, SO2 + BC + OC, no dust
!          AER00D1 -> No aerosol, with dust
!          AER01D1 -> Biomass, SO2 + BC + OC, with dust
!          AER10D1 -> Anthropogenic, SO2 + BC + OC, with dust
!          AER11D1 -> Anthropogenic+Biomass, SO2 + BC + OC, with dust

      character(7) , parameter :: aertyp = 'AER00D0'

! Tracer parameters: number of tracers and bins number for dust

      integer , parameter :: ntr = 10
      integer , parameter :: nbin = 4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of configureation. Below this point things are
!    calculated from above or should be considered as fixed
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

      integer , parameter :: ixsg = ix * nsg
      integer , parameter :: jxsg = jx * nsg
      integer , parameter :: ixm1sg = (ix-1) * nsg
      integer , parameter :: jxm1sg = (jx-1) * nsg
      integer , parameter :: ixm2sg = (ix-2) * nsg
      integer , parameter :: jxm2sg = (jx-2) * nsg
#ifdef MPP1
      integer , parameter :: jxpsg = jxp * nsg
#endif

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
