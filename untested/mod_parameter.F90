      module regcm_param
      implicit none
!
! PARAMETER definitions
!

! Point in X (longitude) direction

      integer , parameter :: ix = 52


#ifdef MPP1

! Number of processor used

      integer , parameter :: nproc = 2

! Point in Y (latitude) direction

      integer , parameter :: mjx = 48

#else

! Point in Y (latitude) direction

      integer , parameter :: jx = 48

#endif

! Point in vertical

      integer , parameter :: kx = 18

! Sub grid decomposition

      integer , parameter :: nsg = 1
      integer , parameter :: nnsg = 1

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

      integer , parameter :: np1 = 2

      character(5) , parameter :: dattyp = 'ERAIN'

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

      integer , parameter :: ilx = ix - 1
      integer , parameter :: ilxm = ix - 2

#ifdef MPP1
      integer , parameter :: jxp = mjx/nproc
      integer , parameter :: jxbb = mjx - 1
#else
      integer , parameter :: jlx = jx - 1
      integer , parameter :: jlxm = jx - 2
#endif

      integer , parameter :: kxm = kx - 1
      integer , parameter :: kxp1 = kx + 1
      integer , parameter :: kxp2 = kx + 2

      integer , parameter :: nspgv = (nspgd+nspgx)*8 + 8
      integer , parameter :: nspgp = nspgx*4
      integer , parameter :: nbmax = ix - 1

#ifdef MPP1
      integer :: myid , iwest , ieast
      integer :: jbegin, jendl, jendx, jendm
#endif

      integer , dimension(289276) :: mdatez

      integer , parameter :: numbat = 21 + 6
      integer , parameter :: numsub = 16

      end module regcm_param
