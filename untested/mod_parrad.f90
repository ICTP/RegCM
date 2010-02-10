      module mod_parrad

      use mod_regcm_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: plon = ix - 1
      integer , parameter :: plev = kx
      integer , parameter :: plat = 1
      integer , parameter :: pcnst = 1
      integer , parameter :: plevp = plev + 1
      integer , parameter :: plond = plon
      integer , parameter :: platd = plat
      integer , parameter :: plevd = plev*(3+pcnst)
      integer , parameter :: plevr = kx
      integer , parameter :: plevrp = plevr + 1

      end module mod_parrad
