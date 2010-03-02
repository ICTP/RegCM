      module mod_block
      implicit none
      integer :: nobs
      integer , parameter :: iter = 2400 , jter = 2400
      integer , parameter :: iblk = (iter*jter)/2
      real(4) , dimension(iblk,2) :: clay , sand
      real(4) , dimension(iblk) :: ht , ht2 , htsd
      real(4) , dimension(iblk) :: xobs , yobs
      real(8) , dimension(iter,jter) :: lnd8
      end module mod_block
