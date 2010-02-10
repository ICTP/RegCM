      module mod_comtim
      implicit none
!
! COMMON /COMTIM/
!
      real(8) :: calday , dtime , twodt
      logical :: doabsems , dolw , dosw
      integer :: mbdate , mbsec , mcdate , mcsec , mdbase , mdcur ,     &
               & msbase , mscur , nelapse , nestep , nnbdat , nnbsec ,  &
               & nndbas , nnsbas , nrstrt , nstep , nstepr , nstop
      end module mod_comtim
