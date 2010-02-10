      module mod_comozp


      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: pnoz = 100 , pozlon = 1
!
! COMMON /COMOZP/
!
      real(8) :: cplol , cplos , ldoyoz , ndoyoz
      real(8) , dimension(1,pnoz) :: ozmix
      real(8) , dimension(pozlon,pnoz,1,2) :: ozmixm
      real(8) , dimension(pnoz) :: pin
!
! COMMON /IOZNCT/
!
      integer :: koz , nyroz
      end module mod_comozp
