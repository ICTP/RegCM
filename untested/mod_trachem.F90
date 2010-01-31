      module trachem

      use regcm_param

      implicit none
!
! COMMON /CHACHEM/
!
      character(5) , dimension(ntr) :: chtrname
!
! COMMON /CHTRPROP/
!
      real(8) , dimension(ntr,2) :: chtrdpv
      real(8) , dimension(nbin,2) :: chtrsize , dustbsiz
      real(8) , dimension(ntr) :: chtrsol
!
! COMMON /IEMOVALT/
!
      integer :: ichcumtra , ichdrdepo , ichremcvc , ichremlsc ,        &
               & ichsursrc
#ifdef MPP1
      integer , dimension(ix,jxp) :: icumbot , icumdwd , icumtop
#else
      integer , dimension(ix,jx) :: icumbot , icumdwd , icumtop
#endif
!
! COMMON /IOXYDANT/
!
      integer :: ibchb , ibchl , iochb , iochl , iso2 , iso4 , mixtype
      integer , dimension(nbin) :: idust
!
! COMMON /MASSFLX/
!
      real(8) , dimension(ix,2) :: mflx
!
! COMMON /OUTAER/
!
#ifdef MPP1
      real(8) , dimension(ix-1,kx,jxp) :: aerasp , aerext , aerssa
      real(8) , dimension(ix-1,jxp) :: aersrrf , aertarf
#else
      real(8) , dimension(ix-1,kx,jx-1) :: aerasp , aerext , aerssa
      real(8) , dimension(ix-1,jx-1) :: aersrrf , aertarf
#endif
!
! COMMON /REMOVALT/
!
#ifdef MPP1
      real(8) , dimension(ix,jxp,ntr) :: cemtr , cemtrac , remdrd
      real(8) , dimension(ix,kx) :: rembc , remrat
      real(8) , dimension(ix,kx,jxp,ntr) :: remcvc , remlsc , rxsaq1 ,  &
           & rxsaq2 , rxsg
#else
      real(8) , dimension(ix,jx,ntr) :: cemtr , cemtrac , remdrd
      real(8) , dimension(ix,kx) :: rembc , remrat
      real(8) , dimension(ix,kx,jx,ntr) :: remcvc , remlsc , rxsaq1 ,   &
           & rxsaq2 , rxsg
#endif
      end module trachem
