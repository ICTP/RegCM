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

#ifdef MPP1
!
! COMMON /TRACEIO/
!
      real(8) , dimension(ix-1,kx,mjx-1) :: aerasp_io , aerext_io ,     &
           & aerssa_io
      real(8) , dimension(ix-1,mjx-1) :: aersrrf_io , aertarf_io
      real(8) , dimension(ix,mjx,ntr) :: cemtrac_io , cemtr_io ,        &
           & wxaq_io , wxsg_io
      real(8) , dimension(ix,mjx) :: dustsotex_io
      real(8) , dimension(ix,kx,mjx,ntr) :: rxsaq1_io , rxsaq2_io ,     &
           & rxsg_io
!
! COMMON /TRACHEMIO/
!
      real(8) , dimension(ix,kx,mjx,ntr) :: remcvc_io , remlsc_io
      real(8) , dimension(ix,mjx,ntr) :: remdrd_io
#endif

      end module trachem
