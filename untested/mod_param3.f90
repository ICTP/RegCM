      module mod_param3

      use mod_regcm_param

      implicit none
!
! COMMON /IPARAM3/
!
      integer :: ispgd , ispgx , julday , julian , jxsex , k700 , kchi ,&
               & kclo , kcmd , kt , kxout , mdate , mdate0 , moutdate , &
               & ncld , ntimax
!
! COMMON /PARAM3/
!
      real(8) , dimension(kx) :: a , anudg , dsigma , qcon
      real(8) :: akht1 , akht2 , alpha , beta , cd , cdsea , ch ,       &
               & chsea , cp , declin , dectim , degrad , deltmx , dpd , &
               & eomeg , g , gmt , gnu , gnuhf , karman , omu , omuhf , &
               & ptop , ptop4 , r , rhos , rovcp , rovg , solcon ,      &
               & stbolt , tauht , thrlh1 , thrlh2
      real(8) , dimension(kxp1) :: sigma
      real(8) , dimension(kx,2) :: twt
      real(8) , dimension(nspgd) :: wgtd
      real(8) , dimension(nspgx) :: wgtx
      end module mod_param3
