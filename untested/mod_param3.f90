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
