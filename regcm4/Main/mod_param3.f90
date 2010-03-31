!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_param3

      use mod_regcm_param

      implicit none
!
      integer :: ispgd , ispgx , jxsex , k700 , kchi , kclo , kcmd ,    &
               & kt , kxout , ncld
!
      real(8) :: ptop , ptop4
      real(8) :: akht1 , akht2

      real(8) , dimension(kz) :: a , anudg , dsigma , qcon
      real(8) , dimension(kzp1) :: sigma
      real(8) , dimension(kz,2) :: twt
      real(8) , dimension(nspgd) :: wgtd
      real(8) , dimension(nspgx) :: wgtx

      end module mod_param3
