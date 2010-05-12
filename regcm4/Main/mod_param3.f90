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

      use mod_dynparam

      implicit none
!
      integer :: ispgd , ispgx , k700 , kchi , kclo , kcmd , kt , ncld
!
      real(8) :: ptop4 , r8pt
      real(8) :: akht1 , akht2

      real(8) , allocatable , dimension(:) :: a , anudg , dsigma , qcon
      real(8) , allocatable , dimension(:) :: sigma
      real(8) , allocatable , dimension(:,:) :: twt
      real(8) , allocatable , dimension(:) :: wgtd
      real(8) , allocatable , dimension(:) :: wgtx

      contains

      subroutine allocate_mod_param3
      implicit none
        allocate(a(kz))
        allocate(anudg(kz))
        allocate(dsigma(kz))
        allocate(qcon(kz))
        allocate(sigma(kzp1))
        allocate(twt(kz,2))
        allocate(wgtd(nspgd))
        allocate(wgtx(nspgx))
      end subroutine allocate_mod_param3

      end module mod_param3
