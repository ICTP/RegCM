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

      module mod_blh_tmp
      use mod_regcm_param
      implicit none
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: cgh , kvc , kvh , kvm , kvq
      real(8) , dimension(ix,jxp) :: hfxv , obklen , th10 , ustr ,      &
                                   & xhfx , xqfx
#else
      real(8) , dimension(ix,kx,jxm1) :: cgh , kvc , kvh , kvm , kvq
      real(8) , dimension(ix,jx) :: hfxv , obklen , th10 , ustr ,      &
                                   & xhfx , xqfx
#endif
      end module mod_blh_tmp
