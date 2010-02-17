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

      module mod_radbuf

      use mod_regcm_param

      implicit none

      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity

#ifdef MPP1
      real(8) , dimension(ixm1,kx,4,jxp) :: absnxt
      real(8) , dimension(ixm1,kxp1,kx + 1,jxp) :: abstot
      real(8) , dimension(ixm1,kxp1,jxp) :: emstot
#else
      real(8) , dimension(ixm1,kx,4,jxm1) :: absnxt
      real(8) , dimension(ixm1,kxp1,kx + 1,jxm1) :: abstot
      real(8) , dimension(ixm1,kxp1,jxm1) :: emstot
#endif

      end module mod_radbuf

