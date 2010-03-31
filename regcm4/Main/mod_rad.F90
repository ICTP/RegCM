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

      module mod_rad

      use mod_regcm_param

      implicit none

#ifdef MPP1
      real(8) , dimension(iym1,kz) :: cldfra , cldlwc
      real(8) , dimension(iym1,kz,jxp) :: heatrt
      real(8) , dimension(iym1,kzp1,jxp) :: o3prof
#else
      real(8) , dimension(iym1,kz) :: cldfra , cldlwc
      real(8) , dimension(iym1,kz,jxm1) :: heatrt
      real(8) , dimension(iym1,kzp1,jxm1) :: o3prof
#endif

      end module mod_rad
