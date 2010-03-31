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

      module mod_outrad
      use mod_regcm_param
      implicit none

#ifdef MPP1
      real(4) , dimension(jxp,iym2,nrad2d) :: frad2d
      real(4) , dimension(jxp,iym2,kz,nrad3d) :: frad3d
#else
      real(4) , dimension(jxm2,iym2,nrad2d) :: frad2d
      real(4) , dimension(jxm2,iym2,kz,nrad3d) :: frad3d
#endif

      end module mod_outrad
