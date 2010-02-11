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

      module mod_outrad
      use mod_regcm_param
      implicit none

#ifdef MPP1
      real(4) , dimension(jxp,ixm2,nrad2d) :: frad2d
      real(4) , dimension(jxp,ixm2,kx,nrad3d) :: frad3d
#else
      real(4) , dimension(jxm2,ixm2,nrad2d) :: frad2d
      real(4) , dimension(jxm2,ixm2,kx,nrad3d) :: frad3d
#endif

      end module mod_outrad
