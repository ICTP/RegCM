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

      module mod_slice

      use mod_regcm_param

      implicit none

#ifdef MPP1
      real(8) , dimension(ix,kx,-1:jxp+2,ntr) :: chib3d
      real(8) , dimension(ix,kx,jxp) :: pb3d , qsb3d , rhb3d , rhob3d , &
                                      & ubx3d , vbx3d
      real(8) , dimension(ix,kx,-1:jxp+2) :: qcb3d , qvb3d , tb3d ,     &
                                      & ubd3d , vbd3d
#else
      real(8) , dimension(ix,kx,jx,ntr) :: chib3d
      real(8) , dimension(ix,kx,jx) :: pb3d , qcb3d , qsb3d , qvb3d ,   &
                                     & rhb3d , rhob3d , tb3d , ubd3d ,  &
                                     & ubx3d , vbd3d , vbx3d
#endif

      end module mod_slice
