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

      module mod_slice

      use mod_dynparam

      implicit none

      real(8) ,allocatable, dimension(:,:,:,:) :: chib3d
      real(8) ,allocatable, dimension(:,:,:) :: pb3d , qsb3d , rhb3d ,  &
                                      & rhob3d , ubx3d , vbx3d
      real(8) ,allocatable, dimension(:,:,:) :: qcb3d , qvb3d , tb3d ,  &
                                      & ubd3d , vbd3d
      contains 
        subroutine allocate_mod_slice
        implicit none   
#ifdef MPP1
        allocate(chib3d(iy,kz,-1:jxp+2,ntr))      
        allocate(pb3d(iy,kz,jxp))
        allocate(qsb3d(iy,kz,jxp))
        allocate(rhb3d(iy,kz,jxp))
        allocate(rhob3d(iy,kz,jxp))
        allocate(ubx3d(iy,kz,jxp))
        allocate(vbx3d(iy,kz,jxp))
        allocate(qcb3d(iy,kz,-1:jxp+2))
        allocate(qvb3d(iy,kz,-1:jxp+2))
        allocate(tb3d(iy,kz,-1:jxp+2))
        allocate(ubd3d(iy,kz,-1:jxp+2))
        allocate(vbd3d(iy,kz,-1:jxp+2))
#else
        allocate(chib3d(iy,kz,jx,ntr))      
        allocate(pb3d(iy,kz,jx))
        allocate(qsb3d(iy,kz,jx))
        allocate(rhb3d(iy,kz,jx))
        allocate(rhob3d(iy,kz,jx))
        allocate(ubx3d(iy,kz,jx))
        allocate(vbx3d(iy,kz,jx))
        allocate(qcb3d(iy,kz,jx))
        allocate(qvb3d(iy,kz,jx))
        allocate(tb3d(iy,kz,jx))
        allocate(ubd3d(iy,kz,jx))
        allocate(vbd3d(iy,kz,jx))

#endif

       end subroutine allocate_mod_slice

      end module mod_slice
