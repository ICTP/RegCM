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

      module mod_bxq

      use mod_regcm_param

      implicit none
!
#ifdef MPP1
      real(8) , dimension(iy,jxp,nsplit) :: ddsum
      real(8) , dimension(iy,jxp,nsplit,3) :: deld
      real(8) , dimension(iy,0:jxp,nsplit,3) :: delh
      real(8) , dimension(iy,0:jxp,nsplit) :: dhsum
      real(8) , dimension(iy,jxp) :: psdot
      real(8) , dimension(iy,jxp,3) :: work
!
      real(8) , dimension(iy,jxp+1) :: uu , vv
!
      real(8) , dimension(iy,kz,jxp+1) :: uuu , vvv
#else
      real(8) , dimension(iy,jx,nsplit) :: ddsum , dhsum
      real(8) , dimension(iy,jx,nsplit,3) :: deld , delh
      real(8) , dimension(iy,jx) :: psdot
      real(8) , dimension(iy,jx,3) :: work
!
      real(8) , dimension(iy,jx) :: uu , vv
!
      real(8) , dimension(iy,kz,jx) :: uuu , vvv
#endif

      end module mod_bxq
