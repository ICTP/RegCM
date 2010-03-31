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

      module mod_mainchem

      use mod_regcm_param

      implicit none
!
#ifdef MPP1
      real(8) , dimension(iy,jxp,12,ntr) :: chemsrc
      real(8) , dimension(iy,kz,-1:jxp+2,ntr) :: chia , chib
      real(8) , dimension(iy,jxp,ntr) :: srclp2
#else
      real(8) , dimension(iy,jx,12,ntr) :: chemsrc
      real(8) , dimension(iy,kz,jx,ntr) :: chia , chib
      real(8) , dimension(iy,jx,ntr) :: srclp2
#endif
!
#ifdef MPP1
      real(8) , dimension(iy,jxp,ntr) :: ddsfc , dtrace , wdcvc ,       &
                                       & wdlsc , wxaq , wxsg
#else
      real(8) , dimension(iy,jx,ntr) :: ddsfc , dtrace , wdcvc , wdlsc ,&
                                      & wxaq , wxsg
#endif

      real(4) , dimension(jxm2,iym2) :: fchem

#ifdef MPP1
      real(8) , dimension(iy,12,ntr,jxp) :: src0
      real(8) , dimension(iy,12,ntr,jx) :: src_0
#endif

      end module mod_mainchem
