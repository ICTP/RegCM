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

      module mod_tmpsav
      use mod_regcm_param
      implicit none

#ifdef MPP1
      real(8) , dimension(iy,kz*4+2,jxp) :: sav0
      real(8) , dimension(iy,kz*4+2,jx) :: sav_0
      real(8) , dimension(iy,kz+nnsg+5,jxp) :: sav0a
      real(8) , dimension(iy,kz+nnsg+5,jx) :: sav_0a
      real(8) , dimension(iy,kzp1,jxp) :: sav0b
      real(8) , dimension(iy,kzp1,jx) :: sav_0b
      real(8) , dimension(iy,kz*2,jxp) :: sav0c
      real(8) , dimension(iy,kz*2,jx) :: sav_0c
      real(8) , dimension(iy,kz,jxp) :: sav0s
      real(8) , dimension(iy,kz,jx) :: sav_0s
      real(8) , dimension(iy,nsplit*2,jxp) :: sav0d
      real(8) , dimension(iy,nsplit*2,jx) :: sav_0d
      real(8) , dimension(iym1,kz*4+(kzp1*kzp2),jxp) :: sav1
      real(8) , dimension(iym1,kz*4+(kzp1*kzp2),jx) :: sav_1
      real(8) , dimension(iym1,nnsg*4+4,jxp) :: sav2
      real(8) , dimension(iym1,nnsg*4+4,jx) :: sav_2
      real(8) , dimension(iym1,nnsg*5+1,jxp) :: sav2a
      real(8) , dimension(iym1,nnsg*5+1,jx) :: sav_2a
      real(8) , dimension(iy,ntr*(kz*4+1),jxp) :: sav4
      real(8) , dimension(iy,ntr*(kz*4+1),jx) :: sav_4
      real(8) , dimension(iym1,7,jxp) :: sav4a
      real(8) , dimension(iym1,7,jx) :: sav_4a
      real(8) , dimension(kz,8,jxp) :: sav6
      real(8) , dimension(kz,8,jx) :: sav_6
#ifdef CLM
      real(8) , dimension(iym1,9,jx) :: sav_clmout
      real(8) , dimension(iym1,9,jxp) :: sav_clmin
#endif

#endif

!---------- DATA init section--------------------------------------------

      end module mod_tmpsav
