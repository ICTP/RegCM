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

      module mod_dust

      use mod_regcm_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nsoil = 152
      integer , parameter :: nats = 12
      integer , parameter :: mode = 3
      integer , parameter :: jsoilm = 1
      integer , parameter :: jfs = 0
      integer , parameter :: ust = 1
!
! COMMON /DUST/
!
#ifdef MPP1
      real(8) , dimension(ix,nats,jxp) :: clay2row2 , sand2row2 ,       &
           & silt2row2
      real(8) , dimension(ix,jxp) :: clayrow2 , dustsotex , sandrow2
      real(8) , dimension(ix,jxp,nsoil) :: srel2d
#else
      real(8) , dimension(ix,nats,jx) :: clay2row2 , sand2row2 ,        &
           & silt2row2
      real(8) , dimension(ix,jx) :: clayrow2 , dustsotex , sandrow2
      real(8) , dimension(ix,jx,nsoil) :: srel2d
#endif
      real(8) , dimension(nsoil) :: dp

      end module mod_dust
