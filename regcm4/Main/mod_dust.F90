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

      module mod_dust

      use mod_dynparam

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
#ifdef MPP1
      real(8) , allocatable, dimension(:,:,:) :: clay2row2 , sand2row2 ,       &
           & silt2row2
      real(8) ,allocatable,  dimension(:,:) :: clayrow2 , dustsotex , sandrow2
      real(8) ,allocatable,  dimension(:,:,:) :: srel2d
#else
      real(8) , dimension(iy,nats,jx) :: clay2row2 , sand2row2 ,        &
           & silt2row2
      real(8) , dimension(iy,jx) :: clayrow2 , dustsotex , sandrow2
      real(8) , dimension(iy,jx,nsoil) :: srel2d
#endif
      real(8) , dimension(nsoil) :: dp

contains 
        subroutine allocate_mod_dust 

#ifdef MPP1
        allocate(clay2row2(iy,nats,jxp) )
        allocate(sand2row2(iy,nats,jxp) )
        allocate(silt2row2(iy,nats,jxp) )

        allocate(clayrow2(iy,jxp) )
        allocate(dustsotex(iy,jxp) )
        allocate(sandrow2(iy,jxp) )

        allocate(srel2d(iy,jxp,nsoil))
        
#endif 

        end subroutine allocate_mod_dust

      end module mod_dust
