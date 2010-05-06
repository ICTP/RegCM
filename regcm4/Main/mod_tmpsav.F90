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
      use mod_dynparam
      implicit none

#ifdef MPP1
      real(8) , allocatable , dimension(:,:,:) :: sav0
      real(8) , allocatable , dimension(:,:,:) :: sav_0
      real(8) , allocatable , dimension(:,:,:) :: sav0a
      real(8) , allocatable , dimension(:,:,:) :: sav_0a
      real(8) , allocatable , dimension(:,:,:) :: sav0b
      real(8) , allocatable , dimension(:,:,:) :: sav_0b
      real(8) , allocatable , dimension(:,:,:) :: sav0c
      real(8) , allocatable , dimension(:,:,:) :: sav_0c
      real(8) , allocatable , dimension(:,:,:) :: sav0s
      real(8) , allocatable , dimension(:,:,:) :: sav_0s
      real(8) , allocatable , dimension(:,:,:) :: sav0d
      real(8) , allocatable , dimension(:,:,:) :: sav_0d
      real(8) , allocatable , dimension(:,:,:) :: sav1
      real(8) , allocatable , dimension(:,:,:) :: sav_1
      real(8) , allocatable , dimension(:,:,:) :: sav2
      real(8) , allocatable , dimension(:,:,:) :: sav_2
      real(8) , allocatable , dimension(:,:,:) :: sav2a
      real(8) , allocatable , dimension(:,:,:) :: sav_2a
      real(8) , allocatable , dimension(:,:,:) :: sav4
      real(8) , allocatable , dimension(:,:,:) :: sav_4
      real(8) , allocatable , dimension(:,:,:) :: sav4a
      real(8) , allocatable , dimension(:,:,:) :: sav_4a
      real(8) , allocatable , dimension(:,:,:) :: sav6
      real(8) , allocatable , dimension(:,:,:) :: sav_6
#ifdef CLM
      real(8) , allocatable , dimension(:,:,:) :: sav_clmout
      real(8) , allocatable , dimension(:,:,:) :: sav_clmin
#endif

#endif

!---------- DATA init section--------------------------------------------

      contains 

        subroutine allocate_mod_tmpsav 
        implicit none
#ifdef MPP1
        allocate(sav0(iy,kz*4+2,jxp))
        allocate(sav_0(iy,kz*4+2,jx))
        allocate(sav0a(iy,kz+nnsg+5,jxp) )
        allocate(sav_0a(iy,kz+nnsg+5,jx) )
        allocate(sav0b(iy,kzp1,jxp))
        allocate(sav_0b(iy,kzp1,jx))
        allocate(sav0c(iy,kz*2,jxp))
        allocate(sav_0c(iy,kz*2,jx))
        allocate(sav0s(iy,kz,jxp))
        allocate(sav_0s(iy,kz,jx))
        allocate(sav0d(iy,nsplit*2,jxp))
        allocate(sav_0d(iy,nsplit*2,jx))
        allocate(sav1(iym1,kz*4+(kzp1*kzp2),jxp))
        allocate(sav_1(iym1,kz*4+(kzp1*kzp2),jx))
        allocate(sav2(iym1,nnsg*4+4,jxp))
        allocate(sav_2(iym1,nnsg*4+4,jx))
        allocate(sav2a(iym1,nnsg*5+1,jxp))
        allocate(sav_2a(iym1,nnsg*5+1,jx))
        allocate(sav4(iy,ntr*(kz*4+1),jxp))
        allocate(sav_4(iy,ntr*(kz*4+1),jx))
        allocate(sav4a(iym1,7,jxp))
        allocate(sav_4a(iym1,7,jx))
        allocate(sav6(kz,8,jxp))
        allocate(sav_6(kz,8,jx))
#ifdef CLM
        allocate(sav_clmout(iym1,9,jx))
        allocate(sav_clmin((iym1,9,jxp))
#endif 

#endif 
        end  subroutine allocate_mod_tmpsav

      end module mod_tmpsav
