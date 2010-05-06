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

      module mod_split

      use mod_regcm_param

      implicit none

      integer , dimension(nsplit) :: m
      real(8) , dimension(kz,kz) :: a , a1 , a2 , a3 , a4 , d1 , d2 ,   &
                                  & e1 , e2 , e3 , g1 , g2 , g3 , s1 ,  &
                                  & s2 , w1 , w2 , x1
      integer , dimension(kz) :: iw2
      real(8) , dimension(kzp1) :: tbarf , thetaf
      real(8) , dimension(kz) :: thetah , tweigh
      real(8) , dimension(kzp1,kz) :: w3
!
      real(8) :: alpha1 , alpha2 , pd , ps , pt , r
      real(8) , dimension(kz) :: cpfac , dsigma , hbar , hweigh , tbarh
      real(8) , dimension(kz,kzp1) :: hydroc , varpa1
      real(8) , dimension(kz,kz) :: hydror , hydros , tau , zmatx ,     &
                                  & zmatxr
      real(8) , dimension(kzp1) :: sigmah
      real(8) , dimension(kzp1,kzp1) :: varpa2
!
      real(8) , dimension(kz,nsplit) :: am
      real(8) , dimension(nsplit) :: an

      real(8) , allocatable, dimension(:,:,:) :: dstor , hstor

contains 
      subroutine allocate_mod_split
#ifdef MPP1
        allocate(dstor(iy,0:jxp+1,nsplit))
        allocate(hstor(iy,0:jxp+1,nsplit))
#else
       allocate(dstor(iy,jx,nsplit))
       allocate(hstor(iy,jx,nsplit))
#endif 
	end subroutine allocate_mod_split
      end module mod_split
