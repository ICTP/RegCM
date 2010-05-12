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

      use mod_dynparam

      implicit none

      integer , allocatable , dimension(:) :: m
      real(8) , allocatable , dimension(:,:) :: a , a1 , a2 , a3 , a4 , &
               & d1 , d2 , e1 , e2 , e3 , g1 , g2 , g3 , s1 , s2 , w1 , &
               & w2 , x1
      integer , allocatable , dimension(:) :: iw2
      real(8) , allocatable , dimension(:) :: tbarf , thetaf
      real(8) , allocatable , dimension(:) :: thetah , tweigh
      real(8) , allocatable , dimension(:,:) :: w3
!
      real(8) :: alpha1 , alpha2 , pd , ps
      real(8) , allocatable , dimension(:) :: cpfac , dsigma , hbar ,   &
               & hweigh , tbarh
      real(8) , allocatable , dimension(:,:) :: hydroc , varpa1
      real(8) , allocatable , dimension(:,:) :: hydror , hydros , tau , &
               & zmatx , zmatxr
      real(8) , allocatable , dimension(:) :: sigmah
      real(8) , allocatable , dimension(:,:) :: varpa2
!
      real(8) , allocatable , dimension(:,:) :: am
      real(8) , allocatable , dimension(:) :: an

      real(8) , allocatable, dimension(:,:,:) :: dstor , hstor

      contains 

      subroutine allocate_mod_split
        implicit none
#ifdef MPP1
        allocate(dstor(iy,0:jxp+1,nsplit))
        allocate(hstor(iy,0:jxp+1,nsplit))
#else
        allocate(dstor(iy,jx,nsplit))
        allocate(hstor(iy,jx,nsplit))
#endif 
        allocate(m(nsplit))
        allocate(a(kz,kz))
        allocate(a1(kz,kz))
        allocate(a2(kz,kz))
        allocate(a3(kz,kz))
        allocate(a4(kz,kz))
        allocate(d1(kz,kz))
        allocate(d2(kz,kz))
        allocate(e1(kz,kz))
        allocate(e2(kz,kz))
        allocate(e3(kz,kz))
        allocate(g1(kz,kz))
        allocate(g2(kz,kz))
        allocate(g3(kz,kz))
        allocate(s1(kz,kz))
        allocate(s2(kz,kz))
        allocate(w1(kz,kz))
        allocate(w2(kz,kz))
        allocate(x1(kz,kz))
        allocate(iw2(kz))
        allocate(thetah(kz))
        allocate(tweigh(kz))
        allocate(tbarf(kzp1))
        allocate(thetaf(kzp1))
        allocate(w3(kzp1,kz))
        allocate(cpfac(kz))
        allocate(dsigma(kz))
        allocate(hbar(kz))
        allocate(hweigh(kz))
        allocate(tbarh(kz))
        allocate(hydroc(kz,kzp1))
        allocate(varpa1(kz,kzp1))
        allocate(hydror(kz,kz))
        allocate(hydros(kz,kz))
        allocate(tau(kz,kz))
        allocate(zmatx(kz,kz))
        allocate(zmatxr(kz,kz))
        allocate(sigmah(kzp1))
        allocate(varpa2(kzp1,kzp1))
        allocate(am(kz,nsplit))
        allocate(an(nsplit))
!
        end subroutine allocate_mod_split
!
      end module mod_split
