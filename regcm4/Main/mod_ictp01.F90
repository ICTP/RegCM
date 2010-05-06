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

      module mod_ictp01

      use mod_dynparam

      implicit none
      real(8) , allocatable , dimension(:,:) :: a , b
      real(8) , allocatable , dimension(:,:) :: ra , rs
      real(8) , allocatable , dimension(:,:) :: cdrd
      real(8) , allocatable , dimension(:,:) :: vpdc
      real(8) , allocatable , dimension(:,:) :: rppq , efe , qsatld
      real(8) , allocatable , dimension(:,:) :: dcd , etrc

      contains

      subroutine allocate_mod_ictp01
      implicit none
      allocate(a(nnsg,iym1))
      allocate(b(nnsg,iym1))
      allocate(ra(nnsg,iym1))
      allocate(rs(nnsg,iym1))
      allocate(cdrd(nnsg,iym1))
      allocate(vpdc(nnsg,iym1))
      allocate(rppq(nnsg,iym1))
      allocate(efe(nnsg,iym1))
      allocate(qsatld(nnsg,iym1))
      allocate(dcd(nnsg,iym1))
      allocate(etrc(nnsg,iym1))
      end subroutine allocate_mod_ictp01

      end module mod_ictp01
