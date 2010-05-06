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

      module mod_grid

      implicit none

      real(4) , allocatable , dimension(:,:) :: coriol , dlat , dlon ,  &
           & msfx , snowcv , topogm , toposdgm , xlandu , xlat , xlon
      real(4) , allocatable , dimension(:,:) :: pa , sst1 , sst2 ,      &
           & tlayer , za , ice1 , ice2
      real(4) , allocatable , dimension(:,:) :: b3pd
      real(4) , allocatable , dimension(:) :: dsigma , sigma2
      real(4) , allocatable , dimension(:) :: sigmaf
      real(4) :: delx , grdfac
      integer :: i0 , i1 , j0 , j1
      real(4) :: lat0 , lat1 , lon0 , lon1

      contains

      subroutine init_grid(iy,jx,kz)
      implicit none
      integer , intent(in) :: iy , jx , kz
      allocate(coriol(jx,iy))
      allocate(dlat(jx,iy))
      allocate(dlon(jx,iy))
      allocate(msfx(jx,iy))
      allocate(snowcv(jx,iy))
      allocate(topogm(jx,iy))
      allocate(toposdgm(jx,iy))
      allocate(xlandu(jx,iy))
      allocate(xlat(jx,iy))
      allocate(xlon(jx,iy))
      allocate(pa(jx,iy))
      allocate(sst1(jx,iy))
      allocate(sst2(jx,iy))
      allocate(tlayer(jx,iy))
      allocate(za(jx,iy))
      allocate(ice1(jx,iy))
      allocate(ice2(jx,iy))
      allocate(b3pd(jx,iy))
      allocate(dsigma(kz))
      allocate(sigma2(kz))
      allocate(sigmaf(kz+1))
      end subroutine init_grid

      subroutine free_grid
      implicit none
      deallocate(coriol)
      deallocate(dlat)
      deallocate(dlon)
      deallocate(msfx)
      deallocate(snowcv)
      deallocate(topogm)
      deallocate(toposdgm)
      deallocate(xlandu)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(pa)
      deallocate(sst1)
      deallocate(sst2)
      deallocate(tlayer)
      deallocate(za)
      deallocate(ice1)
      deallocate(ice2)
      deallocate(b3pd)
      deallocate(dsigma)
      deallocate(sigma2)
      deallocate(sigmaf)
      end subroutine free_grid

      end module mod_grid
