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

      use mod_regcm_param , only : jx , iy , kz

      implicit none

      real(4) , dimension(jx,iy) :: coriol , dlat , dlon , msfx ,       &
                 & snowcv , topogm , toposdgm , xlandu , xlat , xlon
      real(4) , dimension(jx,iy) :: pa , sst1 , sst2 , tlayer , za ,    &
                 & ice1 , ice2
      real(4) , dimension(jx,iy) :: b3pd
      real(4) , dimension(kz) :: dsigma , sigma2
      real(4) , dimension(kz+1) :: sigmaf
      real(4) :: delx , grdfac
      integer :: i0 , i1 , j0 , j1
      real(4) :: lat0 , lat1 , lon0 , lon1

      end module mod_grid
