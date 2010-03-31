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
 
      subroutine radoz2(o3vmr,pint,plol,plos)

!-----------------------------------------------------------------------
!
! Computes the path length integrals to the model interfaces given the
! ozone volume mixing ratio
!
!---------------------------Code history--------------------------------
!
! Original version:     CCM1
! Standardized:         J. Rosinski, June 1992
! Reviewed:             J. Kiehl, B. Briegleb, August 1992
! Mixing ratio version: Bruce Biegleb, September 1992
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_comozp
      implicit none
!
!------------------------------Input arguments--------------------------
!
! o3vmr   - ozone volume mixing ratio
! pint    - Model interface pressures
!
!----------------------------Output arguments---------------------------
!
! plol    - Ozone prs weighted path length (cm)
! plos    - Ozone path length (cm)
!
!
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: o3vmr
      real(8) , dimension(iym1,kzp1) :: pint , plol , plos
      intent (in) o3vmr , pint
      intent (inout) plol , plos
!
!---------------------------Local workspace-----------------------------
!
! i       - longitude index
! k       - level index
!
!-----------------------------------------------------------------------
!
! Local variables
!
      integer :: i , k
!
!     Evaluate the ozone path length integrals to interfaces;
!     factors of .1 and .01 to convert pressures from cgs to mks:
!
!     Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
!
      do i = 1 , iym1
        plos(i,1) = 0.1*cplos*o3vmr(i,1)*pint(i,1)
        plol(i,1) = 0.01*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
      end do
      do k = 2 , kzp1
        do i = 1 , iym1
          plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)              &
                    & *(pint(i,k)-pint(i,k-1))
          plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)             &
                    & *(pint(i,k)*pint(i,k)-pint(i,k-1)*pint(i,k-1))
        end do
      end do
!
      end subroutine radoz2
