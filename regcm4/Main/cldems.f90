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
 
      subroutine cldems(clwp,fice,rei,emis)

!-----------------------------------------------------------------------
!
! Compute cloud emissivity using cloud liquid water path (g/m**2)
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Kiehl
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, J. Kiehl, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      implicit none
!
! PARAMETER definitions
!
!     longwave absorption coeff (m**2/g)
      real(8) , parameter :: kabsl = 0.090361
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! clwp    - cloud liquid water path (g/m**2)
! rei     - ice effective drop size (microns)
! fice    - fractional ice content within cloud
!
!     Output arguments
!
! emis    - cloud emissivity (fraction)
!
!-----------------------------------------------------------------------
!
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: clwp , emis , fice , rei
      intent (in) clwp , fice , rei
      intent (out) emis
!
! Local variables
!
!-----------------------------------------------------------------------
!
! i, k    - longitude, level indices
! kabs    - longwave absorption coeff (m**2/g)
! kabsi   - ice absorption coefficient
!
!-----------------------------------------------------------------------
!
      integer :: i , k
      real(8) :: kabs , kabsi
!
      do k = 1 , kz
        do i = 1 , iym1
          kabsi = 0.005 + 1./rei(i,k)
          kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1. - dexp(-1.66*kabs*clwp(i,k))
        end do
      end do
!
      end subroutine cldems
