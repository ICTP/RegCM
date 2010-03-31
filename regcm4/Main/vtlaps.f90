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
 
      subroutine vtlaps(t,sigma,r,pt,pd,nk)

      use mod_constants , only : rgti , stdp , stdt , lrate
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: tstrat = 218.15 , zstrat = 10769
      real(8) :: p0
!
! Dummy arguments
!
      integer :: nk
      real(8) :: pd , pt , r
      real(8) , dimension(nk) :: sigma , t
      intent (in) nk , pd , pt , r , sigma
      intent (inout) t
!
! Local variables
!
      real(8) :: fac , p , z
      integer :: k
!
!  this routine computes the temperature corresponding to a u. s.
!  standard atmosphere (see text by hess). units of p are cb.
!
      p0 = stdp/1000.D0
      fac = r*lrate*rgti
      do k = 1 , nk
        p = sigma(k)*pd + pt
        t(k) = stdt*((p/p0)**fac)
        z = (stdt-t(k))/lrate
        if ( z.gt.zstrat ) t(k) = tstrat
      end do
!
      end subroutine vtlaps
