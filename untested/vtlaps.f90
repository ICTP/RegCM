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
 
      subroutine vtlaps(t,sigma,r,pt,pd,nk)
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: rate = .0065 , g = 9.8 , tstrat = 218.15 , &
                           & zstrat = 10769 , t0 = 288.15 , p0 = 101.325
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
      fac = r*rate/g
      do k = 1 , nk
        p = sigma(k)*pd + pt
        t(k) = t0*((p/p0)**fac)
        z = (t0-t(k))/rate
        if ( z.gt.zstrat ) t(k) = tstrat
      end do
!
      end subroutine vtlaps
