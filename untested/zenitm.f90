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
!
      subroutine zenitm(coszrs,ivmx,jslc)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subroutine calculates the cosine of the solar zenith angle
!   for all longitude points of the mm42d domain. it needs as inputs
!   the longitude and latitude of the points, the initial date of the
!   simulation and the gmt. all these quantities are specified
!   in the initialization procedure of RegCM
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      use mod_regcm_param
      use mod_param1
      use mod_param3
      use mod_main
      use mod_date
      use mod_constants , only : degrad
      implicit none
!
! Arguments
!
      integer, intent (in) :: ivmx , jslc
      real(kind=8) , intent (inout), dimension(ix) :: coszrs
!
! Local variables
!
      integer :: ill
      real(kind=8) :: omega , tlocap , xt24 , xxlat
!
!***********************************************************************
!
      xt24 = dmod(lhour*60.+xtime,1440.D0)
      do ill = 1 , ivmx
        tlocap = xt24/60. + xlong(ill,jslc)/15.
        tlocap = dmod(tlocap+24.,24.D0)
        omega = 15.*(tlocap-12.)*degrad
        xxlat = xlat(ill,jslc)*degrad
!       coszrs = cosine of solar zenith angle
        coszrs(ill) = dsin(declin)*dsin(xxlat) + dcos(declin)           &
                    & *dcos(xxlat)*dcos(omega)
        coszrs(ill) = dmax1(0.D0,coszrs(ill))
      end do
!
      end subroutine zenitm
