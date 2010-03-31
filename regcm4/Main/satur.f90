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
 
      subroutine satur(qsat,t,p)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****  calculates saturation vapor pressure (eg)
!           and qsat = saturated specific humidity (dimensionless)
!
!           uses tetens formula (1930) (ref. riegel,1974,jam,p606
!                                                 equation 1)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      use mod_regcm_param
      use mod_ictp01
      use mod_bats
      use mod_constants , only : tzero , c1es , c4ies , c4les , c3ies , &
                   &             c3les , ep2
      implicit none
!
! Dummy arguments
!
      real(8) , dimension(nnsg,iym1) :: p , qsat , t
      intent (in) p , t
      intent (out) qsat
!
! Local variables
!
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( t(n,i).le.tzero ) then
            a(n,i) = c3ies
            b(n,i) = c4ies
          else
            a(n,i) = c3les
            b(n,i) = c4les
          end if
          eg(n,i) = c1es*dexp(a(n,i)*(t(n,i)-tzero)/(t(n,i)-b(n,i)))
          qsat(n,i) = ep2*eg(n,i)/(p(n,i)-0.378D0*eg(n,i))
        end do
      end do
 
      end subroutine satur
