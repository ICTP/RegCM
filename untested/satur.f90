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
      implicit none
!
! Dummy arguments
!
      real(8) , dimension(nnsg,nbmax) :: p , qsat , t
      intent (in) p , t
      intent (out) qsat
!
! Local variables
!
      integer :: n , np
! 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( t(n,np).le.c(67) ) then
            a(n,np) = c(70)
            b(n,np) = c(71)
          else
            a(n,np) = c(72)
            b(n,np) = c(73)
          end if
          eg(n,np) = c(74)                                              &
                   & *dexp(a(n,np)*(t(n,np)-c(67))/(t(n,np)-b(n,np)))
          qsat(n,np) = c(75)*eg(n,np)/(p(n,np)-c(76)*eg(n,np))
        end do
      end do
 
      end subroutine satur
