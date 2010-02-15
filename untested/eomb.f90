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
 
      function eomb(x)

!****************************FUNCTION EOMB******************************
!          COMPUTES AIR VAPOR PRESSURE AS A FUNCTION OF TEMP (in K)    *
!***********************************************************************

      use mod_constants , only : stdpmb , tboil , tmelt
      implicit none
!
! Dummy arguments
!
      real(8) :: x
      intent (in) x
      real(8) :: eomb
!
! Local variables
!
      real(8) :: tr1
!
      tr1 = 1.0 - (tboil/(x+tmelt))
      eomb = stdpmb*dexp(13.3185*tr1-1.976*tr1**2-0.6445*tr1**3-       &
           & 0.1299*tr1**4)

      end function eomb
