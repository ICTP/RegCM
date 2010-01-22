C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

*****************************FUNCTION EOMB******************************
*          COMPUTES AIR VAPOR PRESSURE AS A FUNCTION OF TEMP (in mb)   *
************************************************************************
      function eomb(x)
      implicit none
      real(kind=8)  eomb,x
      real(kind=8)  tr1
c
      tr1=1.0-(373.15/(x+273.15))
      eomb=1013.25*dexp(13.3185*tr1-1.976*tr1**2-0.6445*tr1**3
     &   -0.1299*tr1**4)
      return
      end
