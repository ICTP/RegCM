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

      function oned(x,a,b,c,d)
      implicit none
!
! Dummy arguments
!
      real(8) :: a , b , c , d , x
      real(8) :: oned
      intent (in) a , b , c , d , x
!
      oned = 0.
      if ( x==0. ) oned = b
      if ( x==1. ) oned = c
      if ( b*c==0. ) return
      if ( a*d==0. ) then
        oned = b*(1.0-x) + c*x
        if ( a/=0.0 ) oned = b + x*(0.5*(c-a)+x*(0.5*(c+a)-b))
        if ( d/=0.0 ) oned = c + (1.0-x)                                &
                           & *(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
        go to 99999
      end if
      oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))                  &
           & + x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      return
99999 continue
      end function oned
