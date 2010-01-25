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

      subroutine wheneq(n, array, inc, itarg, indx, nval)

      implicit none

      integer n, inc, nval
      integer indx(*), array(*), itarg
C
      integer ina, i
C
      ina  = 1
      nval = 0
      if (inc .lt. 0) ina = (-inc)*(n-1)+1
      do i = 1, n
        if (array(ina) .eq. itarg) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do

      return
      end subroutine wheneq
