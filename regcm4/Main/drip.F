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

      subroutine drip
c
c  excess leaf water is put into rain or snow;
c  leaf water is reset to its maximum value.
c
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"

      integer n,np

      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.0.5) then
            if (sigf(n,np).gt.0.001) then
c  ***********         xrun = leaf drip ; sdrop = snow drop off foliage
              etrrun(n,np) = 0.
              xrun(n,np) = ldew1d(n,np) - dewmx*vegt(n,np)
              sdrop(n,np) = 0.

c  ***********         test on maximum value of dew
              if (xrun(n,np).gt.0.) then
                etrrun(n,np) = xrun(n,np)
                ldew1d(n,np) = dewmx*vegt(n,np)
              end if

c  ***********         below freezing excess leaf water falls as snow
              if ( (xrun(n,np).gt.0.) .and. (tm(n,np) .lt. c(67)) ) then
                etrrun(n,np) = 0.
                sdrop(n,np) = xrun(n,np)
              end if
            end if
          end if
        end do
      end do

      return
      end subroutine drip
