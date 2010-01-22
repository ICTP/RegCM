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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function lenchr(chrstg)
c----------------------------------------------------------------------
c
cl         return position of right-most non-blank, non-null character
cl         in chrstg.
c
c----------------------------------------------------------------------
      implicit none
      character*(*) chrstg
      lenchr = len(chrstg)
      if (chrstg(lenchr:lenchr).eq.' '.or.
     1    chrstg(lenchr:lenchr).eq.char(0)) then
10       continue
         lenchr = lenchr - 1
         if (chrstg(lenchr:lenchr).eq.' '.or.
     1       chrstg(lenchr:lenchr).eq.char(0)) then
            if (lenchr.gt.0) go to 10
         endif
      endif
      return
      end
