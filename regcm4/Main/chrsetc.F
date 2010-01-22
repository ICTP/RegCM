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
      subroutine chrsetc(ichar,len,isetch)
      implicit none
      integer len
c
c****     this routine sets to the first len characters in ichar
c****     to the character variable specified by isetch
c
c****     on input-
c****         len    - number of characters to be set in ichar (integer)
c****         isetch - character to which each character in ichar will
c****                  be set equal to if input is character. (char)
c
c****     on output-
c****         ichar  - len charcters in ichar is set to the character
c****                  specified by either isetch.  (character)
c
      character*(*) ichar
      character*1 isetch
      integer i
c
      do 20 i=1,len
       ichar(i:i) = isetch
20    continue

      return
      end
