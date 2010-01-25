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

      subroutine lshfch(ichar,lchar,ioff,iflg)
c
c****     this routine shifts the characters in ichar ioff characters
c****     to the left.  the right side off ichar is either blank
c****     or null filled depending on iflg.  for example, if ichar
c****     is to start with the 3rd character instead of the 1st ioff
c****     would be set to 2 shifting ichar to the left.
c
c****     on input-
c****         ichar  - character input data. (character)
c****         lchar  - number of characters in ichar. (integer)
c****         ioff   - number of characters ioff is to be offset to
c****                  left. (integer)
c****         iflg   - flag specifying type of fill to be used for
c****                  the ioff rightmost characters in ichar.  if
c****                      0 - blank fill is used
c****                      1 - null fill is used
c
c****     on output-
c****         ichar  - ichar is retured shifted ioff characters to the
c****                  left with the ioff rightmost characters filled
c****                  in accordance with iflg.
c
      implicit none
#ifdef MPP1
#     include "mpiregcm.inc"
      integer ierr
#endif
      integer lchar,ioff,iflg
      character*(*) ichar
      character*1 temp
      integer iblnk,null,ipc,i,ifill
      data iblnk /32/
      data null  /0/
c
      ipc = 1
      do 10 i=ioff+1,lchar
       temp = ichar(i:i)
       ichar(ipc:ipc) = temp
       ipc=ipc+1
10    continue
c
      if (iflg.eq.0) then
       ifill = iblnk
      else if (iflg.eq.1) then
       ifill = null
      else
       write(*,*) ' iflg in lshfch out of range'
#ifdef MPP1
       write(*,*) ' stop in subroutine lshfch'
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
       stop
#endif
      end if
c
      do 20 i=ipc,lchar
       ichar(i:i) = char(ifill)
20    continue

      return
      end
