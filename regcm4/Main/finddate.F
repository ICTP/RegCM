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

      subroutine finddate(npos,idate)
      implicit none
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/message.cb"
      integer         mdatez
      common /datenum/mdatez(289276)
c
      integer npos,idate
      integer i
c
      i=0
 10   continue
      i=i+1
      if (mdatez(i).eq.idate) then
         npos=i
         go to 200
      end if
      if (i.gt.289276) goto 100
      go to 10
 100  continue
      
      write(aline, *) 'error in finddate : ', npos, idate
      call say
      call fatal(__FILE__,__LINE__,'Date not found')
 200  return
      end subroutine finddate
