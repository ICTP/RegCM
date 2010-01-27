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
      subroutine printsp( array, ix, jx, kx, iinc, jinc, kref, aname)
c
      implicit none
      integer ix,jx,kx,iinc,jinc,kref
      real(kind=8)  array(ix,jx,kx)
      character*18 cform
      character*08 aname
      integer i,j,ineg,inum,jnum,jcheck
c
      ineg = - iinc
      inum = (ix-2) / iinc + 1
      jnum = (jx-2) / jinc + 1
c
      jcheck = 11 * jnum
      if(jcheck.gt.127) print *,'too many values printed.'
c
      write(cform,5) jnum
5     format('(1x,i2,2x,',i2,'g11.3)')
c
      print *,'  field = ', aname
      write(6,cform)(i,(array(i,j,kref),j=1,jx-1,jinc),i=ix-1,1,ineg)
c
      return
      end subroutine printsp
