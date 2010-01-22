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

      subroutine vchekt (tbarh,sigmah,sigmaf,xkappa,pt,pd,nk,numerr)
      implicit none
      integer nk,numerr
      real(kind=8)  xkappa,pt,pd
      real(kind=8)  tbarh(nk), sigmah(nk), sigmaf(nk+1)
      logical lstab
      integer k
      real(kind=8)  ds1,ds2,tb,g1,g2
c
c  check if tbar is stable.  this is not the actual stability condition
c  consistent with the model finite differences.
c
      lstab=.true.
      do 1 k=1,nk-1
      ds1=sigmaf(k+1)-sigmaf(k)
      ds2=sigmaf(k+2)-sigmaf(k+1)
      tb=(ds1*tbarh(k)+ds2*tbarh(k+1))/(ds1+ds2)
      g1=xkappa*tb/(sigmaf(k+1)+pt/pd)
      g2=(tbarh(k+1)-tbarh(k))/(sigmah(k+1)-sigmah(k))
      if (g1-g2 .lt. 0.) lstab=.false.
    1 continue
      if (.not.lstab) then
         numerr=numerr+1
         print 2
      endif
    2 format('0 indication that tbarh statically unstable')
c
      return
      end
