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

      subroutine ljchrn(chrv,nlshf)
c
c****     routine left justifies a character variable cnull filling on
c****     the right.  chrv is the character variable.  nlshf tells
c****     how many characters the variable was shifted left
c
      implicit none
      integer nlshf,lnchrv
      character *(*) chrv
      character *1 blank,cnull
      data blank /' '/
c
      cnull = char(0)
      lnchrv = len(chrv)
      nlshf = 0
10    continue
      if (chrv(1:1).eq.blank.or.chrv(1:1).eq.cnull) then
         call lshfch(chrv,lnchrv,1,1)
         nlshf = nlshf+1
         if (nlshf.ne.lnchrv) go to 10
      end if
      return
      end subroutine ljchrn
