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

      subroutine vcheke (er,ei,nk,numerr,aname)
      implicit none
      integer nk,numerr
      real(kind=8)  tol
      parameter (tol=1.e-9)
      character*8 aname
      real(kind=8)  er(nk), ei(nk)
      integer numneg,n,nimag
      real(kind=8)  emax
c
c  check that eigenvalues are real and positive valued.
c
      numneg=0
      emax=0.
      do 1 n=1,nk
      if (er(n).le.0.) numneg=numneg+1
      if (er(n).gt.emax) emax=er(n)
    1 continue
c
      nimag=0
      do 2 n=1,nk
      if ( ei(n)/emax .gt. tol ) nimag=nimag+1
    2 continue
c
      if ( numneg+nimag .eq. 0 ) return
c
      numerr=numerr+1
      print 3,aname,numneg,nimag
    3 format('0 problem with equivalent depths determined from ',a8,/,
     a10x,i3,' depths are nonpositive valued',10x,i3,' are not real')
c
      return
      end subroutine vcheke
