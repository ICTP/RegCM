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

      subroutine elmhes(nm,n,low,igh,a,int)
      implicit none
c
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real(kind=8)  a(nm,n)
      real(kind=8)  x,y
      integer int(igh)
c
c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the input matrix.
c
c     on output
c
c        a contains the hessenberg matrix.  the multipliers
c          which were used in the reduction are stored in the
c          remaining triangle under the hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0e0
         i = m
c
         do 100 j = m, igh
            if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
c
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c     .......... end interchange ..........
  130    if (x .eq. 0.0e0) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0e0) go to 160
            y = y / x
            a(i,mm1) = y
c
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
c
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
c
  160    continue
c
  180 continue
c
  200 return
      end
