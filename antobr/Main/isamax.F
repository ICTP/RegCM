      integer function isamax(n,sx,incx)
      implicit none
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      integer incx,n
      real*8  sx(1+(n-1)*incx),smax
      integer i,ix
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dabs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = dabs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dabs(sx(1))
      do 30 i = 2,n
         if(dabs(sx(i)).le.smax) go to 30
         isamax = i
         smax = dabs(sx(i))
   30 continue
      return
      end
