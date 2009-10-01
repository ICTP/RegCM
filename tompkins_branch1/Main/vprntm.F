c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vprntm(a, n1, n2, nam)
      implicit none
      integer n1,n2
      character nam*8
      real*8  a(n1,n2)
      integer k,l
c
      print 10, nam
      do 1 k=1,n1
    1 print 11,k,(a(k,l),l=1,n2)
   10 format(1h1,a8,/)
   11 format(1x,i3,5x,1p11g11.3,1x,/,9x,1p11g11.3)
      return
      end
