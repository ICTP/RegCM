c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vprntv(a, n, nam)
c
      implicit none
      integer n
      character nam*8
      real*8  a(n)
c
      print 1, nam, a
    1 format('0',a8,1x,1p11g11.3,1x,/,9x,1p11g11.3)
      return
      end
