c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vmultm(a1,a2,a3,nk)
      implicit none
      integer nk
      real*8  a1(nk,nk),a2(nk,nk),a3(nk,nk)
      integer k,l,mm
      real*8  c
c
      do 10 l=1,nk
      do 10 k=1,nk
      a1(k,l)=0.
      do 10 mm=1,nk
   10 a1(k,l)=a1(k,l)+a2(k,mm)*a3(mm,l)
      return
c
      entry vaddms(a1,a2,a3,nk)
      do 20 l=1,nk
      do 20 k=1,nk
   20 a1(k,l)=a2(k,l)+a3(k,l)
      return
c
      entry vsubtm(a1,a2,a3,nk)
      do 30 l=1,nk
      do 30 k=1,nk
   30 a1(k,l)=a2(k,l)-a3(k,l)
      return
c
      entry vmultc(a1,a2,nk,c)
      do 40 l=1,nk
      do 40 k=1,nk
   40 a1(k,l)=a2(k,l)*c
      return
      end
