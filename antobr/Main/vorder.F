c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vorder(z,hbar,wz,wh,nk)
      implicit none
      integer nk
      real*8  z(nk,nk),hbar(nk),wz(nk,nk),wh(nk,2)
      integer k,l,kmax
      real*8  hmax
c
c  this routine orders the components of hbar so they are largest to
c  smallest valued.  the columns of z are reorded so that they
c  correspond to the same (but reordered) components of hbar.
c
      do 1 k=1,nk
      wh(k,1)=hbar(k)
      wh(k,2)=0.
      do 1 l=1,nk
    1 wz(k,l)=z(k,l)
c
      do 4 l=1,nk
      hmax=-1.d100
      do 2 k=1,nk
      if ((wh(k,2).eq.0.).and.(wh(k,1).gt.hmax)) then
         hmax=wh(k,1)
         kmax=k
      endif
    2 continue
c
      hbar(l)=hmax
      wh(kmax,2)=1.
      do 3 k=1,nk
    3 z(k,l)=wz(k,kmax)
    4 continue
c
      return
      end
