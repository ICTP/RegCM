c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vnorml(z,s,nk,nk1)
      implicit none
      integer nk,nk1
      real*8  z(nk,nk),s(nk1)
      integer k,l,kmax
      real*8  zmax,v,a
c
c  this routine normalizes the columns of z such that the component
c  with the largest absolute value is positive, and the sum of the
c  mass-weighted squares equals one.
c
      do 3 l=1,nk
      zmax=-1.
      v=0.
c
      do 1 k=1,nk
      a=dabs(z(k,l))
      if (a.gt.zmax) then
         zmax=a
         kmax=k
      endif
      v=(s(k+1)-s(k))*a*a + v
    1 continue
c
      a=(z(kmax,l)/zmax)/dsqrt(v)
      do 2 k=1,nk
    2 z(k,l)=a*z(k,l)
    3 continue
c
      return
      end
