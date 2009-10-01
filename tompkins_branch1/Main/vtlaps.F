c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vtlaps(t,sigma,r,pt,pd,nk)
      implicit none
      integer nk
      real*8  r,pt,pd
      real*8  rate,g,tstrat,zstrat,t0,p0
      parameter (rate=.0065, g=9.8, tstrat=218.15, zstrat=10769
     a,          t0=288.15, p0=101.325)
      real*8  t(nk),sigma(nk)
      real*8  fac,p,z
      integer k
c
c  this routine computes the temperature corresponding to a u. s.
c  standard atmosphere (see text by hess). units of p are cb.
c
      fac=r*rate/g
      do 1 k=1,nk
      p=sigma(k)*pd+pt
      t(k)=t0*((p/p0)**fac)
      z=(t0-t(k))/rate
      if (z.gt.zstrat) t(k)=tstrat
    1 continue
c
      return
      end
