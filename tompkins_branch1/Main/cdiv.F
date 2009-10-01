
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      implicit none
      real*8  ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      real*8  s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
