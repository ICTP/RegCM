c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine condch
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     dimensional and non-dimensional sensible heat conductances
c               for canopy and soil flux calculations
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'
      real*8  a,b
      common/ab/ a(NNSG,nbmax), b(NNSG,nbmax)
      integer n,np

c     csoilc = constant drag coefficient for soil under canopy;
c                    prescribed in subroutine bconst.
c
c     symbols used for weights are:   wt : weight
c                                      a : air
c                                      l : leaf
c                                      i : inverse
c                                      s : sum
c                                      h : sensible heat
c                                      q : water vapor
c                                      0 : normalized (sums to one)
c                                      g : ground
c
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         uaf(n,np) = vspda(n,np) * dsqrt(cdr(n,np))
         cf(n,np)  = 0.01 * sqrtdi(lveg(n,np)) / dsqrt(uaf(n,np))
         wta(n,np) = sigf(n,np) * cdr(n,np) * vspda(n,np)
         wtlh(n,np) = cf(n,np) * uaf(n,np) * vegt(n,np)
         wtg(n,np) = csoilc * uaf(n,np) * sigf(n,np)
         wtshi(n,np) = 1. / (wta(n,np)+wtlh(n,np)+wtg(n,np))
         wtl0(n,np) = wtlh(n,np) * wtshi(n,np)
         wtg0(n,np) =  wtg(n,np) * wtshi(n,np)
         wtgl(n,np) = wtl0(n,np) +  wtg0(n,np)
         wta0(n,np) = 1. - wtgl(n,np)
         wtga(n,np) = wta0(n,np) + wtg0(n,np)
        endif
       endif
      enddo
90    continue
c
      return
      end
