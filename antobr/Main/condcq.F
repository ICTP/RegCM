
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine condcq
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     dimensional and non-dimensional latent heat conductances
c               for canopy and soil flux calculations
c
c     latent fluxes differ from sensible due to stomatal resistance
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
         rgr(n,np) = gwet1d(n,np)
         wtlq(n,np) = wtlh(n,np) * rpp(n,np)
         wtgq(n,np) = wtg(n,np) * rgr(n,np)
         wtsqi(n,np) = 1. / (wta(n,np) + wtlq(n,np) + wtgq(n,np))
         wtgq0(n,np) = wtgq(n,np) * wtsqi(n,np)
         wtlq0(n,np) = wtlq(n,np) * wtsqi(n,np)
         wtglq(n,np) = wtgq0(n,np) + wtlq0(n,np)
         wtaq0(n,np) = 1. - wtglq(n,np)
         wtgaq(n,np) = wtaq0(n,np) + wtgq0(n,np)
        endif
       endif
      enddo
90    continue

      return
      end
