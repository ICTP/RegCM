
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine vcover
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     provides leaf and stem area parameters;
c     depends on climate through subsoil temperatures.

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'
      integer n,np

      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         seasb(n,np)=dmax1(0.d0
     &                    ,1.-0.0016*dmax1(298.-tgb1d(n,np),0.d0)**2)
        endif
       endif
      enddo
90    continue

      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         xlai(n,np)=xla(lveg(n,np))
         xlai(n,np)=xlai(n,np)
     &             +(xlai0(lveg(n,np))-xlai(n,np))*(1.-seasb(n,np))
         rlai(n,np)  = xlai(n,np) + sai(lveg(n,np))
         xlsai(n,np) = xlai(n,np) + sai(lveg(n,np))
         vegt(n,np)  = sigf(n,np) * xlsai(n,np)
        endif
       endif
      enddo
91    continue
c
      return
      end
