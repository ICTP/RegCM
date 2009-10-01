c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine frawat
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  **** determines fraction of foliage covered by water fwet, and
c  **** the fraction of foliage that is dry transpiring leaf fdry.
c  note: their defns differ - fwet is the fraction of all veg surfaces
c  which are wet because stems can evaporate, fdry is the fraction
c  of lai which is dry because only leaves can transpire
c
c  ldew1d(np) is in kg/m**2/s
c  fwet   = ratio of dew to max value to 2/3 power
c           ( 2/3 power comes from deardorff (1978) )
c              ** keep fwet le 1.0 **
c  dewmxi = inverse of max allowed dew depth on leaf in mm
c
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
         fwet(n,np) = 0.
         if(ldew1d(n,np).gt.0.) then
          fwet(n,np) = ((dewmxi/vegt(n,np))*ldew1d(n,np))**.666666666666
          fwet(n,np) = dmin1(fwet(n,np),1.d0)
         endif
         fdry(n,np) = (1.-fwet(n,np))*xlai(n,np) / xlsai(n,np)
        endif
       endif
      enddo
90    continue

      return
      end
