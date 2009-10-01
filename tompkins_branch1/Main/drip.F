
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine drip
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  excess leaf water is put into rain or snow;
c  leaf water is reset to its maximum value.
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
c  ***********         xrun = leaf drip ; sdrop = snow drop off foliage
         etrrun(n,np) = 0.
         xrun(n,np) = ldew1d(n,np) - dewmx*vegt(n,np)
         sdrop(n,np) = 0.

c  ***********         test on maximum value of dew
         if(xrun(n,np).gt.0.) then
          etrrun(n,np) = xrun(n,np)
          ldew1d(n,np) = dewmx*vegt(n,np)
         endif

c  ***********         below freezing excess leaf water falls as snow
         if( (xrun(n,np).gt.0.) .and. (tm(n,np) .lt. c(67)) ) then
          etrrun(n,np) = 0.
          sdrop(n,np) = xrun(n,np)
         endif
        endif
       endif
      enddo
90    continue

      return
      end
