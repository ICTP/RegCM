      subroutine initlk( veg2d, ix1, jx1 )

      implicit none
      integer ix1,jx1
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'

c ******  lkpts = number of lake points desired

      real*8  t(400), veg2d(ix1,jx1)
      integer ilake, jlake, depth, freeze
      integer lcount,iin,iout,numpts
      common /ldata/ lcount, iin, iout, numpts

      real*8  hi,hice,hsnow,eta
      data hi, hice, hsnow, eta / .01, 0., 0., .5 /
      integer idata,lakeset,n,j


c ******  unit number containing lake pt locations, depths
      idata = 40

      lcount = 0
      iin = 41
      iout = 42
      numpts = lkpts


c read in vegetation type desired for lake points
c  (18 = mixed woodland (no lake model);  14 = inland water (lake model)

      read(idata,89)  lakeset
89    format(i2)
      print *, '*** lake points set to bats surface type ', lakeset



c initialize data for lake model
      do 20 n=1,lkpts
       read(idata,90)  ilake, jlake, depth
90     format(3i4)
       freeze = 1


c ******     Initially set eta to mean values
       if(depth.lt.50) then
        eta = .7
       elseif (depth .gt. 100) then
        eta = .3
       else
          eta = .5
       endif

c ******     Initially set lake points isothermal at 6.0 C (June)
       do 10 j=1,depth
        t(j) = 6.0
10     continue

       veg2d(ilake,jlake) = dble(lakeset)
       write(iin) ilake, jlake, depth, freeze, hi, hice, hsnow
     &            , eta, (t(j),j=1,depth)
       print *, ilake, jlake, depth, freeze

20    continue

      rewind(iin)

      return
      end
