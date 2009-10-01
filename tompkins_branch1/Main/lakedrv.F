      SUBROUTINE LAKEDRV(jslc)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/iunits.cb'
      include 'Commons/bats.cb'
#ifdef MPP1
      include 'Common2/bats.cb2'
#else
      include 'Commons/bats.cb2'
#endif
      include 'Commons/date.cb'

      real*8  dayl,tl,vl,zl,ql,fsw,flw,prec,hsen
      real*8  hlat,tgl,evl,aveice,snow
      integer ilake,jslc,n

      do 7 ilake=1,ix-1
      do n=1,NNSG
        if( veg2d1(n,ilake,jslc).eq.14. ) then
          dayl = (NNNNNN-NSTRT0)/4.+(xtime + dtmin) / 1440.
          tl  =       ts1d(n,ilake)
          vl  = dsqrt( us1d(ilake)**2 + vs1d(ilake)**2 )
          zl  =       z1d(n,ilake)
          ql  =       qs1d(n,ilake)
          fsw =       fsw1d(ilake)
          flw = -1. * flw1d(ilake)
          prec=       prca2d(ilake,jslc)   !  units of prec = mm
          hsen= -1. * sent1d(n,ilake)
          hlat= -1. * evpr1d(n,ilake)
          call lake( iutlak, dayl, dtlake, tl, vl, zl, ql, fsw, flw
     a             , hsen, hlat, tgl, evl, prec, aveice, snow)

          tg1d(n,ilake) = tgl
          tgb1d(n,ilake) = tgl
          if(aveice.le.10.)then
            ldoc1d(n,ilake) = 0.
            sice1d(n,ilake) = 0.
            scv1d(n,ilake) = 0.
            sag1d(n,ilake) = 0.
          else
            ldoc1d(n,ilake) = 2.
            sice1d(n,ilake) = aveice    !  units of ice = mm
            scv1d(n,ilake)  = snow    !  units of snow = mm h2o
          endif

        endif
      enddo
7     continue

      return
      end

