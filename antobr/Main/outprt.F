#ifdef MPP1
      subroutine outprt(iexec)
      implicit none
      integer iexec
      return
      end
#else
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outprt(iexec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine creates the printer output for a quick check.   c
c     the complete analysis is handled by other program "dataflow"    c
c     later.                                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer iexec
      include 'regcm.param'
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
      include 'Commons/iunits.cb'
      include 'Commons/main.cb'
      include 'Commons/cvaria.cb'
      include 'Commons/bdycod.cb'
      include 'Commons/pmoist.cb'
      include 'Commons/date.cb'
c
      real*8  xlatr(jx,ix),xlongr(jx,ix)
      real*8  hscr1(ix,jx),hscr2(ix,jx)                      
      real*8  hscr1r(jx,ix),hscr2r(jx,ix),fr(jx,ix),msfdr(jx,ix)
      real*8  psar(jx,ix),tgar(jx,ix),raincr(jx,ix),rainncr(jx,ix)
      real*8  scrx1(kx,ix),scrx2(kx,ix),scrx1r(ix,kx),scrx2r(ix,kx)
      character c40nam*40
      integer ixn,jxn,kxn
      data    ixn,jxn,kxn/2,2,1/
      integer i,j,k,jcx,kout
      real*8  xth,psabar,psab1,psab2,const,pres,satvp,qvs
c
c----------------------------------------------------------------------
c
      xth = xtime/60.
      jcx = jxsex
      kout = kxout
      print 10,ldatez+xtime/1440.,ktau,jyear
10    format(///1x,'--------------------------------------------------',
     1/1x,'*****',4x,'large domain at t = ',f17.5,
     2' minutes, ktau = ',i7,' in year=',i4,3x,'*****'/1x,
     &'------------------------------------------------------------'///)
c
c-----north-south cross sections:
c
c.....u and v:
c
      if (iexec .eq. 1) then
c
c.....first time :
c
       if (jcx .eq. 1) then
c......for jcx = 1:
        do 50 i=2,ilx
        do 50 k=1,kx
         psabar=0.5*(psa(i,jcx)+psa(i-1,jcx))
         scrx1(k,i)=ua(i,k,jcx)/psabar
         scrx2(k,i)=va(i,k,jcx)/psabar
50      continue
        do 60 k=1,kx
         scrx1(k,1)=ua(1,k,jcx)/psa(1,jcx)
         scrx1(k,ix)=ua(ix,k,jcx)/psa(ilx,jcx)
         scrx2(k,1)=va(1,k,jcx)/psa(1,jcx)
         scrx2(k,ix)=va(ix,k,jcx)/psa(ilx,jcx)
60      continue
       elseif (jcx .eq. jx) then
c......for jcx = jx:
        do 80 i=2,ilx
        do 80 k=1,kx
         psabar=0.5*(psa(i,jlx)+psa(i-1,jlx))
         scrx1(k,i)=ua(i,k,jcx)/psabar
         scrx2(k,i)=va(i,k,jcx)/psabar
80      continue
        do 90 k=1,kx
         scrx1(k,1)=ua(1,k,jcx)/psa(1,jlx)
         scrx1(k,ix)=ua(ix,k,jcx)/psa(ilx,jlx)
         scrx2(k,1)=va(1,k,jcx)/psa(1,jlx)
         scrx2(k,ix)=va(ix,k,jcx)/psa(ilx,jlx)
90      continue
       else
c......interior slice:
        do 20 i=2,ilx
        do 20 k=1,kx
         psabar=0.25*(psa(i,jcx)+psa(i,jcx-1)+psa(i-1,jcx)+
     1          psa(i-1,jcx-1))
         scrx1(k,i)=ua(i,k,jcx)/psabar
         scrx2(k,i)=va(i,k,jcx)/psabar
20      continue
        do 30 k=1,kx
         psab1=0.5*(psa(1,jcx-1)+psa(1,jcx))
         psab2=0.5*(psa(ilx,jcx-1)+psa(ilx,jcx))
         scrx1(k,1)=ua(1,k,jcx)/psab1
         scrx1(k,ix)=ua(ix,k,jcx)/psab2
         scrx2(k,1)=va(1,k,jcx)/psab1
         scrx2(k,ix)=va(ix,k,jcx)/psab2
30      continue
       endif
c
      elseif(iexec.gt.1)then
c
c.....subsequent calls:
c
       if (jcx .eq. 1) then
c......for jcx = 1:
        do 140 i=1,ix
        do 140 k=1,kx
         scrx1(k,i)=uj1(i,k)
         scrx2(k,i)=vj1(i,k)
140     continue
       elseif (jcx .eq. jx) then
c......for jcx = jx:
        do 160 i=1,ix
        do 160 k=1,kx
         scrx1(k,i)=ujl(i,k)
         scrx2(k,i)=vjl(i,k)
160     continue
       else
c......interior slice:
        do 110 i=2,ilx
        do 110 k=1,kx
         psabar=0.25*(psa(i,jcx)+psa(i,jcx-1)+psa(i-1,jcx)+
     1          psa(i-1,jcx-1))
         scrx1(k,i)=ua(i,k,jcx)/psabar
         scrx2(k,i)=va(i,k,jcx)/psabar
110     continue
        do 120 k=1,kx
         scrx1(k,1)=ui1(k,jcx)
         scrx1(k,ix)=uil(k,jcx)
         scrx2(k,1)=vi1(k,jcx)
         scrx2(k,ix)=vil(k,jcx)
120     continue
       endif
c
      endif      !end iexec test
c
      write(c40nam,180)jcx
180   format('  u  cross-section at j = ',i3,8x)
      call mapsmp(scrx1,scrx1r,kx,ix,1,kx,kxn,1,ix,ixn,0.d0,-1,
     &            c40nam,xth)
      write(c40nam,190)jcx
190   format('  v  cross-section at j = ',i3,8x)
      call mapsmp(scrx2,scrx2r,kx,ix,1,kx,kxn,1,ix,ixn,0.d0,-1,
     &            c40nam,xth)
c
c.....t:
c
      do 200 i=1,ilx
      do 200 k=1,kx
       scrx1(k,i)=ta(i,k,jcx)/psa(i,jcx)
200   continue
      write(c40nam,210)jcx
210   format('  t  cross-section at j = ',i3,8x)
      call mapsmp(scrx1,scrx1r,kx,ix,1,kx,kxn,1,ilx,ixn,273.16d0,
     &            -1,c40nam,xth)
c
c.....qv:
c
      do 220 i=1,ilx
      do 220 k=1,kx
       scrx1(k,i)=qva(i,k,jcx)/psa(i,jcx)
220   continue

      write(c40nam,230)jcx
230   format('  qv  cross-section at j = ',i3,7x)
      call mapsmp(scrx1,scrx1r,kx,ix,1,kx,kxn,1,ilx,ixn,0.d0,-1,
     &            c40nam,xth)
c
c....qc and qr:
c
      do 240 i=1,ilx
      do 240 k=1,kx
       scrx1(k,i)=qca(i,k,jcx)/psa(i,jcx)
240   continue
      write(c40nam,250)jcx
250   format('  qc  cross-section at j = ',i3,7x)
      call mapsmp(scrx1,scrx1r,kx,ix,1,kx,kxn,1,ilx,ixn,0.d0,-1,
     &            c40nam,xth)
c
c-----horizontal slices:
c
      if (jyear.eq.jyear0.and.ktau.eq.0) then
       c40nam='       latitude at cross points        '
       call mapsmp(xlat,xlatr,ix,jx,1,ilx,5,1,jlx,5,0.d0,1,c40nam,xth)
       c40nam='       longitude at cross points       '
       call mapsmp(xlong,xlongr,ix,jx,1,ilx,5,1,jlx,5,0.d0,1,c40nam,xth)
       do 295 j=1,jlx
       do 295 i=1,ilx
        hscr1(i,j)=ht(i,j)/g
295    continue
       c40nam='     terrain height at cross points     '
       call mapsmp(hscr1,hscr1r,ix,jx,1,ilx,5,1,jlx,5,0.d0,1,c40nam,xth)
       c40nam='    coriolis parameter at dot points    '
       call mapsmp(f,fr,ix,jx,1,ix,5,1,jx,5,0.d0,1,c40nam,xth)
       c40nam='    map scale factor at dot points     '
       call mapsmp(msfd,msfdr,ix,jx,1,ix,5,1,jx,5,0.d0,1,c40nam,xth)
      endif
c
c....surface pressure
c
      const = -ptop
      c40nam='         surface pressure (cb)         '
      call mapsmp(psa,psar,ix,jx,1,ilx,ixn,1,jlx,jxn,const,1,c40nam,xth)
c
c.....u and v:
c
      do 430 j=2,jlx
      do 430 i=2,ilx
       psabar=0.25*(psa(i,j)+psa(i,j-1)+psa(i-1,j)+psa(i-1,j-1))
       hscr1(i,j)=ua(i,kout,j)/psabar
       hscr2(i,j)=va(i,kout,j)/psabar
430   continue
c
      if (iexec .eq. 1) then
c
c.....first time:
c
       do 440 i=2,ilx
        psab1=0.5*(psa(i,1)+psa(i-1,1))
        psab2=0.5*(psa(i,jlx)+psa(i-1,jlx))
        hscr1(i,1)=ua(i,kout,1)/psab1
        hscr1(i,jx)=ua(i,kout,jx)/psab2
        hscr2(i,1)=va(i,kout,1)/psab1
        hscr2(i,jx)=va(i,kout,jx)/psab2
440    continue
       do 450 j=2,jlx
        psab1=0.5*(psa(1,j)+psa(1,j-1))
        psab2=0.5*(psa(ilx,j)+psa(ilx,j-1))
        hscr1(1,j)=ua(1,kout,j)/psab1
        hscr1(ix,j)=ua(ix,kout,j)/psab2
        hscr2(1,j)=va(1,kout,j)/psab1
        hscr2(ix,j)=va(ix,kout,j)/psab2
450    continue
       hscr1(1,1)=ua(1,kout,1)/psa(1,1)
       hscr1(1,jx)=ua(1,kout,jx)/psa(1,jlx)
       hscr1(ix,1)=ua(ix,kout,1)/psa(ilx,1)
       hscr1(ix,jx)=ua(ix,kout,jx)/psa(ilx,jlx)
       hscr2(1,1)=va(1,kout,1)/psa(1,1)
       hscr2(1,jx)=va(1,kout,jx)/psa(1,jlx)
       hscr2(ix,1)=va(ix,kout,1)/psa(ilx,1)
       hscr2(ix,jx)=va(ix,kout,jx)/psa(ilx,jlx)
c
      elseif(iexec.gt.1)then
c
       do 470 i=1,ix
        hscr1(i,1)=uj1(i,kout)
        hscr1(i,jx)=ujl(i,kout)
        hscr2(i,1)=vj1(i,kout)
        hscr2(i,jx)=vjl(i,kout)
470    continue
       do 480 j=2,jlx
        hscr1(1,j)=ui1(kout,j)
        hscr1(ix,j)=uil(kout,j)
        hscr2(1,j)=vi1(kout,j)
        hscr2(ix,j)=vil(kout,j)
480    continue
c
      endif      !end 2nd iexec test

      write(c40nam,500)kout
500   format('  u   at k = ',i3,24x)
      call mapsmp(hscr1,hscr1r,ix,jx,1,ix,ixn,1,jx,jxn,0.d0,1,
     &            c40nam,xth)
      write(c40nam,510)kout
510   format('  v   at k = ',i3,24x)
      call mapsmp(hscr2,hscr2r,ix,jx,1,ix,ixn,1,jx,jxn,0.d0,1,
     &            c40nam,xth)
c
c.....t:
c
      do 520 j=1,jlx
      do 520 i=1,ilx
       hscr1(i,j)=ta(i,kout,j)/psa(i,j)
520   continue
      write(c40nam,530)kout
530   format('  t (c)  at k = ',i3,21x)
      call mapsmp(hscr1,hscr1r,ix,jx,1,ilx,ixn,1,jlx,jxn,273.16d0,
     &            1,c40nam,xth)
c
c.....relative humidity:
c
      do 540 j=1,jlx
      do 540 i=1,ilx
       pres=a(kout)*psa(i,j)+ptop
       if(hscr1(i,j).gt.273.15)then
c     v8 svp formula
           satvp=svp1*dexp(svp2*(hscr1(i,j)-273.15)/(hscr1(i,j)-svp3))
        else
           satvp=.611*dexp(22.514-6.15e3/hscr1(i,j))
        endif

       qvs=ep2*satvp/(pres-satvp)
       hscr2(i,j)=(qva(i,kout,j)/psa(i,j))/qvs
540   continue
      write(c40nam,550)kout
550   format('  relative humidity at k = ',i3,10x)
      call mapsmp(hscr2,hscr2r,ix,jx,1,ilx,ixn,1,jlx,jxn,0.d0,1,
     1            c40nam,xth)
c
c.....qc and qr:
c
      do 560 j=1,jlx
      do 560 i=1,ilx
       hscr1(i,j)=qca(i,kout,j)/psa(i,j)
560   continue
      write(c40nam,570)kout
570   format('  qc  at k = ',i3,24x)
      call mapsmp(hscr1,hscr1r,ix,jx,1,ilx,ixn,1,jlx,jxn,0.d0,1,
     1            c40nam,xth)
c
c.....surface temperature:
c
      if (ibltyp .ne. 0) then
       const = 273.15
       c40nam =' ground temperature (c)    '
       call mapsmp(tga,tgar,ix,jx,1,ilx,ixn,1,jlx,jxn,const,1,
     1             c40nam,xth)
      endif
c
c.....precipitation:
c
      c40nam='   convective rainfall (mm)          '
      call mapsmp(rainc,raincr,ix,jx,1,ilx,ixn,1,jlx,jxn,0.d0,1,
     1             c40nam,xth)
      c40nam='  nonconvective rainfall (mm)      '
      call mapsmp(rainnc,rainncr,ix,jx,1,ilx,ixn,1,jlx,jxn,0.d0,1,
     1            c40nam,xth)
      do 770 j=1,jlx
      do 770 i=1,ilx
       hscr1(i,j)=rainc(i,j)+rainnc(i,j)
770   continue
      c40nam='  total rainfall (mm)        '
      call mapsmp(hscr1,hscr1r,ix,jx,1,ilx,ixn,1,jlx,jxn,0.d0,1,
     1            c40nam,xth)
c
      return
      end
#endif
