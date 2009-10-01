ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bdyval(xt,dtmn,iexec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine sets the boundary values for p*, p*u, p*v,      c
c     p*t, p*qv, p*qc, and p*qr.                                      c
c                                                                     c
c     ---the boundary values of p*u and p*v are extrapolated from     c
c        the interior points.                                         c
c                                                                     c
c     ---the boundary values of p* and p*t are specified.             c
c                                                                     c
c     ---the boundary values of p*qv, p*qc, and p*qr depend on        c
c        inflow/outflow conditions, if iboudy = 3 or 4.               c
c                                                                     c
c     xt     : is the time in minutes the variables xxa represent.    c
c                                                                     c
c     dtmn  : is the time step in minutes.                            c
c                                                                     c
c     iexec  : = 1 ; represents this subroutine is called for the     c
c                    first time in this forecast run.                 c
c              > 1 ; represents subsequent calls to this subroutine.  c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer iexec
      real*8  xt,dtmn
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
      include 'Commons/iunits.cb'
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
      include 'Common2/bdycod.cb'
      include 'Common2/trachem.cb'
      include 'Common2/cvaria.cb'
      include 'mpiregcm.h'
#else
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
      include 'Commons/bdycod.cb'
      include 'Commons/trachem.cb'
      include 'Commons/cvaria.cb'
#endif
      include 'Commons/date.cb'
      integer i,j,k
      real*8  dtb,qvx1,qvx2,uavg,qvx,vavg,qcx2,qcx
Cchem2
      integer itr
      real*8  chix1,chix2,chix
Cchem2_
c
c*********************************************************************
c*****fill up the boundary value for xxb variables from xxa variables:
c     if this subroutine is called for the first time, this part
c     shall be skipped.
c
      if (iexec .ne. 1) then
c
c-----for p*:
c
       do 10 i=1,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        psb(i,1) = psa(i,1)
#ifdef MPP1
        endif
        if (myid .eq. nproc-1) then
          psb(i,jendx) = psa(i,jendx)
        endif
#else
        psb(i,jlx) = psa(i,jlx)
#endif
10     continue
#ifdef MPP1
       do 20 j=jbegin,jendm
#else
       do 20 j=2,jlxm
#endif
        psb(1,j) = psa(1,j)
        psb(ilx,j) = psa(ilx,j)
20     continue
c
c-----for p*u and p*v:
c
       do 40 k=1,kx
        do 30 i=1,ix
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         ub(i,k,1)=ua(i,k,1)/msfd(i,1)
         vb(i,k,1)=va(i,k,1)/msfd(i,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           ub(i,k,jendl)=ua(i,k,jendl)/msfd(i,jendl)
           vb(i,k,jendl)=va(i,k,jendl)/msfd(i,jendl)
         endif
#else
         ub(i,k,jx)=ua(i,k,jx)/msfd(i,jx)
         vb(i,k,jx)=va(i,k,jx)/msfd(i,jx)
#endif
30      continue
#ifdef MPP1
       do 40 j=jbegin,jendx
#else
       do 40 j=2,jlx
#endif
        ub(1,k,j)=ua(1,k,j)/msfd(1,j)
        ub(ix,k,j)=ua(ix,k,j)/msfd(ix,j)
        vb(1,k,j)=va(1,k,j)/msfd(1,j)
        vb(ix,k,j)=va(ix,k,j)/msfd(ix,j)
40     continue
c
c-----for p*t:
c
       do 60 k=1,kx
        do 50 i=1,ilx
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         tb(i,k,1)=ta(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           tb(i,k,jendx)=ta(i,k,jendx)
         endif
#else
         tb(i,k,jlx)=ta(i,k,jlx)
#endif
50      continue
#ifdef MPP1
       do 60 j=jbegin,jendm
#else
       do 60 j=2,jlxm
#endif
        tb(1,k,j)=ta(1,k,j)
        tb(ilx,k,j)=ta(ilx,k,j)
60     continue
c
c-----for p*qv:
c
       do 80 k=1,kx
        do 70 i=1,ilx
#ifdef MPP1
         if(myid .eq. 0) then
#endif
         qvb(i,k,1)=qva(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           qvb(i,k,jendx)=qva(i,k,jendx)
         endif
#else
         qvb(i,k,jlx)=qva(i,k,jlx)
#endif
70      continue
#ifdef MPP1
       do 80 j=jbegin,jendm
#else
       do 80 j=2,jlxm
#endif
        qvb(1,k,j)=qva(1,k,j)
        qvb(ilx,k,j)=qva(ilx,k,j)
80     continue
c
Cchem2
c
       if(ichem.eq.1) then
c-----for p*chi (tracers)
       do 85 itr = 1,ntr
       do 85 k=1,kx
         do 83 i=1,ilx
#ifdef MPP1
          if (myid .eq. 0) then
#endif
          chib(i,k,1,itr)=chia(i,k,1,itr)
#ifdef MPP1
          endif
          if (myid .eq. nproc-1) then
            chib(i,k,jendx,itr)=chia(i,k,jendx,itr)
          endif
#else
          chib(i,k,jlx,itr)=chia(i,k,jlx,itr)
#endif
83       continue
#ifdef MPP1
        do 85 j=jbegin,jendm
#else
        do 85 j=2,jlxm
#endif
         chib(1,k,j,itr)=chia(1,k,j,itr)
         chib(ilx,k,j,itr)=chia(ilx,k,j,itr)
85      continue
       endif
Cchem2_
c
c-----for p*qc:
c
       do 100 k=1,kx
        do 90 i=1,ilx
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         qcb(i,k,1)=qca(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           qcb(i,k,jendx)=qca(i,k,jendx)
         endif
#else
         qcb(i,k,jlx)=qca(i,k,jlx)
#endif
90      continue
#ifdef MPP1
       do 100 j=jbegin,jendm
#else
       do 100 j=2,jlxm
#endif
        qcb(1,k,j)=qca(1,k,j)
        qcb(ilx,k,j)=qca(ilx,k,j)
100    continue
c
      endif      !end if(iexec.ne.1) test
c**********************************************************************
c*****compute the boundary values for xxa variables:
c
c-----compute the time interval for boundary tendency:
c
      dtb=xt*60.
      if(dabs(xt).lt.0.00001.and.ldatez.gt.IDATE0) dtb=IBDYFRQ*60.*60.
c
c-----set boundary values for p*:
c-----set boundary conditions for p*u and p*v:
c
      if (iexec .eq. 1 .and. ifrest) goto 310
c
      if (iboudy .eq. 0)  then
c.....fixed boundary conditions:
       do 160 i=1,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        psa(i,1)=pwb(i,1)
#ifdef MPP1
        endif
        if (myid .eq. nproc-1) then
          psa(i,jendx)=peb(i,1)
        endif
#else
        psa(i,jlx)=peb(i,1)
#endif
160    continue
#ifdef MPP1
       do 170 j=jbegin,jendm
#else
       do 170 j=2,jlxm
#endif
        psa(1,j)=pss(1,j)
        psa(ilx,j)=pnb(1,j)
170    continue
c
       do 270 k=1,kx
        do 260 i=1,ix
#ifdef MPP1
         if (myid .eq.0) then
#endif
         ua(i,k,1)=uwb(i,k,1)
         va(i,k,1)=vwb(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           ua(i,k,jendl)=ueb(i,k,1)
           va(i,k,jendl)=veb(i,k,1)
         endif
#else
         ua(i,k,jx)=ueb(i,k,1)
         va(i,k,jx)=veb(i,k,1)
#endif
260     continue
#ifdef MPP1
       do 270 j=jbegin,jendx
#else
       do 270 j=2,jlx
#endif
        ua(1,k,j)=usb(1,k,j)
        ua(ix,k,j)=unb(1,k,j)
        va(1,k,j)=vsb(1,k,j)
        va(ix,k,j)=vnb(1,k,j)
270    continue
      endif
c
c.....time-dependent boundary conditions:
c
      do 190 i=1,ilx
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       psa(i,1)=pwb(i,1)+dtb*pwbt(i,1)
#ifdef MPP1
       endif
       if (myid .eq. nproc-1 ) then
         psa(i,jendx)=peb(i,1)+dtb*pebt(i,1)
       endif
#else
       psa(i,jlx)=peb(i,1)+dtb*pebt(i,1)
#endif
190   continue
#ifdef MPP1
      do 200 j=jbegin,jendm
#else
      do 200 j=2,jlxm
#endif
       psa(1,j)=pss(1,j)+dtb*psbt(1,j)
       psa(ilx,j)=pnb(1,j)+dtb*pnbt(1,j)
200   continue
c
      do 300 k=1,kx
       do 290 i=1,ix
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        ua(i,k,1)=uwb(i,k,1)+dtb*uwbt(i,k,1)
        va(i,k,1)=vwb(i,k,1)+dtb*vwbt(i,k,1)
#ifdef MPP1
        endif
        if (myid .eq. nproc-1) then
          ua(i,k,jendl)=ueb(i,k,1)+dtb*uebt(i,k,1)
          va(i,k,jendl)=veb(i,k,1)+dtb*vebt(i,k,1)
        endif
#else
        ua(i,k,jx)=ueb(i,k,1)+dtb*uebt(i,k,1)
        va(i,k,jx)=veb(i,k,1)+dtb*vebt(i,k,1)
#endif
290    continue
#ifdef MPP1
      do 300 j=jbegin,jendx
#else
      do 300 j=2,jlx
#endif
       ua(1,k,j)=usb(1,k,j)+dtb*usbt(1,k,j)
       ua(ix,k,j)=unb(1,k,j)+dtb*unbt(1,k,j)
       va(1,k,j)=vsb(1,k,j)+dtb*vsbt(1,k,j)
       va(ix,k,j)=vnb(1,k,j)+dtb*vnbt(1,k,j)
300   continue
c
c-----get boundary values of u and v:
c
310   call bdyuv(iboudy,dtb)
c
      if (iexec .eq. 1 .and. ifrest) return
c
c-----set boundary values for p*t:
c-----set boundary values for p*qv:
c
      if (iboudy .eq. 0) then
c.....fixed boundary conditions:
       do 340 k=1,kx
        do 330 i=1,ilx
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         ta(i,k,1)=twb(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           ta(i,k,jendx)=teb(i,k,1)
         endif
#else
         ta(i,k,jlx)=teb(i,k,1)
#endif
330     continue
#ifdef MPP1
       do 340 j=jbegin,jendm
#else
       do 340 j=2,jlxm
#endif
        ta(1,k,j)=tsb(1,k,j)
        ta(ilx,k,j)=tnb(1,k,j)
340    continue
       do 410 k=1,kx
        do 400 i=1,ilx
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         qva(i,k,1)=qwb(i,k,1)
#ifdef MPP1
         endif
         if (myid .eq. nproc-1) then
           qva(i,k,jendx)=qeb(i,k,1)
         endif
#else
         qva(i,k,jlx)=qeb(i,k,1)
#endif
400     continue
#ifdef MPP1
       do 410 j=jbegin,jendm
#else
       do 410 j=2,jlxm
#endif
        qva(1,k,j)=qsb(1,k,j)
        qva(ilx,k,j)=qnb(1,k,j)
410    continue
      endif
c
c.....time-dependent boundary conditions:
c
      do 370 k=1,kx
       do 360 i=1,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        ta(i,k,1)=twb(i,k,1)+dtb*twbt(i,k,1)
#ifdef MPP1
        endif
        if (myid.eq. nproc-1) then
          ta(i,k,jendx)=teb(i,k,1)+dtb*tebt(i,k,1)
        endif
#else
        ta(i,k,jlx)=teb(i,k,1)+dtb*tebt(i,k,1)
#endif
360    continue
#ifdef MPP1
      do 370 j=jbegin,jendm
#else
      do 370 j=2,jlxm
#endif
       ta(1,k,j)=tsb(1,k,j)+dtb*tsbt(1,k,j)
       ta(ilx,k,j)=tnb(1,k,j)+dtb*tnbt(1,k,j)
370   continue
      do 440 k=1,kx
       do 430 i=1,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        qva(i,k,1)=qwb(i,k,1)+dtb*qwbt(i,k,1)
#ifdef MPP1
        endif
        if (myid.eq. nproc-1) then
          qva(i,k,jendx)=qeb(i,k,1)+dtb*qebt(i,k,1)
        endif
#else
        qva(i,k,jlx)=qeb(i,k,1)+dtb*qebt(i,k,1)
#endif
430    continue
#ifdef MPP1
      do 440 j=jbegin,jendm
#else
      do 440 j=2,jlxm
#endif
       qva(1,k,j)=qsb(1,k,j)+dtb*qsbt(1,k,j)
       qva(ilx,k,j)=qnb(1,k,j)+dtb*qnbt(1,k,j)
440   continue
c
      if (iboudy .eq. 3. or. iboudy.eq.4) then
c
c-----determine boundary values depends on inflow/outflow:
c
       do 490 k=1,kx
c
c.....west boundary:
c
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        do 450 i=1,ilx
         qvx1 = qva(i,k,1)/psa(i,1)
         qvx2 = qva(i,k,2)/psa(i,2)
         uavg=uj1(i,k)+uj1(i+1,k)+uj2(i,k)+uj2(i+1,k)
         if(uavg.ge.0.)then
          qvx=qvx1
         else
          qvx=qvx2
         endif
         qva(i,k,1)=qvx*psa(i,1)
450     continue
#ifdef MPP1
        endif
#endif
c
c.....east boundary:
c
#ifdef MPP1
         if (myid .eq. nproc-1) then
#endif
        do 460 i=1,ilx
#ifdef MPP1
         qvx1 = qva(i,k,jendx)/psa(i,jendx)
         qvx2 = qva(i,k,jendm)/psa(i,jendm)
#else
         qvx1 = qva(i,k,jlx)/psa(i,jlx)
         qvx2 = qva(i,k,jlxm)/psa(i,jlxm)
#endif
         uavg=ujlx(i,k)+ujlx(i+1,k)+ujl(i,k)+ujl(i+1,k)
         if(uavg.lt.0.)then
          qvx=qvx1
         else
          qvx=qvx2
         endif
#ifdef MPP1
         qva(i,k,jendx)=qvx*psa(i,jendx)
#else
         qva(i,k,jlx)=qvx*psa(i,jlx)
#endif
460     continue
#ifdef MPP1
        endif
#endif
c
c.....south boundary:
c
#ifdef MPP1
        do 470 j=jbegin,jendm
#else
        do 470 j=2,jlxm
#endif
         qvx1 = qva(1,k,j)/psa(1,j)
         qvx2 = qva(2,k,j)/psa(2,j)
         vavg=vi1(k,j)+vi1(k,j+1)+vi2(k,j)+vi2(k,j+1)
         if(vavg.ge.0.)then
          qvx=qvx1
         else
          qvx=qvx2
         endif
         qva(1,k,j) = qvx*psa(1,j)
470     continue
c
c.....north boundary:
c
#ifdef MPP1
        do 480 j=jbegin,jendm
#else
        do 480 j=2,jlxm
#endif
         qvx1 = qva(ilx,k,j)/psa(ilx,j)
         qvx2 = qva(ilxm,k,j)/psa(ilxm,j)
         vavg=vilx(k,j)+vilx(k,j+1)+vil(k,j)+vil(k,j+1)
         if(vavg.lt.0.)then
          qvx=qvx1
         else
          qvx=qvx2
         endif
         qva(ilx,k,j)=qvx*psa(ilx,j)
480     continue
c
490    continue
      endif      !end if(iboudy.eq.3.or.4) test
c
c-----set boundary values for p*qc and p*qr:
c     *** note ***
c         for large domain, we assume the boundary tendencies are not
c         available.
c
c
c-----if the boundary values and tendencies are not available,
c     determine boundary values depends on inflow/outflow:
c        inflow  : set it equal to zero.
c        outflow : get from interior point.
c
      do 730 k=1,kx
c
c.....west boundary:
c
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       do 690 i=1,ilx
        qcx2 = qca(i,k,2)/psa(i,2)
        uavg=uj1(i,k)+uj1(i+1,k)+uj2(i,k)+uj2(i+1,k)
        if(uavg.ge.0.)then
         qcx=0.
        else
         qcx=qcx2
        endif
        qca(i,k,1)=qcx*psa(i,1)
690    continue
#ifdef MPP1
        endif
#endif
c
c.....east boundary:
c
#ifdef MPP1
       if (myid .eq. nproc-1 ) then
#endif
       do 700 i=1,ilx
#ifdef MPP1
        qcx2 = qca(i,k,jendm)/psa(i,jendm)
#else
        qcx2 = qca(i,k,jlxm)/psa(i,jlxm)
#endif
        uavg=ujlx(i,k)+ujlx(i+1,k)+ujl(i,k)+ujl(i+1,k)
        if(uavg.lt.0.)then
         qcx=0.
        else
         qcx=qcx2
        endif
#ifdef MPP1
        qca(i,k,jendx)=qcx*psa(i,jendx)
#else
        qca(i,k,jlx)=qcx*psa(i,jlx)
#endif
700    continue
#ifdef MPP1
        endif
#endif
c
c.....south boundary:
c
#ifdef MPP1
       do 710 j=jbegin,jendm
#else
       do 710 j=2,jlxm
#endif
        qcx2 = qca(2,k,j)/psa(2,j)
        vavg=vi1(k,j)+vi1(k,j+1)+vi2(k,j)+vi2(k,j+1)
        if(vavg.ge.0.)then
         qcx=0.
        else
         qcx=qcx2
        endif
        qca(1,k,j) = qcx*psa(1,j)
710    continue
c
c.....north boundary:
c
#ifdef MPP1
       do 720 j=jbegin,jendm
#else
       do 720 j=2,jlxm
#endif
        qcx2 = qca(ilxm,k,j)/psa(ilxm,j)
        vavg=vilx(k,j)+vilx(k,j+1)+vil(k,j)+vil(k,j+1)
        if(vavg.lt.0.)then
         qcx=0.
        else
         qcx=qcx2
        endif
        qca(ilx,k,j)=qcx*psa(ilx,j)
720    continue
c
730   continue

      if(ichem.eq.1) then
Cchem2

c----add tracer bc's
c
      do 500 itr=1,ntr
      do 500 k=1,kx
c
c.....west boundary:
c
       
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        do 493 i=1,ilx
         chix1 = chia(i,k,1,itr)/psa(i,1)
         chix2 = chia(i,k,2,itr)/psa(i,2)
         uavg=uj1(i,k)+uj1(i+1,k)+uj2(i,k)+uj2(i+1,k)
         if(uavg.ge.0.)then
          chix=chix1
         else
          chix=chix2
         endif
         chia(i,k,1,itr)=chix*psa(i,1)
493     continue
#ifdef MPP1
         endif
#endif
c
c.....east boundary:
c
#ifdef MPP1
        if (myid .eq. nproc-1) then
#endif
        do 495 i=1,ilx
#ifdef MPP1
         chix1 = chia(i,k,jendx,itr)/psa(i,jendx)
         chix2 = chia(i,k,jendm,itr)/psa(i,jendm)
#else
         chix1 = chia(i,k,jlx,itr)/psa(i,jlx)
         chix2 = chia(i,k,jlxm,itr)/psa(i,jlxm)
#endif
         uavg=ujlx(i,k)+ujlx(i+1,k)+ujl(i,k)+ujl(i+1,k)
         if(uavg.lt.0.)then
          chix=chix1
         else
          chix=chix2
         endif
#ifdef MPP1
         chia(i,k,jendx,itr)=chix*psa(i,jendx)
#else
         chia(i,k,jlx,itr)=chix*psa(i,jlx)
#endif
495     continue
#ifdef MPP1
         endif
#endif
c
c.....south boundary:
c
#ifdef MPP1
        do 497 j=jbegin,jendm
#else
        do 497 j=2,jlxm
#endif
         chix1 = chia(1,k,j,itr)/psa(1,j)
         chix2 = chia(2,k,j,itr)/psa(2,j)
         vavg=vi1(k,j)+vi1(k,j+1)+vi2(k,j)+vi2(k,j+1)
         if(vavg.ge.0.)then
          chix=chix1
         else
          chix=chix2
         endif
         chia(1,k,j,itr) = chix*psa(1,j)
497     continue
c
c.....north boundary:
c
#ifdef MPP1
        do 499 j=jbegin,jendm
#else
        do 499 j=2,jlxm
#endif
         chix1 = chia(ilx,k,j,itr)/psa(ilx,j)
         chix2 = chia(ilxm,k,j,itr)/psa(ilxm,j)
         vavg=vilx(k,j)+vilx(k,j+1)+vil(k,j)+vil(k,j+1)
         if(vavg.lt.0.)then
          chix=chix1
         else
          chix=chix2
         endif
         chia(ilx,k,j,itr)=chix*psa(ilx,j)
499     continue
500     continue
Cchem2_
      endif
c
      return
      end
