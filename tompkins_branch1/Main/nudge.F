ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nudge_P(ip,fcoef,gcoef,tbe,xt,pten,c203,j,iboudy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies relaxation boundary conditions to the   c
c     tendency term - pten.                                           c
c                                                                     c
c     ip    : is the number of slices affected by nudging.            c
c                                                                     c
c     tbe   : is the time for the boundary values (peb, pwb, pss,     c
c             and pnb) in minutes.                                    c
c                                                                     c
c     xt    : is the time in minutes for variable "psb".              c
c                                                                     c
c     fcoef : are the coefficients for the newtonian term.            c
c                                                                     c
c     gcoef : are the coefficients for the diffusion term.            c
c                                                                     c
c     pten  : is the tendency calculated from the model.              c
c                                                                     c
c     peb, pwb, pss, pnb : are the observed boundary values at time   c
c             "tbe" on east, west, south, and north boundaries.       c
c                                                                     c
c     pebt, pwbt, psbt, pnbt : are the large-scale or observed        c
c             tendencies at east, west, south, and north boundaries.  c
c                                                                     c
c     psb    : is the variable at tau-1.                               c
c                                                                     c
c     c203  : = 1./(dx*dx), defined in 'param'.                       c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c     iboudy : type of boundary condition relaxation, 1=linear        c
c              5 = exponential                                        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j,iboudy
      real*8  fcoef,gcoef,tbe,xt,c203
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/bdycod.cb'
      include 'Common2/main.cb'
      integer jwb,jeb
#else
      include 'Commons/bdycod.cb'
      include 'Commons/main.cb'
#endif
      include 'Commons/param3.cb'
c
      real*8  pten(ix)
      integer i,ii,jj,jsls,ibeg,iend
      real*8  dtb,fcx,gcx,fls0,fls1,fls2,fls3,fls4
c
      integer mm,kk
      real*8  xfun,xfune
      xfun(mm) = dble(nspgd-mm)/(nspgd-2.)
      xfune(mm,kk)= dexp(-dble(mm-2)/anudg(kk))
c
c----------------------------------------------------------------------
c
      dtb=xt*60.
#ifdef MPP1
      jsls=j+myid*jxp
      jj=mjx-jsls
#else
      jsls=j
      jj=jx-jsls
#endif
      if (jj .le. ip) jsls=jj
#ifdef MPP1
      jwb=jsls
      if(jwb.gt.jxp) jwb=mod(jwb,jxp)
      if(jwb.eq.0) jwb=jxp
      if(myid.eq.nproc-1) then
        jeb=jsls
      else
        jeb=jsls+1
      endif
      if(jeb.gt.jxp) jeb=mod(jeb,jxp)
      if(jeb.eq.0) jeb=jxp
#endif
c-----determine which relaxation method to use:linear/expon.
c
      if(iboudy.eq.1)then
c
c---------use linear method
c
       if (jsls .gt. ip) then
c------interior j slices:
        do 10 i=2,ip
         ii=ilx-i+1
         fcx=fcoef*xfun(i)
         gcx=gcoef*xfun(i)
c.......south boundary:
         fls0=(pss(i,j)+dtb*psbt(i,j))-psb(i,j)
         fls1=(pss(i,j-1)+dtb*psbt(i,j-1))-psb(i,j-1)
         fls2=(pss(i,j+1)+dtb*psbt(i,j+1))-psb(i,j+1)
         fls3=(pss(i-1,j)+dtb*psbt(i-1,j))-psb(i-1,j)
         fls4=(pss(i+1,j)+dtb*psbt(i+1,j))-psb(i+1,j)
         pten(i)=pten(i)+fcx*fls0-gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(pnb(i,j)+dtb*pnbt(i,j))-psb(ii,j)
         fls1=(pnb(i,j-1)+dtb*pnbt(i,j-1))-psb(ii,j-1)
         fls2=(pnb(i,j+1)+dtb*pnbt(i,j+1))-psb(ii,j+1)
         fls3=(pnb(i-1,j)+dtb*pnbt(i-1,j))-psb(ii-1,j)
         fls4=(pnb(i+1,j)+dtb*pnbt(i+1,j))-psb(ii+1,j)
         pten(ii)=pten(ii)+fcx*fls0
     &           -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
10      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 30 i=2,jsls-1
          ii=ilx-i+1
          fcx=fcoef*xfun(i)
          gcx=gcoef*xfun(i)
c........south boundary:
          fls0=(pss(i,j)+dtb*psbt(i,j))-psb(i,j)
          fls1=(pss(i,j-1)+dtb*psbt(i,j-1))-psb(i,j-1)
          fls2=(pss(i,j+1)+dtb*psbt(i,j+1))-psb(i,j+1)
          fls3=(pss(i-1,j)+dtb*psbt(i-1,j))-psb(i-1,j)
          fls4=(pss(i+1,j)+dtb*psbt(i+1,j))-psb(i+1,j)
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(pnb(i,j)+dtb*pnbt(i,j))-psb(ii,j)
          fls1=(pnb(i,j-1)+dtb*pnbt(i,j-1))-psb(ii,j-1)
          fls2=(pnb(i,j+1)+dtb*pnbt(i,j+1))-psb(ii,j+1)
          fls3=(pnb(i-1,j)+dtb*pnbt(i-1,j))-psb(ii-1,j)
          fls4=(pnb(i+1,j)+dtb*pnbt(i+1,j))-psb(ii+1,j)
          pten(ii)=pten(ii)+fcx*fls0
     &                     -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
30       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 50 i=ibeg,iend
#ifdef MPP1
          fls0=(pwb(i,jwb)+dtb*pwbt(i,jwb))-psb(i,j)
          fls1=(pwb(i-1,jwb)+dtb*pwbt(i-1,jwb))-psb(i-1,j)
          fls2=(pwb(i+1,jwb)+dtb*pwbt(i+1,jwb))-psb(i+1,j)
          fls3=(pwb(i,jwb-1)+dtb*pwbt(i,jwb-1))-psb(i,j-1)
          fls4=(pwb(i,jwb+1)+dtb*pwbt(i,jwb+1))-psb(i,j+1)
#else
          fls0=(pwb(i,jsls)+dtb*pwbt(i,jsls))-psb(i,j)
          fls1=(pwb(i-1,jsls)+dtb*pwbt(i-1,jsls))-psb(i-1,j)
          fls2=(pwb(i+1,jsls)+dtb*pwbt(i+1,jsls))-psb(i+1,j)
          fls3=(pwb(i,jsls-1)+dtb*pwbt(i,jsls-1))-psb(i,j-1)
          fls4=(pwb(i,jsls+1)+dtb*pwbt(i,jsls+1))-psb(i,j+1)
#endif
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
50       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 70 i=ibeg,iend
#ifdef MPP1
          fls0=(peb(i,jeb)+dtb*pebt(i,jeb))-psb(i,j)
          fls1=(peb(i-1,jeb)+dtb*pebt(i-1,jeb))-psb(i-1,j)
          fls2=(peb(i+1,jeb)+dtb*pebt(i+1,jeb))-psb(i+1,j)
          fls3=(peb(i,jeb-1)+dtb*pebt(i,jeb-1))-psb(i,j-1)
          fls4=(peb(i,jeb+1)+dtb*pebt(i,jeb+1))-psb(i,j+1)
#else
          fls0=(peb(i,jsls)+dtb*pebt(i,jsls))-psb(i,j)
          fls1=(peb(i-1,jsls)+dtb*pebt(i-1,jsls))-psb(i-1,j)
          fls2=(peb(i+1,jsls)+dtb*pebt(i+1,jsls))-psb(i+1,j)
          fls3=(peb(i,jsls-1)+dtb*pebt(i,jsls-1))-psb(i,j-1)
          fls4=(peb(i,jsls+1)+dtb*pebt(i,jsls+1))-psb(i,j+1)
#endif
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
70       continue
        endif
       endif
c
      elseif(iboudy.eq.5)then

c----------use exponential method

       if (jsls .gt. ip) then
c------interior j slices:
        do 75 i=2,ip
         ii=ilx-i+1
         fcx=fcoef*xfune(i,kx)
         gcx=gcoef*xfune(i,kx)
c........south boundary:
         fls0=(pss(i,j)+dtb*psbt(i,j))-psb(i,j)
         fls1=(pss(i,j-1)+dtb*psbt(i,j-1))-psb(i,j-1)
         fls2=(pss(i,j+1)+dtb*psbt(i,j+1))-psb(i,j+1)
         fls3=(pss(i-1,j)+dtb*psbt(i-1,j))-psb(i-1,j)
         fls4=(pss(i+1,j)+dtb*psbt(i+1,j))-psb(i+1,j)
         pten(i)=pten(i)+fcx*fls0
     &                  -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(pnb(i,j)+dtb*pnbt(i,j))-psb(ii,j)
         fls1=(pnb(i,j-1)+dtb*pnbt(i,j-1))-psb(ii,j-1)
         fls2=(pnb(i,j+1)+dtb*pnbt(i,j+1))-psb(ii,j+1)
         fls3=(pnb(i-1,j)+dtb*pnbt(i-1,j))-psb(ii-1,j)
         fls4=(pnb(i+1,j)+dtb*pnbt(i+1,j))-psb(ii+1,j)
         pten(ii)=pten(ii)+fcx*fls0
     &                    -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
75      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 85 i=2,jsls-1
          ii=ilx-i+1
          fcx=fcoef*xfune(i,kx)
          gcx=gcoef*xfune(i,kx)
c.........south boundary:
          fls0=(pss(i,j)+dtb*psbt(i,j))-psb(i,j)
          fls1=(pss(i,j-1)+dtb*psbt(i,j-1))-psb(i,j-1)
          fls2=(pss(i,j+1)+dtb*psbt(i,j+1))-psb(i,j+1)
          fls3=(pss(i-1,j)+dtb*psbt(i-1,j))-psb(i-1,j)
          fls4=(pss(i+1,j)+dtb*psbt(i+1,j))-psb(i+1,j)
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(pnb(i,j)+dtb*pnbt(i,j))-psb(ii,j)
          fls1=(pnb(i,j-1)+dtb*pnbt(i,j-1))-psb(ii,j-1)
          fls2=(pnb(i,j+1)+dtb*pnbt(i,j+1))-psb(ii,j+1)
          fls3=(pnb(i-1,j)+dtb*pnbt(i-1,j))-psb(ii-1,j)
          fls4=(pnb(i+1,j)+dtb*pnbt(i+1,j))-psb(ii+1,j)
          pten(ii)=pten(ii)+fcx*fls0
     &                     -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
85       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfune(jsls,kx)
         gcx=gcoef*xfune(jsls,kx)
         do 95 i=ibeg,iend
#ifdef MPP1
          fls0=(pwb(i,jwb)+dtb*pwbt(i,jwb))-psb(i,j)
          fls1=(pwb(i-1,jwb)+dtb*pwbt(i-1,jwb))-psb(i-1,j)
          fls2=(pwb(i+1,jwb)+dtb*pwbt(i+1,jwb))-psb(i+1,j)
          fls3=(pwb(i,jwb-1)+dtb*pwbt(i,jwb-1))-psb(i,j-1)
          fls4=(pwb(i,jwb+1)+dtb*pwbt(i,jwb+1))-psb(i,j+1)
#else
          fls0=(pwb(i,jsls)+dtb*pwbt(i,jsls))-psb(i,j)
          fls1=(pwb(i-1,jsls)+dtb*pwbt(i-1,jsls))-psb(i-1,j)
          fls2=(pwb(i+1,jsls)+dtb*pwbt(i+1,jsls))-psb(i+1,j)
          fls3=(pwb(i,jsls-1)+dtb*pwbt(i,jsls-1))-psb(i,j-1)
          fls4=(pwb(i,jsls+1)+dtb*pwbt(i,jsls+1))-psb(i,j+1)
#endif
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
95       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfune(jsls,kx)
         gcx=gcoef*xfune(jsls,kx)
         do 105 i=ibeg,iend
#ifdef MPP1
          fls0=(peb(i,jeb)+dtb*pebt(i,jeb))-psb(i,j)
          fls1=(peb(i-1,jeb)+dtb*pebt(i-1,jeb))-psb(i-1,j)
          fls2=(peb(i+1,jeb)+dtb*pebt(i+1,jeb))-psb(i+1,j)
          fls3=(peb(i,jeb-1)+dtb*pebt(i,jeb-1))-psb(i,j-1)
          fls4=(peb(i,jeb+1)+dtb*pebt(i,jeb+1))-psb(i,j+1)
#else
          fls0=(peb(i,jsls)+dtb*pebt(i,jsls))-psb(i,j)
          fls1=(peb(i-1,jsls)+dtb*pebt(i-1,jsls))-psb(i-1,j)
          fls2=(peb(i+1,jsls)+dtb*pebt(i+1,jsls))-psb(i+1,j)
          fls3=(peb(i,jsls-1)+dtb*pebt(i,jsls-1))-psb(i,j-1)
          fls4=(peb(i,jsls+1)+dtb*pebt(i,jsls+1))-psb(i,j+1)
#endif
          pten(i)=pten(i)+fcx*fls0
     &                   -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
105      continue
        endif
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nudge_T(ip,fcoef,gcoef,tbe,xt,ften,c203,j,iboudy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies relaxation boundary conditions to the   c
c     tendency term - ften.                                           c
c                                                                     c
c     ip    : is the number of slices affected by nudging.            c
c                                                                     c
c     tbe   : is the time for the boundary values (teb, twb, tsb,     c
c             and tnb) in minutes.                                    c
c                                                                     c
c     xt    : is the time in minutes for variable "tb".               c
c                                                                     c
c     fcoef : are the coefficients for the newtonian term.            c
c                                                                     c
c     gcoef : are the coefficients for the diffusion term.            c
c                                                                     c
c     ften  : is the tendency calculated from the model.              c
c                                                                     c
c     teb, twb, tsb, tnb : are the observed boundary values at time   c
c             "tbe" on east, west, south, and north boundaries.       c
c                                                                     c
c     tebt, twbt, tsbt, tnbt : are the large-scale or observed        c
c             tendencies at east, west, south, and north boundaries.  c
c                                                                     c
c     tb    : is the variable at tau-1.                               c
c                                                                     c
c     c203  : = 1./(dx*dx), defined in 'param'.                       c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c     iboudy : type of boundary condition relaxation, 1=linear        c
c              5 = exponential                                        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j,iboudy
      real*8  fcoef,gcoef,tbe,xt,c203
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/main.cb'
      include 'Common2/bdycod.cb'
      integer jwb,jeb
#else
      include 'Commons/main.cb'
      include 'Commons/bdycod.cb'
#endif
      include 'Commons/param3.cb'
c
      real*8  ften(ix,kx)
      integer i,ii,k,jj,jsls,ibeg,iend
      real*8  dtb,fcx,gcx,fls0,fls1,fls2,fls3,fls4
c
      integer mm,kk
      real*8  xfun,xfune
      xfun(mm) = dble(nspgd-mm)/(nspgd-2.)
      xfune(mm,kk)= dexp(-dble(mm-2)/anudg(kk))
c
c----------------------------------------------------------------------
c
      dtb=xt*60.
#ifdef MPP1
      jsls=j+myid*jxp
      jj=mjx-jsls
#else
      jsls=j
      jj=jx-jsls
#endif
      if (jj .le. ip) jsls=jj
#ifdef MPP1
      jwb=jsls
      if(jwb.gt.jxp) jwb=mod(jwb,jxp)
      if(jwb.eq.0) jwb=jxp
      if(myid.eq.nproc-1) then
        jeb=jsls
      else
        jeb=jsls+1
      endif
      if(jeb.gt.jxp) jeb=mod(jeb,jxp)
      if(jeb.eq.0) jeb=jxp
#endif
c-----determine which relaxation method to use:linear/expon.
c
      if(iboudy.eq.1)then
c
c---------use linear method
c
       if (jsls .gt. ip) then
c------interior j slices:
        do 10 i=2,ip
         ii=ilx-i+1
         fcx=fcoef*xfun(i)
         gcx=gcoef*xfun(i)
        do 10 k=1,kx
c.......south boundary:
         fls0=(tsb(i,k,j)+dtb*tsbt(i,k,j))-tb(i,k,j)
         fls1=(tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-tb(i,k,j-1)
         fls2=(tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-tb(i,k,j+1)
         fls3=(tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-tb(i-1,k,j)
         fls4=(tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-tb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(tnb(i,k,j)+dtb*tnbt(i,k,j))-tb(ii,k,j)
         fls1=(tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-tb(ii,k,j-1)
         fls2=(tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-tb(ii,k,j+1)
         fls3=(tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-tb(ii-1,k,j)
         fls4=(tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-tb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
10      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 30 i=2,jsls-1
          ii=ilx-i+1
          fcx=fcoef*xfun(i)
          gcx=gcoef*xfun(i)
         do 30 k=1,kx
c........south boundary:
          fls0=(tsb(i,k,j)+dtb*tsbt(i,k,j))-tb(i,k,j)
          fls1=(tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-tb(i,k,j-1)
          fls2=(tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-tb(i,k,j+1)
          fls3=(tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-tb(i-1,k,j)
          fls4=(tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-tb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(tnb(i,k,j)+dtb*tnbt(i,k,j))-tb(ii,k,j)
          fls1=(tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-tb(ii,k,j-1)
          fls2=(tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-tb(ii,k,j+1)
          fls3=(tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-tb(ii-1,k,j)
          fls4=(tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-tb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
30       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 50 k=1,kx
         do 50 i=ibeg,iend
#ifdef MPP1
          fls0=(twb(i,k,jwb)+dtb*twbt(i,k,jwb))-tb(i,k,j)
          fls1=(twb(i-1,k,jwb)+dtb*twbt(i-1,k,jwb))-tb(i-1,k,j)
          fls2=(twb(i+1,k,jwb)+dtb*twbt(i+1,k,jwb))-tb(i+1,k,j)
          fls3=(twb(i,k,jwb-1)+dtb*twbt(i,k,jwb-1))-tb(i,k,j-1)
          fls4=(twb(i,k,jwb+1)+dtb*twbt(i,k,jwb+1))-tb(i,k,j+1)
#else
          fls0=(twb(i,k,jsls)+dtb*twbt(i,k,jsls))-tb(i,k,j)
          fls1=(twb(i-1,k,jsls)+dtb*twbt(i-1,k,jsls))-tb(i-1,k,j)
          fls2=(twb(i+1,k,jsls)+dtb*twbt(i+1,k,jsls))-tb(i+1,k,j)
          fls3=(twb(i,k,jsls-1)+dtb*twbt(i,k,jsls-1))-tb(i,k,j-1)
          fls4=(twb(i,k,jsls+1)+dtb*twbt(i,k,jsls+1))-tb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
50       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 70 k=1,kx
         do 70 i=ibeg,iend
#ifdef MPP1
          fls0=(teb(i,k,jeb)+dtb*tebt(i,k,jeb))-tb(i,k,j)
          fls1=(teb(i-1,k,jeb)+dtb*tebt(i-1,k,jeb))-tb(i-1,k,j)
          fls2=(teb(i+1,k,jeb)+dtb*tebt(i+1,k,jeb))-tb(i+1,k,j)
          fls3=(teb(i,k,jeb-1)+dtb*tebt(i,k,jeb-1))-tb(i,k,j-1)
          fls4=(teb(i,k,jeb+1)+dtb*tebt(i,k,jeb+1))-tb(i,k,j+1)
#else
          fls0=(teb(i,k,jsls)+dtb*tebt(i,k,jsls))-tb(i,k,j)
          fls1=(teb(i-1,k,jsls)+dtb*tebt(i-1,k,jsls))-tb(i-1,k,j)
          fls2=(teb(i+1,k,jsls)+dtb*tebt(i+1,k,jsls))-tb(i+1,k,j)
          fls3=(teb(i,k,jsls-1)+dtb*tebt(i,k,jsls-1))-tb(i,k,j-1)
          fls4=(teb(i,k,jsls+1)+dtb*tebt(i,k,jsls+1))-tb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
70       continue
        endif
       endif
c
      elseif(iboudy.eq.5)then

c----------use exponential method

       if (jsls .gt. ip) then
c------interior j slices:
        do 75 i=2,ip
         ii=ilx-i+1
        do 75 k=1,kx
         fcx=fcoef*xfune(i,k)
         gcx=gcoef*xfune(i,k)
c........south boundary:
         fls0=(tsb(i,k,j)+dtb*tsbt(i,k,j))-tb(i,k,j)
         fls1=(tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-tb(i,k,j-1)
         fls2=(tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-tb(i,k,j+1)
         fls3=(tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-tb(i-1,k,j)
         fls4=(tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-tb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(tnb(i,k,j)+dtb*tnbt(i,k,j))-tb(ii,k,j)
         fls1=(tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-tb(ii,k,j-1)
         fls2=(tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-tb(ii,k,j+1)
         fls3=(tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-tb(ii-1,k,j)
         fls4=(tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-tb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
75      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 85 i=2,jsls-1
          ii=ilx-i+1
         do 85 k=1,kx
          fcx=fcoef*xfune(i,k)
          gcx=gcoef*xfune(i,k)
c.........south boundary:
          fls0=(tsb(i,k,j)+dtb*tsbt(i,k,j))-tb(i,k,j)
          fls1=(tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-tb(i,k,j-1)
          fls2=(tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-tb(i,k,j+1)
          fls3=(tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-tb(i-1,k,j)
          fls4=(tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-tb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(tnb(i,k,j)+dtb*tnbt(i,k,j))-tb(ii,k,j)
          fls1=(tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-tb(ii,k,j-1)
          fls2=(tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-tb(ii,k,j+1)
          fls3=(tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-tb(ii-1,k,j)
          fls4=(tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-tb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
85       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         do 95 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 95 i=ibeg,iend
#ifdef MPP1
          fls0=(twb(i,k,jwb)+dtb*twbt(i,k,jwb))-tb(i,k,j)
          fls1=(twb(i-1,k,jwb)+dtb*twbt(i-1,k,jwb))-tb(i-1,k,j)
          fls2=(twb(i+1,k,jwb)+dtb*twbt(i+1,k,jwb))-tb(i+1,k,j)
          fls3=(twb(i,k,jwb-1)+dtb*twbt(i,k,jwb-1))-tb(i,k,j-1)
          fls4=(twb(i,k,jwb+1)+dtb*twbt(i,k,jwb+1))-tb(i,k,j+1)
#else
          fls0=(twb(i,k,jsls)+dtb*twbt(i,k,jsls))-tb(i,k,j)
          fls1=(twb(i-1,k,jsls)+dtb*twbt(i-1,k,jsls))-tb(i-1,k,j)
          fls2=(twb(i+1,k,jsls)+dtb*twbt(i+1,k,jsls))-tb(i+1,k,j)
          fls3=(twb(i,k,jsls-1)+dtb*twbt(i,k,jsls-1))-tb(i,k,j-1)
          fls4=(twb(i,k,jsls+1)+dtb*twbt(i,k,jsls+1))-tb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
95       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         do 105 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 105 i=ibeg,iend
#ifdef MPP1
          fls0=(teb(i,k,jeb)+dtb*tebt(i,k,jeb))-tb(i,k,j)
          fls1=(teb(i-1,k,jeb)+dtb*tebt(i-1,k,jeb))-tb(i-1,k,j)
          fls2=(teb(i+1,k,jeb)+dtb*tebt(i+1,k,jeb))-tb(i+1,k,j)
          fls3=(teb(i,k,jeb-1)+dtb*tebt(i,k,jeb-1))-tb(i,k,j-1)
          fls4=(teb(i,k,jeb+1)+dtb*tebt(i,k,jeb+1))-tb(i,k,j+1)
#else
          fls0=(teb(i,k,jsls)+dtb*tebt(i,k,jsls))-tb(i,k,j)
          fls1=(teb(i-1,k,jsls)+dtb*tebt(i-1,k,jsls))-tb(i-1,k,j)
          fls2=(teb(i+1,k,jsls)+dtb*tebt(i+1,k,jsls))-tb(i+1,k,j)
          fls3=(teb(i,k,jsls-1)+dtb*tebt(i,k,jsls-1))-tb(i,k,j-1)
          fls4=(teb(i,k,jsls+1)+dtb*tebt(i,k,jsls+1))-tb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
105      continue
        endif
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nudgeQV(ip,fcoef,gcoef,tbe,xt,ften,c203,j,iboudy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies relaxation boundary conditions to the   c
c     tendency term - ften.                                           c
c                                                                     c
c     ip    : is the number of slices affected by nudging.            c
c                                                                     c
c     tbe   : is the time for the boundary values (qeb, qwb, qsb,     c
c             and qnb) in minutes.                                    c
c                                                                     c
c     xt    : is the time in minutes for variable "qvb".               c
c                                                                     c
c     fcoef : are the coefficients for the newtonian term.            c
c                                                                     c
c     gcoef : are the coefficients for the diffusion term.            c
c                                                                     c
c     ften  : is the tendency calculated from the model.              c
c                                                                     c
c     qeb, qwb, qsb, qnb : are the observed boundary values at time   c
c             "tbe" on east, west, south, and north boundaries.       c
c                                                                     c
c     qebt, qwbt, qsbt, qnbt : are the large-scale or observed        c
c             tendencies at east, west, south, and north boundaries.  c
c                                                                     c
c     tb    : is the variable at tau-1.                               c
c                                                                     c
c     c203  : = 1./(dx*dx), defined in 'param'.                       c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c     iboudy : type of boundary condition relaxation, 1=linear        c
c              5 = exponential                                        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j,iboudy
      real*8  fcoef,gcoef,tbe,xt,c203
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/main.cb'
      include 'Common2/bdycod.cb'
      integer jwb,jeb
#else
      include 'Commons/main.cb'
      include 'Commons/bdycod.cb'
#endif
      include 'Commons/param3.cb'
c
      real*8  ften(ix,kx)
      integer i,ii,k,jj,jsls,ibeg,iend
      real*8  dtb,fcx,gcx,fls0,fls1,fls2,fls3,fls4
c
      integer mm,kk
      real*8  xfun,xfune
      xfun(mm) = dble(nspgd-mm)/(nspgd-2.)
      xfune(mm,kk)= dexp(-dble(mm-2)/anudg(kk))
c
c----------------------------------------------------------------------
c
      dtb=xt*60.
#ifdef MPP1
      jsls=j+myid*jxp
      jj=mjx-jsls
#else
      jsls=j
      jj=jx-jsls
#endif
      if (jj .le. ip) jsls=jj
#ifdef MPP1
      jwb=jsls
      if(jwb.gt.jxp) jwb=mod(jwb,jxp)
      if(jwb.eq.0) jwb=jxp
      if(myid.eq.nproc-1) then
        jeb=jsls
      else
        jeb=jsls+1
      endif
      if(jeb.gt.jxp) jeb=mod(jeb,jxp)
      if(jeb.eq.0) jeb=jxp
#endif
c-----determine which relaxation method to use:linear/expon.
c
      if(iboudy.eq.1)then
c
c---------use linear method
c
       if (jsls .gt. ip) then
c------interior j slices:
        do 10 i=2,ip
         ii=ilx-i+1
         fcx=fcoef*xfun(i)
         gcx=gcoef*xfun(i)
        do 10 k=1,kx
c.......south boundary:
         fls0=(qsb(i,k,j)+dtb*qsbt(i,k,j))-qvb(i,k,j)
         fls1=(qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-qvb(i,k,j-1)
         fls2=(qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-qvb(i,k,j+1)
         fls3=(qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-qvb(i-1,k,j)
         fls4=(qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-qvb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(qnb(i,k,j)+dtb*qnbt(i,k,j))-qvb(ii,k,j)
         fls1=(qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-qvb(ii,k,j-1)
         fls2=(qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-qvb(ii,k,j+1)
         fls3=(qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-qvb(ii-1,k,j)
         fls4=(qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-qvb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
10      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 30 i=2,jsls-1
          ii=ilx-i+1
          fcx=fcoef*xfun(i)
          gcx=gcoef*xfun(i)
         do 30 k=1,kx
c........south boundary:
          fls0=(qsb(i,k,j)+dtb*qsbt(i,k,j))-qvb(i,k,j)
          fls1=(qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-qvb(i,k,j-1)
          fls2=(qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-qvb(i,k,j+1)
          fls3=(qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-qvb(i-1,k,j)
          fls4=(qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-qvb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(qnb(i,k,j)+dtb*qnbt(i,k,j))-qvb(ii,k,j)
          fls1=(qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-qvb(ii,k,j-1)
          fls2=(qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-qvb(ii,k,j+1)
          fls3=(qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-qvb(ii-1,k,j)
          fls4=(qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-qvb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
30       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 50 k=1,kx
         do 50 i=ibeg,iend
#ifdef MPP1
          fls0=(qwb(i,k,jwb)+dtb*qwbt(i,k,jwb))-qvb(i,k,j)
          fls1=(qwb(i-1,k,jwb)+dtb*qwbt(i-1,k,jwb))-qvb(i-1,k,j)
          fls2=(qwb(i+1,k,jwb)+dtb*qwbt(i+1,k,jwb))-qvb(i+1,k,j)
          fls3=(qwb(i,k,jwb-1)+dtb*qwbt(i,k,jwb-1))-qvb(i,k,j-1)
          fls4=(qwb(i,k,jwb+1)+dtb*qwbt(i,k,jwb+1))-qvb(i,k,j+1)
#else
          fls0=(qwb(i,k,jsls)+dtb*qwbt(i,k,jsls))-qvb(i,k,j)
          fls1=(qwb(i-1,k,jsls)+dtb*qwbt(i-1,k,jsls))-qvb(i-1,k,j)
          fls2=(qwb(i+1,k,jsls)+dtb*qwbt(i+1,k,jsls))-qvb(i+1,k,j)
          fls3=(qwb(i,k,jsls-1)+dtb*qwbt(i,k,jsls-1))-qvb(i,k,j-1)
          fls4=(qwb(i,k,jsls+1)+dtb*qwbt(i,k,jsls+1))-qvb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
50       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 70 k=1,kx
         do 70 i=ibeg,iend
#ifdef MPP1
          fls0=(qeb(i,k,jeb)+dtb*qebt(i,k,jeb))-qvb(i,k,j)
          fls1=(qeb(i-1,k,jeb)+dtb*qebt(i-1,k,jeb))-qvb(i-1,k,j)
          fls2=(qeb(i+1,k,jeb)+dtb*qebt(i+1,k,jeb))-qvb(i+1,k,j)
          fls3=(qeb(i,k,jeb-1)+dtb*qebt(i,k,jeb-1))-qvb(i,k,j-1)
          fls4=(qeb(i,k,jeb+1)+dtb*qebt(i,k,jeb+1))-qvb(i,k,j+1)
#else
          fls0=(qeb(i,k,jsls)+dtb*qebt(i,k,jsls))-qvb(i,k,j)
          fls1=(qeb(i-1,k,jsls)+dtb*qebt(i-1,k,jsls))-qvb(i-1,k,j)
          fls2=(qeb(i+1,k,jsls)+dtb*qebt(i+1,k,jsls))-qvb(i+1,k,j)
          fls3=(qeb(i,k,jsls-1)+dtb*qebt(i,k,jsls-1))-qvb(i,k,j-1)
          fls4=(qeb(i,k,jsls+1)+dtb*qebt(i,k,jsls+1))-qvb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
70       continue
        endif
       endif
c
      elseif(iboudy.eq.5)then

c----------use exponential method

       if (jsls .gt. ip) then
c------interior j slices:
        do 75 i=2,ip
         ii=ilx-i+1
        do 75 k=1,kx
         fcx=fcoef*xfune(i,k)
         gcx=gcoef*xfune(i,k)
c........south boundary:
         fls0=(qsb(i,k,j)+dtb*qsbt(i,k,j))-qvb(i,k,j)
         fls1=(qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-qvb(i,k,j-1)
         fls2=(qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-qvb(i,k,j+1)
         fls3=(qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-qvb(i-1,k,j)
         fls4=(qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-qvb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(qnb(i,k,j)+dtb*qnbt(i,k,j))-qvb(ii,k,j)
         fls1=(qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-qvb(ii,k,j-1)
         fls2=(qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-qvb(ii,k,j+1)
         fls3=(qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-qvb(ii-1,k,j)
         fls4=(qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-qvb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
75      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ilx-1
        if (jsls .gt. 2) then
         do 85 i=2,jsls-1
          ii=ilx-i+1
         do 85 k=1,kx
          fcx=fcoef*xfune(i,k)
          gcx=gcoef*xfune(i,k)
c.........south boundary:
          fls0=(qsb(i,k,j)+dtb*qsbt(i,k,j))-qvb(i,k,j)
          fls1=(qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-qvb(i,k,j-1)
          fls2=(qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-qvb(i,k,j+1)
          fls3=(qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-qvb(i-1,k,j)
          fls4=(qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-qvb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(qnb(i,k,j)+dtb*qnbt(i,k,j))-qvb(ii,k,j)
          fls1=(qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-qvb(ii,k,j-1)
          fls2=(qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-qvb(ii,k,j+1)
          fls3=(qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-qvb(ii-1,k,j)
          fls4=(qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-qvb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
85       continue
         ibeg = jsls
         iend = ilx-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         do 95 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 95 i=ibeg,iend
#ifdef MPP1
          fls0=(qwb(i,k,jwb)+dtb*qwbt(i,k,jwb))-qvb(i,k,j)
          fls1=(qwb(i-1,k,jwb)+dtb*qwbt(i-1,k,jwb))-qvb(i-1,k,j)
          fls2=(qwb(i+1,k,jwb)+dtb*qwbt(i+1,k,jwb))-qvb(i+1,k,j)
          fls3=(qwb(i,k,jwb-1)+dtb*qwbt(i,k,jwb-1))-qvb(i,k,j-1)
          fls4=(qwb(i,k,jwb+1)+dtb*qwbt(i,k,jwb+1))-qvb(i,k,j+1)
#else
          fls0=(qwb(i,k,jsls)+dtb*qwbt(i,k,jsls))-qvb(i,k,j)
          fls1=(qwb(i-1,k,jsls)+dtb*qwbt(i-1,k,jsls))-qvb(i-1,k,j)
          fls2=(qwb(i+1,k,jsls)+dtb*qwbt(i+1,k,jsls))-qvb(i+1,k,j)
          fls3=(qwb(i,k,jsls-1)+dtb*qwbt(i,k,jsls-1))-qvb(i,k,j-1)
          fls4=(qwb(i,k,jsls+1)+dtb*qwbt(i,k,jsls+1))-qvb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
95       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         do 105 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 105 i=ibeg,iend
#ifdef MPP1
          fls0=(qeb(i,k,jeb)+dtb*qebt(i,k,jeb))-qvb(i,k,j)
          fls1=(qeb(i-1,k,jeb)+dtb*qebt(i-1,k,jeb))-qvb(i-1,k,j)
          fls2=(qeb(i+1,k,jeb)+dtb*qebt(i+1,k,jeb))-qvb(i+1,k,j)
          fls3=(qeb(i,k,jeb-1)+dtb*qebt(i,k,jeb-1))-qvb(i,k,j-1)
          fls4=(qeb(i,k,jeb+1)+dtb*qebt(i,k,jeb+1))-qvb(i,k,j+1)
#else
          fls0=(qeb(i,k,jsls)+dtb*qebt(i,k,jsls))-qvb(i,k,j)
          fls1=(qeb(i-1,k,jsls)+dtb*qebt(i-1,k,jsls))-qvb(i-1,k,j)
          fls2=(qeb(i+1,k,jsls)+dtb*qebt(i+1,k,jsls))-qvb(i+1,k,j)
          fls3=(qeb(i,k,jsls-1)+dtb*qebt(i,k,jsls-1))-qvb(i,k,j-1)
          fls4=(qeb(i,k,jsls+1)+dtb*qebt(i,k,jsls+1))-qvb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
105      continue
        endif
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nudge_U(ip,fcoef,gcoef,tbe,xt,ften,c203,j,iboudy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies relaxation boundary conditions to the   c
c     tendency term - ften.                                           c
c                                                                     c
c     ip    : is the number of slices affected by nudging.            c
c                                                                     c
c     tbe   : is the time for the boundary values (ueb, uwb, usb,     c
c             and unb) in minutes.                                    c
c                                                                     c
c     xt    : is the time in minutes for variable "ub".               c
c                                                                     c
c     fcoef : are the coefficients for the newtonian term.            c
c                                                                     c
c     gcoef : are the coefficients for the diffusion term.            c
c                                                                     c
c     ften  : is the tendency calculated from the model.              c
c                                                                     c
c     ueb, uwb, usb, unb : are the observed boundary values at time   c
c             "tbe" on east, west, south, and north boundaries.       c
c                                                                     c
c     uebt, uwbt, usbt, unbt : are the large-scale or observed        c
c             tendencies at east, west, south, and north boundaries.  c
c                                                                     c
c     tb    : is the variable at tau-1.                               c
c                                                                     c
c     c203  : = 1./(dx*dx), defined in 'param'.                       c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c     iboudy : type of boundary condition relaxation, 1=linear        c
c              5 = exponential                                        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j,iboudy
      real*8  fcoef,gcoef,tbe,xt,c203
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/main.cb'
      include 'Common2/bdycod.cb'
      integer jew
#else
      include 'Commons/main.cb'
      include 'Commons/bdycod.cb'
#endif
      include 'Commons/param3.cb'
c
      real*8  ften(ix,kx)
      integer i,ii,k,jj,jsls,ibeg,iend
      real*8  dtb,fcx,gcx,fls0,fls1,fls2,fls3,fls4
c
      integer mm,kk
      real*8  xfun,xfune
      xfun(mm) = dble(nspgd-mm)/(nspgd-2.)
      xfune(mm,kk)= dexp(-dble(mm-2)/anudg(kk))
c
c----------------------------------------------------------------------
c
      dtb=xt*60.
#ifdef MPP1
      jsls=j+myid*jxp
      jj=mjx+1-jsls
#else
      jsls=j
      jj=jx+1-jsls
#endif
      if (jj .le. ip) jsls=jj
#ifdef MPP1
      jew=jsls
      if(jew.gt.jxp) jew=mod(jsls,jxp)
      if(jew.eq.0) jew=jxp
#endif
c-----determine which relaxation method to use:linear/expon.
c
      if(iboudy.eq.1)then
c
c---------use linear method
c
       if (jsls .gt. ip) then
c------interior j slices:
        do 10 i=2,ip
         ii=ix-i+1
         fcx=fcoef*xfun(i)
         gcx=gcoef*xfun(i)
        do 10 k=1,kx
c.......south boundary:
         fls0=(usb(i,k,j)+dtb*usbt(i,k,j))-ub(i,k,j)
         fls1=(usb(i,k,j-1)+dtb*usbt(i,k,j-1))-ub(i,k,j-1)
         fls2=(usb(i,k,j+1)+dtb*usbt(i,k,j+1))-ub(i,k,j+1)
         fls3=(usb(i-1,k,j)+dtb*usbt(i-1,k,j))-ub(i-1,k,j)
         fls4=(usb(i+1,k,j)+dtb*usbt(i+1,k,j))-ub(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(unb(i,k,j)+dtb*unbt(i,k,j))-ub(ii,k,j)
         fls1=(unb(i,k,j-1)+dtb*unbt(i,k,j-1))-ub(ii,k,j-1)
         fls2=(unb(i,k,j+1)+dtb*unbt(i,k,j+1))-ub(ii,k,j+1)
         fls3=(unb(i-1,k,j)+dtb*unbt(i-1,k,j))-ub(ii-1,k,j)
         fls4=(unb(i+1,k,j)+dtb*unbt(i+1,k,j))-ub(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
10      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ix-1
        if (jsls .gt. 2) then
         do 30 i=2,jsls-1
          ii=ix-i+1
          fcx=fcoef*xfun(i)
          gcx=gcoef*xfun(i)
         do 30 k=1,kx
c........south boundary:
          fls0=(usb(i,k,j)+dtb*usbt(i,k,j))-ub(i,k,j)
          fls1=(usb(i,k,j-1)+dtb*usbt(i,k,j-1))-ub(i,k,j-1)
          fls2=(usb(i,k,j+1)+dtb*usbt(i,k,j+1))-ub(i,k,j+1)
          fls3=(usb(i-1,k,j)+dtb*usbt(i-1,k,j))-ub(i-1,k,j)
          fls4=(usb(i+1,k,j)+dtb*usbt(i+1,k,j))-ub(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(unb(i,k,j)+dtb*unbt(i,k,j))-ub(ii,k,j)
          fls1=(unb(i,k,j-1)+dtb*unbt(i,k,j-1))-ub(ii,k,j-1)
          fls2=(unb(i,k,j+1)+dtb*unbt(i,k,j+1))-ub(ii,k,j+1)
          fls3=(unb(i-1,k,j)+dtb*unbt(i-1,k,j))-ub(ii-1,k,j)
          fls4=(unb(i+1,k,j)+dtb*unbt(i+1,k,j))-ub(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
30       continue
         ibeg = jsls
         iend = ix-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 50 k=1,kx
         do 50 i=ibeg,iend
#ifdef MPP1
          fls0=(uwb(i,k,jew)+dtb*uwbt(i,k,jew))-ub(i,k,j)
          fls1=(uwb(i-1,k,jew)+dtb*uwbt(i-1,k,jew))-ub(i-1,k,j)
          fls2=(uwb(i+1,k,jew)+dtb*uwbt(i+1,k,jew))-ub(i+1,k,j)
          fls3=(uwb(i,k,jew-1)+dtb*uwbt(i,k,jew-1))-ub(i,k,j-1)
          fls4=(uwb(i,k,jew+1)+dtb*uwbt(i,k,jew+1))-ub(i,k,j+1)
#else
          fls0=(uwb(i,k,jsls)+dtb*uwbt(i,k,jsls))-ub(i,k,j)
          fls1=(uwb(i-1,k,jsls)+dtb*uwbt(i-1,k,jsls))-ub(i-1,k,j)
          fls2=(uwb(i+1,k,jsls)+dtb*uwbt(i+1,k,jsls))-ub(i+1,k,j)
          fls3=(uwb(i,k,jsls-1)+dtb*uwbt(i,k,jsls-1))-ub(i,k,j-1)
          fls4=(uwb(i,k,jsls+1)+dtb*uwbt(i,k,jsls+1))-ub(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
50       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 70 k=1,kx
         do 70 i=ibeg,iend
#ifdef MPP1
          fls0=(ueb(i,k,jew)+dtb*uebt(i,k,jew))-ub(i,k,j)
          fls1=(ueb(i-1,k,jew)+dtb*uebt(i-1,k,jew))-ub(i-1,k,j)
          fls2=(ueb(i+1,k,jew)+dtb*uebt(i+1,k,jew))-ub(i+1,k,j)
          fls3=(ueb(i,k,jew-1)+dtb*uebt(i,k,jew-1))-ub(i,k,j-1)
          fls4=(ueb(i,k,jew+1)+dtb*uebt(i,k,jew+1))-ub(i,k,j+1)
#else
          fls0=(ueb(i,k,jsls)+dtb*uebt(i,k,jsls))-ub(i,k,j)
          fls1=(ueb(i-1,k,jsls)+dtb*uebt(i-1,k,jsls))-ub(i-1,k,j)
          fls2=(ueb(i+1,k,jsls)+dtb*uebt(i+1,k,jsls))-ub(i+1,k,j)
          fls3=(ueb(i,k,jsls-1)+dtb*uebt(i,k,jsls-1))-ub(i,k,j-1)
          fls4=(ueb(i,k,jsls+1)+dtb*uebt(i,k,jsls+1))-ub(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
70       continue
        endif
       endif
c
      elseif(iboudy.eq.5)then

c----------use exponential method

       if (jsls .gt. ip) then
c------interior j slices:
        do 75 i=2,ip
         ii=ix-i+1
        do 75 k=1,kx
         fcx=fcoef*xfune(i,k)
         gcx=gcoef*xfune(i,k)
c........south boundary:
         fls0=(usb(i,k,j)+dtb*usbt(i,k,j))-ub(i,k,j)
         fls1=(usb(i,k,j-1)+dtb*usbt(i,k,j-1))-ub(i,k,j-1)
         fls2=(usb(i,k,j+1)+dtb*usbt(i,k,j+1))-ub(i,k,j+1)
         fls3=(usb(i-1,k,j)+dtb*usbt(i-1,k,j))-ub(i-1,k,j)
         fls4=(usb(i+1,k,j)+dtb*usbt(i+1,k,j))-ub(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(unb(i,k,j)+dtb*unbt(i,k,j))-ub(ii,k,j)
         fls1=(unb(i,k,j-1)+dtb*unbt(i,k,j-1))-ub(ii,k,j-1)
         fls2=(unb(i,k,j+1)+dtb*unbt(i,k,j+1))-ub(ii,k,j+1)
         fls3=(unb(i-1,k,j)+dtb*unbt(i-1,k,j))-ub(ii-1,k,j)
         fls4=(unb(i+1,k,j)+dtb*unbt(i+1,k,j))-ub(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
75      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ix-1
        if (jsls .gt. 2) then
         do 85 i=2,jsls-1
          ii=ix-i+1
         do 85 k=1,kx
          fcx=fcoef*xfune(i,k)
          gcx=gcoef*xfune(i,k)
c.........south boundary:
          fls0=(usb(i,k,j)+dtb*usbt(i,k,j))-ub(i,k,j)
          fls1=(usb(i,k,j-1)+dtb*usbt(i,k,j-1))-ub(i,k,j-1)
          fls2=(usb(i,k,j+1)+dtb*usbt(i,k,j+1))-ub(i,k,j+1)
          fls3=(usb(i-1,k,j)+dtb*usbt(i-1,k,j))-ub(i-1,k,j)
          fls4=(usb(i+1,k,j)+dtb*usbt(i+1,k,j))-ub(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(unb(i,k,j)+dtb*unbt(i,k,j))-ub(ii,k,j)
          fls1=(unb(i,k,j-1)+dtb*unbt(i,k,j-1))-ub(ii,k,j-1)
          fls2=(unb(i,k,j+1)+dtb*unbt(i,k,j+1))-ub(ii,k,j+1)
          fls3=(unb(i-1,k,j)+dtb*unbt(i-1,k,j))-ub(ii-1,k,j)
          fls4=(unb(i+1,k,j)+dtb*unbt(i+1,k,j))-ub(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
85       continue
         ibeg = jsls
         iend = ix-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         do 95 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 95 i=ibeg,iend
#ifdef MPP1
          fls0=(uwb(i,k,jew)+dtb*uwbt(i,k,jew))-ub(i,k,j)
          fls1=(uwb(i-1,k,jew)+dtb*uwbt(i-1,k,jew))-ub(i-1,k,j)
          fls2=(uwb(i+1,k,jew)+dtb*uwbt(i+1,k,jew))-ub(i+1,k,j)
          fls3=(uwb(i,k,jew-1)+dtb*uwbt(i,k,jew-1))-ub(i,k,j-1)
          fls4=(uwb(i,k,jew+1)+dtb*uwbt(i,k,jew+1))-ub(i,k,j+1)
#else
          fls0=(uwb(i,k,jsls)+dtb*uwbt(i,k,jsls))-ub(i,k,j)
          fls1=(uwb(i-1,k,jsls)+dtb*uwbt(i-1,k,jsls))-ub(i-1,k,j)
          fls2=(uwb(i+1,k,jsls)+dtb*uwbt(i+1,k,jsls))-ub(i+1,k,j)
          fls3=(uwb(i,k,jsls-1)+dtb*uwbt(i,k,jsls-1))-ub(i,k,j-1)
          fls4=(uwb(i,k,jsls+1)+dtb*uwbt(i,k,jsls+1))-ub(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
95       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         do 105 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 105 i=ibeg,iend
#ifdef MPP1
          fls0=(ueb(i,k,jew)+dtb*uebt(i,k,jew))-ub(i,k,j)
          fls1=(ueb(i-1,k,jew)+dtb*uebt(i-1,k,jew))-ub(i-1,k,j)
          fls2=(ueb(i+1,k,jew)+dtb*uebt(i+1,k,jew))-ub(i+1,k,j)
          fls3=(ueb(i,k,jew-1)+dtb*uebt(i,k,jew-1))-ub(i,k,j-1)
          fls4=(ueb(i,k,jew+1)+dtb*uebt(i,k,jew+1))-ub(i,k,j+1)
#else
          fls0=(ueb(i,k,jsls)+dtb*uebt(i,k,jsls))-ub(i,k,j)
          fls1=(ueb(i-1,k,jsls)+dtb*uebt(i-1,k,jsls))-ub(i-1,k,j)
          fls2=(ueb(i+1,k,jsls)+dtb*uebt(i+1,k,jsls))-ub(i+1,k,j)
          fls3=(ueb(i,k,jsls-1)+dtb*uebt(i,k,jsls-1))-ub(i,k,j-1)
          fls4=(ueb(i,k,jsls+1)+dtb*uebt(i,k,jsls+1))-ub(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
105      continue
        endif
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nudge_V(ip,fcoef,gcoef,tbe,xt,ften,c203,j,iboudy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies relaxation boundary conditions to the   c
c     tendency term - ften.                                           c
c                                                                     c
c     ip    : is the number of slices affected by nudging.            c
c                                                                     c
c     tbe   : is the time for the boundary values (veb, vwb, vsb,     c
c             and vnb) in minutes.                                    c
c                                                                     c
c     xt    : is the time in minutes for variable "vb".               c
c                                                                     c
c     fcoef : are the coefficients for the newtonian term.            c
c                                                                     c
c     gcoef : are the coefficients for the diffusion term.            c
c                                                                     c
c     ften  : is the tendency calculated from the model.              c
c                                                                     c
c     veb, vwb, vsb, vnb : are the observed boundary values at time   c
c             "tbe" on east, west, south, and north boundaries.       c
c                                                                     c
c     vebt, vwbt, vsbt, vnbt : are the large-scale or observed        c
c             tendencies at east, west, south, and north boundaries.  c
c                                                                     c
c     tb    : is the variable at tau-1.                               c
c                                                                     c
c     c203  : = 1./(dx*dx), defined in 'param'.                       c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c     iboudy : type of boundary condition relaxation, 1=linear        c
c              5 = exponential                                        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j,iboudy
      real*8  fcoef,gcoef,tbe,xt,c203
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/bdycod.cb'
      include 'Common2/main.cb'
      integer jew
#else
      include 'Commons/main.cb'
      include 'Commons/bdycod.cb'
#endif
      include 'Commons/param3.cb'
c
      real*8  ften(ix,kx)
      integer i,ii,k,jj,jsls,ibeg,iend
      real*8  dtb,fcx,gcx,fls0,fls1,fls2,fls3,fls4
c
      integer mm,kk
      real*8  xfun,xfune
      xfun(mm) = dble(nspgd-mm)/(nspgd-2.)
      xfune(mm,kk)= dexp(-dble(mm-2)/anudg(kk))
c
c----------------------------------------------------------------------
c
      dtb=xt*60.
#ifdef MPP1
      jsls=j+myid*jxp
      jj=mjx+1-jsls
#else
      jsls=j
      jj=jx+1-jsls
#endif
      if (jj .le. ip) jsls=jj
#ifdef MPP1
      jew=jsls
      if(jew.gt.jxp) jew=mod(jsls,jxp)
      if(jew.eq.0) jew=jxp
#endif
c-----determine which relaxation method to use:linear/expon.
c
      if(iboudy.eq.1)then
c
c---------use linear method
c
       if (jsls .gt. ip) then
c------interior j slices:
        do 10 i=2,ip
         ii=ix-i+1
         fcx=fcoef*xfun(i)
         gcx=gcoef*xfun(i)
        do 10 k=1,kx
c.......south boundary:
         fls0=(vsb(i,k,j)+dtb*vsbt(i,k,j))-vb(i,k,j)
         fls1=(vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-vb(i,k,j-1)
         fls2=(vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-vb(i,k,j+1)
         fls3=(vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-vb(i-1,k,j)
         fls4=(vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-vb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(vnb(i,k,j)+dtb*vnbt(i,k,j))-vb(ii,k,j)
         fls1=(vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-vb(ii,k,j-1)
         fls2=(vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-vb(ii,k,j+1)
         fls3=(vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-vb(ii-1,k,j)
         fls4=(vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-vb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
10      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ix-1
        if (jsls .gt. 2) then
         do 30 i=2,jsls-1
          ii=ix-i+1
          fcx=fcoef*xfun(i)
          gcx=gcoef*xfun(i)
         do 30 k=1,kx
c........south boundary:
          fls0=(vsb(i,k,j)+dtb*vsbt(i,k,j))-vb(i,k,j)
          fls1=(vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-vb(i,k,j-1)
          fls2=(vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-vb(i,k,j+1)
          fls3=(vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-vb(i-1,k,j)
          fls4=(vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-vb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(vnb(i,k,j)+dtb*vnbt(i,k,j))-vb(ii,k,j)
          fls1=(vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-vb(ii,k,j-1)
          fls2=(vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-vb(ii,k,j+1)
          fls3=(vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-vb(ii-1,k,j)
          fls4=(vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-vb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
30       continue
         ibeg = jsls
         iend = ix-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 50 k=1,kx
         do 50 i=ibeg,iend
#ifdef MPP1
          fls0=(vwb(i,k,jew)+dtb*vwbt(i,k,jew))-vb(i,k,j)
          fls1=(vwb(i-1,k,jew)+dtb*vwbt(i-1,k,jew))-vb(i-1,k,j)
          fls2=(vwb(i+1,k,jew)+dtb*vwbt(i+1,k,jew))-vb(i+1,k,j)
          fls3=(vwb(i,k,jew-1)+dtb*vwbt(i,k,jew-1))-vb(i,k,j-1)
          fls4=(vwb(i,k,jew+1)+dtb*vwbt(i,k,jew+1))-vb(i,k,j+1)
#else
          fls0=(vwb(i,k,jsls)+dtb*vwbt(i,k,jsls))-vb(i,k,j)
          fls1=(vwb(i-1,k,jsls)+dtb*vwbt(i-1,k,jsls))-vb(i-1,k,j)
          fls2=(vwb(i+1,k,jsls)+dtb*vwbt(i+1,k,jsls))-vb(i+1,k,j)
          fls3=(vwb(i,k,jsls-1)+dtb*vwbt(i,k,jsls-1))-vb(i,k,j-1)
          fls4=(vwb(i,k,jsls+1)+dtb*vwbt(i,k,jsls+1))-vb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
50       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         fcx=fcoef*xfun(jsls)
         gcx=gcoef*xfun(jsls)
         do 70 k=1,kx
         do 70 i=ibeg,iend
#ifdef MPP1
          fls0=(veb(i,k,jew)+dtb*vebt(i,k,jew))-vb(i,k,j)
          fls1=(veb(i-1,k,jew)+dtb*vebt(i-1,k,jew))-vb(i-1,k,j)
          fls2=(veb(i+1,k,jew)+dtb*vebt(i+1,k,jew))-vb(i+1,k,j)
          fls3=(veb(i,k,jew-1)+dtb*vebt(i,k,jew-1))-vb(i,k,j-1)
          fls4=(veb(i,k,jew+1)+dtb*vebt(i,k,jew+1))-vb(i,k,j+1)
#else
          fls0=(veb(i,k,jsls)+dtb*vebt(i,k,jsls))-vb(i,k,j)
          fls1=(veb(i-1,k,jsls)+dtb*vebt(i-1,k,jsls))-vb(i-1,k,j)
          fls2=(veb(i+1,k,jsls)+dtb*vebt(i+1,k,jsls))-vb(i+1,k,j)
          fls3=(veb(i,k,jsls-1)+dtb*vebt(i,k,jsls-1))-vb(i,k,j-1)
          fls4=(veb(i,k,jsls+1)+dtb*vebt(i,k,jsls+1))-vb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
70       continue
        endif
       endif
c
      elseif(iboudy.eq.5)then

c----------use exponential method

       if (jsls .gt. ip) then
c------interior j slices:
        do 75 i=2,ip
         ii=ix-i+1
        do 75 k=1,kx
         fcx=fcoef*xfune(i,k)
         gcx=gcoef*xfune(i,k)
c........south boundary:
         fls0=(vsb(i,k,j)+dtb*vsbt(i,k,j))-vb(i,k,j)
         fls1=(vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-vb(i,k,j-1)
         fls2=(vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-vb(i,k,j+1)
         fls3=(vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-vb(i-1,k,j)
         fls4=(vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-vb(i+1,k,j)
         ften(i,k)=ften(i,k)+fcx*fls0
     &                      -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c........north boundary:
         fls0=(vnb(i,k,j)+dtb*vnbt(i,k,j))-vb(ii,k,j)
         fls1=(vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-vb(ii,k,j-1)
         fls2=(vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-vb(ii,k,j+1)
         fls3=(vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-vb(ii-1,k,j)
         fls4=(vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-vb(ii+1,k,j)
         ften(ii,k)=ften(ii,k)+fcx*fls0
     &                        -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
75      continue
c
       elseif(jsls.le.ip)then
c------east or west boundary slices:
        ibeg = 2
        iend = ix-1
        if (jsls .gt. 2) then
         do 85 i=2,jsls-1
          ii=ix-i+1
         do 85 k=1,kx
          fcx=fcoef*xfune(i,k)
          gcx=gcoef*xfune(i,k)
c.........south boundary:
          fls0=(vsb(i,k,j)+dtb*vsbt(i,k,j))-vb(i,k,j)
          fls1=(vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-vb(i,k,j-1)
          fls2=(vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-vb(i,k,j+1)
          fls3=(vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-vb(i-1,k,j)
          fls4=(vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-vb(i+1,k,j)
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
c.........north boundary:
          fls0=(vnb(i,k,j)+dtb*vnbt(i,k,j))-vb(ii,k,j)
          fls1=(vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-vb(ii,k,j-1)
          fls2=(vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-vb(ii,k,j+1)
          fls3=(vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-vb(ii-1,k,j)
          fls4=(vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-vb(ii+1,k,j)
          ften(ii,k)=ften(ii,k)+fcx*fls0
     &                         -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
85       continue
         ibeg = jsls
         iend = ix-jsls+1
        endif
c
        if (jj .gt. ip) then
c-------west-boundary slice:
         do 95 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 95 i=ibeg,iend
#ifdef MPP1
          fls0=(vwb(i,k,jew)+dtb*vwbt(i,k,jew))-vb(i,k,j)
          fls1=(vwb(i-1,k,jew)+dtb*vwbt(i-1,k,jew))-vb(i-1,k,j)
          fls2=(vwb(i+1,k,jew)+dtb*vwbt(i+1,k,jew))-vb(i+1,k,j)
          fls3=(vwb(i,k,jew-1)+dtb*vwbt(i,k,jew-1))-vb(i,k,j-1)
          fls4=(vwb(i,k,jew+1)+dtb*vwbt(i,k,jew+1))-vb(i,k,j+1)
#else
          fls0=(vwb(i,k,jsls)+dtb*vwbt(i,k,jsls))-vb(i,k,j)
          fls1=(vwb(i-1,k,jsls)+dtb*vwbt(i-1,k,jsls))-vb(i-1,k,j)
          fls2=(vwb(i+1,k,jsls)+dtb*vwbt(i+1,k,jsls))-vb(i+1,k,j)
          fls3=(vwb(i,k,jsls-1)+dtb*vwbt(i,k,jsls-1))-vb(i,k,j-1)
          fls4=(vwb(i,k,jsls+1)+dtb*vwbt(i,k,jsls+1))-vb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
95       continue
        elseif(jj.le.ip)then
c-------east-boundary slice:
         do 105 k=1,kx
          fcx=fcoef*xfune(jsls,k)
          gcx=gcoef*xfune(jsls,k)
         do 105 i=ibeg,iend
#ifdef MPP1
          fls0=(veb(i,k,jew)+dtb*vebt(i,k,jew))-vb(i,k,j)
          fls1=(veb(i-1,k,jew)+dtb*vebt(i-1,k,jew))-vb(i-1,k,j)
          fls2=(veb(i+1,k,jew)+dtb*vebt(i+1,k,jew))-vb(i+1,k,j)
          fls3=(veb(i,k,jew-1)+dtb*vebt(i,k,jew-1))-vb(i,k,j-1)
          fls4=(veb(i,k,jew+1)+dtb*vebt(i,k,jew+1))-vb(i,k,j+1)
#else
          fls0=(veb(i,k,jsls)+dtb*vebt(i,k,jsls))-vb(i,k,j)
          fls1=(veb(i-1,k,jsls)+dtb*vebt(i-1,k,jsls))-vb(i-1,k,j)
          fls2=(veb(i+1,k,jsls)+dtb*vebt(i+1,k,jsls))-vb(i+1,k,j)
          fls3=(veb(i,k,jsls-1)+dtb*vebt(i,k,jsls-1))-vb(i,k,j-1)
          fls4=(veb(i,k,jsls+1)+dtb*vebt(i,k,jsls+1))-vb(i,k,j+1)
#endif
          ften(i,k)=ften(i,k)+fcx*fls0
     &                       -gcx*c203*(fls1+fls2+fls3+fls4-4.*fls0)
105      continue
        endif
       endif
c
      endif
      return
      end
