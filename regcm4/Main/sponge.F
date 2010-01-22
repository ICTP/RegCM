C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine sponge_P(ip,wg,ften,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies sponge boundary condition to the        c
c     tendency term - ften.                                           c
c                                                                     c
c     ip   : is the number of slices affected by sponge boundary.     c
c                                                                     c
c     wg   : are the weightings.                                      c
c                                                                     c
c     ften : is the tendency calculated from the model.               c
c                                                                     c
c     pebt, pwbt, pnbt, psbt : are the large-scale or observed        c
c            tendencies at east, west, north, and south boundaries.   c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bdycod.cb"
#ifdef MPP1
      include 'mpiregcm.h'
      integer jwb,jeb
#endif
c
      real(kind=8)  ften(ix),      wg(ip)
      integer i,ii,jj,jsls,ibeg,iend
c
c----------------------------------------------------------------------
c
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
c
      if (jsls .gt. ip) then
c-----interior j slices:
       do 10 i=2,ip
        ii=ix-i
c.......south boundary:
        ften(i)=wg(i)*ften(i)+(1.-wg(i))*psbt(i,j)
c.......north boundary:
        ften(ii)=wg(i)*ften(ii)+(1.-wg(i))*pnbt(i,j)
10     continue
c
      elseif(jsls.le.ip)then
       ibeg = 2
       iend = ilx-1
       if (jsls .gt. 2) then
        do 30 i=2,jsls-1
         ii=ix-i
c........south boundary:
         ften(i)=wg(i)*ften(i)+(1.-wg(i))*psbt(i,j)
c........north boundary:
         ften(ii)=wg(i)*ften(ii)+(1.-wg(i))*pnbt(i,j)
30      continue
        ibeg = jsls
        iend = ix-jsls
       endif
c
       if (jj .gt. ip) then
c------west-boundary slice:
        do 50 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i)=wg(jsls)*ften(i)+(1.-wg(jsls))*pwbt(i,jwb)
#else
         ften(i)=wg(jsls)*ften(i)+(1.-wg(jsls))*pwbt(i,jsls)
#endif
50      continue
       elseif(jj.le.ip)then
c------east-boundary slice:
        do 70 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i)=wg(jsls)*ften(i)+(1.-wg(jsls))*pebt(i,jeb)
#else
         ften(i)=wg(jsls)*ften(i)+(1.-wg(jsls))*pebt(i,jsls)
#endif
70      continue
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sponge_T(ip,wg,ften,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies sponge boundary condition to the        c
c     tendency term - ften.                                           c
c                                                                     c
c     ip   : is the number of slices affected by sponge boundary.     c
c                                                                     c
c     wg   : are the weightings.                                      c
c                                                                     c
c     ften : is the tendency calculated from the model.               c
c                                                                     c
c     tebt, twbt, tnbt, tsbt : are the large-scale or observed        c
c            tendencies at east, west, north, and south boundaries.   c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bdycod.cb"
#ifdef MPP1
      include 'mpiregcm.h'
      integer jwb,jeb
#endif
c
      real(kind=8)  ften(ix,kx),   wg(ip)
      integer i,ii,jj,k,jsls,ibeg,iend
c
c----------------------------------------------------------------------
c
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
c
      if (jsls .gt. ip) then
c-----interior j slices:
       do 10 i=2,ip
        ii=ix-i
       do 10 k=1,kx
c.......south boundary:
        ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*tsbt(i,k,j)
c.......north boundary:
        ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*tnbt(i,k,j)
10     continue
c
      elseif(jsls.le.ip)then
       ibeg = 2
       iend = ilx-1
       if (jsls .gt. 2) then
        do 30 i=2,jsls-1
         ii=ix-i
        do 30 k=1,kx
c........south boundary:
         ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*tsbt(i,k,j)
c........north boundary:
         ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*tnbt(i,k,j)
30      continue
        ibeg = jsls
        iend = ix-jsls
       endif
c
       if (jj .gt. ip) then
c------west-boundary slice:
        do 50 k=1,kx
        do 50 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*twbt(i,k,jwb)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*twbt(i,k,jsls)
#endif
50      continue
       elseif(jj.le.ip)then
c------east-boundary slice:
        do 70 k=1,kx
        do 70 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*tebt(i,k,jeb)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*tebt(i,k,jsls)
#endif
70      continue
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spongeQV(ip,wg,ften,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies sponge boundary condition to the        c
c     tendency term - ften.                                           c
c                                                                     c
c     ip   : is the number of slices affected by sponge boundary.     c
c                                                                     c
c     wg   : are the weightings.                                      c
c                                                                     c
c     ften : is the tendency calculated from the model.               c
c                                                                     c
c     qebt, qwbt, qnbt, qsbt : are the large-scale or observed        c
c            tendencies at east, west, north, and south boundaries.   c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bdycod.cb"
#ifdef MPP1
      include 'mpiregcm.h'
      integer jwb,jeb
#endif
c
      real(kind=8)  ften(ix,kx),   wg(ip)
      integer i,ii,jj,k,jsls,ibeg,iend
c
c----------------------------------------------------------------------
c
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
c
      if (jsls .gt. ip) then
c-----interior j slices:
       do 10 i=2,ip
        ii=ix-i
       do 10 k=1,kx
c.......south boundary:
        ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*qsbt(i,k,j)
c.......north boundary:
        ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*qnbt(i,k,j)
10     continue
c
      elseif(jsls.le.ip)then
       ibeg = 2
       iend = ilx-1
       if (jsls .gt. 2) then
        do 30 i=2,jsls-1
         ii=ix-i
        do 30 k=1,kx
c........south boundary:
         ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*qsbt(i,k,j)
c........north boundary:
         ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*qnbt(i,k,j)
30      continue
        ibeg = jsls
        iend = ix-jsls
       endif
c
       if (jj .gt. ip) then
c------west-boundary slice:
        do 50 k=1,kx
        do 50 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip) 
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*qwbt(i,k,jwb)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*qwbt(i,k,jsls)
#endif
50      continue
       elseif(jj.le.ip)then
c------east-boundary slice:
        do 70 k=1,kx
        do 70 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*qebt(i,k,jeb)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*qebt(i,k,jsls)
#endif
70      continue
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sponge_U(ip,wg,ften,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies sponge boundary condition to the        c
c     tendency term - ften.                                           c
c                                                                     c
c     ip   : is the number of slices affected by sponge boundary.     c
c                                                                     c
c     wg   : are the weightings.                                      c
c                                                                     c
c     ften : is the tendency calculated from the model.               c
c                                                                     c
c     uebt, uwbt, unbt, usbt : are the large-scale or observed        c
c            tendencies at east, west, north, and south boundaries.   c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bdycod.cb"
#ifdef MPP1
      include 'mpiregcm.h'
      integer jew
#endif
c
      real(kind=8)  ften(ix,kx),   wg(ip)
      integer i,ii,jj,k,jsls,ibeg,iend
c
c----------------------------------------------------------------------
c
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
c
      if (jsls .gt. ip) then
c-----interior j slices:
       do 10 i=2,ip
        ii=ix-i+1
       do 10 k=1,kx
c.......south boundary:
        ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*usbt(i,k,j)
c.......north boundary:
        ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*unbt(i,k,j)
10     continue
c
      elseif(jsls.le.ip)then
       ibeg = 2
       iend = ix-1
       if (jsls .gt. 2) then
        do 30 i=2,jsls-1
         ii=ix-i+1
        do 30 k=1,kx
c........south boundary:
         ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*usbt(i,k,j)
c........north boundary:
         ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*unbt(i,k,j)
30      continue
        ibeg = jsls
        iend = ix-jsls+1
       endif
c
       if (jj .gt. ip) then
c------west-boundary slice:
        do 50 k=1,kx
        do 50 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*uwbt(i,k,jew)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*uwbt(i,k,jsls)
#endif
50      continue
       elseif(jj.le.ip)then
c------east-boundary slice:
        do 70 k=1,kx
        do 70 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*uebt(i,k,jew)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*uebt(i,k,jsls)
#endif
70      continue
       endif
c
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sponge_V(ip,wg,ften,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine applies sponge boundary condition to the        c
c     tendency term - ften.                                           c
c                                                                     c
c     ip   : is the number of slices affected by sponge boundary.     c
c                                                                     c
c     wg   : are the weightings.                                      c
c                                                                     c
c     ften : is the tendency calculated from the model.               c
c                                                                     c
c     vebt, vwbt, vnbt, vsbt : are the large-scale or observed        c
c            tendencies at east, west, north, and south boundaries.   c
c                                                                     c
c     ie = ix, je = jx for dot-point variables.                       c
c     ie = ix-1, je = jx-1 for cross-point variables.                 c
c                                                                     c
c     j    : is the j'th slice of the tendency to be adjusted.        c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ip,j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bdycod.cb"
#ifdef MPP1
      include 'mpiregcm.h'
      integer jew
#endif
c
      real(kind=8)  ften(ix,kx),   wg(ip)
      integer i,ii,jj,k,jsls,ibeg,iend
c
c----------------------------------------------------------------------
c
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
c
      if (jsls .gt. ip) then
c-----interior j slices:
       do 10 i=2,ip
        ii=ix-i+1
       do 10 k=1,kx
c.......south boundary:
        ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*vsbt(i,k,j)
c.......north boundary:
        ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*vnbt(i,k,j)
10     continue
c
      elseif(jsls.le.ip)then
       ibeg = 2
       iend = ix-1
       if (jsls .gt. 2) then
        do 30 i=2,jsls-1
         ii=ix-i+1
        do 30 k=1,kx
c........south boundary:
         ften(i,k)=wg(i)*ften(i,k)+(1.-wg(i))*vsbt(i,k,j)
c........north boundary:
         ften(ii,k)=wg(i)*ften(ii,k)+(1.-wg(i))*vnbt(i,k,j)
30      continue
        ibeg = jsls
        iend = ix-jsls+1
       endif
c
       if (jj .gt. ip) then
c------west-boundary slice:
        do 50 k=1,kx
        do 50 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*vwbt(i,k,jew)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*vwbt(i,k,jsls)
#endif
50      continue
       elseif(jj.le.ip)then
c------east-boundary slice:
        do 70 k=1,kx
        do 70 i=ibeg,iend
#ifdef MPP1
         if(jsls.le.ip)
     &   ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*vebt(i,k,jew)
#else
         ften(i,k)=wg(jsls)*ften(i,k)+(1.-wg(jsls))*vebt(i,k,jsls)
#endif
70      continue
       endif
c
      endif
      return
      end
