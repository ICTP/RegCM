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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadv_T(ften,dxx,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfx   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/cvaria.cb"
#     include "include/main.cb"
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)-((ua(i+1,k,j+1)+ua(i,k,j+1))*(t(i,k,j+1)+
     1          t(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))*(t(i,k,j)+
     2          t(i,k,j-1))+(va(i+1,k,j+1)+va(i+1,k,j))*(t(i+1,k,j)+
     3          t(i,k,j))-(va(i,k,j+1)+va(i,k,j))*(t(i-1,k,j)+
     4          t(i,k,j)))/(dxx*msfx(i,j)*msfx(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*t(i,k,j)+fact2*t(i,k,j+1)
        else
         fx2=fact1*t(i,k,j+1)+fact2*t(i,k,j)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*t(i,k,j-1)+fact2*t(i,k,j)
        else
         fx1=fact1*t(i,k,j)+fact2*t(i,k,j-1)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*t(i,k,j)+fact2*t(i+1,k,j)
        else
         fy2=fact1*t(i+1,k,j)+fact2*t(i,k,j)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*t(i-1,k,j)+fact2*t(i,k,j)
        else
         fy1=fact1*t(i,k,j)+fact2*t(i-1,k,j)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfx(i,j)*msfx(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((t(i,k,j+1)+t(i,k,j))*ucmonb-
     &          (t(i,k,j)+t(i,k,j-1))*ucmonc+(t(i+1,k,j)+
     &          t(i,k,j))*vcmonb-(t(i,k,j)+t(i-1,k,j))*
     &          vcmonc)/(dxx*msfx(i,j)*msfx(i,j))
60     continue
c
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadv_U(ften,dxx,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfd   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#     include "include/cvaria.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)-((ua(i+1,k,j+1)+ua(i,k,j+1))*(u(i,k,j+1)+
     1          u(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))*(u(i,k,j)+
     2          u(i,k,j-1))+(va(i+1,k,j+1)+va(i+1,k,j))*(u(i+1,k,j)+
     3          u(i,k,j))-(va(i,k,j+1)+va(i,k,j))*(u(i-1,k,j)+
     4          u(i,k,j)))/(dxx*msfd(i,j)*msfd(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*u(i,k,j)+fact2*u(i,k,j+1)
        else
         fx2=fact1*u(i,k,j+1)+fact2*u(i,k,j)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*u(i,k,j-1)+fact2*u(i,k,j)
        else
         fx1=fact1*u(i,k,j)+fact2*u(i,k,j-1)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*u(i,k,j)+fact2*u(i+1,k,j)
        else
         fy2=fact1*u(i+1,k,j)+fact2*u(i,k,j)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*u(i-1,k,j)+fact2*u(i,k,j)
        else
         fy1=fact1*u(i,k,j)+fact2*u(i-1,k,j)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfd(i,j)*msfd(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((u(i,k,j+1)+u(i,k,j))*ucmonb-
     &          (u(i,k,j)+u(i,k,j-1))*ucmonc+(u(i+1,k,j)+
     &          u(i,k,j))*vcmonb-(u(i,k,j)+u(i-1,k,j))*
     &          vcmonc)/(dxx*msfd(i,j)*msfd(i,j))
60     continue
c
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadv_V(ften,dxx,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfd   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#     include "include/cvaria.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)-((ua(i+1,k,j+1)+ua(i,k,j+1))*(v(i,k,j+1)+
     1          v(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))*(v(i,k,j)+
     2          v(i,k,j-1))+(va(i+1,k,j+1)+va(i+1,k,j))*(v(i+1,k,j)+
     3          v(i,k,j))-(va(i,k,j+1)+va(i,k,j))*(v(i-1,k,j)+
     4          v(i,k,j)))/(dxx*msfd(i,j)*msfd(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*v(i,k,j)+fact2*v(i,k,j+1)
        else
         fx2=fact1*v(i,k,j+1)+fact2*v(i,k,j)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*v(i,k,j-1)+fact2*v(i,k,j)
        else
         fx1=fact1*v(i,k,j)+fact2*v(i,k,j-1)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*v(i,k,j)+fact2*v(i+1,k,j)
        else
         fy2=fact1*v(i+1,k,j)+fact2*v(i,k,j)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*v(i-1,k,j)+fact2*v(i,k,j)
        else
         fy1=fact1*v(i,k,j)+fact2*v(i-1,k,j)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfd(i,j)*msfd(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((v(i,k,j+1)+v(i,k,j))*ucmonb-
     &          (v(i,k,j)+v(i,k,j-1))*ucmonc+(v(i+1,k,j)+
     &          v(i,k,j))*vcmonb-(v(i,k,j)+v(i-1,k,j))*
     &          vcmonc)/(dxx*msfd(i,j)*msfd(i,j))
60     continue
c
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadvQV(ften,dxx,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfx   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#     include "include/cvaria.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)-((ua(i+1,k,j+1)+ua(i,k,j+1))*(qv(i,k,j+1)+
     1          qv(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))*(qv(i,k,j)+
     2          qv(i,k,j-1))+(va(i+1,k,j+1)+va(i+1,k,j))*(qv(i+1,k,j)+
     3          qv(i,k,j))-(va(i,k,j+1)+va(i,k,j))*(qv(i-1,k,j)+
     4          qv(i,k,j)))/(dxx*msfx(i,j)*msfx(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*qv(i,k,j)+fact2*qv(i,k,j+1)
        else
         fx2=fact1*qv(i,k,j+1)+fact2*qv(i,k,j)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*qv(i,k,j-1)+fact2*qv(i,k,j)
        else
         fx1=fact1*qv(i,k,j)+fact2*qv(i,k,j-1)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*qv(i,k,j)+fact2*qv(i+1,k,j)
        else
         fy2=fact1*qv(i+1,k,j)+fact2*qv(i,k,j)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*qv(i-1,k,j)+fact2*qv(i,k,j)
        else
         fy1=fact1*qv(i,k,j)+fact2*qv(i-1,k,j)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfx(i,j)*msfx(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((qv(i,k,j+1)+qv(i,k,j))*ucmonb-
     &          (qv(i,k,j)+qv(i,k,j-1))*ucmonc+(qv(i+1,k,j)+
     &          qv(i,k,j))*vcmonb-(qv(i,k,j)+qv(i-1,k,j))*
     &          vcmonc)/(dxx*msfx(i,j)*msfx(i,j))
60     continue
c
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadvQC(ften,dxx,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfx   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#     include "include/cvaria.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)-((ua(i+1,k,j+1)+ua(i,k,j+1))*(qc(i,k,j+1)+
     1          qc(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))*(qc(i,k,j)+
     2          qc(i,k,j-1))+(va(i+1,k,j+1)+va(i+1,k,j))*(qc(i+1,k,j)+
     3          qc(i,k,j))-(va(i,k,j+1)+va(i,k,j))*(qc(i-1,k,j)+
     4          qc(i,k,j)))/(dxx*msfx(i,j)*msfx(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*qc(i,k,j)+fact2*qc(i,k,j+1)
        else
         fx2=fact1*qc(i,k,j+1)+fact2*qc(i,k,j)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*qc(i,k,j-1)+fact2*qc(i,k,j)
        else
         fx1=fact1*qc(i,k,j)+fact2*qc(i,k,j-1)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*qc(i,k,j)+fact2*qc(i+1,k,j)
        else
         fy2=fact1*qc(i+1,k,j)+fact2*qc(i,k,j)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*qc(i-1,k,j)+fact2*qc(i,k,j)
        else
         fy1=fact1*qc(i,k,j)+fact2*qc(i-1,k,j)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfx(i,j)*msfx(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((qc(i,k,j+1)+qc(i,k,j))*ucmonb-
     &          (qc(i,k,j)+qc(i,k,j-1))*ucmonc+(qc(i+1,k,j)+
     &          qc(i,k,j))*vcmonb-(qc(i,k,j)+qc(i-1,k,j))*
     &          vcmonc)/(dxx*msfx(i,j)*msfx(i,j))
60     continue
c
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hadvCH(ften,dxx,n,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the horizontal flux-divergence terms.  c
c     second-order difference is used.                                c
c                                                                     c
c     ften   : is the tendency for variable 'f'.                      c
c                                                                     c
c     ua, va : are p*u and p*v.                                       c
c                                                                     c
c     msfx   : is the map scale factor.                               c
c                                                                     c
c     dxx    : is the horizontal distance.                            c
c              = dx4  for ind=1.                                      c
c              = dx   for ind=2.                                      c
c              = dx16 for ind=3.                                      c
c                                                                     c
c     j      : is the j'th slice of f anf ften.                       c
c                                                                     c
c     ind = 1 : for t and qv.                                         c
c         = 2 : for qc and qr.                                        c
c         = 3 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer n,j,ind
      real(kind=8)  dxx
#     include "parameter.inc"
#     include "include/cvaria.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/param1.cb"
#     include "include/param3.cb"
c
      real(kind=8)  ften(ix,kx)
      integer i,k,jdx,jdxp1,jdxm1,idx,idxp1,idxm1
      real(kind=8)  fact1,fact2,uavg2,uavg1,fx2,fx1,vavg2,vavg1
      real(kind=8)  fy2,fy1,ucmona,vcmona,ucmonb,vcmonb,ucmonc,vcmonc
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----for t and qv:
c
       do 20 k=1,kx
       do 20 i=2,ilxm
        ften(i,k)=ften(i,k)
     & -((ua(i+1,k,j+1)+ua(i,k,j+1))*(chi(i,k,j+1,n)+chi(i,k,j,n))
     &  -(ua(i+1,k,j)+ua(i,k,j))*(chi(i,k,j,n)+chi(i,k,j-1,n))
     &  +(va(i+1,k,j+1)+va(i+1,k,j))*(chi(i+1,k,j,n)+chi(i,k,j,n))
     &  -(va(i,k,j+1)+va(i,k,j))*(chi(i-1,k,j,n)+chi(i,k,j,n)))
     &        /(dxx*msfx(i,j)*msfx(i,j))
20     continue
c
      else if(ind.eq.2)then
c
c   implement a "relaxed" upstream scheme
c
chy     fact1=0.75
       fact1=0.60
       fact2=1.-fact1
c
c-----for qc and qr:
c     up-wind values of qc and qr are used.
c
       do 40 k=1,kx
       do 40 i=2,ilxm
        uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
        uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
        if(uavg2.ge.0.)then
         fx2=fact1*chi(i,k,j,n)+fact2*chi(i,k,j+1,n)
        else
         fx2=fact1*chi(i,k,j+1,n)+fact2*chi(i,k,j,n)
        end if
        if(uavg1.ge.0.)then
         fx1=fact1*chi(i,k,j-1,n)+fact2*chi(i,k,j,n)
        else
         fx1=fact1*chi(i,k,j,n)+fact2*chi(i,k,j-1,n)
        end if
        vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
        vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
        if(vavg2.ge.0.)then
         fy2=fact1*chi(i,k,j,n)+fact2*chi(i+1,k,j,n)
        else
         fy2=fact1*chi(i+1,k,j,n)+fact2*chi(i,k,j,n)
        end if
        if(vavg1.ge.0.)then
         fy1=fact1*chi(i-1,k,j,n)+fact2*chi(i,k,j,n)
        else
         fy1=fact1*chi(i,k,j,n)+fact2*chi(i-1,k,j,n)
        end if
        ften(i,k)=ften(i,k)-(uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)/
     1          (dxx*msfx(i,j)*msfx(i,j))
40     continue
c
      else if(ind.eq.3)then
c
c-----for u and v:
c
       jdx=j
       jdxp1=j+1
       jdxm1=j-1
#ifdef MPP1
       if (myid .eq. 0 ) then
         jdxm1=max0(jdxm1,2)
       end if
       if (myid .eq. nproc-1) then
         jdxp1=min0(jdxp1,jendl-1)
       end if
#else
       jdxp1=min0(jdxp1,jlx)
       jdxm1=max0(jdxm1,2)
#endif
c
       do 60 k=1,kx
       do 60 i=2,ilx
        idx=i
        idxp1=i+1
        idxp1=min0(idxp1,ilx)
        idxm1=i-1
        idxm1=max0(idxm1,2)
        ucmona=ua(idxp1,k,jdx)+2.*ua(idx,k,jdx)+ua(idxm1,k,jdx)
        vcmona=va(idx,k,jdxp1)+2.*va(idx,k,jdx)+va(idx,k,jdxm1)
        ucmonb=ua(idxp1,k,jdxp1)+2.*ua(idx,k,jdxp1)+ua(idxm1,k,jdxp1)
     &         + ucmona
        vcmonb=va(idxp1,k,jdxp1)+2.*va(idxp1,k,jdx)+va(idxp1,k,jdxm1)
     &         + vcmona
        ucmonc=ua(idxp1,k,jdxm1)+2.*ua(idx,k,jdxm1)+ua(idxm1,k,jdxm1)
     &         + ucmona
        vcmonc=va(idxm1,k,jdxp1)+2.*va(idxm1,k,jdx)+va(idxm1,k,jdxm1)
     &         + vcmona
        ften(i,k)=ften(i,k)-((chi(i,k,j+1,n)+chi(i,k,j,n))*ucmonb
     &                      -(chi(i,k,j,n)+chi(i,k,j-1,n))*ucmonc
     &                      +(chi(i+1,k,j,n)+chi(i,k,j,n))*vcmonb
     &                      -(chi(i,k,j,n)+chi(i-1,k,j,n))*vcmonc)
     &           /(dxx*msfx(i,j)*msfx(i,j))
60     continue
c
      end if
c
      return
      end
