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

      subroutine cuparan(tten,qten,j)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer j
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/bats2.cb"
#     include "include/main.cb"

      real(kind=8)  tten(ix,kx),qten(ix,kx)
c
c--- dimensions needed for little driver routine
c
      real(kind=8)  usk(ix,kx),vsk(ix,kx)
      common /savetmp/usk,vsk
      real(kind=8)  t(ix,kx),q(ix,kx),p(ix,kx),po(ix,kx),vsp(ix,kx)
     &       ,tn(ix,kx),qo(ix,kx),outt(ix,kx),outq(ix,kx),pret(ix)
     &       ,psur(ix),ter11(ix),qcrit(ix)
      integer kdet(ix)
      integer i,k,icut,kk,istart,iend,iconj
      real(kind=8)  dtime,pkdcut,c0,us,vs,pkk,calc,aprdiv,prainx
c
c   zero out radiative clouds
c
      do 2 k=1,kx
      do 2 i=1,ilx
       cldlwc(i,k)=0.
       cldfra(i,k)=0.
2     continue

      icut=0
      dtime=dt
      pkdcut=75.
      c0=0.002
      istart=2+icut
      iend=ilxm-icut
c
c--- prepare input, erase output
c
      do 75 i=istart,iend
       kdet(i)=2
       qcrit(i)=0.
       pret(i)=0.
75    continue
      do 100 k=1,kx
      do 100 i=2+icut,ilxm-icut
       kk=kx-k+1
       us=ua(i,kk,j)/psb(i,j)
       vs=va(i,kk,j)/psb(i,j)
       us=0.25*(ua(i,kk,j  )/psb(i,j  )+ua(i+1,kk,j  )/psb(i+1,j)
     &         +ua(i,kk,j+1)/psb(i,j+1)+ua(i+1,kk,j+1)/psb(i+1,j+1))
       vs=0.25*(va(i,kk,j  )/psb(i,j  )+va(i+1,kk,j  )/psb(i+1,j)
     &         +va(i,kk,j+1)/psb(i,j+1)+va(i+1,kk,j+1)/psb(i+1,j+1))
       t(i,k)=tb(i,kk,j)/psb(i,j)
       q(i,k)=qvb(i,kk,j)/psb(i,j)
       if(q(i,k).lt.1.e-08)q(i,k)=1.e-08
       tn(i,k)=t(i,k)+(tten(i,kk))/psb(i,j)*dtime
       qo(i,k)=q(i,k)+(qten(i,kk))/psb(i,j)*dtime
       p(i,k)=10.*psb(i,j)*a(kk)+10.*ptop
       vsp(i,k)=dsqrt(us**2+vs**2)
       if(qo(i,k).lt.1.e-08)qo(i,k)=1.e-08
c
       po(i,k)=p(i,k)
       psur(i)=10.*psb(i,j)+10.*ptop
       outt(i,k)=0.
       pkk=psur(i)-po(i,k)
       if(pkk.le.pkdcut)kdet(i)=kdet(i)+1
       outq(i,k)=0.
       ter11(i)=ht(i,j)/9.81
       if(ter11(i).le.0.)ter11(i)=1.e-05
       qcrit(i)=qcrit(i)+qten(i,kk)
100   continue
c
c--- call cumulus parameterization
c
      call cup(qcrit,t,q,ter11,tn,qo,po,pret,p,
     &         outt,outq,dtime,psur,vsp,istart,iend,kdet,j)
      do 300 k=1,kx
      do 300 i=istart,iend
       if(pret(i).gt.0.) then
        kk=kx-k+1
        tten(i,kk)=psb(i,j)*outt(i,k)+tten(i,kk)
        qten(i,kk)=psb(i,j)*outq(i,k)+qten(i,kk)
       end if
300   continue
c
c--- rain in cm.
c
      calc=.5
      iconj=0
      do 50 i=istart,iend
       if(pret(i).gt.0.)then
        rainc(i,j)=rainc(i,j)+pret(i)*calc*dt
c       print *,'rainc(',i,j,')=',rainc(i,j)
        iconj=iconj+1
c.....................precipitation rate for bats (mm/s)
        aprdiv=dble(nbatst)
        if (jyear.eq.jyear0.and.ktau.eq.0) aprdiv=1.
        prainx=pret(i)*calc*dt
        pptc(i,j)=pptc(i,j)+prainx/(dtmin*60.)/aprdiv
c.......................................................
       end if
50    continue
      icon(j)=iconj
c
      return
      end subroutine cuparan
