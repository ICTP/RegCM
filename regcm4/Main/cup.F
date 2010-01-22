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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cup(qcrit,t,q,z1,tn,qo,po,pre,
     1  p,outtem,outq,dtime,psur,vsp,istart,iend,kdet,jslc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer istart,iend,jslc
      real(kind=8)  dtime
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/trachem.cb"

c
      real(kind=8)  ae(2),be(2),ht(2)
c
c--- input variables
c
      real(kind=8)  p(ix,kx),t(ix,kx),q(ix,kx),vsp(ix,kx),
     1        po(ix,kx),tn(ix,kx),qo(ix,kx),z1(ix),
     2        qcrit(ix)
c
c--- work variables
c
      real(kind=8)  he(ix,kx),hes(ix,kx),qes(ix,kx),
     1        z(ix,kx),tv(ix,kx),dby(ix,kx),
     2        qc(ix,kx),qrcd(ix,kx),pwd(ix,kx),
     3        pw(ix,kx),edt(ix),edto(ix),edtx(ix),
     4        aa1(ix),aa0(ix),xaa0(ix),hkb(ix),hkbo(ix),
     5        xhkb(ix),qkb(ix),qkbo(ix),xmb(ix)
      real(kind=8)  heo(ix,kx),heso(ix,kx),qeso(ix,kx),
     1        zo(ix,kx),tvo(ix,kx),dbyo(ix,kx),
     2        qco(ix,kx),qrcdo(ix,kx),pwdo(ix,kx),
     3        pwo(ix,kx),pre(ix),psur(ix),
     4        xpwav(ix),xpwev(ix),pwav(ix),pwev(ix),
     5        pwavo(ix),pwevo(ix),hcd(ix),qck(ix),qcko(ix),
     6        hcdo(ix),qcd(ix),qcdo(ix),bu(ix),buo(ix),
     7        xqkb(ix),xhcd(ix),xqcd(ix),xqck(ix)
      real(kind=8)  xhe(ix,kx),xhes(ix,kx),xqes(ix,kx),
     1        xz(ix,kx),xtv(ix,kx),xt(ix,kx),xq(ix,kx),
     2        xdby(ix,kx),xqc(ix,kx),xqrcd(ix,kx),
     3        xpwd(ix,kx),xpw(ix,kx),vshear(ix)
      real(kind=8)  dellah(ix,kx),dellaq(ix,kx),dellat(ix,kx)
      real(kind=8)  dkk(ix,kx)
      integer kdet(ix)
      integer k22(ix),kbcon(ix),kb(ix),jmin(ix),ktop(ix),
     1        kds(ix),kbhe(ix)
      real(kind=8)  f,xk,mbdt
c
c----output variables
c
      real(kind=8)  outtem(ix,kx),outq(ix,kx)
      integer i,k,ierr,iph,ipho,kk,lpt,kbcono,kclth
      real(kind=8)  tcrit,xl,cpinv,tfinv,alsixt,c0,e,eo,tvbar,tvbaro
      real(kind=8)  dh,pbcdif,dz,dz1,gamma,gamma0,qrch,dzo,dz2
      real(kind=8)  gammo,gammo0,qrcho,dq,dt,dp,dv1,dv2,dv3
      real(kind=8)  dv1q,dv2q,dv3q,detdo,detdoq,aup,adw
      real(kind=8)  gamma1,gamma2,dhh,dg,xddd,outtes,akclth
c
      ierr=0
      tcrit=50.
      xl=2.5e06
      rv=461.9

      cpinv = 1./cp
      tfinv = 1./273.
      alsixt = dlog(610.71d0)
      ht(1)= xl * cpinv
      ht(2)= 2.834e6 * cpinv
      be(1)= .622*ht(1) * 3.50
      ae(1)= be(1)*tfinv + alsixt
      be(2)= .622*ht(2) * 3.50
      ae(2)= be(2)*tfinv + alsixt
      mbdt=dtime*5.e-03
      c0=.002
      f=-1.
      xk=-1.
c
c--- environmental conditions, first heights
c
      do 5 k=1,kx
      do 5 i=istart,iend
       dkk(i,k)=1.
       iph=1
       ipho=1
       if(t(i,k).le.tcrit)iph=2
       if(tn(i,k).le.tcrit)ipho=2
       e=dexp(ae(iph)-be(iph)/t(i,k))
       eo=dexp(ae(ipho)-be(ipho)/tn(i,k))
       qes(i,k)=.622*e/(100.*p(i,k)-(1.-.622)*e)
       qeso(i,k)=.622*eo/(100.*po(i,k)-(1.-.622)*eo)
       if(qes(i,k).le.1.e-08)qes(i,k)=1.e-08
       if(q(i,k).gt.qes(i,k))q(i,k)=qes(i,k)
       if(qeso(i,k).le.1.e-08)qeso(i,k)=1.e-08
       if(qo(i,k).gt.qeso(i,k))qo(i,k)=qeso(i,k)
       tv(i,k)=t(i,k)+.608*q(i,k)*t(i,k)
       tvo(i,k)=tn(i,k)+.608*qo(i,k)*tn(i,k)
5     continue
      do i=1,ix
       hkb(i)=0.    ! EES
       qkb(i)=0.
       hkbo(i)=0.
       qkbo(i)=0.
       xhkb(i)=0.
       xqkb(i)=0.
       edt(i)=0.
       edto(i)=0.
       edtx(i)=0.
      enddo
      do 7 i=istart,iend
c      hkb(i)=0.
c      qkb(i)=0.
c      hkbo(i)=0.
c      qkbo(i)=0.
c      xhkb(i)=0.
c      xqkb(i)=0.
       aa1(i)=0.
       aa0(i)=0.
       if(qcrit(i).le.0.)aa0(i)=-1.
       xaa0(i)=0.
       xpwav(i)=0.
       xpwev(i)=0.
       pwav(i)=0.
       pwev(i)=0.
       pwavo(i)=0.
       pwevo(i)=0.
       k22(i)=1
       kbhe(i)=1
       ktop(i)=1
       kbcon(i)=1
       kb(i)=1
       kds(i)=1
       jmin(i)=1
c      edt(i)=0.
c      edto(i)=0.
c      edtx(i)=0.
       xmb(i)=0.
       vshear(i)=0.
       z(i,1)=z1(i)-(dlog(p(i,1))-dlog(psur(i)))*287.*tv(i,1)/9.81
       zo(i,1)=z1(i)-(dlog(po(i,1))-dlog(psur(i)))*287.*tvo(i,1)/9.81
7     continue
      do 10 k=2,kx
      do 10 i=istart,iend
       tvbar=.5*tv(i,k)+.5*tv(i,k-1)
       z(i,k)=z(i,k-1)-(dlog(p(i,k))-dlog(p(i,k-1)))*287.*tvbar/9.81
       tvbaro=.5*tvo(i,k)+.5*tvo(i,k-1)
       zo(i,k)=zo(i,k-1)-(dlog(po(i,k))-dlog(po(i,k-1)))*287.*tvbaro
     1         /9.81
10    continue
c
c--- moist static energy
c
      do 20 k=1,kx
      do 20 i=istart,iend
       cldlwc(i,k)=0.
       cldfra(i,k)=0.
       pw(i,k)=0.
       xpw(i,k)=0.
       pwo(i,k)=0.
       qc(i,k)=0.
       xqc(i,k)=0.
       qco(i,k)=0.
       pwd(i,k)=0.
       pwdo(i,k)=0.
       xpwd(i,k)=0.
       dellah(i,k)=0.
       dellaq(i,k)=0.
       dellat(i,k)=0.
       he(i,k)=9.81*z(i,k)+1004.*t(i,k)+2.5e06*q(i,k)
       hes(i,k)=9.81*z(i,k)+1004.*t(i,k)+2.5e06*qes(i,k)
       if(he(i,k).ge.hes(i,k))he(i,k)=hes(i,k)
       heo(i,k)=9.81*zo(i,k)+1004.*tn(i,k)+2.5e06*qo(i,k)
       heso(i,k)=9.81*zo(i,k)+1004.*tn(i,k)+2.5e06*qeso(i,k)
       if(heo(i,k).ge.heso(i,k))heo(i,k)=heso(i,k)
       xt(i,k)=t(i,k)
       xq(i,k)=q(i,k)
       xhe(i,k)=he(i,k)
       if(k.ne.kx)qrcd(i,k)=.5*(qes(i,k)+qes(i,k+1))
       if(k.ne.kx)qrcdo(i,k)=.5*(qeso(i,k)+qeso(i,k+1))
20    continue
c
c------- determine level with highest moist static energy content.
c
      call maximi(he,ix,kx,1,kbmax2d(i,jslc),k22,istart,iend)
      do 36 i=istart,iend
       if(aa0(i).ge.0.) then
        if(k22(i).ge.kbmax2d(i,jslc))then
         aa0(i)=-1.
         go to 36
        endif
        hkb(i)=he(i,k22(i))
        qkb(i)=q(i,k22(i))
        hkbo(i)=heo(i,k22(i))
        qkbo(i)=qo(i,k22(i))
        qck(i)=qkb(i)
        qcko(i)=qkbo(i)
       endif
36    continue
c
c--- decide for convective cloud base
c
      do 27 i=istart,iend
       if(aa0(i).ge.0.)then
        do 28 k=1,kdet(i)
         kk=kdet(i)-k+1
c        dkk(i,kk)=.75*dkk(i,kk+1)
         dkk(i,k)=1.-dble(kk)/dble(kdet(i))
28      continue
3777    continue
        kb(i)=k22(i)
        kbcon(i)=kb(i)
c----------------------------------
        go to 32
 31     continue
        kbcon(i)=kbcon(i)+1
        if(kbcon(i).gt.kbmax2d(i,jslc))then
         aa0(i)=-1.
         go to 27
        endif
 32     continue
        dh=.5*hes(i,kbcon(i))+.5*hes(i,kbcon(i)+1)
        if(hkb(i).lt.dh)go to 31
c
c--- after large-scale forcing is applied, possible lid should be
c--- removed!!!
c
        kbcono=kb(i)
        go to 332
331     continue
        kbcono=kbcono+1
332     continue
cictp
        if(kbcono.gt.kbmax2d(i,jslc))then
         aa0(i)=-1.
         go to 27
        endif
cictp_
        dh=.5*heso(i,kbcono)+.5*heso(i,kbcono+1)
        if(hkbo(i).lt.dh)go to 331
        pbcdif=-p(i,kbcono)+p(i,kb(i))
c----------------------------- below was commented out
cas uncommenting the following lines for experiment 2/5/95
        if(pbcdif.gt.pbcmax2d(i,jslc))then  !this is where typo was (pbdcdif)
         k22(i)=k22(i)+1
         if(k22(i).ge.kbmax2d(i,jslc))then
          aa0(i)=-1.
          go to 27
         endif
         hkb(i)=he(i,k22(i))
         qkb(i)=q(i,k22(i))
         hkbo(i)=heo(i,k22(i))
         qkbo(i)=qo(i,k22(i))
         qck(i)=qkb(i)
         qcko(i)=qkbo(i)
         go to 3777
        endif
cas
       endif
27    continue

c
c--- downdraft originating level
c
      call minimi(he,ix,kx,kb,kx,jmin,istart,iend)
      call maximi(vsp,ix,kx,1,kx,kds,istart,iend)
c
c**************************** static control
c
c
c--- determine cloud top
c
      do 37 i=istart,iend
       if(aa0(i).ge.0)then
        if(jmin(i).le.3)then
         aa0(i)=-1.
         go to 37
        endif
        if(kds(i).ge.kx)kds(i)=kx-1
        if(kds(i).le.kbcon(i))kds(i)=kbcon(i)
        dby(i,kx)=hkb(i)-hes(i,kx)
        dbyo(i,kx)=hkbo(i)-heso(i,kx)
       endif
37    continue
      do 35 k=1,kx-1
      do 35 i=istart,iend
       if(aa0(i).ne.-1.)then
        dby(i,k)=hkb(i)-.5*(hes(i,k)+hes(i,k+1))
        dbyo(i,k)=hkbo(i)-.5*(heso(i,k)+heso(i,k+1))
       endif
35    continue
      do 42 i=istart,iend
       if(aa0(i).ne.-1.)then
        do 40 k=2,kx-kbcon(i)-1
         kk=kx-k+1
         if(dby(i,kk).ge.0.)then
          ktop(i)=kk+1
          go to 41
         endif
40      continue
        aa0(i)=-1.
        go to 42
41      continue
        if(ktop(i).gt.kx)ktop(i)=kx
        if(p(i,kbcon(i))-p(i,ktop(i)).lt.mincld2d(i,jslc))aa0(i)=-1.
       endif
42    continue
c
   
c------- moisture and cloud work functions
c
      do 55 k=2,kx-1
      do 55 i=istart,iend
       if(aa0(i).ne.-1.)then
        if(k.gt.kbcon(i))then
         if(k.lt.ktop(i))then
          dz=-.5*z(i,k-1)+.5*z(i,k+1)
          dz1=z(i,k)-z(i,k-1)
          gamma=(xl/cp)*(xl/(rv*(t(i,k)**2)))*qes(i,k)
          gamma0=(xl/cp)*(xl/(rv*(t(i,k-1)**2)))*qes(i,k-1)
          qrch=qes(i,k)+(1./xl)*(gamma/(1.+gamma))*dby(i,k)
          qc(i,k)=(qck(i)-qrch)/(1.+c0*dz)+qrch
          pw(i,k)=c0*dz*(qc(i,k)-qrch)
          qck(i)=qc(i,k)
          pwav(i)=pwav(i)+pw(i,k)
          dz1=z(i,k)-z(i,k-1)
          aa0(i)=aa0(i)+dz1*(9.81/(1004.*(.5*
     *          (t(i,k)+t(i,k-1)))))*dby(i,k-1)/
     1       (1.+.5*gamma+.5*gamma0)
          dzo=-.5*zo(i,k-1)+.5*zo(i,k+1)
          dz2=zo(i,k)-zo(i,k-1)
          gammo=(xl/cp)*(xl/(rv*(tn(i,k)**2)))*qeso(i,k)
          gammo0=(xl/cp)*(xl/(rv*(tn(i,k-1)**2)))*qeso(i,k-1)
          qrcho=qeso(i,k)+(1./xl)*(gammo/(1.+gammo))*dbyo(i,k)
          qco(i,k)=(qcko(i)-qrcho)/(1.+c0*dzo)+qrcho
          pwo(i,k)=c0*dzo*(qco(i,k)-qrcho)
          qcko(i)=qco(i,k)
          pwavo(i)=pwavo(i)+pwo(i,k)
          aa1(i)=aa1(i)+dz2*(9.81/(1004.*(.5*
     1       (tn(i,k)+tn(i,k-1)))))*dbyo(i,k-1)/
     1       (1.+.5*gammo+.5*gammo0)
         endif
        endif
       endif
55    continue
c
c
      do 50 i=istart,iend
       if(aa0(i).ne.-1.)then
        k=ktop(i)
        dz=-.5*z(i,k-1)+.5*z(i,k)
        gamma=(xl/cp)*(xl/(rv*(t(i,k)**2)))*qes(i,k)
        qrch=qes(i,k)+(1./xl)*(gamma/(1.+gamma))*dby(i,k)
        qc(i,k)=qes(i,k)
        pw(i,k)=(qrch-qes(i,k))
        pwav(i)=pwav(i)+pw(i,k)
c
        dz=-.5*zo(i,k-1)+.5*zo(i,k)
        gamma=(xl/cp)*(xl/(rv*(tn(i,k)**2)))*qeso(i,k)
        qrcho=qeso(i,k)+(1./xl)*(gamma/(1.+gamma))*dbyo(i,k)
        qco(i,k)=qeso(i,k)
        pwo(i,k)=(qrcho-qeso(i,k))
        pwavo(i)=pwavo(i)+pwo(i,k)
       endif
50    continue
c
c------- downdraft calculations
c
c
c--- determine downdraft strength in terms of windshear
c
      do 62 kk=1,kx/2
      do 62 i=istart,iend
       if(aa0(i).ne.-1.)then
        vshear(i)=vshear(i)+dabs((vsp(i,kk+1)-vsp(i,kk))
     1             /(z(i,kk+1)-z(i,kk)))
       endif
62    continue
      do 72 i=istart,iend
       if(aa0(i).ne.-1.)then
        vshear(i)=1.e3*vshear(i)/dble(kx/2)
        edt(i)=1.-(1.591-.639*vshear(i)+.0953*(vshear(i)**2)
     1      -.00496*(vshear(i)**3))

        if(edt(i).gt.shrmax2d(i,jslc)) edt(i)=shrmax2d(i,jslc)
        if(edt(i).lt.shrmin2d(i,jslc)) edt(i)=shrmin2d(i,jslc)

        edto(i)=edt(i)
        edtx(i)=edt(i)
        qrcd(i,kx)=qes(i,kx)
        hcd(i)=.5*(he(i,jmin(i))+he(i,jmin(i)+1))
        qcd(i)=.5*(q(i,jmin(i))+q(i,jmin(i)+1))
        qrcdo(i,kx)=qeso(i,kx)
        hcdo(i)=heso(i,kx)
        hcdo(i)=.5*(heo(i,jmin(i))+heo(i,jmin(i)+1))
        qcdo(i)=.5*(qo(i,jmin(i))+qo(i,jmin(i)+1))
        bu(i)=0.
        buo(i)=0.
       endif
72    continue
      do 95 k=1,kx-1
      do 95 i=istart,iend
       if(aa0(i).ne.-1.)then
        if(k.lt.jmin(i))then
         kk=jmin(i)-k
         dz=-(z(i,kk)-z(i,kk+2))*.5
         bu(i)=bu(i)+dz*(hcd(i)-.5*(hes(i,kk)+hes(i,kk+1)))
         dq=(qes(i,kk)+qes(i,kk+1))*.5
         dt=(t(i,kk)+t(i,kk+1))*.5
         gamma=(xl/cp)*(xl/(rv*(dt**2)))*dq
         dh=hcd(i)-.5*(hes(i,kk)+hes(i,kk+1))
         qrcd(i,kk)=(dq+(1./xl)*(gamma/(1.+gamma))*dh)
         pwd(i,kk)=dkk(i,kk)*(qcd(i)-qrcd(i,kk))
         qcd(i)=qrcd(i,kk)
         pwev(i)=pwev(i)+pwd(i,kk)
c
         dz=-(zo(i,kk)-zo(i,kk+2))*.5
         buo(i)=buo(i)+dz*(hcdo(i)-.5*(heso(i,kk)+heso(i,kk+1)))
         dq=(qeso(i,kk)+qeso(i,kk+1))*.5
         dt=(tn(i,kk)+tn(i,kk+1))*.5
         gamma=(xl/cp)*(xl/(rv*(dt**2)))*dq
         dh=hcdo(i)-.5*(heso(i,kk)+heso(i,kk+1))
         qrcdo(i,kk)=(dq+(1./xl)*(gamma/(1.+gamma))*dh)
         pwdo(i,kk)=dkk(i,kk)*(qcdo(i)-qrcdo(i,kk))
         qcdo(i)=qrcdo(i,kk)
         pwevo(i)=pwevo(i)+pwdo(i,kk)
        endif
       endif
95    continue
c
      do 96 i=istart,iend
       if(aa0(i).ne.-1.)then
        if(bu(i).ge.0.or.buo(i).ge.0.or.pwev(i).ge.0.or.
     1        pwevo(i).ge.0.)aa0(i)=-1.
        edt(i)=-edt(i)*pwav(i)/pwev(i)
        if(edt(i).gt.edtmax2d(i,jslc))edt(i)=edtmax2d(i,jslc)
        if(edt(i).lt.edtmin2d(i,jslc))edt(i)=edtmin2d(i,jslc)
        edto(i)=-edto(i)*pwavo(i)/pwevo(i)
        if(edto(i).gt.edtmaxo2d(i,jslc))edto(i)=edtmaxo2d(i,jslc)
        if(edto(i).lt.edtmino2d(i,jslc))edto(i)=edtmino2d(i,jslc)
       endif
96    continue
c
c--- what would the change be?
c
      do 98 i=istart,iend
       if(aa0(i).ne.-1.)then
        k=1
        dz=.5*(z(i,2)-z(i,1))
        dp=50.*(psur(i)-p(i,2))
        dellah(i,1)=edt(i)*(dkk(i,1)*hcd(i)-dkk(i,1)*
     *          .5*(he(i,1)+he(i,2)))*g/dp
        dellaq(i,1)=edt(i)*(dkk(i,1)*qrcd(i,1)-dkk(i,1)*
     *          .5*(q(i,1)+q(i,2)))*g/dp
        xhe(i,k)=dellah(i,k)*mbdt+he(i,k)
        xq(i,k)=dellaq(i,k)*mbdt+q(i,k)
        dellat(i,k)=(1./1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))
        xt(i,k)=(mbdt/1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))+t(i,k)
        if(xq(i,k).le.0.)then
         xq(i,k)=1.e-08
        endif
       endif
98    continue
c
      do 100 k=1,kx-1
      do 100 i=istart,iend
       if(aa0(i).ne.-1.)then
        if(k.eq.1.or.k.ge.ktop(i))goto 100
        dv1=.5*(he(i,k)+he(i,k+1))
        dv2=he(i,k)
        dv3=.5*(he(i,k)+he(i,k-1))
        dv1q=.5*(q(i,k)+q(i,k+1))
        dv2q=q(i,k)
        dv3q=.5*(q(i,k)+q(i,k-1))
c
c--- specifiy detrainment of downdraft, has to be consistent
c--- with zd calculations in soundd.
c
        detdo=(1.-dkk(i,k))*(hcd(i)-dv2)
        detdoq=(1.-dkk(i,k))*(qrcd(i,k)-dv2q)
        dz=.5*(z(i,k+1)-z(i,k-1))
c
c--- changed due to subsidence and entrainment
c
        aup=1.
        if(k.le.k22(i))aup=0.
        adw=1.
        if(k.gt.jmin(i))adw=0.
        dp=+50.*(p(i,k-1)-p(i,k+1))
        dellah(i,k)=((aup-adw*edt(i))*(dv1-dv2)
     *            +(aup-adw*edt(i))*(dv2-dv3))*g/dp
     2            +adw*edt(i)*detdo*g/dp
        dellaq(i,k)=((aup-adw*edt(i))*(dv1q-dv2q)
     *             +(aup-adw*edt(i))*(dv2q-dv3q))*g/dp
     2             +adw*edt(i)*detdoq*g/dp
        xhe(i,k)=dellah(i,k)*mbdt+he(i,k)
        xq(i,k)=dellaq(i,k)*mbdt+q(i,k)
        dellat(i,k)=(1./1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))
        xt(i,k)=(mbdt/1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))+t(i,k)
        if(xq(i,k).le.0.)then
         xq(i,k)=1.e-08
        endif
       endif
100   continue
c
c------- cloud top
c
      do 102 i=istart,iend
       if(aa0(i).ne.-1.)then
        lpt=ktop(i)
        dp=100.*(p(i,lpt-1)-p(i,lpt))
        dv1=.5*(he(i,lpt)+he(i,lpt-1))
        dellah(i,lpt)=( hkb(i) -dv1)*g/dp
        dv1=.5*(q(i,lpt)+q(i,lpt-1))
        dellaq(i,lpt)=(qes(i,lpt)-dv1)*g/dp
        k=lpt
        xhe(i,k)=dellah(i,k)*mbdt+he(i,k)
        xq(i,k)=dellaq(i,k)*mbdt+q(i,k)
        dellat(i,k)=(1./1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))
        xt(i,k)=(mbdt/1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))+t(i,k)
        if(xq(i,k).le.0.)then
         xq(i,k)=1.e-08
        endif
        xhkb(i)=dellah(i,kbcon(i))*mbdt+hkb(i)
        xqkb(i)=dellaq(i,kbcon(i))*mbdt+qkb(i)
        if(xqkb(i).le.0.)xqkb(i)=1.e-08
       endif
102   continue
c
c--- environmental conditions, first heights
c
      do 115 k=1,kx
      do 115 i=istart,iend
       if(aa0(i).ne.-1.)then
        iph=1
        if(xt(i,k).le.tcrit)iph=2
        e=dexp(ae(iph)-be(iph)/xt(i,k))
        xqes(i,k)=.622*e/(100.*p(i,k)-(1.-.622)*e)
        if(xqes(i,k).le.1.e-08)xqes(i,k)=1.e-08
        if(xq(i,k).gt.xqes(i,k))xq(i,k)=xqes(i,k)
        xtv(i,k)=xt(i,k)+.608*xq(i,k)*xt(i,k)
       endif
115   continue
c bug fix
       do k=1,kx-1
       do i=istart,iend
        if(aa0(i).ne.-1)xqrcd(i,k)=.5*(xqes(i,k)+xqes(i,k+1))
       enddo
       enddo
c
      do 116 i=istart,iend
       if(aa0(i).ne.-1.)then
        xz(i,1)=z1(i)-(dlog(p(i,1))-dlog(psur(i)))*287.*xtv(i,1)
     1          /9.81
       endif
116   continue
      do 120 k=2,kx
      do 120 i=istart,iend
       if(aa0(i).ne.-1.)then
        tvbar=.5*xtv(i,k)+.5*xtv(i,k-1)
        xz(i,k)=xz(i,k-1)-(dlog(p(i,k))-dlog(p(i,k-1)))*287.*tvbar
     1          /9.81
       endif
120   continue
c
c--- moist static energy
c
      do 130 k=1,kx
      do 130 i=istart,iend
       if(aa0(i).ne.-1.)then
        xhes(i,k)=9.81*xz(i,k)+1004.*xt(i,k)+2.5e06*xqes(i,k)
        if(xhe(i,k).ge.xhes(i,k))xhe(i,k)=xhes(i,k)
       endif
130   continue
c
c
c**************************** static control
c
      do 132 i=istart,iend
       if(aa0(i).ne.-1.)then
        xqck(i)=xqkb(i)
        xdby(i,kx)=xhkb(i)-xhes(i,kx)
       endif
132   continue
c
c------- moisture and cloud work functions
c
      do 140 k=1,kx-1
      do 140 i=istart,iend
       if(aa0(i).ge.0.)then
        xdby(i,k)=xhkb(i)-.5*(xhes(i,k)+xhes(i,k+1))
        if(k.le.kbcon(i).or.k.ge.ktop(i))goto 140
        dz=-.5*xz(i,k-1)+.5*xz(i,k+1)
        dz1=xz(i,k)-xz(i,k-1)
        gamma=(xl/cp)*(xl/(rv*(xt(i,k)**2)))*xqes(i,k)
        gamma0=(xl/cp)*(xl/(rv*(xt(i,k-1)**2)))*xqes(i,k-1)
        qrch=xqes(i,k)+(1./xl)*(gamma/(1.+gamma))*xdby(i,k)
        xqc(i,k)=(xqck(i)-qrch)/(1.+c0*dz)+qrch
        xpw(i,k)=c0*dz*(xqc(i,k)-qrch)
        xqck(i)=xqc(i,k)
        xpwav(i)=xpwav(i)+xpw(i,k)
        xaa0(i)=xaa0(i)+dz1*(9.81/(1004.*(.5*
     1          (xt(i,k)+xt(i,k-1)))))*xdby(i,k-1)/
     1       (1.+.5*gamma+.5*gamma0)
       endif
140   continue
      do 142 i=istart,iend
       if(aa0(i).ge.0.)then
        k=ktop(i)
        dz=-.5*xz(i,k-1)+.5*xz(i,k)
        gamma=(xl/cp)*(xl/(rv*(xt(i,k)**2)))*xqes(i,k)
        qrch=xqes(i,k)+(1./xl)*(gamma/(1.+gamma))*xdby(i,k)
        xqc(i,k)=xqes(i,k)
        xpw(i,k)=(qrch-xqes(i,k))
        xpwav(i)=xpwav(i)+xpw(i,k)
        xqrcd(i,kx)=xqes(i,kx)
        xhcd(i)=.5*(xhe(i,jmin(i))+xhe(i,jmin(i)+1))
        xqcd(i)=.5*(xq(i,jmin(i))+xq(i,jmin(i)+1))
        xpwev(i)=0.
        bu(i)=0.
       endif
142   continue
c
c------- downdraft calculations
c
c
c--- downdraft moisture properties
c
      do 160 k=1,kx-1
      do 160 i=istart,iend
       if(aa0(i).ge.0.)then
        if(k.lt.jmin(i))then
         kk=jmin(i)-k
         dz=-(xz(i,kk)-xz(i,kk+2))*.5
         bu(i)=bu(i)+dz*(xhcd(i)-.5*(xhes(i,kk)+xhes(i,kk+1)))
         dq=(xqes(i,kk)+xqes(i,kk+1))*.5
         dt=(xt(i,kk)+xt(i,kk+1))*.5
         gamma=(xl/cp)*(xl/(rv*(dt**2)))*dq
         dh=xhcd(i)-.5*(xhes(i,kk)+xhes(i,kk+1))
         xqrcd(i,kk)=(dq+(1./xl)*(gamma/(1.+gamma))*dh)
         xpwd(i,kk)=dkk(i,kk)*(xqcd(i)-xqrcd(i,kk))
         xqcd(i)=xqrcd(i,kk)
         xpwev(i)=xpwev(i)+xpwd(i,kk)
        endif
       endif
160   continue
      do 162 i=istart,iend
       if(aa0(i).ge.0.)then
        if(bu(i).ge.0.)then
         aa0(i)=-1.
         go to 162
        endif
        if(xpwev(i).ne.0.)edtx(i)=-edtx(i)*xpwav(i)/xpwev(i)
        if(edtx(i).gt.edtmaxx2d(i,jslc))edtx(i)=edtmaxx2d(i,jslc)
        if(edtx(i).lt.edtminx2d(i,jslc))edtx(i)=edtminx2d(i,jslc)
       endif
162   continue
c
c
c--- downdraft cloudwork functions
c
c
      do 171 k=1,kx-1
      do 171 i=istart,iend
       if(aa0(i).ge.0.)then
        if(k.lt.jmin(i))then
         kk=jmin(i)-k
c
c--- original
c
         gamma1=(xl/cp)*(xl/(rv*(t(i,kk)**2)))*qes(i,kk)
         gamma2=(xl/cp)*(xl/(rv*(t(i,kk+1)**2)))*qes(i,kk+1)
         dhh=hcd(i)
         dt=.5*(t(i,kk)+t(i,kk+1))
         dg=.5*(gamma1+gamma2)
         dh=.5*(hes(i,kk)+hes(i,kk+1))
         dz=(z(i,kk)-z(i,kk+1))*dkk(i,kk)
         xddd=dhh-dh
         aa0(i)=aa0(i)+edt(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))
c
c--- modified by larger scale
c
         gamma1=(xl/cp)*(xl/(rv*(tn(i,kk)**2)))*qeso(i,kk)
         gamma2=(xl/cp)*(xl/(rv*(tn(i,kk+1)**2)))*qeso(i,kk+1)
         dhh=hcdo(i)
         dt=.5*(tn(i,kk)+tn(i,kk+1))
         dg=.5*(gamma1+gamma2)
         dh=.5*(heso(i,kk)+heso(i,kk+1))
         dz=(zo(i,kk)-zo(i,kk+1))*dkk(i,kk)
         aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))
c
c--- modified by cloud
c
         gamma1=(xl/cp)*(xl/(rv*(xt(i,kk)**2)))*xqes(i,kk)
         gamma2=(xl/cp)*(xl/(rv*(xt(i,kk+1)**2)))*xqes(i,kk+1)
         dhh=xhcd(i)
         dt=.5*(xt(i,kk)+xt(i,kk+1))
         dg=.5*(gamma1+gamma2)
         dh=.5*(xhes(i,kk)+xhes(i,kk+1))
         dz=(xz(i,kk)-xz(i,kk+1))*dkk(i,kk)
         xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))
        endif
       endif
171   continue
c
c--- large scale forcing
c
      do 175 i=istart,iend
       if(aa0(i).ge.0.)then
        if (igcc.eq.1) then
          f= (aa1(i)-aa0(i))/dtime     ! Arakawa-Schubert closure
        elseif (igcc.eq.2) then
          f=aa0(i)/dtauc2d(i,jslc)       ! Fritsch-Chappell closure
        end if
        xk=(xaa0(i)-aa0(i))/mbdt
        xmb(i)=-f/xk
        if(f.le.0.or.xk.ge.0.)xmb(i)=0.
       endif
175   continue
Cchem2
      mflx(i,1) = xmb(i)
      mflx(i,2) = xmb(i)*edt(i)
Cchem2_
c
c--- feedback
c
      do 180 k=1,kx
      do 180 i=istart,iend
       if(aa0(i).ge.0.)then
        if(k.le.ktop(i))then
         outtes=dellat(i,k)*xmb(i)*86400.
         if((outtes.gt.htmax2d(i,jslc)).or.
     &      (outtes.lt.htmin2d(i,jslc))) then
           xmb(i)=0.
           aa0(i)=-1.
         else
           outtem(i,k)=outtem(i,k)+dellat(i,k)*xmb(i)
           outq(i,k)=outq(i,k)+dellaq(i,k)*xmb(i)
           pre(i)=pre(i)+(pw(i,k)+edt(i)*pwd(i,k))*xmb(i)
         end if
        endif
       endif
180   continue
c
c    calculate cloud fraction and water content
c
      do 205 i=istart,iend
Cchem2
       icumtop(i,jslc)=0
       icumbot(i,jslc)=0
       icumdwd(i,jslc)=0
Cchem2_
       if(aa0(i).ge.0.)then

        if (ktop(i).gt.1 .and. kbcon(i).gt.1) then
         kclth=ktop(i)-kbcon(i)+1
         akclth=1./dble(kclth)
         do 203 k=kbcon(i),ktop(i)
          kk=kx-k+1
          cldlwc(i,kk)=cllwcv
          cldfra(i,kk)=1.-(1.-clfrcv)**akclth
203      continue
Cchem2
Cchem2 define convection  base and top for tracers 
         if(ichem .eq. 1) then   
          if(ktop(i).gt. 1 .and. k22(i).ge.1 )then 
           icumtop(i,jslc)= kx+1- ktop(i)
           icumbot(i,jslc)= kx+1- k22(i) 
           icumdwd(i,jslc)= kx+1- jmin(i)
          end if
         end if
Cchem2_  
        endif
       endif

205   continue

      return
      end
