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

      subroutine cupara(j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine performs cumulus parameterization scheme.       c
c     the basic method follows anthes and keyser (1979) and           c
c     kuo (1983).                                                     c
c                                                                     c
c     all the other arguments are passed from subroutine "tend" and   c
c     explained in "tend".                                            c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/bats.cb2"
#     include "include/trachem.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif

c
      real(kind=8)  seqt(kx),tmp3(ix,kx)
c
c----------------------------------------------------------------------
c
c     crneg  : critical negative area in the sounding
c     pert   : perturbation temperature
c     perq   : perturbation mixing ratio
c     dlt    : temperature difference used to allow over shooting.
c     cdscld : critical cloud depth in delta sigma.
c
      real(kind=8)  crneg,pert,perq
      data    crneg, pert, perq/1.,1.,1.e-3/
      real(kind=8)  dlt,cdscld
      data    dlt, cdscld/3.0,0.3/
c
      integer i,k,kbase,kk,ktop,kbaseb,kclth
      real(kind=8)  sca,eqtm,ttp,q,psg,t1,eqt,tmax,qmax,pmax,emax,tdmax
      real(kind=8)  dalr,dplr,zlcl,tlcl,tmean,dlnp,plcl,siglcl,es,qs
      real(kind=8)  deqt,dsc,ttsum,sum,sumb,arh,psx,pux,e1,rh,xsav,c301
      real(kind=8)  ttconv,apcnt,eddyf,akclth,prainx,aprdiv,rsht,rswt
c
      do 2 k=1,kx
      do 2 i=1,ilx
       cldlwc(i,k)=0.
       cldfra(i,k)=0.
2     continue
      do 10 k=1,kx
      do 10 i=1,ilx
       qvten(i,k,j) = 0.
10    continue
c
c-----compute the horizontal advection terms:
c
      call hadvQV(qvten(1,1,j),dx4,j,1)
c---------------
Cchem2  
      if(ichem.eq.1) then
c
c     icumtop = top level of cumulus clouds
c     icumbot = bottom level of cumulus clouds
c           (calculated in cupara and stored for tractend)
c  before do 100 put
      do i=2,ilxm
         icumtop(i,j) = 0
         icumbot(i,j) = 0
      end do
      end if
Cchem2__
c
c-----compute the moisture convergence in a column:
c     at this stage, qvten(i,k,j) only includes horizontal advection.
c     sca: is the amount of total moisture convergence
c
      do 100 i=2,ilxm
c
       sca = 0.0
       do 20 k=1,kx
        sca=sca+qvten(i,k,j)*dsigma(k)
20     continue
c
c-----determine if moist convection exists:
c
       if (sca .lt. qdcrit) go to 70
c
c-----check for stability
c
c--1--compute eqt (equivalent potential temperature)
c     between surface and 700 mb, with perturbation temperature
c     and moisture added. the maximum eqt will be regarded
c     as the origin of air parcel that produce cloud.
c
       eqtm = 0.0
       do 210 k = k700, kx
        ttp = ta(i,k,j)/psa(i,j)+pert
        q = qva(i,k,j)/psa(i,j)+perq
        psg = psa(i,j)*a(k)+ptop
        t1 = ttp*(100./psg)**rovcp
        eqt = t1*dexp(xlvocp*q/ttp)
        if (eqt .gt. eqtm) then
         eqtm = eqt
         tmax = ttp
         qmax = q
         pmax = psg
        end if
210    continue
c
c--2--compute lcl, get the sigma and p of lcl
c
       emax = qmax*pmax/(ep2+qmax)
       tdmax = 5418.12/(19.84659-dlog(emax/.611))
       dalr = g/cp
       dplr = (g*tdmax*tdmax)/(ep2*xlv*tmax)
       zlcl = (tmax-tdmax)/(dalr-dplr)
       tlcl = tmax-dalr*zlcl
       tmean = 0.5*(tmax+tlcl)
       dlnp = (g*zlcl)/(r*tmean)
       plcl = pmax*dexp(-dlnp)
       siglcl = (plcl-ptop)/psa(i,j)
c
c--3--compute seqt (saturation equivalent potential temperature)
c     of all the levels that are above the lcl
c
       do 220 k = 1,kx
        if (a(k) .ge. siglcl) go to 230
220    continue
230    continue
       kbase = k
       if (kbase .gt. kx) kbase = kx
c
c.....kbase is the layer where lcl is located.
c
       do 240 k = 1,kbase
        ttp = ta(i,k,j)/psa(i,j)
        psg = psa(i,j)*a(k)+ptop
        es = .611*dexp(19.84659-5418.12/ttp)
        qs = ep2*es/(psg-es)
        t1 = ttp*(100./psg)**rovcp
        seqt(k) = t1*dexp(xlvocp*qs/ttp)
240    continue
c
c--4--when seqt = eqt + dt, cloud top is reached.
c     eqt is the eqt of cloud (same as lcl eqt).
c
       do 260 kk = 1,kbase
        k = kbase+1-kk
        deqt = seqt(k)-eqtm
        if (deqt .gt. dlt) go to 270
260    continue
270    continue
c
c.....cloud top has been reached
c
       ktop = k
c
c--5--check cloud depth
c     if cloud depth is less than critical depth (cdscld = 0.3),
c     the convection is killed
c
       dsc = (siglcl-a(ktop))
       if (dsc .lt. cdscld) go to 70
c
c--6--check negative area
c     if negative area is larger than the positive area
c     convection is killed.
c
       ttsum = 0.
       do 280 k=ktop,kbase
        ttsum =  (eqtm-seqt(k))*dsigma(k) + ttsum
280    continue
       if (ttsum .lt. 0.) go to 70
c
c.....you are here if stability was found.
c
c.....if values dont already exist in array twght,vqflx for this
c     kbase/ktop, then flag it, and set kbase/ktop to standard
c
       if ((kbase .lt. 5) .or. (ktop .gt. kbase-3)) then
        print 300, ktau,jyear,i,j,kbase,ktop
300     format (/,' >>in **cupara**: at ktau=',i8,' in year=',i5,
     &          ' & (i,j)=(',i2,',',i2,'),   ',
     &          ' kbase/ktop are non-standard:',2i3,
     &          '  & will be set to closest standard.')
        if (kbase .lt. 5) kbase = 5
        if (ktop .gt. kbase-3) ktop = kbase-3
       end if
c
c.....convection exist, compute convective flux of water vapor and
c     latent heating
c     icon   : is a counter which keep track the total points where
c              deep convection occurs.
c     c301   : is the 'b' factor in kuo's scheme.
c
       icon(j)=icon(j)+1
       sum = 0.
       sumb = 0.
       arh = 0.
       psx=psa(i,j)
       do 310 k = 1,kx
        qwght(k)=0.0
310    continue
       do 30 k=ktop,kx
        pux=psx*a(k)+ptop
        e1=.611*dexp(19.84659-5418.12/(ta(i,k,j)/psx))
        qs=ep2*e1/(pux-e1)
        rh=qva(i,k,j)/(qs*psx)
        rh=dmin1(rh,1.d0)
        xsav=(1.0-rh)*qs
        qwght(k)=xsav
        sumb=sumb+qs*dsigma(k)
        arh=arh+rh*qs*dsigma(k)
        sum=sum+xsav*dsigma(k)
30     continue
       arh=arh/sumb
       c301=2.0*(1.0-arh)
       if (c301 .lt. 0.0) c301=0.0
       if (c301 .gt. 1.0) c301=1.0
       if (sum .le. 0.0) then
        c301=0.0
        sum=1.0
       end if
       do 50 k=ktop,kx
        qwght(k)=qwght(k)/sum
50     continue
       do 60 k=1,kx
        ttconv = xlvocp*(1.0-c301)*twght(k,kbase,ktop)*sca
        rsheat(i,k,j)=rsheat(i,k,j)+ttconv*dt/2.
cx      if (ttconv*2. .gt. 0.01) write(18,1234) i,j,k,ttconv*2.
c1234    format(1x,'cupara, i=',i4,' j=',i4,' k=',i4,' qteva=',e12.4)
        apcnt=(1.0-c301)*sca/4.3e-3
        eddyf = apcnt*vqflx(k,kbase,ktop)
        qvten(i,k,j)=eddyf
        rswat(i,k,j)=rswat(i,k,j)+c301*qwght(k)*sca*dt/2.
60     continue
c
c  find cloud fractional cover and liquid water content
c
       kbaseb=min0(kbase,kx-2)
       if (ktop .le. kbaseb) then
        kclth=kbaseb-ktop+1
        akclth=1./dble(kclth)
        do 65 k=ktop,kbaseb
         cldlwc(i,k)=cllwcv
         cldfra(i,k)=1.-(1.-clfrcv)**akclth
65      continue
       end if
c.....the unit for rainfall is mm.
       prainx=(1.-c301)*sca*dtmin*60000./g
       rainc(i,j)=rainc(i,j)+prainx
c instantaneous precipitation rate for use in bats (mm/s)
       aprdiv=dble(nbatst)
       if (jyear.eq.jyear0.and.ktau.eq.0) aprdiv=1.
       pptc(i,j)=pptc(i,j)+prainx/(dtmin*60.)/aprdiv
c
Cchem2
       if(ichem.eq.1) then
c  before go to 100 put
         icumtop(i,j) = ktop
         icumbot(i,j) = kbaseb
       end if
Cchem2_

       go to 100
c
c.....convection not exist, compute the vertical advection term:
c
70     continue
       tmp3(i,1) = 0.
       do 80 k=2,kx
        if(qva(i,k,j).lt.1.e-15) then
         tmp3(i,k)=0.0
        else
         tmp3(i,k)=qva(i,k,j)*(qva(i,k-1,j)/qva(i,k,j))**qcon(k)
        end if
80     continue
       qvten(i,1,j)=qvten(i,1,j)-qdot(i,2,j)*tmp3(i,2)/dsigma(1)
       do 90 k=2,kxm
        qvten(i,k,j)=qvten(i,k,j)-(qdot(i,k+1,j)*tmp3(i,k+1)
     &                            -qdot(i,k,j)*tmp3(i,k))/dsigma(k)
90     continue
       qvten(i,kx,j)=qvten(i,kx,j)+qdot(i,kx,j)*tmp3(i,kx)/dsigma(kx)
c
100   continue         !end i=2,ilxm loop
c
      do 105 k=1,kx
      do 105 i=2,ilxm
       rsheat(i,k,j)=dmax1(rsheat(i,k,j),0.d0)
       rswat(i,k,j)=dmax1(rswat(i,k,j),0.d0)
       rsht=rsheat(i,k,j)/tauht
       rswt=rswat(i,k,j)/tauht
       tten(i,k,j)=tten(i,k,j)+rsht
       qvten(i,k,j)=qvten(i,k,j)+rswt
       rsheat(i,k,j)=rsheat(i,k,j)*(1.-dt/(2.*tauht))
       rswat(i,k,j)=rswat(i,k,j)*(1.-dt/(2.*tauht))
Cbxq   if(rsht/psb(i,j).gt..0002)write(18,1222)ktau,jyear,i,j,k,rsht/psb(i,j)
C1222   format (1x,'ktau= ',i7,' jyear= ',i5,' i= ',i5,' j= ',i5,
C     1        ' k= ',i5,' ttconv =',e15.7)
105   continue
c
      return
      end subroutine cupara
