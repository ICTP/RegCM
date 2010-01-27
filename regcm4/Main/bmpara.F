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

      subroutine bmpara(tten, qten, j)
C
C    *** modified by jack kain of penn state to replace the look-up tabl
C        by calculations.
C
C     ******************************************************************
C     *                                                                *
C     *  convective adjustment for deep or shallow convection          *
C     *                                                                *
C     *  references:                                                   *
C     *                                                                *
C     *  betts, a.k., 1986:  a new convective adjustment scheme.       *
C     *    part i: observational and theoretical basis.  quart. j. r.  *
C     *    met. soc., 112, 677-691.                                    *
C     *                                                                *
C     *  betts, a.k., and m.j. miller, 1986:  a new convective         *
C     *    adjustment scheme.  part ii: single column tests using      *
C     *    gate wave, bomex, atex and arctic air mass data sets.       *
C     *    quart. j. r. met. soc., 112, 693-709.                       *
C     *                                                                *
C     *  n.b.  part of the code is scalar.  in global models           *
C     *  convection occurs in less than 30/100 points.  with           *
C     *  simulataneous vector processing for both deep and shallow     *
C     *  convection, there would be a lot of redundant vector          *
C     *  computations.  if vector processing is 10 times faster        *
C     *  than scalar, one might hope that the cpu time will be about   *
C     *  the same for both scalar and vector code.                     *
C     *                                                                *
C     ******************************************************************
C *** warning: this subroutine will not work if kx.lt.12;
C
      implicit none
#     include "parameter.inc"

#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"

#     include "include/main.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/bats.cb2"
#     include "include/trachem.cb"

      integer j
      real(kind=8)  aliq,bliq,cliq,dliq,aice,bice,cice1,dice,xls0,xls1
      common/vapprs/aliq,bliq,cliq,dliq,aice,bice,cice1,dice,xls0,xls1
      logical unis,unil,oct90
      real(kind=8)  lblev
      real(kind=8)  h1,h1d5,h2d5,h3000,h10e5,d00,d125,d50,d608
      real(kind=8)  dm2859,elwv,elivw,row,epsq,a2,a3,a4,t0,t1,pq0,stresh
      real(kind=8)  stabs,stabd,stabfc,dttop,rhf,pmn,epsup,epsdn,epsth
      real(kind=8)  pbm,pqm,pno,pone,psh,pfrz,pshu,zno,zsh,fss,efimn
      real(kind=8)  efmnt,fcc1,fcp,dspbfl,dsp0fl,dsptfl,fsl,dspbfs
      real(kind=8)  dsp0fs
      real(kind=8)  dsptfs,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss
      real(kind=8)  epsntp,efifc,avgefi,dspc,epsp,stefi,slopbl,slop0l
      real(kind=8)  sloptl,slopbs,slop0s,slopts,slope,a23m4l,d273,cporng
      real(kind=8)  elocp,cprlg,rcp
      integer list,la,ksmud,nrow,itb,jtb,lp1,lm1,lscrch,l1,l2,l3,l4
C
C Very long parameter list. Stay with us.
C
      parameter (
     &  h1      = 1.e0,
     &  h1d5    = 1.5e0,
     &  h2d5    = 2.5e0,
     &  h3000   = 3000.e0,
     &  h10e5   = 100000.e0,
     &  d00     = 0.e0,
     &  d125    = 0.125e0,
     &  d50     = 0.5e0,
     &  d608    = 0.608e0,
     &  dm2859  = -287.04/1004.e0,
     &  elwv    = 2.50e6,
     &  elivw   = 2.72e6,
     &  row     = 1.e3,
     &  epsq    = 2.e-12,
     &  a2      = 17.2693882e0,
     &  a3      = 273.16e0,
     &  a4      = 35.86e0,
     &  t0      = 273.16e0,
     &  t1      = 274.16e0,
     &  pq0     = 379.90516e0,
     &  stresh  = 1.10e0,
     &  stabs   = 1.0e0,
     &  stabd   = 0.90e0,
     &  stabfc  = 1.00e0,
     &  dttop   = 1.0e0,
     &  list    = 3,
C RH drop 10%
     &  rhf     = 0.20,
     &  pmn     = 6500.0,
     &  epsup   = 1.00, ! 0.93
     &  epsdn   = 1.05, ! 1.00
     &  epsth   = 6.0,  ! rh drop 10%
     &  pbm     = 30000.0,
     &  pqm     = 20000.0,
     &  pno     = 9000.0,
     &  pone    = 2500.0,
     &  psh     = 29000.0,
     &  pfrz    = 15000.0,
     &  pshu    = 45000.0,
     &  zno     = 750.0,
     &  zsh     = 3999.0,
     &  lblev   = 0.98,
C
C    august `91: scheme having the option of using different fast and
C    slow profiles for sea and for land points; and also the "sea" and
C    the "land" scheme everywhere.  over land profiles depart from the
C    fast (dry) profiles only for precipitation/time step .gt.
C    a prescribed value (currently, in the version #3 done wednesday
C    18 september, 1/4 inch/24 h).  use of various switches as follows.
C
C       the "old" ("hard", =zavisa oct.1990) land scheme with fixed
C       land profiles is run by
C       * setting oct90=.true. in the first executable line following
C            these comments (this reactivates efi=h1 over land if
C            .not.unis, and cldefi(k)=efimn as the leftover cldefi
C            value at swap points);
C       * defining fast land profiles same as fast sea profiles (or by
C            choosing another set of land profiles zavisa used at
C            earlier time);
C       * setting fsl=1.; and
C       * defining stefi (starting efi) equal to avgefi.
C            (the last three points are handled by switching around the
C            "cfm" comments at two places)
C
C       the "old,old" (appr. original betts) scheme is run by
C       * specifying unil=.true.;
C       * setting fsl=1.;
C       * setting oct90=.true.
C            (with these settings fast land profiles only are used).
C
     &  unis    = .false.,
     &  unil    = .true.,
     &  oct90   = .true.,
     &  fss     = 0.60e0,
     &  efimn   = 0.20e0,
     &  efmnt   = 0.70e0,
     &  fcc1    = 0.50,
     &  fcp     = h1-fcc1,
C
C   In this version 3.5, over land and for the fast profiles, dspb
C   is prescribed to be 25 percent drier than the fast sea value
C   (in rough agreement with binder, qj, in press) while dsp0 and
C   dspt are each 20 percent drier than the corresponding fast
C   sea values. With fsl=.875 this makes the average of the fast
C   and the slow land profiles somewhat drier than the oct90 fixed
C   land profiles.
C
C   Values for oct90 = .false. are commented in next lines
C
     &  dspbfl  = -3875.e0, ! -4843.75e0
     &  dsp0fl  = -5875.e0, ! -7050.00e0
     &  dsptfl  = -1875.e0, ! -2250.0e0
     &  fsl     = 1.0e0,    ! .900e0
     &  dspbfs  = -3875.e0,
     &  dsp0fs  = -5875.e0,
     &  dsptfs  = -1875.e0,
     &  dspbsl  = dspbfl*fsl,
     &  dsp0sl  = dsp0fl*fsl,
     &  dsptsl  = dsptfl*fsl,
     &  dspbss  = dspbfs*fss,
     &  dsp0ss  = dsp0fs*fss,
     &  dsptss  = dsptfs*fss,
C*** new convection scheme with crossing dsp profiles ******************
C+-  &, unis=.false.,efimn=.71e0,efmnt=.71,fcc1=0.5,fcp=h1-fcc1
C+-  &, dspbl=-3875.e0,dsp0l=-5875.e0,dsptl=-1875.e0
C+-  &, dspbs=-2875.e0,dsp0s=-5125.e0,dspts=-4875.e0
C+-  &, dspbf=-4375.e0,dsp0f=-4375.e0,dsptf=-1000.e0
C*** betts convection scheme *******************************************
C    &, unis=.false.,efimn=.9999e0,efmnt=.9999e0,fcc1=.50,fcp=h1-fcc1
C    &, dspbl=-3875.e0,dsp0l=-5875.e0,dsptl=-1875.e0
C    &, dspbf=-3875.e0,dsp0f=-5875.e0,dsptf=-1875.e0
C    &, dspbs=-3875.e0,dsp0s=-5875.e0,dspts=-1875.e0
C***********************************************************************
     &  epsntp  = 0.0010e0,
     &  efifc   = 5.0e0,
     &  avgefi  = (efimn+1.e0)*.5e0,
     &  dspc    = -3000.e0,
     &  epsp    = 1.e-7,
C    &  stefi   = 1.e0,
C*** activate the following line and comment out the preceding line if
C*** oct90=.true.
     &  stefi   = avgefi,
     &  slopbl  = (dspbfl-dspbsl)/(h1-efimn),
     &  slop0l  = (dsp0fl-dsp0sl)/(h1-efimn),
     &  sloptl  = (dsptfl-dsptsl)/(h1-efimn),
     &  slopbs  = (dspbfs-dspbss)/(h1-efimn),
     &  slop0s  = (dsp0fs-dsp0ss)/(h1-efimn),
     &  slopts  = (dsptfs-dsptss)/(h1-efimn),
     &  slope   = (h1   -efmnt)/(h1-efimn)
     &          )
      parameter ( 
     &  a23m4l = a2*(a3-a4)*elwv,
     &  d273   = 1./t0,
     &  cporng = 1./dm2859,
     &  elocp  = elwv/1004.,
     &  cprlg  = 1004./(row*9.81*elwv),
     &  rcp    = h1/1004.
     &          )
C
C  Instability for too large lsh
C
      parameter (
     &  la    = 13,
     &  ksmud = 0,
     &  nrow  = 0,
     &  itb   = 100,
     &  jtb   = 150
     &          )
      parameter (
     &  lp1    = kx+1,
     &  lm1    = kx-1,
     &  lscrch = 4*kx+1+la+1,
     &  l1     = la+1,
     &  l2     = la+kx+1,
     &  l3     = la+2*kx+1,
     &  l4     = la+3*kx+1
     &          )
C----------------------------------------------------------------------
      real(kind=8)  pl,thl,rdq,rdth,rdp,rdthe,qs0,sqs
      real(kind=8)  the0,sthe,ptbl,ttbl
      common /tbls/ pl,thl,rdq,rdth,rdp,rdthe,qs0(jtb),sqs(jtb),
     & the0(jtb),sthe(jtb),ptbl(itb,jtb),ttbl(jtb,itb)
C----------------------------------------------------------------------
      real(kind=8)  tten(ix,kx),qten(ix,kx)
      integer kdp(kx)
      real(kind=8)  sm(ix)
C
      real(kind=8)  t(ix,kx), q(ix,kx) ,  z0(ix,kx)
     &       ,dzq(kx), cldhgt(ix), pdp(kx)
C-----------------------------------------------------------------------
      real(kind=8)
     & trefk (kx),qrefk (kx),pk    (kx),apek  (kx),tk    (kx)
     &,thsk  (kx),psk   (kx),apesk (kx),qk    (kx),therk (kx)
     &,thvref(kx),thevrf(kx),thvmod(kx),dift  (kx),difq  (kx)
     &,qsatk (kx),fpk   (kx)
      integer
     & ntopd (kx),nbotd (kx),ntops (kx),nbots (kx)
     &,ndpthd(kx),ndpths(kx)
C
      integer 
     & ltop  (ix),lbot  (ix),kshal (ix),kdeep (ix)
     &,ip300 (ix),ifbuoy(ix),kf(ix),ml(ix)
      real(kind=8)
     & prtop (ix),pbot  (ix)
     &,tthbt (ix)
     &,tthes (ix)
     &,psp   (ix),thbt  (ix),thesp (ix)
     &,p     (ix)
     &,dspb  (ix),dsp0  (ix),dspt  (ix)
     &,ape   (ix,kx)
     &,tref  (ix,kx),tmod  (ix,kx),qqmod  (ix,kx),ths(ix)
      real(kind=8)  tauk,cthrs,cell,tv0,pdiff,pkl,psfck,ee,tdpt,tlcl
      real(kind=8)  ak,es,qs,efi,dentpy,avrgt,preck,dspbk,dsp0k,dsptk
      real(kind=8)  tkl,qkl,apekl,pkb,pkt,pk0,stabdl,rdp0t,dthem,dsp
      real(kind=8)  sumde,sumdp,hcorr,tskl,dhdt,thskl,diftl,difql
      real(kind=8)  avrgtl,ztop,pthrs,drheat,fefi,ptpk,pbtk,thvmkl
      real(kind=8)  pflag,pdiffk,rhdmx,rhl,rhh,thtpk,dpmix,smix,dpkts
      real(kind=8)  thivi,summm,summd,sumdt,rdpsum,tcorr,trfkl,psum
      real(kind=8)  qsum,potsum,qotsum,otsum,dst,fptk,dpkl,rtbar
      real(kind=8)  rotsum,dstq,den,dqref,qrftp,qrfkl,qnew,dtdeta
      integer i,lb,k,kb,ivi,l,khdeep,nswap,n,ltpk,lbtk,ltp1,lbm1
      integer l0,l0m1,iter,lcor,lqm,ndepth,nneg,khshal,ndstn
      integer ndstp,kk,ltsh,lshu,ltp2,nshal,negds,ndeep
      integer icond,iconss
      real(kind=8)  pppk
      integer kbaseb,kclth
      real(kind=8)  prainx,aprdiv,akclth
      real(kind=8)  qu,tpfc
      external tpfc
C-----------------------------------------------------------------------
      do k=1,kx
      do i=1,ilx
       cldlwc(i,k)=0.
       cldfra(i,k)=0.
      end do
      end do
Cchem2
      if(ichem.eq.1) then
C
C     icumtop = top level of cumulus clouds
C     icumbot = bottom level of cumulus clouds
C           (calculated in cupara and stored for tractend)
C  before do 100 put
      do i=2,ilxm
         icumtop(i,j) = 0
         icumbot(i,j) = 0
      end do
      end if
Cchem2_
      icond=0
      iconss=0
C     write(3,12345)
C12345 format(' ---> enter bmpara')
      tauk=dt2/trel
      cthrs=(0.006350/86400.)*dt2/cprlg
C-----------------------------------------------------------------------
C
C... sm is surface mask: =1 water; =0 land
      do 90 i=2,ilxm
        if(veg2d(i,j).ge.0.002)then
          sm(i)=0.
        else
          sm(i)=1.
        end if
   90 continue
      if(jyear.eq.jyear0.and.ktau.eq.0)then
        do 100 i=2,ilxm
          cldefi(i,j)=avgefi*sm(i)+stefi*(h1-sm(i))
  100   continue
C       print 103,zno,zsh
C       print 104,zsh
C103     format(//,'shallow convection occurs when cloud depth is',
C     1      ' .ge. ',f5.0,' m. and .lt. ',f5.0,' m.')
C104     format(//,'deep convection occurs when cloud depth is .ge. ',
C     1      f5.0,' km.')
      end if
C...lb is currently set to kx-1
      lb=kx-1
Cvvvvvvvvvvvvv diagnostics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      do 110 k=1,kx
        ntopd(k)=0
        nbotd(k)=0
        ntops(k)=0
        nbots(k)=0
        ndpths(k)=0
        ndpthd(k)=0
  110 continue
Caaaaaaaaaaaaa diagnostics aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C--------------preparations---------------------------------------------
C...find melting level...
      do 115 i=2,ilxm
  115 ml(i)=kxp1
      do 120 i=2,ilxm
C        if(kf(i).eq.1)goto 120
        do 121 k=1,kx
          t(i,k)=tb(i,k,j)/psb(i,j)
          if(t(i,k).gt.t0.and.ml(i).eq.kxp1)ml(i)=k
          q(i,k)=qvb(i,k,j)/psb(i,j)
          pppk=(a(k)*psb(i,j)+ptop)*1000.
          ape(i,k)=(pppk/h10e5)**dm2859
  121   continue
        lbot(i)=kx
        thesp(i)=d00
        thbt(i)=d00
        psp(i)=9.5e4
C      pdsl (k)=res(k)*pd(k)
        tref(i,1)=t(i,1)
C...ifbuoy = 0 means no positive buoyancy; ifbuoy(i) means yes...
C...ip300 is the highest model level in the lowest 300 mb...
        ifbuoy(i)=0
        ip300(i)=0
        kf(i)=0
C...activate this statement to run only shallow betts-miller with kf and
C...ak deep convection...
C      kf(i)=cvmgt(1,0,(icps.eq.2.and.nca(i,j).gt.0) .or.
C     *                (icps.eq.1.and.nak(i,j).gt.0))
        cell=ptop/psb(i,j)
        do 378 k=1,kx
          tv0=t(i,k)*(1.+d608*q(i,k))
          dzq(k)=r/g*tbase(i,k,j)*
     +           dlog((sigma(k+1)+cell)/(sigma(k)+cell))
C        dzq(k)=(r/g)*tv0*dlog((sigma(k+1)+cell)/(sigma(k)+cell))
  378   continue
        z0(i,kx)=0.5*dzq(kx)
        do 379 k=kx-1,1,-1
          z0(i,k)=z0(i,k+1)+0.5*(dzq(k)+dzq(k+1))
  379   continue
  120 continue
C--------------padding specific humidity if too small-------------------
      do 130 k=1,kx
        do 130 i=2,ilxm
C        if(kf(i).eq.1)goto 130
          if(q(i,k).lt.epsq)q(i,k)=epsq
          pdiff=(1.-a(k))*psb(i,j)
          if(pdiff.lt.30..and.ip300(i).eq.0)ip300(i)=k
  130   continue
C--------------search for maximum buoyancy level------------------------
      do 140 kb=1,kx
C       print *,'i, j, kb, z0, pbl =',
C     *i,j,kb,z0(i,kb),pbl(i,j)
C--------------trial maximum buoyancy level variables-------------------
        do 150 i=2,ilxm
C           if(z0(i,kb).gt.pbl(i,j)+10.)goto 150
C           if(kf(i).eq.1)goto 150
          pkl=(a(kb)*psb(i,j)+ptop)*1000.
          psfck=(a(kx)*psb(i,j)+ptop)*1000.
          if(pkl.lt.psfck-pbm)goto 150
          tthbt(i)=t(i,kb)*ape(i,kb)
          ee=pkl*q(i,kb)/(0.622+q(i,kb))
          tdpt=1./(d273-rv/elwv*dlog(ee/611.))
          tdpt=dmin1(tdpt,t(i,kb))
          tlcl=tdpt-(.212+1.571e-3*(tdpt-t0)-4.36e-4*(t(i,kb)-t0))*(
     +         t(i,kb)-tdpt)
          tthes(i)=tthbt(i)*exp(elocp*q(i,kb)/tlcl)
C--------------check for maximum buoyancy-------------------------------
          if(tthes(i).gt.thesp(i))then
            psp(i)=h10e5*(tthbt(i)/tlcl)**cporng
            thbt(i)=tthbt(i)
            thesp(i)=tthes(i)
          end if
C-----------------------------------------------------------------------
  150   continue
  140 continue
C---------choose cloud base as model level just below psp--------------
      do 240 k=1,lm1
        ak=a(k)
        do 240 i=2,ilxm
C          if(kf(i).eq.1)goto 240
          p(i)=(ak*psb(i,j)+ptop)*1000.
C cloud bottom cannot be above 200 mb
          if(p(i).lt.psp(i).and.p(i).ge.pqm)lbot(i)=k+1
  240   continue
C*** warning: lbot must not be gt kx-1 in shallow convection
C*** make sure the cloud base is at least 25 mb above the surface
      do 250 i=2,ilxm
C       if(kf(i).eq.1)goto 250
        pbot(i)=(a(lbot(i))*psb(i,j)+ptop)*1000.
        psfck=(a(kx)*psb(i,j)+ptop)*1000.
        if(pbot(i).ge.psfck-pone.or.lbot(i).ge.kx)then
C*** cloud bottom is at the surface so recalculate cloud bottom
          do 245 k=1,lm1
            p(i)=(a(kx)*psb(i,j)+ptop)*1000.
            if(p(i).lt.psfck-pone)lbot(i)=k
  245     continue
          pbot(i)=(a(lbot(i))*psb(i,j)+ptop)*1000.
        end if
  250 continue
C--------------cloud top computation------------------------------------
      do 260 i=2,ilxm
        prtop(i)=pbot(i)
  260 ltop(i)=lbot(i)
      do 270 ivi=1,kx
        l=lp1-ivi
C--------------find environmental saturation equiv pot temp...
        do 280 i=2,ilxm
          p(i)=(a(l)*psb(i,j)+ptop)*1000.
          es=aliq*exp((bliq*t(i,l)-cliq)/(t(i,l)-dliq))
          qs=0.622*es/(p(i)-es)
          ths(i)=t(i,l)*ape(i,l)*exp(elocp*qs/t(i,l))
  280   continue
C--------------buoyancy check-------------------------------------------
        do 350 i=2,ilxm
          if(l.le.lbot(i))then
            if(thesp(i).gt.ths(i))ifbuoy(i)=1
            if(thesp(i).gt.ths(i)-1.5.and.ifbuoy(i).eq.1)ltop(i)=l+1
          end if
  350   continue
C------------------------------------------------
  270 continue
C--------------cloud top pressure---------------------------------------
      do 275 i=2,ilxm
C        if(kf(i).eq.1) goto 275
        prtop(i)=(a(ltop(i))*psb(i,j)+ptop)*1000.
  275 continue
C-----------------------------------------------------------------------
C--------------define and smooth dsps and cldefi------------------------
Cvvvvv unified or separate land/sea conv vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      if(unis)then
        do 370 i=2,ilxm
C           if(kf(i).eq.1)goto 370
          efi=cldefi(i,j)
          dspb(i)=(efi-efimn)*slopbs+dspbss
          dsp0(i)=(efi-efimn)*slop0s+dsp0ss
          dspt(i)=(efi-efimn)*slopts+dsptss
  370   continue
      else if(.not.unil)then
        do 371 i=2,ilxm
C           if(kf(i).eq.1)goto 371
          efi=cldefi(i,j)
          dspb(i)=((efi-efimn)*slopbs+dspbss)*sm(i)+((efi-efimn)*slopbl+
     +            dspbsl)*(h1-sm(i))
          dsp0(i)=((efi-efimn)*slop0s+dsp0ss)*sm(i)+((efi-efimn)*slop0l+
     +            dsp0sl)*(h1-sm(i))
          dspt(i)=((efi-efimn)*slopts+dsptss)*sm(i)+((efi-efimn)*sloptl+
     +            dsptsl)*(h1-sm(i))
  371   continue
      else
        do 375 i=2,ilxm
          efi=cldefi(i,j)
          dspb(i)=((efi-efimn)*slopbl+dspbsl)
          dsp0(i)=((efi-efimn)*slop0l+dsp0sl)
  375   dspt(i)=((efi-efimn)*sloptl+dsptsl)
      end if
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C--------------extending sea structures inland along coastline----------
C      if(nrow.gt.0.and..not.unis.and..not.unil)    then
C-----------------------------------------------------------------------
C          do 373 i=2,ilxm
C      wf1(i)=0.
C      wf2(i)=0.
C      wf3(i)=0.
C 373  wf4(i)=0.
C-----------------------------------------------------------------------
C      krow=nrow
C                             do 372 kount=1,krow
C          do 374 i=2,ilxm
C      wf1(i)=(dspb  (i+ksw)+dspb  (i+kse)
C     2       +dspb  (i+knw)+dspb  (i+kne)+dspb  (i)*4.0)*hbm2(i)*0.125
C      wf2(i)=(dsp0  (i+ksw)+dsp0  (i+kse)
C     2       +dsp0  (i+knw)+dsp0  (i+kne)+dsp0  (i)*4.0)*hbm2(i)*0.125
C      wf3(i)=(dspt  (i+ksw)+dspt  (i+kse)
C     2       +dspt  (i+knw)+dspt  (i+kne)+dspt  (i)*4.0)*hbm2(i)*0.125
C 374  wf4(i)=(cldefi(i+ksw)+cldefi(i+kse)
C     2       +cldefi(i+knw)+cldefi(i+kne)+cldefi(i)*4.0)*hbm2(i)*0.125
C          do 372 i=2,ilxm
C      smk =sm(i,j)
C      rsmk=1.0-smk
C      dspb  (i)=dspb  (i)*smk+wf1(i)*rsmk
C      dsp0  (i)=dsp0  (i)*smk+wf2(i)*rsmk
C      dspt  (i)=dspt  (i)*smk+wf3(i)*rsmk
C 372  cldefi(i)=cldefi(i)*smk+wf4(i)*rsmk
C-----------------------------------------------------------------------
C      end if
C--------------initialize changes of t and q due to convection----------
      do 360 k=1,kx
        do 360 i=2,ilxm
          tmod(i,k)=d00
  360   qqmod(i,k)=d00
C--------------clean up and gather deep convection points---------------
      khdeep=0
      nswap=0
      do 380 i=2,ilxm
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     llmh=lmh(k)
C     ubar=(u(k+khs,llmh)+u(k-1  ,llmh)+u(k    ,llmh)+u(k+khn,llmh))*.25
C     vbar=(v(k+khs,llmh)+v(k-1  ,llmh)+v(k    ,llmh)+v(k+khn,llmh))*.25
C     vcnv=sqrt(ubar*ubar+vbar*vbar)
C     if(vcnv.lt.8.            )    ltop(k)=lbot(k)-lsh
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C...activate this statement to run only shallow betts-miller with kf and
C...ak deep convection...
C      if(ltop(i).gt.lbot(i) .or. kf(i).eq.1)  then
        if(ltop(i).gt.lbot(i))then
          ltop(i)=lbot(i)
          prtop(i)=pbot(i)
        end if
C     if(hbm2(i).lt.0.90       )    ltop(i)=lbot(i)
        cldhgt(i)=z0(i,ltop(i))-z0(i,lbot(i))
C  cloud is less than 90 mb deep or less than 3 sigma layers deep
C     if(prtop(i).gt.pbot(i)-pno.or.ltop(i).gt.lbot(i)-2)
C     if(cldhgt(i).lt.zno.or.ltop(i).gt.lbot(i)-2)
        if(cldhgt(i).lt.zno)cldefi(i,j)=avgefi*sm(i)+stefi*(h1-sm(i))
C  cloud has to be at least 290 mb deep
        if(cldhgt(i).ge.zsh)then
C     if(prtop(i).lt.pbot(i)-psh)then
          khdeep=khdeep+1
          kdeep(khdeep)=i
        end if
  380 continue
C************* horizontal loop for deep convection *********************
C***********************************************************************
C     print*,'number of deep convective points=',khdeep
      do 400 n=1,khdeep
        i=kdeep(n)
        dentpy=d00
        avrgt=d00
        preck=d00
C      if(icps.eq.2 .or. icps.eq.1)goto 532
        ltpk=ltop(i)
        lbtk=lbot(i)
Cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
Cdcdcdcdcdcdc  deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
Cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
C     lb   =kx-1
        efi=cldefi(i,j)
        dspbk=dspb(i)
        dsp0k=dsp0(i)
        dsptk=dspt(i)
C--------------initialize variables in the convective column------------
        do 410 k=1,kx
          dift(k)=d00
          difq(k)=d00
          tkl=t(i,k)
          tk(k)=tkl
          trefk(k)=tkl
          qkl=q(i,k)
          qk(k)=qkl
          qrefk(k)=qkl
          pkl=(a(k)*psb(i,j)+ptop)*1000.
***************
          tref(i,k)=tpfc(pkl,thesp(i),t(i,k),d273,elwv,qu,ape(i,k))
****************
          pk(k)=pkl
          psk(k)=pkl
          apekl=ape(i,k)
          apek(k)=apekl
          therk(k)=tref(i,k)*apekl
  410   thvmod(k)=(qkl*d608+h1)*tkl*apekl
C     print *, 'lb,pk(k) = ', lb,pk
C--------------deep convection reference temperature profile------------
        ltp1=ltpk+1
        lbm1=lb-1
        pkb=pk(lb)
        pkt=pk(ltpk)
C--------------temperature reference profile below freezing level-------
        l0=lb
        pk0=pk(lb)
        do 420 l=ltpk,lbm1
          ivi=ltpk+lbm1-l
          if(trefk(ivi+1).le.t1)goto 430
          stabdl=stabd
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C         if(ivi.eq.lbm1.and.sm(i).gt.d00) stabdl=stabd*stabfc
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          trefk(ivi)=((therk(ivi)-therk(ivi+1))*stabdl+trefk(ivi+1)*
     +               apek(ivi+1))/apek(ivi)
          l0=ivi
          pk0=pk(l0)
  420   continue
C--------------freezing level at or above the cloud top-----------------
        l0m1=l0-1
        goto 440
C--------------temperature reference profile above freezing level-------
  430   l0m1=l0-1
        rdp0t=h1/(pk0-pkt)
        dthem=therk(l0)-trefk(l0)*apek(l0)
        do 450 l=ltpk,l0m1
  450   trefk(l)=(therk(l)-(pk(l)-pkt)*dthem*rdp0t)/apek(l)
C--------------deep convection reference humidity profile---------------
  440   do 460 l=ltpk,lb
C--------------saturation pressure difference---------------------------
          if(pkb-pk0.ge.pfrz)then
! p(kx-1) -p(frz) .ge. 150 mb)
            if(l.lt.l0)then
              dsp=((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
            else
              dsp=((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
            end if
          else
            dsp=dspc
          end if
C--------------humidity profile-----------------------------------------
          if(pk(l).gt.pqm)then
! pressure must be below 200 mb
            psk(l)=pk(l)+dsp
            apesk(l)=(psk(l)/h10e5)**dm2859
            thsk(l)=trefk(l)*apek(l)
            qrefk(l)=pq0/psk(l)*
     +               exp(a2*(thsk(l)-a3*apesk(l))/(thsk(l)-a4*apesk(l)))
          else
            qrefk(l)=q(i,l)
          end if
  460   continue
        do 470 l=ltpk,lb
  470   thevrf(l)=trefk(l)*apek(l)*(qrefk(l)*d608+h1)
C--------------enthalpy conservation integral--------------------------
        do 480 iter=1,2
C-----------------------------------------------------------------------
          sumde=d00
          sumdp=d00
          do 490 l=ltpk,lb
            sumde=((tk(l)-trefk(l))*cp+(qk(l)-qrefk(l))*elwv)*dsigma(l)+
     +            sumde
  490     sumdp=sumdp+dsigma(l)
          hcorr=sumde/(sumdp-dsigma(ltpk))
          lcor=ltpk+1
C--------------find lqm-------------------------------------------------
          do 495 l=1,lb
            if(pk(l).le.pqm)lqm=l
  495     continue
C      print 4950, lcor,lqm
C4950  format(1x,' === lcor,lqm = ',2i6)
C--------------above lqm correct temperature only-----------------------
          if(lcor.le.lqm)then
            do 500 l=lcor,lqm
  500       trefk(l)=trefk(l)+hcorr*rcp
            lcor=lqm+1
          end if
C--------------below lqm correct both temperature and moisture----------
          do 510 l=lcor,lb
            tskl=trefk(l)*apek(l)/apesk(l)
            dhdt=qrefk(l)*a23m4l/(tskl-a4)**2+cp
            trefk(l)=hcorr/dhdt+trefk(l)
            thskl=trefk(l)*apek(l)
            qrefk(l)=pq0/psk(l)*
     +               exp(a2*(thskl-a3*apesk(l))/(thskl-a4*apesk(l)))
  510     continue
C-----------------------------------------------------------------------
  480   continue
Cvvvvvvvvvvvvv diagnostics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        do 520 l=1,kx
  520   thvref(l)=trefk(l)*apek(l)*(qrefk(l)*d608+h1)
C--------------heating, moistening, precipitation-----------------------
        do 530 l=ltpk,lb
          tkl=tk(l)
          diftl=(trefk(l)-tkl)*tauk
          difql=(qrefk(l)-qk(l))*tauk
          avrgtl=(tkl+tkl+diftl)
          dentpy=(diftl*cp+difql*elwv)*dsigma(l)/avrgtl+dentpy
          avrgt=avrgtl*dsigma(l)+avrgt
          preck=dsigma(l)*diftl+preck
          dift(l)=diftl
  530   difq(l)=difql
        dentpy=dentpy+dentpy
        avrgt=avrgt/(sumdp+sumdp)
C-------------swap if entropy and/or precip .le. 0 ...------------------
C 532   if(dentpy.lt.epsntp .or. preck.le.d00 .or.
C     *   icps.eq.2)    then
C 532   if(dentpy.lt.epsntp.or.preck.le.d00)then
        if(dentpy.lt.epsntp.or.preck.le.d00)then
          if(oct90)then
            cldefi(i,j)=efimn
          else
            cldefi(i,j)=efimn*sm(i)+stefi*(h1-sm(i))
          end if
          ztop=z0(i,lbot(i))+zsh-0.000001
C         prtop(i)=pbot(i)-psh    ! ptop is now pbot - 290.
          do 535 l=1,lb
            if(z0(i,l).ge.ztop)ltop(i)=l+1
C           if(pk(l).le.prtop(i)) ltop(i)=l+1
  535     continue
          prtop(i)=pk(ltop(i))
C------------cloud must be at least 2 layers thick---------------------
          if(lbot(i)-ltop(i).lt.2)ltop(i)=lbot(i)-2
          prtop(i)=pk(ltop(i))
          cldhgt(i)=z0(i,ltop(i))-z0(i,lbot(i))
          nswap=nswap+1
          goto 400
        end if
C--------------... deep convection otherwise----------------------------
        icond=icond+1
C***  keep the land value of efi equal to 1 until precip surpasses
C***  a threshold value, currently set to 0.25 inches per 24 hrs
        pthrs=cthrs/psb(i,j)
        drheat=(preck*sm(i)+dmax1(epsp,preck-pthrs)*(h1-sm(i)))*cp/avrgt
        efi=efifc*dentpy/drheat
Cvvvvv unified or separate land/sea conv. vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        if(oct90)then
          if(unis)then
            efi=cldefi(i,j)*fcp+efi*fcc1
          else if(.not.unil)then
            efi=(cldefi(i,j)*fcp+efi*fcc1)*sm(i)+h1-sm(i)
          else
            efi=h1
          end if
        else
          efi=cldefi(i,j)*fcp+efi*fcc1
        end if
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        if(efi.gt.h1)efi=h1
        if(efi.lt.efimn)efi=efimn
        cldefi(i,j)=efi
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        fefi=efmnt+slope*(efi-efimn)
C     fefi=dmax1(efi,efmnt)
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        preck=preck*fefi
C--------------update precipitation, temperature & moisture-------------
C         prec  (i,j)=psb(i,j)*preck*cprlg+prec  (i,j)
C correct rainc so that there is no double counting of precipitation (gk
C         rainc(i,j)=((psb(i,j)*1000.*preck*cprlg)*100.)+rainc(i,j)
        prainx=0.5*((psb(i,j)*1000.*preck*cprlg)*100.)
        rainc(i,j)=prainx+rainc(i,j)
C.....................precipitation rate for bats (mm/s)
        aprdiv=dble(nbatst)
        if (jyear.eq.jyear0.and.ktau.eq.0) aprdiv=1.
        pptc(i,j)=pptc(i,j)+prainx/(dtmin*60.)/aprdiv
        do 550 l=ltpk,lb
          tmod(i,l)=dift(l)*fefi/dt2
  550   qqmod(i,l)=difq(l)*fefi/dt2
C      print *, 'after do 550; icond = ', icond
Cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
Cdcdcdcdcdcdc  end of deep convection  dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
Cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
C-----------------------------------------------------------------------
  400 continue
C     print *, 'after loop 400'
Cvvvvvvvvvvvvv diagnostics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ndeep=0
      do 570 i=2,ilxm
        ltpk=ltop(i)
        lbtk=lbot(i)
        ptpk=prtop(i)
        pbtk=pbot(i)
C     lb  =kx-1
        if(cldhgt(i).ge.zsh)then
C     if(ptpk.lt.pbtk-psh)then    ! psh=290 mb ndeep=ndeep+1
          ndeep=ndeep+1
          ndepth=lb-ltpk
          ntopd(ltpk)=ntopd(ltpk)+1
          nbotd(lb)=nbotd(lb)+1
          if(ndepth.gt.0)ndpthd(ndepth)=ndpthd(ndepth)+1
        end if
  570 continue
      nneg=khdeep-ndeep
C      print *,'no. of swapped points, final deep: ',nswap,ndeep
Caaaaaaaaaaaaa diagnostics aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C--------------gather shallow convection points-------------------------
      khshal=0
      ndstn=0
      ndstp=0
      do 560 i=2,ilxm
C  if cloud depth < 90 or < 2 sigma layers no shallow convection
        if(cldhgt(i).lt.zno.or.ltop(i).gt.lbot(i)-2)goto 560
C     if(cldhgt(i).lt.zno)go to 560
C     if(prtop(i).gt.pbot(i)-pno.or.ltop(i).gt.lbot(i)-2)go to 560
        if(cldhgt(i).lt.zsh)then
C     if(prtop(i).ge.pbot(i)-psh)then  ! 90 mb < cld depth < 290 mb
          khshal=khshal+1
          kshal(khshal)=i
        end if
  560 continue
C************* horizontal loop for shallow convection ******************
C***********************************************************************
C      print*,'at j =',j,' number of shallow convective points=',
C     *khshal
Cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
Cscscscscscsc  shallow convection  cscscscscscscscscscscscscscscscscscsc
Cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
C-----------------------------------------------------------------------
      do 600 n=1,khshal
        i=kshal(n)
C         pz0= (psb(i,j)+ptop)*1000.       ! not used!
        do 610 k=1,kx
          tkl=t(i,k)
          tk(k)=tkl
          trefk(k)=tkl
          qkl=q(i,k)
          qk(k)=qkl
          qrefk(k)=qkl
          qsatk(k)=qkl
          pkl=(a(k)*psb(i,j)+ptop)*1000.
          pk(k)=pkl
          apekl=ape(i,k)
Cvvvvvvvvvvvv choose the pressure function vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C         fpk  (l) =dlog(pkl)
C         fpk  (l) =pkl
C         fpk  (l) =-1./pkl
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          apek(k)=apekl
          thvmkl=tkl*apekl*(qkl*d608+h1)
          thvmod(k)=thvmkl
          thvref(k)=thvmkl
C... for kdp calculation
          pdp(k)=pk(k)-pmn
  610   continue
C
C...find kdp...kdp(k) is the model level closest to 65 mb (pmn) above k;
C...this is the depth over which relative humidity drop is measured to
C...estimate shallow cloud top... see do 545...
C
        do 106 kk=kx,1,-1
          pflag=abs(pk(kx)-pdp(kk))
          do 105 k=kx-1,1,-1
            pdiffk=abs(pk(k)-pdp(kk))
            if(pdiffk.lt.pflag)then
              pflag=pdiffk
              if(kk.eq.k)then
                kdp(kk)=k-1
              else
                kdp(kk)=k
              end if
            end if
  105     continue
          kdp(kk)=max(1,kdp(kk))
C         print *,'kk, kdp =',kk,kdp(kk)
  106   continue
C       if (kk.lt.1)
C    &     call fatal(__FILE__,__LINE__, 'KK NOT CONVERGED')
C--------------search for shallow cloud top-----------------------------
        lbtk=lbot(i)
        ltsh=lbtk
        lbm1=lbtk-1
        pbtk=pk(lbtk)
        ztop=z0(i,lbot(i))+zsh-0.000001
C         ptpk=pbtk-psh
C--------------cloud top is level just above pbtk-psh ------------------
        do 612 l=1,kx
C              if (pk(l).le.ptpk) ltpk=l-1
C              if (pk(l).le.ptpk) ltpk=l
          if(z0(i,l).ge.ztop)ltpk=l
  612   continue
        ptpk=pk(ltpk)
C--------------highest level allowed is level just below pshu-----------
        if(ptpk.le.pshu)then
          do 615 l=1,kx
            if(pk(l).le.pshu)lshu=l+1
  615     continue
          ltpk=lshu
          ptpk=pk(ltpk)
        end if
        ltp1=ltpk+1
        ltp2=ltpk+2
C--------------no shallow convection in unstable environment------------
C             do 541 l=ltp2,lbtk
C         if(thvmod(l).gt.thvmod(l-1))    then
C             ltop(k)=lbot(k)
C             go to 600
C         end if
C541      continue
C-----------------------------------------------------------------------
        do 540 l=ltpk,lbtk
C          es=cvmgt(aliq*exp((bliq*tk(l)-cliq)/(tk(l)-dliq)),
C     *   aice*exp((bice*tk(l)-cice1)/(tk(l)-dice)),l.ge.ml(i))
          if(l.ge.ml(i))then
            es=aliq*exp((bliq*tk(l)-cliq)/(tk(l)-dliq))
          else
            es=aice*exp((bice*tk(l)-cice1)/(tk(l)-dice))
          end if
          qsatk(l)=0.622*es/(pk(l)-es)
C          qsatk(l)=pq0/pk(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
Cvvvvvvvvvvvv choose the pressure function vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C         fpk  (l)=qsatk(l)
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
  540   continue
C-----------------------------------------------------------------------
C         rhh=qk(ltpk)/qsatk(ltpk)
C       print *,'before do 545; i, j, lbot, ltop =',i,j,
C     *lbot(i),ltop(i)
C       print *,'cldhgt =',cldhgt(i)
        rhdmx=rhf
        do 545 l=ltp1,lbm1
          rhl=qk(l)/qsatk(l)
C           if(kdp(l).lt.ltpk.or.kdp(l).gt.lbtk)then
C             qsatk(l)=pq0/pk(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
C           end if
          rhh=qk(kdp(l))/qsatk(kdp(l))
          if(rhh+rhf.lt.rhl)ltsh=l
C           if(rhh+rhf.lt.rhl)then
C             rhdmx=dmax1(rhl-rhh,rhdmx)
C             if(rhdmx.eq.rhl-rhh)ltsh=kdp(l)
C           end if
C       print *,'in do 545; l, rhl, kdp, rhh =',
C     *          l,rhl,kdp(l),rhh
  545   continue
C
        ltop(i)=ltsh
        prtop(i)=pk(ltsh)
        ltp1=ltsh
        ltpk=ltsh-1
        cldhgt(i)=z0(i,ltop(i))-z0(i,lbot(i))
C       print *,'after do 545; ltop, cldhgt =',ltop(i),cldhgt(i)
C       if(j.gt.2)
C    &     call fatal(__FILE__,__LINE__, 'KK NOT CONVERGED')
C-----------------------------------------------------------------------
C  if cloud is not at least 90 mb or 3 sigma layers deep, then no cloud
        if(cldhgt(i).lt.zno.or.ltop(i).gt.lbot(i)-2)then
C         if(cldhgt(i).lt.zno)then
C         if(prtop(i).gt.pbot(i)-pno.or.ltop(i).gt.lbot(i)-2)then
          ltop(i)=lbot(i)
          prtop(i)=pbot(i)
C          print *,'cldhgt(i).lt.zn0 .or. ltop.gt.lbot-2...'
          goto 600
        end if
C--------------scaling potential temperature & table index at top-------
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        thtpk=t(i,ltp1)*ape(i,ltp1)
C        thtpk=t(i,ltpk)*ape(i,ltpk)
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        pkl=(a(ltp1)*psb(i,j)+ptop)*1000.
        ee=pkl*q(i,ltp1)/(0.622+q(i,ltp1))
        tdpt=1./(d273-rv/elwv*dlog(ee/611.))
        tdpt=dmin1(tdpt,t(i,ltp1))
        tlcl=tdpt-(.212+1.571e-3*(tdpt-t0)-4.36e-4*(t(i,ltp1)-t0))*(
     +       t(i,ltp1)-tdpt)
        ptpk=h10e5*(thtpk/tlcl)**cporng
C         print *,'t, tdpt, tlcl, p, plcl =',t(i,ltp1),tdpt,tlcl,pkl,ptpk
        dpmix=ptpk-psp(i)
        if(abs(dpmix).lt.h3000)dpmix=-h3000
C--------------temperature propfile slope-------------------------------
        smix=(thtpk-thbt(i))/dpmix*stabs
C        print *,'thtpk, thbt, dpmix, smix  =',
C     *          thtpk,thbt(i),dpmix,smix
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        do 620 l=ltp1,lbtk
          ivi=ltp1+lbtk-l
C         lbm1=lbtk-1
C         trefk(lbtk)=thvmod(lbtk+1)/((qk(lbtk)*d608+1.)*apek(lbtk))
C             do 620 l=ltp1,lbm1
C             ivi=ltp1+lbm1-l
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          trefk(ivi)=((pk(ivi)-pk(ivi+1))*smix+trefk(ivi+1)*apek(ivi+1))
     +               /apek(ivi)
          dpkts=(pk(ivi)-pk(ivi+1))*smix
          thivi=trefk(ivi+1)*apek(ivi+1)
          summm=dpkts+thivi
          summd=summm/apek(ivi)
C       print *,'dpkts, thivip1, sum, sumd =',dpkts,thivi,summm,summd
C       print *,'ivi, pk, pk(ivi+1), trefk(ivi+1), tref(ivi) =',
C     * ivi,pk(ivi),pk(ivi+1),trefk(ivi+1),tref(ivi)
C      print *,'apek(ivi+1), apek(ivi) =',
C     *        apek(ivi+1),apek(ivi)
  620   continue
C--------------temperature reference profile correction-----------------
        sumdt=d00
        sumdp=d00
        do 630 l=ltp1,lbtk
          sumdt=(tk(l)-trefk(l))*dsigma(l)+sumdt
          sumdp=sumdp+dsigma(l)
  630   continue
C
        rdpsum=1./sumdp
        fpk(lbtk)=trefk(lbtk)
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        tcorr=sumdt*rdpsum
        do 640 l=ltp1,lbtk
C         tcorr=sumdt/(sumdp-dsigma(lbtk))
C             do 640 l=ltp1,lbm1
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          trfkl=trefk(l)+tcorr
          trefk(l)=trfkl
          fpk(l)=trfkl
  640   continue
C--------------humidity profile equations-------------------------------
        psum=0.0
        qsum=0.0
        potsum=0.0
        qotsum=0.0
        otsum=0.0
        dst=0.0
        fptk=fpk(ltp1)
        do 642 l=ltp1,lbtk
          dpkl=fpk(l)-fptk
          psum=dpkl*dsigma(l)+psum
          qsum=qk(l)*dsigma(l)+qsum
          rtbar=2./(trefk(l)+tk(l))
          otsum=dsigma(l)*rtbar+otsum
          potsum=dpkl*rtbar*dsigma(l)+potsum
          qotsum=qk(l)*rtbar*dsigma(l)+qotsum
          dst=(trefk(l)-tk(l))*rtbar*dsigma(l)+dst
  642   continue
C
        psum=psum*rdpsum
        qsum=qsum*rdpsum
        rotsum=1./otsum
        potsum=potsum*rotsum
        qotsum=qotsum*rotsum
        dst=dst*rotsum*cp/elwv
C--------------ensure positive entropy change---------------------------
        if(dst.gt.0.)then
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C             dstq=dst*epsup
          prtop(i)=pbot(i)
          ltop(i)=lbot(i)
          ndstp=ndstp+1
C           print *,'dst.gt.0...'
          goto 600
        else
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          dstq=dst*epsdn
        end if
C--------------check for isothermal atmosphere--------------------------
        den=potsum-psum
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        if(-den/psum.lt.0.00005)then
          ltop(i)=lbot(i)
          prtop(i)=pbot(i)
C          print *,'-den/psum .lt. 0.00005'
          goto 600
        else
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C--------------slope of the reference humidity profile------------------
          dqref=(qotsum-dstq-qsum)/den
        end if
C--------------humidity doesn`t increase with height--------------------
        if(dqref.lt.0.0)then
          ltop(i)=lbot(i)
          prtop(i)=pbot(i)
C      print *,'dqref.lt.0....'
          goto 600
        end if
C--------------humidity at the cloud top--------------------------------
        qrftp=qsum-dqref*psum
C--------------humidity profile-----------------------------------------
        do 644 l=ltp1,lbtk
          qrfkl=(fpk(l)-fptk)*dqref+qrftp
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C--------------supersaturation not allowed------------------------------
          qnew=(qrfkl-qk(l))*tauk+qk(l)
          if(qnew.gt.qsatk(l)*stresh)then
            ltop(i)=lbot(i)
            prtop(i)=pbot(i)
C      print *,'qnew.gt.qsatk(l)*stresh...'
            goto 600
          end if
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C-----------------------------------------------------------------------
          thvref(l)=trefk(l)*apek(l)*(qrfkl*d608+h1)
          qrefk(l)=qrfkl
  644   continue
C--------------eliminate impossible slopes (betts, dtheta/dq)-----------
        do 650 l=ltp1,lbtk
          dtdeta=(thvref(l-1)-thvref(l))/(a(l)-a(l-1))
          if(dtdeta.lt.epsth)then
            ltop(i)=lbot(i)
            prtop(i)=pbot(i)
C       print *,'dtdeat.lt.epsth...'
C       print *,'l, thvref, thvref(l-1), a, a(l-1) =',
C     * l,thvref(l),thvref(l-1),a(l),a(l-1)
C       if(j.eq.12.and. i.eq.59)
C    &     call fatal(__FILE__,__LINE__, 'T NOT CONVERGED')
            goto 600
          end if
  650   continue
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
Cvvvvvvvvvvvvv diagnostics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        if(dst.gt.0.)then
          ndstp=ndstp+1
        else
          ndstn=ndstn+1
        end if
        dentpy=d00
        do 660 l=ltp1,lbtk
  660   dentpy=((trefk(l)-tk(l))*cp+(qrefk(l)-qk(l))*elwv)/(tk(l)+
     +         trefk(l))*dsigma(l)+dentpy
Caaaaaaaaaaaaa diagnostics aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C--------------relaxation towards reference profiles--------------------
        iconss=iconss+1
        do 670 l=ltp1,lbtk
          tmod(i,l)=(trefk(l)-tk(l))/trel
          qqmod(i,l)=(qrefk(l)-qk(l))/trel
C         tmod(i,l)=(trefk(l)-tk(l))/dt2
C         qqmod(i,l)=(qrefk(l)-qk(l))/dt2
  670   continue
C       print *,'after do 670; iconss =',iconss
Cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
Cscscscscscsc  end of shallow convection   scscscscscscscscscscscscscscs
Cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
C-----------------------------------------------------------------------
  600 continue
C     print *, 'after loop 600'
Cvvvvvvvvvvvvv diagnostics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      nshal=0
      do 680 i=2,ilxm
        ltpk=ltop(i)
        lbtk=lbot(i)
        ptpk=prtop(i)
        pbtk=pbot(i)
C no shallow convection if cloud is not at least 90 mb or 3 sigma
C layers deep
C     if(cldhgt(i).lt.zno.or.ltop(i).gt.lbot(i)-2)go to 680
        if(cldhgt(i).lt.zno)goto 680
C     if(prtop(i).gt.pbot(i)-pno.or.ltop(i).gt.lbot(i)-2)goto 680
C shallow convection if cloud is .le. 290 mb
        if(cldhgt(i).lt.zsh)then
C     if(ptpk.ge.pbtk-psh)then
          nshal=nshal+1
          ntops(ltpk)=ntops(ltpk)+1
          nbots(lbtk)=nbots(lbtk)+1
          ndepth=lbtk-ltpk
          if(ndepth.gt.0)ndpths(ndepth)=ndpths(ndepth)+1
        end if
C
C  find cloud fractional cover and liquid water content
C
        kbaseb=min0(lbtk,kx-2)
        if(ltpk.le.kbaseb) then
          kclth = kbaseb-ltpk+1
          akclth= 1./dble(kclth)
          do k=ltpk,kbaseb
             cldlwc(i,k)=cllwcv
             cldfra(i,k)=1.-(1.-clfrcv)**akclth
          end do
        end if
Cchem2
        if(ichem.eq.1) then
C  before go to 100 put
          icumtop(i,j) = ltpk
          icumbot(i,j) = kbaseb
        end if
Cchem2_
  680 continue
C     print *,'at j =',j,', final no. shallow conv points =',nshal
      negds=khshal-nshal
C-----------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
         tten(i,k)=tten(i,k)+tmod(i,k)*psb(i,j)
         qten(i,k)=qten(i,k)+qqmod(i,k)*psb(i,j)
      end do
      end do
      icon(j)=icond
      return
      end subroutine bmpara
