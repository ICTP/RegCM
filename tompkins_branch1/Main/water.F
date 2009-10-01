
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine water
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c              update soil moisture and runoff
c
c     new algorithms for three soil layers (dickinson & kennedy 8-88)
c     calculate fluxes through air, surface layer, and root layer faces
c
c          b = b of clapp and hornberger
c       est0 = soil water flux, out top
c      gwatr = net input of water to the soil surface
c       ircp = leaf interception
c     wflux1 = soil water flux, 10 cm
c     wflux2 = soil water flux, 1 m
c     rsubss = soil water flux by grav. drainage thru 10 cm interface
c     rsubsr = soil water flux by grav. drainage thru 1 m interface
c     rsubst = soil water flux by grav. drainage thru 10 m interface
c       rsur = surface runoff
c     rno1d(n,np) = total runoff (mm/day)
c     rnos1d(n,np) = surface runoff (mm/day)
c
c     xkmxr and wflux1 determine flow thru upper/root soil interface
c     evmxt, xkmx1, and xkmx2 determine flow thru lower interfaces
c
c     veg type 10 "irrigated crop" is irrigated through reducing
c          the runoff (rsur), i.e., by adding a negative number
c          if the land isn't at least 60% saturated.
c     veg type 13 and 14 are water covered (lake, swamp, rice paddy);
c          negative runoff keeps this land saturated.
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'

      real*8  rsur(NNSG,nbmax),   rnof(NNSG,nbmax),  gwatr(NNSG,nbmax)
      real*8  wfluxc(NNSG,nbmax), wflux1(NNSG,nbmax), wflux2(NNSG,nbmax)
      real*8  rsubsr(NNSG,nbmax), rsubss(NNSG,nbmax), rsubst(NNSG,nbmax)
      real*8  xkmxr(NNSG,nbmax),  xkmx1(NNSG,nbmax),  xkmx2(NNSG,nbmax)

      integer n,np
      real*8  evmxr,evmxt,b,bfac,bfac2,delwat,xxkb,vakb,wtg2c
      real*8  rap,est0,evmax

c=======================================================================
cl          1.   define soil water fluxes
c=======================================================================
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then

cl          1.1  reduce infiltration for frozen ground
        if(tgb1d(n,np).gt.c(67)) then
         xkmxr(n,np) = xkmx(n,np)
        else
         xkmxr(n,np) = 0.
        endif

cl          1.11 permafrost or ice sheet
        if (lveg(n,np).eq.9 .or. lveg(n,np).eq.12) then
         xkmx1(n,np)=0.
         xkmx2(n,np)=0.
        else
         xkmx1(n,np)=xkmx(n,np)
         xkmx2(n,np)=drain
        endif
cl          1.2  diffusive fluxes
        evmxr = evmx0(n,np)*xkmxr(n,np) / xkmx(n,np)
        evmxt = evmx0(n,np)*xkmx1(n,np) / xkmx(n,np)
        b     = bsw(n,np)
        bfac  = watr(n,np)**(3.+bfc(n,np))*watu(n,np)**(b-bfc(n,np)-1)
        bfac2 = watt(n,np)**(2.+bfc(n,np))*watr(n,np)**(b-bfc(n,np))
        wfluxc(n,np)=evmxr
     &              *(depuv(lveg(n,np))/deprv(lveg(n,np)))**0.4*bfac
        wflux1(n,np)=wfluxc(n,np)*watr(n,np)
        wflux2(n,np)=evmxt*(depuv(lveg(n,np))/deprv(lveg(n,np)))**0.5
     a           * bfac2*(watt(n,np)-watr(n,np))
cl          1.3  gravitational drainage
        rsubss(n,np)=xkmxr(n,np)*watr(n,np)**(b+0.5)*watu(n,np)**(b+2.5)
        rsubsr(n,np)=xkmx1(n,np)*watt(n,np)**(b+0.5)*watr(n,np)**(b+2.5)
        rsubst(n,np)=xkmx2(n,np)*watt(n,np)**(2.*b+3.)
cl             1.32 bog or water
        if((lveg(n,np).ge.13).and.(lveg(n,np).le.15))   then
         rsubst(n,np) = 0.0
         rsubss(n,np) = 0.0
         rsubsr(n,np) = 0.0
        endif
cl          1.4  fluxes through internal surfaces
        wflux1(n,np) = wflux1(n,np)-rsubss(n,np)
        wflux2(n,np) = wflux2(n,np)-rsubsr(n,np)
       endif
      enddo
90    continue

cl          1.5  net flux at air interface
      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then

        gwatr(n,np)=pw(n,np)-evapw(n,np)+sm(n,np)+etrrun(n,np)/c(4)
c=======================================================================
cl          2.   define runoff terms
c=======================================================================
cl          2.1  surface runoff
        wata(n,np) = 0.5*(watu(n,np)+watr(n,np))
cl          2.11 increase surface runoff over frozen ground
        if (tg1d(n,np).lt.c(67)) then
         rsur(n,np) = dmin1(1.d0,wata(n,np)**1)*dmax1(0.d0,gwatr(n,np))
        else
         rsur(n,np) = dmin1(1.d0,wata(n,np)**4)*dmax1(0.d0,gwatr(n,np) )
        endif
cl          2.12 irrigate cropland
        if (lveg(n,np).eq.10.and.watr(n,np).lt.relfc(n,np))
     &  rsur(n,np)=rsur(n,np)+(rsw1d(n,np)-relfc(n,np)*gwmx1(n,np))/c(4)
cl          2.13 saturate swamp or rice paddy
        if((lveg(n,np).ge.13).and.(lveg(n,np).lt.16))
     &  rsur(n,np)=rsur(n,np)+dmin1(0.d0,(rsw1d(n,np)-gwmx1(n,np))/c(4))
cl          2.2  total runoff
        rnof(n,np) = rsur(n,np) + rsubst(n,np)
c=======================================================================
cl          3.   increment soil moisture
c=======================================================================
cl          3.1  update top layer with implicit treatment
cl                    of flux from below
        ssw1d(n,np) = ssw1d(n,np)
     &  +c(4)*(gwatr(n,np)-efpr(n,np)*etr(n,np)-rsur(n,np)+wflux1(n,np))
        ssw1d(n,np)=ssw1d(n,np)/(1.+wfluxc(n,np)*c(4)/gwmx0(n,np))
cl          3.2  update root zone
        rsw1d(n,np)=rsw1d(n,np)
     &             +c(4)*(gwatr(n,np)-etr(n,np)-rsur(n,np)+wflux2(n,np))
cl          3.3  update total water
        tsw1d(n,np)=tsw1d(n,np)+c(4)*(gwatr(n,np)-etr(n,np)-rnof(n,np))
       endif
      enddo
91    continue

      do 92 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then
c=======================================================================
cl          4.   check whether soil water exceeds maximum capacity or
cl               becomes negative (should rarely or never happen)
c=======================================================================
cl          4.1  surface water assumed to move downward into soil
        if(ssw1d(n,np).gt.gwmx0(n,np)) ssw1d(n,np)=gwmx0(n,np)

cl          4.2  excess root layer water assumed to move downward
        if(rsw1d(n,np).gt.gwmx1(n,np)) rsw1d(n,np)=gwmx1(n,np)

cl          4.3  excess total water assumed to go to subsurface runoff
        if(tsw1d(n,np).gt.gwmx2(n,np))then
         delwat=tsw1d(n,np)-gwmx2(n,np)
         tsw1d(n,np)=gwmx2(n,np)
         rsubst(n,np)=rsubst(n,np)+delwat/c(4)
        endif
cl          4.4  check for negative water in top layer
        if(ssw1d(n,np).le.1.e-2) ssw1d(n,np)=1.e-2
c=======================================================================
cl          5.   accumulate leaf interception
c=======================================================================
        ircp1d(n,np) = ircp1d(n,np)
     &      +sigf(n,np)*(c(4)*prcp1d(n,np))-(sdrop(n,np)+etrrun(n,np))
c=======================================================================
cl          6.   evaluate runoff (incremented in ccm)
c=======================================================================
c*  update total runoff
        rnof(n,np) = rsur(n,np) + rsubst(n,np)
        rno1d(n,np)  = rnof(n,np)*tau1
        rnos1d(n,np) = rsur(n,np)*tau1
       else                          ! ocean or sea ice
        rnof(n,np) = 0.
        rno1d(n,np)  = 0.
        rnos1d(n,np) = 0.
       endif
      enddo
92    continue

c=======================================================================
cl          7.   calculate potential evaporation and use to determine
cl              wetness factor, allowing for snow being saturated
c=======================================================================
      do 93 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then
        xxkb = dmin1(rough(lveg(n,np)),1.d0)
        vakb = (1.-sigf(n,np))*vspda(n,np)
     &     + sigf(n,np)*(xxkb*uaf(n,np)+(1.-xxkb)*vspda(n,np))
        wtg2c = (1.-sigf(n,np))*cdrx(n,np)*vakb
        rap  = rhs1d(n,np) * (csoilc*uaf(n,np)*sigf(n,np)*(qg1d(n,np)
     &       +delq1d(n,np)-qs1d(n,np))
     &       + wtg2c*(qg1d(n,np)-qs1d(n,np)))
        bfac = watr(n,np)**(3.+bfc(n,np))
     &        *watu(n,np)**(bsw(n,np)-bfc(n,np)-1)
        est0 = evmx0(n,np)*bfac*watu(n,np)
        evmax = dmax1(est0,0.d0)
        gwet1d(n,np) = dmin1(1.d0,evmax/dmax1(1.d-14,rap))
        gwet1d(n,np) = scvk(n,np) + gwet1d(n,np)*(1.0-scvk(n,np))
       endif
      enddo
93    continue

      return
      end
