c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine initb
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ***  provides constant initial fields to boundary subroutine
c     ***    soil textures are set in soilbc
c     ***    soil   colors are set in albedo
c
c     ***   units are si
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/iunits.cb'
      include 'Commons/bats.cb'
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/bats.cb2'
      include 'mpiregcm.h'
#else
      include 'Commons/main.cb'
      include 'Commons/bats.cb2'
#endif
c
c----------------------------------------------------------------------
c     slmo   : surface moisture availability in fraction of one.
      real*8  slmo(20)                !BATS land types
c     data slmo/0.75,0.60,0.75,0.75,0.75,
c    &          0.80,0.60,0.10,0.90,0.90,
c    &          0.10,0.95,0.90,1.00,1.00,
c    &          0.60,0.75,0.85,0.85,0.85 /
      data slmo/0.65,0.45,0.60,0.60,0.65,
     &          0.65,0.55,0.10,0.90,0.80,
     &          0.20,0.90,0.90,1.00,1.00,
     &          0.50,0.50,0.65,0.60,0.60/

c ****** vgtran transforms mm4 surface types into bats surface types
cMM4  real*8  vgtran(13)
cMM4  data vgtran /5.,1.,2.,5.,3.,18.,0.,13.,11.,9.,12.,6.,7./
cMM4  integer ist
      integer is,i,ill,jll,nlveg,itex,k

c ****** typically resp=1.0 kg/m**2/s and changes by<10% in 10 days
c
#ifdef MPP1
      do 18 jll =1,jendx
#else
      do 18 jll =1,jx-1
#endif
      do 18 ill =1,ix-1
       pptnc(ill,jll) = 0.
       pptc(ill,jll) = 0.
cMM4   ist=nint(satbrt(ill,jll))
cMM4   if(ist.le.13)then
cMM4    ist=ist
cMM4   else
cMM4    ist=7
cMM4   endif
cMM4   veg2d(ill,jll)=vgtran(ist)
ceros  note:  when using bats dataset comment line above
       veg2d(ill,jll)=satbrt(ill,jll)
       if (satbrt(ill,jll).gt.13.9.and.satbrt(ill,jll).lt.15.1)
     &        veg2d(ill,jll) = 0.
18    continue
#ifdef MPP1
      do jll =1,jendx
#else
      do jll =1,jx-1
#endif
      do ill=1,ix-1
      do k=1,NNSG
         veg2d1(k,ill,jll)=satbrt1(k,ill,jll)
         if(satbrt1(k,ill,jll).gt.13.9.and.satbrt1(k,ill,jll).lt.15.1)
     &       veg2d1(k,ill,jll) = 0.
      enddo
      enddo
      enddo

c ******  initialize hostetler lake model
      if (lakemod.eq.1) then
#ifdef MPP1
       call initlk(veg2d1, (ix-1)*NSG, jxp*NSG)
#else
       call initlk(veg2d1, (ix-1)*NSG, (jx-1)*NSG)
#endif
      endif

#ifdef MPP1
      do 20 jll =1,jendx
#else
      do 20 jll =1,jx-1
#endif
      do 20 ill =1,ix-1
      do 20 k=1,NNSG
       if(veg2d1(k,ill,jll).lt.0.5)then
        ocld2d(k,ill,jll)=0.
       else
        ocld2d(k,ill,jll)=1.
       endif
       nlveg=nint(veg2d1(k,ill,jll))
       if(nlveg.eq.0)then
        nlveg=15
       else
        nlveg=nlveg
       endif
csol   itex=int(text2d(k,ill,jll))
       itex=iexsol(nlveg)
       tg2d(k,ill,jll)=tgb(ill,jll)
       tgb2d(k,ill,jll)=tgb(ill,jll)
       taf2d(k,ill,jll)=tgb(ill,jll)
       tlef2d(k,ill,jll)=tgb(ill,jll)

c ******  initialize soil moisture in the 3 layers
       is=nint(satbrt1(k,ill,jll))
       swt2d(k,ill,jll)=deptv(nlveg)*xmopor(itex)*slmo(is)
       srw2d(k,ill,jll)=deprv(nlveg)*xmopor(itex)*slmo(is)
       ssw2d(k,ill,jll)=depuv(nlveg)*xmopor(itex)*slmo(is)

       dew2d(k,ill,jll)=0.
       sag2d(k,ill,jll)=0.
       scv2d(k,ill,jll)=dmax1(snowc(k,ill,jll),0.d0)
       sice2d(k,ill,jll)=0.
       gwet2d(k,ill,jll)=0.5
       sena2d(k,ill,jll)=0.
       evpa2d(k,ill,jll)=0.
       rnos2d(k,ill,jll)=0.
       rno2d(k,ill,jll)=0.
       if(sice2d(k,ill,jll).gt.0.)then
        ocld2d(k,ill,jll)=2.
       endif
       ircp2d(k,ill,jll)=0.
20    continue

#ifdef MPP1
      do 30 jll =1,jendx
#else
      do 30 jll =1,jx-1
#endif
      do 30 ill =1,ix-1
       fsw2d(ill,jll)=0.
       flw2d(ill,jll)=0.
       sabv2d(ill,jll)=0.
       sol2d(ill,jll)=0.
       fswa2d(ill,jll)=0.
       flwa2d(ill,jll)=0.
       prca2d(ill,jll)=0.
       prnca2d(ill,jll)=0.
       svga2d(ill,jll)=0.
       sina2d(ill,jll)=0.
30    continue

      return
      end
