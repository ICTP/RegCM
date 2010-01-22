      implicit none
      include '../Main/regcm.param'
      integer ixb,jxb
      parameter( ixb=ix-1, jxb=jx-1 )
      integer ilx,jlx,kxp1,nsplit
      parameter(ilx=ix-1,jlx=jx-1,kxp1=kx+1,nsplit=2)
      integer mdate0,ktau,ldatez,lyear,lmonth,lday,lhour
      integer*8 ntime
      real*8  xtime
      real*4  UB0(IX,JX,KX),VB0(IX,JX,KX),QB0(IX,JX,KX)
     &       ,TB0(IX,JX,KX),PS0(IX,JX),TS0(IX,JX)
      real*8  ua(ix,jx,kx), ub(ix,jx,kx), va(ix,jx,kx),
     &        vb(ix,jx,kx), ta(ix,jx,kx), tb(ix,jx,kx),
     &        qva(ix,jx,kx),qvb(ix,jx,kx),qca(ix,jx,kx),
     &        qcb(ix,jx,kx)
      real*8  psa(ix,jx),psb(ix,jx),satbrt(ix,jx),f(ix,jx)
      real*8  ht(ix,jx),msfx(ix,jx),msfd(ix,jx),xlat(ix,jx),xlong(ix,jx)
      real*8  tga(ix,jx),tgb(ix,jx),rainc(ix,jx),rainnc(ix,jx)
      real*8  hfx(ix,jx),qfx(ix,jx),snowc(ix,jx),uvdrag(ix,jx)
      real*8  tdini,tdadv,tqini,tqadv,tqeva,tqrai
      real*8  absnxt(ix-1,kx,4,jx-1)
     &      , abstot(ix-1,kx,kx,jx-1),emstot(ix-1,kx,jx-1)
      real*8  fcc(ix,jx,kx),rh0(ix,jx)
      real*8  sol2d(ixb,jxb),solvd2d(ixb,jxb),solvs2d(ixb,jxb)
     &       ,flw2d(ixb,jxb),flwd2d(ixb,jxb),fsw2d(ixb,jxb)
     &       ,sabv2d(ixb,jxb),sinc2d(ixb,jxb)
      real*8  taf2d(ixb,jxb),tlef2d(ixb,jxb),tgbb(ix,jx)
     &       ,ssw2d(ixb,jxb),srw2d(ixb,jxb),tgb2d(ixb,jxb)
     &       ,swt2d(ixb,jxb),scv2d(ixb,jxb),gwet2d(ixb,jxb)
     &       ,veg2d(ixb,jxb),sag2d(ixb,jxb),sice2d(ixb,jxb)
     &       ,dew2d(ixb,jxb),ircp2d(ixb,jxb),text2d(ixb,jxb)
     &       ,col2d(ixb,jxb),ocld2d(ixb,jxb)
     &       ,rsheat(ix,jx,kx),rswat(ix,jx,kx)
     &       ,heatrt(ilx,jlx,kx),o3prof(ilx,jlx,kxp1)
      real*8  pptnc(ixb,jxb),pptc(ixb,jxb),prca2d(ixb,jxb)
     &       ,prnca2d(ixb,jxb)
      real*8  dstor(ix,jx,nsplit),hstor(ix,jx,nsplit)
      real*8  uj1(ix,kx), uj2(ix,kx), ujlx(ix,kx), ujl(ix,kx),
     &        vj1(ix,kx), vj2(ix,kx), vjlx(ix,kx), vjl(ix,kx),
     &        ui1(jx,kx), ui2(jx,kx), uilx(jx,kx), uil(jx,kx),
     &        vi1(jx,kx), vi2(jx,kx), vilx(jx,kx), vil(jx,kx)
      integer iutl
c
      iutl=14
      read (iutl) mdate0
      read (iutl) ktau,xtime,ldatez,lyear,lmonth,lday,lhour,ntime
      read (iutl) UB0,VB0,QB0,TB0,PS0,TS0
      read (iutl) ua
      read (iutl) ub
      read (iutl) va
      read (iutl) vb
      read (iutl) ta
      read (iutl) tb
      read (iutl) qva
      read (iutl) qvb
      read (iutl) qca
      read (iutl) qcb
      read (iutl) psa,psb,satbrt,f
      read (iutl) ht,msfx,msfd,xlat,xlong
      read (iutl) tga,tgb,rainc,rainnc
      read (iutl) hfx,qfx,snowc,uvdrag
      read (iutl) tdini,tdadv,tqini,tqadv,tqeva,tqrai
      read (iutl) absnxt, abstot, emstot
      read (iutl)  fcc,rh0
      read (iutl) sol2d, solvd2d, solvs2d
     &           ,  flw2d, flwd2d, fsw2d
     &           , sabv2d, sinc2d
      read (iutl)  taf2d, tlef2d, tgbb
     b           ,  ssw2d,  srw2d
     c           ,  tgb2d,  swt2d
     d           ,  scv2d, gwet2d,  veg2d
     e           ,  sag2d, sice2d,  dew2d
     &           , ircp2d
     g           , text2d,  col2d, ocld2d
     h           , rsheat,  rswat, heatrt, o3prof
      read (iutl)  pptnc, pptc, prca2d, prnca2d
      read (iutl) dstor
      read (iutl) hstor
      read (iutl) uj1,uj2,ujlx,ujl
      read (iutl) ui1,ui2,uilx,uil
      read (iutl) vj1,vj2,vjlx,vjl
      read (iutl) vi1,vi2,vilx,vil
c
      stop
      end
