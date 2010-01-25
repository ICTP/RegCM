C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    this file is part of regcm model.
C
C    regcm model is free software: you can redistribute it and/or modify
C    it under the terms of the gnu general public license as published by
C    the free software foundation, either version 3 of the license, or
C    (at your option) any later version.
C
C    regcm model is distributed in the hope that it will be useful,
C    but without any warranty; without even the implied warranty of
C    merchantability or fitness for a particular purpose.  see the
C    gnu general public license for more details.
C
C    you should have received a copy of the gnu general public license
C    along with regcm model.  if not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine lutbl(ptop)
C
      implicit none

      real(kind=8)  ptop
      real(kind=8)  a2,a3,a4
      integer itb,jtb
      parameter(itb=100,jtb=150,a2=17.2693882e0,a3=273.16e0,a4=35.86e0)
      real(kind=8)  pq0,eps
      parameter(pq0=379.90516e0,eps=2.e-12)
      real(kind=8)  eliwv,cp,r,rovcp
      parameter(eliwv=2.72e6,cp=1004.e0,r=287.04,rovcp=r/cp)
C
      real(kind=8)  pl,thl,rdq,rdth,rdp,rdthe,qs0,sqs,the0,sthe,ptbl
      real(kind=8)  ttbl
      common /tbls/ pl,thl,rdq,rdth,rdp,rdthe,qs0(jtb),sqs(jtb),
     1             the0(jtb),sthe(jtb),ptbl(itb,jtb),ttbl(jtb,itb)
      real(kind=8)  qsold(jtb),pold(jtb),qsnew(jtb),y2p(jtb),pnew(jtb)
      real(kind=8)  told(jtb),theold(jtb),thenew(jtb),y2t(jtb),tnew(jtb)
      real(kind=8)  app(jtb),aqp(jtb),apt(jtb),aqt(jtb)
      real(kind=8)  thh,ph,dth,dp,th,p,ape,qs0k,sqsk,dqs,qs,the0k,sthek
      real(kind=8)  dthe,pt
      integer kthm,kpm,kthm1,kpm1,kth,kp
C
C--------------coarse look-up table for saturation point----------------
C
      pt = ptop*1000.
! ptop in pascal

      kthm  = jtb
      kpm   = itb
      kthm1 = kthm-1
      kpm1  = kpm-1
C
      thl = 210.
      thh = 385.
      pl  = pt
      ph  = 105000.
C
      dth = (thh-thl)/dble(kthm-1)
      dp  = (ph-pl)/dble(kpm-1)
C
      rdth = 1.0/dth
      rdp  = 1.0/dp
      rdq  = kpm-1
C
      th   = thl-dth

C-----------------------------------------------------------------------

      do kth = 1, kthm
        th = th+dth
        p  = pl-dp
        do kp = 1, kpm
          p         = p+dp
          ape       = (100000./p)**(rovcp)
          qsold(kp) = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          pold(kp)  = p
        end do
C
        qs0k       = qsold(1)
        sqsk       = qsold(kpm)-qsold(1)
        qsold(1)   = 0.0
        qsold(kpm) = 1.0
C
        do kp = 2, kpm1
          qsold(kp) = (qsold(kp)-qs0k)/sqsk
Cwwwwwwwwwwwwww fix due to cyber half prec. limitation wwwwwwwwwwwwwwwww
          if ((qsold(kp)-qsold(kp-1)).lt.eps) qsold(kp)=qsold(kp-1)+eps
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        end do
C
        qs0(kth) = qs0k
        sqs(kth) = sqsk
C-----------------------------------------------------------------------
        qsnew(1)   = 0.0
        qsnew(kpm) = 1.0
        dqs        = 1.0/dble(kpm-1)
C
        do kp = 2, kpm1
          qsnew(kp)=qsnew(kp-1)+dqs
        end do
C
        y2p(1)   = 0.0
        y2p(kpm) = 0.0
C
        call spline(kpm,qsold,pold,y2p,kpm,qsnew,pnew,app,aqp)
C
        do kp = 1, kpm
          ptbl(kp,kth) = pnew(kp)
        end do
C-----------------------------------------------------------------------
      end do
C-----------------------------------------------------------------------

C--------------coarse look-up table for t(p) from constant the----------
      p = pl-dp
      do kp = 1, kpm
        p  = p+dp
        th = thl-dth
        do kth=1,kthm
          th          = th+dth
          ape         = (100000./p)**(rovcp)
          qs          = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          told(kth)   = th/ape
          theold(kth) = th*exp(eliwv*qs/(cp*told(kth)))
        end do
C
        the0k        = theold(1)
        sthek        = theold(kthm)-theold(1)
        theold(1)    = 0.0
        theold(kthm) = 1.0
C
        do kth = 2, kthm1
          theold(kth) = (theold(kth)-the0k)/sthek
Cwwwwwwwwwwwwww fix due to cyber half prec. limitation wwwwwwwwwwwwwwwww
          if ((theold(kth)-theold(kth-1)).lt.eps) theold(kth) =
     +      theold(kth-1)+eps
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        end do 
C
        the0(kp) = the0k
        sthe(kp) = sthek
C-----------------------------------------------------------------------
        thenew(1)    = 0.0
        thenew(kthm) = 1.0
        dthe         = 1.0/dble(kthm-1)
        rdthe        = 1.0/dthe
C
        do kth = 2, kthm1
          thenew(kth) = thenew(kth-1)+dthe
        end do
C
        y2t(1)    = 0.0
        y2t(kthm) = 0.0
C
        call spline(kthm,theold,told,y2t,kthm,thenew,tnew,apt,aqt)
C
        do kth = 1, kthm
          ttbl(kth,kp) = tnew(kth)
        end do
C-----------------------------------------------------------------------
      end do
C-----------------------------------------------------------------------
C
      return
      end subroutine lutbl
