ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cqhy15/07/99 subroutine cumtran(il,jl,kx,ntr,dsigma,icumtop,chia,chib)
      subroutine cumtran
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/trachem.cb'
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
#else
      include 'Commons/trachem.cb'
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
#endif
      include 'Commons/param3.cb'
      real*8  cumfrc,deltas,chiabar,chibbar
      integer i,j,k,n,kcumtop

chy the well-mixing only over cumulus cloud fraction (assuming 0.3)
c add the well-mixing fraction in dt step (5%) - only updraft
chy     cumfrc = 0.3
        cumfrc = 0.015

      do 100 n=1,ntr
#ifdef MPP1
      do 100 j=jbegin,jendx
#else
      do 100 j=2,jx-1
#endif
      do 100 i=2,ix-1
      if (icumtop(i,j).gt.0) then
         deltas=0.
         chiabar=0.
         chibbar=0.
         kcumtop=max0(icumtop(i,j),4)
         do 90 k=kcumtop,kx
         deltas=deltas+dsigma(k)
         chiabar=chiabar+chia(i,k,j,n)*dsigma(k)
         chibbar=chibbar+chib(i,k,j,n)*dsigma(k)
  90     continue
c?         do 95 k=icumtop(i,j),kx      ! yhuang, 12/98
cqhy
         do 95 k=kcumtop,kx
         chia(i,k,j,n)=chia(i,k,j,n)*(1.-cumfrc)+cumfrc*chiabar/deltas
         chib(i,k,j,n)=chib(i,k,j,n)*(1.-cumfrc)+cumfrc*chibbar/deltas
  95     continue
      endif
  100 continue
      return
      end
