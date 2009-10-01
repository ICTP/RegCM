ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vadv(ften,fa,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the vertical flux-divergence terms.    c
c                                                                     c
c     ften   : is the tendency of variable 'f'.                       c
c                                                                     c
c     fa     : is p*f.                                                c
c                                                                     c
c     qdot   : is the vertical sigma-velocity                         c
c                                                                     c
c     f      : is the working space used to store the interlated      c
c              values.                                                c
c                                                                     c
c     psa    : is p* used to interpolate the temperature.             c
c                                                                     c
c     j      : j'th slice of variable fa.                             c
c                                                                     c
c     ind = 1 : for t.                                                c
c           2 : for qv.                                               c
c           3 : for qc and qr.                                        c
c           4 : for u and v.                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/main.cb'
      include 'Common2/cvaria.cb'
#else
      include 'Commons/main.cb'
      include 'Commons/cvaria.cb'
#endif
      include 'Commons/param1.cb'
      include 'Commons/param3.cb'
c
      real*8  ften(ix,kx),fa(ix,kx),fg(ix,kx)
c
      integer i,k
      real*8  f1,f2
c
c----------------------------------------------------------------------
c
      if(ind.eq.1)then
c
c-----vertical advection terms for:
c.....interpolate ta to full sigma levels:
c
       do 20 i=2,ilxm
        fg(i,1) = 0.
20     continue
       do 30 k=2,kx
       do 30 i=2,ilxm
        fg(i,k)=twt(k,1)*fa(i,k)*((psa(i,j)*sigma(k)+ptop)/
     1       (psa(i,j)*a(k)+ptop))**0.287 + twt(k,2)*fa(i,k-1)*
     2       ((psa(i,j)*sigma(k)+ptop)/(psa(i,j)*a(k-1)+ptop))**0.287
30     continue
c......k = 1
       do 35 i=2,ilxm
        ften(i,1)=ften(i,1)-qdot(i,2,j)*fg(i,2)/dsigma(1)
35     continue
c......k = 2,kxm
       do 40 k=2,kxm
       do 40 i=2,ilxm
        ften(i,k)=ften(i,k)-(qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*
     1                     fg(i,k))/dsigma(k)
40     continue
c,.....k = kx
       do 45 i=2,ilxm
        ften(i,kx)=ften(i,kx)+qdot(i,kx,j)*fg(i,kx)/dsigma(kx)
45     continue
c
      elseif(ind.eq.2)then
c
c-----vertical advection term for qv:
c.....interpolate qv to full sigma levels:
c
       do 50 i=2,ilxm
        fg(i,1) = 0.
50     continue
       do 60 k=2,kx
       do 60 i=2,ilxm
c modif !!
        if (fa(i,k).gt.1.e-15 .AND.fa(i,k-1).gt. 1.e-15)then 
         fg(i,k)=fa(i,k)*(fa(i,k-1)/fa(i,k))**qcon(k)
        else
         fg(i,k)= 0.
        end if
60     continue
c......k = 1
       do 80 i=2,ilxm
        ften(i,1)=ften(i,1)-qdot(i,2,j)*fg(i,2)/dsigma(1)
80     continue
c......k = 2,kxm
       do 90 k=2,kxm
       do 90 i=2,ilxm
        ften(i,k)=ften(i,k)-(qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*
     1                     fg(i,k))/dsigma(k)
90     continue
c,.....k = kx
       do 100 i=2,ilxm
        ften(i,kx)=ften(i,kx)+qdot(i,kx,j)*fg(i,kx)/dsigma(kx)
100    continue
c
      elseif(ind.eq.3)then
c
c-----vertical advection terms for qc and qr:
c
c......k = 1
       do 120 i=2,ilxm
        if(qdot(i,2,j).ge.0.)then
         f2=fa(i,1)
        else
         f2=fa(i,2)
        endif
        ften(i,1)=ften(i,1)-qdot(i,2,j)*f2/dsigma(1)
120    continue
c......k = 2,kxm
       do 130 k=2,kxm
       do 130 i=2,ilxm
        if(qdot(i,k+1,j).ge.0.)then
         f2=fa(i,k)
        else
         f2=fa(i,k+1)
        endif
        if(qdot(i,k,j).ge.0.)then
         f1=fa(i,k-1)
        else
         f1=fa(i,k)
        endif
        ften(i,k)=ften(i,k)-(qdot(i,k+1,j)*f2-
     1            qdot(i,k,j)*f1)/dsigma(k)
130    continue
c......k = kx
       do 140 i=2,ilxm
        if(qdot(i,kx,j).ge.0.)then
         f1=fa(i,kxm)
        else
         f1=fa(i,kx)
        endif
        ften(i,kx)=ften(i,kx)+qdot(i,kx,j)*f1/dsigma(kx)
140    continue
c
      elseif(ind.eq.4)then
c
c-----vertical advection terms for u and v:
c.....interpolate ua or va to full sigma levels:
c
       do 160 i=2,ilx
        fg(i,1) = 0.
160    continue
       do 170 k=2,kx
       do 170 i=2,ilx
        fg(i,k)=0.5*(fa(i,k)+fa(i,k-1))/msfd(i,j)
170    continue
c......k = 1
       do 180 i=2,ilx
        ften(i,1)=ften(i,1)-(qdot(i-1,2,j-1)+qdot(i,2,j-1)+qdot(i,2,j)+
     1                     qdot(i-1,2,j))*fg(i,2)/(4.*dsigma(1))
180    continue
c......k = 2,kxm
       do 190 k=2,kxm
       do 190 i=2,ilx
        ften(i,k)=ften(i,k)-((qdot(i,k+1,j-1)+qdot(i-1,k+1,j-1)+
     1          qdot(i,k+1,j)+qdot(i-1,k+1,j))*fg(i,k+1)-(qdot(i,k,j-1)+
     2          qdot(i-1,k,j-1)+qdot(i,k,j)+qdot(i-1,k,j))*fg(i,k))/
     3          (4.*dsigma(k))
190    continue
c......k = kx
       do 200 i=2,ilx
        ften(i,kx)=ften(i,kx)+(qdot(i,kx,j-1)+qdot(i-1,kx,j-1)+
     1             qdot(i,kx,j)+qdot(i-1,kx,j))*fg(i,kx)/
     2             (4.*dsigma(kx))
200    continue
c


      elseif(ind.eq.5)then
       
        do 210  k=2,kx
        do 210  i=2,ilxm
            fg(i,k)= twt(k,1)*fa(i,k) + twt(k,2)*fa(i,k-1)
210     continue        

c......k = 1
       do 220 i=2,ilxm
        ften(i,1)=ften(i,1)-qdot(i,2,j)*fg(i,2)/dsigma(1)
220     continue
c......k = 2,kxm
       do 230 k=2,kxm
       do 230 i=2,ilxm
        ften(i,k)=ften(i,k)-(qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*
     1                     fg(i,k))/dsigma(k)
230     continue
c,.....k = kx
       do 240 i=2,ilxm
        ften(i,kx)=ften(i,kx)+qdot(i,kx,j)*fg(i,kx)/dsigma(kx)
240    continue

      endif
c
      return
      end
