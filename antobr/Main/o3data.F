c----------------------------------------------------------------------
      subroutine o3data
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
      include 'Commons/param3.cb'
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/rad.cb'
      include 'mpiregcm.h'
#else
      include 'Commons/main.cb'
      include 'Commons/rad.cb'
#endif
      real*8  prlevh,ppwrkh,o3wrk,ppwrk,o3sum,ppsum,o3win,ppwin
      real*8  o3ann ,ppann
      common/o3blk/ prlevh(kxp2),ppwrkh(32),
     & o3wrk(31),ppwrk(31),o3sum(31),ppsum(31),o3win(31),ppwin(31),
     & o3ann(31),ppann(31)
c
      integer i,j,k,klevm1,klevp1,klevp2,kj,jj
      real*8  pb1,pb2,pt1,pt2
c
      do 501 k=1,31
       ppann(k)=ppsum(k)
501   continue
      o3ann(1)=0.5*(o3sum(1)+o3win(1))
c
      do 502 k=2,31
       o3ann(k)=o3win(k-1)+(o3win(k)-o3win(k-1))/(ppwin(k)-ppwin(k-1))
     &          * (ppsum(k)-ppwin(k-1))
502   continue
      do 503 k=2,31
       o3ann(k)=0.5*(o3ann(k)+o3sum(k))
503   continue
      do 2 k=1,31
       o3wrk(k)=o3ann(k)
       ppwrk(k)=ppann(k)
2     continue
c
c   calculate half pressure levels for model and data levels
c
      klevm1=kx-1
      klevp1=kx+1
      klevp2=kx+2
c
#ifdef MPP1
      do 101 j=1,jendx
#else
      do 101 j=1,jlx
#endif
       do 100 i=1,ilx
        do 10 k=klevp1,1,-1
         kj=klevp1-k+1
         prlevh(kj)=(sigma(k)*psb(i,j)+ptop)*10.
10      continue
        ppwrkh(1)=1100.
        do 20 k=2,31
         ppwrkh(k)=(ppwrk(k)+ppwrk(k-1))/2.
20      continue
        ppwrkh(32)=0.
        do 30 k=1,kx
         o3prof(i,k,j)=0.
         do 25 jj=1,31
          if((-(prlevh(k)-ppwrkh(jj))).ge.0.)then
           pb1=0.
          else
           pb1=prlevh(k)-ppwrkh(jj)
          endif
          if((-(prlevh(k)-ppwrkh(jj+1))).ge.0.)then
           pb2=0.
          else
           pb2=prlevh(k)-ppwrkh(jj+1)
          endif
          if((-(prlevh(k+1)-ppwrkh(jj))).ge.0.)then
           pt1=0.
          else
           pt1=prlevh(k+1)-ppwrkh(jj)
          endif
          if((-(prlevh(k+1)-ppwrkh(jj+1))).ge.0.)then
           pt2=0.
          else
           pt2=prlevh(k+1)-ppwrkh(jj+1)
          endif
          o3prof(i,k,j)=o3prof(i,k,j)+(pb2-pb1-pt2+pt1)*o3wrk(jj)
25       continue
         o3prof(i,k,j)=o3prof(i,k,j)/(prlevh(k)-prlevh(k+1))
30      continue
100    continue
101   continue
c
      return
      end
