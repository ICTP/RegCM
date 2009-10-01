cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spstep( hbar, dx2, dtau, m)
c
      implicit none
      real*8  dx2
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/main.cb'
      real*8  ddsum(ix,jxp,nsplit), dhsum(ix,0:jxp,nsplit)
      real*8  deld(ix,jxp,nsplit,3), delh(ix,0:jxp,nsplit,3)
      real*8  psdot(ix,jxp)
      real*8  work(ix,jxp,3)
      common /bxq_aaa/ddsum,dhsum,deld,delh,psdot,work
      real*8  uu,vv
      common /bxq_tmp/ uu(ix,jxp+1),vv(ix,jxp+1)
      real*8  wksend(ix*2),wkrecv(ix*2)
#else
      include 'Commons/main.cb'
      real*8  ddsum(ix,jx,nsplit), dhsum(ix,jx,nsplit)
      real*8  deld(ix,jx,nsplit,3), delh(ix,jx,nsplit,3)
      real*8  psdot(ix,jx)
      real*8  work(ix,jx,3)
      common /bxq_aaa/ddsum,dhsum,deld,delh,psdot,work
      real*8  uu,vv
      common /bxq_tmp/ uu(ix,jx),vv(ix,jx)
#endif
      integer m(nsplit)
      real*8  hbar(nsplit), dtau(nsplit)
      integer ns,n0,n1,n2,m2,i,j,nw,n
      real*8  dtau2,fac
c
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        ddsum(i,j,n) = 0.
        dhsum(i,j,n) = 0.
      enddo
      enddo
      enddo
c
      do 100 ns=1,nsplit
c
       n0 = 1
       n1 = 2
       n2 = n0
       m2 = m(ns) * 2
       dtau2 = dtau(ns) * 2.
c
c** below follows madala(1987)
cc     do 101 j=1,jlx
cc     do 101 i=1,ilx
c                            deld, delh: 1,ilx on cross grid
cc      deld(i,j,ns,n1) = deld(i,j,ns,n0)
cc      delh(i,j,ns,n1) = delh(i,j,ns,n0)
cc101  continue
c
#ifdef MPP1
       do 1 j=1,jendx
#else
       do 1 j=1,jlx
#endif
       do 1 i=1,ilx
c                            deld, delh: 1,ilx on cross grid
        ddsum(i,j,ns) = deld(i,j,ns,n0)
        dhsum(i,j,ns) = delh(i,j,ns,n0)
1      continue
c
c** first step, use forward scheme
C=======================================================================
c
c** compute gradient of delh;  output = (work1,work2)
c
#ifdef MPP1
      call MPI_SENDRECV(delh(1,jxp,ns,n0),ix,MPI_DOUBLE_PRECISION,
     &          ieast,1,delh(1,0  ,ns,n0),ix,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          fac = dx2 * msfx(i,j)
          work(i,j,1) = ( delh(i  ,j,ns,n0) + delh(i-1,j  ,ns,n0)
     &                  - delh(i,j-1,ns,n0) - delh(i-1,j-1,ns,n0) )/fac
          work(i,j,2) = ( delh(i  ,j,ns,n0) + delh(i  ,j-1,ns,n0)
     &                  - delh(i-1,j,ns,n0) - delh(i-1,j-1,ns,n0) )/fac
        enddo
      enddo
c
C=======================================================================
       do nw=1,2
#ifdef MPP1
         do j=jbegin,jendx
#else
         do j=2,jlx
#endif
           do i=2,ilx
c                            work: 2,ilx on dot grid
             work(i,j,nw) = work(i,j,nw) * psdot(i,j)
           enddo
         enddo
       enddo
C=======================================================================
c
c** compute divergence z from u and v
c    ( u must be pstar * u ; similarly for v )
c    ( note: map scale factors have been inverted in model (init) )
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          uu(i,j) = work(i,j,1) * msfd(i,j)
          vv(i,j) = work(i,j,2) * msfd(i,j)
        enddo
      enddo
c
#ifdef MPP1
      do i=1,ix
        wksend(i)   =uu(i,1)
        wksend(i+ix)=vv(i,1)
      enddo
      call MPI_SENDRECV(wksend(1),2*ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                  wkrecv(1),2*ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                  MPI_COMM_WORLD,status,ierr)
      do i=1,ix
        uu(i,jxp+1) = wkrecv(i)
        vv(i,jxp+1) = wkrecv(i+ix)
      enddo
#endif
c
#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif
        do i=2,ilxm
          fac = dx2 * msfx(i,j) * msfx(i,j)
          work(i,j,3) = ( -uu(i+1,j)+uu(i+1,j+1)-uu(i,j)+uu(i,j+1)
     &            +vv(i+1,j)+vv(i+1,j+1)-vv(i,j)-vv(i,j+1) ) / fac
        enddo
      enddo
c
C=======================================================================
c
#ifdef MPP1
       do j=jbegin,jendm
#else
       do j=2,jlxm
#endif
         do i=2,ilxm
c                            work3: 2,ilxm on cross grid
           deld(i,j,ns,n1) = deld(i,j,ns,n0)
     a           - dtau(ns) * work(i,j,3) + deld(i,j,ns,3) / m2
           delh(i,j,ns,n1) = delh(i,j,ns,n0)
     a           - dtau(ns) * hbar(ns) * deld(i,j,ns,n0) / psa(i,j)
     b           + delh(i,j,ns,3) / m2
         enddo
       enddo

c** not in madala(1987)
       fac = ( m(ns) - 1. ) / m(ns)
       do i=2,ilxm
#ifdef MPP1
         if (myid .eq. 0) delh(i,  1,ns,n1) = delh(i,  1,ns,n0) * fac
#else
         delh(i,  1,ns,n1) = delh(i,  1,ns,n0) * fac
#endif
#ifdef MPP1
         if (myid .eq. nproc-1)
     &   delh(i,jendx,ns,n1) = delh(i,jendx,ns,n0) * fac
#else
         delh(i,jlx,ns,n1) = delh(i,jlx,ns,n0) * fac
#endif
       enddo
#ifdef MPP1
       do j=1,jendx
#else
       do j=1,jlx
#endif
         delh(  1,j,ns,n1) = delh(  1,j,ns,n0) * fac
         delh(ilx,j,ns,n1) = delh(ilx,j,ns,n0) * fac
       enddo
c
#ifdef MPP1
       do j=1,jendx
#else
       do j=1,jlx
#endif
         do i=1,ilx
           ddsum(i,j,ns) = ddsum(i,j,ns) + deld(i,j,ns,n1)
           dhsum(i,j,ns) = dhsum(i,j,ns) + delh(i,j,ns,n1)
         enddo
       enddo
c
c** subsequent steps, use leapfrog scheme
       do 50 n=2,m2
C=======================================================================
c
c** compute gradient of delh;  output = (work1,work2)
c
#ifdef MPP1
      call MPI_SENDRECV(delh(1,jxp,ns,n1),ix,MPI_DOUBLE_PRECISION,
     &          ieast,1,delh(1,0  ,ns,n1),ix,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          fac = dx2 * msfx(i,j)
          work(i,j,1) = ( delh(i  ,j,ns,n1) + delh(i-1,j  ,ns,n1)
     &                  - delh(i,j-1,ns,n1) - delh(i-1,j-1,ns,n1) )/fac
          work(i,j,2) = ( delh(i  ,j,ns,n1) + delh(i  ,j-1,ns,n1)
     &                  - delh(i-1,j,ns,n1) - delh(i-1,j-1,ns,n1) )/fac
        enddo
      enddo
C=======================================================================
c
        do nw=1,2
#ifdef MPP1
          do j=jbegin,jendx
#else
          do j=2,jlx
#endif
            do i=2,ilx
              work(i,j,nw) = work(i,j,nw) * psdot(i,j)
            enddo
          enddo
        enddo
C=======================================================================
c
c** compute divergence z from u and v
c    ( u must be pstar * u ; similarly for v )
c    ( note: map scale factors have been inverted in model (init) )
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          uu(i,j) = work(i,j,1) * msfd(i,j)
          vv(i,j) = work(i,j,2) * msfd(i,j)
        enddo
      enddo
c
#ifdef MPP1
      do i=1,ix
        wksend(i)   =uu(i,1)
        wksend(i+ix)=vv(i,1)
      enddo
      call MPI_SENDRECV(wksend(1),2*ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                  wkrecv(1),2*ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                  MPI_COMM_WORLD,status,ierr)
      do i=1,ix
        uu(i,jxp+1) = wkrecv(i)
        vv(i,jxp+1) = wkrecv(i+ix)
      enddo
#endif
c
#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif
        do i=2,ilxm
          fac = dx2 * msfx(i,j) * msfx(i,j)
          work(i,j,3) = ( -uu(i+1,j)+uu(i+1,j+1)-uu(i,j)+uu(i,j+1)
     &            +vv(i+1,j)+vv(i+1,j+1)-vv(i,j)-vv(i,j+1) ) / fac
        enddo
      enddo
c
C=======================================================================
c
#ifdef MPP1
        do j=jbegin,jendm
#else
        do j=2,jlxm
#endif
          do i=2,ilxm
            deld(i,j,ns,n2) = deld(i,j,ns,n0)
     a        - dtau2 * work(i,j,3) + deld(i,j,ns,3) / m(ns)
            delh(i,j,ns,n2) = delh(i,j,ns,n0)
     a        - dtau2 * hbar(ns) * deld(i,j,ns,n1) / psa(i,j)
     b        + delh(i,j,ns,3) / m(ns)
          enddo
        enddo
c
c** not in madala(1987)
        do i=2,ilxm
#ifdef MPP1
          if (myid .eq. 0)
     &    delh(i,1,ns,n2)=2.*delh(i,1,ns,n1)-delh(i,1,ns,n0)
#else
          delh(i,1,ns,n2)=2.*delh(i,1,ns,n1)-delh(i,1,ns,n0)
#endif
#ifdef MPP1
          if (myid .eq. nproc-1)
     &    delh(i,jendx,ns,n2)=2.*delh(i,jendx,ns,n1)-delh(i,jendx,ns,n0)
#else
          delh(i,jlx,ns,n2)=2.*delh(i,jlx,ns,n1)-delh(i,jlx,ns,n0)
#endif
        enddo
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          delh(1,j,ns,n2)=2.*delh(1,j,ns,n1)-delh(1,j,ns,n0)
          delh(ilx,j,ns,n2)=2.*delh(ilx,j,ns,n1)-delh(ilx,j,ns,n0)
        enddo
c
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            ddsum(i,j,ns) = ddsum(i,j,ns) + deld(i,j,ns,n2)
            dhsum(i,j,ns) = dhsum(i,j,ns) + delh(i,j,ns,n2)
          enddo
        enddo
c
        n0 = n1
        n1 = n2
        n2 = n0
50     continue
c
100   continue
c
      return
      end
