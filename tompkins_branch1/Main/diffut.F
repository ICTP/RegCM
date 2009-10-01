ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffut_T(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     tb3d   : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable tb3d                            c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/slice.cb'
      include 'Common2/main.cb'
#else
      include 'Commons/slice.cb'
      include 'Commons/main.cb'
#endif
c
      real*8  ften(ix,kx),xkc(ix,kx)
      real*8  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*(tb3d(i,k,j+2)+tb3d(i,k,j-2)+
     &          tb3d(i+2,k,j)+tb3d(i-2,k,j)
     &                               - 4.*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j)) + 12.*tb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
   70  continue
c
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutQV(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     qvb3d  : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable qvb3d                           c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/slice.cb'
      include 'Common2/main.cb'
#else
      include 'Commons/slice.cb'
      include 'Commons/main.cb'
#endif
c
      real*8  ften(ix,kx),xkc(ix,kx)
      real*8  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (qvb3d(i,k,j+2)+qvb3d(i,k,j-2)+
     &       qvb3d(i+2,k,j)+qvb3d(i-2,k,j)
     & - 4.*(qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j)) + 12.*qvb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
   70  continue
c
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutQC(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     qcb3d  : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable qcb3d                           c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/slice.cb'
      include 'Common2/main.cb'
#else
      include 'Commons/slice.cb'
      include 'Commons/main.cb'
#endif
c
      real*8  ften(ix,kx),xkc(ix,kx)
      real*8  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (qcb3d(i,k,j+2)+qcb3d(i,k,j-2)+
     &       qcb3d(i+2,k,j)+qcb3d(i-2,k,j)
     & - 4.*(qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j)) + 12.*qcb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
   70  continue
c
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutCH(ften,xkc,c203,n,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     chib3d : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable chib3d                          c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/slice.cb'
      include 'Common2/main.cb'
#else
      include 'Commons/slice.cb'
      include 'Commons/main.cb'
#endif
c
      real*8  ften(ix,kx),xkc(ix,kx)
      real*8  c203
      integer n,i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (chib3d(i,k,j+2,n)+chib3d(i,k,j-2,n)+
     &       chib3d(i+2,k,j,n)+chib3d(i-2,k,j,n)
     & - 4.*(chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n))
     & + 12.*chib3d(i,k,j,n))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
   70  continue
c
      endif
c
      return
      end
